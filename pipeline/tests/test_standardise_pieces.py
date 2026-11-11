"""Tables standardised in pieces (GWAS_TABLE_CHUNK) and in parallel (GWAS_WORKERS) must give the same output as
the whole table: same VCF records in the same order, the same QC counts and the same step log. Also covers the
step log's reconciliation (every removed row accounted for once) and the prepared-lookup cache.

The synthetic table is seeded (821) and built to cross piece boundaries: duplicate records far apart in input
order (with different effects, so the kept copy is identifiable), palindromic SNPs, rare variants, indels,
alleles that need swapping, complementing or dropping against the FASTA, rows missing required fields, and rows
that cannot be placed on GRCh38 (an rsID absent from the reference, no position or rsID, a failed liftover).
"""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from testlib import chain_block, make_reference, run_standardise, write_chain, write_lookup, write_table  # noqa: I001 (puts pipeline/bin on sys.path)

import gwas_alleles as alleles
from test_standardise_passthrough import write_vcf

BASES = np.array(list("ACGT"))
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
OFFSET = 100  # testlib offset chains: GRCh38 = GRCh37 + 100 on chromosomes 1 and 2
# af_concordance is pooled from moments (compared to 12 places below); build_detection samples the first piece only,
# so its hit counts differ with tiny pieces while the build decision must not.
VOLATILE_QC = {"started_at", "completed_at", "table_pieces", "settings_hash", "provenance", "af_concordance", "build_detection",
               "timings_s"}


def synthetic_case(root: Path, n: int = 240, seed: int = 821) -> tuple[Path, Path]:
    rng = np.random.default_rng(seed)
    genome = {c: "".join(rng.choice(BASES, size=7000)) for c in ("1", "2")}
    reference_rows, table_rows = [], []
    for i in range(n):
        chrom = "1" if i % 2 else "2"
        pos19 = int(rng.integers(1000, 5800))
        ref = genome[chrom][pos19 + OFFSET - 1]
        alt = str(rng.choice([b for b in "ACGT" if b != ref]))
        af = float(np.round(rng.uniform(0.02, 0.98), 3))
        reference_rows.append((chrom, f"rs{i + 1}", pos19, alt, ref, af))
        ea, oa, eaf = alt, ref, af
        kind = i % 10
        if kind == 1:
            ea, oa = ref, alt  # effect allele is the reference base: swapped
            eaf = 1 - af
        elif kind == 2 and COMPLEMENT[ref] != alt:
            ea, oa = COMPLEMENT[alt], COMPLEMENT[ref]  # reverse-strand report: complemented
        elif kind == 3:
            # A palindromic pair of the 2 bases that are neither the reference base nor its complement: its strand
            # resolves against the lookup, then it cannot match the GRCh38 base, so alignment drops it.
            b = next(x for x in "ACGT" if x not in (ref, COMPLEMENT[ref]))
            ea, oa, eaf = b, COMPLEMENT[b], 0.25
            reference_rows[-1] = (chrom, f"rs{i + 1}", pos19, b, COMPLEMENT[b], 0.2)
        elif kind == 4:
            eaf = 0.004  # rare
        elif kind == 5:
            ea = alt + "T"  # indel
        elif kind == 6:
            eaf = None  # filled from the reference
        se = None if i % 37 == 0 else round(float(rng.uniform(0.01, 0.05)), 4)
        beta = round(float(rng.normal(0, 0.05)), 5)
        p = float(rng.uniform(1e-9, 1))
        table_rows.append([chrom, pos19, f"rs{i + 1}", ea, oa, eaf, beta, se, p, 5000])
    # Duplicates of early records, appended at the end with different effects: the first copy must win.
    for j in (0, 7, 8, 14, 20):
        copy = list(table_rows[j])
        copy[6] = 0.999
        table_rows.append(copy)
    # Palindromic SNPs at their own positions, resolvable and ambiguous.
    for k, (study, ref_af) in enumerate([(0.2, 0.25), (0.8, 0.22), (0.48, 0.5), (0.3, 0.3)]):
        pos19 = 5900 + k
        ref_base = genome["1"][pos19 + OFFSET - 1]
        partner = COMPLEMENT[ref_base]
        reference_rows.append(("1", f"rs{9000 + k}", pos19, partner, ref_base, ref_af))
        table_rows.insert(50 * k + 3, ["1", pos19, f"rs{9000 + k}", partner, ref_base, study, 0.01, 0.02, 0.3, 5000])
    # Rows that cannot be placed on GRCh38: an rsID absent from the reference, no position or rsID, no liftover chain.
    table_rows.insert(10, ["", "", "rs99999999", "A", "G", 0.3, 0.01, 0.02, 0.3, 5000])
    table_rows.insert(100, ["", "", "chr_unknown", "A", "G", 0.003, 0.01, 0.02, 0.3, 5000])
    table_rows.insert(200, ["3", 2000, "rs77777777", "A", "G", 0.3, 0.01, 0.02, 0.3, 5000])
    reference = make_reference(root / "ref", lookup_rows=reference_rows, fasta=genome)
    raw = write_table(root / "gwas.tsv", ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"], table_rows)
    return raw, reference


class PiecewiseEquivalenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.raw, self.reference = synthetic_case(self.tmp)

    def run_case(self, name: str, **env: str):
        defaults = {"GWAS_TABLE_CHUNK": "1000000", "GWAS_WORKERS": "1"}  # whole table unless the test overrides
        result = run_standardise(self.tmp, self.raw, self.reference, input_build="hg19",
                                 name=name, env={**os.environ, **defaults, **env})
        self.assertEqual(result.returncode, 0, result.stderr)
        return result

    def assert_same(self, whole, pieces):
        body = lambda r: [line for line in r.text.splitlines() if not line.startswith("##fileDate")]  # noqa: E731
        self.assertEqual(body(pieces), body(whole))
        stable = lambda r: {k: v for k, v in r.qc.items() if k not in VOLATILE_QC}  # noqa: E731
        self.assertEqual(stable(pieces), stable(whole))
        for key in ("n", "frac_abs_diff_gt_0.2"):
            self.assertEqual(pieces.qc["af_concordance"][key], whole.qc["af_concordance"][key])
        self.assertAlmostEqual(pieces.qc["af_concordance"]["pearson_r"], whole.qc["af_concordance"]["pearson_r"], places=12)
        decision = lambda r: {k: r.qc["build_detection"][k] for k in ("build", "method", "declared")}  # noqa: E731
        self.assertEqual(decision(pieces), decision(whole))

    def test_pieces_and_workers_match_the_whole_table(self):
        whole = self.run_case("whole")
        self.assertEqual(whole.qc["table_pieces"]["pieces"], 1)
        self.assertGreater(whole.qc["duplicates_removed"], 0)
        self.assertGreater(whole.qc["palindromic_resolution"]["palindromic_total"], 0)
        self.assertGreater(whole.qc["reference_alignment"]["dropped_mismatch"], 0)
        self.assertGreater(whole.qc["reference_alignment"]["swapped"], 0)
        serial = self.run_case("serial", GWAS_TABLE_CHUNK="7")
        self.assertEqual(serial.qc["table_pieces"]["pieces"], 36)  # 252 rows
        self.assert_same(whole, serial)
        parallel = self.run_case("parallel", GWAS_TABLE_CHUNK="7", GWAS_WORKERS="3")
        self.assertEqual(parallel.qc["table_pieces"]["workers"], 3)
        self.assert_same(whole, parallel)
        self.assert_same(whole, self.run_case("pieces5", GWAS_TABLE_CHUNK="5", GWAS_WORKERS="2"))

    def test_the_first_copy_of_a_duplicate_is_kept_across_pieces(self):
        result = self.run_case("dups", GWAS_TABLE_CHUNK="4")
        self.assertNotIn("0.999", [row["FORMAT"]["ES"] for row in result.rows])
        self.assertNotIn("-0.999", [row["FORMAT"]["ES"] for row in result.rows])


def assert_reconciles(case: unittest.TestCase, qc: dict) -> None:
    steps = qc["steps"]
    case.assertEqual(steps[0]["rows_in"], qc["input_rows"])
    case.assertEqual(steps[-1]["rows_out"], qc["output_rows"])
    for before, after in zip(steps, steps[1:]):
        case.assertEqual(after["rows_in"], before["rows_out"])
    for step in steps:
        case.assertEqual(step["rows_in"] - step["removed"], step["rows_out"])
        case.assertEqual(sum(step["removed_by"].values()), step["removed"])
    case.assertNotIn("does not reconcile", " ".join(qc["warnings"]))


class StepLogTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.raw, self.reference = synthetic_case(self.tmp)

    def run_case(self, name: str):
        result = run_standardise(self.tmp, self.raw, self.reference, input_build="hg19", name=name,
                                 env={**os.environ, "GWAS_TABLE_CHUNK": "1000000", "GWAS_WORKERS": "1"})
        self.assertEqual(result.returncode, 0, result.stderr)
        return result

    def test_every_removed_row_is_counted_once_with_its_cause(self):
        result = self.run_case("steps")
        assert_reconciles(self, result.qc)
        causes = result.qc["incomplete_by_cause"]
        self.assertEqual(causes["rsid_not_in_coordinate_index"]["rows"], 1)
        self.assertEqual(causes["no_position_or_rsid"]["rows"], 1)
        self.assertEqual(causes["no_position_or_rsid"]["study_maf_below_threshold"], 1)
        self.assertEqual(causes["liftover_to_grch38_failed"]["rows"], 1)
        self.assertEqual(sum(c["rows"] for c in causes.values()), result.qc["rows_removed_missing_required"])
        self.assertIn("Alleles and GRCh38 coordinates", result.stderr)  # the summary is printed for the Nextflow task log
        self.assertIn(f"Written: {result.qc['output_rows']:,} of {result.qc['input_rows']:,} input rows", result.stderr)
        self.assertNotIn("UNRECONCILED", result.stderr)

    def test_sample_header_counts_add_up_and_match_the_step_log(self):
        result = self.run_case("sample")
        sample = {k: int(v) for k, v in result.sample.items() if k in ("TotalVariants", "VariantsNotRead", "HarmonisedVariants",
                                                                         "VariantsNotHarmonised", "SwitchedAlleles")}
        self.assertEqual(sample["TotalVariants"], sample["VariantsNotRead"] + sample["VariantsNotHarmonised"] + sample["HarmonisedVariants"])
        self.assertEqual((sample["TotalVariants"], sample["HarmonisedVariants"]), (result.qc["input_rows"], result.qc["output_rows"]))
        self.assertEqual(sample["VariantsNotHarmonised"],
                         result.qc["duplicates_removed"] + result.qc["reference_alignment"]["dropped_mismatch"])
        self.assertEqual(sample["SwitchedAlleles"], sum(r["INFO"]["HM"].endswith("swapped") for r in result.rows))

    def test_step_log_reconciles_for_a_lifted_gwas_vcf(self):
        reference = self.tmp / "vcf_reference"
        write_lookup(reference, [("1", "rs1", 1000, "G", "A", 0.3), ("1", "rs2", 2000, "C", "T", 0.004)])
        write_chain(reference / "genome" / "hg19ToHg38.over.chain.gz", [chain_block("chr1", "chr1", 249250621, 248956422, 0, 100000, 0, 1)])
        raw = write_vcf(self.tmp / "grch37.gwas.vcf", [
            ("1", 1000, "rs1", "A", "G", "0.01:0.01:1:0.3:1000"), ("1", 2000, "rs2", "T", "C", "0.01:0.01:1:0.004:1000"),
            ("1", 1000, "rs1", "A", "G", "0.02:0.01:1:0.3:1000"), ("1", 500000, "rs3", "A", "G", "0.01:0.01:1:0.3:1000"),
            ("GL000192.1", 10, "rs4", "A", "G", "0.01:0.01:1:0.3:1000"), ("2", 3000, "rs5", "A", "G", "0.01:0.01:1:0.3:1000")])
        result = run_standardise(self.tmp, raw, reference, name="vcf_steps", input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        assert_reconciles(self, result.qc)
        self.assertEqual(result.qc["input_rows"], 6)
        coordinates = next(s for s in result.qc["steps"] if s["title"] == "Alleles and GRCh38 coordinates")
        self.assertEqual(coordinates["removed"], 3)  # GL contig, and positions outside the chain on chromosomes 1 and 2
        duplicates = next(s for s in result.qc["steps"] if s["title"] == "Duplicates")
        self.assertEqual(duplicates["removed"], 1)


class PreparedLookupCacheTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.rows = [("1", "rs1", 1000, "A", "G", 0.3), ("2", "rs2", 2000, "C", "T", 0.6), ("X", ".", 3000, "G", "A", 0.1)]

    def test_cache_is_written_and_matches_a_fresh_read(self):
        write_lookup(self.root, self.rows, population="EUR")
        first, pop = alleles.read_reference_lookup(self.root, "EUR")
        self.assertEqual(pop, "EUR")
        cache = self.root / "1000genomes" / "phase3" / "lookup" / "EUR.prepared"
        self.assertTrue((cache / "meta.json").exists())
        cached, _ = alleles.read_reference_lookup(self.root, "EUR")
        pd.testing.assert_frame_equal(cached, first)
        pd.testing.assert_frame_equal(cached, alleles._read_and_prepare(cache.parent / "EUR.tsv"))

    def test_changed_lookup_file_is_read_again(self):
        path = write_lookup(self.root, self.rows, population="EUR")
        alleles.read_reference_lookup(self.root, "EUR")
        write_lookup(self.root, self.rows + [("3", "rs4", 4000, "T", "C", 0.4)], population="EUR")
        os.utime(path, ns=(path.stat().st_atime_ns, path.stat().st_mtime_ns + 10**9))
        table, _ = alleles.read_reference_lookup(self.root, "EUR")
        self.assertEqual(len(table), 4)
        self.assertIn(4, table["RS"].tolist())


if __name__ == "__main__":
    unittest.main()
