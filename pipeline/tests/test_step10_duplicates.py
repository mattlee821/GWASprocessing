"""README step 10, Duplicates: records with the same GRCh38 position and allele pair are duplicates; the first in input
order is kept, except that a record under a current rsID is kept over one under a retired rsID. Tables are standardised
in pieces (GWAS_TABLE_CHUNK) and in parallel (GWAS_WORKERS), so duplicates must be found across pieces and the output
must not depend on the number of pieces or workers.

Fixtures are seeded (821). The piece-equivalence case lifts GRCh37 positions with BCFtools/liftover and is skipped
without it.
"""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from testlib import (bcftools_liftover_available, make_reference, run_standardise, write_fasta,  # noqa: I001
                     write_grch38_references, write_table)

import gwas_alleles as alleles

BASES = np.array(list("ACGT"))
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
OFFSET = 100  # testlib offset chains: GRCh38 = GRCh37 + 100 on chromosomes 1 and 2
HEADER = ["CHR", "POS38", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]
GRCH38_CONTIG = "##contig=<ID=1,length=248956422,assembly=GRCh38>"
# af_concordance is pooled from moments (compared to 12 places); build detection samples the first piece only, so its hit
# counts differ with tiny pieces while the build decision must not.
VOLATILE_QC = {"started_at", "completed_at", "table_pieces", "settings_hash", "provenance", "af_concordance", "build_detection",
               "timings_s"}


def sequence(length: int, fixed: dict[int, str], seed: int = 821) -> str:
    rng = np.random.default_rng(seed)
    bases = list("".join(rng.choice(BASES, size=length)))
    for pos, text in fixed.items():
        bases[pos - 1:pos - 1 + len(text)] = list(text)
    return "".join(bases)


def write_vcf(path: Path, records: list[tuple]) -> Path:
    fmt = "ES:SE:LP:AF:SS"
    formats = [f'##FORMAT=<ID={key},Number=A,Type=Float,Description="{key}">' for key in fmt.split(":")]
    lines = ["##fileformat=VCFv4.2", GRCH38_CONTIG, *formats, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBMI"]
    lines += [f"{c}\t{p}\t{i}\t{r}\t{a}\t.\tPASS\t.\t{fmt}\t{v}" for c, p, i, r, a, v in records]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def synthetic_case(root: Path, n: int = 240, seed: int = 821) -> tuple[Path, Path]:
    """A GRCh37 table built to cross piece boundaries: duplicates far apart in input order (with effect 0.999, so the kept
    copy is identifiable), palindromic SNPs, rare variants, indels, alleles to swap or complement, rows missing required
    fields and rows that cannot be placed on GRCh38."""
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
            ea, oa, eaf = ref, alt, 1 - af
        elif kind == 2 and COMPLEMENT[ref] != alt:
            ea, oa = COMPLEMENT[alt], COMPLEMENT[ref]
        elif kind == 3:
            b = next(x for x in "ACGT" if x not in (ref, COMPLEMENT[ref]))
            ea, oa, eaf = b, COMPLEMENT[b], 0.25
            reference_rows[-1] = (chrom, f"rs{i + 1}", pos19, b, COMPLEMENT[b], 0.2)
        elif kind == 4:
            eaf = 0.004
        elif kind == 5:
            ea = alt + "T"
        elif kind == 6:
            eaf = None
        se = None if i % 37 == 0 else round(float(rng.uniform(0.01, 0.05)), 4)
        beta = round(float(rng.normal(0, 0.05)), 5)
        p = float(rng.uniform(1e-9, 1))
        table_rows.append([chrom, pos19, f"rs{i + 1}", ea, oa, eaf, beta, se, p, 5000])
    for j in (0, 7, 8, 14, 20):  # duplicates of early records, at the end of the input
        copy = list(table_rows[j])
        copy[6] = 0.999
        table_rows.append(copy)
    for k, (study, ref_af) in enumerate([(0.2, 0.25), (0.8, 0.22), (0.48, 0.5), (0.3, 0.3)]):
        pos19 = 5900 + k
        ref_base = genome["1"][pos19 + OFFSET - 1]
        partner = COMPLEMENT[ref_base]
        reference_rows.append(("1", f"rs{9000 + k}", pos19, partner, ref_base, ref_af))
        table_rows.insert(50 * k + 3, ["1", pos19, f"rs{9000 + k}", partner, ref_base, study, 0.01, 0.02, 0.3, 5000])
    table_rows.insert(10, ["", "", "rs99999999", "A", "G", 0.3, 0.01, 0.02, 0.3, 5000])
    table_rows.insert(100, ["", "", "chr_unknown", "A", "G", 0.003, 0.01, 0.02, 0.3, 5000])
    table_rows.insert(200, ["3", 2000, "rs77777777", "A", "G", 0.3, 0.01, 0.02, 0.3, 5000])
    reference = make_reference(root / "ref", lookup_rows=reference_rows, fasta=genome)
    write_fasta(reference / "genome", {c: seq[OFFSET:] for c, seq in genome.items()}, name="GRCh37.fa")
    raw = write_table(root / "gwas.tsv", ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"], table_rows)
    return raw, reference


class DuplicatesStepTest(unittest.TestCase):
    def run_case(self, tmp: Path, raw: Path, reference: Path, name: str, chunk: str = "1000000", workers: str = "1", **record):
        result = run_standardise(tmp, raw, reference, name=name,
                                 env={**os.environ, "GWAS_TABLE_CHUNK": chunk, "GWAS_WORKERS": workers}, **record)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result

    def assert_same(self, whole, pieces) -> None:
        body = lambda r: [line for line in r.text.splitlines() if not line.startswith("##fileDate")]  # noqa: E731
        self.assertEqual(body(pieces), body(whole))
        stable = lambda r: {k: v for k, v in r.qc.items() if k not in VOLATILE_QC}  # noqa: E731
        self.assertEqual(stable(pieces), stable(whole))
        for key in ("n", "frac_abs_diff_gt_0.2"):
            self.assertEqual(pieces.qc["af_concordance"][key], whole.qc["af_concordance"][key])
        if whole.qc["af_concordance"]["pearson_r"] is not None:
            self.assertAlmostEqual(pieces.qc["af_concordance"]["pearson_r"], whole.qc["af_concordance"]["pearson_r"], places=12)
        decision = lambda r: {k: (r.qc.get("build_detection") or {}).get(k) for k in ("build", "method", "declared")}  # noqa: E731
        self.assertEqual(decision(pieces), decision(whole))

    def test_duplicates(self):
        tmp = Path(tempfile.mkdtemp())
        genome = sequence(4000, {1100: "G", 1200: "G", 1300: "G", 1400: "G", 1500: "G", 3001: "GCACACATT"})
        reference = make_reference(tmp / "reference", lookup_rows=[], fasta={"1": genome})
        rows = [
            [1, 1100, "rsA", "A", "G", 0.3, 0.11, 0.02, 1e-6, 5000],
            [1, 1200, "rsB", "A", "G", 0.3, 0.12, 0.02, 1e-6, 5000],
            [1, 1300, "rsC1", "A", "G", 0.3, 0.13, 0.02, 1e-6, 5000],
            [1, 1300, "rsC2", "T", "G", 0.3, 0.14, 0.02, 1e-6, 5000],         # same position, another allele pair
            [1, 3005, "del_right", "A", "ACA", 0.3, 0.15, 0.02, 1e-6, 5000],  # CA deletion written at the right
            [1, 1400, "rsD", "A", "G", 0.3, 0.16, 0.02, 1e-6, 5000],
            [1, 1500, "rsA", "A", "G", 0.3, 0.17, 0.02, 1e-6, 5000],          # the same study ID at another position
            [1, 1100, "rsA_again", "A", "G", 0.3, 0.91, 0.02, 1e-6, 5000],    # repeat under another ID
            [1, 1200, "rsB_swapped", "G", "A", 0.7, -0.92, 0.02, 1e-6, 5000],  # the same variant coded the other way
            [1, 3003, "del_left", "A", "ACA", 0.3, 0.93, 0.02, 1e-6, 5000],   # the same deletion written further left
            [1, 1400, "rsD_reverse", "T", "C", 0.3, 0.94, 0.02, 1e-6, 5000],  # the same variant on the reverse strand
        ]
        raw = write_table(tmp / "dups.tsv", HEADER, rows)
        whole = self.run_case(tmp, raw, reference, "whole")
        kept = {(r["POS"], r["REF"], r["ALT"]): (r["INFO"]["SNPID"], r["FORMAT"]["ES"]) for r in whole.rows}

        with self.subTest("records with the same GRCh38 position and allele pair are duplicates, never keyed on the ID"):
            self.assertEqual(kept[(1100, "G", "A")], ("rsA", "0.11"))
            self.assertEqual(kept[(1500, "G", "A")], ("rsA", "0.17"))
            self.assertEqual(whole.qc["duplicates_removed"], 4)
            self.assertEqual(len(whole.rows), len(rows) - 4)

        with self.subTest("the same position with another allele pair is not a duplicate"):
            self.assertEqual((kept[(1300, "G", "A")], kept[(1300, "G", "T")]), (("rsC1", "0.13"), ("rsC2", "0.14")))

        with self.subTest("the same variant coded both ways, or on the other strand, is kept once"):
            self.assertEqual(kept[(1200, "G", "A")], ("rsB", "0.12"))
            self.assertEqual(kept[(1400, "G", "A")], ("rsD", "0.16"))
            data = pd.DataFrame({"CHR": ["1"] * 3, "POS38": [5000] * 3, "EA": ["A", "G", "A"], "OA": ["G", "A", "C"],
                                 "BETA": [0.1, -0.1, 0.2]})
            unit, removed = alleles.deduplicate(data)
            self.assertEqual((removed, unit["BETA"].tolist()), (1, [0.1, 0.2]))

        with self.subTest("the same indel written 2 ways is recognised after left-alignment"):
            self.assertEqual(kept[(3001, "GCA", "G")], ("del_right", "0.15"))
            self.assertEqual([k for k in kept if k[0] in (3003, 3005)], [])

        with self.subTest("the first in input order is kept"):
            self.assertFalse({"0.91", "-0.92", "0.92", "0.93", "0.94"} & {es for _, es in kept.values()})

        with self.subTest("the step log and the ##SAMPLE header count the duplicates"):
            step = next(s for s in whole.qc["steps"] if s["title"] == "Duplicates")
            self.assertEqual((step["removed"], step["removed_by"]), (4, {"duplicate of an earlier record": 4}))
            self.assertEqual(whole.sample["VariantsNotHarmonised"], "4")

        with self.subTest("duplicates spanning a piece boundary are removed, with the same output for any pieces and workers"):
            serial = self.run_case(tmp, raw, reference, "serial", chunk="2")
            self.assertEqual(serial.qc["table_pieces"]["pieces"], 6)
            self.assert_same(whole, serial)
            parallel = self.run_case(tmp, raw, reference, "parallel", chunk="3", workers="3")
            self.assertEqual((parallel.qc["table_pieces"]["pieces"], parallel.qc["table_pieces"]["workers"]), (4, 3))
            self.assert_same(whole, parallel)

        with self.subTest("duplicates spanning a piece boundary are removed from a GWAS-VCF"):
            # with 2 records per piece, the swapped-allele duplicate at 2000 starts the second piece
            vcf_reference = make_reference(tmp / "vcf_reference", lookup_rows=[],
                                           fasta={"1": sequence(3500, {1000: "A", 2000: "C", 3000: "G"})})
            vcf = write_vcf(tmp / "chunks.gwas.vcf", [
                ("1", 1000, "rs1", "A", "G", "0.1:0.02:6:0.3:5000"), ("1", 2000, "rs2", "C", "T", "0.1:0.02:6:0.3:5000"),
                ("1", 2000, "rs2", "T", "C", "-0.1:0.02:6:0.7:5000"), ("1", 3000, "rs3", "G", "A", "0.1:0.02:6:0.3:5000")])
            result = self.run_case(tmp, vcf, vcf_reference, "vcf", chunk="2")
            self.assertEqual(result.qc["table_pieces"]["pieces"], 2)
            self.assertEqual([(r["POS"], r["REF"], r["ALT"], r["FORMAT"]["ES"]) for r in result.rows],
                             [(1000, "A", "G", "0.1"), (2000, "C", "T", "0.1"), (3000, "G", "A", "0.1")])
            self.assertEqual(result.qc["duplicates_removed"], 1)

        with self.subTest("a record under a current rsID is kept over an earlier one under a retired rsID"):
            merged_reference = make_reference(tmp / "merged_reference", lookup_rows=[])
            write_grch38_references(merged_reference, [("1", 1000, 5004, "G", "A", 0.3), ("1", 2000, 5005, "C", "T", 0.3)],
                                    merges={9003: 5004})
            merged_rows = [["rs9003", "A", "G", 0.3, 0.5, 0.02, 1e-6, 4000],   # retired rsID, first in input order
                           ["rs5005", "T", "C", 0.3, 0.2, 0.02, 1e-6, 5000],
                           ["rs5004", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]]   # the same variant under its current rsID
            merged_raw = write_table(tmp / "merged.tsv", ["SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"], merged_rows)
            for name, chunk, workers in (("merged_whole", "1000000", "1"), ("merged_pieces", "1", "2")):
                merged = self.run_case(tmp, merged_raw, merged_reference, name, chunk=chunk, workers=workers)
                at_1000 = [r for r in merged.rows if r["POS"] == 1000]
                self.assertEqual([(r["ID"], r["INFO"]["SNPID"], r["FORMAT"]["ES"]) for r in at_1000], [("rs5004", "rs5004", "0.1")])
                self.assertEqual((merged.qc["duplicates_removed"], merged.qc["reference_fills"]["POS38_merged_rsid"]), (1, 1))

        with self.subTest("a seeded GRCh37 table gives the same output whatever the number of pieces and workers"):
            if not bcftools_liftover_available():
                self.skipTest("bcftools with the +liftover plugin is required")
            big_raw, big_reference = synthetic_case(tmp / "synthetic")
            big = self.run_case(tmp, big_raw, big_reference, "big_whole", input_build="hg19")
            self.assertEqual(big.qc["table_pieces"]["pieces"], 1)
            self.assertGreater(big.qc["duplicates_removed"], 0)
            self.assertGreater(big.qc["palindromic_resolution"]["palindromic_total"], 0)
            self.assertGreater(big.qc["reference_alignment"]["swapped"], 0)
            big_serial = self.run_case(tmp, big_raw, big_reference, "big_serial", chunk="7", input_build="hg19")
            self.assertEqual(big_serial.qc["table_pieces"]["pieces"], 36)  # 252 rows
            self.assert_same(big, big_serial)
            self.assertFalse({"0.999", "-0.999"} & {r["FORMAT"]["ES"] for r in big_serial.rows})  # first copies kept
            parallel = self.run_case(tmp, big_raw, big_reference, "big_parallel", chunk="5", workers="3", input_build="hg19")
            self.assertEqual(parallel.qc["table_pieces"]["workers"], 3)
            self.assert_same(big, parallel)


if __name__ == "__main__":
    unittest.main()
