"""README step 5, GRCh38 coordinates: rsID-only tables take GRCh38 positions from dbSNP (retired rsIDs through the
merge map, rsIDs at more than 1 position left out), and GRCh37 positions are lifted with BCFtools/liftover, following the
effect allele through reference changes and minus-strand chain blocks and counting the rows that cannot be lifted by
reason.

The liftover fixture is seeded (821): GRCh37 chromosomes 1 and 2 of 4,000 bases, GRCh38 of 4,100, and a UCSC chain with
a plus-strand block (GRCh37 chr1 [0, 1000) to GRCh38 chr1 [100, 1100)), a minus-strand block (GRCh37 chr1 [1000, 2000) to
the reverse strand of GRCh38 chr1, so GRCh37 position P is GRCh38 position 4101 - P), a block to another chromosome
(GRCh37 chr1 [2000, 3000) to GRCh38 chr2 [3100, 4100)) and chromosome 2 shifted by 100; GRCh37 chr1 [3000, 4000) has no
chain.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from testlib import (SRC, chain_block, example_reference, make_reference, needs_liftover, run_standardise, write_chain,  # noqa: I001
                     write_fasta, write_grch38_references, write_table)

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
HEADER = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P"]
RS_HEADER = ["SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]
GRCH38_CONTIG = "##contig=<ID=1,length=248956422,assembly=GRCh38>"
L37, L38 = 4000, 4100


def other_base(ref: str) -> str:
    """An allele that is neither the base nor its complement, so the pair is not palindromic."""
    return next(b for b in "ACGT" if b not in (ref, COMPLEMENT[ref]))


def write_vcf(path: Path, records: list[tuple], contig: str) -> Path:
    """A minimal GWAS-VCF (ES:SE:LP:AF:SS) with 1 contig line."""
    lines = ["##fileformat=VCFv4.2", contig,
             '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">',
             '##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error">',
             '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value">',
             '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">',
             '##FORMAT=<ID=SS,Number=A,Type=Float,Description="Sample size">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBMI"]
    lines += [f"{c}\t{p}\t{i}\t{r}\t{a}\t.\tPASS\t.\tES:SE:LP:AF:SS\t{v}" for c, p, i, r, a, v in records]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def liftover_fixture(root: Path):
    """Reference root with the chains and genomes described in the module docstring; returns (root, GRCh37, GRCh38)."""
    rng = np.random.default_rng(821)
    g37 = {c: list("".join(rng.choice(list("ACGT"), size=L37))) for c in ("1", "2")}
    g37["2"][2001 - 1] = "C"  # the base before the indel is not the inserted base, so left-alignment leaves it in place
    g38 = {"1": list("".join(rng.choice(list("ACGT"), size=L38))), "2": list("".join(rng.choice(list("ACGT"), size=L38)))}
    g38["1"][100:1100] = g37["1"][0:1000]
    for p0 in range(1000, 2000):
        g38["1"][L38 - 1 - p0] = COMPLEMENT[g37["1"][p0]]
    g38["2"][100:3100] = g37["2"][0:3000]
    g38["2"][3100:4100] = g37["1"][2000:3000]
    reference = make_reference(root, lookup_rows=[], offset_chains=False)
    write_chain(reference / "genome" / "hg19ToHg38.over.chain.gz", [
        chain_block("chr1", "chr1", L37, L38, 0, 1000, 100, 1),
        chain_block("chr1", "chr1", L37, L38, 1000, 1000, 1000, 2, dst_strand="-"),
        chain_block("chr1", "chr2", L37, L38, 2000, 1000, 3100, 3),
        chain_block("chr2", "chr2", L37, L38, 0, 3000, 100, 4)])
    return reference, g37, g38


def write_genomes(reference: Path, g37: dict, g38: dict) -> None:
    write_fasta(reference / "genome", {c: "".join(s) for c, s in g37.items()}, name="GRCh37.fa")
    write_fasta(reference / "genome", {c: "".join(s) for c, s in g38.items()}, name="GRCh38.fa")


def build_dbsnp_index(reference: Path, records: list[tuple]) -> None:
    """Build dbsnp/GRCh38 with src/reference_data/build_dbsnp_index.py from a tiny dbSNP-style VCF (RefSeq contigs)."""
    work = reference / "dbsnp_source"
    work.mkdir(parents=True, exist_ok=True)
    report = work / "assembly_report.txt"
    report.write_text("# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\t"
                      "Relationship\tRefSeq-Accn\n1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\tNC_000001.11\n"
                      "2\tassembled-molecule\t2\tChromosome\tCM000664.2\t=\tNC_000002.12\n", encoding="utf-8")
    lines = ["##fileformat=VCFv4.2", "##dbSNP_BUILD_ID=157", "##contig=<ID=NC_000001.11>", "##contig=<ID=NC_000002.12>",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    accession = {"1": "NC_000001.11", "2": "NC_000002.12"}
    lines += [f"{accession[c]}\t{p}\trs{rs}\t{r}\t{a}\t.\t.\t." for c, p, rs, r, a in sorted(records, key=lambda x: (x[0], x[1]))]
    vcf = work / "dbsnp.vcf"
    vcf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    env = dict(os.environ)
    tool = os.environ.get("GWAS_BCFTOOLS") or shutil.which("bcftools")
    if tool:
        env["PATH"] = str(Path(tool).parent) + os.pathsep + env.get("PATH", "")
    proc = subprocess.run([sys.executable, str(SRC / "reference_data" / "build_dbsnp_index.py"), "--dbsnp", str(vcf),
                           "--assembly-report", str(report), "--out", str(reference / "dbsnp" / "GRCh38"), "--force"],
                          capture_output=True, text=True, env=env, timeout=300)
    if proc.returncode != 0:
        raise AssertionError(f"build_dbsnp_index.py exited {proc.returncode}\n{proc.stderr}")


@needs_liftover
class Step05GRCh38CoordinatesTest(unittest.TestCase):
    def test_step05_grch38_coordinates(self):
        tmp = Path(tempfile.mkdtemp())

        # --- rsID-only table: positions from dbSNP, retired rsIDs through the merge map -------------------------------
        variants = [("1", 1000 + 33 * i, 5000 + i, "G", "A", 0.3) for i in range(5)]
        merged_ref = write_grch38_references(tmp / "merged_reference", variants,
                                             merges={9000: 5002, 9002: 999998, 9003: 5004})
        rows = [["rs9003", "A", "G", 0.3, 0.5, 0.02, 1e-6, 4000],  # retired rsID of rs5004, reported first
                ["rs5000", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],
                ["rs9000", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],  # retired rsID of rs5002
                ["rs9002", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],  # retired into an rsID that is not in the index
                ["rs999999", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],  # not in dbSNP
                ["rs5004", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]]
        result = run_standardise(tmp, write_table(tmp / "rsid.tsv", RS_HEADER, rows), merged_ref, ancestry="EUR", name="rsid")
        self.assertEqual(result.returncode, 0, result.stderr)
        by_id = result.by_id()

        with self.subTest("files with rsIDs and no positions take GRCh38 positions from dbSNP, without liftover"):
            self.assertEqual((by_id["rs5000"]["CHROM"], by_id["rs5000"]["POS"]), ("1", 1000))
            self.assertEqual(result.qc["coordinate_source"], "dbsnp")
            self.assertNotIn("hg19_to_hg38", result.qc.get("liftover", {}))
            self.assertEqual(result.metadata["coordinate_build"], "GRCh38")
            self.assertEqual(result.qc["reference_fills"]["POS38"], 4)

        with self.subTest("an rsID absent from dbSNP gets no position and is counted"):
            self.assertEqual(result.qc["incomplete_by_cause"]["rsid_not_in_coordinate_index"]["rows"], 2)

        with self.subTest("a retired rsID takes the position of the current rsID from the merge map"):
            moved = by_id["rs5002"]
            self.assertEqual((moved["CHROM"], moved["POS"], moved["FORMAT"]["ID"]), ("1", 1066, "rs5002"))
            self.assertEqual(result.qc["reference_fills"]["POS38_merged_rsid"], 2)

        with self.subTest("INFO/SNPID keeps the study's rsID"):
            self.assertEqual(by_id["rs5002"]["INFO"]["SNPID"], "rs9000")

        with self.subTest("a row under the current rsID is kept over one under a retired rsID"):
            self.assertEqual((by_id["rs5004"]["INFO"]["SNPID"], by_id["rs5004"]["FORMAT"]["ES"]), ("rs5004", "0.1"))
            self.assertEqual(result.qc["duplicates_removed"], 1)
            self.assertEqual(len(result.rows), 3)

        with self.subTest("rsIDs at more than 1 position are not used"):
            multi_ref = write_grch38_references(tmp / "multi_reference", [])
            build_dbsnp_index(multi_ref, [("1", 1500, 7000, "G", "A"), ("2", 2500, 7000, "G", "A"), ("1", 1600, 7001, "C", "T")])
            rows = [["rs7000", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000], ["rs7001", "T", "C", 0.3, 0.1, 0.02, 1e-6, 5000]]
            multi = run_standardise(tmp, write_table(tmp / "multi.tsv", RS_HEADER, rows), multi_ref, ancestry="EUR", name="multi")
            self.assertEqual(multi.returncode, 0, multi.stderr)
            self.assertEqual([(r["INFO"]["SNPID"], r["CHROM"], r["POS"]) for r in multi.rows], [("rs7001", "1", 1600)])
            self.assertEqual(multi.qc["incomplete_by_cause"]["rsid_not_in_coordinate_index"]["rows"], 1)

        # --- GRCh37 table lifted with BCFtools/liftover ------------------------------------------------------------
        reference, g37, g38 = liftover_fixture(tmp / "liftover_reference")
        b37 = lambda c, p: g37[c][p - 1]  # noqa: E731
        sites = {}
        r = b37("1", 201); sites["plus"] = (r, other_base(r))
        r = b37("1", 301); sites["swap"] = (r, other_base(r))
        g38["1"][401 - 1] = sites["swap"][1]  # GRCh38 carries the GRCh37 ALT base: the reference allele changed
        r = b37("1", 701); sites["third"] = (r, other_base(r))
        g38["1"][801 - 1] = COMPLEMENT[sites["third"][1]]  # neither allele in GRCh38: a third allele at the site
        r = b37("1", 1501); sites["minus"] = (r, other_base(r))
        r = b37("1", 1601); sites["minus_rev"] = (r, other_base(r))
        r = b37("2", 1000); sites["plus_rev"] = (r, other_base(r))
        r = b37("1", 501); sites["mismatch"] = (r, other_base(r))
        r = b37("1", 3501); sites["unchained"] = (r, other_base(r))
        r = b37("1", 2501); sites["other_chrom"] = (r, other_base(r))
        write_genomes(reference, g37, g38)
        rows = [
            [1, 201, "plus", sites["plus"][1], sites["plus"][0], 0.3, 0.1, 0.02, 1e-6],
            [1, 301, "swap", sites["swap"][1], sites["swap"][0], 0.3, 0.1, 0.02, 1e-6],
            [1, 1501, "minus", sites["minus"][1], sites["minus"][0], 0.3, 0.1, 0.02, 1e-6],
            ["chr2", 1000, "plus_rev", COMPLEMENT[sites["plus_rev"][1]], COMPLEMENT[sites["plus_rev"][0]], 0.3, 0.1, 0.02, 1e-6],
            [1, 1601, "minus_rev", COMPLEMENT[sites["minus_rev"][1]], COMPLEMENT[sites["minus_rev"][0]], 0.3, 0.1, 0.02, 1e-6],
            [2, 2001, "indel", "CA", "C", 0.3, 0.1, 0.02, 1e-6],
            # alleles matching neither strand of GRCh37 (the 2 bases that are neither the reference base nor its complement)
            [1, 501, "mismatch", other_base(sites["mismatch"][0]), COMPLEMENT[other_base(sites["mismatch"][0])], 0.3, 0.1, 0.02, 1e-6],
            [1, 3501, "unchained", sites["unchained"][1], sites["unchained"][0], 0.3, 0.1, 0.02, 1e-6],
            [1, 2501, "other_chrom", sites["other_chrom"][1], sites["other_chrom"][0], 0.3, 0.1, 0.02, 1e-6],
            [1, 701, "third", sites["third"][1], sites["third"][0], 0.3, 0.1, 0.02, 1e-6],
            [1, 901, "chr1:901:D", "D", "I", 0.3, 0.1, 0.02, 1e-6],
        ]
        lifted = run_standardise(tmp, write_table(tmp / "b37.tsv", HEADER, rows), reference, input_build="hg19",
                                 ancestry="EUR", name="liftover")
        self.assertEqual(lifted.returncode, 0, lifted.stderr)
        out = lifted.by_snpid()
        lift = lifted.qc["liftover"]["hg19_to_hg38"]
        causes = lifted.qc["incomplete_by_cause"]

        with self.subTest("GRCh37 positions are lifted with BCFtools/liftover"):
            self.assertEqual((lift["tool"], lift["ok"]), ("bcftools +liftover", 6))
            plus = out["plus"]
            self.assertEqual((plus["CHROM"], plus["POS"], plus["REF"], plus["ALT"], plus["FORMAT"]["ES"], plus["INFO"]["HM"]),
                             ("1", 301, sites["plus"][0], sites["plus"][1], "0.1", "forward"))

        with self.subTest("the effect allele is followed through a reference change between builds"):
            swap = out["swap"]
            self.assertEqual(lift["swapped_between_builds"], 1)
            # the effect allele is the GRCh38 reference base, so it is written as REF and ES is negated relative to ALT
            self.assertEqual((swap["POS"], swap["REF"], swap["ALT"], swap["FORMAT"]["ES"]),
                             (401, sites["swap"][1], sites["swap"][0], "-0.1"))
            self.assertEqual(swap["INFO"]["HM"], "swapped")

        with self.subTest("a minus-strand chain block complements the alleles and HM records it (STRAND_CHANGED)"):
            minus = out["minus"]
            self.assertEqual((minus["POS"], minus["REF"], minus["ALT"], minus["FORMAT"]["ES"]),
                             (4101 - 1501, COMPLEMENT[sites["minus"][0]], COMPLEMENT[sites["minus"][1]], "0.1"))
            self.assertEqual(minus["INFO"]["HM"], "complemented")
            self.assertEqual(lift["minus_strand"], 2)

        with self.subTest("alleles are first oriented to the GRCh37 genome: a reverse-strand report reads complemented"):
            plus_rev = out["plus_rev"]
            self.assertEqual((plus_rev["CHROM"], plus_rev["POS"], plus_rev["REF"], plus_rev["ALT"], plus_rev["INFO"]["HM"]),
                             ("2", 1100, sites["plus_rev"][0], sites["plus_rev"][1], "complemented"))

        with self.subTest("a reverse-strand report lifted through a minus-strand block cancels to forward"):
            minus_rev = out["minus_rev"]
            self.assertEqual((minus_rev["POS"], minus_rev["REF"], minus_rev["ALT"], minus_rev["INFO"]["HM"]),
                             (4101 - 1601, COMPLEMENT[sites["minus_rev"][0]], COMPLEMENT[sites["minus_rev"][1]], "forward"))

        with self.subTest("chr-prefixed contigs are renamed"):
            self.assertEqual(out["plus_rev"]["CHROM"], "2")

        with self.subTest("indels are lifted"):
            self.assertEqual((out["indel"]["CHROM"], out["indel"]["POS"], out["indel"]["REF"], out["indel"]["ALT"]),
                             ("2", 2101, "C", "CA"))

        with self.subTest("rows that cannot be lifted are counted by reason"):
            self.assertEqual(lift["grch37_allele_mismatch"], 1)
            self.assertEqual(causes["alleles_do_not_match_grch37"]["rows"], 1)
            self.assertEqual(sum(lift["rejected"].values()), 1)  # the tool's reject reason, as reported
            self.assertNotIn("rejected", lift["rejected"])
            self.assertEqual(lift["other_chrom"], 1)
            self.assertEqual(lift["reference_allele_added"], 1)
            self.assertEqual(causes["liftover_to_grch38_failed"]["rows"], 3)
            self.assertFalse({"mismatch", "unchained", "other_chrom", "third"} & set(out))

        with self.subTest("I/D-coded indels with GRCh37 positions count as invalid alleles, not liftover failures"):
            self.assertEqual(causes["invalid_allele"]["rows"], 1)
            self.assertNotIn("chr1:901:D", out)
            self.assertEqual(sorted(out), ["indel", "minus", "minus_rev", "plus", "plus_rev", "swap"])

        # --- GWAS-VCF inputs ----------------------------------------------------------------------------------------
        with self.subTest("a GWAS-VCF with chr-prefixed GRCh37 contigs is renamed and lifted"):
            raw = write_vcf(tmp / "chr37.gwas.vcf", [("chr1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")],
                            contig="##contig=<ID=chr1,length=249250621,assembly=GRCh37>")
            vcf = run_standardise(tmp, raw, example_reference(genomes=True), name="vcf37")
            self.assertEqual(vcf.returncode, 0, vcf.stderr)
            row = vcf.rows[0]
            self.assertEqual((row["CHROM"], row["POS"], row["REF"], row["ALT"], row["INFO"]["HM"]), ("1", 1100, "G", "A", "forward"))
            self.assertEqual(vcf.qc["liftover"]["hg19_to_hg38"]["ok"], 1)
            self.assertFalse(any("GRCh37" in h for h in vcf.header if h.startswith("##contig")))

        with self.subTest("a GWAS-VCF with chr-prefixed GRCh38 contigs is renamed and not lifted"):
            raw = write_vcf(tmp / "chr38.gwas.vcf", [("chr1", 2000, "rs2", "C", "T", "0.1:0.02:6:0.3:5000")],
                            contig="##contig=<ID=chr1,length=248956422,assembly=GRCh38>")
            vcf = run_standardise(tmp, raw, example_reference(), name="vcf38")
            self.assertEqual(vcf.returncode, 0, vcf.stderr)
            self.assertEqual((vcf.rows[0]["CHROM"], vcf.rows[0]["POS"]), ("1", 2000))
            self.assertNotIn("hg19_to_hg38", vcf.qc.get("liftover", {}))

        with self.subTest("no liftable records is a data error"):
            raw = write_vcf(tmp / "far.gwas.vcf", [("1", 5_000_000, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")],
                            contig="##contig=<ID=1,length=249250621,assembly=GRCh37>")
            far = run_standardise(tmp, raw, example_reference(genomes=True), name="far")
            self.assertEqual(far.returncode, 65)
            self.assertIn("no variants", far.stderr.lower())

        with self.subTest("GRCh38 chr:pos identifiers never borrow GRCh37 annotations"):
            # 30 GRCh38 variants, plus 3 whose GRCh38 position coincides with the GRCh37 position of an unrelated variant
            true_rows = [("1", f"rs{9000 + i}", 9000 + 1000 * i, "A", "G", 0.3) for i in range(30)]
            decoys = [("1", f"rs{7000 + j}", 200000 + 1000 * j, "A", "G", 0.3) for j in range(3)]
            ref38 = make_reference(tmp / "decoy_reference", lookup_rows=true_rows + decoys)
            positions = [9100 + 1000 * i for i in range(30)] + [200000 + 1000 * j for j in range(3)]
            raw = write_table(tmp / "chrpos.tsv", ["CHR", "POS", "MarkerName", "EA", "OA", "EAF", "BETA", "SE", "P", "N"],
                              [[1, p, f"1:{p}", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000] for p in positions])
            decoy = run_standardise(tmp, raw, ref38, input_build="hg38", name="decoy")
            self.assertEqual(decoy.returncode, 0, decoy.stderr)
            by_pos = {r["POS"]: r for r in decoy.rows}
            self.assertEqual((by_pos[9100]["ID"], by_pos[200000]["ID"]), ("rs9000", "."))
            self.assertNotIn("hg19_to_hg38", decoy.qc.get("liftover", {}))


if __name__ == "__main__":
    unittest.main()
