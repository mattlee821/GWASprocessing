"""End-to-end tests for GWAS-VCF inputs (e.g. OpenGWAS) in standardise_gwas.py. Records on GRCh37 need BCFtools/liftover
(needs_liftover) and the example reference with consistent GRCh37 and GRCh38 genomes."""

import os
import tempfile
import unittest
from pathlib import Path

from testlib import example_reference, make_reference, needs_liftover, run_standardise



def write_vcf(path, records, contig="##contig=<ID=1,length=249250621,assembly=HG19/GRCh37>", extra_header=(),
              fmt="ES:SE:LP:AF:SS", sample="BMI", formats=None):
    formats = formats or [
        '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">',
        '##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error">',
        '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">',
        '##FORMAT=<ID=SS,Number=A,Type=Float,Description="Sample size">',
    ]
    lines = ["##fileformat=VCFv4.2", *([contig] if contig else []), *extra_header, *formats,
             f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"]
    for chrom, pos, rsid, ref, alt, values, *info in records:
        lines.append(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t{info[0] if info else '.'}\t{fmt}\t{values}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class AdaptedPassthroughTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    @needs_liftover
    def test_existing_gwas_vcf_gets_manifest_sex_metadata_and_declared_build(self):
        raw = write_vcf(self.tmp / "existing_gwas.vcf", [(1, 1000, "rs1", "G", "A", "0.2:0.05:3:0.1:1000", "POS19=900;SNPID=1:1000")],
                        contig=None, extra_header=["##source=test"])
        undeclared = run_standardise(self.tmp, raw, example_reference(genomes=True), sex="male", ancestry="EUR", name="undeclared")
        self.assertEqual(undeclared.returncode, 78)  # 1 record is too few to detect a build, and none is declared
        result = run_standardise(self.tmp, raw, example_reference(genomes=True), sex="male", ancestry="EUR", input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.metadata["sex"], "male")
        self.assertEqual(result.sample["Sex"], "male")
        row = result.rows[0]
        self.assertEqual(row["POS"], 1100)
        self.assertEqual((result.qc["build_detection"]["build"], result.qc["build_detection"]["method"]), ("hg19", "declared"))

    def test_opengwas_lp_999_sentinel_is_recovered_from_the_z_score(self):
        raw = write_vcf(self.tmp / "p0.gwas.vcf", [("1", 1100, "rs1", "G", "A", "0.5:0.01:999:0.3:1000")], contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, example_reference(), input_build="hg38")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertAlmostEqual(float(result.rows[0]["FORMAT"]["LP"]), 544.2, delta=1.0)  # z = 50
        self.assertEqual(result.qc["p_value_source"]["from_z"], 1)

    @needs_liftover
    def test_passthrough_hg19_lifted_to_grch38(self):
        raw = write_vcf(self.tmp / "grch37.gwas.vcf", [("1", 1000, "rs1", "G", "A", "0.01:0.001:2:0.1:1000")])
        result = run_standardise(self.tmp, raw, example_reference(genomes=True))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.metadata["coordinate_build"], "GRCh38")
        self.assertEqual(result.qc["build_detection"]["declared_source"].split()[0].lower()[:3], "vcf")
        row = result.rows[0]
        self.assertEqual((row["POS"], row["ID"], row["INFO"]["HM"]), (1100, "rs1", "forward"))
        self.assertEqual(result.sample["Ancestry"], "ALL")
        self.assertEqual(result.sample["Sex"], "combined")

    def test_passthrough_output_is_position_sorted(self):
        raw = write_vcf(self.tmp / "unsorted.gwas.vcf", [("1", 5000, "rsB", "A", "G", "0.02:0.01:1:0.3:1000"),
                                                        ("1", 2000, "rsA", "C", "T", "0.01:0.01:1:0.3:1000"),
                                                        ("1", 9000, "rsC", "G", "A", "0.03:0.01:1:0.3:1000")],
                        contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, example_reference())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r["POS"] for r in result.rows], [2000, 5000, 9000])

    def test_passthrough_palindromic_strand_resolution(self):
        # GRCh37 900, 1900, 2900 are GRCh38 1000, 2000, 3000 in the fixture
        reference = make_reference(self.tmp / "reference", lookup_rows=[("1", "rsP1", 900, "A", "T", 0.20), ("1", "rsP2", 1900, "C", "G", 0.15),
                                                                        ("1", "rsP3", 2900, "A", "T", 0.47)])
        raw = write_vcf(self.tmp / "pass.gwas.vcf", [
            ("1", 1000, "rsP1", "T", "A", "0.1:0.01:1:0.20:1000"), ("1", 2000, "rsP2", "G", "C", "0.1:0.01:1:0.85:1000"),
            ("1", 3000, "rsP3", "T", "A", "0.1:0.01:1:0.47:1000"), ("1", 4000, "rsN1", "G", "A", "0.1:0.01:1:0.30:1000")],
            contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, reference, source_type="opengwas")
        self.assertEqual(result.returncode, 0, result.stderr)
        variants = {r["INFO"]["SNPID"]: (r["REF"], r["ALT"], r["INFO"]["HM"]) for r in result.rows}
        self.assertEqual(variants.get("rsP1")[:2], ("T", "A"))
        self.assertEqual(variants.get("rsP2")[:2], ("C", "G"))
        self.assertTrue(variants["rsP2"][2].startswith("palindrome_complemented"))
        self.assertTrue(variants["rsP3"][2].startswith("palindrome_ambiguous"))
        self.assertIn("rsN1", variants)
        detail = result.qc["palindromic_resolution"]
        self.assertEqual((detail["flipped"], detail["kept_same_strand"], detail["ambiguous_frequency_near_half"]), (1, 1, 1))

class PassthroughRegressionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    @needs_liftover
    def test_chr_prefixed_grch37_contigs_are_renamed_and_lifted(self):
        raw = write_vcf(self.tmp / "chr.gwas.vcf", [("chr1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")],
                        contig="##contig=<ID=chr1,length=249250621,assembly=GRCh37>")
        result = run_standardise(self.tmp, raw, example_reference(genomes=True))
        self.assertEqual(result.returncode, 0, result.stderr)
        row = result.rows[0]
        self.assertEqual((row["CHROM"], row["POS"]), ("1", 1100))

    def test_chr_prefixed_grch38_contigs_are_renamed(self):
        raw = write_vcf(self.tmp / "chr38.gwas.vcf", [("chr1", 2000, "rs2", "C", "T", "0.1:0.02:6:0.3:5000")],
                        contig="##contig=<ID=chr1,length=248956422,assembly=GRCh38>")
        row = run_standardise(self.tmp, raw, example_reference()).rows[0]
        self.assertEqual((row["CHROM"], row["POS"]), ("1", 2000))

    @needs_liftover
    def test_sample_size_is_added_from_opengwas_metadata(self):
        formats = ['##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect">', '##FORMAT=<ID=SE,Number=A,Type=Float,Description="SE">',
                   '##FORMAT=<ID=LP,Number=A,Type=Float,Description="LP">', '##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">',
                   '##FORMAT=<ID=ID,Number=1,Type=String,Description="Study variant identifier">']
        raw = write_vcf(self.tmp / "noss.gwas.vcf", [("1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:rs1")], fmt="ES:SE:LP:AF:ID",
                        formats=formats, sample="ieu-b-1")
        result = run_standardise(self.tmp, raw, example_reference(genomes=True), source_type="opengwas",
                                 metadata={"gwasinfo": [{"id": "ieu-b-1", "sample_size": 46468, "build": "HG19/GRCh37"}]})
        self.assertEqual(result.returncode, 0, result.stderr)
        row = result.rows[0]
        self.assertEqual(row["FORMAT"]["ID"], "rs1")
        self.assertNotIn("SS", row["FORMAT"])  # per-variant SS only from the data; the study total goes in ##SAMPLE
        self.assertEqual((result.sample["TotalControls"], result.sample["StudyType"]), ("46468", "Continuous"))

    def test_sample_size_is_added_to_records_that_lack_it_even_when_the_header_declares_it(self):
        # OpenGWAS headers declare FORMAT/SS although many of their records carry no SS value.
        raw = write_vcf(self.tmp / "declared.gwas.vcf", [("1", 1100, "rs1", "G", "A", "0.1:0.02:6:0.3:rs1"),
                                                        ("1", 1101, "rs2", "G", "A", "0.1:0.02:6:0.3:777:rs2")],
                        fmt="ES:SE:LP:AF:ID", contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        text = raw.read_text(encoding="utf-8").replace("ES:SE:LP:AF:ID\t0.1:0.02:6:0.3:777:rs2", "ES:SE:LP:AF:SS:ID\t0.1:0.02:6:0.3:777:rs2")
        raw.write_text(text, encoding="utf-8")
        result = run_standardise(self.tmp, raw, example_reference(), source_type="opengwas",
                                 metadata={"gwasinfo": [{"id": "ieu-b-1", "sample_size": 46468}]})
        self.assertEqual(result.returncode, 0, result.stderr)
        sizes = {r["INFO"]["SNPID"]: r["FORMAT"].get("SS") for r in result.rows}
        self.assertEqual(sizes, {"rs1": ".", "rs2": "777"})
        self.assertEqual(sum(line.startswith("##FORMAT=<ID=SS,") for line in result.header), 1)
        self.assertEqual(result.sample["TotalControls"], "46468")

    @needs_liftover
    def test_source_build_headers_are_replaced_by_grch38_declarations(self):
        raw = write_vcf(self.tmp / "b37.gwas.vcf", [("1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")],
                        extra_header=["##reference=file:///data/b37/human_g1k_v37.fasta", "##study_id=stale", "##ancestry=stale"])
        result = run_standardise(self.tmp, raw, example_reference(genomes=True))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([h for h in result.header if h.startswith("##reference")], ["##reference=GRCh38.fa"])  # the GRCh38 FASTA used
        self.assertFalse(any("GRCh37" in h for h in result.header if h.startswith("##contig")))
        self.assertEqual(result.sample["StudyID"], "test_2026_1")
        self.assertTrue(all("assembly=GRCh38" in h for h in result.header if h.startswith("##contig")))

    def test_undeclared_build_that_cannot_be_detected_fails_loudly(self):
        raw = write_vcf(self.tmp / "nobuild.gwas.vcf", [("1", 123457, "rsX", "G", "A", "0.1:0.02:6:0.3:5000")], contig="")
        result = run_standardise(self.tmp, raw, example_reference())
        self.assertEqual(result.returncode, 78)
        self.assertIn("input_build", result.stderr)

    def test_alignment_to_the_grch38_reference_reorients_effects(self):
        sequence = ["N"] * 1500
        sequence[1100 - 1] = "A"  # GRCh38 base equals the source ALT allele
        reference = make_reference(self.tmp / "reference", lookup_rows=[], fasta={"1": "".join(sequence)})
        formats = ['##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
                   '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect">', '##FORMAT=<ID=SE,Number=A,Type=Float,Description="SE">',
                   '##FORMAT=<ID=LP,Number=A,Type=Float,Description="LP">', '##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">',
                   '##FORMAT=<ID=EZ,Number=A,Type=Float,Description="Z">', '##FORMAT=<ID=SS,Number=A,Type=Float,Description="SS">']
        raw = write_vcf(self.tmp / "align.gwas.vcf", [("1", 1100, "rs1", "G", "A", "0.1:0.02:6:0.3:5:5000", "AF=0.3")],
                        fmt="ES:SE:LP:AF:EZ:SS", formats=formats, contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, reference)
        self.assertEqual(result.returncode, 0, result.stderr)
        row = result.rows[0]
        self.assertEqual((row["REF"], row["ALT"]), ("A", "G"))
        self.assertEqual((row["FORMAT"]["ES"], row["FORMAT"]["AF"], row["FORMAT"]["EZ"], row["FORMAT"]["SE"]), ("-0.1", "0.7", "-5", "0.02"))
        self.assertEqual(result.qc["reference_alignment"]["swapped"], 1)

    @needs_liftover
    def test_no_liftable_records_is_a_data_error(self):
        raw = write_vcf(self.tmp / "far.gwas.vcf", [("1", 5_000_000, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")])
        result = run_standardise(self.tmp, raw, example_reference(genomes=True))
        self.assertEqual(result.returncode, 65)
        self.assertIn("no variants", result.stderr.lower())

    def test_genomic_inflation_is_recorded_for_passthrough_outputs(self):
        # LP values 2, 1.30103 and 0.69897 are P = 0.01, 0.05 and 0.2: median 0.05, so lambda = 8.44394650653 (R).
        raw = write_vcf(self.tmp / "lambda.gwas.vcf", [("1", 1000 + i, f"rs{i}", "G", "A", f"0.1:0.02:{lp}:0.3:5000")
                                                      for i, lp in enumerate([2, 1.30103, 0.69897])],
                        contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, example_reference())
        self.assertAlmostEqual(result.qc["lambda_gc"], 8.44394650653, places=3)

    def test_duplicates_spanning_a_chunk_boundary_are_removed(self):
        # With 2 records per chunk, the swapped-allele duplicate at 2000 starts the second chunk.
        raw = write_vcf(self.tmp / "chunks.gwas.vcf", [
            ("1", 1000, "rs1", "A", "G", "0.1:0.02:6:0.3:5000"), ("1", 2000, "rs2", "C", "T", "0.1:0.02:6:0.3:5000"),
            ("1", 2000, "rs2", "T", "C", "-0.1:0.02:6:0.7:5000"), ("1", 3000, "rs3", "G", "A", "0.1:0.02:6:0.3:5000")],
            contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        env = {**os.environ, "GWAS_PASSTHROUGH_CHUNK": "2"}
        result = run_standardise(self.tmp, raw, example_reference(), env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([(r["POS"], r["REF"], r["ALT"]) for r in result.rows], [(1000, "A", "G"), (2000, "C", "T"), (3000, "G", "A")])
        self.assertEqual(result.qc["duplicates_removed"], 1)

    def test_alignment_counts_add_up_across_chunks(self):
        sequence = ["N"] * 3500
        for position in (1000, 2000, 3000):
            sequence[position - 1] = "A"  # the source ALT is the GRCh38 base at every record
        reference = make_reference(self.tmp / "reference", lookup_rows=[], fasta={"1": "".join(sequence)})
        raw = write_vcf(self.tmp / "align_chunks.gwas.vcf", [("1", pos, f"rs{pos}", "G", "A", "0.1:0.02:6:0.3:5000") for pos in (1000, 2000, 3000)],
                        contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, reference, env={**os.environ, "GWAS_PASSTHROUGH_CHUNK": "1"})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.qc["reference_alignment"]["swapped"], 3)
        self.assertEqual({(r["REF"], r["ALT"], r["FORMAT"]["ES"]) for r in result.rows}, {("A", "G", "-0.1")})

    def test_no_intermediate_files_are_left_in_the_output_or_temporary_folders(self):
        scratch = self.tmp / "scratch"
        scratch.mkdir()
        raw = write_vcf(self.tmp / "clean.gwas.vcf", [("1", 1100, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")], contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, example_reference(), env={**os.environ, "GWAS_TMPDIR": str(scratch)})
        self.assertEqual(result.returncode, 0, result.stderr)
        outputs = sorted(p.name for p in Path(result.record["standard_vcf"]).parent.iterdir())
        self.assertEqual(outputs, ["standardised.gwas.vcf.gz", "standardised.gwas.vcf.gz.tbi", "standardised.metadata.json",
                                   "standardised.qc.json"])
        self.assertEqual(list(scratch.iterdir()), [])

    def test_chromosomes_are_written_in_numeric_order(self):
        reference = make_reference(self.tmp / "reference")
        raw = write_vcf(self.tmp / "order.gwas.vcf", [("10", 50, "rs10", "G", "A", "0.1:0.02:6:0.3:5000"),
                                                     ("2", 50, "rs2", "G", "A", "0.1:0.02:6:0.3:5000"),
                                                     ("1", 50, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")],
                        contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, reference)
        self.assertEqual([r["CHROM"] for r in result.rows], ["1", "2", "10"])

if __name__ == "__main__":
    unittest.main()
