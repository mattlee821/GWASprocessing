"""End-to-end tests for the passthrough path of standardise_gwas.py (pre-built GWAS-VCF inputs, e.g. OpenGWAS)."""

import os
import tempfile
import unittest
from pathlib import Path

from testlib import EXAMPLES, chain_block, make_reference, run_standardise, write_chain, write_lookup

import standardise_gwas


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

    def test_existing_gwas_vcf_gets_manifest_sex_metadata_and_detected_build(self):
        raw = write_vcf(self.tmp / "existing_gwas.vcf", [(1, 1000, "rs1", "G", "A", "0.2:0.05:3:0.1:1000", "POS19=900;SNPID=1:1000")],
                        contig=None, extra_header=["##source=test"])
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference", sex="male", ancestry="EUR")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.metadata["sex"], "male")
        self.assertIn("##sex=male", result.header)
        row = result.rows[0]
        self.assertEqual((row["POS"], row["INFO"]["POS19"]), (1100, "1000"))
        self.assertEqual(result.qc["source_build"], "GRCh37")

    def test_passthrough_hg19_lifted_to_grch38(self):
        reference = self.tmp / "reference"
        write_chain(reference / "genome" / "hg19ToHg38.over.chain.gz",
                    [chain_block("chr1", "chr1", 249250621, 248956422, 700000, 100000, 764620, 1)])
        raw = write_vcf(self.tmp / "grch37.gwas.vcf", [("1", 721290, "rs12565286", "G", "C", "0.01:0.001:2:0.1:1000")])
        result = run_standardise(self.tmp, raw, reference)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.metadata["coordinate_build"], "GRCh38")
        row = result.rows[0]
        self.assertEqual((row["POS"], row["ID"], row["INFO"]["POS19"]), (785910, "rs12565286", "721290"))
        self.assertIn("##ancestry=ALL", result.header)
        self.assertIn("##sex=combined", result.header)

    def test_passthrough_output_is_position_sorted(self):
        raw = write_vcf(self.tmp / "unsorted.gwas.vcf", [("1", 5000, "rsB", "A", "G", "0.02:0.01:1:0.3:1000"),
                                                        ("1", 2000, "rsA", "C", "T", "0.01:0.01:1:0.3:1000"),
                                                        ("1", 9000, "rsC", "G", "A", "0.03:0.01:1:0.3:1000")],
                        contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r["POS"] for r in result.rows], [2000, 5000, 9000])
        self.assertEqual(result.rows[0]["INFO"]["POS19"], "1900")

    def test_passthrough_palindromic_strand_resolution(self):
        reference = self.tmp / "reference"
        write_lookup(reference, [("1", "rsP1", 1000, "A", "T", 0.20), ("1", "rsP2", 2000, "C", "G", 0.15), ("1", "rsP3", 3000, "A", "T", 0.47)])
        write_chain(reference / "genome" / "hg38ToHg19.over.chain.gz",
                    [chain_block("chr1", "chr1", 300000000, 300000000, 0, 300000000, 0, 1)])
        raw = write_vcf(self.tmp / "pass.gwas.vcf", [
            ("1", 1000, "rsP1", "T", "A", "0.1:0.01:1:0.20:1000"), ("1", 2000, "rsP2", "G", "C", "0.1:0.01:1:0.85:1000"),
            ("1", 3000, "rsP3", "T", "A", "0.1:0.01:1:0.47:1000"), ("1", 4000, "rsN1", "G", "A", "0.1:0.01:1:0.30:1000")],
            contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, reference, analysis_qc=True, source_type="opengwas")
        self.assertEqual(result.returncode, 0, result.stderr)
        variants = {r["ID"]: (r["REF"], r["ALT"]) for r in result.rows}
        self.assertEqual(variants.get("rsP1"), ("T", "A"))
        self.assertEqual(variants.get("rsP2"), ("C", "G"))
        self.assertNotIn("rsP3", variants)
        self.assertIn("rsN1", variants)
        detail = result.qc["palindromic_resolution"]
        self.assertEqual((detail["flipped"], detail["kept_same_strand"], detail["dropped_ambiguous"]), (1, 1, 1))


class PassthroughRegressionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_chr_prefixed_grch37_contigs_are_renamed_and_lifted(self):
        raw = write_vcf(self.tmp / "chr.gwas.vcf", [("chr1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")],
                        contig="##contig=<ID=chr1,length=249250621,assembly=GRCh37>")
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference")
        self.assertEqual(result.returncode, 0, result.stderr)
        row = result.rows[0]
        self.assertEqual((row["CHROM"], row["POS"], row["INFO"]["POS19"]), ("1", 1100, "1000"))

    def test_chr_prefixed_grch38_contigs_are_renamed(self):
        raw = write_vcf(self.tmp / "chr38.gwas.vcf", [("chr1", 2000, "rs2", "C", "T", "0.1:0.02:6:0.3:5000")],
                        contig="##contig=<ID=chr1,length=248956422,assembly=GRCh38>")
        row = run_standardise(self.tmp, raw, EXAMPLES / "reference").rows[0]
        self.assertEqual((row["CHROM"], row["POS"], row["INFO"]["POS19"]), ("1", 2000, "1900"))

    def test_minus_strand_liftover_complements_ref_and_alt(self):
        reference = self.tmp / "reference"
        write_chain(reference / "genome" / "hg19ToHg38.over.chain.gz",
                    [chain_block("chr1", "chr1", 249250621, 248956422, 7000, 1000, 10000, 1, dst_strand="-")])
        raw = write_vcf(self.tmp / "minus.gwas.vcf", [("1", 7500, "rs51", "C", "A", "0.1:0.02:6:0.3:5000")])
        row = run_standardise(self.tmp, raw, reference).rows[0]
        self.assertEqual((row["POS"], row["REF"], row["ALT"], row["FORMAT"]["ES"]), (248945923, "G", "T", "0.1"))

    def test_sample_size_is_added_from_opengwas_metadata(self):
        formats = ['##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect">', '##FORMAT=<ID=SE,Number=A,Type=Float,Description="SE">',
                   '##FORMAT=<ID=LP,Number=A,Type=Float,Description="LP">', '##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">',
                   '##FORMAT=<ID=ID,Number=1,Type=String,Description="Study variant identifier">']
        raw = write_vcf(self.tmp / "noss.gwas.vcf", [("1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:rs1")], fmt="ES:SE:LP:AF:ID",
                        formats=formats, sample="ieu-b-1")
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference", source_type="opengwas",
                                 metadata={"gwasinfo": [{"id": "ieu-b-1", "sample_size": 46468, "build": "HG19/GRCh37"}]})
        self.assertEqual(result.returncode, 0, result.stderr)
        row = result.rows[0]
        self.assertEqual((row["FORMAT"]["SS"], row["FORMAT"]["ID"]), ("46468", "rs1"))
        self.assertTrue(any(line.startswith("##FORMAT=<ID=SS,") for line in result.header))
        self.assertEqual(result.qc["sample_size_source"], "gwasinfo")

    def test_sample_size_is_added_to_records_that_lack_it_even_when_the_header_declares_it(self):
        # OpenGWAS headers declare FORMAT/SS although many of their records carry no SS value.
        raw = write_vcf(self.tmp / "declared.gwas.vcf", [("1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:rs1"),
                                                        ("1", 1001, "rs2", "G", "A", "0.1:0.02:6:0.3:777:rs2")],
                        fmt="ES:SE:LP:AF:ID")
        text = raw.read_text(encoding="utf-8").replace("ES:SE:LP:AF:ID\t0.1:0.02:6:0.3:777:rs2", "ES:SE:LP:AF:SS:ID\t0.1:0.02:6:0.3:777:rs2")
        raw.write_text(text, encoding="utf-8")
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference", source_type="opengwas",
                                 metadata={"gwasinfo": [{"id": "ieu-b-1", "sample_size": 46468}]})
        self.assertEqual(result.returncode, 0, result.stderr)
        sizes = {r["ID"]: r["FORMAT"].get("SS") for r in result.rows}
        self.assertEqual(sizes, {"rs1": "46468", "rs2": "777"})
        self.assertEqual(sum(line.startswith("##FORMAT=<ID=SS,") for line in result.header), 1)
        self.assertEqual(result.qc["sample_size_filled_records"], 1)

    def test_source_build_headers_are_replaced_by_grch38_declarations(self):
        raw = write_vcf(self.tmp / "b37.gwas.vcf", [("1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")],
                        extra_header=["##reference=file:///data/b37/human_g1k_v37.fasta", "##study_id=stale", "##ancestry=stale"])
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([h for h in result.header if h.startswith("##reference")], ["##reference=GRCh38"])
        self.assertFalse(any("GRCh37" in h for h in result.header if h.startswith("##contig")))
        self.assertEqual([h for h in result.header if h.startswith("##study_id")], ["##study_id=test_2026_1"])
        out_vcf = Path(result.record["standard_vcf"])
        self.assertEqual(standardise_gwas.detect_vcf_build(out_vcf), "GRCh38")

    def test_undeclared_build_that_cannot_be_detected_fails_loudly(self):
        raw = write_vcf(self.tmp / "nobuild.gwas.vcf", [("1", 123457, "rsX", "G", "A", "0.1:0.02:6:0.3:5000")], contig="")
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference")
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
        raw = write_vcf(self.tmp / "align.gwas.vcf", [("1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:5:5000", "AF=0.3")],
                        fmt="ES:SE:LP:AF:EZ:SS", formats=formats)
        result = run_standardise(self.tmp, raw, reference)
        self.assertEqual(result.returncode, 0, result.stderr)
        row = result.rows[0]
        self.assertEqual((row["REF"], row["ALT"]), ("A", "G"))
        self.assertEqual((row["FORMAT"]["ES"], row["FORMAT"]["AF"], row["FORMAT"]["EZ"], row["FORMAT"]["SE"]), ("-0.1", "0.7", "-5", "0.02"))
        self.assertEqual(row["INFO"]["AF"], "0.7")
        self.assertEqual(result.qc["reference_alignment"]["swapped"], 1)

    def test_no_liftable_records_is_a_data_error(self):
        raw = write_vcf(self.tmp / "far.gwas.vcf", [("1", 5_000_000, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")])
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference")
        self.assertEqual(result.returncode, 65)
        self.assertIn("no variants", result.stderr.lower())

    def test_genomic_inflation_is_recorded_for_passthrough_outputs(self):
        # LP values 2, 1.30103 and 0.69897 are P = 0.01, 0.05 and 0.2: median 0.05, so lambda = 8.44394650653 (R).
        raw = write_vcf(self.tmp / "lambda.gwas.vcf", [("1", 1000 + i, f"rs{i}", "G", "A", f"0.1:0.02:{lp}:0.3:5000")
                                                      for i, lp in enumerate([2, 1.30103, 0.69897])],
                        contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference")
        self.assertAlmostEqual(result.qc["lambda_gc"], 8.44394650653, places=3)

    def test_duplicates_spanning_a_chunk_boundary_are_removed(self):
        # With 2 records per chunk, the swapped-allele duplicate at 2000 starts the second chunk.
        raw = write_vcf(self.tmp / "chunks.gwas.vcf", [
            ("1", 1000, "rs1", "A", "G", "0.1:0.02:6:0.3:5000"), ("1", 2000, "rs2", "C", "T", "0.1:0.02:6:0.3:5000"),
            ("1", 2000, "rs2", "T", "C", "-0.1:0.02:6:0.7:5000"), ("1", 3000, "rs3", "G", "A", "0.1:0.02:6:0.3:5000")],
            contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        env = {**os.environ, "GWAS_PASSTHROUGH_CHUNK": "2"}
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference", analysis_qc=True, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([(r["POS"], r["REF"], r["ALT"]) for r in result.rows], [(1000, "A", "G"), (2000, "C", "T"), (3000, "G", "A")])
        self.assertEqual(result.qc["qc_removed"]["duplicate"], 1)

    def test_alignment_counts_add_up_across_chunks(self):
        sequence = ["N"] * 3500
        for position in (1000, 2000, 3000):
            sequence[position - 1] = "A"  # the source ALT is the GRCh38 base at every record
        reference = make_reference(self.tmp / "reference", lookup_rows=[], fasta={"1": "".join(sequence)})
        write_chain(reference / "genome" / "hg38ToHg19.over.chain.gz", [chain_block("chr1", "chr1", 300000000, 300000000, 0, 300000000, 0, 1)])
        raw = write_vcf(self.tmp / "align_chunks.gwas.vcf", [("1", pos, f"rs{pos}", "G", "A", "0.1:0.02:6:0.3:5000") for pos in (1000, 2000, 3000)],
                        contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, reference, env={**os.environ, "GWAS_PASSTHROUGH_CHUNK": "1"})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.qc["reference_alignment"]["swapped"], 3)
        self.assertEqual({(r["REF"], r["ALT"], r["FORMAT"]["ES"]) for r in result.rows}, {("A", "G", "-0.1")})

    def test_no_intermediate_files_are_left_in_the_output_or_temporary_folders(self):
        scratch = self.tmp / "scratch"
        scratch.mkdir()
        raw = write_vcf(self.tmp / "clean.gwas.vcf", [("1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")])
        result = run_standardise(self.tmp, raw, EXAMPLES / "reference", analysis_qc=True, env={**os.environ, "GWAS_TMPDIR": str(scratch)})
        self.assertEqual(result.returncode, 0, result.stderr)
        outputs = sorted(p.name for p in Path(result.record["standard_vcf"]).parent.iterdir())
        self.assertEqual(outputs, ["standardised.gwas.vcf.gz", "standardised.gwas.vcf.gz.tbi", "standardised.metadata.json",
                                   "standardised.qc.json", "standardised.steps.txt"])
        self.assertEqual(list(scratch.iterdir()), [])

    def test_chromosomes_are_written_in_numeric_order(self):
        reference = self.tmp / "reference"
        write_chain(reference / "genome" / "hg38ToHg19.over.chain.gz",
                    [chain_block(f"chr{c}", f"chr{c}", 300000000, 300000000, 0, 300000000, 0, i) for i, c in enumerate(["1", "2", "10"], 1)])
        raw = write_vcf(self.tmp / "order.gwas.vcf", [("10", 50, "rs10", "G", "A", "0.1:0.02:6:0.3:5000"),
                                                     ("2", 50, "rs2", "G", "A", "0.1:0.02:6:0.3:5000"),
                                                     ("1", 50, "rs1", "G", "A", "0.1:0.02:6:0.3:5000")],
                        contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
        result = run_standardise(self.tmp, raw, reference)
        self.assertEqual([r["CHROM"] for r in result.rows], ["1", "2", "10"])


if __name__ == "__main__":
    unittest.main()
