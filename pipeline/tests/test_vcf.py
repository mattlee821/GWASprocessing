"""Tests for gwas_vcf: GWAS-VCF 1.0.0 header, body formatting, ordering and bgzip/tabix output."""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from testlib import parse_vcf_rows, read_vcf_text

import gwas_vcf

RECORD = {"study_id": "test_2026_1", "output_id": "local_test 1", "phenotype": "Body mass index\t(BMI)\nline2",
          "ancestry": "EUR", "sex": "", "PMID": "1"}
KEYS = gwas_vcf.format_keys(["N"])
FIELDS = {"TotalVariants": 10, "VariantsNotRead": 2, "HarmonisedVariants": 7, "VariantsNotHarmonised": 1, "SwitchedAlleles": 3,
          "AFFilled": 1, "TotalControls": 5000, "StudyType": "Continuous", "StudyID": "test_2026_1",
          "Phenotype": "Body mass index (BMI)", "Ancestry": "EUR", "Sex": "", "PMID": "1"}


def table(**columns):
    n = len(next(iter(columns.values())))
    base = {"CHR": ["1"] * n, "POS38": [1100 + i for i in range(n)], "SNP": [f"rs{i}" for i in range(n)], "SNPID": [pd.NA] * n,
            "EA": ["A"] * n, "OA": ["G"] * n, "EAF": [0.3] * n, "EAF_SOURCE": ["input"] * n, "BETA": [0.1] * n,
            "SE": [0.02] * n, "LP": [6.0] * n, "N": [5000] * n, "HM": ["forward"] * n}
    base.update(columns)
    return pd.DataFrame(base)


class HeaderTests(unittest.TestCase):
    def setUp(self):
        self.header = gwas_vcf.table_header(RECORD, [], KEYS, FIELDS, "Homo_sapiens.GRCh38.dna.primary_assembly.fa", "abc123")

    def test_header_declares_grch38_contigs_reference_and_settings(self):
        self.assertIn("##reference=Homo_sapiens.GRCh38.dna.primary_assembly.fa", self.header)
        self.assertIn("##settings_hash=abc123", self.header)
        self.assertIn("##GWASVCF=1.0", self.header)
        contigs = [line for line in self.header if line.startswith("##contig=")]
        self.assertEqual(len(contigs), 25)
        self.assertIn("##contig=<ID=1,length=248956422,assembly=GRCh38>", contigs)
        self.assertIn("##contig=<ID=MT,length=16569,assembly=GRCh38>", contigs)

    def test_sample_is_the_output_id_with_a_sample_line_and_meta_for_every_key(self):
        self.assertTrue(self.header[-1].endswith("\tFORMAT\tlocal_test_1"))
        sample = [line for line in self.header if line.startswith("##SAMPLE=")]
        self.assertEqual(len(sample), 1)
        self.assertTrue(sample[0].startswith("##SAMPLE=<ID=local_test_1,TotalVariants=10,VariantsNotRead=2,"))
        self.assertIn('Phenotype="Body mass index (BMI)"', sample[0])
        self.assertIn("Sex=NA", sample[0])
        meta = {line.split(",")[0].split("=")[-1] for line in self.header if line.startswith("##META=")}
        self.assertEqual(meta, set(FIELDS))

    def test_format_and_info_are_declared(self):
        formats = [line.split(",")[0].split("=")[-1] for line in self.header if line.startswith("##FORMAT=")]
        self.assertEqual(formats, ["ES", "SE", "LP", "AF", "AFS", "SS", "ID"])
        infos = [line.split(",")[0].split("=")[-1] for line in self.header if line.startswith("##INFO=")]
        self.assertEqual(infos, ["SNPID", "HM"])
        self.assertFalse(any("POS19" in line or "EAF_SOURCE" in line for line in self.header))

    def test_optional_format_keys_follow_the_columns_present(self):
        self.assertEqual(gwas_vcf.format_keys([]), ["ES", "SE", "LP", "AF", "AFS", "ID"])
        self.assertEqual(gwas_vcf.format_keys(["N", "NCASE", "IMPINFO"]), ["ES", "SE", "LP", "AF", "AFS", "SS", "NC", "SI", "ID"])

    def test_extras_are_declared_as_info_strings(self):
        header = gwas_vcf.table_header(RECORD, ["INFO", "9score", "HM"], KEYS, FIELDS)
        self.assertIn('##INFO=<ID=INFO,Number=1,Type=String,Description="Source column INFO">', header)
        self.assertIn('##INFO=<ID=X_9score,Number=1,Type=String,Description="Source column 9score">', header)
        self.assertIn('##INFO=<ID=SRC_HM,Number=1,Type=String,Description="Source column HM">', header)


class BodyTests(unittest.TestCase):
    def test_lines_carry_af_source_rsid_format_id_and_sanitised_info(self):
        df = table(EAF=[0.3, 0.8, np.nan], EAF_SOURCE=["input", "reference", ""], LP=[544.665305866, 6.0, np.nan],
                   N=[5000.4, np.nan, 10], SNPID=["1:100;x", pd.NA, ""], SNP=["rs1", pd.NA, "1:200"],
                   HM=["forward", "palindrome_ambiguous_swapped", "swapped"], INFO=["0.98", "0.5 ; x=y", np.nan])
        rows = parse_vcf_rows("\n".join(gwas_vcf.table_body_lines(df, ["INFO"], KEYS, "1000G_EUR")))
        self.assertEqual([r["ID"] for r in rows], ["rs1", ".", "."])
        self.assertEqual(rows[0]["INFO"], {"SNPID": "1:100_x", "HM": "forward", "INFO": "0.98"})
        self.assertEqual(rows[1]["INFO"]["HM"], "palindrome_ambiguous_swapped")
        self.assertEqual(rows[1]["INFO"]["INFO"], "0.5___x_y")
        self.assertEqual(rows[2]["INFO"]["SNPID"], ".")
        self.assertEqual(rows[0]["FORMAT"], {"ES": "0.1", "SE": "0.02", "LP": "544.66531", "AF": "0.3", "AFS": "study",
                                             "SS": "5000.4", "ID": "rs1"})
        self.assertEqual((rows[1]["FORMAT"]["AF"], rows[1]["FORMAT"]["AFS"], rows[1]["FORMAT"]["SS"]), ("0.8", "1000G_EUR", "."))
        self.assertEqual((rows[2]["FORMAT"]["LP"], rows[2]["FORMAT"]["AF"], rows[2]["FORMAT"]["AFS"], rows[2]["FORMAT"]["ID"]),
                         (".", ".", ".", "."))

    def test_filled_frequencies_are_written_with_6_significant_digits_and_study_ones_unchanged(self):
        # 1000 Genomes frequencies arrive as float32, so 0.8744980 reads back as 0.87449801
        filled = [float(np.float32(0.874498)), float(np.float32(0.3)), float(np.float32(0.0592369))]
        df = table(EAF=filled + [0.123456789], EAF_SOURCE=["reference", "reference", "reference", "input"])
        rows = parse_vcf_rows("\n".join(gwas_vcf.table_body_lines(df, [], KEYS, "1000G_EUR")))
        self.assertEqual([r["FORMAT"]["AF"] for r in rows], ["0.874498", "0.3", "0.0592369", "0.12345679"])
        self.assertEqual([r["FORMAT"]["AFS"] for r in rows], ["1000G_EUR", "1000G_EUR", "1000G_EUR", "study"])
        self.assertEqual((rows[0]["REF"], rows[0]["ALT"]), ("G", "A"))

    def test_rows_sort_in_numeric_chromosome_order(self):
        df = table(CHR=["10", "2", "X", "1", "MT", "2"], POS38=[5, 9, 1, 7, 3, 4])
        ordered = gwas_vcf.sort_frame(df)
        self.assertEqual(list(zip(ordered["CHR"], ordered["POS38"])), [("1", 7), ("2", 4), ("2", 9), ("10", 5), ("X", 1), ("MT", 3)])


@unittest.skipUnless(shutil.which("bgzip") and shutil.which("tabix") and shutil.which("bcftools"), "htslib tools required")
class WriterTests(unittest.TestCase):
    def test_header_and_body_are_bgzipped_indexed_and_valid(self):
        tmp = Path(tempfile.mkdtemp())
        vcf, header, body = tmp / "out.gwas.vcf.gz", tmp / "header.vcf", tmp / "body.vcf"
        header.write_text("\n".join(gwas_vcf.table_header(RECORD, [], KEYS, FIELDS)) + "\n", encoding="utf-8")
        ordered = gwas_vcf.sort_frame(table(CHR=["2", "1", "1"], POS38=[50, 70, 60]))
        body.write_text("\n".join(gwas_vcf.table_body_lines(ordered, [], KEYS)) + "\n", encoding="utf-8")
        gwas_vcf.bgzip_and_index([header, body], vcf)
        self.assertTrue(Path(str(vcf) + ".tbi").exists())
        self.assertFalse(header.exists() or body.exists())
        listed = subprocess.run(["tabix", "-l", str(vcf)], capture_output=True, text=True, check=True).stdout.split()
        self.assertEqual(listed, ["1", "2"])
        check = subprocess.run(["bcftools", "view", "-Ou", "-o", os.devnull, str(vcf)], capture_output=True, text=True)
        self.assertEqual((check.returncode, check.stderr.strip()), (0, ""))
        positions = [(r["CHROM"], r["POS"]) for r in parse_vcf_rows(read_vcf_text(vcf))]
        self.assertEqual(positions, [("1", 60), ("1", 70), ("2", 50)])

    def test_missing_htslib_is_a_configuration_error(self):
        path = os.environ["PATH"]
        try:
            os.environ["PATH"] = os.pathsep.join(p for p in path.split(os.pathsep) if not Path(p, "bgzip").exists())
            with self.assertRaises(SystemExit) as caught:
                gwas_vcf.require_htslib()
            self.assertEqual(caught.exception.code, 78)
        finally:
            os.environ["PATH"] = path


if __name__ == "__main__":
    unittest.main()
