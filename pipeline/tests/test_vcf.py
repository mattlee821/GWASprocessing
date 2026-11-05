"""Tests for gwas_vcf: header contract, body formatting, ordering and bgzip/tabix output."""

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

RECORD = {"study_id": "test_2026_1", "phenotype": "Body mass index\t(BMI)\nline2", "ancestry": "EUR", "sex": "", "PMID": "1"}


def table(**columns):
    n = len(next(iter(columns.values())))
    base = {"CHR": ["1"] * n, "POS38": [1100 + i for i in range(n)], "POS19": [1000 + i for i in range(n)],
            "SNP": [f"rs{i}" for i in range(n)], "SNPID": [pd.NA] * n, "EA": ["A"] * n, "OA": ["G"] * n,
            "EAF": [0.3] * n, "EAF_SOURCE": ["input"] * n, "BETA": [0.1] * n, "SE": [0.02] * n, "LP": [6.0] * n, "N": [5000] * n}
    base.update(columns)
    return pd.DataFrame(base)


class HeaderTests(unittest.TestCase):
    def test_header_declares_grch38_contigs_and_sanitised_metadata(self):
        header = gwas_vcf.table_header(RECORD, extras=[])
        self.assertIn("##reference=GRCh38", header)
        contigs = [line for line in header if line.startswith("##contig=")]
        self.assertEqual(len(contigs), 25)
        self.assertIn("##contig=<ID=1,length=248956422,assembly=GRCh38>", contigs)
        self.assertIn("##contig=<ID=MT,length=16569,assembly=GRCh38>", contigs)
        self.assertIn("##phenotype=Body mass index (BMI) line2", header)
        self.assertIn("##sex=NA", header)
        self.assertTrue(header[-1].endswith("\tFORMAT\tBody_mass_index_(BMI)_line2"))

    def test_extras_are_declared_as_info_strings(self):
        header = gwas_vcf.table_header(RECORD, extras=["INFO", "9score"])
        self.assertIn('##INFO=<ID=INFO,Number=1,Type=String,Description="Source column INFO">', header)
        self.assertIn('##INFO=<ID=X_9score,Number=1,Type=String,Description="Source column 9score">', header)


class BodyTests(unittest.TestCase):
    def test_lines_carry_per_row_eaf_source_lp_and_sanitised_info(self):
        df = table(EAF=[0.3, 0.8, np.nan], EAF_SOURCE=["input", "reference", ""], LP=[544.665305866, 6.0, np.nan],
                   N=[5000.4, np.nan, 10], SNPID=["1:100;x", pd.NA, ""], SNP=["rs1", pd.NA, "1:200"],
                   INFO=["0.98", "0.5 ; x=y", np.nan])
        rows = parse_vcf_rows("\n".join(gwas_vcf.table_body_lines(df, extras=["INFO"])))
        self.assertEqual([r["ID"] for r in rows], ["rs1", ".", "."])
        self.assertEqual(rows[0]["INFO"], {"POS19": "1000", "SNPID": "1:100_x", "EAF_SOURCE": "input", "INFO": "0.98"})
        self.assertEqual(rows[1]["INFO"]["EAF_SOURCE"], "reference")
        self.assertEqual(rows[1]["INFO"]["INFO"], "0.5___x_y")
        self.assertEqual(rows[2]["INFO"]["EAF_SOURCE"], ".")
        self.assertEqual(rows[2]["INFO"]["SNPID"], ".")
        self.assertEqual(rows[0]["FORMAT"], {"ES": "0.1", "SE": "0.02", "LP": "544.66531", "AF": "0.3", "SS": "5000"})
        self.assertEqual((rows[1]["FORMAT"]["SS"], rows[2]["FORMAT"]["LP"], rows[2]["FORMAT"]["AF"]), (".", ".", "."))
        self.assertEqual((rows[0]["REF"], rows[0]["ALT"]), ("G", "A"))

    def test_rows_sort_in_numeric_chromosome_order(self):
        df = table(CHR=["10", "2", "X", "1", "MT", "2"], POS38=[5, 9, 1, 7, 3, 4])
        ordered = gwas_vcf.sort_frame(df)
        self.assertEqual(list(zip(ordered["CHR"], ordered["POS38"])), [("1", 7), ("2", 4), ("2", 9), ("10", 5), ("X", 1), ("MT", 3)])


@unittest.skipUnless(shutil.which("bgzip") and shutil.which("tabix") and shutil.which("bcftools"), "htslib tools required")
class WriterTests(unittest.TestCase):
    def test_written_file_is_valid_bgzipped_and_indexed(self):
        tmp = Path(tempfile.mkdtemp())
        vcf = tmp / "out.gwas.vcf.gz"
        gwas_vcf.write_table_vcf(table(CHR=["2", "1", "1"], POS38=[50, 70, 60]), vcf, RECORD, extras=[])
        self.assertTrue(Path(str(vcf) + ".tbi").exists())
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

    def test_plain_vcf_body_is_sorted_numerically_in_place(self):
        tmp = Path(tempfile.mkdtemp())
        plain = tmp / "body.vcf"
        plain.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n10\t5\nX\t1\n2\t30\n2\t4\n1\t9\n", encoding="utf-8")
        gwas_vcf.sort_vcf_body(plain)
        lines = plain.read_text(encoding="utf-8").splitlines()
        self.assertEqual(lines[:2], ["##fileformat=VCFv4.2", "#CHROM\tPOS"])
        self.assertEqual(lines[2:], ["1\t9", "2\t4", "2\t30", "10\t5", "X\t1"])


if __name__ == "__main__":
    unittest.main()
