"""Tests for report.py: -log10 P extraction, QQ thinning, genomic inflation and plot regeneration."""

import json
import os
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np

from testlib import EXAMPLES, run_script, run_standardise, write_table

import report


class LpExtractionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_table_and_vcf_values_are_read_in_full(self):
        table = write_table(self.tmp / "p.tsv", ["SNP", "P"], [[f"rs{i}", (i + 1) / 1000] for i in range(50)])
        self.assertEqual(len(report.extract_lp(table)), 50)
        vcf = self.tmp / "p.gwas.vcf"
        body = "".join(f"1\t{1000 + i}\trs{i}\tA\tG\t.\tPASS\t.\tES:SE:LP:AF:SS\t0.1:0.01:{(i % 5) + 0.1}:0.3:1000\n" for i in range(30))
        vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBMI\n" + body, encoding="utf-8")
        np.testing.assert_allclose(report.extract_lp(vcf)[:3], [0.1, 1.1, 2.1])

    def test_raw_table_uses_the_standardisers_columns_and_keeps_extreme_values(self):
        table = write_table(self.tmp / "raw.tsv", ["SNP", "P", "mlog10p_extreme"], [["rs1", 0.5, 450.0], ["rs2", 0.5, 12.5]])
        lp = report.extract_lp(table, {"LOG10P": "mlog10p_extreme"}, "\t")
        np.testing.assert_allclose(lp, [450.0, 12.5])

    def test_underflowed_p_is_recovered_from_beta_and_se(self):
        table = write_table(self.tmp / "raw.tsv", ["SNP", "P", "BETA", "SE"], [["rs1", 0, 0.5, 0.01]])
        lp = report.extract_lp(table, {"P": "P", "BETA": "BETA", "SE": "SE"}, "\t")
        self.assertAlmostEqual(lp[0], 544.665305866, delta=1e-3)


class QqThinningTests(unittest.TestCase):
    def test_every_point_above_lp_two_is_kept_and_the_bulk_is_thinned(self):
        rng = np.random.default_rng(821)
        lp = -np.log10(rng.uniform(size=100_000))
        expected, observed = report.qq_points(lp)
        self.assertEqual(int((observed >= 2).sum()), int((lp >= 2).sum()))
        self.assertLess(len(observed), 5_000)
        self.assertAlmostEqual(float(observed.max()), float(lp.max()))
        self.assertAlmostEqual(float(expected.max()), -np.log10(0.5 / 100_000))


class PlotTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        result = run_standardise(self.tmp, EXAMPLES / "tables/generic.tsv", EXAMPLES / "reference", input_build="hg19")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.record = Path(result.record["standard_vcf"]).parent.parent / "record.json"

    def plot(self, *extra):
        out = self.tmp / "plot.json"
        run_script("report.py", "--standardise-record", self.record, "--output", out, *extra)
        return json.loads(out.read_text(encoding="utf-8"))

    def test_self_contained_html_report_is_written_with_steps_and_qq_plots(self):
        result = self.plot()
        self.assertTrue(result["report_created"])
        text = Path(result["report_html"]).read_text(encoding="utf-8")
        self.assertIn("data:image/png;base64,", text)
        self.assertIn("Palindromic (A/T, C/G) SNPs", text)
        self.assertIn("TotalVariants", text)
        self.assertNotIn("http://", text.replace("http://www.w3.org", ""))  # no external links or files
        self.assertNotIn("<script src", text)
        self.assertIn("lambda_gc_raw", result)
        self.assertIn("lambda_gc_standardised", result)

    def test_existing_report_is_kept_only_while_it_is_newer_than_the_vcf(self):
        first = self.plot()
        self.assertFalse(self.plot()["report_created"])
        later = time.time() + 5
        os.utime(first["standard_vcf"], (later, later))
        self.assertTrue(self.plot()["report_created"])

    def test_a_report_can_be_made_after_processing_from_the_output_folder(self):
        folder = Path(json.loads(self.record.read_text(encoding="utf-8"))["standard_vcf"]).parent
        run_script("report.py", "--dir", folder)
        self.assertIn("<h2>Steps</h2>", (folder / "standardised.report.html").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
