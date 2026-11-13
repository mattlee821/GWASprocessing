"""Step 3, Parse and derive (README "Pipeline in detail", 3. Parse and derive).

BETA is log OR for binary traits. SE comes from a confidence interval, (log upper - log lower) / (2 x 1.96) for odds
ratios or (upper - lower) / (2 x 1.96) for betas, or else from abs(BETA / qnorm(P / 2)). LP (-log10 P) is read from a
LOG10P column or computed from P, and recovered from BETA/SE when P is 0 or missing.

Expected values were computed independently in R:
    log(1.2); (log(1.3091) - log(1.1)) / (2 * qnorm(0.975)); (0.15 - 0.05) / (2 * qnorm(0.975)); 0.1 / qnorm(1 - 1e-6 / 2)
    lp <- function(z) -(log(2) + pnorm(-abs(z), log.p = TRUE)) / log(10)
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from testlib import make_reference, run_standardise, write_table

import gwas_stats


class Step03ParseTest(unittest.TestCase):
    def test_parse_and_derive(self):
        tmp = Path(tempfile.mkdtemp(prefix="step03_"))
        reference = make_reference(tmp / "reference", lookup_rows=[])

        # 1 row per derivation, positions given in GRCh38
        header = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "OR", "L95", "U95", "BETA", "SE", "P", "N"]
        derive = run_standardise(tmp, write_table(tmp / "derive.tsv", header, [
            [1, 1100, "rs_or", "A", "G", 0.3, 1.2, 1.1, 1.3091, None, None, 1e-4, 5000],
            [1, 1200, "rs_ci_beta", "A", "G", 0.3, None, 0.05, 0.15, 0.1, None, 1e-3, 5000],
            [1, 1300, "rs_se_p", "A", "G", 0.3, None, None, None, -0.1, None, 1e-6, 5000],
            [1, 1400, "rs_p", "A", "G", 0.3, None, None, None, 0.1, 0.02, 1e-6, 5000],
            [1, 1500, "rs_p0", "A", "G", 0.3, None, None, None, 0.5, 0.01, 0, 5000],
            [1, 1600, "rs_p_missing", "A", "G", 0.3, None, None, None, 0.5, 0.01, None, 5000],
        ]), reference, input_build="hg38", name="derive")
        self.assertEqual(derive.returncode, 0, derive.stderr)
        rows = {k: v["FORMAT"] for k, v in derive.by_snpid().items()}

        with self.subTest("BETA is log OR for binary traits"):
            self.assertAlmostEqual(float(rows["rs_or"]["ES"]), 0.182321556794, places=7)
            self.assertEqual(derive.qc["derivations"]["beta_from_or"], 1)

        with self.subTest("SE from a confidence interval for odds ratios: (log upper - log lower) / (2 x 1.96)"):
            self.assertAlmostEqual(float(rows["rs_or"]["SE"]), 0.044396147014, places=7)

        with self.subTest("SE from a confidence interval for betas: (upper - lower) / (2 x 1.96)"):
            self.assertAlmostEqual(float(rows["rs_ci_beta"]["SE"]), 0.0255106728462, places=7)
            self.assertEqual(derive.qc["derivations"]["se_from_ci"], 2)

        with self.subTest("SE otherwise from abs(BETA / qnorm(P / 2))"):
            self.assertAlmostEqual(float(rows["rs_se_p"]["SE"]), 0.0204430479678, places=8)
            self.assertEqual(rows["rs_se_p"]["ES"], "-0.1")
            self.assertEqual(derive.qc["derivations"]["se_from_p"], 1)

        with self.subTest("LP is computed from P"):
            self.assertEqual(rows["rs_p"]["LP"], "6")
            self.assertEqual(rows["rs_or"]["LP"], "4")
            lp, info = gwas_stats.compute_lp(p=np.array([0.05, 1e-8, 1.0]))
            np.testing.assert_allclose(lp, [1.30102999566, 8.0, 0.0], rtol=1e-9)
            self.assertEqual(info["from_p"], 3)

        with self.subTest("LP is recovered from BETA/SE when P is 0"):
            self.assertEqual(rows["rs_p0"]["LP"], "544.66531")

        with self.subTest("LP is recovered from BETA/SE when P is missing"):
            self.assertEqual(rows["rs_p_missing"]["LP"], "544.66531")
            self.assertEqual(derive.qc["p_value_source"]["from_z"], 2)

        with self.subTest("LP from BETA/SE keeps its precision beyond double precision, matching R"):
            cases = {1.959963984540054: 1.30102999566, 5.0: 6.24161567673, 36.9: 297.33525943,
                     37.5: 307.035707079, 50.0: 544.665305866, 200.0: 8688.28873886}
            for z, expected, observed in zip(cases, cases.values(), gwas_stats.lp_from_z(np.array(list(cases)))):
                self.assertAlmostEqual(observed, expected, delta=1e-6 * max(1.0, expected), msg=f"z = {z}")
            np.testing.assert_allclose(gwas_stats.lp_from_z(np.array([-5.0, 0.0])), [6.24161567673, 0.0], rtol=1e-9)
            self.assertTrue(np.isnan(gwas_stats.lp_from_z(np.array([np.inf, np.nan]))).all())

        with self.subTest("P outside [0, 1] and a non-positive or missing SE give no LP"):
            lp, info = gwas_stats.compute_lp(p=np.array([1.5, -0.1, 0.0, 0.0]), beta=np.array([0.1, 0.1, 0.1, 0.1]),
                                             se=np.array([0.01, 0.01, 0.0, np.nan]))
            self.assertTrue(np.isnan(lp).all())
            self.assertEqual(info["invalid"], 4)

        # LOG10P column
        header = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "LOG10P", "N"]
        log10p = run_standardise(tmp, write_table(tmp / "log10p.tsv", header, [
            [1, 2100, "rs_450", "A", "G", 0.3, 0.5, 0.01, None, 450, 5000],
            [1, 2200, "rs_12", "A", "G", 0.3, 0.1, 0.02, None, 12.5, 5000],
            [1, 2300, "rs_underflow", "A", "G", 0.3, 0.1, 0.1, 0, 400, 5000],
            [1, 2400, "rs_no_log10p", "A", "G", 0.3, 0.1, 0.02, 0.01, None, 5000],
        ]), reference, input_build="hg38", name="log10p")
        self.assertEqual(log10p.returncode, 0, log10p.stderr)
        rows = {k: v["FORMAT"] for k, v in log10p.by_snpid().items()}

        with self.subTest("LP is read from a LOG10P column, beyond double precision"):
            self.assertEqual((rows["rs_450"]["LP"], rows["rs_12"]["LP"]), ("450", "12.5"))
            self.assertEqual(log10p.qc["source_columns"]["LOG10P"], "LOG10P")

        with self.subTest("A LOG10P value wins over a P column that underflowed to 0"):
            self.assertEqual(rows["rs_underflow"]["LP"], "400")
            self.assertEqual(log10p.qc["p_value_source"]["from_log10p"], 3)

        with self.subTest("Rows without a LOG10P value fall back to P"):
            self.assertEqual(rows["rs_no_log10p"]["LP"], "2")
            self.assertEqual(log10p.qc["p_value_source"]["from_p"], 1)

        with self.subTest("A LOG10P column holding log10 P (negative values) is negated"):
            lp, info = gwas_stats.compute_lp(log10p=np.array([-2.0, -8.0, -0.5]))
            np.testing.assert_allclose(lp, [2.0, 8.0, 0.5])
            self.assertTrue(info["log10p_negated"])

        with self.subTest("The OpenGWAS LP 999 sentinel (P = 0) is recovered from BETA/SE"):
            vcf = tmp / "p0.gwas.vcf"
            vcf.write_text("\n".join([
                "##fileformat=VCFv4.2", "##contig=<ID=1,length=248956422,assembly=GRCh38>",
                *[f'##FORMAT=<ID={k},Number=A,Type=Float,Description="{k}">' for k in ("ES", "SE", "LP", "AF", "SS")],
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGWAS",
                "1\t3100\trs1\tG\tA\t.\tPASS\t.\tES:SE:LP:AF:SS\t0.5:0.01:999:0.3:1000",
                "1\t3200\trs2\tG\tA\t.\tPASS\t.\tES:SE:LP:AF:SS\t0.1:0.02:6:0.3:1000",
            ]) + "\n", encoding="utf-8")
            result = run_standardise(tmp, vcf, reference, input_build="hg38", name="sentinel")
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = result.by_snpid()
            self.assertAlmostEqual(float(rows["rs1"]["FORMAT"]["LP"]), 544.665305866, places=3)  # z = 50
            self.assertEqual(rows["rs2"]["FORMAT"]["LP"], "6")
            self.assertEqual(result.qc["p_value_source"]["from_z"], 1)


if __name__ == "__main__":
    unittest.main()
