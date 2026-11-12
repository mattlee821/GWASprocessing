"""Tests for gwas_stats: -log10 P derivation and genomic inflation.

Expected values were computed independently in R:
    lp <- function(z) -(log(2) + pnorm(-abs(z), log.p = TRUE)) / log(10)
    qchisq(0.95, 1) / qchisq(0.5, 1)
"""

import unittest

import numpy as np

import testlib  # noqa: F401  (puts pipeline/bin on sys.path)

assert testlib.BIN.is_dir()

import gwas_stats


class LpFromZTests(unittest.TestCase):
    def test_matches_r_across_the_double_precision_boundary(self):
        cases = {1.959963984540054: 1.30102999566, 5.0: 6.24161567673, 36.9: 297.33525943,
                 37.5: 307.035707079, 50.0: 544.665305866, 200.0: 8688.28873886}
        z = np.array(list(cases))
        got = gwas_stats.lp_from_z(z)
        for value, expected, observed in zip(cases, cases.values(), got):
            with self.subTest(z=value):
                self.assertAlmostEqual(observed, expected, delta=1e-6 * max(1.0, expected))

    def test_sign_of_z_does_not_matter_and_zero_gives_zero(self):
        np.testing.assert_allclose(gwas_stats.lp_from_z(np.array([-5.0, 0.0])), [6.24161567673, 0.0], rtol=1e-9)

    def test_non_finite_z_gives_nan(self):
        self.assertTrue(np.isnan(gwas_stats.lp_from_z(np.array([np.inf, np.nan]))).all())


class ComputeLpTests(unittest.TestCase):
    def test_p_values_convert_directly(self):
        lp, info = gwas_stats.compute_lp(p=np.array([0.05, 1e-8, 1.0]))
        np.testing.assert_allclose(lp, [1.30102999566, 8.0, 0.0], rtol=1e-9)
        self.assertEqual(info["from_p"], 3)

    def test_p_of_zero_is_recovered_from_beta_and_se(self):
        lp, info = gwas_stats.compute_lp(p=np.array([0.0]), beta=np.array([0.5]), se=np.array([0.01]))
        self.assertAlmostEqual(lp[0], 544.665305866, delta=1e-3)
        self.assertEqual(info["from_z"], 1)

    def test_log10p_column_is_used_as_is_even_beyond_double_precision(self):
        lp, info = gwas_stats.compute_lp(p=np.array([np.nan, np.nan]), log10p=np.array([450.0, 12.5]))
        np.testing.assert_allclose(lp, [450.0, 12.5])
        self.assertEqual(info["from_log10p"], 2)

    def test_log10p_wins_over_an_underflowed_p_column(self):
        lp, _ = gwas_stats.compute_lp(p=np.array([0.0]), log10p=np.array([400.0]), beta=np.array([0.1]), se=np.array([0.1]))
        self.assertEqual(lp[0], 400.0)

    def test_negative_log10p_column_is_treated_as_log10_p(self):
        lp, info = gwas_stats.compute_lp(log10p=np.array([-2.0, -8.0, -0.5]))
        np.testing.assert_allclose(lp, [2.0, 8.0, 0.5])
        self.assertTrue(info["log10p_negated"])

    def test_rows_missing_log10p_fall_back_to_p(self):
        lp, _ = gwas_stats.compute_lp(p=np.array([0.01, 0.1]), log10p=np.array([np.nan, 3.0]))
        np.testing.assert_allclose(lp, [2.0, 3.0])

    def test_out_of_range_p_and_bad_se_are_invalid(self):
        lp, info = gwas_stats.compute_lp(p=np.array([1.5, -0.1, 0.0, 0.0]), beta=np.array([0.1, 0.1, 0.1, 0.1]),
                                         se=np.array([0.01, 0.01, 0.0, np.nan]))
        self.assertTrue(np.isnan(lp).all())
        self.assertEqual(info["invalid"], 4)


class LambdaGcTests(unittest.TestCase):
    def test_uniform_median_gives_one(self):
        lp = -np.log10(np.array([0.1, 0.5, 0.9]))
        self.assertAlmostEqual(gwas_stats.lambda_gc(lp), 1.0, places=9)

    def test_inflated_median_matches_chi_square_ratio(self):
        lp = -np.log10(np.array([0.01, 0.05, 0.2]))
        self.assertAlmostEqual(gwas_stats.lambda_gc(lp), 8.44394650653, places=6)

    def test_no_finite_values_gives_none(self):
        self.assertIsNone(gwas_stats.lambda_gc(np.array([np.nan])))


if __name__ == "__main__":
    unittest.main()
