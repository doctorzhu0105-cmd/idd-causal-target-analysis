import math
import sys
import unittest
from pathlib import Path

from scipy.stats import t as student_t


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py"))

from idd_genetics_utils import (  # noqa: E402
    benjamini_hochberg,
    coloc_posteriors,
    effective_binary_n,
    harmonize_effects,
    ivw_fixed,
    ivw_multiplicative_random_effects,
    mr_egger,
    observed_r2_from_z,
    standardized_effect_from_z,
    wald_ratio,
    wakefield_log_abf,
    weighted_median,
)


class HarmonizationTests(unittest.TestCase):
    def test_swapped_outcome_is_flipped(self):
        exposure = {
            "effect_allele": "A",
            "other_allele": "C",
            "beta": 0.2,
            "se": 0.05,
            "eaf": 0.2,
        }
        outcome = {
            "effect_allele": "C",
            "other_allele": "A",
            "beta": 0.1,
            "se": 0.02,
            "eaf": 0.8,
        }
        aligned, reason = harmonize_effects(exposure, outcome)
        self.assertEqual(reason, "ok")
        self.assertEqual(aligned.action, "swapped")
        self.assertAlmostEqual(aligned.beta_outcome, -0.1)

    def test_ambiguous_palindrome_is_removed(self):
        exposure = {
            "effect_allele": "A",
            "other_allele": "T",
            "beta": 0.2,
            "se": 0.05,
            "eaf": 0.49,
        }
        outcome = {
            "effect_allele": "A",
            "other_allele": "T",
            "beta": 0.1,
            "se": 0.02,
            "eaf": 0.48,
        }
        aligned, reason = harmonize_effects(exposure, outcome)
        self.assertIsNone(aligned)
        self.assertEqual(reason, "palindromic_ambiguous")


class StatisticsTests(unittest.TestCase):
    def test_standardized_effect_reconstructs_z(self):
        beta, se = standardized_effect_from_z(5.0, 0.25, 30000)
        self.assertAlmostEqual(beta / se, 5.0)

    def test_wald_ratio(self):
        result = wald_ratio(0.2, 0.01, 0.1, 0.02)
        self.assertAlmostEqual(result["beta"], 0.5)
        self.assertGreater(result["se"], 0)

    def test_large_z_pvalue_does_not_round_to_zero(self):
        result = wald_ratio(1.0, 0.001, 1.0, 0.1)
        self.assertGreater(result["pvalue"], 0.0)
        self.assertLess(result["pvalue"], 1e-20)

    def test_ivw_fixed(self):
        result = ivw_fixed(
            [
                {"beta_exposure": 0.2, "beta_outcome": 0.1, "se_outcome": 0.02},
                {"beta_exposure": 0.4, "beta_outcome": 0.2, "se_outcome": 0.04},
            ]
        )
        self.assertAlmostEqual(result["beta"], 0.5)
        self.assertAlmostEqual(result["q"], 0.0)

    def test_multiplicative_random_effects_inflates_se_under_heterogeneity(self):
        rows = [
            {"beta_exposure": 0.2, "beta_outcome": 0.1, "se_outcome": 0.01},
            {"beta_exposure": 0.4, "beta_outcome": -0.1, "se_outcome": 0.01},
            {"beta_exposure": 0.3, "beta_outcome": 0.2, "se_outcome": 0.01},
        ]
        fixed = ivw_fixed(rows)
        random_effects = ivw_multiplicative_random_effects(rows)
        self.assertGreater(random_effects["dispersion"], 1.0)
        self.assertGreater(random_effects["se"], fixed["se"])

    def test_r2_and_effective_sample_size(self):
        self.assertAlmostEqual(effective_binary_n(100, 300), 300.0)
        self.assertAlmostEqual(observed_r2_from_z(2.0, 100), 4.0 / 102.0)

    def test_mr_egger_is_invariant_to_joint_allele_recoding(self):
        rows = [
            {"beta_exposure": 0.10, "beta_outcome": 0.04, "se_outcome": 0.01},
            {"beta_exposure": 0.20, "beta_outcome": 0.09, "se_outcome": 0.02},
            {"beta_exposure": 0.35, "beta_outcome": 0.13, "se_outcome": 0.02},
        ]
        recoded = [dict(row) for row in rows]
        recoded[1]["beta_exposure"] *= -1
        recoded[1]["beta_outcome"] *= -1
        original = mr_egger(rows)
        changed = mr_egger(recoded)
        self.assertAlmostEqual(original["slope_beta"], changed["slope_beta"])
        self.assertAlmostEqual(original["intercept"], changed["intercept"])

    def test_mr_egger_uses_k_minus_2_student_t_inference(self):
        rows = [
            {"beta_exposure": 0.10, "beta_outcome": 0.04, "se_outcome": 0.01},
            {"beta_exposure": 0.20, "beta_outcome": 0.09, "se_outcome": 0.02},
            {"beta_exposure": 0.35, "beta_outcome": 0.13, "se_outcome": 0.02},
        ]
        result = mr_egger(rows)
        slope_t = result["slope_beta"] / result["slope_se"]
        intercept_t = result["intercept"] / result["intercept_se"]
        self.assertAlmostEqual(result["slope_pvalue"], 2.0 * student_t.sf(abs(slope_t), 1))
        self.assertAlmostEqual(result["intercept_pvalue"], 2.0 * student_t.sf(abs(intercept_t), 1))

    def test_coloc_posteriors_sum_to_one(self):
        trait1 = [wakefield_log_abf(0.2, 0.01, 0.0225), wakefield_log_abf(0.01, 0.01, 0.0225)]
        trait2 = [wakefield_log_abf(0.25, 0.01, 0.04), wakefield_log_abf(0.01, 0.01, 0.04)]
        posterior = coloc_posteriors(trait1, trait2)
        self.assertAlmostEqual(sum(posterior.values()), 1.0)
        self.assertTrue(all(0.0 <= value <= 1.0 for value in posterior.values()))

    def test_weighted_median_and_egger_recover_linear_effect(self):
        rows = [
            {
                "beta_exposure": bx,
                "se_exposure": 0.005,
                "beta_outcome": 0.5 * bx,
                "se_outcome": 0.01,
            }
            for bx in (0.1, 0.2, 0.3, 0.4)
        ]
        median = weighted_median(rows, seed=1, bootstrap_iterations=200)
        egger = mr_egger(rows)
        self.assertAlmostEqual(median["beta"], 0.5)
        self.assertAlmostEqual(egger["slope_beta"], 0.5)
        self.assertAlmostEqual(egger["intercept"], 0.0)

    def test_benjamini_hochberg_preserves_missing(self):
        adjusted = benjamini_hochberg([0.01, 0.02, 0.2, None])
        self.assertAlmostEqual(adjusted[0], 0.03)
        self.assertAlmostEqual(adjusted[1], 0.03)
        self.assertAlmostEqual(adjusted[2], 0.2)
        self.assertIsNone(adjusted[3])


if __name__ == "__main__":
    unittest.main()
