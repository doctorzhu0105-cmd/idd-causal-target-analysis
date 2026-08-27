from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis_p0_correlated_mr_20260814" / "scripts"))

from correlated_mr import (  # noqa: E402
    _omega_inverse,
    egger_correlated,
    ivw_correlated,
    leave_one_out,
    orient_to_exposure_increasing,
)


class Master14CorrelatedMRTests(unittest.TestCase):
    def setUp(self) -> None:
        self.bx = np.array([0.22, -0.16, 0.31, -0.27, 0.19])
        self.by = np.array([0.020, 0.018, 0.041, -0.006, 0.029])
        self.se_y = np.array([0.014, 0.018, 0.016, 0.015, 0.017])
        self.rho = np.array(
            [
                [1.00, 0.18, -0.08, 0.04, 0.02],
                [0.18, 1.00, 0.12, -0.06, 0.03],
                [-0.08, 0.12, 1.00, 0.16, -0.05],
                [0.04, -0.06, 0.16, 1.00, 0.10],
                [0.02, 0.03, -0.05, 0.10, 1.00],
            ]
        )

    def test_orientation_changes_egger_intercept_but_not_ivw(self) -> None:
        omega_inv, _ = _omega_inverse(self.se_y, self.rho)
        design = np.column_stack([np.ones(len(self.bx)), self.bx])
        unoriented = np.linalg.pinv(design.T @ omega_inv @ design) @ design.T @ omega_inv @ self.by
        _, oriented_intercept = egger_correlated(self.bx, self.by, self.se_y, self.rho)
        self.assertNotAlmostEqual(float(unoriented[0]), oriented_intercept.beta, places=10)

        bx, by, rho = orient_to_exposure_increasing(self.bx, self.by, self.rho)
        before = ivw_correlated(self.bx, self.by, self.se_y, self.rho)
        after = ivw_correlated(bx, by, self.se_y, rho)
        self.assertAlmostEqual(before.beta, after.beta, places=14)
        self.assertAlmostEqual(before.se, after.se, places=14)
        self.assertAlmostEqual(before.pvalue, after.pvalue, places=14)

    def test_generalized_q_and_k_minus_one_df(self) -> None:
        result = ivw_correlated(self.bx, self.by, self.se_y, self.rho)
        omega_inv, _ = _omega_inverse(self.se_y, self.rho)
        residual = self.by - result.beta * self.bx
        expected_q = float(residual @ omega_inv @ residual)
        self.assertAlmostEqual(result.q_stat, expected_q, places=13)
        self.assertAlmostEqual(result.overdispersion, max(1.0, expected_q / (len(self.bx) - 1)), places=13)

    def test_leave_one_out_uses_reduced_omega_without_repruning(self) -> None:
        rows = leave_one_out(self.bx, self.by, self.se_y, self.rho)
        self.assertEqual(len(rows), len(self.bx))
        for index, row in enumerate(rows):
            keep = np.arange(len(self.bx)) != index
            direct = ivw_correlated(
                self.bx[keep], self.by[keep], self.se_y[keep], self.rho[np.ix_(keep, keep)]
            )
            self.assertEqual(row["dropped_index"], index)
            self.assertAlmostEqual(row["beta"], direct.beta, places=14)
            self.assertAlmostEqual(row["se"], direct.se, places=14)

    def test_pseudoinverse_handles_singular_omega(self) -> None:
        rho = np.ones((3, 3), dtype=float)
        bx = np.array([0.1, 0.2, 0.3])
        by = np.array([0.02, 0.04, 0.06])
        se_y = np.array([0.01, 0.01, 0.01])
        result = ivw_correlated(bx, by, se_y, rho, random_effects=False)
        self.assertTrue(np.isfinite(result.beta))
        self.assertAlmostEqual(result.beta, 0.2, places=12)
        self.assertGreater(result.condition_number, 1e12)


if __name__ == "__main__":
    unittest.main()
