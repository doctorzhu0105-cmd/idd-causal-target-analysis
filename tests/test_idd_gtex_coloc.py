import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py"))

from idd_gtex_coloc import run_coloc, split_rsids  # noqa: E402


class GtexColocTests(unittest.TestCase):
    def test_split_rsids_accepts_registered_delimiters(self):
        self.assertEqual(split_rsids("rs1,rs2;rs3|rs4"), {"rs1", "rs2", "rs3", "rs4"})

    def test_coloc_requires_minimum_shared_variants(self):
        eqtl = [{"rsid": "rs1", "beta": 0.1, "se": 0.02, "pvalue": 1e-6}]
        outcome = [{"rsid": "rs1", "beta": 0.05, "se": 0.01, "pvalue": 1e-5}]
        result = run_coloc(eqtl, outcome, 50)
        self.assertEqual(result["status"], "not_estimable")
        self.assertEqual(result["shared_variants"], 1)

    def test_coloc_posteriors_are_emitted_for_sufficient_overlap(self):
        eqtl = []
        outcome = []
        for index in range(50):
            eqtl.append(
                {
                    "rsid": f"rs{index}",
                    "beta": 0.2 if index == 0 else 0.001,
                    "se": 0.02,
                    "pvalue": 1e-12 if index == 0 else 0.9,
                }
            )
            outcome.append(
                {
                    "rsid": f"rs{index}",
                    "beta": 0.15 if index == 0 else 0.001,
                    "se": 0.02,
                    "pvalue": 1e-10 if index == 0 else 0.9,
                }
            )
        result = run_coloc(eqtl, outcome, 50)
        self.assertEqual(result["status"], "estimable")
        self.assertAlmostEqual(sum(result[f"pp{i}"] for i in range(5)), 1.0)


if __name__ == "__main__":
    unittest.main()
