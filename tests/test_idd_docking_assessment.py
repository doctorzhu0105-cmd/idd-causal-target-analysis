import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py"))

from idd_docking_assessment import assess_pair  # noqa: E402


class DockingAssessmentTests(unittest.TestCase):
    def test_pair_requires_structure_smiles_and_aligned_direction(self):
        priority = {
            "gene_symbol": "GENE",
            "top_modulator_smiles": "CCO",
            "top_modulator_direction_aligned": "True",
        }
        self.assertTrue(assess_pair(priority, {"pdb_ids": ["1ABC"]})["docking_applicable"])
        self.assertFalse(assess_pair(priority, {"pdb_ids": []})["docking_applicable"])
        self.assertFalse(
            assess_pair({**priority, "top_modulator_smiles": ""}, {"pdb_ids": ["1ABC"]})[
                "docking_applicable"
            ]
        )
        self.assertFalse(
            assess_pair(
                {**priority, "top_modulator_direction_aligned": "False"},
                {"pdb_ids": ["1ABC"]},
            )["docking_applicable"]
        )


if __name__ == "__main__":
    unittest.main()
