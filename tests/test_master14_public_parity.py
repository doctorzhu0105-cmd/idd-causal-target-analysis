from __future__ import annotations

import csv
import hashlib
import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "py"))

from generate_master14_assets import validate_table  # noqa: E402
from master14_sensitivity import bootstrap_ld, effective_k, multiplicity_sensitivity, shrink_ld  # noqa: E402


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Master14PublicParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = ROOT / "results" / "master_1_4"
        cls.version = json.loads((ROOT / "provenance" / "ANALYSIS_VERSION.json").read_text(encoding="utf-8"))

    def test_master_identity_is_locked_1_4_0(self) -> None:
        self.assertEqual(self.version["master_version"], "1.4.0")
        self.assertEqual(
            self.version["master_hash"],
            "d5bb9062efbd170abd7b9380673a423e5cc54bd94de6f50c48092b503a512006",
        )
        self.assertEqual(self.version["provenance_unresolved"], 0)

    def test_multiplicity_expected_counts(self) -> None:
        data = pd.read_csv(self.results / "multiplicity_summary.tsv", sep="\t")
        global_row = data.loc[data.family.eq("global_all_estimable")].iloc[0]
        self.assertEqual(int(global_row.number_of_tests), 7284)
        self.assertEqual(int(global_row.BH_pass_0_05), 26)
        self.assertEqual(int(global_row.BY_pass_0_05), 5)
        k2 = data.loc[data.family.eq("K>=2")].iloc[0]
        self.assertEqual((int(k2.number_of_tests), int(k2.BH_pass_0_05)), (6210, 13))

    def test_p2ry13_is_not_rescued_by_current_families(self) -> None:
        row = pd.read_csv(self.results / "p2ry13_multiplicity_path.tsv", sep="\t").iloc[0]
        self.assertEqual(row.gene, "P2RY13")
        self.assertEqual(int(row.K), 21)
        self.assertLess(float(row.historical_15_row_q), 0.05)
        for column in ("global_BH", "global_BY", "endpoint_BH", "endpoint_BY", "K_ge_2_BH", "K_ge_3_BH", "K_ge_4_BH"):
            self.assertGreater(float(row[column]), 0.05)

    def test_public_sensitivity_primitives_are_deterministic(self) -> None:
        frame = pd.DataFrame({
            "pvalue": [0.01, 0.02, 0.40, 0.90],
            "outcome": ["A", "A", "B", "B"],
            "k": [1, 2, 3, 4],
        })
        tables = multiplicity_sensitivity(frame)
        self.assertEqual(len(tables["global"]), 4)
        self.assertEqual(len(tables["endpoint"]), 4)
        self.assertEqual(set(tables["k_family"].family_definition), {"all_estimable", "K>=2", "K>=3", "K>=4"})
        self.assertAlmostEqual(effective_k(np.eye(3)), 3.0)

        dosages = np.array([[0.0, 1.0, 2.0, 1.0], [2.0, 1.0, 0.0, 1.0]])
        first = bootstrap_ld(dosages, replicates=3, seed=17)
        second = bootstrap_ld(dosages, replicates=3, seed=17)
        for left, right in zip(first, second):
            np.testing.assert_allclose(left, right)
        np.testing.assert_allclose(shrink_ld(np.eye(2), 0.5), np.eye(2))

    def test_coloc_prior_registered_rows_and_range(self) -> None:
        data = read_tsv(self.results / "coloc_prior_sensitivity_summary.tsv")
        self.assertTrue(data.p12.astype(float).between(1e-7, 1e-4).all())
        registered = data.loc[data.registered_p12.astype(str).str.lower().eq("true")]
        self.assertTrue((registered.p12.astype(float) == 1e-5).all())
        self.assertGreater(len(registered), 0)

    def test_bootstrap_identity_matches_public_provenance(self) -> None:
        expected = {}
        with (ROOT / "provenance" / "OUTPUT_HASHES.tsv").open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                expected[row["public_path"]] = row["sha256"]
        relative = "results/master_1_4/ld_reference_bootstrap_summary.tsv"
        self.assertEqual(sha256(ROOT / relative), expected[relative])
        data = read_tsv(ROOT / relative)
        self.assertEqual(len(data), 10)
        self.assertEqual(int(data.descriptive_material_instability_flag.astype(str).str.lower().eq("true").sum()), 0)

    def test_current_evidence_matrix_has_no_tier_or_grade_encoding(self) -> None:
        table = pd.read_csv(self.results / "table3_evidence_matrix.tsv", sep="\t", dtype=str)
        validate_table(table)
        self.assertFalse(any("tier" in column.lower() or "grade" in column.lower() for column in table.columns))
        self.assertFalse(table.astype(str).apply(lambda column: column.str.contains(r"\bTier [ABCD]\b", regex=True).any()).any())

    def test_current_figure1_label_and_current_docs(self) -> None:
        generator = (ROOT / "analysis_p0_correlated_mr_20260814" / "scripts" / "make_current_figure1.py").read_text(encoding="utf-8")
        self.assertIn("historical candidates", generator)
        self.assertNotIn("formal candidates", generator)
        current_docs = [
            ROOT / "README.md", ROOT / "README.zh-CN.md",
            ROOT / "docs/REPRODUCIBILITY.md", ROOT / "provenance/REPRODUCIBILITY_SCOPE.md",
        ]
        prohibited = ("four novel targets", "validated targets", "primary 13-row discovery family")
        for path in current_docs:
            text = path.read_text(encoding="utf-8").lower()
            for phrase in prohibited:
                self.assertNotIn(phrase.lower(), text, f"{phrase!r} in {path}")


if __name__ == "__main__":
    unittest.main()
