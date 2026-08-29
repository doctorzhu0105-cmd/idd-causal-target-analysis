from __future__ import annotations

import csv
import hashlib
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "master_1_5"
SCRIPTS = ROOT / "analysis_p0_correlated_mr_20260814" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from catalogue_compute import effective_k, joint_information  # noqa: E402
from catalogue_extension import endpoint_associations, exposure_pairs, matched_associations, precision_pairs  # noqa: E402


def read_tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / name, sep="\t", comment="#")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Master15CatalogueExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.effective = read_tsv("CATALOGUE_EFFECTIVE_K.tsv")
        cls.joint = read_tsv("CATALOGUE_JOINT_EXPOSURE_INFORMATION.tsv")
        cls.givw = read_tsv("CATALOGUE_GIVW_LD_LADDER.tsv")
        cls.exposure = read_tsv("CATALOGUE_EXPOSURE_PAIRED_CHANGES.tsv")
        cls.precision = read_tsv("CATALOGUE_PRECISION_PAIRED_CHANGES.tsv")

    def test_public_output_bytes_match_governed_provenance(self) -> None:
        with (ROOT / "provenance" / "OUTPUT_HASHES.tsv").open(encoding="utf-8", newline="") as handle:
            rows = [row for row in csv.DictReader(handle, delimiter="\t")
                    if row["public_path"].startswith("results/master_1_5/")]
        self.assertEqual(len(rows), 17)
        for row in rows:
            path = ROOT / row["public_path"]
            self.assertEqual(path.stat().st_size, int(row["bytes"]))
            self.assertEqual(sha256(path), row["sha256"])
            self.assertEqual(row["sha256"], row["upstream_sha256"])
            self.assertEqual(row["relationship"], "byte_identical_copy")

    def test_catalogue_row_complete_identity(self) -> None:
        self.assertEqual(len(self.effective), 2479 * 5)
        self.assertEqual(len(self.joint), 2479 * 5)
        self.assertEqual(len(self.givw), 2479 * 3 * 5)
        self.assertEqual(int(self.effective.status.eq("estimable").sum()), 2438 * 5)
        self.assertEqual(int(self.givw.status.eq("estimable").sum()), 2428 * 3 * 5)

    def test_effective_k_and_joint_information_primitives(self) -> None:
        rho = np.array([[1.0, 0.25], [0.25, 1.0]])
        beta = np.array([0.2, -0.1])
        se = np.array([0.05, 0.08])
        self.assertAlmostEqual(effective_k(np.eye(3)), 3.0)
        expected = float(beta @ np.linalg.pinv(np.outer(se, se) * rho) @ beta)
        self.assertAlmostEqual(joint_information(beta, se, rho), expected)

    def test_joint_information_identity(self) -> None:
        estimable = self.joint[self.joint.status.eq("estimable")]
        error = (estimable.conditional_strength_analogue * estimable.K - estimable.joint_exposure_chi_square).abs()
        self.assertLessEqual(float(error.max()), 1e-11)

    def test_primary_pair_regeneration(self) -> None:
        regenerated_exposure = exposure_pairs(self.effective, self.joint)
        regenerated_precision = precision_pairs(self.givw)
        self.assertEqual((len(regenerated_exposure), len(regenerated_precision)), (4876, 14568))
        expected_e = self.exposure.sort_values(["gene", "contrast"]).reset_index(drop=True)
        actual_e = regenerated_exposure.sort_values(["gene", "contrast"]).reset_index(drop=True)
        expected_p = self.precision.sort_values(["gene", "endpoint", "contrast"]).reset_index(drop=True)
        actual_p = regenerated_precision.sort_values(["gene", "endpoint", "contrast"]).reset_index(drop=True)
        np.testing.assert_allclose(actual_e.select_dtypes("number"), expected_e.select_dtypes("number"), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(actual_p.select_dtypes("number"), expected_p.select_dtypes("number"), rtol=1e-12, atol=1e-12)

    def test_historical_75_row_ladder_is_unchanged(self) -> None:
        old_e = pd.read_csv(ROOT / "results/master_1_4/effective_instrument_count.tsv", sep="\t")
        old_j = pd.read_csv(ROOT / "results/master_1_4/conditional_strength.tsv", sep="\t")
        genes = set(old_e.gene)
        new = self.givw[self.givw.gene.isin(genes)].copy()
        merged = new.merge(old_e, left_on=["gene", "ensembl_gene_id", "endpoint", "threshold"],
                           right_on=["gene", "ensembl_gene_id", "outcome", "r2_threshold"], validate="one_to_one")
        merged = merged.merge(old_j, left_on=["gene", "endpoint", "threshold"],
                              right_on=["gene", "outcome", "r2_threshold"], suffixes=("", "_joint"), validate="one_to_one")
        self.assertEqual(len(merged), 75)
        self.assertTrue((merged.K == merged.k).all())
        self.assertTrue((merged.K == merged.k_joint).all())
        np.testing.assert_allclose(merged.effective_K, merged.effective_independent_k, rtol=0, atol=1e-12)
        np.testing.assert_allclose(merged.SE, merged.se, rtol=0, atol=1e-12)
        np.testing.assert_allclose(merged.P, merged.pvalue, rtol=0, atol=1e-12)
        np.testing.assert_allclose(merged.joint_exposure_chi_square, merged.joint_chi2, rtol=0, atol=1e-12)

    def test_matched_and_endpoint_associations(self) -> None:
        matched = matched_associations(self.precision)
        expected = read_tsv("INFORMATION_PRECISION_MATCHED_ASSOCIATION.tsv")
        np.testing.assert_allclose(matched[["spearman_rho", "P"]], expected[["spearman_rho", "P"]], rtol=1e-12, atol=1e-15)
        endpoints = endpoint_associations(self.precision)
        expected_endpoints = read_tsv("ENDPOINT_STRATIFIED_INFORMATION_PRECISION.tsv")
        np.testing.assert_allclose(
            endpoints[["nominal_K_relative_rho", "effective_K_relative_rho", "joint_information_relative_rho"]],
            expected_endpoints[["nominal_K_relative_rho", "effective_K_relative_rho", "joint_information_relative_rho"]],
            rtol=1e-12, atol=1e-15,
        )
        expected_values = {
            "M13_INTERVERTEB": (0.12281333740090211, 0.11496901995825172, 0.4024438027966377),
            "M13_LOWBACKPAIN": (0.18403059886219386, 0.17747954384627393, 0.4483619486697231),
            "M13_SCIATICA": (0.22784014421617463, 0.22294165628973073, 0.4669617130564762),
        }
        for endpoint, values in expected_values.items():
            row = endpoints[endpoints.endpoint.eq(endpoint)].iloc[0]
            np.testing.assert_allclose(row[["nominal_K_relative_rho", "effective_K_relative_rho", "joint_information_relative_rho"]].astype(float), values, rtol=0, atol=1e-12)

    def test_k1_second_order_and_rank_propagation(self) -> None:
        data = read_tsv("K1_SECOND_ORDER_MULTIPLICITY_REPLICATION.tsv")
        self.assertEqual((len(data), int(data.k.eq(1).sum())), (7284, 1074))
        self.assertEqual((int(data.original_global_BH_pass.sum()), int(data.second_order_global_BH_pass.sum())), (26, 18))
        self.assertEqual((int(data.original_global_BY_pass.sum()), int(data.second_order_global_BY_pass.sum())), (5, 2))
        losses = sorted(data.loc[(data.k >= 2) & data.original_global_BH_pass & ~data.second_order_global_BH_pass, "gene"])
        self.assertEqual(losses, ["BMP6", "FBLN2", "LYG1", "RIOK3"])

    def test_figure1_zero_gain_rows_are_retained(self) -> None:
        exposure = self.exposure[self.exposure.contrast.eq("0.001_to_0.1")]
        precision = self.precision[self.precision.contrast.eq("0.001_to_0.1")]
        self.assertEqual((len(exposure), int(exposure.K_percent_change.eq(0).sum()),
                          int(exposure.joint_chi_square_percent_change.eq(0).sum())), (2438, 392, 377))
        self.assertEqual((len(precision), int(precision.joint_chi_square_percent_change.eq(0).sum()),
                          int(precision.precision_gain_percent.eq(0).sum())), (7284, 1125, 1125))

    def test_no_restricted_raw_input_is_bundled(self) -> None:
        prohibited = {".bed", ".bim", ".fam", ".gz", ".pkl", ".rds", ".h5ad"}
        files = [path for path in ROOT.rglob("*") if path.is_file() and ".git" not in path.parts]
        self.assertFalse([path for path in files if path.suffix.lower() in prohibited])
        headers = "\n".join(path.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
                            for path in RESULTS.glob("*.tsv"))
        for field in ("variant_id", "rsid", "effect_allele", "other_allele", "eaf"):
            self.assertNotIn(field, headers.lower())


if __name__ == "__main__":
    unittest.main()
