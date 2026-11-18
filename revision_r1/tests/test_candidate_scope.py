from pathlib import Path
import csv
import json
import unittest

ROOT = Path(__file__).resolve().parents[2]

class CandidateScopeTests(unittest.TestCase):
    def test_version_metadata_and_dated_preparation_observation(self):
        status = json.loads((ROOT / "RELEASE_STATUS.json").read_text(encoding="utf-8"))
        self.assertEqual((status["version"], status["update_type"]), ("1.3.1", "maintenance"))
        snapshot = status["preparation_snapshot"]
        self.assertEqual(snapshot["as_of"], "2026-11-18")
        for field in ("git_tag", "github_release", "zenodo_version_doi", "published_at"):
            self.assertIsNone(snapshot[field])  # None had been assigned at preparation.
        records = status["publication_records"]
        self.assertEqual(records["github_releases"], "https://github.com/doctorzhu0105-cmd/idd-causal-target-analysis/releases")
        previous = status["historical_releases"]
        r1 = previous["v1.3.0"]
        self.assertEqual(r1["git_tag"], "v1.3.0")
        self.assertEqual(r1["commit"], "6ada8319292387c63e27346d22ffd847708ef71e")
        self.assertEqual(r1["published_at"], "2026-11-17T13:59:36Z")
        self.assertEqual(r1["zenodo_version_doi"], "10.5281/zenodo.24312891")
        self.assertEqual(records["r1_zenodo_series"], "https://doi.org/" + r1["zenodo_concept_doi"])
        self.assertEqual(previous["v1.2.0"]["zenodo_concept_doi"], "10.5281/zenodo.21486870")
        self.assertNotEqual(r1["zenodo_concept_doi"], previous["v1.2.0"]["zenodo_concept_doi"])

    def test_conditional_evidence_tiers_and_separate_populations(self):
        p = ROOT / "revision_r1/results/K1_CONDITIONAL_EVIDENCE_COUNTS.tsv"
        with p.open(encoding="utf-8", newline="") as f:
            rows = {r["measure"]: r for r in csv.DictReader(f, delimiter="\t")}
        n = lambda k: int(rows[k]["n"])
        self.assertEqual(n("historically_saved_matrix_complete") + n("conditional_singleton_reconstruction"), 7155)
        self.assertEqual(n("combined_conditional_evaluable") + n("incomplete"), 7284)
        self.assertEqual(n("empty_common_set") + n("nonempty_high_threshold_cache_missing") + n("nonempty_low_threshold_cache_missing"), 129)
        self.assertIn("Separate", rows["missing_membership_configurations"]["boundary"])
        self.assertIn("6087", rows["RE_direction_changes"]["boundary"])
        self.assertIn("6087", rows["FE_direction_changes"]["boundary"])

    def test_verifiers_present_and_restricted_rows_absent(self):
        for name in ("verify_ab.py", "verify_rho.py", "verify_r37.py", "verify_psd.py"):
            self.assertTrue((ROOT / "revision_r1/code/verification" / name).is_file())
            cfg = ROOT / "revision_r1/config" / Path(name).with_suffix(".json")
            self.assertTrue(cfg.is_file())
            self.assertIsInstance(json.loads(cfg.read_text(encoding="utf-8")), dict)
        self.assertFalse((ROOT / "revision_r1/code/k1_conditional_forensic_generator.py").exists())
        self.assertFalse(any(ROOT.glob("revision_r1/**/*1068*.tsv")))

if __name__ == "__main__":
    unittest.main()
