from pathlib import Path
import csv
import json
import unittest

ROOT = Path(__file__).resolve().parents[2]

class CandidateScopeTests(unittest.TestCase):
    def test_release_status_is_not_published(self):
        status = json.loads((ROOT / "RELEASE_STATUS.json").read_text(encoding="utf-8"))
        self.assertEqual(status["target_version"], "1.3.1")
        self.assertEqual(status["release_state"], "RELEASE_CANDIDATE")
        self.assertEqual(status["publication_status"], "NOT_YET_PUBLISHED")
        self.assertIsNone(status["git_tag"])
        self.assertIsNone(status["zenodo_version_doi"])
        self.assertEqual(status["last_verified_public_release"]["git_tag"], "v1.3.0")

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
