from __future__ import annotations

import csv
import hashlib
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from build_manifest import canonical_bytes  # noqa: E402


class ReleaseManifestV120Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tsv_path = ROOT / "provenance/RELEASE_MANIFEST_v1.2.0.tsv"
        cls.json_path = ROOT / "provenance/RELEASE_MANIFEST_v1.2.0.json"
        with cls.tsv_path.open(encoding="utf-8", newline="") as handle:
            cls.tsv = list(csv.DictReader(handle, delimiter="\t"))
        cls.payload = json.loads(cls.json_path.read_text(encoding="utf-8"))

    def test_json_and_tsv_are_identical(self) -> None:
        expected = [{key: str(value) for key, value in row.items()} for row in self.payload["files"]]
        self.assertEqual(self.tsv, expected)

    def test_membership_and_hashes(self) -> None:
        excluded = set(self.payload["excluded_from_membership"])
        actual = {path.relative_to(ROOT).as_posix() for path in ROOT.rglob("*")
                  if path.is_file() and ".git" not in path.relative_to(ROOT).parts
                  and path.relative_to(ROOT).as_posix() not in excluded}
        self.assertEqual({row["relative_path"] for row in self.tsv}, actual)
        for row in self.tsv:
            data = canonical_bytes(ROOT / row["relative_path"])
            self.assertEqual(len(data), int(row["bytes"]))
            self.assertEqual(hashlib.sha256(data).hexdigest(), row["sha256"])

    def test_master15_outputs_have_derived_redistribution_class(self) -> None:
        rows = [row for row in self.tsv if row["relative_path"].startswith("results/master_1_5/")]
        self.assertEqual(len(rows), 17)
        self.assertTrue(all(row["redistribution_class"] == "REDISTRIBUTABLE_DERIVED_OUTPUT" for row in rows))


if __name__ == "__main__":
    unittest.main()
