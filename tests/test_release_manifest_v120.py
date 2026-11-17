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

    def test_historical_paths_and_preserved_hashes(self) -> None:
        # The v1.2.0 manifest remains an immutable release snapshot, not the
        # membership manifest for a later forward extension. Existing results,
        # scientific code and provenance must remain unchanged.
        explicitly_updated = {
            "README.md", "README.zh-CN.md", "RELEASE_STATUS.json",
            "UPLOAD_GUIDE.zh-CN.md", "docs/CODE_INVENTORY.tsv",
            "docs/REPRODUCIBILITY.md", "docs/DATA_ACCESS.md",
            "provenance/REPRODUCIBILITY_SCOPE.md", "tools/verify_package.py",
            "tests/test_release_manifest_v120.py",
        }
        historical = {row["relative_path"] for row in self.tsv}
        self.assertTrue(explicitly_updated <= historical)
        for row in self.tsv:
            relative = row["relative_path"]
            path = ROOT / relative
            self.assertTrue(path.is_file(), relative)
            if relative in explicitly_updated:
                continue
            data = canonical_bytes(path)
            self.assertEqual(len(data), int(row["bytes"]), relative)
            self.assertEqual(hashlib.sha256(data).hexdigest(), row["sha256"], relative)

    def test_master15_outputs_have_derived_redistribution_class(self) -> None:
        rows = [row for row in self.tsv if row["relative_path"].startswith("results/master_1_5/")]
        self.assertEqual(len(rows), 17)
        self.assertTrue(all(row["redistribution_class"] == "REDISTRIBUTABLE_DERIVED_OUTPUT" for row in rows))


if __name__ == "__main__":
    unittest.main()
