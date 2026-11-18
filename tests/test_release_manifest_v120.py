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
        # Pin both generations for the three deliberately maintained files.
        # This also works in an exported source archive without Git history.
        maintenance_versions = {
            "Makefile": (
                (3019, "a4954e7274e795c9f9a47d649edac30545d1cedfe29d5fc0073f08cbbb7cd67d"),
                (3137, "d74aad70fe128e10e0601274d87e8df0f28d0bcd8736eb298c9452feb83df510"),
            ),
            "py/idd_manuscript_verify.py": (
                (28703, "fe50e31a91d92dece638cfbb84e4b6e895970579af1113bd7537b2e4f4f8460b"),
                (28701, "e0ce8ff26fb36bca880d530947b3f3cb047d6ff8fde28afd0db4c97224c7f2cf"),
            ),
            "py/idd_ukb_replication.py": (
                (24831, "395bb105a7132d07d52d9c27d8b18fe519f05f8f82d197ac0658fbe3f9772c46"),
                (24825, "bcc216412c67aa6484a3892cf8d8bdd8b125b8d993153c18e04e682213b56d6c"),
            ),
        }
        self.assertTrue(maintenance_versions.keys() <= historical)
        for row in self.tsv:
            relative = row["relative_path"]
            path = ROOT / relative
            self.assertTrue(path.is_file(), relative)
            if relative in explicitly_updated:
                continue
            data = canonical_bytes(path)
            if relative in maintenance_versions:
                frozen, candidate = maintenance_versions[relative]
                self.assertEqual((int(row["bytes"]), row["sha256"]), frozen, relative)
                self.assertEqual((len(data), hashlib.sha256(data).hexdigest()), candidate, relative)
                continue
            self.assertEqual(len(data), int(row["bytes"]), relative)
            self.assertEqual(hashlib.sha256(data).hexdigest(), row["sha256"], relative)

    def test_master15_outputs_have_derived_redistribution_class(self) -> None:
        rows = [row for row in self.tsv if row["relative_path"].startswith("results/master_1_5/")]
        self.assertEqual(len(rows), 17)
        self.assertTrue(all(row["redistribution_class"] == "REDISTRIBUTABLE_DERIVED_OUTPUT" for row in rows))


if __name__ == "__main__":
    unittest.main()
