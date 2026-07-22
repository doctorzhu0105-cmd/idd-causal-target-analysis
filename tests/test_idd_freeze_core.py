import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py"))

from idd_freeze_core import (  # noqa: E402
    FROZEN_FILES,
    MANIFEST_FIELDS,
    MANIFEST_NAME,
    FreezeError,
    freeze_core,
    sha256_file,
    verify_core,
)


class IddCoreFreezeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self._write_complete_inputs()

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _write_complete_inputs(self, omit: str | None = None) -> None:
        for index, item in enumerate(FROZEN_FILES):
            if item.source_path == omit:
                continue
            path = self.root / Path(item.source_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.suffix == ".json":
                path.write_text(json.dumps({"source": item.source_path}, sort_keys=True) + "\n", encoding="utf-8")
            else:
                path.write_bytes(f"source-{index}\n".encode("utf-8"))

    def test_freeze_writes_exact_manifest_and_is_verifiable(self) -> None:
        rows = freeze_core(self.root)
        core = self.root / "artifacts" / "core"
        self.assertEqual(len(rows), len(FROZEN_FILES))
        expected_rows = sorted(
            ({key: str(value) for key, value in row.items()} for row in rows),
            key=lambda row: row["destination_path"],
        )
        self.assertEqual(verify_core(self.root), expected_rows)
        manifest = core / MANIFEST_NAME
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, list(MANIFEST_FIELDS))
            manifest_rows = list(reader)
        self.assertEqual(len(manifest_rows), len(FROZEN_FILES))
        self.assertNotIn(MANIFEST_NAME, {row["destination_path"] for row in manifest_rows})
        for row in manifest_rows:
            destination = core / Path(row["destination_path"])
            self.assertEqual(int(row["bytes"]), destination.stat().st_size)
            self.assertEqual(row["sha256"], sha256_file(destination))

    def test_missing_required_input_fails_before_touching_existing_core(self) -> None:
        freeze_core(self.root)
        core = self.root / "artifacts" / "core"
        marker = core / "targets" / "priority_targets.tsv"
        original = marker.read_bytes()
        missing = FROZEN_FILES[-1].source_path
        (self.root / Path(missing)).unlink()
        with self.assertRaisesRegex(FreezeError, "publication aborted"):
            freeze_core(self.root)
        self.assertEqual(marker.read_bytes(), original)
        self.assertTrue((core / MANIFEST_NAME).is_file())
        verify_core(self.root)

    def test_verify_rejects_tampering_and_extra_files(self) -> None:
        freeze_core(self.root)
        core = self.root / "artifacts" / "core"
        target = core / "bulk" / "gse70362_deg.tsv"
        target.write_bytes(b"tampered\n")
        with self.assertRaisesRegex(FreezeError, "SHA-256 differs"):
            verify_core(self.root)
        # Restore the release, then test exact-inventory rejection separately.
        freeze_core(self.root)
        (core / "unexpected.txt").write_text("drift\n", encoding="utf-8")
        with self.assertRaisesRegex(FreezeError, "inventory differs"):
            verify_core(self.root)

    def test_manifest_cannot_claim_non_whitelisted_paths(self) -> None:
        freeze_core(self.root)
        manifest = self.root / "artifacts" / "core" / MANIFEST_NAME
        lines = manifest.read_text(encoding="utf-8").splitlines()
        lines[1] = lines[1].replace("bulk/", "../", 1)
        manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
        with self.assertRaises(FreezeError):
            verify_core(self.root)


if __name__ == "__main__":
    unittest.main()
