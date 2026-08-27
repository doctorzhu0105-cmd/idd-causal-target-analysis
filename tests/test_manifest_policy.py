from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from build_manifest import canonical_bytes, iter_package_files  # noqa: E402


class ManifestPolicyTests(unittest.TestCase):
    def test_text_newlines_are_canonical(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "example.md"
            path.write_bytes(b"a\r\nb\rc\n")
            self.assertEqual(canonical_bytes(path), b"a\nb\nc\n")

    def test_git_metadata_is_outside_package_membership(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.tsv"
            manifest.write_text("", encoding="utf-8")
            (root / "kept.txt").write_text("kept\n", encoding="utf-8")
            (root / ".git").mkdir()
            (root / ".git" / "config").write_text("secret-like metadata\n", encoding="utf-8")
            files = [path.relative_to(root).as_posix() for path in iter_package_files(root, manifest)]
            self.assertEqual(files, ["kept.txt"])


if __name__ == "__main__":
    unittest.main()
