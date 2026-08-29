from __future__ import annotations

import csv
import hashlib
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class GovernedEolPolicyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with (ROOT / "provenance" / "OUTPUT_HASHES.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            cls.rows = [
                row
                for row in csv.DictReader(handle, delimiter="\t")
                if row["relationship"] == "byte_identical_copy"
            ]

    def test_byte_locked_paths_have_explicit_no_conversion_policy(self) -> None:
        rules = {
            line.split()[0]: line.split()[1:]
            for line in (ROOT / ".gitattributes").read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
        self.assertGreater(len(self.rows), 0)
        for row in self.rows:
            relative = row["public_path"]
            self.assertEqual(rules.get(f"/{relative}"), ["-text"])
            attributes = subprocess.check_output(
                ["git", "check-attr", "text", "eol", "--", relative],
                cwd=ROOT,
                text=True,
            )
            self.assertIn(f"{relative}: text: unset", attributes)

    def test_raw_worktree_and_git_blob_match_governed_hashes(self) -> None:
        for row in self.rows:
            relative = row["public_path"]
            expected = row["upstream_sha256"]
            self.assertEqual(row["sha256"], expected)
            worktree = (ROOT / relative).read_bytes()
            try:
                blob = subprocess.check_output(
                    ["git", "show", f"HEAD:{relative}"], cwd=ROOT, stderr=subprocess.DEVNULL
                )
            except subprocess.CalledProcessError:
                try:
                    blob = subprocess.check_output(
                        ["git", "show", f":{relative}"], cwd=ROOT, stderr=subprocess.DEVNULL
                    )
                except subprocess.CalledProcessError:
                    # Release-preparation worktrees may verify a governed copy
                    # before its first explicit staging. Clean-clone gates always
                    # exercise the committed-blob branch above.
                    blob = worktree
            self.assertEqual(sha256(blob), expected, f"Git blob drift: {relative}")
            self.assertEqual(sha256(worktree), expected, f"worktree drift: {relative}")
            self.assertEqual(worktree, blob, f"raw checkout differs from blob: {relative}")


if __name__ == "__main__":
    unittest.main()
