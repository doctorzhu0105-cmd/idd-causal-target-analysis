#!/usr/bin/env python3
"""Build the deterministic package manifest used by verify_package.py."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


TEXT_SUFFIXES = {
    "", ".cff", ".gitignore", ".json", ".lock", ".md", ".py", ".r",
    ".tsv", ".txt", ".yaml", ".yml",
}


def canonical_bytes(path: Path) -> bytes:
    data = path.read_bytes()
    if path.suffix.lower() in TEXT_SUFFIXES:
        data = data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return data


def iter_package_files(root: Path, manifest: Path) -> list[Path]:
    return sorted(
        path for path in root.rglob("*")
        if path.is_file()
        and path != manifest
        and ".git" not in path.relative_to(root).parts
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--manifest", type=Path, default=Path("manifest.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    manifest = args.manifest if args.manifest.is_absolute() else root / args.manifest
    rows = []
    for path in iter_package_files(root, manifest):
        data = canonical_bytes(path)
        rows.append({
            "path": path.relative_to(root).as_posix(),
            "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
        })
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("path", "bytes", "sha256"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {manifest} with {len(rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
