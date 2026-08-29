#!/usr/bin/env python3
"""Build the non-recursive v1.2.0 public release manifests."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from build_manifest import canonical_bytes


ROOT = Path(__file__).resolve().parents[1]
TSV = ROOT / "provenance" / "RELEASE_MANIFEST_v1.2.0.tsv"
JSON = ROOT / "provenance" / "RELEASE_MANIFEST_v1.2.0.json"
EXCLUDED = {"manifest.tsv", TSV.relative_to(ROOT).as_posix(), JSON.relative_to(ROOT).as_posix()}


def classify(relative: str) -> tuple[str, str, str]:
    if relative.startswith("results/master_1_5/"):
        return "governed_master_1_5_output", "MASTER_1.5.0_SCIENTIFIC_CANDIDATE", "REDISTRIBUTABLE_DERIVED_OUTPUT"
    if relative.startswith("results/master_1_4/"):
        return "preserved_v1_1_output", "MASTER_1.4.0_PUBLIC_LAYER", "REDISTRIBUTABLE_DERIVED_OUTPUT"
    if relative.endswith((".py", ".R")) or relative.startswith(("py/", "R/", "tools/", "tests/")):
        return "executable_code_or_test", "PUBLIC_REPOSITORY", "MIT_CODE"
    if relative.startswith(("docs/", "manuscript/", "provenance/")) or relative.endswith((".md", ".json", ".tsv")):
        return "documentation_or_provenance", "PUBLIC_REPOSITORY", "PUBLIC_DOCUMENTATION_OR_METADATA"
    if relative.startswith(("environment/", "config/")):
        return "configuration_or_environment", "PUBLIC_REPOSITORY", "PUBLIC_CONFIGURATION"
    return "repository_support_file", "PUBLIC_REPOSITORY", "PUBLIC_REPOSITORY_CONTENT"


def rows() -> list[dict[str, object]]:
    files = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file()
        and ".git" not in path.relative_to(ROOT).parts
        and path.relative_to(ROOT).as_posix() not in EXCLUDED
    )
    result = []
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        data = canonical_bytes(path)
        role, authority, redistribution = classify(relative)
        result.append({
            "relative_path": relative,
            "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "role": role,
            "source_authority": authority,
            "redistribution_class": redistribution,
        })
    return result


def main() -> int:
    entries = rows()
    TSV.parent.mkdir(parents=True, exist_ok=True)
    with TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=entries[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(entries)
    JSON.write_text(json.dumps({
        "release": "v1.2.0",
        "hash_mode": "LF-canonicalized text and raw bytes for non-text files",
        "excluded_from_membership": sorted(EXCLUDED),
        "files": entries,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"wrote v1.2.0 release manifests with {len(entries)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
