"""Build the v1.3.1 maintenance inventory without touching prior release manifests.

File sizes and SHA-256 correspond to canonical LF text, or to raw binary bytes,
using the same policy as the root package manifest. Exclude only Git metadata,
the enclosing root manifest, and this manifest itself to avoid self-reference.
"""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

from build_manifest import canonical_bytes, iter_package_files


MEMBERSHIP = "provenance/RELEASE_MANIFEST_v1.3.1.tsv"
FIELDS = ("path", "bytes", "sha256", "payload_role")


def role(name: str) -> str:
    if name.startswith(("results/", "revision_r1/results/")):
        return "derived_result"
    if name.startswith("provenance/"):
        return "provenance"
    if name.endswith((".py", ".R")):
        return "software_or_test"
    if name.startswith(("environment/", "config/", "revision_r1/config/")):
        return "environment_or_configuration"
    return "documentation_or_project_file"


def manifest_rows(root: Path):
    paths = iter_package_files(root, root / "manifest.tsv")
    for path in sorted(paths, key=lambda p: p.relative_to(root).as_posix()):
        name = path.relative_to(root).as_posix()
        if name == MEMBERSHIP:
            continue
        data = canonical_bytes(path)
        yield dict(path=name, bytes=str(len(data)),
                   sha256=hashlib.sha256(data).hexdigest(), payload_role=role(name))


def main():
    root = Path(__file__).resolve().parents[1]
    dest = root / MEMBERSHIP
    with dest.open("w", encoding="utf-8", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows(root))
    print(f"Wrote {MEMBERSHIP}; rebuild manifest.tsv next.")


if __name__ == "__main__":
    main()
