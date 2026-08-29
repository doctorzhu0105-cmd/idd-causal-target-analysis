#!/usr/bin/env python3
"""Verify the public reproducibility package without network access."""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import re
from pathlib import Path

from build_manifest import TEXT_SUFFIXES, canonical_bytes, iter_package_files


FORBIDDEN_SUFFIXES = {
    ".bed", ".bim", ".docx", ".exe", ".fam", ".gz", ".h5ad", ".pdf",
    ".pkl", ".png", ".rds", ".tiff", ".xlsx",
}
SENSITIVE_PATTERNS = {
    "personal_email": re.compile(r"\b[A-Z0-9._%+-]+@(?!example\.)[A-Z0-9.-]+\.[A-Z]{2,}\b", re.I),
    "local_user_path": re.compile(r"(?:[A-Z]:\\Users\\|D:\\研究生\\)", re.I),
    "private_key": re.compile(r"BEGIN (?:RSA |OPENSSH )?PRIVATE KEY", re.I),
    "credential_assignment": re.compile(
        r"(?:api[_-]?key|access[_-]?token|client[_-]?secret|password)\s*[:=]\s*['\"][^'\"]+", re.I
    ),
}


def canonical_sha256(path: Path) -> str:
    return hashlib.sha256(canonical_bytes(path)).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    return parser.parse_args()


def main() -> int:
    root = parse_args().root.resolve()
    manifest_path = root / "manifest.tsv"
    failures: list[str] = []
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing manifest: {manifest_path}")

    with manifest_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required_columns = {"path", "bytes", "sha256"}
    if not rows or not required_columns.issubset(rows[0]):
        failures.append("manifest_schema_invalid")

    listed = {row["path"] for row in rows}
    actual_paths = iter_package_files(root, manifest_path)
    actual = {path.relative_to(root).as_posix() for path in actual_paths}
    if listed != actual:
        failures.append(f"manifest_membership_mismatch:missing={sorted(listed - actual)}:unlisted={sorted(actual - listed)}")

    for row in rows:
        path = root / Path(row["path"])
        if not path.is_file():
            continue
        data = canonical_bytes(path)
        if len(data) != int(row["bytes"]):
            failures.append(f"byte_mismatch:{row['path']}")
        if hashlib.sha256(data).hexdigest().lower() != row["sha256"].lower():
            failures.append(f"sha256_mismatch:{row['path']}")

    for relative in sorted(actual):
        path = root / relative
        lowered = relative.lower()
        if lowered.startswith(("data/raw/", "data/processed/", "artifacts/")):
            failures.append(f"excluded_data_present:{relative}")
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            failures.append(f"forbidden_binary_present:{relative}")
        if "api_cache" in lowered or "__pycache__" in lowered:
            failures.append(f"cache_present:{relative}")
        if path.stat().st_size > 50 * 1024 * 1024:
            failures.append(f"file_exceeds_50_mib:{relative}")
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = canonical_bytes(path).decode("utf-8", errors="replace")
        for label, pattern in SENSITIVE_PATTERNS.items():
            if label == "personal_email" and relative == "environment/renv.lock":
                continue
            if pattern.search(text):
                failures.append(f"{label}:{relative}")

    for relative in sorted(path for path in actual if path.lower().endswith(".json")):
        try:
            json.loads((root / relative).read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            failures.append(f"invalid_json:{relative}:{exc}")

    for relative in sorted(path for path in actual if path.lower().endswith(".py")):
        try:
            ast.parse((root / relative).read_text(encoding="utf-8"), filename=relative)
        except (SyntaxError, UnicodeDecodeError) as exc:
            failures.append(f"invalid_python:{relative}:{exc}")

    local_modules = {path.stem for path in (root / "py").glob("*.py")}
    for python_path in sorted((root / "py").glob("*.py")):
        tree = ast.parse(python_path.read_text(encoding="utf-8"), filename=str(python_path))
        imported_modules: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported_modules.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported_modules.add(node.module.split(".")[0])
        missing_local = sorted(
            module for module in imported_modules
            if module.startswith(("idd_", "mr_", "master14_")) and module not in local_modules
        )
        if missing_local:
            failures.append(f"missing_local_python_dependency:{python_path.relative_to(root).as_posix()}:{missing_local}")

    makefile_text = (root / "Makefile").read_text(encoding="utf-8")
    referenced_code = set(re.findall(r"(?:py|R|tools|analysis_p0_correlated_mr_20260814/scripts)/[A-Za-z0-9_.-]+\.(?:py|R)", makefile_text))
    missing_makefile_code = sorted(path for path in referenced_code if not (root / path).is_file())
    if missing_makefile_code:
        failures.append(f"missing_makefile_code:{missing_makefile_code}")

    inventory_path = root / "docs" / "CODE_INVENTORY.tsv"
    with inventory_path.open(encoding="utf-8", newline="") as handle:
        inventory_rows = list(csv.DictReader(handle, delimiter="\t"))
    inventory_paths = {row["path"] for row in inventory_rows}
    inventory_dirs = (root / "py", root / "R", root / "tests", root / "analysis_p0_correlated_mr_20260814" / "scripts")
    expected_inventory_paths = {
        path.relative_to(root).as_posix()
        for directory in inventory_dirs
        for path in directory.glob("*")
        if path.is_file() and path.suffix.lower() in {".py", ".r"}
    }
    if inventory_paths != expected_inventory_paths:
        failures.append(f"code_inventory_mismatch:missing={sorted(expected_inventory_paths - inventory_paths)}:extra={sorted(inventory_paths - expected_inventory_paths)}")

    expected_roots = {
        "R", "analysis_p0_correlated_mr_20260814", "config", "data", "docs",
        "environment", "manuscript", "provenance", "py", "results", "tests", "tools",
    }
    observed_roots = {path.name for path in root.iterdir() if path.is_dir() and path.name != ".git"}
    missing_roots = sorted(expected_roots - observed_roots)
    if missing_roots:
        failures.append(f"missing_directories:{missing_roots}")

    release = json.loads((root / "RELEASE_STATUS.json").read_text(encoding="utf-8"))
    expected_release = {
        "identity": "V1.2.0_CATALOGUE_EXTENSION_RELEASE_CONTENT",
        "master_version": "1.5.0",
        "release_state": "VERSIONED_RELEASE_CONTENT",
        "git_tag": "v1.2.0",
        "github_release": "https://github.com/doctorzhu0105-cmd/idd-causal-target-analysis/releases/tag/v1.2.0",
        "zenodo_version": "RESOLVED_BY_PUBLICATION_RECORD",
    }
    for key, expected in expected_release.items():
        if release.get(key) != expected:
            failures.append(f"release_status_mismatch:{key}:{release.get(key)!r}")

    required_current = {
        "results/master_1_4/table3_evidence_matrix.tsv",
        "results/master_1_4/supplementary_figure_s2_source.tsv",
        "results/master_1_4/multiplicity_summary.tsv",
        "results/master_1_4/p2ry13_multiplicity_path.tsv",
        "results/master_1_4/ld_reference_bootstrap_summary.tsv",
        "provenance/ANALYSIS_VERSION.json",
        "provenance/OUTPUT_HASHES.tsv",
        "provenance/RELEASE_MANIFEST.tsv",
        "provenance/REPRODUCIBILITY_SCOPE.md",
        "provenance/RELEASE_MANIFEST_v1.2.0.tsv",
        "provenance/RELEASE_MANIFEST_v1.2.0.json",
        "results/master_1_5/CATALOGUE_EFFECTIVE_K.tsv",
        "results/master_1_5/CATALOGUE_GIVW_LD_LADDER.tsv",
        "results/master_1_5/CATALOGUE_JOINT_EXPOSURE_INFORMATION.tsv",
        "results/master_1_5/CATALOGUE_PRECISION_PAIRED_CHANGES.tsv",
        "results/master_1_5/K1_SECOND_ORDER_MULTIPLICITY_REPLICATION.tsv",
    }
    for relative in sorted(required_current - actual):
        failures.append(f"missing_current_release_surface:{relative}")

    status = {
        "status": "PASS" if not failures else "FAIL",
        "file_count": len(actual),
        "manifest_rows": len(rows),
        "failure_count": len(failures),
        "failures": failures,
        "manifest_hash_policy": "LF-normalized text bytes; raw bytes for non-text files; .git excluded",
    }
    print(json.dumps(status, ensure_ascii=False))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
