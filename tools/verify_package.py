#!/usr/bin/env python3
"""Verify a code-only GitHub upload package without network access."""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import re
from pathlib import Path


TEXT_SUFFIXES = {
    "",
    ".cff",
    ".gitignore",
    ".json",
    ".lock",
    ".md",
    ".py",
    ".r",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}
FORBIDDEN_SUFFIXES = {
    ".bed",
    ".bim",
    ".docx",
    ".exe",
    ".fam",
    ".gz",
    ".h5ad",
    ".pdf",
    ".pkl",
    ".png",
    ".rds",
    ".tiff",
    ".xlsx",
}
SENSITIVE_PATTERNS = {
    "personal_email": re.compile(r"\b[A-Z0-9._%+-]+@(?!example\.)[A-Z0-9.-]+\.[A-Z]{2,}\b", re.I),
    "local_user_path": re.compile(r"(?:[A-Z]:\\Users\\|D:\\研究生\\)", re.I),
    "private_key": re.compile(r"BEGIN (?:RSA |OPENSSH )?PRIVATE KEY", re.I),
    "credential_assignment": re.compile(
        r"(?:api[_-]?key|access[_-]?token|client[_-]?secret|password)\s*[:=]\s*['\"][^'\"]+",
        re.I,
    ),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    actual = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path != manifest_path
    }
    if listed != actual:
        failures.append(
            f"manifest_membership_mismatch:missing={sorted(listed - actual)}:unlisted={sorted(actual - listed)}"
        )

    for row in rows:
        path = root / Path(row["path"])
        if not path.is_file():
            continue
        if path.stat().st_size != int(row["bytes"]):
            failures.append(f"byte_mismatch:{row['path']}")
        if sha256_file(path).lower() != row["sha256"].lower():
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
        text = path.read_text(encoding="utf-8", errors="replace")
        for label, pattern in SENSITIVE_PATTERNS.items():
            if label == "personal_email" and relative == "environment/renv.lock":
                continue
            if pattern.search(text):
                failures.append(f"{label}:{relative}")

    for json_path in sorted(root.rglob("*.json")):
        try:
            json.loads(json_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            failures.append(f"invalid_json:{json_path.relative_to(root).as_posix()}:{exc}")

    for python_path in sorted(root.rglob("*.py")):
        try:
            ast.parse(python_path.read_text(encoding="utf-8"), filename=str(python_path))
        except (SyntaxError, UnicodeDecodeError) as exc:
            failures.append(f"invalid_python:{python_path.relative_to(root).as_posix()}:{exc}")

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
            module
            for module in imported_modules
            if module.startswith(("idd_", "mr_")) and module not in local_modules
        )
        if missing_local:
            failures.append(
                f"missing_local_python_dependency:{python_path.relative_to(root).as_posix()}:{missing_local}"
            )

    makefile_text = (root / "Makefile").read_text(encoding="utf-8")
    referenced_code = set(re.findall(r"(?:py|R|tools)/[A-Za-z0-9_.-]+\.(?:py|R)", makefile_text))
    missing_makefile_code = sorted(path for path in referenced_code if not (root / path).is_file())
    if missing_makefile_code:
        failures.append(f"missing_makefile_code:{missing_makefile_code}")

    inventory_path = root / "docs" / "CODE_INVENTORY.tsv"
    with inventory_path.open(encoding="utf-8", newline="") as handle:
        inventory_rows = list(csv.DictReader(handle, delimiter="\t"))
    inventory_paths = {row["path"] for row in inventory_rows}
    expected_inventory_paths = {
        path.relative_to(root).as_posix()
        for directory in ("py", "R", "tests")
        for path in (root / directory).glob("*")
        if path.is_file() and path.suffix.lower() in {".py", ".r"}
    }
    if inventory_paths != expected_inventory_paths:
        failures.append(
            "code_inventory_mismatch:"
            f"missing={sorted(expected_inventory_paths - inventory_paths)}:"
            f"extra={sorted(inventory_paths - expected_inventory_paths)}"
        )

    expected_roots = {"R", "config", "data", "docs", "environment", "manuscript", "py", "tests", "tools"}
    observed_roots = {path.name for path in root.iterdir() if path.is_dir()}
    missing_roots = sorted(expected_roots - observed_roots)
    if missing_roots:
        failures.append(f"missing_directories:{missing_roots}")

    status = {
        "status": "PASS" if not failures else "FAIL",
        "file_count": len(actual),
        "manifest_rows": len(rows),
        "failure_count": len(failures),
        "failures": failures,
    }
    print(json.dumps(status, ensure_ascii=False))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
