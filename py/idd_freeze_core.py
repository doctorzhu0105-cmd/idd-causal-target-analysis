#!/usr/bin/env python3
"""Freeze the IDD target dossier into an auditable, hash-locked core tree.

The freeze contract is deliberately explicit.  Only files listed in
``FROZEN_FILES`` may enter ``artifacts/core``; missing or unsafe inputs abort
before the existing core is touched.  Publication uses a staging directory
and directory renames, then verifies the published tree before removing the
previous release.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_NAME = "manifest.tsv"
MANIFEST_FIELDS = ("source_path", "destination_path", "bytes", "sha256")
CHUNK_SIZE = 1024 * 1024


@dataclass(frozen=True)
class FrozenFile:
    """A source path and its stable destination inside the core release."""

    source_path: str
    destination_path: str


# This is the complete release contract.  Keep paths relative and POSIX-style
# so the manifest is portable across Windows and Unix workstations.
FROZEN_FILES: tuple[FrozenFile, ...] = (
    # Bulk discovery and WGCNA.
    FrozenFile("artifacts/bulk/gse70362_sample_metadata.tsv", "bulk/gse70362_sample_metadata.tsv"),
    FrozenFile("artifacts/bulk/gse70362_deg.tsv", "bulk/gse70362_deg.tsv"),
    FrozenFile("artifacts/bulk/gse70362_wgcna_modules.tsv", "bulk/gse70362_wgcna_modules.tsv"),
    FrozenFile("artifacts/bulk/gse70362_wgcna_genes.tsv", "bulk/gse70362_wgcna_genes.tsv"),
    FrozenFile("artifacts/bulk/gse70362_wgcna_soft_threshold.tsv", "bulk/gse70362_wgcna_soft_threshold.tsv"),
    FrozenFile("artifacts/bulk/druggable_bulk_candidates.tsv", "bulk/druggable_bulk_candidates.tsv"),
    FrozenFile("logs/idd/bulk_summary.json", "bulk/bulk_summary.json"),
    FrozenFile("logs/idd/gate0_summary.json", "bulk/gate0_summary.json"),
    # Forward MR, colocalization, SMR/HEIDI, Steiger, reverse MR and source
    # validation.  GTEx proxy-tissue colocalization is kept beside the genetic
    # results so its tissue/source labels cannot be separated from the result.
    FrozenFile("artifacts/mr/formal_candidate_genes.tsv", "genetics/formal_candidate_genes.tsv"),
    FrozenFile("artifacts/mr/druggable_mr_screen.tsv", "genetics/druggable_mr_screen.tsv"),
    FrozenFile("artifacts/mr/druggable_lead_instruments.tsv", "genetics/druggable_lead_instruments.tsv"),
    FrozenFile("artifacts/mr/candidate_extraction_summary.tsv", "genetics/candidate_extraction_summary.tsv"),
    FrozenFile("artifacts/mr/instrument_candidates.tsv", "genetics/instrument_candidates.tsv"),
    FrozenFile("artifacts/mr/harmonized_instruments.tsv", "genetics/harmonized_instruments.tsv"),
    FrozenFile("artifacts/mr/instrument_screen_summary.tsv", "genetics/instrument_screen_summary.tsv"),
    FrozenFile("artifacts/mr/mr_results.tsv", "genetics/mr_results.tsv"),
    FrozenFile("artifacts/mr/mr_sensitivity_results.tsv", "genetics/mr_sensitivity_results.tsv"),
    FrozenFile("artifacts/mr/leave_one_out_results.tsv", "genetics/leave_one_out_results.tsv"),
    FrozenFile("artifacts/mr/coloc_results.tsv", "genetics/coloc_results.tsv"),
    FrozenFile("artifacts/mr/smr_heidi_results.tsv", "genetics/smr_heidi_results.tsv"),
    FrozenFile("artifacts/mr/smr_run_metadata.json", "genetics/smr_run_metadata.json"),
    FrozenFile("artifacts/mr/steiger_results.tsv", "genetics/steiger_results.tsv"),
    FrozenFile("artifacts/mr/reverse_mr_results.tsv", "genetics/reverse_mr_results.tsv"),
    FrozenFile("artifacts/mr/genetics_run_summary.json", "genetics/genetics_run_summary.json"),
    FrozenFile("artifacts/mr/source_checks_eqtl_ld.tsv", "genetics/source_checks_eqtl_ld.tsv"),
    FrozenFile("artifacts/mr/source_validation.tsv", "genetics/source_validation.tsv"),
    FrozenFile("artifacts/mr/finngen_metadata.json", "genetics/finngen_metadata.json"),
    FrozenFile("artifacts/mr/gtex_eqtl_summary.tsv", "genetics/gtex_eqtl_summary.tsv"),
    FrozenFile("artifacts/mr/gtex_coloc_results.tsv", "genetics/gtex_coloc_results.tsv"),
    FrozenFile("artifacts/mr/gtex_coloc_summary.json", "genetics/gtex_coloc_summary.json"),
    # Human disc single-cell localization.  The UMAP is included as a visual
    # audit artifact; it does not authorize a disease-effect claim.
    FrozenFile("artifacts/scrna/gse230809_sample_manifest.tsv", "scrna/gse230809_sample_manifest.tsv"),
    FrozenFile("artifacts/scrna/gse230809_qc.tsv", "scrna/gse230809_qc.tsv"),
    FrozenFile("artifacts/scrna/gse230809_cluster_annotations.tsv", "scrna/gse230809_cluster_annotations.tsv"),
    FrozenFile("artifacts/scrna/target_localization.tsv", "scrna/target_localization.tsv"),
    FrozenFile("logs/idd/scrna_summary.json", "scrna/scrna_summary.json"),
    FrozenFile("figures/idd/umap_idd_cell_types.png", "scrna/umap_idd_cell_types.png"),
    # Drug, delivery and docking evidence.  docking_assessment.tsv is a
    # required status/result contract: a truthful NOT_APPLICABLE row is valid,
    # whereas silently omitting docking would make the dossier incomplete.
    FrozenFile("artifacts/immune_drug/open_targets_summary.tsv", "drug_delivery/open_targets_summary.tsv"),
    FrozenFile("artifacts/immune_drug/chembl_summary.tsv", "drug_delivery/chembl_summary.tsv"),
    FrozenFile("artifacts/immune_drug/drug_interactions.tsv", "drug_delivery/drug_interactions.tsv"),
    FrozenFile("artifacts/immune_drug/delivery_assessment.tsv", "drug_delivery/delivery_assessment.tsv"),
    FrozenFile("artifacts/immune_drug/docking_tool_metadata.json", "drug_delivery/docking_tool_metadata.json"),
    FrozenFile("artifacts/immune_drug/docking_assessment.tsv", "drug_delivery/docking_assessment.tsv"),
    FrozenFile("artifacts/immune_drug/docking_summary.json", "drug_delivery/docking_summary.json"),
    FrozenFile("logs/idd/drug_evidence_summary.json", "drug_delivery/drug_evidence_summary.json"),
    # Integrated target ranking and priority nomination.
    FrozenFile("artifacts/targets/target_evidence_matrix.tsv", "targets/target_evidence_matrix.tsv"),
    FrozenFile("artifacts/targets/target_ranking.tsv", "targets/target_ranking.tsv"),
    FrozenFile("artifacts/targets/priority_targets.tsv", "targets/priority_targets.tsv"),
    FrozenFile("artifacts/targets/locus_alternatives.tsv", "targets/locus_alternatives.tsv"),
    FrozenFile("artifacts/targets/gtex_targets.tsv", "targets/gtex_targets.tsv"),
    FrozenFile("artifacts/targets/ranking_summary.json", "targets/ranking_summary.json"),
    # Analysis contracts and source/environment metadata.
    FrozenFile("config/idd_analysis.json", "metadata/idd_analysis.json"),
    FrozenFile("config/idd_sources.json", "metadata/idd_sources.json"),
    FrozenFile("config/mr_sources.json", "metadata/mr_sources.json"),
    FrozenFile("config/idd_docking_structures.json", "metadata/idd_docking_structures.json"),
    FrozenFile("env/environment.yml", "metadata/environment.yml"),
    FrozenFile("logs/idd/scrna_source_manifest.json", "metadata/scrna_source_manifest.json"),
    FrozenFile("HANDOFF/IDD_TARGET_DOSSIER.md", "reports/IDD_TARGET_DOSSIER.md"),
)


class FreezeError(RuntimeError):
    """Raised when the release contract cannot be satisfied."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="verify the existing artifacts/core tree without changing it",
    )
    return parser.parse_args(argv)


def _relative_path(value: str, *, label: str) -> Path:
    """Parse a safe repository-relative path and reject traversal."""

    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise FreezeError(f"Unsafe {label}: {value}")
    if not value or path == Path("."):
        raise FreezeError(f"Empty {label}")
    return path


def _normalise_rel(path: Path) -> str:
    return path.as_posix()


def _validate_contract() -> None:
    sources = [_normalise_rel(_relative_path(item.source_path, label="source path")) for item in FROZEN_FILES]
    destinations = [
        _normalise_rel(_relative_path(item.destination_path, label="destination path"))
        for item in FROZEN_FILES
    ]
    if len(sources) != len(set(sources)):
        raise FreezeError("FROZEN_FILES contains duplicate source paths")
    if len(destinations) != len(set(destinations)):
        raise FreezeError("FROZEN_FILES contains duplicate destination paths")
    if MANIFEST_NAME in destinations:
        raise FreezeError("The manifest cannot be a frozen source destination")


def _reject_symlink_components(root: Path, relative: Path, *, label: str) -> None:
    current = root
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise FreezeError(f"Symlink is not allowed in {label}: {relative.as_posix()}")


def _source_path(repo_root: Path, relative: str) -> Path:
    rel = _relative_path(relative, label="source path")
    _reject_symlink_components(repo_root, rel, label="source")
    candidate = (repo_root / rel).resolve(strict=False)
    root_resolved = repo_root.resolve()
    try:
        candidate.relative_to(root_resolved)
    except ValueError as exc:
        raise FreezeError(f"Source escapes repository root: {relative}") from exc
    if not candidate.is_file():
        raise FreezeError(f"Required freeze input missing or not a file: {relative}")
    return candidate


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(CHUNK_SIZE)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _preflight_sources(repo_root: Path) -> list[tuple[FrozenFile, Path]]:
    _validate_contract()
    resolved: list[tuple[FrozenFile, Path]] = []
    missing: list[str] = []
    for item in FROZEN_FILES:
        try:
            resolved.append((item, _source_path(repo_root, item.source_path)))
        except FreezeError:
            missing.append(item.source_path)
    if missing:
        raise FreezeError(
            "Required freeze inputs are unavailable; publication aborted: "
            + ", ".join(missing)
        )
    return resolved


def _write_manifest(path: Path, rows: Iterable[dict[str, object]]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    ordered = sorted(rows, key=lambda row: str(row["destination_path"]))
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(MANIFEST_FIELDS),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(ordered)
    temporary.replace(path)


def _read_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise FreezeError(f"Manifest is missing or unsafe: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(MANIFEST_FIELDS):
            raise FreezeError(
                f"Manifest header mismatch: expected {list(MANIFEST_FIELDS)}, got {reader.fieldnames}"
            )
        rows = list(reader)
    if not rows:
        raise FreezeError("Manifest is empty")
    return rows


def _expected_maps() -> tuple[dict[str, FrozenFile], dict[str, FrozenFile]]:
    by_source = {_normalise_rel(_relative_path(item.source_path, label="source path")): item for item in FROZEN_FILES}
    by_destination = {
        _normalise_rel(_relative_path(item.destination_path, label="destination path")): item
        for item in FROZEN_FILES
    }
    return by_source, by_destination


def _validate_inventory(core_dir: Path, rows: list[dict[str, str]]) -> None:
    """Validate manifest rows and every regular file below ``core_dir``."""

    by_source, by_destination = _expected_maps()
    seen_sources: set[str] = set()
    seen_destinations: set[str] = set()
    for row in rows:
        source = _normalise_rel(_relative_path(row.get("source_path", ""), label="manifest source path"))
        destination = _normalise_rel(
            _relative_path(row.get("destination_path", ""), label="manifest destination path")
        )
        if source in seen_sources:
            raise FreezeError(f"Manifest repeats source path: {source}")
        if destination in seen_destinations:
            raise FreezeError(f"Manifest repeats destination path: {destination}")
        seen_sources.add(source)
        seen_destinations.add(destination)
        expected_source = by_source.get(source)
        expected_destination = by_destination.get(destination)
        if expected_source is None or expected_destination is None or expected_source != expected_destination:
            raise FreezeError(f"Manifest path is outside the freeze whitelist: {source} -> {destination}")
        try:
            byte_count = int(row.get("bytes", ""))
        except ValueError as exc:
            raise FreezeError(f"Manifest byte count is invalid: {destination}") from exc
        digest = row.get("sha256", "")
        if byte_count < 0 or len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            raise FreezeError(f"Manifest hash or byte count is invalid: {destination}")
        destination_path = core_dir / _relative_path(destination, label="manifest destination path")
        _reject_symlink_components(core_dir, _relative_path(destination, label="manifest destination path"), label="core")
        if not destination_path.is_file():
            raise FreezeError(f"Frozen file is missing: {destination}")
        if destination_path.stat().st_size != byte_count:
            raise FreezeError(f"Frozen file byte count differs: {destination}")
        if sha256_file(destination_path) != digest:
            raise FreezeError(f"Frozen file SHA-256 differs: {destination}")

    expected_sources = set(by_source)
    expected_destinations = set(by_destination)
    if seen_sources != expected_sources or seen_destinations != expected_destinations:
        raise FreezeError(
            "Manifest inventory differs from whitelist: "
            f"missing_sources={sorted(expected_sources - seen_sources)}, "
            f"missing_destinations={sorted(expected_destinations - seen_destinations)}"
        )

    actual_files: set[str] = set()
    for candidate in core_dir.rglob("*"):
        relative = candidate.relative_to(core_dir)
        if candidate.is_symlink():
            raise FreezeError(f"Symlink is not allowed in frozen core: {relative.as_posix()}")
        if candidate.is_file():
            actual_files.add(_normalise_rel(relative))
    expected_files = expected_destinations | {MANIFEST_NAME}
    if actual_files != expected_files:
        raise FreezeError(
            "Frozen core file inventory differs: "
            f"missing={sorted(expected_files - actual_files)}, "
            f"extra={sorted(actual_files - expected_files)}"
        )
    # Empty or unrelated directories are also drift: every directory must be a
    # parent needed by an expected file (the root itself is allowed).
    expected_directories = {Path(".")}
    for relative in expected_files:
        parent = Path(relative).parent
        while parent != Path("."):
            expected_directories.add(parent)
            parent = parent.parent
    for candidate in core_dir.rglob("*"):
        if candidate.is_dir() and not candidate.is_symlink():
            relative = candidate.relative_to(core_dir)
            if relative not in expected_directories:
                raise FreezeError(f"Unexpected directory in frozen core: {relative.as_posix()}")


def verify_core(repo_root: Path, core_dir: Path | None = None) -> list[dict[str, str]]:
    """Verify the existing core and return its manifest rows."""

    root = repo_root.resolve()
    requested = core_dir or root / "artifacts" / "core"
    if requested.is_symlink() or not requested.is_dir():
        raise FreezeError(f"Frozen core directory is missing or unsafe: {requested}")
    target = requested.resolve()
    try:
        target.relative_to(root)
    except ValueError as exc:
        raise FreezeError(f"Frozen core escapes repository root: {requested}") from exc
    rows = _read_manifest(target / MANIFEST_NAME)
    _validate_inventory(target, rows)
    return rows


def _copy_to_staging(
    repo_root: Path,
    stage_dir: Path,
    resolved: Iterable[tuple[FrozenFile, Path]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item, source in resolved:
        destination = _relative_path(item.destination_path, label="destination path")
        target = stage_dir / destination
        target.parent.mkdir(parents=True, exist_ok=True)
        source_hash_before = sha256_file(source)
        source_bytes_before = source.stat().st_size
        shutil.copyfile(source, target)
        source_hash_after = sha256_file(source)
        source_bytes_after = source.stat().st_size
        if (source_bytes_before, source_hash_before) != (source_bytes_after, source_hash_after):
            raise FreezeError(f"Source changed while staging: {item.source_path}")
        target_bytes = target.stat().st_size
        target_hash = sha256_file(target)
        if (target_bytes, target_hash) != (source_bytes_after, source_hash_after):
            raise FreezeError(f"Staged copy differs from source: {item.source_path}")
        rows.append(
            {
                "source_path": _normalise_rel(_relative_path(item.source_path, label="source path")),
                "destination_path": _normalise_rel(destination),
                "bytes": target_bytes,
                "sha256": target_hash,
            }
        )
    _write_manifest(stage_dir / MANIFEST_NAME, rows)
    return rows


def _remove_directory(path: Path) -> None:
    if path.exists() or path.is_symlink():
        if path.is_symlink() or path.is_file():
            path.unlink()
        else:
            shutil.rmtree(path)


def _publish_verified(stage_dir: Path, core_dir: Path, rows: list[dict[str, object]]) -> None:
    """Replace ``core_dir`` transactionally and retain rollback on verify failure."""

    core_dir.parent.mkdir(parents=True, exist_ok=True)
    if core_dir.is_symlink() or (core_dir.exists() and not core_dir.is_dir()):
        raise FreezeError(f"Existing core path is unsafe: {core_dir}")
    backup: Path | None = None
    if core_dir.exists():
        backup = core_dir.parent / f".core-backup-{uuid.uuid4().hex}"
        os.replace(core_dir, backup)
    try:
        os.replace(stage_dir, core_dir)
        _validate_inventory(core_dir, [dict(row) for row in rows])
    except Exception:
        _remove_directory(core_dir)
        if backup is not None and backup.exists():
            os.replace(backup, core_dir)
        raise
    if backup is not None and backup.exists():
        _remove_directory(backup)


def freeze_core(repo_root: Path) -> list[dict[str, object]]:
    """Build and atomically publish a complete core release."""

    root = repo_root.resolve()
    resolved = _preflight_sources(root)
    artifacts_dir = root / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    stage_dir = Path(tempfile.mkdtemp(prefix=".core-staging-", dir=str(artifacts_dir)))
    try:
        rows = _copy_to_staging(root, stage_dir, resolved)
        _validate_inventory(stage_dir, [dict(row) for row in rows])
        core_dir = artifacts_dir / "core"
        _publish_verified(stage_dir, core_dir, rows)
        return rows
    except Exception:
        if stage_dir.exists():
            _remove_directory(stage_dir)
        raise


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        root = args.repo_root.resolve()
        if args.verify_only:
            rows = verify_core(root)
            print(json.dumps({"status": "verified", "files": len(rows)}, sort_keys=True))
        else:
            rows = freeze_core(root)
            print(json.dumps({"status": "published", "files": len(rows)}, sort_keys=True))
    except (OSError, FreezeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
