from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path


class InputError(ValueError):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists() or not path.is_file():
        raise InputError(f"required input does not exist: {path}")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise InputError(f"input has no header: {path}")
        rows = list(reader)
    return list(reader.fieldnames), rows


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def audit_attachment(path: Path, output_dir: Path) -> dict:
    fields, rows = read_table(path)
    key_candidates = [key for key in ("ensembl_gene_id", "endpoint", "outcome", "scenario", "r2_threshold") if key in fields]
    duplicate_rows = 0
    if key_candidates:
        seen: set[tuple[str, ...]] = set()
        for row in rows:
            key = tuple(row.get(field, "") for field in key_candidates)
            if key in seen:
                duplicate_rows += 1
            seen.add(key)
    payload = {
        "analysis": "attachment-audit",
        "input": str(path),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
        "rows": len(rows),
        "columns": len(fields),
        "key_fields": key_candidates,
        "duplicate_key_rows": duplicate_rows,
        "status_counts": dict(Counter(row.get("status", row.get("match_status", "")) for row in rows)),
        "p_q_bh_fields": [field for field in fields if any(token in field.lower() for token in ("pvalue", "q_pvalue", "bh", "by", "fdr"))],
        "mode": "ROW_LEVEL_AUDIT",
    }
    write_json(output_dir / "attachment_audit.json", payload)
    return payload


def k_transition(path: Path, output_dir: Path) -> dict:
    fields, rows = read_table(path)
    required = {"K_from", "K_to", "K_transition"}
    if not required.issubset(fields):
        raise InputError(f"k-transition input missing fields {sorted(required - set(fields))}: {path}")
    status_field = "match_status" if "match_status" in fields else "status"
    matched = [row for row in rows if row.get(status_field) == "MATCHED"]
    derived = Counter()
    disagreements = 0
    for row in matched:
        try:
            k_from = float(row["K_from"])
            k_to = float(row["K_to"])
        except (TypeError, ValueError) as exc:
            raise InputError(f"non-numeric K values in row {row}") from exc
        if k_from == 1 and k_to == 1:
            label = "K1_to_K1"
        elif k_from == 1 and k_to >= 2:
            label = "K1_to_multi"
        elif k_from >= 2 and k_to >= 2:
            label = "multi_to_multi"
        else:
            label = "UNCLASSIFIED"
        derived[label] += 1
        if row.get("K_transition") != label:
            disagreements += 1
    payload = {
        "analysis": "k-transition",
        "input": str(path),
        "sha256": sha256(path),
        "rows": len(rows),
        "matched_rows": len(matched),
        "transition_counts": dict(derived),
        "serialized_label_disagreements": disagreements,
        "status_field": status_field,
        "mode": "ROW_LEVEL_RECONSTRUCTION",
    }
    write_json(output_dir / "k_transition.json", payload)
    return payload


def historical_multiplicity_audit(path: Path, output_dir: Path) -> dict:
    fields, rows = read_table(path)
    expected = {"pvalue", "q_pvalue", "BH_q", "BY_q", "official_K", "bh_membership_changed", "by_membership_changed"}
    missing = sorted(expected - set(fields))
    if missing:
        raise InputError(f"multiplicity input missing required historical fields {missing}")
    k1 = [row for row in rows if row.get("official_K") == "1"]
    payload = {
        "analysis": "historical-multiplicity-audit",
        "input": str(path),
        "sha256": sha256(path),
        "rows": len(rows),
        "columns": len(fields),
        "k1_rows": len(k1),
        "bh_changed_rows": sum(row.get("bh_membership_changed") == "True" for row in rows),
        "by_changed_rows": sum(row.get("by_membership_changed") == "True" for row in rows),
        "p_q_bh_fields": sorted(expected),
        "relationship": "historical_attachment_separate_from_S6_S16",
        "regeneration_source": "recovered Stage10R1 propagation; run_scientific_analysis.py k1-propagation; upstream inputs required",
        "mode": "HISTORICAL_ATTACHMENT_AUDIT",
    }
    write_json(output_dir / "historical_multiplicity_audit.json", payload)
    return payload


def attenuation_source_audit(path: Path, output_dir: Path) -> dict:
    fields, rows = read_table(path)
    required = {"gene", "endpoint", "threshold"}
    missing = sorted(required - set(fields))
    payload = {
        "analysis": "attenuation-source-audit",
        "input": str(path),
        "sha256": sha256(path),
        "rows": len(rows),
        "columns": len(fields),
        "required_fields": sorted(required),
        "missing_required_fields": missing,
        "source_class": "historical_C_M0_or_planned_configuration_table",
        "calculation_status": "SOURCE_SUMMARY_ONLY",
        "full_raw_to_result_source": "historical Batch14 source recovered as inert reference; portable full orchestration not verified",
        "mode": "SOURCE_CHAIN_AUDIT",
    }
    write_json(output_dir / "attenuation_source_audit.json", payload)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Parameterized, output-isolated release analyses")
    parser.add_argument("--analysis", choices=("attachment-audit", "k-transition", "historical-multiplicity-audit", "attenuation-source-audit"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        source = args.input.resolve()
        destination = args.output_dir.resolve()
        if not source.is_file():
            raise InputError("required input file does not exist")
        if destination.exists():
            raise InputError("output directory must be new")
        if destination == source.parent or source.parent in destination.parents or destination in source.parents:
            raise InputError("output directory must be disjoint from input directory")
        code_root = Path(__file__).resolve().parent
        if destination == code_root or code_root in destination.parents:
            raise InputError("output cannot be inside the code directory")
        destination.mkdir(parents=True)
        if args.analysis == "attachment-audit":
            result = audit_attachment(args.input, args.output_dir)
        elif args.analysis == "k-transition":
            result = k_transition(args.input, args.output_dir)
        elif args.analysis == "historical-multiplicity-audit":
            result = historical_multiplicity_audit(args.input, args.output_dir)
        else:
            result = attenuation_source_audit(args.input, args.output_dir)
    except (InputError, OSError, ValueError) as exc:
        print(f"INPUT_ERROR: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
