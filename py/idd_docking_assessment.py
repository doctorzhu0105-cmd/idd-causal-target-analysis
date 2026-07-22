#!/usr/bin/env python3
"""Apply a fail-closed applicability gate before molecular docking."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Mapping


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    return parser.parse_args()


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def assess_pair(priority: Mapping[str, str], structure: Mapping[str, object]) -> dict[str, object]:
    pdb_ids = [str(value) for value in structure.get("pdb_ids", [])]
    reasons = []
    if not pdb_ids:
        reasons.append("no_human_experimental_pdb")
    if not str(priority.get("top_modulator_smiles", "")).strip():
        reasons.append("selected_modulator_has_no_canonical_smiles")
    if not truthy(priority.get("top_modulator_direction_aligned")):
        reasons.append("selected_modulator_direction_conflicts_with_mr")
    eligible = not reasons
    return {
        "priority_role": priority.get("priority_role", ""),
        "gene_symbol": priority.get("gene_symbol", ""),
        "desired_intervention": priority.get("desired_intervention", ""),
        "selected_modulator": priority.get("top_modulator_name", ""),
        "selected_modulator_id": priority.get("top_modulator_id", ""),
        "uniprot_accession": structure.get("uniprot_accession", ""),
        "experimental_pdb_count": len(pdb_ids),
        "experimental_pdb_ids": ";".join(pdb_ids),
        "canonical_smiles_available": bool(str(priority.get("top_modulator_smiles", "")).strip()),
        "direction_aligned": truthy(priority.get("top_modulator_direction_aligned")),
        "docking_applicable": eligible,
        "docking_status": "ELIGIBLE" if eligible else "NOT_APPLICABLE",
        "reason": ";".join(reasons),
        "vina_affinity_kcal_mol": "not_estimable",
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[Mapping[str, object]]) -> None:
    fields = list(rows[0]) if rows else []
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> int:
    args = parse_args()
    root = args.repo_root.resolve()
    priority = read_tsv(root / "artifacts" / "targets" / "priority_targets.tsv")
    inventory = json.loads(
        (root / "config" / "idd_docking_structures.json").read_text(encoding="utf-8")
    )
    rows = []
    for target in priority:
        symbol = target["gene_symbol"]
        if symbol not in inventory["targets"]:
            raise ValueError(f"Missing structure assessment for priority target: {symbol}")
        rows.append(assess_pair(target, inventory["targets"][symbol]))
    eligible = [row for row in rows if row["docking_applicable"]]
    output_dir = root / "artifacts" / "immune_drug"
    write_tsv(output_dir / "docking_assessment.tsv", rows)
    summary = {
        "status": "ELIGIBLE_PAIR_REQUIRES_VINA" if eligible else "NOT_APPLICABLE",
        "priority_target_count": len(rows),
        "eligible_pair_count": len(eligible),
        "vina_run": False,
        "reason": (
            "No priority target-modulator pair simultaneously has a human experimental PDB, "
            "canonical ligand structure, and MR-aligned intervention direction."
            if not eligible
            else "At least one pair passed the applicability gate; Vina must run before freezing."
        ),
        "structure_source": inventory["source"],
        "structure_checked_at": inventory["checked_at"],
    }
    (output_dir / "docking_summary.json").write_text(
        json.dumps(summary, ensure_ascii=True, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, sort_keys=True))
    return 2 if eligible else 0


if __name__ == "__main__":
    raise SystemExit(main())
