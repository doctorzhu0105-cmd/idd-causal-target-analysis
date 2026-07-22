#!/usr/bin/env python3
"""Integrate IDD evidence and rank targets for local biomaterial delivery."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    return parser.parse_args()


def read_tsv(path: Path, required: Iterable[str] = ()) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    columns = set(rows[0]) if rows else set()
    missing = set(required) - columns
    if missing:
        raise ValueError(f"{path} lacks required columns: {sorted(missing)}")
    return rows


def grouped(rows: Iterable[Mapping[str, str]], key: str = "gene_symbol") -> dict[str, list[Mapping[str, str]]]:
    result: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        result[str(row.get(key, ""))].append(row)
    return dict(result)


def optional_float(value: object) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def optional_int(value: object) -> int | None:
    parsed = optional_float(value)
    return int(parsed) if parsed is not None else None


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def value_or_blank(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return value


def best_numeric(rows: Iterable[Mapping[str, str]], column: str, *, maximum: bool = False) -> Mapping[str, str] | None:
    values = [(optional_float(row.get(column)), row) for row in rows]
    values = [(value, row) for value, row in values if value is not None]
    if not values:
        return None
    chooser = max if maximum else min
    return chooser(values, key=lambda item: item[0])[1]


def index_outcome(rows: Iterable[Mapping[str, str]]) -> dict[tuple[str, str], Mapping[str, str]]:
    result = {}
    for row in rows:
        result[(str(row["gene_symbol"]), str(row["outcome_id"]))] = row
    return result


def select_distinct_loci(
    rows: Iterable[Mapping[str, object]],
    count: int,
    separation_bp: int,
) -> tuple[list[Mapping[str, object]], list[dict[str, object]]]:
    selected: list[Mapping[str, object]] = []
    alternatives: list[dict[str, object]] = []
    for row in rows:
        chromosome = str(row.get("chr_b37", ""))
        start = optional_int(row.get("start_b37"))
        end = optional_int(row.get("end_b37"))
        conflict = None
        if chromosome and start is not None and end is not None:
            for chosen in selected:
                chosen_start = optional_int(chosen.get("start_b37"))
                chosen_end = optional_int(chosen.get("end_b37"))
                if (
                    chromosome == str(chosen.get("chr_b37", ""))
                    and chosen_start is not None
                    and chosen_end is not None
                    and start <= chosen_end + separation_bp
                    and chosen_start <= end + separation_bp
                ):
                    conflict = chosen
                    break
        if conflict is not None:
            alternatives.append(
                {
                    "gene_symbol": row.get("gene_symbol", ""),
                    "alternative_to": conflict.get("gene_symbol", ""),
                    "chr_b37": chromosome,
                    "start_b37": start,
                    "end_b37": end,
                    "reason": f"within_{separation_bp}_bp_priority_locus",
                }
            )
            continue
        selected.append(row)
        if len(selected) == count:
            break
    return selected, alternatives


def pass_smr(row: Mapping[str, str] | None) -> bool:
    if row is None or row.get("status") != "estimable":
        return False
    fdr = optional_float(row.get("p_smr_fdr"))
    heidi = optional_float(row.get("p_heidi"))
    n_heidi = optional_int(row.get("n_snp_heidi"))
    return bool(fdr is not None and fdr < 0.05 and heidi is not None and heidi > 0.01 and n_heidi is not None and n_heidi >= 10)


def delivery_fraction(label: str) -> float:
    if label.startswith("high_"):
        return 1.0
    if label.startswith("moderate_"):
        return 0.5
    return 0.0


def directional_delivery_fraction(
    delivery_row: Mapping[str, str] | None,
    desired_intervention: str,
    direction_aligned: bool,
) -> float:
    """Score payload feasibility for the MR-implied direction, not binding alone."""
    if not delivery_row or desired_intervention == "not_estimable":
        return 0.0
    payloads = str(delivery_row.get("preferred_local_payload_options", "")).lower()
    surface = truthy(delivery_row.get("extracellular_or_surface_location"))
    if desired_intervention == "inhibit_or_reduce_expression":
        return 1.0 if direction_aligned or "sirna" in payloads else 0.5
    if desired_intervention == "increase_or_replace_expression":
        if surface and ("protein" in payloads or "antibody" in payloads or "sirna" in payloads):
            return 1.0
        return 0.5 if "sirna" in payloads else 0.0
    return 0.0


def recommended_local_payload(
    delivery_row: Mapping[str, str] | None,
    desired_intervention: str,
    modulator: Mapping[str, str] | None,
    direction_aligned: bool,
) -> str:
    if not delivery_row or desired_intervention == "not_estimable":
        return "not_estimable"
    surface = truthy(delivery_row.get("extracellular_or_surface_location"))
    has_smiles = bool(modulator and str(modulator.get("canonical_smiles", "")).strip())
    if desired_intervention == "increase_or_replace_expression":
        if surface:
            return "mrna_nanoparticle_in_hydrogel;recombinant_protein_in_hydrogel"
        return "mrna_nanoparticle_in_hydrogel"
    if direction_aligned and has_smiles:
        return "small_molecule_in_hydrogel;sirna_nanoparticle_in_hydrogel"
    return "sirna_nanoparticle_in_hydrogel"


def choose_localization(rows: Iterable[Mapping[str, str]], minimum_fraction: float) -> tuple[Mapping[str, str] | None, bool]:
    localized = [row for row in rows if row.get("status") == "localized"]
    if not localized:
        return None, False
    eligible = [
        row
        for row in localized
        if (optional_float(row.get("cell_fraction_expressing")) or 0.0) >= minimum_fraction
        and (optional_float(row.get("donor_detection_fraction")) or 0.0) >= 0.5
        and (optional_int(row.get("donors")) or 0) >= 2
        and (optional_int(row.get("cells")) or 0) >= 50
    ]
    pool = eligible or localized
    best = max(
        pool,
        key=lambda row: (
            optional_float(row.get("cell_fraction_expressing")) or -1.0,
            optional_float(row.get("mean_donor_expression")) or -1.0,
            str(row.get("broad_cell_type", "")),
        ),
    )
    return best, bool(eligible)


def choose_modulator(
    rows: Iterable[Mapping[str, str]],
    desired_intervention: str,
) -> Mapping[str, str] | None:
    values = [row for row in rows if is_credible_modulator(row)]
    if not values:
        return None

    def key(row: Mapping[str, str]) -> tuple[int, int, int, float, str]:
        source = str(row.get("source", ""))
        has_structure = bool(str(row.get("canonical_smiles", "")).strip())
        score = optional_float(row.get("evidence_score")) or -1.0
        source_rank = 2 if source == "ChEMBL" else 1 if source.startswith("Open Targets") else 0
        return (
            int(modulator_direction_aligned(row, desired_intervention)),
            int(has_structure),
            source_rank,
            score,
            str(row.get("drug_name", "")),
        )

    return max(values, key=key)


def is_credible_modulator(row: Mapping[str, str]) -> bool:
    source = str(row.get("source", ""))
    if source == "ChEMBL" or source.startswith("Open Targets"):
        return True
    if source != "DGIdb":
        return False
    interaction_type = str(row.get("interaction_type", "")).strip()
    directionality = str(row.get("directionality", "")).strip().upper()
    return bool(interaction_type or directionality in {"ACTIVATING", "INHIBITORY"})


def modulator_direction_aligned(row: Mapping[str, str], desired_intervention: str) -> bool:
    directionality = str(row.get("directionality", "")).strip().upper()
    interaction_type = str(row.get("interaction_type", "")).strip().lower()
    if desired_intervention == "inhibit_or_reduce_expression":
        return directionality == "INHIBITORY" or "inhib" in interaction_type or "antagon" in interaction_type
    if desired_intervention == "increase_or_replace_expression":
        return directionality == "ACTIVATING" or "agon" in interaction_type or "activat" in interaction_type
    return False


def evidence_tier(row: Mapping[str, object]) -> str:
    if row["causal_complete"] and row["single_cell_pass"] and row["known_modulator"]:
        return "A"
    if row["mr_fdr_pass"] and row["coloc_pass"] and row["single_cell_pass"]:
        return "B"
    if int(row["evidence_modalities"]) >= 4:
        return "C"
    return "D"


def main() -> int:
    args = parse_args()
    root = args.repo_root.resolve()
    config = json.loads((root / "config" / "idd_analysis.json").read_text(encoding="utf-8"))
    weights = config["ranking"]["weights"]
    minimum_fraction = float(config["single_cell"]["minimum_cells_expressing_fraction"])

    candidates = read_tsv(
        root / "artifacts" / "mr" / "formal_candidate_genes.tsv",
        ("gene_symbol", "ensembl_gene_id"),
    )
    gene_locations = {
        row["hgnc_names"]: row
        for row in read_tsv(root / "data" / "processed" / "mr" / "druggable_genes.tsv")
        if row.get("hgnc_names")
    }
    bulk = grouped(read_tsv(root / "artifacts" / "bulk" / "druggable_bulk_candidates.tsv"))
    screen = grouped(read_tsv(root / "artifacts" / "mr" / "druggable_mr_screen.tsv"))
    forward = index_outcome(read_tsv(root / "artifacts" / "mr" / "mr_results.tsv"))
    coloc = index_outcome(read_tsv(root / "artifacts" / "mr" / "coloc_results.tsv"))
    smr = index_outcome(read_tsv(root / "artifacts" / "mr" / "smr_heidi_results.tsv"))
    steiger = index_outcome(read_tsv(root / "artifacts" / "mr" / "steiger_results.tsv"))
    reverse_rows = read_tsv(root / "artifacts" / "mr" / "reverse_mr_results.tsv")
    reverse = index_outcome(reverse_rows)
    localization = grouped(read_tsv(root / "artifacts" / "scrna" / "target_localization.tsv"))
    delivery = grouped(read_tsv(root / "artifacts" / "immune_drug" / "delivery_assessment.tsv"))
    interactions = grouped(read_tsv(root / "artifacts" / "immune_drug" / "drug_interactions.tsv"))
    gtex_path = root / "artifacts" / "mr" / "gtex_coloc_results.tsv"
    gtex = grouped(read_tsv(gtex_path)) if gtex_path.exists() else {}

    evidence_rows: list[dict[str, object]] = []
    for candidate in candidates:
        symbol = candidate["gene_symbol"]
        bulk_row = bulk.get(symbol, [None])[0]
        delivery_row = delivery.get(symbol, [None])[0]
        local_row, single_cell_pass = choose_localization(localization.get(symbol, []), minimum_fraction)
        outcome_ids = sorted(
            {
                str(row["outcome_id"])
                for row in screen.get(symbol, [])
                if row.get("outcome_id")
            }
        )
        outcome_summaries = []
        for outcome_id in outcome_ids:
            screen_row = next(
                (row for row in screen.get(symbol, []) if row.get("outcome_id") == outcome_id),
                None,
            )
            mr_row = forward.get((symbol, outcome_id))
            coloc_row = coloc.get((symbol, outcome_id))
            smr_row = smr.get((symbol, outcome_id))
            steiger_row = steiger.get((symbol, outcome_id))
            reverse_row = reverse.get((symbol, outcome_id))
            screen_fdr = optional_float(screen_row.get("fdr")) if screen_row else None
            pp4 = optional_float(coloc_row.get("pp4")) if coloc_row else None
            mr_pass = bool(screen_fdr is not None and screen_fdr < 0.05)
            coloc_pass = bool(pp4 is not None and pp4 >= float(config["genetics"]["coloc_pp4_min"]))
            smr_pass = pass_smr(smr_row)
            steiger_pass = bool(
                steiger_row
                and steiger_row.get("status") == "estimable"
                and truthy(steiger_row.get("direction_supported"))
            )
            reverse_fdr = optional_float(reverse_row.get("fdr")) if reverse_row else None
            reverse_no_signal = bool(
                reverse_row
                and reverse_row.get("status") == "estimable"
                and reverse_fdr is not None
                and reverse_fdr >= 0.05
            )
            outcome_summaries.append(
                {
                    "outcome_id": outcome_id,
                    "screen": screen_row,
                    "mr": mr_row,
                    "coloc": coloc_row,
                    "smr": smr_row,
                    "steiger": steiger_row,
                    "reverse": reverse_row,
                    "screen_fdr": screen_fdr,
                    "pp4": pp4,
                    "mr_pass": mr_pass,
                    "coloc_pass": coloc_pass,
                    "smr_pass": smr_pass,
                    "steiger_pass": steiger_pass,
                    "reverse_no_signal": reverse_no_signal,
                    "integrated": mr_pass and coloc_pass and smr_pass,
                }
            )
        outcome_summaries.sort(
            key=lambda item: (
                int(item["integrated"]),
                int(item["mr_pass"]) + int(item["coloc_pass"]) + int(item["smr_pass"]),
                int(item["steiger_pass"]),
                -(item["screen_fdr"] if item["screen_fdr"] is not None else 1.0),
                item["pp4"] if item["pp4"] is not None else -1.0,
                item["outcome_id"],
            ),
            reverse=True,
        )
        best = outcome_summaries[0] if outcome_summaries else {}
        mr_pass_outcomes = [item["outcome_id"] for item in outcome_summaries if item["mr_pass"]]
        coloc_pass_outcomes = [item["outcome_id"] for item in outcome_summaries if item["coloc_pass"]]
        smr_pass_outcomes = [item["outcome_id"] for item in outcome_summaries if item["smr_pass"]]
        integrated_outcomes = [item["outcome_id"] for item in outcome_summaries if item["integrated"]]
        steiger_pass = any(item["steiger_pass"] and item["mr_pass"] for item in outcome_summaries)
        reverse_no_signal = any(item["reverse_no_signal"] and item["mr_pass"] for item in outcome_summaries)
        reverse_signal_rows = [
            row
            for row in reverse_rows
            if row.get("gene_symbol") == symbol
            if row.get("status") == "estimable"
            and optional_float(row.get("fdr")) is not None
            and optional_float(row.get("fdr")) < 0.05
        ]
        reverse_signal_rows.sort(key=lambda row: optional_float(row.get("fdr")) or 1.0)
        reverse_signal_detected = bool(reverse_signal_rows)

        mr_row = best.get("mr") or {}
        coloc_row = best.get("coloc") or {}
        smr_row = best.get("smr") or {}
        reverse_row = best.get("reverse") or {}
        mr_beta = optional_float(mr_row.get("beta"))
        if mr_beta is None:
            desired_intervention = "not_estimable"
        elif mr_beta > 0:
            desired_intervention = "inhibit_or_reduce_expression"
        elif mr_beta < 0:
            desired_intervention = "increase_or_replace_expression"
        else:
            desired_intervention = "not_estimable"
        modulator = choose_modulator(interactions.get(symbol, []), desired_intervention)
        modulator_aligned = bool(modulator and modulator_direction_aligned(modulator, desired_intervention))
        directional_delivery_score = directional_delivery_fraction(
            delivery_row, desired_intervention, modulator_aligned
        )
        recommended_payload = recommended_local_payload(
            delivery_row, desired_intervention, modulator, modulator_aligned
        )
        gtex_rows = [
            row
            for row in gtex.get(symbol, [])
            if row.get("status") == "estimable" and optional_float(row.get("pp4")) is not None
        ]
        best_gtex = max(
            gtex_rows,
            key=lambda row: optional_float(row.get("pp4")) or -1.0,
            default=None,
        )
        gtex_pp4 = optional_float(best_gtex.get("pp4")) if best_gtex else None

        bulk_deg = bool(bulk_row and truthy(bulk_row.get("pass_primary")))
        wgcna_module = bool(bulk_row and truthy(bulk_row.get("in_selected_module")))
        known_modulator = bool(delivery_row and truthy(delivery_row.get("known_modulator")))
        delivery_label = str(delivery_row.get("local_delivery_feasibility", "")) if delivery_row else ""
        mr_fdr_pass = bool(mr_pass_outcomes)
        coloc_pass = bool(coloc_pass_outcomes)
        smr_heidi_pass = bool(smr_pass_outcomes)
        causal_complete = bool(integrated_outcomes)
        modality_flags = [
            bulk_deg,
            wgcna_module,
            mr_fdr_pass,
            coloc_pass,
            smr_heidi_pass,
            steiger_pass,
            reverse_no_signal,
            single_cell_pass,
            known_modulator,
            directional_delivery_score > 0,
        ]
        score_components = {
            "score_bulk_deg": float(weights["bulk_deg"]) * bulk_deg,
            "score_wgcna_module": float(weights["wgcna_module"]) * wgcna_module,
            "score_mr_fdr": float(weights["mr_fdr"]) * mr_fdr_pass,
            "score_coloc_pp4": float(weights["coloc_pp4"]) * coloc_pass,
            "score_smr_heidi": float(weights["smr_heidi"]) * smr_heidi_pass,
            "score_steiger": float(weights["steiger"]) * steiger_pass,
            "score_reverse_mr_no_evidence": float(weights["reverse_mr_no_evidence"]) * reverse_no_signal,
            "score_single_cell_localization": float(weights["single_cell_localization"]) * single_cell_pass,
            "score_known_modulator": float(weights["known_modulator"]) * known_modulator,
            "score_local_delivery_feasibility": float(weights["local_delivery_feasibility"]) * directional_delivery_score,
        }
        row: dict[str, object] = {
            **candidate,
            "chr_b37": gene_locations.get(symbol, {}).get("chr_b37", ""),
            "start_b37": gene_locations.get(symbol, {}).get("start_b37", ""),
            "end_b37": gene_locations.get(symbol, {}).get("end_b37", ""),
            "score": sum(score_components.values()),
            "evidence_modalities": sum(modality_flags),
            "causal_complete": causal_complete,
            "integrated_genetic_outcomes": ";".join(integrated_outcomes),
            "mr_fdr_pass": mr_fdr_pass,
            "mr_pass_outcomes": ";".join(mr_pass_outcomes),
            "coloc_pass": coloc_pass,
            "coloc_pass_outcomes": ";".join(coloc_pass_outcomes),
            "smr_heidi_pass": smr_heidi_pass,
            "smr_pass_outcomes": ";".join(smr_pass_outcomes),
            "steiger_pass": steiger_pass,
            "reverse_mr_no_signal_detected": reverse_no_signal,
            "reverse_mr_conflict": reverse_signal_detected,
            "reverse_mr_signal_outcomes": ";".join(str(row.get("outcome_id", "")) for row in reverse_signal_rows),
            "reverse_mr_min_fdr": reverse_signal_rows[0].get("fdr", "") if reverse_signal_rows else "",
            "bulk_deg_pass": bulk_deg,
            "wgcna_module_pass": wgcna_module,
            "single_cell_pass": single_cell_pass,
            "known_modulator": known_modulator,
            "local_delivery_feasibility": delivery_label or "not_assessed",
            "preferred_local_payload_options": delivery_row.get("preferred_local_payload_options", "") if delivery_row else "",
            "directional_delivery_feasibility": directional_delivery_score,
            "recommended_local_payload": recommended_payload,
            "desired_intervention": desired_intervention,
            "best_outcome_id": best.get("outcome_id", ""),
            "screen_mr_fdr": value_or_blank(best.get("screen_fdr")),
            "mr_method": mr_row.get("method", "not_estimable"),
            "mr_n_instruments": mr_row.get("n_instruments", ""),
            "mr_beta": mr_row.get("beta", ""),
            "mr_ci_lower": mr_row.get("ci_lower", ""),
            "mr_ci_upper": mr_row.get("ci_upper", ""),
            "mr_odds_ratio": mr_row.get("odds_ratio", ""),
            "mr_or_ci_lower": mr_row.get("or_ci_lower", ""),
            "mr_or_ci_upper": mr_row.get("or_ci_upper", ""),
            "mr_pvalue": mr_row.get("pvalue", ""),
            "mr_formal_set_fdr": mr_row.get("fdr", ""),
            "coloc_shared_variants": coloc_row.get("shared_variants", ""),
            "coloc_pp4": coloc_row.get("pp4", ""),
            "smr_beta": smr_row.get("beta_smr", ""),
            "smr_se": smr_row.get("se_smr", ""),
            "smr_pvalue": smr_row.get("p_smr", ""),
            "smr_fdr": smr_row.get("p_smr_fdr", ""),
            "heidi_pvalue": smr_row.get("p_heidi", ""),
            "heidi_n_snp": smr_row.get("n_snp_heidi", ""),
            "reverse_mr_status": reverse_row.get("status", "not_estimable"),
            "reverse_mr_fdr": reverse_row.get("fdr", ""),
            "eqtl_source": "eQTLGen whole blood full cis",
            "gtex_coloc_pass": bool(gtex_pp4 is not None and gtex_pp4 >= 0.8),
            "gtex_coloc_best_tissue": best_gtex.get("tissue", "") if best_gtex else "not_estimable",
            "gtex_coloc_best_outcome_id": best_gtex.get("outcome_id", "") if best_gtex else "not_estimable",
            "gtex_coloc_best_pp4": gtex_pp4 if gtex_pp4 is not None else "",
            "scrna_cell_type": local_row.get("broad_cell_type", "") if local_row else "not_detected",
            "scrna_tissue": local_row.get("tissue", "") if local_row else "",
            "scrna_cell_fraction_expressing": local_row.get("cell_fraction_expressing", "") if local_row else "",
            "scrna_donors": local_row.get("donors", "") if local_row else "",
            "top_modulator_name": modulator.get("drug_name", "") if modulator else "",
            "top_modulator_id": modulator.get("drug_id", "") if modulator else "",
            "top_modulator_source": modulator.get("source", "") if modulator else "",
            "top_modulator_activity": modulator.get("evidence_score", "") if modulator else "",
            "top_modulator_smiles": modulator.get("canonical_smiles", "") if modulator else "",
            "top_modulator_direction_aligned": modulator_aligned,
            **score_components,
        }
        row["evidence_tier"] = evidence_tier(row)
        evidence_rows.append(row)

    evidence_rows.sort(key=lambda row: str(row["gene_symbol"]))
    evidence_rows.sort(
        key=lambda row: (
            row["evidence_tier"] == "A",
            row["causal_complete"],
            float(row["score"]),
            row["single_cell_pass"],
            len(str(row["integrated_genetic_outcomes"]).split(";")) if row["integrated_genetic_outcomes"] else 0,
            float(row["coloc_pp4"] or -1),
            row["known_modulator"],
            row["directional_delivery_feasibility"],
            row["top_modulator_direction_aligned"],
            row["gtex_coloc_pass"],
        ),
        reverse=True,
    )
    for rank, row in enumerate(evidence_rows, start=1):
        row["rank"] = rank

    priority_count = 1 + int(config["ranking"]["required_backup_targets_max"])
    selected_rows, locus_alternatives = select_distinct_loci(
        evidence_rows,
        priority_count,
        int(config["ranking"]["priority_locus_separation_bp"]),
    )
    priority_rows = []
    for index, row in enumerate(selected_rows):
        priority_rows.append(
            {
                "priority_role": "primary" if index == 0 else f"backup_{index}",
                **row,
            }
        )

    output_dir = root / "artifacts" / "targets"
    write_tsv(output_dir / "target_evidence_matrix.tsv", evidence_rows)
    write_tsv(output_dir / "target_ranking.tsv", evidence_rows)
    write_tsv(output_dir / "priority_targets.tsv", priority_rows)
    write_tsv(output_dir / "locus_alternatives.tsv", locus_alternatives)
    alternative_symbols = {str(row["gene_symbol"]) for row in locus_alternatives}
    gtex_rows = list(priority_rows)
    gtex_rows.extend(
        {"priority_role": "locus_alternative", **row}
        for row in evidence_rows
        if str(row["gene_symbol"]) in alternative_symbols
    )
    write_tsv(output_dir / "gtex_targets.tsv", gtex_rows)
    summary = {
        "candidate_count": len(evidence_rows),
        "tier_counts": {
            tier: sum(row["evidence_tier"] == tier for row in evidence_rows)
            for tier in ("A", "B", "C", "D")
        },
        "causal_complete_count": sum(bool(row["causal_complete"]) for row in evidence_rows),
        "priority_genes": [row["gene_symbol"] for row in priority_rows],
        "locus_alternatives": locus_alternatives,
        "score_definition": "config/idd_analysis.json ranking weights; moderate delivery receives half weight",
        "absence_rule": "not_estimable fields receive no support score",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "ranking_summary.json").write_text(
        json.dumps(summary, ensure_ascii=True, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


def write_tsv(path: Path, rows: list[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


if __name__ == "__main__":
    raise SystemExit(main())
