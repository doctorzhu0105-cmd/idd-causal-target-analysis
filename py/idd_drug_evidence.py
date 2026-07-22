#!/usr/bin/env python3
"""Collect database-backed tractability and drug evidence for IDD candidates."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


ROOT = Path(__file__).resolve().parents[1]
OPEN_TARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
DGIDB_URL = "https://dgidb.org/api/graphql"
CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data"

OPEN_TARGETS_QUERY = """
query Target($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    targetClass { label level }
    tractability { label modality value }
    subcellularLocations { location source }
    chemicalProbes {
      id
      drugId
      drugFromSourceId
      isHighQuality
      mechanismOfAction
      probesDrugsScore
    }
    drugAndClinicalCandidates {
      count
      rows {
        id
        maxClinicalStage
        drug { id name drugType maximumClinicalStage description }
      }
    }
  }
}
"""

DGIDB_QUERY = """
query Interactions($genes: [String!]) {
  interactions(geneNames: $genes, first: 500) {
    totalCount
    nodes {
      id
      interactionScore
      gene { name conceptId }
      drug { name conceptId approved }
      interactionTypes { type directionality }
      publications { pmid }
      sources { sourceDbName }
    }
  }
}
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidates",
        type=Path,
        default=ROOT / "artifacts" / "mr" / "formal_candidate_genes.tsv",
    )
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--reuse-existing", action="store_true")
    return parser.parse_args()


def session() -> requests.Session:
    retry = Retry(
        total=4,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "POST"),
    )
    result = requests.Session()
    result.mount("https://", HTTPAdapter(max_retries=retry))
    result.headers.update({"User-Agent": "op-bioinfo-idd-target-dossier/1.0"})
    return result


def read_candidates(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"gene_symbol", "ensembl_gene_id"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Candidate file lacks required columns: {path}")
    result = []
    seen = set()
    for row in rows:
        gene = row["ensembl_gene_id"].split(".")[0]
        symbol = row["gene_symbol"].strip()
        if gene and symbol and gene not in seen:
            seen.add(gene)
            result.append(
                {
                    "gene_symbol": symbol,
                    "ensembl_gene_id": gene,
                    "druggability_tier": row.get("druggability_tier", ""),
                    "selection_reason": row.get("selection_reason", ""),
                }
            )
    return result


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def post_graphql(client: requests.Session, url: str, query: str, variables: dict[str, Any]) -> dict[str, Any]:
    response = client.post(url, json={"query": query, "variables": variables}, timeout=(20, 60))
    response.raise_for_status()
    payload = response.json()
    if payload.get("errors"):
        raise RuntimeError(json.dumps(payload["errors"], ensure_ascii=False))
    return payload["data"]


def fetch_open_targets(
    client: requests.Session,
    candidates: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summaries = []
    interactions = []
    for candidate in candidates:
        base = dict(candidate)
        try:
            target = post_graphql(
                client,
                OPEN_TARGETS_URL,
                OPEN_TARGETS_QUERY,
                {"ensemblId": candidate["ensembl_gene_id"]},
            ).get("target")
            if target is None:
                summaries.append({**base, "source": "Open Targets", "status": "target_not_found"})
                continue
            tractability = target.get("tractability") or []
            locations = target.get("subcellularLocations") or []
            clinical = target.get("drugAndClinicalCandidates") or {"count": 0, "rows": []}
            probes = target.get("chemicalProbes") or []
            summaries.append(
                {
                    **base,
                    "source": "Open Targets",
                    "status": "retrieved",
                    "approved_name": target.get("approvedName", ""),
                    "target_classes": join_values(item.get("label") for item in target.get("targetClass") or []),
                    "subcellular_locations": join_values(item.get("location") for item in locations),
                    "small_molecule_tractability_labels": join_values(
                        item.get("label") for item in tractability if item.get("modality") == "SM" and item.get("value")
                    ),
                    "antibody_tractability_labels": join_values(
                        item.get("label") for item in tractability if item.get("modality") == "AB" and item.get("value")
                    ),
                    "protein_degrader_tractability_labels": join_values(
                        item.get("label") for item in tractability if item.get("modality") == "PR" and item.get("value")
                    ),
                    "clinical_candidate_count": int(clinical.get("count") or 0),
                    "chemical_probe_count": len(probes),
                    "high_quality_chemical_probe_count": sum(bool(item.get("isHighQuality")) for item in probes),
                }
            )
            for row in clinical.get("rows") or []:
                drug = row.get("drug") or {}
                interactions.append(
                    {
                        **base,
                        "source": "Open Targets",
                        "drug_id": drug.get("id", ""),
                        "drug_name": drug.get("name", ""),
                        "drug_type": drug.get("drugType", ""),
                        "interaction_type": "clinical_candidate",
                        "directionality": "",
                        "approved": "",
                        "maximum_clinical_stage": row.get("maxClinicalStage", ""),
                        "evidence_score": "",
                        "references": "",
                    }
                )
            for probe in probes:
                interactions.append(
                    {
                        **base,
                        "source": "Open Targets chemical probe",
                        "drug_id": probe.get("drugId") or probe.get("drugFromSourceId") or probe.get("id", ""),
                        "drug_name": probe.get("drugFromSourceId") or probe.get("id", ""),
                        "drug_type": "chemical_probe",
                        "interaction_type": join_values(probe.get("mechanismOfAction") or []),
                        "directionality": "",
                        "approved": False,
                        "maximum_clinical_stage": "",
                        "evidence_score": probe.get("probesDrugsScore", ""),
                        "references": "",
                    }
                )
        except Exception as error:  # source failure must remain visible, not become a zero count
            summaries.append({**base, "source": "Open Targets", "status": "unavailable", "reason": str(error)})
    return summaries, interactions


def fetch_dgidb(client: requests.Session, candidates: list[dict[str, str]]) -> tuple[list[dict[str, Any]], dict[str, int]]:
    by_symbol = {row["gene_symbol"].upper(): row for row in candidates}
    rows = []
    counts = defaultdict_int()
    for start in range(0, len(candidates), 10):
        symbols = [row["gene_symbol"] for row in candidates[start : start + 10]]
        data = post_graphql(client, DGIDB_URL, DGIDB_QUERY, {"genes": symbols})
        connection = data.get("interactions") or {"nodes": []}
        for item in connection.get("nodes") or []:
            gene = (item.get("gene") or {}).get("name", "").upper()
            candidate = by_symbol.get(gene)
            if candidate is None:
                continue
            counts[candidate["gene_symbol"]] += 1
            interaction_types = item.get("interactionTypes") or []
            drug = item.get("drug") or {}
            rows.append(
                {
                    **candidate,
                    "source": "DGIdb",
                    "drug_id": drug.get("conceptId", ""),
                    "drug_name": drug.get("name", ""),
                    "drug_type": "",
                    "interaction_type": join_values(value.get("type") for value in interaction_types),
                    "directionality": join_values(value.get("directionality") for value in interaction_types),
                    "approved": drug.get("approved", ""),
                    "maximum_clinical_stage": "",
                    "evidence_score": item.get("interactionScore", ""),
                    "references": join_values(
                        [str(value.get("pmid")) for value in item.get("publications") or []]
                        + [value.get("sourceDbName") for value in item.get("sources") or []]
                    ),
                }
            )
    return rows, dict(counts)


def exact_chembl_targets(payload: dict[str, Any], symbol: str) -> list[dict[str, Any]]:
    result = []
    for target in payload.get("targets") or []:
        if target.get("organism") != "Homo sapiens" or target.get("target_type") != "SINGLE PROTEIN":
            continue
        synonyms = {
            synonym.get("component_synonym", "").upper()
            for component in target.get("target_components") or []
            for synonym in component.get("target_component_synonyms") or []
        }
        if symbol.upper() in synonyms:
            result.append(target)
    return result


def fetch_chembl(
    client: requests.Session,
    candidates: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summaries = []
    rows = []
    for candidate in candidates:
        try:
            response = client.get(
                f"{CHEMBL_URL}/target/search.json",
                params={"q": candidate["gene_symbol"], "limit": 20},
                timeout=(20, 60),
            )
            response.raise_for_status()
            targets = exact_chembl_targets(response.json(), candidate["gene_symbol"])
            chembl_ids = [target["target_chembl_id"] for target in targets]
            candidate_rows = []
            for target_id in chembl_ids:
                activity_response = client.get(
                    f"{CHEMBL_URL}/activity.json",
                    params={
                        "target_chembl_id": target_id,
                        "pchembl_value__gte": 6,
                        "limit": 100,
                        "order_by": "-pchembl_value",
                    },
                    timeout=(20, 90),
                )
                activity_response.raise_for_status()
                candidate_rows.extend(activity_response.json().get("activities") or [])
            candidate_rows.sort(
                key=lambda row: float(row.get("pchembl_value") or "-inf"),
                reverse=True,
            )
            for activity in candidate_rows[:20]:
                molecule_id = activity.get("molecule_chembl_id", "")
                rows.append(
                    {
                        **candidate,
                        "source": "ChEMBL",
                        "drug_id": molecule_id,
                        "drug_name": activity.get("molecule_pref_name") or molecule_id,
                        "drug_type": activity.get("modality", ""),
                        "interaction_type": activity.get("standard_type", ""),
                        "directionality": "",
                        "approved": "",
                        "maximum_clinical_stage": "",
                        "evidence_score": activity.get("pchembl_value", ""),
                        "references": activity.get("document_chembl_id", ""),
                        "target_chembl_id": activity.get("target_chembl_id", ""),
                        "assay_chembl_id": activity.get("assay_chembl_id", ""),
                        "standard_value": activity.get("standard_value", ""),
                        "standard_units": activity.get("standard_units", ""),
                        "canonical_smiles": activity.get("canonical_smiles", ""),
                        "standard_inchi_key": "",
                    }
                )
            summaries.append(
                {
                    **candidate,
                    "source": "ChEMBL",
                    "status": "retrieved",
                    "target_chembl_ids": join_values(chembl_ids),
                    "activity_count_pchembl_gte_6": len(candidate_rows),
                    "top_pchembl_value": candidate_rows[0].get("pchembl_value", "") if candidate_rows else "",
                }
            )
        except Exception as error:
            summaries.append({**candidate, "source": "ChEMBL", "status": "unavailable", "reason": str(error)})
    return summaries, rows


def aggregate_delivery(
    candidates: list[dict[str, str]],
    open_targets: list[dict[str, Any]],
    interactions: list[dict[str, Any]],
    dgidb_counts: Mapping[str, int],
) -> list[dict[str, Any]]:
    target_summary = {
        row["gene_symbol"]: row
        for row in open_targets
        if row.get("status") == "retrieved"
    }
    rows = []
    for candidate in candidates:
        symbol = candidate["gene_symbol"]
        target = target_summary.get(symbol, {})
        gene_interactions = [row for row in interactions if row.get("gene_symbol") == symbol]
        chembl_rows = [row for row in gene_interactions if row.get("source") == "ChEMBL"]
        small_molecule = bool(target.get("small_molecule_tractability_labels")) or bool(chembl_rows)
        antibody = bool(target.get("antibody_tractability_labels"))
        locations = str(target.get("subcellular_locations", ""))
        extracellular = any(token in locations.lower() for token in ("extracellular", "secreted", "cell surface", "plasma membrane"))
        known_modulator = any(is_credible_modulator(row) for row in gene_interactions)
        if small_molecule:
            feasibility = "high_for_local_small_molecule_or_siRNA_formulation"
            payloads = "small_molecule;siRNA"
        elif extracellular or antibody:
            feasibility = "moderate_for_local_protein_antibody_or_siRNA_formulation"
            payloads = "protein_or_antibody;siRNA"
        else:
            feasibility = "moderate_preclinical_for_siRNA_nanoparticle_hydrogel"
            payloads = "siRNA"
        rows.append(
            {
                **candidate,
                "known_modulator": known_modulator,
                "small_molecule_tractable": small_molecule,
                "antibody_tractable": antibody,
                "extracellular_or_surface_location": extracellular,
                "preferred_local_payload_options": payloads,
                "local_delivery_feasibility": feasibility,
                "open_targets_clinical_candidate_count": target.get("clinical_candidate_count", ""),
                "open_targets_chemical_probe_count": target.get("chemical_probe_count", ""),
                "dgidb_interaction_count": dgidb_counts.get(symbol, 0),
                "chembl_high_activity_rows_retained": len(chembl_rows),
                "delivery_claim_scope": "formulation feasibility only; efficacy and disc retention require wet-lab validation",
            }
        )
    return rows


def is_credible_modulator(row: Mapping[str, Any]) -> bool:
    source = str(row.get("source", ""))
    if source == "ChEMBL" or source.startswith("Open Targets"):
        return True
    if source != "DGIdb":
        return False
    interaction_type = str(row.get("interaction_type", "")).strip()
    directionality = str(row.get("directionality", "")).strip().upper()
    return bool(interaction_type or directionality in {"ACTIVATING", "INHIBITORY"})


def defaultdict_int() -> dict[str, int]:
    from collections import defaultdict

    return defaultdict(int)


def join_values(values: Iterable[Any]) -> str:
    unique = []
    seen = set()
    for value in values:
        if value is None:
            continue
        text = str(value).strip()
        if text and text not in seen:
            seen.add(text)
            unique.append(text)
    return ";".join(unique)


def write_tsv(path: Path, rows: list[Mapping[str, Any]]) -> None:
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


def main() -> None:
    args = parse_args()
    root = args.repo_root.resolve()
    candidate_path = args.candidates if args.candidates.is_absolute() else root / args.candidates
    candidates = read_candidates(candidate_path)
    output_dir = root / "artifacts" / "immune_drug"
    if args.reuse_existing:
        open_target_summary = read_tsv(output_dir / "open_targets_summary.tsv")
        chembl_summary = read_tsv(output_dir / "chembl_summary.tsv")
        interactions = read_tsv(output_dir / "drug_interactions.tsv")
        dgidb_counts = defaultdict_int()
        for row in interactions:
            if row.get("source") == "DGIdb":
                dgidb_counts[row["gene_symbol"]] += 1
        dgidb_status = "reused_existing_snapshot"
        dgidb_reason = ""
    else:
        client = session()
        open_target_summary, open_target_interactions = fetch_open_targets(client, candidates)
        try:
            dgidb_interactions, dgidb_counts = fetch_dgidb(client, candidates)
            dgidb_status = "retrieved"
            dgidb_reason = ""
        except Exception as error:
            dgidb_interactions, dgidb_counts = [], {}
            dgidb_status = "unavailable_fallback_used"
            dgidb_reason = str(error)
        chembl_summary, chembl_interactions = fetch_chembl(client, candidates)
        interactions = open_target_interactions + dgidb_interactions + chembl_interactions
    delivery = aggregate_delivery(candidates, open_target_summary, interactions, dgidb_counts)
    if not args.reuse_existing:
        write_tsv(output_dir / "open_targets_summary.tsv", open_target_summary)
        write_tsv(output_dir / "chembl_summary.tsv", chembl_summary)
        write_tsv(output_dir / "drug_interactions.tsv", interactions)
    write_tsv(output_dir / "delivery_assessment.tsv", delivery)
    previous_summary = {}
    summary_path = root / "logs" / "idd" / "drug_evidence_summary.json"
    if args.reuse_existing and summary_path.exists():
        previous_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary = {
        "generated_at_utc": previous_summary.get(
            "generated_at_utc", "2026-07-16T03:00:00+00:00"
        ),
        "candidate_count": len(candidates),
        "sources": {
            "Open Targets": OPEN_TARGETS_URL,
            "DGIdb": DGIDB_URL,
            "ChEMBL": CHEMBL_URL,
        },
        "dgidb_status": dgidb_status,
        "dgidb_reason": dgidb_reason,
        "interaction_rows": len(interactions),
        "known_modulator_candidates": sum(bool(row["known_modulator"]) for row in delivery),
        "delivery_rule_scope": "local formulation feasibility, not therapeutic efficacy",
    }
    log_path = summary_path
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
