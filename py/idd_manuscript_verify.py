#!/usr/bin/env python3
"""Verify manuscript citations, figures, derived results, and frozen-core integrity."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from PIL import Image


CORE_MANIFEST_SHA256 = "85610e4e52fbad961c8e6e8055ce32ec0a188dc093bdc7033a8671c4cf5fc7da"
REQUIRED_FILES = [
    "manuscript/idd_mr_manuscript.md",
    "manuscript/refs.bib",
    "manuscript/FIGURE_LEGENDS.md",
    "manuscript/SUPPLEMENTARY_MATERIALS.md",
    "manuscript/STROBE_MR_CHECKLIST.md",
    "manuscript/analysis_summary.json",
    "manuscript/tables/table_1_data_sources.tsv",
    "manuscript/tables/table_2_key_genetic_evidence.tsv",
    "manuscript/tables/table_3_translational_prioritization.tsv",
    "manuscript/supplement/mr_presso_results.tsv",
    "manuscript/supplement/mr_presso_summary.tsv",
    "manuscript/supplement/mr_egger_student_t_results.tsv",
    "manuscript/supplement/multiplicative_random_effects_ivw.tsv",
    "manuscript/supplement/selection_flow.tsv",
    "manuscript/supplement/enrichment_results.tsv",
    "manuscript/supplement/enrichment_queries.tsv",
    "manuscript/supplement/ppi_mapping.tsv",
    "manuscript/supplement/ppi_physical_edges.tsv",
    "manuscript/supplement/ppi_summary.tsv",
    "manuscript/supplement/external_analysis_metadata.json",
    "manuscript/supplement/wgcna_signed_modules.tsv",
    "manuscript/supplement/wgcna_signed_genes.tsv",
    "manuscript/supplement/wgcna_signed_soft_threshold.tsv",
    "manuscript/supplement/wgcna_signed_summary.json",
    "manuscript/supplement/ukb_saige_722_mr.tsv",
    "manuscript/supplement/ukb_saige_722_metadata.json",
    "manuscript/figures/figure_manifest.tsv",
    "HANDOFF/PAPER_READINESS.md",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--skip-online", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_doi(value: str) -> str:
    return value.strip().lower().replace("https://doi.org/", "").replace("http://dx.doi.org/", "")


def parse_bib_entries(text: str) -> list[dict[str, str]]:
    entries = []
    starts = list(re.finditer(r"(?m)^@(\w+)\{([^,]+),", text))
    for index, match in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        block = text[match.start() : end]
        doi_match = re.search(r"(?mi)^\s*doi\s*=\s*\{([^}]+)\}", block)
        title_match = re.search(r"(?mi)^\s*title\s*=\s*\{([^\n]+)\}", block)
        entries.append(
            {
                "entry_type": match.group(1),
                "citekey": match.group(2).strip(),
                "doi": normalize_doi(doi_match.group(1)) if doi_match else "",
                "title": title_match.group(1).rstrip("},").strip() if title_match else "",
            }
        )
    return entries


def audit_citations(repo: Path, skip_online: bool) -> tuple[list[dict[str, Any]], list[str]]:
    manuscript_paths = [
        repo / "manuscript" / "idd_mr_manuscript.md",
        repo / "manuscript" / "STROBE_MR_CHECKLIST.md",
    ]
    manuscript_text = "\n".join(path.read_text(encoding="utf-8") for path in manuscript_paths)
    bib_text = (repo / "manuscript" / "refs.bib").read_text(encoding="utf-8")
    entries = parse_bib_entries(bib_text)
    citekeys = {entry["citekey"] for entry in entries}
    used = set(re.findall(r"@([A-Za-z0-9_:-]+)", manuscript_text))
    failures = []
    missing = sorted(used - citekeys)
    if missing:
        failures.append(f"Missing bibliography citekeys: {missing}")
    duplicate_keys = sorted(
        key for key in citekeys if sum(entry["citekey"] == key for entry in entries) > 1
    )
    if duplicate_keys:
        failures.append(f"Duplicate bibliography citekeys: {duplicate_keys}")
    unresolved_patterns = [r"\[@NEW:", r"\[UNVERIFIED", r"\bTODO\b", r"\bTBD\b"]
    for pattern in unresolved_patterns:
        if re.search(pattern, manuscript_text, flags=re.IGNORECASE):
            failures.append(f"Unresolved manuscript marker matched {pattern}")

    audit_rows = []
    session = requests.Session()
    session.headers["User-Agent"] = "IDD manuscript DOI audit"
    for entry in entries:
        row: dict[str, Any] = {
            **entry,
            "cited_in_manuscript": entry["citekey"] in used,
            "crossref_status": "not_applicable" if not entry["doi"] else "not_checked",
            "crossref_title": "",
            "reason": "",
        }
        if entry["doi"] and not skip_online:
            try:
                response = session.get(
                    "https://api.crossref.org/works/" + requests.utils.quote(entry["doi"], safe=""),
                    timeout=45,
                )
                if response.status_code == 200:
                    message = response.json()["message"]
                    row["crossref_status"] = "verified"
                    row["crossref_title"] = (message.get("title") or [""])[0]
                else:
                    row["crossref_status"] = "failed"
                    row["reason"] = f"HTTP {response.status_code}"
                    failures.append(f"Crossref failed for {entry['citekey']}: HTTP {response.status_code}")
            except (requests.RequestException, ValueError, KeyError) as exc:
                row["crossref_status"] = "failed"
                row["reason"] = str(exc)
                failures.append(f"Crossref failed for {entry['citekey']}: {exc}")
        audit_rows.append(row)
    return audit_rows, failures


def audit_figures(repo: Path) -> tuple[list[dict[str, Any]], list[str]]:
    manifest_path = repo / "manuscript" / "figures" / "figure_manifest.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t")
    rows = []
    failures = []
    for item in manifest.to_dict(orient="records"):
        path = repo / item["path"]
        row = {**item, "exists": path.exists(), "bytes": path.stat().st_size if path.exists() else 0}
        if not path.exists() or path.stat().st_size == 0:
            failures.append(f"Missing or empty figure: {item['path']}")
            rows.append(row)
            continue
        suffix = path.suffix.lower()
        if suffix in (".png", ".tif", ".tiff"):
            with Image.open(path) as image:
                dpi = image.info.get("dpi", (0, 0))
                row["width_px"] = image.width
                row["height_px"] = image.height
                row["dpi_x"] = float(dpi[0]) if dpi else 0.0
                row["dpi_y"] = float(dpi[1]) if dpi else 0.0
                if min(row["dpi_x"], row["dpi_y"]) < 295:
                    failures.append(f"Raster figure below 300 dpi tolerance: {item['path']} ({dpi})")
        rows.append(row)
    expected_main = {"Figure 1", "Figure 2", "Figure 2B", "Figure 3", "Figure 4", "Figure 5"}
    observed = set(manifest["figure_id"].astype(str))
    missing_main = sorted(expected_main - observed)
    if missing_main:
        failures.append(f"Figure manifest missing main entries: {missing_main}")
    return rows, failures


def check_close(actual: float, expected: float, tolerance: float, label: str, failures: list[str]) -> None:
    if abs(actual - expected) > tolerance:
        failures.append(f"{label}: expected {expected}, observed {actual}")


def audit_results(repo: Path) -> tuple[dict[str, Any], list[str]]:
    failures = []
    core = repo / "artifacts" / "core"
    manifest_hash = sha256(core / "manifest.tsv")
    if manifest_hash != CORE_MANIFEST_SHA256:
        failures.append(f"Frozen core manifest drift: {manifest_hash}")
    manifest_rows = len(pd.read_csv(core / "manifest.tsv", sep="\t"))
    if manifest_rows != 57:
        failures.append(f"Frozen core manifest row count: expected 57, observed {manifest_rows}")

    genetic = pd.read_csv(repo / "manuscript" / "tables" / "table_2_key_genetic_evidence.tsv", sep="\t")
    akr = genetic[
        genetic["gene_symbol"].eq("AKR1C2")
        & genetic["outcome_id"].eq("M13_INTERVERTEB")
    ].iloc[0]
    check_close(float(akr["odds_ratio"]), 0.934064704476, 1e-10, "AKR1C2 FinnGen OR", failures)
    check_close(float(akr["screen_fdr"]), 0.00332746319302, 1e-12, "AKR1C2 screen FDR", failures)
    check_close(float(akr["pp4"]), 0.941246566595, 1e-10, "AKR1C2 FinnGen PP4", failures)
    check_close(
        float(akr["p_smr_bonferroni"]),
        0.00130692611,
        1e-12,
        "AKR1C2 SMR Bonferroni P",
        failures,
    )
    if int(akr["smr_bonferroni_test_count"]) != 71:
        failures.append("AKR1C2 SMR Bonferroni denominator is not 71 tested probes")
    check_close(float(akr["p_heidi"]), 0.4465942, 1e-8, "AKR1C2 HEIDI P", failures)
    nt5e_direct = genetic[
        genetic["gene_symbol"].eq("NT5E")
        & genetic["outcome_id"].eq("M13_INTERVERTEB")
    ].iloc[0]
    check_close(
        float(nt5e_direct["p_smr_bonferroni"]),
        0.0653940388,
        1e-11,
        "NT5E direct-outcome SMR Bonferroni P",
        failures,
    )
    observed_counts = set(
        pd.to_numeric(genetic["smr_bonferroni_test_count"], errors="coerce").dropna().astype(int)
    )
    if observed_counts != {71}:
        failures.append(f"Unexpected SMR Bonferroni test counts: {sorted(observed_counts)}")

    corrected_egger = pd.read_csv(
        repo / "manuscript" / "supplement" / "mr_egger_student_t_results.tsv", sep="\t"
    )
    akr_egger = corrected_egger[
        corrected_egger["gene_symbol"].eq("AKR1C2")
        & corrected_egger["outcome_id"].eq("M13_INTERVERTEB")
    ].iloc[0]
    check_close(float(akr_egger["pvalue"]), 0.180998087664, 1e-11, "AKR1C2 Egger t P", failures)
    p2_egger = corrected_egger[
        corrected_egger["gene_symbol"].eq("P2RY13")
        & corrected_egger["outcome_id"].eq("M13_LOWBACKPAIN")
    ].iloc[0]
    check_close(
        float(p2_egger["egger_intercept_pvalue"]),
        0.254146937846,
        1e-11,
        "P2RY13 Egger intercept t P",
        failures,
    )

    random_effects = pd.read_csv(
        repo / "manuscript" / "supplement" / "multiplicative_random_effects_ivw.tsv",
        sep="\t",
    )
    akr1c1_mre = random_effects[
        random_effects["gene_symbol"].eq("AKR1C1")
        & random_effects["outcome_id"].eq("M13_INTERVERTEB")
    ].iloc[0]
    check_close(
        float(akr1c1_mre["pvalue"]),
        0.180267029521,
        1e-11,
        "AKR1C1 multiplicative random-effects P",
        failures,
    )

    selection = pd.read_csv(repo / "manuscript" / "supplement" / "selection_flow.tsv", sep="\t")
    selection_counts = dict(zip(selection["selection_step"], selection["gene_count"]))
    expected_selection = {
        "primary_threshold_af_deg": 40,
        "druggable_af_deg_entering_bulk_route": 12,
        "corrected_signed_af_non_grey_modules": 31,
        "corrected_signed_af_fdr_modules": 0,
        "top_corrected_signed_af_module_genes": 71,
        "bulk_route_genes_in_top_signed_af_module": 2,
        "unique_genes_passing_lead_eqtl_mr_screen": 64,
        "overlap_mr_screen_and_bulk_route": 0,
        "frozen_formal_candidate_set": 76,
    }
    if selection_counts != expected_selection:
        failures.append(f"Selection-flow counts differ: {selection_counts}")
    translation = pd.read_csv(
        repo / "manuscript" / "tables" / "table_3_translational_prioritization.tsv",
        sep="\t",
    ).set_index("gene_symbol")
    if translation.loc["AKR1C2", "preferred_local_payload_options"] != "mrna_nanoparticle_in_hydrogel":
        failures.append("AKR1C2 direction-matched payload is not mRNA nanoparticle in hydrogel")
    if translation.loc["NT5E", "preferred_local_payload_options"] != "mrna_or_recombinant_protein_in_hydrogel":
        failures.append("NT5E direction-matched payload is not mRNA or recombinant protein")
    signed_wgcna = pd.read_csv(
        repo / "manuscript" / "supplement" / "wgcna_signed_modules.tsv",
        sep="\t",
    )
    if set(signed_wgcna["network_type"].astype(str)) != {"signed"}:
        failures.append("Corrected WGCNA does not report signed network type")
    if set(signed_wgcna["tom_type"].astype(str)) != {"signed"}:
        failures.append("Corrected WGCNA does not report signed TOM type")
    signed_genes = pd.read_csv(
        repo / "manuscript" / "supplement" / "wgcna_signed_genes.tsv",
        sep="\t",
    )
    if set(signed_genes["network_type"].astype(str)) != {"signed"}:
        failures.append("Corrected WGCNA gene assignments do not report signed network type")
    signed_soft_threshold = pd.read_csv(
        repo / "manuscript" / "supplement" / "wgcna_signed_soft_threshold.tsv",
        sep="\t",
    )
    selected_soft_thresholds = {"AF": 16, "NP": 12}
    for tissue, power in selected_soft_thresholds.items():
        selected_row = signed_soft_threshold[
            signed_soft_threshold["tissue"].eq(tissue)
            & pd.to_numeric(signed_soft_threshold["Power"], errors="coerce").eq(power)
        ]
        if len(selected_row) != 1:
            failures.append(
                f"Corrected WGCNA soft-threshold table lacks one {tissue} power-{power} row"
            )
    af_signed = signed_wgcna[
        signed_wgcna["tissue"].eq("AF") & signed_wgcna["module"].astype(str).ne("0")
    ]
    if int((af_signed["non_grey_bh_fdr"] < 0.05).sum()) != 0:
        failures.append("Corrected signed AF WGCNA unexpectedly has an FDR-significant module")
    af_top_module = af_signed.sort_values("adjusted_p_value").iloc[0]
    check_close(
        float(af_top_module["non_grey_bh_fdr"]),
        0.0518434054151,
        1e-11,
        "Top corrected signed AF module BH-FDR",
        failures,
    )
    if str(af_top_module["module"]) != "23":
        failures.append("Top corrected signed AF module is not M23")
    af_top_genes = signed_genes[
        signed_genes["tissue"].eq("AF")
        & signed_genes["module"].astype(str).eq(str(af_top_module["module"]))
    ]
    if af_top_genes["gene_symbol"].nunique() != 71:
        failures.append("Top corrected signed AF module does not contain 71 unique genes")

    signed_summary_path = repo / "manuscript" / "supplement" / "wgcna_signed_summary.json"
    signed_summary_file = json.loads(signed_summary_path.read_text(encoding="utf-8"))
    if signed_summary_file.get("network_type") != "signed" or signed_summary_file.get("tom_type") != "signed":
        failures.append("Corrected WGCNA summary does not report signed network and TOM")
    if signed_summary_file.get("af_soft_threshold_power") != 16:
        failures.append("Corrected WGCNA summary AF soft-threshold power is not 16")
    if signed_summary_file.get("np_soft_threshold_power") != 12:
        failures.append("Corrected WGCNA summary NP soft-threshold power is not 12")

    supplementary_text = (
        repo / "manuscript" / "SUPPLEMENTARY_MATERIALS.md"
    ).read_text(encoding="utf-8")
    expected_wgcna_index = {
        "Table S21": "supplement/wgcna_signed_modules.tsv",
        "Table S22": "supplement/wgcna_signed_genes.tsv",
        "Table S23": "supplement/wgcna_signed_soft_threshold.tsv",
    }
    for table_id, relative_path in expected_wgcna_index.items():
        expected_entry = f"| {table_id} | `{relative_path}` |"
        if expected_entry not in supplementary_text:
            failures.append(f"Supplementary index lacks {table_id} signed WGCNA entry")
    legacy_wgcna_name = "wgcna_multiple_testing_sensitivity.tsv"
    if legacy_wgcna_name in supplementary_text:
        failures.append("Supplementary index still references the legacy WGCNA sensitivity table")
    if (repo / "manuscript" / "supplement" / legacy_wgcna_name).exists():
        failures.append("Legacy WGCNA sensitivity table remains in the supplementary bundle")

    ukb = pd.read_csv(repo / "manuscript" / "supplement" / "ukb_saige_722_mr.tsv", sep="\t")
    ukb_akr = ukb[
        ukb["gene_symbol"].eq("AKR1C2")
        & ukb["analysis_mode"].eq("primary_eaf_qc")
    ].iloc[0]
    check_close(float(ukb_akr["odds_ratio"]), 0.876802737798, 1e-10, "AKR1C2 UKB OR", failures)
    check_close(float(ukb_akr["fdr"]), 0.0107106169251, 1e-10, "AKR1C2 UKB FDR", failures)
    ukb_sensitivity = pd.read_csv(
        repo / "manuscript" / "supplement" / "ukb_saige_722_sensitivity.tsv", sep="\t"
    )
    ukb_egger = ukb_sensitivity[
        ukb_sensitivity["gene_symbol"].eq("AKR1C2")
        & ukb_sensitivity["analysis_mode"].eq("primary_eaf_qc")
        & ukb_sensitivity["method"].eq("mr_egger")
    ].iloc[0]
    check_close(float(ukb_egger["pvalue"]), 0.293340972814, 1e-11, "AKR1C2 UKB Egger t P", failures)
    akr_weighted = ukb_sensitivity[
        ukb_sensitivity["gene_symbol"].eq("AKR1C2")
        & ukb_sensitivity["method"].eq("weighted_median")
    ]
    if akr_weighted["pvalue"].astype(float).nunique() != 1:
        failures.append("Identical AKR1C2 UKB instrument sets produced different weighted-median P values")
    check_close(
        float(akr_weighted.iloc[0]["pvalue"]),
        0.00495777217819,
        1e-12,
        "AKR1C2 UKB weighted-median P",
        failures,
    )
    ukb_random_effects = ukb_sensitivity[
        ukb_sensitivity["analysis_mode"].eq("primary_eaf_qc")
        & ukb_sensitivity["method"].eq("multiplicative_random_effects_ivw")
    ].set_index("gene_symbol")
    check_close(
        float(ukb_random_effects.loc["AKR1C1", "pvalue"]),
        0.0559596266245,
        1e-12,
        "AKR1C1 UKB multiplicative random-effects P",
        failures,
    )
    check_close(
        float(ukb_random_effects.loc["AKR1C2", "pvalue"]),
        0.00345766786923,
        1e-12,
        "AKR1C2 UKB multiplicative random-effects P",
        failures,
    )
    instruments = pd.read_csv(
        repo / "manuscript" / "supplement" / "ukb_saige_722_instruments.tsv", sep="\t"
    )
    if len(instruments) != 14 or instruments["primary_included"].astype(bool).sum() != 13:
        failures.append("UKB instrument counts do not equal 14 queried and 13 primary")
    rs470 = instruments[instruments["snp"].eq("rs4707219")].iloc[0]
    if rs470["status"] != "sensitivity_only":
        failures.append("rs4707219 is not marked sensitivity_only")

    presso = pd.read_csv(repo / "manuscript" / "supplement" / "mr_presso_summary.tsv", sep="\t")
    if int(presso.iloc[0]["estimable_combinations"]) != 12:
        failures.append("MR-PRESSO estimable count is not 12")
    if int(presso.iloc[0]["bh_significant_global_tests"]) != 0:
        failures.append("MR-PRESSO significant global test count is not 0")

    enrichment = pd.read_csv(
        repo / "manuscript" / "supplement" / "enrichment_queries.tsv", sep="\t"
    ).set_index("query_name")
    if set(enrichment.index.astype(str)) != {"AF_DEG", "AF_WGCNA_SIGNED"}:
        failures.append(f"Unexpected enrichment query set: {sorted(enrichment.index.astype(str))}")
    if int(float(enrichment.loc["AF_DEG", "foreground_size"])) != 40:
        failures.append("AF DEG enrichment foreground is not 40")
    if int(float(enrichment.loc["AF_DEG", "mapped_foreground_size"])) != 39:
        failures.append("AF DEG enrichment mapped foreground is not 39")
    if int(float(enrichment.loc["AF_WGCNA_SIGNED", "foreground_size"])) != 0:
        failures.append("Corrected signed WGCNA enrichment foreground is not empty")
    if enrichment.loc["AF_WGCNA_SIGNED", "status"] != "not_estimable":
        failures.append("Corrected signed WGCNA enrichment is not marked not_estimable")
    expected_wgcna_ora_reason = (
        "no_non_grey_module_passed_bh_fdr_0.05_in_corrected_signed_reanalysis"
    )
    if enrichment.loc["AF_WGCNA_SIGNED", "reason"] != expected_wgcna_ora_reason:
        failures.append("Corrected signed WGCNA enrichment reason is unexpected")
    enrichment_results = pd.read_csv(
        repo / "manuscript" / "supplement" / "enrichment_results.tsv", sep="\t"
    )
    enrichment_result_names = set(
        enrichment_results.get("query_name", pd.Series(dtype=str)).astype(str)
    )
    if not enrichment_result_names <= {"AF_DEG"}:
        failures.append(
            f"Non-DEG enrichment results remain: {sorted(enrichment_result_names - {'AF_DEG'})}"
        )

    ppi_summary = pd.read_csv(
        repo / "manuscript" / "supplement" / "ppi_summary.tsv", sep="\t"
    ).set_index("query_name")
    if set(ppi_summary.index.astype(str)) != {"AF_DEG", "CAUSAL_COMPLETE"}:
        failures.append(f"Unexpected PPI query set: {sorted(ppi_summary.index.astype(str))}")
    expected_ppi_inputs = {"AF_DEG": 40, "CAUSAL_COMPLETE": 11}
    for query_name, expected_count in expected_ppi_inputs.items():
        if int(ppi_summary.loc[query_name, "input_gene_count"]) != expected_count:
            failures.append(f"{query_name} PPI input count is not {expected_count}")
    for relative_path in ("ppi_mapping.tsv", "ppi_physical_edges.tsv"):
        ppi_frame = pd.read_csv(repo / "manuscript" / "supplement" / relative_path, sep="\t")
        ppi_query_names = set(ppi_frame.get("query_name", pd.Series(dtype=str)).astype(str))
        unexpected_ppi_names = ppi_query_names - {"AF_DEG", "CAUSAL_COMPLETE"}
        if unexpected_ppi_names:
            failures.append(
                f"Unexpected PPI queries remain in {relative_path}: {sorted(unexpected_ppi_names)}"
            )
    external_metadata = json.loads(
        (
            repo / "manuscript" / "supplement" / "external_analysis_metadata.json"
        ).read_text(encoding="utf-8")
    )
    metadata_enrichment_names = {
        str(row.get("query_name"))
        for row in external_metadata.get("gprofiler", {}).get("query_status", [])
    }
    if metadata_enrichment_names != {"AF_DEG", "AF_WGCNA_SIGNED"}:
        failures.append("External metadata enrichment query set is unexpected")
    metadata_ppi_names = {
        str(row.get("query_name"))
        for row in external_metadata.get("string", {}).get("query_status", [])
    }
    if metadata_ppi_names != {"AF_DEG", "CAUSAL_COMPLETE"}:
        failures.append("External metadata PPI query set is unexpected")

    summary = json.loads((repo / "manuscript" / "analysis_summary.json").read_text(encoding="utf-8"))
    if summary.get("core_manifest_sha256") != CORE_MANIFEST_SHA256:
        failures.append("Manuscript analysis summary core hash differs")
    smr_counts = summary.get("genetics", {}).get("smr_bonferroni_test_counts_by_outcome", {})
    if set(smr_counts.values()) != {71} or len(smr_counts) != 3:
        failures.append(f"Analysis summary SMR Bonferroni counts are unexpected: {smr_counts}")
    if summary.get("genetics", {}).get("heidi_estimable") != 207:
        failures.append("Analysis summary HEIDI estimable count is not 207")
    signed_summary = summary.get("wgcna_signed_reanalysis", {})
    if signed_summary.get("network_type") != "signed":
        failures.append("Analysis summary corrected WGCNA network type is not signed")
    if signed_summary.get("af_fdr_significant_module_count") != 0:
        failures.append("Analysis summary corrected AF WGCNA module count is not zero")
    check_close(
        float(signed_summary.get("af_top_module_p", float("nan"))),
        0.00360339384516422,
        1e-15,
        "Analysis summary corrected AF top-module P",
        failures,
    )
    check_close(
        float(signed_summary.get("af_top_module_fdr", float("nan"))),
        0.0518434054151057,
        1e-15,
        "Analysis summary corrected AF top-module FDR",
        failures,
    )
    bulk_summary = summary.get("bulk", {})
    stale_bulk_wgcna_fields = {"selected_af_modules", "module_selection_rule", "wgcna_power"}
    retained_stale_fields = sorted(stale_bulk_wgcna_fields & set(bulk_summary))
    if retained_stale_fields:
        failures.append(
            "Analysis summary retains superseded frozen WGCNA fields: "
            + ", ".join(retained_stale_fields)
        )
    if "corrected WGCNA is reported separately" not in str(bulk_summary.get("scope", "")):
        failures.append("Analysis summary bulk scope does not separate corrected WGCNA")
    summary_enrichment_names = {
        str(row.get("query_name")) for row in summary.get("enrichment_queries", [])
    }
    if summary_enrichment_names != {"AF_DEG", "AF_WGCNA_SIGNED"}:
        failures.append("Analysis summary enrichment query set is unexpected")
    summary_ppi_names = {str(row.get("query_name")) for row in summary.get("ppi_queries", [])}
    if summary_ppi_names != {"AF_DEG", "CAUSAL_COMPLETE"}:
        failures.append("Analysis summary PPI query set is unexpected")
    for source_name, rows in (
        ("enrichment", summary.get("enrichment_queries", [])),
        ("PPI", summary.get("ppi_queries", [])),
    ):
        for row in rows:
            if row.get("status") == "estimable" and row.get("reason") != "none":
                failures.append(
                    f"Analysis summary {source_name} query has estimable status with reason "
                    f"{row.get('reason')!r}"
                )
    if (repo / "manuscript" / "figures" / "figure_s6_source_umap.png").exists():
        failures.append("Excluded frozen source UMAP remains in the manuscript figure package")
    report_text = "\n".join(
        (repo / relative).read_text(encoding="utf-8")
        for relative in (
            "manuscript/idd_mr_manuscript.md",
            "manuscript/FIGURE_LEGENDS.md",
            "manuscript/SUPPLEMENTARY_MATERIALS.md",
            "HANDOFF/PAPER_READINESS.md",
        )
    )
    for stale_label in ("SMR FDR", "SMR -log10 FDR"):
        if stale_label in report_text:
            failures.append(f"Stale SMR multiple-testing label remains: {stale_label}")
    for stale_term in ("reverse MR", "reverse-MR conflict", "MR FDR |"):
        if stale_term in report_text:
            failures.append(f"Stale or overbroad manuscript term remains: {stale_term}")
    for stale_wgcna_pattern in (
        r"AF_WGCNA_MODULE_4",
        r"0\.0674",
        r"0\.0379",
        r"module[- ]4 (?:comprised|contained|produced|was retained|member|genes?)",
        r"identifies AF module 4",
    ):
        if re.search(stale_wgcna_pattern, report_text, flags=re.I):
            failures.append(
                f"Stale WGCNA module-4 inference remains: {stale_wgcna_pattern}"
            )
    return {
        "core_manifest_sha256": manifest_hash,
        "core_manifest_rows": manifest_rows,
        "akr1c2_finngen_or": float(akr["odds_ratio"]),
        "akr1c2_smr_bonferroni_p": float(akr["p_smr_bonferroni"]),
        "akr1c2_ukb_or": float(ukb_akr["odds_ratio"]),
        "mr_presso_estimable": int(presso.iloc[0]["estimable_combinations"]),
    }, failures


def main() -> None:
    args = parse_args()
    repo = Path(args.repo_root).resolve()
    output = repo / "manuscript" / "qc"
    output.mkdir(parents=True, exist_ok=True)
    failures = []
    for relative in REQUIRED_FILES:
        path = repo / relative
        if not path.exists() or path.stat().st_size == 0:
            failures.append(f"Missing required file: {relative}")

    citation_rows, citation_failures = audit_citations(repo, args.skip_online)
    figure_rows, figure_failures = audit_figures(repo)
    result_summary, result_failures = audit_results(repo)
    failures.extend(citation_failures)
    failures.extend(figure_failures)
    failures.extend(result_failures)
    pd.DataFrame(citation_rows).to_csv(output / "citation_audit.tsv", sep="\t", index=False)
    pd.DataFrame(figure_rows).to_csv(output / "figure_audit.tsv", sep="\t", index=False)

    manuscript = (repo / "manuscript" / "idd_mr_manuscript.md").read_text(encoding="utf-8")
    word_count = len(re.findall(r"\b[\w'-]+\b", re.sub(r"```.*?```", "", manuscript, flags=re.S)))
    payload = {
        "status": "verified" if not failures else "failed",
        "failure_count": len(failures),
        "failures": failures,
        "manuscript_word_count_including_tables_and_declarations": word_count,
        "bibliography_entries": len(citation_rows),
        "cited_entries": sum(bool(row["cited_in_manuscript"]) for row in citation_rows),
        "figure_files": len(figure_rows),
        "result_checks": result_summary,
    }
    (output / "verification.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
