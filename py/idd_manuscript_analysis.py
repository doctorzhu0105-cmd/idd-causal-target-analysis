#!/usr/bin/env python3
"""Build manuscript tables and figures from the frozen IDD evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import networkx as nx
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from PIL import Image

from idd_genetics_utils import ivw_multiplicative_random_effects, mr_egger
from idd_ukb_replication import plot_replication_axis


SEED = 20260716
CORE_MANIFEST_SHA256 = "85610e4e52fbad961c8e6e8055ce32ec0a188dc093bdc7033a8671c4cf5fc7da"
FULL_WIDTH_MM = 170.0
FULL_WIDTH_IN = FULL_WIDTH_MM / 25.4
OUTPUT_DPI = 300
OUTCOME_LABELS = {
    "M13_INTERVERTEB": "Intervertebral disc disorders",
    "M13_LOWBACKPAIN": "Low back pain",
    "M13_SCIATICA": "Sciatica with lumbago",
}
PAPER_GENE_ORDER = ["AKR1C2", "NT5E", "CCL28", "P2RY13", "AKR1C1"]
ENRICHMENT_QUERY_NAMES = {"AF_DEG", "AF_WGCNA_SIGNED"}
PPI_QUERY_NAMES = {"AF_DEG", "CAUSAL_COMPLETE"}
LEGACY_WGCNA_SUPPLEMENT = "wgcna_multiple_testing_sensitivity.tsv"
PAPER_ROLES = {
    "AKR1C2": "genetic lead for the primary disc-disorder phenotype",
    "NT5E": "low-back-pain-associated candidate; no independent outcome-cohort support for the primary disc-disorder phenotype",
    "CCL28": "single-instrument candidate for the primary disc-disorder phenotype",
    "P2RY13": "low-back-pain-associated candidate",
    "AKR1C1": "AKR1C2-locus alternative",
}
DIRECTION_MATCHED_PAYLOADS = {
    "AKR1C2": "mrna_nanoparticle_in_hydrogel",
    "NT5E": "mrna_or_recombinant_protein_in_hydrogel",
    "CCL28": "sirna_nanoparticle_in_hydrogel",
    "P2RY13": "sirna_or_direction_aligned_small_molecule_in_hydrogel",
    "AKR1C1": "mrna_nanoparticle_in_hydrogel",
}
PALETTE = {
    "navy": "#0072B2",
    "blue": "#56B4E9",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "gray": "#6B6B6B",
    "light": "#E6E6E6",
    "ink": "#222222",
}

OUTCOME_COLORS = {
    "M13_INTERVERTEB": PALETTE["navy"],
    "M13_LOWBACKPAIN": PALETTE["orange"],
    "M13_SCIATICA": PALETTE["green"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--refresh-external", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="not_estimable")


def json_records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    return frame.astype(object).where(pd.notna(frame), None).to_dict(orient="records")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_figure(fig: plt.Figure, stem: Path) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix in (".pdf", ".png", ".tiff"):
        path = stem.with_suffix(suffix)
        kwargs: dict[str, Any] = {"bbox_inches": "tight", "facecolor": "white"}
        if suffix != ".pdf":
            kwargs["dpi"] = OUTPUT_DPI
        if suffix == ".tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def set_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.titleweight": "bold",
        }
    )


def add_panel_label(ax: plt.Axes, label: str, *, x: float = -0.14, y: float = 1.05) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )


def format_compact_value(value: float) -> str:
    if pd.isna(value):
        return "NE"
    if value == 0:
        return "0"
    if abs(value) < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def load_inputs(repo: Path) -> dict[str, Any]:
    core = repo / "artifacts" / "core"
    manifest = core / "manifest.tsv"
    actual_hash = sha256(manifest)
    if actual_hash != CORE_MANIFEST_SHA256:
        raise RuntimeError(
            f"Frozen core manifest drifted: expected {CORE_MANIFEST_SHA256}, got {actual_hash}"
        )
    files = {
        "manifest": read_tsv(manifest),
        "deg": read_tsv(core / "bulk" / "gse70362_deg.tsv"),
        "wgcna_genes": read_tsv(
            repo / "manuscript" / "supplement" / "wgcna_signed_genes.tsv"
        ),
        "wgcna_modules": read_tsv(
            repo / "manuscript" / "supplement" / "wgcna_signed_modules.tsv"
        ),
        "bulk_candidates": read_tsv(core / "bulk" / "druggable_bulk_candidates.tsv"),
        "sample_metadata": read_tsv(core / "bulk" / "gse70362_sample_metadata.tsv"),
        "priority": read_tsv(core / "targets" / "priority_targets.tsv"),
        "evidence": read_tsv(core / "targets" / "target_evidence_matrix.tsv"),
        "locus": read_tsv(core / "targets" / "locus_alternatives.tsv"),
        "mr": read_tsv(core / "genetics" / "mr_results.tsv"),
        "mr_screen": read_tsv(core / "genetics" / "druggable_mr_screen.tsv"),
        "formal_candidates": read_tsv(core / "genetics" / "formal_candidate_genes.tsv"),
        "mr_sensitivity": read_tsv(core / "genetics" / "mr_sensitivity_results.tsv"),
        "harmonized": read_tsv(core / "genetics" / "harmonized_instruments.tsv"),
        "loo": read_tsv(core / "genetics" / "leave_one_out_results.tsv"),
        "coloc": read_tsv(core / "genetics" / "coloc_results.tsv"),
        "gtex_coloc": read_tsv(core / "genetics" / "gtex_coloc_results.tsv"),
        "smr": read_tsv(core / "genetics" / "smr_heidi_results.tsv"),
        "steiger": read_tsv(core / "genetics" / "steiger_results.tsv"),
        "reverse": read_tsv(core / "genetics" / "reverse_mr_results.tsv"),
        "scrna": read_tsv(core / "scrna" / "target_localization.tsv"),
        "clusters": read_tsv(core / "scrna" / "gse230809_cluster_annotations.tsv"),
        "drug": read_tsv(core / "drug_delivery" / "drug_interactions.tsv"),
        "delivery": read_tsv(core / "drug_delivery" / "delivery_assessment.tsv"),
        "docking": read_tsv(core / "drug_delivery" / "docking_assessment.tsv"),
        "open_targets": read_tsv(core / "drug_delivery" / "open_targets_summary.tsv"),
        "chembl": read_tsv(core / "drug_delivery" / "chembl_summary.tsv"),
    }
    for name in (
        "bulk_summary",
        "scrna_summary",
        "genetics_summary",
        "gtex_summary",
        "ranking_summary",
        "drug_summary",
        "docking_summary",
    ):
        paths = {
            "bulk_summary": core / "bulk" / "bulk_summary.json",
            "scrna_summary": core / "scrna" / "scrna_summary.json",
            "genetics_summary": core / "genetics" / "genetics_run_summary.json",
            "gtex_summary": core / "genetics" / "gtex_coloc_summary.json",
            "ranking_summary": core / "targets" / "ranking_summary.json",
            "drug_summary": core / "drug_delivery" / "drug_evidence_summary.json",
            "docking_summary": core / "drug_delivery" / "docking_summary.json",
        }
        files[name] = json.loads(paths[name].read_text(encoding="utf-8"))
    files["idd_config"] = json.loads(
        (repo / "config" / "idd_analysis.json").read_text(encoding="utf-8")
    )
    files["source_config"] = json.loads(
        (repo / "config" / "idd_sources.json").read_text(encoding="utf-8")
    )
    files["mr_source_config"] = json.loads(
        (repo / "config" / "mr_sources.json").read_text(encoding="utf-8")
    )
    ukb_metadata = repo / "manuscript" / "supplement" / "ukb_saige_722_metadata.json"
    files["ukb_metadata"] = (
        json.loads(ukb_metadata.read_text(encoding="utf-8")) if ukb_metadata.exists() else None
    )
    ukb_mr = repo / "manuscript" / "supplement" / "ukb_saige_722_mr.tsv"
    if not ukb_mr.exists():
        raise FileNotFoundError(
            "UKB outcome-cohort table is required before building the composite Figure 2"
        )
    files["ukb_mr"] = read_tsv(ukb_mr)
    files["wgcna_signed_summary"] = json.loads(
        (repo / "manuscript" / "supplement" / "wgcna_signed_summary.json").read_text(
            encoding="utf-8"
        )
    )
    return files


def build_source_table(data: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    source_config = data["source_config"]
    for item in source_config["bulk"]:
        rows.append(
            {
                "source": item["accession"],
                "data_type": "bulk expression",
                "role": item["role"],
                "sample_size": item["samples"],
                "cases": np.nan,
                "controls": np.nan,
                "platform_or_build": item["platform"],
                "analysis_status": "primary analysis"
                if item["accession"] == "GSE70362"
                else "availability-verified supporting source",
            }
        )
    sc = source_config["single_cell"]
    rows.append(
        {
            "source": sc["accession"],
            "data_type": "single-cell RNA sequencing",
            "role": sc["role"],
            "sample_size": data["scrna_summary"]["samples"],
            "cases": np.nan,
            "controls": np.nan,
            "platform_or_build": "human disc; healthy and diseased subseries",
            "analysis_status": sc["status"],
        }
    )
    for item in source_config["gwas"]:
        rows.append(
            {
                "source": item["outcome_id"],
                "data_type": "GWAS summary statistics",
                "role": "genetic outcome",
                "sample_size": item["cases"] + item["controls"],
                "cases": item["cases"],
                "controls": item["controls"],
                "platform_or_build": item["genome_build"],
                "analysis_status": item["status"],
            }
        )
    mr_sources = {item["source_id"]: item for item in data["mr_source_config"]["sources"]}
    finan = mr_sources["finan_2017_druggable_genome"]
    rows.extend(
        [
            {
                "source": "Finan 2017 Table S1",
                "data_type": "druggable-gene universe",
                "role": "gene universe",
                "sample_size": finan["expected_rows_including_header"] - 1,
                "cases": np.nan,
                "controls": np.nan,
                "platform_or_build": finan["genome_build"],
                "analysis_status": finan["status"],
            },
            {
                "source": "eQTLGen significant and full cis",
                "data_type": "whole-blood cis-eQTL summary statistics",
                "role": "lead-eQTL screen, exposure effects, colocalization, and SMR/HEIDI",
                "sample_size": int(
                    pd.to_numeric(data["harmonized"]["exposure_n"], errors="coerce").max()
                ),
                "cases": np.nan,
                "controls": np.nan,
                "platform_or_build": mr_sources["eqtlgen_full_cis"]["genome_build"],
                "analysis_status": "downloaded_verified; variant-specific maximum N",
            },
        ]
    )
    for source_id, tissue in (
        ("gtex_v8_blood_eqtl_catalogue", "GTEx v8 whole blood"),
        ("gtex_v8_muscle_eqtl_catalogue", "GTEx v8 skeletal muscle"),
        ("gtex_v8_fibroblast_eqtl_catalogue", "GTEx v8 cultured fibroblasts"),
    ):
        item = mr_sources[source_id]
        rows.append(
            {
                "source": f"{tissue} / eQTL Catalogue {item['dataset_id']}",
                "data_type": "proxy-tissue cis-eQTL summary statistics",
                "role": "proxy-tissue colocalization",
                "sample_size": item["sample_size"],
                "cases": np.nan,
                "controls": np.nan,
                "platform_or_build": item["genome_build"],
                "analysis_status": item["status"],
            }
        )
    ld_source = mr_sources["one_thousand_genomes_phase3_eur"]
    ld_sample_size = next(
        item["expected_rows"] for item in ld_source["files"] if item["name"].endswith(".fam")
    )
    rows.append(
        {
            "source": "1000 Genomes Phase 3 EUR",
            "data_type": "LD reference",
            "role": "clumping and SMR/HEIDI LD",
            "sample_size": ld_sample_size,
            "cases": np.nan,
            "controls": np.nan,
            "platform_or_build": ld_source["genome_build"],
            "analysis_status": ld_source["status"],
        }
    )
    ukb = data["ukb_metadata"]
    if ukb:
        rows.append(
            {
                "source": "UKB-SAIGE phecode 722",
                "data_type": "GWAS summary statistics",
                "role": "independent outcome-cohort evaluation",
                "sample_size": ukb["sample_size"],
                "cases": ukb["cases"],
                "controls": ukb["controls"],
                "platform_or_build": ukb["genome_build"],
                "analysis_status": (
                    "public exact-position API; source application 24460; "
                    "no participant-level access in this study"
                ),
            }
        )
    return pd.DataFrame(rows)


def build_selection_flow(data: dict[str, Any]) -> pd.DataFrame:
    deg = data["deg"]
    primary_af_deg = deg[deg["tissue"].eq("AF") & deg["pass_primary"].astype(bool)]
    bulk = data["bulk_candidates"]
    af_modules = data["wgcna_modules"][
        data["wgcna_modules"]["tissue"].eq("AF")
        & data["wgcna_modules"]["module"].astype(str).ne("0")
    ].copy()
    top_module = str(af_modules.sort_values("adjusted_p_value").iloc[0]["module"])
    top_module_genes = set(
        data["wgcna_genes"].loc[
            data["wgcna_genes"]["tissue"].eq("AF")
            & data["wgcna_genes"]["module"].astype(str).eq(top_module),
            "gene_symbol",
        ]
        .dropna()
        .astype(str)
    )
    mr_hits = data["mr_screen"][
        data["mr_screen"]["status"].eq("estimable")
        & (pd.to_numeric(data["mr_screen"]["fdr"], errors="coerce") < 0.05)
    ]
    mr_genes = set(mr_hits["gene_symbol"].dropna().astype(str))
    bulk_genes = set(bulk["gene_symbol"].dropna().astype(str))
    rows = [
        ("primary_threshold_af_deg", primary_af_deg["gene_symbol"].nunique(), "AF feature-level BH-FDR <0.05 and absolute log2 fold-change >=0.5"),
        ("druggable_af_deg_entering_bulk_route", len(bulk_genes), "Primary AF DEG, Finan-druggable, and present in the WGCNA input"),
        ("corrected_signed_af_non_grey_modules", len(af_modules), "Non-grey modules from explicit signed adjacency and signed TOM"),
        ("corrected_signed_af_fdr_modules", int((af_modules["non_grey_bh_fdr"] < 0.05).sum()), "Signed AF modules passing BH-FDR <0.05 across non-grey modules"),
        ("top_corrected_signed_af_module_genes", len(top_module_genes), f"Genes in exploratory signed AF module {top_module}"),
        ("bulk_route_genes_in_top_signed_af_module", len(bulk_genes & top_module_genes), "Overlap of DEG-defined bulk entries with the top exploratory signed AF module"),
        ("unique_genes_passing_lead_eqtl_mr_screen", len(mr_genes), "Unique genes with BH-FDR <0.05 in at least one outcome-wide lead-eQTL screen"),
        ("overlap_mr_screen_and_bulk_route", len(mr_genes & bulk_genes), "Overlap between MR-screen hits and DEG-defined bulk entries"),
        ("frozen_formal_candidate_set", data["formal_candidates"]["gene_symbol"].nunique(), "Union of MR-screen hits and DEG-defined bulk entries"),
    ]
    return pd.DataFrame(rows, columns=["selection_step", "gene_count", "definition"])


def build_genetic_table(data: dict[str, Any]) -> pd.DataFrame:
    genes = PAPER_GENE_ORDER
    mr = data["mr"]
    mr = mr[
        mr["gene_symbol"].isin(genes)
        & (mr["direction"] == "expression_to_disease")
        & (mr["status"] == "estimable")
    ].copy()
    screen = data["mr_screen"][
        data["mr_screen"]["gene_symbol"].isin(genes)
        & data["mr_screen"]["status"].eq("estimable")
    ].copy()
    screen = screen[
        ["gene_symbol", "outcome_id", "snp", "beta", "pvalue", "fdr"]
    ].rename(
        columns={
            "snp": "screen_lead_snp",
            "beta": "screen_beta",
            "pvalue": "screen_pvalue",
            "fdr": "screen_fdr",
        }
    )
    coloc = data["coloc"][data["coloc"]["gene_symbol"].isin(genes)].copy()
    smr = data["smr"].copy()
    smr_p = pd.to_numeric(smr["p_smr"], errors="coerce")
    smr_estimable = smr["status"].eq("estimable") & smr_p.notna()
    smr_test_counts = smr.loc[smr_estimable].groupby("outcome_id").size()
    smr["smr_bonferroni_test_count"] = smr["outcome_id"].map(smr_test_counts)
    smr["p_smr_bonferroni"] = (
        smr_p * pd.to_numeric(smr["smr_bonferroni_test_count"], errors="coerce")
    ).clip(upper=1.0)
    smr.loc[~smr_estimable, ["smr_bonferroni_test_count", "p_smr_bonferroni"]] = np.nan
    smr = smr[smr["gene_symbol"].isin(genes)].copy()
    steiger = data["steiger"][data["steiger"]["gene_symbol"].isin(genes)].copy()
    reverse = data["reverse"][
        data["reverse"]["gene_symbol"].isin(genes)
        & (data["reverse"]["direction"] == "disease_to_expression")
    ].copy()
    gtex = data["gtex_coloc"][
        data["gtex_coloc"]["gene_symbol"].isin(genes)
        & (data["gtex_coloc"]["status"] == "estimable")
    ].copy()
    gtex["pp4"] = pd.to_numeric(gtex["pp4"], errors="coerce")
    gtex = (
        gtex.sort_values("pp4", ascending=False)
        .drop_duplicates(["gene_symbol", "outcome_id"])
        .loc[:, ["gene_symbol", "outcome_id", "tissue", "pp4"]]
        .rename(columns={"tissue": "gtex_best_tissue", "pp4": "gtex_best_pp4"})
    )
    columns = [
        "gene_symbol",
        "ensembl_gene_id",
        "outcome_id",
        "n_instruments",
        "method",
        "beta",
        "se",
        "odds_ratio",
        "or_ci_lower",
        "or_ci_upper",
        "pvalue",
        "fdr",
        "q_pvalue",
    ]
    table = mr[columns].copy()
    table["formal_set_fdr"] = table["fdr"]
    table = table.merge(screen, on=["gene_symbol", "outcome_id"], how="left")
    table = table.merge(
        coloc[["gene_symbol", "outcome_id", "shared_variants", "pp4", "status"]].rename(
            columns={"status": "coloc_status"}
        ),
        on=["gene_symbol", "outcome_id"],
        how="left",
    )
    table = table.merge(
        smr[
            [
                "gene_symbol",
                "outcome_id",
                "beta_smr",
                "se_smr",
                "p_smr",
                "p_smr_fdr",
                "smr_bonferroni_test_count",
                "p_smr_bonferroni",
                "p_heidi",
                "n_snp_heidi",
                "heidi_pass_p_gt_0_01",
                "status",
            ]
        ].rename(columns={"status": "smr_status"}),
        on=["gene_symbol", "outcome_id"],
        how="left",
    )
    table = table.merge(
        steiger[
            ["gene_symbol", "outcome_id", "direction_supported", "status"]
        ].rename(columns={"status": "steiger_status"}),
        on=["gene_symbol", "outcome_id"],
        how="left",
    )
    table = table.merge(
        reverse[
            ["gene_symbol", "outcome_id", "n_instruments", "method", "pvalue", "fdr", "status"]
        ].rename(
            columns={
                "n_instruments": "regional_reverse_n_instruments",
                "method": "regional_reverse_method",
                "pvalue": "regional_reverse_pvalue",
                "fdr": "regional_reverse_fdr",
                "status": "regional_reverse_status",
            }
        ),
        on=["gene_symbol", "outcome_id"],
        how="left",
    )
    table = table.merge(gtex, on=["gene_symbol", "outcome_id"], how="left")
    table.insert(1, "paper_role", table["gene_symbol"].map(PAPER_ROLES))
    table.insert(3, "outcome_label", table["outcome_id"].map(OUTCOME_LABELS))
    table["direct_disc_endpoint"] = table["outcome_id"].eq("M13_INTERVERTEB")
    order = {gene: index for index, gene in enumerate(genes)}
    outcome_order = {key: index for index, key in enumerate(OUTCOME_LABELS)}
    table["_gene_order"] = table["gene_symbol"].map(order)
    table["_outcome_order"] = table["outcome_id"].map(outcome_order)
    return table.sort_values(["_gene_order", "_outcome_order"]).drop(
        columns=["_gene_order", "_outcome_order"]
    )


def recalculate_mr_egger(data: dict[str, Any]) -> pd.DataFrame:
    harmonized = data["harmonized"]
    harmonized = harmonized[harmonized["direction"].eq("expression_to_disease")].copy()
    grouped = {
        key: group
        for key, group in harmonized.groupby(["gene_symbol", "outcome_id"], sort=False)
    }
    rows: list[dict[str, Any]] = []
    formal = data["mr"][data["mr"]["direction"].eq("expression_to_disease")]
    for gene, outcome in formal[["gene_symbol", "outcome_id"]].drop_duplicates().itertuples(index=False):
        group = grouped.get((gene, outcome), pd.DataFrame())
        values = (
            group[["beta_exposure", "se_exposure", "beta_outcome", "se_outcome"]]
            .to_dict(orient="records")
            if not group.empty
            else []
        )
        row: dict[str, Any] = {
            "gene_symbol": gene,
            "outcome_id": outcome,
            "n_instruments": len(values),
            "method": "mr_egger_student_t",
            "inference_distribution": "student_t",
            "degrees_of_freedom": len(values) - 2 if len(values) >= 3 else np.nan,
            "status": "not_estimable",
            "reason": "fewer_than_3_instruments",
        }
        if len(values) >= 3:
            try:
                result = mr_egger(values)
                row.update(
                    {
                        "beta": result["slope_beta"],
                        "se": result["slope_se"],
                        "pvalue": result["slope_pvalue"],
                        "ci_lower": result["slope_ci_lower"],
                        "ci_upper": result["slope_ci_upper"],
                        "odds_ratio": math.exp(result["slope_beta"]),
                        "or_ci_lower": math.exp(result["slope_ci_lower"]),
                        "or_ci_upper": math.exp(result["slope_ci_upper"]),
                        "egger_intercept": result["intercept"],
                        "egger_intercept_se": result["intercept_se"],
                        "egger_intercept_pvalue": result["intercept_pvalue"],
                        "q": result["q"],
                        "q_df": result["q_df"],
                        "q_pvalue": result["q_pvalue"],
                        "status": "estimable",
                        "reason": "",
                    }
                )
            except ValueError as exc:
                row["reason"] = str(exc)
        rows.append(row)
    return pd.DataFrame(rows)


def recalculate_multiplicative_random_effects(data: dict[str, Any]) -> pd.DataFrame:
    harmonized = data["harmonized"]
    harmonized = harmonized[harmonized["direction"].eq("expression_to_disease")].copy()
    grouped = {
        key: group
        for key, group in harmonized.groupby(["gene_symbol", "outcome_id"], sort=False)
    }
    rows: list[dict[str, Any]] = []
    formal = data["mr"][data["mr"]["direction"].eq("expression_to_disease")]
    for gene, outcome in formal[["gene_symbol", "outcome_id"]].drop_duplicates().itertuples(index=False):
        group = grouped.get((gene, outcome), pd.DataFrame())
        values = (
            group[["beta_exposure", "se_exposure", "beta_outcome", "se_outcome"]]
            .to_dict(orient="records")
            if not group.empty
            else []
        )
        row: dict[str, Any] = {
            "gene_symbol": gene,
            "outcome_id": outcome,
            "n_instruments": len(values),
            "method": "multiplicative_random_effects_ivw",
            "status": "not_estimable",
            "reason": "fewer_than_2_instruments",
        }
        if len(values) >= 2:
            result = ivw_multiplicative_random_effects(values)
            row.update(
                {
                    "beta": result["beta"],
                    "se": result["se"],
                    "pvalue": result["pvalue"],
                    "ci_lower": result["ci_lower"],
                    "ci_upper": result["ci_upper"],
                    "odds_ratio": result["odds_ratio"],
                    "or_ci_lower": result["or_ci_lower"],
                    "or_ci_upper": result["or_ci_upper"],
                    "q": result["q"],
                    "q_df": result["q_df"],
                    "q_pvalue": result["q_pvalue"],
                    "dispersion": result["dispersion"],
                    "status": "estimable",
                    "reason": "",
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def build_sensitivity_table(
    data: dict[str, Any],
    genetic: pd.DataFrame,
    corrected_egger: pd.DataFrame,
    random_effects: pd.DataFrame,
) -> pd.DataFrame:
    sensitivity = data["mr_sensitivity"]
    sensitivity = sensitivity[
        sensitivity["gene_symbol"].isin(PAPER_GENE_ORDER)
        & (sensitivity["direction"] == "expression_to_disease")
        & (sensitivity["status"] == "estimable")
        & (sensitivity["method"] != "mr_egger")
    ].copy()
    pivot = sensitivity.pivot_table(
        index=["gene_symbol", "outcome_id"],
        columns="method",
        values=["beta", "pvalue", "egger_intercept", "egger_intercept_pvalue"],
        aggfunc="first",
    )
    pivot.columns = [f"{method}_{metric}" for metric, method in pivot.columns]
    pivot = pivot.reset_index()
    egger = corrected_egger[
        corrected_egger["gene_symbol"].isin(PAPER_GENE_ORDER)
        & corrected_egger["status"].eq("estimable")
    ][
        [
            "gene_symbol",
            "outcome_id",
            "beta",
            "pvalue",
            "egger_intercept",
            "egger_intercept_pvalue",
        ]
    ].rename(
        columns={
            "beta": "mr_egger_beta",
            "pvalue": "mr_egger_pvalue",
            "egger_intercept": "mr_egger_egger_intercept",
            "egger_intercept_pvalue": "mr_egger_egger_intercept_pvalue",
        }
    )
    mre = random_effects[
        random_effects["gene_symbol"].isin(PAPER_GENE_ORDER)
        & random_effects["status"].eq("estimable")
    ][
        [
            "gene_symbol",
            "outcome_id",
            "beta",
            "se",
            "pvalue",
            "or_ci_lower",
            "or_ci_upper",
            "dispersion",
        ]
    ].rename(
        columns={
            "beta": "multiplicative_random_effects_beta",
            "se": "multiplicative_random_effects_se",
            "pvalue": "multiplicative_random_effects_pvalue",
            "or_ci_lower": "multiplicative_random_effects_or_ci_lower",
            "or_ci_upper": "multiplicative_random_effects_or_ci_upper",
            "dispersion": "multiplicative_random_effects_dispersion",
        }
    )
    loo = data["loo"]
    loo = loo[
        loo["gene_symbol"].isin(PAPER_GENE_ORDER)
        & (loo["direction"] == "expression_to_disease")
        & (loo["status"] == "estimable")
    ].copy()
    full_beta = genetic.set_index(["gene_symbol", "outcome_id"])["beta"].to_dict()
    loo_rows = []
    for (gene, outcome), group in loo.groupby(["gene_symbol", "outcome_id"], sort=False):
        reference = full_beta.get((gene, outcome), np.nan)
        signs = np.sign(pd.to_numeric(group["beta"], errors="coerce"))
        loo_rows.append(
            {
                "gene_symbol": gene,
                "outcome_id": outcome,
                "loo_estimable_runs": len(group),
                "loo_sign_reversal": bool(np.isfinite(reference) and (signs != np.sign(reference)).any()),
                "loo_any_ci_crosses_null": bool(
                    ((group["ci_lower"] <= 0) & (group["ci_upper"] >= 0)).any()
                ),
                "loo_min_pvalue": group["pvalue"].min(),
                "loo_max_pvalue": group["pvalue"].max(),
            }
        )
    loo_summary = pd.DataFrame(loo_rows)
    base = genetic[
        ["gene_symbol", "outcome_id", "n_instruments", "beta", "pvalue", "q_pvalue"]
    ].rename(
        columns={
            "beta": "primary_beta",
            "pvalue": "primary_pvalue",
            "q_pvalue": "cochran_q_pvalue",
        }
    )
    return (
        base.merge(pivot, on=["gene_symbol", "outcome_id"], how="left")
        .merge(egger, on=["gene_symbol", "outcome_id"], how="left")
        .merge(mre, on=["gene_symbol", "outcome_id"], how="left")
        .merge(loo_summary, on=["gene_symbol", "outcome_id"], how="left")
    )


def aggregate_scrna(data: dict[str, Any]) -> pd.DataFrame:
    frame = data["scrna"][data["scrna"]["gene_symbol"].isin(PAPER_GENE_ORDER)].copy()
    frame["cells"] = pd.to_numeric(frame["cells"], errors="coerce")
    for column in ("cell_fraction_expressing", "mean_log1p_expression"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    rows = []
    for (gene, cell_type), group in frame.groupby(["gene_symbol", "broad_cell_type"], sort=False):
        weights = group["cells"].fillna(0).to_numpy()
        total = weights.sum()
        rows.append(
            {
                "gene_symbol": gene,
                "broad_cell_type": cell_type,
                "cells": total,
                "cell_fraction_expressing": np.average(
                    group["cell_fraction_expressing"], weights=weights
                )
                if total > 0
                else np.nan,
                "mean_log1p_expression": np.average(
                    group["mean_log1p_expression"], weights=weights
                )
                if total > 0
                else np.nan,
                "aggregation_scope": "weighted across tissue and disease-state strata; descriptive only",
            }
        )
    return pd.DataFrame(rows)


def build_translation_table(data: dict[str, Any]) -> pd.DataFrame:
    evidence = data["evidence"].set_index("gene_symbol")
    delivery = data["delivery"].set_index("gene_symbol")
    docking = data["docking"].set_index("gene_symbol")
    rows = []
    for gene in PAPER_GENE_ORDER:
        ev = evidence.loc[gene]
        de = delivery.loc[gene]
        dock = docking.loc[gene] if gene in docking.index else None
        rows.append(
            {
                "gene_symbol": gene,
                "paper_role": PAPER_ROLES[gene],
                "druggability_tier": ev["druggability_tier"],
                "desired_intervention": ev["desired_intervention"],
                "known_modulator": bool(de["known_modulator"]),
                "database_example": ev.get("top_modulator_name", np.nan),
                "database_example_id": ev.get("top_modulator_id", np.nan),
                "database_source": ev.get("top_modulator_source", np.nan),
                "direction_aligned": bool(ev.get("top_modulator_direction_aligned", False)),
                "preferred_local_payload_options": DIRECTION_MATCHED_PAYLOADS[gene],
                "direction_matched_payload_status": "plausible_formulation_hypothesis",
                "generic_formulation_feasibility": de["local_delivery_feasibility"],
                "docking_status": dock["docking_status"] if dock is not None else "not_assessed",
                "docking_reason": dock["reason"] if dock is not None else "locus_alternative_not_prioritized",
                "claim_scope": de["delivery_claim_scope"],
            }
        )
    return pd.DataFrame(rows)


def external_snapshot_paths(out: Path) -> dict[str, Path]:
    supplement = out / "supplement"
    return {
        "enrichment": supplement / "enrichment_results.tsv",
        "enrichment_queries": supplement / "enrichment_queries.tsv",
        "ppi_mapping": supplement / "ppi_mapping.tsv",
        "ppi_edges": supplement / "ppi_physical_edges.tsv",
        "ppi_summary": supplement / "ppi_summary.tsv",
        "metadata": supplement / "external_analysis_metadata.json",
    }


def run_enrichment(
    data: dict[str, Any], paths: dict[str, Path], refresh: bool
) -> tuple[pd.DataFrame, pd.DataFrame]:
    enrichment_path = paths["enrichment"]
    query_path = paths["enrichment_queries"]
    af_signed_modules = data["wgcna_modules"][
        data["wgcna_modules"]["tissue"].eq("AF")
        & data["wgcna_modules"]["module"].astype(str).ne("0")
    ]
    af_signed_fdr_module_count = int(
        (pd.to_numeric(af_signed_modules["non_grey_bh_fdr"], errors="coerce") < 0.05).sum()
    )
    if enrichment_path.exists() and query_path.exists() and not refresh:
        cached_results = read_tsv(enrichment_path)
        cached_queries = read_tsv(query_path)
        cached_query_names = set(cached_queries.get("query_name", pd.Series(dtype=str)).astype(str))
        cached_result_names = set(cached_results.get("query_name", pd.Series(dtype=str)).astype(str))
        cached_wgcna = cached_queries[cached_queries["query_name"].eq("AF_WGCNA_SIGNED")]
        expected_wgcna_status = "not_estimable" if af_signed_fdr_module_count == 0 else "not_run"
        wgcna_cache_valid = (
            len(cached_wgcna) == 1
            and float(cached_wgcna.iloc[0]["foreground_size"]) == 0
            and cached_wgcna.iloc[0]["status"] == expected_wgcna_status
        )
        if (
            cached_query_names == ENRICHMENT_QUERY_NAMES
            and cached_result_names <= {"AF_DEG"}
            and wgcna_cache_valid
        ):
            return cached_results, cached_queries
    deg = data["deg"]
    deg = deg[deg["tissue"].eq("AF")].copy()
    deg_foreground = sorted(
        deg.loc[deg["pass_primary"].astype(bool), "gene_symbol"].dropna().astype(str).unique()
    )
    deg_background = sorted(deg["gene_symbol"].dropna().astype(str).unique())
    queries = [("AF_DEG", deg_foreground, deg_background)]
    result_rows: list[dict[str, Any]] = []
    query_rows: list[dict[str, Any]] = []
    for name, foreground, background in queries:
        payload = {
            "organism": "hsapiens",
            "query": foreground,
            "background": background,
            "sources": ["GO:BP", "GO:CC", "GO:MF", "REAC"],
            "user_threshold": 0.05,
            "significance_threshold_method": "fdr",
            "domain_scope": "custom",
            "no_evidences": True,
        }
        status = "not_estimable"
        reason = "none"
        try:
            response = requests.post(
                "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
                json=payload,
                timeout=90,
            )
            response.raise_for_status()
            result = response.json().get("result", [])
            status = "estimable" if result else "no_significant_terms"
            for item in result:
                intersections = item.get("intersections", [])
                if intersections and isinstance(intersections[0], list):
                    intersection_genes = [
                        gene
                        for gene, evidence in zip(foreground, intersections)
                        if evidence
                    ]
                else:
                    intersection_genes = item.get("intersection", [])
                result_rows.append(
                    {
                        "query_name": name,
                        "source": item.get("source"),
                        "term_id": item.get("native"),
                        "term_name": item.get("name"),
                        "adjusted_pvalue": item.get("p_value"),
                        "term_size": item.get("term_size"),
                        "query_size": item.get("query_size"),
                        "intersection_size": item.get("intersection_size"),
                        "effective_domain_size": item.get("effective_domain_size"),
                        "intersection_genes": ";".join(intersection_genes),
                    }
                )
        except (requests.RequestException, ValueError) as exc:
            reason = str(exc)
        query_rows.append(
            {
                "query_name": name,
                "foreground_size": len(foreground),
                "background_size": len(background),
                "significant_term_count": sum(
                    row["query_name"] == name for row in result_rows
                ),
                "status": status,
                "reason": reason,
                "method": "g:Profiler ORA; custom measured-gene background; FDR <= 0.05",
            }
        )
    query_rows.append(
        {
            "query_name": "AF_WGCNA_SIGNED",
            "foreground_size": 0,
            "background_size": int(
                data["wgcna_genes"].loc[data["wgcna_genes"]["tissue"].eq("AF"), "gene_symbol"].nunique()
            ),
            "significant_term_count": 0,
            "status": "not_estimable" if af_signed_fdr_module_count == 0 else "not_run",
            "reason": (
                "no_non_grey_module_passed_bh_fdr_0.05_in_corrected_signed_reanalysis"
                if af_signed_fdr_module_count == 0
                else "wgcna_ora_outside_registered_manuscript_scope"
            ),
            "method": "WGCNA ORA not run; manuscript ORA is restricted to primary AF DEGs",
        }
    )
    result_columns = [
        "query_name",
        "source",
        "term_id",
        "term_name",
        "adjusted_pvalue",
        "term_size",
        "query_size",
        "intersection_size",
        "effective_domain_size",
        "intersection_genes",
    ]
    results = pd.DataFrame(result_rows, columns=result_columns)
    query_frame = pd.DataFrame(query_rows)
    write_tsv(results, enrichment_path)
    write_tsv(query_frame, query_path)
    return results, query_frame


def request_string_mapping(genes: list[str]) -> pd.DataFrame:
    response = requests.post(
        "https://string-db.org/api/tsv/get_string_ids",
        data={"identifiers": "\r".join(genes), "species": 9606, "limit": 1},
        timeout=90,
    )
    response.raise_for_status()
    lines = response.text.strip()
    if not lines:
        return pd.DataFrame()
    from io import StringIO

    return pd.read_csv(StringIO(lines), sep="\t")


def request_string_edges(genes: list[str]) -> pd.DataFrame:
    response = requests.post(
        "https://string-db.org/api/tsv/network",
        data={
            "identifiers": "\r".join(genes),
            "species": 9606,
            "required_score": 700,
            "network_type": "physical",
            "add_nodes": 0,
        },
        timeout=90,
    )
    response.raise_for_status()
    lines = response.text.strip()
    if not lines:
        return pd.DataFrame()
    from io import StringIO

    return pd.read_csv(StringIO(lines), sep="\t")


def run_ppi(
    data: dict[str, Any], paths: dict[str, Path], refresh: bool
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mapping_path = paths["ppi_mapping"]
    edge_path = paths["ppi_edges"]
    summary_path = paths["ppi_summary"]
    deg = data["deg"]
    evidence = data["evidence"]
    gene_sets = {
        "AF_DEG": sorted(
            deg.loc[
                deg["tissue"].eq("AF") & deg["pass_primary"].astype(bool), "gene_symbol"
            ]
            .dropna()
            .astype(str)
            .unique()
        ),
        "CAUSAL_COMPLETE": sorted(
            evidence.loc[evidence["causal_complete"].astype(bool), "gene_symbol"]
            .dropna()
            .astype(str)
            .unique()
        ),
    }
    if all(path.exists() for path in (mapping_path, edge_path, summary_path)) and not refresh:
        mappings = read_tsv(mapping_path)
        edges = read_tsv(edge_path)
        summary = read_tsv(summary_path)
        summary_query_names = set(summary.get("query_name", pd.Series(dtype=str)).astype(str))
        mapping_query_names = set(mappings.get("query_name", pd.Series(dtype=str)).astype(str))
        edge_query_names = set(edges.get("query_name", pd.Series(dtype=str)).astype(str))
        cached_input_counts = {
            str(row.query_name): int(float(row.input_gene_count))
            for row in summary[["query_name", "input_gene_count"]].itertuples(index=False)
        }
        expected_input_counts = {name: len(genes) for name, genes in gene_sets.items()}
        if (
            summary_query_names == PPI_QUERY_NAMES
            and mapping_query_names <= PPI_QUERY_NAMES
            and edge_query_names <= PPI_QUERY_NAMES
            and cached_input_counts == expected_input_counts
        ):
            mapping_key = "queryItem" if "queryItem" in mappings.columns else "queryIndex"
            if mapping_key in mappings.columns:
                mapped_counts = mappings.groupby("query_name")[mapping_key].nunique()
                summary["mapped_gene_count"] = (
                    summary["query_name"].map(mapped_counts).fillna(0).astype(int)
                )
                write_tsv(summary, summary_path)
            return mappings, edges, summary
    mapping_frames = []
    edge_frames = []
    summary_rows = []
    for name, genes in gene_sets.items():
        status = "not_estimable"
        reason = "none"
        mapping = pd.DataFrame()
        edges = pd.DataFrame()
        try:
            mapping = request_string_mapping(genes)
            edges = request_string_edges(genes)
            status = "estimable"
        except (requests.RequestException, ValueError) as exc:
            reason = str(exc)
        if not mapping.empty:
            mapping.insert(0, "query_name", name)
            mapping_frames.append(mapping)
        if not edges.empty:
            edges.insert(0, "query_name", name)
            edge_frames.append(edges)
        summary_rows.append(
            {
                "query_name": name,
                "input_gene_count": len(genes),
                "mapped_gene_count": (
                    mapping["queryItem"].nunique()
                    if "queryItem" in mapping.columns
                    else mapping["queryIndex"].nunique()
                    if "queryIndex" in mapping.columns
                    else 0
                ),
                "physical_edge_count": len(edges),
                "connected_node_count": len(
                    set(edges.get("preferredName_A", [])).union(
                        set(edges.get("preferredName_B", []))
                    )
                ),
                "status": status,
                "reason": reason,
                "method": "STRING physical network; score >= 0.700; no added nodes",
            }
        )
    mappings = (
        pd.concat(mapping_frames, ignore_index=True)
        if mapping_frames
        else pd.DataFrame(columns=["query_name"])
    )
    edges = (
        pd.concat(edge_frames, ignore_index=True)
        if edge_frames
        else pd.DataFrame(columns=["query_name"])
    )
    summary = pd.DataFrame(summary_rows)
    write_tsv(mappings, mapping_path)
    write_tsv(edges, edge_path)
    write_tsv(summary, summary_path)
    return mappings, edges, summary


def make_study_design_figure(data: dict[str, Any], output: Path) -> list[Path]:
    bulk = data["bulk_summary"]
    genetics = data["genetics_summary"]
    sc = data["scrna_summary"]
    drug = data["drug_summary"]
    docking = data["docking_summary"]
    mr_estimable = int((data["mr"]["status"] == "estimable").sum())
    coloc_pass = int(
        ((data["coloc"]["status"] == "estimable") & (data["coloc"]["pp4"] >= 0.8)).sum()
    )
    fig, ax = plt.subplots(figsize=(FULL_WIDTH_IN, 7.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    columns = [
        (
            "1. Tissue\ndiscovery",
            "GSE70362; 48 profiles / 18 donors\n"
            f"AF: {bulk['deg_counts'][0]['fdr_effect']} primary DEGs\n"
            "Signed WGCNA\n"
            f"{data['wgcna_signed_summary']['af_fdr_significant_module_count']} FDR-significant modules",
            "#005A8D",
        ),
        (
            "2. Genetic\ntriangulation",
            f"{genetics['candidate_count']} formal candidates; 3 FinnGen outcomes\n"
            f"{mr_estimable} MR estimates; {coloc_pass} coloc PP4 \N{GREATER-THAN OR EQUAL TO} 0.80\n"
            "SMR / HEIDI / Steiger\n"
            "Regional direction analysis",
            "#805400",
        ),
        (
            "3. Cell\nlocalization",
            f"GSE230809; {sc['samples']} samples; {sc['cells_after_qc']:,} QC cells\n"
            f"{sc['clustered_cells']:,} graph-subset cells; {sc['clusters']} graph clusters\n"
            "5 marker-score broad strata\n"
            "Descriptive localization only",
            "#006B4F",
        ),
        (
            "4. Translation\nscreen",
            f"{drug['interaction_rows']:,} drug-target interactions\n"
            "Local payload feasibility\n"
            f"{docking['eligible_pair_count']} docking-eligible pairs",
            "#983700",
        ),
    ]
    y_positions = [0.70, 0.51, 0.32, 0.13]
    box_height = 0.145
    for index, ((title, body, color), y) in enumerate(zip(columns, y_positions)):
        box = FancyBboxPatch(
            (0.06, y),
            0.88,
            box_height,
            boxstyle="round,pad=0.008,rounding_size=0.012",
            linewidth=1.5,
            edgecolor=color,
            facecolor="white",
        )
        ax.add_patch(box)
        ax.add_patch(
            FancyBboxPatch(
                (0.06, y),
                0.26,
                box_height,
                boxstyle="round,pad=0.008,rounding_size=0.012",
                linewidth=0,
                facecolor=color,
            )
        )
        ax.text(
            0.19,
            y + box_height / 2,
            title,
            ha="center",
            va="center",
            color="white",
            weight="bold",
            fontsize=8.5,
        )
        ax.text(
            0.63,
            y + box_height / 2,
            body,
            ha="center",
            va="center",
            color=PALETTE["ink"],
            fontsize=8.5,
            linespacing=1.15,
        )
        if index < len(columns) - 1:
            ax.annotate(
                "",
                xy=(0.5, y_positions[index + 1] + box_height + 0.008),
                xytext=(0.5, y - 0.008),
                arrowprops={"arrowstyle": "-|>", "lw": 1.4, "color": PALETTE["gray"]},
            )
    ax.text(
        0.5,
        0.94,
        "Multi-layer prioritization of druggable targets for intervertebral disc disorders",
        ha="center",
        va="center",
        fontsize=12,
        weight="bold",
        color=PALETTE["ink"],
    )
    ax.text(
        0.5,
        0.055,
        "Interpretation is endpoint-specific: primary disc-disorder evidence is separated\n"
        "from low-back-pain and sciatica evidence.",
        ha="center",
        va="center",
        fontsize=8.5,
        color=PALETTE["gray"],
    )
    return save_figure(fig, output)


def plot_mr_axis(ax: plt.Axes, genetic: pd.DataFrame, *, show_legend: bool = True) -> None:
    frame = genetic.copy()
    frame["label"] = frame["gene_symbol"] + "  |  " + frame["outcome_label"]
    frame = frame.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(frame))
    for index, row in frame.iterrows():
        color = OUTCOME_COLORS[row["outcome_id"]]
        significant = bool(row["screen_fdr"] < 0.05)
        ax.plot([row["or_ci_lower"], row["or_ci_upper"]], [index, index], color=color, lw=1.6)
        marker = "D" if row["method"] == "wald_ratio" else "o"
        ax.scatter(
            row["odds_ratio"],
            index,
            facecolor=color if significant else "white",
            edgecolor=color,
            marker=marker,
            s=42,
            zorder=3,
            linewidth=1.2,
        )
        ax.text(
            1.43,
            index,
            f"{row['odds_ratio']:.3f} ({row['or_ci_lower']:.3f}\N{EN DASH}{row['or_ci_upper']:.3f})",
            va="center",
            fontsize=8.5,
        )
    ax.axvline(1.0, color=PALETTE["gray"], ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xlim(0.80, 1.60)
    ax.set_xticks([0.8, 0.9, 1.0, 1.1, 1.25, 1.5])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks(y)
    ax.set_yticklabels(frame["label"])
    ax.set_xlabel("Odds ratio per genetically predicted unit increase in expression")
    ax.set_title("Forward cis-eQTL Mendelian randomization across prioritized genes")
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    handles = [
        plt.Line2D([0], [0], color=OUTCOME_COLORS[key], lw=2, label=label)
        for key, label in OUTCOME_LABELS.items()
    ]
    handles.extend(
        [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="white",
                markeredgecolor=PALETTE["gray"],
                label="IVW",
                markersize=6,
            ),
            plt.Line2D(
                [0],
                [0],
                marker="D",
                color="none",
                markerfacecolor="white",
                markeredgecolor=PALETTE["gray"],
                label="Wald ratio",
                markersize=6,
            ),
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color=PALETTE["ink"],
                markerfacecolor=PALETTE["ink"],
                label="Screen FDR < 0.05",
                markersize=6,
            ),
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color=PALETTE["ink"],
                markerfacecolor="white",
                label="Screen FDR \N{GREATER-THAN OR EQUAL TO} 0.05",
                markersize=6,
            ),
        ]
    )
    if show_legend:
        ax.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.08),
            frameon=False,
            ncol=4,
        )


def make_mr_forest(genetic: pd.DataFrame, output: Path) -> list[Path]:
    fig, ax = plt.subplots(figsize=(FULL_WIDTH_IN, 6.8))
    plot_mr_axis(ax, genetic)
    return save_figure(fig, output)


def make_figure_2_composite(
    genetic: pd.DataFrame,
    ukb_mr: pd.DataFrame,
    output: Path,
) -> list[Path]:
    fig = plt.figure(figsize=(FULL_WIDTH_IN, 9.4))
    grid = fig.add_gridspec(2, 1, height_ratios=[2.45, 1.0], hspace=0.46)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[1, 0])
    plot_mr_axis(ax_a, genetic)
    plot_replication_axis(ax_b, ukb_mr, show_legend=False)
    add_panel_label(ax_a, "A", x=-0.16, y=1.04)
    add_panel_label(ax_b, "B", x=-0.16, y=1.08)
    return save_figure(fig, output)


def heatmap_matrix(frame: pd.DataFrame, value: str) -> pd.DataFrame:
    matrix = frame.pivot(index="gene_symbol", columns="outcome_id", values=value)
    return matrix.reindex(index=PAPER_GENE_ORDER, columns=list(OUTCOME_LABELS))


def make_genetic_heatmap(genetic: pd.DataFrame, output: Path) -> list[Path]:
    frame = genetic.copy()
    frame["mr_strength"] = -np.log10(
        pd.to_numeric(frame["screen_fdr"], errors="coerce").clip(lower=1e-12)
    )
    frame["smr_strength"] = -np.log10(
        pd.to_numeric(frame["p_smr_bonferroni"], errors="coerce").clip(lower=1e-12)
    )
    strength_threshold = -math.log10(0.05)
    matrices = [
        (
            heatmap_matrix(frame, "mr_strength"),
            "Lead-eQTL screen: -log10(FDR)\nOutlined: FDR < 0.05",
            "YlOrRd",
            0,
            5,
            lambda value: value > strength_threshold,
        ),
        (
            heatmap_matrix(frame, "pp4"),
            "eQTLGen coloc: PP4\nOutlined: PP4 \N{GREATER-THAN OR EQUAL TO} 0.80",
            "Blues",
            0,
            1,
            lambda value: value >= 0.80,
        ),
        (
            heatmap_matrix(frame, "smr_strength"),
            "SMR: -log10(Bonferroni P)\nOutlined: adjusted P < 0.05",
            "YlGn",
            0,
            5,
            lambda value: value > strength_threshold,
        ),
        (
            heatmap_matrix(frame, "p_heidi"),
            "HEIDI P value\nOutlined: P > 0.01",
            "Purples",
            0,
            1,
            lambda value: value > 0.01,
        ),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(FULL_WIDTH_IN, 7.4))
    axes = axes.ravel()
    labels = [OUTCOME_LABELS[key].replace("Intervertebral disc disorders", "Disc disorders") for key in OUTCOME_LABELS]
    for panel, (ax, (matrix, title, cmap, vmin, vmax, passes_threshold)) in enumerate(
        zip(axes, matrices)
    ):
        colorbar_breaks = np.linspace(vmin, vmax, 6)
        annotation = matrix.map(
            lambda value: "NE"
            if pd.isna(value)
            else "0.00"
            if abs(value) < 0.005
            else f"{value:.2f}"
        )
        sns.heatmap(
            matrix,
            ax=ax,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            annot=annotation,
            fmt="",
            linewidths=0.8,
            linecolor="white",
            cbar_kws={
                "shrink": 0.65,
                "boundaries": colorbar_breaks,
                "ticks": colorbar_breaks,
                "drawedges": True,
            },
            annot_kws={"fontsize": 8.5},
        )
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_yticklabels(PAPER_GENE_ORDER, rotation=0)
        for row_index in range(matrix.shape[0]):
            for column_index in range(matrix.shape[1]):
                value = matrix.iloc[row_index, column_index]
                if pd.notna(value) and passes_threshold(float(value)):
                    ax.add_patch(
                        Rectangle(
                            (column_index, row_index),
                            1,
                            1,
                            fill=False,
                            edgecolor=PALETTE["ink"],
                            linewidth=1.5,
                        )
                    )
        add_panel_label(ax, chr(ord("A") + panel), x=-0.30, y=1.14)
    fig.suptitle("Endpoint-specific genetic evidence", fontsize=12, weight="bold", y=0.995)
    fig.text(
        0.5,
        0.012,
        "NE, not estimable. HEIDI outlines indicate no detected heterogeneity, not causal proof.",
        ha="center",
        va="bottom",
        fontsize=8,
        color=PALETTE["gray"],
    )
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.19, top=0.91, hspace=0.68, wspace=0.38)
    return save_figure(fig, output)


def make_scrna_heatmap(scrna: pd.DataFrame, output: Path) -> list[Path]:
    cell_order = ["nucleus_pulposus", "chondrocyte_like", "annulus_fibrosus", "immune", "endothelial"]
    fraction = scrna.pivot(index="gene_symbol", columns="broad_cell_type", values="cell_fraction_expressing").reindex(
        index=PAPER_GENE_ORDER, columns=cell_order
    )
    expression = scrna.pivot(index="gene_symbol", columns="broad_cell_type", values="mean_log1p_expression").reindex(
        index=PAPER_GENE_ORDER, columns=cell_order
    )
    fig, axes = plt.subplots(2, 1, figsize=(FULL_WIDTH_IN, 7.4), sharex=True)
    fraction_max = max(0.01, float(np.nanmax(fraction.to_numpy())))
    expression_max = max(0.01, float(np.nanmax(expression.to_numpy())))
    fraction_breaks = np.linspace(0, fraction_max, 6)
    expression_breaks = np.linspace(0, expression_max, 6)
    sns.heatmap(
        fraction,
        ax=axes[0],
        cmap="Blues",
        vmin=0,
        vmax=fraction_max,
        annot=fraction.map(format_compact_value),
        fmt="",
        linewidths=0.8,
        linecolor="white",
        cbar_kws={
            "label": "Fraction expressing",
            "boundaries": fraction_breaks,
            "ticks": fraction_breaks,
            "drawedges": True,
            "format": "%.2f",
        },
        annot_kws={"fontsize": 8.5},
    )
    sns.heatmap(
        expression,
        ax=axes[1],
        cmap="Oranges",
        vmin=0,
        vmax=expression_max,
        annot=expression.map(format_compact_value),
        fmt="",
        linewidths=0.8,
        linecolor="white",
        cbar_kws={
            "label": "Mean log1p expression",
            "boundaries": expression_breaks,
            "ticks": expression_breaks,
            "drawedges": True,
            "format": "%.2f",
        },
        annot_kws={"fontsize": 8.5},
    )
    for ax, matrix, maximum in zip(
        axes,
        (fraction, expression),
        (fraction_max, expression_max),
    ):
        for text, value in zip(ax.texts, matrix.to_numpy().ravel()):
            if pd.notna(value) and value >= maximum * 0.55:
                text.set_color("white")
    pretty = [
        "Nucleus pulposus-like",
        "Chondrocyte-like",
        "Annulus fibrosus-like",
        "Immune-like",
        "Endothelial-like",
    ]
    for index, (ax, title) in enumerate(
        zip(axes, ["Cell fraction expressing", "Mean normalized expression"])
    ):
        ax.set_title(title)
        ax.set_ylabel("")
        ax.set_yticklabels(PAPER_GENE_ORDER, rotation=0)
        if index == 0:
            ax.set_xlabel("")
            ax.tick_params(axis="x", labelbottom=False)
        else:
            ax.set_xlabel("Marker-score broad cell stratum")
            ax.set_xticklabels(pretty, rotation=35, ha="right")
            ax.tick_params(axis="x", labelbottom=True)
        add_panel_label(ax, chr(ord("A") + index), x=-0.13, y=1.05)
    fig.suptitle("Descriptive single-cell localization in GSE230809", fontsize=12, weight="bold", y=0.995)
    fig.subplots_adjust(left=0.13, right=0.92, bottom=0.12, top=0.93, hspace=0.32)
    return save_figure(fig, output)


def make_translation_figure(translation: pd.DataFrame, genetic: pd.DataFrame, output: Path) -> list[Path]:
    direct = genetic[genetic["outcome_id"].eq("M13_INTERVERTEB")].set_index("gene_symbol")
    symptom = genetic[genetic["outcome_id"].isin(["M13_LOWBACKPAIN", "M13_SCIATICA"])].copy()
    symptom_support = symptom.groupby("gene_symbol").apply(
        lambda group: bool(((group["screen_fdr"] < 0.05) & (group["pp4"] >= 0.8)).any()),
        include_groups=False,
    )
    matrix = pd.DataFrame(index=PAPER_GENE_ORDER)
    matrix["Primary disc\nMR + coloc"] = [
        "Supported"
        if bool(
            gene in direct.index
            and direct.loc[gene, "screen_fdr"] < 0.05
            and direct.loc[gene, "pp4"] >= 0.8
        )
        else "Not supported"
        for gene in PAPER_GENE_ORDER
    ]
    matrix["Symptom\nMR + coloc"] = [
        "Supported" if bool(symptom_support.get(gene, False)) else "Not supported"
        for gene in PAPER_GENE_ORDER
    ]
    trans = translation.set_index("gene_symbol")
    matrix["Database\nmodulator"] = [
        "Database hit" if bool(trans.loc[gene, "known_modulator"]) else "Not found"
        for gene in PAPER_GENE_ORDER
    ]
    matrix["Direction\nalignment"] = [
        "N/A"
        if not bool(trans.loc[gene, "known_modulator"])
        else "Aligned"
        if bool(trans.loc[gene, "direction_aligned"])
        else "Not aligned"
        for gene in PAPER_GENE_ORDER
    ]
    matrix["Local payload"] = [
        "Formulation hypothesis"
        if trans.loc[gene, "direction_matched_payload_status"]
        == "plausible_formulation_hypothesis"
        else "Not assessed"
        for gene in PAPER_GENE_ORDER
    ]
    docking_labels = {
        "COMPLETED": "Completed",
        "NOT_APPLICABLE": "N/A",
        "not_assessed": "Not assessed",
    }
    matrix["Docking"] = [
        docking_labels.get(str(trans.loc[gene, "docking_status"]), "Not assessed")
        for gene in PAPER_GENE_ORDER
    ]

    cell_styles = {
        "Supported": ("#DCEAF7", PALETTE["navy"], ""),
        "Not supported": ("#EFF1F3", PALETTE["gray"], ""),
        "Database hit": ("#DCEAF7", PALETTE["navy"], ""),
        "Not found": ("white", PALETTE["gray"], "///"),
        "Aligned": ("#DDEFE8", PALETTE["green"], ""),
        "Not aligned": ("#F9E7CF", "#8A5A00", ""),
        "Formulation hypothesis": ("#EFE6F2", PALETTE["purple"], ""),
        "Completed": ("#DDEFE8", PALETTE["green"], ""),
        "N/A": ("white", PALETTE["gray"], "///"),
        "Not assessed": ("#EFF1F3", PALETTE["gray"], ""),
    }
    display_labels = {
        "Not supported": "Not\nsupported",
        "Database hit": "Database\nhit",
        "Not aligned": "Not\naligned",
        "Formulation hypothesis": "Formulation\nhypothesis",
        "Not assessed": "Not\nassessed",
    }
    row_labels = {
        "CCL28": "CCL28\n(1 instrument)",
        "AKR1C1": "AKR1C1\n(locus alternative)",
    }
    group_headers = [
        (0, 2, "Genetic evidence", "#005A8D"),
        (2, 2, "Modulator evidence", "#805400"),
        (4, 2, "Delivery / structure", "#6F4A78"),
    ]

    fig, ax = plt.subplots(figsize=(FULL_WIDTH_IN, 5.2))
    row_count, column_count = matrix.shape
    ax.set_xlim(-0.02, column_count)
    ax.set_ylim(row_count + 0.62, -1.35)
    for start, width, label, color in group_headers:
        ax.add_patch(
            Rectangle(
                (start + 0.02, -1.28),
                width - 0.04,
                0.38,
                facecolor=color,
                edgecolor="none",
            )
        )
        ax.text(
            start + width / 2,
            -1.09,
            label,
            ha="center",
            va="center",
            color="white",
            fontsize=8.5,
            fontweight="bold",
        )
    for column_index, column in enumerate(matrix.columns):
        ax.text(
            column_index + 0.5,
            -0.43,
            column,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
            color=PALETTE["ink"],
        )
    for row_index, gene in enumerate(matrix.index):
        ax.text(
            -0.10,
            row_index + 0.5,
            row_labels.get(gene, gene),
            ha="right",
            va="center",
            fontsize=9,
            fontstyle="italic",
            color=PALETTE["ink"],
            clip_on=False,
        )
        for column_index, value in enumerate(matrix.loc[gene]):
            facecolor, edgecolor, hatch = cell_styles[value]
            ax.add_patch(
                Rectangle(
                    (column_index + 0.03, row_index + 0.03),
                    0.94,
                    0.94,
                    facecolor=facecolor,
                    edgecolor=edgecolor,
                    linewidth=0.9,
                    hatch=hatch,
                )
            )
            ax.text(
                column_index + 0.5,
                row_index + 0.5,
                display_labels.get(value, value),
                ha="center",
                va="center",
                fontsize=7.7,
                color=PALETTE["ink"],
            )
    ax.text(
        0,
        row_count + 0.34,
        "Hatched cells denote not found or not applicable.",
        ha="left",
        va="center",
        fontsize=7.8,
        color=PALETTE["gray"],
    )
    ax.axis("off")
    fig.suptitle(
        "Genetic evidence and translational feasibility are distinct dimensions",
        fontsize=12,
        fontweight="bold",
        y=0.99,
    )
    fig.subplots_adjust(left=0.22, right=0.99, bottom=0.08, top=0.90)
    return save_figure(fig, output)


def make_bulk_volcano(data: dict[str, Any], output: Path) -> list[Path]:
    frame = data["deg"][data["deg"]["tissue"].eq("AF")].copy()
    frame["minus_log10_fdr"] = -np.log10(pd.to_numeric(frame["adj.P.Val"], errors="coerce").clip(lower=1e-300))
    frame["category"] = np.where(frame["pass_primary"].astype(bool), "FDR and effect threshold", "Other")
    fig, ax = plt.subplots(figsize=(FULL_WIDTH_IN, 5.2))
    for category, color, size, alpha in [
        ("Other", "#B9C0C7", 8, 0.45),
        ("FDR and effect threshold", PALETTE["red"], 17, 0.85),
    ]:
        group = frame[frame["category"].eq(category)]
        ax.scatter(group["logFC"], group["minus_log10_fdr"], s=size, alpha=alpha, color=color, label=category, linewidth=0)
    ax.axvline(-0.5, color=PALETTE["gray"], ls="--", lw=0.8)
    ax.axvline(0.5, color=PALETTE["gray"], ls="--", lw=0.8)
    ax.axhline(-math.log10(0.05), color=PALETTE["gray"], ls="--", lw=0.8)
    labels = (
        frame[frame["pass_primary"].astype(bool)]
        .nsmallest(6, "adj.P.Val")
        .sort_values("logFC")
    )
    label_x_positions = np.linspace(-2.55, 2.55, len(labels))
    label_y = float(frame["minus_log10_fdr"].max()) * 1.10
    for (_, row), label_x in zip(labels.iterrows(), label_x_positions):
        ax.annotate(
            str(row["gene_symbol"]),
            xy=(row["logFC"], row["minus_log10_fdr"]),
            xytext=(label_x, label_y),
            textcoords="data",
            fontsize=8.5,
            ha="center",
            va="bottom",
            arrowprops={"arrowstyle": "-", "color": PALETTE["gray"], "lw": 0.5},
        )
    ax.set_ylim(-0.08, float(frame["minus_log10_fdr"].max()) * 1.30)
    ax.set_xlabel("log2 fold change (high-grade vs low-grade AF)")
    ax.set_ylabel("-log10 BH-FDR")
    fdr_count = int((pd.to_numeric(frame["adj.P.Val"], errors="coerce") < 0.05).sum())
    primary_count = int(frame["pass_primary"].astype(bool).sum())
    ax.set_title("GSE70362 annulus fibrosus differential expression", pad=30)
    ax.text(
        0.5,
        1.025,
        f"BH-FDR < 0.05: {fdr_count} features; FDR + |log2FC| \N{GREATER-THAN OR EQUAL TO} 0.5: {primary_count} features",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8.5,
        color=PALETTE["gray"],
    )
    ax.legend(frameon=False, loc="lower right")
    return save_figure(fig, output)


def make_wgcna_figure(wgcna_modules: pd.DataFrame, output: Path) -> list[Path]:
    frame = wgcna_modules[wgcna_modules["module"].astype(str).ne("0")].copy()
    frame["label"] = "M" + frame["module"].astype(str)
    af_top_module = str(
        frame[frame["tissue"].eq("AF")].sort_values("adjusted_p_value").iloc[0]["module"]
    )
    fig, axes = plt.subplots(1, 2, figsize=(FULL_WIDTH_IN, 8.0), sharex=True)
    for ax, tissue in zip(axes, ("AF", "NP")):
        group = (
            frame[frame["tissue"].eq(tissue)]
            .sort_values("adjusted_log2fc")
            .reset_index(drop=True)
        )
        highlighted = (
            group["module"].astype(str).eq(af_top_module)
            if tissue == "AF"
            else pd.Series(False, index=group.index)
        )
        colors = np.where(highlighted, PALETTE["orange"], "#AAB2B9")
        y = np.arange(len(group))
        ax.barh(y, group["adjusted_log2fc"], color=colors, height=0.72)
        ax.axvline(0, color=PALETTE["ink"], lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(group["label"], fontsize=8.5)
        ax.grid(axis="x", color="#D9D9D9", linewidth=0.6)
        ax.grid(axis="y", visible=False)
        if tissue == "AF":
            for tick, module in zip(ax.get_yticklabels(), group["module"].astype(str)):
                if module == af_top_module:
                    tick.set_color(PALETTE["orange"])
                    tick.set_weight("bold")
            m23 = group[group["module"].astype(str).eq(af_top_module)].iloc[0]
            ax.set_title(
                f"AF ({len(group)} non-gray modules)\n"
                f"M23: FDR = {float(m23['non_grey_bh_fdr']):.4f}; "
                f"P = {float(m23['adjusted_p_value']):.5f}; {int(m23['size'])} genes",
                fontsize=9,
                pad=6,
            )
        else:
            ax.set_title(f"NP ({len(group)} non-gray modules)", fontsize=9, pad=6)
    fig.suptitle("Signed WGCNA module effects", fontsize=12, weight="bold", y=0.99)
    fig.supxlabel(
        "Adjusted module eigengene difference (high vs low grade)",
        fontsize=9,
        y=0.02,
    )
    fig.subplots_adjust(left=0.11, right=0.99, bottom=0.07, top=0.90, wspace=0.28)
    return save_figure(fig, output)


def make_loo_figure(data: dict[str, Any], output: Path) -> list[Path]:
    combinations = [
        ("AKR1C2", "M13_INTERVERTEB"),
        ("NT5E", "M13_LOWBACKPAIN"),
        ("P2RY13", "M13_LOWBACKPAIN"),
    ]
    loo = data["loo"]
    mr = data["mr"].set_index(["gene_symbol", "outcome_id", "direction"])
    fig, axes = plt.subplots(3, 1, figsize=(FULL_WIDTH_IN, 8.4), sharex=True)
    for ax, (gene, outcome) in zip(axes, combinations):
        group = loo[
            loo["gene_symbol"].eq(gene)
            & loo["outcome_id"].eq(outcome)
            & loo["direction"].eq("expression_to_disease")
            & loo["status"].eq("estimable")
        ].copy()
        group = group.sort_values("excluded_snp").reset_index(drop=True)
        full = mr.loc[(gene, outcome, "expression_to_disease")]
        rows = pd.concat(
            [
                pd.DataFrame(
                    {
                        "excluded_snp": ["All instruments"],
                        "odds_ratio": [full["odds_ratio"]],
                        "or_ci_lower": [full["or_ci_lower"]],
                        "or_ci_upper": [full["or_ci_upper"]],
                    }
                ),
                group[["excluded_snp", "odds_ratio", "or_ci_lower", "or_ci_upper"]],
            ],
            ignore_index=True,
        ).iloc[::-1]
        y = np.arange(len(rows))
        for index, row in enumerate(rows.itertuples(index=False)):
            color = OUTCOME_COLORS[outcome]
            significant = bool(row.or_ci_lower > 1 or row.or_ci_upper < 1)
            marker = "D" if row.excluded_snp == "All instruments" else "o"
            ax.plot([row.or_ci_lower, row.or_ci_upper], [index, index], color=color, lw=1.5)
            ax.scatter(
                row.odds_ratio,
                index,
                facecolor=color if significant else "white",
                edgecolor=color,
                linewidth=1.1,
                marker=marker,
                s=38,
                zorder=3,
            )
        ax.axvline(1.0, color=PALETTE["gray"], ls="--", lw=0.9)
        ax.set_xscale("log")
        ax.set_xlim(0.70, 1.20)
        ax.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
        plot_labels = rows["excluded_snp"].map(
            lambda value: value if value == "All instruments" else f"Excluding {value}"
        )
        ax.set_yticks(y)
        ax.set_yticklabels(plot_labels)
        ax.set_title(f"{gene}: {OUTCOME_LABELS[outcome]}")
        ax.grid(axis="x", color="#D9D9D9", linewidth=0.6)
        ax.grid(axis="y", visible=False)
    fig.suptitle("Leave-one-out sensitivity analyses", fontsize=12, weight="bold", y=0.98)
    fig.supxlabel("Odds ratio", fontsize=9, y=0.025)
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="D",
            color=PALETTE["gray"],
            markerfacecolor="white",
            label="All instruments",
            markersize=6,
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color=PALETTE["gray"],
            markerfacecolor="white",
            label="One instrument excluded",
            markersize=6,
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color=PALETTE["ink"],
            markerfacecolor=PALETTE["ink"],
            label="95% CI excludes 1",
            markersize=6,
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color=PALETTE["ink"],
            markerfacecolor="white",
            label="95% CI includes 1",
            markersize=6,
        ),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.055), ncol=4, frameon=False)
    fig.subplots_adjust(left=0.21, right=0.98, bottom=0.13, top=0.91, hspace=0.58)
    return save_figure(fig, output)


def make_enrichment_figure(enrichment: pd.DataFrame, output: Path) -> list[Path]:
    if enrichment.empty:
        fig, ax = plt.subplots(figsize=(FULL_WIDTH_IN, 2.8))
        ax.axis("off")
        ax.text(0.5, 0.5, "No FDR-significant enrichment terms", ha="center", va="center")
        return save_figure(fig, output)
    frame = enrichment.copy()
    frame["strength"] = -np.log10(pd.to_numeric(frame["adjusted_pvalue"], errors="coerce").clip(lower=1e-300))
    frame = (
        frame.sort_values(["query_name", "adjusted_pvalue"])
        .groupby("query_name", group_keys=False)
        .head(10)
        .sort_values("strength")
    )
    frame["label"] = frame["term_name"].map(
        lambda value: "\n".join(textwrap.wrap(str(value), width=42))
    )
    fig, ax = plt.subplots(figsize=(FULL_WIDTH_IN, max(2.8, len(frame) * 0.55 + 2.1)))
    y = np.arange(len(frame))
    for index, row in frame.reset_index(drop=True).iterrows():
        ax.hlines(index, 0, row["strength"], color=PALETTE["orange"], linewidth=2)
        ax.scatter(
            row["strength"],
            index,
            s=35 + 4 * float(row["intersection_size"]),
            facecolor=PALETTE["orange"],
            edgecolor="#805400",
            linewidth=0.8,
            zorder=3,
        )
        ax.text(
            float(row["strength"]) + 0.035,
            index,
            f"FDR = {float(row['adjusted_pvalue']):.3f}; "
            f"{int(row['intersection_size'])}/{int(row['query_size'])} mapped genes",
            ha="left",
            va="center",
            fontsize=8.5,
            color=PALETTE["ink"],
        )
    fdr_threshold = -math.log10(0.05)
    ax.axvline(fdr_threshold, color=PALETTE["gray"], ls="--", lw=1)
    ax.text(
        fdr_threshold,
        0.92,
        "FDR = 0.05",
        transform=ax.get_xaxis_transform(),
        ha="right",
        va="top",
        fontsize=8,
        color=PALETTE["gray"],
    )
    x_max = max(float(frame["strength"].max()) * 1.55, fdr_threshold * 1.30)
    ax.set_xlim(0, x_max)
    ax.set_yticks(y)
    ax.set_yticklabels(frame.reset_index(drop=True)["label"])
    ax.set_ylim(-0.7, len(frame) - 0.3)
    ax.set_xlabel("-log10 BH-FDR")
    ax.set_ylabel("")
    ax.set_title("Measured-background functional enrichment")
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    return save_figure(fig, output)


def make_ppi_figure(edges: pd.DataFrame, output: Path) -> list[Path]:
    frame = edges[edges.get("query_name", pd.Series(dtype=str)).eq("CAUSAL_COMPLETE")].copy()
    fig, ax = plt.subplots(figsize=(FULL_WIDTH_IN, 5.0))
    graph = nx.Graph()
    queried_genes = data_gene_order()
    graph.add_nodes_from(queried_genes)
    for row in frame.itertuples(index=False):
        graph.add_edge(row.preferredName_A, row.preferredName_B, weight=float(row.score))

    positions = {
        "AKR1C1": (0.34, 0.76),
        "AKR1C2": (0.66, 0.76),
        "BMP6": (0.18, 0.48),
        "CCL28": (0.39, 0.48),
        "DLK1": (0.61, 0.48),
        "GPX1": (0.82, 0.48),
        "NR6A1": (0.18, 0.22),
        "NT5E": (0.39, 0.22),
        "P2RY13": (0.61, 0.22),
        "TSSK6": (0.82, 0.22),
        "UTP11L": (0.50, 0.05),
    }
    connected_nodes = sorted({node for edge in graph.edges for node in edge})
    isolated_nodes = [gene for gene in queried_genes if gene not in connected_nodes]
    nx.draw_networkx_nodes(
        graph,
        positions,
        nodelist=isolated_nodes,
        node_color="#E5E7E9",
        edgecolors=PALETTE["gray"],
        linewidths=1,
        node_size=1050,
        ax=ax,
    )
    if connected_nodes:
        nx.draw_networkx_nodes(
            graph,
            positions,
            nodelist=connected_nodes,
            node_color="#F9E7CF",
            edgecolors="#805400",
            linewidths=1.4,
            node_size=1250,
            ax=ax,
        )
        nx.draw_networkx_edges(
            graph,
            positions,
            edge_color=PALETTE["navy"],
            width=2.2,
            ax=ax,
        )
        if graph.number_of_edges() == 1:
            edge = next(iter(graph.edges))
            midpoint_x = (positions[edge[0]][0] + positions[edge[1]][0]) / 2
            midpoint_y = (positions[edge[0]][1] + positions[edge[1]][1]) / 2
            ax.text(
                midpoint_x,
                midpoint_y + 0.085,
                f"Physical STRING score = {float(graph.edges[edge]['weight']):.3f}",
                ha="center",
                va="center",
                fontsize=8.5,
                color=PALETTE["navy"],
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
            )
    nx.draw_networkx_labels(
        graph,
        positions,
        font_size=8.5,
        font_weight="bold",
        font_color=PALETTE["ink"],
        ax=ax,
    )
    ax.text(
        0.5,
        -0.09,
        "Gray nodes had no physical STRING edge at score \N{GREATER-THAN OR EQUAL TO} 0.700 within the queried set.",
        ha="center",
        va="center",
        fontsize=8,
        color=PALETTE["gray"],
    )
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.14, 0.96)
    ax.set_title("Integrated genetic set: high-confidence physical STRING network")
    return save_figure(fig, output)


def data_gene_order() -> list[str]:
    return ["AKR1C1", "AKR1C2", "BMP6", "CCL28", "DLK1", "GPX1", "NR6A1", "NT5E", "P2RY13", "TSSK6", "UTP11L"]


def prepare_submission_tiffs(repo: Path, figures: Path) -> pd.DataFrame:
    source_by_target = {
        "Figure_1.tiff": figures / "figure_1_study_design.tiff",
        "Figure_2.tiff": figures / "figure_2_composite.tiff",
        "Figure_3.tiff": figures / "figure_3_genetic_evidence.tiff",
        "Figure_4.tiff": figures / "figure_4_single_cell_localization.tiff",
        "Figure_5.tiff": figures / "figure_5_translational_evidence.tiff",
        "Figure_S1.tiff": figures / "figure_s1_af_volcano.tiff",
        "Figure_S2.tiff": figures / "figure_s2_wgcna_modules.tiff",
        "Figure_S3.tiff": figures / "figure_s3_leave_one_out.tiff",
        "Figure_S4.tiff": figures / "figure_s4_enrichment.tiff",
        "Figure_S5.tiff": figures / "figure_s5_ppi.tiff",
        "Figure_2B_UKB_outcome_cohort.tiff": figures / "figure_2b_ukb_outcome_cohort.tiff",
    }
    target_width_px = round(FULL_WIDTH_IN * OUTPUT_DPI)
    rows = []
    for target_name, source in source_by_target.items():
        if not source.exists():
            raise FileNotFoundError(f"Missing source figure: {source}")
        target = repo / "manuscript" / target_name
        with Image.open(source) as opened:
            image = opened.convert("RGB")
            target_height_px = round(image.height * target_width_px / image.width)
            if image.size != (target_width_px, target_height_px):
                image = image.resize(
                    (target_width_px, target_height_px),
                    Image.Resampling.LANCZOS,
                )
            image.save(
                target,
                format="TIFF",
                compression="tiff_lzw",
                dpi=(OUTPUT_DPI, OUTPUT_DPI),
            )
        with Image.open(target) as check:
            compression_tag = int(check.tag_v2.get(259, -1))
            width_px, height_px = check.size
        rows.append(
            {
                "submission_file": target.relative_to(repo).as_posix(),
                "source_file": source.relative_to(repo).as_posix(),
                "width_px": width_px,
                "height_px": height_px,
                "width_mm": width_px / OUTPUT_DPI * 25.4,
                "height_mm": height_px / OUTPUT_DPI * 25.4,
                "dpi": OUTPUT_DPI,
                "compression": "LZW" if compression_tag == 5 else f"TIFF tag {compression_tag}",
                "size_bytes": target.stat().st_size,
                "sha256": sha256(target),
            }
        )
    return pd.DataFrame(rows)


def write_external_metadata(paths: dict[str, Path], enrichment_queries: pd.DataFrame, ppi_summary: pd.DataFrame) -> None:
    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "gprofiler": {
            "endpoint": "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
            "organism": "hsapiens",
            "correction": "FDR",
            "query_status": json_records(enrichment_queries),
        },
        "string": {
            "version_label": "STRING 12.0 API endpoint",
            "mapping_endpoint": "https://string-db.org/api/tsv/get_string_ids",
            "network_endpoint": "https://string-db.org/api/tsv/network",
            "species": 9606,
            "required_score": 700,
            "network_type": "physical",
            "query_status": json_records(ppi_summary),
        },
    }
    paths["metadata"].write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    args = parse_args()
    repo = Path(args.repo_root).resolve()
    out = repo / "manuscript"
    figures = out / "figures"
    tables = out / "tables"
    supplement = out / "supplement"
    for directory in (figures, tables, supplement, out / "qc"):
        directory.mkdir(parents=True, exist_ok=True)
    legacy_wgcna_supplement = supplement / LEGACY_WGCNA_SUPPLEMENT
    if legacy_wgcna_supplement.exists():
        legacy_wgcna_supplement.unlink()
    np.random.seed(SEED)
    set_style()
    data = load_inputs(repo)

    source_table = build_source_table(data)
    selection_flow = build_selection_flow(data)
    wgcna_modules = data["wgcna_modules"].copy()
    genetic_table = build_genetic_table(data)
    corrected_egger = recalculate_mr_egger(data)
    random_effects = recalculate_multiplicative_random_effects(data)
    sensitivity_table = build_sensitivity_table(
        data, genetic_table, corrected_egger, random_effects
    )
    scrna_table = aggregate_scrna(data)
    translation_table = build_translation_table(data)
    write_tsv(source_table, tables / "table_1_data_sources.tsv")
    write_tsv(selection_flow, supplement / "selection_flow.tsv")
    write_tsv(genetic_table, tables / "table_2_key_genetic_evidence.tsv")
    write_tsv(translation_table, tables / "table_3_translational_prioritization.tsv")
    write_tsv(sensitivity_table, supplement / "sensitivity_summary.tsv")
    write_tsv(corrected_egger, supplement / "mr_egger_student_t_results.tsv")
    write_tsv(random_effects, supplement / "multiplicative_random_effects_ivw.tsv")
    write_tsv(scrna_table, supplement / "single_cell_localization_summary.tsv")
    write_tsv(
        data["mr"][
            data["mr"]["gene_symbol"].isin(PAPER_GENE_ORDER)
            & data["mr"]["direction"].eq("expression_to_disease")
        ],
        supplement / "forward_mr_all_priority_endpoints.tsv",
    )
    write_tsv(
        data["loo"][
            data["loo"]["gene_symbol"].isin(PAPER_GENE_ORDER)
            & data["loo"]["direction"].eq("expression_to_disease")
        ],
        supplement / "leave_one_out_priority_endpoints.tsv",
    )

    external_paths = external_snapshot_paths(out)
    enrichment, enrichment_queries = run_enrichment(data, external_paths, args.refresh_external)
    if not enrichment.empty and "query_size" in enrichment.columns:
        mapped_sizes = (
            enrichment.groupby("query_name", sort=False)["query_size"].first().to_dict()
        )
        enrichment_queries["mapped_foreground_size"] = enrichment_queries["query_name"].map(
            mapped_sizes
        )
    write_tsv(enrichment_queries, external_paths["enrichment_queries"])
    _, ppi_edges, ppi_summary = run_ppi(data, external_paths, args.refresh_external)
    write_external_metadata(external_paths, enrichment_queries, ppi_summary)

    figure_entries = []
    figure_jobs = [
        ("Figure 1", "Study design and evidence flow", make_study_design_figure(data, figures / "figure_1_study_design")),
        ("Figure 2A", "Forward MR estimates", make_mr_forest(genetic_table, figures / "figure_2_mr_forest")),
        ("Figure 3", "Endpoint-specific genetic evidence", make_genetic_heatmap(genetic_table, figures / "figure_3_genetic_evidence")),
        ("Figure 4", "Single-cell localization", make_scrna_heatmap(scrna_table, figures / "figure_4_single_cell_localization")),
        ("Figure 5", "Translational evidence", make_translation_figure(translation_table, genetic_table, figures / "figure_5_translational_evidence")),
        ("Figure S1", "AF differential expression", make_bulk_volcano(data, figures / "figure_s1_af_volcano")),
        ("Figure S2", "Signed WGCNA module effects", make_wgcna_figure(wgcna_modules, figures / "figure_s2_wgcna_modules")),
        ("Figure S3", "Leave-one-out sensitivity", make_loo_figure(data, figures / "figure_s3_leave_one_out")),
        ("Figure S4", "Functional enrichment", make_enrichment_figure(enrichment, figures / "figure_s4_enrichment")),
        ("Figure S5", "Sparse physical PPI network", make_ppi_figure(ppi_edges, figures / "figure_s5_ppi")),
    ]
    composite_paths = make_figure_2_composite(
        genetic_table,
        data["ukb_mr"],
        figures / "figure_2_composite",
    )
    figure_jobs.append(
        ("Figure 2", "Forward MR and independent outcome-cohort evaluation", composite_paths)
    )
    for figure_id, title, paths in figure_jobs:
        for path in paths:
            figure_entries.append(
                {
                    "figure_id": figure_id,
                    "title": title,
                    "path": path.relative_to(repo).as_posix(),
                    "format": path.suffix.lstrip("."),
                    "dpi": "vector" if path.suffix == ".pdf" else 300,
                    "source_scope": "frozen artifacts/core plus manuscript supplementary snapshots",
                }
            )
    for suffix in (".pdf", ".png", ".tiff"):
        ukb_figure = figures / f"figure_2b_ukb_outcome_cohort{suffix}"
        if ukb_figure.exists():
            figure_entries.append(
                {
                    "figure_id": "Figure 2B",
                    "title": "Independent outcome evaluation using public UKB-SAIGE summary statistics",
                    "path": ukb_figure.relative_to(repo).as_posix(),
                    "format": suffix.lstrip("."),
                    "dpi": "vector" if suffix == ".pdf" else 300,
                    "source_scope": "UKB-SAIGE phecode 722 public summary-statistics API",
                }
            )
    derived_umap = figures / "figure_s6_source_umap.png"
    if derived_umap.exists():
        derived_umap.unlink()
    figure_manifest = pd.DataFrame(figure_entries)
    write_tsv(figure_manifest, figures / "figure_manifest.tsv")
    submission_manifest = prepare_submission_tiffs(repo, figures)
    write_tsv(submission_manifest, out / "submission_figure_manifest.tsv")

    summary = {
        "seed": SEED,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "core_manifest_sha256": sha256(repo / "artifacts" / "core" / "manifest.tsv"),
        "core_manifest_rows": len(data["manifest"]),
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "submission_figures": {
            "count": len(submission_manifest),
            "width_mm": FULL_WIDTH_MM,
            "dpi": OUTPUT_DPI,
            "compression": "TIFF LZW",
        },
        "bulk": {
            "scope": "frozen DEG and preprocessing fields only; corrected WGCNA is reported separately",
            "deg_counts": data["bulk_summary"]["deg_counts"],
            "druggable_bulk_candidate_count": data["bulk_summary"][
                "druggable_bulk_candidate_count"
            ],
            "repeat_combat_performed": data["bulk_summary"]["repeat_combat_performed"],
        },
        "wgcna_signed_reanalysis": data["wgcna_signed_summary"],
        "genetics": {
            "formal_candidates": data["genetics_summary"]["candidate_count"],
            "screen_estimable_by_outcome": {
                str(outcome_id): int(count)
                for outcome_id, count in data["mr_screen"]
                .loc[data["mr_screen"]["status"].eq("estimable")]
                .groupby("outcome_id")
                .size()
                .items()
            },
            "mr_estimable": int((data["mr"]["status"] == "estimable").sum()),
            "eqtlgen_coloc_pp4_ge_0_8": int(
                ((data["coloc"]["status"] == "estimable") & (data["coloc"]["pp4"] >= 0.8)).sum()
            ),
            "smr_estimable": int((data["smr"]["status"] == "estimable").sum()),
            "heidi_estimable": int(
                pd.to_numeric(data["smr"]["p_heidi"], errors="coerce").notna().sum()
            ),
            "smr_bonferroni_test_counts_by_outcome": {
                str(outcome_id): int(count)
                for outcome_id, count in data["smr"]
                .loc[data["smr"]["status"].eq("estimable")]
                .groupby("outcome_id")
                .size()
                .items()
            },
            "mr_egger_student_t_estimable": int(corrected_egger["status"].eq("estimable").sum()),
            "multiplicative_random_effects_estimable": int(
                random_effects["status"].eq("estimable").sum()
            ),
            "gtex_coloc_estimable": data["gtex_summary"]["coloc_estimable"],
            "gtex_coloc_pp4_ge_0_8": data["gtex_summary"]["coloc_pp4_at_least_0_8"],
        },
        "single_cell": data["scrna_summary"],
        "drug_delivery": data["drug_summary"],
        "docking": data["docking_summary"],
        "enrichment_queries": json_records(enrichment_queries),
        "ppi_queries": json_records(ppi_summary),
        "figure_count": figure_manifest["figure_id"].nunique(),
    }
    ukb_metadata_path = supplement / "ukb_saige_722_metadata.json"
    if ukb_metadata_path.exists():
        summary["ukb_replication"] = json.loads(ukb_metadata_path.read_text(encoding="utf-8"))
    (out / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({"status": "completed", "output": str(out), "figure_count": summary["figure_count"]}))


if __name__ == "__main__":
    main()
