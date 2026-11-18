#!/usr/bin/env python3
"""Replicate prioritized IDD target MR estimates in UKB-SAIGE phecode 722."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests


SEED = 20260716
GENES = ["AKR1C2", "NT5E", "CCL28", "P2RY13", "AKR1C1"]
OUTCOME_ID = "UKB_SAIGE_722"
OUTCOME_LABEL = "UKB-SAIGE intervertebral disc disorders"
PHENOTYPE_URL = "https://pheweb.org/UKB-SAIGE/pheno/722"
API_URL = "https://pheweb.org/UKB-SAIGE/api/region/722/lz-results/"
CASES = 9241
CONTROLS = 391917
SAMPLE_SIZE = CASES + CONTROLS
EAF_DIFFERENCE_MAX = 0.05
CORE_MANIFEST_SHA256 = "85610e4e52fbad961c8e6e8055ce32ec0a188dc093bdc7033a8671c4cf5fc7da"
FULL_WIDTH_IN = 170.0 / 25.4
SIGNIFICANT_COLOR = "#0072B2"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--refresh", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> pd.Series:
    pvalues = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    adjusted = np.full(len(pvalues), np.nan)
    valid = np.where(np.isfinite(pvalues))[0]
    if not len(valid):
        return pd.Series(adjusted, index=values.index)
    order = valid[np.argsort(pvalues[valid])]
    ranked = pvalues[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.minimum(ranked, 1.0)
    return pd.Series(adjusted, index=values.index)


def save_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".png", ".tiff"):
        path = stem.with_suffix(suffix)
        kwargs: dict[str, Any] = {"bbox_inches": "tight", "facecolor": "white"}
        if suffix != ".pdf":
            kwargs["dpi"] = 300
        if suffix == ".tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, **kwargs)
    plt.close(fig)


def set_figure_style() -> None:
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


def flatten_api_data(payload: dict[str, Any]) -> list[dict[str, Any]]:
    data = payload.get("data", {})
    if not data:
        return []
    lengths = {len(value) for value in data.values() if isinstance(value, list)}
    if not lengths:
        return []
    if len(lengths) != 1:
        raise RuntimeError("UKB-SAIGE API returned non-parallel arrays")
    count = lengths.pop()
    return [
        {key: value[index] for key, value in data.items() if isinstance(value, list)}
        for index in range(count)
    ]


def load_exposure_instruments(repo: Path) -> pd.DataFrame:
    core = pd.read_csv(
        repo / "artifacts" / "core" / "genetics" / "harmonized_instruments.tsv",
        sep="\t",
    )
    core = core[
        core["direction"].eq("expression_to_disease") & core["gene_symbol"].isin(GENES)
    ].copy()
    consistency_columns = [
        "beta_exposure",
        "se_exposure",
        "effect_allele",
        "other_allele",
    ]
    uniqueness = core.groupby(["gene_symbol", "snp"])[consistency_columns].nunique(dropna=False)
    if (uniqueness > 1).any().any():
        raise RuntimeError("Exposure fields differ across repeated FinnGen outcome rows")
    core = core.drop_duplicates(["gene_symbol", "snp"])[
        [
            "gene_symbol",
            "ensembl_gene_id",
            "snp",
            "beta_exposure",
            "se_exposure",
            "effect_allele",
            "other_allele",
        ]
    ]
    exposure = pd.read_csv(
        repo / "data" / "processed" / "idd" / "genetics" / "candidate_significant_eqtl.tsv.gz",
        sep="\t",
    )
    exposure = exposure[
        exposure["gene_symbol"].isin(GENES)
        & exposure.set_index(["gene_symbol", "snp"]).index.isin(
            core.set_index(["gene_symbol", "snp"]).index
        )
    ].copy()
    exposure = exposure.drop_duplicates(["gene_symbol", "snp"])
    enriched = core.merge(
        exposure[
            [
                "gene_symbol",
                "snp",
                "chr_b37",
                "pos_b37",
                "eaf",
                "beta",
                "se",
                "f_statistic",
                "n",
                "effect_allele",
                "other_allele",
            ]
        ].rename(
            columns={
                "beta": "source_beta_exposure",
                "se": "source_se_exposure",
                "effect_allele": "source_effect_allele",
                "other_allele": "source_other_allele",
            }
        ),
        on=["gene_symbol", "snp"],
        how="left",
        validate="one_to_one",
    )
    if enriched[["chr_b37", "pos_b37", "eaf"]].isna().any().any():
        raise RuntimeError("One or more exposure instruments lack position or EAF")
    if not np.allclose(enriched["beta_exposure"], enriched["source_beta_exposure"], rtol=0, atol=1e-12):
        raise RuntimeError("Core and processed exposure beta fields differ")
    if not np.allclose(enriched["se_exposure"], enriched["source_se_exposure"], rtol=0, atol=1e-12):
        raise RuntimeError("Core and processed exposure SE fields differ")
    for allele in ("effect_allele", "other_allele"):
        if not enriched[allele].str.upper().equals(enriched[f"source_{allele}"].str.upper()):
            raise RuntimeError(f"Core and processed exposure {allele} fields differ")
    return enriched.sort_values(["gene_symbol", "snp"]).reset_index(drop=True)


def query_outcomes(instruments: pd.DataFrame, cache_path: Path, refresh: bool) -> dict[str, Any]:
    if cache_path.exists() and not refresh:
        return json.loads(cache_path.read_text(encoding="utf-8"))
    cache: dict[str, Any] = {
        "accessed_at_utc": datetime.now(timezone.utc).isoformat(),
        "api_url": API_URL,
        "records": {},
    }
    session = requests.Session()
    session.headers["User-Agent"] = "IDD manuscript outcome-cohort evaluation"
    for row in instruments.itertuples(index=False):
        chromosome = str(int(row.chr_b37))
        position = int(row.pos_b37)
        filter_expression = (
            f"chromosome in '{chromosome}' and position ge {position} and position le {position}"
        )
        response = session.get(API_URL, params={"filter": filter_expression}, timeout=90)
        response.raise_for_status()
        cache["records"][f"{row.gene_symbol}|{row.snp}"] = {
            "request_url": response.url,
            "response": response.json(),
        }
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cache, indent=2), encoding="utf-8")
    return cache


def harmonize(instruments: pd.DataFrame, cache: dict[str, Any]) -> pd.DataFrame:
    from idd_genetics_utils import harmonize_effects

    rows = []
    for exposure in instruments.to_dict(orient="records"):
        key = f"{exposure['gene_symbol']}|{exposure['snp']}"
        records = flatten_api_data(cache["records"][key]["response"])
        matches = [record for record in records if str(record.get("rsid", "")) == exposure["snp"]]
        row = {
            **exposure,
            "outcome_id": OUTCOME_ID,
            "outcome_label": OUTCOME_LABEL,
            "outcome_cases": CASES,
            "outcome_controls": CONTROLS,
            "status": "not_estimable",
            "reason": "rsid_not_found_at_exact_position",
        }
        if len(matches) != 1:
            rows.append(row)
            continue
        outcome_raw = matches[0]
        chromosome_matches = str(outcome_raw["chr"]) == str(int(exposure["chr_b37"]))
        position_matches = int(outcome_raw["position"]) == int(exposure["pos_b37"])
        if not chromosome_matches or not position_matches:
            row["reason"] = "coordinate_mismatch"
            rows.append(row)
            continue
        outcome = {
            "effect_allele": outcome_raw["alt"],
            "other_allele": outcome_raw["ref"],
            "eaf": outcome_raw["af"],
            "beta": outcome_raw["beta"],
            "se": outcome_raw["sebeta"],
        }
        exposure_for_harmonization = {
            "effect_allele": exposure["effect_allele"],
            "other_allele": exposure["other_allele"],
            "eaf": exposure["eaf"],
            "beta": exposure["beta_exposure"],
            "se": exposure["se_exposure"],
        }
        aligned, reason = harmonize_effects(exposure_for_harmonization, outcome)
        row.update(
            {
                "outcome_ref": outcome_raw["ref"],
                "outcome_alt": outcome_raw["alt"],
                "outcome_eaf": outcome_raw["af"],
                "outcome_beta_raw": outcome_raw["beta"],
                "outcome_se": outcome_raw["sebeta"],
                "outcome_pvalue_raw": outcome_raw["pvalue"],
            }
        )
        if aligned is None:
            row["reason"] = reason
            rows.append(row)
            continue
        swapped = aligned.action in ("swapped", "strand_flip_and_swap")
        aligned_outcome_eaf = 1.0 - float(outcome_raw["af"]) if swapped else float(outcome_raw["af"])
        eaf_difference = abs(float(exposure["eaf"]) - aligned_outcome_eaf)
        row.update(
            {
                "beta_outcome": aligned.beta_outcome,
                "se_outcome": aligned.se_outcome,
                "harmonization_action": aligned.action,
                "aligned_outcome_eaf": aligned_outcome_eaf,
                "eaf_difference": eaf_difference,
                "sensitivity_included": True,
                "primary_included": eaf_difference <= EAF_DIFFERENCE_MAX,
                "status": "included"
                if eaf_difference <= EAF_DIFFERENCE_MAX
                else "sensitivity_only",
                "reason": ""
                if eaf_difference <= EAF_DIFFERENCE_MAX
                else "eaf_difference_gt_0.05",
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_effect(effect: dict[str, float]) -> dict[str, float]:
    return {
        **effect,
        "odds_ratio": math.exp(effect["beta"]),
        "or_ci_lower": math.exp(effect["ci_lower"]),
        "or_ci_upper": math.exp(effect["ci_upper"]),
    }


def run_mr(harmonized: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    from idd_genetics_utils import (
        effective_binary_n,
        ivw_fixed,
        ivw_multiplicative_random_effects,
        mr_egger,
        observed_r2_from_z,
        wald_ratio,
        weighted_median,
    )

    modes = {"primary_eaf_qc": "primary_included", "all_allele_aligned": "sensitivity_included"}
    mr_rows = []
    sensitivity_rows = []
    loo_rows = []
    for mode, column in modes.items():
        for gene in GENES:
            group = harmonized[
                harmonized["gene_symbol"].eq(gene) & harmonized[column].fillna(False).astype(bool)
            ].copy()
            values = group[
                ["beta_exposure", "se_exposure", "beta_outcome", "se_outcome"]
            ].to_dict(orient="records")
            base = {
                "gene_symbol": gene,
                "outcome_id": OUTCOME_ID,
                "outcome_label": OUTCOME_LABEL,
                "analysis_mode": mode,
                "n_instruments": len(values),
                "method": "not_estimable",
                "status": "not_estimable",
                "reason": "no_instruments_after_harmonization",
            }
            if values:
                if len(values) == 1:
                    effect = wald_ratio(
                        values[0]["beta_exposure"],
                        values[0]["se_exposure"],
                        values[0]["beta_outcome"],
                        values[0]["se_outcome"],
                    )
                    method = "wald_ratio"
                else:
                    effect = ivw_fixed(values)
                    method = "ivw_fixed"
                exposure_r2 = sum(
                    observed_r2_from_z(float(row.beta_exposure / row.se_exposure), float(row.n))
                    for row in group.itertuples(index=False)
                )
                effective_n = effective_binary_n(CASES, CONTROLS)
                outcome_r2 = sum(
                    observed_r2_from_z(float(row.beta_outcome / row.se_outcome), effective_n)
                    for row in group.itertuples(index=False)
                )
                base.update(
                    {
                        **summarize_effect(effect),
                        "method": method,
                        "status": "estimable",
                        "reason": "",
                        "exposure_r2": exposure_r2,
                        "outcome_r2": outcome_r2,
                        "steiger_direction_supported": exposure_r2 > outcome_r2,
                    }
                )
            mr_rows.append(base)

            for method in (
                "weighted_median",
                "mr_egger",
                "multiplicative_random_effects_ivw",
            ):
                minimum_instruments = 2 if method == "multiplicative_random_effects_ivw" else 3
                sensitivity = {
                    "gene_symbol": gene,
                    "outcome_id": OUTCOME_ID,
                    "analysis_mode": mode,
                    "n_instruments": len(values),
                    "method": method,
                    "status": "not_estimable",
                    "reason": f"fewer_than_{minimum_instruments}_instruments",
                }
                if len(values) >= minimum_instruments:
                    try:
                        if method == "weighted_median":
                            instrument_signature = "|".join(sorted(group["snp"].astype(str)))
                            seed_material = f"{gene}|{OUTCOME_ID}|{instrument_signature}|{SEED}"
                            seed = int(hashlib.sha256(seed_material.encode("ascii")).hexdigest()[:8], 16)
                            result = weighted_median(values, seed=seed)
                            sensitivity.update(summarize_effect(result))
                        elif method == "mr_egger":
                            result = mr_egger(values)
                            slope = {
                                key.replace("slope_", ""): value
                                for key, value in result.items()
                                if key.startswith("slope_")
                            }
                            sensitivity.update(summarize_effect(slope))
                            sensitivity.update(
                                {
                                    "egger_intercept": result["intercept"],
                                    "egger_intercept_se": result["intercept_se"],
                                    "egger_intercept_pvalue": result["intercept_pvalue"],
                                    "q": result["q"],
                                    "q_df": result["q_df"],
                                    "q_pvalue": result["q_pvalue"],
                                }
                            )
                        else:
                            result = ivw_multiplicative_random_effects(values)
                            sensitivity.update(summarize_effect(result))
                            sensitivity.update(
                                {
                                    "q": result["q"],
                                    "q_df": result["q_df"],
                                    "q_pvalue": result["q_pvalue"],
                                    "dispersion": result["dispersion"],
                                }
                            )
                        sensitivity.update({"status": "estimable", "reason": ""})
                    except ValueError as exc:
                        sensitivity["reason"] = str(exc)
                sensitivity_rows.append(sensitivity)

            if len(values) >= 2:
                for excluded in group["snp"]:
                    subset = group[group["snp"].ne(excluded)]
                    subset_values = subset[
                        ["beta_exposure", "se_exposure", "beta_outcome", "se_outcome"]
                    ].to_dict(orient="records")
                    if len(subset_values) == 1:
                        effect = wald_ratio(
                            subset_values[0]["beta_exposure"],
                            subset_values[0]["se_exposure"],
                            subset_values[0]["beta_outcome"],
                            subset_values[0]["se_outcome"],
                        )
                        method = "wald_ratio"
                    else:
                        effect = ivw_fixed(subset_values)
                        method = "ivw_fixed"
                    loo_rows.append(
                        {
                            "gene_symbol": gene,
                            "outcome_id": OUTCOME_ID,
                            "analysis_mode": mode,
                            "excluded_snp": excluded,
                            "remaining_instruments": len(subset_values),
                            "method": method,
                            **summarize_effect(effect),
                        }
                    )
    mr = pd.DataFrame(mr_rows)
    mr["fdr"] = np.nan
    for _, index in mr[mr["status"].eq("estimable")].groupby("analysis_mode").groups.items():
        mr.loc[index, "fdr"] = bh_adjust(mr.loc[index, "pvalue"])
    return mr, pd.DataFrame(sensitivity_rows), pd.DataFrame(loo_rows)


def plot_replication_axis(
    ax: plt.Axes,
    mr: pd.DataFrame,
    *,
    show_legend: bool = True,
) -> None:
    frame = mr[
        mr["analysis_mode"].eq("primary_eaf_qc") & mr["status"].eq("estimable")
    ].copy()
    order = {gene: index for index, gene in enumerate(GENES)}
    frame["order"] = frame["gene_symbol"].map(order)
    frame = frame.sort_values("order", ascending=False).reset_index(drop=True)
    y = np.arange(len(frame))
    for index, row in frame.iterrows():
        significant = bool(row["fdr"] < 0.05)
        ax.plot(
            [row["or_ci_lower"], row["or_ci_upper"]],
            [index, index],
            color=SIGNIFICANT_COLOR,
            lw=1.8,
        )
        marker = "D" if row["method"] == "wald_ratio" else "o"
        ax.scatter(
            row["odds_ratio"],
            index,
            facecolor=SIGNIFICANT_COLOR if significant else "white",
            edgecolor=SIGNIFICANT_COLOR,
            linewidth=1.3,
            s=50,
            marker=marker,
            zorder=3,
        )
        ax.text(
            1.56,
            index,
            f"{row['odds_ratio']:.3f} ({row['or_ci_lower']:.3f}\N{EN DASH}{row['or_ci_upper']:.3f})",
            va="center",
            fontsize=9,
        )
    ax.axvline(1, color="#777777", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(frame["gene_symbol"])
    ax.set_xlim(0.65, 1.75)
    ax.set_xticks([0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 1.7])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("Odds ratio per genetically predicted unit increase in expression")
    ax.set_title("Independent outcome cohort: UKB-SAIGE phecode 722", loc="left", pad=10)
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color=SIGNIFICANT_COLOR,
            markerfacecolor=SIGNIFICANT_COLOR,
            markeredgecolor=SIGNIFICANT_COLOR,
            label="BH-FDR < 0.05",
            markersize=6,
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color=SIGNIFICANT_COLOR,
            markerfacecolor="white",
            markeredgecolor=SIGNIFICANT_COLOR,
            label="BH-FDR >= 0.05",
            markersize=6,
        ),
    ]
    if show_legend:
        ax.legend(
            handles=handles,
            loc="upper right",
            bbox_to_anchor=(1.0, 1.18),
            ncol=2,
            frameon=False,
            borderaxespad=0,
        )


def make_replication_figure(mr: pd.DataFrame, stem: Path) -> None:
    set_figure_style()
    fig, ax = plt.subplots(figsize=(FULL_WIDTH_IN, 3.8))
    plot_replication_axis(ax, mr)
    save_figure(fig, stem)


def main() -> None:
    args = parse_args()
    repo = Path(args.repo_root).resolve()
    sys.path.insert(0, str(repo / "py"))
    manifest = repo / "artifacts" / "core" / "manifest.tsv"
    if sha256(manifest) != CORE_MANIFEST_SHA256:
        raise RuntimeError("Frozen core manifest drifted before UKB outcome-cohort analysis")
    supplement = repo / "manuscript" / "supplement"
    figures = repo / "manuscript" / "figures"
    supplement.mkdir(parents=True, exist_ok=True)
    instruments = load_exposure_instruments(repo)
    if len(instruments) != 14:
        raise RuntimeError(f"Expected 14 unique priority-gene instruments, found {len(instruments)}")
    cache_path = supplement / "ukb_saige_722_api_cache.json"
    cache = query_outcomes(instruments, cache_path, args.refresh)
    harmonized = harmonize(instruments, cache)
    mr, sensitivity, loo = run_mr(harmonized)
    harmonized.to_csv(
        supplement / "ukb_saige_722_instruments.tsv", sep="\t", index=False, na_rep="not_estimable"
    )
    mr.to_csv(
        supplement / "ukb_saige_722_mr.tsv", sep="\t", index=False, na_rep="not_estimable"
    )
    sensitivity.to_csv(
        supplement / "ukb_saige_722_sensitivity.tsv", sep="\t", index=False, na_rep="not_estimable"
    )
    loo.to_csv(
        supplement / "ukb_saige_722_leave_one_out.tsv", sep="\t", index=False, na_rep="not_estimable"
    )
    make_replication_figure(mr, figures / "figure_2b_ukb_outcome_cohort")
    metadata = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED,
        "phenotype": "Phecode 722: intervertebral disc disorders",
        "phenotype_url": PHENOTYPE_URL,
        "api_url": API_URL,
        "access_route": "public UKB-SAIGE PheWeb aggregate-summary-statistics API",
        "source_ukb_application": "24460",
        "present_study_participant_level_ukb_access": False,
        "api_cache_distribution": "internal provenance only; excluded from external submission package",
        "population": "White British",
        "gwas_method": "SAIGE",
        "genome_build": "GRCh37",
        "cases": CASES,
        "controls": CONTROLS,
        "sample_size": SAMPLE_SIZE,
        "effective_sample_size": 4.0 / (1.0 / CASES + 1.0 / CONTROLS),
        "effect_allele": "ALT",
        "primary_eaf_difference_max": EAF_DIFFERENCE_MAX,
        "primary_instrument_count": int(harmonized["primary_included"].fillna(False).sum()),
        "sensitivity_instrument_count": int(harmonized["sensitivity_included"].fillna(False).sum()),
        "core_manifest_sha256": sha256(manifest),
        "interpretation_limit": "Independent outcome-cohort evaluation using the same eQTLGen exposure instruments; Phecode 722 is broader than imaging-defined IDD.",
    }
    (supplement / "ukb_saige_722_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": "completed",
                "instruments": len(harmonized),
                "primary_instruments": metadata["primary_instrument_count"],
                "estimable_genes": int(
                    (
                        mr["analysis_mode"].eq("primary_eaf_qc")
                        & mr["status"].eq("estimable")
                    ).sum()
                ),
            }
        )
    )


if __name__ == "__main__":
    main()
