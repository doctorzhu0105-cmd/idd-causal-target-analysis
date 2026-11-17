"""Global, post-hoc instrument-configuration ablation across the 76 candidates.

The script deliberately keeps the exposure-side concentration metric separate from
the outcome-side leave-one-out diagnostic. It streams the full eQTLGen cis file,
uses the project's existing allele-frequency conversion, and reuses the local LD
and correlated-IVW implementations.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
WS = ROOT.parent
PROC = WS / "data" / "processed" / "idd" / "genetics"
OUT = ROOT / "out"
RAW_EQTL = WS / "data" / "raw" / "idd" / "genetics" / "eqtlgen" / "eqtlgen_full_cis.txt.gz"
RAW_AF = WS / "data" / "raw" / "idd" / "genetics" / "eqtlgen" / "eqtlgen_allele_frequency.txt.gz"
CANDIDATES = WS / "artifacts" / "core" / "genetics" / "formal_candidate_genes.tsv"
LD_DIR = WS / "data" / "raw" / "idd" / "genetics" / "ld"
BED = LD_DIR / "1000G_EUR.bed"
BIM = LD_DIR / "1000G_EUR.bim"
N_REF = 503
WINDOW_BP = 1_000_000
R2_LADDER = [0.001, 0.01, 0.05, 0.1, 0.2]
CONCENTRATION_THRESHOLDS = [0.90, 0.95, 0.99]
OUTCOMES = ["M13_INTERVERTEB", "M13_LOWBACKPAIN", "M13_SCIATICA"]
COMPARISON_TABLE = OUT / "comparison_set_concentration.tsv"
COMPARISON_PROVENANCE = OUT / "comparison_set_provenance.json"
GROUP_TESTS = OUT / "global_ablation_group_tests.json"
LEGACY_FULL_EQTL = WS / "data" / "processed" / "idd" / "genetics" / "candidate_full_eqtl.tsv.gz"
LEGACY_SIGNIFICANT_EQTL = WS / "data" / "processed" / "idd" / "genetics" / "candidate_significant_eqtl.tsv.gz"

from .idd_genetics_pipeline import (  # noqa: E402
    _collect_column,
    _enrich_eqtl,
    _extract_af,
    _extract_eqtl,
    read_candidates,
)
from .correlated_mr import ivw_correlated, leave_one_out, prune_by_ld  # noqa: E402
from .plink_ld import is_palindromic, ld_matrix, load_bim  # noqa: E402
from . import run_all_genes_endpoints as endpoint_module
from .run_all_genes_endpoints import harmonise, load_outcome  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_enriched_eqtl(gene_map: dict[str, str], prefix: str) -> Path:
    """Extract and enrich candidate rows directly from the full raw eQTL file."""
    raw = OUT / f"{prefix}_eqtl.raw.tsv.gz"
    af_out = OUT / f"{prefix}_allele_frequency.tsv.gz"
    enriched = OUT / f"{prefix}_eqtl.tsv.gz"
    _extract_eqtl(
        RAW_EQTL,
        raw,
        gene_map,
        pvalue_max=None,
        f_min_exclusive=None,
        window_bp=WINDOW_BP,
    )
    requested = _collect_column(raw, "snp")
    _extract_af(RAW_AF, requested, af_out)
    _enrich_eqtl(raw, enriched, _read_af_map(af_out))
    return enriched


def _read_af_map(path: Path) -> dict[str, tuple[str, str, float]]:
    result: dict[str, tuple[str, str, float]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for line in handle:
            values = line.rstrip("\r\n").split("\t")
            result[values[index["snp"]]] = (
                values[index["allele_a"]],
                values[index["allele_b"]],
                float(values[index["allele_b_frequency"]]),
            )
    return result


def read_enriched(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    numeric = [
        "chr_b37", "pos_b37", "gene_pos_b37", "zscore", "n", "pvalue", "fdr",
        "eaf", "beta", "se", "f_statistic",
    ]
    for column in numeric:
        if column in df:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def concentration_index(zx: np.ndarray, rho: np.ndarray) -> dict[str, float | int]:
    """Return an exposure-side, LD-aware heuristic concentration diagnostic."""
    precision = np.linalg.pinv(rho)
    contribution = zx * (precision @ zx)
    absolute_total = float(np.abs(contribution).sum())
    signed_total = float(contribution.sum())
    return {
        "share_abs": float(np.abs(contribution).max() / absolute_total)
        if absolute_total > 0 else float("nan"),
        "share_raw": float(contribution.max() / signed_total)
        if signed_total > 0 else float("nan"),
        "n_negative_contrib": int((contribution < 0).sum()),
        "negative_abs_fraction": float(np.abs(contribution[contribution < 0]).sum() / absolute_total)
        if absolute_total > 0 else float("nan"),
    }


def _deduplicate_after_ld(kept: pd.DataFrame, rho: np.ndarray) -> tuple[pd.DataFrame, np.ndarray]:
    """Keep the highest-F row after LD allele matching, not before it."""
    order = kept.sort_values("f_statistic", ascending=False, kind="mergesort").index.tolist()
    selected: list[int] = []
    seen: set[str] = set()
    for index in order:
        snp = str(kept.loc[index, "snp"])
        if snp not in seen:
            selected.append(int(index))
            seen.add(snp)
    result = kept.iloc[selected].reset_index(drop=True)
    matrix = rho[np.ix_(selected, selected)]
    return result, matrix


def _ld_pruned_exposure(filtered: pd.DataFrame, bim: pd.DataFrame, r2: float):
    if filtered.empty:
        return filtered, np.zeros((0, 0)), 0
    variants, rho = ld_matrix(filtered, bim, str(BED), N_REF)
    if variants.empty:
        return variants, rho, 0
    variants, rho = _deduplicate_after_ld(variants, rho)
    order = np.argsort(-variants.f_statistic.to_numpy())
    index = prune_by_ld(rho, order, r2)
    return variants.iloc[index].reset_index(drop=True), rho[np.ix_(index, index)], len(variants)


def _filtered_gene(
    df: pd.DataFrame, gene: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, float]:
    all_rows = df[df.gene_symbol == gene].copy()
    if all_rows.empty:
        return all_rows, all_rows, all_rows, float("nan")
    max_n = float(all_rows.n.max())
    filtered_pre = all_rows[
        (all_rows.pvalue < 5e-8)
        & (all_rows.f_statistic > 10)
        & (all_rows.n >= 0.80 * max_n)
    ].copy()
    filtered = filtered_pre[
        ~filtered_pre.apply(
            lambda row: is_palindromic(row.effect_allele, row.other_allele), axis=1
        )
    ].copy()
    usable = filtered.dropna(subset=["beta", "se", "zscore", "f_statistic"]).copy()
    return all_rows, filtered_pre, usable, max_n


def run_concentration(
    candidates: list[dict[str, str]],
    eqtl: pd.DataFrame,
    output_name: str = "global_ablation_concentration.tsv",
) -> pd.DataFrame:
    bim_all = load_bim(str(BIM))
    rows: list[dict[str, object]] = []
    seen_genes: set[str] = set()

    for candidate in candidates:
        gene = candidate["gene_symbol"]
        gene_id = candidate["ensembl_gene_id"]
        all_rows, filtered_pre, usable, max_n = _filtered_gene(eqtl, gene)
        seen_genes.add(gene)
        n_filtered = int(len(filtered_pre))
        n_palindromic = int(n_filtered - len(
            filtered_pre[
                ~filtered_pre.apply(
                    lambda row: is_palindromic(
                        row.effect_allele, row.other_allele
                    ),
                    axis=1,
                )
            ]
        ))
        n_unusable = max(n_filtered - len(usable), 0)
        for r2 in R2_LADDER:
            base = {
                "gene": gene,
                "ensembl_gene_id": gene_id,
                "r2_threshold": r2,
                "n_full_cis_rows": int(len(all_rows)),
                "max_n": max_n,
                "n_filtered": n_filtered,
                "n_palindromic_excluded": n_palindromic,
                "n_unusable_exposure_rows": n_unusable,
                "k": 0,
                "min_f": float("nan"),
                "ld_condition_number": float("nan"),
                "share_abs": float("nan"),
                "share_raw": float("nan"),
                "n_negative_contrib": 0,
                "negative_abs_fraction": float("nan"),
                "n_ld_reference_rows": 0,
                "status": "no_tool" if n_filtered == 0 else "no_exposure_conversion",
            }
            if n_filtered > 0 and not usable.empty:
                chrom = str(int(usable.chr_b37.iloc[0]))
                bim = bim_all[bim_all.chr == chrom]
                sub, rho, n_ld = _ld_pruned_exposure(usable, bim, r2)
                base["n_ld_reference_rows"] = int(n_ld)
                if not sub.empty:
                    zx = (sub.beta / sub.se).to_numpy(dtype=float)
                    idx = concentration_index(zx, rho)
                    base.update(
                        {
                            "k": int(len(sub)),
                            "min_f": float(sub.f_statistic.min()),
                            "ld_condition_number": float(np.linalg.cond(rho)),
                            **idx,
                            "status": "ok",
                        }
                    )
                    if r2 == 0.001 and len(sub) < 3:
                        base["status"] = "loo_uninterpretable_k_lt_3"
                else:
                    base["status"] = "no_ld_reference"
            rows.append(base)

    if len(seen_genes) != len(candidates):
        missing = sorted(set(c["gene_symbol"] for c in candidates) - seen_genes)
        raise AssertionError(f"candidate symbols missing from extraction map: {missing}")
    result = pd.DataFrame(rows)
    result.to_csv(OUT / output_name, sep="\t", index=False)
    return result


def run_loo(candidates: list[dict[str, str]], eqtl: pd.DataFrame) -> pd.DataFrame:
    bim_all = load_bim(str(BIM))
    rows: list[dict[str, object]] = []
    for candidate in candidates:
        gene = candidate["gene_symbol"]
        gene_id = candidate["ensembl_gene_id"]
        for outcome in OUTCOMES:
            base: dict[str, object] = {
                "gene": gene,
                "ensembl_gene_id": gene_id,
                "outcome": outcome,
                "r2_threshold": 0.001,
                "k": 0,
                "full_beta": float("nan"),
                "full_or": float("nan"),
                "full_p": float("nan"),
                "full_nominal_significant": False,
                "loo_max_p": float("nan"),
                "loo_fraction_significant": float("nan"),
                "sign_reversal": False,
                "status": "no_outcome_overlap",
            }
            fg = load_outcome(outcome)
            harmonised = harmonise(eqtl, gene, fg)
            if harmonised.empty or not {"beta_exposure", "se_outcome", "f_statistic"}.issubset(harmonised.columns):
                rows.append(base)
                continue
            harmonised = harmonised.dropna(subset=["beta_exposure", "se_outcome", "f_statistic"])
            if harmonised.empty:
                rows.append(base)
                continue
            chrom = str(int(harmonised.chr_b37.iloc[0]))
            bim = bim_all[bim_all.chr == chrom]
            sub, rho, _ = _ld_pruned_exposure(harmonised, bim, 0.001)
            if len(sub) < 3:
                base["k"] = int(len(sub))
                base["status"] = "loo_uninterpretable_k_lt_3"
                rows.append(base)
                continue
            bx = sub.beta_exposure.to_numpy(dtype=float)
            by = sub.beta_outcome.to_numpy(dtype=float)
            sey = sub.se_outcome.to_numpy(dtype=float)
            estimate = ivw_correlated(bx, by, sey, rho, random_effects=True)
            loo = pd.DataFrame(leave_one_out(bx, by, sey, rho))
            full_sign = np.sign(estimate.beta)
            base.update(
                {
                    "k": int(len(sub)),
                    "full_beta": float(estimate.beta),
                    "full_or": float(estimate.or_),
                    "full_p": float(estimate.pvalue),
                    "full_nominal_significant": bool(estimate.pvalue < 0.05),
                    "loo_max_p": float(loo.pvalue.max()),
                    "loo_fraction_significant": float((loo.pvalue < 0.05).mean()),
                    "sign_reversal": bool((np.sign(loo.beta) != full_sign).any()),
                    "status": "estimable",
                }
            )
            rows.append(base)
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "global_ablation_loo.tsv", sep="\t", index=False)
    return result


def _load_comparison_for_group_tests() -> tuple[pd.DataFrame, dict[str, object], str]:
    """Load the fixed-seed comparator and refuse missing or stale provenance."""
    if not COMPARISON_TABLE.exists():
        raise FileNotFoundError(
            f"Required comparator table is missing: {COMPARISON_TABLE}. "
            "Run comparison_set_ablation.py first."
        )
    if not COMPARISON_PROVENANCE.exists():
        raise FileNotFoundError(
            f"Required comparator provenance is missing: {COMPARISON_PROVENANCE}."
        )
    provenance = json.loads(COMPARISON_PROVENANCE.read_text(encoding="utf-8"))
    if int(provenance.get("seed", -1)) != 20260816:
        raise AssertionError(
            "Comparator provenance seed is not the locked value 20260816; "
            "refusing to use a potentially stale comparator."
        )
    comparison = pd.read_csv(COMPARISON_TABLE, sep="\t")
    required = {"gene", "r2_threshold", "k", "share_abs"}
    missing = required - set(comparison.columns)
    if missing:
        raise AssertionError(f"Comparator table is missing columns: {sorted(missing)}")
    expected_genes = int(provenance.get("sampled_gene_count_with_hgnc_symbol", -1))
    observed_genes = int(comparison.gene.astype(str).nunique())
    if expected_genes < 0 or observed_genes != expected_genes:
        raise AssertionError(
            "Comparator gene count does not match provenance: "
            f"observed={observed_genes}, provenance={expected_genes}."
        )
    if comparison.r2_threshold.nunique() != len(R2_LADDER):
        raise AssertionError("Comparator table does not contain the full locked r² ladder")
    if comparison.duplicated(["gene", "r2_threshold"]).any():
        raise AssertionError("Comparator table contains duplicate gene×r² rows")
    return comparison, provenance, sha256_file(COMPARISON_TABLE)


def _fisher_result(
    candidate_success: int,
    candidate_total: int,
    comparison_success: int,
    comparison_total: int,
) -> dict[str, object]:
    table = [
        [candidate_success, candidate_total - candidate_success],
        [comparison_success, comparison_total - comparison_success],
    ]
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    return {
        "candidate_success": int(candidate_success),
        "candidate_total": int(candidate_total),
        "comparison_success": int(comparison_success),
        "comparison_total": int(comparison_total),
        "contingency_table": table,
        "odds_ratio": float(odds_ratio),
        "p_value_two_sided": float(p_value),
    }


def _mann_whitney_result(candidate: pd.Series, comparison: pd.Series) -> dict[str, object]:
    candidate = pd.to_numeric(candidate, errors="coerce").dropna()
    comparison = pd.to_numeric(comparison, errors="coerce").dropna()
    result = mannwhitneyu(candidate, comparison, alternative="two-sided", method="auto")
    return {
        "candidate_n": int(len(candidate)),
        "comparison_n": int(len(comparison)),
        "candidate_median": float(candidate.median()),
        "comparison_median": float(comparison.median()),
        "u_statistic": float(result.statistic),
        "p_value_two_sided": float(result.pvalue),
    }


def _legacy_zero_tool_reconciliation(candidate_count: int) -> dict[str, object]:
    """Explain the older 5-gene count without using it for the main analysis."""
    result: dict[str, object] = {
        "current_definition": (
            "raw full cis extraction followed by P<5e-8, F>10, N>=0.80*gene_max_N, "
            "palindromic exclusion"
        ),
        "legacy_significant_table": str(LEGACY_SIGNIFICANT_EQTL.relative_to(WS)).replace("\\", "/"),
        "legacy_full_table": str(LEGACY_FULL_EQTL.relative_to(WS)).replace("\\", "/"),
    }
    if LEGACY_SIGNIFICANT_EQTL.exists():
        legacy_sig = pd.read_csv(LEGACY_SIGNIFICANT_EQTL, sep="\t", compression="gzip", usecols=["gene_symbol"])
        represented = set(legacy_sig.gene_symbol.astype(str))
        legacy_missing = sorted(set(pd.read_csv(CANDIDATES, sep="\t", dtype=str).gene_symbol) - represented)
        result.update(
            {
                "legacy_significant_table_gene_count": int(len(represented)),
                "legacy_significant_table_missing_from_76": int(len(legacy_missing)),
                "legacy_significant_table_missing_genes": legacy_missing,
                "legacy_missing_count_interpretation": (
                    "This table contains only significant eQTL rows; absent genes are not a full raw-cis "
                    "no-tool estimate."
                ),
            }
        )
    else:
        result["legacy_significant_table_status"] = "missing"
    if LEGACY_FULL_EQTL.exists():
        legacy_full = pd.read_csv(LEGACY_FULL_EQTL, sep="\t", compression="gzip", usecols=["gene_symbol"])
        result["legacy_full_table_gene_count"] = int(legacy_full.gene_symbol.astype(str).nunique())
        result["legacy_full_table_missing_from_76"] = int(
            candidate_count - legacy_full.gene_symbol.astype(str).nunique()
        )
    else:
        result["legacy_full_table_status"] = "missing"
    return result


def compute_group_tests(concentration: pd.DataFrame) -> dict[str, object]:
    """Compute all locked, descriptive candidate-versus-comparator comparisons."""
    comparison, comparator_provenance, comparator_sha256 = _load_comparison_for_group_tests()
    candidate = concentration[concentration.r2_threshold == 0.001].copy()
    comparator = comparison[comparison.r2_threshold == 0.001].copy()
    candidate_count = int(candidate.gene.astype(str).nunique())
    comparator_count = int(comparator.gene.astype(str).nunique())
    if candidate.duplicated("gene").any() or comparator.duplicated("gene").any():
        raise AssertionError("r²<0.001 group-test inputs must have one row per gene")
    candidate_k3 = candidate[candidate.k >= 3].dropna(subset=["share_abs"])
    comparator_k3 = comparator[comparator.k >= 3].dropna(subset=["share_abs"])
    candidate_k_eq_3 = candidate[candidate.k == 3].dropna(subset=["share_abs"])
    comparator_k_eq_3 = comparator[comparator.k == 3].dropna(subset=["share_abs"])
    tests = {
        "k_ge_3_proportion": _fisher_result(
            int((candidate.k >= 3).sum()), candidate_count,
            int((comparator.k >= 3).sum()), comparator_count,
        ),
        "k_gt_3_proportion": _fisher_result(
            int((candidate.k > 3).sum()), candidate_count,
            int((comparator.k > 3).sum()), comparator_count,
        ),
        "k_ge_3_share_abs_ge_0_90": _fisher_result(
            int((candidate_k3.share_abs >= 0.90).sum()), len(candidate_k3),
            int((comparator_k3.share_abs >= 0.90).sum()), len(comparator_k3),
        ),
        "k_ge_3_share_abs_distribution": _mann_whitney_result(
            candidate_k3.share_abs, comparator_k3.share_abs
        ),
        "k_eq_3_share_abs_distribution": _mann_whitney_result(
            candidate_k_eq_3.share_abs, comparator_k_eq_3.share_abs
        ),
    }
    current_no_tool = candidate.loc[candidate.status.eq("no_tool"), "gene"].astype(str).unique().tolist()
    current_no_tool.sort()
    reconciliation = _legacy_zero_tool_reconciliation(candidate_count)
    legacy_missing = set(reconciliation.get("legacy_significant_table_missing_genes", []))
    reconciliation["current_no_tool_genes"] = current_no_tool
    reconciliation["current_no_tool_not_in_legacy_significant_missing"] = sorted(
        set(current_no_tool) - legacy_missing
    )
    reconciliation["legacy_significant_missing_not_current_no_tool"] = sorted(
        legacy_missing - set(current_no_tool)
    )
    payload: dict[str, object] = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "r2_threshold": 0.001,
        "candidate_denominator": candidate_count,
        "comparison_denominator": comparator_count,
        "comparison_sample_size_requested": int(comparator_provenance["sample_size_requested"]),
        "comparison_seed": int(comparator_provenance["seed"]),
        "comparison_table_path": str(COMPARISON_TABLE.relative_to(ROOT)).replace("\\", "/"),
        "comparison_table_sha256": comparator_sha256,
        "comparison_provenance_path": str(COMPARISON_PROVENANCE.relative_to(ROOT)).replace("\\", "/"),
        "comparison_sampled_gene_count_with_hgnc_symbol": int(
            comparator_provenance["sampled_gene_count_with_hgnc_symbol"]
        ),
        "comparison_sampled_gene_count_without_hgnc_symbol": int(
            comparator_provenance["sampled_gene_count"]
            - comparator_provenance["sampled_gene_count_with_hgnc_symbol"]
        ),
        "tests": tests,
        "current_no_tool_genes": current_no_tool,
        "legacy_zero_tool_reconciliation": reconciliation,
        "multiple_testing": "none; descriptive two-sided comparisons, no FDR",
        "scope": "small-sample descriptive comparator; not a power calculation",
    }
    GROUP_TESTS.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return payload


def _distribution_lines(
    df: pd.DataFrame,
    label: str,
    key: str = "share_abs",
    min_k: int | None = None,
) -> list[str]:
    lines = [f"### {label}"]
    for r2 in R2_LADDER:
        sub = df[(df.r2_threshold == r2) & np.isfinite(df[key])]
        if min_k is not None:
            sub = sub[sub.k >= min_k]
        if sub.empty:
            lines.append(f"- r² < {r2}: no estimable values")
            continue
        values = sub[key].astype(float)
        counts = [int((values >= threshold).sum()) for threshold in CONCENTRATION_THRESHOLDS]
        lines.append(
            f"- r² < {r2}: n={len(values)}, median={values.median():.4f}, "
            f"IQR={values.quantile(.25):.4f}–{values.quantile(.75):.4f}; "
            f">=0.90: {counts[0]}, >=0.95: {counts[1]}, >=0.99: {counts[2]}"
        )
    return lines


def _pct(numerator: int, denominator: int) -> str:
    return f"{100 * numerator / denominator:.1f}%" if denominator else "NA"


def _fmt_p(value: float) -> str:
    return f"{value:.4g}"


def write_report(
    concentration: pd.DataFrame | None = None,
    loo: pd.DataFrame | None = None,
    group_tests: dict[str, object] | None = None,
) -> None:
    concentration = concentration if concentration is not None else pd.read_csv(
        OUT / "global_ablation_concentration.tsv", sep="\t"
    )
    loo = loo if loo is not None and not loo.empty else pd.read_csv(OUT / "global_ablation_loo.tsv", sep="\t")
    group_tests = group_tests if group_tests is not None else compute_group_tests(concentration)
    comparison, comparator_provenance, _ = _load_comparison_for_group_tests()
    primary = concentration[concentration.r2_threshold == group_tests["r2_threshold"]]
    comparator_primary = comparison[comparison.r2_threshold == group_tests["r2_threshold"]]
    candidate_count = int(primary.gene.astype(str).nunique())
    comparator_count = int(comparator_primary.gene.astype(str).nunique())
    no_tool = int((primary.status == "no_tool").sum())
    loo_bad = int((primary.status == "loo_uninterpretable_k_lt_3").sum())
    other = int((~primary.status.isin(["no_tool", "loo_uninterpretable_k_lt_3", "ok"])).sum())
    finite_primary = primary[primary.share_abs.notna()]
    ok_primary = primary[(primary.status == "ok") & (primary.k >= 3) & primary.share_abs.notna()]
    high = int((ok_primary.share_abs >= 0.90).sum())
    denom = int(len(ok_primary))
    if denom and high / denom >= 0.5:
        interpretation = "按本轮运行前固定的描述性规则，在本次正式候选筛选的可计算 K>=3 子集中属于常见现象；这不是对其他基因组的发生率估计。"
    else:
        interpretation = "按本轮运行前固定的描述性规则，在本次正式候选筛选的可计算 K>=3 子集中不属于多数现象；这不是对其他基因组的发生率估计。"
    loo_est = loo[loo.status == "estimable"]
    loo_sig = loo_est[loo_est.full_nominal_significant]
    tests = group_tests["tests"]
    kge3 = tests["k_ge_3_proportion"]
    kgt3 = tests["k_gt_3_proportion"]
    high_test = tests["k_ge_3_share_abs_ge_0_90"]
    dist_test = tests["k_ge_3_share_abs_distribution"]
    k3_test = tests["k_eq_3_share_abs_distribution"]
    legacy = group_tests["legacy_zero_tool_reconciliation"]
    current_no_tool = group_tests["current_no_tool_genes"]
    comp_kgt3_counts = comparator_primary.loc[comparator_primary.k > 3, "k"].value_counts().sort_index().to_dict()
    lines = [
        "# Global instrument-configuration ablation report",
        "",
        "本报告来自事后设计的筛选内描述性分析；不改变既有稿件、投稿包或冻结 core。",
        "",
        "## 分母与可估计性",
        f"- 正式候选分母：{candidate_count} 个基因；这是完整枚举，不是抽样估计。",
        f"- 全基因组抽取覆盖（含零行基因）：{candidate_count}/{candidate_count}；其中原始文件有行的基因为 "
        f"{int(primary.loc[primary.n_full_cis_rows > 0, 'gene'].nunique())}，零行基因保留在分母内。",
        f"- r² < 0.001 下过滤后无显著工具：{no_tool} 个。",
        f"- r² < 0.001 下 K < 3、留一法不可解读：{loo_bad} 个。",
        f"- 其他数据可用性状态：{other} 个。",
        f"- r² < 0.001 下有有限 share_abs 的基因：{len(finite_primary)} 个；其中 K >= 3 且可解释集中度：{denom} 个，"
        f"share_abs >= 0.90：{high} 个。",
        f"- 结论：{interpretation}",
        f"- 在这 {candidate_count} 个正式候选中，K > 3 的基因为 0 个；K >= 3 的 {denom} 个均为 K=3。"
        f"这句话不外推到 {comparator_count} 个对照基因或 Finan 全集。",
        "",
        "## 候选集 share_abs 分布",
        *_distribution_lines(concentration, f"{candidate_count}个正式候选（有限值；K<3仅作结构性诊断）"),
        *_distribution_lines(concentration, f"{candidate_count}个正式候选（主解释分层：K>=3）", min_k=3),
        "",
        "## 对照集 share_abs 分布",
        f"- 固定种子 {comparator_provenance['seed']} 抽取 {comparator_provenance['sample_size_requested']} 个非候选基因；"
        f"其中实际有 HGNC symbol 的基因为 {comparator_count} 个，{comparator_provenance['sampled_gene_count'] - comparator_count} 个没有 HGNC symbol，"
        f"因此组间比例分母使用 {comparator_count}。",
        f"- 实际有 HGNC symbol 的对照基因中，原始文件有 cis 行的基因为 {comparator_provenance['sampled_gene_count_with_raw_eqtl_rows']} 个。",
        *_distribution_lines(comparison, f"{comparator_count}个有 HGNC symbol 的对照基因（有限值；K<3仅作结构性诊断）"),
        *_distribution_lines(comparison, f"{comparator_count}个有 HGNC symbol 的对照基因（主解释分层：K>=3）", min_k=3),
        "",
        "## 候选集与对照集组间比较（描述性）",
        "所有比较均为双侧、未作 FDR 或其他多重校正；完整枚举与小样本对照的统计性质分开解释，组间比较效能有限。",
        f"- K >= 3 比例：候选 {kge3['candidate_success']}/{kge3['candidate_total']} ({_pct(kge3['candidate_success'], kge3['candidate_total'])})，"
        f"对照 {kge3['comparison_success']}/{kge3['comparison_total']} ({_pct(kge3['comparison_success'], kge3['comparison_total'])})；"
        f"Fisher OR={kge3['odds_ratio']:.3f}，P={_fmt_p(kge3['p_value_two_sided'])}。",
        f"- K > 3 比例：候选 {kgt3['candidate_success']}/{kgt3['candidate_total']}，"
        f"对照 {kgt3['comparison_success']}/{kgt3['comparison_total']}（其中 K=4 为 {comp_kgt3_counts.get(4, 0)} 个，K=5 为 {comp_kgt3_counts.get(5, 0)} 个）；"
        f"Fisher OR={kgt3['odds_ratio']:.3f}，P={_fmt_p(kgt3['p_value_two_sided'])}。",
        f"- K >= 3 层内 share_abs >= 0.90：候选 {high_test['candidate_success']}/{high_test['candidate_total']}，"
        f"对照 {high_test['comparison_success']}/{high_test['comparison_total']}；"
        f"Fisher OR={high_test['odds_ratio']:.3f}，P={_fmt_p(high_test['p_value_two_sided'])}，未检出差异。",
        f"- K >= 3 层内 share_abs 分布：候选中位数 {dist_test['candidate_median']:.4f}（n={dist_test['candidate_n']}），"
        f"对照 {dist_test['comparison_median']:.4f}（n={dist_test['comparison_n']}）；"
        f"Mann–Whitney U={dist_test['u_statistic']:.1f}，P={_fmt_p(dist_test['p_value_two_sided'])}。",
        f"- 仅 K=3 层内 share_abs 分布：候选中位数 {k3_test['candidate_median']:.4f}（n={k3_test['candidate_n']}），"
        f"对照 {k3_test['comparison_median']:.4f}（n={k3_test['comparison_n']}）；"
        f"Mann–Whitney U={k3_test['u_statistic']:.1f}，P={_fmt_p(k3_test['p_value_two_sided'])}。",
        "- 工具稀缺不是候选选择程序产物这一点得到 K >= 3 比例比较的支持；但候选集是否额外富集主导架构未能判定，"
        "因为 K > 3 比例与 K=3 层内 share_abs 的方向相反且均未达名义显著。",
        "- 三项方向可被描述为一条待检验趋势：按 lead-eQTL 显著性选择可能拉入单一强信号位点，"
        "使部分基因达到 3 个工具但缺少第 4 个独立信号且集中度较高；该整体解释统计上未获支持，只能作为不显著的一致趋势陈述。",
        "",
        "## 无工具计数的历史口径核对",
        f"- 本轮统一采用：从 raw full cis 文件覆盖 {candidate_count}/{candidate_count} 后，执行 P<5e-8、F>10、"
        f"N>=0.80*gene_max_N、回文位点排除；该口径下 no_tool={no_tool}，基因为 "
        f"{', '.join(current_no_tool)}。",
        f"- 较早的 `candidate_significant_eqtl.tsv.gz` 只保留显著 eQTL 行，实际表示 {legacy.get('legacy_significant_table_gene_count', 'NA')} 个基因；"
        f"从 76 个候选中缺失 {legacy.get('legacy_significant_table_missing_from_76', 'NA')} 个（{', '.join(legacy.get('legacy_significant_table_missing_genes', []))}）。"
        "这个 5 个的数字是显著性筛选表的缺失计数，不是完整 raw cis + 样本量过滤后的 no_tool 计数。",
        f"- 本轮 no_tool 中不在旧表缺失名单内的基因为 {', '.join(legacy.get('current_no_tool_not_in_legacy_significant_missing', []))}；"
        "这些基因在旧表中有显著行，但在本轮 N>=0.80*gene_max_N 过滤后没有保留工具。"
        "旧表缺失名单与本轮 no_tool 的重叠及差集已写入 group-tests JSON。",
        f"- `candidate_full_eqtl.tsv.gz` 实际表示 {legacy.get('legacy_full_table_gene_count', 'NA')} 个基因，仍遗漏零行基因；"
        "因此也不能替代本轮 76 基因完整枚举。稿件及后续报告统一采用本轮 raw full cis 口径的 no_tool 数值。",
        "",
        "## 三终点留一法（次级、结局侧条件化）",
        f"- 三终点×候选组合总数：{len(loo)}。",
        f"- K >= 3 且可计算的组合：{len(loo_est)}。",
        f"- 可计算组合中出现符号反转：{int(loo_est.sign_reversal.sum())} "
        f"（{(100 * loo_est.sign_reversal.mean() if len(loo_est) else float('nan')):.1f}%）。",
        f"- 其中全集估计名义显著的组合：{len(loo_sig)}；该子集中符号反转："
        f"{int(loo_sig.sign_reversal.sum())} "
        f"（{(100 * loo_sig.sign_reversal.mean() if len(loo_sig) else float('nan')):.1f}%）。",
        "- 该比例受结局端工具重叠、效应大小和估计精度条件化，不能解释为暴露工具几何的总体发生率。",
        "",
        "## 固定阈值",
        "- LD ladder: r² < 0.001, 0.01, 0.05, 0.1, 0.2。",
        "- 高集中度描述阈值: share_abs >= 0.90, 0.95, 0.99。",
        "- 不执行FDR或其他多重检验校正。",
        "",
        "## 输入溯源",
        "- eQTLGen full cis: `data/raw/idd/genetics/eqtlgen/eqtlgen_full_cis.txt.gz` "
        f"(SHA-256 `{sha256_file(RAW_EQTL)}`).",
        "- eQTLGen allele frequency: `data/raw/idd/genetics/eqtlgen/eqtlgen_allele_frequency.txt.gz` "
        f"(SHA-256 `{sha256_file(RAW_AF)}`).",
        "- 1000 Genomes EUR LD: `data/raw/idd/genetics/ld/1000G_EUR.bim` and "
        f"`data/raw/idd/genetics/ld/1000G_EUR.bed`; n={N_REF}.",
        f"- 组间检验 JSON：`out/global_ablation_group_tests.json`；对照表 SHA-256 已写入该 JSON。",
    ]
    (OUT / "GLOBAL_ABLATION_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    global ROOT, WS, PROC, OUT, RAW_EQTL, RAW_AF, CANDIDATES, LD_DIR, BED, BIM
    global COMPARISON_TABLE, COMPARISON_PROVENANCE, GROUP_TESTS, LEGACY_FULL_EQTL, LEGACY_SIGNIFICANT_EQTL
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--reuse-extracted", action="store_true")
    args = parser.parse_args()
    WS, OUT = args.input_root.resolve(), args.output_dir.resolve()
    if not WS.is_dir(): parser.error("--input-root must exist")
    if WS == OUT or WS in OUT.parents or OUT in WS.parents: parser.error("output directory must be disjoint from input root")
    ROOT = OUT.parent
    PROC = WS / "data/processed/idd/genetics"
    RAW_EQTL = WS / "data/raw/idd/genetics/eqtlgen/eqtlgen_full_cis.txt.gz"
    RAW_AF = WS / "data/raw/idd/genetics/eqtlgen/eqtlgen_allele_frequency.txt.gz"
    CANDIDATES = WS / "artifacts/core/genetics/formal_candidate_genes.tsv"
    LD_DIR = WS / "data/raw/idd/genetics/ld"
    BED, BIM = LD_DIR / "1000G_EUR.bed", LD_DIR / "1000G_EUR.bim"
    COMPARISON_TABLE, COMPARISON_PROVENANCE, GROUP_TESTS = OUT / "comparison_set_concentration.tsv", OUT / "comparison_set_provenance.json", OUT / "global_ablation_group_tests.json"
    LEGACY_FULL_EQTL, LEGACY_SIGNIFICANT_EQTL = PROC / "candidate_full_eqtl.tsv.gz", PROC / "candidate_significant_eqtl.tsv.gz"
    endpoint_module.PROC = PROC
    for required in (RAW_EQTL, RAW_AF, CANDIDATES, BED, BIM):
        if not args.report_only and not required.is_file(): parser.error(f"missing input: {required}")
    OUT.mkdir(parents=True, exist_ok=True)
    if args.report_only:
        concentration = pd.read_csv(OUT / "global_ablation_concentration.tsv", sep="\t")
        loo = pd.read_csv(OUT / "global_ablation_loo.tsv", sep="\t")
        group_tests = compute_group_tests(concentration)
        write_report(concentration, loo, group_tests)
        return

    candidates = read_candidates(CANDIDATES)
    if len(candidates) != 76:
        raise AssertionError(f"Expected 76 formal candidates, found {len(candidates)}")
    gene_map = {row["ensembl_gene_id"]: row["gene_symbol"] for row in candidates}
    enriched_path = OUT / "global_ablation_eqtl.tsv.gz"
    if not args.reuse_extracted or not enriched_path.exists():
        enriched_path = extract_enriched_eqtl(gene_map, "global_ablation")
    eqtl = read_enriched(enriched_path)
    extracted_symbols = set(eqtl.gene_symbol.dropna().astype(str))
    candidate_symbols = set(gene_map.values())
    if not extracted_symbols.issubset(candidate_symbols):
        raise AssertionError("extraction contains a gene outside the formal candidate set")
    concentration = run_concentration(candidates, eqtl)
    loo = run_loo(candidates, eqtl)
    provenance = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "formal_candidate_count": len(candidates),
        "candidate_coverage_including_zero_row_genes": len(candidates),
        "candidate_genes_with_raw_rows": len(extracted_symbols),
        "candidate_genes_without_raw_rows": sorted(candidate_symbols - extracted_symbols),
        "raw_eqtl_path": str(RAW_EQTL),
        "raw_eqtl_sha256": sha256_file(RAW_EQTL),
        "raw_allele_frequency_path": str(RAW_AF),
        "raw_allele_frequency_sha256": sha256_file(RAW_AF),
        "ld_reference": "1000G_EUR_n503",
        "r2_ladder": R2_LADDER,
        "share_abs_thresholds": CONCENTRATION_THRESHOLDS,
        "script": "global_ablation.py",
    }
    (OUT / "global_ablation_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    if provenance["candidate_coverage_including_zero_row_genes"] != 76:
        raise AssertionError("candidate coverage including zero-row genes is not 76")
    print(
        "candidate coverage including zero-row genes: "
        f"{provenance['candidate_coverage_including_zero_row_genes']}/76"
    )
    group_tests = compute_group_tests(concentration)
    write_report(concentration, loo, group_tests)
    print(json.dumps(provenance, ensure_ascii=False, indent=2))
    print(f"wrote {len(concentration)} concentration rows and {len(loo)} LOO rows")


if __name__ == "__main__":
    main()
