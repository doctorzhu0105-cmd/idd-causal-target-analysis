#!/usr/bin/env python3
"""Recompute the five-threshold catalogue from explicitly supplied provider inputs.

This public entry point contains no manuscript, package, submission, or audit-tree
dependencies. The required eQTLGen, FinnGen, and 1000 Genomes inputs are excluded
from the repository and must be obtained from their providers under their terms.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


R2_LADDER = (0.001, 0.01, 0.05, 0.1, 0.2)
OUTCOMES = ("M13_INTERVERTEB", "M13_LOWBACKPAIN", "M13_SCIATICA")
PINV_RCOND = 1e-15


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def effective_k(rho: np.ndarray) -> float:
    eigenvalues = np.clip(np.linalg.eigvalsh(rho), 0.0, None)
    denominator = float(np.square(eigenvalues).sum())
    return float(np.square(eigenvalues.sum()) / denominator) if denominator > 0 else float("nan")


def matrix_metrics(matrix: np.ndarray) -> dict[str, object]:
    if matrix.size == 0:
        return {
            "numerical_rank": 0,
            "minimum_eigenvalue": np.nan,
            "condition_number": np.nan,
            "pinv_singular_direction_discard_count": 0,
        }
    eigenvalues = np.linalg.eigvalsh(matrix)
    singular = np.linalg.svd(matrix, compute_uv=False)
    cutoff = PINV_RCOND * float(singular.max()) if len(singular) else 0.0
    return {
        "numerical_rank": int(np.linalg.matrix_rank(matrix)),
        "minimum_eigenvalue": float(eigenvalues.min()),
        "condition_number": float(np.linalg.cond(matrix)),
        "pinv_singular_direction_discard_count": int((singular <= cutoff).sum()),
    }


def mean_pairwise_r2(rho: np.ndarray) -> float:
    if len(rho) < 2:
        return float("nan")
    upper = rho[np.triu_indices(len(rho), 1)]
    return float(np.mean(np.square(upper)))


def joint_information(beta_x: np.ndarray, se_x: np.ndarray, rho: np.ndarray) -> float:
    sigma_x = np.outer(se_x, se_x) * rho
    return float(beta_x @ np.linalg.pinv(sigma_x) @ beta_x)


def hat_diagonal(x: np.ndarray, covariance: np.ndarray) -> np.ndarray:
    design = np.asarray(x, dtype=float).reshape(-1, 1)
    precision = np.linalg.pinv(covariance)
    scalar = float((design.T @ precision @ design)[0, 0])
    if not np.isfinite(scalar) or scalar <= 0:
        return np.full(len(design), np.nan)
    return np.diag(design @ np.array([[1.0 / scalar]]) @ design.T @ precision)


def failure_rows(gene: str, gene_id: str, endpoint: str | None, status: str, reason: str = "") -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for threshold in R2_LADDER:
        row: dict[str, object] = {
            "gene": gene,
            "ensembl_gene_id": gene_id,
            "threshold": threshold,
            "K": 0,
            "effective_K": np.nan,
            "joint_exposure_chi_square": np.nan,
            "conditional_strength_analogue": np.nan,
            "status": status,
            "reason": reason,
        }
        if endpoint is not None:
            row.update({"endpoint": endpoint, "beta": np.nan, "SE": np.nan, "P": np.nan})
        rows.append(row)
    return rows


def exposure_rows(gene_df: pd.DataFrame, gene: str, gene_id: str, bim: pd.DataFrame, bed: Path, n_ref: int,
                  ld_matrix, deduplicate_after_ld, prune_by_ld) -> list[dict[str, object]]:
    usable = gene_df.dropna(subset=["beta", "se", "f_statistic"]).copy()
    if usable.empty:
        return failure_rows(gene, gene_id, None, "no_exposure_conversion")
    try:
        chrom = str(int(usable.chr_b37.iloc[0]))
        variants, rho = ld_matrix(usable, bim[bim.chr == chrom], str(bed), n_ref)
        if variants.empty:
            return failure_rows(gene, gene_id, None, "no_ld_reference")
        variants, rho = deduplicate_after_ld(variants, rho)
        order = np.argsort(-variants.f_statistic.to_numpy())
        rows: list[dict[str, object]] = []
        for threshold in R2_LADDER:
            index = prune_by_ld(rho, order, threshold)
            if len(index) == 0:
                rows.append(failure_rows(gene, gene_id, None, "no_pruned_instruments")[0] | {"threshold": threshold})
                continue
            subset = variants.iloc[index].reset_index(drop=True)
            sub_rho = rho[np.ix_(index, index)]
            beta_x = subset.beta.to_numpy(float)
            se_x = subset.se.to_numpy(float)
            k = len(subset)
            info = joint_information(beta_x, se_x, sub_rho)
            metrics = matrix_metrics(sub_rho)
            covariance = np.outer(se_x, se_x) * sub_rho
            leverage = hat_diagonal(beta_x, covariance)
            rows.append({
                "gene": gene, "ensembl_gene_id": gene_id, "threshold": threshold, "K": k,
                "effective_K": effective_k(sub_rho), "effective_K_over_K": effective_k(sub_rho) / k,
                "K_minus_effective_K": k - effective_k(sub_rho), "mean_pairwise_r2": mean_pairwise_r2(sub_rho),
                "LD_condition_number": metrics["condition_number"], "numerical_rank": metrics["numerical_rank"],
                "minimum_eigenvalue": metrics["minimum_eigenvalue"],
                "pinv_singular_direction_discard_count": metrics["pinv_singular_direction_discard_count"],
                "joint_exposure_chi_square": info, "conditional_strength_analogue": info / k,
                "max_generalized_leverage": float(np.nanmax(leverage)), "leverage_sum": float(np.nansum(leverage)),
                "status": "estimable", "reason": "",
            })
        return rows
    except Exception as exc:
        return failure_rows(gene, gene_id, None, f"exception_{type(exc).__name__}", str(exc))


def fit_rows(eqtl: pd.DataFrame, gene: str, gene_id: str, endpoint: str, outcome: pd.DataFrame,
             bim: pd.DataFrame, bed: Path, n_ref: int, ld_matrix, harmonise, prune_by_ld, ivw_correlated) -> list[dict[str, object]]:
    try:
        harmonized = harmonise(eqtl, gene, outcome).dropna(subset=["beta_exposure", "se_outcome", "f_statistic"])
        if harmonized.empty:
            return failure_rows(gene, gene_id, endpoint, "no_harmonized_instruments")
        chrom = str(int(harmonized.chr_b37.iloc[0]))
        kept, rho = ld_matrix(harmonized, bim[bim.chr == chrom], str(bed), n_ref)
        if kept.empty:
            return failure_rows(gene, gene_id, endpoint, "no_ld_reference")
        kept = kept.reset_index(drop=True)
        order = np.argsort(-kept.f_statistic.to_numpy())
        rows: list[dict[str, object]] = []
        for threshold in R2_LADDER:
            index = prune_by_ld(rho, order, threshold)
            if len(index) == 0:
                rows.append(failure_rows(gene, gene_id, endpoint, "no_pruned_instruments")[0] | {"threshold": threshold})
                continue
            subset = kept.iloc[index].reset_index(drop=True)
            sub_rho = rho[np.ix_(index, index)]
            beta_x = subset.beta_exposure.to_numpy(float)
            beta_y = subset.beta_outcome.to_numpy(float)
            se_y = subset.se_outcome.to_numpy(float)
            se_x = subset.se.to_numpy(float)
            estimate = ivw_correlated(beta_x, beta_y, se_y, sub_rho, random_effects=True)
            k = len(subset)
            metrics = matrix_metrics(sub_rho)
            omega = np.outer(se_y, se_y) * sub_rho
            info = joint_information(beta_x, se_x, sub_rho)
            leverage = hat_diagonal(beta_x, omega)
            rows.append({
                "gene": gene, "ensembl_gene_id": gene_id, "endpoint": endpoint, "threshold": threshold,
                "K": k, "n_harmonized": len(harmonized), "beta": estimate.beta, "SE": estimate.se,
                "OR": estimate.or_, "CI95_low": estimate.or_low, "CI95_high": estimate.or_high,
                "P": estimate.pvalue, "Q": estimate.q_stat, "Q_P": estimate.q_pvalue,
                "random_effects_inflation_factor": estimate.overdispersion,
                "numerical_rank": metrics["numerical_rank"], "condition_number": metrics["condition_number"],
                "omega_condition_number": float(np.linalg.cond(omega)),
                "minimum_eigenvalue": metrics["minimum_eigenvalue"],
                "pinv_singular_direction_discard_count": metrics["pinv_singular_direction_discard_count"],
                "effective_K": effective_k(sub_rho), "joint_exposure_chi_square": info,
                "conditional_strength_analogue": info / k,
                "max_outcome_weighted_leverage": float(np.nanmax(leverage)),
                "leverage_sum": float(np.nansum(leverage)), "status": "estimable", "reason": "",
            })
        return rows
    except Exception as exc:
        return failure_rows(gene, gene_id, endpoint, f"exception_{type(exc).__name__}", str(exc))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--feasibility-json", type=Path, required=True)
    parser.add_argument("--eqtl-batch-dir", type=Path, required=True)
    parser.add_argument("--outcome-dir", type=Path, required=True)
    parser.add_argument("--bim", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--start-batch", type=int, default=0)
    parser.add_argument("--stop-batch", type=int)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--n-ref", type=int, default=503)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo = args.repo_root.resolve()
    script_dir = repo / "analysis_p0_correlated_mr_20260814" / "scripts"
    sys.path.insert(0, str(script_dir))
    from correlated_mr import ivw_correlated, prune_by_ld
    from global_ablation import _deduplicate_after_ld
    from plink_ld import ld_matrix, load_bim
    from run_all_genes_endpoints import harmonise

    payload = json.loads(args.feasibility_json.read_text(encoding="utf-8"))
    gene_ids = sorted(payload["finan_universe"]["eligible_gene_ids"])
    total_batches = math.ceil(len(gene_ids) / args.batch_size)
    stop = total_batches if args.stop_batch is None else min(args.stop_batch, total_batches)
    bim = load_bim(str(args.bim))
    outcomes = {name: read_tsv(args.outcome_dir / f"full_universe_{name}.tsv.gz") for name in OUTCOMES}
    args.work_dir.mkdir(parents=True, exist_ok=True)
    for batch in range(args.start_batch, stop):
        batch_ids = gene_ids[batch * args.batch_size:(batch + 1) * args.batch_size]
        eqtl = read_tsv(args.eqtl_batch_dir / f"eqtl_batch_{batch:03d}.tsv.gz")
        symbols = eqtl.groupby("ensembl_gene_id")["gene_symbol"].first().to_dict()
        exposure, fits = [], []
        for gene_id in batch_ids:
            gene = str(symbols.get(gene_id, gene_id))
            gene_frame = eqtl[eqtl.ensembl_gene_id == gene_id].copy()
            exposure.extend(exposure_rows(gene_frame, gene, gene_id, bim, args.bed, args.n_ref,
                                          ld_matrix, _deduplicate_after_ld, prune_by_ld))
            for endpoint, outcome in outcomes.items():
                fits.extend(fit_rows(eqtl, gene, gene_id, endpoint, outcome, bim, args.bed, args.n_ref,
                                     ld_matrix, harmonise, prune_by_ld, ivw_correlated))
        write_tsv(pd.DataFrame(exposure), args.work_dir / f"exposure_batch_{batch:03d}.tsv")
        write_tsv(pd.DataFrame(fits), args.work_dir / f"givw_batch_{batch:03d}.tsv")
        print(f"batch {batch + 1}/{total_batches}: exposure={len(exposure)}, fits={len(fits)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
