"""Reusable MASTER 1.4 correlated-cis-MR sensitivity primitives.

The functions are public-facing extractions of the governed Stage 3 Version-B
implementation. They do not select biological targets and do not change the
formal 7,284-row multiplicity reference. LD bootstrap helpers quantify
same-panel reference-sample uncertainty; they are not independent replication.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def adjust_pvalues(values: pd.Series | np.ndarray, method: str) -> np.ndarray:
    """Return Benjamini-Hochberg or Benjamini-Yekutieli adjusted P values."""
    raw = np.asarray(values, dtype=float)
    if method not in {"BH", "BY"}:
        raise ValueError("method must be BH or BY")
    if np.any(~np.isfinite(raw)) or np.any((raw < 0) | (raw > 1)):
        raise ValueError("values must be finite probabilities")
    order = np.argsort(raw, kind="mergesort")
    ranks = np.arange(1, len(raw) + 1, dtype=float)
    factor = 1.0 if method == "BH" else float(np.sum(1.0 / ranks))
    adjusted_sorted = raw[order] * len(raw) * factor / ranks
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    adjusted = np.empty(len(raw), dtype=float)
    adjusted[order] = np.clip(adjusted_sorted, 0.0, 1.0)
    return adjusted


def multiplicity_sensitivity(
    frame: pd.DataFrame,
    *,
    p_column: str = "pvalue",
    outcome_column: str = "outcome",
    k_column: str = "k",
) -> dict[str, pd.DataFrame]:
    """Build formal global and prespecified sensitivity-family corrections."""
    required = {p_column, outcome_column, k_column}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"missing columns: {sorted(missing)}")
    data = frame.loc[np.isfinite(pd.to_numeric(frame[p_column], errors="coerce"))].copy()
    data["raw_P"] = pd.to_numeric(data[p_column], errors="raise")

    global_frame = data.copy()
    global_frame["BH_q"] = adjust_pvalues(global_frame.raw_P, "BH")
    global_frame["BY_q"] = adjust_pvalues(global_frame.raw_P, "BY")

    endpoint_parts = []
    for endpoint, part in data.groupby(outcome_column, sort=True):
        part = part.copy()
        part["endpoint"] = endpoint
        part["endpoint_BH_q"] = adjust_pvalues(part.raw_P, "BH")
        part["endpoint_BY_q"] = adjust_pvalues(part.raw_P, "BY")
        endpoint_parts.append(part)

    k_parts = []
    for label, minimum in (("all_estimable", 0), ("K>=2", 2), ("K>=3", 3), ("K>=4", 4)):
        part = data.loc[pd.to_numeric(data[k_column], errors="raise") >= minimum].copy()
        part["family_definition"] = label
        part["adjusted_q"] = adjust_pvalues(part.raw_P, "BH")
        part["BY_adjusted_q"] = adjust_pvalues(part.raw_P, "BY")
        k_parts.append(part)

    return {
        "global": global_frame,
        "endpoint": pd.concat(endpoint_parts, ignore_index=True),
        "k_family": pd.concat(k_parts, ignore_index=True),
    }


def effective_k(rho: np.ndarray) -> float:
    """Eigenvalue participation ratio used as an effective-dimension diagnostic."""
    eigenvalues = np.clip(np.linalg.eigvalsh(np.asarray(rho, dtype=float)), 0.0, None)
    denominator = float(np.square(eigenvalues).sum())
    return float(np.square(eigenvalues.sum()) / denominator) if denominator > 0 else float("nan")


def maximum_hat_diagonal(bx: np.ndarray, se_y: np.ndarray, rho: np.ndarray) -> float:
    """Maximum diagonal of the one-predictor generalized-IVW hat matrix."""
    bx = np.asarray(bx, dtype=float).reshape(-1, 1)
    se_y = np.asarray(se_y, dtype=float)
    rho = np.asarray(rho, dtype=float)
    omega = np.outer(se_y, se_y) * rho
    precision = np.linalg.pinv(omega)
    scalar = float((bx.T @ precision @ bx)[0, 0])
    if not np.isfinite(scalar) or scalar <= 0:
        return float("nan")
    hat = bx @ np.array([[1.0 / scalar]]) @ bx.T @ precision
    return float(np.max(np.diag(hat)))


def mean_impute_dosages(dosages: np.ndarray) -> np.ndarray:
    """Mean-impute each SNP once before same-panel individual resampling."""
    values = np.asarray(dosages, dtype=float).copy()
    if values.ndim != 2:
        raise ValueError("dosages must have shape variants x samples")
    for row in values:
        missing = ~np.isfinite(row)
        if missing.any():
            observed = row[~missing]
            if observed.size == 0:
                raise ValueError("a variant has no observed dosages")
            row[missing] = observed.mean()
    return values


def correlation_from_resample(dosages: np.ndarray, sample_indices: np.ndarray) -> np.ndarray:
    """Compute signed SNP correlation after individual-level panel resampling."""
    sampled = np.asarray(dosages, dtype=float)[:, np.asarray(sample_indices, dtype=int)]
    centered = sampled - sampled.mean(axis=1, keepdims=True)
    scale = np.sqrt(np.sum(centered * centered, axis=1, keepdims=True))
    standardized = np.divide(centered, scale, out=np.zeros_like(centered), where=scale > 0)
    rho = standardized @ standardized.T
    rho[~np.isfinite(rho)] = 0.0
    np.fill_diagonal(rho, 1.0)
    return rho


def bootstrap_ld(dosages: np.ndarray, *, replicates: int = 200, seed: int = 20260826) -> list[np.ndarray]:
    """Return fixed-instrument LD matrices from the same reference panel."""
    values = mean_impute_dosages(dosages)
    rng = np.random.default_rng(seed)
    n_samples = values.shape[1]
    return [
        correlation_from_resample(values, rng.integers(0, n_samples, size=n_samples))
        for _ in range(replicates)
    ]


def shrink_ld(rho: np.ndarray, value: float) -> np.ndarray:
    """Apply the registered deterministic shrinkage rho_lambda=(1-lambda)rho+lambda I."""
    if not 0.0 <= value <= 1.0:
        raise ValueError("value must be in [0, 1]")
    rho = np.asarray(rho, dtype=float)
    return (1.0 - value) * rho + value * np.eye(len(rho))
