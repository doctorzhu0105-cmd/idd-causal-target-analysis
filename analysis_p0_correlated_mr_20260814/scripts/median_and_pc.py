"""Weighted median and principal-component IVW for correlated cis instruments.

Two different problems, two different tools:

Weighted median (Bowden et al., Genet Epidemiol 2016) assumes *independent*
instruments. At r2 < 0.1 that assumption fails, so it is reported here both on a
near-independent subset (valid, low power) and on the full set (higher power,
assumption violated) with the status recorded on each row. It is never presented
as if the correlated version were assumption-clean.

PC-IVW (Burgess et al., Genet Epidemiol 2017) is the estimator actually designed
for many correlated cis variants: project the instruments onto the leading
principal components of the LD-and-precision-weighted matrix, then run
generalized IVW in that basis. It is the appropriate robustness check on the
main generalized-IVW estimate.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from correlated_mr import _finalise, MREstimate


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cum = np.cumsum(w) - 0.5 * w
    cum /= w.sum()
    return float(np.interp(0.5, cum, v))


def weighted_median(bx, by, se_x, se_y, n_boot: int = 1000, seed: int = 20260814,
                    assumption: str = "independent") -> MREstimate:
    """Weighted median ratio estimate with a parametric bootstrap standard error."""
    bx, by, se_x, se_y = map(np.asarray, (bx, by, se_x, se_y))
    k = len(bx)
    if k < 3:
        return _finalise("weighted_median", np.nan, np.nan, k,
                         note="需要至少 3 个工具")

    ratios = by / bx
    # first-order delta-method variance of each ratio
    weights = 1.0 / ((se_y / np.abs(bx)) ** 2)
    est = _weighted_median(ratios, weights)

    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        bx_b = rng.normal(bx, se_x)
        by_b = rng.normal(by, se_y)
        # guard against a near-zero exposure draw inverting a ratio
        safe = np.abs(bx_b) > 1e-8
        if safe.sum() < 3:
            draws[i] = np.nan
            continue
        draws[i] = _weighted_median(
            by_b[safe] / bx_b[safe],
            1.0 / ((se_y[safe] / np.abs(bx_b[safe])) ** 2),
        )
    se = float(np.nanstd(draws, ddof=1))

    out = _finalise("weighted_median", float(est), se, k)
    out.note = f"LD 假设: {assumption}"
    return out


def pc_ivw(bx, by, se_y, rho, variance_kept: float = 0.99) -> MREstimate:
    """Principal-component IVW for correlated instruments (Burgess et al. 2017).

    Instruments are projected onto the leading eigenvectors of the
    precision-and-LD weighted matrix, which keeps the informative dimensions
    while dropping the near-collinear ones that make the raw LD matrix unstable.
    """
    bx, by, se_y = map(np.asarray, (bx, by, se_y))
    k = len(bx)
    if k < 3:
        return _finalise("pc_ivw", np.nan, np.nan, k, note="需要至少 3 个工具")

    scale = bx / se_y
    phi = np.outer(scale, scale) * rho
    phi = (phi + phi.T) / 2.0

    eigvals, eigvecs = np.linalg.eigh(phi)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]

    positive = eigvals > 0
    eigvals, eigvecs = eigvals[positive], eigvecs[:, positive]
    if eigvals.size == 0:
        return _finalise("pc_ivw", np.nan, np.nan, k, note="无正特征值")

    share = np.cumsum(eigvals) / eigvals.sum()
    n_pc = int(np.searchsorted(share, variance_kept) + 1)
    n_pc = max(1, min(n_pc, eigvals.size))
    lam = eigvecs[:, :n_pc]

    bx_pc = lam.T @ bx
    by_pc = lam.T @ by
    omega_pc = lam.T @ (np.outer(se_y, se_y) * rho) @ lam
    omega_inv = np.linalg.pinv(omega_pc)

    denom = float(bx_pc @ omega_inv @ bx_pc)
    if denom <= 0:
        return _finalise("pc_ivw", np.nan, np.nan, k, note="非正定")
    beta = float(bx_pc @ omega_inv @ by_pc) / denom
    se_fixed = float(np.sqrt(1.0 / denom))

    resid = by_pc - beta * bx_pc
    q = float(resid @ omega_inv @ resid)
    df = max(n_pc - 1, 1)
    phi_over = max(1.0, q / df)
    se = se_fixed * np.sqrt(phi_over) if n_pc > 1 else se_fixed

    out = _finalise("pc_ivw", beta, se, k,
                    q_stat=q, q_pvalue=float(stats.chi2.sf(q, df)),
                    overdispersion=phi_over)
    out.note = f"保留 {n_pc}/{k} 个主成分（{variance_kept:.0%} 方差）"
    return out
