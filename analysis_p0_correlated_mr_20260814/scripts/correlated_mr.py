"""Generalized (LD-aware) IVW and MR-Egger for cis instruments in linkage disequilibrium.

The frozen analysis clumped at r2 < 0.001, which left three instruments for
AKR1C2 and made every pleiotropy-robust estimator inestimable. Cis drug-target
MR normally keeps correlated instruments and accounts for the LD structure in
the weight matrix instead (Burgess, Dudbridge & Thompson, Stat Med 2016).

Weight matrix: Omega = diag(se_Y) %*% rho %*% diag(se_Y)
    IVW    beta = (bX' Omega^-1 bX)^-1 bX' Omega^-1 bY
    Egger  [intercept, slope] = (X' Omega^-1 X)^-1 X' Omega^-1 bY,  X = [1, bX]
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from scipy import stats


@dataclass
class MREstimate:
    method: str
    n_snp: int
    beta: float
    se: float
    pvalue: float
    ci_low: float
    ci_high: float
    or_: float
    or_low: float
    or_high: float
    q_stat: float = float("nan")
    q_pvalue: float = float("nan")
    overdispersion: float = float("nan")
    condition_number: float = float("nan")
    note: str = ""

    def as_dict(self):
        d = asdict(self)
        d["or"] = d.pop("or_")
        return d


def _finalise(method, beta, se, k, **extra) -> MREstimate:
    z = beta / se if se > 0 else np.nan
    p = 2 * stats.norm.sf(abs(z)) if np.isfinite(z) else np.nan
    lo, hi = beta - 1.96 * se, beta + 1.96 * se
    return MREstimate(
        method=method,
        n_snp=k,
        beta=beta,
        se=se,
        pvalue=p,
        ci_low=lo,
        ci_high=hi,
        or_=float(np.exp(beta)),
        or_low=float(np.exp(lo)),
        or_high=float(np.exp(hi)),
        **extra,
    )


def _omega_inverse(se_y: np.ndarray, rho: np.ndarray):
    omega = np.outer(se_y, se_y) * rho
    cond = float(np.linalg.cond(omega))
    # pinv keeps near-singular LD blocks from producing garbage rather than an error
    return np.linalg.pinv(omega), cond


def orient_to_exposure_increasing(bx, by, rho):
    """Orient every instrument to an exposure-increasing effect allele.

    MR-Egger contains an intercept and is not invariant to arbitrary allele
    coding.  If D = diag(sign(bx)), the standard univariable convention is
    bx* = |bx|, by* = D by, and rho* = D rho D.  Outcome standard errors are
    unchanged by the sign flip.
    """
    bx = np.asarray(bx, dtype=float)
    by = np.asarray(by, dtype=float)
    rho = np.asarray(rho, dtype=float)
    if bx.ndim != 1 or by.ndim != 1 or bx.shape != by.shape:
        raise ValueError("bx and by must be one-dimensional arrays of equal length")
    if rho.shape != (len(bx), len(bx)):
        raise ValueError("rho must be a square matrix matching bx")
    signs = np.where(bx < 0.0, -1.0, 1.0)
    return np.abs(bx), by * signs, rho * np.outer(signs, signs)


def ivw_correlated(bx, by, se_y, rho, random_effects=True) -> MREstimate:
    """Generalized IVW allowing correlated instruments."""
    bx, by, se_y = map(np.asarray, (bx, by, se_y))
    k = len(bx)
    if k == 0:
        return _finalise("ivw_correlated", np.nan, np.nan, 0)

    omega_inv, cond = _omega_inverse(se_y, rho)
    denom = float(bx @ omega_inv @ bx)
    beta = float(bx @ omega_inv @ by) / denom
    se_fixed = float(np.sqrt(1.0 / denom))

    resid = by - beta * bx
    q = float(resid @ omega_inv @ resid)
    df = max(k - 1, 1)
    q_p = float(stats.chi2.sf(q, df)) if k > 1 else np.nan
    phi = max(1.0, q / df)

    se = se_fixed * np.sqrt(phi) if (random_effects and k > 1) else se_fixed
    method = "ivw_correlated_random" if (random_effects and k > 1) else "ivw_correlated_fixed"

    return _finalise(
        method, beta, se, k,
        q_stat=q, q_pvalue=q_p, overdispersion=phi, condition_number=cond,
    )


def egger_correlated(bx, by, se_y, rho) -> tuple[MREstimate, MREstimate]:
    """LD-aware MR-Egger. Returns (slope, intercept) estimates.

    Uses a t distribution on K-2 degrees of freedom, matching the correction the
    project already applied to the independent-instrument Egger analysis.
    """
    bx, by, se_y = map(np.asarray, (bx, by, se_y))
    k = len(bx)
    if k < 3:
        nan = _finalise("egger_correlated_slope", np.nan, np.nan, k, note="需要至少 3 个工具")
        nan_i = _finalise("egger_correlated_intercept", np.nan, np.nan, k, note="需要至少 3 个工具")
        return nan, nan_i

    bx, by, rho = orient_to_exposure_increasing(bx, by, rho)
    omega_inv, cond = _omega_inverse(se_y, rho)
    design = np.column_stack([np.ones(k), bx])
    xtwx = design.T @ omega_inv @ design
    xtwx_inv = np.linalg.pinv(xtwx)
    coef = xtwx_inv @ design.T @ omega_inv @ by

    resid = by - design @ coef
    df = k - 2
    phi = max(1.0, float(resid @ omega_inv @ resid) / df)
    se = np.sqrt(np.diag(xtwx_inv) * phi)

    def t_estimate(name, est, err):
        t = est / err if err > 0 else np.nan
        p = 2 * stats.t.sf(abs(t), df) if np.isfinite(t) else np.nan
        crit = stats.t.ppf(0.975, df)
        lo, hi = est - crit * err, est + crit * err
        return MREstimate(
            method=name, n_snp=k, beta=float(est), se=float(err), pvalue=float(p),
            ci_low=float(lo), ci_high=float(hi),
            or_=float(np.exp(est)), or_low=float(np.exp(lo)), or_high=float(np.exp(hi)),
            overdispersion=phi, condition_number=cond,
        )

    slope = t_estimate("egger_correlated_slope", coef[1], se[1])
    intercept = t_estimate("egger_correlated_intercept", coef[0], se[0])
    return slope, intercept


def leave_one_out(bx, by, se_y, rho) -> list[dict]:
    """Refit the generalized IVW with each instrument dropped in turn."""
    bx, by, se_y = map(np.asarray, (bx, by, se_y))
    k = len(bx)
    if k < 3:
        return []

    out = []
    for i in range(k):
        keep = np.arange(k) != i
        est = ivw_correlated(bx[keep], by[keep], se_y[keep], rho[np.ix_(keep, keep)])
        row = est.as_dict()
        row["dropped_index"] = i
        out.append(row)
    return out


def prune_by_ld(rho: np.ndarray, order: np.ndarray, r2_threshold: float) -> np.ndarray:
    """Greedy LD pruning: walk `order` (best instrument first), keeping a variant
    only when its r2 against every already-kept variant is below the threshold."""
    kept: list[int] = []
    for idx in order:
        if all(rho[idx, j] ** 2 <= r2_threshold for j in kept):
            kept.append(int(idx))
    return np.array(sorted(kept), dtype=int)
