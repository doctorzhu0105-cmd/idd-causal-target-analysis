"""Shared statistical and harmonization helpers for the IDD genetics pipeline."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Iterable, Mapping

from scipy.stats import t as student_t


COMPLEMENTS = {"A": "T", "T": "A", "C": "G", "G": "C"}


@dataclass(frozen=True)
class HarmonizedEffect:
    beta_exposure: float
    se_exposure: float
    beta_outcome: float
    se_outcome: float
    effect_allele: str
    other_allele: str
    action: str


def is_palindromic(a1: str, a2: str) -> bool:
    alleles = {a1.upper(), a2.upper()}
    return alleles in ({"A", "T"}, {"C", "G"})


def harmonize_effects(
    exposure: Mapping[str, object],
    outcome: Mapping[str, object],
    *,
    palindromic_maf_max: float = 0.42,
) -> tuple[HarmonizedEffect | None, str]:
    """Align outcome effects to the exposure assessed allele.

    Strand complements are accepted for non-palindromic SNPs. Palindromic SNPs
    are retained only when both data sets provide unambiguous allele frequency.
    """
    ea_exp = str(exposure["effect_allele"]).upper()
    oa_exp = str(exposure["other_allele"]).upper()
    ea_out = str(outcome["effect_allele"]).upper()
    oa_out = str(outcome["other_allele"]).upper()
    valid = {"A", "C", "G", "T"}
    if any(a not in valid for a in (ea_exp, oa_exp, ea_out, oa_out)):
        return None, "non_snv_or_missing_allele"

    beta_out = float(outcome["beta"])
    action = "aligned"
    if (ea_exp, oa_exp) == (ea_out, oa_out):
        pass
    elif (ea_exp, oa_exp) == (oa_out, ea_out):
        beta_out = -beta_out
        action = "swapped"
    elif (COMPLEMENTS[ea_exp], COMPLEMENTS[oa_exp]) == (ea_out, oa_out):
        action = "strand_flip"
    elif (COMPLEMENTS[ea_exp], COMPLEMENTS[oa_exp]) == (oa_out, ea_out):
        beta_out = -beta_out
        action = "strand_flip_and_swap"
    else:
        return None, "allele_mismatch"

    if is_palindromic(ea_exp, oa_exp):
        eaf_exp = _optional_float(exposure.get("eaf"))
        eaf_out = _optional_float(outcome.get("eaf"))
        if eaf_exp is None or eaf_out is None:
            return None, "palindromic_frequency_missing"
        maf_exp = min(eaf_exp, 1.0 - eaf_exp)
        maf_out = min(eaf_out, 1.0 - eaf_out)
        if maf_exp > palindromic_maf_max or maf_out > palindromic_maf_max:
            return None, "palindromic_ambiguous"
        aligned_out_eaf = eaf_out if beta_out == float(outcome["beta"]) else 1.0 - eaf_out
        if abs(eaf_exp - aligned_out_eaf) > 0.1:
            return None, "palindromic_frequency_mismatch"

    return (
        HarmonizedEffect(
            beta_exposure=float(exposure["beta"]),
            se_exposure=float(exposure["se"]),
            beta_outcome=beta_out,
            se_outcome=float(outcome["se"]),
            effect_allele=ea_exp,
            other_allele=oa_exp,
            action=action,
        ),
        "ok",
    )


def standardized_effect_from_z(z: float, eaf: float, n: float) -> tuple[float, float]:
    """Recover a standardized additive effect and SE from Z, EAF, and N.

    This is the summary-data conversion used by SMR for a standardized trait:
    SE = 1 / sqrt(2*p*(1-p)*(N + Z^2)); beta = Z*SE.
    """
    if not (0.0 < eaf < 1.0):
        raise ValueError("eaf must be between zero and one")
    if n <= 0:
        raise ValueError("n must be positive")
    denominator = 2.0 * eaf * (1.0 - eaf) * (n + z * z)
    se = 1.0 / math.sqrt(denominator)
    return z * se, se


def wald_ratio(beta_x: float, se_x: float, beta_y: float, se_y: float) -> dict[str, float]:
    if beta_x == 0:
        raise ValueError("exposure effect must be nonzero")
    beta = beta_y / beta_x
    se = math.sqrt((se_y * se_y) / (beta_x * beta_x) + (beta_y * beta_y * se_x * se_x) / (beta_x**4))
    z = beta / se
    return _effect_summary(beta, se, _two_sided_p(z))


def ivw_fixed(rows: Iterable[Mapping[str, float]]) -> dict[str, float]:
    values = list(rows)
    if not values:
        raise ValueError("at least one instrument is required")
    bx = [float(r["beta_exposure"]) for r in values]
    by = [float(r["beta_outcome"]) for r in values]
    sy = [float(r["se_outcome"]) for r in values]
    weights = [1.0 / (s * s) for s in sy]
    denominator = sum(w * x * x for w, x in zip(weights, bx))
    if denominator <= 0:
        raise ValueError("invalid IVW denominator")
    beta = sum(w * x * y for w, x, y in zip(weights, bx, by)) / denominator
    se = math.sqrt(1.0 / denominator)
    residuals = [y - beta * x for x, y in zip(bx, by)]
    q = sum(w * r * r for w, r in zip(weights, residuals))
    q_df = len(values) - 1
    q_p = chi_square_survival(q, q_df) if q_df > 0 else math.nan
    result = _effect_summary(beta, se, _two_sided_p(beta / se))
    result.update({"q": q, "q_df": q_df, "q_pvalue": q_p})
    return result


def ivw_multiplicative_random_effects(
    rows: Iterable[Mapping[str, float]],
) -> dict[str, float]:
    values = list(rows)
    if len(values) < 2:
        raise ValueError("multiplicative random-effects IVW requires at least two instruments")
    fixed = ivw_fixed(values)
    dispersion = max(1.0, fixed["q"] / fixed["q_df"])
    se = fixed["se"] * math.sqrt(dispersion)
    result = _effect_summary(fixed["beta"], se, _two_sided_p(fixed["beta"] / se))
    result.update(
        {
            "q": fixed["q"],
            "q_df": fixed["q_df"],
            "q_pvalue": fixed["q_pvalue"],
            "dispersion": dispersion,
        }
    )
    return result


def weighted_median(
    rows: Iterable[Mapping[str, float]],
    *,
    seed: int,
    bootstrap_iterations: int = 2000,
) -> dict[str, float]:
    values = list(rows)
    if len(values) < 3:
        raise ValueError("weighted median requires at least three instruments")
    ratios = [float(row["beta_outcome"]) / float(row["beta_exposure"]) for row in values]
    weights = [
        (float(row["beta_exposure"]) ** 2) / (float(row["se_outcome"]) ** 2)
        for row in values
    ]
    beta = _weighted_median_value(ratios, weights)
    rng = random.Random(seed)
    bootstrap = []
    for _ in range(bootstrap_iterations):
        simulated_ratios = []
        simulated_weights = []
        valid = True
        for row in values:
            bx = rng.gauss(float(row["beta_exposure"]), float(row["se_exposure"]))
            by = rng.gauss(float(row["beta_outcome"]), float(row["se_outcome"]))
            if abs(bx) < 1e-12:
                valid = False
                break
            simulated_ratios.append(by / bx)
            simulated_weights.append((bx * bx) / (float(row["se_outcome"]) ** 2))
        if valid:
            bootstrap.append(_weighted_median_value(simulated_ratios, simulated_weights))
    if len(bootstrap) < max(100, bootstrap_iterations // 2):
        raise ValueError("weighted median bootstrap produced too few valid replicates")
    mean = sum(bootstrap) / len(bootstrap)
    se = math.sqrt(sum((value - mean) ** 2 for value in bootstrap) / (len(bootstrap) - 1))
    if se <= 0:
        raise ValueError("weighted median bootstrap standard error is zero")
    return _effect_summary(beta, se, _two_sided_p(beta / se))


def mr_egger(rows: Iterable[Mapping[str, float]]) -> dict[str, float]:
    values = list(rows)
    if len(values) < 3:
        raise ValueError("MR-Egger requires at least three instruments")
    orientations = [1.0 if float(row["beta_exposure"]) >= 0 else -1.0 for row in values]
    x = [abs(float(row["beta_exposure"])) for row in values]
    y = [orientation * float(row["beta_outcome"]) for orientation, row in zip(orientations, values)]
    weights = [1.0 / (float(row["se_outcome"]) ** 2) for row in values]
    sw = sum(weights)
    swx = sum(w * value for w, value in zip(weights, x))
    swy = sum(w * value for w, value in zip(weights, y))
    swxx = sum(w * value * value for w, value in zip(weights, x))
    swxy = sum(w * vx * vy for w, vx, vy in zip(weights, x, y))
    denominator = sw * swxx - swx * swx
    if denominator <= 1e-20:
        raise ValueError("MR-Egger exposure effects have insufficient spread")
    slope = (sw * swxy - swx * swy) / denominator
    intercept = (swy - slope * swx) / sw
    residuals = [vy - intercept - slope * vx for vx, vy in zip(x, y)]
    q = sum(w * residual * residual for w, residual in zip(weights, residuals))
    df = len(values) - 2
    dispersion = max(1.0, q / df)
    slope_se = math.sqrt(dispersion * sw / denominator)
    intercept_se = math.sqrt(dispersion * swxx / denominator)
    critical = float(student_t.ppf(0.975, df))
    slope_pvalue = float(2.0 * student_t.sf(abs(slope / slope_se), df))
    result = {
        "slope_beta": slope,
        "slope_se": slope_se,
        "slope_pvalue": slope_pvalue,
        "slope_ci_lower": slope - critical * slope_se,
        "slope_ci_upper": slope + critical * slope_se,
        "slope_odds_ratio": math.exp(slope),
        "slope_or_ci_lower": math.exp(slope - critical * slope_se),
        "slope_or_ci_upper": math.exp(slope + critical * slope_se),
        "intercept": intercept,
        "intercept_se": intercept_se,
        "intercept_pvalue": float(2.0 * student_t.sf(abs(intercept / intercept_se), df)),
        "q": q,
        "q_df": df,
        "q_pvalue": chi_square_survival(q, df),
        "dispersion": dispersion,
    }
    return result


def observed_r2_from_z(z: float, n: float) -> float:
    if n <= 2:
        raise ValueError("n must exceed two")
    return (z * z) / (z * z + n - 2.0)


def effective_binary_n(cases: float, controls: float) -> float:
    if cases <= 0 or controls <= 0:
        raise ValueError("case and control counts must be positive")
    return 4.0 / (1.0 / cases + 1.0 / controls)


def wakefield_log_abf(beta: float, variance: float, prior_variance: float) -> float:
    if variance <= 0 or prior_variance <= 0:
        raise ValueError("variances must be positive")
    r = prior_variance / (prior_variance + variance)
    z2 = (beta * beta) / variance
    return 0.5 * (math.log1p(-r) + r * z2)


def coloc_posteriors(
    log_abf_trait1: Iterable[float],
    log_abf_trait2: Iterable[float],
    *,
    p1: float = 1e-4,
    p2: float = 1e-4,
    p12: float = 1e-5,
) -> dict[str, float]:
    """Compute coloc.abf hypothesis posteriors from per-variant log ABFs."""
    a = list(log_abf_trait1)
    b = list(log_abf_trait2)
    if len(a) != len(b) or not a:
        raise ValueError("ABF vectors must have equal nonzero length")
    l1 = _logsumexp(a)
    l2 = _logsumexp(b)
    l12_same = _logsumexp([x + y for x, y in zip(a, b)])
    # Sum over distinct causal variants: sum(BF1)*sum(BF2)-sum(BF1*BF2).
    log_product = l1 + l2
    ratio = math.exp(min(0.0, l12_same - log_product))
    log_distinct = log_product + math.log1p(-ratio) if ratio < 1.0 else -math.inf
    log_h3 = math.log(p1) + math.log(p2) + log_distinct
    log_weights = [
        0.0,
        math.log(p1) + l1,
        math.log(p2) + l2,
        log_h3,
        math.log(p12) + l12_same,
    ]
    normalizer = _logsumexp(log_weights)
    return {f"pp{i}": math.exp(value - normalizer) for i, value in enumerate(log_weights)}


def chi_square_survival(value: float, degrees_freedom: int) -> float:
    """Chi-square survival function for integer df via regularized gamma Q."""
    if degrees_freedom <= 0:
        return math.nan
    return _regularized_gamma_q(degrees_freedom / 2.0, value / 2.0)


def _effect_summary(beta: float, se: float, pvalue: float) -> dict[str, float]:
    return {
        "beta": beta,
        "se": se,
        "pvalue": pvalue,
        "ci_lower": beta - 1.959963984540054 * se,
        "ci_upper": beta + 1.959963984540054 * se,
        "odds_ratio": math.exp(beta),
        "or_ci_lower": math.exp(beta - 1.959963984540054 * se),
        "or_ci_upper": math.exp(beta + 1.959963984540054 * se),
    }


def _weighted_median_value(values: list[float], weights: list[float]) -> float:
    if len(values) != len(weights) or not values:
        raise ValueError("weighted median inputs must have equal nonzero length")
    if any(weight < 0 or not math.isfinite(weight) for weight in weights):
        raise ValueError("weighted median weights must be finite and nonnegative")
    total = sum(weights)
    if total <= 0:
        raise ValueError("weighted median weights sum to zero")
    ordered = sorted(zip(values, weights), key=lambda item: item[0])
    cumulative = 0.0
    for value, weight in ordered:
        cumulative += weight
        if cumulative >= total / 2.0:
            return value
    return ordered[-1][0]


def benjamini_hochberg(pvalues: Iterable[float | None]) -> list[float | None]:
    values = list(pvalues)
    indexed = [
        (index, float(value))
        for index, value in enumerate(values)
        if value is not None and math.isfinite(float(value))
    ]
    indexed.sort(key=lambda item: item[1])
    adjusted: list[float | None] = [None] * len(values)
    running = 1.0
    m = len(indexed)
    for rank_index in range(m - 1, -1, -1):
        original_index, pvalue = indexed[rank_index]
        rank = rank_index + 1
        running = min(running, pvalue * m / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def _optional_float(value: object) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _two_sided_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def _logsumexp(values: Iterable[float]) -> float:
    vals = list(values)
    maximum = max(vals)
    if maximum == -math.inf:
        return -math.inf
    return maximum + math.log(sum(math.exp(value - maximum) for value in vals))


def _regularized_gamma_q(a: float, x: float) -> float:
    # Numerical Recipes series/continued-fraction implementation.
    if x < 0 or a <= 0:
        raise ValueError("invalid gamma arguments")
    if x == 0:
        return 1.0
    eps = 3e-14
    tiny = 1e-300
    if x < a + 1.0:
        term = 1.0 / a
        total = term
        ap = a
        for _ in range(10000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * eps:
                break
        p = total * math.exp(-x + a * math.log(x) - math.lgamma(a))
        return max(0.0, min(1.0, 1.0 - p))
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    q = math.exp(-x + a * math.log(x) - math.lgamma(a)) * h
    return max(0.0, min(1.0, q))
