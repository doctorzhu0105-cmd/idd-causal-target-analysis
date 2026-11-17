"""Selected P1/P2/negative-control computations extracted from frozen R3-9 source.

The original six numerical function bodies are retained byte-for-byte at source
text level. Only I/O, selection plumbing and output presentation are adapted.
This is not end-to-end GWAS processing or verification of all S17-S19 cells.
"""
import math
import numpy as np
import pandas as pd
from scipy.stats import norm

def rank_spearman(x, y):
    from scipy.stats import rankdata
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return float("nan")
    rx, ry = rankdata(x[m], method="average"), rankdata(y[m], method="average")
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])

def blocks_for(q, size_mb=5):
    chrom = q["chr"].astype(str).str.replace(r"^chr", "", regex=True)
    pos = pd.to_numeric(q["pos_lead"], errors="coerce")
    return chrom + ":" + (np.floor(pos / (size_mb * 1_000_000)).astype("Int64").astype(str))

def jackknife(metric_fn, q, size_mb=5):
    labels = blocks_for(q, size_mb).to_numpy()
    unique = pd.unique(labels[pd.notna(labels)])
    full = float(metric_fn(q))
    vals = []
    for b in unique:
        keep = labels != b
        if keep.sum() < 3:
            continue
        v = metric_fn(q.loc[keep])
        if np.isfinite(v):
            vals.append(float(v))
    if len(vals) < 2:
        return full, float("nan"), len(unique), len(vals)
    vals = np.asarray(vals)
    se = math.sqrt((len(vals) - 1) / len(vals) * np.sum((vals - vals.mean()) ** 2))
    return full, float(se), len(unique), len(vals)

def ci(est, se, confidence):
    if not np.isfinite(est) or not np.isfinite(se):
        return float("nan"), float("nan")
    z = norm.ppf(0.5 + confidence / 2.0)
    return float(est - z * se), float(est + z * se)

def p1_fn(q): return rank_spearman(q["delta_theta_fg"], q["delta_theta_uk"])

def p2_fn(q): return rank_spearman(q["dlogphi_fg"], q["dlogphi_uk"])

def compute(source):
    df = pd.read_csv(source, sep="\t", low_memory=False)
    required = ["finngen_endpoint", "external_outcome", "contrast", "chr", "pos_lead",
                "delta_theta_fg", "delta_theta_uk", "dlogphi_fg", "dlogphi_uk",
                "K_lo_fg", "K_hi_fg", "K_lo_uk", "K_hi_uk"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError("Missing required source fields: " + ", ".join(missing))
    sel = lambda endpoint, outcome: df[(df.finngen_endpoint == endpoint) &
                                        (df.external_outcome == outcome) &
                                        (df.contrast == "r2 0.001->0.1 (lo->hi)")].copy()
    p1 = sel("M13_INTERVERTEB", "PanUKB_EUR_icd10_M51")
    neg = sel("M13_INTERVERTEB", "PanUKB_EUR_phecode_366_NEGATIVE_CONTROL")
    if len(p1) != 2427 or len(neg) != 2427:
        raise ValueError(f"Unexpected primary/negative n: {len(p1)}, {len(neg)}")
    p1["k2_all"] = ((p1.K_lo_fg >= 2) & (p1.K_hi_fg >= 2) &
                    (p1.K_lo_uk >= 2) & (p1.K_hi_uk >= 2))
    p2 = p1[p1.k2_all].copy()
    result = {}
    for key, frame, fn in [("A-1", p1, p1_fn), ("A-2", p2, p2_fn), ("C-1", neg, p1_fn)]:
        estimate, se, n_blocks, used = jackknife(fn, frame)
        low, high = ci(estimate, se, 0.975)
        result[key] = {"n": len(frame), "estimate": estimate, "jackknife_se": se,
                       "ci_low": low, "ci_high": high, "blocks": n_blocks,
                       "jackknife_blocks_used": used, "confidence": 0.975}
    return result
