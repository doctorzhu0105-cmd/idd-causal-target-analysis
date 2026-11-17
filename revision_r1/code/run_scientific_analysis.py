"""Run recovered numerical components on caller-supplied, isolated inputs.

The K=1 propagation follows Stage10R1 k1_replication/p_adjust. The attenuation
component follows Batch14 sensitivity_worker.py with frozen membership, F,
SE and signed LD. These entries do not generate the upstream MR summary data.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from project_modules.correlated_mr import ivw_correlated


def p_adjust(values: np.ndarray, method: str) -> np.ndarray:
    """Stage10R1 p_adjust, with explicit validation of its two supported methods."""
    if method not in {"BH", "BY"}:
        raise ValueError("method must be BH or BY")
    p = np.asarray(values, dtype=float)
    if len(p) == 0 or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError("p-values must be a nonempty finite vector in [0,1]")
    m = len(p)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    factor = float(m)
    if method == "BY":
        factor *= float(np.sum(1.0 / np.arange(1, m + 1)))
    raw = ranked * factor / np.arange(1, m + 1)
    adjusted_sorted = np.minimum.accumulate(raw[::-1])[::-1]
    adjusted_sorted = np.clip(adjusted_sorted, 0.0, 1.0)
    result = np.empty(m, dtype=float)
    result[order] = adjusted_sorted
    return result


def propagate(fdr: pd.DataFrame, official: pd.DataFrame) -> pd.DataFrame:
    required_fdr = {"gene", "outcome", "status", "k", "se", "beta", "min_f", "pvalue"}
    if not required_fdr.issubset(fdr.columns):
        raise ValueError(f"FDR input missing fields: {sorted(required_fdr - set(fdr.columns))}")
    official = official.rename(columns={"endpoint": "outcome", "K": "official_K"})
    required_official = {"gene", "outcome", "official_K", "BH_q", "BY_q"}
    if not required_official.issubset(official.columns):
        raise ValueError(f"official-family input missing fields: {sorted(required_official - set(official.columns))}")
    fdr = fdr[fdr.status == "estimable"].copy().reset_index(drop=True)
    if fdr.empty:
        raise ValueError("FDR input has no estimable rows")
    fdr = fdr.merge(official[list(required_official)], on=["gene", "outcome"], validate="one_to_one", how="left", indicator=True)
    if not fdr["_merge"].eq("both").all():
        raise ValueError("official-family input is missing one or more gene-outcome keys")
    fdr = fdr.drop(columns="_merge")
    for col in ["k", "se", "beta", "min_f", "pvalue", "official_K", "BH_q", "BY_q"]:
        fdr[col] = pd.to_numeric(fdr[col], errors="raise")
    if not np.isfinite(fdr[["k", "se", "beta", "min_f", "pvalue", "official_K", "BH_q", "BY_q"]].to_numpy()).all() or (fdr.se <= 0).any() or (fdr.min_f <= 0).any() or not (fdr.k == fdr.official_K).all():
        raise ValueError("invalid estimates, F values, or K disagreement")
    fdr["reproduced_BH_q"] = p_adjust(fdr.pvalue.to_numpy(), "BH")
    fdr["reproduced_BY_q"] = p_adjust(fdr.pvalue.to_numpy(), "BY")
    fdr["second_order_se"] = fdr.se
    fdr["second_order_p"] = fdr.pvalue
    k1 = fdr.k == 1
    fdr.loc[k1, "second_order_se"] = np.sqrt(np.square(fdr.loc[k1, "se"]) + np.square(fdr.loc[k1, "beta"]) / fdr.loc[k1, "min_f"])
    fdr.loc[k1, "second_order_p"] = 2 * stats.norm.sf(np.abs(fdr.loc[k1, "beta"] / fdr.loc[k1, "second_order_se"]))
    fdr["second_order_global_BH_q"] = p_adjust(fdr.second_order_p.to_numpy(), "BH")
    fdr["second_order_global_BY_q"] = p_adjust(fdr.second_order_p.to_numpy(), "BY")
    fdr["original_global_BH_pass"] = fdr.BH_q < .05
    fdr["original_global_BY_pass"] = fdr.BY_q < .05
    fdr["second_order_global_BH_pass"] = fdr.second_order_global_BH_q < .05
    fdr["second_order_global_BY_pass"] = fdr.second_order_global_BY_q < .05
    fdr["bh_membership_changed"] = fdr.original_global_BH_pass != fdr.second_order_global_BH_pass
    fdr["by_membership_changed"] = fdr.original_global_BY_pass != fdr.second_order_global_BY_pass
    fdr["se_percent_change"] = 100 * (fdr.second_order_se / fdr.se - 1)
    return fdr


def attenuation(data: dict, family: str, lam: float) -> dict:
    if family not in {"A", "B"} or lam < 0:
        raise ValueError("candidate must be A or B and lambda nonnegative")
    x = np.asarray(data["beta_x"], dtype=float)
    sx = np.asarray(data["se_x"], dtype=float)
    f = np.asarray(data["original_f"], dtype=float)
    rho = np.asarray(data["signed_ld"], dtype=float)
    n = len(x)
    if n == 0 or sx.shape != (n,) or f.shape != (n,) or rho.shape != (n, n) or not all(np.isfinite(v).all() for v in (x, sx, f, rho)) or (sx <= 0).any() or (f <= 0).any():
        raise ValueError("invalid vector/matrix dimensions or values")
    if not np.allclose(f, (x / sx) ** 2, rtol=1e-10, atol=1e-10):
        raise ValueError("original F must equal the frozen (beta_x/se_x)^2")
    covariance = np.outer(sx, sx) * rho
    if not np.allclose(covariance, covariance.T, rtol=1e-10, atol=1e-12):
        raise ValueError("exposure covariance is not symmetric")
    singular = np.linalg.svd(covariance, compute_uv=False)
    if np.linalg.eigvalsh(covariance).min() < -1e-10 * singular[0] or np.any(singular <= 1e-15 * singular[0]):
        raise ValueError("Batch14 covariance PSD/rank gate failed")
    inv_x = np.linalg.pinv(covariance, rcond=1e-15)
    factor = np.maximum(.5, 1 - lam / f) if family == "A" else f / (f + lam)
    adjusted = x * factor
    result = {"candidate": family, "lambda": lam, "K": n, "factor": factor.tolist(), "information": float(adjusted @ inv_x @ adjusted), "original_frozen_f": f.tolist(), "rcond": 1e-15}
    if "beta_y" in data or "se_y" in data:
        y = np.asarray(data["beta_y"], dtype=float); sy = np.asarray(data["se_y"], dtype=float)
        if y.shape != (n,) or sy.shape != (n,) or not np.isfinite(y).all() or not np.isfinite(sy).all() or (sy <= 0).any():
            raise ValueError("invalid outcome vectors")
        omega = np.outer(sy, sy) * rho
        singular_y = np.linalg.svd(omega, compute_uv=False)
        if np.linalg.eigvalsh(omega).min() < -1e-10 * singular_y[0] or np.any(singular_y <= 1e-15 * singular_y[0]):
            raise ValueError("Batch14 outcome covariance PSD/rank gate failed")
        w = np.linalg.pinv(omega, rcond=1e-15)
        d = float(adjusted @ w @ adjusted)
        if d <= 0 or not np.isfinite(d): raise ValueError("nonpositive GIVW denominator")
        theta = float(adjusted @ w @ y) / d
        residual = y - theta * adjusted
        q = float(residual @ w @ residual)
        df = n - 1 if n > 1 else None
        phi = max(1.0, q / df) if df is not None else 1.0
        fe = float(np.sqrt(1 / d))
        result.update({"D": d, "theta": theta, "Q": q, "df": df, "phi": phi, "SE_FE": fe, "SE_RE": fe * np.sqrt(phi)})
    return result


def require_isolated_output(output: Path, inputs: list[Path]) -> None:
    output = output.resolve()
    if output.exists():
        raise ValueError("output directory must be new")
    for item in inputs:
        item = item.resolve()
        if not item.is_file(): raise ValueError(f"input file missing: {item.name}")
        if output == item.parent or item.parent in output.parents or output in item.parents:
            raise ValueError("output directory must be disjoint from each input directory")
    code_root = Path(__file__).resolve().parent
    if output == code_root or code_root in output.parents:
        raise ValueError("output cannot be inside the code directory")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="analysis", required=True)
    m = sub.add_parser("k1-propagation"); m.add_argument("--fdr", type=Path, required=True); m.add_argument("--official", type=Path, required=True); m.add_argument("--output-dir", type=Path, required=True)
    a = sub.add_parser("attenuation"); a.add_argument("--input", type=Path, required=True); a.add_argument("--candidate", choices=["A", "B"], required=True); a.add_argument("--lambda", dest="lam", type=float, required=True); a.add_argument("--output-dir", type=Path, required=True)
    i = sub.add_parser("givw"); i.add_argument("--input", type=Path, required=True); i.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        inputs = [args.fdr, args.official] if args.analysis == "k1-propagation" else [args.input]
        require_isolated_output(args.output_dir, inputs)
        if args.analysis == "k1-propagation":
            fdr = pd.read_csv(args.fdr, sep="\t", comment="#")
            official = pd.read_csv(args.official, sep="\t", comment="#")
            result = propagate(fdr, official)
            args.output_dir.mkdir(parents=True)
            result.to_csv(args.output_dir / "k1_second_order_multiplicity.tsv", sep="\t", index=False)
            summary = {"rows": len(result), "k1_rows": int((result.k == 1).sum()), "bh_changes": int(result.bh_membership_changed.sum()), "by_changes": int(result.by_membership_changed.sum()), "family": "all supplied estimable rows"}
        else:
            data = json.loads(args.input.read_text(encoding="utf-8-sig"))
            if args.analysis == "attenuation":
                summary = attenuation(data, args.candidate, args.lam)
            else:
                result = ivw_correlated(np.asarray(data["beta_x"]), np.asarray(data["beta_y"]), np.asarray(data["se_y"]), np.asarray(data["signed_ld"]))
                summary = {key: (None if isinstance(value, float) and not np.isfinite(value) else value) for key, value in result.as_dict().items()}
            args.output_dir.mkdir(parents=True)
        (args.output_dir / "result.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(json.dumps(summary, ensure_ascii=False)); return 0
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"INPUT_ERROR: {exc}", file=sys.stderr); return 2


if __name__ == "__main__":
    raise SystemExit(main())
