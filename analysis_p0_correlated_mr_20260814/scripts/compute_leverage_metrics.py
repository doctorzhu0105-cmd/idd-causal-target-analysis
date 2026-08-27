"""Compute bounded generalized-IVW hat-matrix leverage diagnostics.

This is a new reporting metric.  It does not overwrite the historical
``share_abs``/``share_raw`` outputs.  For a one-predictor GLS design, the hat
matrix is H = X (X' W X)^-1 X' W and its diagonal sums to one.  The reported
summary is the largest diagonal element, with the sum retained as a numerical
check.  Exposure-side census rows use W = pinv(rho); the AKR1C2 ablation uses
the same outcome-weighted Omega definition as correlated MR.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
WS = ROOT.parent
OUT = ROOT / "out"
EQTL = OUT / "full_universe_eqtl.tsv.gz"
ADEQUACY = OUT / "full_universe_adequacy.tsv"
LD_DIR = WS / "data" / "raw" / "idd" / "genetics" / "ld"
BED = LD_DIR / "1000G_EUR.bed"
BIM = LD_DIR / "1000G_EUR.bim"
N_REF = 503
R2_LADDER = [0.001, 0.01, 0.05, 0.1, 0.2]

sys.path.insert(0, str(SCRIPT_DIR))
from global_ablation import _filtered_gene, _ld_pruned_exposure  # noqa: E402
from plink_ld import load_bim  # noqa: E402
from correlated_mr import prune_by_ld  # noqa: E402
from run_all_genes_endpoints import harmonise, load_outcome  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def hat_diagonal(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float).reshape(-1, 1)
    precision = np.linalg.pinv(w)
    scalar = float((x.T @ precision @ x)[0, 0])
    if not np.isfinite(scalar) or scalar <= 0:
        return np.full(len(x), np.nan)
    hat = x @ np.array([[1.0 / scalar]]) @ x.T @ precision
    return np.diag(hat)


def census() -> pd.DataFrame:
    adequacy = pd.read_csv(ADEQUACY, sep="\t")
    eqtl = pd.read_csv(EQTL, sep="\t", compression="gzip")
    bim_all = load_bim(str(BIM))
    rows: list[dict[str, object]] = []
    for gene, gene_rows in adequacy.groupby("gene", sort=True):
        gene_id = str(gene_rows["ensembl_gene_id"].iloc[0])
        _, _, usable, _ = _filtered_gene(eqtl, gene)
        if usable.empty:
            for r2 in R2_LADDER:
                rows.append({"gene": gene, "ensembl_gene_id": gene_id, "r2_threshold": r2, "k": 0, "max_leverage": np.nan, "leverage_sum": np.nan, "min_leverage": np.nan, "status": "no_filtered_instruments"})
            continue
        chrom = str(int(usable.chr_b37.iloc[0]))
        bim = bim_all[bim_all.chr == chrom]
        for r2 in R2_LADDER:
            sub, rho, _ = _ld_pruned_exposure(usable, bim, r2)
            if sub.empty:
                rows.append({"gene": gene, "ensembl_gene_id": gene_id, "r2_threshold": r2, "k": 0, "max_leverage": np.nan, "leverage_sum": np.nan, "min_leverage": np.nan, "status": "no_ld_reference"})
                continue
            z = (sub.beta / sub.se).to_numpy(float)
            leverage = hat_diagonal(z, rho)
            rows.append({"gene": gene, "ensembl_gene_id": gene_id, "r2_threshold": r2, "k": int(len(sub)), "max_leverage": float(np.nanmax(leverage)), "leverage_sum": float(np.nansum(leverage)), "min_leverage": float(np.nanmin(leverage)), "status": "ok"})
    result = pd.DataFrame(rows)
    if len(result) != 2479 * len(R2_LADDER):
        raise AssertionError(f"Unexpected census rows: {len(result)}")
    result.to_csv(OUT / "full_universe_leverage.tsv", sep="\t", index=False)
    return result


def ablation() -> pd.DataFrame:
    eqtl = pd.read_csv(OUT / "all5_full_cis.tsv", sep="\t")
    bim_all = load_bim(str(BIM))
    harm = harmonise(eqtl, "AKR1C2", load_outcome("M13_INTERVERTEB"))
    if harm.empty:
        raise AssertionError("AKR1C2 ablation harmonisation is empty")
    chrom = str(int(harm.chr_b37.iloc[0]))
    kept, rho = __import__("plink_ld").ld_matrix(harm, bim_all[bim_all.chr == chrom], str(BED), N_REF)
    order = np.argsort(-kept.f_statistic.to_numpy())
    rows: list[dict[str, object]] = []
    for r2 in R2_LADDER:
        idx = prune_by_ld(rho, order, r2)
        sub = kept.iloc[idx].reset_index(drop=True)
        sr = rho[np.ix_(idx, idx)]
        if sub.empty:
            rows.append({"gene": "AKR1C2", "outcome": "M13_INTERVERTEB", "r2_threshold": r2, "k": 0, "max_leverage": np.nan, "leverage_sum": np.nan, "min_leverage": np.nan, "status": "no_ld_reference"})
            continue
        bx = sub.beta_exposure.to_numpy(float)
        sey = sub.se_outcome.to_numpy(float)
        omega = np.outer(sey, sey) * sr
        leverage = hat_diagonal(bx, omega)
        rows.append({"gene": "AKR1C2", "outcome": "M13_INTERVERTEB", "r2_threshold": r2, "k": int(len(sub)), "max_leverage": float(np.nanmax(leverage)), "leverage_sum": float(np.nansum(leverage)), "min_leverage": float(np.nanmin(leverage)), "status": "ok"})
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "ablation_leverage.tsv", sep="\t", index=False)
    return result


def main() -> None:
    census_frame = census()
    ablation_frame = ablation()
    provenance = {
        "schema_version": "generalized_ivw_leverage_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "metric": "diagonal of one-predictor generalized-IVW hat matrix; summary is maximum diagonal leverage",
        "formula": "H = X (X' W X)^-1 X' W; diagonal(H) sums to one",
        "census_weight": "W = pinv(rho), with X = beta_exposure / se_exposure",
        "ablation_weight": "W = pinv(outer(se_outcome, se_outcome) * rho), with X = beta_exposure",
        "pseudoinverse": "numpy.linalg.pinv with default rcond; no explicit eigenvalue truncation or regularization",
        "inputs": {"full_universe_eqtl": str(EQTL.relative_to(ROOT)), "full_universe_eqtl_sha256": sha256(EQTL), "adequacy": str(ADEQUACY.relative_to(ROOT)), "adequacy_sha256": sha256(ADEQUACY), "ld_panel": str(BED)},
        "outputs": ["out/full_universe_leverage.tsv", "out/ablation_leverage.tsv"],
        "census_rows": int(len(census_frame)),
        "census_max_abs_sum_error": float(np.nanmax(np.abs(census_frame.loc[census_frame.status == "ok", "leverage_sum"] - 1.0))),
        "ablation_max_abs_sum_error": float(np.nanmax(np.abs(ablation_frame.loc[ablation_frame.status == "ok", "leverage_sum"] - 1.0))),
    }
    (OUT / "leverage_metrics_provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(provenance, indent=2, ensure_ascii=False))
    print("census max leverage summary:")
    print(census_frame[census_frame.r2_threshold == 0.001].groupby("status").max(numeric_only=True)[["k", "max_leverage", "leverage_sum"]])
    print("ablation:")
    print(ablation_frame.to_string(index=False))


if __name__ == "__main__":
    main()
