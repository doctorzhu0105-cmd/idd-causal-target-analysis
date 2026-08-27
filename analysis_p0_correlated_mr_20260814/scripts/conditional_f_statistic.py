"""Exposure-side joint-strength diagnostics for the five main genes.

For one endogenous exposure, the Sanderson--Windmeijer conditional-F idea
reduces to the joint first-stage strength of the exposure effects after
accounting for their covariance.  Here the covariance is
Sigma_x = diag(SE_x) rho diag(SE_x), and the reported statistic is

    F_joint = (b_x' Sigma_x^{-1} b_x) / K.

This is reported as a correlated-instrument joint-strength analogue (not a
literal multi-endogenous-variable conditional F), not as a new causal result.
The implementation uses the same allele harmonisation,
1000G EUR LD matrix, F-ordering and greedy pruning as the full-universe screen.
The statistic is motivated by Sanderson & Windmeijer, Journal of Econometrics
(2016), with the one-endogenous-exposure reduction made explicit here.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
WS = ROOT.parent
OUT = ROOT / "out"
LD_DIR = WS / "data" / "raw" / "idd" / "genetics" / "ld"
EQTL = OUT / "full_universe_eqtl.tsv.gz"
MAIN = ["NT5E", "AKR1C1", "AKR1C2", "P2RY13", "CCL28"]
OUTCOMES = ["M13_INTERVERTEB", "M13_LOWBACKPAIN", "M13_SCIATICA"]
R2_LADDER = [0.001, 0.01, 0.05, 0.1, 0.2]
N_REF = 503

sys.path.insert(0, str(SCRIPT_DIR))
from correlated_mr import prune_by_ld  # noqa: E402
from plink_ld import ld_matrix, load_bim  # noqa: E402
from run_all_genes_endpoints import harmonise  # noqa: E402


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    eqtl = pd.read_csv(EQTL, sep="\t")
    bim_all = load_bim(str(LD_DIR / "1000G_EUR.bim"))
    rows: list[dict[str, object]] = []
    for outcome in OUTCOMES:
        fg = pd.read_csv(OUT / f"full_universe_{outcome}.tsv.gz", sep="\t")
        for gene in MAIN:
            harm = harmonise(eqtl, gene, fg)
            if harm.empty:
                for r2 in R2_LADDER:
                    rows.append({"gene": gene, "outcome": outcome, "r2_threshold": r2, "k": 0, "conditional_f": np.nan, "joint_chi2": np.nan, "min_marginal_f": np.nan, "ld_condition_number": np.nan, "status": "no_harmonized_instruments"})
                continue
            chrom = str(int(harm.chr_b37.iloc[0]))
            kept, rho = ld_matrix(harm, bim_all[bim_all.chr == chrom], str(LD_DIR / "1000G_EUR.bed"), N_REF)
            kept = kept.reset_index(drop=True)
            order = np.argsort(-kept.f_statistic.to_numpy())
            for r2 in R2_LADDER:
                idx = prune_by_ld(rho, order, r2)
                sub = kept.iloc[idx].reset_index(drop=True)
                sr = rho[np.ix_(idx, idx)]
                k = len(sub)
                if k == 0:
                    rows.append({"gene": gene, "outcome": outcome, "r2_threshold": r2, "k": 0, "conditional_f": np.nan, "joint_chi2": np.nan, "min_marginal_f": np.nan, "ld_condition_number": np.nan, "status": "no_ld_reference"})
                    continue
                bx = sub.beta_exposure.to_numpy(float)
                # ``harmonise`` renames beta to beta_exposure but keeps the
                # eQTL standard error column as ``se``.
                sex = sub.se.to_numpy(float)
                sigma_x = np.outer(sex, sex) * sr
                sigma_inv = np.linalg.pinv(sigma_x)
                joint_chi2 = float(bx @ sigma_inv @ bx)
                conditional_f = joint_chi2 / k
                rows.append({
                    "gene": gene, "outcome": outcome, "r2_threshold": r2,
                    "k": int(k), "conditional_f": conditional_f,
                    "joint_chi2": joint_chi2,
                    "min_marginal_f": float(sub.f_statistic.min()),
                    "ld_condition_number": float(np.linalg.cond(sr)),
                    "status": "estimable",
                })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "conditional_f_statistic.tsv", sep="\t", index=False)
    provenance = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": "conditional_f_statistic.py",
        "method": "joint correlated-instrument first-stage strength: b_x' Sigma_x^-1 b_x / K",
        "interpretation": "correlated-instrument joint-strength analogue motivated by the one-exposure reduction of the Sanderson-Windmeijer framework; diagnostic only, not a literal multi-endogenous-variable conditional F",
        "sigma_x": "diag(SE_x) * rho * diag(SE_x)",
        "matrix_inverse": "numpy.linalg.pinv with NumPy default rcond; no explicit eigenvalue truncation or regularization",
        "reference": "Sanderson and Windmeijer, Journal of Econometrics (2016), conditional F-statistics for IV models",
        "inputs": {"eqtl": str(EQTL), "eqtl_sha256": sha256_file(EQTL), "outcomes": OUTCOMES, "r2_ladder": R2_LADDER, "ld_panel": str(LD_DIR / "1000G_EUR.bed")},
        "rows": int(len(out)),
    }
    (OUT / "conditional_f_statistic_provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
