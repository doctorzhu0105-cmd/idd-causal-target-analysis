"""QWEN independent numerical verifier, adapted for isolated inputs and outputs."""
import json, os
import numpy as np
import pandas as pd

def run(paired, output_dir):
    gz = str(paired)
    OUT = str(output_dir)
    df = pd.read_csv(gz, sep="\t")
    sub = df[(df["K_transition"] == "K1_to_multi")].copy()
    N = len(sub)
    TOL = 1e-15

    v0f = sub["R1_Vdelta0_from"].to_numpy()
    v0t = sub["R1_Vdelta0_to"].to_numpy()
    clf = sub["R1_C_LD_from"].to_numpy()
    clt = sub["R1_C_LD_to"].to_numpy()
    psf = sub["R1_V_upper_PSD_from"].to_numpy()
    pst = sub["R1_V_upper_PSD_to"].to_numpy()

    # 1) rho=0 point decrease
    n_rho0 = int(np.sum(v0t < v0f - TOL))

    # 2) simultaneous universal interval: (v0t-v0f) + 2*r*|clt-clf| < -TOL
    def simul(r):
        gap = (v0t - v0f) + 2.0 * r * np.abs(clt - clf)
        return int(np.sum(gap < -TOL))

    n_0177 = simul(0.0177)
    n_010 = simul(0.10)
    n_1 = simul(1.0)

    # 3) distribution-free PSD Cauchy-Schwarz envelope
    n_psd = int(np.sum(pst < psf - TOL))

    # signed-grid point counts (cross-check disposition table)
    def point(r):
        vf = sub[f"R1_Vdelta_{'neg' if r<0 else 'pos'}{abs(r)*10000:05.0f}_from"] if False else None
        return None
    # use stored signed columns where available
    signed = {}
    for label, r in [("neg00177", -0.0177), ("pos00177", 0.0177),
                     ("neg01000", -0.10), ("pos01000", 0.10),
                     ("neg10000", -1.0), ("pos10000", 1.0)]:
        cf = f"R1_Vdelta_{label}_from"; ct = f"R1_Vdelta_{label}_to"
        if cf in sub.columns:
            signed[r] = int(np.sum(sub[ct].to_numpy() < sub[cf].to_numpy() - TOL))

    res = {
        "scope": "R1_EUR_503_BASELINE / K1_to_multi",
        "N_pairs": N,
        "tolerance": TOL,
        "rho0_point_decrease": {"n": n_rho0, "pct": 100.0 * n_rho0 / N, "claim_pct": 98.0625},
        "simul_abs_rho_le_0.0177": {"n": n_0177, "pct": 100.0 * n_0177 / N, "claim_pct": 97.9439},
        "simul_abs_rho_le_0.10": {"n": n_010, "pct": 100.0 * n_010 / N, "claim_pct": 97.5880},
        "simul_abs_rho_le_1.0": {"n": n_1, "pct": 100.0 * n_1 / N, "claim_pct": 88.8098},
        "psd_envelope": {"n": n_psd, "pct": 100.0 * n_psd / N, "claim_pct": 68.7228},
        "signed_grid_point_counts": signed,
        "monotonicity_check": n_rho0 >= n_0177 >= n_010 >= n_1 >= n_psd,
    }
    with open(os.path.join(OUT, "qwen_wp3_verification_results.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)
    print(json.dumps(res, indent=2, ensure_ascii=False))