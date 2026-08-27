"""Run registered-prior coloc sensitivity for the current raw-region A tier."""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[1]
OUT = ROOT / "out"
RAW_TIERS = OUT / "coloc_smr_authoritative_raw_tiers.tsv"
RAW_CONFIG = ROOT / "config" / "coloc_smr_authoritative_raw.json"
PAIRS = OUT / "coloc_prior_sensitivity_pairs.tsv"
OUTPUT = OUT / "coloc_prior_sensitivity.tsv"
PROV = OUT / "coloc_prior_sensitivity_provenance.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    x = pd.read_csv(RAW_TIERS, sep="\t")
    top = x[x.tier.astype(str).str.startswith("A:")][["gene", "ensembl_gene_id", "outcome", "tier", "pp4"]].copy()
    top = top.rename(columns={"gene": "gene_symbol", "outcome": "outcome_id"})
    if len(top) != 4:
        raise RuntimeError(f"expected four current raw-region A-tier combinations, found {len(top)}")
    top.to_csv(PAIRS, sep="\t", index=False)
    rscript = shutil.which("Rscript") or "Rscript"
    r_root = Path("D:/idd_ws_tmp/analysis_p0_correlated_mr_20260814") if Path("D:/idd_ws_tmp/analysis_p0_correlated_mr_20260814").exists() else ROOT
    r = subprocess.run([rscript, str(SCRIPT.with_suffix(".R")), "--repo-root", r_root.as_posix(), "--pairs", (r_root / "out/coloc_prior_sensitivity_pairs.tsv").as_posix(), "--output", (r_root / "out/coloc_prior_sensitivity.tsv").as_posix()], cwd=str(ROOT), capture_output=True, text=True, timeout=1800)
    (OUT / "coloc_prior_sensitivity.stdout.log").write_text(r.stdout, encoding="utf-8")
    (OUT / "coloc_prior_sensitivity.stderr.log").write_text(r.stderr, encoding="utf-8")
    if r.returncode != 0:
        raise RuntimeError(f"coloc prior sensitivity failed ({r.returncode})")
    result = pd.read_csv(OUTPUT, sep="\t")
    registered = result[result["p12"].astype(float) == 1e-5].copy()
    check = registered.merge(top[["gene_symbol", "outcome_id", "pp4"]], on=["gene_symbol", "outcome_id"], suffixes=("_rerun", "_authority"))
    if len(check) != 4 or (check["pp4_rerun"].astype(float) - check["pp4_authority"].astype(float)).abs().max() > 1e-8:
        raise AssertionError("registered-p12 rerun does not reproduce the raw-region authority")
    intervals = []
    for (g, o), d in result.groupby(["gene_symbol", "outcome_id"]):
        ok = d.loc[d.p4_ge_0_80 == True, "p12"]
        intervals.append({"gene_symbol": g, "outcome_id": o, "pp4_ge_0_80_p12_min": float(ok.min()) if len(ok) else None, "pp4_ge_0_80_p12_max": float(ok.max()) if len(ok) else None, "n_grid_values": int(len(d))})
    prov = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(), "script": str(SCRIPT),
        "registered_p12": 1e-5, "registered_p1": 1e-4, "registered_p2": 1e-4,
        "p12_grid": [1e-7,3e-7,1e-6,3e-6,1e-5,3e-5,1e-4],
        "highest_tier_pairs": top.to_dict(orient="records"), "robustness_intervals": intervals,
        "selection_rule": "all and only rows whose current raw authoritative tier starts with A:",
        "stale_selection_removed": "integrated_evidence.tsv and P2RY13-specific historical selection are not read",
        "inputs": {
            "raw_tiers": str(RAW_TIERS), "raw_tiers_sha256": sha256(RAW_TIERS),
            "raw_config": str(RAW_CONFIG), "raw_config_sha256": sha256(RAW_CONFIG),
            "pairs": str(PAIRS), "pairs_sha256": sha256(PAIRS),
        },
        "registered_p12_reproduces_authority_max_abs_pp4_difference": float((check["pp4_rerun"].astype(float) - check["pp4_authority"].astype(float)).abs().max()),
        "output_sha256": sha256(OUTPUT),
    }
    PROV.write_text(json.dumps(prov, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(result.to_string(index=False)); print(pd.DataFrame(intervals).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
