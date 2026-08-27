"""Leave-one-out summaries for the 13-row full-universe primary screen.

This add-on reproduces the existing all-gene endpoint LOO implementation at
the locked r2 < 0.1 threshold. It writes new files only; historical LOO
outputs and manuscript files are not modified.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[1]
WS = ROOT.parent
OUT = ROOT / "out"
R2_PRIMARY = 0.1
N_REF = 503
OUTCOMES = ["M13_INTERVERTEB", "M13_LOWBACKPAIN", "M13_SCIATICA"]

sys.path.insert(0, str(SCRIPT.parent))
from correlated_mr import ivw_correlated, leave_one_out, prune_by_ld  # noqa: E402
from plink_ld import _alleles_match, is_palindromic, ld_matrix, load_bim  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_outcome(outcome: str) -> pd.DataFrame:
    source = OUT / f"full_universe_{outcome}.tsv.gz"
    tab = pd.read_csv(source, sep="\t", dtype={"snp": str})
    return tab[["snp", "fg_ea", "fg_oa", "fg_beta", "fg_se"]].copy()


def harmonise(eqtl_all: pd.DataFrame, gene: str, fg: pd.DataFrame) -> pd.DataFrame:
    eqtl = eqtl_all[eqtl_all.gene_symbol == gene].copy()
    if eqtl.empty:
        return eqtl
    max_n = float(eqtl.n.max())
    sig = eqtl[
        (eqtl.pvalue < 5e-8)
        & (eqtl.f_statistic > 10)
        & (eqtl.n >= 0.8 * max_n)
    ].copy()
    merged = sig.merge(fg, on="snp", how="inner")
    if merged.empty:
        return merged
    merged = merged[
        ~merged.apply(
            lambda row: is_palindromic(row.effect_allele, row.other_allele),
            axis=1,
        )
    ].copy()
    if merged.empty:
        return merged
    orientations = [
        _alleles_match(row.effect_allele, row.other_allele, row.fg_ea, row.fg_oa)
        for _, row in merged.iterrows()
    ]
    merged["orientation"] = orientations
    merged = merged[merged.orientation.notna()].copy()
    if merged.empty:
        return merged
    merged["beta_outcome"] = merged.fg_beta * merged.orientation
    merged["se_outcome"] = merged.fg_se
    # Allele matching precedes deduplication, as in the existing endpoint run.
    return merged.rename(columns={"beta": "beta_exposure"}).drop_duplicates("snp")


def run_one(
    gene: str,
    outcome: str,
    eqtl_all: pd.DataFrame,
    fg: pd.DataFrame,
    bim_all: pd.DataFrame,
    bed: str,
) -> tuple[dict, list[dict]]:
    harmonised = harmonise(eqtl_all, gene, fg)
    base = {
        "gene": gene,
        "outcome": outcome,
        "r2_threshold": R2_PRIMARY,
        "n_harmonised": int(len(harmonised)),
        "k": 0,
        "loo_max_p": np.nan,
        "loo_frac_sig": np.nan,
        "loo_sign_flip": np.nan,
        "status": "no_harmonised_variants",
        "reason": "",
    }
    if harmonised.empty:
        return base, []
    chrom = str(int(harmonised.chr_b37.iloc[0]))
    bim = bim_all[bim_all.chr == chrom]
    kept, rho = ld_matrix(harmonised, bim, bed, N_REF)
    if kept.empty:
        base["status"] = "no_ld_reference"
        base["reason"] = "no harmonised variant matched the LD reference panel"
        return base, []
    kept = kept.reset_index(drop=True)
    order = np.argsort(-kept.f_statistic.to_numpy())
    idx = prune_by_ld(rho, order, R2_PRIMARY)
    sub = kept.iloc[idx].reset_index(drop=True)
    sr = rho[np.ix_(idx, idx)]
    bx = sub.beta_exposure.to_numpy()
    by = sub.beta_outcome.to_numpy()
    sey = sub.se_outcome.to_numpy()
    base["k"] = int(len(sub))
    if len(sub) < 3:
        base["status"] = "loo_uninterpretable_k_lt_3"
        base["reason"] = "existing LOO implementation requires at least three retained instruments"
        return base, []
    full = ivw_correlated(bx, by, sey, sr, random_effects=True)
    loo = pd.DataFrame(leave_one_out(bx, by, sey, sr))
    loo["snp_dropped"] = sub.snp.to_numpy()[loo.dropped_index.astype(int)]
    loo["gene"] = gene
    loo["outcome"] = outcome
    loo["r2_threshold"] = R2_PRIMARY
    loo["full_beta"] = full.beta
    loo["full_se"] = full.se
    loo["full_pvalue"] = full.pvalue
    loo["full_or"] = full.or_
    loo["k"] = len(sub)
    base["status"] = "estimable"
    base["loo_max_p"] = float(loo.pvalue.max())
    base["loo_frac_sig"] = float((loo.pvalue < 0.05).mean())
    base["loo_sign_flip"] = bool((np.sign(loo.beta) != np.sign(full.beta)).any())
    return base, loo.to_dict(orient="records")


def update_tiers(loo_summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    downstream = pd.read_csv(OUT / "primary13_downstream.tsv", sep="\t")
    old_loo_cols = ["loo_max_p", "loo_frac_sig", "loo_sign_flip"]
    downstream = downstream.drop(columns=[c for c in old_loo_cols if c in downstream])
    merged = downstream.merge(loo_summary, on=["gene", "outcome"], how="left")
    if "k_x" in merged:
        merged["k"] = merged["k_x"]
        merged = merged.drop(columns=[c for c in ["k_x", "k_y"] if c in merged])
    merged["mr_pass"] = (
        (merged.bh_fdr < 0.05)
        & (merged.q_pvalue > 0.05)
        & (merged.egger_intercept_p > 0.05)
        & (merged.loo_frac_sig == 1.0)
    )
    merged["coloc_pass"] = merged.pp4 >= 0.80
    merged["smr_pass"] = (merged.p_smr_bonferroni < 0.05) & (merged.p_heidi > 0.01)
    merged["tier"] = np.select(
        [
            merged.mr_pass & merged.coloc_pass & merged.smr_pass,
            merged.mr_pass & merged.smr_pass & ~merged.coloc_pass,
            merged.mr_pass & ~merged.smr_pass,
            merged.bh_fdr < 0.05,
        ],
        [
            "A: MR+coloc+SMR three-modal concordance",
            "B: MR+SMR without coloc support",
            "C: MR only",
            "D: significant but sensitivity flagged",
        ],
        default="not assigned",
    )
    loo_unavailable = merged.loo_frac_sig.isna() & (merged.bh_fdr < 0.05)
    merged.loc[loo_unavailable, "tier"] = "D: significant; LOO not interpretable with K < 3"
    tier_cols = [
        "gene", "outcome", "k", "bh_fdr", "q_pvalue", "pp0", "pp1", "pp2", "pp3", "pp4",
        "p_smr", "p_smr_bonferroni", "p_heidi", "loo_max_p", "loo_frac_sig", "loo_sign_flip",
        "mr_pass", "coloc_pass", "smr_pass", "tier", "outside_original_76",
    ]
    return merged, merged[tier_cols].copy()


def main() -> int:
    primary = pd.read_csv(OUT / "full_universe_discoveries_primary.tsv", sep="\t")
    if len(primary) != 13:
        raise AssertionError(f"expected 13 primary rows, got {len(primary)}")
    eqtl = pd.read_csv(OUT / "full_universe_eqtl.tsv.gz", sep="\t")
    bim = load_bim(str(WS / "data/raw/idd/genetics/ld/1000G_EUR.bim"))
    bed = str(WS / "data/raw/idd/genetics/ld/1000G_EUR.bed")
    outcomes = {o: load_outcome(o) for o in OUTCOMES}
    summaries, details = [], []
    for row in primary[["gene", "outcome"]].itertuples(index=False):
        summary, detail = run_one(row.gene, row.outcome, eqtl, outcomes[row.outcome], bim, bed)
        summaries.append(summary)
        details.extend(detail)
    loo_summary = pd.DataFrame(summaries)
    loo_summary.to_csv(OUT / "primary13_loo.tsv", sep="\t", index=False)
    pd.DataFrame(details).to_csv(OUT / "primary13_loo_details.tsv", sep="\t", index=False)
    merged, tiers = update_tiers(loo_summary)
    merged.to_csv(OUT / "primary13_downstream_with_loo.tsv", sep="\t", index=False)
    tiers.to_csv(OUT / "primary13_tiers_with_loo.tsv", sep="\t", index=False)
    provenance = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(SCRIPT.relative_to(ROOT)).replace("\\", "/"),
        "input_primary_rows": 13,
        "r2_threshold": R2_PRIMARY,
        "ld_reference": "1000G EUR",
        "n_reference": N_REF,
        "implementation": "same harmonisation, F-ordering, greedy LD pruning, correlated random-effects IVW and leave_one_out functions as run_all_genes_endpoints.py",
        "inputs": {
            "primary": str(OUT / "full_universe_discoveries_primary.tsv"),
            "primary_sha256": sha256(OUT / "full_universe_discoveries_primary.tsv"),
            "eqtl": str(OUT / "full_universe_eqtl.tsv.gz"),
            "eqtl_sha256": sha256(OUT / "full_universe_eqtl.tsv.gz"),
        },
        "outputs": {
            "summary": str(OUT / "primary13_loo.tsv"),
            "details": str(OUT / "primary13_loo_details.tsv"),
            "downstream_with_loo": str(OUT / "primary13_downstream_with_loo.tsv"),
            "tiers_with_loo": str(OUT / "primary13_tiers_with_loo.tsv"),
        },
    }
    (OUT / "primary13_loo_provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(tiers[["gene", "outcome", "k", "loo_max_p", "loo_frac_sig", "loo_sign_flip", "tier"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
