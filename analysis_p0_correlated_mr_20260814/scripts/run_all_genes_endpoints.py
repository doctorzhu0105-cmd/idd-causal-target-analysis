"""LD-aware cis-MR across all five priority genes and all three FinnGen endpoints.

Extends the AKR1C-only P0 analysis so endpoint specificity can be re-checked
under adequate instrument counts. The frozen analysis clumped at r2 < 0.001 and
ended up with 1-4 instruments per gene-outcome pair, which is what made
weighted median, MR-Egger and MR-PRESSO inestimable for most of them.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plink_ld import load_bim, ld_matrix, is_palindromic, _alleles_match  # noqa: E402
from correlated_mr import (  # noqa: E402
    ivw_correlated, egger_correlated, leave_one_out, prune_by_ld,
)

ROOT = Path(__file__).resolve().parents[1]
WS = ROOT.parent
LD_DIR = WS / "data/raw/idd/genetics/ld"
PROC = WS / "data/processed/idd/genetics"
OUT = ROOT / "out"
OUT.mkdir(exist_ok=True)

N_REF = 503
GENES = ["AKR1C2", "AKR1C1", "NT5E", "CCL28", "P2RY13"]
OUTCOMES = ["M13_INTERVERTEB", "M13_LOWBACKPAIN", "M13_SCIATICA"]
# Primary threshold. NOT pre-registered: chosen after inspecting the r2 ladder on
# this data (r2<0.2 introduced heterogeneity and Egger-intercept signal). Treated
# as an analyst degree of freedom — see ANALYSIS_PLAN.md.
R2_PRIMARY = 0.1
R2_LADDER = [0.001, 0.05, 0.1, 0.2]
MIN_N_FRACTION = 0.8


def load_outcome(outcome: str) -> pd.DataFrame:
    """All outcome rows, deliberately NOT deduplicated by rsID.

    1,327 rsIDs in this file are multi-allelic (2,660 rows). Dropping duplicates
    here would keep an arbitrary alt allele and could discard the row that
    actually matches the exposure alleles. Deduplication happens in harmonise(),
    after allele matching has eliminated the non-matching rows.
    """
    fg = pd.read_csv(PROC / f"{outcome}.candidate_variants.tsv.gz", sep="\t")
    fg = fg[["snp", "effect_allele", "other_allele", "beta", "se"]]
    fg.columns = ["snp", "fg_ea", "fg_oa", "fg_beta", "fg_se"]
    return fg


def harmonise(eqtl_all: pd.DataFrame, gene: str, fg: pd.DataFrame) -> pd.DataFrame:
    eqtl = eqtl_all[eqtl_all.gene_symbol == gene]
    max_n = eqtl.n.max()
    sig = eqtl[
        (eqtl.pvalue < 5e-8)
        & (eqtl.f_statistic > 10)
        & (eqtl.n >= MIN_N_FRACTION * max_n)
    ].copy()

    m = sig.merge(fg, on="snp", how="inner")
    if m.empty:
        return m
    m = m[~m.apply(lambda r: is_palindromic(r.effect_allele, r.other_allele), axis=1)]
    if m.empty:
        return m

    orient = [
        _alleles_match(r.effect_allele, r.other_allele, r.fg_ea, r.fg_oa)
        for _, r in m.iterrows()
    ]
    m["o"] = orient
    # allele matching first: at a multi-allelic site this drops the alt alleles
    # that do not correspond to the exposure, leaving the correct row to survive
    # the deduplication below
    m = m[m.o.notna()].copy()
    m["beta_outcome"] = m.fg_beta * m.o
    m["se_outcome"] = m.fg_se
    return m.rename(columns={"beta": "beta_exposure"}).drop_duplicates("snp")


def main() -> None:
    bim_all = load_bim(str(LD_DIR / "1000G_EUR.bim"))
    bed = str(LD_DIR / "1000G_EUR.bed")
    eqtl_all = pd.read_csv(OUT / "all5_full_cis.tsv", sep="\t")

    rows, loo_rows = [], []
    for outcome in OUTCOMES:
        fg = load_outcome(outcome)
        for gene in GENES:
            harm = harmonise(eqtl_all, gene, fg)
            if harm.empty:
                continue
            chrom = str(int(harm.chr_b37.iloc[0]))
            bim = bim_all[bim_all.chr == chrom]
            kept, rho = ld_matrix(harm, bim, bed, N_REF)
            if len(kept) == 0:
                continue
            kept = kept.reset_index(drop=True)
            order = np.argsort(-kept.f_statistic.to_numpy())

            for r2 in R2_LADDER:
                idx = prune_by_ld(rho, order, r2)
                sub = kept.iloc[idx]
                sr = rho[np.ix_(idx, idx)]
                bx = sub.beta_exposure.to_numpy()
                by = sub.beta_outcome.to_numpy()
                sey = sub.se_outcome.to_numpy()

                est = ivw_correlated(bx, by, sey, sr, random_effects=True)
                slope, intercept = egger_correlated(bx, by, sey, sr)

                rows.append(dict(
                    gene=gene, outcome=outcome, r2=r2, k=len(idx),
                    beta=est.beta, se=est.se, pvalue=est.pvalue,
                    or_=est.or_, or_low=est.or_low, or_high=est.or_high,
                    q_pvalue=est.q_pvalue, cond_rho=float(np.linalg.cond(sr)),
                    min_f=float(sub.f_statistic.min()),
                    egger_slope_or=slope.or_, egger_slope_p=slope.pvalue,
                    egger_intercept=intercept.beta,
                    egger_intercept_p=intercept.pvalue,
                ))

                if r2 == R2_PRIMARY:
                    loo = leave_one_out(bx, by, sey, sr)
                    if loo:
                        df = pd.DataFrame(loo)
                        df["snp_dropped"] = sub.snp.to_numpy()[df.dropped_index]
                        df["gene"], df["outcome"] = gene, outcome
                        loo_rows.append(df)
                    sub.assign(outcome=outcome).to_csv(
                        OUT / f"instruments_{gene}_{outcome}_r2_{r2}.tsv",
                        sep="\t", index=False,
                    )

    res = pd.DataFrame(rows).rename(columns={"or_": "or"})
    res.to_csv(OUT / "all_genes_endpoints_mr.tsv", sep="\t", index=False)
    if loo_rows:
        pd.concat(loo_rows).to_csv(OUT / "all_genes_endpoints_loo.tsv", sep="\t", index=False)

    # LOO stability summary at the primary threshold
    summary = []
    if loo_rows:
        for (g, o), d in pd.concat(loo_rows).groupby(["gene", "outcome"]):
            main_beta = res[(res.gene == g) & (res.outcome == o) & (res.r2 == R2_PRIMARY)].beta
            sign = np.sign(main_beta.iloc[0]) if len(main_beta) else np.nan
            summary.append(dict(
                gene=g, outcome=o, k=len(d),
                loo_max_p=float(d.pvalue.max()),
                loo_frac_sig=float((d.pvalue < 0.05).mean()),
                loo_sign_flip=bool((np.sign(d.beta) != sign).any()),
            ))
    pd.DataFrame(summary).to_csv(OUT / "all_genes_endpoints_loo_summary.tsv",
                                 sep="\t", index=False)

    pd.set_option("display.width", 240)
    prim = res[res.r2 == R2_PRIMARY].copy()
    prim["sig"] = np.where(prim.pvalue < 0.05, "*", "")
    print(f"=== LD-aware IVW (random effects), r2 < {R2_PRIMARY} ===")
    print(prim[["gene", "outcome", "k", "or", "or_low", "or_high", "pvalue",
                "q_pvalue", "egger_intercept_p", "min_f", "sig"]].to_string(index=False))
    if summary:
        print("\n=== leave-one-out stability at primary threshold ===")
        print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
