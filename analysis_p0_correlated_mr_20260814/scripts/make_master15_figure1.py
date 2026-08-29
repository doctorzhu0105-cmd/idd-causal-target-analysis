#!/usr/bin/env python3
"""Render the MASTER 1.5 catalogue information/precision figure from public TSVs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {"black": "#000000", "sky": "#56B4E9", "blue": "#0072B2",
          "grey": "#6B6B6B", "vermillion": "#D55E00", "purple": "#8E44AD"}
THRESHOLDS = (0.001, 0.01, 0.05, 0.1, 0.2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = args.repo_root.resolve() / "results" / "master_1_5"
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    ladder = pd.read_csv(source / "CATALOGUE_EFFECTIVE_K.tsv", sep="\t")
    exposure = pd.read_csv(source / "CATALOGUE_EXPOSURE_PAIRED_CHANGES.tsv", sep="\t")
    precision = pd.read_csv(source / "CATALOGUE_PRECISION_PAIRED_CHANGES.tsv", sep="\t")
    ladder = ladder[ladder.status.eq("estimable")]
    exposure = exposure[exposure.contrast.eq("0.001_to_0.1")]
    precision = precision[precision.contrast.eq("0.001_to_0.1")]
    if (len(ladder), len(exposure), len(precision)) != (12190, 2438, 7284):
        raise ValueError("Figure 1 source row-count identity failed")
    if not np.isfinite(exposure[["K_percent_change", "joint_chi_square_percent_change"]]).all().all():
        raise ValueError("Panel B contains non-finite coordinates")
    if not np.isfinite(precision[["joint_chi_square_percent_change", "precision_gain_percent"]]).all().all():
        raise ValueError("Panel C contains non-finite coordinates")

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 9,
                         "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.2,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(8.4, 2.9), constrained_layout=True)
    for axis, label in zip(axes, "ABC"):
        axis.text(-0.17, 1.08, label, transform=axis.transAxes, fontsize=11, fontweight="bold", va="top")
    values = [np.log10(ladder.loc[np.isclose(ladder.threshold, value), "K"].to_numpy(float)) for value in THRESHOLDS]
    parts = axes[0].violinplot(values, positions=np.arange(5), widths=.78, showmedians=True, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor(COLORS["sky"]); body.set_edgecolor(COLORS["blue"]); body.set_alpha(.72)
    parts["cmedians"].set_color(COLORS["black"])
    axes[0].set_xticks(np.arange(5), ["0.001", "0.01", "0.05", "0.1", "0.2"])
    axes[0].set_yticks(np.log10([1, 2, 5, 10, 20, 50]), ["1", "2", "5", "10", "20", "50"])
    axes[0].set(xlabel="LD pruning threshold (r²)", ylabel="Retained instruments (K)",
                title="Full catalogue\nLD-matched genes")
    axes[0].title.set_ha("left"); axes[0].title.set_position((0, 1)); axes[0].grid(axis="y", color="#DDDDDD", linewidth=.45)
    hb = axes[1].hexbin(exposure.K_percent_change, exposure.joint_chi_square_percent_change,
                        gridsize=38, bins="log", mincnt=1, cmap="Blues", linewidths=.15)
    axes[1].axline((0, 0), slope=1, color=COLORS["grey"], linestyle="--", linewidth=.8)
    axes[1].set(xlabel="Nominal K gain (%)", ylabel="Joint information gain (%)",
                title="2,438 paired genes\n0.001 to 0.1")
    axes[1].title.set_ha("left"); axes[1].title.set_position((0, 1))
    axes[1].set_xscale("symlog", linthresh=25); axes[1].set_yscale("symlog", linthresh=10)
    fig.colorbar(hb, ax=axes[1], label="log10 bin count", fraction=.05, pad=.02)
    hb2 = axes[2].hexbin(precision.joint_chi_square_percent_change, precision.precision_gain_percent,
                         gridsize=38, bins="log", mincnt=1, cmap="Greens", linewidths=.15)
    for gene, color in (("AKR1C1", COLORS["vermillion"]), ("AKR1C2", COLORS["purple"])):
        row = precision[(precision.gene.eq(gene)) & (precision.endpoint.eq("M13_INTERVERTEB"))].iloc[0]
        axes[2].scatter(row.joint_chi_square_percent_change, row.precision_gain_percent,
                        s=32, color=color, edgecolor="white", linewidth=.8, zorder=5)
        axes[2].annotate(gene, (row.joint_chi_square_percent_change, row.precision_gain_percent),
                         xytext=(4, 4), textcoords="offset points", fontsize=7, color=color)
    axes[2].axhline(0, color=COLORS["grey"], linestyle="--", linewidth=.8)
    axes[2].set(xlabel="Joint information gain (%)", ylabel="Precision gain (-SE change, %)",
                title="7,284 paired\ngene-endpoint rows")
    axes[2].title.set_ha("left"); axes[2].title.set_position((0, 1))
    axes[2].set_xscale("symlog", linthresh=10); axes[2].set_yscale("symlog", linthresh=10)
    fig.colorbar(hb2, ax=axes[2], label="log10 bin count", fraction=.05, pad=.02)
    stem = output / "Figure_1_catalogue_information_precision"
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=.06)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=.06)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
