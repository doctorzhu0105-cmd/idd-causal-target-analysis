"""Render current Figure 1 from the redistributable MASTER 1.4 source table."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


COLORS = {
    "black": "#000000",
    "sky": "#56B4E9",
    "green": "#009E73",
    "orange": "#E69F00",
    "vermillion": "#D55E00",
    "blue": "#0072B2",
}
THRESHOLDS = [0.001, 0.01, 0.05, 0.1, 0.2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.repo_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(root / "results/master_1_4/figure1_source.tsv", sep="\t")
    if len(data) != 2479 * len(THRESHOLDS):
        raise ValueError(f"expected 12,395 rows, found {len(data)}")
    historical = data.loc[data.group.eq("historical_candidate"), "gene"].nunique()
    if historical != 70:
        raise ValueError(f"expected 70 data-eligible historical candidates, found {historical}")

    bins = ("LD reference unavailable", "K=1", "K=2", "K=3", "K≥4")
    bin_colors = (COLORS["black"], COLORS["sky"], COLORS["green"], COLORS["orange"], COLORS["vermillion"])
    counts, historical_counts, leverage, historical_leverage = [], [], [], []
    for threshold in THRESHOLDS:
        part = data.loc[np.isclose(data.r2_threshold, threshold)]
        subset = part.loc[part.group.eq("historical_candidate")]

        def count_rows(frame: pd.DataFrame) -> list[int]:
            return [
                int(frame.status.eq("no_ld_reference").sum()),
                int(frame.k.eq(1).sum()), int(frame.k.eq(2).sum()),
                int(frame.k.eq(3).sum()), int(frame.k.ge(4).sum()),
            ]

        counts.append(count_rows(part))
        historical_counts.append(count_rows(subset))
        leverage.append(part.loc[part.status.eq("ok") & part.k.ge(3), "max_leverage"].dropna().tolist())
        historical_leverage.append(subset.loc[subset.status.eq("ok") & subset.k.ge(3), "max_leverage"].dropna().tolist())

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.55), constrained_layout=False)
    fig.subplots_adjust(top=0.90, bottom=0.25, left=0.10, right=0.98, wspace=0.30)
    for index, label in enumerate(("A", "B")):
        axes[index].text(-0.16, 1.08, label, transform=axes[index].transAxes, fontsize=12, fontweight="bold", va="top")

    x = np.arange(len(THRESHOLDS))
    bottom = np.zeros(len(THRESHOLDS))
    historical_bottom = np.zeros(len(THRESHOLDS))
    for index, color in enumerate(bin_colors):
        values = np.array([row[index] for row in counts], dtype=float)
        axes[0].bar(x - 0.18, values, bottom=bottom, color=color, edgecolor="white", linewidth=0.45, width=0.32)
        bottom += values
        values = np.array([row[index] for row in historical_counts], dtype=float)
        axes[0].bar(x + 0.18, values, bottom=historical_bottom, color=color, edgecolor="#333333", linewidth=0.35, width=0.32, hatch="//", alpha=0.72)
        historical_bottom += values
    labels = ["<0.001", "<0.01", "<0.05", "<0.1", "<0.2"]
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels)
    axes[0].set_xlabel("LD threshold (r²)")
    axes[0].set_ylabel("Data-eligible genes")
    axes[0].set_ylim(0, 2600)

    positions = np.arange(1, len(THRESHOLDS) + 1)
    axes[1].boxplot(
        leverage, positions=positions - 0.11, widths=0.18, patch_artist=True, showfliers=False,
        medianprops={"color": COLORS["black"], "linewidth": 1.1},
        whiskerprops={"color": "#555555", "linewidth": 0.8},
        capprops={"color": "#555555", "linewidth": 0.8},
        boxprops={"facecolor": COLORS["sky"], "edgecolor": "#555555", "linewidth": 0.8},
    )
    axes[1].boxplot(
        historical_leverage, positions=positions + 0.11, widths=0.18, patch_artist=True, showfliers=False,
        medianprops={"color": COLORS["black"], "linewidth": 1.1},
        whiskerprops={"color": "#555555", "linewidth": 0.8},
        capprops={"color": "#555555", "linewidth": 0.8},
        boxprops={"facecolor": COLORS["orange"], "edgecolor": "#555555", "linewidth": 0.8},
    )
    for position, values in zip(positions - 0.11, leverage):
        offsets = np.linspace(-0.05, 0.05, len(values)) if values else np.array([])
        axes[1].scatter(np.full(len(values), position) + offsets, values, s=7,
                        color=COLORS["blue"], alpha=0.32, linewidth=0, zorder=3)
    for position, values in zip(positions + 0.11, historical_leverage):
        offsets = np.linspace(-0.05, 0.05, len(values)) if values else np.array([])
        axes[1].scatter(np.full(len(values), position) + offsets, values, s=9,
                        color=COLORS["vermillion"], alpha=0.55, linewidth=0, zorder=3)
    axes[1].set_xticks(positions)
    axes[1].set_xticklabels(labels)
    axes[1].set_xlabel("LD threshold (r²)")
    axes[1].set_ylabel("Maximum generalized-IVW leverage")
    maximum = max(max(values) for values in leverage + historical_leverage if values)
    axes[1].set_ylim(0, max(1.05, maximum * 1.08))
    axes[1].set_title("Genes retaining at least three instruments", loc="left", fontsize=8)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#333333")
        ax.spines["bottom"].set_color("#333333")
        ax.tick_params(direction="out", length=3, width=0.7, labelsize=8)
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.45, alpha=0.65)

    fig.legend([Patch(facecolor=color, edgecolor="white") for color in bin_colors], bins,
               frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.115), ncol=5,
               handlelength=1.1, columnspacing=0.9, borderaxespad=0.0, fontsize=7)
    fig.legend([
        Patch(facecolor=COLORS["sky"], edgecolor="#555555"),
        Patch(facecolor=COLORS["orange"], edgecolor="#555555"),
    ], ["Panel B: all data-eligible genes", "Panel B: 70 data-eligible historical candidates"],
       frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.025), ncol=2,
       handlelength=1.1, columnspacing=1.2, borderaxespad=0.0, fontsize=7)
    stem = output / "Figure_1_instrument_adequacy"
    fig.savefig(stem.with_suffix(".png"), dpi=610, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
