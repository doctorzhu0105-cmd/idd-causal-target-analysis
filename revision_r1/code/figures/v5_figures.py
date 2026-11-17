from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch
from matplotlib.scale import SymmetricalLogTransform
from matplotlib.ticker import FixedLocator, FixedFormatter, LogFormatterMathtext, LogLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

CAT_DIR = FIG2_DATA = S2_DATA = OUT_DIR = None

COLORS = {
    "black": "#111111",
    "sky": "#56B4E9",
    "blue": "#0072B2",
    "grey": "#555555",
    "vermillion": "#D55E00",
    "orange": "#E69F00",
    "green": "#009E73",
    "purple": "#8E44AD",
}
THRESHOLDS = [0.001, 0.01, 0.05, 0.1, 0.2]


def configure_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 9.0,
        "xtick.labelsize": 7.8,
        "ytick.labelsize": 7.8,
        "legend.fontsize": 7.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "savefig.dpi": 600,
    })


def save_both(fig: plt.Figure, stem: str) -> None:
    for ext in ("png", "tiff"):
        p = OUT_DIR / f"{stem}.{ext}"
        kwargs = {"format": ext, "dpi": 600, "bbox_inches": "tight", "pad_inches": 0.08}
        if ext == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(p, **kwargs)
    plt.close(fig)


def make_figure1() -> None:
    configure_style()
    ladder = pd.read_csv(CAT_DIR / "CATALOGUE_EFFECTIVE_K.tsv", sep="\t")
    exposure = pd.read_csv(CAT_DIR / "CATALOGUE_EXPOSURE_PAIRED_CHANGES.tsv", sep="\t")
    precision = pd.read_csv(CAT_DIR / "CATALOGUE_PRECISION_PAIRED_CHANGES.tsv", sep="\t")

    ladder = ladder[ladder["status"].eq("estimable")].copy()
    exposure = exposure[exposure["contrast"].eq("0.001_to_0.1")].copy()
    precision = precision[precision["contrast"].eq("0.001_to_0.1")].copy()
    assert (len(ladder), len(exposure), len(precision)) == (12190, 2438, 7284)

    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.35), gridspec_kw={"wspace": 0.42})

    # Panel A: Violin plot of log10(K)
    ax0 = axes[0]
    ax0.text(-0.18, 1.09, "A", transform=ax0.transAxes, fontsize=11.5, fontweight="bold", va="top")
    vals = [np.log10(ladder.loc[np.isclose(ladder["threshold"], t), "K"].to_numpy(float)) for t in THRESHOLDS]
    parts = ax0.violinplot(vals, positions=np.arange(5), widths=0.76, showmedians=True, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor(COLORS["sky"])
        body.set_edgecolor(COLORS["blue"])
        body.set_linewidth(0.8)
        body.set_alpha(0.75)
    parts["cmedians"].set_color(COLORS["black"])
    parts["cmedians"].set_linewidth(1.5)
    ax0.set_xticks(np.arange(5), ["0.001", "0.01", "0.05", "0.1", "0.2"])
    yticks_k = [1, 2, 5, 10, 20, 50, 100]
    ax0.set_yticks(np.log10(yticks_k), [str(k) for k in yticks_k])
    ax0.set_ylim(-0.08, np.log10(130))
    ax0.set_xlabel("LD pruning threshold ($\\tau$)")
    ax0.set_ylabel("Retained instruments (K)")
    ax0.set_title("Full catalogue (2,438 LD-matched genes)", loc="left", fontsize=8.8)
    ax0.grid(axis="y", color="#DDDDDD", linewidth=0.45)
    ax0.set_axisbelow(True)

    # Panel B: Symlog-space fine hexagonal binning for Exposure-side I_X gain vs Nominal K gain
    ax1 = axes[1]
    ax1.text(-0.20, 1.09, "B", transform=ax1.transAxes, fontsize=11.5, fontweight="bold", va="top")
    tx_b = SymmetricalLogTransform(base=10, linthresh=25.0, linscale=1.0)
    ty_b = SymmetricalLogTransform(base=10, linthresh=10.0, linscale=1.0)
    xb = tx_b.transform(exposure["K_percent_change"].to_numpy(float))
    yb = ty_b.transform(exposure["joint_chi_square_percent_change"].to_numpy(float))

    x_pad_b = 0.08 * (xb.max() - xb.min())
    y_pad_b = 0.08 * (yb.max() - yb.min())
    xlim_b = (xb.min() - x_pad_b, xb.max() + x_pad_b)
    ylim_b = (yb.min() - y_pad_b, yb.max() + y_pad_b)

    hb1 = ax1.hexbin(xb, yb, gridsize=34, bins="log", mincnt=1, cmap="Blues",
                     extent=(xlim_b[0], xlim_b[1], ylim_b[0], ylim_b[1]),
                     linewidths=0.15, edgecolors="#EBF3FA", zorder=2)

    # Descriptive y = x equality line mapped into transformed symlog coordinates
    u_eq = np.linspace(0.0, 350.0, 600)
    ux = tx_b.transform(u_eq)
    uy = ty_b.transform(u_eq)
    mask_eq = (ux >= xlim_b[0]) & (ux <= xlim_b[1]) & (uy >= ylim_b[0]) & (uy <= ylim_b[1])
    ax1.plot(ux[mask_eq], uy[mask_eq], color=COLORS["grey"], linestyle="--", linewidth=0.95, zorder=4)

    ax1.set_xlim(xlim_b)
    ax1.set_ylim(ylim_b)
    xt_b = [0, 10, 100, 1000]
    yt_b = [0, 10, 100]
    ax1.set_xticks(tx_b.transform(xt_b), ["0", "10", "$10^2$", "$10^3$"])
    ax1.set_yticks(ty_b.transform(yt_b), ["0", "10", "$10^2$"])
    ax1.set_xlabel("Nominal K gain (%)")
    ax1.set_ylabel("Exposure-side joint info. gain (%)")
    ax1.set_title("2,438 paired genes ($\\tau$: 0.001 to 0.1)", loc="left", fontsize=8.8)

    div1 = make_axes_locatable(ax1)
    cax1 = div1.append_axes("right", size="5%", pad=0.07)
    cb1 = fig.colorbar(hb1, cax=cax1)
    cb1.set_label("Bin count (log scale)", fontsize=7.5)
    cb1.ax.tick_params(labelsize=7.0)

    # Panel C: Symlog-space fine hexagonal binning for Random-effects precision gain vs Endpoint-matched I_X gain
    ax2 = axes[2]
    ax2.text(-0.22, 1.09, "C", transform=ax2.transAxes, fontsize=11.5, fontweight="bold", va="top")
    tx_c = SymmetricalLogTransform(base=10, linthresh=10.0, linscale=1.0)
    ty_c = SymmetricalLogTransform(base=10, linthresh=10.0, linscale=1.0)
    xc = tx_c.transform(precision["joint_chi_square_percent_change"].to_numpy(float))
    yc = ty_c.transform(precision["precision_gain_percent"].to_numpy(float))

    x_pad_c = 0.08 * (xc.max() - xc.min())
    y_pad_c = 0.08 * (yc.max() - yc.min())
    xlim_c = (xc.min() - x_pad_c, xc.max() + x_pad_c)
    ylim_c = (yc.min() - y_pad_c, yc.max() + y_pad_c)

    hb2 = ax2.hexbin(xc, yc, gridsize=36, bins="log", mincnt=1, cmap="Greens",
                     extent=(xlim_c[0], xlim_c[1], ylim_c[0], ylim_c[1]),
                     linewidths=0.15, edgecolors="#ECF7F0", zorder=2)

    ax2.axhline(ty_c.transform([0.0])[0], color=COLORS["grey"], linestyle="--", linewidth=0.95, zorder=3)

    for gene, color, offset, ha_val in (
        ("AKR1C1", COLORS["vermillion"], (-8, -5), "right"),
        ("AKR1C2", COLORS["purple"], (7, 5), "left"),
    ):
        row = precision[(precision["gene"].eq(gene)) & (precision["endpoint"].eq("M13_INTERVERTEB"))].iloc[0]
        gx = float(tx_c.transform([row["joint_chi_square_percent_change"]])[0])
        gy = float(ty_c.transform([row["precision_gain_percent"]])[0])
        ax2.scatter(gx, gy, s=38, color=color, edgecolor="white", linewidth=0.9, zorder=5)
        ax2.annotate(
            f"{gene}\n(disc)",
            (gx, gy),
            xytext=offset,
            textcoords="offset points",
            fontsize=7.0,
            fontweight="bold",
            color=color,
            ha=ha_val,
            va="top" if gene == "AKR1C1" else "bottom",
            bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="none", alpha=0.82),
            zorder=6,
        )

    ax2.set_xlim(xlim_c)
    ax2.set_ylim(ylim_c)
    xt_c = [-10, 0, 10, 100]
    yt_c = [-100, -10, 0, 10, 50]
    ax2.set_xticks(tx_c.transform(xt_c), ["$-10$", "0", "10", "$10^2$"])
    ax2.set_yticks(ty_c.transform(yt_c), ["$-10^2$", "$-10$", "0", "10", "50"])
    ax2.set_xlabel("Endpoint-matched info. gain (%)")
    ax2.set_ylabel("Random-effects precision gain ($-$SE change, %)")
    ax2.set_title("7,284 paired gene–endpoint rows", loc="left", fontsize=8.8)

    div2 = make_axes_locatable(ax2)
    cax2 = div2.append_axes("right", size="5%", pad=0.07)
    cb2 = fig.colorbar(hb2, cax=cax2)
    cb2.set_label("Bin count (log scale)", fontsize=7.5)
    cb2.ax.tick_params(labelsize=7.0)

    fig.subplots_adjust(left=0.065, right=0.96, top=0.86, bottom=0.16, wspace=0.44)
    save_both(fig, "Figure_1")


def make_figure2() -> None:
    configure_style()
    df = pd.read_csv(FIG2_DATA, sep="\t")
    assert len(df) == 5
    x = df["threshold"].to_numpy(float)
    k = df["K"].to_numpy(int)
    odds = df["odds_ratio"].to_numpy(float)
    se = df["se_log_odds"].to_numpy(float)
    loo = df["max_loo_p"].to_numpy(float)
    flip = df["loo_sign_reversal"].astype(str).str.lower().eq("true").to_numpy()
    lev = df["max_generalized_leverage"].to_numpy(float)

    lo = np.exp(np.log(odds) - 1.96 * se)
    hi = np.exp(np.log(odds) + 1.96 * se)

    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.15), constrained_layout=True)
    fig.suptitle("AKR1C2 (M13_INTERVERTEB) — non-representative high-count, limited-information configuration",
                 fontsize=9.2, fontweight="bold", y=1.03)

    for ax, label in zip(axes, "ABC"):
        ax.text(-0.18, 1.09, label, transform=ax.transAxes, fontsize=11.5, fontweight="bold", va="top")
        ax.set_xscale("log")
        ax.set_xlim(6.5e-4, 0.30)
        ax.set_xlabel("LD pruning threshold ($\\tau$)")
        ax.xaxis.set_major_locator(LogLocator(base=10))
        ax.xaxis.set_major_formatter(LogFormatterMathtext())
        ax.grid(axis="y", color="#E2E2E2", linewidth=0.45)
        ax.set_axisbelow(True)

    # Panel A: Odds ratio + 95% CI, K labels placed cleanly above upper CI caps
    ax0 = axes[0]
    ax0.errorbar(x, odds, yerr=[odds - lo, hi - odds], marker="o", markersize=4.8,
                 color=COLORS["blue"], ecolor=COLORS["blue"], capsize=2.5, linewidth=1.15, zorder=3)
    ax0.axhline(1.0, color=COLORS["grey"], linestyle="--", linewidth=0.85, zorder=1)
    for xx, hiy, kk in zip(x, hi, k):
        ax0.annotate(f"K={kk}", (xx, hiy), xytext=(0, 3.5), textcoords="offset points",
                     ha="center", va="bottom", fontsize=7.0, color="#222222")
    ax0.set_ylim(0.922, 1.006)
    ax0.set_ylabel("Odds ratio (95% CI)")
    ax0.set_title("Estimate and retained K", loc="left", fontsize=8.8)

    # Panel B: Max LOO P with P=0.05 line clearly explained
    ax1 = axes[1]
    ax1.set_yscale("log")
    ax1.plot(x, loo, color=COLORS["orange"], linewidth=1.2, zorder=2)
    for xx, yy, ff in zip(x, loo, flip):
        ax1.scatter(xx, yy, marker="D" if ff else "o", s=34 if ff else 30,
                    facecolor=COLORS["orange"],
                    edgecolor=COLORS["vermillion"] if ff else "white",
                    linewidth=1.0 if ff else 0.7, zorder=4)
    ax1.axhline(0.05, color=COLORS["grey"], linestyle="--", linewidth=0.9, zorder=1)
    ax1.text(0.24, 0.05 * 1.18, "P = 0.05", ha="right", va="bottom", fontsize=7.0, color=COLORS["grey"])
    rev_idx = int(np.flatnonzero(flip)[0])
    ax1.annotate("Sign reversal", (x[rev_idx], loo[rev_idx]), xytext=(8, -14),
                 textcoords="offset points", fontsize=7.0, color=COLORS["vermillion"],
                 fontweight="bold", ha="left", va="top",
                 arrowprops=dict(arrowstyle="-", color=COLORS["vermillion"], lw=0.75))
    ax1.set_ylim(1.5e-4, 2.2)
    ax1.set_ylabel("Maximum leave-one-out P")
    ax1.set_title("Single-variant deletion", loc="left", fontsize=8.8)

    # Panel C: Endpoint-specific leverage
    ax2 = axes[2]
    ax2.plot(x, lev, marker="o", markersize=4.8, color=COLORS["green"], linewidth=1.2, zorder=3)
    for xx, yy in zip(x, lev):
        ax2.annotate(f"{yy:.3f}", (xx, yy), xytext=(0, 5.5), textcoords="offset points",
                     ha="center", va="bottom", fontsize=6.8, color="#111111")
    ax2.set_ylim(0.905, 1.008)
    ax2.set_ylabel("Maximum endpoint-specific leverage")
    ax2.set_title("Endpoint-specific GIVW leverage", loc="left", fontsize=8.8)

    save_both(fig, "Figure_2")


def make_supp_s1() -> None:
    configure_style()
    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    def draw_box(cx, cy, w, h, text, edge, fc="white", fs=8.5, bold=False):
        x0, y0 = cx - w / 2, cy - h / 2
        box = FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0.012,rounding_size=0.015",
                             facecolor=fc, edgecolor=edge, linewidth=1.4)
        ax.add_patch(box)
        ax.text(cx, cy, text, ha="center", va="center", fontsize=fs,
                fontweight="bold" if bold else "normal", color="#111111", linespacing=1.28)

    def arrow(p1, p2):
        ax.annotate("", xy=p2, xytext=p1,
                    arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.25",
                                    color="#444444", lw=1.2, shrinkA=0, shrinkB=0))

    # Top flow row (y = 0.78)
    w_box, h_box = 0.195, 0.165
    y_top = 0.78
    xs = [0.115, 0.370, 0.625, 0.880]
    draw_box(xs[0], y_top, w_box, h_box, "Finan druggable\ncatalogue\n4,479 genes", COLORS["blue"], bold=True)
    draw_box(xs[1], y_top, w_box, h_box, "Raw cis-eQTL\nrecords present\n2,902 genes", COLORS["sky"], bold=True)
    draw_box(xs[2], y_top, w_box, h_box, "Data-eligible\ngenes\n2,479 genes", COLORS["green"], bold=True)
    draw_box(xs[3], y_top, w_box, h_box, "LD-reference\nmatched\n2,438 genes", COLORS["orange"], bold=True)

    # Horizontal arrows between top boxes (with gap for clean arrowhead)
    for i in range(3):
        arrow((xs[i] + w_box / 2 + 0.014, y_top), (xs[i + 1] - w_box / 2 - 0.014, y_top))

    # Upstream exclusion branch below 2,902 -> 2,479
    y_mid = 0.44
    draw_box(xs[1], y_mid, 0.235, 0.17,
             "Upstream exclusions\n423 genes\n(351 P; 71 rel.-N; 1 palindr.)", COLORS["grey"], fs=8.0)
    arrow((xs[1], y_top - h_box / 2 - 0.014), (xs[1], y_mid + 0.17 / 2 + 0.014))

    # LD unavailable branch cleanly below 2,438 box, routed vertically from Data-eligible? Wait: it branches from 2,479!
    # Route cleanly: orthogonal drop from right side of Data-eligible or dedicated diagonal that never touches box text
    draw_box(xs[3], y_mid, 0.205, 0.16,
             "LD reference\nunavailable\n41 genes", COLORS["grey"], fs=8.3)
    # Orthogonal L-shaped connector from bottom-right of Data-eligible (2,479) into left edge of LD unavailable
    elbow_x = (xs[2] + w_box / 2 + xs[3] - 0.205 / 2) / 2
    ax.plot([xs[2] + 0.065, xs[2] + 0.065, xs[3] - 0.205 / 2 - 0.014],
            [y_top - h_box / 2 - 0.012, y_mid, y_mid], color="#444444", lw=1.2)
    arrow((xs[3] - 0.205 / 2 - 0.025, y_mid), (xs[3] - 0.205 / 2 - 0.014, y_mid))

    # Bottom partition from Data-eligible (2,479)
    y_bot = 0.16
    draw_box(0.505, y_bot, 0.22, 0.155, "Historical candidate list\n70 data-eligible genes", COLORS["purple"], fs=8.2)
    draw_box(0.765, y_bot, 0.22, 0.155, "Other data-eligible genes\n2,409 genes", COLORS["sky"], fs=8.2)

    # Clean vertical stem down from 2,479 then horizontal split to 70 and 2,409
    y_split = 0.30
    ax.plot([xs[2], xs[2]], [y_top - h_box / 2 - 0.012, y_split], color="#444444", lw=1.2)
    ax.plot([0.505, 0.765], [y_split, y_split], color="#444444", lw=1.2)
    arrow((0.505, y_split), (0.505, y_bot + 0.155 / 2 + 0.014))
    arrow((0.765, y_split), (0.765, y_bot + 0.155 / 2 + 0.014))

    ax.text(0.50, 0.02,
            "Note: The 70 and 2,409 groups partition the 2,479 data-eligible genes prior to LD-reference matching.",
            ha="center", va="bottom", fontsize=7.8, color="#444444")

    save_both(fig, "Supplementary_Figure_S1")


def make_supp_s2() -> None:
    configure_style()
    s2 = pd.read_csv(S2_DATA, sep="\t")
    assert len(s2) == 13
    s2["Endpoint"] = s2["Endpoint"].str.replace("Sciatica+with lumbago", "Sciatica with lumbago", regex=False)
    s2 = s2.iloc[::-1].reset_index(drop=True)

    ep_colors = {
        "M13_INTERVERTEB": "#2166ac",
        "M13_LOWBACKPAIN": "#b2182b",
        "M13_SCIATICA": "#4d9221",
    }
    ep_markers = {
        "M13_INTERVERTEB": "o",
        "M13_LOWBACKPAIN": "s",
        "M13_SCIATICA": "^",
    }
    ep_labels = {
        "M13_INTERVERTEB": "Other intervertebral disc disorders",
        "M13_LOWBACKPAIN": "Low back pain",
        "M13_SCIATICA": "Sciatica with lumbago",
    }

    fig, ax = plt.subplots(figsize=(8.6, 6.4), constrained_layout=True)
    for idx, row in s2.iterrows():
        code = row["Endpoint_code"]
        ax.errorbar(
            row["OR"], idx,
            xerr=[[row["OR"] - row["CI_low"]], [row["CI_high"] - row["OR"]]],
            fmt=ep_markers[code], color=ep_colors[code], ecolor=ep_colors[code],
            markersize=5.6, capsize=2.3, linewidth=1.25, zorder=3,
        )
    ax.axvline(1.0, color="#555555", linestyle="--", linewidth=1.0, zorder=1)
    ax.set_xlim(0.48, 1.33)
    ax.set_yticks(range(len(s2)))
    ax.set_yticklabels([f"{g} | {e}" for g, e in zip(s2["Gene"], s2["Endpoint"])], fontsize=8.2)
    ax.set_xlabel("Odds ratio per standardized additive eQTL effect")
    ax.set_title("Historical 13-row downstream subset ($\\tau = 0.1$)", loc="left", fontsize=9.5)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.6)
    ax.set_axisbelow(True)

    handles = [
        plt.Line2D([0], [0], marker=ep_markers[c], color="none",
                   markerfacecolor=ep_colors[c], markeredgecolor=ep_colors[c],
                   label=lbl, markersize=6.2)
        for c, lbl in ep_labels.items()
    ]
    ax.legend(handles=handles, loc="upper right", frameon=True, facecolor="white",
              edgecolor="#DDDDDD", framealpha=0.92, fontsize=8.0)

    save_both(fig, "Supplementary_Figure_S2")


def make_supp_s3() -> None:
    configure_style()
    lev = pd.read_csv(CAT_DIR / "CATALOGUE_LEVERAGE_PAIRED_CHANGES.tsv", sep="\t")
    lev = lev[lev["contrast"].eq("0.001_to_0.1")].copy()
    assert len(lev) == 2438
    med = float(lev["max_leverage_change"].median())

    fig, ax = plt.subplots(figsize=(6.4, 3.8), constrained_layout=True)
    ax.hist(lev["max_leverage_change"], bins=70, color=COLORS["sky"], edgecolor="white", linewidth=0.3)
    ax.axvline(med, color=COLORS["vermillion"], linewidth=1.4, label=f"Median = {med:.3f}")
    ax.axvline(0.0, color=COLORS["grey"], linestyle="--", linewidth=1.0, label="No change (0)")
    ax.set_xlabel("Change in maximum exposure-side generalized leverage ($\\tau$: 0.001 to 0.1)")
    ax.set_ylabel("Genes")
    ax.set_title("Catalogue-wide paired exposure-side influence diagnostic (n = 2,438)", loc="left", fontsize=9.0)
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="y", color="#E2E2E2", linewidth=0.45)
    ax.set_axisbelow(True)

    save_both(fig, "Supplementary_Figure_S3")


def sha_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def generate(name: str, files: dict, output: Path) -> None:
    """Run one archived figure function with caller-supplied, validated inputs."""
    global CAT_DIR, FIG2_DATA, S2_DATA, OUT_DIR
    OUT_DIR = output
    if name == "figure-1":
        expected = {"effective_k":"CATALOGUE_EFFECTIVE_K.tsv", "exposure":"CATALOGUE_EXPOSURE_PAIRED_CHANGES.tsv", "precision":"CATALOGUE_PRECISION_PAIRED_CHANGES.tsv"}
        assert all(files[k].name == n for k, n in expected.items()), "Figure 1 input names must match the documented schema"
        assert len({files[k].parent for k in expected}) == 1, "Figure 1 inputs must share a directory"
        CAT_DIR = files["effective_k"].parent
        make_figure1()
    elif name == "figure-2":
        FIG2_DATA = files["figure2"]
        make_figure2()
    elif name == "figure-s1":
        make_supp_s1()
    elif name == "figure-s2":
        S2_DATA = files["figure_s2"]
        make_supp_s2()
    elif name == "figure-s3":
        assert files["leverage"].name == "CATALOGUE_LEVERAGE_PAIRED_CHANGES.tsv"
        CAT_DIR = files["leverage"].parent
        make_supp_s3()
    else:
        raise ValueError("Unknown figure")
