"""Validate and render the public MASTER 1.4 evidence matrix and neutral forest."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ENDPOINT_COLORS = {
    "Other intervertebral disc disorders": "#2166ac",
    "Low back pain": "#b2182b",
    "Sciatica+with lumbago": "#4d9221",
}
ENDPOINT_MARKERS = {
    "Other intervertebral disc disorders": "o",
    "Low back pain": "s",
    "Sciatica+with lumbago": "^",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def validate_table(frame: pd.DataFrame) -> None:
    required = {
        "Gene", "Endpoint", "K", "OR", "CI_low", "CI_high", "Raw_MR_P",
        "Global_BH_q_7284_row_family", "Generalized_Q_P",
        "Oriented_GLS_Egger_intercept_P", "LOO_status",
        "PP3_registered_p12_1e-5", "PP4_registered_p12_1e-5",
        "Locus_attribution_status",
    }
    missing = required - set(frame.columns)
    if len(frame) != 13 or missing:
        raise ValueError(f"invalid evidence matrix: rows={len(frame)}, missing={sorted(missing)}")
    prohibited = {column for column in frame.columns if "tier" in column.lower() or "grade" in column.lower()}
    if prohibited:
        raise ValueError(f"current evidence matrix must not encode tiers/grades: {sorted(prohibited)}")


def main() -> int:
    args = parse_args()
    root = args.repo_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = root / "results" / "master_1_4" / "table3_evidence_matrix.tsv"
    expected_forest = root / "results" / "master_1_4" / "supplementary_figure_s2_source.tsv"
    table = pd.read_csv(source, sep="\t", dtype=str)
    validate_table(table)

    forest = table[["Gene", "Endpoint", "OR", "CI_low", "CI_high"]].copy()
    forest.insert(1, "Source_gene", forest["Gene"].str.extract(r"^([^ ]+)", expand=False))
    forest["Source_gene"] = forest["Source_gene"].replace({"UTP11": "UTP11L"})
    code = {
        "Other intervertebral disc disorders": "M13_INTERVERTEB",
        "Low back pain": "M13_LOWBACKPAIN",
        "Sciatica+with lumbago": "M13_SCIATICA",
    }
    forest.insert(3, "Endpoint_code", forest["Endpoint"].map(code))
    expected = pd.read_csv(expected_forest, sep="\t", dtype=str)
    pd.testing.assert_frame_equal(forest.reset_index(drop=True), expected.reset_index(drop=True))

    shutil.copyfile(source, output / "Table_3_evidence_matrix.tsv")
    forest.to_csv(output / "Supplementary_Figure_S2_Source_Data.tsv", sep="\t", index=False, lineterminator="\n")

    plot = forest.iloc[::-1].reset_index(drop=True)
    for column in ("OR", "CI_low", "CI_high"):
        plot[column] = plot[column].astype(float)
    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    for index, row in plot.iterrows():
        color = ENDPOINT_COLORS[row.Endpoint]
        ax.errorbar(
            row.OR, index,
            xerr=[[row.OR - row.CI_low], [row.CI_high - row.OR]],
            fmt=ENDPOINT_MARKERS[row.Endpoint], color=color, ecolor=color,
            markersize=5.5, capsize=2, linewidth=1.2,
        )
    ax.axvline(1.0, color="#555555", linestyle="--", linewidth=1)
    ax.set_xlim(0.48, 1.32)
    ax.set_yticks(range(len(plot)))
    ax.set_yticklabels([f"{g} | {e}" for g, e in zip(plot.Gene, plot.Endpoint)], fontsize=8)
    ax.set_xlabel("Odds ratio per standardized additive eQTL effect")
    ax.set_title("13-row K>=2 sensitivity subset with downstream evidence")
    ax.grid(axis="x", color="#dddddd", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output / "Supplementary_Figure_S2.png", dpi=610, bbox_inches="tight")
    fig.savefig(output / "Supplementary_Figure_S2.pdf", bbox_inches="tight")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
