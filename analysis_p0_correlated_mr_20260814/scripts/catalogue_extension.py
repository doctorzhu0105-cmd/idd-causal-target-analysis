#!/usr/bin/env python3
"""Validate and regenerate redistributable MASTER 1.5 catalogue derivatives."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


CONTRASTS = ((0.001, 0.1, "0.001_to_0.1"), (0.001, 0.2, "0.001_to_0.2"))


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def exposure_pairs(effective: pd.DataFrame, joint: pd.DataFrame) -> pd.DataFrame:
    source = effective.merge(
        joint[["gene", "ensembl_gene_id", "threshold", "joint_exposure_chi_square"]],
        on=["gene", "ensembl_gene_id", "threshold"], validate="one_to_one",
    )
    source = source[source.status.eq("estimable")]
    rows = []
    for low, high, label in CONTRASTS:
        left = source[np.isclose(source.threshold, low)]
        right = source[np.isclose(source.threshold, high)]
        pair = left.merge(right, on=["gene", "ensembl_gene_id"], suffixes=("_from", "_to"), validate="one_to_one")
        pair["contrast"] = label
        pair["K_ratio"] = pair.K_to / pair.K_from
        pair["K_absolute_change"] = pair.K_to - pair.K_from
        pair["K_percent_change"] = 100 * (pair.K_ratio - 1)
        pair["effective_K_ratio"] = pair.effective_K_to / pair.effective_K_from
        pair["effective_K_absolute_change"] = pair.effective_K_to - pair.effective_K_from
        pair["effective_K_percent_change"] = 100 * (pair.effective_K_ratio - 1)
        pair["joint_chi_square_ratio"] = pair.joint_exposure_chi_square_to / pair.joint_exposure_chi_square_from
        pair["joint_chi_square_percent_change"] = 100 * (pair.joint_chi_square_ratio - 1)
        pair["effective_K_over_K_change"] = pair.effective_K_over_K_to - pair.effective_K_over_K_from
        rows.append(pair[[
            "gene", "ensembl_gene_id", "contrast", "threshold_from", "threshold_to",
            "K_from", "K_to", "K_ratio", "K_absolute_change", "K_percent_change",
            "effective_K_from", "effective_K_to", "effective_K_ratio", "effective_K_absolute_change",
            "effective_K_percent_change", "joint_exposure_chi_square_from", "joint_exposure_chi_square_to",
            "joint_chi_square_ratio", "joint_chi_square_percent_change", "effective_K_over_K_from",
            "effective_K_over_K_to", "effective_K_over_K_change",
        ]])
    return pd.concat(rows, ignore_index=True)


def precision_pairs(givw: pd.DataFrame) -> pd.DataFrame:
    source = givw[givw.status.eq("estimable")]
    rows = []
    for low, high, label in CONTRASTS:
        left = source[np.isclose(source.threshold, low)]
        right = source[np.isclose(source.threshold, high)]
        pair = left.merge(right, on=["gene", "ensembl_gene_id", "endpoint"], suffixes=("_from", "_to"), validate="one_to_one")
        pair["contrast"] = label
        pair["K_change"] = pair.K_to - pair.K_from
        pair["K_ratio"] = pair.K_to / pair.K_from
        pair["K_percent_change"] = 100 * (pair.K_ratio - 1)
        pair["effective_K_change"] = pair.effective_K_to - pair.effective_K_from
        pair["effective_K_ratio"] = pair.effective_K_to / pair.effective_K_from
        pair["effective_K_percent_change"] = 100 * (pair.effective_K_ratio - 1)
        pair["joint_chi_square_change"] = pair.joint_exposure_chi_square_to - pair.joint_exposure_chi_square_from
        pair["joint_chi_square_ratio"] = pair.joint_exposure_chi_square_to / pair.joint_exposure_chi_square_from
        pair["joint_chi_square_percent_change"] = 100 * (pair.joint_chi_square_ratio - 1)
        pair["SE_ratio"] = pair.SE_to / pair.SE_from
        pair["SE_percent_change"] = 100 * (pair.SE_ratio - 1)
        pair["precision_gain_percent"] = -pair.SE_percent_change
        pair["beta_change"] = pair.beta_to - pair.beta_from
        pair["direction_consistency"] = np.sign(pair.beta_from) == np.sign(pair.beta_to)
        pair["P_change"] = pair.P_to - pair.P_from
        rows.append(pair[[
            "gene", "ensembl_gene_id", "endpoint", "contrast", "threshold_from", "threshold_to",
            "K_from", "K_to", "K_change", "K_ratio", "K_percent_change", "effective_K_from", "effective_K_to",
            "effective_K_change", "effective_K_ratio", "effective_K_percent_change",
            "joint_exposure_chi_square_from", "joint_exposure_chi_square_to", "joint_chi_square_change",
            "joint_chi_square_ratio", "joint_chi_square_percent_change", "SE_from", "SE_to", "SE_ratio",
            "SE_percent_change", "precision_gain_percent", "beta_from", "beta_to", "beta_change",
            "direction_consistency", "P_from", "P_to", "P_change",
        ]])
    return pd.concat(rows, ignore_index=True)


def distribution_summary(frame: pd.DataFrame, groups: list[str], metrics: list[str]) -> pd.DataFrame:
    rows = []
    for key, group in frame.groupby(groups, dropna=False):
        keys = key if isinstance(key, tuple) else (key,)
        identity = dict(zip(groups, keys))
        for metric in metrics:
            values = pd.to_numeric(group[metric], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            rows.append(identity | {"metric": metric, "N": len(values), "median": values.median(),
                                    "p05": values.quantile(.05), "p25": values.quantile(.25),
                                    "p75": values.quantile(.75), "p95": values.quantile(.95),
                                    "IQR": values.quantile(.75) - values.quantile(.25)})
    return pd.DataFrame(rows)


def matched_associations(precision: pd.DataFrame) -> pd.DataFrame:
    primary = precision[precision.contrast.eq("0.001_to_0.1")]
    specs = {
        "nominal_K": {"absolute_gain": primary.K_change, "relative_percent_gain": primary.K_percent_change,
                      "fold_gain": primary.K_ratio, "log_ratio": np.log(primary.K_ratio)},
        "effective_K": {"absolute_gain": primary.effective_K_change,
                        "relative_percent_gain": primary.effective_K_percent_change,
                        "fold_gain": primary.effective_K_ratio, "log_ratio": np.log(primary.effective_K_ratio)},
        "joint_exposure_chi_square": {"absolute_gain": primary.joint_chi_square_change,
                                      "relative_percent_gain": primary.joint_chi_square_percent_change,
                                      "fold_gain": primary.joint_chi_square_ratio,
                                      "log_ratio": np.log(primary.joint_chi_square_ratio)},
    }
    rows = []
    for predictor, parameters in specs.items():
        for parameterization, values in parameters.items():
            valid = np.isfinite(values.to_numpy(float)) & np.isfinite(primary.precision_gain_percent.to_numpy(float))
            result = stats.spearmanr(values.to_numpy(float)[valid], primary.precision_gain_percent.to_numpy(float)[valid])
            rows.append({"contrast": "r2_0.001_to_0.1", "predictor": predictor,
                         "parameterization": parameterization,
                         "precision_outcome": "precision_gain_percent=-SE_percent_change", "N": int(valid.sum()),
                         "spearman_rho": float(result.statistic), "P": float(result.pvalue),
                         "predictor_median": float(np.median(values.to_numpy(float)[valid])),
                         "precision_gain_median": float(np.median(primary.precision_gain_percent.to_numpy(float)[valid])),
                         "status": "MATCHED_PARAMETERIZATION"})
    return pd.DataFrame(rows)


def endpoint_associations(precision: pd.DataFrame) -> pd.DataFrame:
    primary = precision[precision.contrast.eq("0.001_to_0.1")]
    labels = {"M13_INTERVERTEB": "Other intervertebral disc disorders",
              "M13_LOWBACKPAIN": "Low back pain", "M13_SCIATICA": "Sciatica with lumbago"}
    rows = []
    for endpoint, frame in primary.groupby("endpoint", sort=True):
        nominal = stats.spearmanr(frame.K_percent_change, frame.precision_gain_percent).statistic
        effective = stats.spearmanr(frame.effective_K_percent_change, frame.precision_gain_percent).statistic
        joint = stats.spearmanr(frame.joint_chi_square_percent_change, frame.precision_gain_percent).statistic
        rows.append({"endpoint": endpoint, "endpoint_label": labels[endpoint], "N": len(frame),
                     "nominal_K_relative_rho": nominal, "effective_K_relative_rho": effective,
                     "joint_information_relative_rho": joint, "joint_gt_nominal": bool(joint > nominal)})
    return pd.DataFrame(rows)


def regenerate(results: Path, output: Path) -> None:
    effective = read_tsv(results / "CATALOGUE_EFFECTIVE_K.tsv")
    joint = read_tsv(results / "CATALOGUE_JOINT_EXPOSURE_INFORMATION.tsv")
    givw = read_tsv(results / "CATALOGUE_GIVW_LD_LADDER.tsv")
    exposure = exposure_pairs(effective, joint)
    precision = precision_pairs(givw)
    write_tsv(exposure, output / "CATALOGUE_EXPOSURE_PAIRED_CHANGES.tsv")
    write_tsv(precision, output / "CATALOGUE_PRECISION_PAIRED_CHANGES.tsv")
    write_tsv(distribution_summary(exposure, ["contrast"], [
        "K_ratio", "K_absolute_change", "K_percent_change", "effective_K_ratio",
        "effective_K_absolute_change", "effective_K_percent_change", "joint_chi_square_ratio",
        "joint_chi_square_percent_change", "effective_K_over_K_change",
    ]), output / "CATALOGUE_EXPOSURE_PAIRED_CHANGE_SUMMARY.tsv")
    all_precision = precision.copy()
    all_precision["endpoint"] = "ALL"
    metrics = ["K_change", "K_ratio", "effective_K_change", "joint_chi_square_percent_change",
               "SE_ratio", "SE_percent_change", "beta_change", "P_change"]
    summary = pd.concat([distribution_summary(precision, ["contrast", "endpoint"], metrics),
                         distribution_summary(all_precision, ["contrast", "endpoint"], metrics)], ignore_index=True)
    write_tsv(summary, output / "CATALOGUE_PRECISION_PAIRED_CHANGE_SUMMARY.tsv")
    write_tsv(matched_associations(precision), output / "INFORMATION_PRECISION_MATCHED_ASSOCIATION.tsv")
    write_tsv(endpoint_associations(precision), output / "ENDPOINT_STRATIFIED_INFORMATION_PRECISION.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    regenerate(args.results_dir.resolve(), args.output_dir.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
