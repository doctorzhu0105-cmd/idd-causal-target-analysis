#!/usr/bin/env python3
"""Localize registered IDD targets in broad human disc cell populations."""

from __future__ import annotations

import argparse
import gc
import gzip
import hashlib
import json
import tarfile
from pathlib import Path

import anndata as ad
import h5py
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from scipy import sparse

matplotlib.use("Agg")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument(
        "--targets",
        type=Path,
        default=Path("artifacts/bulk/druggable_bulk_candidates.tsv"),
    )
    parser.add_argument("--reuse-qc", action="store_true")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_soft(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)
                current = {"sample_id": line.split(" = ", 1)[1].strip()}
                continue
            if current is None:
                continue
            if line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 = "):
                current["tissue"] = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1].strip()
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key.lower().replace(" ", "_")] = item
        if current is not None:
            records.append(current)
    frame = pd.DataFrame.from_records(records)
    expected = {"sample_id", "tissue", "disease_state", "patient_id", "age", "sex"}
    missing = expected.difference(frame.columns)
    if missing:
        raise ValueError(f"GSE230809 SOFT is missing fields: {sorted(missing)}")
    if len(frame) != 24 or frame["sample_id"].nunique() != 24:
        raise ValueError(f"Expected 24 unique GSE230809 samples, observed {len(frame)}")
    frame["age"] = pd.to_numeric(frame["age"], errors="raise")
    frame["tissue"] = frame["tissue"].map(
        {"annulus fibrosus": "AF", "nucleus pulposus": "NP"}
    )
    if frame["tissue"].isna().any():
        raise ValueError("Unrecognized tissue in GSE230809 SOFT")
    return frame


def extract_tar(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    if len(list(destination.glob("*_matrix.mtx.gz"))) == 24:
        return
    root = destination.resolve()
    with tarfile.open(archive, "r") as handle:
        for member in handle.getmembers():
            target = (destination / member.name).resolve()
            if root not in target.parents:
                raise ValueError(f"Unsafe archive member: {member.name}")
        handle.extractall(destination, filter="data")
    matrices = list(destination.glob("*_matrix.mtx.gz"))
    if len(matrices) != 24:
        raise ValueError(f"Expected 24 matrix files, observed {len(matrices)}")


def read_10x_mtx_streaming(raw_dir: Path, prefix: str) -> ad.AnnData:
    """Read a column-sorted 10x Matrix Market file without a large COO copy."""
    matrix_path = raw_dir / f"{prefix}matrix.mtx.gz"
    features_path = raw_dir / f"{prefix}features.tsv.gz"
    barcodes_path = raw_dir / f"{prefix}barcodes.tsv.gz"
    for path in (matrix_path, features_path, barcodes_path):
        if not path.exists():
            raise FileNotFoundError(path)

    with gzip.open(matrix_path, "rb") as handle:
        banner = handle.readline().strip()
        if banner != b"%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"Unsupported Matrix Market banner in {matrix_path.name}: {banner!r}")
        dimensions = handle.readline()
        while dimensions.startswith(b"%"):
            dimensions = handle.readline()
        n_features, n_barcodes, expected_nnz = map(int, dimensions.split())
        values = np.empty(expected_nnz, dtype=np.float32)
        feature_indices = np.empty(expected_nnz, dtype=np.int32)
        column_counts = np.zeros(n_barcodes, dtype=np.int64)
        offset = 0
        remainder = b""
        previous_column = -1

        while block := handle.read(16 * 1024 * 1024):
            block = remainder + block
            final_newline = block.rfind(b"\n")
            if final_newline < 0:
                remainder = block
                continue
            complete = block[: final_newline + 1]
            remainder = block[final_newline + 1 :]
            parsed = np.fromstring(complete, dtype=np.int64, sep=" ")
            if parsed.size % 3:
                raise ValueError(f"Malformed Matrix Market coordinates in {matrix_path.name}")
            triples = parsed.reshape(-1, 3)
            rows = triples[:, 0] - 1
            columns = triples[:, 1] - 1
            count = len(triples)
            if offset + count > expected_nnz:
                raise ValueError(f"Matrix Market nnz exceeds header in {matrix_path.name}")
            if count:
                if rows.min() < 0 or rows.max() >= n_features:
                    raise ValueError(f"Feature index outside declared dimensions in {matrix_path.name}")
                if columns.min() < 0 or columns.max() >= n_barcodes:
                    raise ValueError(f"Barcode index outside declared dimensions in {matrix_path.name}")
                if columns[0] < previous_column or np.any(columns[1:] < columns[:-1]):
                    raise ValueError(f"Matrix Market entries are not column sorted in {matrix_path.name}")
                previous_column = int(columns[-1])
            feature_indices[offset : offset + count] = rows
            values[offset : offset + count] = triples[:, 2]
            column_counts += np.bincount(columns, minlength=n_barcodes)
            offset += count

        if remainder.strip():
            parsed = np.fromstring(remainder, dtype=np.int64, sep=" ")
            if parsed.size != 3:
                raise ValueError(f"Malformed final Matrix Market row in {matrix_path.name}")
            row, column, value = parsed
            row -= 1
            column -= 1
            if not (0 <= row < n_features and 0 <= column < n_barcodes):
                raise ValueError(f"Final Matrix Market index outside dimensions in {matrix_path.name}")
            if column < previous_column:
                raise ValueError(f"Matrix Market entries are not column sorted in {matrix_path.name}")
            feature_indices[offset] = row
            values[offset] = value
            column_counts[column] += 1
            offset += 1

    if offset != expected_nnz:
        raise ValueError(
            f"Matrix Market nnz mismatch in {matrix_path.name}: expected {expected_nnz}, observed {offset}"
        )
    indptr = np.empty(n_barcodes + 1, dtype=np.int64)
    indptr[0] = 0
    np.cumsum(column_counts, out=indptr[1:])
    matrix = sparse.csc_matrix(
        (values, feature_indices, indptr),
        shape=(n_features, n_barcodes),
    ).transpose(copy=False)

    features = pd.read_csv(features_path, sep="\t", header=None, dtype=str)
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, dtype=str)
    if len(features) != n_features or len(barcodes) != n_barcodes:
        raise ValueError(
            f"10x dimension mismatch for {prefix}: matrix={n_features}x{n_barcodes}, "
            f"features={len(features)}, barcodes={len(barcodes)}"
        )
    var = pd.DataFrame(index=pd.Index(features.iloc[:, 1].astype(str), name="gene_symbol"))
    var["gene_id"] = features.iloc[:, 0].to_numpy()
    if features.shape[1] >= 3:
        var["feature_type"] = features.iloc[:, 2].to_numpy()
    obs = pd.DataFrame(index=pd.Index(barcodes.iloc[:, 0].astype(str), name="barcode"))
    return ad.AnnData(X=matrix, obs=obs, var=var)


def load_sample(
    raw_dir: Path,
    record: pd.Series,
    qc: dict[str, float],
) -> tuple[ad.AnnData, dict[str, int | float | str]]:
    sample_id = record["sample_id"]
    matches = list(raw_dir.glob(f"{sample_id}_*_matrix.mtx.gz"))
    if len(matches) != 1:
        raise ValueError(f"Expected one matrix for {sample_id}, observed {len(matches)}")
    prefix = matches[0].name[: -len("matrix.mtx.gz")]
    data = read_10x_mtx_streaming(raw_dir, prefix)
    data.var_names_make_unique()
    data.obs_names = [f"{sample_id}_{barcode}" for barcode in data.obs_names]
    data.var["mt"] = data.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(data, qc_vars=["mt"], inplace=True, percent_top=None)
    before = data.n_obs
    keep = (
        (data.obs["n_genes_by_counts"] >= int(qc["minimum_genes"]))
        & (data.obs["total_counts"] >= int(qc["minimum_counts"]))
        & (data.obs["n_genes_by_counts"] <= int(qc["maximum_genes"]))
        & (data.obs["pct_counts_mt"] <= float(qc["maximum_mitochondrial_percent"]))
    )
    data = data[keep].copy()
    for field in ("sample_id", "tissue", "disease_state", "patient_id", "sex"):
        data.obs[field] = str(record[field])
    data.obs["age"] = int(record["age"])
    summary = {
        "sample_id": sample_id,
        "tissue": record["tissue"],
        "disease_state": record["disease_state"],
        "patient_id": record["patient_id"],
        "age": int(record["age"]),
        "cells_before_qc": before,
        "cells_after_qc": data.n_obs,
    }
    return data, summary


def marker_fraction(data: ad.AnnData, genes: list[str]) -> np.ndarray:
    available = [gene for gene in genes if gene in data.var_names]
    if not available:
        return np.zeros(data.n_obs, dtype=float)
    matrix = data[:, available].X
    if sparse.issparse(matrix):
        return np.asarray((matrix > 0).mean(axis=1)).ravel()
    return (matrix > 0).mean(axis=1)


def annotate_clusters(data: ad.AnnData, marker_sets: dict[str, list[str]], seed: int) -> pd.DataFrame:
    score_columns = []
    for label, genes in marker_sets.items():
        available = [gene for gene in genes if gene in data.var_names]
        if len(available) < 2:
            raise ValueError(f"Too few marker genes available for {label}: {available}")
        score_name = f"score_{label}"
        sc.tl.score_genes(data, available, score_name=score_name, random_state=seed)
        data.obs[f"fraction_{label}"] = marker_fraction(data, available)
        score_columns.append(score_name)
    rows = []
    for cluster, frame in data.obs.groupby("leiden", observed=True):
        mean_scores = frame[score_columns].mean()
        best = mean_scores.idxmax().removeprefix("score_")
        ordered = mean_scores.sort_values(ascending=False)
        margin = float(ordered.iloc[0] - ordered.iloc[1]) if len(ordered) > 1 else np.nan
        rows.append(
            {
                "cluster": str(cluster),
                "broad_cell_type": best,
                "annotation_margin": margin,
                "annotation_confidence": "higher" if margin >= 0.1 else "lower",
                "cells": len(frame),
                "af_fraction": float((frame["tissue"] == "AF").mean()),
                "np_fraction": float((frame["tissue"] == "NP").mean()),
                **{column: float(mean_scores[column]) for column in score_columns},
            }
        )
    result = pd.DataFrame(rows).sort_values("cluster")
    mapping = result.set_index("cluster")["broad_cell_type"]
    data.obs["broad_cell_type"] = data.obs["leiden"].astype(str).map(mapping)
    return result


def select_clustering_cells(data: ad.AnnData, max_per_sample: int, seed: int) -> np.ndarray:
    """Select a deterministic, sample-balanced subset for expensive graph operations."""
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    sample_values = data.obs["sample_id"].astype(str).to_numpy()
    for sample_id in sorted(set(sample_values)):
        positions = np.flatnonzero(sample_values == sample_id)
        take = min(len(positions), max_per_sample)
        selected.append(np.sort(rng.choice(positions, size=take, replace=False)))
    if not selected:
        raise ValueError("No cells available for clustering subset")
    return np.sort(np.concatenate(selected))


def annotate_cells_by_markers(data: ad.AnnData, marker_sets: dict[str, list[str]], seed: int) -> None:
    """Assign broad marker-based labels to all cells after subset clustering."""
    score_columns = []
    for label, genes in marker_sets.items():
        available = [gene for gene in genes if gene in data.var_names]
        if len(available) < 2:
            raise ValueError(f"Too few marker genes available for {label}: {available}")
        score_name = f"full_score_{label}"
        sc.tl.score_genes(data, available, score_name=score_name, random_state=seed)
        score_columns.append(score_name)
    scores = data.obs[score_columns].to_numpy()
    labels = [column.removeprefix("full_score_") for column in score_columns]
    data.obs["broad_cell_type"] = [labels[index] for index in np.argmax(scores, axis=1)]


def target_localization(data: ad.AnnData, targets: list[str]) -> pd.DataFrame:
    rows = []
    group_columns = ["broad_cell_type", "tissue", "disease_state"]
    for target in targets:
        if target not in data.var_names:
            rows.append(
                {
                    "gene_symbol": target,
                    "status": "gene_not_detected_in_feature_space",
                }
            )
            continue
        values = data[:, target].X
        values = values.toarray().ravel() if sparse.issparse(values) else np.asarray(values).ravel()
        frame = data.obs[group_columns + ["patient_id", "sample_id"]].copy()
        frame["expression"] = values
        for keys, subset in frame.groupby(group_columns, observed=True):
            donor = subset.groupby("patient_id", observed=True)["expression"].agg(
                donor_mean="mean",
                donor_fraction=lambda value: float((value > 0).mean()),
            )
            rows.append(
                {
                    "gene_symbol": target,
                    "status": "localized",
                    "broad_cell_type": keys[0],
                    "tissue": keys[1],
                    "disease_state": keys[2],
                    "cells": len(subset),
                    "donors": subset["patient_id"].nunique(),
                    "samples": subset["sample_id"].nunique(),
                    "cell_fraction_expressing": float((subset["expression"] > 0).mean()),
                    "mean_log1p_expression": float(subset["expression"].mean()),
                    "mean_donor_expression": float(donor["donor_mean"].mean()),
                    "donor_detection_fraction": float((donor["donor_fraction"] >= 0.05).mean()),
                }
            )
    return pd.DataFrame(rows)


def matrix_barcode_count(raw_dir: Path, sample_id: str) -> int:
    matches = list(raw_dir.glob(f"{sample_id}_*_matrix.mtx.gz"))
    if len(matches) != 1:
        raise ValueError(f"Expected one matrix for {sample_id}, observed {len(matches)}")
    with gzip.open(matches[0], "rt", encoding="ascii") as handle:
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        _, n_barcodes, _ = map(int, line.split())
    return n_barcodes


def reused_qc_rows(sample_manifest: pd.DataFrame, raw_dir: Path, sample_qc_dir: Path) -> tuple[list[Path], list[dict[str, object]]]:
    sample_paths: list[Path] = []
    qc_rows: list[dict[str, object]] = []
    for _, record in sample_manifest.iterrows():
        sample_path = sample_qc_dir / f"{record['sample_id']}.h5ad"
        if not sample_path.exists():
            raise FileNotFoundError(f"Missing reusable QC object: {sample_path}")
        with h5py.File(sample_path, "r") as handle:
            cells_after = len(handle["obs"]["_index"])
        cells_before = matrix_barcode_count(raw_dir, str(record["sample_id"]))
        sample_paths.append(sample_path)
        qc_rows.append(
            {
                "sample_id": record["sample_id"],
                "tissue": record["tissue"],
                "disease_state": record["disease_state"],
                "patient_id": record["patient_id"],
                "age": int(record["age"]),
                "cells_before_qc": cells_before,
                "cells_after_qc": cells_after,
            }
        )
    return sample_paths, qc_rows


def main() -> None:
    args = parse_args()
    root = args.repo_root.resolve()
    target_path = args.targets if args.targets.is_absolute() else root / args.targets
    config = json.loads((root / "config" / "idd_analysis.json").read_text(encoding="utf-8"))
    single = config["single_cell"]
    seed = int(config["seed"])
    np.random.seed(seed)
    sc.settings.seed = seed
    raw_root = root / "data" / "raw" / "idd" / "GSE230809"
    archive = raw_root / "GSE230809_RAW.tar"
    soft = raw_root / "GSE230809_family.soft.gz"
    if not archive.exists() or not soft.exists():
        raise FileNotFoundError("Run py/idd_acquire_scrna.py before localization")
    sample_manifest = parse_soft(soft)
    extracted = root / "data" / "processed" / "idd" / "scrna" / "GSE230809_raw"
    extract_tar(archive, extracted)

    sample_qc_dir = root / "data" / "processed" / "idd" / "scrna" / "sample_qc"
    sample_qc_dir.mkdir(parents=True, exist_ok=True)
    combined_path = root / "data" / "processed" / "idd" / "scrna" / "gse230809_qc_combined.h5ad"
    if args.reuse_qc:
        if not combined_path.exists():
            raise FileNotFoundError(f"Missing reusable combined QC object: {combined_path}")
        sample_paths, qc_rows = reused_qc_rows(sample_manifest, extracted, sample_qc_dir)
    else:
        sample_paths = []
        qc_rows = []
        for _, record in sample_manifest.iterrows():
            data, summary = load_sample(extracted, record, single["qc"])
            sample_path = sample_qc_dir / f"{record['sample_id']}.h5ad"
            temporary_sample = sample_path.with_name(f"{sample_path.stem}.tmp.h5ad")
            data.write_h5ad(temporary_sample, compression="lzf")
            temporary_sample.replace(sample_path)
            sample_paths.append(sample_path)
            qc_rows.append(summary)
            print(
                f"processed {record['sample_id']}: {summary['cells_after_qc']}/{summary['cells_before_qc']} cells",
                flush=True,
            )
            del data
            gc.collect()
        temporary_combined = combined_path.with_name(f"{combined_path.stem}.tmp.h5ad")
        if temporary_combined.exists():
            temporary_combined.unlink()
        ad.experimental.concat_on_disk(
            sample_paths,
            temporary_combined,
            join="inner",
            merge="same",
            index_unique=None,
            max_loaded_elems=20_000_000,
        )
        temporary_combined.replace(combined_path)
    combined = ad.read_h5ad(combined_path)
    gc.collect()
    combined.obs_names_make_unique()
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    cluster_indices = select_clustering_cells(
        combined,
        int(single["clustering"]["max_cells_per_sample"]),
        seed,
    )
    cluster_data = combined[cluster_indices].copy()
    sc.pp.highly_variable_genes(
        cluster_data,
        n_top_genes=int(single["clustering"]["highly_variable_genes"]),
        flavor="seurat",
        batch_key="sample_id",
    )
    sc.pp.pca(
        cluster_data,
        n_comps=int(single["clustering"]["principal_components"]),
        mask_var="highly_variable",
        random_state=seed,
    )
    sce.pp.harmony_integrate(
        cluster_data,
        key="sample_id",
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
        random_state=seed,
        max_iter_harmony=20,
    )
    sc.pp.neighbors(
        cluster_data,
        n_neighbors=int(single["clustering"]["neighbors"]),
        use_rep="X_pca_harmony",
        random_state=seed,
    )
    sc.tl.umap(cluster_data, random_state=seed)
    sc.tl.leiden(
        cluster_data,
        resolution=float(single["clustering"]["leiden_resolution"]),
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added="leiden",
    )
    cluster_count = int(cluster_data.obs["leiden"].nunique())
    cluster_table = annotate_clusters(cluster_data, single["marker_sets"], seed)
    annotate_cells_by_markers(combined, single["marker_sets"], seed)
    target_table = pd.read_csv(target_path, sep="\t")
    targets = sorted(target_table["gene_symbol"].dropna().astype(str).unique())
    localization = target_localization(combined, targets)

    output_dir = root / "artifacts" / "scrna"
    figure_dir = root / "figures" / "idd"
    log_dir = root / "logs" / "idd"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    sample_manifest.to_csv(output_dir / "gse230809_sample_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(qc_rows).to_csv(output_dir / "gse230809_qc.tsv", sep="\t", index=False)
    cluster_table.to_csv(output_dir / "gse230809_cluster_annotations.tsv", sep="\t", index=False)
    localization.to_csv(output_dir / "target_localization.tsv", sep="\t", index=False)

    sc.settings.figdir = figure_dir
    sc.pl.umap(
        cluster_data,
        color=["broad_cell_type", "tissue"],
        show=False,
        save="_idd_cell_types.png",
    )
    del cluster_data
    gc.collect()
    summary = {
        "seed": seed,
        "dataset": "GSE230809",
        "raw_tar_bytes": archive.stat().st_size,
        "raw_tar_sha256": sha256_file(archive),
        "samples": len(sample_manifest),
        "cells_before_qc": int(sum(row["cells_before_qc"] for row in qc_rows)),
        "cells_after_qc": int(combined.n_obs),
        "clustered_cells": int(len(cluster_indices)),
        "clustering_mode": "sample_balanced_subset_for_graph_and_umap; localization_on_all_qc_cells",
        "genes": int(combined.n_vars),
        "clusters": cluster_count,
        "targets_requested": len(targets),
        "targets_detected": int(
            localization.loc[localization["status"] == "localized", "gene_symbol"].nunique()
        ),
        "disease_effect_claim_allowed": False,
        "doublet_detection": "not_run_dependency_unavailable",
    }
    (log_dir / "scrna_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
