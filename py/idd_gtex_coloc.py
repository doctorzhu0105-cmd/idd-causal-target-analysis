#!/usr/bin/env python3
"""Retrieve free GTEx cis-eQTL records and run proxy-tissue colocalization."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from idd_genetics_utils import coloc_posteriors, wakefield_log_abf


ROOT = Path(__file__).resolve().parents[1]
API_ROOT = "https://www.ebi.ac.uk/eqtl/api/v2"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument(
        "--targets",
        type=Path,
        default=Path("artifacts/targets/priority_targets.tsv"),
    )
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument("--minimum-shared-variants", type=int, default=50)
    return parser.parse_args()


def session() -> requests.Session:
    retry = Retry(
        total=5,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
    )
    result = requests.Session()
    result.mount("https://", HTTPAdapter(max_retries=retry))
    result.headers.update({"User-Agent": "op-bioinfo-idd-gtex-coloc/1.0"})
    return result


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fetch_dataset_metadata(client: requests.Session, dataset_id: str) -> dict[str, Any]:
    response = client.get(f"{API_ROOT}/datasets/{dataset_id}", timeout=(20, 90))
    response.raise_for_status()
    return response.json()


def fetch_gene_associations(
    client: requests.Session,
    dataset_id: str,
    gene_id: str,
    page_size: int,
) -> tuple[list[dict[str, Any]], str]:
    rows = []
    start = 0
    while True:
        response = client.get(
            f"{API_ROOT}/datasets/{dataset_id}/associations",
            params={"gene_id": gene_id, "start": start, "size": page_size},
            timeout=(20, 120),
        )
        if response.status_code == 400 and response.json().get("message") == "No results":
            return [], response.url
        response.raise_for_status()
        page = response.json()
        if not isinstance(page, list):
            raise ValueError(f"Unexpected eQTL Catalogue response for {dataset_id}/{gene_id}")
        rows.extend(page)
        if len(page) < page_size:
            return rows, response.url
        start += len(page)
        if start > 200_000:
            raise ValueError(f"Unexpectedly large association result for {dataset_id}/{gene_id}")


def normalize_associations(
    rows: Iterable[Mapping[str, Any]],
    *,
    dataset_id: str,
    tissue: str,
    gene_id: str,
    gene_symbol: str,
) -> list[dict[str, object]]:
    normalized = []
    for row in rows:
        rsid = str(row.get("rsid") or "").strip()
        beta = optional_float(row.get("beta"))
        se = optional_float(row.get("se"))
        pvalue = optional_float(row.get("pvalue"))
        if not rsid.startswith("rs") or beta is None or se is None or se <= 0 or pvalue is None:
            continue
        normalized.append(
            {
                "dataset_id": dataset_id,
                "tissue": tissue,
                "gene_symbol": gene_symbol,
                "ensembl_gene_id": gene_id,
                "rsid": rsid,
                "variant": row.get("variant", ""),
                "chromosome": row.get("chromosome", ""),
                "position": row.get("position", ""),
                "ref": row.get("ref", ""),
                "alt": row.get("alt", ""),
                "beta": beta,
                "se": se,
                "pvalue": pvalue,
                "maf": value_or_blank(optional_float(row.get("maf"))),
            }
        )
    normalized.sort(key=lambda row: (str(row["chromosome"]), int(row["position"]), str(row["rsid"])))
    return normalized


def extract_finngen_variants(
    source: Path,
    requested_rsids: set[str],
    destination: Path,
    outcome_id: str,
) -> list[dict[str, object]]:
    rows = []
    with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for record in reader:
            matched = split_rsids(record["rsids"]) & requested_rsids
            for rsid in sorted(matched):
                rows.append(
                    {
                        "outcome_id": outcome_id,
                        "rsid": rsid,
                        "chromosome": record["#chrom"],
                        "position": record["pos"],
                        "ref": record["ref"],
                        "alt": record["alt"],
                        "beta": record["beta"],
                        "se": record["sebeta"],
                        "pvalue": record["pval"],
                        "eaf": record["af_alt"],
                    }
                )
    rows.sort(key=lambda row: (str(row["chromosome"]), int(row["position"]), str(row["rsid"])))
    write_tsv_gzip(destination, rows)
    return rows


def best_by_rsid(rows: Iterable[Mapping[str, object]]) -> dict[str, Mapping[str, object]]:
    result = {}
    for row in rows:
        rsid = str(row["rsid"])
        if rsid not in result or float(row["pvalue"]) < float(result[rsid]["pvalue"]):
            result[rsid] = row
    return result


def run_coloc(
    eqtl_rows: Iterable[Mapping[str, object]],
    outcome_rows: Iterable[Mapping[str, object]],
    minimum_shared: int,
) -> dict[str, object]:
    eqtl = best_by_rsid(eqtl_rows)
    outcome = best_by_rsid(outcome_rows)
    shared = sorted(eqtl.keys() & outcome.keys())
    base: dict[str, object] = {
        "shared_variants": len(shared),
        "pp0": "",
        "pp1": "",
        "pp2": "",
        "pp3": "",
        "pp4": "",
        "status": "not_estimable",
        "reason": "insufficient_shared_variants",
    }
    if len(shared) < minimum_shared:
        return base
    eqtl_abf = [
        wakefield_log_abf(float(eqtl[rsid]["beta"]), float(eqtl[rsid]["se"]) ** 2, 0.15**2)
        for rsid in shared
    ]
    outcome_abf = [
        wakefield_log_abf(float(outcome[rsid]["beta"]), float(outcome[rsid]["se"]) ** 2, 0.20**2)
        for rsid in shared
    ]
    base.update(coloc_posteriors(eqtl_abf, outcome_abf))
    base.update({"status": "estimable", "reason": ""})
    return base


def main() -> int:
    args = parse_args()
    root = args.repo_root.resolve()
    target_path = args.targets if args.targets.is_absolute() else root / args.targets
    targets = read_tsv(target_path)
    if not targets:
        raise ValueError("Priority target file is empty")
    source_config = json.loads((root / "config" / "mr_sources.json").read_text(encoding="utf-8"))
    gtex_sources = [
        source
        for source in source_config["sources"]
        if source.get("role") == "proxy_tissue_colocalization"
    ]
    client = session()
    processed_dir = root / "data" / "processed" / "idd" / "genetics" / "gtex"
    association_rows: dict[tuple[str, str], list[dict[str, object]]] = {}
    summary_rows = []
    requested_rsids: set[str] = set()

    for source in gtex_sources:
        dataset_id = source["dataset_id"]
        metadata = fetch_dataset_metadata(client, dataset_id)
        tissue = str(metadata["tissue_label"])
        if int(metadata["sample_size"]) != int(source["sample_size"]):
            raise ValueError(f"eQTL Catalogue sample-size mismatch for {dataset_id}")
        for target in targets:
            gene_id = target["ensembl_gene_id"].split(".")[0]
            gene_symbol = target["gene_symbol"]
            raw_rows, final_url = fetch_gene_associations(
                client,
                dataset_id,
                gene_id,
                args.page_size,
            )
            rows = normalize_associations(
                raw_rows,
                dataset_id=dataset_id,
                tissue=tissue,
                gene_id=gene_id,
                gene_symbol=gene_symbol,
            )
            association_rows[(dataset_id, gene_symbol)] = rows
            requested_rsids.update(str(row["rsid"]) for row in rows)
            path = processed_dir / f"{dataset_id}_{gene_symbol}.tsv.gz"
            write_tsv_gzip(path, rows)
            summary_rows.append(
                {
                    "dataset_id": dataset_id,
                    "study_id": metadata["study_id"],
                    "study_label": metadata["study_label"],
                    "tissue": tissue,
                    "sample_size": metadata["sample_size"],
                    "gene_symbol": gene_symbol,
                    "ensembl_gene_id": gene_id,
                    "association_rows": len(rows),
                    "status": "retrieved" if rows else "gene_not_available_in_dataset",
                    "api_url": final_url,
                    "local_path": path.relative_to(root).as_posix(),
                }
            )

    outcome_rows = {}
    outcome_config = json.loads((root / "config" / "idd_sources.json").read_text(encoding="utf-8"))["gwas"]
    for source in outcome_config:
        outcome_id = source["outcome_id"].removeprefix("FinnGen_R13_")
        path = root / "data" / "raw" / "idd" / "genetics" / "finngen" / f"finngen_R13_{outcome_id}.gz"
        destination = processed_dir / f"{outcome_id}_gtex_variants.tsv.gz"
        outcome_rows[outcome_id] = extract_finngen_variants(path, requested_rsids, destination, outcome_id)

    coloc_rows = []
    for source in gtex_sources:
        dataset_id = source["dataset_id"]
        metadata = fetch_dataset_metadata(client, dataset_id)
        for target in targets:
            symbol = target["gene_symbol"]
            for outcome_id, rows in outcome_rows.items():
                result = run_coloc(
                    association_rows[(dataset_id, symbol)],
                    rows,
                    args.minimum_shared_variants,
                )
                coloc_rows.append(
                    {
                        "gene_symbol": symbol,
                        "ensembl_gene_id": target["ensembl_gene_id"],
                        "dataset_id": dataset_id,
                        "tissue": metadata["tissue_label"],
                        "sample_size": metadata["sample_size"],
                        "outcome_id": outcome_id,
                        **result,
                    }
                )

    artifact_dir = root / "artifacts" / "mr"
    write_tsv(artifact_dir / "gtex_eqtl_summary.tsv", summary_rows)
    write_tsv(artifact_dir / "gtex_coloc_results.tsv", coloc_rows)
    summary = {
        "target_count": len(targets),
        "dataset_count": len(gtex_sources),
        "requested_unique_rsids": len(requested_rsids),
        "coloc_estimable": sum(row["status"] == "estimable" for row in coloc_rows),
        "coloc_pp4_at_least_0_8": sum(
            optional_float(row["pp4"]) is not None and float(row["pp4"]) >= 0.8
            for row in coloc_rows
        ),
        "source_api": API_ROOT,
    }
    (artifact_dir / "gtex_coloc_summary.json").write_text(
        json.dumps(summary, ensure_ascii=True, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


def split_rsids(value: str) -> set[str]:
    return {
        token
        for token in value.replace(";", ",").replace("|", ",").split(",")
        if token.startswith("rs")
    }


def optional_float(value: object) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math_isfinite(parsed) else None


def math_isfinite(value: float) -> bool:
    return value == value and value not in (float("inf"), float("-inf"))


def value_or_blank(value: object) -> object:
    return "" if value is None else value


def write_tsv(path: Path, rows: list[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = collect_fields(rows)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def write_tsv_gzip(path: Path, rows: list[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = collect_fields(rows)
    temporary = path.with_name(path.name + ".tmp")
    with gzip.GzipFile(filename=str(temporary), mode="wb", mtime=0) as compressed:
        with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    temporary.replace(path)


def collect_fields(rows: Iterable[Mapping[str, object]]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    return fields


if __name__ == "__main__":
    raise SystemExit(main())
