#!/usr/bin/env python3
"""Validate registered FinnGen IDD outcomes and write deterministic metadata."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import BinaryIO


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_HEADER = [
    "#chrom",
    "pos",
    "ref",
    "alt",
    "rsids",
    "nearest_genes",
    "pval",
    "mlogp",
    "beta",
    "sebeta",
    "af_alt",
    "af_alt_cases",
    "af_alt_controls",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    return parser.parse_args()


def file_hashes(path: Path) -> tuple[str, str]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        update_hashes(handle, md5, sha256)
    return md5.hexdigest(), sha256.hexdigest()


def update_hashes(handle: BinaryIO, *digests: object) -> None:
    for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
        for digest in digests:
            digest.update(block)


def validate_bgzf(path: Path) -> tuple[list[str], list[str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        first_row = handle.readline().rstrip("\r\n").split("\t")
        while handle.read(8 * 1024 * 1024):
            pass
    if header != EXPECTED_HEADER:
        raise ValueError(f"Unexpected FinnGen header in {path}: {header}")
    if len(first_row) != len(EXPECTED_HEADER):
        raise ValueError(f"Unexpected first-row width in {path}: {len(first_row)}")
    return header, first_row


def main() -> int:
    args = parse_args()
    root = args.repo_root.resolve()
    config = json.loads((root / "config" / "idd_sources.json").read_text(encoding="utf-8"))
    raw_dir = root / "data" / "raw" / "idd" / "genetics" / "finngen"
    manifest_path = raw_dir / "finngen_R13_manifest.tsv"
    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        manifest_rows = {
            row["phenocode"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    manifest_md5, manifest_sha256 = file_hashes(manifest_path)
    outcomes = []
    for source in config["gwas"]:
        phenocode = source["outcome_id"].removeprefix("FinnGen_R13_")
        manifest = manifest_rows.get(phenocode)
        if manifest is None:
            raise ValueError(f"FinnGen manifest lacks {phenocode}")
        if int(manifest["num_cases"]) != int(source["cases"]):
            raise ValueError(f"FinnGen case count mismatch for {phenocode}")
        if int(manifest["num_controls"]) != int(source["controls"]):
            raise ValueError(f"FinnGen control count mismatch for {phenocode}")
        if manifest["path_https"] != source["url"]:
            raise ValueError(f"FinnGen URL mismatch for {phenocode}")

        summary_path = raw_dir / f"finngen_R13_{phenocode}.gz"
        index_path = raw_dir / f"finngen_R13_{phenocode}.gz.tbi"
        if summary_path.stat().st_size != int(source["expected_content_length"]):
            raise ValueError(f"FinnGen summary length mismatch for {phenocode}")
        if index_path.stat().st_size != int(source["index_expected_content_length"]):
            raise ValueError(f"FinnGen index length mismatch for {phenocode}")
        summary_md5, summary_sha256 = file_hashes(summary_path)
        index_md5, index_sha256 = file_hashes(index_path)
        if summary_md5 != source["expected_md5"]:
            raise ValueError(f"FinnGen summary MD5 mismatch for {phenocode}")
        if index_md5 != source["index_expected_md5"]:
            raise ValueError(f"FinnGen index MD5 mismatch for {phenocode}")
        header, first_row = validate_bgzf(summary_path)
        cases = int(source["cases"])
        controls = int(source["controls"])
        outcomes.append(
            {
                "phenocode": phenocode,
                "phenotype": manifest["phenotype"],
                "category": manifest["category"],
                "cases": cases,
                "controls": controls,
                "analyzed_sample_size": cases + controls,
                "effective_sample_size": 4 * cases * controls / (cases + controls),
                "summary_stats_url": source["url"],
                "local_path": summary_path.relative_to(root).as_posix(),
                "content_length_bytes": summary_path.stat().st_size,
                "md5_hex": summary_md5,
                "sha256": summary_sha256,
                "index_url": source["url"] + ".tbi",
                "index_local_path": index_path.relative_to(root).as_posix(),
                "index_content_length_bytes": index_path.stat().st_size,
                "index_md5_hex": index_md5,
                "index_sha256": index_sha256,
                "header_check": {
                    "decompression": "passed",
                    "header_column_count": len(header),
                    "first_row_column_count": len(first_row),
                    "first_data_row": "\t".join(first_row),
                },
                "verification_status": "verified_complete",
            }
        )

    payload = {
        "schema_version": "1.1.0",
        "release": {"name": "FinnGen R13", "genome_build": "GRCh38"},
        "effective_sample_size": {
            "formula": "4 * cases * controls / (cases + controls)",
            "purpose": "Binary-trait effective sample size for summary-statistic analyses",
        },
        "manifest": {
            "source_url": "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv",
            "local_path": manifest_path.relative_to(root).as_posix(),
            "content_length_bytes": manifest_path.stat().st_size,
            "md5_hex": manifest_md5,
            "sha256": manifest_sha256,
        },
        "file_format": {
            "compression": "BGZF",
            "delimiter": "tab",
            "header": EXPECTED_HEADER,
            "effect_allele": "alt",
        },
        "outcomes": outcomes,
    }
    output = root / "artifacts" / "mr" / "finngen_metadata.json"
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=True, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(json.dumps({"outcomes": len(outcomes), "status": "verified_complete"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
