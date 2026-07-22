#!/usr/bin/env python3
"""Extract and validate the Finan 2017 druggable gene universe."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import openpyxl


EXPECTED_COLUMNS = [
    "ensembl_gene_id",
    "druggability_tier",
    "hgnc_names",
    "chr_b37",
    "start_b37",
    "end_b37",
    "strand",
    "description",
    "no_of_gwas_regions",
    "small_mol_druggable",
    "bio_druggable",
    "adme_gene",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    source_path = repo_root / "data" / "raw" / "mr" / "finan_2017_table_s1.xlsx"
    output_path = repo_root / "data" / "processed" / "mr" / "druggable_genes.tsv"
    summary_path = repo_root / "logs" / "mr_druggable_gene_summary.json"

    workbook = openpyxl.load_workbook(source_path, read_only=True, data_only=True)
    if workbook.sheetnames != ["Data", "Column Descriptions", "Legend"]:
        raise ValueError(f"Unexpected workbook sheets: {workbook.sheetnames}")

    rows = workbook["Data"].iter_rows(values_only=True)
    header = list(next(rows))
    if header != EXPECTED_COLUMNS:
        raise ValueError(f"Unexpected Data header: {header}")

    records: list[list[object]] = []
    seen_gene_ids: set[str] = set()
    tier_counts: dict[str, int] = {}
    missing_symbol_count = 0
    non_empty_symbols: set[str] = set()
    flag_counts = {
        "small_mol_druggable_y": 0,
        "bio_druggable_y": 0,
        "adme_gene_y": 0,
    }
    for row_number, row in enumerate(rows, start=2):
        record = list(row)
        gene_id = str(record[0]).strip() if record[0] is not None else ""
        tier = str(record[1]).strip() if record[1] is not None else ""
        symbol = str(record[2]).strip() if record[2] is not None else ""
        if not gene_id.startswith("ENSG"):
            raise ValueError(f"Invalid Ensembl gene ID at row {row_number}: {gene_id!r}")
        if gene_id in seen_gene_ids:
            raise ValueError(f"Duplicate Ensembl gene ID at row {row_number}: {gene_id}")
        if tier not in {"Tier 1", "Tier 2", "Tier 3A", "Tier 3B"}:
            raise ValueError(f"Unexpected druggability tier at row {row_number}: {tier!r}")
        if not symbol:
            missing_symbol_count += 1
        else:
            non_empty_symbols.add(symbol)
        for output_key, column_index in (
            ("small_mol_druggable_y", 9),
            ("bio_druggable_y", 10),
            ("adme_gene_y", 11),
        ):
            if str(record[column_index]).strip() == "Y":
                flag_counts[output_key] += 1
        seen_gene_ids.add(gene_id)
        tier_counts[tier] = tier_counts.get(tier, 0) + 1
        records.append(record)
    workbook.close()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(EXPECTED_COLUMNS) + "\n")
        for record in records:
            handle.write("\t".join("" if value is None else str(value) for value in record) + "\n")

    summary = {
        "source_path": source_path.relative_to(repo_root).as_posix(),
        "source_bytes": source_path.stat().st_size,
        "source_sha256": sha256_file(source_path),
        "output_path": output_path.relative_to(repo_root).as_posix(),
        "output_bytes": output_path.stat().st_size,
        "output_sha256": sha256_file(output_path),
        "gene_count": len(records),
        "unique_ensembl_gene_count": len(seen_gene_ids),
        "unique_non_empty_hgnc_symbol_count": len(non_empty_symbols),
        "missing_hgnc_symbol_count": missing_symbol_count,
        "tier_counts": dict(sorted(tier_counts.items())),
        "flag_counts": flag_counts,
        "ensembl_version_note": (
            "The workbook column description says Ensembl version 75; the associated paper text "
            "mentions version 73. This discrepancy is retained rather than resolved by assumption."
        ),
        "gene_count_claim_scope": "Finan 2017 Table S1 rows after exact Ensembl ID deduplication",
    }
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
