"""Primary cis-eQTL MR screen across the Finan 2017 druggable genome."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd

from idd_genetics_pipeline import (
    _collect_column,
    _enrich_eqtl,
    _extract_af,
    _extract_eqtl,
    _extract_outcome,
    _load_finngen_manifest,
    _require_complete_gzip,
    best_by,
    load_rows,
    write_json,
    write_tsv,
)
from idd_genetics_utils import effective_binary_n, harmonize_effects, observed_r2_from_z, wald_ratio


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--finan", type=Path, default=ROOT / "data" / "raw" / "mr" / "finan_2017_table_s1.xlsx")
    parser.add_argument("--bulk-candidates", type=Path, default=ROOT / "artifacts" / "bulk" / "druggable_bulk_candidates.tsv")
    parser.add_argument("--raw-dir", type=Path, default=ROOT / "data" / "raw" / "idd" / "genetics")
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data" / "processed" / "idd" / "genetics")
    parser.add_argument("--artifacts-dir", type=Path, default=ROOT / "artifacts" / "mr")
    parser.add_argument("--ld-bim", type=Path, default=ROOT / "data" / "raw" / "idd" / "genetics" / "ld" / "1000G_EUR.bim")
    parser.add_argument("--instrument-p", type=float, default=5e-8)
    parser.add_argument("--f-min-exclusive", type=float, default=10.0)
    parser.add_argument("--mr-fdr", type=float, default=0.05)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--reuse-prepared-eqtl", action="store_true")
    parser.add_argument("--reuse-outcomes", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.artifacts_dir.mkdir(parents=True, exist_ok=True)
    universe = read_finan_universe(args.finan)
    gene_to_symbol = {row["ensembl_gene_id"]: row["gene_symbol"] for row in universe}
    eqtl_dir = args.raw_dir / "eqtlgen"
    significant_source = eqtl_dir / "eqtlgen_significant_cis.txt.gz"
    af_source = eqtl_dir / "eqtlgen_allele_frequency.txt.gz"
    _require_complete_gzip(significant_source)
    _require_complete_gzip(af_source)

    raw_extract = args.processed_dir / "druggable_significant_eqtl.raw.tsv.gz"
    enriched_extract = args.processed_dir / "druggable_significant_eqtl.tsv.gz"
    if not args.reuse_prepared_eqtl:
        _extract_eqtl(
            significant_source,
            raw_extract,
            gene_to_symbol,
            pvalue_max=args.instrument_p,
            f_min_exclusive=args.f_min_exclusive,
            window_bp=1_000_000,
        )
        snps = _collect_column(raw_extract, "snp")
        af = _extract_af(
            af_source,
            snps,
            args.processed_dir / "druggable_instrument_allele_frequency.tsv.gz",
        )
        _enrich_eqtl(raw_extract, enriched_extract, af)
    elif not enriched_extract.exists():
        raise FileNotFoundError(enriched_extract)
    snps = _collect_column(enriched_extract, "snp")
    enriched_rows = load_rows(enriched_extract)
    lead_candidates = choose_lead_instruments(enriched_rows)
    if not args.ld_bim.exists():
        raise FileNotFoundError(f"1000G EUR BIM is required for the primary screen: {args.ld_bim}")
    present_snps = read_reference_snps(args.ld_bim, {row["snp"] for row in enriched_rows})
    lead_by_gene = choose_lead_instruments(row for row in enriched_rows if row["snp"] in present_snps)
    lead_table = [
        {
            **row,
            "strongest_overall_snp": lead_candidates[gene]["snp"],
            "strongest_overall_present_in_1000g_eur": lead_candidates[gene]["snp"] in present_snps,
        }
        for gene, row in lead_by_gene.items()
    ]
    write_tsv(args.artifacts_dir / "druggable_lead_instruments.tsv", lead_table)
    if args.prepare_only:
        write_json(
            args.artifacts_dir / "druggable_eqtl_screen_summary.json",
            {
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "finan_genes_with_hgnc": len(universe),
                "genes_with_significant_cis_eqtl": len(lead_candidates),
                "genes_with_lead_in_1000g_eur": len(lead_by_gene),
                "instrument_unique_snps": len(snps),
                "status": "eqtl_prepared_outcomes_pending",
            },
        )
        return 0

    lead_snps = {row["snp"] for row in lead_by_gene.values()}
    outcomes = _load_finngen_manifest(args.raw_dir / "finngen" / "finngen_R13_manifest.tsv")
    outcome_extracts: dict[str, Path] = {}
    for outcome_id in outcomes:
        source = args.raw_dir / "finngen" / f"finngen_R13_{outcome_id}.gz"
        destination = args.processed_dir / f"{outcome_id}.druggable_leads.tsv.gz"
        if not args.reuse_outcomes:
            _require_complete_gzip(source)
            _extract_outcome(source, destination, outcome_id, lead_snps)
        elif not destination.exists():
            raise FileNotFoundError(destination)
        outcome_extracts[outcome_id] = destination

    screen_rows: list[dict[str, object]] = []
    for outcome_id, metadata in outcomes.items():
        outcome_by_snp = best_by(load_rows(outcome_extracts[outcome_id]), "snp", "pvalue")
        n_effective = effective_binary_n(float(metadata["cases"]), float(metadata["controls"]))
        outcome_rows = []
        for gene in universe:
            lead = lead_by_gene.get(gene["ensembl_gene_id"])
            base: dict[str, object] = {
                **gene,
                "outcome_id": outcome_id,
                "snp": lead["snp"] if lead else "",
                "beta": "",
                "se": "",
                "ci_lower": "",
                "ci_upper": "",
                "odds_ratio": "",
                "or_ci_lower": "",
                "or_ci_upper": "",
                "pvalue": "",
                "fdr": "",
                "exposure_r2": "",
                "outcome_r2": "",
                "steiger_direction_supported": "",
                "status": "not_estimable",
                "reason": "no_genome_wide_significant_cis_eqtl",
            }
            if lead is None:
                if gene["ensembl_gene_id"] in lead_candidates:
                    base["reason"] = "no_significant_cis_eqtl_present_in_1000g_eur"
                outcome_rows.append(base)
                continue
            outcome = outcome_by_snp.get(lead["snp"])
            if outcome is None:
                base["reason"] = "lead_snp_missing_from_outcome"
                outcome_rows.append(base)
                continue
            if not lead.get("beta") or not lead.get("se"):
                base["reason"] = f"eqtl_effect_unavailable:{lead.get('af_status', 'unknown')}"
                outcome_rows.append(base)
                continue
            aligned, reason = harmonize_effects(lead, outcome)
            if aligned is None:
                base["reason"] = f"harmonization:{reason}"
                outcome_rows.append(base)
                continue
            effect = wald_ratio(
                aligned.beta_exposure,
                aligned.se_exposure,
                aligned.beta_outcome,
                aligned.se_outcome,
            )
            exposure_z = aligned.beta_exposure / aligned.se_exposure
            outcome_z = aligned.beta_outcome / aligned.se_outcome
            exposure_r2 = observed_r2_from_z(exposure_z, float(lead["n"]))
            outcome_r2 = observed_r2_from_z(outcome_z, n_effective)
            base.update(effect)
            base.update(
                {
                    "exposure_r2": exposure_r2,
                    "outcome_r2": outcome_r2,
                    "steiger_direction_supported": exposure_r2 > outcome_r2,
                    "status": "estimable",
                    "reason": "",
                }
            )
            outcome_rows.append(base)
        apply_bh(outcome_rows, "pvalue", "fdr")
        screen_rows.extend(outcome_rows)

    write_tsv(args.artifacts_dir / "druggable_mr_screen.tsv", screen_rows)
    bulk = read_tsv_candidates(args.bulk_candidates)
    hits = {
        row["ensembl_gene_id"]
        for row in screen_rows
        if row["status"] == "estimable" and float(row["fdr"]) < args.mr_fdr
    }
    formal_by_gene = {row["ensembl_gene_id"]: {**row, "selection_reason": "mr_fdr"} for row in universe if row["ensembl_gene_id"] in hits}
    for row in bulk:
        current = formal_by_gene.get(row["ensembl_gene_id"])
        if current:
            current["selection_reason"] = "mr_fdr_and_bulk_candidate"
        else:
            formal_by_gene[row["ensembl_gene_id"]] = {**row, "selection_reason": "bulk_candidate"}
    formal = sorted(formal_by_gene.values(), key=lambda row: (row["gene_symbol"], row["ensembl_gene_id"]))
    write_tsv(args.artifacts_dir / "formal_candidate_genes.tsv", formal)
    write_json(
        args.artifacts_dir / "druggable_screen_summary.json",
        {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "finan_genes_with_hgnc": len(universe),
            "genes_with_significant_cis_eqtl": len(lead_candidates),
            "genes_with_lead_in_1000g_eur": len(lead_by_gene),
            "instrument_unique_snps": len(snps),
            "mr_fdr_threshold": args.mr_fdr,
            "mr_fdr_candidate_genes": len(hits),
            "bulk_candidate_genes": len(bulk),
            "formal_candidate_genes": len(formal),
            "primary_screen_method": "single_strongest_cis_eqtl_wald_ratio",
            "ld_note": "One lead variant per gene is LD-independent by construction; formal candidates undergo 1000G EUR multi-variant clumping.",
        },
    )
    return 0


def read_finan_universe(path: Path) -> list[dict[str, str]]:
    frame = pd.read_excel(path, sheet_name="Data", dtype=str)
    required = {"ensembl_gene_id", "hgnc_names", "druggability_tier"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Finan table lacks columns: {sorted(required - set(frame.columns))}")
    result = []
    seen = set()
    for row in frame.to_dict(orient="records"):
        gene = clean_string(row.get("ensembl_gene_id")).split(".")[0]
        symbol = clean_string(row.get("hgnc_names"))
        if not gene or not symbol or gene in seen:
            continue
        seen.add(gene)
        result.append(
            {
                "gene_symbol": symbol,
                "ensembl_gene_id": gene,
                "druggability_tier": clean_string(row.get("druggability_tier")),
            }
        )
    return result


def read_tsv_candidates(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = []
    seen = set()
    for row in rows:
        gene = row["ensembl_gene_id"].split(".")[0]
        if gene and gene not in seen:
            seen.add(gene)
            result.append(
                {
                    "gene_symbol": row["gene_symbol"],
                    "ensembl_gene_id": gene,
                    "druggability_tier": row.get("druggability_tier", ""),
                }
            )
    return result


def choose_lead_instruments(rows: Iterable[Mapping[str, str]]) -> dict[str, Mapping[str, str]]:
    result: dict[str, Mapping[str, str]] = {}
    for row in rows:
        if not row.get("beta") or not row.get("se"):
            continue
        gene = row["ensembl_gene_id"]
        if gene not in result or float(row["pvalue"]) < float(result[gene]["pvalue"]):
            result[gene] = row
    return result


def read_reference_snps(path: Path, requested: set[str]) -> set[str]:
    present = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2 and fields[1] in requested:
                present.add(fields[1])
                if len(present) == len(requested):
                    break
    return present


def apply_bh(rows: list[dict[str, object]], p_column: str, output_column: str) -> None:
    estimable = [(index, float(row[p_column])) for index, row in enumerate(rows) if row.get(p_column) not in (None, "")]
    estimable.sort(key=lambda item: item[1])
    m = len(estimable)
    adjusted = [1.0] * m
    running = 1.0
    for rank_index in range(m - 1, -1, -1):
        _, pvalue = estimable[rank_index]
        rank = rank_index + 1
        running = min(running, pvalue * m / rank)
        adjusted[rank_index] = min(1.0, running)
    for (row_index, _), value in zip(estimable, adjusted):
        rows[row_index][output_column] = value


def clean_string(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


if __name__ == "__main__":
    raise SystemExit(main())
