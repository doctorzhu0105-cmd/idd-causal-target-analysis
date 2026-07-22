"""Run the eQTL-based causal genetics workflow for IDD/LDH candidates.

The pipeline intentionally streams the multi-gigabyte source files. It writes
candidate-only extracts before statistical analysis, so an interrupted later
stage never requires another full-source scan.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Mapping

from idd_genetics_utils import (
    benjamini_hochberg,
    coloc_posteriors,
    effective_binary_n,
    harmonize_effects,
    ivw_fixed,
    mr_egger,
    observed_r2_from_z,
    standardized_effect_from_z,
    wald_ratio,
    wakefield_log_abf,
    weighted_median,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RAW = ROOT / "data" / "raw" / "idd" / "genetics"
DEFAULT_PROCESSED = ROOT / "data" / "processed" / "idd" / "genetics"
DEFAULT_ARTIFACTS = ROOT / "artifacts" / "mr"
EQTL_COLUMNS = [
    "gene_symbol",
    "ensembl_gene_id",
    "snp",
    "chr_b37",
    "pos_b37",
    "gene_pos_b37",
    "effect_allele",
    "other_allele",
    "zscore",
    "n",
    "pvalue",
    "fdr",
]
ENRICHED_EQTL_COLUMNS = EQTL_COLUMNS + ["eaf", "beta", "se", "f_statistic", "af_status"]
OUTCOME_COLUMNS = [
    "outcome_id",
    "snp",
    "chr_b38",
    "pos_b38",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "eaf",
    "pvalue",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("extract", "analyze", "validate", "run"))
    parser.add_argument("--candidates", type=Path, default=ROOT / "artifacts" / "mr" / "formal_candidate_genes.tsv")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--processed-dir", type=Path, default=DEFAULT_PROCESSED)
    parser.add_argument("--artifacts-dir", type=Path, default=DEFAULT_ARTIFACTS)
    parser.add_argument("--plink", type=Path, default=ROOT / "tools" / "idd" / "plink" / "plink.exe")
    parser.add_argument("--ld-prefix", type=Path, default=DEFAULT_RAW / "ld" / "1000G_EUR")
    parser.add_argument("--instrument-p", type=float, default=5e-8)
    parser.add_argument("--f-min-exclusive", type=float, default=10.0)
    parser.add_argument("--clump-r2", type=float, default=0.001)
    parser.add_argument("--clump-kb", type=int, default=10000)
    parser.add_argument("--coloc-window-bp", type=int, default=1_000_000)
    parser.add_argument("--coloc-min-shared", type=int, default=50)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.artifacts_dir.mkdir(parents=True, exist_ok=True)
    if args.command in ("extract", "run"):
        extract(args)
    if args.command in ("analyze", "run"):
        analyze(args)
    if args.command in ("validate", "run"):
        validate_sources(args)
    return 0


def read_candidates(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"gene_symbol", "ensembl_gene_id"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"candidate file lacks required columns: {path}")
    seen: set[str] = set()
    result: list[dict[str, str]] = []
    for row in rows:
        gene = row["ensembl_gene_id"].split(".")[0]
        if gene and gene not in seen:
            seen.add(gene)
            result.append({"ensembl_gene_id": gene, "gene_symbol": row["gene_symbol"]})
    return result


def extract(args: argparse.Namespace) -> None:
    candidates = read_candidates(args.candidates)
    gene_to_symbol = {row["ensembl_gene_id"]: row["gene_symbol"] for row in candidates}
    eqtl_dir = args.raw_dir / "eqtlgen"
    significant_source = eqtl_dir / "eqtlgen_significant_cis.txt.gz"
    full_source = eqtl_dir / "eqtlgen_full_cis.txt.gz"
    af_source = eqtl_dir / "eqtlgen_allele_frequency.txt.gz"
    significant_raw = args.processed_dir / "candidate_significant_eqtl.raw.tsv.gz"
    full_raw = args.processed_dir / "candidate_full_eqtl.raw.tsv.gz"
    _require_complete_gzip(significant_source)
    _require_complete_gzip(full_source)
    _require_complete_gzip(af_source)

    _extract_eqtl(
        significant_source,
        significant_raw,
        gene_to_symbol,
        pvalue_max=args.instrument_p,
        f_min_exclusive=args.f_min_exclusive,
        window_bp=args.coloc_window_bp,
    )
    _extract_eqtl(
        full_source,
        full_raw,
        gene_to_symbol,
        pvalue_max=None,
        f_min_exclusive=None,
        window_bp=args.coloc_window_bp,
    )

    requested_snps = _collect_column(full_raw, "snp") | _collect_column(significant_raw, "snp")
    af = _extract_af(af_source, requested_snps, args.processed_dir / "candidate_allele_frequency.tsv.gz")
    _enrich_eqtl(significant_raw, args.processed_dir / "candidate_significant_eqtl.tsv.gz", af)
    _enrich_eqtl(full_raw, args.processed_dir / "candidate_full_eqtl.tsv.gz", af)

    outcome_meta = _load_finngen_manifest(args.raw_dir / "finngen" / "finngen_R13_manifest.tsv")
    for outcome_id in outcome_meta:
        source = args.raw_dir / "finngen" / f"finngen_R13_{outcome_id}.gz"
        _require_complete_gzip(source)
        destination = args.processed_dir / f"{outcome_id}.candidate_variants.tsv.gz"
        _extract_outcome(source, destination, outcome_id, requested_snps)

    summary_rows = []
    sig_counts = _count_by(args.processed_dir / "candidate_significant_eqtl.tsv.gz", "ensembl_gene_id")
    full_counts = _count_by(args.processed_dir / "candidate_full_eqtl.tsv.gz", "ensembl_gene_id")
    for candidate in candidates:
        gene = candidate["ensembl_gene_id"]
        summary_rows.append(
            {
                **candidate,
                "significant_instrument_candidates": sig_counts.get(gene, 0),
                "full_cis_variants": full_counts.get(gene, 0),
            }
        )
    write_tsv(args.artifacts_dir / "candidate_extraction_summary.tsv", summary_rows)


def _extract_eqtl(
    source: Path,
    destination: Path,
    gene_to_symbol: Mapping[str, str],
    *,
    pvalue_max: float | None,
    f_min_exclusive: float | None,
    window_bp: int,
) -> None:
    genes = set(gene_to_symbol)
    kept = 0
    scanned = 0
    with gzip.open(source, "rt", encoding="utf-8", newline="") as input_handle, atomic_gzip_writer(destination) as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=EQTL_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        header = input_handle.readline().rstrip("\r\n").split("\t")
        index = {name: position for position, name in enumerate(header)}
        required = {
            "Gene", "Pvalue", "Zscore", "SNPPos", "GenePos", "SNP", "SNPChr",
            "AssessedAllele", "OtherAllele", "NrSamples", "FDR",
        }
        if not required.issubset(index):
            raise ValueError(f"unexpected eQTLGen header in {source}")
        for line in input_handle:
            scanned += 1
            fields = line.rstrip("\r\n").split("\t")
            gene = fields[index["Gene"]].split(".")[0]
            if gene not in genes:
                continue
            pvalue = float(fields[index["Pvalue"]])
            zscore = float(fields[index["Zscore"]])
            if pvalue_max is not None and pvalue > pvalue_max:
                continue
            if f_min_exclusive is not None and zscore * zscore <= f_min_exclusive:
                continue
            snp_pos = int(fields[index["SNPPos"]])
            gene_pos = int(fields[index["GenePos"]])
            if abs(snp_pos - gene_pos) > window_bp:
                continue
            writer.writerow(
                {
                    "gene_symbol": gene_to_symbol[gene],
                    "ensembl_gene_id": gene,
                    "snp": fields[index["SNP"]],
                    "chr_b37": fields[index["SNPChr"]],
                    "pos_b37": snp_pos,
                    "gene_pos_b37": gene_pos,
                    "effect_allele": fields[index["AssessedAllele"]].upper(),
                    "other_allele": fields[index["OtherAllele"]].upper(),
                    "zscore": fields[index["Zscore"]],
                    "n": fields[index["NrSamples"]],
                    "pvalue": fields[index["Pvalue"]],
                    "fdr": fields[index["FDR"]],
                }
            )
            kept += 1
            if scanned % 10_000_000 == 0:
                print(f"{source.name}: scanned={scanned:,}, kept={kept:,}", file=sys.stderr)
    print(f"{source.name}: scanned={scanned:,}, kept={kept:,}", file=sys.stderr)


def _extract_af(source: Path, requested_snps: set[str], destination: Path) -> dict[str, tuple[str, str, float]]:
    result: dict[str, tuple[str, str, float]] = {}
    with gzip.open(source, "rt", encoding="utf-8", newline="") as input_handle, atomic_gzip_writer(destination) as output_handle:
        fields = ["snp", "allele_a", "allele_b", "allele_b_frequency"]
        writer = csv.DictWriter(output_handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        header = input_handle.readline().rstrip("\r\n").split("\t")
        index = {name: position for position, name in enumerate(header)}
        for line in input_handle:
            values = line.rstrip("\r\n").split("\t")
            snp = values[index["SNP"]]
            if snp not in requested_snps:
                continue
            value = (
                values[index["AlleleA"]].upper(),
                values[index["AlleleB"]].upper(),
                float(values[index["AlleleB_all"]]),
            )
            result[snp] = value
            writer.writerow(dict(zip(fields, (snp, *value))))
    return result


def _enrich_eqtl(source: Path, destination: Path, af: Mapping[str, tuple[str, str, float]]) -> None:
    with gzip.open(source, "rt", encoding="utf-8", newline="") as input_handle, atomic_gzip_writer(destination) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        writer = csv.DictWriter(output_handle, fieldnames=ENRICHED_EQTL_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            af_record = af.get(row["snp"])
            eaf: float | None = None
            af_status = "missing"
            if af_record:
                allele_a, allele_b, allele_b_frequency = af_record
                if row["effect_allele"] == allele_b and row["other_allele"] == allele_a:
                    eaf = allele_b_frequency
                    af_status = "matched_allele_b"
                elif row["effect_allele"] == allele_a and row["other_allele"] == allele_b:
                    eaf = 1.0 - allele_b_frequency
                    af_status = "matched_allele_a"
                else:
                    af_status = "allele_mismatch"
            beta = se = None
            if eaf is not None and 0 < eaf < 1:
                beta, se = standardized_effect_from_z(float(row["zscore"]), eaf, float(row["n"]))
            writer.writerow(
                {
                    **row,
                    "eaf": format_value(eaf),
                    "beta": format_value(beta),
                    "se": format_value(se),
                    "f_statistic": format_value(float(row["zscore"]) ** 2),
                    "af_status": af_status,
                }
            )


def _extract_outcome(source: Path, destination: Path, outcome_id: str, requested_snps: set[str]) -> None:
    with gzip.open(source, "rt", encoding="utf-8", newline="") as input_handle, atomic_gzip_writer(destination) as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=OUTCOME_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        header = input_handle.readline().rstrip("\r\n").split("\t")
        index = {name: position for position, name in enumerate(header)}
        for line in input_handle:
            values = line.rstrip("\r\n").split("\t")
            rsid_value = values[index["rsids"]]
            if "rs" not in rsid_value:
                continue
            matched = split_rsids(rsid_value) & requested_snps
            for snp in sorted(matched):
                writer.writerow(
                    {
                        "outcome_id": outcome_id,
                        "snp": snp,
                        "chr_b38": values[index["#chrom"]],
                        "pos_b38": values[index["pos"]],
                        "effect_allele": values[index["alt"]].upper(),
                        "other_allele": values[index["ref"]].upper(),
                        "beta": values[index["beta"]],
                        "se": values[index["sebeta"]],
                        "eaf": values[index["af_alt"]],
                        "pvalue": values[index["pval"]],
                    }
                )


def analyze(args: argparse.Namespace) -> None:
    candidates = read_candidates(args.candidates)
    outcome_meta = _load_finngen_manifest(args.raw_dir / "finngen" / "finngen_R13_manifest.tsv")
    sig_rows = load_rows(args.processed_dir / "candidate_significant_eqtl.tsv.gz")
    full_rows = load_rows(args.processed_dir / "candidate_full_eqtl.tsv.gz")
    sig_by_gene = group_best(sig_rows, "ensembl_gene_id", "snp", "pvalue")
    full_by_gene = group_best(full_rows, "ensembl_gene_id", "snp", "pvalue")

    clump_dir = args.processed_dir / "clump"
    clump_dir.mkdir(parents=True, exist_ok=True)
    selected_forward: dict[str, set[str]] = {}
    forward_clump_status: dict[str, str] = {}
    instrument_summary: list[dict[str, object]] = []
    instrument_rows: list[dict[str, object]] = []
    for candidate in candidates:
        gene = candidate["ensembl_gene_id"]
        rows = list(sig_by_gene.get(gene, {}).values())
        selected, clump_status = run_clump(
            rows,
            p_column="pvalue",
            label=f"forward_{gene}",
            args=args,
            clump_dir=clump_dir,
        )
        selected_forward[gene] = selected
        forward_clump_status[gene] = clump_status
        instrument_summary.append(
            {
                **candidate,
                "screened_variants": len(rows),
                "clumped_instruments": len(selected),
                "status": "estimable" if selected else "not_estimable",
                "reason": "" if selected else clump_status,
            }
        )
        for row in rows:
            instrument_rows.append(
                {
                    **candidate,
                    "snp": row["snp"],
                    "pvalue": row["pvalue"],
                    "f_statistic": row["f_statistic"],
                    "selected_after_clump": row["snp"] in selected,
                    "af_status": row["af_status"],
                }
            )

    write_tsv(args.artifacts_dir / "instrument_screen_summary.tsv", instrument_summary)
    write_tsv(args.artifacts_dir / "instrument_candidates.tsv", instrument_rows)

    mr_results: list[dict[str, object]] = []
    steiger_results: list[dict[str, object]] = []
    reverse_results: list[dict[str, object]] = []
    coloc_results: list[dict[str, object]] = []
    harmonized_details: list[dict[str, object]] = []
    sensitivity_results: list[dict[str, object]] = []
    leave_one_out_results: list[dict[str, object]] = []

    for outcome_id, metadata in outcome_meta.items():
        outcome_rows = load_rows(args.processed_dir / f"{outcome_id}.candidate_variants.tsv.gz")
        outcome_by_snp = best_by(outcome_rows, "snp", "pvalue")
        n_effective = effective_binary_n(metadata["cases"], metadata["controls"])
        for candidate in candidates:
            gene = candidate["ensembl_gene_id"]
            gene_symbol = candidate["gene_symbol"]
            forward = []
            exclusion_counts: dict[str, int] = defaultdict(int)
            for snp in sorted(selected_forward.get(gene, set())):
                exposure = sig_by_gene.get(gene, {}).get(snp)
                outcome = outcome_by_snp.get(snp)
                if exposure is None or outcome is None:
                    exclusion_counts["missing_outcome"] += 1
                    continue
                if not exposure.get("beta") or not exposure.get("se"):
                    exclusion_counts["missing_exposure_effect"] += 1
                    continue
                aligned, reason = harmonize_effects(exposure, outcome)
                if aligned is None:
                    exclusion_counts[reason] += 1
                    continue
                detail = {
                    "gene_symbol": gene_symbol,
                    "ensembl_gene_id": gene,
                    "outcome_id": outcome_id,
                    "direction": "expression_to_disease",
                    "snp": snp,
                    "beta_exposure": aligned.beta_exposure,
                    "se_exposure": aligned.se_exposure,
                    "exposure_n": float(exposure["n"]),
                    "beta_outcome": aligned.beta_outcome,
                    "se_outcome": aligned.se_outcome,
                    "effect_allele": aligned.effect_allele,
                    "other_allele": aligned.other_allele,
                    "harmonization_action": aligned.action,
                }
                forward.append(detail)
                harmonized_details.append(detail)
            if not forward and not exclusion_counts:
                exclusion_counts[forward_clump_status.get(gene, "no_harmonized_instruments")] += 1
            mr_results.append(_mr_result(candidate, outcome_id, "expression_to_disease", forward, exclusion_counts))
            sensitivity_results.extend(
                _sensitivity_results(candidate, outcome_id, "expression_to_disease", forward)
            )
            leave_one_out_results.extend(
                _leave_one_out_results(candidate, outcome_id, "expression_to_disease", forward)
            )
            steiger_results.append(_steiger_result(candidate, outcome_id, forward, n_effective))

            full_gene = full_by_gene.get(gene, {})
            reverse_candidates = [
                outcome_by_snp[snp]
                for snp in full_gene.keys() & outcome_by_snp.keys()
                if float(outcome_by_snp[snp]["pvalue"]) <= args.instrument_p
            ]
            reverse_selected, reverse_clump_status = run_clump(
                reverse_candidates,
                p_column="pvalue",
                label=f"reverse_{outcome_id}_{gene}",
                args=args,
                clump_dir=clump_dir,
            )
            reverse = []
            reverse_exclusions: dict[str, int] = defaultdict(int)
            for snp in sorted(reverse_selected):
                disease = outcome_by_snp[snp]
                expression = full_gene[snp]
                if not expression.get("beta") or not expression.get("se"):
                    reverse_exclusions["missing_expression_effect"] += 1
                    continue
                aligned, reason = harmonize_effects(disease, expression)
                if aligned is None:
                    reverse_exclusions[reason] += 1
                    continue
                detail = {
                    "gene_symbol": gene_symbol,
                    "ensembl_gene_id": gene,
                    "outcome_id": outcome_id,
                    "direction": "disease_to_expression",
                    "snp": snp,
                    "beta_exposure": aligned.beta_exposure,
                    "se_exposure": aligned.se_exposure,
                    "beta_outcome": aligned.beta_outcome,
                    "se_outcome": aligned.se_outcome,
                    "effect_allele": aligned.effect_allele,
                    "other_allele": aligned.other_allele,
                    "harmonization_action": aligned.action,
                }
                reverse.append(detail)
                harmonized_details.append(detail)
            if not reverse and not reverse_exclusions:
                reverse_exclusions[reverse_clump_status] += 1
            reverse_results.append(_mr_result(candidate, outcome_id, "disease_to_expression", reverse, reverse_exclusions))
            sensitivity_results.extend(
                _sensitivity_results(candidate, outcome_id, "disease_to_expression", reverse)
            )
            leave_one_out_results.extend(
                _leave_one_out_results(candidate, outcome_id, "disease_to_expression", reverse)
            )
            coloc_results.append(
                _coloc_result(
                    candidate,
                    outcome_id,
                    full_gene,
                    outcome_by_snp,
                    args.coloc_min_shared,
                )
            )

    apply_fdr(mr_results, "pvalue", "fdr", ("outcome_id", "direction"))
    apply_fdr(mr_results, "q_pvalue", "q_fdr", ("outcome_id", "direction"))
    apply_fdr(reverse_results, "pvalue", "fdr", ("outcome_id", "direction"))
    apply_fdr(reverse_results, "q_pvalue", "q_fdr", ("outcome_id", "direction"))
    apply_fdr(sensitivity_results, "pvalue", "fdr", ("outcome_id", "direction", "method"))
    apply_fdr(
        sensitivity_results,
        "egger_intercept_pvalue",
        "egger_intercept_fdr",
        ("outcome_id", "direction"),
    )
    apply_fdr(sensitivity_results, "q_pvalue", "q_fdr", ("outcome_id", "direction", "method"))
    apply_fdr(leave_one_out_results, "pvalue", "fdr", ("outcome_id", "direction"))
    for row in [*mr_results, *reverse_results]:
        row["fdr_scope"] = "formal_candidate_set_by_outcome_direction"
        row["q_fdr_scope"] = "formal_candidate_set_by_outcome_direction"
    for row in coloc_results:
        row["fdr"] = ""
        row["fdr_scope"] = "not_applicable"
        row["fdr_reason"] = "coloc_pp_are_posterior_probabilities_not_pvalues"

    write_tsv(args.artifacts_dir / "harmonized_instruments.tsv", harmonized_details)
    write_tsv(args.artifacts_dir / "mr_results.tsv", mr_results)
    write_tsv(args.artifacts_dir / "mr_sensitivity_results.tsv", sensitivity_results)
    write_tsv(args.artifacts_dir / "leave_one_out_results.tsv", leave_one_out_results)
    write_tsv(args.artifacts_dir / "steiger_results.tsv", steiger_results)
    write_tsv(args.artifacts_dir / "reverse_mr_results.tsv", reverse_results)
    write_tsv(args.artifacts_dir / "coloc_results.tsv", coloc_results)
    analysis_config = json.loads((ROOT / "config" / "idd_analysis.json").read_text(encoding="utf-8"))
    run_summary = {
        "analysis_registered_at": analysis_config["registered_at"],
        "candidate_count": len(candidates),
        "outcomes": outcome_meta,
        "instrument_pvalue_max": args.instrument_p,
        "instrument_f_min_exclusive": args.f_min_exclusive,
        "clump_r2_max": args.clump_r2,
        "clump_window_kb": args.clump_kb,
        "coloc_priors": {"p1": 1e-4, "p2": 1e-4, "p12": 1e-5},
        "steiger_scale": "observed_summary_r2_using_effective_binary_sample_size",
        "genome_build_alignment": "rsid_and_alleles; no coordinate inference between GRCh37 and GRCh38",
    }
    write_json(args.artifacts_dir / "genetics_run_summary.json", run_summary)


def _mr_result(
    candidate: Mapping[str, str],
    outcome_id: str,
    direction: str,
    rows: list[dict[str, object]],
    exclusions: Mapping[str, int],
) -> dict[str, object]:
    base: dict[str, object] = {
        **candidate,
        "outcome_id": outcome_id,
        "direction": direction,
        "n_instruments": len(rows),
        "method": "not_estimable",
        "beta": "",
        "se": "",
        "ci_lower": "",
        "ci_upper": "",
        "odds_ratio": "",
        "or_ci_lower": "",
        "or_ci_upper": "",
        "pvalue": "",
        "q": "",
        "q_df": "",
        "q_pvalue": "",
        "status": "not_estimable",
        "reason": json.dumps(dict(exclusions), sort_keys=True),
    }
    if not rows:
        if not exclusions:
            base["reason"] = "no_harmonized_instruments"
        return base
    if len(rows) == 1:
        result = wald_ratio(
            float(rows[0]["beta_exposure"]),
            float(rows[0]["se_exposure"]),
            float(rows[0]["beta_outcome"]),
            float(rows[0]["se_outcome"]),
        )
        method = "wald_ratio"
    else:
        result = ivw_fixed(rows)
        method = "ivw_fixed"
    base.update(result)
    if direction == "disease_to_expression":
        base.update({"odds_ratio": "", "or_ci_lower": "", "or_ci_upper": ""})
    base.update({"method": method, "status": "estimable", "reason": ""})
    return base


def _sensitivity_results(
    candidate: Mapping[str, str],
    outcome_id: str,
    direction: str,
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    common = {
        **candidate,
        "outcome_id": outcome_id,
        "direction": direction,
        "n_instruments": len(rows),
        "beta": "",
        "se": "",
        "ci_lower": "",
        "ci_upper": "",
        "odds_ratio": "",
        "or_ci_lower": "",
        "or_ci_upper": "",
        "pvalue": "",
        "fdr": "",
        "q": "",
        "q_df": "",
        "q_pvalue": "",
        "q_fdr": "",
        "egger_intercept": "",
        "egger_intercept_se": "",
        "egger_intercept_pvalue": "",
        "egger_intercept_fdr": "",
        "status": "not_estimable",
        "reason": "fewer_than_3_instruments",
        "fdr_scope": "formal_candidate_set_by_outcome_direction_method",
    }
    if len(rows) < 3:
        results.append({**common, "method": "weighted_median"})
        results.append({**common, "method": "mr_egger"})
        return results

    seed_material = f"{candidate['ensembl_gene_id']}|{outcome_id}|{direction}|20260716"
    seed = int(hashlib.sha256(seed_material.encode("ascii")).hexdigest()[:8], 16)
    try:
        median = weighted_median(rows, seed=seed)
        median_row = {**common, "method": "weighted_median", "status": "estimable", "reason": ""}
        median_row.update(median)
        if direction == "disease_to_expression":
            median_row.update({"odds_ratio": "", "or_ci_lower": "", "or_ci_upper": ""})
        results.append(median_row)
    except ValueError as error:
        results.append({**common, "method": "weighted_median", "reason": str(error)})

    try:
        egger = mr_egger(rows)
        egger_row = {
            **common,
            "method": "mr_egger",
            "beta": egger["slope_beta"],
            "se": egger["slope_se"],
            "ci_lower": egger["slope_ci_lower"],
            "ci_upper": egger["slope_ci_upper"],
            "odds_ratio": egger["slope_odds_ratio"],
            "or_ci_lower": egger["slope_or_ci_lower"],
            "or_ci_upper": egger["slope_or_ci_upper"],
            "pvalue": egger["slope_pvalue"],
            "q": egger["q"],
            "q_df": egger["q_df"],
            "q_pvalue": egger["q_pvalue"],
            "egger_intercept": egger["intercept"],
            "egger_intercept_se": egger["intercept_se"],
            "egger_intercept_pvalue": egger["intercept_pvalue"],
            "status": "estimable",
            "reason": "",
        }
        if direction == "disease_to_expression":
            egger_row.update({"odds_ratio": "", "or_ci_lower": "", "or_ci_upper": ""})
        results.append(egger_row)
    except ValueError as error:
        results.append({**common, "method": "mr_egger", "reason": str(error)})
    return results


def _leave_one_out_results(
    candidate: Mapping[str, str],
    outcome_id: str,
    direction: str,
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    common = {
        **candidate,
        "outcome_id": outcome_id,
        "direction": direction,
        "n_instruments_full": len(rows),
        "excluded_snp": "",
        "remaining_instruments": max(0, len(rows) - 1),
        "method": "not_estimable",
        "beta": "",
        "se": "",
        "ci_lower": "",
        "ci_upper": "",
        "odds_ratio": "",
        "or_ci_lower": "",
        "or_ci_upper": "",
        "pvalue": "",
        "fdr": "",
        "status": "not_estimable",
        "reason": "fewer_than_2_instruments",
        "fdr_scope": "formal_candidate_leave_one_out_by_outcome_direction",
    }
    if len(rows) < 2:
        return [common]
    results = []
    for excluded_index, excluded in enumerate(rows):
        remaining = rows[:excluded_index] + rows[excluded_index + 1 :]
        if len(remaining) == 1:
            item = remaining[0]
            effect = wald_ratio(
                float(item["beta_exposure"]),
                float(item["se_exposure"]),
                float(item["beta_outcome"]),
                float(item["se_outcome"]),
            )
            method = "wald_ratio"
        else:
            effect = ivw_fixed(remaining)
            method = "ivw_fixed"
        row = {
            **common,
            "excluded_snp": excluded["snp"],
            "method": method,
            "status": "estimable",
            "reason": "",
        }
        row.update(effect)
        if direction == "disease_to_expression":
            row.update({"odds_ratio": "", "or_ci_lower": "", "or_ci_upper": ""})
        results.append(row)
    return results


def _steiger_result(
    candidate: Mapping[str, str],
    outcome_id: str,
    rows: list[dict[str, object]],
    outcome_effective_n: float,
) -> dict[str, object]:
    base: dict[str, object] = {
        **candidate,
        "outcome_id": outcome_id,
        "n_instruments": len(rows),
        "exposure_r2": "",
        "outcome_r2": "",
        "direction_supported": "",
        "status": "not_estimable",
        "reason": "no_harmonized_instruments",
    }
    if not rows:
        return base
    exposure_r2 = 0.0
    outcome_r2 = 0.0
    for row in rows:
        exposure_z = float(row["beta_exposure"]) / float(row["se_exposure"])
        outcome_z = float(row["beta_outcome"]) / float(row["se_outcome"])
        exposure_source_n = float(row["exposure_n"])
        exposure_r2 += observed_r2_from_z(exposure_z, exposure_source_n)
        outcome_r2 += observed_r2_from_z(outcome_z, outcome_effective_n)
    base.update(
        {
            "exposure_r2": exposure_r2,
            "outcome_r2": outcome_r2,
            "direction_supported": exposure_r2 > outcome_r2,
            "status": "estimable",
            "reason": "",
        }
    )
    return base


def _coloc_result(
    candidate: Mapping[str, str],
    outcome_id: str,
    eqtl_by_snp: Mapping[str, dict[str, str]],
    outcome_by_snp: Mapping[str, dict[str, str]],
    min_shared: int,
) -> dict[str, object]:
    shared = sorted(eqtl_by_snp.keys() & outcome_by_snp.keys())
    abf_eqtl: list[float] = []
    abf_outcome: list[float] = []
    for snp in shared:
        eqtl = eqtl_by_snp[snp]
        outcome = outcome_by_snp[snp]
        if not eqtl.get("beta") or not eqtl.get("se"):
            continue
        try:
            abf_eqtl.append(wakefield_log_abf(float(eqtl["beta"]), float(eqtl["se"]) ** 2, 0.15**2))
            abf_outcome.append(wakefield_log_abf(float(outcome["beta"]), float(outcome["se"]) ** 2, 0.20**2))
        except (ValueError, OverflowError):
            continue
    base: dict[str, object] = {
        **candidate,
        "outcome_id": outcome_id,
        "shared_variants": len(abf_eqtl),
        "pp0": "",
        "pp1": "",
        "pp2": "",
        "pp3": "",
        "pp4": "",
        "status": "not_estimable",
        "reason": "insufficient_shared_variants",
    }
    if len(abf_eqtl) < min_shared:
        return base
    posterior = coloc_posteriors(abf_eqtl, abf_outcome)
    base.update(posterior)
    base.update({"status": "estimable", "reason": ""})
    return base


def run_clump(
    rows: Iterable[Mapping[str, str]],
    *,
    p_column: str,
    label: str,
    args: argparse.Namespace,
    clump_dir: Path,
) -> tuple[set[str], str]:
    best = best_by(rows, "snp", p_column)
    if not best:
        return set(), "no_genome_wide_significant_variants"
    if not args.plink.exists():
        return set(), "plink_missing"
    if not all(args.ld_prefix.with_suffix(ext).exists() for ext in (".bed", ".bim", ".fam")):
        return set(), "ld_reference_missing"
    input_path = clump_dir / f"{label}.tsv"
    write_tsv(input_path, [{"SNP": snp, "P": row[p_column]} for snp, row in best.items()])
    prefix = clump_dir / label
    command = [
        str(args.plink),
        "--bfile",
        str(args.ld_prefix),
        "--clump",
        str(input_path),
        "--clump-snp-field",
        "SNP",
        "--clump-field",
        "P",
        "--clump-p1",
        str(args.instrument_p),
        "--clump-r2",
        str(args.clump_r2),
        "--clump-kb",
        str(args.clump_kb),
        "--allow-no-sex",
        "--out",
        str(prefix),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if completed.returncode != 0:
        return set(), f"plink_failed_{completed.returncode}"
    clumped = prefix.with_suffix(".clumped")
    if not clumped.exists() or not clumped.read_text(encoding="utf-8", errors="replace").strip():
        return set(), "no_variants_present_in_ld_reference"
    selected: set[str] = set()
    lines = [line.split() for line in clumped.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip()]
    header = lines[0]
    snp_index = header.index("SNP")
    for values in lines[1:]:
        if len(values) > snp_index:
            selected.add(values[snp_index])
    return selected, "ok" if selected else "no_variants_present_in_ld_reference"


def apply_fdr(
    rows: list[dict[str, object]],
    p_column: str,
    output_column: str,
    group_columns: tuple[str, ...],
) -> None:
    groups: dict[tuple[object, ...], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        groups[tuple(row.get(column, "") for column in group_columns)].append(index)
    for indices in groups.values():
        pvalues = []
        for index in indices:
            value = rows[index].get(p_column)
            try:
                parsed = float(value) if value not in (None, "") else None
            except (TypeError, ValueError):
                parsed = None
            pvalues.append(parsed if parsed is not None and math.isfinite(parsed) else None)
        adjusted = benjamini_hochberg(pvalues)
        for index, value in zip(indices, adjusted):
            rows[index][output_column] = "" if value is None else value


def validate_sources(args: argparse.Namespace) -> None:
    expected = {
        args.raw_dir / "eqtlgen" / "eqtlgen_significant_cis.txt.gz": (322_775_879, None),
        args.raw_dir / "eqtlgen" / "eqtlgen_full_cis.txt.gz": (3_880_187_026, None),
        args.raw_dir / "eqtlgen" / "eqtlgen_smr_besd.tar.gz": (1_372_512_049, None),
        args.raw_dir / "eqtlgen" / "eqtlgen_allele_frequency.txt.gz": (240_045_342, None),
        args.raw_dir / "finngen" / "finngen_R13_M13_INTERVERTEB.gz": (810_165_287, None),
        args.raw_dir / "finngen" / "finngen_R13_M13_LOWBACKPAIN.gz": (809_321_364, None),
        args.raw_dir / "finngen" / "finngen_R13_M13_SCIATICA.gz": (808_334_919, None),
        args.ld_prefix.with_suffix(".bed"): (231_387_159, "a163c74ea92676c9ada76fe5eb3b2ec9"),
        args.ld_prefix.with_suffix(".bim"): (51_537_337, "81b1ee40a05ef5be98ad3e496c71d708"),
        args.ld_prefix.with_suffix(".fam"): (9_557, "669a4260fda7a9e1dd7df374aff294ea"),
        args.ld_prefix.with_suffix(".afreq"): (56_496_658, "7746ab38e42c6abd5980c8ce57172327"),
    }
    rows = []
    for path, (expected_size, expected_md5) in expected.items():
        exists = path.exists()
        size = path.stat().st_size if exists else None
        size_matches = size == expected_size
        sha256 = sha256_file(path) if size_matches else ""
        observed_md5 = hash_file(path, "md5") if size_matches and expected_md5 else ""
        publisher_checksum_matches = observed_md5 == expected_md5 if expected_md5 else ""
        rows.append(
            {
                "path": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                "expected_size": expected_size,
                "observed_size": format_value(size),
                "size_matches": size_matches,
                "sha256": sha256,
                "publisher_md5": expected_md5 or "",
                "observed_md5": observed_md5,
                "publisher_checksum_matches": publisher_checksum_matches,
                "status": (
                    "verified"
                    if size_matches and sha256 and (not expected_md5 or publisher_checksum_matches)
                    else "missing_or_checksum_mismatch"
                ),
            }
        )
    write_tsv(args.artifacts_dir / "source_validation.tsv", rows)


def _load_finngen_manifest(path: Path) -> dict[str, dict[str, object]]:
    wanted = {"M13_INTERVERTEB", "M13_LOWBACKPAIN", "M13_SCIATICA"}
    result: dict[str, dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["phenocode"] in wanted:
                result[row["phenocode"]] = {
                    "phenotype": row["phenotype"],
                    "cases": int(row["num_cases"]),
                    "controls": int(row["num_controls"]),
                    "url": row["path_https"],
                    "genome_build": "GRCh38",
                }
    missing = wanted - result.keys()
    if missing:
        raise ValueError(f"FinnGen manifest lacks outcomes: {sorted(missing)}")
    return result


def split_rsids(value: str) -> set[str]:
    normalized = value.replace(";", ",").replace("|", ",").replace(" ", ",")
    return {token for token in normalized.split(",") if token.startswith("rs")}


def load_rows(path: Path) -> list[dict[str, str]]:
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def group_best(
    rows: Iterable[dict[str, str]], group_column: str, key_column: str, p_column: str
) -> dict[str, dict[str, dict[str, str]]]:
    result: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        group = row[group_column]
        key = row[key_column]
        current = result[group].get(key)
        if current is None or float(row[p_column]) < float(current[p_column]):
            result[group][key] = row
    return result


def best_by(rows: Iterable[Mapping[str, str]], key_column: str, p_column: str) -> dict[str, Mapping[str, str]]:
    result: dict[str, Mapping[str, str]] = {}
    for row in rows:
        key = str(row[key_column])
        if key not in result or float(row[p_column]) < float(result[key][p_column]):
            result[key] = row
    return result


def _collect_column(path: Path, column: str) -> set[str]:
    with open_text(path, "rt") as handle:
        return {row[column] for row in csv.DictReader(handle, delimiter="\t")}


def _count_by(path: Path, column: str) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    with open_text(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts[row[column]] += 1
    return counts


def _require_complete_gzip(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    with path.open("rb") as handle:
        if handle.read(2) != b"\x1f\x8b":
            raise ValueError(f"not a gzip stream: {path}")


class atomic_gzip_writer:
    def __init__(self, path: Path):
        self.path = path
        self.temp = path.with_name(path.name + ".tmp")
        self.handle = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = gzip.open(self.temp, "wt", encoding="utf-8", newline="")
        return self.handle

    def __exit__(self, exc_type, exc, traceback):
        assert self.handle is not None
        self.handle.close()
        if exc_type is None:
            self.temp.replace(self.path)
        elif self.temp.exists():
            self.temp.unlink()


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]]) -> None:
    data = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not data:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in data:
        for key in row:
            if key not in fields:
                fields.append(key)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in data:
            writer.writerow({key: format_value(row.get(key)) for key in fields})
    temp.replace(path)


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n", encoding="utf-8")
    temp.replace(path)


def format_value(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return format(value, ".12g")
    if isinstance(value, bool):
        return "true" if value else "false"
    return value


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    return hash_file(path, "sha256", block_size)


def hash_file(path: Path, algorithm: str, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
