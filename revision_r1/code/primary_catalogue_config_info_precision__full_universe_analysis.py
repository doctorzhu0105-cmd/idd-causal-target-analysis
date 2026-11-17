"""Full Finan-universe LD-aware MR, adequacy census, and locked diagnostics.

This post-hoc extension is deliberately isolated from the frozen core and all
submission packages.  It uses the raw eQTLGen full-cis file, the raw FinnGen
R13 files, the existing allele harmonisation, LD reader, pruning, and correlated
IVW modules.  It writes only files whose names begin with ``full_universe`` (or
the explicitly named M5/M12 outputs) and supports resumable gene batches.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
WS = ROOT.parent
OUT = ROOT / "out"
RAW_DIR = WS / "data" / "raw" / "idd" / "genetics"
RAW_EQTL = RAW_DIR / "eqtlgen" / "eqtlgen_full_cis.txt.gz"
RAW_AF = RAW_DIR / "eqtlgen" / "eqtlgen_allele_frequency.txt.gz"
FINAN = WS / "data" / "raw" / "mr" / "finan_2017_table_s1.xlsx"
FINNGEN_DIR = RAW_DIR / "finngen"
LD_DIR = RAW_DIR / "ld"
BED = LD_DIR / "1000G_EUR.bed"
BIM = LD_DIR / "1000G_EUR.bim"
FEASIBILITY = OUT / "feasibility_gate_full_finan_20260816.json"
PREP_PROVENANCE = OUT / "full_universe_prepare_provenance.json"
MR_CHECKPOINT = OUT / "full_universe_mr_checkpoint.json"
MR_BATCH_DIR = OUT / "full_universe_mr_batches"
EQTL_BATCH_DIR = OUT / "full_universe_eqtl_batches"
ADEQUACY_BATCH_DIR = OUT / "full_universe_adequacy_batches"

P_THRESHOLD = 5e-8
F_THRESHOLD = 10.0
N_FRACTION = 0.80
R2_PRIMARY = 0.1
R2_LADDER = [0.001, 0.01, 0.05, 0.1, 0.2]
OUTCOMES = ["M13_INTERVERTEB", "M13_LOWBACKPAIN", "M13_SCIATICA"]
MAIN_GENES = ["NT5E", "AKR1C1", "AKR1C2", "P2RY13", "CCL28"]
N_REF = 503
BATCH_SIZE = 100
SMR_BONFERRONI_N = 71

EQTL_COLUMNS = [
    "gene_symbol", "ensembl_gene_id", "snp", "chr_b37", "pos_b37",
    "gene_pos_b37", "effect_allele", "other_allele", "zscore", "n",
    "pvalue", "fdr",
]
ENRICHED_COLUMNS = EQTL_COLUMNS + ["eaf", "beta", "se", "f_statistic", "af_status"]
OUTCOME_COLUMNS = ["snp", "fg_ea", "fg_oa", "fg_beta", "fg_se"]

from project_modules.correlated_mr import egger_correlated, ivw_correlated, prune_by_ld  # noqa: E402
from project_modules.idd_genetics_pipeline import split_rsids  # noqa: E402
from project_modules.idd_genetics_utils import standardized_effect_from_z  # noqa: E402
from project_modules.plink_ld import is_palindromic, ld_matrix, load_bim  # noqa: E402
from project_modules.run_all_genes_endpoints import harmonise  # noqa: E402
from project_modules import global_ablation
from project_modules.global_ablation import concentration_index, _ld_pruned_exposure  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_eligible() -> tuple[list[str], dict[str, dict[str, str]]]:
    if not FEASIBILITY.exists():
        raise FileNotFoundError(f"Missing feasibility gate: {FEASIBILITY}")
    payload = json.loads(FEASIBILITY.read_text(encoding="utf-8"))
    eligible = list(payload["finan_universe"]["eligible_gene_ids"])
    if len(eligible) != int(payload["finan_universe"]["eligible_genes_with_at_least_one_filtered_instrument"]):
        raise AssertionError("feasibility eligible list/count mismatch")
    from openpyxl import load_workbook

    workbook = load_workbook(FINAN, read_only=True, data_only=True)
    sheet = workbook["Data"]
    iterator = sheet.iter_rows(values_only=True)
    header = list(next(iterator))
    i_id = header.index("ensembl_gene_id")
    i_symbol = header.index("hgnc_names")
    all_genes: dict[str, dict[str, str]] = {}
    for row in iterator:
        gene_id = str(row[i_id]).split(".")[0]
        all_genes[gene_id] = {
            "ensembl_gene_id": gene_id,
            "finan_symbol": "" if row[i_symbol] is None else str(row[i_symbol]).strip(),
        }
    workbook.close()
    if len(all_genes) != 4479:
        raise AssertionError(f"expected 4,479 Finan rows, found {len(all_genes)}")
    if not set(eligible).issubset(all_genes):
        raise AssertionError("eligible list contains an Ensembl ID absent from Finan")
    return sorted(eligible), all_genes


def batch_map(eligible: list[str]) -> dict[str, int]:
    return {gene: i // BATCH_SIZE for i, gene in enumerate(eligible)}


def n_batches(eligible: list[str]) -> int:
    return (len(eligible) + BATCH_SIZE - 1) // BATCH_SIZE


def _raw_indices(header: list[str]) -> dict[str, int]:
    index = {name: i for i, name in enumerate(header)}
    required = {
        "Gene", "GeneSymbol", "Pvalue", "Zscore", "SNPPos", "GenePos", "SNP",
        "SNPChr", "AssessedAllele", "OtherAllele", "NrSamples", "FDR",
    }
    missing = required - set(index)
    if missing:
        raise ValueError(f"raw eQTLGen header missing {sorted(missing)}")
    return index


def extract_raw_subset(eligible: set[str], destination: Path) -> tuple[dict[str, float], int, int]:
    max_n: dict[str, float] = {}
    rows = 0
    genes: set[str] = set()
    partial = destination.with_suffix(destination.suffix + ".partial")
    with gzip.open(RAW_EQTL, "rt", encoding="utf-8", newline="") as source, gzip.open(partial, "wt", encoding="utf-8", newline="") as out:
        header = source.readline().rstrip("\r\n").split("\t")
        index = _raw_indices(header)
        writer = csv.DictWriter(out, fieldnames=EQTL_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line in source:
            values = line.rstrip("\r\n").split("\t")
            gene_id = values[index["Gene"]].split(".")[0]
            if gene_id not in eligible:
                continue
            n = float(values[index["NrSamples"]])
            max_n[gene_id] = max(max_n.get(gene_id, float("-inf")), n)
            genes.add(gene_id)
            rows += 1
            writer.writerow({
                "gene_symbol": values[index["GeneSymbol"]] or gene_id, "ensembl_gene_id": gene_id,
                "snp": values[index["SNP"]], "chr_b37": values[index["SNPChr"]],
                "pos_b37": values[index["SNPPos"]], "gene_pos_b37": values[index["GenePos"]],
                "effect_allele": values[index["AssessedAllele"]].upper(),
                "other_allele": values[index["OtherAllele"]].upper(),
                "zscore": values[index["Zscore"]], "n": values[index["NrSamples"]],
                "pvalue": values[index["Pvalue"]], "fdr": values[index["FDR"]],
            })
    partial.replace(destination)
    return max_n, rows, len(genes)


def filter_raw_subset(source_path: Path, max_n: dict[str, float], destination: Path) -> tuple[set[str], dict[str, int]]:
    requested: set[str] = set()
    counts: dict[str, int] = {}
    partial = destination.with_suffix(destination.suffix + ".partial")
    with gzip.open(source_path, "rt", encoding="utf-8", newline="") as source, gzip.open(partial, "wt", encoding="utf-8", newline="") as out:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(out, fieldnames=EQTL_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            gene = row["ensembl_gene_id"]
            p = float(row["pvalue"])
            z = float(row["zscore"])
            n = float(row["n"])
            if not (p < P_THRESHOLD and z * z > F_THRESHOLD and n >= N_FRACTION * max_n[gene]):
                continue
            if is_palindromic(row["effect_allele"], row["other_allele"]):
                continue
            writer.writerow(row)
            requested.add(row["snp"])
            counts[gene] = counts.get(gene, 0) + 1
    partial.replace(destination)
    return requested, counts


def read_af_subset(requested: set[str], destination: Path) -> dict[str, tuple[str, str, float]]:
    af_map: dict[str, tuple[str, str, float]] = {}
    partial = destination.with_suffix(destination.suffix + ".partial")
    fields = ["snp", "allele_a", "allele_b", "allele_b_frequency"]
    with gzip.open(RAW_AF, "rt", encoding="utf-8", newline="") as source, gzip.open(partial, "wt", encoding="utf-8", newline="") as out:
        header = source.readline().rstrip("\r\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line in source:
            values = line.rstrip("\r\n").split("\t")
            snp = values[index["SNP"]]
            if snp not in requested:
                continue
            record = {
                "snp": snp, "allele_a": values[index["AlleleA"]].upper(),
                "allele_b": values[index["AlleleB"]].upper(),
                "allele_b_frequency": values[index["AlleleB_all"]],
            }
            writer.writerow(record)
            af_map[snp] = (record["allele_a"], record["allele_b"], float(record["allele_b_frequency"]))
    partial.replace(destination)
    return af_map


def enrich_and_split(
    filtered_path: Path,
    af_map: dict[str, tuple[str, str, float]],
    gene_to_batch: dict[str, int],
    batch_count: int,
) -> tuple[int, int]:
    enriched_path = OUT / "full_universe_eqtl.tsv.gz"
    partial = enriched_path.with_suffix(enriched_path.suffix + ".partial")
    EQTL_BATCH_DIR.mkdir(parents=True, exist_ok=True)
    handles = []
    writers = []
    try:
        for batch in range(batch_count):
            path = EQTL_BATCH_DIR / f"eqtl_batch_{batch:03d}.tsv.gz"
            handle = gzip.open(path.with_suffix(path.suffix + ".partial"), "wt", encoding="utf-8", newline="")
            handles.append((path, handle))
            writer = csv.DictWriter(handle, fieldnames=ENRICHED_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writers.append(writer)
        total = converted = 0
        with gzip.open(filtered_path, "rt", encoding="utf-8", newline="") as source, gzip.open(partial, "wt", encoding="utf-8", newline="") as out:
            reader = csv.DictReader(source, delimiter="\t")
            writer = csv.DictWriter(out, fieldnames=ENRICHED_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                total += 1
                eaf = beta = se = None
                status = "missing"
                af = af_map.get(row["snp"])
                if af:
                    allele_a, allele_b, freq_b = af
                    if row["effect_allele"] == allele_b and row["other_allele"] == allele_a:
                        eaf = freq_b; status = "matched_allele_b"
                    elif row["effect_allele"] == allele_a and row["other_allele"] == allele_b:
                        eaf = 1.0 - freq_b; status = "matched_allele_a"
                    else:
                        status = "allele_mismatch"
                if eaf is not None and 0 < eaf < 1:
                    beta, se = standardized_effect_from_z(float(row["zscore"]), float(eaf), float(row["n"]))
                    converted += 1
                enriched = {
                    **row, "eaf": "" if eaf is None else eaf,
                    "beta": "" if beta is None else beta, "se": "" if se is None else se,
                    "f_statistic": float(row["zscore"]) ** 2, "af_status": status,
                }
                writer.writerow(enriched)
                writers[gene_to_batch[row["ensembl_gene_id"]]].writerow(enriched)
        return total, converted
    finally:
        for path, handle in handles:
            handle.close()
            partial_path = path.with_suffix(path.suffix + ".partial")
            partial_path.replace(path)
        partial.replace(enriched_path)


def extract_outcomes(requested: set[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for outcome in OUTCOMES:
        destination = OUT / f"full_universe_{outcome}.tsv.gz"
        partial = destination.with_suffix(destination.suffix + ".partial")
        count = 0
        with gzip.open(FINNGEN_DIR / f"finngen_R13_{outcome}.gz", "rt", encoding="utf-8", newline="") as source, gzip.open(partial, "wt", encoding="utf-8", newline="") as out:
            header = source.readline().rstrip("\r\n").split("\t")
            index = {name: i for i, name in enumerate(header)}
            writer = csv.DictWriter(out, fieldnames=OUTCOME_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for line in source:
                values = line.rstrip("\r\n").split("\t")
                for snp in split_rsids(values[index["rsids"]]) & requested:
                    writer.writerow({
                        "snp": snp, "fg_ea": values[index["alt"]].upper(),
                        "fg_oa": values[index["ref"]].upper(),
                        "fg_beta": values[index["beta"]], "fg_se": values[index["sebeta"]],
                    })
                    count += 1
        partial.replace(destination)
        counts[outcome] = count
    return counts


def prepare() -> None:
    eligible, finan = load_eligible()
    gene_to_batch = batch_map(eligible)
    batch_count = n_batches(eligible)
    raw_path = OUT / "full_universe_eqtl.raw.tsv.gz"
    filtered_path = OUT / "full_universe_eqtl_filtered.raw.tsv.gz"
    af_path = OUT / "full_universe_allele_frequency.tsv.gz"
    print(f"Preparing {len(eligible):,} eligible genes in {batch_count} batches", flush=True)
    max_n, raw_rows, raw_gene_count = extract_raw_subset(set(eligible), raw_path)
    requested, counts = filter_raw_subset(raw_path, max_n, filtered_path)
    af_map = read_af_subset(requested, af_path)
    total, converted = enrich_and_split(filtered_path, af_map, gene_to_batch, batch_count)
    outcome_counts = extract_outcomes(requested)
    if len(counts) != len(eligible):
        missing = sorted(set(eligible) - set(counts))
        raise AssertionError(f"eligible genes missing filtered rows: {missing[:10]}")
    payload = {
        "analysis": "full_finan_universe_prepare",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "formal_denominator": len(eligible), "eligible_gene_ids": eligible,
        "batch_size": BATCH_SIZE, "batch_count": batch_count,
        "raw_eqtl_path": str(RAW_EQTL), "processed_candidate_extract_not_used": True,
        "raw_rows_extracted": raw_rows, "raw_genes_with_rows": raw_gene_count,
        "filtered_rows": total, "filtered_genes": len(counts),
        "af_requested_snps": len(requested), "af_found_snps": len(af_map),
        "exposure_rows_with_beta_se": converted, "outcome_rows": outcome_counts,
        "filter": {"p_lt": P_THRESHOLD, "f_gt": F_THRESHOLD, "n_fraction": N_FRACTION, "exclude_palindromic": True},
        "outputs": {"eqtl_batch_dir": str(EQTL_BATCH_DIR), "mr_batch_dir": str(MR_BATCH_DIR)},
    }
    PREP_PROVENANCE.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"prepared {total:,} filtered rows; outcome rows={outcome_counts}", flush=True)


def load_outcome_extract(outcome: str) -> pd.DataFrame:
    path = OUT / f"full_universe_{outcome}.tsv.gz"
    df = pd.read_csv(path, sep="\t", compression="gzip")
    return df


def run_gene_endpoint(eqtl: pd.DataFrame, gene: str, outcome: str, fg: pd.DataFrame, bim_all: pd.DataFrame) -> dict[str, object]:
    base: dict[str, object] = {
        "gene": gene, "outcome": outcome, "r2_threshold": R2_PRIMARY,
        "k": 0, "n_harmonised": 0, "or": np.nan, "or_low": np.nan, "or_high": np.nan,
        "beta": np.nan, "se": np.nan, "pvalue": np.nan, "q_stat": np.nan, "q_pvalue": np.nan,
        "egger_intercept": np.nan, "egger_intercept_p": np.nan, "min_f": np.nan,
        "ld_condition_number": np.nan, "status": "no_harmonized_instruments", "reason": "",
    }
    harm = harmonise(eqtl, gene, fg)
    if harm.empty:
        return base
    harm = harm.dropna(subset=["beta_exposure", "se_outcome", "f_statistic"])
    base["n_harmonised"] = int(len(harm))
    if harm.empty:
        base["reason"] = "missing_beta_or_se"
        return base
    chrom = str(int(harm.chr_b37.iloc[0]))
    kept, rho = ld_matrix(harm, bim_all[bim_all.chr == chrom], str(BED), N_REF)
    if kept.empty:
        base["status"] = "no_ld_reference"
        return base
    order = np.argsort(-kept.f_statistic.to_numpy())
    idx = prune_by_ld(rho, order, R2_PRIMARY)
    if len(idx) == 0:
        base["status"] = "no_pruned_instruments"
        return base
    sub = kept.iloc[idx].reset_index(drop=True)
    rho_sub = rho[np.ix_(idx, idx)]
    bx = sub.beta_exposure.to_numpy(dtype=float)
    by = sub.beta_outcome.to_numpy(dtype=float)
    sey = sub.se_outcome.to_numpy(dtype=float)
    estimate = ivw_correlated(bx, by, sey, rho_sub, random_effects=True)
    _, intercept = egger_correlated(bx, by, sey, rho_sub)
    base.update({
        "k": int(len(sub)), "or": estimate.or_, "or_low": estimate.or_low,
        "or_high": estimate.or_high, "beta": estimate.beta, "se": estimate.se,
        "pvalue": estimate.pvalue, "q_stat": estimate.q_stat, "q_pvalue": estimate.q_pvalue,
        "egger_intercept": intercept.beta, "egger_intercept_p": intercept.pvalue,
        "min_f": float(sub.f_statistic.min()), "ld_condition_number": float(np.linalg.cond(rho_sub)),
        "status": "estimable", "reason": "",
    })
    return base


def run_mr() -> None:
    if not PREP_PROVENANCE.exists():
        raise FileNotFoundError("Run --stage prepare first")
    eligible, _ = load_eligible()
    batches = n_batches(eligible)
    MR_BATCH_DIR.mkdir(parents=True, exist_ok=True)
    if MR_CHECKPOINT.exists():
        checkpoint = json.loads(MR_CHECKPOINT.read_text(encoding="utf-8"))
        completed = set(int(x) for x in checkpoint.get("completed_batches", []))
    else:
        completed = set()
    fg = {outcome: load_outcome_extract(outcome) for outcome in OUTCOMES}
    bim_all = load_bim(str(BIM))
    for batch in range(batches):
        batch_path = MR_BATCH_DIR / f"mr_batch_{batch:03d}.tsv"
        if batch in completed:
            if not batch_path.exists():
                raise AssertionError(f"checkpoint says batch {batch} complete but file is missing")
            continue
        eqtl_path = EQTL_BATCH_DIR / f"eqtl_batch_{batch:03d}.tsv.gz"
        eqtl = pd.read_csv(eqtl_path, sep="\t", compression="gzip")
        gene_ids = eligible[batch * BATCH_SIZE:(batch + 1) * BATCH_SIZE]
        symbols = eqtl.groupby("ensembl_gene_id")["gene_symbol"].first().to_dict()
        rows: list[dict[str, object]] = []
        for gene_id in gene_ids:
            symbol = str(symbols.get(gene_id, gene_id))
            for outcome in OUTCOMES:
                rows.append(run_gene_endpoint(eqtl, symbol, outcome, fg[outcome], bim_all) | {
                    "ensembl_gene_id": gene_id,
                })
        pd.DataFrame(rows).to_csv(batch_path, sep="\t", index=False)
        completed.add(batch)
        MR_CHECKPOINT.write_text(json.dumps({
            "analysis": "full_universe_mr", "formal_denominator": len(eligible),
            "r2_threshold": R2_PRIMARY, "completed_batches": sorted(completed),
            "total_batches": batches, "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        }, indent=2) + "\n", encoding="utf-8")
        print(f"MR batch {batch + 1}/{batches} complete", flush=True)
    frames = [pd.read_csv(MR_BATCH_DIR / f"mr_batch_{batch:03d}.tsv", sep="\t") for batch in range(batches)]
    result = pd.concat(frames, ignore_index=True)
    expected = len(eligible) * len(OUTCOMES)
    if len(result) != expected or result.duplicated(["ensembl_gene_id", "outcome"]).any():
        raise AssertionError(f"MR row/key check failed: rows={len(result)}, expected={expected}")
    result.to_csv(OUT / "full_universe_mr.tsv", sep="\t", index=False)
    print(f"wrote full_universe_mr.tsv: {len(result):,} rows", flush=True)


def bh_adjust(pvalues: pd.Series) -> tuple[pd.Series, pd.Series]:
    values = pd.to_numeric(pvalues, errors="coerce")
    finite = values.notna() & np.isfinite(values)
    order = np.argsort(values[finite].to_numpy())
    sorted_p = values[finite].to_numpy()[order]
    m = len(sorted_p)
    q_sorted = np.minimum.accumulate((sorted_p * m / np.arange(1, m + 1))[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    q = pd.Series(np.nan, index=values.index, dtype=float)
    rank = pd.Series(pd.NA, index=values.index, dtype="Int64")
    finite_indices = values[finite].index.to_numpy()[order]
    q.loc[finite_indices] = q_sorted
    rank.loc[finite_indices] = np.arange(1, m + 1)
    return q, rank


def run_fdr() -> None:
    mr = pd.read_csv(OUT / "full_universe_mr.tsv", sep="\t")
    mr["bh_fdr"], mr["bh_rank"] = bh_adjust(mr["pvalue"].where(mr["status"] == "estimable"))
    mr.to_csv(OUT / "full_universe_fdr.tsv", sep="\t", index=False)
    estimable = mr[mr.status == "estimable"].copy()
    main = mr[mr.gene.isin(MAIN_GENES)].copy()
    old = pd.read_csv(OUT / "integrated_evidence.tsv", sep="\t")
    old = old.rename(columns={"gene_symbol": "gene", "outcome_id": "outcome"})
    old = old[["gene", "outcome", "bh_fdr"]].rename(columns={"bh_fdr": "bh_fdr_15_combinations"})
    main = main.merge(old, on=["gene", "outcome"], how="left")
    main = main[["gene", "ensembl_gene_id", "outcome", "pvalue", "bh_fdr", "bh_rank", "bh_fdr_15_combinations", "status"]]
    main.to_csv(OUT / "full_universe_main_gene_fdr.tsv", sep="\t", index=False)
    provenance = {
        "analysis": "full_universe_bh_fdr",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "formal_denominator_genes": int(mr.ensembl_gene_id.nunique()),
        "outcomes": OUTCOMES,
        "total_gene_endpoint_rows": int(len(mr)),
        "m": int(len(estimable)),
        "m_definition": "number of rows with status=estimable and finite primary r2<0.1 MR P values among the 2,479×3 full-universe rows",
        "bh_method": "Benjamini-Hochberg over this complete Finan-universe estimable family; non-estimable rows are not part of m",
        "main_gene_fdr_output": "out/full_universe_main_gene_fdr.tsv",
    }
    (OUT / "full_universe_fdr_provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"full-universe BH: m={len(estimable)}", flush=True)


def run_adequacy() -> None:
    eligible, _ = load_eligible()
    batches = n_batches(eligible)
    bim_all = load_bim(str(BIM))
    ADEQUACY_BATCH_DIR.mkdir(parents=True, exist_ok=True)
    for batch in range(batches):
        eqtl_path = EQTL_BATCH_DIR / f"eqtl_batch_{batch:03d}.tsv.gz"
        eqtl = pd.read_csv(eqtl_path, sep="\t", compression="gzip")
        rows: list[dict[str, object]] = []
        gene_ids = eligible[batch * BATCH_SIZE:(batch + 1) * BATCH_SIZE]
        for gene_id in gene_ids:
            gene_df = eqtl[eqtl.ensembl_gene_id == gene_id].copy()
            symbol = str(gene_df.gene_symbol.iloc[0]) if not gene_df.empty else gene_id
            for r2 in R2_LADDER:
                base = {
                    "gene": symbol, "ensembl_gene_id": gene_id, "r2_threshold": r2,
                    "n_filtered_rows": int(len(gene_df)), "k": 0, "min_f": np.nan,
                    "ld_condition_number": np.nan, "share_abs": np.nan, "share_raw": np.nan,
                    "n_negative_contrib": 0, "status": "no_exposure_conversion",
                }
                usable = gene_df.dropna(subset=["beta", "se", "f_statistic"])
                if usable.empty:
                    rows.append(base); continue
                chrom = str(int(usable.chr_b37.iloc[0]))
                sub, rho, _ = _ld_pruned_exposure(usable, bim_all[bim_all.chr == chrom], r2)
                if sub.empty:
                    base["status"] = "no_ld_reference"; rows.append(base); continue
                idx = concentration_index((sub.beta / sub.se).to_numpy(float), rho)
                base.update({"k": int(len(sub)), "min_f": float(sub.f_statistic.min()),
                             "ld_condition_number": float(np.linalg.cond(rho)), **idx, "status": "ok"})
                rows.append(base)
        pd.DataFrame(rows).to_csv(ADEQUACY_BATCH_DIR / f"adequacy_batch_{batch:03d}.tsv", sep="\t", index=False)
        print(f"adequacy batch {batch + 1}/{batches} complete", flush=True)
    frames = [pd.read_csv(ADEQUACY_BATCH_DIR / f"adequacy_batch_{batch:03d}.tsv", sep="\t") for batch in range(batches)]
    result = pd.concat(frames, ignore_index=True)
    result.to_csv(OUT / "full_universe_adequacy.tsv", sep="\t", index=False)
    distributions: dict[str, object] = {}
    for r2 in R2_LADDER:
        sub = result[result.r2_threshold == r2]
        distributions[str(r2)] = {
            "n_genes": int(sub.ensembl_gene_id.nunique()),
            "k_0": int((sub.k == 0).sum()), "k_1": int((sub.k == 1).sum()),
            "k_2": int((sub.k == 2).sum()), "k_3": int((sub.k == 3).sum()),
            "k_gt_3": int((sub.k > 3).sum()),
            "share_abs_median_k_ge_3": float(sub.loc[sub.k >= 3, "share_abs"].median()),
            "share_abs_q1_k_ge_3": float(sub.loc[sub.k >= 3, "share_abs"].quantile(.25)),
            "share_abs_q3_k_ge_3": float(sub.loc[sub.k >= 3, "share_abs"].quantile(.75)),
        }
    legacy_group = OUT / "global_ablation_group_tests.json"
    old_note = {"status": "not_available"}
    if legacy_group.exists():
        old = json.loads(legacy_group.read_text(encoding="utf-8"))
        old_note = {
            "status": "retained_but_logically_superseded",
            "old_candidate_denominator": old.get("candidate_denominator"),
            "old_comparison_denominator": old.get("comparison_denominator"),
            "old_tests": old.get("tests"),
            "comparison": "not directly comparable: the full-universe enumeration replaces the sampled comparator rather than reproducing its sampling contrast",
        }
    provenance = {
        "analysis": "full_universe_instrument_adequacy_census",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "formal_denominator": len(eligible), "complete_enumeration": True,
        "r2_ladder": R2_LADDER, "primary_metric": "share_abs",
        "metric_definition": "max absolute LD-aware quadratic-form contribution divided by the sum of absolute contributions; heuristic, not a standard IVW weight or variance share",
        "diagnostic_metrics": {"share_raw": "signed diagnostic; may be negative and is not a share decomposition", "n_negative_contrib": "count of negative LD-aware contributions"},
        "distributions": distributions,
        "legacy_sampled_comparator": old_note,
    }
    (OUT / "full_universe_adequacy_provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def effective_k(rho: np.ndarray) -> float:
    eigenvalues = np.linalg.eigvalsh(rho)
    eigenvalues = np.clip(eigenvalues, 0.0, None)
    denom = float(np.square(eigenvalues).sum())
    return float(np.square(eigenvalues.sum()) / denom) if denom > 0 else float("nan")


def run_m5() -> None:
    eligible, _ = load_eligible()
    batch_index = batch_map(eligible)
    bim_all = load_bim(str(BIM))
    fg = {outcome: load_outcome_extract(outcome) for outcome in OUTCOMES}
    rows: list[dict[str, object]] = []
    for gene in MAIN_GENES:
        matches = [gene_id for gene_id in eligible if gene_id in batch_index and gene_id]
        # Resolve the Ensembl ID from the batch files, rather than assuming a hard-coded mapping.
        gene_id = None
        for batch in sorted(set(batch_index.values())):
            path = EQTL_BATCH_DIR / f"eqtl_batch_{batch:03d}.tsv.gz"
            probe = pd.read_csv(path, sep="\t", compression="gzip", usecols=["ensembl_gene_id", "gene_symbol"])
            found = probe.loc[probe.gene_symbol == gene, "ensembl_gene_id"].unique().tolist()
            if found:
                gene_id = found[0]; break
        if gene_id is None:
            raise AssertionError(f"main gene absent from prepared eligible set: {gene}")
        eqtl = pd.read_csv(EQTL_BATCH_DIR / f"eqtl_batch_{batch_index[gene_id]:03d}.tsv.gz", sep="\t", compression="gzip")
        for outcome in OUTCOMES:
            harm = harmonise(eqtl, gene, fg[outcome])
            if harm.empty:
                for r2 in R2_LADDER:
                    rows.append({"gene": gene, "ensembl_gene_id": gene_id, "outcome": outcome, "r2_threshold": r2, "k": 0, "effective_independent_k": np.nan, "se": np.nan, "pvalue": np.nan, "status": "no_harmonized_instruments"})
                continue
            harm = harm.dropna(subset=["beta_exposure", "se_outcome", "f_statistic"])
            chrom = str(int(harm.chr_b37.iloc[0]))
            kept, rho = ld_matrix(harm, bim_all[bim_all.chr == chrom], str(BED), N_REF)
            order = np.argsort(-kept.f_statistic.to_numpy()) if not kept.empty else np.array([], dtype=int)
            for r2 in R2_LADDER:
                if kept.empty:
                    rows.append({"gene": gene, "ensembl_gene_id": gene_id, "outcome": outcome, "r2_threshold": r2, "k": 0, "effective_independent_k": np.nan, "se": np.nan, "pvalue": np.nan, "status": "no_ld_reference"}); continue
                idx = prune_by_ld(rho, order, r2)
                sub = kept.iloc[idx]
                sr = rho[np.ix_(idx, idx)]
                est = ivw_correlated(sub.beta_exposure.to_numpy(float), sub.beta_outcome.to_numpy(float), sub.se_outcome.to_numpy(float), sr, random_effects=True)
                rows.append({"gene": gene, "ensembl_gene_id": gene_id, "outcome": outcome, "r2_threshold": r2, "k": len(sub), "effective_independent_k": effective_k(sr), "se": est.se, "pvalue": est.pvalue, "or": est.or_, "ld_condition_number": float(np.linalg.cond(sr)), "status": "estimable"})
    pd.DataFrame(rows).to_csv(OUT / "effective_instrument_count.tsv", sep="\t", index=False)
    # The method is explicitly recorded in the file-level JSON so that the count is not mistaken for a physical SNP count.
    (OUT / "effective_instrument_count_provenance.json").write_text(json.dumps({
        "analysis": "effective_independent_instrument_count",
        "method": "eigenvalue participation ratio (sum eigenvalues)^2 / sum eigenvalues^2 after LD-aware greedy pruning",
        "eigenvalues": "non-negative eigenvalues of the signed 1000G EUR LD correlation submatrix; numerical negative values clipped to zero",
        "literature_basis": "effective-number-of-tests/effective-dimension eigenvalue participation ratio; reported as a diagnostic, not as a replacement for K",
        "outcomes": OUTCOMES, "r2_ladder": R2_LADDER,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def run_m12() -> None:
    core = WS / "artifacts" / "core" / "genetics"
    coloc = pd.read_csv(core / "coloc_results.tsv", sep="\t")
    smr = pd.read_csv(core / "smr_heidi_results.tsv", sep="\t")
    coloc = coloc[coloc.gene_symbol.isin(MAIN_GENES)].copy()
    if len(coloc) != 15:
        raise AssertionError(f"expected 15 coloc rows for the five main genes, found {len(coloc)}")
    smr = smr[smr.gene_symbol.isin(MAIN_GENES)][["gene_symbol", "outcome_id", "p_smr", "p_heidi"]]
    out = coloc[["ensembl_gene_id", "gene_symbol", "outcome_id", "shared_variants", "pp0", "pp1", "pp2", "pp3", "pp4", "status", "reason"]].merge(
        smr, left_on=["gene_symbol", "outcome_id"], right_on=["gene_symbol", "outcome_id"], how="left"
    )
    out["smr_bonferroni_n"] = SMR_BONFERRONI_N
    out["smr_bonferroni_family"] = "71 estimable probes per outcome in the 76-candidate screen; not 5 genes or 15 combinations"
    out["p_smr_bonferroni"] = (out.p_smr * SMR_BONFERRONI_N).clip(upper=1.0)
    out["smr_bonferroni_truncated_at_1"] = (out.p_smr * SMR_BONFERRONI_N) >= 1.0
    out.to_csv(OUT / "coloc_full_posteriors.tsv", sep="\t", index=False)
    (OUT / "coloc_smr_provenance.json").write_text(json.dumps({
        "analysis": "full_coloc_posteriors_and_smr_correction_audit",
        "coloc_source": str(core / "coloc_results.tsv"), "rows": 15,
        "posterior_columns": ["pp0", "pp1", "pp2", "pp3", "pp4"],
        "smr_correction_family": "71 estimable probes per outcome in the original 76-candidate screen",
        "smr_bonferroni_n": SMR_BONFERRONI_N,
        "not_5_or_15_combinations": True,
        "n_rows_with_bonferroni_value_clipped_to_1": int(out.smr_bonferroni_truncated_at_1.sum()),
        "p_equals_1_is_truncated": True,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["prepare", "mr", "fdr", "adequacy", "m5", "m12", "all"], required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--feasibility", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    global ROOT, WS, OUT, RAW_DIR, RAW_EQTL, RAW_AF, FINAN, FINNGEN_DIR, LD_DIR, BED, BIM, FEASIBILITY
    global PREP_PROVENANCE, MR_CHECKPOINT, MR_BATCH_DIR, EQTL_BATCH_DIR, ADEQUACY_BATCH_DIR
    args = parse_args()
    WS, OUT, FEASIBILITY = args.input_root.resolve(), args.output_dir.resolve(), args.feasibility.resolve()
    if not WS.is_dir() or not FEASIBILITY.is_file(): raise SystemExit("input root or feasibility file does not exist")
    if OUT == WS or WS in OUT.parents or OUT in WS.parents or OUT == FEASIBILITY.parent or OUT in FEASIBILITY.parents or FEASIBILITY.parent in OUT.parents:
        raise SystemExit("output directory must be disjoint from inputs")
    code_root = Path(__file__).resolve().parent
    if OUT == code_root or code_root in OUT.parents:
        raise SystemExit("output directory must be outside the code directory")
    ROOT = OUT.parent
    RAW_DIR = WS / "data/raw/idd/genetics"
    RAW_EQTL, RAW_AF = RAW_DIR / "eqtlgen/eqtlgen_full_cis.txt.gz", RAW_DIR / "eqtlgen/eqtlgen_allele_frequency.txt.gz"
    FINAN, FINNGEN_DIR, LD_DIR = WS / "data/raw/mr/finan_2017_table_s1.xlsx", RAW_DIR / "finngen", RAW_DIR / "ld"
    BED, BIM = LD_DIR / "1000G_EUR.bed", LD_DIR / "1000G_EUR.bim"
    PREP_PROVENANCE, MR_CHECKPOINT = OUT / "full_universe_prepare_provenance.json", OUT / "full_universe_mr_checkpoint.json"
    MR_BATCH_DIR, EQTL_BATCH_DIR, ADEQUACY_BATCH_DIR = OUT / "full_universe_mr_batches", OUT / "full_universe_eqtl_batches", OUT / "full_universe_adequacy_batches"
    global_ablation.BED, global_ablation.BIM = BED, BIM
    required = [FINAN]
    if args.stage in {"prepare", "all"}: required += [RAW_EQTL, RAW_AF]
    if args.stage in {"mr", "adequacy", "m5", "all"}: required += [BED, BIM]
    missing = [str(p) for p in required if not p.is_file()]
    if missing: raise SystemExit("missing required inputs: " + ", ".join(missing))
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage in {"prepare", "all"}:
        prepare()
    if args.stage in {"mr", "all"}:
        run_mr()
    if args.stage in {"fdr", "all"}:
        run_fdr()
    if args.stage in {"adequacy", "all"}:
        run_adequacy()
    if args.stage in {"m5", "all"}:
        run_m5()
    if args.stage in {"m12", "all"}:
        run_m12()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
