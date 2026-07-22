"""Run official SMR/HEIDI for the IDD candidate genes and FinnGen outcomes."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import shutil
import subprocess
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd

from idd_genetics_utils import benjamini_hochberg


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=ROOT / "artifacts" / "mr" / "formal_candidate_genes.tsv")
    parser.add_argument("--manifest", type=Path, default=ROOT / "data" / "raw" / "idd" / "genetics" / "finngen" / "finngen_R13_manifest.tsv")
    parser.add_argument("--candidate-dir", type=Path, default=ROOT / "data" / "processed" / "idd" / "genetics")
    parser.add_argument("--besd-archive", type=Path, default=ROOT / "data" / "raw" / "idd" / "genetics" / "eqtlgen" / "eqtlgen_smr_besd.tar.gz")
    parser.add_argument("--besd-dir", type=Path, default=ROOT / "data" / "processed" / "idd" / "genetics" / "smr_besd")
    parser.add_argument("--ld-prefix", type=Path, default=ROOT / "data" / "raw" / "idd" / "genetics" / "ld" / "1000G_EUR")
    parser.add_argument("--smr", type=Path, default=ROOT / "tools" / "idd" / "smr" / "smr-1.3.1-win-x86_64" / "smr-1.3.1-win.exe")
    parser.add_argument("--work-dir", type=Path, default=ROOT / "data" / "processed" / "idd" / "genetics" / "smr_runs")
    parser.add_argument("--output", type=Path, default=ROOT / "artifacts" / "mr" / "smr_heidi_results.tsv")
    parser.add_argument("--peqtl-smr", type=float, default=5e-8)
    parser.add_argument("--heidi-min-m", type=int, default=10)
    parser.add_argument("--reuse-ma", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    candidates = read_candidates(args.candidates)
    outcomes = read_manifest(args.manifest)
    require_inputs(args)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    besd_prefix = ensure_besd(args.besd_archive, args.besd_dir)
    probe_map = map_candidate_probes(prefix_file(besd_prefix, ".epi"), candidates)
    results: list[dict[str, object]] = []
    command_records = []

    for outcome_id, metadata in outcomes.items():
        extract_path = args.candidate_dir / f"{outcome_id}.candidate_variants.tsv.gz"
        ma_path = args.work_dir / f"{outcome_id}.ma"
        if not args.reuse_ma or not ma_path.exists() or ma_path.stat().st_size == 0:
            write_ma(extract_path, ma_path, int(metadata["cases"]) + int(metadata["controls"]))
        probes_path = args.work_dir / f"{outcome_id}.probes.txt"
        exact_probes = sorted(set(probe_map.values()))
        probes_path.write_text("\n".join(exact_probes) + ("\n" if exact_probes else ""), encoding="ascii")
        output_prefix = args.work_dir / outcome_id
        previous_output = output_prefix.with_suffix(".smr")
        if previous_output.exists():
            previous_output.unlink()
        command = [
            str(args.smr),
            "--bfile",
            str(args.ld_prefix),
            "--gwas-summary",
            str(ma_path),
            "--beqtl-summary",
            str(besd_prefix),
            "--extract-probe",
            str(probes_path),
            "--smr",
            "--peqtl-smr",
            str(args.peqtl_smr),
            "--heidi-min-m",
            str(args.heidi_min_m),
            "--out",
            str(output_prefix),
        ]
        completed = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace")
        (args.work_dir / f"{outcome_id}.stdout.log").write_text(completed.stdout, encoding="utf-8")
        (args.work_dir / f"{outcome_id}.stderr.log").write_text(completed.stderr, encoding="utf-8")
        command_records.append(
            {
                "outcome_id": outcome_id,
                "command": command,
                "returncode": completed.returncode,
            }
        )
        parsed = parse_smr_output(output_prefix.with_suffix(".smr")) if completed.returncode == 0 else {}
        for candidate in candidates:
            gene = candidate["ensembl_gene_id"]
            exact_probe = probe_map.get(gene)
            base: dict[str, object] = {
                **candidate,
                "outcome_id": outcome_id,
                "probe_id": exact_probe or "",
                "top_snp": "",
                "effect_allele": "",
                "other_allele": "",
                "eaf": "",
                "beta_gwas": "",
                "se_gwas": "",
                "p_gwas": "",
                "beta_eqtl": "",
                "se_eqtl": "",
                "p_eqtl": "",
                "beta_smr": "",
                "se_smr": "",
                "p_smr": "",
                "p_smr_fdr": "",
                "p_smr_fdr_scope": "formal_candidate_cis_regions_by_outcome",
                "p_heidi": "",
                "p_heidi_fdr": "",
                "p_heidi_fdr_reason": "HEIDI_is_a_heterogeneity_null_not_a_target_discovery_pvalue",
                "n_snp_heidi": "",
                "heidi_pass_p_gt_0_01": "",
                "status": "not_estimable",
                "reason": "",
            }
            if exact_probe is None:
                base["reason"] = "probe_absent_from_official_besd"
            elif completed.returncode != 0:
                base["reason"] = f"smr_failed_{completed.returncode}"
            elif exact_probe not in parsed:
                base["reason"] = "no_smr_result_after_official_tool_run"
            else:
                record = parsed[exact_probe]
                base.update(record)
                p_heidi = optional_float(record.get("p_heidi"))
                base["heidi_pass_p_gt_0_01"] = p_heidi > 0.01 if p_heidi is not None else ""
                base["status"] = "estimable"
                base["reason"] = ""
            results.append(base)

    apply_smr_fdr(results)
    write_tsv(args.output, results)
    analysis_config = json.loads((ROOT / "config" / "idd_analysis.json").read_text(encoding="utf-8"))
    portable_commands = [
        {
            **record,
            "command": [portable_path_token(token) for token in record["command"]],
        }
        for record in command_records
    ]
    metadata_path = args.output.with_name("smr_run_metadata.json")
    metadata_path.write_text(
        json.dumps(
            {
                "analysis_registered_at": analysis_config["registered_at"],
                "smr_executable": portable_path_token(str(args.smr)),
                "besd_prefix": portable_path_token(str(besd_prefix)),
                "ld_prefix": portable_path_token(str(args.ld_prefix)),
                "peqtl_smr": args.peqtl_smr,
                "heidi_min_m": args.heidi_min_m,
                "gwas_input_scope": "formal_candidate_cis_variants_matched_by_rsid_to_eqtlgen_full_cis",
                "commands": portable_commands,
            },
            indent=2,
            ensure_ascii=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


def portable_path_token(value: str) -> str:
    path = Path(value)
    if not path.is_absolute():
        return value
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def read_candidates(path: Path) -> list[dict[str, str]]:
    if path.suffix.lower() == ".xlsx":
        frame = pd.read_excel(path, sheet_name="Data", dtype=str)
        rows = []
        seen = set()
        for row in frame.to_dict(orient="records"):
            gene = clean_string(row.get("ensembl_gene_id")).split(".")[0]
            symbol = clean_string(row.get("hgnc_names"))
            if gene and symbol and gene not in seen:
                seen.add(gene)
                rows.append({"gene_symbol": symbol, "ensembl_gene_id": gene})
        return rows
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = []
    seen = set()
    for row in rows:
        gene = row["ensembl_gene_id"].split(".")[0]
        if gene not in seen:
            seen.add(gene)
            result.append({"gene_symbol": row["gene_symbol"], "ensembl_gene_id": gene})
    return result


def read_manifest(path: Path) -> dict[str, dict[str, object]]:
    wanted = {"M13_INTERVERTEB", "M13_LOWBACKPAIN", "M13_SCIATICA"}
    result = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["phenocode"] in wanted:
                result[row["phenocode"]] = {
                    "cases": int(row["num_cases"]),
                    "controls": int(row["num_controls"]),
                    "phenotype": row["phenotype"],
                }
    if result.keys() != wanted:
        raise ValueError(f"manifest outcomes mismatch: {sorted(result)}")
    return result


def require_inputs(args: argparse.Namespace) -> None:
    required = [args.smr, args.besd_archive]
    required.extend(args.ld_prefix.with_suffix(ext) for ext in (".bed", ".bim", ".fam"))
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("missing SMR inputs: " + ", ".join(missing))


def ensure_besd(archive: Path, destination: Path) -> Path:
    destination.mkdir(parents=True, exist_ok=True)
    prefixes = discover_besd_prefixes(destination)
    if not prefixes:
        with tarfile.open(archive, "r:gz") as handle:
            handle.extractall(destination)
        for compressed in destination.rglob("*.gz"):
            if compressed.name.endswith((".besd.gz", ".esi.gz", ".epi.gz")):
                decompress_gzip(compressed, compressed.with_suffix(""))
        prefixes = discover_besd_prefixes(destination)
    if not prefixes:
        raise FileNotFoundError("official archive did not contain a complete .besd/.esi/.epi prefix")
    if len(prefixes) > 1:
        prefixes.sort(key=lambda path: prefix_file(path, ".besd").stat().st_size, reverse=True)
    return prefixes[0]


def decompress_gzip(source: Path, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size > 0:
        return
    temp = destination.with_name(destination.name + ".tmp")
    with gzip.open(source, "rb") as input_handle, temp.open("wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle, length=8 * 1024 * 1024)
    temp.replace(destination)


def discover_besd_prefixes(root: Path) -> list[Path]:
    prefixes = []
    for besd in root.rglob("*.besd"):
        prefix = Path(str(besd)[: -len(".besd")])
        if prefix_file(prefix, ".esi").exists() and prefix_file(prefix, ".epi").exists():
            prefixes.append(prefix)
    return prefixes


def prefix_file(prefix: Path, extension: str) -> Path:
    return Path(str(prefix) + extension)


def map_candidate_probes(epi_path: Path, candidates: Iterable[Mapping[str, str]]) -> dict[str, str]:
    candidates = list(candidates)
    wanted = {row["ensembl_gene_id"] for row in candidates}
    token_to_gene = {row["ensembl_gene_id"]: row["ensembl_gene_id"] for row in candidates}
    token_to_gene.update({row["gene_symbol"].upper(): row["ensembl_gene_id"] for row in candidates})
    result: dict[str, str] = {}
    with epi_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                continue
            probe = fields[1]
            tokens = [probe, *fields[4:]]
            for token in tokens:
                versionless = token.split(".")[0]
                gene = token_to_gene.get(versionless) or token_to_gene.get(versionless.upper())
                if gene in wanted and gene not in result:
                    result[gene] = probe
    return result


def write_ma(source: Path, destination: Path, total_n: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    best: dict[str, dict[str, str]] = {}
    with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            snp = row["snp"]
            if snp not in best or float(row["pvalue"]) < float(best[snp]["pvalue"]):
                best[snp] = row
    temp = destination.with_name(destination.name + ".tmp")
    with temp.open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "A1", "A2", "freq", "b", "se", "p", "N"])
        for snp in sorted(best):
            row = best[snp]
            writer.writerow(
                [
                    snp,
                    row["effect_allele"],
                    row["other_allele"],
                    row["eaf"],
                    row["beta"],
                    row["se"],
                    row["pvalue"],
                    total_n,
                ]
            )
    temp.replace(destination)


def write_finngen_ma(source: Path, destination: Path, total_n: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temp = destination.with_name(destination.name + ".tmp")
    with gzip.open(source, "rt", encoding="utf-8", newline="") as input_handle, temp.open(
        "w", encoding="ascii", newline=""
    ) as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "A1", "A2", "freq", "b", "se", "p", "N"])
        header = input_handle.readline().rstrip("\r\n").split("\t")
        index = {name: position for position, name in enumerate(header)}
        seen = set()
        for line in input_handle:
            row = line.rstrip("\r\n").split("\t")
            rsids = [token for token in row[index["rsids"]].replace(";", ",").split(",") if token.startswith("rs")]
            for snp in rsids:
                if snp in seen:
                    continue
                seen.add(snp)
                writer.writerow(
                    [
                        snp,
                        row[index["alt"]],
                        row[index["ref"]],
                        row[index["af_alt"]],
                        row[index["beta"]],
                        row[index["sebeta"]],
                        row[index["pval"]],
                        total_n,
                    ]
                )
    temp.replace(destination)


def parse_smr_output(path: Path) -> dict[str, dict[str, object]]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {}
    aliases = {
        "top_snp": ("topSNP", "top_SNP"),
        "effect_allele": ("A1",),
        "other_allele": ("A2",),
        "eaf": ("Freq", "freq"),
        "beta_gwas": ("b_GWAS",),
        "se_gwas": ("se_GWAS",),
        "p_gwas": ("p_GWAS",),
        "beta_eqtl": ("b_eQTL",),
        "se_eqtl": ("se_eQTL",),
        "p_eqtl": ("p_eQTL",),
        "beta_smr": ("b_SMR",),
        "se_smr": ("se_SMR",),
        "p_smr": ("p_SMR",),
        "p_heidi": ("p_HEIDI",),
        "n_snp_heidi": ("nsnp_HEIDI",),
    }
    for row in rows:
        probe = row.get("ProbeID") or row.get("probeID") or row.get("Probe")
        if not probe:
            continue
        normalized = {}
        for target, source_names in aliases.items():
            normalized[target] = next((row[name] for name in source_names if name in row), "")
        result[probe] = normalized
    return result


def optional_float(value: object) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def apply_smr_fdr(rows: list[dict[str, object]]) -> None:
    outcomes = sorted({str(row["outcome_id"]) for row in rows})
    for outcome in outcomes:
        indices = [index for index, row in enumerate(rows) if row["outcome_id"] == outcome]
        pvalues = [optional_float(rows[index].get("p_smr")) for index in indices]
        for index, adjusted in zip(indices, benjamini_hochberg(pvalues)):
            rows[index]["p_smr_fdr"] = "" if adjusted is None else adjusted


def clean_string(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def write_tsv(path: Path, rows: list[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    for row in rows[1:]:
        for key in row:
            if key not in fields:
                fields.append(key)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


if __name__ == "__main__":
    raise SystemExit(main())
