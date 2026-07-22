# Reproducibility guide

## 1. Environments

Create the Python environment:

```bash
conda env create -f environment/environment.yml
conda activate idd-causal-target-analysis
```

Restore R packages from the recorded lockfile:

```r
install.packages("renv")
renv::restore(lockfile = "environment/renv.lock", prompt = FALSE)
```

The original external tools were PLINK v1.9.0-b.7.11 (19 Aug 2025) and SMR v1.3.1. They are not redistributed. The Python commands accept explicit `--plink` or `--smr` paths where supported; some diagnostic R scripts retain the original Windows directory contract.

## 2. Offline checks

```bash
make test
make package-verify
```

These checks use synthetic fixtures and package metadata. They do not claim that absent source data have been reproduced.

## 3. Prepared-input analysis stages

After placing the registered inputs at the paths in `docs/DATA_ACCESS.md` and the JSON source registries, run stages explicitly:

```bash
make prepare-druggable
make gate0
make bulk
make genetics-validate
make genetics-screen
make smr
make scrna
make drug
make rank
make gtex
make docking
make handoff
```

The repository deliberately does not expose a misleading one-command fresh-clone target. Several inputs are multi-gigabyte, source-specific, or require external executables.

## 4. Derived and sensitivity analyses

With the upstream results present, the corrected or article-facing derived analyses are available as separate targets:

```bash
make wgcna-signed
make ukb-outcome
make mr-presso
make official-coloc
make coloc-susie
make tables-figures
```

`make article-verify` becomes applicable only after the manuscript assets are added. They are intentionally absent from this code-only upload candidate.
