# IDD causal target analysis

This repository contains the code-only reproducibility package for a translational bioinformatics workflow that prioritizes druggable and locally deliverable targets for intervertebral disc disorders.

## Scope

Included code covers bulk transcriptomics, exploratory WGCNA, druggable-genome cis-eQTL Mendelian randomization, colocalization, official SMR/HEIDI, Steiger and regional direction analyses, single-cell localization, drug-target evidence, docking applicability, target ranking, and deterministic table/figure generation.

This snapshot intentionally excludes raw data, processed source extracts, participant-level data, API-response caches, external executables, derived result tables, manuscript text, author contact information, and submission files. Source accessions, URLs, expected file metadata, analysis parameters, and environment locks are retained in `config/` and `environment/`.

## Quick start

```bash
conda env create -f environment/environment.yml
conda activate idd-causal-target-analysis
python -m unittest discover -s tests -p "test_*.py"
python tools/verify_package.py --root .
```

Restore the R environment from `environment/renv.lock` before running R analyses. See `docs/REPRODUCIBILITY.md` for staged commands and `docs/DATA_ACCESS.md` for the required external inputs.

## Important boundaries

- The package is code-only; a fresh clone cannot reproduce the numerical results until the registered public inputs and external tools are installed at the documented paths.
- pQTL analyses were deferred and are not represented as completed.
- The single-cell workflow supports descriptive localization and does not, by itself, establish a disease effect.
- Docking is gated by target identity, intervention direction, human experimental structure, and ligand availability; `NOT_APPLICABLE` is an allowed outcome.
- The coloc-SuSiE diagnostic preserves non-convergence rather than fabricating posterior probabilities.

## Repository status

The code is released under the MIT License. The stable repository URL, release tag, archival DOI, and article citation remain pending until real deposition; they must not be invented in advance.
