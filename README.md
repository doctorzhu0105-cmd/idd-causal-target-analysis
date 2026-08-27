# Correlated cis-MR analysis for intervertebral disc disease

This repository is the public reproducibility and code layer synchronized to the locked scientific MASTER 1.4.0. It implements configuration-dependent cis-Mendelian randomization with correlated instruments and uses the intervertebral-disc-disease screen as an empirical test bed. The manuscript and its internal governance tree remain authoritative outside this repository.

## Current reproducibility layer

The MASTER 1.4 code and public outputs cover:

- signed-LD generalized IVW with NumPy's Moore-Penrose pseudoinverse;
- generalized Cochran Q and multiplicative random-effects standard errors;
- exposure-increasing GLS MR-Egger orientation;
- leave-one-out refits using the reduced covariance matrix without re-pruning;
- LD-threshold sweeps, effective dimension, leverage and joint exposure information;
- global 7,284-row BH as the formal multiplicity reference, with endpoint, BY and K-family sensitivities;
- the P2RY13 multiplicity-family example, which is not rescued by the current formal families;
- 200-replicate same-panel LD-reference bootstrap, AKR1C2 re-pruning, matrix/missingness QC and shrinkage sensitivity;
- colocalization-prior sensitivity;
- the current neutral 13-row evidence matrix, Supplementary Figure S2 source and Figure 1 source/generator.

Small, redistributable outputs are in `results/master_1_4/`. Their identities and upstream source hashes are recorded in `provenance/`.

## Boundaries

- This is not a manuscript authority tree and does not assign discovery tiers or validated-target labels.
- The 13-row table is a K>=2 sensitivity subset with downstream evidence, not the formal multiplicity family.
- The LD bootstrap resamples the same 503-person 1000 Genomes EUR panel. It measures reference-sample uncertainty conditional on the available panel, selected instruments and summary statistics; it is not independent-panel replication.
- A known-truth simulation was deferred and is not represented as completed.
- Raw eQTLGen, FinnGen, 1000 Genomes PLINK files, UKB caches, participant-level data and redistribution-uncertain extracts are not included.
- Historical v1.0.0 modules remain for reproducibility of the earlier public workflow. They are secondary to the current correlated cis-MR layer and have not been rewritten or deleted.

## Quick start

```bash
conda env create -f environment/environment.yml
conda activate idd-causal-target-analysis
python -B -m unittest discover -s tests -p "test_*.py"
python tools/verify_package.py --root .
```

To render the current public assets outside the repository tree:

```bash
python analysis_p0_correlated_mr_20260814/scripts/make_current_figure1.py --repo-root . --output-dir ../idd-master14-render
python py/generate_master14_assets.py --repo-root . --output-dir ../idd-master14-render
```

Full numerical reruns require the external inputs documented in `docs/DATA_ACCESS.md`. See `docs/REPRODUCIBILITY.md` for execution order and explicit reproducibility limits.

## Candidate status

This branch is an unreleased MASTER 1.4.0 synchronization candidate. No tag, GitHub Release or new Zenodo version has been created. The prior public baseline remains commit `67da2c560d01dc24f1b9d3769b30ede505df931f` / tag `v1.0.0`; this candidate is an ordinary forward update and does not rewrite that history.
