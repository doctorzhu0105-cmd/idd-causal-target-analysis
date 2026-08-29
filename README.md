# Correlated cis-MR analysis for intervertebral disc disease

This repository is the public code and redistributable-results layer for the scientific MASTER 1.5.0 candidate. It evaluates how LD-pruning configuration changes nominal instrument count, LD effective dimensionality, joint exposure information, generalized-IVW precision, and leverage across a 2,479-gene catalogue and three FinnGen endpoints. The manuscript governance tree remains authoritative outside this repository.

## Release layers

### v1.2.0

The current release adds the catalogue-wide configuration, information, and precision extension:

- Effective K and joint exposure information at five LD thresholds;
- a five-threshold generalized-IVW ladder;
- paired exposure-side, precision, and generalized-leverage changes;
- matched and endpoint-stratified information-precision associations;
- AKR1C1/AKR1C2 catalogue positioning;
- corrected five-threshold historical-candidate comparisons;
- K=1 second-order standard-error sensitivity with global BH/BY rank propagation;
- deterministic source data and code for the revised Figure 1.

These are descriptive configuration analyses. They are not a new estimator, estimator validation, independent replication, or therapeutic-target validation.

### v1.1.0

The previous release added the correlated cis-MR reproducibility layer: signed-LD generalized IVW, generalized Q, exposure-increasing GLS MR-Egger orientation, reduced-covariance leave-one-out fits, multiplicity-family sensitivities, same-panel LD-reference bootstrap, shrinkage, colocalization-prior sensitivity, and the neutral historical downstream evidence subset.

All v1.1.0 material remains in place. The v1.2.0 release is an ordinary forward extension and does not rewrite previous tags or releases.

## Boundaries

- This is not a manuscript authority tree and does not assign discovery tiers or validated-target labels.
- The 7,284-row global BH family remains the formal multiplicity reference.
- The catalogue extension uses provider-controlled eQTLGen, FinnGen, and 1000 Genomes inputs. Raw source files are not redistributed.
- Included gene-level tables are derived model outputs. They contain no SNP identifiers, alleles, coordinates, source-level association rows, EAF values, or participant records.
- The K=1 second-order analysis is a sensitivity analysis; it does not replace the primary first-order estimator.
- The null-construction results are diagnostics, not known-truth simulation.
- Historical v1.0.0 modules remain for earlier workflow reproducibility and are secondary to the current layer.

## Quick start

```bash
conda env create -f environment/environment.yml
conda activate idd-causal-target-analysis
python -B -m unittest discover -s tests -p "test_*.py"
python tools/verify_package.py --root .
```

Regenerate the redistributable v1.2 derivatives and Figure 1 outside the repository tree:

```bash
python analysis_p0_correlated_mr_20260814/scripts/catalogue_extension.py \
  --results-dir results/master_1_5 \
  --output-dir ../idd-master15-regenerated
python analysis_p0_correlated_mr_20260814/scripts/make_master15_figure1.py \
  --repo-root . \
  --output-dir ../idd-master15-render
```

Full numerical recomputation requires the provider inputs documented in [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md). See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the explicit input interface and verification order.
