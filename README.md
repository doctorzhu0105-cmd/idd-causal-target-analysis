# Correlated cis-MR analysis for intervertebral disc disease

The public repository contains code and redistributable-results layers for correlated cis-MR. MASTER 1.5.0 is the historical label for the v1.2.0 catalogue layer, not the R1 revision release number. It evaluates how LD-pruning configuration changes nominal instrument count, LD effective dimensionality, joint exposure information, generalized-IVW precision, and leverage across a 2,479-gene catalogue and three FinnGen endpoints. The manuscript governance tree remains authoritative outside this repository.

## R1 software and the v1.3.1 maintenance update

The R1 software was published on 17 November 2026 as [v1.3.0](https://github.com/doctorzhu0105-cmd/idd-causal-target-analysis/releases/tag/v1.3.0) (commit `6ada8319292387c63e27346d22ffd847708ef71e`; [version DOI](https://doi.org/10.5281/zenodo.24312891), [R1 concept DOI](https://doi.org/10.5281/zenodo.24312890)). The v1.3.1 maintenance update aligns verifier names, command-line entry points and release documentation. It does not change estimators or reported scientific results. For the release date, consult [GitHub Releases](https://github.com/doctorzhu0105-cmd/idd-causal-target-analysis/releases); for the version DOI, use the record for the exact software version in the [R1 archive on Zenodo](https://doi.org/10.5281/zenodo.24312890).

The R1 matched-set comparison has 6,087 historically complete matrix-based comparisons and 1,068 additional scalar comparisons conditionally reconstructed from archived C3 inputs. The C3 state remains COMPUTED_UNVERIFIED. Together they yield 7,155 of 7,284 evaluable pairs in a combined historical and conditional framework, not 7,155 historical matrix outputs. Of the remaining 129 incomplete pairs (matched-set precision UNKNOWN), 54 have empty common sets, 51 have non-empty sets with a missing high-threshold cache, and 24 have non-empty sets with a missing low-threshold cache. The 150 missing membership configurations are a separate population. The 25 random-effects and 6 fixed-effect direction changes refer only to the historical 6,087.

The four programs in [revision_r1/code/verification](revision_r1/code/verification) are independent numerical validators of previously saved inputs, not original analysis generators. R1 code, an aggregate evidence census and a synthetic fixture are included; third-party inputs, restricted matrices, internal reports and the 1,068-row derived table are not redistributed. The latter remains pending a file-level redistribution review. See revision_r1/README.md and docs/R1_DATA_AND_CODE.md for exact scope. No full end-to-end historical reproduction is claimed.

The root [manifest.tsv](manifest.tsv) records current package files using LF-normalized text hashes. See also the [v1.3.1 release manifest](provenance/RELEASE_MANIFEST_v1.3.1.tsv) and [maintenance notes](docs/RELEASE_NOTES_v1.3.1.md). The historical provenance/RELEASE_MANIFEST.tsv refers to the earlier MASTER 1.4.0 layer, and provenance/RELEASE_MANIFEST_v1.2.0.tsv and its JSON companion are immutable v1.2.0 snapshots; none is the current full-tree membership manifest.

## Historical release layers

### v1.2.0

The v1.2.0 release added the catalogue-wide configuration, information, and precision extension:

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
python -B -m unittest discover -s revision_r1/tests -p "test_*.py"
python -B tools/verify_package.py --root .
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
