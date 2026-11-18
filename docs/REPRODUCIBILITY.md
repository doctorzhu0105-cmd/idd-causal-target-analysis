# Reproducibility guide

## Environment and offline gates

```bash
conda env create -f environment/environment.yml
conda activate idd-causal-target-analysis
python -B -m unittest discover -s tests -p "test_*.py"
python -B -m unittest discover -s revision_r1/tests -p "test_*.py"
python -B tools/verify_package.py --root .
```

The Python environment records Python, NumPy, pandas, SciPy, and Matplotlib versions. Restore R dependencies from `environment/renv.lock`. PLINK and SMR are external executables and are not redistributed.

## Included v1.2.0 results

`results/master_1_5/` contains byte-governed derived outputs for:

- five-threshold Effective K and joint exposure information;
- the row-complete generalized-IVW ladder, including explicit failure-status rows;
- exposure, precision, and leverage paired changes and summaries;
- matched and endpoint-stratified information-precision associations;
- the corrected historical-candidate comparison;
- K=1 second-order uncertainty and multiplicity propagation;
- AKR1C1/AKR1C2 positioning and null-construction diagnostics.

`provenance/OUTPUT_HASHES.tsv` records byte identity to the governed scientific candidate. `provenance/RELEASE_MANIFEST_v1.2.0.tsv` and its JSON equivalent record every release file, role, authority, and redistribution class.

## Regenerate included derivatives

The deterministic postprocessor derives the exposure/precision pair tables, summaries, matched associations, and endpoint associations from the included five-threshold source tables:

```bash
python analysis_p0_correlated_mr_20260814/scripts/catalogue_extension.py \
  --results-dir results/master_1_5 \
  --output-dir ../idd-master15-regenerated
```

Render Figure 1 without filtering zero-gain observations:

```bash
python analysis_p0_correlated_mr_20260814/scripts/make_master15_figure1.py \
  --repo-root . \
  --output-dir ../idd-master15-render
```

Panel B must receive 2,438 rows, including 392 zero nominal-K gains and 377 zero joint-information gains. Panel C must receive 7,284 rows, including 1,125 zero joint-information and precision gains. The symlog axes preserve zero coordinates.

## Full catalogue recomputation

The public `catalogue_compute.py` accepts explicit paths and contains no manuscript, submission, package, handoff, or internal audit dependency. Its required arguments identify:

- the feasibility JSON containing the 2,479 eligible Ensembl IDs;
- prepared eQTLGen batches;
- prepared FinnGen endpoint tables;
- the 1000 Genomes EUR PLINK BIM/BED files;
- an external work directory.

Example interface:

```bash
python analysis_p0_correlated_mr_20260814/scripts/catalogue_compute.py \
  --repo-root . \
  --feasibility-json /authorized/feasibility.json \
  --eqtl-batch-dir /authorized/eqtl_batches \
  --outcome-dir /authorized/finngen_outcomes \
  --bim /authorized/1000G_EUR.bim \
  --bed /authorized/1000G_EUR.bed \
  --work-dir ../idd-master15-work
```

Inputs must first be prepared with the v1.1.0 public scripts and provider-specific conventions. Users are responsible for access and terms. The public release does not claim a full rerun without those external inputs.

## Verification order

1. Run all unit tests.
2. Run `tools/verify_package.py`.
3. Regenerate derivatives into a new external directory.
4. Compare keys and numeric values against `results/master_1_5/`.
5. Render Figure 1 and confirm the source row counts and zero-preservation assertions.
6. For a full external-input rerun, compare the five-threshold GIVW, Effective K, and joint-information tables to the governed output hashes.

The K=1 second-order calculation is `sqrt(SE^2 + beta^2 / min_F)` for K=1 rows only, followed by global BH/BY recomputation over the unchanged 7,284-row family. The observed BH membership changes from 26 to 18 and BY from 5 to 2; four K>=2 BH losses arise through rank propagation.

## Limits

The associations between information changes and precision changes are descriptive. No scalar adequacy score was optimized. The same-panel bootstrap is not independent replication. The included null diagnostics are not known-truth simulation. No raw or participant-level third-party data are bundled.


## R1 components and v1.3.1 maintenance candidate

The components under revision_r1 are bounded code and verification tools for authorized supplied inputs. The R1 singleton census separates 6,087 historical complete comparisons from 1,068 conditionally reconstructed scalar comparisons, whose C3 source state remains COMPUTED_UNVERIFIED. The keyed numeric derivative and underlying third-party inputs are not distributed pending file-level rights review. Four independent numerical programs verify saved inputs; they are not original historical generators. No complete original-data-to-manuscript replay is claimed.


For the current CLI choices, config names and authorized-input interface, see [the R1 guide](../revision_r1/README.md). Release identity and citation details are in [RELEASE_STATUS.json](../RELEASE_STATUS.json).
