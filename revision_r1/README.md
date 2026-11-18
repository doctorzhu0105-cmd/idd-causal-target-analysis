# R1 code and evidence scope

This directory contains bounded R1 code components and an aggregate matched-set census. Existing v1.2.0 repository contents remain in their original locations. These are not a complete historical raw-data-to-result pipeline.

## Running bounded components

The dispatcher requires a named analysis, an authorized input directory, a JSON config and a new output directory outside both inputs and this repository. Inputs and outputs cannot overlap. All paths in the config are relative to the input directory.

    python -B revision_r1/code/run_analysis.py --help

A fully artificial fixture can be used for a single aligned configuration:

    python -B revision_r1/code/run_analysis.py --analysis givw --input-dir revision_r1/fixtures --config revision_r1/config/givw_fixture.json --output-dir ../r1-fixture-output

The fixture is synthetic, and this command does not validate the catalogue or regenerate historical MR/GWAS data.

The dispatcher exposes single-configuration GIVW and attenuation, saved-input K=1 multiplicity propagation, K-transition counting, selected R3-9 summary calculations, and figure adapters. Saved-input comparisons previously matched historical outputs where documented; required upstream tables and real figure inputs are not included. The figure adapters were syntax-checked and were not rendered in this build. The selected external outcome comparison shares exposure and LD inputs. It therefore does not provide an end-to-end independent replication. The earlier independent audit flagged 22 historical supplementary claims: 15 have since been reproduced and one checked against metadata; six remain traceable only to historical summaries. Those six do not supply the final S17 co-primary outcomes or their confidence intervals. The source-adapted full-universe entry has not been rerun end-to-end; it requires external authorized inputs and software. Its presence does not establish that the historical pipeline can be fully regenerated.

## Four independent numerical verification programs

- verify-ab checks distributions and directions in supplied saved A/B results; it does not create those historical scenarios.
- verify-rho recounts fixed-W rho/PSD-envelope criteria in a supplied saved K1-to-multivariate paired table.
- verify-r37 checks supplied K1-to-multivariate FE/RE, phi and regime values; project convention phi=1 for K=1 is not evidence that heterogeneity is absent.
- verify-psd checks at most 40 eligible supplied matrices against their SHA-256 identities and evaluates selected joint PSD and Schur conditions. The matrix files are not included; such checks do not establish empirical overlap or catalogue-wide verification.

Run from the repository root, with `inputs` pointing to authorized inputs and a fresh external output directory for each call:

```bash
python -B revision_r1/code/run_analysis.py --analysis verify-ab  --config revision_r1/config/verify_ab.json  --input-dir "$inputs" --output-dir ../verify-ab-output
python -B revision_r1/code/run_analysis.py --analysis verify-psd --config revision_r1/config/verify_psd.json --input-dir "$inputs" --output-dir ../verify-psd-output
python -B revision_r1/code/run_analysis.py --analysis verify-r37 --config revision_r1/config/verify_r37.json --input-dir "$inputs" --output-dir ../verify-r37-output
python -B revision_r1/code/run_analysis.py --analysis verify-rho --config revision_r1/config/verify_rho.json --input-dir "$inputs" --output-dir ../verify-rho-output
```

The supplied configs use the input basenames listed in each JSON. For saved inputs held under different filenames, make a local copy of the relevant config and enter the relative input name; the file contents and table schema must match. In particular, the supplied R37 config expects `R37_RECOMPUTED_2529.tsv`, the rho config expects `R2M1_RHO_GRID_PAIRED_7284.tsv.gz`, and the external outcome config expects `R39_GENE_LEVEL_REPRODUCIBILITY.tsv`. The AB and rho reports are written to `ab_verification_results.json` and `rho_verification_results.json`, respectively.

All four are independent saved-input validators. Historical selected outputs matched in an earlier local verification; the source inputs and derived verification files are not distributed. A local syntax/import or help check cannot replace execution against authorized inputs.

## R1 matched-set evidence

The 7,284 matched gene-endpoint pairs are not by themselves proof of identical SNP membership. The relevant census is results/K1_CONDITIONAL_EVIDENCE_COUNTS.tsv. It distinguishes 6,087 historical complete matrix comparisons from 1,068 current singleton scalar reconstructions, labelled CONDITIONALLY_RECONSTRUCTED_FROM_ARCHIVED_C3_INPUTS with upstream C3 state COMPUTED_UNVERIFIED. The combined 7,155 is conditional. The 129 incomplete pairs have UNKNOWN matched-set precision and are distinct from the 150 missing membership configurations; six K1 pairs remain UNKNOWN. The 25 RE and 6 FE direction changes refer only to the 6,087 historical complete comparisons.

The 1,068 keyed derived numeric rows, source-level data, allele and SNP records, and LD matrices are not included while redistribution conditions remain unverified. The local forensic generator and internal review reports are not part of this repository. See docs/R1_DATA_AND_CODE.md.
