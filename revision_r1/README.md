# R1 code and evidence scope

This directory contains bounded R1 code components and an aggregate matched-set census. Existing v1.2.0 repository contents remain in their original locations. These are not a complete historical raw-data-to-result pipeline.

## Running bounded components

The dispatcher requires a named analysis, an authorized input directory, a JSON config and a new output directory outside both inputs and this repository. Inputs and outputs cannot overlap. All paths in the config are relative to the input directory.

    python revision_r1/code/run_analysis.py --help

A fully artificial fixture can be used for a single aligned configuration:

    python revision_r1/code/run_analysis.py --analysis givw --input-dir revision_r1/fixtures --config revision_r1/config/givw_fixture.json --output-dir ../r1-fixture-output

The fixture is synthetic, and this command does not validate the catalogue or regenerate historical MR/GWAS data.

The dispatcher exposes single-configuration GIVW and attenuation, saved-input K=1 multiplicity propagation, K-transition counting, selected R3-9 summary calculations, and figure adapters. Saved-input comparisons previously matched historical outputs where documented; required upstream tables and real figure inputs are not included. The figure adapters were syntax-checked and were not rendered in this build. The selected R3-9 comparison shares exposure and LD inputs, and 22 supplementary cells were not verified in the earlier audit; it is not an end-to-end independent replication. The source-adapted full-universe entry has not been rerun end-to-end; it requires external authorized inputs and software. Its presence does not establish that the historical pipeline can be fully regenerated.

## Four independent numerical verification programs

- qwen-ab checks distributions and directions in supplied saved A/B results; it does not create those historical scenarios.
- qwen-rho recounts fixed-W rho/PSD-envelope criteria in a supplied saved K1-to-multivariate paired table.
- qwen-r37 checks supplied K1-to-multivariate FE/RE, phi and regime values; project convention phi=1 for K=1 is not evidence that heterogeneity is absent.
- qwen-psd checks at most 40 eligible supplied matrices against their SHA-256 identities and evaluates selected joint PSD and Schur conditions. The matrix files are not included; such checks do not establish empirical overlap or catalogue-wide verification.

All four are independent saved-input validators. Historical selected outputs matched in an earlier local verification; the source inputs and derived verification files are not distributed. A local syntax/import or help check cannot replace execution against authorized inputs.

## R1 matched-set evidence

The 7,284 matched gene-endpoint pairs are not by themselves proof of identical SNP membership. The relevant census is results/K1_CONDITIONAL_EVIDENCE_COUNTS.tsv. It distinguishes 6,087 historical complete matrix comparisons from 1,068 current singleton scalar reconstructions, labelled CONDITIONALLY_RECONSTRUCTED_FROM_ARCHIVED_C3_INPUTS with upstream C3 state COMPUTED_UNVERIFIED. The combined 7,155 is conditional. The 129 incomplete pairs have UNKNOWN matched-set precision and are distinct from the 150 missing membership configurations; six K1 pairs remain UNKNOWN. The 25 RE and 6 FE direction changes refer only to the 6,087 historical complete comparisons.

The 1,068 keyed derived numeric rows, source-level data, allele and SNP records, and LD matrices are not included while redistribution conditions remain unverified. The local forensic generator and internal review reports are not part of this repository candidate. See docs/R1_DATA_AND_CODE.md.
