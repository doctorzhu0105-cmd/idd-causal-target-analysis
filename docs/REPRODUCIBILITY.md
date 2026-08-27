# Reproducibility guide

## 1. Environment

```bash
conda env create -f environment/environment.yml
conda activate idd-causal-target-analysis
python -B -m unittest discover -s tests -p "test_*.py"
python tools/verify_package.py --root .
```

The Python lock records Python 3.12.10, NumPy 2.4.6, pandas 2.3.3, SciPy 1.18.0 and Matplotlib 3.11.0. Restore R packages from `environment/renv.lock`. PLINK 1.9 and SMR 1.3.1 are external executables and are not redistributed.

## 2. Current public assets

`results/master_1_4/` contains small redistributable derivatives. Generate the current Figure 1 and the neutral evidence assets outside the repository tree:

```bash
python analysis_p0_correlated_mr_20260814/scripts/make_current_figure1.py --repo-root . --output-dir ../idd-master14-render
python py/generate_master14_assets.py --repo-root . --output-dir ../idd-master14-render
```

The evidence-asset command verifies that the public Table 3 has 13 rows, contains the registered fields and has no tier/grade encoding. It verifies that its derived forest source is value-identical to the promoted Supplementary Figure S2 source.

## 3. Prepared-input order

After obtaining the external inputs and staging them at the paths in `docs/DATA_ACCESS.md`:

1. Run `full_universe_analysis.py prepare` to produce filtered exposure and endpoint extracts.
2. Run `full_universe_analysis.py mr`, then `fdr`, `adequacy` and `m5` as applicable.
3. Run `compute_leverage_metrics.py` and `conditional_f_statistic.py`.
4. Run the correlated MR and LOO paths with `correlated_mr.py`, `run_all_genes_endpoints.py` and `primary13_loo.py`.
5. Use `py/master14_sensitivity.py` for BH/BY families, effective dimension, leverage, same-panel LD bootstrap primitives and deterministic shrinkage.
6. Run `coloc_prior_sensitivity.py` with the registered raw-region configuration.
7. Compare newly generated summaries with `provenance/OUTPUT_HASHES.tsv` and the included MASTER 1.4 derivatives.

The original Stage 3 governed runner is not published wholesale because it embeds private MASTER paths, handoff text and authority-lock checks. The public module extracts its scientific functions; exact upstream source hashes are retained in `provenance/ANALYSIS_VERSION.json`.

## 4. Limits

A fresh clone supports syntax, synthetic-regression, manifest and included-output parity checks. It does not claim full end-to-end numerical reproduction without the external eQTL, outcome and LD inputs. The bootstrap uses the same reference panel and is not independent replication. Simulation remains deferred.

The older `py/`, `R/` and Makefile targets are retained as the historical v1.0.0 workflow. They are not the authority for current manuscript-facing correlated cis-MR results.
