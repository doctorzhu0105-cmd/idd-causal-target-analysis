PYTHON ?= python
RSCRIPT ?= Rscript
PLINK ?= tools/idd/plink/plink.exe
SMR ?= tools/idd/smr/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe
RENDER_OUT ?= ../idd-master14-render
MASTER15_OUT ?= ../idd-master15-regenerated
MASTER15_RENDER_OUT ?= ../idd-master15-render

.DEFAULT_GOAL := help
.PHONY: help test package-verify release-manifest master14-assets master14-figure1 master15-derivatives master15-figure1 prepare-druggable gate0 bulk genetics-validate genetics-screen smr scrna drug rank gtex docking handoff freeze-verify wgcna-signed ukb-outcome mr-presso official-coloc coloc-susie tables-figures article-verify

help:
	@echo "v1.3.1 maintenance candidate: run 'make test' and 'make package-verify'. Full reruns require docs/DATA_ACCESS.md."

test:
	$(PYTHON) -B -m unittest discover -s tests -p "test_*.py"
	$(PYTHON) -B -m unittest discover -s revision_r1/tests -p "test_*.py"

package-verify:
	$(PYTHON) -B tools/verify_package.py --root .

release-manifest:
	$(PYTHON) -B tools/build_r1_release_manifest.py
	$(PYTHON) -B tools/build_manifest.py

master14-assets:
	$(PYTHON) py/generate_master14_assets.py --repo-root . --output-dir $(RENDER_OUT)

master14-figure1:
	$(PYTHON) analysis_p0_correlated_mr_20260814/scripts/make_current_figure1.py --repo-root . --output-dir $(RENDER_OUT)

master15-derivatives:
	$(PYTHON) analysis_p0_correlated_mr_20260814/scripts/catalogue_extension.py --results-dir results/master_1_5 --output-dir $(MASTER15_OUT)

master15-figure1:
	$(PYTHON) analysis_p0_correlated_mr_20260814/scripts/make_master15_figure1.py --repo-root . --output-dir $(MASTER15_RENDER_OUT)

prepare-druggable:
	$(PYTHON) py/mr_prepare_druggable_genes.py --repo-root .

gate0:
	$(RSCRIPT) --vanilla R/idd_gate0_precheck.R --repo-root .

bulk: prepare-druggable gate0
	$(RSCRIPT) --vanilla R/idd_bulk_analysis.R --repo-root .

genetics-validate:
	$(PYTHON) py/idd_genetics_pipeline.py validate --plink $(PLINK)

genetics-screen:
	$(PYTHON) py/idd_genetics_screen.py --reuse-prepared-eqtl --reuse-outcomes

smr:
	$(PYTHON) py/idd_run_smr.py --smr $(SMR) --reuse-ma

scrna:
	$(PYTHON) py/idd_scrna_localization.py --repo-root . --targets artifacts/mr/formal_candidate_genes.tsv --reuse-qc

drug:
	$(PYTHON) py/idd_drug_evidence.py --repo-root . --reuse-existing

rank:
	$(PYTHON) py/idd_rank_targets.py --repo-root .

gtex:
	$(PYTHON) py/idd_gtex_coloc.py --repo-root . --targets artifacts/targets/gtex_targets.tsv

docking:
	$(PYTHON) py/idd_docking_assessment.py --repo-root .

handoff:
	$(PYTHON) py/idd_make_handoff.py --repo-root .

freeze-verify:
	$(PYTHON) py/idd_freeze_core.py --repo-root . --verify-only

wgcna-signed:
	$(RSCRIPT) --vanilla R/idd_manuscript_wgcna.R --repo-root .

ukb-outcome:
	$(PYTHON) py/idd_ukb_replication.py --repo-root .

mr-presso:
	$(RSCRIPT) --vanilla R/idd_mr_presso.R --repo-root . --nb-distribution 10000

official-coloc:
	$(RSCRIPT) --vanilla R/idd_official_coloc.R --repo-root .

coloc-susie:
	$(RSCRIPT) --vanilla R/idd_coloc_susie.R --repo-root .

tables-figures:
	$(PYTHON) py/idd_manuscript_analysis.py --repo-root .

article-verify:
	$(PYTHON) py/idd_manuscript_verify.py --repo-root .
