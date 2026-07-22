PYTHON ?= python
RSCRIPT ?= Rscript
PLINK ?= tools/idd/plink/plink.exe
SMR ?= tools/idd/smr/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe

.DEFAULT_GOAL := help
.PHONY: help test package-verify prepare-druggable gate0 bulk genetics-validate genetics-screen smr scrna drug rank gtex docking handoff freeze-verify wgcna-signed ukb-outcome mr-presso official-coloc coloc-susie tables-figures article-verify

help:
	@echo "Code-only release: run 'make test' first. Full analyses require the registered external inputs documented in docs/DATA_ACCESS.md."

test:
	$(PYTHON) -B -m unittest discover -s tests -p "test_*.py"

package-verify:
	$(PYTHON) tools/verify_package.py --root .

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
