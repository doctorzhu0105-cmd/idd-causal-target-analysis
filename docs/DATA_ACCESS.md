# Data access and redistribution boundaries

No source data are bundled. The exact registered URLs, expected sizes or checksums, genome builds, and local paths are stored in `config/idd_sources.json` and `config/mr_sources.json`.

## Registered inputs

- Bulk GEO series: GSE70362, GSE56081, GSE23130, GSE153761, GSE124272, and GSE150408.
- Human disc single-cell source: GSE230809, including GSE229711 and GSE230808 subseries.
- FinnGen R13 outcomes: M13_INTERVERTEB, M13_LOWBACKPAIN, and M13_SCIATICA.
- eQTLGen: significant cis-eQTL, full cis-eQTL, allele-frequency, and official SMR BESD releases.
- eQTL Catalogue GTEx v8 datasets: QTD000356, QTD000281, and QTD000216.
- LD reference: 1000 Genomes Phase 3 EUR, GRCh37, 503 participants.
- Independent outcome-cohort evaluation: public UKB-SAIGE phecode 722 aggregate summary statistics.
- Druggable-gene universe: Finan 2017 supplementary Table S1.
- Drug-target evidence: public Open Targets, DGIdb, and ChEMBL endpoints used by the scripted queries.

## Not included

- UK Biobank participant-level data.
- Raw UKB-SAIGE API-response caches.
- deCODE or UKB-PPP account-dependent pQTL data.
- GEO matrices, scRNA archives, FinnGen files, eQTLGen files, GTEx summary files, and 1000G genotype files.
- PLINK, SMR, or other third-party executables.

Users must obtain each input from its registered provider and comply with that provider's terms. A public URL does not automatically grant redistribution rights for copied source files or large derived extracts.
