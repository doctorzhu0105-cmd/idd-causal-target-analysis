# External data access

The repository deliberately excludes source and participant-level data. Users must obtain each input from its provider and comply with its terms.

## Current correlated cis-MR inputs

Stage these inputs beneath `data/raw/` before a full rerun:

- eQTLGen full cis-eQTL summary statistics and allele-frequency file;
- FinnGen R13 summary statistics for `M13_INTERVERTEB`, `M13_LOWBACKPAIN` and `M13_SCIATICA`;
- Finan et al. druggable-genome Table S1;
- 1000 Genomes Phase 3 EUR PLINK `.bed/.bim/.fam` files, preserving signed allele orientation;
- registered raw-region coloc inputs and SMR/HEIDI source files;
- any UKB-SAIGE summary data required by the historical five-gene secondary analysis.

The exact source registry and configuration conventions are in `config/`. File hashes in `provenance/` identify redistributable derivatives only; they do not grant redistribution rights for upstream inputs.

## Not redistributed

No UKB cache, participant-level genotype/expression data, API response cache, credentials, third-party executable, raw eQTL/FinnGen extract or source file with uncertain redistribution permission is included. The public code accepts local paths; never commit those inputs to this repository.
