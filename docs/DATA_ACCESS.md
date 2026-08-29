# External data access

The repository excludes source and participant-level data. Users must obtain each input from its provider and comply with the applicable terms.

## Required provider inputs

- eQTLGen Phase I cis-eQTL summary statistics and allele-frequency material;
- FinnGen R13 summary statistics for `M13_INTERVERTEB`, `M13_LOWBACKPAIN`, and `M13_SCIATICA`;
- Finan et al. druggable-genome Table S1;
- 1000 Genomes Phase 3 EUR PLINK `.bed/.bim/.fam` files, preserving signed allele orientation;
- registered raw-region colocalization and SMR/HEIDI inputs when those secondary analyses are rerun;
- UKB-SAIGE summary data only for the historical secondary analysis and only under its access terms.

Configuration and provider URLs are recorded under `config/`. The new catalogue entry point requires explicit local paths and never assumes a private project directory.

## Redistribution boundary

The v1.2.0 gene-level output tables are analytical derivatives. They do not contain SNP identifiers, genomic coordinates, alleles, EAF, source-level beta/SE rows, participant identifiers, or participant records. Their inclusion does not grant rights to the upstream inputs.

Not redistributed:

- eQTLGen or FinnGen source extracts;
- 1000 Genomes PLINK files;
- UKB caches or controlled-access content;
- GEO matrices or API caches;
- participant-level data;
- credentials, cookies, or third-party executables.

Users must regenerate restricted-input-dependent components locally after authorized access.
