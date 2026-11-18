# Public reproducibility scope

The historical v1.2.0 repository contains the executable code and redistributable provenance layer associated with the MASTER 1.5.0 catalogue. It is not the manuscript authority and does not publish internal governance records.

## Reproducible from a fresh clone

- package syntax, imports, privacy checks, and deterministic manifests;
- correlated generalized-IVW and related synthetic regression tests retained from v1.1.0;
- row, hash, formula, and cross-version identity checks for the v1.2.0 catalogue;
- deterministic exposure/precision derivatives and matched associations;
- K=1 second-order sensitivity membership counts and K>=2 rank-propagation losses;
- Figure 1 source identities, zero-gain preservation, and rendering.

## Requires authorized external inputs

Full eQTLGen/FinnGen catalogue recomputation, signed 1000 Genomes EUR LD construction, SMR/HEIDI, raw-region colocalization, and UKB secondary analysis require provider-controlled inputs described in `docs/DATA_ACCESS.md`.

## Public analytical derivatives

The v1.2.0 tables contain gene-level model outputs and diagnostics. They contain no SNP identifiers, alleles, coordinates, EAF, source-level association rows, participant identifiers, or participant records. Raw third-party data and redistribution-uncertain caches remain excluded.

The information-precision analyses are descriptive, the K=1 calculation is a sensitivity analysis, and the null diagnostics are not known-truth simulation. Previous v1.0.0 and v1.1.0 reproducibility material remains intact.


## Scope of the published R1 software and maintenance candidate

The revision_r1 code and aggregate census are a bounded R1 extension. The four numerical verification programs verify previously saved inputs; they do not generate the original catalogue. The 1,068 C3-derived singleton values are conditional, their underlying C3 state remains COMPUTED_UNVERIFIED, and the keyed data and source inputs are not redistributed in this candidate. This addition does not establish end-to-end historical reproduction.
