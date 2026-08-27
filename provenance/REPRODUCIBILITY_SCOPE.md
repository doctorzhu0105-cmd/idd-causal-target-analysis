# Public reproducibility scope

Repository A is the executable public code and redistributable provenance layer for MASTER 1.4.0. It is not the manuscript authority and does not publish internal governance records.

## Reproducible from a fresh clone

- package syntax, imports and deterministic manifest checks;
- synthetic orientation, generalized-Q, pseudoinverse and reduced-Omega LOO regression tests;
- identity/schema checks for included multiplicity, P2RY13, prior-sensitivity, bootstrap and evidence-matrix outputs;
- deterministic rendering of Figure 1 and the neutral Supplementary Figure S2 from included source tables.

## Requires external inputs

Full eQTLGen/FinnGen screening, signed 1000G EUR LD construction, SMR/HEIDI, raw-region colocalization, individual-resample LD bootstrap and re-pruning require the provider-controlled inputs described in `docs/DATA_ACCESS.md`.

## Intentionally not public

MASTER governance files, manuscripts, internal audits/handoffs, author-approval evidence, restricted caches, participant-level data, credentials and redistribution-uncertain extracts are excluded. Replicate-level bootstrap and full multiplicity tables are omitted where included summaries and executable functions provide the required public audit surface.

The bootstrap conditions on the same 503-person 1000G EUR panel and is not independent replication. Simulation is deferred. Historical v1.0.0 analyses remain available but are secondary to the current correlated cis-MR layer.

Some byte-identical Stage 3 tables retain an embedded `MASTER 1.3.0` origin header. Their scientific values were promoted unchanged into MASTER 1.4.0; `ANALYSIS_VERSION.json` and `OUTPUT_HASHES.tsv` record the current promotion identity without rewriting those source bytes.
