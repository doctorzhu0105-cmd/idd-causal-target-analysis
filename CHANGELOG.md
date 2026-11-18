# Changelog

## v1.3.1 — maintenance update

- Use `verify-ab`, `verify-psd`, `verify-r37` and `verify-rho` for saved-input verification, with matching script/config names and aligned output labels.
- Synchronize entry points, documentation, tests and manifests. Estimators, thresholds, scientific results and data-access limits are unchanged.
- For citation, use the date on [GitHub Releases](https://github.com/doctorzhu0105-cmd/idd-causal-target-analysis/releases) and the version DOI on the matching record in the [R1 archive on Zenodo](https://doi.org/10.5281/zenodo.24312890).

## v1.3.0 — R1 software, published 17 November 2026

- Preserve the v1.2.0 code, data derivatives and history; add bounded R1 code and a synthetic fixture under revision_r1.
- Include four independent numerical verification components for saved A/B, fixed-W rho, K1-to-multivariate FE/RE and PSD checks; they are verifiers, not original analysis generators.
- State the distinction between 6,087 historical complete comparisons and 1,068 conditional singleton reconstructions, with C3 COMPUTED_UNVERIFIED. The 1,068-row derivative is withheld pending provider and file-level rights review.
- Record the 129 incomplete comparisons and the separate 150 missing membership configurations.
- Archived at [version DOI 10.5281/zenodo.24312891](https://doi.org/10.5281/zenodo.24312891); the associated R1 concept DOI is [10.5281/zenodo.24312890](https://doi.org/10.5281/zenodo.24312890). The release does not contain third-party source inputs or establish full historical reproduction.

## v1.2.0 — historical catalogue extension

The historical catalogue-wide extension, outputs, manifest and release are retained at the existing tag and commit. See docs/RELEASE_NOTES_v1.2.0.md.
