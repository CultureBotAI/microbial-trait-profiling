# Harmonized composition outputs

This branch contains the corrected composition producers and the controlled
denominator comparison for issue #174. It does not import the entire
`clustered-heatmap-generation` branch or its independent downstream figures.

- [Measured before/after report](denominator-harmonization/README.md)
- [Methods, reproduction, and captions](../methods/composition_denominator_harmonization.md)
- [Per-trait, per-stratum comparison](denominator-harmonization/trait_before_after.tsv)
- [Correction families](denominator-harmonization/correction_families.tsv)

Run `just composition-denominator-audit` to regenerate the bundle and the
corrected conventional `study_trait_profile_*` producer caches here. The
statistical tables and infographic are in `trait-profile-investigation/`.
The active derived matrix in `data/derived/` is not overwritten.

The `before/` and `after/` audit subdirectories intentionally preserve both
versions. Only `after/` outputs are harmonized. Other branch-specific
renderers must use matching matrix, metadata, and transformed-cache versions.
