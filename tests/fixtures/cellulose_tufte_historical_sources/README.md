# Historical sources for the cellulose plate-number tests

These are verbatim copies of committed repository files, not regenerated analyses
or numeric snapshots taken from the plates. Tests read the source values and
compare the plates' displayed values with them. They do not need Git history in CI.

| Fixture | Source commit | Original repository path |
|---|---|---|
| `study_modal_7970cba.json` | `7970cba2ddd40a636470ce767e22ed0b805bee33` | `data/analysis/cellulose_decr/cellulose_decr_summary.json` |
| `integrated_evidence_a3c32cf.md` | `a3c32cf9e32e46ff55b2077ed596b498e2c3e1cb` | `data/derived/integrated_evidence/integrated_evidence.md` |
| `methods_integrated_evidence_a3c32cf.md` | `a3c32cf9e32e46ff55b2077ed596b498e2c3e1cb` | `methods/integrated_evidence.md` |

The JSON tests select `label_source_comparison.study_modal`: the reconstructed
study-modal labels that the dated README compares with experiment labels. The
other contents are preserved for byte identity; they are not permission to mix
the two populations. The integrated-evidence snapshots describe the historical
grid, not the current grid or evidence of tested physiological differences.

The extraction command is recorded in
`scripts/adhoc/extract_plate_historical_sources_1551.sh`. It refuses to overwrite
fixtures, caps each source at 256 KiB, and checks each output against the original
Git blob ID. Missing or conflicting sources are reported separately rather than
replacing them with the plate's own literals.
