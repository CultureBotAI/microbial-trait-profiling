# Historical sources for the cellulose plate-number tests

These are verbatim copies of committed repository files, not regenerated analyses
or numeric snapshots taken from the plates. Tests read selected source values and
compare selected displayed values with them. They do not need Git history in CI.

| Fixture | Source commit | Original repository path | Git blob ID |
|---|---|---|---|
| `study_modal_7970cba.json` | `7970cba2ddd40a636470ce767e22ed0b805bee33` | `data/analysis/cellulose_decr/cellulose_decr_summary.json` | `ce39f2e47808aededbda73fdb5e68649bccce301` |
| `integrated_evidence_a3c32cf.md` | `a3c32cf9e32e46ff55b2077ed596b498e2c3e1cb` | `data/derived/integrated_evidence/integrated_evidence.md` | `e74b90f6f9eed341b44f3364123a5c43f6889b50` |
| `methods_integrated_evidence_a3c32cf.md` | `a3c32cf9e32e46ff55b2077ed596b498e2c3e1cb` | `methods/integrated_evidence.md` | `c72893d0d9c541c5579af188b556fac9d57e8037` |
| `table1_log_ratio_validation_1e0d0eef.md` | `1e0d0eef52095d3032b14e0291e52a8417ec21d8` | `methods/table1_log_ratio_validation.md` | `b6557f185d6ea628ce58bfeb88454a256354e05a` |

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


## Plate V document version (#1599)

Plate V names `methods/table1_log_ratio_validation.md` but does not name a
revision. The snapshot here is the document present at the first tracked Plate V
producer, commit `1e0d0eef52095d3032b14e0291e52a8417ec21d8` (#1125), rather than a
claim that the plate quoted that commit. Its blob is identical to the document at
its last prior edit, `f21a44467be67388bcde5831f26bf23197b3877b`. The analysis source
revision named inside the document is a separate provenance field. Later edits
to the live methods note must not update this historical fixture or the plate.

The three earlier fixtures were extracted by the repository script above. The
new document was read with `git show <source commit>:<original repository path>`
using this table, and compared byte-for-byte with that Git object. All four
fixtures are checked in CI without Git history: SHA-1 of `blob <byte size>` plus
a NUL byte plus the contents must equal the documented Git blob ID. The manifest
also covers every preserved source file in this folder. These are content
identity checks, not cryptographic authentication of the historical claims.

## Exact scope of the source pins (#1596)

These tests bind selected displayed values and labels to identified historical
sources. They do not pin every numeric token or every repeated occurrence,
validate confidence-interval geometry, rerun statistical analyses, or establish
scientific validity. Fixture byte identity protects the source, not all claims
made by a plate. No broader numeric-token coverage is claimed by this change.

| Test in `test_cellulose_tufte_statistics_sources.py` | Values pinned |
|---|---|
| `test_plate_five_shared_corpus_and_cascade_match_the_pinned_methods` | Plate V subtitle and cascade profile, experiment, study and retained-trait counts. |
| `test_plate_five_fraction_labels_match_the_tracked_validation_rows` | Six selected fraction summaries: proportions, displayed sample counts, q-values and stars; CEL's summary has no sample count. |
| `test_plate_five_ratio_labels_match_the_tracked_validation_rows` | Four selected ratio summaries (FA/CEL, FA/COC, CEL/SPH, CEL/ARS), plus three displayed effects; no confidence-interval geometry. |
| `test_plate_five_retained_trait_family_matches_the_matrix_source` | Matrix column counts and the displayed six-trait/fifteen-ratio family sizes. |
| `test_plate_five_ora_numbers_match_the_historical_export_audit` | Selected pooled ORA profile, test and hit totals, using the separately dated tracked export audit. |
| `test_plate_six_cell_counts_match_the_historical_inventory` | Four cell-state counts and percentages. |
| `test_plate_six_total_rows_match_the_historical_methods_snapshot` | Selected overall, method and grid-row totals. |
| `test_plate_six_method_counts_match_the_historical_methods_snapshot` | Fraction-direction and log-ratio-direction row counts and two q-threshold counts. |

`test_cellulose_tufte_overlap_sources.py` covers selected set and direction labels.
`test_cellulose_tufte_historical_labels.py` checks historical/nominal-p wording and
current producer image routes; it does not establish comprehensive numerical
provenance. The historical #1573 review's token counts are not a current coverage
measurement. Its examples outside these pins include repeated union/disease/set
size sentences in fig7 and fig9, narrative duplicates in Plates I/III/IV, Plate V's
sampling-unit panel pair count, and Plate VI's free-form fraction summary. They
must not be described as covered merely because the source fixtures are intact.
Retirement notices or labels likewise do not turn those values into verified
current-run results. Broader source pins remain a separate task.
