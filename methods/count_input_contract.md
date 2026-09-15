# Count inputs and profile annotation coverage (#184)

The taxon-by-trait and study-profile TSVs use the shared `count-input-v1`
contract implemented in `src/bsdbng/count_inputs.py`. The contract is enforced
before filtering profiles, calculating denominators or backgrounds, selecting
traits, or deciding whether a statistical result is testable.

## Count-table contract

- Taxon matrices start with `taxon_id`, `taxon_name`, and `taxonomic_rank`.
  `taxon_id` is the unique row key. Profile matrices start with `study_id`,
  `experiment_id`, and `inc_dec`; all three together are the unique row key.
- Headers must be nonempty and unique, include at least one trait column, and
  retain the declared leading-column order. Every physical record must have
  exactly the header width. Blank identifiers and duplicate keys are errors.
- Every count must be an ASCII, finite, nonnegative integer. Explicit `0` is an
  observed numeric zero. Blank, negative, fractional, nonnumeric, `NaN`, and
  infinite tokens are errors carrying the file, physical row, and column name.
  Readers never use floating-point conversion, truncate values, overwrite a
  duplicate key, pad a short row, or turn a malformed token into zero.
- `complete` is the default and publication format. The only supported legacy
  alternative is the explicitly requested `legacy_sparse_zero` mode. It treats
  an in-width blank as a sparse zero, still rejects ragged or malformed rows,
  counts every converted blank, and records that mode/count in the generated
  provenance. Legacy mode must not be inferred from a bad cell.
- Both plain and gzip-compressed TSVs pass through the same streaming reader.
  `just validate-profile-counts [path] [mode]` and
  `just validate-taxa-counts [path] [mode]` provide standalone preflight checks.

`read_taxa_trait_matrix()` and the shared `load_profiles()` used by the guild,
Table 1, and composition code paths call this validator directly. Statistical
method branches must consume the same reader or run the same preflight before
their own parsing; a method-specific permissive fallback is not a valid input
mode. The shared loader validates every count, including rows and columns later
removed by host, pairing, MediaDive, or metadata filters, so filtering cannot
hide corrupt input.

## Profile coverage sidecar

`just study-trait-profiles` writes three coordinated artifacts:

1. the numeric profile matrix;
2. `<profile-stem>_coverage.tsv`;
3. `<profile-stem>_provenance.json`.

The sidecar is keyed by `study_id`, `experiment_id`, and `inc_dec`. For each
profile it records:

| Field | Meaning |
| --- | --- |
| `reported_distinct_taxa` | Distinct valid taxon identifiers collected from the source signatures. |
| `lookup_matched_taxa` / `lookup_unmatched_taxa` | Reported taxa present or absent in the taxon-by-trait lookup. |
| `matched_taxa_with_annotations` / `matched_taxa_without_annotations` | Matched taxa with at least one positive annotation, or none, under the declared trait-column filter. |
| `unmatched_taxon_ids` | Sorted JSON list preserving identifiers absent from the lookup. |
| `matched_without_annotation_taxon_ids` | Sorted JSON list preserving matched taxa with no annotation under the filter. |
| `trait_filter_id` | Column count and digest identifying the filter used for the coverage classification. |
| `coverage_status` | `no_reported_taxa`, `no_lookup_matches`, `partial_lookup_coverage`, `matched_without_annotations`, `partial_annotation_coverage`, or `complete_annotation_coverage`. |

The reader verifies that matched plus unmatched equals reported, annotated plus
unannotated equals matched, identifier-list lengths agree with their counts, and
the status follows from those counts. When a statistical consumer is given a
coverage sidecar, its profile keys must match the numeric matrix exactly. Coverage
is attached to the profile as source metadata; it is not converted into a count,
an exclusion, a p/q value, or a display mask.

The provenance JSON stores input/output SHA-256 hashes, count mode, any legacy
blank conversions, trait-filter identity, and sums of every coverage field. Those
sums reconcile the per-profile sidecar with the run summary. A profile whose
numeric vector is all zero can therefore still be distinguished as fully mapped
but unannotated, partly mapped, or entirely absent from the lookup. None of these
states means biological absence or statistical nonsignificance.

## Validation

Fixtures cover short and long rows, duplicate columns and keys, blanks versus
explicit zeros, negative/fractional/nonnumeric tokens, gzip input, legacy sparse
mode, coverage-key reconciliation, and the distinction between a matched but
unannotated taxon and a missing-lookup taxon. Complete valid fixtures retain their
previous numeric profile and statistical results.

On 2026-09-15, `just validate-profile-counts` validated the compressed profile
artifact committed on this branch: 12,108 rows, 2,651 count columns, no blank
counts, and no contract violations. The taxon matrix is generated and gitignored;
validate it with `just validate-taxa-counts` whenever that upstream artifact is
available.
