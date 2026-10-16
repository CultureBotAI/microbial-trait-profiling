# Reviewed integration closeout, 2026-10-16

The initial claim that all listed work was already resolved on `452426c` was
partly correct. The clustering, oral figures, contributor figures, release
recovery instructions, and production denominator wiring had reached `main`.
PR #197 still contained unique historical audit evidence, two stale source
hashes, and work that must remain distinguishable from current results.
Adversarial review also found the calibration and download defects below.

## Corrections in this change

- **#534:** Composition calibration normalizes against the full retained
  non-MediaDive universe before selecting the tested traits. Selecting one
  trait previously made every positive fraction one. Regression tests compare
  subset and full-universe scores and raw inference, including numeric metadata
  and zero-signal profiles. The tested family still controls its own BH correction.
- **#535:** A matching pinned download URL or HTTP 304 does not establish local
  file integrity. Cached BugSigDB dump bytes now pass the pinned size/checksum
  verification before reuse. Corruption raises `IntegrityError`; an explicit
  forced refresh still verifies the replacement. Regression cases exercise
  missing validators, ETag, Last-Modified, HTTP 304, and forced replacement.
- **#536:** The ecological reporter note now distinguishes a positive difference
  score from the two-axis sign criterion. A difference can be positive when both
  component contrasts are positive. The maintained note is corrected; its prose
  is not emitted by a generator, and no figure or statistical result is changed.
- **#202:** Composition and taxon-influence calculations share the scalar paired
  log2-ratio change, including the four-positive-count rule. Taking logs before
  subtraction also avoids overflow in intermediate ratios.
- **#528:** Caller regressions retained from the superseded PR #526 cover trait
  selection, eligibility counts, cached-total rejection and valid-cache reuse.
  An additional case covers the follow-up loader that #526 left untested.
- **#537:** The secret scanner reviews the archive manifest normally. Its exact
  public hashes are audited in the existing baseline; there is no new path
  exclusion or inline suppression. `just secrets-scan` now accepts the scanner's
  arguments so the supported baseline update and audit workflow is available.

## Resolution evidence for the requested issues

This table records code and artifact evidence, not an assertion that all
historical outputs have been regenerated after every later change.

| Issue | Resolution evidence |
| --- | --- |
| #527 | #533 filters MediaDive before both totals and trait counts in the study and condition-module loaders. Eligibility thresholds and profile rankings can therefore change on rerun; they now use the same column universe. New caller regressions check this. |
| #531 | #533 rejects `--include-media-traits` in both composition heatmaps, avoiding media numerators with non-media denominators. |
| #532 | The #526 body correctly states 5,828 of 5,829 metadata rows match the old all-column totals. The remaining stale row is acknowledged; this text correction does not refresh the cache. |
| #202 | #503/#533 separate `composition_denominators.py` from `trait_composition.py`, retain the canonical contribution-input module and MediaDive marker, and this change completes the shared log-ratio arithmetic. |
| #207 | #443 corrected LF Python hashes and named permanent baseline tags. This change preserves the final #197 head and publishes corrected final-head Justfile/lock hashes in the separate archive manifest below. |
| #198 | #509, integrated through #516, supplies clustered overviews, readable pages, tree tables, and comparison diagnostics. |
| #220 | Canonical label sorting precedes linkage; rows with no defined neighbor are excluded into a separate block. Row-permutation and missing-overlap tests cover both rules. |
| #221 | BH borders, explicit missingness, twelve-significant-digit exports, and five-trait comparison eligibility replace the earlier ambiguous diagnostics. Superseded branch-label tables are explicitly historical. |
| #258 | #266/#516 supply the inclusive oral-site OR oral-diagnosis selection, readable figures, audit tables, source provenance, and unchanged source inference. |
| #259 | The secret hook allows the exact generated oral provenance path while retaining scanning elsewhere. |
| #261 | Oral rendering rejects missing/conflicting/noncanonical row keys and manifests without six distinct nonempty trait IDs. |
| #262 | The figure and row audit label support as the maximum per-trait study count; unequal-support regression coverage verifies the calculation. |
| #273 | #440/#516 supply explicit nominal fields and documented compatibility aliases. Overall and oral borders use BH selection; a nominal-positive/BH-negative oral cell remains unbordered. |
| #277 | #278/#516 preserve grey unresolved mass and full-cell directional shares; direction-specific co-carriage exports display only changed pairs. |
| #285 | #516 includes only conditions with a six-trait individual interval excluding zero and sufficient paired studies. Composite scores cannot qualify conditions. Existing selection tests and committed panel outputs agree. |
| #353 | Abandoned #315 is closed without merging its unverified download path. Main retains the pinned URL and integrity checks; #535 closes the separate cache-verification gap found in this review. |
| #91 | README's raw-data release recovery instructions and CONTRIBUTING's pointer explain restoring `data/raw/` and running ingest without downloading again. |

The clustering audit was run against the pre-reconciliation `7b7abf2` source
table using `just audit-six-trait-diagnostics`. All 3,126 contrast records retain
their effects, p/q values, support fields, flags and missingness within the
declared export precision. The oral subset retains 381 records and has zero
BH-selected cells and one nominal-positive cell. This is a preservation check,
not a new calibration or an endorsement of causal interpretation.

## Preserved PR history and provenance correction

| PR | Disposition | Exact preserved head |
| --- | --- | --- |
| #197 | Production implementation superseded by #503/#533; unique audit scripts, methods, tests and results preserved under [`archive/pr-197-20261016`](https://github.com/CultureBotAI/microbial-trait-profiling/tree/archive/pr-197-20261016). | `865b15fcdc2be0bc9f0be976a37dbd0a2ed9f33d` |
| #526 | Production implementation superseded by #533; missing caller regressions retained in this change. Original three commits preserved under [`archive/pr-526-20261016`](https://github.com/CultureBotAI/microbial-trait-profiling/tree/archive/pr-526-20261016). | `cf4679c5a8becc638f88d749386a40fa3599aa24` |
| #248 | Already merged into the clustering stack; its remote branch was already deleted when this review started. Current lockfile retains the fixed AnyIO version. | `e3083b1a4071de7091d2247e5dee42cccbcdf016` |

The #197 archive retains 215 files (60,023,623 bytes) under
`trait-visualization/denominator-harmonization/`, including its original
provenance. Those historical files are preserved unchanged. They are not
current-main results and were not rerun in this closeout. Adoption of the
historical composition results remains a separate task (#205).

The [archive manifest](2026-10-16-composition-archive-manifest.json) measures the
final archived implementation's LF-normalized SHA256 values and records both
old and corrected hashes for `Justfile` and `uv.lock`. Later CI edits made those
two original fields stale; the six Python implementation hashes and two replay
source hashes still match. The new manifest attests to archived source identity,
not to which later environment generated the original numerical results.

Recreate the manifest from the exact snapshots without running an analysis:

```bash
git fetch origin tag archive/pr-197-20261016 tag archive/pr-526-20261016
just audit-retired-composition-prs
```

The helper verifies each tar's git commit ID and inventories the historical
result files. It writes the same manifest on repeat runs. The original
provenance is never silently overwritten. Use its
`archived_implementation_sha256_lf` fields to identify the final archived
implementation; retain the old values when interpreting the historical record.

To inspect or reproduce the historical audit in an isolated checkout, fetch its
explicit baseline dependencies as well:

```bash
git fetch origin tag issue174-baseline tag issue174-ora tag issue174-matrix
git worktree add --detach ../mtp-pr197-archive archive/pr-197-20261016
```

The three dependency tags resolve respectively to `a2b2df98`, `3ae7ac84`, and
`c75744aa`; their full IDs and original `composition-*` recipes are recorded in
the archived Justfile. Historical runs also require the pinned inputs described
in that archive's methods note. This closeout downloads no release or KG data.

## Validation and limits

Validation uses `just check`, `just precommit`, the archive audit recipe, and
independent adversarial review. New tests exercise the failures above rather
than only repeating the implementation. Four existing real-data-dependent
tests require input bundles unavailable in the isolated checkout; they remain
explicit skips, and the integration-marked test is excluded by project policy.
No statistical pipeline or publication bundle is represented as freshly rerun.
