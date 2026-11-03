# Deferred generated-output deletions from PR #1219

**Decision: retain all 522 candidate files in the #1344 port.** This inventory
does not authorize deletion. Existing code still reads some candidates, and
others hold provenance or explanatory text. A regeneration recipe alone does
not establish that it recreates the same historical evidence.

## Identity and search scope

This is the final pre-gate 2026-11-02 working-tree inventory (grouped scan a03),
captured after the #1344 source and test authors froze their files.
The deleted-path input is the complete result of
`git diff --diff-filter=D --name-only ab9156f...origin/pr-1219-readonly`:

- Source PR #1219 base: `ab9156f4429e91e1c2278e171b9728bbf4edc4da`.
- Source PR #1219 head: `1eda7969abb8a4e22e8a0b7fdb08012236b59c3f`.
- Review checkout starting commit: `cd720eef19b28d41dbb551b31c02d1bbcb91e848`.

The helper's `scan.json`, retained with the execution evidence outside the
repository, binds the deleted-path list, helper, and every scanned source file
by SHA256, with file sizes and a canonical inventory digest. The command below
recreates that record with each inventory. The starting commit is not a claim
that the working tree was clean. Candidate output files were checked for
presence, not independently hashed or scientifically revalidated.

The helper enumerates the filesystem directly, including hidden and ignored
files. It scans first-party text under `docs`, `methods`, `tests`, `scripts`,
`src`, `.github`, `schema`, and `workflows`; documentation under
`trait-visualization` and `data/analysis`; and root-level text documentation
and `Justfile`. It does not consult Git ignore rules. All 738 discovered text
files were read; ten Python bytecode-cache directories were excluded. The helper
declares the text extensions; the manifest records the cache/vendor exclusions
and 10 MB per-text-file limit. The helper and its report artifacts are excluded.

This is a bounded text search, not a search of every possible consumer:
binary decks/images, other output data, external projects, separately stored
runs, Git history and names assembled at runtime are not exhaustively parsed.
No filename match therefore means only no lexical match in this scan.

## Complete inventories

- [Candidate files](issue1344_deferred_generated_outputs.files.tsv): all 522
  paths, directory, format, checkout presence, reference counts and the
  `deferred_no_deletion` decision.
- [Filename references](issue1344_deferred_generated_outputs.references.tsv):
  3,464 path/basename matches grouped into 2,115 source/target/token records.
- [Directory references](issue1344_deferred_generated_outputs.directories.tsv):
  8,105 directory path/basename matches grouped into 1,135 records, useful for
  locating constructed names, output-directory defaults and glob readers.
- [Deleted-path input](issue1344_deferred_generated_outputs.deleted-paths.txt):
  the unchanged Git-diff list used by the helper.

Both reference tables retain every matched line in `line_numbers`; grouping
does not discard occurrences.

All 522 candidates are present. There are 333 TSV, 126 PNG, 25 SVG, 24 Markdown,
11 JSON, and one each of PDF, PowerPoint and HTML. These are deletion
candidates from generated-output directories; the inventory does not establish
that every Markdown or provenance file can be regenerated.

336 candidates have at least one filename match in a source file that the PR
does not delete. 114 have no filename match at all. Neither number is a count
of runtime dependencies: common basenames such as `README.md`, directory
tokens such as `data`, fixture paths, historical citations, and output writers
produce lexical matches too. The TSVs explicitly flag references originating
in another deletion candidate and identify ambiguous tokens.

| Candidate directory | Files |
|---|---:|
| `data/analysis/cellulose_decr` | 19 |
| `trait-visualization` | 18 |
| `trait-visualization/bsdb-taxa-only` | 4 |
| `trait-visualization/trait-profile-investigation` | 22 |
| `…/trait-profile-investigation/PD-cellulose` | 78 |
| `…/trait-profile-investigation/aerobe-enrichment-v2` | 46 |
| `…/trait-profile-investigation/all-bugsigdb-condition-body-site-maps` | 24 |
| `…/trait-profile-investigation/anaerobe-depletion-v2` | 41 |
| `…/trait-profile-investigation/enrichment-infographic-followup` | 36 |
| `…/trait-profile-investigation/enrichment-infographic-followup/my-top3-enriched-traits` | 12 |
| `…/trait-profile-investigation/parkinsons-disease` | 15 |
| `…/trait-profile-investigation/profile-condition-bars` | 12 |
| `…/trait-profile-investigation/six-shared-traits` | 152 |
| `…/trait-profile-investigation/six-shared-traits/archive/2026-09-23-branch-labels` | 5 |
| `…/trait-profile-investigation/six-shared-traits/oral` | 11 |
| `…/trait-profile-investigation/taxa-coverage` | 19 |
| `…/trait-profile-investigation/taxa-oxygen-switches` | 8 |
| **Total** | **522** |

The ellipses abbreviate `trait-visualization`; the complete TSV has no
abbreviated paths.

## Opened references requiring retention or migration

These are representative, manually inspected consumers at the captured
source bytes, not a claim that every lexical match is executable.

| Consumer | Observed use and deletion risk |
|---|---|
| `scripts/build_slide_deck.py:446–449` | Constructs `trait_profile_paired_contrast_trait_stats_{stratum}.tsv` for pooled, stool, stool-excluded and saliva strata. All four saved tables occur in the deletion list. A full-path-only search would miss the construction. |
| `scripts/build_slide_deck.py:452–499` | Reads saved oxygen-switch scores, study metadata and condition-trait metadata. The coverage-summary read is conditional on existence; deletion can change the available evidence rather than always causing an immediate exception. |
| `scripts/build_slide_deck.py:116–133` | Maps deleted saved PNG paths to their producer scripts for deck provenance. This is a provenance reference, not by itself an image read. |
| `scripts/build_slide_deck.py:2657–2670` | Recursively collects PNGs under `trait-visualization` for the appendix, then applies exclusions and duplicate checks. Removing a candidate can silently remove an appendix entry without a literal filename reference. |
| `scripts/generate_tufte_plates.py:379–393,641–696,1137` | Reads saved cellulose summary/interaction tables, paired-statistic tables and coverage tables; one read constructs the stratum filename. These are consumers, separate from the recipes that can write outputs. |
| `src/mtp/cellulose_decr/figures.py:125–155` | Declares and loads nine saved JSON/TSV inputs in `FIGURE_INPUTS`; all nine are deletion candidates. |
| `tests/test_deck_builders.py:676–681` | `_same_classifier` reads the actual checkout's `cellulose_decr_provenance.json` to construct a test fixture. This is an actual saved-file dependency; many other test matches merely create temporary fixtures and must not be counted the same way. |
| `Justfile:288–295` | The oral-six-trait recipe records the source revision for saved direction-contrast and trait-selection tables, then invokes the oral-only consumer. Both saved TSVs are deletion candidates. |
| `docs/reviews/2026-10-23-final-rerun-record.md` | Separately describes tracked downstream outputs and core-verified runs, including provenance and input distinctions. Retiring saved outputs requires preserving that distinction and the evidence needed by their citations. |

The deletion list also includes 24 Markdown documents (methods, interpretation
notes and README files) and 11 JSON records (including manifests and
provenance). They require content-level retention decisions rather than a
blanket image/cache cleanup. Conversely, a commit-pinned GitHub citation
continues to refer to its historical Git object: deleting a current checkout
copy does not itself invalidate that historical URL. Relative local reads and
unpinned links need separate review.

## Requirements for a later deletion proposal

Keep this inventory separate from the #1344 source port. A later proposal
should identify the replacement or retained archive for each affected
consumer, preserve provenance and explanatory documents, and update direct
readers, constructed filenames and glob behavior. It should demonstrate the
relevant deck/figure/test paths against the replacement and compare the
resulting artifacts and evidence, not rely on a successful generation command
alone. The 114 no-filename-match candidates receive no automatic exemption.

No pipeline, deck, figure or scientific calculation was executed for this
inventory. No generated file was removed or rewritten, and `.gitignore` was
not changed by this task.

## Reproduce the inventory

From the review checkout, choose a fresh durable output directory in the
`INVENTORY_OUTPUT` shell variable, then run the committed helper:

```bash
uv run python scripts/adhoc/issue1344_inventory_deferred_deletions.py \
  --repository . \
  --deleted-list docs/reviews/issue1344_deferred_generated_outputs.deleted-paths.txt \
  --source-base ab9156f4429e91e1c2278e171b9728bbf4edc4da \
  --source-head 1eda7969abb8a4e22e8a0b7fdb08012236b59c3f \
  --review-base cd720eef19b28d41dbb551b31c02d1bbcb91e848 \
  --output "$INVENTORY_OUTPUT"
```

The helper refuses an existing output directory. Earlier scans remain in the
execution evidence. Any later source change needs a fresh scan: a03 binds the
observed pre-gate working bytes, not a future commit or merged tree.
