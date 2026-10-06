# Plan: three quick fixes (nested test, #295, #254 code)

Base: `enrichment-analysis-heatmap-cluster` @ `5ed546f`. Its CI is green, and it includes #300 and #296.
Each fix gets its own issue, issue-linked branch and PR into `enrichment-analysis-heatmap-cluster`. Fixes are merged one at a time, in order.
Nothing is pushed, merged, closed or commented on without your approval at that step.

```mermaid
flowchart LR
    A[Read issue bodies and comments] --> B[Fix 1: nested test]
    B --> C[Fix 2: 295 calibrate]
    C --> D[Fix 3: 254 complete export code]
    D --> E[Check in with user]
```

## Fix 1: never-collected test (new issue)

Facts:
- `tests/test_enrichment_pipeline.py:542` defines `test_malformed_taxa_count_strict_rejection` inside `test_repeated_taxon_fixture_matrix_wiring`. Pytest never collects it.
- Line 720 duplicates an assert from the outer test.

Steps:
1. Open an issue: the test is nested and never runs. It was added in `85dea78` (refs #131 and #184).
2. Run `gh issue develop` from `5ed546f`.
3. Move the function out to the top level and remove the duplicated assert on line 720.
4. Run it. If a case fails because the analyzer check it covers is broken, fix the analyzer in the same PR and record which check was broken.
5. Run `just check`, then open the PR with `Closes #N`.

Done when pytest reports one more test collected and `just check` passes.

## Fix 2: #295 `just calibrate` exits successfully without running calibration

Facts:
- `ba97ae8` is already on the branch. It added `EXIT_SUCCESS`/`EXIT_FAILURE`/`EXIT_PARTIAL`, input validation, a failure summary JSON, removal of stale outputs, and `test_standard_module_command_executes_benchmark`.
- That commit is not linked to #295, and nobody has posted evidence on the issue.

Steps:
1. Read #295 and its comments. List each acceptance criterion and match it against `ba97ae8` and the current `harness.py`.
2. Close any gaps found. Likely candidates:
   - a test that the `just calibrate` recipe arguments parse;
   - a test that a partial run exits 2;
   - the summary keeping its "not calibrated" wording.
3. Run `just calibrate 2 42 all` on a synthetic fixture, not real data, and capture the exit code, `calibration_summary.json` and the output file list.
4. If gaps were found, open a linked branch and PR. If not, ask whether to post only an evidence comment that references `ba97ae8`.

Done when every criterion in #295 maps to code, a test, or recorded evidence. This does not mean the model is calibrated; that work stays in #263 and #182.

## Fix 3: #254 complete unfiltered export (code part only)

Facts:
- Both `trait_enrichment.py` and `trait_depletion.py` build `all_hits`, which holds every test performed, but they only write `significant_hits`.
- The complete export needs these columns, per the reviewer's #296 comment:
  - raw identifiers;
  - p and q;
  - effects (count, expected, fold);
  - correction family (`fdr_scope` plus a family id);
  - model denominators M, N, K and k;
  - eligibility and provenance: `is_significant`, the thresholds, and the background mode.

Steps:
1. Read #254 and its comments, and lock the column list and file names with you before writing code.
2. Add a `--all-results-output` flag to both analyzers. It writes one row per test, including nonsignificant ones. Its columns are declared in `enrichment_results.py`, and the metadata carries the same `Result-Status` and counts as #296.
3. Wire it into `just analyze` and `just analyze-depletion`, writing `trait_enrichment_all_results.tsv` and `trait_depletion_all_results.tsv`. The existing selected-results file and the exporter stay unchanged.
4. Add tests:
   - the complete-file row count equals `Total-Tests-Performed`;
   - the selected rows are a subset of the complete rows;
   - hand-computed M, N, K and k values match on a small fixture;
   - a zero-discovery run still writes every tested row.
5. Run `just check` and open the PR with `Refs #254`, not Closes.

Done when the code and tests are merged. #254 stays open until real-data output files are produced and checked. That needs a large download, which I will ask about first.

## Out of scope for this pass

#263, #264, #191 (PR #239 is open) and #190 (PR #194 and related PRs are open).
