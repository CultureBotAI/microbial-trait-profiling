# ORA provenance refresh

Date: 2026-10-19. Issue: [#720](https://github.com/CultureBotAI/microbial-trait-profiling/issues/720).
Implementation: [PR #722](https://github.com/CultureBotAI/microbial-trait-profiling/pull/722).

The user requested regeneration of the historical ORA exports after #705 added
required input fingerprints. The execution uses committed analysis revision
`b284b29af92b7e28c0956e66acf195e2272a94c2` and existing verified core
`run-49a0c712ad3d4fc59d10804951fe511a`. Historical exports remain unchanged.

The output root is:

```text
/Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/microbial-trait-profiling/data/derived/ora_provenance_refresh/issue720
```

## Input identity and execution

The core manifest is
`/private/tmp/mtp-workflow-629/.tmp/core-production-final/run-49a0c712ad3d4fc59d10804951fe511a/run.json`.
The measured input identities are retained in `input-sha256.txt` under the output
root. The committed [audit manifest](2026-10-19-ora-provenance-refresh-inputs.json)
identifies every historical/current result pair and the actual input files.
Stratified fingerprints refer to generated partition files; their membership is
determined using the core experiment table, whose hash is also retained.
The refresh reuses the core and native inputs at the paths listed here. Full
future revalidation requires retaining those inputs; the new output directory
is not a self-contained archival copy of the core and other analysis arms.

The production commands use these bindings (`OUTPUT` must be a fresh directory):

```bash
CORE=/private/tmp/mtp-workflow-629/.tmp/core-production-final/run-49a0c712ad3d4fc59d10804951fe511a
NATIVE=/private/tmp/mtp-issue551/.tmp/issue551-native-core
NOEL=/private/tmp/mtp-workflow-629/.tmp/issue658-noel-production
OUTPUT=/absolute/path/to/a/new/ora-refresh

just core-ora "$CORE/run.json" "$OUTPUT/enrichment" upper
just core-ora "$CORE/run.json" "$OUTPUT/depletion" lower
just analyze-strata-inputs \
  --profiles "$CORE/steps/ora_profiles/derived/study_trait_profiles_ora_wide.tsv" \
  --experiments "$CORE/steps/export_tables/derived/experiments.tsv" \
  --taxa-matrix "$CORE/steps/taxa_matrix/derived/taxa_trait_matrix.tsv" \
  --output-dir "$OUTPUT/strata" \
  --p-threshold 0.001 --min-fold 1.5 --fdr-scope global

just integrated-evidence \
  --core-run "$CORE/run.json" \
  --arm1-pooled-stats "$NATIVE/trait_profile_paired_contrast_trait_stats_pooled.tsv" \
  --arm2-enrichment-results "$OUTPUT/enrichment/all_results.tsv" \
  --arm2-depletion-results "$OUTPUT/depletion/all_results.tsv" \
  --full-screen "$NOEL/trait_direction_screen/trait_direction_screen.tsv" \
  --full-screen-provenance "$NOEL/trait_direction_screen/trait_screen_provenance.json" \
  --influence-status "$NOEL/table1_taxon_influence/native_manifest.json" \
  --calibration-status "$NOEL/calibration/calibration_status.json" \
  --orientation-provenance "$NOEL/contrast_orientation/contrast_orientation_provenance.json" \
  --output-dir "$OUTPUT/integrated_with_assessments"

just core-verify-consumer "$OUTPUT/enrichment/report.core_consumption.json"
just core-verify-consumer "$OUTPUT/depletion/report.core_consumption.json"
just core-verify-consumer "$OUTPUT/integrated_with_assessments/core_consumption.json"

# The committed manifest contains this execution's absolute paths.
# Update its current/partition paths if reproducing into another fresh directory.
just audit-ora-refresh \
  docs/reviews/2026-10-19-ora-provenance-refresh-inputs.json \
  "$OUTPUT/audit.json"
```

Global BH uses q ≤ 0.001, enrichment fold ≥ 1.5 and depletion fold ≤ 0.5.
The audit compares every table-body byte, all scientific metadata, actual input
fingerprints, row counts and selected-hit counts. It does not fabricate input
fingerprints for the historical files. Native composition, directional screen
and influence evidence retain their earlier source identities when imported.
The calibration and orientation assessments are also imported from the matching
earlier execution, with their input fingerprints checked. The calibration used
one replicate per scenario and retains its FAIL/insufficient-budget findings;
it does not validate native ORA or composition. The orientation review retains
58/60 agreement and two unresolved disagreements. This refresh does not rerun
those assessments. The earlier `integrated/` directory is a diagnostic import
without those two assessment inputs; use `integrated_with_assessments/` for the
integrated result that preserves all supplied assessments.

## Verification results

The [five-export audit](2026-10-19-ora-provenance-refresh-audit.json) passed.
Every result-table body is byte-identical to its historical counterpart; all
scientific metadata is unchanged. The only metadata additions are the profile
and taxa-matrix fingerprints, both checked against actual input bytes. Declared
test counts and selected-hit counts match the observed complete rows.

| Export | Tested profiles | Complete tests | Discoveries |
| --- | ---: | ---: | ---: |
| Pooled enrichment | 4,670 | 1,480,390 | 66,240 |
| Pooled depletion | 4,670 | 1,480,390 | 23 |
| Stratified pooled enrichment | 4,670 | 1,480,390 | 66,240 |
| Feces enrichment | 2,950 | 935,150 | 36,562 |
| Feces-excluded enrichment | 1,720 | 545,240 | 29,498 |

The three complete-trait summaries, annotation-density deciles and stratification
sensitivity TSV are also byte-identical to the corrected #664 baseline under
`/private/tmp/mtp-issue664/.tmp/issue664-strata/` (checked with `cmp`).

The completed pooled producers wrote these new core-bound consumer records:

- Enrichment: `consumer:6e55f7ba-46c0-49ab-8417-7e587a754737`, at
  `enrichment/report.core_consumption.json` under the output root.
- Depletion: `consumer:43d0542f-43e0-424c-ba85-f8cf03a09d58`, at
  `depletion/report.core_consumption.json` under the output root.

Both bind core manifest SHA-256
`7fb72fa862a63c8297189aef0f9e4c7ea81403f56cabf26a3d777c6c5faed572`.
The verification commands above independently recheck the full recorded chain.
PR #722 records their final status and the integrated consumer record after
execution; a completed export alone is not a claim that integration or later
verification has already succeeded.

The implementation's full `just check` passed 1,456 tests, with four existing
missing-data skips, one intentional deselection and one upstream warning. All
precommit hooks and CI passed on execution revision `b284b29`. Ten focused
actual-producer regressions cover successful refresh, input mismatch, truncated
results, changed metadata, missing fingerprints and output/input aliases.
Independent review also checked the five production file hashes and headers,
the core/experiment fingerprints, and the additional stratification comparisons.
The PR carries final-head validation status for this report and audit snapshot.
The secret scanner flagged 17 distinct public checksum values in the audit
snapshot. Each was reviewed as a file, table-body or manifest SHA-256 and added
to the existing `.secrets.baseline` as a false positive; scanner rules and
coverage were not changed.

These outputs describe the representation of organism-level trait annotations
among taxa reported increased or decreased in BugSigDB. Regeneration does not
establish in vivo physiological activity, causality, calibration or a
procedure-level influence verdict.
