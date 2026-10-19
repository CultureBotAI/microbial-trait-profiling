# Analysis limitation fixes and execution evidence

Date: 2026-10-19. Issues:
[#551](https://github.com/CultureBotAI/microbial-trait-profiling/issues/551),
[#664](https://github.com/CultureBotAI/microbial-trait-profiling/issues/664), and
[#671](https://github.com/CultureBotAI/microbial-trait-profiling/issues/671).

## Native evidence and calibration (#551)

The native composition producer and integrated evidence command now accept an
executed core manifest. They resolve the core inputs, reject conflicting paths,
protect source files from output collisions, and write validated LinkML consumer
records. Native composition derives its human paired cohort from the core
experiment table; it no longer needs a pre-existing figure metadata table in
this mode. Historical zero-total pairs remain included and explicitly counted.
The native study-level t-test and its support/variance-filtered BH family are
preserved. Log-ratio BY15 results and the separately corrected full directional
screen retain their own estimands and families.

Core-bound composition imports require matching profile and experiment hashes.
Checking experiment identity matters because changing host assignments can
change the included cohort without changing the all-host count matrix. Complete
ORA exports are checked across their entire correction family, including traits
outside Table 1; malformed or truncated non-selected rows cannot silently pass.
Profile hit counts remain summaries without invented aggregate p/q values.

The adversarial follow-up (#705) additionally requires the complete ORA exports'
profile and taxa-matrix input SHA-256 values to match the verified core products.
Profile identifiers alone cannot establish a shared KG background or annotation
counts. Rerun the updated ORA producer before core integration; archived exports
without fingerprints remain usable only outside core mode, where their input
identity is explicitly unverified. Historical outputs are not assigned new
fingerprints after the fact.
The lower-tail export is auxiliary to a complete upper-tail export; supplying
depletion alone is rejected instead of silently discarding that evidence (#709).

Native composition results with too few studies, all-zero study contrasts, or
insufficient study variance remain untested (#706). The native producer records
these states explicitly; the importer retains descriptive study effects and any
separately testable pair-level result without inventing primary p-values, q-values
or confidence intervals. Their original nominal BH family remains recorded.

Calibration rejects malformed, negative, fractional and infinite counts before
binarizing the taxon matrix. A failed or partial execution cannot report a
successful harness check. Regression fixtures compare the calibrated family
and q/fold decisions with the production analyzer, including MediaDive removal,
zero-background traits, universal traits, unpaired profiles and zero-N profiles.

Use fresh output directories and the same completed core run for these steps:

```bash
just trait-profile-pooled-stats \
  --core-run "$CORE_RUN" --output-dir "$NATIVE_OUTPUT"

just integrated-evidence \
  --core-run "$CORE_RUN" \
  --arm1-pooled-stats "$NATIVE_OUTPUT/trait_profile_paired_contrast_trait_stats_pooled.tsv" \
  --arm2-enrichment-results "$ILYAS_OUTPUT/trait_enrichment_all_results.tsv" \
  --arm2-depletion-results "$ILYAS_OUTPUT/trait_depletion_all_results.tsv" \
  --full-screen "$NOEL_OUTPUT/trait_direction_screen/trait_direction_screen.tsv" \
  --full-screen-provenance "$NOEL_OUTPUT/trait_direction_screen/trait_screen_provenance.json" \
  --output-dir "$INTEGRATED_OUTPUT"

just calibrate-from \
  "$ILYAS_OUTPUT/inputs/study_trait_profiles_ora_wide.tsv" \
  "$ILYAS_OUTPUT/inputs/taxa_trait_matrix.tsv" \
  "$CALIBRATION_OUTPUT" 10 42 all
```

The explicit calibration recipe uses supplied inputs without rebuilding legacy
checkout data. Core integration does not silently import checkout assessment
files; supply assessment paths explicitly when they apply to this run. After
#671 is available, the native influence manifest can be supplied with
`--influence-status "$NOEL_OUTPUT/table1_taxon_influence/native_manifest.json"`.

These commands are an explicit post-arm integration step. Once these fixes and
the owner rerun changes in PR #660 are present together, Ilyas's rerun inherits
the calibration and stratification fixes, and Noel's inherits the native
influence adapter. An individual owner rerun does not import the other owners'
outputs. PR #660 does not run the native pooled-statistics producer or pass the
full-screen paths into its integrated table; use the commands above to obtain
the complete three-arm evidence table.

The current ten-replicate calibration run completed 10/10 replicates with no
failed or invalid replicates, using 317 biological traits, 4,670 paired profiles
and 1,480,390 tests per replicate. Its complete replicate table is byte-identical
to the 2026-10-17 implementation's result. The refreshed summary records the
actual 2026-10-19 execution time and supplied core-derived input paths; input
checksums and numerical results are unchanged. The original report remains a
historical account of that earlier generation.

At q=0.05 and q=0.001, respectively, complete-null empirical FDR remains zero,
mixed declared-spike FDP means are 0.222857 and 0.075, and TPR is 0.000428266 at
both thresholds. These results do not establish mixed-signal FDR control.

The mixed benchmark's declared truth omits some correlated traits shifted by
whole-taxon spiking. Its FDP is therefore a declared-spike benchmark, and the
legacy null-only `is_calibrated` flag is not a general scientific verdict. See
the [original report](2026-10-17-prompt-1-implementation-report.md) for numerical
results, uncertainty and interpretation limits.

The production pooled-statistics output was regenerated from core run
`run-49a0c712ad3d4fc59d10804951fe511a`. Its metadata starts with 7,329 all-host
count profiles, excludes 1,500 nonhuman and 1,159 unpaired profiles, and retains
4,670 profiles in 2,335 human pairs. The native population contains 959 studies,
including 94 pairs with a zero total on at least one side; 2,241 pairs in 939
studies have positive totals on both sides. The retained denominator universe
has 317 traits and the tested native correction family has 315.

Its completed LinkML consumer record is
`consumer:78276952-788f-46c2-b6d1-5faf6351bbb6`, bound to core manifest SHA-256
`7fb72fa862a63c8297189aef0f9e4c7ea81403f56cabf26a3d777c6c5faed572`.
The bounded validation destinations in the issue checkout are
`.tmp/issue551-native-core/` and `.tmp/issue551-integrated-core/`.
They retain the run-specific outputs and consumer records. The committed native
composition and integrated snapshots retain their earlier generation identity;
their provenance is not relabeled as a current-core execution.

## Complete tested denominators (#664)

Stratification now exports a complete trait summary separately from selective
discovery rankings. Tested denominators, cross-stratum support and annotation
density diagnostics use that complete summary. Zero discoveries, zero eligible
tests and unavailable analyses remain distinct. A legacy text-only report with
no known tested denominator reports it as unavailable, never as an infinite
prevalence.

The bounded production regeneration left all three full analyzer-result TSVs
byte-identical. Each stratum retains all 317 eligible traits. In the pooled
example, `METPO:1000487` has 363 significant profiles out of 4,670 tested
(approximately 0.07773), with background taxon count K=394. The prior summary
used 363/363 and replaced K with a profile count sum. The corrected density
summary includes 11 zero-hit and eight one-hit traits omitted from the earlier
298-trait selected ranking. These remain descriptive diagnostics, not a test
of calibration or a biological mechanism.

## Native influence standing (#671)

The adapter reads native execution manifests or compact published review
bundles, validates their corpus and artifact identities, and preserves the
scenario/scope/background inventory and scoped interpretation flags.
Completed execution with no procedure assessment reports
`native influence results reported; no procedure-level verdict supplied #190`.
Only an explicit versioned assessment with a stated procedure-wide basis can
supply a procedure verdict. Missing evidence does not establish that the
analysis never ran.

The bounded production check reused the completed native influence run,
retaining 25 scenarios, four scopes, 200 ORA-family inventory rows and both
existing interpretation flags. The procedure-verdict mapping remains empty.
Completion does not establish robustness, calibrated error rates, or causality.

All three analyses concern trait annotations represented among taxa reported
increased or decreased in the literature. They do not measure microbial
physiology, cell abundance or metabolite flux.

## Validation and delivery

The #551 implementation is delivered in
[PR #553](https://github.com/CultureBotAI/microbial-trait-profiling/pull/553).
Its final `just check` passed 1,317 tests with four existing missing-data skips
and one intentional integration-test deselection. Generation checks, strict
typing, lint, workflow/schema validation, security and dependency audits passed.
The suite reports one upstream seaborn pending-deprecation warning in an
existing heatmap test. The focused native/core tests include an actual executed
core, the production ORA analyzer, native composition output, integrated imports
and consumer-record verification.

Independent review identified and resolved the missing experiment-hash check
and the legacy report-only denominator bug. The stratification final suite
passed 1,238 tests, with four existing missing-data skips and one deselection.
The native influence change in
[PR #697](https://github.com/CultureBotAI/microbial-trait-profiling/pull/697)
passed 1,262 tests with the same skip/deselection counts, plus its full precommit
checks and bounded production integration. Each PR records its current CI and
review standing. Copilot's reviewer-registration problem remains tracked in
[#267](https://github.com/CultureBotAI/microbial-trait-profiling/issues/267);
an absent review does not count as approval.
