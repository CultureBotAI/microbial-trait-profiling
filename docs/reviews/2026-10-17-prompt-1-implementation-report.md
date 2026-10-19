# Prompt 1 implementation and evidence handoff

Date: 2026-10-17 (America/Los_Angeles). Tracking issue: [#551](https://github.com/CultureBotAI/microbial-trait-profiling/issues/551).

This report describes the original implementation and generation. The
[2026-10-19 follow-up](2026-10-19-analysis-limitations-followup.md) records
current-core compatibility, additional review fixes and refreshed execution
evidence. The retained calibration artifacts now record that later execution;
its replicate results are byte-identical to those reported here.

## Implemented scope

The integrated evidence pipeline now imports the native pooled composition tests and complete enrichment/depletion exports, retains the existing BY-corrected 15 log-ratios, and attaches the compatible full-trait directional screen. Native source readers reject incompatible test families, populations, duplicate profile tests, inconsistent backgrounds, and mismatched input hashes where supplied.

Historical composition denominators are explicitly reported rather than changing the original estimand: 2,335 pairs and 959 studies, including 94 zero-total pairs; 2,241 pairs and 939 studies have positive totals on both sides. Zero-total pairs remain in this historical native analysis. The integrated directional analysis uses S1 and excludes them.

ORA calibration now excludes MediaDive before simulation and spike selection and reports FDR and TPR at q = 0.05 and q = 0.001 on the same simulated replicates. Automatic enrichment spike selection requires a positive shift relative to background prevalence. Declared truth is restricted to testable profile/trait pairs, with deterministic spike eligibility independent of random taxon identities. The CLI import warning was corrected at its source.

## Corrections to the supplied prompt

- The canonical coccus and sphere IDs are `METPO:1000668` and `METPO:1000683`.
- Native composition uses two-sided one-sample **t-tests**, with separate study and experiment results. Its support/variance-filtered correction family has **315 traits**, although the retained denominator universe has 317. The integration uses study-level p/q/effect/interval values together; experiment-level p/q remain explicitly secondary.
- The current 15 log-ratio rows already had BY15 q-values. They are preserved, with 3 significant at 0.05 and none at 0.001.
- ORA ranked files omit traits with fewer than two significant hits. Complete result exports are required to establish zero discoveries and actual tested denominators.
- The production depletion rule is q ≤ 0.001 and fold ≤ **0.5**, in its own lower-tail global BH family.
- An ORA profile-hit summary has no aggregate trait-level or cell-level p/q. Imported profile summaries are explicitly `tested_profiles_summarized`; descriptive aggregated support remains `untested`. This is not missing evidence.
- Full-screen corrections must match the S1 historical population. Disease-versus-healthy results and approximate p × 317 calculations were not substituted for the actual full-family q-values.

## Resulting evidence

The full integrated output has 1,912 rows: 40 non-cell method records and 1,872 descriptive oxygen-guild cell records. All 18 previously absent native source rows are supplied. All 15 log-ratio rows are tested and carry BY15 q-values. Native composition has six study-level results significant at the integration's 0.05 threshold; that cutoff is explicitly distinguished from the original native output, which publishes q without imposing a decision threshold.

The full-screen values below are recorded alongside, without replacing, the within-six BH values. BY is the full screen's primary correction.

| Trait | Full-family BH317 q | Full-family BY317 q |
| --- | ---: | ---: |
| Facultatively aerobic | 0.000189459 | 0.00120073 |
| Cellulose catabolic process | 0.0046464 | 0.0294475 |
| Coccus shaped | 0.0366373 | 0.232196 |
| Sphere shaped | 0.0141507 | 0.0896829 |
| Arsenate reduction | 0.0213112 | 0.135064 |
| Moderately halophilic | 0.00157528 | 0.00998364 |

All six meet BH317 q ≤ 0.05; only facultatively aerobic meets BH317 q ≤ 0.001. Under the primary BY317 correction, three meet 0.05 and none meet 0.001. Same-corpus selection still precludes a claim of independent confirmation.

Each Table 1 trait has 4,670 native ORA profile tests in the global 1,480,390-test enrichment family:

| Trait | Enriched increased profiles | Enriched decreased profiles |
| --- | ---: | ---: |
| Facultatively aerobic | 187 | 150 |
| Cellulose catabolic process | 181 | 186 |
| Coccus shaped | 243 | 200 |
| Sphere shaped | 235 | 180 |
| Arsenate reduction | 130 | 74 |
| Moderately halophilic | 222 | 179 |

None of these six traits has a significant lower-tail depletion hit. Best/median q and median/maximum fold summaries are conditional on significant enrichment hits. Enrichment among taxa reported decreased is different from lower-tail ORA depletion.

## Calibration result and interpretation

The production-aligned run completed 10 of 10 replicates, with no failed or invalid replicates. It uses 317 retained/testable biological traits, 4,670 paired profiles, 1,480,390 tests per replicate, and minimum enrichment fold 1.5. Its ORA input contains 2,333 excluded MediaDive columns; this differs from the 2,290 in the raw annotation-count input and is recorded separately.

| Discovery q | Complete-null empirical FDR | Mixed declared-spike FDR | Mixed TPR | Mean true discoveries / true signals |
| --- | ---: | ---: | ---: | ---: |
| 0.05 | 0 | 0.222857 | 0.000428266 | 1 / 2,335 |
| 0.001 | 0 | 0.075000 | 0.000428266 | 1 / 2,335 |

**These mixed-signal results do not establish FDR control at either threshold.** Whole-taxon spiking can also shift correlated traits outside the declared spike/identical-twin truth set. The reported mixed FDR is therefore a declared-spike benchmark with that limitation, not a definitive false-discovery rate over all induced alternatives. Ten replicates also provide limited precision; the replicate-level results and dispersion are retained. The legacy `is_calibrated` boolean is explicitly scoped to the complete-null raw-p check and is not a mixed-FDR verdict.

TPR measures the complete q/fold decision rule against directionally shifted declared signals, including shifts below the 1.5-fold discovery requirement. It is not power conditional on every spiked alternative exceeding that effect cutoff.

These analyses measure annotation representation in literature-reported taxon signatures. They do not measure microbial activity, cell abundance, or metabolite flux.

## Reproduction and artifacts

From the issue-linked checkout with the same local input corpus:

```bash
just sync
just trait-profile-pooled-stats
just calibrate
just integrated-evidence
just check
```

The regeneration recipes write:

- [Native pooled statistics](../../data/derived/trait_profile_paired_contrast_trait_stats_pooled.tsv), including test identity, family size, zero-total counts and input checksums.
- [Integrated method results](../../data/derived/integrated_evidence/integrated_evidence_method_results.tsv), [readable evidence table](../../data/derived/integrated_evidence/integrated_evidence.md), and [provenance](../../data/derived/integrated_evidence/integrated_evidence_provenance.json). The primary-family inventory and assessment/feature sidecars are retained beside them.
- [Calibration summary](../../data/derived/ora_calibration/calibration_summary.json), [replicate metrics](../../data/derived/ora_calibration/calibration_replicates.tsv), and [null ECDF](../../data/derived/ora_calibration/calibration_ecdf.tsv).

The complete 1,912-row `integrated_evidence.tsv` and raw calibration null p-values are regenerated locally and remain ignored because they exceed the repository's 500 KB policy for newly committed derived files. The committed method-results table contains all 40 non-cell records; the Markdown view preserves counts of every cell state. Small cited artifacts are opted into version control individually.

The native historical fallback can read the legacy pooled TSV when the regenerated source is absent; legacy sources without input hashes explicitly lack that provenance and should be regenerated for publication. ORA exports do not contain original input hashes: their exact file hashes and profile identity checks are retained, but this does not prove the upstream taxon background was built from a newly changed corpus. Full-screen joins require both input hashes and matching study tallies/contrasts.

The first full verification exposed pre-existing warnings. Minimal repairs replaced deprecated colormap mutation, handled exact constant t-test and zero-probability Fisher limits without changing the ordinary statistical path, and modernized LinkML checks with a 1.11 minimum (also addressing [#500](https://github.com/CultureBotAI/microbial-trait-profiling/issues/500)). The standard `just lock` invocation migrated lockfile metadata to the installed pinned uv format; dependency package versions did not change. These validation repairs did not regenerate parallel figure or coupling outputs.

Work was isolated in `/private/tmp/mtp-pipeline-20261017` on the issue-linked branch `551-integrate-native-evidence-and-ora-calibration`. The parallel cellulose, methods, slide files, and existing main-checkout edits were preserved. No ad hoc helper scripts were needed.

## Validation and delivery

Final full test run: **862 passed, 1 intentionally deselected, no warnings**, in 113.83 seconds. The 36 native-import regressions, 30 calibration tests, and exact-limit numerical tests are included. Separately, 47 coupling/heatmap tests passed with warnings treated as errors.

**`just check` and `just precommit` both passed.** This includes generated-model checks, agent-doc checks, Ruff, mypy, schema validation/lint, Bandit, dependency checks, secret scanning and package vulnerability audits.

The secret scanner initially flagged the `sha256` values for the profile input and taxa matrix at lines 150 and 155 of `data/derived/ora_calibration/calibration_summary.json`. Both were independently recomputed with `shasum -a 256` and match the actual input files.

After explicit user approval, exactly those two reviewed input-checksum entries were added to `.secrets.baseline`. The complete scanner then passed; no inline suppression was added.

Implementation and reviewed artifacts are delivered through the branch linked to [#551](https://github.com/CultureBotAI/microbial-trait-profiling/issues/551), with local verification complete before the first commit. Remote CI and Copilot review are separate publication checks, whose current status is recorded on the PR. The existing Copilot registration problem is tracked in [#267](https://github.com/CultureBotAI/microbial-trait-profiling/issues/267); an absent review must not be treated as approval. This work does not authorize merging the PR.

Generated evidence records the base revision `402ef5c25b7105c9aa73398085956dcf1f372bb4+dirty` and the exact integrated-analysis source checksum. This describes the implementation state at generation time, before the publication commit.
