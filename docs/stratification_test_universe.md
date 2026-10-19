# Stratification test universe and density diagnostics

The stratification stage runs the ORA analyzer separately for each profile stratum
against the supplied taxa-matrix reference population. It changes the eligible
profiles and the multiple-testing family; it does not rebuild a different
taxa-matrix background for each body site.

Issue #664 corrected selection errors in the ranking denominators and density
diagnostic. The analyzer's `trait_enrichment_all_results_<stratum>.tsv` is now the
source for both selected rankings and the complete trait summary. The exporter
checks the declared complete test count, distinct profile–trait tests, eligible
states, discovery flags, per-trait profile coverage, and background counts.

Outputs have different purposes:

- `trait_enrichment_traits_overall_ranked_<stratum>.tsv` remains a discovery
  ranking requiring at least two selected profiles per trait. Its prevalence is
  selected profiles divided by **all eligible tested profiles**, including
  nonsignificant tests. Its fold summaries describe the selected tests.
- `trait_enrichment_trait_study_pairs_ranked_<stratum>.tsv` remains discovery-only.
- `trait_enrichment_traits_tested_<stratum>.tsv` includes every eligible trait,
  including traits with zero or one discovery. Its minimum p/q values and mean
  folds summarize **all eligible tests**. Traits present in the input profile
  schema but with no eligible tests are recorded separately.
- `annotation_density_deciles.tsv` bins all eligible tested traits by the
  analyzer's reference-population `background_k`, exposed as `trait_total`.
  Repeated counts of annotated taxa across profiles are not substituted for K.
  Each decile reports the fraction of traits with at least one discovery under
  the analyzer's declared p/q **and fold** thresholds. It does not silently use
  a different q threshold.
- `stratification_sensitivity_summary.tsv` compares complete trait summaries.
  Discovery proportions are descriptive, not effect sizes; the reported mean
  fold contrasts now cover all eligible tests, avoiding selection on discoveries.

`tested_nonsignificant` means tests were performed and none passed the analyzer's
selection rule. It is distinct from `zero_eligible_tests`, `unavailable_trait`,
and `unavailable_analysis`. Untested or unavailable results have missing
statistics and prevalence rather than a fabricated zero or negative finding.
They are excluded from density denominators and retain explicit support states.
No fallback to a discovery-only ranking is allowed for the diagnostic.
Legacy text-report rankings without an eligible-profile count retain missing
denominators and prevalence (for example, `2 / NA`), not an infinite proportion.

Standalone exports using `--tested-output` must supply `--profiles` from the
analyzed input. The exporter reads its trait schema strictly before writing any
output, including traits with no eligible tests. Missing, unreadable or malformed
schemas and tested traits absent from that schema fail instead of silently
publishing zero untested traits (#708). The normal stratification runner already
supplies this input. The schema establishes the declared trait vocabulary; it is
not a substitute for core provenance or an input-content checksum.

These quantities describe the relationship between reported taxon sets and
organism-level trait annotations. They do not measure physiological activity or
establish calibration bias. The diagnostic is descriptive; see
[the biological interpretation guide](../methods/biological_interpretation.md).

To rerun only stratification using existing core products, without rebuilding
default data products:

```bash
just analyze-strata-inputs \
  --profiles /path/to/study_trait_profiles_ora_wide.tsv \
  --experiments /path/to/experiments.tsv \
  --taxa-matrix /path/to/taxa_trait_matrix.tsv \
  --output-dir /path/to/new-stratification-output
```

`just test-stratification` exercises the actual analyzer, exporter, and diagnostic
with discoveries, a single-hit trait, tested null traits, unsupported traits,
zero-draw profiles, unpaired profiles, and an unavailable stratum.

## Production verification on 2026-10-19

Only the deterministic stratification stage was rerun against the preserved,
verified Ilyas inputs from core run `run-49a0c712ad3d4fc59d10804951fe511a`.
Fresh output is in `/private/tmp/mtp-issue664/.tmp/issue664-strata`.
The comparison can be reproduced with
`just audit-stratification-universe BEFORE_STRATA_DIR AFTER_STRATA_DIR`; its
implementation is tracked in
`scripts/adhoc/issue664_compare_stratification.py`.

| Quantity | Pooled | Feces | Feces excluded |
| --- | ---: | ---: | ---: |
| Eligible tests | 1,480,390 | 935,150 | 545,240 |
| Eligible traits | 317 | 317 | 317 |
| Traits with no discoveries | 11 | 14 | 17 |
| Traits with one discovery | 8 | 14 | 11 |
| Traits in the selected ranking | 298 | 289 | 289 |
| Corrected `METPO:1000487` prevalence | 363 / 4,670 | 187 / 2,950 | 176 / 1,720 |

Previously the example prevalence was 100% in every stratum because its selected
profile count was also used as its denominator. Its reference K is 394 in all
strata; the old profile-count substitutions were 9,835, 5,657, and 4,178.

The pooled diagnostic now includes all 317 eligible traits, including the 11
with no discoveries and eight with only one discovery. Its lowest-K decile has
21 of 32 traits with discoveries (65.625%); the other nine deciles remain 100%.
Those high descriptive discovery rates do not demonstrate calibrated inference
or establish biological effects. Calibration is a separate analysis.
