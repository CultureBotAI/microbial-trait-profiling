# Oral microbiomes and oral diseases: six-trait heatmap

[Main figure](oral_six_trait_heatmap.png) · [All 80 rows, vector PDF](oral_six_trait_heatmap_all_rows.pdf)

This is an oral-focused view of Jed's saved individual-condition signed contrasts
from [PR #249](https://github.com/CultureBotAI/microbial-trait-profiling/pull/249),
artifact commit `fb7afae18d37881ae071c256cbc93c8adb5cb8b8`. It retains the source's
six traits, effect values, support rules, Fisher p-values, and BH q-values.

## Inclusion and coverage

A condition/site row qualifies if its source body-site group is **Oral cavity /
dental OR its condition contains an explicitly oral diagnosis**. Diagnosis
components are split on commas and matched exactly, ignoring case. The complete
diagnosis list and source table SHA-256 are recorded in
[the provenance JSON](oral_six_trait_provenance.json). Broad grouped condition
categories are excluded; individual condition labels are preserved.

Nonspecific mucositis, Sjogren syndrome, and head/neck cancer are not assumed to be
oral diseases. They can qualify through oral sampling. This is an explicit label
selection for the saved source snapshot, not an ontology-wide disease classifier.

The selection contains 80 condition/site rows and 381 existing source trait-cell
records: 67 rows qualify through oral sampling alone and 13 through both criteria.
No additional nonoral-site rows qualify through oral disease in this snapshot.
The implementation also accepts oral diagnoses sampled at nonoral sites.

The main figure includes all 13 oral-disease rows, including eight with no
supported cells, plus 16 other rows with supported cells: 29 rows in total.
All 80 rows remain in the PDF, tables, and readable PNG pages:
[1](oral_six_trait_heatmap_page_01.png),
[2](oral_six_trait_heatmap_page_02.png),
[3](oral_six_trait_heatmap_page_03.png),
[4](oral_six_trait_heatmap_page_04.png).

## Values and significance

Each displayed value is the fraction of increased-direction profiles with a
positive trait z-score minus the corresponding fraction of decreased-direction
profiles. Colors use the same fixed -1 to +1 scale as the source. These are
profile fractions, not study-weighted estimates. The source Fisher tests are
profile-level diagnostics; within-study dependence and denominator sensitivity
are not resolved by this focused display. The separate methods review in
[issue #256](https://github.com/CultureBotAI/microbial-trait-profiling/issues/256)
remains applicable.

Dark borders indicate source BH q <= 0.05. **None of the 79 supported oral cells
meet that threshold.** The correction across the original 2,379 individual-view
cells is retained; filtering does not refit the tests or recompute BH. The six
traits were selected by a separate pooled analysis; their pooled significance
does not imply significance within these oral rows.

Source support requires at least three positive-z profiles in either direction
and at least two studies. Gray cells represent missing source records or cells
below the source display requirements; they do not establish a null effect.
The source is sparse: 381 recorded cells out of 80 x 6 possible cells. Row `max n`
is the largest source study count across that row's recorded traits, not a count
of independent observations in the Fisher test or the union of studies across
traits. The row audit names it `max_trait_study_count`. Study counts must not be summed
across overlapping condition rows.

[Cell audit](oral_six_trait_cells.tsv) preserves source counts, p/q values, flags,
inclusion reasons, and figure membership. The legacy `significant_in_heatmap`
column refers to source raw-p screening; borders use `selected_after_bh` instead.
[Row audit](oral_six_trait_rows.tsv) records support and main-figure membership.

## Regenerate

From the repository root:

```bash
just oral-six-trait-heatmap
```

The recipe reads the committed source contrasts and trait selection in the parent
directory, records the source artifact commit, and writes this directory's figures
and audit tables. It does not download data or rerun the upstream analysis.
