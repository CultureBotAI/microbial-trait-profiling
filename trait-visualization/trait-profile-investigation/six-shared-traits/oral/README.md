# Oral microbiomes and oral diseases: six-trait heatmap

[Main figure](oral_six_trait_heatmap.png) · [All 117 rows, vector PDF](oral_six_trait_heatmap_all_rows.pdf)

This is an oral-focused view of the saved individual-condition signed contrasts
in the parent directory. It retains the source's six traits, effect values,
support and display rules, Fisher p-values, and BH q-values. The
[provenance JSON](oral_six_trait_provenance.json) records the source artifact
commit and the SHA-256 hashes of the exact contrast and trait-selection files
used for this rendering.

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

The selection contains 117 condition/site rows and 522 existing source trait-cell
records: 103 rows qualify through oral sampling alone and 14 through both criteria.
No additional nonoral-site rows qualify through oral disease in this snapshot.
The implementation also accepts oral diagnoses sampled at nonoral sites.

The main figure includes all 14 oral-disease rows, including nine with no
displayed cells, plus 19 other rows with displayed cells: 33 rows in total.
All 117 rows remain in the PDF, tables, and readable PNG pages:
[1](oral_six_trait_heatmap_page_01.png),
[2](oral_six_trait_heatmap_page_02.png),
[3](oral_six_trait_heatmap_page_03.png),
[4](oral_six_trait_heatmap_page_04.png),
[5](oral_six_trait_heatmap_page_05.png).

## Values and significance

Each displayed value is the fraction of increased-direction profiles with a
positive trait z-score minus the corresponding fraction of decreased-direction
profiles. Colors use the same fixed -1 to +1 scale as the source. These are
profile fractions, not study-weighted estimates. The source Fisher tests are
profile-level diagnostics; within-study dependence and denominator sensitivity
are not resolved by this focused display. The separate methods review in
[issue #256](https://github.com/CultureBotAI/microbial-trait-profiling/issues/256)
remains applicable.

Dark borders indicate BH-adjusted q ≤ 0.05. The oral subset has 83 displayed
cells and no BH-selected cells. **None of the 83 displayed oral cells meet that
threshold.** The correction across the full 3,144 individual-view cells is
retained; filtering does not refit the tests or recompute BH. The single
BH-selected cell in the full individual-condition view is fecal.

| Cohort | Displayed cells | Nominal Fisher p < 0.01 | BH-adjusted q ≤ 0.05 |
| --- | ---: | ---: | ---: |
| All body sites (full individual view) | 632 | 29 | 1 |
| Oral cavity / dental | 83 | 2 | 0 |

Two oral cells meet nominal Fisher p < 0.01:
`Colorectal cancer | Oral cavity / dental` × `coccus shaped`
(METPO:1000668; effect +0.714, p = 0.003394, q = 0.587927) and
`Parkinson's disease | Oral cavity / dental` × `moderately halophilic`
(METPO:1000623; effect +1.0, p = 0.007937, q = 0.881907). Their nominal flag
`nominal_p_value_selected` is true, but neither has a dark BH border because
q > 0.05.

### Terminology and analysis distinctions

- **Signed contrast significance borders (Approach 4)**: A dark border strictly marks BH-adjusted q ≤ 0.05 (`selected_after_bh`) evaluated over the full 3,144-cell family. Nominal Fisher p < 0.01 is audit-only and does not trigger a border.
- **Trait-pattern display filtering (Approach 2)**: In trait-pattern figures, nominal Fisher p < 0.01, n ≥ 3, and s ≥ 2 are exploratory display filters for coloring cells. This display filtering rule must not be conflated with the signed-contrast BH border rule.
- **Pooled trait selection**: The six traits were prioritized by a separate pooled across-disease association analysis; their pooled significance does not imply significance within these individual condition rows or the oral subset.

Source support requires at least three positive-z profiles in either direction
and at least two studies. A displayed cell must also have a finite Fisher test
and an absolute increased-minus-decreased positive-z fraction strictly greater
than 0.01. Passing support alone does not make a cell visible. Gray cells
represent missing source records, insufficient support, an unavailable test,
or an effect at or below that cutoff; they do not establish a null effect.
The source is sparse: 522 recorded cells out of 117 x 6 possible cells. Row `max n`
is the largest source study count across that row's recorded traits, not a count
of independent observations in the Fisher test or the union of studies across
traits. The row audit names it `max_trait_study_count`. Study counts must not be summed
across overlapping condition rows.

[Cell audit](oral_six_trait_cells.tsv) preserves source counts, p/q values, flags,
inclusion reasons, and figure membership. `nominal_p_value_selected` explicitly
records nominal Fisher p < 0.01; `significant_in_heatmap` is retained as a
deprecated compatibility alias. Borders use `selected_after_bh`.
[Row audit](oral_six_trait_rows.tsv) records displayed-trait counts and
main-figure membership; its legacy `supported_traits` column counts displayed
traits after the effect cutoff.
[Provenance JSON](oral_six_trait_provenance.json) tracks `nominal_p_value_cell_count`
alongside deprecated `significant_cell_count`, and `bh_selected_cell_count`
alongside `bh_significant_displayed_cells`.

## Regenerate

From the repository root:

```bash
just oral-six-trait-heatmap
```

The recipe reads the saved source contrasts and trait selection in the parent
directory, records their hashes and source artifact commit, and writes this
directory's figures and audit tables. Modified source artifacts are marked
`uncommitted` until committed and regenerated. It does not download data or rerun
the upstream analysis.
