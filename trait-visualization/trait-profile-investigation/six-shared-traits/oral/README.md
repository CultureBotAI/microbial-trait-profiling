# Oral microbiomes and oral diseases: six-trait heatmap

[Clustered figure](oral_six_trait_heatmap_clustered.png) ·
[Clustered vector PDF](oral_six_trait_heatmap_clustered.pdf) ·
[All 80 rows, clustered PDF](oral_six_trait_heatmap_clustered_all_rows.pdf) ·
[Original ordering](oral_six_trait_heatmap.png)

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

## Clustering both axes

The clustered view orders both condition/site rows and trait columns using the
source workflow's average-linkage hierarchical clustering on `1 - Spearman rho`,
with optimal leaf ordering. Distances use only jointly supported finite values,
with at least three shared values and nonconstant ranks in both vectors. For
rows those shared values are traits; for columns they are condition/site rows.
Clustering uses all 80 oral rows; the main figure retains the same 29-row subset.
Blank cells remain missing and are never filled with biological zeros.

All six traits and 16 rows have at least one defined comparison and enter their
respective trees. Of the main figure's 29 rows, 13 have insufficient overlap or
constant values and remain below the tree in alphabetical order. The full PDF
retains all 64 unclustered rows. Thus a displayed effect alone does not guarantee
enough information to cluster that row.

Among the 120 pairs within the row tree, nine lack a defined correlation and
receive the source method's fallback distance of 1. All 15 trait-pair correlations
are defined. That fallback is a layout convention, not measured similarity.
The trees summarize exploratory rank patterns; they do not establish biological
groups or new statistical significance. Neither effects nor source p/q values
are changed by reordering.

[Axis order and membership](oral_six_trait_cluster_order.tsv) records each label's
order, supported-value count, and number of defined peers. The
[pairwise comparison audit](oral_six_trait_cluster_pairs.tsv) records every pair
used within each tree, its overlap, whether its correlation is defined, and the
distance used. Provenance also records the method, leaf orders, and fallback
counts. Vectors with no defined peers are excluded from these pair tables and
identified explicitly in the membership table.

## Regenerate

From the repository root:

```bash
just oral-six-trait-heatmap
```

The recipe reads the committed source contrasts and trait selection in the parent
directory, records the source artifact commit, and writes this directory's figures
and audit tables. It does not download data or rerun the upstream analysis.
