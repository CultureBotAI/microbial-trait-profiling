# Profile-Condition Enrichment Bars

These plots are profile-data companions to Ilyas's ORA enrichment infographic from `data/derived/enrichment_infographic.png` on the `enrichment-analysis-heatmap-cluster` branch.

They do not rerun the hypergeometric ORA. Instead, they use the selected profile heatmap z-score matrix. A selected trait is counted when its z-score is positive, meaning that trait's fraction is above the dataset-average fraction for the selected trait column.

## Inputs

- Matrix: `trait-visualization/study_trait_profile_clustered_heatmap_matrix.tsv`
- Metadata: `trait-visualization/study_trait_profile_study_metadata.tsv`
- Minimum condition profiles: 10
- Minimum condition studies: 2
- Broad metadata labels included: False
- Body-site filter: none
- Focus direction: decreased profiles
- Focus-condition sort: depletion-bias
- Minimum focus-direction profiles: 3

## Outputs

- `condition_depleted_profile_signal_diagram.png`: two-panel diagram combining the between-direction trait-signal difference and coverage.
- `condition_depleted_signal_vs_coverage_scatter.png`: uncapped scatterplot of the between-direction signal difference versus positive-signal profile coverage.
- `condition_depleted_signal_vs_coverage_scatter.tsv`: full data table used by the scatterplot.
- `condition_depleted_trait_signal_heatmap.png`: uncapped condition x selected-trait heatmap of each trait's positive-signal fraction difference between directions.
- `condition_depleted_trait_signal_heatmap_matrix.tsv`: full condition x trait matrix used by the heatmap.
- `condition_profile_trait_summary.tsv`: profile-data summary table used by the plots.

## Top Conditions By Profile Count

| condition | positive profiles | increased | decreased | studies |
| --- | ---: | ---: | ---: | ---: |
| Parkinson's disease | 271 | 141 | 130 | 57 |
| Colorectal cancer | 197 | 105 | 92 | 55 |
| Obesity | 171 | 89 | 82 | 37 |
| COVID-19 | 122 | 61 | 61 | 31 |
| Irritable bowel syndrome | 98 | 50 | 48 | 16 |

## Top Conditions By Average Positive Trait Depth

| condition | decreased mean | increased mean | direction-mean total | profiles | studies |
| --- | ---: | ---: | ---: | ---: | ---: |
| Chronic obstructive pulmonary disease | 46.00 | 48.00 | 94.00 | 13 | 4 |
| Glaucoma | 44.30 | 47.00 | 91.30 | 23 | 5 |
| Response to vaccine | 43.20 | 46.67 | 89.87 | 11 | 2 |
| Human immunodeficiency virus | 46.62 | 39.36 | 85.98 | 24 | 4 |
| Alzheimer's disease | 43.29 | 41.86 | 85.14 | 14 | 5 |

## Top Depleted-Profile Conditions

| condition | direction delta | focus mean | comparator mean | focus profiles | comparator profiles | studies |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Acute lymphoblastic leukemia | 18.85 | 45.00 | 26.15 | 6 | 13 | 5 |
| Diarrhea | 16.67 | 37.00 | 20.33 | 12 | 12 | 3 |
| Response to allogeneic hematopoietic stem cell transplant | 13.68 | 36.45 | 22.78 | 11 | 9 | 5 |
| Lung cancer | 11.01 | 42.82 | 31.81 | 11 | 16 | 8 |
| Endometrial cancer | 10.90 | 45.40 | 34.50 | 5 | 6 | 3 |

## Uncapped Follow-Up Figures

The scatterplot and condition x trait heatmap include all 47 focus-eligible condition rows after the shared filters, minimum focus-direction profile filter, and positive direction-delta filter. The heatmap includes all 60 selected trait columns from the input matrix.

The scatterplot is intended to separate large but sparsely supported direction deltas from conditions with broader profile/study support. `coverage_weighted_signal` in the scatter TSV is a descriptive product of the direction delta, square-root positive-signal profile count, and square-root study count; it is not used as a filter.

The heatmap decomposes the focus-direction signal by trait. Each cell is the focus-direction fraction of profiles with a positive selected-trait z-score minus the comparator-direction fraction for the same condition and trait.

## Interpretation Note

Use these as a profile-level analogue to the ORA condition panel. The two condition tables above show where profile coverage is deepest and where profiles have the broadest above-average selected trait signal. Because the metric is z-score based, it should be described as `positive selected trait signal` unless it is later replaced with formal ORA q-value enrichment counts.

The depleted-condition diagram ranks conditions by the mean number of selected traits with z > 0 in decreased profiles minus increased profiles while keeping the comparator visible as a guardrail. A positive direction delta means the focus direction has a higher mean positive selected-trait count than the comparator direction. Hatched rows have limited study support.

The diagram displays the top 15 conditions after filtering. The scatterplot and heatmap are not top-N capped.
