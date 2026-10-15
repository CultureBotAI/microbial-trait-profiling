# Profile-Condition Enrichment Bars

These plots are profile-data companions to Ilyas's ORA enrichment infographic from `data/derived/enrichment_infographic.png` on the `enrichment-analysis-heatmap-cluster` branch.

They do not rerun the hypergeometric ORA. Instead, they use the selected profile heatmap z-score matrix. A selected trait is counted when its z-score is positive, meaning that trait's fraction is above the dataset-average fraction for the selected trait column.

## Inputs

- Matrix: `trait-visualization\study_trait_profile_clustered_heatmap_matrix.tsv`
- Metadata: `trait-visualization\study_trait_profile_study_metadata.tsv`
- Minimum condition profiles: 10
- Minimum condition studies: 2
- Broad metadata labels included: False
- Body-site filter: none
- Focus direction: decreased profiles
- Focus-condition sort: depletion-bias
- Minimum focus-direction profiles: 3

## Outputs

- `condition_positive_profile_count_bar.png`: condition x number of profiles with at least one positive selected trait signal, stacked by increased/decreased BugSigDB profile direction.
- `condition_average_positive_traits_bar.png`: condition x average number of positive selected traits, stacked as decreased-profile mean plus increased-profile mean in the same style as the ORA infographic condition panel.
- `condition_depleted_trait_signal_bar.png`: focus-condition bars ranked by the selected sort option; signal is `decreased_minus_increased_mean_positive_selected_trait_count`.
- `condition_depleted_profile_coverage_bar.png`: profile coverage for the same focus conditions, stacked by profile direction.
- `condition_depleted_profile_signal_diagram.png`: two-panel diagram combining direction-specific trait signal and coverage.
- `condition_depleted_signal_vs_coverage_scatter.png`: uncapped scatterplot of direction-specific signal difference versus positive-signal profile coverage.
- `condition_depleted_signal_vs_coverage_scatter.tsv`: full data table used by the scatterplot.
- `condition_depleted_trait_signal_heatmap.png`: uncapped condition x selected-trait heatmap showing which traits drive the direction-specific signal.
- `condition_depleted_trait_signal_heatmap_matrix.tsv`: full condition x trait matrix used by the heatmap.
- `condition_profile_trait_summary.tsv`: profile-data summary table used by the plots.

## Top Conditions By Profile Count

| condition | positive profiles | increased | decreased | studies |
| --- | ---: | ---: | ---: | ---: |
| Parkinson's disease | 373 | 192 | 181 | 76 |
| Colorectal cancer | 290 | 159 | 131 | 67 |
| Obesity | 251 | 134 | 117 | 56 |
| COVID-19 | 192 | 99 | 93 | 45 |
| Response to transplant | 170 | 105 | 65 | 18 |

## Top Conditions By Average Positive Trait Depth

| condition | decreased mean | increased mean | direction-mean total | profiles | studies |
| --- | ---: | ---: | ---: | ---: | ---: |
| Glaucoma | 42.64 | 45.33 | 87.98 | 32 | 6 |
| Human immunodeficiency virus | 48.27 | 38.90 | 87.17 | 21 | 2 |
| Diabetic retinopathy | 46.00 | 39.60 | 85.60 | 10 | 2 |
| Epilepsy | 39.25 | 42.44 | 81.69 | 24 | 5 |
| Attention deficit-hyperactivity disorder | 43.00 | 37.67 | 80.67 | 22 | 3 |

## Top Depleted-Profile Conditions

| condition | direction delta | focus mean | comparator mean | focus profiles | comparator profiles | studies |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Bronchiectasis | 20.00 | 42.67 | 22.67 | 6 | 6 | 2 |
| Arthritis | 18.60 | 40.60 | 22.00 | 5 | 5 | 6 |
| Acute lymphoblastic leukemia | 13.61 | 38.80 | 25.19 | 10 | 16 | 7 |
| Gastric adenocarcinoma | 13.06 | 36.86 | 23.80 | 7 | 10 | 3 |
| Chlamydia trachomatis infectious disease | 12.37 | 23.17 | 10.80 | 6 | 10 | 2 |

## Uncapped Follow-Up Figures

The scatterplot and condition x trait heatmap include all 75 focus-eligible condition rows after the shared filters, minimum focus-direction profile filter, and positive direction-delta filter. The heatmap includes all 60 selected trait columns from the input matrix.

The scatterplot is intended to separate large but sparsely supported direction deltas from conditions with broader profile/study support. `coverage_weighted_signal` in the scatter TSV is a descriptive product of the direction delta, square-root positive-signal profile count, and square-root study count; it is not used as a filter.

The heatmap decomposes the focus-direction signal by trait. Each cell is the focus-direction fraction of profiles with a positive selected-trait z-score minus the comparator-direction fraction for the same condition and trait.

## Interpretation Note

Use these as a profile-level analogue to the ORA condition panel. The count plot shows where profile coverage is deepest; the average-depth plot shows where profiles have the broadest above-average selected trait signal. Because the metric is z-score based, it should be described as `positive selected trait signal` unless it is later replaced with formal ORA q-value enrichment counts.

The depleted-condition diagram ranks conditions by how much stronger selected traits are in decreased profiles than increased profiles while keeping the comparator visible as a guardrail. A positive direction delta means the focus direction has a higher mean positive selected-trait count than the comparator direction. Hatched rows have limited study support.

The bar plots display the top 15 conditions after filtering. The scatterplot and heatmap are not top-N capped.
