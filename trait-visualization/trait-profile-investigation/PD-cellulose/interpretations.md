# Parkinson's Disease and Cellulose Catabolic Process

## Summary

This folder focuses on the direction-specific connection between Parkinson's disease (PD), fecal profiles, and the cellulose catabolic process trait from the six-shared-traits analysis.

The main result is that cellulose catabolic process is not a broad PD marker across directions. It is specifically stronger in PD decreased fecal profiles.

## Direction-Specific Cytoscape Signal

In the individual-condition Cytoscape display edges, PD feces has `34` cellulose-positive increased profiles from `22` studies, compared with `69` cellulose-positive decreased profiles from `34` studies.

The decreased edge is stronger: decreased PD feces has a focus-minus-comparator z contrast of `1.610`, versus `0.659` for increased PD feces.

See `pd_cellulose_direction_specificity.png` and `pd_cellulose_cytoscape_highlight.png`.

## PD Six-Trait Context

Among the six traits in PD feces, cellulose has an increased-minus-decreased mean z difference of `-0.701`. Negative values indicate decreased-profile skew.

The companion lollipop plot shows cellulose as the clearest decreased-skewed PD trait, while moderately halophilic, facultatively aerobic, and arsenate reduction are increased-skewed.

See `pd_six_trait_direction_lollipop.png`.

## PD vs Other Fecal Profiles

Compared with other fecal profiles, PD decreased feces has a higher cellulose-positive z-rate (`58.0%` for PD vs `45.9%` for other feces in the same direction).

The decreased-direction PD-vs-other Fisher test gives an odds ratio of `1.63` and `p=0.0131`. In increased profiles, the PD cellulose-positive z-rate is `26.6%` versus `36.8%` for other feces.

The feces-only linear model includes PD, decreased direction, and the PD-by-decreased interaction. The cellulose z-score interaction term is `0.487` with study-clustered `p=0.0422`. This tests whether the decreased direction has an extra PD-associated cellulose shift beyond the main PD and direction effects.

See `pd_cellulose_enrichment_test.png` and `pd_cellulose_interaction_model.tsv`.

## Taxon-Level Drivers

Observed PD feces taxa were intersected with the KG-Microbe taxa-trait matrix, then limited to taxa carrying the cellulose catabolic process annotation. Counts are profile-level occurrences after deduplicating each taxon within each profile.

Top cellulose-carrier taxon labels in PD feces:

| taxon | rank | decreased profiles | increased profiles |
| --- | --- | ---: | ---: |
| Lachnospiraceae | family | 40 | 9 |
| Roseburia | genus | 33 | 4 |
| Oscillospiraceae | family | 12 | 15 |
| Bacillota | phylum | 15 | 5 |
| Blautia | genus | 16 | 2 |
| Pseudomonadota | phylum | 1 | 16 |
| Clostridiaceae | family | 13 | 2 |
| Eubacteriales | order | 7 | 7 |

The top decreased-skewed cellulose-carrier labels include Lachnospiraceae, Roseburia, Blautia, Clostridiaceae, Clostridium, Ruminococcus, and Lachnospira. These include a mix of genus/family/higher-rank reported labels, so they should be treated as taxonomic-unit drivers from BugSigDB curation rather than always as species-resolved drivers.

The taxon-trait heatmap checks whether these cellulose-carrier taxa also carry the other five shared traits. This helps identify whether the PD-cellulose signal is isolated or part of a broader multi-trait microbial phenotype pattern.

See `pd_cellulose_taxon_driver_bar.png`, `pd_cellulose_taxon_six_trait_heatmap.png`, and `pd_cellulose_taxon_drivers.tsv`.

## Generated Files

- `pd_cellulose_direction_specificity.png`: two-panel direction-specific summary.
- `pd_six_trait_direction_lollipop.png`: all-six-trait PD feces direction contrast.
- `pd_cellulose_cytoscape_highlight.png`: focused Cytoscape-style PD-cellulose highlight.
- `pd_cellulose_taxon_driver_bar.png`: observed cellulose-carrier taxa by direction.
- `pd_cellulose_taxon_six_trait_heatmap.png`: six-trait annotation context for top taxa.
- `pd_cellulose_enrichment_test.png`: feces-only model and group mean summary.
- `pd_cellulose_display_edge_summary.tsv`: PD display-edge audit data.
- `pd_cellulose_all_edge_summary.tsv`: PD all-edge audit data.
- `pd_six_trait_direction_summary.tsv`: PD feces six-trait z-score summary.
- `pd_cellulose_vs_other_feces_summary.tsv`: PD vs other feces comparison.
- `pd_cellulose_interaction_model.tsv`: study-clustered OLS model table.
- `pd_cellulose_model_group_summary.tsv`: model input group summaries.
- `pd_cellulose_taxon_drivers.tsv`: taxon-driver audit table.
- `pd_cellulose_taxon_trait_presence.tsv`: long-form top-taxon trait matrix.
