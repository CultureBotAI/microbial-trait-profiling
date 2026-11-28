# Parkinson's Disease and Cellulose Catabolic Process

## Summary

This folder follows the cellulose catabolic process trait (`METPO:2000103|GO:0030245`) of the six-shared-traits analysis in Parkinson's disease (PD) fecal profiles. Each profile is one curated BugSigDB signature. Its trait values summarize KG-Microbe annotations of the reported taxa, including annotations rolled up from their descendants; they are not measurements of cellulose degradation, gene expression or metabolites.

In PD fecal profiles the mean cellulose annotation z-score is lower in increased than in decreased profiles (increased minus decreased `-0.706`). A study-clustered model asks whether the decreased-minus-increased difference is larger in PD than in other fecal profiles: the extra PD effect in decreased profiles is `+0.494` cellulose z (95% CI `+0.018` to `+0.970`, `nominal p=0.0419`, 664 study clusters), and the interval excludes zero. This is the interaction's nominal p, one of 3 non-intercept model terms tested in the figure (8 terms across 2 outcomes in `pd_cellulose_interaction_model.tsv`), not corrected for the other conditions and traits examined; the corrected screen of every eligible fecal condition-trait contrast is described in `pd_condition_taxon_trait_network_methods.md`.

## Cellulose Edges by Signature Direction

In the individual-condition Cytoscape display edges, PD feces has `33` increased profiles with cellulose z > 0 from `22` studies, and `68` decreased profiles with cellulose z > 0 from `34` studies.

The focus-minus-comparator z contrast is `1.597` for the decreased edge and `0.666` for the increased edge.

See `pd_cellulose_direction_specificity.png` and `pd_cellulose_cytoscape_highlight.png`.

## PD Six-Trait Context

Among the six traits in PD feces, cellulose has an increased-minus-decreased mean z difference of `-0.706`. Negative values indicate decreased-profile skew.

In the companion lollipop plot the decreased-skewed traits, most negative first, are cellulose catabolic process, sphere shaped and coccus shaped; the increased-skewed traits, most positive first, are moderately halophilic, facultative (aerobe or anaerobe) and arsenate reduction.

See `pd_six_trait_direction_lollipop.png`.

## PD vs Other Fecal Profiles

Among decreased fecal profiles, cellulose z > 0 in `58.6%` of PD profiles and `46.1%` of other profiles; among increased fecal profiles, in `26.2%` and `37.0%`.

The decreased-direction PD-vs-other Fisher test gives an odds ratio of `1.66` and `nominal p=0.00928` (two-sided Mann-Whitney `nominal p=0.00819`; uncorrected across 2 directions x 2 tests = 4 nominal tests in `pd_cellulose_vs_other_feces_summary.tsv`). It treats the profiles of one study as independent, so the study-clustered model below is the better guide.

The feces-only linear model includes PD, decreased direction, and the PD-by-decreased interaction, with standard errors clustered by study. The cellulose z-score interaction term is `0.494` (95% CI `0.018` to `0.970`, `nominal p=0.0419`, uncorrected across 3 non-intercept model terms plotted in `pd_cellulose_enrichment_test.png` and 8 terms across 2 outcomes in `pd_cellulose_interaction_model.tsv`). It tests whether the decreased direction has an extra PD-associated shift in cellulose annotation beyond the main PD and direction effects.

See `pd_cellulose_enrichment_test.png` and `pd_cellulose_interaction_model.tsv`.

## Reported Taxa With the Cellulose Annotation

Observed PD feces taxa were intersected with the KG-Microbe taxa-trait matrix, then limited to taxa annotated with cellulose catabolic process, themselves or through annotated descendants. Counts are profile-level occurrences after deduplicating each taxon within each profile.

Cellulose-annotated taxon labels most often reported in PD feces:

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

Of the 18 labels plotted, those reported in more decreased than increased PD fecal profiles are Lachnospiraceae, Roseburia, Bacillota, Blautia, Clostridiaceae, Clostridium, Ruminococcus, Lachnospira, Bacteroides, Roseburia faecis and Bacteroidia. The labels mix genus, family and higher ranks, so they are BugSigDB reporting units rather than species-resolved organisms, and an annotation reached through descendants does not mean every member of the taxon carries it. The counts describe how often PD signatures report these labels; they do not show that the taxa produce the cellulose contrast.

The taxon-trait heatmap shows which of the other five shared traits the same taxa are annotated with. A higher-rank label can collect different annotations from different descendants, so co-annotation does not show that the traits occur together in one organism.

See `pd_cellulose_taxon_driver_bar.png`, `pd_cellulose_taxon_six_trait_heatmap.png`, and `pd_cellulose_taxon_drivers.tsv`.

## Generated Files

`just pd-cellulose-analysis` writes these files, after `just six-shared-trait-sankeys` has written the individual-condition Sankey edge tables it reads. The file names predate the current figure titles; the titles describe what each figure shows.

- `pd_cellulose_direction_specificity.png`: PD fecal profiles with cellulose z > 0 by signature direction, and PD vs other feces.
- `pd_six_trait_direction_lollipop.png`: all-six-trait PD feces direction contrast.
- `pd_cellulose_cytoscape_highlight.png`: focused Cytoscape-style PD-cellulose highlight.
- `pd_cellulose_taxon_driver_bar.png`: PD fecal profiles reporting each cellulose-annotated taxon, by direction.
- `pd_cellulose_taxon_six_trait_heatmap.png`: six-trait annotations of the most often reported cellulose-annotated taxa.
- `pd_cellulose_enrichment_test.png`: feces-only model and group mean summary.
- `pd_cellulose_display_edge_summary.tsv`: PD display-edge audit data.
- `pd_cellulose_all_edge_summary.tsv`: PD all-edge audit data.
- `pd_six_trait_direction_summary.tsv`: PD feces six-trait z-score summary.
- `pd_cellulose_vs_other_feces_summary.tsv`: PD vs other feces comparison.
- `pd_cellulose_interaction_model.tsv`: study-clustered OLS model table.
- `pd_cellulose_model_group_summary.tsv`: model input group summaries.
- `pd_cellulose_taxon_drivers.tsv`: audit table of the cellulose-annotated taxa reported in PD feces.
- `pd_cellulose_taxon_trait_presence.tsv`: long-form top-taxon trait matrix.
