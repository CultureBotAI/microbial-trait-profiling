# PD Trait Contributor Dot Heatmap

![Contributor heatmap](pd_trait_contributors_species_family_dot_heatmap.png)

The background row replaces the separate unresolved-share bars. Within each trait, paired red/blue circle area compares the increased/decreased propagated count; the grey outer area is the portion unresolved or ambiguous at the displayed rank. Contributor circle area uses the full cell denominator, so unresolved mass already reduces the colored contributor share rather than being discarded.

The heatmap retains `11` of `53` candidate taxa. A row must have at least one absolute increased-versus-decreased full-cell contributor share difference of `5%`. This removes rows containing only small directional changes.

The direction-difference table retains 149 of 318 taxon-trait pairs whose full-cell contributor share changes between directions. Equal pairs, including pairs absent in both directions, are omitted.

The right strip colors taxa by the profile-pool-normalized prevalence difference (`PD prevalence - other-feces prevalence`) in percentage points. This finite difference is used instead of a count ratio, so a small or zero other-feces count cannot create an infinite ratio. Raw count displays with a zero denominator, such as `1/0`, are rendered as `---`. The metric remains descriptive and does not account for study clustering.

## Direction-changing co-carriage pairs

![Trait co-carriage network](pd_trait_contributors_trait_cocarriage_network.png)

Only 3 of 15 trait pairs have different selected-taxon co-carriage counts between directions and are displayed:

| trait pair | increased | decreased | inc - dec | increased-only taxa | decreased-only taxa |
| --- | ---: | ---: | ---: | --- | --- |
| coccus shaped + moderately halophilic | 19 | 18 | +1 | Bifidobacteriaceae (31953) | none |
| facultatively aerobic + coccus shaped | 28 | 27 | +1 | Bifidobacteriaceae (31953) | none |
| facultatively aerobic + moderately halophilic | 28 | 27 | +1 | Bifidobacteriaceae (31953) | none |

## Aggregate direction patterns

The contributor counts are not valid independent observations for a significance test. As external context, the existing paired, study-weighted aggregate analysis has 95% intervals excluding zero for these traits:

| trait | mean increased - decreased z | 95% interval | pattern |
| --- | ---: | ---: | --- |
| moderately halophilic | 0.473 | [0.208, 0.738] | inc |
| facultatively aerobic | 0.412 | [0.045, 0.778] | inc |
| arsenate reduction | 0.516 | [0.043, 0.989] | inc |
| cellulose catabolic process | -0.698 | [-1.200, -0.196] | dec |

These intervals are aggregate directional context, not contributor-level tests, and the table does not add a new multiple-testing correction.

## Largest unresolved or ambiguous shares

| rank | trait | direction | unresolved share |
| --- | --- | --- | ---: |
| species | sphere shaped | inc | 80.3% |
| species | sphere shaped | dec | 70.6% |
| species | arsenate reduction | dec | 65.7% |
| species | moderately halophilic | inc | 63.9% |
| species | moderately halophilic | dec | 57.7% |
| species | cellulose catabolic process | inc | 45.2% |

High unresolved shares limit rank-specific attribution. They do not create extra colored contributor mass because every contributor share is divided by the full trait-direction count.

## PD-versus-other-feces occurrence context

Top profile-pool-normalized reported taxon labels (descriptive):

| rank | taxon | PD/other profiles | PD/other prevalence | prevalence difference |
| --- | --- | ---: | ---: | ---: |
| family | Lachnospiraceae (186803) | 49/307 | 19.37%/8.98% | +10.39 pp |
| family | Oscillospiraceae (216572) | 28/222 | 11.07%/6.49% | +4.57 pp |
| family | Bifidobacteriaceae (31953) | 16/106 | 6.32%/3.10% | +3.22 pp |
| family | Lactobacillaceae (33958) | 17/128 | 6.72%/3.74% | +2.98 pp |
| family | Streptococcaceae (1300) | 14/107 | 5.53%/3.13% | +2.40 pp |

## Counting and selection bias

- Counts are propagated taxon-trait assertions, not independent patients, samples, or original source assertions. Repeated experiments and descendant propagation can increase them.
- The displayed taxa are the union of the top five resolved contributors per trait-direction cell. Broadly annotated or frequently reported taxa therefore have more opportunities to enter and remain in the selected set.
- Species and family panels are separate resolution views. Parent and descendant labels must not be added as if they were disjoint contributors.
- Co-carriage means that the selected label contributes to both annotation cells; it does not establish within-sample expression, mechanism, or causal PD specificity.

Audit files:

- `pd_trait_contributors_cells.tsv`: long-form source dot cells.
- `pd_trait_contributors_direction_differences.tsv`: change-only taxon-trait direction pairs.
- `pd_trait_contributors_heatmap_taxa.tsv`: taxa passing the large-difference row filter.
- `pd_trait_contributors_selected_taxa.tsv`: selected taxa with normalized specificity context.
- `pd_trait_contributors_taxon_specificity.tsv`: PD-vs-other-feces counts and prevalence metrics.
- `pd_trait_contributors_trait_cocarriage_taxa.tsv`: direction-specific taxon-by-trait presence.
- `pd_trait_contributors_trait_cocarriage_all_edges.tsv`: all pairwise direction counts.
- `pd_trait_contributors_trait_cocarriage_edges.tsv`: only pairs that change by direction.
- `pd_trait_contributors_species_family_dot_heatmap.svg`: vector contributor heatmap.
- `pd_trait_contributors_trait_cocarriage_network.svg`: vector change-only network.
