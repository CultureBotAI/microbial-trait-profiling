# Parkinson's Disease and Cellulose Catabolic Process

## Summary

This folder focuses on the direction-specific connection between Parkinson's disease (PD), fecal profiles, and the cellulose catabolic process trait from the six-shared-traits analysis.

The main cellulose result is that the cellulose catabolic process signal is concentrated in PD decreased fecal profiles rather than acting as a broad PD marker across directions. The within-direction PD-versus-other-feces comparisons are nominally significant, but the study-clustered PD-by-direction interactions do not cross the conventional `p<0.05` threshold. In the separate global FDR screen across 56 conditions and all six traits, PD-cellulose does not pass `q<=0.05`; only PD-moderately-halophilic does. The result therefore supports a decreased-profile concentration of the PD cellulose signal, not a definitive PD-specific interaction or an FDR-confirmed PD-cellulose biomarker.

## Direction-Specific Cytoscape Signal

In the individual-condition Cytoscape display edges, PD feces has `35` cellulose-positive increased profiles from `22` studies, compared with `67` cellulose-positive decreased profiles from `32` studies.

The decreased edge is stronger: decreased PD feces has a focus-minus-comparator z contrast of `1.581`, versus `0.752` for increased PD feces.

See `pd_cellulose_direction_specificity.png` and `pd_cellulose_cytoscape_highlight.png`.

## PD Six-Trait Context

Among the six traits in PD feces, cellulose has a profile-level increased-minus-decreased mean z difference of `-0.609`. Negative values indicate decreased-profile skew. The species-influence workflow uses a different estimand: paired experiments are averaged within study and studies are weighted equally.

| trait | study-weighted inc-dec z | 95% CI | aggregate direction |
| --- | ---: | ---: | --- |
| arsenate reduction | `+0.516` | `0.043` to `0.989` | increased, CI excludes zero |
| moderately halophilic | `+0.473` | `0.208` to `0.738` | increased, CI excludes zero |
| facultatively aerobic | `+0.412` | `0.045` to `0.778` | increased, CI excludes zero |
| coccus shaped | `+0.138` | `-0.260` to `0.536` | positive estimate, uncertain |
| sphere shaped | `+0.133` | `-0.230` to `0.495` | positive estimate, uncertain |
| cellulose catabolic process | `-0.698` | `-1.200` to `-0.196` | decreased, CI excludes zero |

The study-weighted summary therefore supports increased-profile skew for moderately halophilic, facultatively aerobic, and arsenate reduction, and decreased-profile skew for cellulose catabolic process. Coccus shaped and sphere shaped remain directionally uncertain. These are aggregate trait contrasts, not species-level contribution tests, and they do not replace the profile-level `-0.609` cellulose result.

See `pd_six_trait_direction_lollipop.png` and `pd_species_six_traits_influence_heatmap_traits.tsv`.

## Cross-Condition FDR Evidence Landscape

The condition-taxon-trait evidence landscape shows every eligible result while preserving the global multiple-testing boundary. It includes `56` human fecal conditions with at least three studies containing both increased and decreased signatures, `264` eligible condition-taxon comparisons, and `336` condition-trait comparisons. Genus and species IDs remain separate.

Condition-taxon direction is tested with an exact McNemar/binomial test on discordant paired studies. Condition-trait contrasts are tested with a two-sided Wilcoxon signed-rank test across paired study-level increased-minus-decreased z scores. Benjamini-Hochberg correction is applied separately across the complete eligible taxon family and the complete trait family.

The figure uses three evidence tiers:

- Filled points are FDR-significant at `q<=0.05`.
- Hollow amber points have `0.05<q<=0.20` and are exploratory, not significant.
- Grey points have `q>0.20` and show the remaining tested background.

Only one result passes the global FDR threshold: Parkinson's disease to moderately halophilic (`mean contrast=+0.530`, raw `p=0.000143`, `q=0.048`, `n=45` paired studies). Four results are retained only as exploratory follow-up hypotheses:

| evidence tier | condition | taxon or trait | effect | raw p | global q | study support |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| significant | Parkinson's disease | moderately halophilic | `+0.530` | `0.000143` | `0.048` | `45` |
| exploratory | Parkinson's disease | *Roseburia* genus | `-0.333` directional study fraction | `0.000729` | `0.192` | `20` reporting studies |
| exploratory | Parkinson's disease | arsenate reduction | `+0.506` | `0.00159` | `0.164` | `45` |
| exploratory | COVID-19 | cellulose catabolic process | `-1.307` | `0.00195` | `0.164` | `11` |
| exploratory | colorectal cancer | cellulose catabolic process | `-0.821` | `0.00194` | `0.164` | `24` |

No condition-taxon comparison passes `q<=0.05`; consequently, no genus, species, or KG-Microbe taxon-to-trait provenance link is presented as an inferential discovery. The significant-only node and edge tables contain only PD to moderately halophilic. The raw p-values and earlier confidence intervals provide descriptive context, but they do not supersede the global q-values.

See `pd_condition_taxon_trait_network.png`, `pd_condition_taxon_trait_network_methods.md`, `pd_condition_taxon_trait_network_condition_taxon_tests.tsv`, and `pd_condition_taxon_trait_network_condition_trait_tests.tsv`.

## PD vs Other Fecal Profiles

Compared with other fecal profiles, PD decreased feces has a higher cellulose-positive z-rate (`54.5%` for PD vs `44.8%` for other feces in the same direction).

The decreased-direction PD-vs-other Fisher test gives an odds ratio of `1.48` and nominal `p=0.0391`; the corresponding raw-profile z-score distribution comparison gives Mann-Whitney `p=0.0232`. In increased profiles, the PD cellulose-positive z-rate is `26.9%` versus `36.1%` for other feces.

The feces-only linear model includes PD, decreased direction, and the PD-by-decreased interaction. The cellulose z-score interaction term is `0.400` with study-clustered `p=0.086`; the positive-z-rate interaction is `0.189` with `p=0.053`. These terms test whether the decreased direction has an extra PD-associated cellulose shift beyond the main PD and direction effects. Both are suggestive, but neither meets the conventional `p<0.05` threshold.

See `pd_cellulose_enrichment_test.png` and `pd_cellulose_interaction_model.tsv`.

## Taxon-Level Associations

Observed PD feces taxa were intersected with the KG-Microbe taxa-trait matrix, then limited to taxa carrying the cellulose catabolic process annotation. Counts are profile-level occurrences after deduplicating each taxon within each profile.

Top cellulose-carrier taxon labels in PD feces:

| taxon | rank | decreased profiles | increased profiles |
| --- | --- | ---: | ---: |
| Lachnospiraceae | family | 39 | 10 |
| Roseburia | genus | 31 | 4 |
| Oscillospiraceae | family | 13 | 15 |
| Bacillota | phylum | 14 | 5 |
| Blautia | genus | 15 | 2 |
| Pseudomonadota | phylum | 1 | 15 |
| Clostridiaceae | family | 13 | 1 |
| Eubacteriales | order | 7 | 7 |

The top decreased-skewed cellulose-carrier labels include Lachnospiraceae, Roseburia, Blautia, Clostridiaceae, Clostridium, Ruminococcus, and Lachnospira. These include a mix of genus, family, and higher-rank reported labels, so they should be treated as recurrent taxonomic-unit associations from BugSigDB curation rather than as species-resolved or causal drivers.

The taxon-trait heatmap checks whether these cellulose-carrier taxa also carry the other five shared traits. This helps identify whether the PD-cellulose signal is isolated or part of a broader multi-trait microbial phenotype pattern.

See `pd_cellulose_taxon_driver_bar.png`, `pd_cellulose_taxon_six_trait_heatmap.png`, and `pd_cellulose_taxon_drivers.tsv`.

### Direction-Exclusive Reports

Among `60` cellulose-annotated taxon labels observed in PD feces, `30` are reported in only one direction. Nine of those direction-exclusive labels occur in at least two studies; the other `21` are supported by one study each.

*Roseburia faecis* has the strongest repeated-study species-level decreased-only pattern in this shortlist: `6/123` decreased PD profiles from `5/46` decreased-direction studies, compared with no increased PD profiles. Its decreased reporting rate is `4.9%` in PD feces versus `0.7%` in the non-PD fecal background. This is a selected descriptive candidate, not a formally established biomarker; formal inference would need study clustering and multiple-testing control across the full candidate set.

The separate global direction screen tests the *Roseburia* genus, not *R. faecis* specifically. The genus favors decreased reporting in `17` discordant studies versus `2` favoring increased, but its global `q=0.192`; it is therefore exploratory rather than significant. This genus result must not be transferred to every *Roseburia* species.

See `pd_cellulose_direction_exclusive_taxa_context.md`.

## Species-by-Six-Trait Influence Diagnostic

The influence workflow evaluates reported species across all six PD lollipop traits. It removes each species' complete non-MediaDive annotation vector from profiles where that species was reported, recomputes trait fractions against the unchanged human-corpus reference, forms paired increased-minus-decreased experiment contrasts, averages experiments within study, and weights studies equally. Species whose six cells all have fewer than three supporting studies are omitted, *Mediterraneibacter torques* is explicitly excluded, and species with no mean influence greater than `+0.01` or less than `-0.01` are omitted. Five species remain, with `30` audited species-trait cells and nine cells beyond the display threshold.

All `30` audited cells for the five retained species have intervals that cross zero and are marked `mixed`; none receives a supported `inc` or `dec` classification. Thus, the aggregate trait directions above are not directionally attributable to any single displayed species in this diagnostic. The cell TSV retains the exact value for every retained species-trait combination and uses `displayed_in_heatmap` to identify the nine threshold-passing cells.

Rows are sorted from the strongest BugSigDB increased skew to the strongest decreased skew using the PD log2 increased/decreased profile-count ratio. Infinite-ratio ties are ordered by the signed profile-count difference, then total direction support and taxon name. This ordering controls display only; it does not change influence values or literature comparisons.

The first right-hand annotation shows BugSigDB increased/decreased recurrence. The second uses red fill for literature-increased and blue fill for literature-decreased. Directional literature and BugSigDB recurrence are compared with one-position Hamming distance: `0` means similar and receives a green border; `1` means different and receives a red border. Four of the five retained species have distance `0`. *Escherichia coli* has distance `1` because its BugSigDB recurrence is decreased (`1/2`) while the matched literature direction is increased.

The x-axis `baseline Δz` is the full-data, study-weighted paired PD increased-minus-decreased trait contrast before taxon removal. Each heatmap cell is that baseline minus the contrast recomputed after taxon removal. Positive cells make the positive contrast larger; negative cells make the negative contrast larger. The values do not sum to the baseline because trait fractions are nonlinear and taxon annotations overlap.

Outlined cells include direct removal of the displayed trait; unoutlined cells show the denominator-mediated effect of removing a species that lacks a direct annotation for that trait. For example, the *E. coli*–arsenate reduction cell is directly annotated and has influence `+0.1012`, so removing *E. coli* weakens the positive arsenate contrast by that amount in this diagnostic. Its 95% interval crosses zero (`-0.1372` to `+0.3395`), so the direction remains `mixed`. This is contribution/influence evidence, not causal attribution or a strain-level measurement.

See `pd_species_six_traits_influence_heatmap.png`, `pd_species_six_traits_influence_heatmap.md`, and `pd_species_six_traits_influence_heatmap_cells.tsv`.

### Species + Genus + Family View

The expanded view applies the same strict magnitude display filter and excludes *M. torques*, but lowers the exploratory support threshold to two studies. It selects up to 30 literature-matched or best-supported rows, then groups the y-axis into species, genus, and family blocks separated by dark horizontal rules. The current figure contains eight species, 14 genera, and eight families, with `71` visible cells beyond `±0.01`. Genus and family rows can mix multiple descendant species and remain hypothesis-generating taxonomic roll-ups.

See `pd_exploratory_taxa_six_traits_influence_heatmap.png`, `pd_exploratory_taxa_six_traits_influence_heatmap.md`, and `pd_exploratory_taxa_six_traits_influence_heatmap_cells.tsv`.

### Agreement With Published PD Directions

| displayed species | BugSigDB inc/dec | literature value | comparison level | Hamming result |
| --- | ---: | --- | --- | --- |
| *Streptococcus salivarius* | `4/0` increased | increased | *Streptococcus* genus proxy | `0`, similar |
| *Methanobrevibacter smithii* | `3/0` increased | increased | exact species term in a genus reference | `0`, similar |
| *Escherichia coli* | `1/2` decreased | increased | *Escherichia-Shigella* genus-group proxy | `1`, different |
| *Roseburia faecis* | `0/6` decreased | decreased | genus plus species-level evidence | `0`, similar |
| *Faecalibacterium prausnitzii* | `0/13` decreased | decreased | genus/species evidence | `0`, similar |

The *E. coli* comparison uses broader genus-group literature rather than a pooled species-specific estimate, so its Hamming distance is a descriptive sensitivity finding rather than evidence that either source is wrong.

The DOI-linked evidence statements and taxon-by-taxon interpretations are in `pd_species_six_traits_reference_interpretation.md` and `pd_species_six_traits_reference_table.tsv`.

### Methods: Literature Review and Direction Comparison

The direction reference is a targeted narrative review updated `2026-09-22`, not a systematic review or meta-analysis. PubMed, PubMed Central, and publisher records were searched using Parkinson disease, gut microbiome, each taxon name or synonym, and the terms increased, decreased, enriched, or depleted. Eligible evidence comprised human PD-versus-control fecal or intestinal microbiome studies and cross-cohort human meta-analyses with an explicit taxon direction. Animal-only, intervention-only, non-gut, and severity-only findings were excluded from increased/decreased calls.

An exact species result was preferred. When unavailable, a genus or family result was retained as an explicitly labeled parent-level proxy; conflicting eligible evidence at the matched level is classified as `variable`, and association-only evidence is not treated as directional. BugSigDB direction is defined from curated PD fecal profile recurrence: increased when increased-profile count exceeds decreased-profile count, decreased for the reverse, and balanced for a tie. Directional literature and BugSigDB calls are encoded as one binary categorical position and compared with Hamming distance (`0` similar, `1` different). Balanced, variable, associated, and unlisted cases have no defined distance. This is not a statistical agreement test and does not resolve cohort, medication, geography, sequencing, or taxonomy-version heterogeneity.

## Contributor Specificity and Trait Co-Carriage

The contributor analysis selects `28` species and `25` families across the six PD traits. To remove rows containing only small dots, the heatmap retains `11/53` taxa—`3` species and `8` families—with at least one absolute increased-versus-decreased full-cell contributor-share difference of `5` percentage points. The complete direction-difference audit still retains `149/318` taxon-trait pairs whose share changes between directions. Equal pairs, including the `169` pairs absent in both directions, are omitted. Its first row uses paired red/blue circles for the directional background and a grey outer portion for unresolved or ambiguous rank mass; contributor circles remain fractions of the full count, so unresolved mass is not silently renormalized away.

The specificity strip uses `253` PD fecal profiles and `3,419` other fecal profiles as the respective denominators and colors by the normalized prevalence difference (`PD prevalence - other-feces prevalence`) rather than a count ratio. The color scale is clipped at `+/-5` percentage points so one extreme row does not compress the remaining colors. Raw count displays with a zero denominator, such as `1/0`, are rendered as `---`; thus an undefined or infinite ratio cannot dominate the color scale. The prevalence difference remains descriptive rather than inferential because profiles and experiments repeat within studies.

Only `3/15` trait co-carriage pairs change across directions, and every change is one additional increased-direction family: Bifidobacteriaceae. Facultatively aerobic plus coccus shaped changes from `27` decreased to `28` increased shared taxa, facultatively aerobic plus moderately halophilic also changes from `27` to `28`, and coccus shaped plus moderately halophilic changes from `18` to `19`. The other `12` pairs are retained in the all-edges audit but omitted from the network.

No contributor-level or co-carriage-edge significance test is supported by these propagated counts: they are neither independent samples nor original assertions, and the taxa were selected using the same counts displayed. In the separate paired, study-weighted aggregate analysis, 95% intervals exclude zero for moderately halophilic (`0.473`, interval `0.208` to `0.738`), facultatively aerobic (`0.412`, `0.045` to `0.778`), and arsenate reduction (`0.516`, `0.043` to `0.989`) toward increased profiles, and for cellulose catabolic process (`-0.698`, `-1.200` to `-0.196`) toward decreased profiles. Those intervals are trait-level contextual patterns, not contributor-level tests or globally adjusted discoveries. In the complete cross-condition FDR screen, only PD-moderately-halophilic passes `q<=0.05`.

Species-level attribution is particularly limited where unresolved or ambiguous mass is high: sphere shaped is `80.3%` unresolved in increased and `70.6%` in decreased profiles; decreased arsenate reduction is `65.7%` unresolved. Count selection also favors frequently reported or broadly annotated taxa, and species/family parent-descendant labels cannot be summed as disjoint contributors.

See `pd_trait_contributors_species_family_dot_heatmap.png`, `pd_trait_contributors_trait_cocarriage_network.png`, and `pd_trait_contributors.md`.

## Interpretation Boundaries

- Profile counts are curated differential-abundance signature occurrences, not patient counts or direct prevalence measurements.
- Exact reported labels span multiple taxonomic ranks. Counts across ranks must not be added, and broad annotations can reflect trait-bearing descendants rather than every member of a higher-rank taxon.
- A zero count in one direction means that the exact label was not reported in that curated direction; it does not establish biological absence.
- Direction-exclusive and contributor tables are selected descriptive views. They are not substitutes for study-clustered, multiple-testing-adjusted inference.
- Species-removal estimates change both trait numerators and the nonlinear denominator. They measure model influence under the local annotation structure, not biological causation.

## Generated Files

- `pd_cellulose_direction_specificity.png`: two-panel direction-specific summary.
- `pd_six_trait_direction_lollipop.png`: all-six-trait PD feces direction contrast.
- `pd_condition_taxon_trait_network.png`: all-test FDR-tiered condition-taxon-trait evidence landscape.
- `pd_condition_taxon_trait_network.svg`: vector evidence landscape.
- `pd_condition_taxon_trait_network_methods.md`: statistical universe, thresholds, result counts, and interpretation boundary.
- `pd_condition_taxon_trait_network_condition_eligibility.tsv`: all conditions and eligibility or exclusion reasons.
- `pd_condition_taxon_trait_network_condition_taxon_tests.tsv`: complete condition-taxon test audit.
- `pd_condition_taxon_trait_network_condition_trait_tests.tsv`: complete condition-trait test audit.
- `pd_condition_taxon_trait_network_study_trait_contrasts.tsv`: study-level trait contrasts used by the Wilcoxon tests.
- `pd_condition_taxon_trait_network_nodes.tsv`: significant-only node table.
- `pd_condition_taxon_trait_network_edges.tsv`: significant-only edge table.
- `pd_condition_taxon_trait_network_taxon_selection.tsv`: FDR-selected taxon table, currently empty.
- `pd_cellulose_cytoscape_highlight.png`: focused Cytoscape-style PD-cellulose highlight.
- `pd_cellulose_taxon_driver_bar.png`: observed cellulose-carrier taxa by direction.
- `pd_cellulose_taxon_six_trait_heatmap.png`: six-trait annotation context for top taxa.
- `pd_cellulose_direction_exclusive_taxa_context.md`: direction-exclusive cellulose-carrier shortlist with PD and non-PD fecal denominators.
- `pd_trait_taxa_lookup.md`: reported carrier lookup for all six PD lollipop traits.
- `pd_species_six_traits_influence_heatmap.png`: study-weighted species-by-six-trait influence diagnostic.
- `pd_species_six_traits_influence_heatmap_cells.tsv`: long-form influence estimates, intervals, support, and classifications.
- `pd_species_six_traits_influence_heatmap_taxa.tsv`: heatmap row annotations and occurrence summaries.
- `pd_species_six_traits_influence_heatmap_specificity.tsv`: exact-label PD-versus-other-feces occurrence context.
- `pd_species_six_traits_influence_heatmap_traits.tsv`: study-weighted overall trait contrasts.
- `pd_species_six_traits_reference_table.tsv`: BugSigDB/literature direction, match scope, six-trait annotations, and strongest influence per displayed species.
- `pd_species_six_traits_reference_interpretation.md`: DOI-linked agreement summary, literature-review methods, and taxon details.
- `pd_taxa_literature_reference.tsv`: curated literature directions, evidence statements, DOI citations, and match terms.
- `pd_exploratory_taxa_six_traits_influence_heatmap.png`: thresholded influence view grouped into species, genus, and family blocks.
- `pd_exploratory_taxa_six_traits_influence_heatmap_cells.tsv`: rank-expanded influence estimates and display-filter flags.
- `pd_exploratory_taxa_six_traits_influence_heatmap_taxa.tsv`: selected rank-grouped rows, occurrence summaries, and literature comparisons.
- `pd_trait_contributors_species_family_dot_heatmap.png`: selected species/family contributor view with specificity strip.
- `pd_trait_contributors_direction_differences.tsv`: change-only taxon-trait direction pairs.
- `pd_trait_contributors_heatmap_taxa.tsv`: taxa retained by the large direction-difference filter.
- `pd_trait_contributors_taxon_specificity.tsv`: raw occurrence counts, display values, and normalized prevalence metrics.
- `pd_trait_contributors_trait_cocarriage_all_edges.tsv`: all 15 direction-specific co-carriage pairs.
- `pd_trait_contributors_trait_cocarriage_network.png`: direction-changing co-carriage pairs among selected contributors.
- `pd_trait_contributors.md`: contributor-analysis interpretation and audit-file index.
- `pd_cellulose_enrichment_test.png`: feces-only model and group mean summary.
- `pd_cellulose_display_edge_summary.tsv`: PD display-edge audit data.
- `pd_cellulose_all_edge_summary.tsv`: PD all-edge audit data.
- `pd_six_trait_direction_summary.tsv`: PD feces six-trait z-score summary.
- `pd_cellulose_vs_other_feces_summary.tsv`: PD vs other feces comparison.
- `pd_cellulose_interaction_model.tsv`: study-clustered OLS model table.
- `pd_cellulose_model_group_summary.tsv`: model input group summaries.
- `pd_cellulose_taxon_drivers.tsv`: taxon-driver audit table.
- `pd_cellulose_taxon_trait_presence.tsv`: long-form top-taxon trait matrix.
