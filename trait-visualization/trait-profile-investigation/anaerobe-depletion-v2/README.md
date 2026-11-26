# Anaerobe v2 Depletion Investigation

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices. Guild-trait z-scores are computed from the retained counts of data/derived/study_trait_profiles_wide.tsv on the heatmap's scale, so a selected trait outside the heatmap's top traits is still scored (#2399).

## Scope

- Selected traits: anaerobic (METPO:1000603), obligately anaerobic (METPO:1000607)
- Target taxa: Faecalibacterium prausnitzii, Anaerobutyricum hallii, Lachnospira eligens, Agathobacter rectalis, Roseburia intestinalis, Ruminococcus bromii, Roseburia inulinivorans
- Formal KG-Microbe taxon universe: 1,908 taxa with at least one core oxygen-tolerance annotation among anaerobic, obligately anaerobic
- Trait selection audit: `guild_trait_selection.tsv` records the selected trait IDs, short rationale, and KG-Microbe taxon counts for each selected trait.
- Context: all body sites, with body-site and condition summaries to separate gut, oral, airway, urogenital, and other patterns.

## Plot Guild Specificity

| plot | guild-specific? | basis |
| --- | --- | --- |
| `condition_decreased_guild_signal.png` | yes | v2 trait-enrichment score by condition |
| `guild_trait_direction_shift.png` | yes | the two selected v2 traits |
| `trait_enriched_profile_histogram.png` | yes | trait counts inside profile-level v2 trait-enrichment filtering |
| `condition_enriched_profile_histogram.png` | yes | condition counts after profile-level v2 trait-enrichment filtering |
| `anatomical_sites_enriched_profile_diagram.png` | yes | anatomical body-site counts after profile-level v2 trait-enrichment filtering |
| `anatomical_sites_increased_profile_diagram.png` | yes | increased-profile v2 trait-enrichment score by body-site group |
| `anatomical_sites_diagram.png` | yes | v2 trait-enrichment score by body-site group |
| `candidate_taxa_depletion_bias.png` | no | fixed target-taxon direction counts |
| `condition_target_taxa_stacked_bar.png` | no | fixed target-taxon occurrence counts stacked by taxon |
| `target_taxa_condition_breakdown_bar.png` | no | fixed target-taxon occurrence counts |

## Main Signal

The selected anaerobe-only oxygen-tolerance traits are consistently enriched in decreased profiles relative to increased profiles. This is the decreased-profile counterpart to the aerobe v2 comparison, not a fermentation-defined gut trait set.

Guardrail: the increased-profile anaerobe score is retained in every condition and body-site summary. It checks whether the anaerobe signal is decreased-profile-associated rather than simply high in both directions.

| trait | decreased - increased z |
| --- | ---: |
| obligately anaerobic | 0.100 |
| anaerobic | 0.080 |

## Strongest Decreased-Profile-Biased Taxa

| taxon | decreased / increased profiles | decreased-profile bias |
| --- | ---: | ---: |
| Agathobacter rectalis | 63 / 17 | 0.575 |
| Roseburia intestinalis | 42 / 18 | 0.400 |
| Lachnospira eligens | 60 / 26 | 0.395 |
| Roseburia inulinivorans | 40 / 18 | 0.379 |
| Faecalibacterium prausnitzii | 99 / 59 | 0.253 |
| Ruminococcus bromii | 40 / 26 | 0.212 |
| Anaerobutyricum hallii | 49 / 44 | 0.054 |

## Recurrent Target-Taxon Decreased Profiles By Condition

| taxon | condition | body site | decreased / increased profiles | decreased studies |
| --- | --- | --- | ---: | ---: |
| Faecalibacterium prausnitzii | Parkinson's disease | Feces | 12 / 0 | 7 |
| Roseburia intestinalis | Parkinson's disease | Feces | 12 / 0 | 7 |
| Faecalibacterium prausnitzii | Colorectal cancer | Feces | 5 / 0 | 5 |
| Lachnospira eligens | Colorectal cancer | Feces | 8 / 0 | 4 |
| Faecalibacterium prausnitzii | COVID-19 | Feces | 7 / 0 | 4 |
| Anaerobutyricum hallii | Colorectal cancer | Feces | 4 / 0 | 4 |
| Agathobacter rectalis | COVID-19 | Feces | 5 / 0 | 3 |
| Ruminococcus bromii | COVID-19 | Feces | 5 / 0 | 3 |
| Agathobacter rectalis | Parkinson's disease | Feces | 3 / 0 | 3 |
| Roseburia intestinalis | Colorectal cancer | Feces | 3 / 0 | 3 |

Target-taxon occurrence views (`target_taxa_condition_breakdown_bar.png`, `condition_target_taxa_stacked_bar.png`, `candidate_taxa_profile_occurrences.tsv` and `candidate_taxa_condition_direction.tsv`) count human-host 16S and shotgun profiles only, the corpus of the profile metadata. 97 occurrences (60 decreased, 37 increased) in 68 profiles from 25 studies lie outside the corpus: 60 come from experiments whose host is not recorded as Homo sapiens (Mus musculus 40, Sus scrofa domesticus 7, Hoolock tianxing 3, Rattus norvegicus 3, Theropithecus gelada 2, Bos taurus 1, Felis catus 1, Gallus gallus 1, Nomascus hainanus 1, Ovis aries 1), and 37 from Homo sapiens experiments with a sequencing method other than 16S or WMS, or none recorded (PCR 26, not recorded 11). They have no profile metadata, so they are not counted under any condition; `candidate_taxa_excluded_host_occurrences.tsv` lists them.

## Conditions With Anaerobe v2 Trait Enrichment in Decreased Profiles

| condition | profiles | studies | decreased - increased z | paired d | nominal paired p | matched pairs | decreased mean trait z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Respiratory tract infectious disease | 10 | 2 | 1.327 |  |  | 4 | 0.707 |
| Graft versus host disease | 14 | 5 | 1.290 |  |  | 4 | 0.907 |
| Diarrhea | 24 | 3 | 1.264 | 1.107 | 0.00277 | 12 | 0.695 |
| Age-related macular degeneration | 22 | 5 | 1.249 | 1.835 | 0.000118 | 11 | 0.770 |
| Autism | 14 | 4 | 1.171 |  |  | 4 | 0.656 |
| Response to antibiotic | 9 | 3 | 1.154 |  |  | 3 | 0.602 |
| Colorectal carcinoma | 13 | 3 | 1.036 | 1.246 | 0.0283 | 6 | 0.450 |
| Response to antiviral drug | 32 | 6 | 0.961 | 1.222 | 0.000525 | 14 | 0.718 |

Reported paired p values are uncorrected two-sided paired t-test nominal p-values across 99 tested conditions after excluding broad study-design labels in `condition_guild_summary.tsv` (shown above: top 8 eligible conditions by guild-z difference, 4 of them tested). Paired d and nominal paired p are blank for conditions with fewer than 5 matched decreased/increased profile pairs, or with zero variance among the paired differences.

## Condition Enriched-Profile Histogram

The enrichment-filtered histogram counts non-excluded conditions (`is_excluded_condition_label`) by how many profiles remain after dropping profiles with zero positive anaerobe v2 selected-trait z-scores. The dropped profiles are still retained in the main trait-enrichment summaries and are reported here by profile direction across condition-filtered profiles.

| direction | conditions with retained profiles | retained profiles | zero-enriched profiles | retained / total profiles |
| --- | ---: | ---: | ---: | ---: |
| decreased profiles | 347 | 1631 | 752 | 0.684 |
| increased profiles | 342 | 1645 | 982 | 0.626 |

## Trait Enriched-Profile Histogram

The trait-side companion histogram counts selected heatmap traits by how many retained enriched profiles have a positive z-score for each trait across all human-host profiles with at least one positive anaerobe v2 selected-trait z-score (before `is_excluded_condition_label` condition-label filtering used in the condition histogram above).

| direction | trait | positive retained profiles | retained profile fraction |
| --- | --- | ---: | ---: |
| increased profiles | anaerobic | 1694 | 0.902 |
| decreased profiles | anaerobic | 1684 | 0.915 |
| decreased profiles | bacillus shaped | 1420 | 0.772 |
| decreased profiles | GC high | 1376 | 0.748 |
| increased profiles | rod shaped | 1349 | 0.718 |
| decreased profiles | GC low | 1348 | 0.733 |
| increased profiles | non motile | 1345 | 0.716 |
| increased profiles | GC high | 1343 | 0.715 |
| decreased profiles | rod shaped | 1339 | 0.728 |
| increased profiles | gram positive | 1338 | 0.712 |

## Enrichment-Filtered Anatomical Sites

The anatomical companion keeps body-site groups with at least 10 profiles containing a positive anaerobe v2 selected-trait z-score. Counts below are retained profiles after zero-enriched profiles drop out.

| body_site_group | retained profiles | increased retained | decreased retained | zero-enriched profiles |
| --- | ---: | ---: | ---: | ---: |
| Feces | 2534 | 1218 | 1316 | 959 |
| Saliva | 202 | 109 | 93 | 168 |
| Vagina | 91 | 53 | 38 | 42 |
| Uterine cervix | 52 | 29 | 23 | 29 |
| Subgingival dental plaque | 47 | 27 | 20 | 28 |
| Mouth | 46 | 26 | 20 | 30 |
| Oral cavity | 40 | 28 | 12 | 30 |
| Throat | 33 | 13 | 20 | 31 |

## Data-Inferred Guild Evidence

For the top condition signals, this table uses all taxa observed in the decreased profiles for that condition that also belong to the KG-Microbe-wide anaerobe candidate taxon universe, then reports the profile v2 traits that are enriched in those same decreased profiles and the universe taxa that carry those trait annotations. This is a data-derived guild summary, not a fixed focal-taxon list.

| condition | decreased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Age-related macular degeneration | 11 | 32 / 34 | obligately anaerobic (+1.58); anaerobic (+0.92) | Prevotella (2 traits, 4 profiles); Lachnospiraceae (2 traits, 3 profiles); Clostridia (2 traits, 2 profiles); Bacteroides (2 traits, 2 profiles); Bacteroides eggerthii (2 traits, 2 profiles); Blautia (2 traits, 2 profiles); Bacillota (2 traits, 1 profile); Bacteroidales (2 traits, 1 profile); Prevotellaceae (2 traits, 1 profile); Agathobacter (2 traits, 1 profile); Methanobrevibacter (2 traits, 1 profile); Lachnospira (2 traits, 1 profile) |
| Autism | 6 | 37 / 37 | anaerobic (+1.23); obligately anaerobic (+1.12) | Bacillota (2 traits, 2 profiles); Streptococcus (2 traits, 2 profiles); Ruminococcus (2 traits, 1 profile); Clostridium (2 traits, 1 profile); Peptoniphilus (2 traits, 1 profile); Actinomyces (2 traits, 1 profile); Bifidobacterium (2 traits, 1 profile); Bacteroidales (2 traits, 1 profile); Lachnospiraceae (2 traits, 1 profile); Actinomycetota (2 traits, 1 profile); Gemmiger (2 traits, 1 profile); Ethanoligenens (2 traits, 1 profile) |
| Colorectal carcinoma | 7 | 48 / 56 | anaerobic (+1.11); obligately anaerobic (+0.96) | Ruminococcus (2 traits, 4 profiles); Bacillota (2 traits, 3 profiles); Lachnospira eligens (2 traits, 3 profiles); Eubacterium (2 traits, 3 profiles); Oscillospiraceae (2 traits, 2 profiles); Bacteroides (2 traits, 2 profiles); Bacteroidales (2 traits, 1 profile); Agathobacter (2 traits, 1 profile); Clostridia (2 traits, 1 profile); Eubacteriales (2 traits, 1 profile); Campylobacter (2 traits, 1 profile); Bacteroidia (2 traits, 1 profile) |
| Diarrhea | 12 | 45 / 46 | obligately anaerobic (+1.53); anaerobic (+0.99) | Bacteroides (2 traits, 8 profiles); Eubacterium (2 traits, 7 profiles); Bifidobacterium (2 traits, 6 profiles); Ruminococcus (2 traits, 6 profiles); Gemmiger (2 traits, 5 profiles); Blautia (2 traits, 4 profiles); Coriobacteriaceae (2 traits, 4 profiles); Gemmiger formicilis (2 traits, 3 profiles); Bacteroides fragilis (2 traits, 3 profiles); Streptococcus (2 traits, 2 profiles); Lachnospiraceae (2 traits, 2 profiles); Erysipelotrichaceae (2 traits, 1 profile) |
| Graft versus host disease | 7 | 19 / 22 | obligately anaerobic (+1.33); anaerobic (+1.25) | Lachnospiraceae (2 traits, 5 profiles); Blautia (2 traits, 4 profiles); Clostridia (2 traits, 2 profiles); Eubacteriales (2 traits, 2 profiles); Flavonifractor (2 traits, 2 profiles); Erysipelotrichaceae (2 traits, 1 profile); Peptostreptococcaceae (2 traits, 1 profile); Eubacteriaceae (2 traits, 1 profile); Bacteroidia (2 traits, 1 profile); Alphaproteobacteria (2 traits, 1 profile); Lachnoclostridium (1 trait, 3 profiles); Sellimonas (1 trait, 2 profiles) |
| Respiratory tract infectious disease | 4 | 18 / 20 | anaerobic (+1.33); obligately anaerobic (+1.32) | Erysipelotrichaceae (2 traits, 1 profile); Eubacterium sp. (2 traits, 1 profile); Eubacterium (2 traits, 1 profile); Clostridia (2 traits, 1 profile); Eubacteriales (2 traits, 1 profile); Desulfobacterales (2 traits, 1 profile); Desulfobulbaceae (2 traits, 1 profile); Oscillospiraceae (2 traits, 1 profile); Dethiosulfovibrio (2 traits, 1 profile); Synergistota (2 traits, 1 profile); Erysipelotrichales (2 traits, 1 profile); Verrucomicrobiota (2 traits, 1 profile) |
| Response to antibiotic | 4 | 19 / 23 | anaerobic (+1.36); obligately anaerobic (+0.95) | Bifidobacterium (2 traits, 2 profiles); Bifidobacterium adolescentis (2 traits, 2 profiles); Bifidobacterium breve (2 traits, 2 profiles); Bacillati (2 traits, 1 profile); Oscillospiraceae (2 traits, 1 profile); Clostridiaceae (2 traits, 1 profile); Parabacteroides (2 traits, 1 profile); Butyricicoccus (2 traits, 1 profile); Bacteroides (2 traits, 1 profile); Phocaeicola (2 traits, 1 profile); Bifidobacterium pseudocatenulatum (1 trait, 3 profiles); Bifidobacterium angulatum (1 trait, 2 profiles) |
| Response to antiviral drug | 18 | 33 / 37 | obligately anaerobic (+1.21); anaerobic (+0.71) | Clostridium (2 traits, 5 profiles); Rikenellaceae (2 traits, 3 profiles); Blautia (2 traits, 3 profiles); Lachnospira (2 traits, 3 profiles); Bacteroides (2 traits, 3 profiles); Prevotella (2 traits, 2 profiles); Bifidobacterium (2 traits, 1 profile); Bacteroidales (2 traits, 1 profile); Clostridia (2 traits, 1 profile); Eubacteriales (2 traits, 1 profile); Veillonellaceae (2 traits, 1 profile); Parabacteroides (2 traits, 1 profile) |

## Body-Site Check

| body_site_group | profiles | increased mean trait z | decreased mean trait z | decreased - increased z |
| --- | ---: | ---: | ---: | ---: |
| Feces | 3493 | 0.076 | 0.267 | 0.191 |
| Throat | 64 | -0.336 | -0.078 | 0.258 |
| Vagina | 133 | 0.028 | -0.167 | -0.194 |
| Oral cavity | 70 | -0.100 | -0.231 | -0.131 |
| Uterine cervix | 81 | -0.040 | -0.263 | -0.222 |
| Saliva | 370 | -0.259 | -0.270 | -0.011 |
| Mouth | 76 | -0.192 | -0.283 | -0.091 |
| Other | 1267 | -0.285 | -0.385 | -0.100 |
| Subgingival dental plaque | 75 | 0.089 | -0.388 | -0.477 |

## Top Decreased Profiles To Inspect

| profile | condition | body site | trait z | title |
| --- | --- | --- | ---: | --- |
| bsdb:30669509\|bsdb:30669509-2\|dec | Constipation | Feces | 2.605 | The Effect of Psyllium Husk on Intestinal Microbiota in Constipated Patients and Health... |
| bsdb:28222761\|bsdb:28222761-4\|dec | Autism | Feces | 2.203 | New evidences on the altered gut microbiota in autism spectrum disorders |
| bsdb:36621696\|bsdb:36621696-2\|dec | Psoriasis vulgaris | Feces | 1.847 | Altered gut microbiome composition in nontreated plaque psoriasis patients |
| bsdb:39976263\|bsdb:39976263-10\|dec | Obesity | Feces | 1.847 | <i>Roseburia hominis</i> improves host metabolism in diet-induced obesity |
| bsdb:31171880\|bsdb:31171880-2\|dec | Colorectal cancer | Feces | 1.685 | Metagenomic and metabolomic analyses reveal distinct stage-specific phenotypes of the g... |
| bsdb:41011244\|bsdb:41011244-4\|dec | Therapeutic Procedure | Feces | 1.626 | Metagenomic Characterization of Gut Microbiota in Children with Autism Spectrum Disorde... |
| bsdb:38384431\|bsdb:38384431-8\|dec | Diet measurement | Feces | 1.626 | Short-term pectin-enriched smoothie consumption has beneficial effects on the gut micro... |
| bsdb:38384431\|bsdb:38384431-10\|dec | Diet measurement | Feces | 1.626 | Short-term pectin-enriched smoothie consumption has beneficial effects on the gut micro... |
| bsdb:26416813\|bsdb:26416813-1\|dec | Diet | Feces | 1.626 | High-level adherence to a Mediterranean diet beneficially impacts the gut microbiota an... |
| bsdb:35847083\|bsdb:35847083-7\|dec | Acne | Feces | 1.626 | Alterations of Gut Microbiome and Fecal Fatty Acids in Patients With Polycystic Ovary S... |

## Suggested Follow-Up

1. Treat target-taxon recurrence as an audit view and the selected oxygen-tolerance traits as the formal trait-enrichment score.
2. Prioritize conditions where decreased profiles have high trait-enrichment scores and target taxa are repeatedly reported in decreased profiles across multiple studies.
3. Compare this anaerobe-only v2 signal against the original five-trait set to separate oxygen-tolerance signal from metabolic-capability signal.

## Files

- `guild_trait_direction_shift.tsv`
- `guild_trait_selection.tsv`
- `body_site_guild_summary.tsv`
- `body_site_enriched_profile_summary.tsv`
- `condition_guild_summary.tsv`
- `trait_enriched_profile_counts.tsv`
- `condition_enriched_profile_counts.tsv`
- `candidate_taxa_direction_summary.tsv`
- `candidate_taxa_trait_presence.tsv`
- `candidate_taxa_condition_direction.tsv`
- `candidate_taxa_profile_occurrences.tsv`
- `candidate_taxa_excluded_host_occurrences.tsv`
- `kg_microbe_candidate_taxon_universe.tsv`
- `condition_inferred_guild_summary.tsv`
- `condition_inferred_guild_trait_evidence.tsv`
- `condition_inferred_guild_taxa_evidence.tsv`
- `condition_body_site_guild_map_edges.tsv`
- `condition_body_site_guild_map_sampled_body_sites.tsv`
- `condition_body_site_guild_map_v2_edges.tsv`
- `condition_body_site_guild_map_v2_sampled_body_sites.tsv`
- `condition_body_site_guild_map_v3_individual_edges.tsv`
- `condition_body_site_guild_map_v3_grouped_edges.tsv`
- `condition_body_site_guild_map_v3_grouped_condition_rollups.tsv`
- `condition_body_site_guild_map_v3_grouped_body_site_rollups.tsv`
- `top_decreased_guild_profiles.tsv`
- Methods are consolidated in `../../figure_methods.md`.
- `candidate_taxa_depletion_bias.png`
- `condition_decreased_guild_signal.png`
- `condition_body_site_guild_map.png`
- `condition_body_site_guild_map_v2.png`
- `condition_body_site_guild_map_v3_individual.png`
- `condition_body_site_guild_map_v3_grouped.png`
- `trait_enriched_profile_histogram.png`
- `condition_enriched_profile_histogram.png`
- `condition_target_taxa_stacked_bar.png`
- `guild_trait_direction_shift.png`
- `anatomical_sites_diagram.png`
- `anatomical_sites_increased_profile_diagram.png`
- `anatomical_sites_enriched_profile_diagram.png`
- `target_taxa_condition_breakdown_bar.png`
