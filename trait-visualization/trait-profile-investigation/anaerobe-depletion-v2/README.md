# Anaerobe v2 Depletion Investigation

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

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
| anaerobic | 0.174 |
| obligately anaerobic | 0.099 |

## Strongest Decreased-Profile-Biased Taxa

| taxon | decreased / increased profiles | decreased-profile bias |
| --- | ---: | ---: |
| Agathobacter rectalis | 21 / 2 | 0.826 |
| Roseburia inulinivorans | 22 / 6 | 0.571 |
| Roseburia intestinalis | 22 / 6 | 0.571 |
| Faecalibacterium prausnitzii | 41 / 16 | 0.439 |
| Lachnospira eligens | 30 / 13 | 0.395 |
| Anaerobutyricum hallii | 22 / 13 | 0.257 |
| Ruminococcus bromii | 12 / 11 | 0.043 |

## Recurrent Target-Taxon Decreased Profiles By Condition

| taxon | condition | body site | decreased / increased profiles | decreased studies |
| --- | --- | --- | ---: | ---: |
| Roseburia intestinalis | Parkinson's disease | Feces | 11 / 0 | 6 |
| Faecalibacterium prausnitzii | Parkinson's disease | Feces | 9 / 0 | 5 |
| Anaerobutyricum hallii | Colorectal cancer | Feces | 4 / 0 | 4 |
| Faecalibacterium prausnitzii | Colorectal cancer | Feces | 3 / 0 | 3 |
| Lachnospira eligens | Colorectal cancer | Feces | 3 / 0 | 3 |
| Roseburia intestinalis | Colorectal cancer | Feces | 3 / 0 | 3 |
| Anaerobutyricum hallii | Colorectal carcinoma | Feces | 3 / 0 | 2 |
| Faecalibacterium prausnitzii | Polycystic ovary syndrome | Feces | 3 / 0 | 2 |
| Lachnospira eligens | Colorectal carcinoma | Feces | 3 / 0 | 2 |
| Roseburia inulinivorans | Parkinson's disease | Feces | 7 / 1 | 4 |

Target-taxon occurrence views (`target_taxa_condition_breakdown_bar.png`, `condition_target_taxa_stacked_bar.png`, `candidate_taxa_profile_occurrences.tsv` and `candidate_taxa_condition_direction.tsv`) count human-host 16S and shotgun profiles only, the corpus of the profile metadata. 97 occurrences (60 decreased, 37 increased) in 68 profiles from 25 studies lie outside the corpus: 60 come from experiments whose host is not recorded as Homo sapiens (Mus musculus 40, Sus scrofa domesticus 7, Hoolock tianxing 3, Rattus norvegicus 3, Theropithecus gelada 2, Bos taurus 1, Felis catus 1, Gallus gallus 1, Nomascus hainanus 1, Ovis aries 1), and 37 from Homo sapiens experiments with a sequencing method other than 16S or WMS, or none recorded (PCR 26, not recorded 11). They have no profile metadata, so they are not counted under any condition; `candidate_taxa_excluded_host_occurrences.tsv` lists them.

## Conditions With Anaerobe v2 Trait Enrichment in Decreased Profiles

| condition | profiles | studies | decreased - increased z | paired d | paired p | matched pairs | decreased mean trait z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Diarrhea | 18 | 3 | 1.652 | 1.879 | 0.000489 | 9 | 0.836 |
| Colorectal carcinoma | 6 | 2 | 1.357 |  |  | 3 | 0.462 |
| Depressive disorder | 6 | 2 | 1.313 |  |  | 3 | 1.107 |
| Endometrial cancer | 6 | 2 | 1.219 |  |  | 3 | 0.334 |
| Graft versus host disease | 6 | 2 | 1.188 |  |  | 3 | 0.797 |
| Premature birth | 8 | 2 | 1.138 |  |  | 4 | 0.240 |
| Autism | 8 | 4 | 0.869 |  |  | 4 | 0.511 |
| Gastric adenocarcinoma | 10 | 2 | 0.821 | 0.515 | 0.313 | 5 | 0.003 |

Paired d and paired p are blank for conditions with fewer than 5 matched decreased/increased profile pairs, or with zero variance among the paired differences.

## Condition Enriched-Profile Histogram

The enrichment-filtered histogram counts non-excluded conditions (`is_excluded_condition_label`) by how many profiles remain after dropping profiles with zero positive anaerobe v2 selected-trait z-scores. The dropped profiles are still retained in the main trait-enrichment summaries and are reported here by profile direction across condition-filtered profiles.

| direction | conditions with retained profiles | retained profiles | zero-enriched profiles | retained / total profiles |
| --- | ---: | ---: | ---: | ---: |
| decreased profiles | 171 | 639 | 293 | 0.686 |
| increased profiles | 154 | 554 | 378 | 0.594 |

## Trait Enriched-Profile Histogram

The trait-side companion histogram counts selected heatmap traits by how many retained enriched profiles have a positive z-score for each trait across all human-host profiles with at least one positive anaerobe v2 selected-trait z-score (before `is_excluded_condition_label` condition-label filtering used in the condition histogram above).

| direction | trait | positive retained profiles | retained profile fraction |
| --- | --- | ---: | ---: |
| decreased profiles | anaerobic | 590 | 0.922 |
| decreased profiles | GC high | 503 | 0.786 |
| increased profiles | anaerobic | 495 | 0.889 |
| decreased profiles | GC low | 495 | 0.773 |
| decreased profiles | rod shaped | 489 | 0.764 |
| decreased profiles | gram positive | 468 | 0.731 |
| decreased profiles | temperature optimum mid4 | 468 | 0.731 |
| decreased profiles | obligately anaerobic | 467 | 0.730 |
| decreased profiles | Fermentation | 459 | 0.717 |
| decreased profiles | pH optimum mid1 | 455 | 0.711 |

## Enrichment-Filtered Anatomical Sites

The anatomical companion keeps body-site groups with at least 10 profiles containing a positive anaerobe v2 selected-trait z-score. Counts below are retained profiles after zero-enriched profiles drop out.

| body_site_group | retained profiles | increased retained | decreased retained | zero-enriched profiles |
| --- | ---: | ---: | ---: | ---: |
| Feces | 835 | 358 | 477 | 355 |
| Saliva | 74 | 39 | 35 | 62 |
| Vagina | 25 | 12 | 13 | 13 |
| Subgingival dental plaque | 23 | 14 | 9 | 17 |
| Uterine cervix | 20 | 11 | 9 | 12 |
| Oral cavity | 14 | 7 | 7 | 10 |
| Throat | 12 | 7 | 5 | 12 |

## Data-Inferred Guild Evidence

For the top condition signals, this table uses all taxa observed in the decreased profiles for that condition that also belong to the KG-Microbe-wide anaerobe candidate taxon universe, then reports the profile v2 traits that are enriched in those same decreased profiles and the universe taxa that carry those trait annotations. This is a data-derived guild summary, not a fixed focal-taxon list.

| condition | decreased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Autism | 4 | 35 / 35 | anaerobic (+1.17); obligately anaerobic (+0.57) | Bacillota (2 traits, 2 profiles); Streptococcus (2 traits, 2 profiles); Ruminococcus (2 traits, 1 profile); Clostridium (2 traits, 1 profile); Peptoniphilus (2 traits, 1 profile); Actinomyces (2 traits, 1 profile); Bifidobacterium (2 traits, 1 profile); Bacteroidales (2 traits, 1 profile); Lachnospiraceae (2 traits, 1 profile); Actinomycetota (2 traits, 1 profile); Ethanoligenens (2 traits, 1 profile); Veillonellaceae (2 traits, 1 profile) |
| Colorectal carcinoma | 3 | 39 / 47 | anaerobic (+1.48); obligately anaerobic (+1.23) | Ruminococcus (2 traits, 3 profiles); Lachnospira eligens (2 traits, 3 profiles); Bacillota (2 traits, 2 profiles); Eubacterium (2 traits, 2 profiles); Agathobacter (2 traits, 1 profile); Campylobacter (2 traits, 1 profile); Oscillospiraceae (2 traits, 1 profile); Lachnospira (2 traits, 1 profile); Ruminococcus bromii (2 traits, 1 profile); Butyricicoccus (2 traits, 1 profile); Bacteroides (2 traits, 1 profile); Roseburia (2 traits, 1 profile) |
| Depressive disorder | 3 | 29 / 46 | obligately anaerobic (+1.94); anaerobic (+0.69) | Faecalibacillus intestinalis (2 traits, 2 profiles); Ruminococcus bromii (2 traits, 2 profiles); Gemmiger formicilis (2 traits, 2 profiles); Phocaeicola massiliensis (2 traits, 1 profile); Roseburia hominis (2 traits, 1 profile); Lachnospira eligens (2 traits, 1 profile); Synergistota (2 traits, 1 profile); Bacteroides uniformis (2 traits, 1 profile); Phocaeicola vulgatus (2 traits, 1 profile); Agathobaculum butyriciproducens (1 trait, 2 profiles); Blautia intestinalis (1 trait, 2 profiles); Eubacterium coprostanoligenes (1 trait, 2 profiles) |
| Diarrhea | 9 | 43 / 44 | obligately anaerobic (+1.87); anaerobic (+1.44) | Bacteroides (2 traits, 6 profiles); Bifidobacterium (2 traits, 5 profiles); Ruminococcus (2 traits, 5 profiles); Eubacterium (2 traits, 5 profiles); Gemmiger (2 traits, 4 profiles); Blautia (2 traits, 4 profiles); Gemmiger formicilis (2 traits, 3 profiles); Bacteroides fragilis (2 traits, 3 profiles); Coriobacteriaceae (2 traits, 3 profiles); Lachnospiraceae (2 traits, 2 profiles); Streptococcus (2 traits, 1 profile); Clostridium (2 traits, 1 profile) |
| Endometrial cancer | 3 | 12 / 12 | obligately anaerobic (+1.40); anaerobic (+1.04) | Bacillota (2 traits, 3 profiles); Clostridia (2 traits, 2 profiles); Eubacteriales (2 traits, 2 profiles); Oscillospiraceae (2 traits, 2 profiles); Bacteroidales (2 traits, 1 profile); Bacteroidia (2 traits, 1 profile); Gemmiger (2 traits, 1 profile); Gemmiger formicilis (2 traits, 1 profile); Roseburia (2 traits, 1 profile); Bacteroidota (2 traits, 1 profile); Faecalibacterium (1 trait, 2 profiles); Faecalibacterium prausnitzii (1 trait, 1 profile) |
| Gastric adenocarcinoma | 5 | 14 / 15 | anaerobic (+0.94) | Lachnoclostridium (1 trait, 2 profiles); Clostridia (1 trait, 2 profiles); Eubacteriales (1 trait, 2 profiles); Lachnospiraceae (1 trait, 2 profiles); Lachnospira (1 trait, 2 profiles); Agathobacter rectalis (1 trait, 2 profiles); Roseburia (1 trait, 2 profiles); Bacillota (1 trait, 1 profile); Faecalibacterium (1 trait, 1 profile); Betaproteobacteria (1 trait, 1 profile); Bifidobacteriaceae (1 trait, 1 profile); Oscillibacter (1 trait, 1 profile) |
| Graft versus host disease | 3 | 16 / 19 | anaerobic (+1.19); obligately anaerobic (+1.18) | Lachnospiraceae (2 traits, 2 profiles); Blautia (2 traits, 2 profiles); Flavonifractor (2 traits, 2 profiles); Clostridia (2 traits, 1 profile); Eubacteriales (2 traits, 1 profile); Bacteroidia (2 traits, 1 profile); Alphaproteobacteria (2 traits, 1 profile); Lachnoclostridium (1 trait, 2 profiles); Sellimonas (1 trait, 2 profiles); Gemella (1 trait, 1 profile); Intestinimonas (1 trait, 1 profile); Ruminiclostridium (1 trait, 1 profile) |
| Premature birth | 4 | 35 / 49 | obligately anaerobic (+1.60) | Bacteroidaceae (1 trait, 3 profiles); Porphyromonadaceae (1 trait, 2 profiles); Clostridiaceae (1 trait, 2 profiles); Pasteurellaceae (1 trait, 2 profiles); Streptococcus (1 trait, 1 profile); Actinomyces (1 trait, 1 profile); Campylobacter (1 trait, 1 profile); Actinomycetota (1 trait, 1 profile); Bifidobacteriaceae (1 trait, 1 profile); Veillonellaceae (1 trait, 1 profile); Campylobacter gracilis (1 trait, 1 profile); Porphyromonas (1 trait, 1 profile) |

## Body-Site Check

| body_site_group | profiles | increased mean trait z | decreased mean trait z | decreased - increased z |
| --- | ---: | ---: | ---: | ---: |
| Feces | 1190 | 0.011 | 0.335 | 0.324 |
| Oral cavity | 24 | -0.291 | -0.159 | 0.133 |
| Vagina | 38 | 0.073 | -0.175 | -0.247 |
| Saliva | 136 | -0.203 | -0.280 | -0.078 |
| Other | 364 | -0.259 | -0.414 | -0.155 |
| Uterine cervix | 32 | -0.008 | -0.447 | -0.439 |
| Throat | 24 | -0.197 | -0.462 | -0.265 |
| Nasopharynx | 22 | -0.570 | -0.696 | -0.126 |
| Subgingival dental plaque | 40 | 0.078 | -0.763 | -0.842 |

## Top Decreased Profiles To Inspect

| profile | condition | body site | trait z | title |
| --- | --- | --- | ---: | --- |
| bsdb:36621696\|bsdb:36621696-2\|dec | Psoriasis vulgaris | Feces | 1.948 | Altered gut microbiome composition in nontreated plaque psoriasis patients |
| bsdb:36621696\|bsdb:36621696-3\|dec | Psoriasis vulgaris | Feces | 1.634 | Altered gut microbiome composition in nontreated plaque psoriasis patients |
| bsdb:36621696\|bsdb:36621696-4\|dec | Psoriasis vulgaris | Feces | 1.634 | Altered gut microbiome composition in nontreated plaque psoriasis patients |
| bsdb:39773949\|bsdb:39773949-1\|dec | Malnutrition | Feces | 1.598 | Interconnected pathways link faecal microbiota plasma lipids and brain activity to chil... |
| bsdb:30669509\|bsdb:30669509-3\|dec | Constipation | Feces | 1.594 | The Effect of Psyllium Husk on Intestinal Microbiota in Constipated Patients and Health... |
| bsdb:39656003\|bsdb:39656003-1\|dec | Psoriasis | Feces | 1.499 | Gut dysbiosis narrative in psoriasis: matched-pair approach identifies only subtle shif... |
| bsdb:32280686\|bsdb:32280686-1\|dec | Colorectal cancer | Feces | 1.463 | Alterations of the Predominant Fecal Microbiota and Disruption of the Gut Mucosal Barri... |
| bsdb:28852182\|bsdb:28852182-2\|dec | Chronic constipation | Feces | 1.413 | Unveiling the gut microbiota composition and functionality associated with constipation... |
| bsdb:24013136\|bsdb:24013136-4\|dec | Diabetes mellitus | Feces | 1.412 | Human gut microbiota changes reveal the progression of glucose intolerance |
| bsdb:24013136\|bsdb:24013136-5\|dec | Prediabetes syndrome | Feces | 1.412 | Human gut microbiota changes reveal the progression of glucose intolerance |

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
