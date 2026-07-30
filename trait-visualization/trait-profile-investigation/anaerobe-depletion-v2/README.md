# Anaerobe v2 Depletion Investigation

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Guild traits: anaerobic (METPO:1000603), obligately anaerobic (METPO:1000607)
- Target taxa: Faecalibacterium prausnitzii, Anaerobutyricum hallii, Lachnospira eligens, Agathobacter rectalis, Roseburia intestinalis, Ruminococcus bromii, Roseburia inulinivorans
- Formal KG-Microbe taxon universe: 1,959 taxa with at least one core oxygen-tolerance annotation among anaerobic, obligately anaerobic
- Guild trait selection audit: `guild_trait_selection.tsv` records the selected trait IDs, short rationale, and KG-Microbe taxon counts for each selected trait.
- Context: all body sites, with body-site and condition summaries to separate gut, oral, airway, urogenital, and other patterns.

## Plot Guild Specificity

| plot | guild-specific? | basis |
| --- | --- | --- |
| `condition_decreased_guild_signal.png` | yes | v2 guild score by condition |
| `guild_trait_direction_shift.png` | yes | the two selected v2 guild traits |
| `trait_enriched_profile_histogram.png` | yes | trait counts inside profile-level v2 trait-enrichment filtering |
| `condition_enriched_profile_histogram.png` | yes | condition counts after profile-level v2 trait-enrichment filtering |
| `anatomical_sites_enriched_profile_diagram.png` | yes | anatomical body-site counts after profile-level v2 trait-enrichment filtering |
| `anatomical_sites_diagram.png` | yes | v2 guild score by body-site group |
| `candidate_taxa_depletion_bias.png` | no | fixed target-taxon direction counts |
| `condition_target_taxa_stacked_bar.png` | no | fixed target-taxon occurrence counts stacked by taxon |
| `target_taxa_condition_breakdown_bar.png` | no | fixed target-taxon occurrence counts |

## Main Signal

The selected anaerobe-only oxygen-tolerance traits are consistently higher in decreased profiles than in increased profiles. This is the decreased-profile counterpart to the facultative-aerobe v2 comparison, not a fermentation-defined gut guild.

Guardrail: the increased-profile anaerobe score is retained in every condition and body-site summary. It checks whether the anaerobe signal is depletion-associated rather than simply high in both directions.

| trait | decreased - increased z |
| --- | ---: |
| anaerobic | 0.077 |
| obligately anaerobic | 0.074 |

## Strongest Depletion-Biased Taxa

| taxon | decreased / increased profiles | depletion bias |
| --- | ---: | ---: |
| Agathobacter rectalis | 92 / 33 | 0.472 |
| Lachnospira eligens | 101 / 39 | 0.443 |
| Roseburia intestinalis | 76 / 31 | 0.421 |
| Roseburia inulinivorans | 62 / 29 | 0.363 |
| Ruminococcus bromii | 61 / 38 | 0.232 |
| Faecalibacterium prausnitzii | 177 / 111 | 0.229 |
| Anaerobutyricum hallii | 86 / 63 | 0.154 |

## Recurrent Target-Taxon Depletions By Condition

| taxon | condition | body site | decreased / increased profiles | decreased studies |
| --- | --- | --- | ---: | ---: |
| Faecalibacterium prausnitzii | Parkinson's disease | Feces | 16 / 0 | 11 |
| Roseburia intestinalis | Parkinson's disease | Feces | 12 / 0 | 7 |
| Faecalibacterium prausnitzii | COVID-19 | Feces | 11 / 0 | 6 |
| Lachnospira eligens | Colorectal cancer | Feces | 9 / 0 | 5 |
| Ruminococcus bromii | COVID-19 | Feces | 9 / 0 | 5 |
| Agathobacter rectalis | COVID-19 | Feces | 7 / 0 | 5 |
| Anaerobutyricum hallii | Colorectal cancer | Feces | 5 / 0 | 5 |
| Roseburia intestinalis | Colorectal cancer | Feces | 5 / 0 | 5 |
| Agathobacter rectalis | Colorectal cancer | Feces | 4 / 0 | 3 |
| Agathobacter rectalis | Parkinson's disease | Feces | 3 / 0 | 3 |

## Conditions With High Decreased-Profile Guild Signal

| condition | profiles | studies | decreased - increased z | paired d | paired p | matched pairs | decreased mean guild z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Age-related macular degeneration | 27 | 5 | 1.269 | 2.101 | 6.53e-06 | 13 | 0.767 |
| Graft versus host disease | 13 | 4 | 1.254 | 1.480 | 0.0595 | 4 | 0.848 |
| Acute graft vs. host disease | 11 | 2 | 1.080 | 0.675 | 0.27 | 4 | 0.569 |
| Lung transplantation | 7 | 2 | 1.049 |  |  | 2 | 0.479 |
| Head and neck carcinoma | 7 | 2 | 0.828 | 0.966 | 0.236 | 3 | 0.132 |
| Acute pancreatitis | 20 | 2 | 0.818 | 0.977 | 0.0129 | 10 | 0.394 |
| Graves ophthalmopathy | 13 | 2 | 0.818 | 0.887 | 0.118 | 5 | 0.023 |
| Diarrhea | 46 | 5 | 0.817 | 0.676 | 0.0046 | 22 | 0.523 |

## Condition Enriched-Profile Histogram

The enrichment-filtered histogram counts conditions by how many profiles remain after dropping profiles with zero positive anaerobe v2 guild-trait z-scores. The dropped profiles are still retained in the main guild-score summaries and are reported here by profile direction.

| direction | conditions with retained profiles | retained profiles | zero-enriched profiles | retained / total profiles |
| --- | ---: | ---: | ---: | ---: |
| Decreased | 430 | 3414 | 1569 | 0.685 |
| Increased | 438 | 3666 | 2146 | 0.631 |

## Trait Enriched-Profile Histogram

The trait-side companion histogram counts selected heatmap traits by how many retained enriched profiles have a positive z-score for each trait. The retained-profile universe is the same one used by the condition histogram: profiles with at least one positive anaerobe v2 guild-trait z-score.

| direction | trait | positive retained profiles | retained profile fraction |
| --- | --- | ---: | ---: |
| Increased | anaerobic | 3824 | 0.927 |
| Decreased | anaerobic | 3579 | 0.936 |
| Increased | non motile | 3126 | 0.758 |
| Increased | GC high | 2999 | 0.727 |
| Increased | GC low | 2986 | 0.724 |
| Increased | bacillus shaped | 2964 | 0.719 |
| Increased | gram positive | 2920 | 0.708 |
| Decreased | GC low | 2872 | 0.751 |
| Decreased | non motile | 2870 | 0.751 |
| Decreased | GC high | 2863 | 0.749 |

## Enrichment-Filtered Anatomical Sites

The anatomical companion keeps body-site groups with at least 10 profiles containing a positive anaerobe v2 guild-trait z-score. Counts below are retained profiles after zero-enriched profiles drop out.

| body_site_group | retained profiles | inc retained | dec retained | zero-enriched profiles |
| --- | ---: | ---: | ---: | ---: |
| Feces | 5459 | 2750 | 2709 | 2136 |
| Saliva | 384 | 202 | 182 | 254 |
| Subgingival dental plaque | 127 | 76 | 51 | 84 |
| Caecum | 124 | 73 | 51 | 41 |
| Vagina | 123 | 77 | 46 | 73 |
| Colon | 105 | 55 | 50 | 34 |
| Mouth | 79 | 39 | 40 | 91 |
| Nasopharynx | 70 | 51 | 19 | 56 |

## Data-Inferred Guild Evidence

For the top condition signals, this table uses all taxa observed in the decreased profiles for that condition that also belong to the KG-Microbe-wide anaerobe candidate taxon universe, then reports the profile v2 guild traits that are enriched in those same decreased profiles and the universe taxa that carry those trait annotations. This is a data-derived guild summary, not a fixed focal-taxon list.

| condition | decreased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Acute graft vs. host disease | 7 | 28 / 30 | obligately anaerobic (+1.09); anaerobic (+1.07) | taxon_841 (2 traits, 4 profiles); taxon_1263 (2 traits, 2 profiles); taxon_171549 (2 traits, 2 profiles); taxon_200643 (2 traits, 2 profiles); taxon_28050 (2 traits, 2 profiles); taxon_816 (2 traits, 2 profiles); taxon_976 (2 traits, 2 profiles); taxon_1224 (2 traits, 1 profile); taxon_1236 (2 traits, 1 profile); taxon_1678 (2 traits, 1 profile); taxon_171550 (2 traits, 1 profile); taxon_186801 (2 traits, 1 profile) |
| Acute pancreatitis | 10 | 49 / 77 | anaerobic (+0.92); obligately anaerobic (+0.72) | taxon_200940 (2 traits, 3 profiles); taxon_208479 (2 traits, 3 profiles); taxon_1496 (2 traits, 2 profiles); taxon_28889 (2 traits, 2 profiles); taxon_292800 (2 traits, 2 profiles); taxon_872 (2 traits, 2 profiles); taxon_1224 (2 traits, 1 profile); taxon_1239 (2 traits, 1 profile); taxon_1350 (2 traits, 1 profile); taxon_201174 (2 traits, 1 profile); taxon_246787 (2 traits, 1 profile); taxon_28111 (2 traits, 1 profile) |
| Age-related macular degeneration | 14 | 36 / 39 | obligately anaerobic (+1.57); anaerobic (+0.97) | taxon_838 (2 traits, 6 profiles); taxon_186803 (2 traits, 3 profiles); taxon_186801 (2 traits, 2 profiles); taxon_816 (2 traits, 2 profiles); taxon_28111 (2 traits, 2 profiles); taxon_572511 (2 traits, 2 profiles); taxon_1239 (2 traits, 1 profile); taxon_171549 (2 traits, 1 profile); taxon_171552 (2 traits, 1 profile); taxon_2172 (2 traits, 1 profile); taxon_28050 (2 traits, 1 profile); taxon_3085636 (2 traits, 1 profile) |
| Diarrhea | 22 | 86 / 92 | obligately anaerobic (+0.92); anaerobic (+0.72) | taxon_816 (2 traits, 9 profiles); taxon_84107 (2 traits, 6 profiles); taxon_1678 (2 traits, 6 profiles); taxon_572511 (2 traits, 6 profiles); taxon_1263 (2 traits, 6 profiles); taxon_817 (2 traits, 4 profiles); taxon_1239 (2 traits, 4 profiles); taxon_186803 (2 traits, 3 profiles); taxon_201174 (2 traits, 3 profiles); taxon_91061 (2 traits, 2 profiles); taxon_1485 (2 traits, 1 profile); taxon_1643824 (2 traits, 1 profile) |
| Graft versus host disease | 6 | 19 / 22 | obligately anaerobic (+1.47); anaerobic (+1.04) | taxon_186803 (2 traits, 5 profiles); taxon_572511 (2 traits, 3 profiles); taxon_186801 (2 traits, 2 profiles); taxon_186802 (2 traits, 2 profiles); taxon_946234 (2 traits, 2 profiles); taxon_186804 (2 traits, 1 profile); taxon_200643 (2 traits, 1 profile); taxon_1506553 (1 trait, 3 profiles); taxon_1769710 (1 trait, 2 profiles); taxon_128827 (1 trait, 1 profile); taxon_1378 (1 trait, 1 profile); taxon_1392389 (1 trait, 1 profile) |
| Graves ophthalmopathy | 5 | 14 / 15 | anaerobic (+1.38) | taxon_29523 (1 trait, 2 profiles); taxon_976 (1 trait, 2 profiles); taxon_102106 (1 trait, 1 profile); taxon_1239 (1 trait, 1 profile); taxon_1263 (1 trait, 1 profile); taxon_1407607 (1 trait, 1 profile); taxon_204516 (1 trait, 1 profile); taxon_207244 (1 trait, 1 profile); taxon_328814 (1 trait, 1 profile); taxon_39488 (1 trait, 1 profile); taxon_39496 (1 trait, 1 profile); taxon_41200 (1 trait, 1 profile) |
| Head and neck carcinoma | 4 | 53 / 62 | anaerobic (+0.88); obligately anaerobic (+0.78) | taxon_836 (2 traits, 3 profiles); taxon_186804 (2 traits, 2 profiles); taxon_203490 (2 traits, 2 profiles); taxon_203491 (2 traits, 2 profiles); taxon_203492 (2 traits, 2 profiles); taxon_32066 (2 traits, 2 profiles); taxon_838 (2 traits, 2 profiles); taxon_841 (2 traits, 2 profiles); taxon_136 (2 traits, 1 profile); taxon_137 (2 traits, 1 profile); taxon_157 (2 traits, 1 profile); taxon_1654 (2 traits, 1 profile) |
| Lung transplantation | 3 | 25 / 25 | obligately anaerobic (+1.06); anaerobic (+1.04) | taxon_171549 (2 traits, 2 profiles); taxon_171552 (2 traits, 2 profiles); taxon_186802 (2 traits, 2 profiles); taxon_200643 (2 traits, 2 profiles); taxon_203490 (2 traits, 2 profiles); taxon_203491 (2 traits, 2 profiles); taxon_203492 (2 traits, 2 profiles); taxon_216572 (2 traits, 2 profiles); taxon_310298 (2 traits, 2 profiles); taxon_32066 (2 traits, 2 profiles); taxon_815 (2 traits, 2 profiles); taxon_816 (2 traits, 2 profiles) |

## Body-Site Check

| body_site_group | profiles | decreased mean guild z | decreased - increased z |
| --- | ---: | ---: | ---: |
| Colon | 139 | 0.445 | 0.366 |
| Feces | 7595 | 0.231 | 0.152 |
| Caecum | 165 | 0.100 | -0.069 |
| Saliva | 638 | -0.153 | -0.006 |
| Vagina | 196 | -0.193 | -0.218 |
| Subgingival dental plaque | 211 | -0.247 | -0.261 |
| Other | 2829 | -0.344 | -0.066 |
| Nasopharynx | 126 | -0.373 | 0.036 |
| Mouth | 170 | -0.519 | 0.011 |
| unknown | 39 | -1.427 | -0.147 |

## Top Decreased Profiles To Inspect

| profile | condition | body site | guild z | title |
| --- | --- | --- | ---: | --- |
| bsdb:36621696\|bsdb:36621696-2\|dec | Psoriasis vulgaris | Feces | 1.921 | Altered gut microbiome composition in nontreated plaque psoriasis patients |
| bsdb:39976263\|bsdb:39976263-10\|dec | Alanine measurement | Feces | 1.921 | <i>Roseburia hominis</i> improves host metabolism in diet-induced obesity |
| bsdb:31171880\|bsdb:31171880-2\|dec | Colorectal cancer | Feces | 1.891 | Metagenomic and metabolomic analyses reveal distinct stage-specific phenotypes of the g... |
| bsdb:34698096\|bsdb:34698096-3\|dec | Sensory perception of bitter taste | Internal cheek pouch | 1.866 | Association of Bitter Taste Receptor T2R38 Polymorphisms, Oral Microbiota, and Rheumato... |
| bsdb:39417117\|bsdb:39417117-24\|dec | Metabolic toxicity | Feces | 1.815 | Toxic Metals Impact Gut Microbiota and Metabolic Risk in Five African-Origin Populations |
| bsdb:40811957\|bsdb:40811957-7\|dec | Peri-Implantitis | Dental plaque | 1.803 | Microbiome Shifts in Peri-Implantitis: Longitudinal Characterization of Dysbiosis and R... |
| bsdb:39191760\|bsdb:39191760-11\|dec | Response to transplant | Feces | 1.743 | Fecal microbiota transplantation influences microbiota without connection to symptom re... |
| bsdb:36566315\|bsdb:36566315-1\|dec | Fasting blood glucose measurement | Feces | 1.721 | The gut microbiome of obese postpartum women with and without previous gestational diab... |
| bsdb:40520642\|bsdb:40520642-1\|dec | Type II diabetes mellitus | Feces | 1.721 | Fecal microbiome analysis in patients with metabolic syndrome and type 2 diabetes |
| bsdb:35140064\|bsdb:35140064-15\|dec | Response to vaccine | Feces | 1.721 | Gut microbiota composition is associated with SARS-CoV-2 vaccine immunogenicity and adv... |

## Suggested Follow-Up

1. Treat target-taxon recurrence as an audit view and the selected oxygen-tolerance traits as the formal guild score.
2. Prioritize conditions where decreased profiles have high guild scores and target taxa are repeatedly decreased across multiple studies.
3. Compare this anaerobe-only v2 signal against the original five-trait guild to separate oxygen-tolerance signal from metabolic-capability signal.

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
- `kg_microbe_candidate_taxon_universe.tsv`
- `condition_inferred_guild_summary.tsv`
- `condition_inferred_guild_trait_evidence.tsv`
- `condition_inferred_guild_taxa_evidence.tsv`
- `top_decreased_guild_profiles.tsv`
- `figure_methods.md`
- `candidate_taxa_depletion_bias.png`
- `condition_decreased_guild_signal.png`
- `trait_enriched_profile_histogram.png`
- `condition_enriched_profile_histogram.png`
- `condition_target_taxa_stacked_bar.png`
- `guild_trait_direction_shift.png`
- `anatomical_sites_diagram.png`
- `anatomical_sites_enriched_profile_diagram.png`
- `target_taxa_condition_breakdown_bar.png`
