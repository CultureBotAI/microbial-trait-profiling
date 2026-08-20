# Anaerobe v2 Depletion Investigation

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Selected traits: anaerobic (METPO:1000603), obligately anaerobic (METPO:1000607)
- Target taxa: Faecalibacterium prausnitzii, Anaerobutyricum hallii, Lachnospira eligens, Agathobacter rectalis, Roseburia intestinalis, Ruminococcus bromii, Roseburia inulinivorans
- Formal KG-Microbe taxon universe: 1,959 taxa with at least one core oxygen-tolerance annotation among anaerobic, obligately anaerobic
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

The selected anaerobe-only oxygen-tolerance traits are consistently enriched in decreased profiles relative to increased profiles. This is the decreased-profile counterpart to the facultative-aerobe v2 comparison, not a fermentation-defined gut trait set.

Guardrail: the increased-profile anaerobe score is retained in every condition and body-site summary. It checks whether the anaerobe signal is decreased-profile-associated rather than simply high in both directions.

| trait | decreased - increased z |
| --- | ---: |
| obligately anaerobic | 0.084 |
| anaerobic | 0.067 |

## Strongest Decreased-Profile-Biased Taxa

| taxon | decreased / increased profiles | decreased-profile bias |
| --- | ---: | ---: |
| Agathobacter rectalis | 92 / 33 | 0.472 |
| Lachnospira eligens | 101 / 39 | 0.443 |
| Roseburia intestinalis | 76 / 31 | 0.421 |
| Roseburia inulinivorans | 62 / 29 | 0.363 |
| Ruminococcus bromii | 61 / 38 | 0.232 |
| Faecalibacterium prausnitzii | 177 / 111 | 0.229 |
| Anaerobutyricum hallii | 86 / 63 | 0.154 |

## Recurrent Target-Taxon Decreased Profiles By Condition

| taxon | condition | body site | decreased / increased profiles | decreased studies |
| --- | --- | --- | ---: | ---: |
| Faecalibacterium prausnitzii | Parkinson's disease | Feces | 16 / 0 | 11 |
| Roseburia intestinalis | Parkinson's disease | Feces | 12 / 0 | 7 |
| Faecalibacterium prausnitzii | COVID-19 | Feces | 11 / 0 | 6 |
| Faecalibacterium prausnitzii | Colorectal cancer | Feces | 6 / 0 | 6 |
| Lachnospira eligens | Colorectal cancer | Feces | 9 / 0 | 5 |
| Ruminococcus bromii | COVID-19 | Feces | 9 / 0 | 5 |
| Agathobacter rectalis | COVID-19 | Feces | 7 / 0 | 5 |
| Anaerobutyricum hallii | Colorectal cancer | Feces | 5 / 0 | 5 |
| Roseburia intestinalis | Colorectal cancer | Feces | 5 / 0 | 5 |
| Agathobacter rectalis | Colorectal cancer | Feces | 4 / 0 | 3 |

## Conditions With Anaerobe v2 Trait Enrichment in Decreased Profiles

| condition | profiles | studies | decreased - increased z | paired d | paired p | matched pairs | decreased mean trait z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Acute graft vs. host disease | 6 | 2 | 2.216 |  |  | 2 | 0.939 |
| Age-related macular degeneration | 27 | 5 | 1.269 | 2.101 | 6.53e-06 | 13 | 0.769 |
| Graft versus host disease | 13 | 4 | 1.253 | 1.480 | 0.0596 | 4 | 0.850 |
| Diarrhea | 32 | 4 | 1.233 | 1.131 | 0.000629 | 15 | 0.775 |
| Lung transplantation | 7 | 2 | 1.049 |  |  | 2 | 0.481 |
| Depressive disorder | 15 | 3 | 0.876 | 0.980 | 0.0617 | 6 | 0.931 |
| Head and neck carcinoma | 7 | 2 | 0.828 | 0.966 | 0.236 | 3 | 0.134 |
| Graves ophthalmopathy | 13 | 2 | 0.819 | 0.886 | 0.119 | 5 | 0.026 |

## Condition Enriched-Profile Histogram

The enrichment-filtered histogram counts conditions by how many profiles remain after dropping profiles with zero positive anaerobe v2 selected-trait z-scores. The dropped profiles are still retained in the main trait-enrichment summaries and are reported here by profile direction.

| direction | conditions with retained profiles | retained profiles | zero-enriched profiles | retained / total profiles |
| --- | ---: | ---: | ---: | ---: |
| decreased profiles | 373 | 2692 | 1241 | 0.684 |
| increased profiles | 378 | 2755 | 1684 | 0.621 |

## Trait Enriched-Profile Histogram

The trait-side companion histogram counts selected heatmap traits by how many retained enriched profiles have a positive z-score for each trait. The retained-profile universe is the same one used by the condition histogram: profiles with at least one positive anaerobe v2 selected-trait z-score.

| direction | trait | positive retained profiles | retained profile fraction |
| --- | --- | ---: | ---: |
| increased profiles | anaerobic | 2866 | 0.927 |
| decreased profiles | anaerobic | 2805 | 0.936 |
| increased profiles | non motile | 2307 | 0.746 |
| decreased profiles | bacillus shaped | 2245 | 0.749 |
| decreased profiles | GC low | 2243 | 0.748 |
| decreased profiles | GC high | 2238 | 0.747 |
| increased profiles | GC high | 2235 | 0.723 |
| decreased profiles | non motile | 2223 | 0.742 |
| increased profiles | bacillus shaped | 2196 | 0.710 |
| increased profiles | gram positive | 2185 | 0.706 |

## Enrichment-Filtered Anatomical Sites

The anatomical companion keeps body-site groups with at least 10 profiles containing a positive anaerobe v2 selected-trait z-score. Counts below are retained profiles after zero-enriched profiles drop out.

| body_site_group | retained profiles | increased retained | decreased retained | zero-enriched profiles |
| --- | ---: | ---: | ---: | ---: |
| Feces | 4116 | 2012 | 2104 | 1653 |
| Saliva | 372 | 197 | 175 | 247 |
| Subgingival dental plaque | 126 | 75 | 51 | 84 |
| Vagina | 119 | 77 | 42 | 69 |
| Nasopharynx | 66 | 47 | 19 | 49 |
| Oral cavity | 55 | 31 | 24 | 39 |
| Uterine cervix | 49 | 27 | 22 | 34 |
| Sputum | 34 | 16 | 18 | 48 |

## Data-Inferred Guild Evidence

For the top condition signals, this table uses all taxa observed in the decreased profiles for that condition that also belong to the KG-Microbe-wide anaerobe candidate taxon universe, then reports the profile v2 traits that are enriched in those same decreased profiles and the universe taxa that carry those trait annotations. This is a data-derived guild summary, not a fixed focal-taxon list.

| condition | decreased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Acute graft vs. host disease | 4 | 19 / 20 | anaerobic (+2.50); obligately anaerobic (+1.93) | taxon_841 (2 traits, 4 profiles); taxon_28050 (2 traits, 2 profiles); taxon_816 (2 traits, 2 profiles); taxon_1263 (2 traits, 1 profile); taxon_1678 (2 traits, 1 profile); taxon_171549 (2 traits, 1 profile); taxon_171550 (2 traits, 1 profile); taxon_186801 (2 traits, 1 profile); taxon_186802 (2 traits, 1 profile); taxon_186803 (2 traits, 1 profile); taxon_200643 (2 traits, 1 profile); taxon_216572 (2 traits, 1 profile) |
| Age-related macular degeneration | 14 | 36 / 39 | obligately anaerobic (+1.57); anaerobic (+0.97) | taxon_838 (2 traits, 6 profiles); taxon_186803 (2 traits, 3 profiles); taxon_186801 (2 traits, 2 profiles); taxon_816 (2 traits, 2 profiles); taxon_28111 (2 traits, 2 profiles); taxon_572511 (2 traits, 2 profiles); taxon_1239 (2 traits, 1 profile); taxon_171549 (2 traits, 1 profile); taxon_171552 (2 traits, 1 profile); taxon_2172 (2 traits, 1 profile); taxon_28050 (2 traits, 1 profile); taxon_3085636 (2 traits, 1 profile) |
| Depressive disorder | 7 | 48 / 105 | obligately anaerobic (+1.03); anaerobic (+0.72) | taxon_31979 (2 traits, 3 profiles); taxon_39485 (2 traits, 2 profiles); taxon_135625 (2 traits, 2 profiles); taxon_1485 (2 traits, 2 profiles); taxon_186804 (2 traits, 2 profiles); taxon_40518 (2 traits, 2 profiles); taxon_712 (2 traits, 2 profiles); taxon_841 (2 traits, 2 profiles); taxon_1239 (2 traits, 1 profile); taxon_1263 (2 traits, 1 profile); taxon_1502 (2 traits, 1 profile); taxon_186801 (2 traits, 1 profile) |
| Diarrhea | 15 | 46 / 48 | obligately anaerobic (+1.64); anaerobic (+0.83) | taxon_816 (2 traits, 9 profiles); taxon_1678 (2 traits, 6 profiles); taxon_1263 (2 traits, 6 profiles); taxon_84107 (2 traits, 5 profiles); taxon_817 (2 traits, 4 profiles); taxon_572511 (2 traits, 4 profiles); taxon_186803 (2 traits, 2 profiles); taxon_1485 (2 traits, 1 profile); taxon_171552 (2 traits, 1 profile); taxon_186802 (2 traits, 1 profile); taxon_201174 (2 traits, 1 profile); taxon_31953 (2 traits, 1 profile) |
| Graft versus host disease | 6 | 19 / 22 | obligately anaerobic (+1.47); anaerobic (+1.04) | taxon_186803 (2 traits, 5 profiles); taxon_572511 (2 traits, 3 profiles); taxon_186801 (2 traits, 2 profiles); taxon_186802 (2 traits, 2 profiles); taxon_946234 (2 traits, 2 profiles); taxon_186804 (2 traits, 1 profile); taxon_200643 (2 traits, 1 profile); taxon_1506553 (1 trait, 3 profiles); taxon_1769710 (1 trait, 2 profiles); taxon_128827 (1 trait, 1 profile); taxon_1378 (1 trait, 1 profile); taxon_1392389 (1 trait, 1 profile) |
| Graves ophthalmopathy | 5 | 14 / 15 | anaerobic (+1.38) | taxon_29523 (1 trait, 2 profiles); taxon_976 (1 trait, 2 profiles); taxon_102106 (1 trait, 1 profile); taxon_1239 (1 trait, 1 profile); taxon_1263 (1 trait, 1 profile); taxon_1407607 (1 trait, 1 profile); taxon_204516 (1 trait, 1 profile); taxon_207244 (1 trait, 1 profile); taxon_328814 (1 trait, 1 profile); taxon_39488 (1 trait, 1 profile); taxon_39496 (1 trait, 1 profile); taxon_41200 (1 trait, 1 profile) |
| Head and neck carcinoma | 4 | 53 / 62 | anaerobic (+0.88); obligately anaerobic (+0.78) | taxon_836 (2 traits, 3 profiles); taxon_186804 (2 traits, 2 profiles); taxon_203490 (2 traits, 2 profiles); taxon_203491 (2 traits, 2 profiles); taxon_203492 (2 traits, 2 profiles); taxon_32066 (2 traits, 2 profiles); taxon_838 (2 traits, 2 profiles); taxon_841 (2 traits, 2 profiles); taxon_136 (2 traits, 1 profile); taxon_137 (2 traits, 1 profile); taxon_157 (2 traits, 1 profile); taxon_1654 (2 traits, 1 profile) |
| Lung transplantation | 3 | 25 / 25 | obligately anaerobic (+1.05); anaerobic (+1.04) | taxon_171549 (2 traits, 2 profiles); taxon_171552 (2 traits, 2 profiles); taxon_186802 (2 traits, 2 profiles); taxon_200643 (2 traits, 2 profiles); taxon_203490 (2 traits, 2 profiles); taxon_203491 (2 traits, 2 profiles); taxon_203492 (2 traits, 2 profiles); taxon_216572 (2 traits, 2 profiles); taxon_310298 (2 traits, 2 profiles); taxon_32066 (2 traits, 2 profiles); taxon_815 (2 traits, 2 profiles); taxon_816 (2 traits, 2 profiles) |

## Body-Site Check

| body_site_group | profiles | increased mean trait z | decreased mean trait z | decreased - increased z |
| --- | ---: | ---: | ---: | ---: |
| Feces | 5769 | 0.079 | 0.239 | 0.160 |
| Saliva | 619 | -0.129 | -0.159 | -0.029 |
| Vagina | 188 | 0.043 | -0.190 | -0.233 |
| Subgingival dental plaque | 210 | 0.016 | -0.244 | -0.261 |
| Uterine cervix | 83 | -0.023 | -0.272 | -0.249 |
| Oral cavity | 94 | -0.175 | -0.274 | -0.099 |
| Nasopharynx | 115 | -0.331 | -0.336 | -0.004 |
| Other | 2184 | -0.272 | -0.336 | -0.064 |
| Sputum | 82 | -0.598 | -0.391 | 0.208 |

## Top Decreased Profiles To Inspect

| profile | condition | body site | trait z | title |
| --- | --- | --- | ---: | --- |
| bsdb:36621696\|bsdb:36621696-2\|dec | Psoriasis vulgaris | Feces | 1.922 | Altered gut microbiome composition in nontreated plaque psoriasis patients |
| bsdb:39976263\|bsdb:39976263-10\|dec | Alanine measurement | Feces | 1.922 | <i>Roseburia hominis</i> improves host metabolism in diet-induced obesity |
| bsdb:31171880\|bsdb:31171880-2\|dec | Colorectal cancer | Feces | 1.892 | Metagenomic and metabolomic analyses reveal distinct stage-specific phenotypes of the g... |
| bsdb:34698096\|bsdb:34698096-3\|dec | Sensory perception of bitter taste | Internal cheek pouch | 1.867 | Association of Bitter Taste Receptor T2R38 Polymorphisms, Oral Microbiota, and Rheumato... |
| bsdb:39417117\|bsdb:39417117-24\|dec | Metabolic toxicity | Feces | 1.815 | Toxic Metals Impact Gut Microbiota and Metabolic Risk in Five African-Origin Populations |
| bsdb:40811957\|bsdb:40811957-7\|dec | Peri-Implantitis | Dental plaque | 1.803 | Microbiome Shifts in Peri-Implantitis: Longitudinal Characterization of Dysbiosis and R... |
| bsdb:39191760\|bsdb:39191760-11\|dec | Response to transplant | Feces | 1.744 | Fecal microbiota transplantation influences microbiota without connection to symptom re... |
| bsdb:36566315\|bsdb:36566315-1\|dec | Fasting blood glucose measurement | Feces | 1.722 | The gut microbiome of obese postpartum women with and without previous gestational diab... |
| bsdb:40520642\|bsdb:40520642-1\|dec | Type II diabetes mellitus | Feces | 1.722 | Fecal microbiome analysis in patients with metabolic syndrome and type 2 diabetes |
| bsdb:35140064\|bsdb:35140064-15\|dec | Response to vaccine | Feces | 1.722 | Gut microbiota composition is associated with SARS-CoV-2 vaccine immunogenicity and adv... |

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
- `anatomical_sites_increased_profile_diagram.png`
- `anatomical_sites_enriched_profile_diagram.png`
- `target_taxa_condition_breakdown_bar.png`
