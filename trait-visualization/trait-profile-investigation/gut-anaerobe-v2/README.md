# Gut Anaerobe v2 Depletion Investigation

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Guild traits: anaerobic (METPO:1000603), obligately anaerobic (METPO:1000607)
- Target taxa: Faecalibacterium prausnitzii, Anaerobutyricum hallii, Lachnospira eligens, Agathobacter rectalis, Roseburia intestinalis, Ruminococcus bromii, Roseburia inulinivorans
- Gut context: body-site labels matching feces/faeces/stool/colon/caecum/cecum/gut/intestin/digestive, plus matching grouped body-site labels.

## Plot Guild Specificity

| plot | guild-specific? | basis |
| --- | --- | --- |
| `condition_decreased_guild_signal.png` | yes | v2 guild score by condition |
| `guild_trait_direction_shift.png` | yes | the two selected v2 guild traits |
| `anatomical_sites_diagram.png` | yes | v2 guild score by body-site group |
| `candidate_taxa_depletion_bias.png` | no | fixed target-taxon direction counts |
| `target_taxa_condition_breakdown_bar.png` | no | fixed target-taxon occurrence counts |

## Main Signal

The selected anaerobe-only guild traits are consistently higher in decreased gut-context profiles than in increased gut-context profiles, which supports the gut anaerobe v2 depletion hypothesis.

| trait | decreased - increased z |
| --- | ---: |
| anaerobic | 0.138 |
| obligately anaerobic | 0.127 |

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

## Gut Contexts With High Decreased-Profile Guild Signal

| condition | profiles | studies | decreased - increased z | paired d | paired p | matched pairs | decreased mean guild z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Graft versus host disease | 13 | 4 | 1.254 | 1.480 | 0.0595 | 4 | 0.848 |
| Age-related macular degeneration | 10 | 3 | 1.169 | 2.298 | 0.0068 | 5 | 0.935 |
| Acute graft vs. host disease | 11 | 2 | 1.080 | 0.675 | 0.27 | 4 | 0.569 |
| Health study participation | 28 | 4 | 1.012 | 0.483 | 0.185 | 9 | 0.194 |
| Timepoint | 5 | 2 | 0.946 |  |  | 1 | 0.349 |
| Aging | 43 | 4 | 0.825 | 0.873 | 0.00176 | 18 | 0.284 |
| Acute pancreatitis | 20 | 2 | 0.818 | 0.977 | 0.0129 | 10 | 0.394 |
| Graves ophthalmopathy | 13 | 2 | 0.818 | 0.887 | 0.118 | 5 | 0.023 |

## Body-Site Check

| body_site_group | profiles | decreased mean guild z | decreased - increased z |
| --- | ---: | ---: | ---: |
| Colon | 139 | 0.445 | 0.366 |
| Feces | 7595 | 0.231 | 0.152 |
| Caecum | 165 | 0.100 | -0.069 |
| Other | 391 | -0.184 | -0.229 |

## Top Decreased Profiles To Inspect

| profile | condition | body site | guild z | title |
| --- | --- | --- | ---: | --- |
| bsdb:36621696\|bsdb:36621696-2\|dec | Psoriasis vulgaris | Feces | 1.921 | Altered gut microbiome composition in nontreated plaque psoriasis patients |
| bsdb:39976263\|bsdb:39976263-10\|dec | Alanine measurement | Feces | 1.921 | <i>Roseburia hominis</i> improves host metabolism in diet-induced obesity |
| bsdb:31171880\|bsdb:31171880-2\|dec | Colorectal cancer | Feces | 1.891 | Metagenomic and metabolomic analyses reveal distinct stage-specific phenotypes of the g... |
| bsdb:39417117\|bsdb:39417117-24\|dec | Metabolic toxicity | Feces | 1.815 | Toxic Metals Impact Gut Microbiota and Metabolic Risk in Five African-Origin Populations |
| bsdb:39191760\|bsdb:39191760-11\|dec | Response to transplant | Feces | 1.743 | Fecal microbiota transplantation influences microbiota without connection to symptom re... |
| bsdb:36566315\|bsdb:36566315-1\|dec | Fasting blood glucose measurement | Feces | 1.721 | The gut microbiome of obese postpartum women with and without previous gestational diab... |
| bsdb:40520642\|bsdb:40520642-1\|dec | Type II diabetes mellitus | Feces | 1.721 | Fecal microbiome analysis in patients with metabolic syndrome and type 2 diabetes |
| bsdb:35140064\|bsdb:35140064-15\|dec | Response to vaccine | Feces | 1.721 | Gut microbiota composition is associated with SARS-CoV-2 vaccine immunogenicity and adv... |
| bsdb:37894412\|bsdb:37894412-2\|dec | Intestinal polyp | Colonic mucosa | 1.721 | Mucosal Microbiome in Patients with Early Bowel Polyps: Inferences from Short-Read and ... |
| bsdb:31530647\|bsdb:31530647-2\|dec | Colorectal cancer | Feces | 1.682 | Altered Fecal Small RNA Profiles in Colorectal Cancer Reflect Gut Microbiome Compositio... |

## Suggested Follow-Up

1. Treat the target taxa as a depletion guild and test enrichment by condition/body site rather than analyzing each taxon in isolation.
2. Prioritize conditions where decreased profiles have high guild scores and target taxa are repeatedly decreased across multiple studies.
3. Compare this anaerobe-only v2 signal against the original five-trait guild to separate oxygen-tolerance signal from metabolic-capability signal.

## Files

- `guild_trait_direction_shift.tsv`
- `gut_body_site_guild_summary.tsv`
- `gut_condition_guild_summary.tsv`
- `candidate_taxa_direction_summary.tsv`
- `candidate_taxa_trait_presence.tsv`
- `candidate_taxa_condition_direction.tsv`
- `candidate_taxa_profile_occurrences.tsv`
- `top_decreased_guild_profiles.tsv`
- `figure_methods.md`
- `candidate_taxa_depletion_bias.png`
- `condition_decreased_guild_signal.png`
- `guild_trait_direction_shift.png`
- `anatomical_sites_diagram.png`
- `target_taxa_condition_breakdown_bar.png`
