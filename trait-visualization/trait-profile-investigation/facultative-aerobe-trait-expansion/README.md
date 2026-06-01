# Facultative-Aerobe Trait Expansion Investigation

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Guild traits: facultatively aerobic (METPO:1000608), aerobic (METPO:1000602), gram negative (METPO:1000699), motile (METPO:1000702), microaerophilic (METPO:1000604), nitrate_reduction (METPO:2000103|pathways:nitrate_reduction)
- Target taxa: Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, Fusobacterium nucleatum, Enterocloster bolteae, Eggerthella lenta, Prevotella intermedia, Porphyromonas gingivalis, Veillonella dispar
- Context: all body sites, with body-site and condition summaries to separate gut, oral, airway, and other patterns.

## Main Signal

The selected guild traits are higher in increased profiles than in decreased profiles, supporting a facultative-aerobe trait expansion signal.

| trait | increased - decreased z |
| --- | ---: |
| facultatively aerobic | 0.127 |
| aerobic | 0.085 |
| gram negative | 0.066 |
| microaerophilic | 0.060 |
| motile | 0.056 |
| nitrate_reduction | 0.029 |

## Strongest Expansion-Biased Taxa

| taxon | increased / decreased profiles | expansion bias |
| --- | ---: | ---: |
| Parvimonas micra | 55 / 17 | 0.528 |
| Klebsiella pneumoniae | 64 / 22 | 0.488 |
| Streptococcus anginosus | 54 / 19 | 0.479 |
| Escherichia coli | 185 / 66 | 0.474 |
| Fusobacterium nucleatum | 76 / 28 | 0.462 |
| Enterocloster bolteae | 56 / 21 | 0.455 |
| Eggerthella lenta | 51 / 20 | 0.437 |
| Prevotella intermedia | 48 / 23 | 0.352 |
| Porphyromonas gingivalis | 50 / 24 | 0.351 |
| Veillonella dispar | 67 / 34 | 0.327 |

## Recurrent Target-Taxon Expansions By Condition

| taxon | condition | body site | increased / decreased profiles | increased studies |
| --- | --- | --- | ---: | ---: |
| Fusobacterium nucleatum | Colorectal cancer | Feces | 16 / 0 | 10 |
| Parvimonas micra | Colorectal cancer | Feces | 14 / 0 | 9 |
| Prevotella intermedia | Colorectal cancer | Feces | 11 / 0 | 7 |
| Escherichia coli | Colorectal cancer | Feces | 6 / 0 | 5 |
| Enterocloster bolteae | COVID-19 | Feces | 7 / 0 | 4 |
| Fusobacterium nucleatum | Smoking behavior | Subgingival dental plaque | 5 / 0 | 3 |
| Klebsiella pneumoniae | Parkinson's disease | Feces | 4 / 0 | 3 |
| Streptococcus anginosus | Rheumatoid arthritis | Subgingival dental plaque | 4 / 0 | 3 |
| Enterocloster bolteae | Colorectal cancer | Feces | 3 / 0 | 3 |
| Enterocloster bolteae | Type II diabetes mellitus | Feces | 14 / 0 | 2 |
| Escherichia coli | Type II diabetes mellitus | Feces | 14 / 0 | 2 |
| Porphyromonas gingivalis | Colorectal cancer | Saliva | 6 / 0 | 2 |

## Contexts With High Increased-Profile Guild Signal

| condition | profiles | studies | increased mean guild z | increased - decreased z |
| --- | ---: | ---: | ---: | ---: |
| Lung transplantation | 7 | 2 | 0.453 | 0.845 |
| Female infertility | 22 | 3 | 0.284 | 0.694 |
| Multiple myeloma | 11 | 2 | 0.515 | 0.681 |
| Tuberculosis | 8 | 3 | 0.336 | 0.678 |
| Autism | 14 | 4 | 0.530 | 0.597 |
| Graft versus host disease | 13 | 4 | 0.324 | 0.596 |
| Colorectal carcinoma | 15 | 3 | 0.492 | 0.570 |
| Vulvovaginitis | 5 | 2 | 0.307 | 0.543 |

## Body-Site Check

| body_site_group | profiles | increased mean guild z | increased - decreased z |
| --- | ---: | ---: | ---: |
| Nasopharynx | 126 | 0.349 | 0.221 |
| Saliva | 638 | 0.259 | -0.068 |
| Subgingival dental plaque | 211 | 0.200 | -0.159 |
| Other | 2829 | 0.177 | 0.035 |
| Caecum | 165 | 0.050 | 0.111 |
| Mouth | 170 | 0.026 | -0.164 |
| Feces | 7595 | -0.042 | 0.105 |
| Vagina | 196 | -0.085 | 0.014 |
| Colon | 139 | -0.110 | 0.112 |
| unknown | 39 | -0.602 | -0.114 |

## Top Increased Profiles To Inspect

| profile | condition | body site | guild z | title |
| --- | --- | --- | ---: | --- |
| bsdb:40232948\|bsdb:40232948-1\|inc | Type II diabetes mellitus | Buccal mucosa | 1.349 | Consensus approach to differential abundance analysis detects few differences in the or... |
| bsdb:40232948\|bsdb:40232948-2\|inc | Type II diabetes mellitus | Buccal mucosa | 1.349 | Consensus approach to differential abundance analysis detects few differences in the or... |
| bsdb:36203175\|bsdb:36203175-1\|inc | Down syndrome | Oral cavity | 1.314 | Salivary microbiome in children with Down syndrome: a case-control study |
| bsdb:40827923\|bsdb:40827923-1\|inc | Benign prostatic hyperplasia,Periodontitis | Gingival groove | 1.300 | Profiles of oral microbiota and metabolites in periodontitis and benign prostatic hyper... |
| bsdb:39325441\|bsdb:39325441-2\|inc | Head and neck squamous cell carcinoma | Oral cavity | 1.271 | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer |
| bsdb:36090083\|bsdb:36090083-8\|inc | Small for gestational age | Feces | 1.268 | Characteristics of gut microbiota of term small gestational age infants within 1 week a... |
| bsdb:40495530\|bsdb:40495530-1\|inc | Parkinson's disease | Blood | 1.246 | Investigating the Blood Microbiome in Parkinson's Disease, Schizophrenia, and Posttraum... |
| bsdb:35093093\|bsdb:35093093-3\|inc | Chronic obstructive pulmonary disease | Saliva | 1.233 | Saliva as a non-invasive specimen for COPD assessment |
| bsdb:37936242\|bsdb:37936242-7\|inc | Response to surgery | Feces | 1.214 | Surgery-induced gut microbial dysbiosis promotes cognitive impairment via regulation of... |
| bsdb:34248903\|bsdb:34248903-2\|inc | Dental caries | Saliva,Supragingival dental plaque | 1.213 | Characterization of Supragingival Plaque and Oral Swab Microbiomes in Children With Sev... |

## Suggested Follow-Up

1. Split the signal by body site before interpreting it as facultative-aerobe trait expansion; oral and airway profiles can dominate aerobic/microaerophilic traits.
2. In gut-only analyses, test the target taxa against inflammation markers, nitrate respiration, oxygen tolerance, virulence genes, and resistome features.
3. Compare expansion-biased taxa against the anaerobe-depletion guild to see whether the same conditions show reciprocal loss and expansion patterns.

## Files

- `guild_trait_direction_shift.tsv`
- `body_site_guild_summary.tsv`
- `condition_guild_summary.tsv`
- `candidate_taxa_direction_summary.tsv`
- `candidate_taxa_trait_presence.tsv`
- `candidate_taxa_condition_direction.tsv`
- `candidate_taxa_profile_occurrences.tsv`
- `top_increased_guild_profiles.tsv`
- `figure_methods.md`
- `candidate_taxa_expansion_bias.png`
- `condition_increased_guild_signal.png`
- `guild_trait_direction_shift.png`
- `anatomical_sites_diagram.png`
- `target_taxa_condition_breakdown_pie.png`
- `target_taxa_condition_breakdown_histogram.png`
