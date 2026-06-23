# Facultative-Aerobe v2 Increased-Profile Signal Investigation

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Guild traits: facultatively aerobic (METPO:1000608), aerobic (METPO:1000602), microaerophilic (METPO:1000604)
- Target taxa: Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, Fusobacterium nucleatum, Enterocloster bolteae, Eggerthella lenta, Prevotella intermedia, Porphyromonas gingivalis, Veillonella dispar
- Formal KG-Microbe taxon universe: 2,667 taxa with at least one core oxygen-tolerance annotation among facultatively aerobic, aerobic, microaerophilic
- Guild trait selection audit: `guild_trait_selection.tsv` records the selected trait IDs, short rationale, and KG-Microbe taxon counts for each selected trait.
- Context: all body sites, with body-site and condition summaries to separate gut, oral, airway, and other patterns.

## Plot Guild Specificity

| plot | guild-specific? | basis |
| --- | --- | --- |
| `condition_increased_guild_signal.png` | yes | v2 guild score by condition |
| `guild_trait_direction_shift.png` | yes | the three selected v2 guild traits |
| `anatomical_sites_diagram.png` | yes | v2 guild score by body-site group |
| `candidate_taxa_increased_profile_tendency.png` | no | fixed target-taxon direction counts |
| `target_taxa_condition_breakdown_bar.png` | no | fixed target-taxon occurrence counts |

## Main Signal

The selected guild traits are higher in increased profiles than in decreased profiles, supporting a facultative-aerobe increased-profile signal.

| trait | increased - decreased z |
| --- | ---: |
| facultatively aerobic | 0.127 |
| aerobic | 0.085 |
| microaerophilic | 0.060 |

## Target Taxa Most Often Reported In Increased Profiles

| taxon | increased / decreased profiles | increased-profile tendency |
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

## Recurrent Target-Taxon Increased Profiles By Condition

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

| condition | profiles | studies | increased - decreased z | paired d | paired p | matched pairs | increased mean guild z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Lung transplantation | 7 | 2 | 1.213 |  |  | 2 | 0.630 |
| Age-related macular degeneration | 27 | 5 | 1.019 | 1.189 | 0.00106 | 13 | 0.504 |
| Multiple myeloma | 11 | 2 | 1.010 | 1.255 | 0.0485 | 5 | 0.387 |
| Tuberculosis | 8 | 3 | 0.964 |  |  | 1 | 0.364 |
| Sarcopenia | 26 | 2 | 0.903 | 0.298 | 0.428 | 8 | 0.172 |
| Graft versus host disease | 13 | 4 | 0.870 | 1.333 | 0.0759 | 4 | 0.163 |
| Acute graft vs. host disease | 11 | 2 | 0.782 | 0.658 | 0.28 | 4 | 0.111 |
| Autism | 14 | 4 | 0.755 | 0.538 | 0.361 | 4 | 0.459 |

Broad labels such as `Exposure`, `Population`, and labels containing `measurement` come from BugSigDB condition metadata but usually describe study design, cohort, exposure, or measurement context rather than a disease/condition context. They are excluded from displayed condition bar charts and condition summaries.

Oral-site increased-profile signals are expected in this analysis because the oral cavity is oxygen-exposed and several focal taxa are common oral or airway-associated aerotolerant, facultatively aerobic, or microaerophilic lineages. These signals should be interpreted as body-site-specific oral microbial ecology, not as a gut increased-profile signal.

## Data-Inferred Guild Evidence

For the top condition signals, this table uses all taxa observed in the increased profiles for that condition that also belong to the KG-Microbe-wide candidate taxon universe, then reports the profile v2 guild traits that are enriched in those same increased profiles and the universe taxa that carry those trait annotations. This is a data-derived guild summary, not a literature-defined facultative-aerobe/aerotolerant guild.

| condition | increased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Acute graft vs. host disease | 4 | 15 / 16 | aerobic (+0.92); facultatively aerobic (+0.71) | taxon_106589 (2 traits, 1 profile); taxon_119060 (2 traits, 1 profile); taxon_1239 (2 traits, 1 profile); taxon_1350 (2 traits, 1 profile); taxon_186801 (2 traits, 1 profile); taxon_203494 (2 traits, 1 profile); taxon_48461 (2 traits, 1 profile); taxon_74201 (2 traits, 1 profile); taxon_815 (2 traits, 1 profile); taxon_816 (2 traits, 1 profile); taxon_1647988 (1 trait, 1 profile); taxon_186803 (1 trait, 1 profile) |
| Age-related macular degeneration | 13 | 31 / 45 | facultatively aerobic (+1.32); microaerophilic (+1.19) | taxon_1301 (2 traits, 5 profiles); taxon_1378 (2 traits, 4 profiles); taxon_29317 (2 traits, 2 profiles); taxon_1224 (2 traits, 1 profile); taxon_1236 (2 traits, 1 profile); taxon_1239 (2 traits, 1 profile); taxon_1300 (2 traits, 1 profile); taxon_1306 (2 traits, 1 profile); taxon_1350 (2 traits, 1 profile); taxon_135621 (2 traits, 1 profile); taxon_1578 (2 traits, 1 profile); taxon_171549 (2 traits, 1 profile) |
| Autism | 8 | 24 / 36 | aerobic (+1.05); facultatively aerobic (+0.80); microaerophilic (+0.42) | taxon_1224 (3 traits, 3 profiles); taxon_1350 (3 traits, 3 profiles); taxon_31953 (3 traits, 3 profiles); taxon_85023 (3 traits, 3 profiles); taxon_1653 (3 traits, 2 profiles); taxon_976 (3 traits, 2 profiles); taxon_1678 (3 traits, 2 profiles); taxon_128827 (3 traits, 1 profile); taxon_1578 (3 traits, 1 profile); taxon_1716 (3 traits, 1 profile); taxon_1760 (3 traits, 1 profile); taxon_186817 (3 traits, 1 profile) |
| Graft versus host disease | 7 | 12 / 15 | aerobic (+1.24); facultatively aerobic (+1.01) | taxon_1239 (2 traits, 2 profiles); taxon_40323 (2 traits, 2 profiles); taxon_543 (2 traits, 2 profiles); taxon_1224 (2 traits, 1 profile); taxon_1236 (2 traits, 1 profile); taxon_135614 (2 traits, 1 profile); taxon_1357 (2 traits, 1 profile); taxon_32033 (2 traits, 1 profile); taxon_91347 (2 traits, 1 profile); taxon_171551 (1 trait, 1 profile); taxon_216851 (1 trait, 1 profile); taxon_47670 (1 trait, 1 profile) |
| Lung transplantation | 4 | 30 / 34 | microaerophilic (+1.37); facultatively aerobic (+1.27); aerobic (+1.00) | taxon_1129771 (3 traits, 2 profiles); taxon_1268 (3 traits, 2 profiles); taxon_1300 (3 traits, 2 profiles); taxon_1301 (3 traits, 2 profiles); taxon_1350 (3 traits, 2 profiles); taxon_1386 (3 traits, 2 profiles); taxon_1654 (3 traits, 2 profiles); taxon_186826 (3 traits, 2 profiles); taxon_186827 (3 traits, 2 profiles); taxon_201174 (3 traits, 2 profiles); taxon_203490 (3 traits, 2 profiles); taxon_203491 (3 traits, 2 profiles) |
| Multiple myeloma | 6 | 24 / 41 | facultatively aerobic (+1.48); aerobic (+0.96) | taxon_1224 (2 traits, 2 profiles); taxon_1302 (2 traits, 2 profiles); taxon_1304 (2 traits, 2 profiles); taxon_1313 (2 traits, 2 profiles); taxon_244366 (2 traits, 2 profiles); taxon_28037 (2 traits, 2 profiles); taxon_54291 (2 traits, 2 profiles); taxon_546 (2 traits, 2 profiles); taxon_548 (2 traits, 2 profiles); taxon_550 (2 traits, 2 profiles); taxon_573 (2 traits, 2 profiles); taxon_816 (2 traits, 2 profiles) |
| Sarcopenia | 13 | 59 / 155 | aerobic (+1.01); microaerophilic (+0.89); facultatively aerobic (+0.81) | taxon_1350 (3 traits, 4 profiles); taxon_547 (3 traits, 3 profiles); taxon_3031852 (3 traits, 2 profiles); taxon_1352 (3 traits, 1 profile); taxon_1624 (3 traits, 1 profile); taxon_186801 (3 traits, 1 profile); taxon_186802 (3 traits, 1 profile); taxon_1903414 (3 traits, 1 profile); taxon_204441 (3 traits, 1 profile); taxon_213849 (3 traits, 1 profile); taxon_28211 (3 traits, 1 profile); taxon_287 (3 traits, 1 profile) |
| Tuberculosis | 7 | 36 / 58 | facultatively aerobic (+1.24); aerobic (+1.08); microaerophilic (+0.58) | taxon_1301 (3 traits, 2 profiles); taxon_186826 (3 traits, 2 profiles); taxon_543 (3 traits, 2 profiles); taxon_91061 (3 traits, 2 profiles); taxon_91347 (3 traits, 2 profiles); taxon_118968 (3 traits, 1 profile); taxon_118969 (3 traits, 1 profile); taxon_119060 (3 traits, 1 profile); taxon_1236 (3 traits, 1 profile); taxon_1253 (3 traits, 1 profile); taxon_128827 (3 traits, 1 profile); taxon_1300 (3 traits, 1 profile) |

## Body-Site Check

| body_site_group | profiles | increased mean guild z | increased - decreased z |
| --- | ---: | ---: | ---: |
| Nasopharynx | 126 | 0.536 | 0.114 |
| Saliva | 638 | 0.383 | -0.017 |
| Mouth | 170 | 0.296 | -0.202 |
| Other | 2829 | 0.228 | 0.022 |
| Subgingival dental plaque | 211 | 0.219 | -0.304 |
| Vagina | 196 | 0.060 | -0.135 |
| Caecum | 165 | 0.028 | 0.091 |
| Feces | 7595 | -0.069 | 0.146 |
| Colon | 139 | -0.227 | 0.170 |
| unknown | 39 | -0.430 | 0.060 |

## Top Increased Profiles To Inspect

| profile | condition | body site | guild z | title |
| --- | --- | --- | ---: | --- |
| bsdb:32115690\|bsdb:32115690-1\|inc | Acute lymphoblastic leukemia | Feces | 1.822 | Gastrointestinal toxicity during induction treatment for childhood acute lymphoblastic ... |
| bsdb:37596536\|bsdb:37596536-10\|inc | Population | Saliva | 1.677 | Oral microbiome homogeneity across diverse human groups from southern Africa: first res... |
| bsdb:39738315\|bsdb:39738315-2\|inc | Response to diet | Feces | 1.623 | Modulating the human gut microbiome and health markers through kombucha consumption: a ... |
| bsdb:39738315\|bsdb:39738315-3\|inc | Response to diet | Feces | 1.623 | Modulating the human gut microbiome and health markers through kombucha consumption: a ... |
| bsdb:37721699\|bsdb:37721699-10\|inc | Chlamydia trachomatis infectious disease | Vagina | 1.503 | Changes in the Vaginal Microbiome During Pregnancy and the Postpartum Period in South A... |
| bsdb:38535538\|bsdb:38535538-4\|inc | Alzheimer's disease | Oral cavity | 1.486 | Oral Microbiome Stamp in Alzheimer's Disease |
| bsdb:39642873\|bsdb:39642873-1\|inc | Respiratory Syncytial Virus Infection | Nasopharynx | 1.476 | The respiratory microbiome is linked to the severity of RSV infections and the persiste... |
| bsdb:39642873\|bsdb:39642873-4\|inc | Respiratory Syncytial Virus Infection | Nasopharynx | 1.476 | The respiratory microbiome is linked to the severity of RSV infections and the persiste... |
| bsdb:36203175\|bsdb:36203175-1\|inc | Down syndrome | Oral cavity | 1.469 | Salivary microbiome in children with Down syndrome: a case-control study |
| bsdb:36171758\|bsdb:36171758-10\|inc | Transport | Nasopharynx | 1.447 | Geography, niches, and transportation influence bovine respiratory microbiome and health |

## Suggested Follow-Up

1. Split the signal by body site before interpreting it as a facultative-aerobe increased-profile signal; oral and airway profiles can dominate aerobic/microaerophilic traits.
2. In gut-only analyses, test the target taxa against inflammation markers, oxygen tolerance, virulence genes, and resistome features.
3. Compare increased-profile taxa against the anaerobe-depletion guild to see whether the same conditions show reciprocal decreased- and increased-profile patterns.

## Files

- `guild_trait_direction_shift.tsv`
- `guild_trait_selection.tsv`
- `body_site_guild_summary.tsv`
- `condition_guild_summary.tsv`
- `candidate_taxa_direction_summary.tsv`
- `candidate_taxa_trait_presence.tsv`
- `candidate_taxa_condition_direction.tsv`
- `candidate_taxa_profile_occurrences.tsv`
- `kg_microbe_candidate_taxon_universe.tsv`
- `condition_inferred_guild_summary.tsv`
- `condition_inferred_guild_trait_evidence.tsv`
- `condition_inferred_guild_taxa_evidence.tsv`
- `top_increased_guild_profiles.tsv`
- `figure_methods.md`
- `candidate_taxa_increased_profile_tendency.png`
- `condition_increased_guild_signal.png`
- `guild_trait_direction_shift.png`
- `anatomical_sites_diagram.png`
- `target_taxa_condition_breakdown_bar.png`
