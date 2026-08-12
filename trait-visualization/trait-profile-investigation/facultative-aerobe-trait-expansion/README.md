# Facultative-Aerobe Increased-Profile Signal Investigation

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Guild traits: facultatively aerobic (METPO:1000608), aerobic (METPO:1000602), gram negative (METPO:1000699), motile (METPO:1000702), microaerophilic (METPO:1000604), nitrate_reduction (METPO:2000103|pathways:nitrate_reduction)
- Target taxa: Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, Fusobacterium nucleatum, Enterocloster bolteae, Eggerthella lenta, Prevotella intermedia, Porphyromonas gingivalis, Veillonella dispar
- Formal KG-Microbe taxon universe: 2,779 taxa with at least one core oxygen/respiration criterion annotation facultatively aerobic, aerobic, microaerophilic, nitrate_reduction
- Context: all body sites, with body-site and condition summaries to separate gut, oral, airway, and other patterns.

## Main Signal

The selected traits are enriched in increased profiles relative to decreased profiles, supporting a facultative-aerobe increased-profile signal.

| trait | increased - decreased z |
| --- | ---: |
| facultatively aerobic | 0.127 |
| aerobic | 0.085 |
| gram negative | 0.066 |
| microaerophilic | 0.060 |
| motile | 0.056 |
| nitrate_reduction | 0.029 |

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
| Lung transplantation | 7 | 2 | 0.845 |  |  | 2 | 0.453 |
| Female infertility | 22 | 3 | 0.694 | 1.172 | 0.101 | 4 | 0.284 |
| Multiple myeloma | 11 | 2 | 0.681 | 2.610 | 0.0043 | 5 | 0.515 |
| Tuberculosis | 8 | 3 | 0.678 |  |  | 1 | 0.336 |
| Autism | 14 | 4 | 0.597 | 0.745 | 0.233 | 4 | 0.530 |
| Graft versus host disease | 13 | 4 | 0.596 | 1.431 | 0.0645 | 4 | 0.324 |
| Colorectal carcinoma | 15 | 3 | 0.570 | 0.500 | 0.234 | 7 | 0.492 |
| Vulvovaginitis | 5 | 2 | 0.543 |  |  | 2 | 0.307 |

Broad labels such as `Exposure`, `Population`, and labels containing `measurement` come from BugSigDB condition metadata but usually describe study design, cohort, exposure, or measurement context rather than a disease/condition context. They are excluded from displayed condition bar charts and condition summaries.

Oral-site increased-profile signals are expected in this analysis because the oral cavity is oxygen-exposed and several focal taxa are common oral or airway-associated aerotolerant, facultatively aerobic, or microaerophilic lineages. These signals should be interpreted as body-site-specific oral microbial ecology, not as a gut increased-profile signal.

## Data-Inferred Guild Evidence

For the top condition signals, this table uses all taxa observed in the increased profiles for that condition that also belong to the KG-Microbe-wide candidate taxon universe, then reports the profile traits that are enriched in those same increased profiles and the universe taxa that carry those trait annotations. This is a data-derived guild summary, not a literature-defined facultative-aerobe/aerotolerant guild.

| condition | increased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Autism | 8 | 27 / 36 | GC mid1 (+1.12); aerobic (+1.05); GC mid2 (+1.00); bacillus shaped (+0.83); obligately aerobic (+0.82) | taxon_1224 (12 traits, 3 profiles); taxon_85023 (12 traits, 3 profiles); taxon_1653 (12 traits, 2 profiles); taxon_976 (12 traits, 2 profiles); taxon_1716 (12 traits, 1 profile); taxon_1760 (12 traits, 1 profile); taxon_186817 (12 traits, 1 profile); taxon_201174 (12 traits, 1 profile); taxon_200643 (11 traits, 1 profile); taxon_543 (11 traits, 1 profile); taxon_561 (10 traits, 1 profile); taxon_186824 (9 traits, 1 profile) |
| Colorectal carcinoma | 7 | 53 / 91 | nitrate_reduction (+1.00); gram negative (+0.96); GC mid1 (+0.94); filament shaped (+0.92); microaerophilic (+0.88) | taxon_1224 (12 traits, 3 profiles); taxon_1716 (11 traits, 1 profile); taxon_80864 (11 traits, 1 profile); taxon_186826 (10 traits, 1 profile); taxon_41275 (10 traits, 1 profile); taxon_543 (10 traits, 1 profile); taxon_13687 (9 traits, 1 profile); taxon_1663 (9 traits, 1 profile); taxon_1827 (9 traits, 1 profile); taxon_203490 (8 traits, 3 profiles); taxon_1300 (8 traits, 1 profile); taxon_1301 (8 traits, 1 profile) |
| Female infertility | 13 | 99 / 125 | temperature range low (+1.27); obligately aerobic (+1.25); temperature range mid2 (+1.13); temperature range mid3 (+1.12); temperature range mid1 (+1.12) | taxon_135621 (12 traits, 3 profiles); taxon_286 (12 traits, 3 profiles); taxon_1386 (12 traits, 2 profiles); taxon_469 (12 traits, 2 profiles); taxon_1224 (12 traits, 1 profile); taxon_1236 (12 traits, 1 profile); taxon_1239 (12 traits, 1 profile); taxon_135623 (12 traits, 1 profile); taxon_1385 (12 traits, 1 profile); taxon_171549 (12 traits, 1 profile); taxon_186801 (12 traits, 1 profile); taxon_186802 (12 traits, 1 profile) |
| Graft versus host disease | 7 | 12 / 15 | aerobic (+1.24); facultatively aerobic (+1.01); obligately aerobic (+0.92); GC mid1 (+0.53); nitrate_reduction (+0.43) | taxon_1239 (12 traits, 2 profiles); taxon_40323 (12 traits, 2 profiles); taxon_1224 (12 traits, 1 profile); taxon_1236 (12 traits, 1 profile); taxon_135614 (12 traits, 1 profile); taxon_32033 (12 traits, 1 profile); taxon_91347 (12 traits, 1 profile); taxon_543 (11 traits, 2 profiles); taxon_171551 (7 traits, 1 profile); taxon_1357 (4 traits, 1 profile); taxon_216851 (3 traits, 1 profile); taxon_47670 (2 traits, 1 profile) |
| Lung transplantation | 4 | 30 / 34 | GC mid1 (+1.66); sphere shaped (+1.45); microaerophilic (+1.37); temperature optimum mid2 (+1.33); facultatively aerobic (+1.27) | taxon_1268 (12 traits, 2 profiles); taxon_201174 (12 traits, 2 profiles); taxon_2037 (12 traits, 2 profiles); taxon_2049 (12 traits, 2 profiles); taxon_85006 (12 traits, 2 profiles); taxon_91061 (12 traits, 2 profiles); taxon_186802 (12 traits, 1 profile); taxon_32207 (12 traits, 1 profile); taxon_91347 (12 traits, 1 profile); taxon_1386 (11 traits, 2 profiles); taxon_1654 (11 traits, 2 profiles); taxon_186826 (11 traits, 2 profiles) |
| Multiple myeloma | 6 | 25 / 41 | facultatively aerobic (+1.48); gram negative (+1.12); non motile (+1.04); aerobic (+0.96); coccus shaped (+0.95) | taxon_1224 (12 traits, 2 profiles); taxon_815 (11 traits, 1 profile); taxon_287 (10 traits, 1 profile); taxon_1313 (8 traits, 2 profiles); taxon_816 (8 traits, 2 profiles); taxon_1328 (8 traits, 1 profile); taxon_573 (7 traits, 2 profiles); taxon_28141 (7 traits, 1 profile); taxon_29465 (7 traits, 1 profile); taxon_413496 (7 traits, 1 profile); taxon_1303 (6 traits, 2 profiles); taxon_54291 (6 traits, 2 profiles) |
| Tuberculosis | 7 | 37 / 58 | motile (+2.02); temperature range mid4 (+1.78); trehalose (+1.72); temperature range mid3 (+1.72); galactose (+1.61) | taxon_1301 (12 traits, 2 profiles); taxon_186826 (12 traits, 2 profiles); taxon_543 (12 traits, 2 profiles); taxon_91061 (12 traits, 2 profiles); taxon_91347 (12 traits, 2 profiles); taxon_186803 (12 traits, 2 profiles); taxon_1236 (12 traits, 1 profile); taxon_1300 (12 traits, 1 profile); taxon_186804 (12 traits, 1 profile); taxon_31953 (12 traits, 1 profile); taxon_33958 (12 traits, 1 profile); taxon_85004 (12 traits, 1 profile) |
| Vulvovaginitis | 3 | 28 / 39 | bacillus shaped (+1.90); gram negative (+1.84); cell width very small (+1.80); cell length very small (+1.70); obligately anaerobic (+1.58) | taxon_1224 (12 traits, 1 profile); taxon_171549 (12 traits, 1 profile); taxon_171552 (12 traits, 1 profile); taxon_186801 (12 traits, 1 profile); taxon_186802 (12 traits, 1 profile); taxon_200643 (12 traits, 1 profile); taxon_201174 (12 traits, 1 profile); taxon_2037 (12 traits, 1 profile); taxon_2049 (12 traits, 1 profile); taxon_31953 (12 traits, 1 profile); taxon_85004 (12 traits, 1 profile); taxon_976 (12 traits, 1 profile) |

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

1. Split the signal by body site before interpreting it as a facultative-aerobe increased-profile signal; oral and airway profiles can dominate aerobic/microaerophilic traits.
2. In gut-only analyses, test the target taxa against inflammation markers, nitrate respiration, oxygen tolerance, virulence genes, and resistome features.
3. Compare increased-profile taxa against the decreased-profile anaerobe trait set to see whether the same conditions show reciprocal decreased- and increased-profile patterns.

## Files

- `guild_trait_direction_shift.tsv`
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
