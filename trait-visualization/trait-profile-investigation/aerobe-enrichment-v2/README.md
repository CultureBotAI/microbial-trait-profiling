# Aerobe v2 Trait Enrichment in Increased Profiles

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Selected traits: facultatively aerobic (METPO:1000608), aerobic (METPO:1000602), microaerophilic (METPO:1000604)
- Target taxa: Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, Fusobacterium nucleatum, Enterocloster bolteae, Eggerthella lenta, Prevotella intermedia, Porphyromonas gingivalis, Veillonella dispar (the 9 aerobe-v2 discovery taxa in `candidate_taxa_direction_summary.tsv` plus Porphyromonas gingivalis NCBITaxon:837 included as an anaerobic comparator in the BugSigDB occurrence breakdown charts and tables)
- Formal KG-Microbe taxon universe: 1,992 taxa with at least one core oxygen-tolerance annotation among facultatively aerobic, aerobic, microaerophilic
- Trait selection audit: `guild_trait_selection.tsv` records the selected trait IDs, short rationale, and KG-Microbe taxon counts for each selected trait.
- Context: all body sites, with body-site and condition summaries to separate gut, oral, airway, and other patterns.

## Plot Guild Specificity

| plot | guild-specific? | basis |
| --- | --- | --- |
| `condition_increased_guild_signal.png` | yes | v2 trait-enrichment score by condition |
| `guild_trait_direction_shift.png` | yes | the three selected v2 traits |
| `anatomical_sites_diagram.png` | yes | v2 trait-enrichment score by body-site group |
| `candidate_taxa_increased_profile_tendency.png` | no | fixed target-taxon direction counts |
| `condition_target_taxa_stacked_bar.png` | no | fixed target-taxon occurrence counts stacked by taxon |
| `target_taxa_condition_breakdown_bar.png` | no | fixed target-taxon occurrence counts |

## Main Signal

The selected traits are enriched in increased profiles relative to decreased profiles, supporting an aerobe v2 increased-profile trait-enrichment signal.

| trait | increased - decreased z |
| --- | ---: |
| facultatively aerobic | 0.199 |
| aerobic | 0.110 |
| microaerophilic | 0.106 |

## Target Taxa Most Often Reported In Increased Profiles

| taxon | increased / decreased profiles | increased-profile tendency |
| --- | ---: | ---: |
| Eggerthella lenta | 36 / 8 | 0.636 |
| Enterocloster bolteae | 27 / 6 | 0.636 |
| Fusobacterium nucleatum | 47 / 11 | 0.621 |
| Parvimonas micra | 32 / 9 | 0.561 |
| Klebsiella pneumoniae | 39 / 12 | 0.529 |
| Streptococcus anginosus | 32 / 14 | 0.391 |
| Escherichia coli | 82 / 37 | 0.378 |
| Prevotella intermedia | 23 / 12 | 0.314 |
| Veillonella dispar | 38 / 22 | 0.267 |

## Recurrent Target-Taxon Increased Profiles By Condition

| taxon | condition | body site | increased / decreased profiles | increased studies |
| --- | --- | --- | ---: | ---: |
| Fusobacterium nucleatum | Colorectal cancer | Feces | 11 / 0 | 7 |
| Parvimonas micra | Colorectal cancer | Feces | 7 / 0 | 5 |
| Escherichia coli | Colorectal cancer | Feces | 5 / 0 | 4 |
| Prevotella intermedia | Colorectal cancer | Feces | 4 / 0 | 4 |
| Enterocloster bolteae | COVID-19 | Feces | 4 / 0 | 3 |
| Klebsiella pneumoniae | Parkinson's disease | Feces | 4 / 0 | 3 |
| Streptococcus anginosus | Rheumatoid arthritis | Subgingival dental plaque | 4 / 0 | 3 |
| Veillonella dispar | Type II diabetes mellitus | Feces | 4 / 0 | 2 |
| Fusobacterium nucleatum | Ovarian cancer | Other | 4 / 1 | 2 |
| Escherichia coli | Obesity | Feces | 4 / 10 | 2 |

Target-taxon occurrence views (`target_taxa_condition_breakdown_bar.png`, `condition_target_taxa_stacked_bar.png`, `candidate_taxa_profile_occurrences.tsv` and `candidate_taxa_condition_direction.tsv`) count human-host profiles only, the corpus of the profile metadata. 72 occurrences (18 decreased, 54 increased) in 63 profiles from 19 studies come from experiments whose host is not recorded as Homo sapiens: Mus musculus 56, Rattus norvegicus 7, Sus scrofa domesticus 4, Bos taurus 2, Canis lupus familiaris 2, Elephas maximus 1. They have no profile metadata, so they are not counted under any condition; `candidate_taxa_excluded_host_occurrences.tsv` lists them.

## Contexts With Aerobe v2 Trait Enrichment in Increased Profiles

| condition | profiles | studies | increased - decreased z | paired d | paired p | matched pairs | increased mean trait z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Multiple myeloma | 9 | 2 | 1.159 |  |  | 4 | 0.353 |
| Response to antibiotic | 9 | 3 | 1.005 |  |  | 3 | 0.627 |
| Age-related macular degeneration | 24 | 5 | 0.992 | 1.105 | 0.00281 | 12 | 0.477 |
| Graft versus host disease | 14 | 5 | 0.923 |  |  | 4 | 0.182 |
| Epilepsy | 16 | 5 | 0.891 | 1.586 | 0.0116 | 6 | 0.233 |
| Cognitive impairment | 15 | 4 | 0.889 | 0.815 | 0.142 | 5 | 0.493 |
| Diarrhea | 24 | 3 | 0.806 | 0.885 | 0.0107 | 12 | 0.661 |
| Bipolar disorder | 18 | 4 | 0.796 |  |  | 4 | 0.010 |

Paired d and paired p are blank for conditions with fewer than 5 matched increased/decreased profile pairs, or with zero variance among the paired differences.

Broad labels such as `Exposure`, `Population`, `Time`, `Timepoint`, and labels containing `measurement` come from BugSigDB condition metadata but usually describe study design, cohort, sampling time, exposure, or measurement context rather than a disease/condition context. They are excluded from displayed condition bar charts and condition summaries.

Oral-site increased-profile trait signals are expected in this analysis because the oral cavity is oxygen-exposed and several focal taxa are common oral or airway-associated aerotolerant, facultatively aerobic, or microaerophilic lineages. These signals should be interpreted as body-site-specific oral microbial ecology, not as a gut increased-profile trait-enrichment signal.

## Data-Inferred Guild Evidence

For the top condition signals, this table uses all taxa observed in the increased profiles for that condition that also belong to the KG-Microbe-wide candidate taxon universe, then reports the profile v2 traits that are enriched in those same increased profiles and the universe taxa that carry those trait annotations. This is a data-derived guild summary, not a literature-defined aerobe/aerotolerant guild.

| condition | increased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Age-related macular degeneration | 12 | 29 / 42 | facultatively aerobic (+1.29); microaerophilic (+1.14); aerobic (+0.55) | Streptococcus (3 traits, 5 profiles); Pseudomonadota (3 traits, 1 profile); Gammaproteobacteria (3 traits, 1 profile); Bacillota (3 traits, 1 profile); Streptococcaceae (3 traits, 1 profile); Enterococcus (3 traits, 1 profile); Pseudomonadaceae (3 traits, 1 profile); Lactobacillus (3 traits, 1 profile); Bacteroidales (3 traits, 1 profile); Corynebacterium (3 traits, 1 profile); Lactobacillales (3 traits, 1 profile); Bacteroidia (3 traits, 1 profile) |
| Bipolar disorder | 11 | 11 / 15 | microaerophilic (+1.28) | Lactobacillus (1 trait, 4 profiles); Lactobacillaceae (1 trait, 4 profiles); Actinomycetota (1 trait, 2 profiles); Clostridiaceae (1 trait, 2 profiles); Coriobacteriia (1 trait, 2 profiles); Streptococcaceae (1 trait, 1 profile); Streptococcus (1 trait, 1 profile); Lachnospiraceae (1 trait, 1 profile); Lactobacillales (1 trait, 1 profile); Enterobacteriaceae (1 trait, 1 profile); Bacilli (1 trait, 1 profile) |
| Cognitive impairment | 8 | 11 / 25 | microaerophilic (+1.10); facultatively aerobic (+1.10); aerobic (+0.47) | Streptococcus (3 traits, 2 profiles); Lactobacillus (3 traits, 2 profiles); Lactobacillaceae (3 traits, 2 profiles); Pseudomonadota (3 traits, 1 profile); Streptococcaceae (3 traits, 1 profile); Lactococcus (3 traits, 1 profile); Lactobacillales (3 traits, 1 profile); Enterobacteriaceae (3 traits, 1 profile); Verrucomicrobiota (2 traits, 1 profile); Prevotellaceae (1 trait, 1 profile); Faecalibacterium (1 trait, 1 profile) |
| Diarrhea | 12 | 17 / 30 | facultatively aerobic (+1.03); microaerophilic (+0.83); aerobic (+0.55) | Streptococcus (3 traits, 7 profiles); Enterococcus (3 traits, 4 profiles); Escherichia (3 traits, 4 profiles); Campylobacter (3 traits, 3 profiles); Bifidobacterium (3 traits, 2 profiles); Lactobacillus (3 traits, 1 profile); Bacteroidales (3 traits, 1 profile); Bacteroidia (3 traits, 1 profile); Betaproteobacteria (3 traits, 1 profile); Enterobacteriaceae (3 traits, 1 profile); Burkholderiales (3 traits, 1 profile); Bacteroides (3 traits, 1 profile) |
| Epilepsy | 10 | 50 / 85 | facultatively aerobic (+1.27); microaerophilic (+0.83); aerobic (+0.57) | Lactobacillales (3 traits, 2 profiles); Enterobacteriaceae (3 traits, 2 profiles); Escherichia (3 traits, 2 profiles); Bacilli (3 traits, 2 profiles); Enterobacterales (3 traits, 2 profiles); Actinomycetota (3 traits, 2 profiles); Pseudomonadota (3 traits, 1 profile); Gammaproteobacteria (3 traits, 1 profile); Bacillota (3 traits, 1 profile); Erysipelotrichaceae (3 traits, 1 profile); Streptococcus (3 traits, 1 profile); Enterococcus (3 traits, 1 profile) |
| Graft versus host disease | 7 | 12 / 15 | aerobic (+1.30); facultatively aerobic (+1.07) | Bacillota (2 traits, 2 profiles); Stenotrophomonas (2 traits, 2 profiles); Enterobacteriaceae (2 traits, 2 profiles); Pseudomonadota (2 traits, 1 profile); Gammaproteobacteria (2 traits, 1 profile); Lysobacterales (2 traits, 1 profile); Lactococcus (2 traits, 1 profile); Lysobacteraceae (2 traits, 1 profile); Enterobacterales (2 traits, 1 profile); Porphyromonadaceae (1 trait, 1 profile); Faecalibacterium (1 trait, 1 profile); Lautropia (1 trait, 1 profile) |
| Multiple myeloma | 5 | 20 / 37 | facultatively aerobic (+1.60); aerobic (+1.20) | Streptococcus gordonii (2 traits, 2 profiles); Streptococcus salivarius (2 traits, 2 profiles); Streptococcus pneumoniae (2 traits, 2 profiles); Klebsiella variicola (2 traits, 2 profiles); Streptococcus mitis (2 traits, 2 profiles); Klebsiella ornithinolytica (2 traits, 2 profiles); Citrobacter freundii (2 traits, 2 profiles); Klebsiella aerogenes (2 traits, 2 profiles); Enterobacter cloacae (2 traits, 2 profiles); Klebsiella pneumoniae (2 traits, 2 profiles); Pseudomonadota (2 traits, 1 profile); Streptococcus thermophilus (2 traits, 1 profile) |
| Response to antibiotic | 5 | 42 / 65 | aerobic (+1.28); microaerophilic (+1.27); facultatively aerobic (+0.46) | Enterobacter (3 traits, 2 profiles); Streptococcus (3 traits, 1 profile); Enterococcus (3 traits, 1 profile); Lactococcus (3 traits, 1 profile); Clostridia (3 traits, 1 profile); Eubacteriales (3 traits, 1 profile); Rothia (3 traits, 1 profile); Methylobacterium (3 traits, 1 profile); Mycoplasmatota (3 traits, 1 profile); Escherichia (3 traits, 1 profile); Klebsiella (3 traits, 1 profile); Klebsiella pneumoniae (3 traits, 1 profile) |

## Body-Site Check

| body_site_group | profiles | increased mean trait z | increased - decreased z |
| --- | ---: | ---: | ---: |
| Oral cavity | 75 | 0.493 | -0.057 |
| Saliva | 374 | 0.390 | -0.029 |
| Mouth | 80 | 0.294 | -0.131 |
| Uterine cervix | 81 | 0.241 | -0.062 |
| Throat | 66 | 0.217 | 0.068 |
| Other | 1339 | 0.213 | -0.003 |
| Subgingival dental plaque | 77 | 0.107 | -0.321 |
| Vagina | 133 | -0.006 | -0.230 |
| Feces | 3604 | -0.046 | 0.233 |

## Top Increased Profiles To Inspect

| profile | condition | body site | trait z | title |
| --- | --- | --- | ---: | --- |
| bsdb:32115690\|bsdb:32115690-1\|inc | Acute lymphoblastic leukemia | Feces | 1.992 | Gastrointestinal toxicity during induction treatment for childhood acute lymphoblastic ... |
| bsdb:37596536\|bsdb:37596536-10\|inc | Population | Saliva | 1.691 | Oral microbiome homogeneity across diverse human groups from southern Africa: first res... |
| bsdb:39738315\|bsdb:39738315-2\|inc | Response to diet | Feces | 1.584 | Modulating the human gut microbiome and health markers through kombucha consumption: a ... |
| bsdb:39738315\|bsdb:39738315-3\|inc | Response to diet | Feces | 1.584 | Modulating the human gut microbiome and health markers through kombucha consumption: a ... |
| bsdb:31619731\|bsdb:31619731-4\|inc | Methamphetamine dependence | Rectum | 1.556 | Alterations to the Gastrointestinal Microbiome Associated with Methamphetamine Use amon... |
| bsdb:39940045\|bsdb:39940045-9\|inc | Diet | Feces | 1.505 | The impact of regular sauerkraut consumption on the human gut microbiota: a crossover i... |
| bsdb:39940045\|bsdb:39940045-8\|inc | Diet | Feces | 1.505 | The impact of regular sauerkraut consumption on the human gut microbiota: a crossover i... |
| bsdb:39940045\|bsdb:39940045-6\|inc | Diet | Feces | 1.505 | The impact of regular sauerkraut consumption on the human gut microbiota: a crossover i... |
| bsdb:39940045\|bsdb:39940045-5\|inc | Diet | Feces | 1.505 | The impact of regular sauerkraut consumption on the human gut microbiota: a crossover i... |
| bsdb:39940045\|bsdb:39940045-4\|inc | Diet | Feces | 1.505 | The impact of regular sauerkraut consumption on the human gut microbiota: a crossover i... |

## Suggested Follow-Up

1. Split the signal by body site before interpreting it as an aerobe v2 increased-profile trait-enrichment signal; oral and airway profiles can dominate aerobic/microaerophilic traits.
2. In gut-only analyses, test the target taxa against inflammation markers, oxygen tolerance, virulence genes, and resistome features.
3. Compare increased-profile taxa against the anaerobe v2 trait set to see whether the same conditions show reciprocal decreased- and increased-profile patterns.

## Files

- `guild_trait_direction_shift.tsv`
- `guild_trait_selection.tsv`
- `body_site_guild_summary.tsv`
- `condition_guild_summary.tsv`
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
- `top_increased_guild_profiles.tsv`
- Methods are consolidated in `../../figure_methods.md`.
- `candidate_taxa_increased_profile_tendency.png`
- `condition_increased_guild_signal.png`
- `condition_body_site_guild_map.png`
- `condition_body_site_guild_map_v2.png`
- `condition_body_site_guild_map_v3_individual.png`
- `condition_body_site_guild_map_v3_grouped.png`
- `condition_target_taxa_stacked_bar.png`
- `guild_trait_direction_shift.png`
- `anatomical_sites_diagram.png`
- `target_taxa_condition_breakdown_bar.png`
