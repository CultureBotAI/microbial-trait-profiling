# Aerobe v2 Trait Enrichment in Increased Profiles

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Selected traits: facultatively aerobic (METPO:1000608), aerobic (METPO:1000602), microaerophilic (METPO:1000604)
- Target taxa: Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, Fusobacterium nucleatum, Enterocloster bolteae, Eggerthella lenta, Prevotella intermedia, Porphyromonas gingivalis, Veillonella dispar
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
| facultatively aerobic | 0.206 |
| aerobic | 0.114 |
| microaerophilic | 0.112 |

## Target Taxa Most Often Reported In Increased Profiles

| taxon | increased / decreased profiles | increased-profile tendency |
| --- | ---: | ---: |
| Enterocloster bolteae | 27 / 6 | 0.636 |
| Eggerthella lenta | 36 / 8 | 0.636 |
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
| Fusobacterium nucleatum | Colorectal cancer | Feces | 12 / 0 | 7 |
| Parvimonas micra | Colorectal cancer | Feces | 7 / 0 | 5 |
| Escherichia coli | Colorectal cancer | Feces | 5 / 0 | 4 |
| Prevotella intermedia | Colorectal cancer | Feces | 5 / 0 | 4 |
| Enterocloster bolteae | COVID-19 | Feces | 4 / 0 | 3 |
| Klebsiella pneumoniae | Parkinson's disease | Feces | 4 / 0 | 3 |
| Streptococcus anginosus | Rheumatoid arthritis | Subgingival dental plaque | 4 / 0 | 3 |
| Porphyromonas gingivalis | Colorectal cancer | Saliva | 6 / 0 | 2 |
| Escherichia coli | Ulcerative colitis | Feces | 3 / 0 | 2 |
| Fusobacterium nucleatum | Ovarian cancer | Other | 4 / 2 | 2 |
| Escherichia coli | Obesity | Feces | 4 / 10 | 2 |

## Contexts With Aerobe v2 Trait Enrichment in Increased Profiles

| condition | profiles | studies | increased - decreased z | paired d | paired p | matched pairs | increased mean trait z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Multiple myeloma | 9 | 2 | 1.177 | 1.422 | 0.0654 | 4 | 0.370 |
| Age-related macular degeneration | 27 | 5 | 1.010 | 1.188 | 0.00106 | 13 | 0.531 |
| Response to antibiotic | 7 | 2 | 0.987 | 0.590 | 0.415 | 3 | 0.579 |
| Graft versus host disease | 13 | 4 | 0.862 | 1.336 | 0.0756 | 4 | 0.192 |
| Diarrhea | 25 | 3 | 0.833 | 0.910 | 0.00922 | 12 | 0.657 |
| Bipolar disorder | 18 | 4 | 0.816 | 0.964 | 0.149 | 4 | 0.010 |
| Stroke | 44 | 5 | 0.785 | 0.953 | 0.00243 | 15 | 0.231 |
| High fat diet | 5 | 2 | 0.773 |  |  | 2 | 0.369 |

Broad labels such as `Exposure`, `Population`, and labels containing `measurement` come from BugSigDB condition metadata but usually describe study design, cohort, exposure, or measurement context rather than a disease/condition context. They are excluded from displayed condition bar charts and condition summaries.

Oral-site increased-profile trait signals are expected in this analysis because the oral cavity is oxygen-exposed and several focal taxa are common oral or airway-associated aerotolerant, facultatively aerobic, or microaerophilic lineages. These signals should be interpreted as body-site-specific oral microbial ecology, not as a gut increased-profile trait-enrichment signal.

## Data-Inferred Guild Evidence

For the top condition signals, this table uses all taxa observed in the increased profiles for that condition that also belong to the KG-Microbe-wide candidate taxon universe, then reports the profile v2 traits that are enriched in those same increased profiles and the universe taxa that carry those trait annotations. This is a data-derived guild summary, not a literature-defined aerobe/aerotolerant guild.

| condition | increased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Age-related macular degeneration | 13 | 31 / 45 | facultatively aerobic (+1.30); microaerophilic (+1.19) | Streptococcus (2 traits, 5 profiles); Gemella (2 traits, 4 profiles); Actinomyces sp. (2 traits, 2 profiles); Pseudomonadota (2 traits, 1 profile); Gammaproteobacteria (2 traits, 1 profile); Bacillota (2 traits, 1 profile); Streptococcaceae (2 traits, 1 profile); Streptococcus sp. (2 traits, 1 profile); Enterococcus (2 traits, 1 profile); Pseudomonadaceae (2 traits, 1 profile); Lactobacillus (2 traits, 1 profile); Bacteroidales (2 traits, 1 profile) |
| Bipolar disorder | 11 | 11 / 15 | microaerophilic (+1.27) | Lactobacillus (1 trait, 4 profiles); Lactobacillaceae (1 trait, 4 profiles); Actinomycetota (1 trait, 2 profiles); Clostridiaceae (1 trait, 2 profiles); Coriobacteriia (1 trait, 2 profiles); Streptococcaceae (1 trait, 1 profile); Streptococcus (1 trait, 1 profile); Lachnospiraceae (1 trait, 1 profile); Lactobacillales (1 trait, 1 profile); Enterobacteriaceae (1 trait, 1 profile); Bacilli (1 trait, 1 profile) |
| Diarrhea | 13 | 23 / 36 | facultatively aerobic (+0.97); microaerophilic (+0.85); aerobic (+0.67) | Streptococcus (3 traits, 7 profiles); Enterococcus (3 traits, 4 profiles); Escherichia (3 traits, 4 profiles); Campylobacter (3 traits, 3 profiles); Bifidobacterium (3 traits, 2 profiles); Lactobacillus (3 traits, 1 profile); Bacteroidales (3 traits, 1 profile); Bacteroidia (3 traits, 1 profile); Actinomycetota (3 traits, 1 profile); Betaproteobacteria (3 traits, 1 profile); Bifidobacteriaceae (3 traits, 1 profile); Enterobacteriaceae (3 traits, 1 profile) |
| Graft versus host disease | 7 | 12 / 15 | aerobic (+1.24); facultatively aerobic (+0.99) | Bacillota (2 traits, 2 profiles); Stenotrophomonas (2 traits, 2 profiles); Enterobacteriaceae (2 traits, 2 profiles); Pseudomonadota (2 traits, 1 profile); Gammaproteobacteria (2 traits, 1 profile); Lysobacterales (2 traits, 1 profile); Lactococcus (2 traits, 1 profile); Lysobacteraceae (2 traits, 1 profile); Enterobacterales (2 traits, 1 profile); Porphyromonadaceae (1 trait, 1 profile); Faecalibacterium (1 trait, 1 profile); Lautropia (1 trait, 1 profile) |
| High fat diet | 3 | 6 / 9 | facultatively aerobic (+1.55); microaerophilic (+0.89) | Streptococcus (2 traits, 1 profile); Enterococcus (2 traits, 1 profile); Lactococcus (2 traits, 1 profile); Escherichia (2 traits, 1 profile); Desulfovibrio (2 traits, 1 profile); Granulicatella (1 trait, 1 profile) |
| Multiple myeloma | 5 | 20 / 37 | facultatively aerobic (+1.62); aerobic (+1.22) | Streptococcus gordonii (2 traits, 2 profiles); Streptococcus salivarius (2 traits, 2 profiles); Streptococcus pneumoniae (2 traits, 2 profiles); Klebsiella variicola (2 traits, 2 profiles); Streptococcus mitis (2 traits, 2 profiles); Klebsiella ornithinolytica (2 traits, 2 profiles); Citrobacter freundii (2 traits, 2 profiles); Klebsiella aerogenes (2 traits, 2 profiles); Enterobacter cloacae (2 traits, 2 profiles); Klebsiella pneumoniae (2 traits, 2 profiles); Pseudomonadota (2 traits, 1 profile); Streptococcus thermophilus (2 traits, 1 profile) |
| Response to antibiotic | 3 | 38 / 61 | facultatively aerobic (+1.20); microaerophilic (+0.97); aerobic (+0.80) | Enterobacter (3 traits, 2 profiles); Streptococcus (3 traits, 1 profile); Enterococcus (3 traits, 1 profile); Lactococcus (3 traits, 1 profile); Clostridia (3 traits, 1 profile); Eubacteriales (3 traits, 1 profile); Methylobacterium (3 traits, 1 profile); Mycoplasmatota (3 traits, 1 profile); Escherichia (3 traits, 1 profile); Klebsiella (3 traits, 1 profile); Klebsiella pneumoniae (3 traits, 1 profile); Enterococcus raffinosus (3 traits, 1 profile) |
| Stroke | 29 | 68 / 125 | facultatively aerobic (+0.94); microaerophilic (+0.75); aerobic (+0.67) | Enterobacteriaceae (3 traits, 8 profiles); Enterobacterales (3 traits, 8 profiles); Pseudomonadota (3 traits, 7 profiles); Enterococcus (3 traits, 7 profiles); Gammaproteobacteria (3 traits, 6 profiles); Lactobacillales (3 traits, 6 profiles); Enterococcaceae (3 traits, 6 profiles); Lactobacillus (3 traits, 5 profiles); Bacilli (3 traits, 5 profiles); Lactobacillaceae (3 traits, 3 profiles); Fusobacteriia (3 traits, 2 profiles); Fusobacteriales (3 traits, 2 profiles) |

## Body-Site Check

| body_site_group | profiles | increased mean trait z | increased - decreased z |
| --- | ---: | ---: | ---: |
| Nasopharynx | 74 | 0.518 | 0.079 |
| Oral cavity | 82 | 0.498 | -0.060 |
| Saliva | 299 | 0.470 | 0.033 |
| Mouth | 74 | 0.314 | -0.113 |
| Other | 1352 | 0.215 | 0.008 |
| Uterine cervix | 82 | 0.192 | -0.030 |
| Subgingival dental plaque | 70 | 0.067 | -0.441 |
| Vagina | 124 | -0.006 | -0.276 |
| Feces | 3672 | -0.047 | 0.229 |

## Top Increased Profiles To Inspect

| profile | condition | body site | trait z | title |
| --- | --- | --- | ---: | --- |
| bsdb:32115690\|bsdb:32115690-1\|inc | Acute lymphoblastic leukemia | Feces | 1.842 | Gastrointestinal toxicity during induction treatment for childhood acute lymphoblastic ... |
| bsdb:37596536\|bsdb:37596536-10\|inc | Population | Saliva | 1.698 | Oral microbiome homogeneity across diverse human groups from southern Africa: first res... |
| bsdb:39738315\|bsdb:39738315-2\|inc | Response to diet | Feces | 1.643 | Modulating the human gut microbiome and health markers through kombucha consumption: a ... |
| bsdb:39738315\|bsdb:39738315-3\|inc | Response to diet | Feces | 1.643 | Modulating the human gut microbiome and health markers through kombucha consumption: a ... |
| bsdb:37721699\|bsdb:37721699-10\|inc | Chlamydia trachomatis infectious disease | Vagina | 1.523 | Changes in the Vaginal Microbiome During Pregnancy and the Postpartum Period in South A... |
| bsdb:39642873\|bsdb:39642873-1\|inc | Respiratory Syncytial Virus Infection | Nasopharynx | 1.498 | The respiratory microbiome is linked to the severity of RSV infections and the persiste... |
| bsdb:39642873\|bsdb:39642873-4\|inc | Respiratory Syncytial Virus Infection | Nasopharynx | 1.498 | The respiratory microbiome is linked to the severity of RSV infections and the persiste... |
| bsdb:36203175\|bsdb:36203175-1\|inc | Down syndrome | Oral cavity | 1.494 | Salivary microbiome in children with Down syndrome: a case-control study |
| bsdb:31339887\|bsdb:31339887-28\|inc | Waist circumference | Feces | 1.453 | The human microbiota is associated with cardiometabolic risk across the epidemiologic t... |
| bsdb:37587195\|bsdb:37587195-2\|inc | Health study participation | Feces | 1.453 | A comparison of the infant gut microbiome before versus after the start of the covid-19... |

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
