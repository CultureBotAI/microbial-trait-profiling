# Aerobe v2 Trait Enrichment in Increased Profiles

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Selected traits: facultatively aerobic (METPO:1000608), aerobic (METPO:1000602), microaerophilic (METPO:1000604)
- Target taxa: Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, Fusobacterium nucleatum, Enterocloster bolteae, Eggerthella lenta, Prevotella intermedia, Porphyromonas gingivalis, Veillonella dispar (the 9 aerobe-v2 discovery taxa in `candidate_taxa_direction_summary.tsv` plus Porphyromonas gingivalis NCBITaxon:837 included as an anaerobic comparator in the BugSigDB occurrence breakdown charts and tables)
- Formal KG-Microbe taxon universe: 2,093 taxa with at least one core oxygen-tolerance annotation among facultatively aerobic, aerobic, microaerophilic
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
| facultatively aerobic | 0.219 |
| aerobic | 0.172 |
| microaerophilic | 0.093 |

## Target Taxa Most Often Reported In Increased Profiles

| taxon | increased / decreased profiles | increased-profile tendency |
| --- | ---: | ---: |
| Eggerthella lenta | 20 / 1 | 0.905 |
| Fusobacterium nucleatum | 16 / 3 | 0.684 |
| Parvimonas micra | 14 / 3 | 0.647 |
| Escherichia coli | 30 / 8 | 0.579 |
| Prevotella intermedia | 14 / 4 | 0.556 |
| Streptococcus anginosus | 15 / 5 | 0.500 |
| Veillonella dispar | 18 / 7 | 0.440 |

## Recurrent Target-Taxon Increased Profiles By Condition

| taxon | condition | body site | increased / decreased profiles | increased studies |
| --- | --- | --- | ---: | ---: |
| Escherichia coli | Colorectal cancer | Feces | 4 / 0 | 4 |
| Fusobacterium nucleatum | Colorectal cancer | Feces | 4 / 0 | 4 |
| Parvimonas micra | Colorectal cancer | Feces | 3 / 0 | 3 |
| Klebsiella pneumoniae | Parkinson's disease | Feces | 3 / 0 | 2 |
| Streptococcus anginosus | Rheumatoid arthritis | Subgingival dental plaque | 3 / 0 | 2 |

Target-taxon occurrence views (`target_taxa_condition_breakdown_bar.png`, `condition_target_taxa_stacked_bar.png`, `candidate_taxa_profile_occurrences.tsv` and `candidate_taxa_condition_direction.tsv`) count human-host 16S and shotgun profiles only, the corpus of the profile metadata. 82 occurrences (19 decreased, 63 increased) in 70 profiles from 24 studies lie outside the corpus: 72 come from experiments whose host is not recorded as Homo sapiens (Mus musculus 56, Rattus norvegicus 7, Sus scrofa domesticus 4, Bos taurus 2, Canis lupus familiaris 2, Elephas maximus 1), and 10 from Homo sapiens experiments with a sequencing method other than 16S or WMS, or none recorded (not recorded 7, PCR 3). They have no profile metadata, so they are not counted under any condition; `candidate_taxa_excluded_host_occurrences.tsv` lists them.

## Contexts With Aerobe v2 Trait Enrichment in Increased Profiles

| condition | profiles | studies | increased - decreased z | paired d | paired p | matched pairs | increased mean trait z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Multiple myeloma | 8 | 2 | 1.207 |  |  | 4 | 0.236 |
| Diarrhea | 18 | 3 | 1.041 | 1.450 | 0.00245 | 9 | 0.729 |
| Esophageal squamous cell carcinoma | 6 | 2 | 1.033 |  |  | 3 | 0.390 |
| Chronic kidney disease | 20 | 9 | 0.864 | 0.929 | 0.0166 | 10 | 0.333 |
| Polycystic ovary syndrome | 34 | 9 | 0.834 | 0.741 | 0.00758 | 17 | 0.075 |
| Gastric adenocarcinoma | 10 | 2 | 0.689 | 1.066 | 0.0756 | 5 | 0.548 |
| Chronic obstructive pulmonary disease | 6 | 2 | 0.671 |  |  | 3 | 0.780 |
| Autism | 8 | 4 | 0.656 |  |  | 4 | 0.475 |

Paired d and paired p are blank for conditions with fewer than 5 matched increased/decreased profile pairs, or with zero variance among the paired differences.

Broad labels such as `Exposure`, `Population`, `Time`, `Timepoint`, and labels containing `measurement` come from BugSigDB condition metadata but usually describe study design, cohort, sampling time, exposure, or measurement context rather than a disease/condition context. They are excluded from displayed condition bar charts and condition summaries.

Oral-site increased-profile trait signals are expected in this analysis because the oral cavity is oxygen-exposed and several focal taxa are common oral or airway-associated aerotolerant, facultatively aerobic, or microaerophilic lineages. These signals should be interpreted as body-site-specific oral microbial ecology, not as a gut increased-profile trait-enrichment signal.

## Data-Inferred Guild Evidence

For the top condition signals, this table uses all taxa observed in the increased profiles for that condition that also belong to the KG-Microbe-wide candidate taxon universe, then reports the profile v2 traits that are enriched in those same increased profiles and the universe taxa that carry those trait annotations. This is a data-derived guild summary, not a literature-defined aerobe/aerotolerant guild.

| condition | increased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Autism | 4 | 17 / 27 | aerobic (+1.24); facultatively aerobic (+0.38); microaerophilic (+0.35) | Corynebacteriaceae (3 traits, 2 profiles); Bacteroidota (3 traits, 2 profiles); Pseudomonadota (3 traits, 1 profile); Erysipelotrichaceae (3 traits, 1 profile); Lactobacillus (3 traits, 1 profile); Corynebacterium (3 traits, 1 profile); Bacillaceae (3 traits, 1 profile); Bacteroidia (3 traits, 1 profile); Bifidobacteriaceae (3 traits, 1 profile); Enterobacteriaceae (3 traits, 1 profile); Bacteroides (3 traits, 1 profile); Enterococcaceae (3 traits, 1 profile) |
| Chronic kidney disease | 10 | 93 / 153 | facultatively aerobic (+1.06); aerobic (+0.82); microaerophilic (+0.72) | Streptococcus (3 traits, 4 profiles); Gammaproteobacteria (3 traits, 3 profiles); Actinomycetota (3 traits, 3 profiles); Enterobacteriaceae (3 traits, 3 profiles); Enterobacterales (3 traits, 3 profiles); Pseudomonadota (3 traits, 2 profiles); Erysipelotrichaceae (3 traits, 2 profiles); Streptococcaceae (3 traits, 2 profiles); Enterococcus (3 traits, 2 profiles); Lactobacillus (3 traits, 2 profiles); Actinomyces (3 traits, 2 profiles); Corynebacterium (3 traits, 2 profiles) |
| Chronic obstructive pulmonary disease | 3 | 44 / 69 | microaerophilic (+1.18); facultatively aerobic (+1.12) | Actinomyces (2 traits, 3 profiles); Actinomycetaceae (2 traits, 3 profiles); Micrococcaceae (2 traits, 2 profiles); Lactobacillales (2 traits, 2 profiles); Rothia (2 traits, 2 profiles); Bacillota (2 traits, 2 profiles); Enterococcus (2 traits, 2 profiles); Actinomycetales (2 traits, 2 profiles); Enterococcaceae (2 traits, 2 profiles); Capnocytophaga (2 traits, 1 profile); Leptotrichiaceae (2 traits, 1 profile); Erysipelotrichaceae (2 traits, 1 profile) |
| Diarrhea | 9 | 14 / 24 | facultatively aerobic (+1.54); microaerophilic (+0.94); aerobic (+0.64) | Streptococcus (3 traits, 6 profiles); Enterococcus (3 traits, 4 profiles); Escherichia (3 traits, 4 profiles); Campylobacter (3 traits, 3 profiles); Haemophilus parainfluenzae (3 traits, 2 profiles); Bifidobacterium (3 traits, 1 profile); Bacteroidales (3 traits, 1 profile); Bacteroidia (3 traits, 1 profile); Betaproteobacteria (3 traits, 1 profile); Enterobacteriaceae (3 traits, 1 profile); Burkholderiales (3 traits, 1 profile); Ligilactobacillus ruminis (2 traits, 2 profiles) |
| Esophageal squamous cell carcinoma | 3 | 12 / 16 | facultatively aerobic (+1.54); microaerophilic (+1.46) | Actinomyces (2 traits, 2 profiles); Leptotrichiaceae (2 traits, 1 profile); Fusobacteriia (2 traits, 1 profile); Fusobacteriales (2 traits, 1 profile); Actinomycetaceae (2 traits, 1 profile); Fusobacteriota (2 traits, 1 profile); Leptotrichia (2 traits, 1 profile); Cardiobacterium (1 trait, 1 profile); Parvimonas (1 trait, 1 profile); Coriobacteriia (1 trait, 1 profile) |
| Gastric adenocarcinoma | 5 | 21 / 33 | facultatively aerobic (+0.93); aerobic (+0.65); microaerophilic (+0.48) | Enterobacteriaceae (3 traits, 4 profiles); Bifidobacterium (3 traits, 3 profiles); Bifidobacteriaceae (3 traits, 3 profiles); Lactobacillaceae (3 traits, 3 profiles); Gammaproteobacteria (3 traits, 2 profiles); Streptococcaceae (3 traits, 2 profiles); Streptococcus (3 traits, 2 profiles); Lactobacillus (3 traits, 2 profiles); Actinomycetes (3 traits, 2 profiles); Lactobacillales (3 traits, 2 profiles); Actinomycetota (3 traits, 2 profiles); Escherichia (3 traits, 2 profiles) |
| Multiple myeloma | 4 | 19 / 35 | facultatively aerobic (+1.84); aerobic (+0.95) | Streptococcus gordonii (2 traits, 2 profiles); Streptococcus oralis (2 traits, 2 profiles); Streptococcus salivarius (2 traits, 2 profiles); Streptococcus pneumoniae (2 traits, 2 profiles); Klebsiella variicola (2 traits, 2 profiles); Streptococcus mitis (2 traits, 2 profiles); Klebsiella ornithinolytica (2 traits, 2 profiles); Citrobacter freundii (2 traits, 2 profiles); Klebsiella aerogenes (2 traits, 2 profiles); Enterobacter cloacae (2 traits, 2 profiles); Klebsiella pneumoniae (2 traits, 2 profiles); Pseudomonadota (2 traits, 1 profile) |
| Polycystic ovary syndrome | 17 | 76 / 161 | facultatively aerobic (+0.99); aerobic (+0.80) | Enterococcus (2 traits, 3 profiles); Actinomycetota (2 traits, 3 profiles); Enterobacteriaceae (2 traits, 3 profiles); Citrobacter (2 traits, 3 profiles); Bacteroidaceae (2 traits, 3 profiles); Bacteroides (2 traits, 3 profiles); Gammaproteobacteria (2 traits, 2 profiles); Comamonas (2 traits, 2 profiles); Pseudomonas (2 traits, 2 profiles); Enterobacter (2 traits, 2 profiles); Acidimicrobiia (2 traits, 2 profiles); Bifidobacteriaceae (2 traits, 2 profiles) |

## Body-Site Check

| body_site_group | profiles | increased mean trait z | increased - decreased z |
| --- | ---: | ---: | ---: |
| Oral cavity | 24 | 0.901 | 0.276 |
| Saliva | 136 | 0.369 | -0.023 |
| Vagina | 38 | 0.334 | -0.151 |
| Nasopharynx | 22 | 0.284 | -0.287 |
| Other | 364 | 0.187 | 0.034 |
| Subgingival dental plaque | 40 | 0.083 | -0.552 |
| Uterine cervix | 32 | 0.075 | -0.544 |
| Throat | 24 | 0.011 | -0.528 |
| Feces | 1190 | -0.012 | 0.294 |

## Top Increased Profiles To Inspect

| profile | condition | body site | trait z | title |
| --- | --- | --- | ---: | --- |
| bsdb:39642873\|bsdb:39642873-1\|inc | Respiratory Syncytial Virus Infection | Nasopharynx | 1.389 | The respiratory microbiome is linked to the severity of RSV infections and the persiste... |
| bsdb:39642873\|bsdb:39642873-4\|inc | Respiratory Syncytial Virus Infection | Nasopharynx | 1.389 | The respiratory microbiome is linked to the severity of RSV infections and the persiste... |
| bsdb:39999841\|bsdb:39999841-2\|inc | Nasopharyngeal squamous cell carcinoma | Saliva | 1.370 | Unveiling familial aggregation of nasopharyngeal carcinoma: Insights from oral microbio... |
| bsdb:36203175\|bsdb:36203175-1\|inc | Down syndrome | Oral cavity | 1.367 | Salivary microbiome in children with Down syndrome: a case-control study |
| bsdb:40155303\|bsdb:40155303-1\|inc | Pruritus | Feces | 1.351 | Gut microbial signatures associated with uremic pruritus in hemodialysis patients |
| bsdb:29313867\|bsdb:29313867-1\|inc | Chlamydia trachomatis infectious disease | Endocervix | 1.247 | HPV/Chlamydia trachomatis co-infection: metagenomic analysis of cervical microbiota in ... |
| bsdb:23459324\|bsdb:23459324-3\|inc | Obesity | Feces | 1.182 | Correlation between body mass index and gut concentrations of Lactobacillus reuteri, Bi... |
| bsdb:30463989\|bsdb:30463989-1\|inc | Human papilloma virus infection | Uterine cervix | 1.176 | Evaluation of the Associations Between Cervical Microbiota and HPV Infection, Clearance... |
| bsdb:28038683\|bsdb:28038683-7\|inc | Colorectal adenoma | Feces | 1.159 | The gut microbiota in conventional and serrated precursors of colorectal cancer |
| bsdb:38745393\|bsdb:38745393-5\|inc | Phenylketonuria | Saliva | 1.144 | Oral microbiota of patients with phenylketonuria: A nation-based cross-sectional study |

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
