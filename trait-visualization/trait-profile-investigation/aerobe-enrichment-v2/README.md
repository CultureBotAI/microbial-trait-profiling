# Aerobe v2 Trait Enrichment in Increased Profiles

This investigation uses the existing study trait profile clustered heatmap matrix, study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.

## Scope

- Selected traits: aerobic (METPO:1000602), microaerophilic (METPO:1000604)
- Target taxa: Parvimonas micra, Klebsiella pneumoniae, Streptococcus anginosus, Escherichia coli, Fusobacterium nucleatum, Enterocloster bolteae, Eggerthella lenta, Prevotella intermedia, Porphyromonas gingivalis, Veillonella dispar (the 9 aerobe-v2 discovery taxa in `candidate_taxa_direction_summary.tsv` plus Porphyromonas gingivalis NCBITaxon:837 included as an anaerobic comparator in the BugSigDB occurrence breakdown charts and tables)
- Formal KG-Microbe taxon universe: 1,978 taxa with at least one core oxygen-tolerance annotation among aerobic, microaerophilic
- Trait selection audit: `guild_trait_selection.tsv` records the selected trait IDs, short rationale, and KG-Microbe taxon counts for each selected trait.
- Context: all body sites, with body-site and condition summaries to separate gut, oral, airway, and other patterns.

## Plot Guild Specificity

| plot | guild-specific? | basis |
| --- | --- | --- |
| `condition_increased_guild_signal.png` | yes | v2 trait-enrichment score by condition |
| `guild_trait_direction_shift.png` | yes | the two selected v2 traits |
| `anatomical_sites_diagram.png` | yes | v2 trait-enrichment score by body-site group |
| `candidate_taxa_increased_profile_tendency.png` | no | fixed target-taxon direction counts |
| `condition_target_taxa_stacked_bar.png` | no | fixed target-taxon occurrence counts stacked by taxon |
| `target_taxa_condition_breakdown_bar.png` | no | fixed target-taxon occurrence counts |

## Main Signal

The selected traits are enriched in increased profiles relative to decreased profiles, supporting an aerobe v2 increased-profile trait-enrichment signal.

| trait | increased - decreased z |
| --- | ---: |
| microaerophilic | 0.124 |
| aerobic | 0.106 |

## Target Taxa Most Often Reported In Increased Profiles

| taxon | increased / decreased profiles | increased-profile tendency |
| --- | ---: | ---: |
| Eggerthella lenta | 36 / 8 | 0.636 |
| Fusobacterium nucleatum | 46 / 11 | 0.614 |
| Parvimonas micra | 32 / 9 | 0.561 |
| Klebsiella pneumoniae | 39 / 12 | 0.529 |
| Streptococcus anginosus | 31 / 14 | 0.378 |
| Escherichia coli | 78 / 36 | 0.368 |
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
| Veillonella dispar | Type II diabetes mellitus | Feces | 4 / 0 | 2 |
| Streptococcus anginosus | Rheumatoid arthritis | Subgingival dental plaque | 3 / 0 | 2 |
| Escherichia coli | Obesity | Feces | 4 / 10 | 2 |

Target-taxon occurrence views (`target_taxa_condition_breakdown_bar.png`, `condition_target_taxa_stacked_bar.png`, `candidate_taxa_profile_occurrences.tsv` and `candidate_taxa_condition_direction.tsv`) count human-host 16S and shotgun profiles only, the corpus of the profile metadata. 82 occurrences (19 decreased, 63 increased) in 70 profiles from 24 studies lie outside the corpus: 72 come from experiments whose host is not recorded as Homo sapiens (Mus musculus 56, Rattus norvegicus 7, Sus scrofa domesticus 4, Bos taurus 2, Canis lupus familiaris 2, Elephas maximus 1), and 10 from Homo sapiens experiments with a sequencing method other than 16S or WMS, or none recorded (not recorded 7, PCR 3). They have no profile metadata, so they are not counted under any condition; `candidate_taxa_excluded_host_occurrences.tsv` lists them.

## Contexts With Aerobe v2 Trait Enrichment in Increased Profiles

| condition | profiles | studies | increased - decreased z | paired d | nominal paired p | matched pairs | increased mean trait z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Response to antibiotic | 9 | 3 | 1.178 |  |  | 3 | 0.819 |
| Response to transplant | 26 | 5 | 0.973 | 0.998 | 0.0079 | 11 | 0.100 |
| Multiple myeloma | 9 | 2 | 0.944 |  |  | 4 | 0.175 |
| Graft versus host disease | 14 | 5 | 0.869 |  |  | 4 | 0.035 |
| Age-related macular degeneration | 22 | 5 | 0.817 | 0.950 | 0.0103 | 11 | 0.323 |
| Urinary tract infection | 13 | 4 | 0.752 |  |  | 4 | 0.160 |
| Cognitive impairment | 15 | 4 | 0.744 | 0.709 | 0.188 | 5 | 0.408 |
| Epilepsy | 16 | 5 | 0.731 | 1.100 | 0.0431 | 6 | 0.192 |

Reported paired p values are uncorrected two-sided paired t-test nominal p-values across 99 tested conditions after excluding broad study-design labels in `condition_guild_summary.tsv` (shown above: top 8 eligible conditions by guild-z difference, 4 of them tested). Paired d and nominal paired p are blank for conditions with fewer than 5 matched increased/decreased profile pairs, or with zero variance among the paired differences.

Broad labels such as `Exposure`, `Population`, `Time`, `Timepoint`, and labels containing `measurement` come from BugSigDB condition metadata but usually describe study design, cohort, sampling time, exposure, or measurement context rather than a disease/condition context. They are excluded from displayed condition bar charts and condition summaries.

Oral-site increased-profile trait signals are expected in this analysis because the oral cavity is oxygen-exposed and several focal taxa are common oral or airway-associated aerotolerant, facultatively aerobic, or microaerophilic lineages. These signals should be interpreted as body-site-specific oral microbial ecology, not as a gut increased-profile trait-enrichment signal.

## Data-Inferred Guild Evidence

For the top condition signals, this table uses all taxa observed in the increased profiles for that condition that also belong to the KG-Microbe-wide candidate taxon universe, then reports the profile v2 traits that are enriched in those same increased profiles and the universe taxa that carry those trait annotations. This is a data-derived guild summary, not a literature-defined aerobe/aerotolerant guild.

| condition | increased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Age-related macular degeneration | 11 | 28 / 42 | microaerophilic (+1.09) | Streptococcus (1 trait, 4 profiles); Gemella (1 trait, 4 profiles); Actinomyces sp. (1 trait, 2 profiles); Pseudomonadota (1 trait, 1 profile); Gammaproteobacteria (1 trait, 1 profile); Bacillota (1 trait, 1 profile); Streptococcaceae (1 trait, 1 profile); Streptococcus sp. (1 trait, 1 profile); Enterococcus (1 trait, 1 profile); Pseudomonadaceae (1 trait, 1 profile); Lactobacillus (1 trait, 1 profile); Peptoniphilus (1 trait, 1 profile) |
| Cognitive impairment | 8 | 10 / 25 | microaerophilic (+1.13); aerobic (+0.36) | Streptococcus (2 traits, 2 profiles); Lactobacillus (2 traits, 2 profiles); Lactobacillaceae (2 traits, 2 profiles); Pseudomonadota (2 traits, 1 profile); Streptococcaceae (2 traits, 1 profile); Lactococcus (2 traits, 1 profile); Lactobacillales (2 traits, 1 profile); Enterobacteriaceae (2 traits, 1 profile); Prevotellaceae (1 trait, 1 profile); Verrucomicrobiota (1 trait, 1 profile) |
| Epilepsy | 10 | 49 / 85 | microaerophilic (+0.90); aerobic (+0.56) | Lactobacillales (2 traits, 2 profiles); Enterobacteriaceae (2 traits, 2 profiles); Escherichia (2 traits, 2 profiles); Negativicutes (2 traits, 2 profiles); Bacilli (2 traits, 2 profiles); Enterobacterales (2 traits, 2 profiles); Veillonellales (2 traits, 2 profiles); Actinomycetota (2 traits, 2 profiles); Veillonellaceae (2 traits, 2 profiles); Pseudomonadota (2 traits, 1 profile); Gammaproteobacteria (2 traits, 1 profile); Bacillota (2 traits, 1 profile) |
| Graft versus host disease | 7 | 11 / 15 | aerobic (+1.31) | Bacillota (1 trait, 2 profiles); Stenotrophomonas (1 trait, 2 profiles); Enterobacteriaceae (1 trait, 2 profiles); Pseudomonadota (1 trait, 1 profile); Gammaproteobacteria (1 trait, 1 profile); Lysobacterales (1 trait, 1 profile); Lactococcus (1 trait, 1 profile); Porphyromonadaceae (1 trait, 1 profile); Lysobacteraceae (1 trait, 1 profile); Lautropia (1 trait, 1 profile); Enterobacterales (1 trait, 1 profile) |
| Multiple myeloma | 5 | 19 / 37 | aerobic (+1.19) | Streptococcus gordonii (1 trait, 2 profiles); Streptococcus oralis (1 trait, 2 profiles); Streptococcus salivarius (1 trait, 2 profiles); Streptococcus pneumoniae (1 trait, 2 profiles); Klebsiella variicola (1 trait, 2 profiles); Streptococcus mitis (1 trait, 2 profiles); Klebsiella ornithinolytica (1 trait, 2 profiles); Citrobacter freundii (1 trait, 2 profiles); Klebsiella aerogenes (1 trait, 2 profiles); Enterobacter cloacae (1 trait, 2 profiles); Klebsiella pneumoniae (1 trait, 2 profiles); Pseudomonadota (1 trait, 1 profile) |
| Response to antibiotic | 5 | 41 / 65 | microaerophilic (+1.25); aerobic (+1.11) | Moraxella (2 traits, 2 profiles); Enterobacter (2 traits, 2 profiles); Haemophilus (2 traits, 2 profiles); Streptococcus (2 traits, 1 profile); Enterococcus (2 traits, 1 profile); Lactococcus (2 traits, 1 profile); Clostridium (2 traits, 1 profile); Clostridia (2 traits, 1 profile); Eubacteriales (2 traits, 1 profile); Desulfovibrionaceae (2 traits, 1 profile); Neisseria subflava (2 traits, 1 profile); Clostridiaceae (2 traits, 1 profile) |
| Response to transplant | 15 | 74 / 201 | microaerophilic (+1.18) | Streptococcus (1 trait, 3 profiles); Lactobacillus (1 trait, 3 profiles); Lactobacillales (1 trait, 2 profiles); Helicobacter (1 trait, 2 profiles); Enterococcus hermanniensis (1 trait, 2 profiles); Lactobacillaceae (1 trait, 2 profiles); Neisseria (1 trait, 2 profiles); Capnocytophaga (1 trait, 1 profile); Granulicatella (1 trait, 1 profile); Sphingobacteriia (1 trait, 1 profile); Gammaproteobacteria (1 trait, 1 profile); Erysipelotrichaceae (1 trait, 1 profile) |
| Urinary tract infection | 7 | 47 / 62 | aerobic (+0.91) | Klebsiella pneumoniae (1 trait, 2 profiles); Burkholderiaceae (1 trait, 1 profile); Pseudomonadota (1 trait, 1 profile); Staphylococcus epidermidis (1 trait, 1 profile); Staphylococcus hominis (1 trait, 1 profile); Streptococcus (1 trait, 1 profile); Enterococcus faecalis (1 trait, 1 profile); Aerococcus (1 trait, 1 profile); Aerococcus urinae (1 trait, 1 profile); Lactobacillus (1 trait, 1 profile); Corynebacteriaceae (1 trait, 1 profile); Prevotellaceae (1 trait, 1 profile) |

## Body-Site Check

| body_site_group | profiles | increased mean trait z | increased - decreased z |
| --- | ---: | ---: | ---: |
| Oral cavity | 70 | 0.518 | -0.253 |
| Mouth | 76 | 0.427 | -0.118 |
| Saliva | 370 | 0.384 | -0.121 |
| Throat | 64 | 0.275 | -0.013 |
| Other | 1267 | 0.272 | -0.035 |
| Uterine cervix | 81 | 0.245 | 0.017 |
| Subgingival dental plaque | 75 | 0.137 | -0.426 |
| Vagina | 133 | -0.076 | -0.224 |
| Feces | 3493 | -0.087 | 0.219 |

## Top Increased Profiles To Inspect

| profile | condition | body site | trait z | title |
| --- | --- | --- | ---: | --- |
| bsdb:23631345\|bsdb:23631345-1\|inc | Obesity | Feces | 2.272 | Differences in gut microbiota composition between obese and lean children: a cross-sect... |
| bsdb:26551842\|bsdb:26551842-4\|inc | Obesity | Feces | 2.272 | Correlation between body mass index and faecal microbiota from children |
| bsdb:32115690\|bsdb:32115690-1\|inc | Acute lymphoblastic leukemia | Feces | 2.044 | Gastrointestinal toxicity during induction treatment for childhood acute lymphoblastic ... |
| bsdb:32332850\|bsdb:32332850-4\|inc | Cervical glandular intraepithelial neoplasia | Uterine cervix | 1.935 | The vaginal microbiota associates with the regression of untreated cervical intraepithe... |
| bsdb:32719372\|bsdb:32719372-8\|inc | Sampling site | Vagina | 1.935 | Changes in the vaginal microbiota across a gradient of urbanization |
| bsdb:36203175\|bsdb:36203175-1\|inc | Down syndrome | Oral cavity | 1.725 | Salivary microbiome in children with Down syndrome: a case-control study |
| bsdb:37596536\|bsdb:37596536-10\|inc | Population | Saliva | 1.699 | Oral microbiome homogeneity across diverse human groups from southern Africa: first res... |
| bsdb:39642873\|bsdb:39642873-3\|inc | Respiratory Syncytial Virus Infection | Nasopharynx | 1.679 | The respiratory microbiome is linked to the severity of RSV infections and the persiste... |
| bsdb:39642873\|bsdb:39642873-2\|inc | Respiratory Syncytial Virus Infection | Nasopharynx | 1.679 | The respiratory microbiome is linked to the severity of RSV infections and the persiste... |
| bsdb:36803868\|bsdb:36803868-5\|inc | Response to antibiotic | Sputum | 1.679 | Sputum bacterial load and bacterial composition correlate with lung function and are al... |

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
