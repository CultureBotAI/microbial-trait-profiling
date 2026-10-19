# Enrichment Infographic Follow-Up

Source figure: https://github.com/CultureBotAI/bsdbng/blob/enrichment-analysis-heatmap-cluster/data/derived/enrichment_infographic.png

The target conditions for this run are selected from the infographic top three condition-level outliers. The local follow-up uses the full `data/derived/study_trait_profiles_wide.tsv` matrix to list titles, find similar profiles, and connect enriched traits to organisms. Profiles carry their own experiment's condition, body site and host from `trait-visualization/study_trait_profile_study_metadata.tsv`; experiment details and reported taxa come from the pinned `data/raw/full_dump.csv`, and only profiles in this corpus are counted.

## Target Conditions

| condition | infographic avg enrichment depth | local profiles | studies | experiments | local mean nonzero traits/profile | local percentile | local outlier? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Sleep Disorder | 216.00 | 7 | 2 | 5 | 107.0 | 48.9% | no |
| Head and neck squamous cell carcinoma | 210.33 | 25 | 4 | 13 | 60.0 | 21.6% | no |

Target conditions with fewer than 2 profiles in this corpus are not summarized or queried: Thalassemia (0 profiles).

Interpretation: condition support matters because profile counts, study counts, and body-site mix vary across exact BugSigDB condition labels.

## Study And Experiment Titles

| condition | study | year | experiment | profiles |
| --- | --- | ---: | --- | ---: |
| Head and neck squamous cell carcinoma | Microbiomic differences in tumor and paired-normal tissue in head and neck squamous cell carcinomas (bsdb:28173873) | 2017 | Experiment 1: normal adjacent tissues vs squamous cell carcinoma tissues | 2 |
| Head and neck squamous cell carcinoma | Microbiomic differences in tumor and paired-normal tissue in head and neck squamous cell carcinomas (bsdb:28173873) | 2017 | Experiment 2: T0-T2 low tumor stages patients vs T3-T4 high tumor stages patients | 2 |
| Head and neck squamous cell carcinoma | Association of Oral Microbiome With Risk for Incident Head and Neck Squamous Cell Cancer (bsdb:29327043) | 2018 | Experiment 1: controls vs HNSCC | 2 |
| Head and neck squamous cell carcinoma | Gut microbiome predicts gastrointestinal toxicity outcomes from chemoradiation therapy in patients with head and neck squamous cell carcinoma (bsdb:38006691) | 2023 | Experiment 1: Oral mucositis CTCAE grade 1-2 (mild) vs Oral mucositis CTCAE grade 3 (severe) | 2 |
| Head and neck squamous cell carcinoma | Gut microbiome predicts gastrointestinal toxicity outcomes from chemoradiation therapy in patients with head and neck squamous cell carcinoma (bsdb:38006691) | 2023 | Experiment 2: No enteral feeding plus mild mucositis vs Enteral feeding plus mucositis of any grade | 2 |
| Head and neck squamous cell carcinoma | Gut microbiome predicts gastrointestinal toxicity outcomes from chemoradiation therapy in patients with head and neck squamous cell carcinoma (bsdb:38006691) | 2023 | Experiment 3: No Enteral feeding tube placement vs Enteral feeding tube placement | 2 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 1: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 2 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 2: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 2 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 3: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 2 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 4: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 2 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 5: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 1 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 6: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 2 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 7: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 2 |
| Sleep Disorder | Short-Chain Fatty Acid-Producing Gut Microbiota Is Decreased in Parkinson's Disease but Not in Rapid-Eye-Movement Sleep Behavior Disorder (bsdb:33293403) | 2020 | Experiment 1: Healthy Controls vs Rapid-eye-movement sleep behavior disorder (iRBD) patients | 1 |
| Sleep Disorder | Short-Chain Fatty Acid-Producing Gut Microbiota Is Decreased in Parkinson's Disease but Not in Rapid-Eye-Movement Sleep Behavior Disorder (bsdb:33293403) | 2020 | Experiment 2: Healthy Controls vs Rapid-eye-movement sleep behavior disorder (iRBD) patients | 1 |
| Sleep Disorder | Short-Chain Fatty Acid-Producing Gut Microbiota Is Decreased in Parkinson's Disease but Not in Rapid-Eye-Movement Sleep Behavior Disorder (bsdb:33293403) | 2020 | Experiment 3: Healthy Controls vs Rapid-eye-movement sleep behavior disorder (iRBD) patients | 1 |
| Sleep Disorder | Preliminary analysis of salivary microbiota in catathrenia (nocturnal groaning) using machine learning algorithms (bsdb:40247863) | 2025 | Experiment 1: Healthy control group vs Catathrenia group | 2 |
| Sleep Disorder | Preliminary analysis of salivary microbiota in catathrenia (nocturnal groaning) using machine learning algorithms (bsdb:40247863) | 2025 | Experiment 2: Pre treatment of Mandibular Advancement Devices (MAD) group vs Post treatment of Mandibular Advancement Devices (MAD) group | 2 |

## Autism Label Check

| exact label | EFO IDs | profiles | studies | experiments | interpretation |
| --- | --- | ---: | ---: | ---: | --- |
| Autism | EFO:0003758 | 14 | 4 | 10 | Distinct exact BugSigDB label; analyze separately unless an explicit autism/ASD grouping is declared. |
| Autism spectrum disorder | EFO:0003756 | 7 | 3 | 5 | Distinct exact BugSigDB label; analyze separately unless an explicit autism/ASD grouping is declared. |

Autism and autism spectrum disorder should not be silently treated as the same exact condition label. For a clinical/domain rollup, they can be grouped as `Autism / autism spectrum disorder` if the merge is declared, but exact label and EFO-level analyses should keep them separate.

## Similar Profiles

Profile similarity uses cosine similarity over log1p-transformed full non-MediaDive trait vectors. For each target condition, the query vector is the centroid of its active profiles, and the top five individual profiles exclude the same exact condition label.

| query condition | rank | neighbor condition | similarity | study title | experiment |
| --- | ---: | --- | ---: | --- | --- |
| Sleep Disorder | 1 | Schizophrenia | 0.935 | Alterations in fecal bacteriome virome interplay and microbiota-derived dysfunction in patients with schizophrenia | Experiment 1: Healthy control group (HC) vs Schizophrenia group (SCZ) |
| Sleep Disorder | 2 | Pancreatic carcinoma | 0.925 | Fecal microbiome signatures of pancreatic cancer patients | Experiment 4: Healthy control vs NBO PC |
| Sleep Disorder | 3 | Endometriosis | 0.924 | Associations Between Endometriosis and Gut Microbiota | Experiment 1: Healthy controls taken from a larger cohort study vs Endometriosis patients |
| Sleep Disorder | 4 | Endometriosis | 0.922 | Associations Between Endometriosis and Gut Microbiota | Experiment 2: Healthy controls taken from a larger cohort study vs Endometriosis patients |
| Sleep Disorder | 5 | Obesity | 0.922 | Metatranscriptomic analysis to define the Secrebiome, and 16S rRNA profiling of the gut microbiome in obesity and metabolic syndrome of Mexican children | Experiment 2: Normal weight (NW) vs Obesity with metabolic syndrome (OMS) |
| Head and neck squamous cell carcinoma | 1 | COVID-19 | 0.938 | Metagenomic assessment of gut microbial communities and risk of severe COVID-19 | Experiment 1: Moderate COVID-19 vs Severe COVID-19 |
| Head and neck squamous cell carcinoma | 2 | Multiple myeloma | 0.936 | Alterations of gut microbiome accelerate multiple myeloma progression by increasing the relative abundances of nitrogen-recycling bacteria | Experiment 1: healthy control vs multiple myeloma patients |
| Head and neck squamous cell carcinoma | 3 | Gingivitis | 0.935 | Quantitative analysis of the effects of essential oil mouthrinses on clinical plaque microbiome: a parallel-group, randomized trial | Experiment 1: Healthy subjects vs Gingivitis subjects |
| Head and neck squamous cell carcinoma | 4 | Psoriasis | 0.934 | Multiomic Analysis of the Gut Microbiome in Psoriasis Reveals Distinct Host‒Microbe Associations | Experiment 1: Healthy controls vs Psoriasis patients (PSO) |
| Head and neck squamous cell carcinoma | 5 | Depressive disorder | 0.933 | Gut Microbiota Alterations and Their Functional Differences in Depression According to Enterotypes in Asian Individuals | Experiment 2: Healthy group vs Depressive group |

### Cross-Target Check

| query condition | any other target in top 5 profiles? | matched target conditions |
| --- | --- | --- |
| Sleep Disorder | no | NA |
| Head and neck squamous cell carcinoma | no | NA |

Top five condition-centroid neighbors are also written to `target_condition_centroid_neighbors.tsv`; this checks whether the three conditions find each other at condition rather than individual-profile level.

| query condition | rank | centroid neighbor | similarity | other target? |
| --- | ---: | --- | ---: | --- |
| Sleep Disorder | 1 | Type I diabetes mellitus | 0.959 | no |
| Sleep Disorder | 2 | HIV-1 infection | 0.957 | no |
| Sleep Disorder | 3 | Response to antiviral drug | 0.956 | no |
| Sleep Disorder | 4 | Health study participation | 0.952 | no |
| Sleep Disorder | 5 | Milk allergic reaction | 0.947 | no |
| Head and neck squamous cell carcinoma | 1 | Response to diet | 0.980 | no |
| Head and neck squamous cell carcinoma | 2 | Lifestyle measurement | 0.974 | no |
| Head and neck squamous cell carcinoma | 3 | Crohn's disease | 0.972 | no |
| Head and neck squamous cell carcinoma | 4 | Waist circumference | 0.971 | no |
| Head and neck squamous cell carcinoma | 5 | Squamous cell carcinoma | 0.969 | no |

## How Trait Sets Connect To Profiles

A profile is one BugSigDB study/experiment/direction row. Its trait vector counts KG-Microbe trait annotations among organisms reported as increased or decreased in that experiment. A trait set is a named group of biologically related traits, such as `aerobe_v2` or `anaerobe_v2_depletion`. A profile connects to a trait set when the selected trait columns carry signal in that profile; organism support is then checked by intersecting organisms observed in the same profile with organism-trait annotations from `data/derived/taxa_trait_matrix.tsv`.

| trait set | condition | profiles | mean trait-set signal | mean relative trait-set signal |
| --- | --- | ---: | ---: | ---: |
| aerobe_v2 | Head and neck squamous cell carcinoma | 25 | 181.360 | 0.02076 |
| aerobe_v2 | Sleep Disorder | 7 | 564.571 | 0.01147 |
| anaerobe_v2_depletion | Head and neck squamous cell carcinoma | 25 | 45.640 | 0.03019 |
| anaerobe_v2_depletion | Sleep Disorder | 7 | 230.714 | 0.04637 |

## Companion Traits Vs Organisms Heatmap

`traits_vs_organisms_heatmap.png` plots 45 organisms by 34 traits. Rows are organisms observed in the target-condition BugSigDB profiles; columns are binary KG-Microbe trait annotations selected from the infographic top traits plus high-signal target-condition traits.

## Files

- `outlier_condition_study_experiment_titles.tsv`
- `outlier_condition_depth_summary.tsv`
- `autism_condition_label_check.tsv`
- `outlier_condition_profile_neighbors.tsv`
- `outlier_condition_neighbor_crosscheck.tsv`
- `target_condition_centroid_neighbors.tsv`
- `guild_profile_connection_summary.tsv`
- `traits_vs_organisms_heatmap_matrix.tsv`
- `traits_vs_organisms_taxa_metadata.tsv`
- `traits_vs_organisms_trait_metadata.tsv`
- `traits_vs_organisms_heatmap.png`
