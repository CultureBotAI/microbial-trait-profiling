# Enrichment Infographic Follow-Up

Source figure: https://github.com/CultureBotAI/bsdbng/blob/enrichment-analysis-heatmap-cluster/data/derived/enrichment_infographic.png

The referenced infographic ranks the following three conditions as the top condition-level enrichment-depth outliers. The local follow-up uses the full `data/derived/study_trait_profiles_wide.tsv` matrix and raw BugSigDB experiment metadata to list titles, find similar profiles, and connect enriched traits to organisms.

## Target Conditions

| condition | infographic avg enrichment depth | local profiles | studies | experiments | local mean nonzero traits/profile | local percentile | local outlier? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Thalassemia | 225.00 | 4 | 1 | 2 | 1822.2 | 100.0% | yes |
| Sleep Disorder | 216.00 | 7 | 2 | 5 | 249.3 | 54.7% | no |
| Head and neck squamous cell carcinoma | 210.33 | 25 | 4 | 13 | 116.8 | 27.1% | no |

Interpretation: the infographic itself identifies these as high average enrichment-depth conditions. In this local profile-trait-count proxy, condition support matters: Thalassemia is concentrated in one study and two experiments, Sleep Disorder spans two studies and five experiments, and head and neck squamous cell carcinoma spans more experiments but has a mix of saliva/oral-cavity/fecal comparisons.

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
| Thalassemia | Exploring alterations of gut/blood microbes in addressing iron overload-induced gut dysbiosis and cognitive impairment in thalassemia patients (bsdb:39438708) | 2024 | Experiment 1: Control vs Thalasemia | 2 |
| Thalassemia | Exploring alterations of gut/blood microbes in addressing iron overload-induced gut dysbiosis and cognitive impairment in thalassemia patients (bsdb:39438708) | 2024 | Experiment 2: Non-transfusion dependent thalassemia vs Transfusion dependent thalasemia | 2 |

## Autism Label Check

| exact label | EFO IDs | profiles | studies | experiments | interpretation |
| --- | --- | ---: | ---: | ---: | --- |
| Autism | EFO:0003758 | 14 | 4 | 10 | Distinct exact BugSigDB label; analyze separately unless an explicit autism/ASD grouping is declared. |
| Autism spectrum disorder | EFO:0003756 | 62 | 15 | 38 | Distinct exact BugSigDB label; analyze separately unless an explicit autism/ASD grouping is declared. |

Autism and autism spectrum disorder should not be silently treated as the same exact condition label. For a clinical/domain rollup, they can be grouped as `Autism / autism spectrum disorder` if the merge is declared, but exact label and EFO-level analyses should keep them separate.

## Similar Profiles

Profile similarity uses cosine similarity over log1p-transformed full trait vectors. For each target condition, the query vector is the centroid of its active profiles, and the top five individual profiles exclude the same exact condition label.

| query condition | rank | neighbor condition | similarity | study title | experiment |
| --- | ---: | --- | ---: | --- | --- |
| Thalassemia | 1 | Cognitive impairment | 0.976 | Exploring alterations of gut/blood microbes in addressing iron overload-induced gut dysbiosis and cognitive impairment in thalassemia patients | Experiment 3: Normal vs Mild cognitive impairment |
| Thalassemia | 2 | Colorectal adenoma | 0.967 | Human oral microbiome dysbiosis as a novel non-invasive biomarker in detection of colorectal cancer | Experiment 2: Healthy Individuals vs CRA |
| Thalassemia | 3 | Chronic obstructive pulmonary disease | 0.965 | Sputum Microbiome and Chronic Obstructive Pulmonary Disease in a Rural Ugandan Cohort of Well-Controlled HIV Infection | Experiment 1: HIV-negative individuals with COPD vs HIV-positive individuals diagnosed with COPD |
| Thalassemia | 4 | Parkinson's disease | 0.964 | Altered gut microbiota in Parkinson's disease patients/healthy spouses and its association with clinical features | Experiment 1: Healthy Controls vs Parkinson's Disease Patients |
| Thalassemia | 5 | Lifestyle measurement | 0.964 | Gut microbiome profiling of a rural and urban South African cohort reveals biomarkers of a population in lifestyle transition | Experiment 2: Lean female subjects in Bushbuckridge vs Lean female subjects in Soweto |
| Sleep Disorder | 1 | Aging | 0.915 | The microbiome is associated with obesity-related metabolome signature in the process of aging | Experiment 1: Young mice (8 weeks old) vs Aged mice (18 months old) |
| Sleep Disorder | 2 | Schizophrenia | 0.914 | Alterations in fecal bacteriome virome interplay and microbiota-derived dysfunction in patients with schizophrenia | Experiment 1: Healthy control group (HC) vs Schizophrenia group (SCZ) |
| Sleep Disorder | 3 | Pancreatic carcinoma | 0.913 | Fecal microbiome signatures of pancreatic cancer patients | Experiment 4: Healthy control vs NBO PC |
| Sleep Disorder | 4 | Air pollution | 0.910 | Effects of inhaled air pollution on markers of integrity, inflammation, and microbiota profiles of the intestines in Apolipoprotein E knockout mice | Experiment 2: mice exposed to Filtered Air vs mice exposed to WS |
| Sleep Disorder | 5 | Obesity | 0.910 | Metatranscriptomic analysis to define the Secrebiome, and 16S rRNA profiling of the gut microbiome in obesity and metabolic syndrome of Mexican children | Experiment 2: Normal weight (NW) vs Obesity with metabolic syndrome (OMS) |
| Head and neck squamous cell carcinoma | 1 | Cognitive impairment | 0.921 | Large-scale metagenomic analysis of oral microbiomes reveals markers for autism spectrum disorders | Experiment 5: Low Intelligence Quotients(Low IQ ≤70) by diet vs High Intelligence Quotients(High IQ > 85) by diet |
| Head and neck squamous cell carcinoma | 2 | Gingivitis | 0.920 | Quantitative analysis of the effects of essential oil mouthrinses on clinical plaque microbiome: a parallel-group, randomized trial | Experiment 1: Healthy subjects vs Gingivitis subjects |
| Head and neck squamous cell carcinoma | 3 | Response to transplant | 0.916 | Exclusive enteral nutrition initiates individual protective microbiome changes to induce remission in pediatric Crohn's disease | Experiment 4: Inactive group vs Active group |
| Head and neck squamous cell carcinoma | 4 | Cognitive impairment | 0.915 | Large-scale metagenomic analysis of oral microbiomes reveals markers for autism spectrum disorders | Experiment 4: Low Intelligence Quotients(Low IQ ≤70) vs High Intelligence Quotients(High IQ > 85) |
| Head and neck squamous cell carcinoma | 5 | COVID-19 | 0.913 | Gut microbiota composition reflects disease severity and dysfunctional immune responses in patients with COVID-19 | Experiment 3: Healthy controls vs Recovered COVID-19 cases |

### Cross-Target Check

| query condition | any other target in top 5 profiles? | matched target conditions |
| --- | --- | --- |
| Thalassemia | no |  |
| Sleep Disorder | no |  |
| Head and neck squamous cell carcinoma | no |  |

Top five condition-centroid neighbors are also written to `target_condition_centroid_neighbors.tsv`; this checks whether the three conditions find each other at condition rather than individual-profile level.

| query condition | rank | centroid neighbor | similarity | other target? |
| --- | ---: | --- | ---: | --- |
| Thalassemia | 1 | Seasonal gut microbiome measurement | 0.962 | no |
| Thalassemia | 2 | Schistosomiasis | 0.959 | no |
| Thalassemia | 3 | Gastric carcinoma | 0.957 | no |
| Thalassemia | 4 | Carrier status | 0.953 | no |
| Thalassemia | 5 | Renal transplant outcome measurement | 0.953 | no |
| Sleep Disorder | 1 | Non-alcoholic steatohepatitis | 0.950 | no |
| Sleep Disorder | 2 | Behcet's syndrome | 0.944 | no |
| Sleep Disorder | 3 | Type I diabetes mellitus | 0.941 | no |
| Sleep Disorder | 4 | Peripheral neuron | 0.941 | no |
| Sleep Disorder | 5 | Urolithiasis,Nephrolithiasis | 0.937 | no |
| Head and neck squamous cell carcinoma | 1 | Response to antibiotic | 0.970 | no |
| Head and neck squamous cell carcinoma | 2 | Periodontitis | 0.969 | no |
| Head and neck squamous cell carcinoma | 3 | Mycobacterium tuberculosis | 0.967 | no |
| Head and neck squamous cell carcinoma | 4 | Asthma | 0.966 | no |
| Head and neck squamous cell carcinoma | 5 | Cervical glandular intraepithelial neoplasia | 0.966 | no |

## How Trait Sets Connect To Profiles

A profile is one BugSigDB study/experiment/direction row. Its trait vector counts KG-Microbe trait annotations among organisms reported as increased or decreased in that experiment. A trait set is a named group of biologically related traits, such as `aerobe_v2` or `anaerobe_v2_depletion`. A profile connects to a trait set when the selected trait columns carry signal in that profile; organism support is then checked by intersecting organisms observed in the same profile with organism-trait annotations from `data/derived/taxa_trait_matrix.tsv`.

| trait set | condition | profiles | mean trait-set signal | mean relative trait-set signal |
| --- | --- | ---: | ---: | ---: |
| aerobe_v2 | Head and neck squamous cell carcinoma | 25 | 181.360 | 0.01820 |
| aerobe_v2 | Sleep Disorder | 7 | 564.571 | 0.00985 |
| aerobe_v2 | Thalassemia | 4 | 9083.834 | 0.02288 |
| anaerobe_v2_depletion | Head and neck squamous cell carcinoma | 25 | 45.640 | 0.02535 |
| anaerobe_v2_depletion | Sleep Disorder | 7 | 230.714 | 0.03673 |
| anaerobe_v2_depletion | Thalassemia | 4 | 6618.375 | 0.01886 |

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
