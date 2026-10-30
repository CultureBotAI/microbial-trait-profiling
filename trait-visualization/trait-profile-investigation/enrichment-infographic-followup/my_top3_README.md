# Enrichment Infographic Follow-Up

Source figure: https://github.com/CultureBotAI/bsdbng/blob/enrichment-analysis-heatmap-cluster/data/derived/enrichment_infographic.png

The target conditions for this run are selected from the top 3 paired-profile condition contrast rows from `trait-visualization\trait-profile-investigation\trait_profile_paired_contrast_condition_summary.tsv`. The local follow-up uses the full `data/derived/study_trait_profiles_wide.tsv` matrix to list titles, find similar profiles, and connect enriched traits to organisms. Profiles carry their own experiment's condition, body site and host from `trait-visualization/study_trait_profile_study_metadata.tsv`; experiment details and reported taxa come from the pinned `data/raw/full_dump.csv`, and only profiles in this corpus are counted.

## Target Conditions

| condition | paired absolute contrast load | local profiles | studies | experiments | local mean nonzero traits/profile | local percentile | local outlier? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Multiple sclerosis | 0.87 | 22 | 6 | 11 | 113.0 | 30.8% | no |
| Eczema | 0.70 | 20 | 4 | 10 | 110.0 | 27.5% | no |
| Pancreatic carcinoma | 0.70 | 24 | 5 | 12 | 163.9 | 65.4% | no |

Interpretation: condition support matters because profile counts, study counts, and body-site mix vary across exact BugSigDB condition labels.

## Study And Experiment Titles

| condition | study | year | experiment | profiles |
| --- | --- | ---: | --- | ---: |
| Eczema | Microarray analysis reveals marked intestinal microbiota aberrancy in infants having eczema compared to healthy children in at-risk for atopic disease (bsdb:23339708) | 2013 | Experiment 1: healthy control vs child with eczema | 2 |
| Eczema | Whole metagenome profiling reveals skin microbiome-dependent susceptibility to atopic dermatitis flare (bsdb:27562258) | 2016 | Experiment 1: Healthy controls vs Atopic dermatitis case group | 2 |
| Eczema | Whole metagenome profiling reveals skin microbiome-dependent susceptibility to atopic dermatitis flare (bsdb:27562258) | 2016 | Experiment 3: Males (control) vs Males (case) | 2 |
| Eczema | Whole metagenome profiling reveals skin microbiome-dependent susceptibility to atopic dermatitis flare (bsdb:27562258) | 2016 | Experiment 4: Normal control vs AD cases with no reported cream usage | 2 |
| Eczema | Altered Gut Microbiota Composition Associated with Eczema in Infants (bsdb:27812181) | 2016 | Experiment 1: healthy control infants vs infants with eczema | 2 |
| Eczema | Variations in early gut microbiome are associated with childhood eczema (bsdb:30860574) | 2019 | Experiment 1: healthy control vs child with eczema | 2 |
| Eczema | Variations in early gut microbiome are associated with childhood eczema (bsdb:30860574) | 2019 | Experiment 2: healthy control vs child with eczema | 2 |
| Eczema | Variations in early gut microbiome are associated with childhood eczema (bsdb:30860574) | 2019 | Experiment 3: healthy control vs child with eczema | 2 |
| Eczema | Variations in early gut microbiome are associated with childhood eczema (bsdb:30860574) | 2019 | Experiment 4: healthy control vs child with eczema | 2 |
| Eczema | Variations in early gut microbiome are associated with childhood eczema (bsdb:30860574) | 2019 | Experiment 5: healthy control vs child with eczema | 2 |
| Multiple sclerosis | Dysbiosis in the Gut Microbiota of Patients with Multiple Sclerosis, with a Striking Depletion of Species Belonging to Clostridia XIVa and IV Clusters (bsdb:26367776) | 2015 | Experiment 1: Healthy subjects (controls) vs Multiple sclerosis (MS), | 2 |
| Multiple sclerosis | Dysbiosis in the Gut Microbiota of Patients with Multiple Sclerosis, with a Striking Depletion of Species Belonging to Clostridia XIVa and IV Clusters (bsdb:26367776) | 2015 | Experiment 2: Healthy subjects (controls) vs Multiple sclerosis (MS), | 2 |
| Multiple sclerosis | Gut microbiota in early pediatric multiple sclerosis: a case-control study (bsdb:27176462) | 2016 | Experiment 1: healthy controls vs pediatric multiple sclerosis | 2 |
| Multiple sclerosis | Multiple sclerosis patients have a distinct gut microbiota compared to healthy controls (bsdb:27346372) | 2016 | Experiment 1: Healthy controls vs Multiple sclerosis | 2 |
| Multiple sclerosis | Alterations of the human gut microbiome in multiple sclerosis (bsdb:27352007) | 2016 | Experiment 1: healthy controls vs multiple sclerosis | 2 |
| Multiple sclerosis | Alterations of the human gut microbiome in multiple sclerosis (bsdb:27352007) | 2016 | Experiment 2: controls vs multiple sclerosis patients | 2 |
| Multiple sclerosis | Alterations of the human gut microbiome in multiple sclerosis (bsdb:27352007) | 2016 | Experiment 3: controls vs untreated multiple sclerosis | 2 |
| Multiple sclerosis | Alterations of the human gut microbiome in multiple sclerosis (bsdb:27352007) | 2016 | Experiment 4: controls vs untreated multiple sclerosis | 2 |
| Multiple sclerosis | A comparative study of the gut microbiota in immune-mediated inflammatory diseases-does a common dysbiosis exist? (bsdb:30545401) | 2018 | Experiment 3: healthy controls vs multiple sceloris | 2 |
| Multiple sclerosis | Multiple sclerosis and gut microbiota: Lachnospiraceae from the ileum of MS twins trigger MS-like disease in germfree transgenic mice-An unbiased functional study (bsdb:40258140) | 2025 | Experiment 1: Unaffected twins vs Affected twins | 2 |
| Multiple sclerosis | Multiple sclerosis and gut microbiota: Lachnospiraceae from the ileum of MS twins trigger MS-like disease in germfree transgenic mice-An unbiased functional study (bsdb:40258140) | 2025 | Experiment 2: Unaffected twins vs Affected twins | 2 |
| Pancreatic carcinoma | Human oral microbiome and prospective risk for pancreatic cancer: a population-based nested case-control study (bsdb:27742762) | 2018 | Experiment 1: controls vs pancreatic cancer cases | 2 |
| Pancreatic carcinoma | Gut microbial profile analysis by MiSeq sequencing of pancreatic carcinoma patients in China (bsdb:29221120) | 2017 | Experiment 1: healthy controls vs pancreatic cancer | 2 |
| Pancreatic carcinoma | Fecal microbiome signatures of pancreatic cancer patients (bsdb:31727922) | 2019 | Experiment 1: Healthy control vs Pancreatic cancer patients | 2 |
| Pancreatic carcinoma | Fecal microbiome signatures of pancreatic cancer patients (bsdb:31727922) | 2019 | Experiment 4: Healthy control vs NBO PC | 2 |
| Pancreatic carcinoma | Microbiome Markers of Pancreatic Cancer Based on Bacteria-Derived Extracellular Vesicles Acquired from Blood Samples: A Retrospective Propensity Score Matching Analysis (bsdb:33805810) | 2021 | Experiment 1: Healthy controls vs Pancreatic cancer patients | 2 |
| Pancreatic carcinoma | Microbiome Markers of Pancreatic Cancer Based on Bacteria-Derived Extracellular Vesicles Acquired from Blood Samples: A Retrospective Propensity Score Matching Analysis (bsdb:33805810) | 2021 | Experiment 3: Healthy controls vs Pancreatic cancer patients | 2 |
| Pancreatic carcinoma | Microbiome Markers of Pancreatic Cancer Based on Bacteria-Derived Extracellular Vesicles Acquired from Blood Samples: A Retrospective Propensity Score Matching Analysis (bsdb:33805810) | 2021 | Experiment 4: Healthy controls vs Pancreatic cancer patients | 2 |
| Pancreatic carcinoma | Microbiome Markers of Pancreatic Cancer Based on Bacteria-Derived Extracellular Vesicles Acquired from Blood Samples: A Retrospective Propensity Score Matching Analysis (bsdb:33805810) | 2021 | Experiment 6: Healthy controls vs Pancreatic cancer patients | 2 |
| Pancreatic carcinoma | Microbiome Markers of Pancreatic Cancer Based on Bacteria-Derived Extracellular Vesicles Acquired from Blood Samples: A Retrospective Propensity Score Matching Analysis (bsdb:33805810) | 2021 | Experiment 7: Healthy controls vs Pancreatic cancer patients | 2 |
| Pancreatic carcinoma | Microbiome Markers of Pancreatic Cancer Based on Bacteria-Derived Extracellular Vesicles Acquired from Blood Samples: A Retrospective Propensity Score Matching Analysis (bsdb:33805810) | 2021 | Experiment 8: Healthy controls vs Pancreatic cancer patients | 2 |
| Pancreatic carcinoma | Oral and fecal microbiome alterations in pancreatic cancer: insights into potential diagnostic biomarkers (bsdb:41039266) | 2025 | Experiment 1: Healthy control group vs Pancreatic cancer patients | 2 |
| Pancreatic carcinoma | Oral and fecal microbiome alterations in pancreatic cancer: insights into potential diagnostic biomarkers (bsdb:41039266) | 2025 | Experiment 2: Healthy control group vs Pancreatic cancer patients | 2 |

## Autism Label Check

| exact label | EFO IDs | profiles | studies | experiments | interpretation |
| --- | --- | ---: | ---: | ---: | --- |
| Autism | EFO:0003758 | 8 | 4 | 4 | Distinct exact BugSigDB label; analyze separately unless an explicit autism/ASD grouping is declared. |
| Autism spectrum disorder | EFO:0003756 | 4 | 2 | 2 | Distinct exact BugSigDB label; analyze separately unless an explicit autism/ASD grouping is declared. |

Autism and autism spectrum disorder should not be silently treated as the same exact condition label. For a clinical/domain rollup, they can be grouped as `Autism / autism spectrum disorder` if the merge is declared, but exact label and EFO-level analyses should keep them separate.

## Similar Profiles

Profile similarity uses cosine similarity over log1p-transformed full non-MediaDive trait vectors. For each target condition, the query vector is the centroid of its active profiles, and the top five individual profiles exclude the same exact condition label.

| query condition | rank | neighbor condition | similarity | study title | experiment |
| --- | ---: | --- | ---: | --- | --- |
| Multiple sclerosis | 1 | Gestational diabetes | 0.930 | Gestational diabetes is associated with change in the gut microbiota composition in third trimester of pregnancy and postpartum | Experiment 4: Previous normoglycaemic  pregnancy vs Previous gestational diabetes mellitus (GDM) pregnancy |
| Multiple sclerosis | 2 | Chlamydia trachomatis | 0.930 | Cervicovaginal microbiome and natural history of Chlamydia trachomatis in adolescents and young women | Experiment 1: Chlamydia trichomonas Controls (t0) vs Chlamydia trichomonas Cases (t0) |
| Multiple sclerosis | 3 | Rheumatoid arthritis | 0.923 | Self-Balance of Intestinal Flora in Spouses of Patients With Rheumatoid Arthritis | Experiment 1: Healthy controls vs Rheumatoid arthritis patient |
| Multiple sclerosis | 4 | Endometriosis | 0.922 | The Endobiota Study: Comparison of Vaginal, Cervical and Gut Microbiota Between Women with Stage 3/4 Endometriosis and Healthy Controls | Experiment 1: Healthy controls vs Stage 3/4 endometriosis patients |
| Multiple sclerosis | 5 | Rheumatoid arthritis | 0.922 | Dysbiotic Subgingival Microbial Communities in Periodontally Healthy Patients With Rheumatoid Arthritis | Experiment 1: healthy controls vs Rheumatoid Arthritis |
| Eczema | 1 | Milk allergic reaction | 0.948 | Lactobacillus rhamnosus GG-supplemented formula expands butyrate-producing bacterial strains in food allergic infants | Experiment 1: healthy controls vs IgE mediated cow's milk allergy |
| Eczema | 2 | Colorectal cancer | 0.947 | Dysbiosis of human gut microbiome in young-onset colorectal cancer | Experiment 2: oControls vs oCRC |
| Eczema | 3 | Human papilloma virus infection | 0.945 | The Alterations of Vaginal Microbiome in HPV16 Infection as Identified by Shotgun Metagenomic Sequencing | Experiment 1: Healthy control vs HPV-16+ |
| Eczema | 4 | Multiple chemical sensitivity | 0.945 | Metagenomic gut microbiome analysis of Japanese patients with multiple chemical sensitivity/idiopathic environmental intolerance | Experiment 1: Healthy Control (HC) vs Multiple Chemical Sensitivity (MCS) |
| Eczema | 5 | COVID-19 | 0.943 | Limited Impact of SARS-CoV-2 on the Human Naso-Oropharyngeal Microbiota in Hospitalized Patients | Experiment 6: local controls vs Hospitalized patients |
| Pancreatic carcinoma | 1 | Cervical glandular intraepithelial neoplasia,Cervical cancer | 0.970 | Taxonomic and Functional Differences in Cervical Microbiome Associated with Cervical Cancer Development | Experiment 1: Healthy controls vs CIN2/3-CC |
| Pancreatic carcinoma | 2 | Autism spectrum disorder | 0.968 | Altered Gut Microbiota and Short-chain Fatty Acids in Chinese Children with Constipated Autism Spectrum Disorder | Experiment 1: Typically developing (TD) group vs Constipated autism spectrum disorder (C-ASD) group |
| Pancreatic carcinoma | 3 | Helminthiasis | 0.967 | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia | Experiment 15: Uninfected (Zero STH infection) vs Heavy Infected (Necator) |
| Pancreatic carcinoma | 4 | Gestational diabetes | 0.966 | Gestational diabetes-related gut microbiome dysbiosis is not influenced by different Asian ethnicities and dietary interventions: a pilot study | Experiment 5: Non-Gestational diabetes mellitus (non-GDM) vs Gestational diabetes mellitus (GDM) |
| Pancreatic carcinoma | 5 | Endometriosis | 0.965 | Gut microbiota imbalance and its correlations with hormone and inflammatory factors in patients with stage 3/4 endometriosis | Experiment 1: Healthy controls vs Stage 3/4 endometriosis patients |

### Cross-Target Check

| query condition | any other target in top 5 profiles? | matched target conditions |
| --- | --- | --- |
| Multiple sclerosis | no | NA |
| Eczema | no | NA |
| Pancreatic carcinoma | no | NA |

Top five condition-centroid neighbors are also written to `target_condition_centroid_neighbors.tsv`; this checks whether the three conditions find each other at condition rather than individual-profile level.

| query condition | rank | centroid neighbor | similarity | other target? |
| --- | ---: | --- | ---: | --- |
| Multiple sclerosis | 1 | Polycystic ovary syndrome | 0.974 | no |
| Multiple sclerosis | 2 | Obesity | 0.972 | no |
| Multiple sclerosis | 3 | Psoriasis | 0.972 | no |
| Multiple sclerosis | 4 | Chronic kidney disease | 0.972 | no |
| Multiple sclerosis | 5 | Breast cancer | 0.971 | no |
| Eczema | 1 | Breast cancer | 0.987 | no |
| Eczema | 2 | Schizophrenia | 0.983 | no |
| Eczema | 3 | Colorectal cancer | 0.982 | no |
| Eczema | 4 | Parkinson's disease | 0.981 | no |
| Eczema | 5 | Major depressive disorder | 0.980 | no |
| Pancreatic carcinoma | 1 | Chronic constipation | 0.990 | no |
| Pancreatic carcinoma | 2 | Gestational diabetes | 0.990 | no |
| Pancreatic carcinoma | 3 | HIV infection | 0.988 | no |
| Pancreatic carcinoma | 4 | Atopic asthma | 0.987 | no |
| Pancreatic carcinoma | 5 | Premature birth | 0.986 | no |

## How Trait Sets Connect To Profiles

A profile is one BugSigDB study/experiment/direction row. Its trait vector counts KG-Microbe trait annotations among organisms reported as increased or decreased in that experiment. A trait set is a named group of biologically related traits, such as `aerobe_v2` or `anaerobe_v2_depletion`. A profile connects to a trait set when the selected trait columns carry signal in that profile; organism support is then checked by intersecting organisms observed in the same profile with organism-trait annotations from `data/derived/taxa_trait_matrix.tsv`.

| trait set | condition | profiles | mean trait-set signal | mean relative trait-set signal |
| --- | --- | ---: | ---: | ---: |
| aerobe_v2 | Eczema | 20 | 487.133 | 0.02527 |
| aerobe_v2 | Multiple sclerosis | 22 | 423.606 | 0.01540 |
| aerobe_v2 | Pancreatic carcinoma | 24 | 846.625 | 0.02302 |
| anaerobe_v2_depletion | Eczema | 20 | 279.825 | 0.03104 |
| anaerobe_v2_depletion | Multiple sclerosis | 22 | 248.545 | 0.05107 |
| anaerobe_v2_depletion | Pancreatic carcinoma | 24 | 241.688 | 0.02874 |

## Companion Traits Vs Organisms Heatmap

`my_top3_traits_vs_organisms_heatmap.png` plots 45 organisms by 34 traits. Rows are organisms observed in the target-condition BugSigDB profiles; columns are binary KG-Microbe trait annotations selected from the top 3 paired-profile contrast traits from `trait-visualization\trait-profile-investigation\trait_profile_paired_contrast_trait_stats.tsv` plus high-signal target-condition traits selected with the original heatmap ranking logic.

## Files

- `my_top3_outlier_condition_study_experiment_titles.tsv`
- `my_top3_outlier_condition_depth_summary.tsv`
- `my_top3_autism_condition_label_check.tsv`
- `my_top3_outlier_condition_profile_neighbors.tsv`
- `my_top3_outlier_condition_neighbor_crosscheck.tsv`
- `my_top3_target_condition_centroid_neighbors.tsv`
- `my_top3_guild_profile_connection_summary.tsv`
- `my_top3_traits_vs_organisms_heatmap_matrix.tsv`
- `my_top3_traits_vs_organisms_taxa_metadata.tsv`
- `my_top3_traits_vs_organisms_trait_metadata.tsv`
- `my_top3_traits_vs_organisms_heatmap.png`
