# Enrichment Infographic Follow-Up

Source figure: https://github.com/CultureBotAI/bsdbng/blob/enrichment-analysis-heatmap-cluster/data/derived/enrichment_infographic.png

The target conditions for this run are selected from the top 3 paired-profile condition contrast rows from `trait-visualization\trait-profile-investigation\trait_profile_paired_contrast_condition_summary.tsv`. The local follow-up uses the full `data/derived/study_trait_profiles_wide.tsv` matrix and raw BugSigDB experiment metadata to list titles, find similar profiles, and connect enriched traits to organisms.

## Target Conditions

| condition | paired absolute contrast load | local profiles | studies | experiments | local mean nonzero traits/profile | local percentile | local outlier? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Gingivitis | 1.11 | 53 | 8 | 36 | 111.2 | 24.9% | no |
| Gastroesophageal reflux disease | 1.08 | 19 | 1 | 10 | 366.9 | 73.1% | no |
| Lung carcinoma | 1.08 | 12 | 1 | 6 | 359.2 | 72.1% | no |

Interpretation: condition support matters because profile counts, study counts, and body-site mix vary across exact BugSigDB condition labels.

## Study And Experiment Titles

| condition | study | year | experiment | profiles |
| --- | --- | ---: | --- | ---: |
| Gastroesophageal reflux disease | Signatures within the esophageal microbiome are associated with host genetics, age, and disease (bsdb:30558669) | 2018 | Experiment 13: Normal (Normal esophagus) vs GERD (Gastroesophageal reflux disease) | 1 |
| Gastroesophageal reflux disease | Signatures within the esophageal microbiome are associated with host genetics, age, and disease (bsdb:30558669) | 2018 | Experiment 14: Normal (Normal esophagus) vs GERD (Gastroesophageal reflux disease)within subjects classified as PPI NO | 2 |
| Gastroesophageal reflux disease | Signatures within the esophageal microbiome are associated with host genetics, age, and disease (bsdb:30558669) | 2018 | Experiment 18: Normal (Cluster 1) vs GERD(Gastroesophageal reflux disease) (Cluster 1) | 2 |
| Gastroesophageal reflux disease | Signatures within the esophageal microbiome are associated with host genetics, age, and disease (bsdb:30558669) | 2018 | Experiment 19: Normal (Cluster 2) vs GERD(Gastroesophageal reflux disease) (Cluster 2) | 2 |
| Gastroesophageal reflux disease | Signatures within the esophageal microbiome are associated with host genetics, age, and disease (bsdb:30558669) | 2018 | Experiment 20: Normal (Cluster 3) vs GERD(Gastroesophageal reflux disease) (Cluster 3) | 2 |
| Gastroesophageal reflux disease | Signatures within the esophageal microbiome are associated with host genetics, age, and disease (bsdb:30558669) | 2018 | Experiment 21: Normal (Cluster 1) vs GERD(Gastroesophageal reflux disease) (Cluster 1) | 2 |
| Gastroesophageal reflux disease | Signatures within the esophageal microbiome are associated with host genetics, age, and disease (bsdb:30558669) | 2018 | Experiment 22: Normal (Cluster 2) vs GERD(Gastroesophageal reflux disease) (Cluster 2) | 2 |
| Gastroesophageal reflux disease | Signatures within the esophageal microbiome are associated with host genetics, age, and disease (bsdb:30558669) | 2018 | Experiment 23: Normal (Cluster 3) vs GERD(Gastroesophageal reflux disease) (Cluster 3) | 2 |
| Gastroesophageal reflux disease | Signatures within the esophageal microbiome are associated with host genetics, age, and disease (bsdb:30558669) | 2018 | Experiment 7: Normal (Normal esophagus) vs GERD (Gastroesophageal reflux disease) | 2 |
| Gastroesophageal reflux disease | Signatures within the esophageal microbiome are associated with host genetics, age, and disease (bsdb:30558669) | 2018 | Experiment 8: Normal vs Gastroesophageal reflux disease GERD(within subjects classified as PPI NO) | 2 |
| Gingivitis | Subgingival Microbiome of Gingivitis in Chinese Undergraduates (bsdb:28808698) | 2017 | Experiment 1: healthy controls vs gingivitis | 2 |
| Gingivitis | Subgingival Microbiota during Healthy Pregnancy and Pregnancy Gingivitis (bsdb:32777190) | 2021 | Experiment 2: pregnant women without gingivitis vs pregnant women with gingivitis | 1 |
| Gingivitis | Subgingival Microbiota during Healthy Pregnancy and Pregnancy Gingivitis (bsdb:32777190) | 2021 | Experiment 3: pregnant women without gingivitis vs pregnant women with gingivitis | 2 |
| Gingivitis | Microbiomes in supragingival biofilms and saliva of adolescents with gingivitis and gingival health (bsdb:33876475) | 2022 | Experiment 1: Healthy controls vs Induced gingivitis patients | 2 |
| Gingivitis | Microbiomes in supragingival biofilms and saliva of adolescents with gingivitis and gingival health (bsdb:33876475) | 2022 | Experiment 2: Spontaneous gingivitis patients vs Induced gingivitis patients | 2 |
| Gingivitis | Microbiomes in supragingival biofilms and saliva of adolescents with gingivitis and gingival health (bsdb:33876475) | 2022 | Experiment 3: Healthy controls vs Spontaneous gingivitis patients | 2 |
| Gingivitis | The microbiome of dental and peri-implant subgingival plaque during peri-implant mucositis therapy: A randomized clinical trial (bsdb:34664294) | 2022 | Experiment 1: Healthy controls vs Gingivitis cases | 2 |
| Gingivitis | The microbiome of dental and peri-implant subgingival plaque during peri-implant mucositis therapy: A randomized clinical trial (bsdb:34664294) | 2022 | Experiment 16: Placebo vs Delmopinol | 2 |
| Gingivitis | The microbiome of dental and peri-implant subgingival plaque during peri-implant mucositis therapy: A randomized clinical trial (bsdb:34664294) | 2022 | Experiment 17: Placebo vs Chlorhexidine | 2 |
| Gingivitis | Interdental oral hygiene interventions elicit varying compositional microbiome changes in naturally occurring gingivitis: Secondary data analysis from a clinical trial (bsdb:38088457) | 2024 | Experiment 1: Oral-B Water Flosser + manual toothbrushing vs Waterpik Water Flosser + manual toothbrushing | 2 |
| Gingivitis | Interdental oral hygiene interventions elicit varying compositional microbiome changes in naturally occurring gingivitis: Secondary data analysis from a clinical trial (bsdb:38088457) | 2024 | Experiment 2: Dental floss + manual toothbrushing vs Waterpik Water Flosser + manual toothbrushing | 2 |
| Gingivitis | Interdental oral hygiene interventions elicit varying compositional microbiome changes in naturally occurring gingivitis: Secondary data analysis from a clinical trial (bsdb:38088457) | 2024 | Experiment 3: Dental floss + manual toothbrushing vs Oral-B Water Flosser + manual toothbrushing | 2 |
| Gingivitis | Interdental oral hygiene interventions elicit varying compositional microbiome changes in naturally occurring gingivitis: Secondary data analysis from a clinical trial (bsdb:38088457) | 2024 | Experiment 4: Waterpik Water Flosser + manual toothbrushing vs Dental floss + manual toothbrushing | 2 |
| Gingivitis | Interdental oral hygiene interventions elicit varying compositional microbiome changes in naturally occurring gingivitis: Secondary data analysis from a clinical trial (bsdb:38088457) | 2024 | Experiment 5: Waterpik Water Flosser + manual toothbrushing vs Oral-B Water Flosser + manual toothbrushing | 2 |
| Gingivitis | Interdental oral hygiene interventions elicit varying compositional microbiome changes in naturally occurring gingivitis: Secondary data analysis from a clinical trial (bsdb:38088457) | 2024 | Experiment 6: Dental floss + manual toothbrushing vs Oral-B Water Flosser + manual toothbrushing | 2 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 1: Healthy subjects vs Gingivitis subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 10: Brushing and Flossing (BF) subjects vs Brushing and rinsing with Listerine® Cool Mint® Antiseptic (BA) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 11: Brushing and Flossing (BF) subjects vs Brushing and rinsing with Listerine® Cool Mint® Zero Alcohol group (BZ) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 12: Brushing and Flossing (BF) subjects vs Brushing, flossing, and rinsing with Listerine® Cool® Mint Zero Alcohol group (BFZ) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 13: Brushing and rinsing with Listerine® Cool Mint® Zero Alcohol group (BZ) subjects vs Brushing, flossing, and rinsing with Listerine® Cool® Mint Zero Alcohol group (BFZ) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 2: Brushing Control (B) vs Brushing and Flossing (BF) subjects | 2 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 3: Brushing Control (B) vs Brushing and rinsing with Listerine® Cool Mint® Antiseptic (BA) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 4: Brushing Control (B) vs Brushing and rinsing with Listerine® Cool Mint® Zero Alcohol group (BZ) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 5: Brushing Control (B) vs Brushing, flossing, and rinsing with Listerine® Cool® Mint Zero Alcohol group (BFZ) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 6: Brushing and Flossing (BF) subjects vs Brushing, flossing, and rinsing with Listerine® Cool® Mint Zero Alcohol group (BFZ) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 7: Brushing and rinsing with Listerine® Cool Mint® Zero Alcohol group (BZ) subjects vs Brushing, flossing, and rinsing with Listerine® Cool® Mint Zero Alcohol group (BFZ) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 8: Brushing Control (B) vs Brushing and Flossing (BF) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of brushing, flossing, and mouthrinsing on supragingival and subgingival plaque microbiota: 12-week clinical trial (bsdb:38760758) | 2024 | Experiment 9: Brushing and Flossing (BF) subjects vs Brushing and rinsing with Listerine® Cool Mint® Zero Alcohol group (BZ) subjects | 1 |
| Gingivitis | Quantitative analysis of the effects of essential oil mouthrinses on clinical plaque microbiome: a parallel-group, randomized trial (bsdb:38762482) | 2024 | Experiment 1: Healthy subjects vs Gingivitis subjects | 2 |
| Gingivitis | Quantitative analysis of the effects of essential oil mouthrinses on clinical plaque microbiome: a parallel-group, randomized trial (bsdb:38762482) | 2024 | Experiment 2: Hydroalcohol negative control (HA) week 4 vs Listerine Cool Mint Antiseptic mouthrinse (LCM) week 4 | 1 |
| Gingivitis | Quantitative analysis of the effects of essential oil mouthrinses on clinical plaque microbiome: a parallel-group, randomized trial (bsdb:38762482) | 2024 | Experiment 3: Hydroalcohol negative control (HA) week 4 vs Alcohol-containing prototype mouthrinse (ACPM) week 4 | 1 |
| Gingivitis | Quantitative analysis of the effects of essential oil mouthrinses on clinical plaque microbiome: a parallel-group, randomized trial (bsdb:38762482) | 2024 | Experiment 4: Hydroalcohol negative control (HA) week 4 vs Alcohol-free prototype mouthrinse (AFPM) week 4 | 1 |
| Gingivitis | Quantitative analysis of the effects of essential oil mouthrinses on clinical plaque microbiome: a parallel-group, randomized trial (bsdb:38762482) | 2024 | Experiment 5: Hydroalcohol negative control (HA) week 6 vs Listerine Cool Mint Antiseptic mouthrinse (LCM) week 6 | 1 |
| Gingivitis | Quantitative analysis of the effects of essential oil mouthrinses on clinical plaque microbiome: a parallel-group, randomized trial (bsdb:38762482) | 2024 | Experiment 6: Hydroalcohol negative control (HA) week 6 vs Alcohol-containing prototype mouthrinse (ACPM) week 6 | 1 |
| Gingivitis | Quantitative analysis of the effects of essential oil mouthrinses on clinical plaque microbiome: a parallel-group, randomized trial (bsdb:38762482) | 2024 | Experiment 7: Hydroalcohol negative control (HA) week 6 vs Alcohol-free prototype mouthrinse (AFPM) week 6 | 1 |
| Gingivitis | Potential biomarkers for early periodontal inflammation: investigating CD5<sup>+</sup> B cells, salivary cytokines and oral microbiome (bsdb:41708693) | 2026 | Experiment 2: Healthy controls vs Gingivitis | 2 |
| Lung carcinoma | Characteristics of gut and lung microbiota in patients with lung masses and their relationship with clinical features (bsdb:40849619) | 2025 | Experiment 1: Benign mass vs Malignant mass | 2 |
| Lung carcinoma | Characteristics of gut and lung microbiota in patients with lung masses and their relationship with clinical features (bsdb:40849619) | 2025 | Experiment 2: Benign mass vs Malignant mass | 2 |
| Lung carcinoma | Characteristics of gut and lung microbiota in patients with lung masses and their relationship with clinical features (bsdb:40849619) | 2025 | Experiment 3: Benign mass vs Malignant mass | 2 |
| Lung carcinoma | Characteristics of gut and lung microbiota in patients with lung masses and their relationship with clinical features (bsdb:40849619) | 2025 | Experiment 4: Benign mass vs Malignant mass | 2 |
| Lung carcinoma | Characteristics of gut and lung microbiota in patients with lung masses and their relationship with clinical features (bsdb:40849619) | 2025 | Experiment 5: Small cell lung cancer (SCLC) vs Non-small cell lung cancer (NSCLC) | 2 |
| Lung carcinoma | Characteristics of gut and lung microbiota in patients with lung masses and their relationship with clinical features (bsdb:40849619) | 2025 | Experiment 6: Small cell lung cancer (NSCLC) vs Non-small cell lung cancer (NSCLC) | 2 |

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
| Gingivitis | 1 | Autism spectrum disorder | 0.926 | Large-scale metagenomic analysis of oral microbiomes reveals markers for autism spectrum disorders | Experiment 3: Low Microbial Load vs High Microbial Load |
| Gingivitis | 2 | Mycobacterium tuberculosis | 0.923 | Full-length 16S rDNA sequencing based on Oxford Nanopore Technologies revealed the association between gut-pharyngeal microbiota and tuberculosis in cynomolgus macaques | Experiment 16: Tuberculosis Control (TB[-]control) vs Tuberculosis Latent (TB[+]Latent) |
| Gingivitis | 3 | Cognitive impairment | 0.920 | Large-scale metagenomic analysis of oral microbiomes reveals markers for autism spectrum disorders | Experiment 8: Low Intelligence Quotients(Low IQ ≤70) vs High Intelligence Quotients(High IQ > 85) |
| Gingivitis | 4 | Cognitive function measurement | 0.920 | Large-scale metagenomic analysis of oral microbiomes reveals markers for autism spectrum disorders | Experiment 6: Lower Social Communication Questionnaire (Low SCQ) scale scale vs Upper Social Communication Questionnaire (Upper SCQ) scale |
| Gingivitis | 5 | Head and neck squamous cell carcinoma | 0.917 | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer | Experiment 4: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants |
| Gastroesophageal reflux disease | 1 | Gastric cancer | 0.943 | Mucosa-Associated Microbiota in Gastric Cancer Tissues Compared With Non-cancer Tissues | Experiment 3: adjacent normal tissue vs cancer tissue |
| Gastroesophageal reflux disease | 2 | COVID-19 | 0.941 | Exploring nasopharyngeal microbiota profile in children affected by SARS-CoV-2 infection | Experiment 1: healthy controls vs COVID-19 patients |
| Gastroesophageal reflux disease | 3 | Fungal infectious disease | 0.940 | Associations of Batrachochytrium dendrobatidis with skin bacteria and fungi on Asian amphibian hosts | Experiment 9: Theloderma rhododiscus frog uninfected with Batrachocytrium dendrobatidis (Bd-) vs Theloderma rhododiscus frog infected with Batrachocytrium dendrobatidis (Bd+) |
| Gastroesophageal reflux disease | 4 | Colorectal cancer | 0.940 | Distinct intratumoral microbiome of young-onset and average-onset colorectal cancer | Experiment 1: Average-Onset Colorectal Cancer(aoCRC)(tumor tissue) vs Young-Onset Colorectal Cancer(yoCRC)(tumor tissue) |
| Gastroesophageal reflux disease | 5 | Environmental factor | 0.939 | Shared environments can facilitate microbial transmission and alter metabolic outcomes | Experiment 5: 4-weeks of SAS (single vinyl isolator with a shared air supply) Donor pair-1 (dp1) Thailand (Thai) colonized mice vs 4-weeks of SAS (single vinyl isolator with a shared air supply) Donor pair-1 (dp1) United States (US) colonized mice |
| Lung carcinoma | 1 | Species design | 0.902 | Captive environments reshape the compositions of carbohydrate active enzymes and virulence factors in wolf gut microbiome | Experiment 3: Red fox vs Wild wolf |
| Lung carcinoma | 2 | Breed | 0.901 | Dysbiosis of skin microbiome and gut microbiome in melanoma progression | Experiment 8: Crossbreed between MeLiM and white piglets vs Homogenous breed MeLiM piglets |
| Lung carcinoma | 3 | Chronic constipation | 0.901 | Enteric nervous system damage caused by abnormal intestinal butyrate metabolism may lead to functional constipation | Experiment 2: Group 2 vs Group 3 |
| Lung carcinoma | 4 | Endometriosis | 0.898 | Associations Between Endometriosis and Gut Microbiota | Experiment 1: Healthy controls taken from a larger cohort study vs Endometriosis patients |
| Lung carcinoma | 5 | HIV infection | 0.898 | Impact of Acute HIV Infection and Early Antiretroviral Therapy on the Human Gut Microbiome | Experiment 1: Human immunodeficiency virus negative (HIV-)Uninfected Controls vs Human immunodeficiency virus positive(HIV+) pre - antiretroviral therapy(ART) |

### Cross-Target Check

| query condition | any other target in top 5 profiles? | matched target conditions |
| --- | --- | --- |
| Gingivitis | no | NA |
| Gastroesophageal reflux disease | no | NA |
| Lung carcinoma | no | NA |

Top five condition-centroid neighbors are also written to `target_condition_centroid_neighbors.tsv`; this checks whether the three conditions find each other at condition rather than individual-profile level.

| query condition | rank | centroid neighbor | similarity | other target? |
| --- | ---: | --- | ---: | --- |
| Gingivitis | 1 | Periodontitis | 0.975 | no |
| Gingivitis | 2 | Esophageal squamous cell carcinoma | 0.973 | no |
| Gingivitis | 3 | Squamous cell carcinoma | 0.969 | no |
| Gingivitis | 4 | Chronic periodontitis | 0.968 | no |
| Gingivitis | 5 | SARS coronavirus | 0.964 | no |
| Gastroesophageal reflux disease | 1 | Gastric cancer | 0.989 | no |
| Gastroesophageal reflux disease | 2 | Health study participation | 0.988 | no |
| Gastroesophageal reflux disease | 3 | Lung cancer | 0.987 | no |
| Gastroesophageal reflux disease | 4 | Breast cancer | 0.985 | no |
| Gastroesophageal reflux disease | 5 | Sampling site | 0.984 | no |
| Lung carcinoma | 1 | Crohn's disease | 0.979 | no |
| Lung carcinoma | 2 | Alzheimer's disease | 0.978 | no |
| Lung carcinoma | 3 | Colorectal cancer | 0.975 | no |
| Lung carcinoma | 4 | Chronic kidney disease | 0.974 | no |
| Lung carcinoma | 5 | Chronic hepatitis B virus infection | 0.974 | no |

## How Guilds Connect To Profiles

A profile is one BugSigDB study/experiment/direction row. Its trait vector counts KG-Microbe trait annotations among organisms reported as increased or decreased in that experiment. A guild is a named set of biologically related traits, such as `aerobe_v2` or `anaerobe_v2_depletion`. A profile connects to a guild when the guild's trait columns carry signal in that profile; organism support is then checked by intersecting organisms observed in the same profile with organism-trait annotations from `data/derived/taxa_trait_matrix.tsv`.

| guild | condition | profiles | mean guild signal | mean relative guild signal |
| --- | --- | ---: | ---: | ---: |
| aerobe_v2 | Gingivitis | 53 | 101.409 | 0.02880 |
| aerobe_v2 | Gastroesophageal reflux disease | 19 | 1010.351 | 0.02852 |
| aerobe_v2 | Lung carcinoma | 12 | 664.306 | 0.01281 |
| anaerobe_v2_depletion | Gingivitis | 53 | 62.311 | 0.03084 |
| anaerobe_v2_depletion | Gastroesophageal reflux disease | 19 | 400.079 | 0.02179 |
| anaerobe_v2_depletion | Lung carcinoma | 12 | 404.458 | 0.03006 |

## Companion Traits Vs Organisms Heatmap

`my_top3_conditions_relevant_traits_traits_vs_organisms_heatmap.png` plots 45 organisms by 34 traits. Rows are organisms observed in the target-condition BugSigDB profiles; columns are binary KG-Microbe trait annotations selected from the top 3 paired-profile contrast traits from `trait-visualization\trait-profile-investigation\trait_profile_paired_contrast_trait_stats.tsv` plus high-signal target-condition traits selected with the original heatmap ranking logic.

## Files

- `my_top3_conditions_relevant_traits_outlier_condition_study_experiment_titles.tsv`
- `my_top3_conditions_relevant_traits_outlier_condition_depth_summary.tsv`
- `my_top3_conditions_relevant_traits_autism_condition_label_check.tsv`
- `my_top3_conditions_relevant_traits_outlier_condition_profile_neighbors.tsv`
- `my_top3_conditions_relevant_traits_outlier_condition_neighbor_crosscheck.tsv`
- `my_top3_conditions_relevant_traits_target_condition_centroid_neighbors.tsv`
- `my_top3_conditions_relevant_traits_guild_profile_connection_summary.tsv`
- `my_top3_conditions_relevant_traits_traits_vs_organisms_heatmap_matrix.tsv`
- `my_top3_conditions_relevant_traits_traits_vs_organisms_taxa_metadata.tsv`
- `my_top3_conditions_relevant_traits_traits_vs_organisms_trait_metadata.tsv`
- `my_top3_conditions_relevant_traits_traits_vs_organisms_heatmap.png`
