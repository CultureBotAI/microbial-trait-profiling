# Enrichment Infographic Follow-Up

Source figure: https://github.com/CultureBotAI/bsdbng/blob/enrichment-analysis-heatmap-cluster/data/derived/enrichment_infographic.png

The target conditions for this run are selected from the top 3 paired-profile condition contrast rows from `trait-visualization/trait-profile-investigation/trait_profile_paired_contrast_condition_summary.tsv`. The local follow-up uses the full `data/derived/study_trait_profiles_wide.tsv` matrix to list titles, find similar profiles, and connect enriched traits to organisms. Profiles carry their own experiment's condition, body site and host from `trait-visualization/study_trait_profile_study_metadata.tsv`; experiment details and reported taxa come from the pinned `data/raw/full_dump.csv`, and only profiles in this corpus are counted.

## Target Conditions

| condition | paired absolute contrast load | local profiles | studies | experiments | local mean nonzero traits/profile | local percentile | local outlier? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Multiple sclerosis | 1.00 | 28 | 6 | 14 | 97.1 | 40.3% | no |
| Stroke | 0.92 | 43 | 5 | 28 | 134.0 | 68.0% | no |
| Constipation | 0.76 | 23 | 7 | 12 | 104.0 | 45.7% | no |

Interpretation: condition support matters because profile counts, study counts, and body-site mix vary across exact BugSigDB condition labels.

## Study And Experiment Titles

| condition | study | year | experiment | profiles |
| --- | --- | ---: | --- | ---: |
| Constipation | Structural changes in the gut microbiome of constipated patients (bsdb:25073603) | 2014 | Experiment 1: Control vs Constipation | 2 |
| Constipation | The Effect of Psyllium Husk on Intestinal Microbiota in Constipated Patients and Healthy Controls (bsdb:30669509) | 2019 | Experiment 1: healthy participants at  baseline vs healthy participants following high dose of Psyllium supplementation | 2 |
| Constipation | The Effect of Psyllium Husk on Intestinal Microbiota in Constipated Patients and Healthy Controls (bsdb:30669509) | 2019 | Experiment 2: constipated patients at baseline vs constipated participants following Psyllium supplementation | 2 |
| Constipation | The Effect of Psyllium Husk on Intestinal Microbiota in Constipated Patients and Healthy Controls (bsdb:30669509) | 2019 | Experiment 3: healthy participants at baseline vs constipated participants at baseline | 2 |
| Constipation | Fecal Microbiota Transplantation Modulates the Gut Flora Favoring Patients With Functional Constipation (bsdb:34690948) | 2021 | Experiment 2: bef (pre-Fecal microbiota transplantation) vs aft (post-Fecal microbiota transplantation) | 2 |
| Constipation | Fecal Microbiota Transplantation Modulates the Gut Flora Favoring Patients With Functional Constipation (bsdb:34690948) | 2021 | Experiment 5: Low BSFS (Bristol stool form scale) score vs High BSFS (Bristol stool form scale) score | 1 |
| Constipation | Protective effect of L-pipecolic acid on constipation in C57BL/6 mice based on gut microbiome and serum metabolomic (bsdb:37210496) | 2023 | Experiment 1: Healthy Controls vs Functional Constipation(FC) | 2 |
| Constipation | The gut microbiome and metabolites are altered and interrelated in patients with functional constipation (bsdb:38125567) | 2023 | Experiment 5: Low level of straining during defecation vs High level of straining during defecation | 2 |
| Constipation | Altered gut microbial profile is associated with differentially expressed fecal microRNAs in patients with functional constipation (bsdb:38274754) | 2023 | Experiment 5: Low levels of straining during defecation vs High levels of straining during defecation | 2 |
| Constipation | Intestinal microbiome changes and mechanisms of maintenance hemodialysis patients with constipation (bsdb:39588509) | 2024 | Experiment 2: Healthy Control Group (HCG) vs Maintenance Hemodialysis with Constipation (MHDCG) | 2 |
| Constipation | Intestinal microbiome changes and mechanisms of maintenance hemodialysis patients with constipation (bsdb:39588509) | 2024 | Experiment 4: Maintenance Hemodialysis without Constipation (MHDNCG) vs Maintenance Hemodialysis with Constipation (MHDCG) | 2 |
| Constipation | Intestinal microbiome changes and mechanisms of maintenance hemodialysis patients with constipation (bsdb:39588509) | 2024 | Experiment 5: Maintenance Hemodialysis without Constipation (MHDNCG) vs Maintenance Hemodialysis with Constipation (MHDCG) | 2 |
| Multiple sclerosis | Dysbiosis in the Gut Microbiota of Patients with Multiple Sclerosis, with a Striking Depletion of Species Belonging to Clostridia XIVa and IV Clusters (bsdb:26367776) | 2015 | Experiment 1: Healthy subjects (controls) vs Multiple sclerosis (MS), | 2 |
| Multiple sclerosis | Dysbiosis in the Gut Microbiota of Patients with Multiple Sclerosis, with a Striking Depletion of Species Belonging to Clostridia XIVa and IV Clusters (bsdb:26367776) | 2015 | Experiment 2: Healthy subjects (controls) vs Multiple sclerosis (MS), | 2 |
| Multiple sclerosis | Gut microbiota in early pediatric multiple sclerosis: a case-control study (bsdb:27176462) | 2016 | Experiment 1: healthy controls vs pediatric multiple sclerosis | 2 |
| Multiple sclerosis | Gut microbiota in early pediatric multiple sclerosis: a case-control study (bsdb:27176462) | 2016 | Experiment 2: IMD naive cases vs IMD exposure multiple sclerosis cases | 2 |
| Multiple sclerosis | Multiple sclerosis patients have a distinct gut microbiota compared to healthy controls (bsdb:27346372) | 2016 | Experiment 1: Healthy controls vs Multiple sclerosis | 2 |
| Multiple sclerosis | Alterations of the human gut microbiome in multiple sclerosis (bsdb:27352007) | 2016 | Experiment 1: healthy controls vs multiple sclerosis | 2 |
| Multiple sclerosis | Alterations of the human gut microbiome in multiple sclerosis (bsdb:27352007) | 2016 | Experiment 2: controls vs multiple sclerosis patients | 2 |
| Multiple sclerosis | Alterations of the human gut microbiome in multiple sclerosis (bsdb:27352007) | 2016 | Experiment 3: controls vs untreated multiple sclerosis | 2 |
| Multiple sclerosis | Alterations of the human gut microbiome in multiple sclerosis (bsdb:27352007) | 2016 | Experiment 4: controls vs untreated multiple sclerosis | 2 |
| Multiple sclerosis | Alterations of the human gut microbiome in multiple sclerosis (bsdb:27352007) | 2016 | Experiment 5: untreated multiple sclerosis vs treated multiple sclerosis | 2 |
| Multiple sclerosis | Alterations of the human gut microbiome in multiple sclerosis (bsdb:27352007) | 2016 | Experiment 6: untreated multiple sclerosis vs treated multiple sclerosis | 2 |
| Multiple sclerosis | A comparative study of the gut microbiota in immune-mediated inflammatory diseases-does a common dysbiosis exist? (bsdb:30545401) | 2018 | Experiment 3: healthy controls vs multiple sceloris | 2 |
| Multiple sclerosis | Multiple sclerosis and gut microbiota: Lachnospiraceae from the ileum of MS twins trigger MS-like disease in germfree transgenic mice-An unbiased functional study (bsdb:40258140) | 2025 | Experiment 1: Unaffected twins vs Affected twins | 2 |
| Multiple sclerosis | Multiple sclerosis and gut microbiota: Lachnospiraceae from the ileum of MS twins trigger MS-like disease in germfree transgenic mice-An unbiased functional study (bsdb:40258140) | 2025 | Experiment 2: Unaffected twins vs Affected twins | 2 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 1: Low Risk (LR) and Medium Risk (MR) vs High Risk (HR) | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 10: Low Risk (LR) and Medium Risk (MR) vs High Risk (HR) | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 11: Low Risk (LR) and High Risk (HR) vs Medium Risk (MR) | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 13: Low Risk (LR) vs High Risk (HR) | 2 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 14: Low Risk (LR) vs High Risk (HR) | 2 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 15: Low Risk (LR) vs Medium Risk (MR) | 2 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 16: Medium Risk (LR) and High Risk (HR) vs Low Risk | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 2: Low Risk (LR) and High Risk (HR) vs Medium Risk (MR) | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 3: High Risk (HR) and Medium Risk (MR) vs Low Risk (LR) | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 4: Low Risk (LR) and Medium Risk (MR) Without Medication use vs High Risk (HR) Without Medication use | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 5: Low Risk (LR) and High Risk (HR) Without Medication use vs Medium Risk (MR) Without Medication use | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 6: High Risk (HR) and Medium Risk (MR) Without Medication use vs Low Risk (LR) Without Medication use | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 7: Low Risk (LR) and Medium Risk (MR) With Medication use vs High Risk (HR) With  Medication use | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 8: Low Risk (LR) and High Risk (HR) With Medication use vs Medium Risk (MR) With Medication use | 1 |
| Stroke | Higher Risk of Stroke Is Correlated With Increased Opportunistic Pathogen Load and Reduced Levels of Butyrate-Producing Bacteria in the Gut (bsdb:30778376) | 2019 | Experiment 9: High Risk (HR) and Medium Risk (MR) With Medication use vs Low Risk (LR) With Medication use | 1 |
| Stroke | Dysbiosis of Gut Microbiota Is an Independent Risk Factor of Stroke-Associated Pneumonia: A Chinese Pilot Study (bsdb:34414134) | 2021 | Experiment 1: Non-Stroke Associated Pneumonia (non-SAP) Training cohort vs Stroke Associated Pneumonia (SAP) Training cohort | 2 |
| Stroke | Dysbiosis of Gut Microbiota Is an Independent Risk Factor of Stroke-Associated Pneumonia: A Chinese Pilot Study (bsdb:34414134) | 2021 | Experiment 2: Non-Stroke Associated Pneumonia (non-SAP) independent cohort vs Stroke Associated Pneumonia (SAP) independent cohort | 2 |
| Stroke | Dysbiosis of Gut Microbiota Is an Independent Risk Factor of Stroke-Associated Pneumonia: A Chinese Pilot Study (bsdb:34414134) | 2021 | Experiment 3: Non-Stroke Associated Pneumonia (non-SAP) subgroup in the training cohort vs Stroke Associated Pneumonia (SAP) subgroup in the training cohort | 2 |
| Stroke | Dysbiosis of Gut Microbiota Is an Independent Risk Factor of Stroke-Associated Pneumonia: A Chinese Pilot Study (bsdb:34414134) | 2021 | Experiment 4: Non-Stroke Associated Pneumonia (non-SAP) Training cohort vs Stroke Associated Pneumonia (SAP) Training cohort | 2 |
| Stroke | Dysbiosis of Gut Microbiota Is an Independent Risk Factor of Stroke-Associated Pneumonia: A Chinese Pilot Study (bsdb:34414134) | 2021 | Experiment 5: Non-Stroke Associated Pneumonia (non-SAP) subgroup in the training cohort vs Stroke Associated Pneumonia (SAP) subgroup in the training cohort | 2 |
| Stroke | Dysbiosis of Gut Microbiota Is an Independent Risk Factor of Stroke-Associated Pneumonia: A Chinese Pilot Study (bsdb:34414134) | 2021 | Experiment 6: Non-Stroke Associated Pneumonia (non-SAP) vs Stroke Associated Pneumonia (non-SAP) | 1 |
| Stroke | Dysbiosis of Gut Microbiota Is an Independent Risk Factor of Stroke-Associated Pneumonia: A Chinese Pilot Study (bsdb:34414134) | 2021 | Experiment 8: Non-Stroke Associated Pneumonia (non-SAP) Training cohort vs Stroke Associated Pneumonia (SAP) Training cohort | 2 |
| Stroke | Dysbiosis of Gut Microbiota Is an Independent Risk Factor of Stroke-Associated Pneumonia: A Chinese Pilot Study (bsdb:34414134) | 2021 | Experiment 9: Non-Stroke Associated Pneumonia (non-SAP) validation cohort vs Stroke Associated Pneumonia (SAP) validation cohort | 2 |
| Stroke | Distinctive Gut Microbiota Alteration Is Associated with Poststroke Functional Recovery: Results from a Prospective Cohort Study (bsdb:34917142) | 2021 | Experiment 1: healthy controls vs stroke patients | 2 |
| Stroke | Distinctive Gut Microbiota Alteration Is Associated with Poststroke Functional Recovery: Results from a Prospective Cohort Study (bsdb:34917142) | 2021 | Experiment 2: subacute vs chronic patients | 2 |
| Stroke | Pivotal interplays between fecal metabolome and gut microbiome reveal functional signatures in cerebral ischemic stroke (bsdb:36209079) | 2022 | Experiment 1: Healthy individuals vs Cerebral ischemic stroke patients | 2 |
| Stroke | Imbalance of Microbacterial Diversity Is Associated with Functional Prognosis of Stroke (bsdb:37197229) | 2023 | Experiment 1: Healthy Control vs Post Stroke Patients | 2 |
| Stroke | Imbalance of Microbacterial Diversity Is Associated with Functional Prognosis of Stroke (bsdb:37197229) | 2023 | Experiment 2: Healthy Control vs Post Stroke Patients | 2 |

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
| Multiple sclerosis | 1 | HIV infection | 0.935 | Longitudinal gut microbiota composition of South African and Nigerian infants in relation to tetanus vaccine responses | Experiment 3: Infants who are HIV unexposed and uninfected (iHUU) at week 15 vs Infants who are HIV exposed but uninfected (iHEU) at week 15 |
| Multiple sclerosis | 2 | Chlamydia trachomatis | 0.931 | Cervicovaginal microbiome and natural history of Chlamydia trachomatis in adolescents and young women | Experiment 1: Chlamydia trichomonas Controls (t0) vs Chlamydia trichomonas Cases (t0) |
| Multiple sclerosis | 3 | Adverse effect | 0.931 | Intestinal microbiota influences clinical outcome and side effects of early breast cancer treatment | Experiment 13: No - increased BMI (Post-Chemo) vs Yes - increased BMI (Post-Chemo) |
| Multiple sclerosis | 4 | Root caries | 0.929 | Comparative analysis of microbiome in coronal and root caries | Experiment 1: Coronal Caries group vs Root Caries group |
| Multiple sclerosis | 5 | Age | 0.927 | Faecal microbiota of schoolchildren is associated with nutritional status and markers of inflammation: a double-blinded cluster-randomized controlled trial using multi-micronutrient fortified rice | Experiment 8: 10-14 years old vs 6-9 years old |
| Stroke | 1 | Response to drug | 0.959 | The Parkinson's disease drug entacapone disrupts gut microbiome homoeostasis via iron sequestration | Experiment 10: No drug - Dimethyl sulfoxide only (DMSO only) vs Entacapone High Concentration(ENT-High) after 24 hrs of incubation - Donor 3 |
| Stroke | 2 | Response to drug | 0.959 | The Parkinson's disease drug entacapone disrupts gut microbiome homoeostasis via iron sequestration | Experiment 9: No drug - Dimethyl sulfoxide only (DMSO only) vs Entacapone High Concentration(ENT-High) after 24 hrs of incubation - Donor 2 |
| Stroke | 3 | HIV infection | 0.957 | Impact of Acute HIV Infection and Early Antiretroviral Therapy on the Human Gut Microbiome | Experiment 1: Human immunodeficiency virus negative (HIV-)Uninfected Controls vs Human immunodeficiency virus positive(HIV+) pre - antiretroviral therapy(ART) |
| Stroke | 4 | Parkinson's disease | 0.956 | Gut microbiota in Parkinson's disease: Temporal stability and relations to disease progression | Experiment 4: Tremor-dominant phenotype (TD) at baseline vs Postural instability and gait difficulty (PIGD) phenotype at baseline |
| Stroke | 5 | Antimicrobial agent | 0.956 | Association of prenatal antibiotics with measures of infant adiposity and the gut microbiome | Experiment 3: 12-month-old infants not exposed to antibiotics during the second trimester vs 12-month-old infants exposed to antibiotics during the second trimester |
| Constipation | 1 | COVID-19 | 0.953 | The gut microbiome of COVID-19 recovered patients returns to uninfected status in a minority-dominated United States cohort | Experiment 1: Non-infected patients with unrelated respiratory medical conditions vs SARS-CoV-2 infected patients |
| Constipation | 2 | Gestational diabetes | 0.950 | Gestational diabetes is associated with change in the gut microbiota composition in third trimester of pregnancy and postpartum | Experiment 4: Previous normoglycaemic  pregnancy vs Previous gestational diabetes mellitus (GDM) pregnancy |
| Constipation | 3 | Response to drug | 0.948 | The Parkinson's disease drug entacapone disrupts gut microbiome homoeostasis via iron sequestration | Experiment 10: No drug - Dimethyl sulfoxide only (DMSO only) vs Entacapone High Concentration(ENT-High) after 24 hrs of incubation - Donor 3 |
| Constipation | 4 | Response to drug | 0.948 | The Parkinson's disease drug entacapone disrupts gut microbiome homoeostasis via iron sequestration | Experiment 9: No drug - Dimethyl sulfoxide only (DMSO only) vs Entacapone High Concentration(ENT-High) after 24 hrs of incubation - Donor 2 |
| Constipation | 5 | Population | 0.948 | Oral microbiome homogeneity across diverse human groups from southern Africa: first results from southwestern Angola and Zimbabwe | Experiment 6: Kuvale vs Tshwa |

### Cross-Target Check

| query condition | any other target in top 5 profiles? | matched target conditions |
| --- | --- | --- |
| Multiple sclerosis | no | NA |
| Stroke | no | NA |
| Constipation | no | NA |

Top five condition-centroid neighbors are also written to `target_condition_centroid_neighbors.tsv`; this checks whether the three conditions find each other at condition rather than individual-profile level.

| query condition | rank | centroid neighbor | similarity | other target? |
| --- | ---: | --- | ---: | --- |
| Multiple sclerosis | 1 | Ulcerative colitis | 0.979 | no |
| Multiple sclerosis | 2 | Physical activity | 0.974 | no |
| Multiple sclerosis | 3 | Chronic kidney disease | 0.974 | no |
| Multiple sclerosis | 4 | Obesity | 0.973 | no |
| Multiple sclerosis | 5 | Socioeconomic status | 0.973 | no |
| Stroke | 1 | Age | 0.989 | no |
| Stroke | 2 | Response to immunochemotherapy | 0.986 | no |
| Stroke | 3 | Inflammatory bowel disease | 0.986 | no |
| Stroke | 4 | Cervical cancer | 0.985 | no |
| Stroke | 5 | Polycystic ovary syndrome | 0.985 | no |
| Constipation | 1 | Schizophrenia | 0.991 | no |
| Constipation | 2 | Colorectal cancer | 0.988 | no |
| Constipation | 3 | HIV infection | 0.988 | no |
| Constipation | 4 | Polycystic ovary syndrome | 0.987 | no |
| Constipation | 5 | Obesity | 0.987 | no |

## How Trait Sets Connect To Profiles

A profile is one BugSigDB study/experiment/direction row. Its trait vector counts KG-Microbe trait annotations among organisms reported as increased or decreased in that experiment. A trait set is a named group of biologically related traits, such as `aerobe_v2` or `anaerobe_v2_depletion`. A profile connects to a trait set when the selected trait columns carry signal in that profile; organism support is then checked by intersecting organisms observed in the same profile with organism-trait annotations from `data/derived/taxa_trait_matrix.tsv`.

| trait set | condition | profiles | mean trait-set signal | mean relative trait-set signal |
| --- | --- | ---: | ---: | ---: |
| aerobe_v2 | Constipation | 23 | 104.507 | 0.01764 |
| aerobe_v2 | Multiple sclerosis | 28 | 318.560 | 0.01480 |
| aerobe_v2 | Stroke | 43 | 1259.124 | 0.01902 |
| anaerobe_v2_depletion | Constipation | 23 | 178.196 | 0.06034 |
| anaerobe_v2_depletion | Multiple sclerosis | 28 | 182.714 | 0.06272 |
| anaerobe_v2_depletion | Stroke | 43 | 398.965 | 0.04019 |

## Companion Traits Vs Organisms Heatmap

`traits_vs_organisms_heatmap.png` plots 45 organisms by 34 traits. Rows are organisms observed in the target-condition BugSigDB profiles; columns are binary KG-Microbe trait annotations selected from the top 3 paired-profile contrast traits from `trait-visualization/trait-profile-investigation/trait_profile_paired_contrast_trait_stats.tsv` plus high-signal target-condition traits selected with the original heatmap ranking logic.

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
