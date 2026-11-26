# Enrichment Infographic Follow-Up

Historical source figure: `data/derived/enrichment_infographic.png` in
`CultureBotAI/microbial-trait-profiling`. This is a historical repository-relative
reference, not a claim that the historical figure is present in the current checkout.

## Report variants

Paths below are relative to `trait-visualization/trait-profile-investigation/enrichment-infographic-followup/`.
These saved runs are not interchangeable; retain their results and read each run's
condition-selection and heatmap-selection paragraphs before comparing them.

- `README.md`: unprefixed report in this directory; its saved run uses infographic
  condition selection and infographic traits plus high-signal traits.
- `my_top3_README.md`: paired-condition and paired-trait run with `my_top3_`
  prefixed artifacts.
- `my_top3_conditions_relevant_traits_README.md`: paired-condition and paired-trait
  run with `my_top3_conditions_relevant_traits_` prefixed artifacts.
- `my-top3-enriched-traits/README.md`: paired-condition and paired-trait run in a
  separate output directory with unprefixed artifacts.

The two prefixed reports and the subdirectory report preserve separately named
outputs of similar selections; they are not additional independent evidence.
The generator uses this common provenance format; run-specific selections and
artifact names are recorded below.

The target conditions for this run are selected from the top 3 paired-profile condition contrast rows from `trait-visualization/trait-profile-investigation/trait_profile_paired_contrast_condition_summary.tsv`. The local follow-up uses the full `data/derived/study_trait_profiles_wide.tsv` matrix to list titles, find similar profiles, and connect enriched traits to organisms. Profiles carry their own experiment's condition, body site and host from `trait-visualization/study_trait_profile_study_metadata.tsv`; experiment details and reported taxa come from the pinned `data/raw/full_dump.csv`, and only profiles in this corpus are counted.

## Target Conditions

| condition | paired absolute contrast load | local profiles | studies | experiments | local mean nonzero traits/profile | local percentile | local outlier? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Multiple sclerosis | 0.94 | 28 | 6 | 14 | 101.5 | 42.6% | no |
| Stroke | 0.92 | 43 | 5 | 28 | 140.1 | 69.3% | no |
| Helminthiasis | 0.76 | 27 | 1 | 15 | 223.7 | 91.4% | no |

Interpretation: condition support matters because profile counts, study counts, and body-site mix vary across exact BugSigDB condition labels.

## Study And Experiment Titles

| condition | study | year | experiment | profiles |
| --- | --- | ---: | --- | ---: |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 1: Zero Soil Transmitted Helminth infection (Zero STH infection) vs High Soil Transmitted Helminth infection (High STH infection) | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 10: 2008 - Uninfected vs 2010 - Uninfected | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 11: 2010 - Treated with albendazole vs 2008 - Infected | 1 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 12: 2010 - Treated with albendazole and dewormed vs 2008 - Infected | 1 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 13: 2010 - Self cleared individuals vs 2008 - Infected | 1 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 14: Uninfected (Zero STH infection) vs Heavy Infected (Ascaris) | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 15: Uninfected (Zero STH infection) vs Heavy Infected (Necator) | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 16: Uninfected (Zero STH infection) vs Heavy Infected (Trichuris) | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 2: Zero Soil Transmitted Helminth infection (Zero STH infection) vs High Soil Transmitted Helminth infection (High STH infection) | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 3: Zero Soil Transmitted Helminth infection (Zero STH infection) vs High Soil Transmitted Helminth infection (High STH infection) | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 5: Zero Soil Transmitted Helminth infection (Zero STH infection) vs High Soil Transmitted Helminth infection (High STH infection) | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 6: Zero Soil Transmitted Helminth infection (Zero STH infection) vs High Soil Transmitted Helminth infection (High STH infection) | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 7: 2008 - Uninfected vs 2010 - Uninfected | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 8: 2008 - Uninfected vs 2010 - Uninfected | 2 |
| Helminthiasis | Differential human gut microbiome assemblages during soil-transmitted helminth infections in Indonesia and Liberia (bsdb:29486796) | 2018 | Experiment 9: 2008 - Uninfected vs 2010 - Uninfected | 2 |
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
| Multiple sclerosis | 1 | Chlamydia trachomatis | 0.935 | Cervicovaginal microbiome and natural history of Chlamydia trachomatis in adolescents and young women | Experiment 1: Chlamydia trichomonas Controls (t0) vs Chlamydia trichomonas Cases (t0) |
| Multiple sclerosis | 2 | HIV infection | 0.932 | Longitudinal gut microbiota composition of South African and Nigerian infants in relation to tetanus vaccine responses | Experiment 3: Infants who are HIV unexposed and uninfected (iHUU) at week 15 vs Infants who are HIV exposed but uninfected (iHEU) at week 15 |
| Multiple sclerosis | 3 | Root caries | 0.928 | Comparative analysis of microbiome in coronal and root caries | Experiment 1: Coronal Caries group vs Root Caries group |
| Multiple sclerosis | 4 | Endometriosis | 0.928 | The Endobiota Study: Comparison of Vaginal, Cervical and Gut Microbiota Between Women with Stage 3/4 Endometriosis and Healthy Controls | Experiment 1: Healthy controls vs Stage 3/4 endometriosis patients |
| Multiple sclerosis | 5 | Adverse effect | 0.927 | Intestinal microbiota influences clinical outcome and side effects of early breast cancer treatment | Experiment 13: No - increased BMI (Post-Chemo) vs Yes - increased BMI (Post-Chemo) |
| Stroke | 1 | Response to drug | 0.958 | The Parkinson's disease drug entacapone disrupts gut microbiome homoeostasis via iron sequestration | Experiment 10: No drug - Dimethyl sulfoxide only (DMSO only) vs Entacapone High Concentration(ENT-High) after 24 hrs of incubation - Donor 3 |
| Stroke | 2 | Response to drug | 0.958 | The Parkinson's disease drug entacapone disrupts gut microbiome homoeostasis via iron sequestration | Experiment 9: No drug - Dimethyl sulfoxide only (DMSO only) vs Entacapone High Concentration(ENT-High) after 24 hrs of incubation - Donor 2 |
| Stroke | 3 | HIV infection | 0.955 | Impact of Acute HIV Infection and Early Antiretroviral Therapy on the Human Gut Microbiome | Experiment 1: Human immunodeficiency virus negative (HIV-)Uninfected Controls vs Human immunodeficiency virus positive(HIV+) pre - antiretroviral therapy(ART) |
| Stroke | 4 | Chronic hepatitis B virus infection | 0.955 | Gut Microbial Dysbiosis Is Associated with Altered Hepatic Functions and Serum Metabolites in Chronic Hepatitis B Patients | Experiment 1: healthy controls vs chronic hepatitis B patients |
| Stroke | 5 | Parkinson's disease | 0.955 | Gut microbiota in Parkinson's disease: Temporal stability and relations to disease progression | Experiment 4: Tremor-dominant phenotype (TD) at baseline vs Postural instability and gait difficulty (PIGD) phenotype at baseline |
| Helminthiasis | 1 | Cardiovascular disease | 0.992 | Sex-specific differences in intestinal microbiota associated with cardiovascular diseases | Experiment 4: Women without cardiovascular disease vs Women with coronary heart disease (CHD) |
| Helminthiasis | 2 | Chronic constipation | 0.990 | Mucosa-associated gut microbiome in Japanese patients with functional constipation | Experiment 2: Healthy controls vs Patients with Functional constipation (FC) |
| Helminthiasis | 3 | HIV infection | 0.990 | The semen microbiome and its relationship with local immunology and viral load in HIV infection | Experiment 1: HIV-uninfected men vs HIV-infected men |
| Helminthiasis | 4 | Age-related macular degeneration | 0.988 | Beneficial Effects of Micronutrient Supplementation in Restoring the Altered Microbiota and Gut-Retina Axis in Patients with Neovascular Age-Related Macular Degeneration-A Randomized Clinical Trial | Experiment 1: Healthy controls vs nAMD (neovascular Age-related macular degeneration) patients |
| Helminthiasis | 5 | Gastric cancer | 0.988 | Association Between Gut Microbiota and <i>Helicobacter pylori</i>-Related Gastric Lesions in a High-Risk Population of Gastric Cancer | Experiment 2: past H.pylori infection vs current H.pylori infection patients |

### Cross-Target Check

| query condition | any other target in top 5 profiles? | matched target conditions |
| --- | --- | --- |
| Multiple sclerosis | no | NA |
| Stroke | no | NA |
| Helminthiasis | no | NA |

Top five condition-centroid neighbors are also written to `target_condition_centroid_neighbors.tsv`; this checks whether the three conditions find each other at condition rather than individual-profile level.

| query condition | rank | centroid neighbor | similarity | other target? |
| --- | ---: | --- | ---: | --- |
| Multiple sclerosis | 1 | Ulcerative colitis | 0.981 | no |
| Multiple sclerosis | 2 | Psoriasis | 0.977 | no |
| Multiple sclerosis | 3 | Chronic kidney disease | 0.977 | no |
| Multiple sclerosis | 4 | Physical activity | 0.976 | no |
| Multiple sclerosis | 5 | Obesity | 0.976 | no |
| Stroke | 1 | Age | 0.990 | no |
| Stroke | 2 | Cervical cancer | 0.987 | no |
| Stroke | 3 | Response to immunochemotherapy | 0.987 | no |
| Stroke | 4 | Parkinson's disease | 0.986 | no |
| Stroke | 5 | Polycystic ovary syndrome | 0.986 | no |
| Helminthiasis | 1 | Air pollution | 0.994 | no |
| Helminthiasis | 2 | Thyroid carcinoma | 0.993 | no |
| Helminthiasis | 3 | Chronic kidney disease,Stage 5 chronic kidney disease | 0.993 | no |
| Helminthiasis | 4 | Peutz-Jeghers Polyp | 0.993 | no |
| Helminthiasis | 5 | Constipation,Parkinson's disease | 0.993 | no |

## How Trait Sets Connect To Profiles

A profile is one BugSigDB study/experiment/direction row. Its trait vector counts KG-Microbe trait annotations among organisms reported as increased or decreased in that experiment. A trait set is a named group of biologically related traits, such as `aerobe_v2` or `anaerobe_v2_depletion`. A profile connects to a trait set when the selected trait columns carry signal in that profile; organism support is then checked by intersecting organisms observed in the same profile with organism-trait annotations from `data/derived/taxa_trait_matrix.tsv`.

| trait set | condition | profiles | mean trait-set signal | mean relative trait-set signal |
| --- | --- | ---: | ---: | ---: |
| aerobe_v2 | Multiple sclerosis | 28 | 421.607 | 0.01390 |
| aerobe_v2 | Helminthiasis | 27 | 1861.500 | 0.03166 |
| aerobe_v2 | Stroke | 43 | 1154.360 | 0.01637 |
| anaerobe_v2_depletion | Multiple sclerosis | 28 | 207.893 | 0.05752 |
| anaerobe_v2_depletion | Helminthiasis | 27 | 909.667 | 0.01614 |
| anaerobe_v2_depletion | Stroke | 43 | 444.012 | 0.04248 |

## Companion Traits Vs Organisms Heatmap

`my_top3_traits_vs_organisms_heatmap.png` plots 45 organisms by 34 traits. Rows are organisms observed in the target-condition BugSigDB profiles; columns are binary KG-Microbe trait annotations selected from the top 3 paired-profile contrast traits from `trait-visualization/trait-profile-investigation/trait_profile_paired_contrast_trait_stats.tsv` plus high-signal target-condition traits selected with the original heatmap ranking logic.

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
