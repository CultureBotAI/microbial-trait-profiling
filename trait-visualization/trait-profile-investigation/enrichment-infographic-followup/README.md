# Enrichment Infographic Follow-Up

Source figure: https://github.com/CultureBotAI/bsdbng/blob/enrichment-analysis-heatmap-cluster/data/derived/enrichment_infographic.png

The target conditions for this run are selected from the infographic top three condition-level outliers. The local follow-up uses the full `data/derived/study_trait_profiles_wide.tsv` matrix to list titles, find similar profiles, and connect enriched traits to organisms. Profiles carry their own experiment's condition, body site and host from `trait-visualization/study_trait_profile_study_metadata.tsv`; experiment details and reported taxa come from the pinned `data/raw/full_dump.csv`, and only profiles in this corpus are counted.

## Target Conditions

| condition | infographic avg enrichment depth | local profiles | studies | experiments | local mean nonzero traits/profile | local percentile | local outlier? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Sleep Disorder | 216.00 | 2 | 1 | 1 | 248.0 | 90.0% | no |
| Head and neck squamous cell carcinoma | 210.33 | 10 | 2 | 5 | 82.9 | 15.2% | no |

Target conditions with fewer than 2 profiles in this corpus are not summarized or queried: Thalassemia (0 profiles).

Interpretation: condition support matters because profile counts, study counts, and body-site mix vary across exact BugSigDB condition labels.

## Study And Experiment Titles

| condition | study | year | experiment | profiles |
| --- | --- | ---: | --- | ---: |
| Head and neck squamous cell carcinoma | Association of Oral Microbiome With Risk for Incident Head and Neck Squamous Cell Cancer (bsdb:29327043) | 2018 | Experiment 1: controls vs HNSCC | 2 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 1: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 2 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 2: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 2 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 3: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 2 |
| Head and neck squamous cell carcinoma | Oral Microbiome and Subsequent Risk of Head and Neck Squamous Cell Cancer (bsdb:39325441) | 2024 | Experiment 4: Healthy Controls vs Head and neck squamous cell cancer (HNSCC) participants | 2 |
| Sleep Disorder | Preliminary analysis of salivary microbiota in catathrenia (nocturnal groaning) using machine learning algorithms (bsdb:40247863) | 2025 | Experiment 1: Healthy control group vs Catathrenia group | 2 |

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
| Sleep Disorder | 1 | HIV infection | 0.993 | The semen microbiome and its relationship with local immunology and viral load in HIV infection | Experiment 1: HIV-uninfected men vs HIV-infected men |
| Sleep Disorder | 2 | Cardiovascular disease | 0.992 | Sex-specific differences in intestinal microbiota associated with cardiovascular diseases | Experiment 4: Women without cardiovascular disease vs Women with coronary heart disease (CHD) |
| Sleep Disorder | 3 | Helicobacter pylori | 0.991 | Microbial Diversity and Composition in Six Different Gastrointestinal Sites among Participants Undergoing Upper Gastrointestinal Endoscopy in Henan, China | Experiment 7: H. pylori-negative cardia group (NC) vs H. pylori-positive cardia group (PC) |
| Sleep Disorder | 4 | Chronic constipation | 0.991 | Mucosa-associated gut microbiome in Japanese patients with functional constipation | Experiment 2: Healthy controls vs Patients with Functional constipation (FC) |
| Sleep Disorder | 5 | Colorectal cancer | 0.991 | Alteration of fecal tryptophan metabolism correlates with shifted microbiota and may be involved in pathogenesis of colorectal cancer | Experiment 1: Healthy Controls (HC) vs CRC |
| Head and neck squamous cell carcinoma | 1 | Gingivitis | 0.946 | Quantitative analysis of the effects of essential oil mouthrinses on clinical plaque microbiome: a parallel-group, randomized trial | Experiment 1: Healthy subjects vs Gingivitis subjects |
| Head and neck squamous cell carcinoma | 2 | Multiple myeloma | 0.940 | Alterations of gut microbiome accelerate multiple myeloma progression by increasing the relative abundances of nitrogen-recycling bacteria | Experiment 1: healthy control vs multiple myeloma patients |
| Head and neck squamous cell carcinoma | 3 | Colorectal cancer | 0.939 | The oral microbiota in colorectal cancer is distinctive and predictive | Experiment 1: Healthy Controls vs Individuals with CRC |
| Head and neck squamous cell carcinoma | 4 | Chlamydia trachomatis | 0.937 | Cervicovaginal microbiome and natural history of Chlamydia trachomatis in adolescents and young women | Experiment 1: Chlamydia trichomonas Controls (t0) vs Chlamydia trichomonas Cases (t0) |
| Head and neck squamous cell carcinoma | 5 | Psoriasis | 0.929 | Multiomic Analysis of the Gut Microbiome in Psoriasis Reveals Distinct Host‒Microbe Associations | Experiment 1: Healthy controls vs Psoriasis patients (PSO) |

### Cross-Target Check

| query condition | any other target in top 5 profiles? | matched target conditions |
| --- | --- | --- |
| Sleep Disorder | no | NA |
| Head and neck squamous cell carcinoma | no | NA |

Top five condition-centroid neighbors are also written to `target_condition_centroid_neighbors.tsv`; this checks whether the three conditions find each other at condition rather than individual-profile level.

| query condition | rank | centroid neighbor | similarity | other target? |
| --- | ---: | --- | ---: | --- |
| Sleep Disorder | 1 | Autoimmune type 1 diabetes | 0.996 | no |
| Sleep Disorder | 2 | Peptic esophagitis | 0.995 | no |
| Sleep Disorder | 3 | Thyroid carcinoma | 0.995 | no |
| Sleep Disorder | 4 | Hepatitis | 0.994 | no |
| Sleep Disorder | 5 | Blastocystis hominis infectious disease | 0.992 | no |
| Head and neck squamous cell carcinoma | 1 | Squamous cell carcinoma | 0.984 | no |
| Head and neck squamous cell carcinoma | 2 | Asthma | 0.973 | no |
| Head and neck squamous cell carcinoma | 3 | Oral squamous cell carcinoma | 0.964 | no |
| Head and neck squamous cell carcinoma | 4 | Phenylketonuria | 0.960 | no |
| Head and neck squamous cell carcinoma | 5 | Chronic periodontitis | 0.957 | no |

## How Trait Sets Connect To Profiles

A profile is one BugSigDB study/experiment/direction row. Its trait vector counts KG-Microbe trait annotations among organisms reported as increased or decreased in that experiment. A trait set is a named group of biologically related traits, such as `aerobe_v2` or `anaerobe_v2_depletion`. A profile connects to a trait set when the selected trait columns carry signal in that profile; organism support is then checked by intersecting organisms observed in the same profile with organism-trait annotations from `data/derived/taxa_trait_matrix.tsv`.

| trait set | condition | profiles | mean trait-set signal | mean relative trait-set signal |
| --- | --- | ---: | ---: | ---: |
| aerobe_v2 | Head and neck squamous cell carcinoma | 10 | 186.367 | 0.03821 |
| aerobe_v2 | Sleep Disorder | 2 | 2601.167 | 0.02128 |
| anaerobe_v2_depletion | Head and neck squamous cell carcinoma | 10 | 41.850 | 0.02807 |
| anaerobe_v2_depletion | Sleep Disorder | 2 | 860.500 | 0.01533 |

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
