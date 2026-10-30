# Taxa Oxygen Switch Investigation

## Scope

- Human-host disease-like profile universe: 1,860 profiles, 545 studies.
- Genus/species taxon-profile occurrences after metadata join: 14,445.
- Oxygen classes are based on direct KG-Microbe trait columns: aerobic side = aerobic, microaerophilic, obligately aerobic, facultatively aerobic; anaerobic side = anaerobic, obligately anaerobic, strictly anaerobic.
- Disease-like filtering excludes Healthy/reference rows, broad cohort/design labels, measurements, exposures, diet/intervention response labels, and other non-disease metadata labels.

## Global Direction Balance by Oxygen Class

| oxygen_class | taxon_profile_count | taxon_count | study_count | increased_taxon_profile_count | decreased_taxon_profile_count | direction_balance |
| --- | --- | --- | --- | --- | --- | --- |
| anaerobic_only | 6488 | 736 | 489 | 3290 | 3198 | 0.014 |
| mixed_oxygen_annotation | 4387 | 367 | 501 | 2325 | 2062 | 0.060 |
| aerobic_only | 1818 | 733 | 280 | 938 | 880 | 0.032 |
| unannotated | 1752 | 972 | 288 | 883 | 869 | 0.008 |

## Strongest Condition-Level Oxygen Switches

Positive scores mean aerobic-only taxa skew increased while anaerobic-only taxa skew decreased.

| condition | body_site_display | aerobe_up_anaerobe_down_score | aerobic_only_direction_balance | anaerobic_only_direction_balance | aerobic_only_taxon_profile_count | anaerobic_only_taxon_profile_count |
| --- | --- | --- | --- | --- | --- | --- |
| Crohn's disease | Feces | 1.108 | 1.000 | -0.108 | 21.0 | 65.0 |
| Polycystic ovary syndrome | Feces | 1.096 | 0.933 | -0.163 | 30.0 | 160.0 |
| Irritable bowel syndrome | Feces | 1.053 | 0.889 | -0.164 | 18.0 | 292.0 |
| Colorectal cancer | Feces | 0.709 | 0.700 | -0.009 | 20.0 | 319.0 |
| COVID-19 | Feces | 0.673 | 0.267 | -0.406 | 30.0 | 155.0 |
| Breast cancer | Feces | 0.615 | 0.538 | -0.077 | 13.0 | 104.0 |
| Alzheimer's disease | Feces | 0.352 | 0.385 | 0.032 | 13.0 | 62.0 |
| Parkinson's disease | Feces | 0.334 | 0.429 | 0.095 | 14.0 | 475.0 |

Negative scores mean the reciprocal pattern: anaerobic-only taxa skew increased while aerobic-only taxa skew decreased.

| condition | body_site_display | aerobe_up_anaerobe_down_score | aerobic_only_direction_balance | anaerobic_only_direction_balance | aerobic_only_taxon_profile_count | anaerobic_only_taxon_profile_count |
| --- | --- | --- | --- | --- | --- | --- |
| Rheumatoid arthritis | Oral cavity / dental | -1.109 | -0.625 | 0.484 | 32.0 | 93.0 |
| Type II diabetes mellitus | Oral cavity / dental | -0.842 | -0.579 | 0.263 | 19.0 | 19.0 |
| Glaucoma | Ocular / eye | -0.804 | -0.200 | 0.604 | 50.0 | 101.0 |
| Periodontitis | Oral cavity / dental | -0.750 | -0.375 | 0.375 | 16.0 | 32.0 |
| Gestational diabetes | Feces | -0.745 | -0.778 | -0.033 | 72.0 | 244.0 |
| Head and neck squamous cell carcinoma | Oral cavity / dental | -0.630 | -0.600 | 0.030 | 15.0 | 33.0 |
| Glaucoma | Oral cavity / dental | -0.625 | -0.833 | -0.208 | 48.0 | 48.0 |
| Human papilloma virus infection | Vaginal / reproductive | -0.404 | 0.267 | 0.671 | 30.0 | 85.0 |

## Top Cross-Disease Taxon Switches

These taxa have at least one disease condition where they are predominantly reported increased and at least one where they are predominantly reported decreased.

| taxon_name | taxonomic_rank | oxygen_class | switch_score | increased_dominant_condition_count | decreased_dominant_condition_count | top_increased_conditions | top_decreased_conditions |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Bacteroides | genus | mixed_oxygen_annotation | 90.66 | 7 | 9 | Chronic constipation (7 profiles, 4 studies); Obesity (7 profiles, 3 studies); Type II diabetes mellitus (5 profiles, 2 studies); Irritable bowel syndrome (4 profiles, 2 studies) | HIV infection (7 profiles, 7 studies); Parkinson's disease (6 profiles, 6 studies); Colorectal cancer (6 profiles, 4 studies); Diarrhea (6 profiles, 2 studies) |
| Prevotella | genus | anaerobic_only | 88.72 | 7 | 8 | HIV infection (6 profiles, 6 studies); Human papilloma virus infection (6 profiles, 4 studies); Obesity (5 profiles, 5 studies); Alzheimer's disease (4 profiles, 3 studies) | Irritable bowel syndrome (9 profiles, 4 studies); COVID-19 (4 profiles, 4 studies); Colorectal cancer (4 profiles, 4 studies); Gestational diabetes (4 profiles, 2 studies) |
| Bifidobacterium | genus | mixed_oxygen_annotation | 62.17 | 3 | 6 | Parkinson's disease (19 profiles, 11 studies); Gestational diabetes (5 profiles, 2 studies); Gastric adenocarcinoma (3 profiles, 2 studies) | Irritable bowel syndrome (11 profiles, 6 studies); HIV infection (5 profiles, 5 studies); Diarrhea (5 profiles, 2 studies); Colorectal cancer (4 profiles, 4 studies) |
| Lactobacillus | genus | mixed_oxygen_annotation | 45.56 | 6 | 4 | Parkinson's disease (10 profiles, 7 studies); Obesity (6 profiles, 3 studies); COVID-19 (3 profiles, 3 studies); HIV infection (3 profiles, 3 studies) | Human papilloma virus infection (5 profiles, 4 studies); Gestational diabetes (5 profiles, 2 studies); Human immunodeficiency virus (5 profiles, 2 studies); Chronic constipation (4 profiles, 3 studies) |
| Ruminococcus | genus | anaerobic_only | 44.73 | 4 | 7 | Obesity (6 profiles, 2 studies); Gestational diabetes (5 profiles, 2 studies); Colorectal cancer (4 profiles, 4 studies); Type II diabetes mellitus (3 profiles, 2 studies) | Irritable bowel syndrome (8 profiles, 2 studies); COVID-19 (5 profiles, 4 studies); Parkinson's disease (5 profiles, 4 studies); HIV infection (4 profiles, 4 studies) |
| Veillonella | genus | anaerobic_only | 35.97 | 7 | 3 | COVID-19 (10 profiles, 5 studies); Schizophrenia (5 profiles, 4 studies); Type II diabetes mellitus (5 profiles, 2 studies); Irritable bowel syndrome (4 profiles, 3 studies) | Gestational diabetes (8 profiles, 3 studies); Food allergy (4 profiles, 3 studies); Colorectal cancer (3 profiles, 3 studies) |
| Faecalibacterium | genus | mixed_oxygen_annotation | 32.50 | 3 | 11 | Glaucoma (5 profiles, 3 studies); Chronic constipation (4 profiles, 2 studies); Obesity (3 profiles, 2 studies) | Parkinson's disease (25 profiles, 12 studies); Irritable bowel syndrome (12 profiles, 5 studies); Colorectal cancer (7 profiles, 7 studies); HIV infection (6 profiles, 5 studies) |
| Blautia | genus | anaerobic_only | 29.93 | 3 | 6 | Gestational diabetes (6 profiles, 2 studies); Glaucoma (4 profiles, 3 studies); Polycystic ovary syndrome (3 profiles, 3 studies) | Parkinson's disease (12 profiles, 10 studies); Colorectal cancer (8 profiles, 7 studies); COVID-19 (5 profiles, 3 studies); Schizophrenia (3 profiles, 3 studies) |
| Roseburia | genus | anaerobic_only | 29.03 | 2 | 11 | Gestational diabetes (6 profiles, 2 studies); Glaucoma (5 profiles, 2 studies) | Parkinson's disease (29 profiles, 17 studies); Colorectal cancer (6 profiles, 6 studies); COVID-19 (5 profiles, 3 studies); Chronic kidney disease (3 profiles, 3 studies) |
| Lachnospira | genus | anaerobic_only | 28.56 | 3 | 5 | Polycystic ovary syndrome (6 profiles, 4 studies); Obesity (4 profiles, 2 studies); Chronic constipation (3 profiles, 3 studies) | COVID-19 (5 profiles, 3 studies); Parkinson's disease (5 profiles, 3 studies); Colorectal cancer (4 profiles, 4 studies); HIV infection (3 profiles, 3 studies) |

## Within-Condition Direction Conflicts

34 condition/taxon rows have at least 3 increased and 3 decreased profile occurrences. These are not necessarily contradictions; many combine different experiments, cohorts, or body sites under the same condition label.

## Files

- `global_oxygen_direction_summary.tsv`
- `condition_oxygen_direction_summary.tsv`
- `condition_oxygen_switch_scores.tsv`
- `taxon_condition_direction_summary.tsv`
- `taxa_cross_disease_switch_summary.tsv`
- `taxa_within_condition_direction_conflicts.tsv`
- `condition_oxygen_switch_scores.png`
