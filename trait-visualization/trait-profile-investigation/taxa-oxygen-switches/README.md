# Taxa Oxygen Switch Investigation

## Scope

- Human-host disease-like profile universe: 4,102 profiles, 860 studies.
- Genus/species taxon-profile occurrences after metadata join: 27,057.
- Oxygen classes are based on direct KG-Microbe trait columns: aerobic side = aerobic, microaerophilic, obligately aerobic, facultatively aerobic; anaerobic side = anaerobic, obligately anaerobic, strictly anaerobic.
- Disease-like filtering excludes Healthy/reference rows, broad cohort/design labels, measurements, exposures, diet/intervention response labels, and other non-disease metadata labels.

## Global Direction Balance by Oxygen Class

| oxygen_class | taxon_profile_count | taxon_count | study_count | increased_taxon_profile_count | decreased_taxon_profile_count | direction_balance |
| --- | --- | --- | --- | --- | --- | --- |
| anaerobic_only | 11717 | 853 | 755 | 6019 | 5698 | 0.027 |
| mixed_oxygen_annotation | 7852 | 359 | 772 | 4293 | 3559 | 0.093 |
| aerobic_only | 3951 | 1026 | 491 | 2037 | 1914 | 0.031 |
| unannotated | 3537 | 1579 | 465 | 1837 | 1700 | 0.039 |

## Strongest Condition-Level Oxygen Switches

Positive scores mean aerobic-only taxa skew increased while anaerobic-only taxa skew decreased.

| condition | body_site_display | aerobe_up_anaerobe_down_score | aerobic_only_direction_balance | anaerobic_only_direction_balance | aerobic_only_taxon_profile_count | anaerobic_only_taxon_profile_count |
| --- | --- | --- | --- | --- | --- | --- |
| Multiple myeloma | Feces | 1.071 | 0.833 | -0.238 | 12.0 | 42.0 |
| Polycystic ovary syndrome | Feces | 0.987 | 0.892 | -0.095 | 37.0 | 190.0 |
| COVID-19 | Lower airway / lung | 0.967 | 0.795 | -0.172 | 39.0 | 29.0 |
| Crohn's disease | Feces | 0.903 | 0.586 | -0.317 | 29.0 | 120.0 |
| Irritable bowel syndrome | Feces | 0.798 | 0.652 | -0.146 | 23.0 | 349.0 |
| COVID-19 | Feces | 0.532 | 0.214 | -0.318 | 56.0 | 343.0 |
| Colorectal cancer | Feces | 0.418 | 0.455 | 0.037 | 44.0 | 517.0 |
| Alzheimer's disease | Feces | 0.401 | 0.385 | -0.016 | 13.0 | 61.0 |

Negative scores mean the reciprocal pattern: anaerobic-only taxa skew increased while aerobic-only taxa skew decreased.

| condition | body_site_display | aerobe_up_anaerobe_down_score | aerobic_only_direction_balance | anaerobic_only_direction_balance | aerobic_only_taxon_profile_count | anaerobic_only_taxon_profile_count |
| --- | --- | --- | --- | --- | --- | --- |
| Colorectal cancer | Intestinal / GI mucosa | -1.048 | -0.714 | 0.333 | 21.0 | 12.0 |
| Colorectal cancer | Lower GI tract | -0.976 | -0.538 | 0.437 | 26.0 | 167.0 |
| Rheumatoid arthritis | Oral cavity / dental | -0.862 | -0.491 | 0.371 | 55.0 | 124.0 |
| Human immunodeficiency virus | Feces | -0.795 | -0.875 | -0.080 | 16.0 | 113.0 |
| Gestational diabetes | Feces | -0.761 | -0.806 | -0.044 | 72.0 | 272.0 |
| Periodontitis | Oral cavity / dental | -0.717 | -0.250 | 0.467 | 32.0 | 75.0 |
| Glaucoma | Ocular / eye | -0.709 | -0.075 | 0.635 | 67.0 | 104.0 |
| Type II diabetes mellitus | Oral cavity / dental | -0.500 | -0.417 | 0.083 | 24.0 | 24.0 |

## Top Cross-Disease Taxon Switches

These taxa have at least one disease condition where they are predominantly reported increased and at least one where they are predominantly reported decreased.

| taxon_name | taxonomic_rank | oxygen_class | switch_score | increased_dominant_condition_count | decreased_dominant_condition_count | top_increased_conditions | top_decreased_conditions |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Bacteroides | genus | mixed_oxygen_annotation | 179.28 | 10 | 11 | Obesity (15 profiles, 5 studies); Chronic constipation (12 profiles, 5 studies); Gestational diabetes (7 profiles, 3 studies); Type II diabetes mellitus (5 profiles, 2 studies) | HIV infection (18 profiles, 10 studies); Diarrhea (8 profiles, 2 studies); Parkinson's disease (6 profiles, 6 studies); Gastric cancer (5 profiles, 3 studies) |
| Prevotella | genus | anaerobic_only | 126.73 | 13 | 8 | HIV infection (12 profiles, 10 studies); Colorectal cancer (10 profiles, 7 studies); Obesity (9 profiles, 5 studies); Rheumatoid arthritis (8 profiles, 4 studies) | Irritable bowel syndrome (10 profiles, 5 studies); COVID-19 (9 profiles, 9 studies); Ulcerative colitis (5 profiles, 4 studies); Gestational diabetes (4 profiles, 2 studies) |
| Lactobacillus | genus | mixed_oxygen_annotation | 110.50 | 9 | 7 | Parkinson's disease (16 profiles, 13 studies); Gastric cancer (8 profiles, 3 studies); Endometriosis (6 profiles, 5 studies); HIV infection (6 profiles, 4 studies) | Human papilloma virus infection (8 profiles, 5 studies); Chronic constipation (7 profiles, 4 studies); Cervical glandular intraepithelial neoplasia (6 profiles, 4 studies); Human immunodeficiency virus (6 profiles, 2 studies) |
| Bifidobacterium | genus | mixed_oxygen_annotation | 105.16 | 5 | 7 | Parkinson's disease (26 profiles, 14 studies); Gestational diabetes (6 profiles, 3 studies); Gastric cancer (3 profiles, 3 studies); Ulcerative colitis (3 profiles, 3 studies) | Irritable bowel syndrome (11 profiles, 6 studies); Eczema (8 profiles, 2 studies); HIV infection (6 profiles, 5 studies); COVID-19 (6 profiles, 3 studies) |
| Streptococcus | genus | mixed_oxygen_annotation | 74.94 | 20 | 5 | Gastric cancer (14 profiles, 7 studies); Colorectal cancer (12 profiles, 9 studies); Parkinson's disease (10 profiles, 9 studies); Obesity (10 profiles, 4 studies) | COVID-19 (8 profiles, 7 studies); Chronic constipation (4 profiles, 3 studies); Eczema (4 profiles, 3 studies); Chlamydia trachomatis (4 profiles, 2 studies) |
| Ruminococcus | genus | anaerobic_only | 68.62 | 4 | 9 | Obesity (14 profiles, 3 studies); Gestational diabetes (5 profiles, 2 studies); Chronic constipation (4 profiles, 4 studies); Type II diabetes mellitus (3 profiles, 2 studies) | COVID-19 (9 profiles, 8 studies); Parkinson's disease (9 profiles, 6 studies); Irritable bowel syndrome (9 profiles, 2 studies); Colorectal cancer (5 profiles, 5 studies) |
| Veillonella | genus | anaerobic_only | 59.58 | 10 | 4 | COVID-19 (12 profiles, 6 studies); HIV infection (5 profiles, 5 studies); Parkinson's disease (5 profiles, 4 studies); Schizophrenia (5 profiles, 4 studies) | Gestational diabetes (8 profiles, 3 studies); Colorectal cancer (7 profiles, 7 studies); Food allergy (4 profiles, 3 studies); Human immunodeficiency virus (3 profiles, 2 studies) |
| Blautia | genus | anaerobic_only | 50.14 | 5 | 8 | Gestational diabetes (6 profiles, 2 studies); Glaucoma (4 profiles, 3 studies); Endometriosis (3 profiles, 3 studies); HIV infection (3 profiles, 3 studies) | Parkinson's disease (17 profiles, 14 studies); Colorectal cancer (10 profiles, 8 studies); COVID-19 (9 profiles, 7 studies); Graft versus host disease (4 profiles, 3 studies) |
| Roseburia | genus | anaerobic_only | 50.06 | 3 | 15 | Gestational diabetes (8 profiles, 4 studies); Glaucoma (5 profiles, 2 studies); Irritable bowel syndrome (4 profiles, 2 studies) | Parkinson's disease (33 profiles, 20 studies); Stroke (11 profiles, 3 studies); COVID-19 (9 profiles, 7 studies); Colorectal cancer (8 profiles, 8 studies) |
| Clostridium | genus | mixed_oxygen_annotation | 45.56 | 6 | 4 | Cervical glandular intraepithelial neoplasia (4 profiles, 2 studies); Chronic constipation (3 profiles, 3 studies); Colorectal adenoma (3 profiles, 2 studies); Polycystic ovary syndrome (3 profiles, 2 studies) | Colorectal cancer (11 profiles, 7 studies); Parkinson's disease (10 profiles, 7 studies); COVID-19 (8 profiles, 4 studies); HIV infection (7 profiles, 5 studies) |

## Within-Condition Direction Conflicts

200 condition/taxon rows have at least 3 increased and 3 decreased profile occurrences. These are not necessarily contradictions; many combine different experiments, cohorts, or body sites under the same condition label.

## Files

- `global_oxygen_direction_summary.tsv`
- `condition_oxygen_direction_summary.tsv`
- `condition_oxygen_switch_scores.tsv`
- `taxon_condition_direction_summary.tsv`
- `taxa_cross_disease_switch_summary.tsv`
- `taxa_within_condition_direction_conflicts.tsv`
- `condition_oxygen_switch_scores.png`
