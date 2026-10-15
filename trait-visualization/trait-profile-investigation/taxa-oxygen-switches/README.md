# Taxa Oxygen Switch Investigation

## Scope

- Human-host disease-like profile universe: 4,373 profiles, 882 studies.
- Genus/species taxon-profile occurrences after metadata join: 28,325.
- Oxygen classes are based on direct KG-Microbe trait columns: aerobic side = aerobic, microaerophilic, obligately aerobic, facultatively aerobic; anaerobic side = anaerobic, obligately anaerobic, strictly anaerobic.
- Disease-like filtering excludes Healthy/reference rows, broad cohort/design labels, measurements, exposures, diet/intervention response labels, and other non-disease metadata labels.

## Global Direction Balance by Oxygen Class

| oxygen_class | taxon_profile_count | taxon_count | study_count | increased_taxon_profile_count | decreased_taxon_profile_count | direction_balance |
| --- | --- | --- | --- | --- | --- | --- |
| anaerobic_only | 11863 | 847 | 774 | 6040 | 5823 | 0.018 |
| mixed_oxygen_annotation | 8126 | 366 | 789 | 4457 | 3669 | 0.097 |
| unannotated | 4214 | 1903 | 482 | 2125 | 2089 | 0.009 |
| aerobic_only | 4122 | 1059 | 506 | 2158 | 1964 | 0.047 |

## Strongest Condition-Level Oxygen Switches

Positive scores mean aerobic-only taxa skew increased while anaerobic-only taxa skew decreased.

| condition | body_site_display | aerobe_up_anaerobe_down_score | aerobic_only_direction_balance | anaerobic_only_direction_balance | aerobic_only_taxon_profile_count | anaerobic_only_taxon_profile_count |
| --- | --- | --- | --- | --- | --- | --- |
| Polycystic ovary syndrome | Feces | 1.183 | 0.943 | -0.240 | 35.0 | 179.0 |
| Ulcerative colitis | Feces | 1.128 | 0.730 | -0.398 | 37.0 | 133.0 |
| Multiple myeloma | Feces | 1.071 | 0.833 | -0.238 | 12.0 | 42.0 |
| COVID-19 | Lower airway / lung | 0.967 | 0.795 | -0.172 | 39.0 | 29.0 |
| Irritable bowel syndrome | Feces | 0.806 | 0.652 | -0.153 | 23.0 | 352.0 |
| Colorectal cancer | Feces | 0.683 | 0.783 | 0.101 | 120.0 | 536.0 |
| COVID-19 | Feces | 0.495 | 0.132 | -0.363 | 53.0 | 336.0 |
| HIV infection | Feces | 0.400 | 0.524 | 0.124 | 21.0 | 226.0 |

Negative scores mean the reciprocal pattern: anaerobic-only taxa skew increased while aerobic-only taxa skew decreased.

| condition | body_site_display | aerobe_up_anaerobe_down_score | aerobic_only_direction_balance | anaerobic_only_direction_balance | aerobic_only_taxon_profile_count | anaerobic_only_taxon_profile_count |
| --- | --- | --- | --- | --- | --- | --- |
| Rheumatoid arthritis | Oral cavity / dental | -1.012 | -0.571 | 0.441 | 56.0 | 118.0 |
| Colorectal cancer | Lower GI tract | -0.793 | -0.429 | 0.365 | 28.0 | 192.0 |
| Glaucoma | Ocular / eye | -0.788 | -0.184 | 0.604 | 49.0 | 91.0 |
| Gestational diabetes | Feces | -0.737 | -0.806 | -0.068 | 72.0 | 279.0 |
| Chronic constipation | Feces | -0.509 | -0.471 | 0.038 | 34.0 | 443.0 |
| Chronic obstructive pulmonary disease | Lower airway / lung | -0.486 | -0.200 | 0.286 | 25.0 | 56.0 |
| Periodontitis | Oral cavity / dental | -0.472 | -0.250 | 0.222 | 48.0 | 108.0 |
| Type II diabetes mellitus | Oral cavity / dental | -0.463 | -0.263 | 0.200 | 19.0 | 20.0 |

## Top Cross-Disease Taxon Switches

These taxa have at least one disease condition where they are predominantly reported increased and at least one where they are predominantly reported decreased.

| taxon_name | taxonomic_rank | oxygen_class | switch_score | increased_dominant_condition_count | decreased_dominant_condition_count | top_increased_conditions | top_decreased_conditions |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Bacteroides | genus | mixed_oxygen_annotation | 191.27 | 10 | 12 | Obesity (15 profiles, 7 studies); Chronic constipation (13 profiles, 7 studies); Gestational diabetes (7 profiles, 3 studies); Type II diabetes mellitus (6 profiles, 3 studies) | HIV infection (21 profiles, 10 studies); Diarrhea (8 profiles, 2 studies); Parkinson's disease (6 profiles, 7 studies); Gastric cancer (5 profiles, 4 studies) |
| Prevotella | genus | anaerobic_only | 182.45 | 14 | 11 | HIV infection (14 profiles, 14 studies); Colorectal cancer (12 profiles, 13 studies); Obesity (9 profiles, 9 studies); Rheumatoid arthritis (7 profiles, 3 studies) | Irritable bowel syndrome (10 profiles, 5 studies); COVID-19 (9 profiles, 12 studies); Crohn's disease (7 profiles, 3 studies); Age-related macular degeneration (6 profiles, 2 studies) |
| Lactobacillus | genus | mixed_oxygen_annotation | 145.28 | 12 | 9 | Parkinson's disease (16 profiles, 17 studies); Gastric cancer (8 profiles, 5 studies); Obesity (6 profiles, 3 studies); Colorectal cancer (5 profiles, 7 studies) | Cervical glandular intraepithelial neoplasia (10 profiles, 4 studies); Human papilloma virus infection (7 profiles, 7 studies); Chronic constipation (6 profiles, 3 studies); Human immunodeficiency virus (6 profiles, 2 studies) |
| Bifidobacterium | genus | mixed_oxygen_annotation | 105.61 | 5 | 9 | Parkinson's disease (23 profiles, 14 studies); Gestational diabetes (6 profiles, 3 studies); Crohn's disease (4 profiles, 3 studies); Gastric cancer (3 profiles, 3 studies) | Irritable bowel syndrome (13 profiles, 9 studies); Eczema (8 profiles, 2 studies); COVID-19 (7 profiles, 5 studies); Chronic constipation (6 profiles, 5 studies) |
| Ruminococcus | genus | anaerobic_only | 79.17 | 5 | 8 | Obesity (14 profiles, 5 studies); Gestational diabetes (5 profiles, 3 studies); Chronic constipation (4 profiles, 4 studies); Diabetes mellitus (4 profiles, 2 studies) | COVID-19 (10 profiles, 8 studies); Parkinson's disease (9 profiles, 8 studies); Irritable bowel syndrome (9 profiles, 4 studies); Colorectal cancer (5 profiles, 9 studies) |
| Streptococcus | genus | mixed_oxygen_annotation | 78.23 | 24 | 5 | Colorectal cancer (15 profiles, 14 studies); Gastric cancer (12 profiles, 5 studies); Parkinson's disease (11 profiles, 14 studies); Obesity (10 profiles, 5 studies) | COVID-19 (8 profiles, 10 studies); Chronic constipation (4 profiles, 5 studies); Eczema (4 profiles, 3 studies); Chlamydia trachomatis (4 profiles, 2 studies) |
| Blautia | genus | anaerobic_only | 67.70 | 5 | 9 | Gestational diabetes (6 profiles, 3 studies); HIV infection (5 profiles, 6 studies); Glaucoma (5 profiles, 3 studies); Type II diabetes mellitus (5 profiles, 2 studies) | Parkinson's disease (17 profiles, 15 studies); Colorectal cancer (10 profiles, 11 studies); COVID-19 (10 profiles, 7 studies); Diabetes mellitus (5 profiles, 2 studies) |
| Parabacteroides | genus | anaerobic_only | 67.09 | 5 | 6 | Chronic constipation (11 profiles, 6 studies); Parkinson's disease (9 profiles, 8 studies); Obesity (7 profiles, 2 studies); Atopic eczema (3 profiles, 2 studies) | Irritable bowel syndrome (8 profiles, 4 studies); HIV infection (5 profiles, 4 studies); COVID-19 (4 profiles, 6 studies); Gestational diabetes (4 profiles, 2 studies) |
| Clostridium | genus | mixed_oxygen_annotation | 66.69 | 8 | 4 | Cervical glandular intraepithelial neoplasia (5 profiles, 2 studies); Chronic constipation (3 profiles, 3 studies); Colorectal adenoma (3 profiles, 3 studies); Schizophrenia (3 profiles, 3 studies) | Colorectal cancer (11 profiles, 11 studies); Parkinson's disease (10 profiles, 10 studies); COVID-19 (8 profiles, 4 studies); HIV infection (7 profiles, 6 studies) |
| Veillonella | genus | anaerobic_only | 65.16 | 12 | 4 | COVID-19 (13 profiles, 11 studies); Parkinson's disease (5 profiles, 6 studies); Schizophrenia (5 profiles, 5 studies); Type II diabetes mellitus (5 profiles, 2 studies) | Colorectal cancer (8 profiles, 12 studies); Gestational diabetes (8 profiles, 3 studies); Food allergy (4 profiles, 4 studies); Human immunodeficiency virus (3 profiles, 2 studies) |

## Within-Condition Direction Conflicts

225 condition/taxon rows have at least 3 increased and 3 decreased profile occurrences. These are not necessarily contradictions; many combine different experiments, cohorts, or body sites under the same condition label.

## Files

- `global_oxygen_direction_summary.tsv`
- `condition_oxygen_direction_summary.tsv`
- `condition_oxygen_switch_scores.tsv`
- `taxon_condition_direction_summary.tsv`
- `taxa_cross_disease_switch_summary.tsv`
- `taxa_within_condition_direction_conflicts.tsv`
- `condition_oxygen_switch_scores.png`
