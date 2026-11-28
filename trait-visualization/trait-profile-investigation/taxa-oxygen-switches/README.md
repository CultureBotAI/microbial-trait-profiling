# Taxa Oxygen Switch Investigation

## Scope

- Human-host disease-like profile universe: 4,102 profiles, 860 studies.
- Genus/species taxon-profile occurrences after metadata join: 27,057.
- Oxygen classes are based on direct KG-Microbe trait columns: aerobic side = aerobic, microaerophilic, obligately aerobic, facultatively aerobic; anaerobic side = anaerobic, obligately anaerobic, strictly anaerobic.
- Note: this taxon oxygen-switch classification (`OXYGEN_SWITCH_AEROBIC_TERMS`) is its own construct and is not the shared oxygen guild (`OXYGEN_GUILDS`); its aerobic side includes facultative (`METPO:1000608`), microaerophilic (`METPO:1000604`), and obligately aerobic (`METPO:1000606`) assertions by design.
- Disease-like filtering excludes Healthy/reference rows, broad cohort/design labels, measurements, exposures, diet/intervention response labels, and other non-disease metadata labels.

## Global Direction Balance by Oxygen Class

| oxygen_class | taxon_profile_count | taxon_count | study_count | increased_taxon_profile_count | decreased_taxon_profile_count | direction_balance |
| --- | --- | --- | --- | --- | --- | --- |
| anaerobic_only | 11500 | 905 | 747 | 5862 | 5638 | 0.019 |
| mixed_oxygen_annotation | 8937 | 443 | 783 | 4856 | 4081 | 0.087 |
| aerobic_only | 3487 | 1024 | 444 | 1852 | 1635 | 0.062 |
| unannotated | 3133 | 1445 | 435 | 1616 | 1517 | 0.032 |

## Strongest Condition-Level Oxygen Switches

Positive scores mean aerobic-only taxa skew increased while anaerobic-only taxa skew decreased.

| condition | body_site_display | aerobe_up_anaerobe_down_score | aerobic_only_direction_balance | anaerobic_only_direction_balance | aerobic_only_taxon_profile_count | anaerobic_only_taxon_profile_count |
| --- | --- | --- | --- | --- | --- | --- |
| Polycystic ovary syndrome | Feces | 0.993 | 0.886 | -0.108 | 35.0 | 195.0 |
| Irritable bowel syndrome | Feces | 0.981 | 0.810 | -0.171 | 21.0 | 333.0 |
| COVID-19 | Lower airway / lung | 0.946 | 0.946 | 0.000 | 37.0 | 26.0 |
| Crohn's disease | Feces | 0.933 | 0.600 | -0.333 | 30.0 | 117.0 |
| COVID-19 | Feces | 0.619 | 0.280 | -0.339 | 50.0 | 336.0 |
| Breast cancer | Feces | 0.598 | 0.538 | -0.060 | 13.0 | 117.0 |
| Colorectal cancer | Feces | 0.592 | 0.625 | 0.033 | 32.0 | 523.0 |
| HIV infection | Feces | 0.422 | 0.545 | 0.124 | 22.0 | 226.0 |

Negative scores mean the reciprocal pattern: anaerobic-only taxa skew increased while aerobic-only taxa skew decreased.

| condition | body_site_display | aerobe_up_anaerobe_down_score | aerobic_only_direction_balance | anaerobic_only_direction_balance | aerobic_only_taxon_profile_count | anaerobic_only_taxon_profile_count |
| --- | --- | --- | --- | --- | --- | --- |
| Colorectal cancer | Intestinal / GI mucosa | -1.033 | -0.700 | 0.333 | 20.0 | 12.0 |
| Colorectal cancer | Lower GI tract | -0.884 | -0.407 | 0.476 | 27.0 | 168.0 |
| Rheumatoid arthritis | Oral cavity / dental | -0.861 | -0.476 | 0.385 | 42.0 | 117.0 |
| Human immunodeficiency virus | Feces | -0.796 | -0.833 | -0.037 | 12.0 | 108.0 |
| Gestational diabetes | Feces | -0.759 | -0.778 | -0.019 | 72.0 | 267.0 |
| Periodontitis | Oral cavity / dental | -0.731 | -0.259 | 0.472 | 27.0 | 72.0 |
| Type II diabetes mellitus | Oral cavity / dental | -0.672 | -0.619 | 0.053 | 21.0 | 19.0 |
| Glaucoma | Ocular / eye | -0.667 | 0.000 | 0.667 | 60.0 | 108.0 |

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
| Veillonella | genus | mixed_oxygen_annotation | 59.58 | 10 | 4 | COVID-19 (12 profiles, 6 studies); HIV infection (5 profiles, 5 studies); Parkinson's disease (5 profiles, 4 studies); Schizophrenia (5 profiles, 4 studies) | Gestational diabetes (8 profiles, 3 studies); Colorectal cancer (7 profiles, 7 studies); Food allergy (4 profiles, 3 studies); Human immunodeficiency virus (3 profiles, 2 studies) |
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
