# PD Cellulose Taxa: Reports by Direction and Background

Prepared from the local analysis inputs on 2026-09-11. Counts updated on 2026-10-18 after profiles took their own experiment's condition and body site instead of the study majority (#555), and on 2026-10-20 for the `count_input_16s_wms_v1` profile policy (#742).

## Scope and Selection

This table expands the earlier nine-taxon shortlist. Taxa must have a cellulose catabolic process annotation in the local taxa-trait matrix, be reported in PD fecal profiles in exactly one direction, and occur in at least two distinct studies in that direction. Among 58 cellulose-annotated taxon labels observed in PD feces, 29 have reports in only one direction; eight meet the two-study threshold. The other 21 are supported by one study each. Micrococcaceae met the threshold before that relabelling; it now has one decreased PD fecal profile from one study.

Counts refer to exact reported taxon IDs. A profile is a study-experiment-direction combination; repeated mentions of the same taxon within a profile count once. These are curated differential-abundance signatures, not counts of patients or measurements of taxon prevalence in patients. The taxon table includes all PD fecal profiles, regardless of the profile's cellulose z-score.

## Background Denominators

| Dataset | Increased profiles | Increased studies | Decreased profiles | Decreased studies | Total profiles | Unique studies across directions |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PD feces | 126 | 51 | 116 | 47 | 242 | 53 |
| Other feces | 1,677 | 594 | 1,574 | 584 | 3,251 | 617 |
| All feces | 1,803 | 640 | 1,690 | 629 | 3,493 | 664 |
| All human profiles | 2,948 | 1,016 | 2,681 | 987 | 5,629 | 1,047 |

Other feces means non-PD-condition fecal profiles in the same human-host metadata universe. It is not a healthy-control cohort. All feces combines PD and other fecal profiles; all human profiles additionally includes other body sites. Direction-specific study counts overlap: 45 PD studies contribute both directions, so 51 + 47 - 45 = 53 unique studies. Do not add increased and decreased study counts.

## Shortlist Within PD Feces

Each profile cell is the number of profiles reporting the taxon / all PD fecal profiles in that direction (percentage). Study cells use the number of studies with at least one such report / all PD fecal studies in the supporting direction. Opposite-direction study counts are zero for every shortlisted taxon.

| Taxon | Rank | Reported PD direction | Increased profiles / 126 | Decreased profiles / 116 | Supporting studies / direction total |
| --- | --- | --- | ---: | ---: | ---: |
| Roseburia faecis | species | decreased-only | 0/126 (0.0%) | 6/116 (5.2%) | 5/47 (10.6%) |
| Serratia | genus | increased-only | 6/126 (4.8%) | 0/116 (0.0%) | 3/51 (5.9%) |
| Blautia obeum | species | increased-only | 5/126 (4.0%) | 0/116 (0.0%) | 2/51 (3.9%) |
| Actinomycetota | phylum | increased-only | 8/126 (6.3%) | 0/116 (0.0%) | 7/51 (13.7%) |
| Bacteroidota | phylum | decreased-only | 0/126 (0.0%) | 3/116 (2.6%) | 3/47 (6.4%) |
| Enterobacterales | order | increased-only | 3/126 (2.4%) | 0/116 (0.0%) | 2/51 (3.9%) |
| Burkholderiales | order | increased-only | 2/126 (1.6%) | 0/116 (0.0%) | 2/51 (3.9%) |
| Comamonadaceae | family | decreased-only | 0/126 (0.0%) | 2/116 (1.7%) | 2/47 (4.3%) |

## Same Taxa in the Non-PD Fecal Background

These are exact-ID reporting counts in all non-PD fecal profiles, using the same deduplication as the PD table. Profile and study denominators differ and are shown explicitly.

| Taxon | Increased profiles / 1,677 | Decreased profiles / 1,574 | Increased studies / 594 | Decreased studies / 584 |
| --- | ---: | ---: | ---: | ---: |
| Roseburia faecis | 7/1,677 (0.4%) | 10/1,574 (0.6%) | 7/594 (1.2%) | 10/584 (1.7%) |
| Serratia | 2/1,677 (0.1%) | 5/1,574 (0.3%) | 2/594 (0.3%) | 5/584 (0.9%) |
| Blautia obeum | 15/1,677 (0.9%) | 38/1,574 (2.4%) | 13/594 (2.2%) | 19/584 (3.3%) |
| Actinomycetota | 82/1,677 (4.9%) | 82/1,574 (5.2%) | 59/594 (9.9%) | 54/584 (9.2%) |
| Bacteroidota | 88/1,677 (5.2%) | 90/1,574 (5.7%) | 60/594 (10.1%) | 60/584 (10.3%) |
| Enterobacterales | 79/1,677 (4.7%) | 17/1,574 (1.1%) | 55/594 (9.3%) | 14/584 (2.4%) |
| Burkholderiales | 24/1,677 (1.4%) | 23/1,574 (1.5%) | 18/594 (3.0%) | 18/584 (3.1%) |
| Comamonadaceae | 8/1,677 (0.5%) | 9/1,574 (0.6%) | 4/594 (0.7%) | 7/584 (1.2%) |

## Relation to the Cellulose Profile Plot

The plot counts profiles whose normalized cellulose trait z-score is above zero. That is a different filter from taxon occurrence. There are 34/126 (27.0%) positive-z increased PD fecal profiles, supported by 22 studies, and 69/116 (59.5%) positive-z decreased PD fecal profiles, supported by 34 studies. Do not use 34 and 69 as denominators for the taxon-occurrence table: its numerators were counted across all 126 and 116 profiles.

## Interpretation

- Roseburia faecis has the strongest repeated-study, species-level decreased-only support in this shortlist: 6 profiles from 5 studies. Its decreased reporting rate is 5.2% in PD fecal profiles versus 0.6% in the non-PD fecal background.
- Serratia has increased-only reports in 6 PD profiles from 3 studies. The genus is reported in both directions outside PD (2 increased and 5 decreased profiles), so PD exclusivity should not be generalized across conditions.
- Blautia obeum has increased-only reports in 5 PD profiles, but these come from 2 studies. The broader Blautia genus has mixed-direction reports (2 increased and 16 decreased PD profiles); exact taxonomic resolution matters.
- Actinomycetota has the largest supporting study count (7), but it is a phylum-level label. Higher-rank labels are less specific than species/genus candidates and can overlap taxonomically with other rows.
- Zero opposite-direction reports means no curated differential-abundance signature for that exact label in the opposite direction. It does not establish biological absence or consistent measurement across every study.
- These are descriptive candidates, not newly established significant biomarkers. Profiles within a study are not independent, studies vary in what they report, and the shortlist was selected for directional purity. Formal testing would need to account for study clustering and multiple testing across the full tested candidate set, rather than treating this selected table as an unbiased test set.
- Cellulose annotations in the local matrix can be inherited from trait-bearing descendant taxa. A broad taxon annotation does not establish the trait in every member, and counts must not be summed across taxonomic ranks.

## Taxon IDs

| Taxon | NCBI taxon ID |
| --- | --- |
| Roseburia faecis | NCBITaxon:301302 |
| Serratia | NCBITaxon:613 |
| Blautia obeum | NCBITaxon:40520 |
| Actinomycetota | NCBITaxon:201174 |
| Bacteroidota | NCBITaxon:976 |
| Enterobacterales | NCBITaxon:91347 |
| Burkholderiales | NCBITaxon:80840 |
| Comamonadaceae | NCBITaxon:80864 |

## Sources and Verification

- [Original PD cellulose taxon counts](pd_cellulose_taxon_drivers.tsv). All eight PD profile and study counts were recomputed from the raw occurrence table and matched this source.
- [PD and other-feces summary](pd_cellulose_vs_other_feces_summary.tsv). The four PD/background direction-specific denominators were recomputed from metadata and matched this source.
- [Cellulose profile support](pd_cellulose_display_edge_summary.tsv). Source of the 34/69 positive-z profile counts and 22/34 supporting-study counts.
- [Profile metadata](../../study_trait_profile_study_metadata.tsv). Source of condition/body-site restrictions and deduplicated cohort study totals.
- [Taxon occurrence records](../../../data/derived/taxa.tsv). Joined to profile metadata by study, experiment, and direction; deduplicated by taxon ID and profile.
- [Taxa-trait counting implementation](../../../src/bsdbng/taxa_trait_matrix.py). See compute_trait_counts for annotation inheritance.
