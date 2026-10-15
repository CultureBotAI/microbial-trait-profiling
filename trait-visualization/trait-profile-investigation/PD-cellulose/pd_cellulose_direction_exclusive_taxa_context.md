# PD Cellulose Taxa: Direction-Specific Reports and Background

Prepared from the local analysis inputs on 2026-09-11.

## Scope and Selection

This table expands the earlier nine-taxon shortlist. Taxa must have a cellulose catabolic process annotation in the local taxa-trait matrix, be reported in PD fecal profiles in exactly one direction, and occur in at least two distinct studies in that direction. Among 60 cellulose-annotated taxon labels observed in PD feces, 30 have reports in only one direction; nine meet the two-study threshold. The other 21 are supported by one study each.

Counts refer to exact reported taxon IDs. A profile is a study-experiment-direction combination; repeated mentions of the same taxon within a profile count once. These are curated differential-abundance signatures, not counts of patients or measurements of taxon prevalence in patients. The taxon table includes all PD fecal profiles, regardless of the profile's cellulose z-score.

## Background Denominators

| Dataset | Increased profiles | Increased studies | Decreased profiles | Decreased studies | Total profiles | Unique studies across directions |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PD feces | 130 | 50 | 123 | 46 | 253 | 51 |
| Other feces | 1,766 | 591 | 1,653 | 590 | 3,419 | 615 |
| All feces | 1,896 | 641 | 1,776 | 636 | 3,672 | 666 |
| All human profiles | 3,050 | 1,045 | 2,779 | 1,020 | 5,829 | 1,080 |

Other feces means non-PD-condition fecal profiles in the same human-host metadata universe. It is not a healthy-control cohort. All feces combines PD and other fecal profiles; all human profiles additionally includes other body sites. Direction-specific study counts overlap: 45 PD studies contribute both directions, so 50 + 46 - 45 = 51 unique studies. Do not add increased and decreased study counts.

## Shortlist Within PD Feces

Each profile cell is the number of profiles reporting the taxon / all PD fecal profiles in that direction (percentage). Study cells use the number of studies with at least one such report / all PD fecal studies in the supporting direction. Opposite-direction study counts are zero for every shortlisted taxon.

| Taxon | Rank | Reported PD direction | Increased profiles / 130 | Decreased profiles / 123 | Supporting studies / direction total |
| --- | --- | --- | ---: | ---: | ---: |
| Roseburia faecis | species | decreased-only | 0/130 (0.0%) | 6/123 (4.9%) | 5/46 (10.9%) |
| Serratia | genus | increased-only | 6/130 (4.6%) | 0/123 (0.0%) | 3/50 (6.0%) |
| Blautia obeum | species | increased-only | 5/130 (3.8%) | 0/123 (0.0%) | 2/50 (4.0%) |
| Actinomycetota | phylum | increased-only | 8/130 (6.2%) | 0/123 (0.0%) | 7/50 (14.0%) |
| Bacteroidota | phylum | decreased-only | 0/130 (0.0%) | 3/123 (2.4%) | 3/46 (6.5%) |
| Enterobacterales | order | increased-only | 3/130 (2.3%) | 0/123 (0.0%) | 2/50 (4.0%) |
| Burkholderiales | order | increased-only | 2/130 (1.5%) | 0/123 (0.0%) | 2/50 (4.0%) |
| Comamonadaceae | family | decreased-only | 0/130 (0.0%) | 2/123 (1.6%) | 2/46 (4.3%) |
| Micrococcaceae | family | decreased-only | 0/130 (0.0%) | 2/123 (1.6%) | 2/46 (4.3%) |

## Same Taxa in the Non-PD Fecal Background

These are exact-ID reporting counts in all non-PD fecal profiles, using the same deduplication as the PD table. Profile and study denominators differ and are shown explicitly.

| Taxon | Increased profiles / 1,766 | Decreased profiles / 1,653 | Increased studies / 591 | Decreased studies / 590 |
| --- | ---: | ---: | ---: | ---: |
| Roseburia faecis | 7/1,766 (0.4%) | 11/1,653 (0.7%) | 7/591 (1.2%) | 11/590 (1.9%) |
| Serratia | 2/1,766 (0.1%) | 6/1,653 (0.4%) | 2/591 (0.3%) | 6/590 (1.0%) |
| Blautia obeum | 18/1,766 (1.0%) | 39/1,653 (2.4%) | 14/591 (2.4%) | 20/590 (3.4%) |
| Actinomycetota | 84/1,766 (4.8%) | 82/1,653 (5.0%) | 61/591 (10.3%) | 54/590 (9.2%) |
| Bacteroidota | 92/1,766 (5.2%) | 92/1,653 (5.6%) | 63/591 (10.7%) | 61/590 (10.3%) |
| Enterobacterales | 78/1,766 (4.4%) | 17/1,653 (1.0%) | 54/591 (9.1%) | 14/590 (2.4%) |
| Burkholderiales | 26/1,766 (1.5%) | 22/1,653 (1.3%) | 20/591 (3.4%) | 17/590 (2.9%) |
| Comamonadaceae | 8/1,766 (0.5%) | 9/1,653 (0.5%) | 4/591 (0.7%) | 7/590 (1.2%) |
| Micrococcaceae | 24/1,766 (1.4%) | 7/1,653 (0.4%) | 20/591 (3.4%) | 6/590 (1.0%) |

## Relation to the Cellulose Profile Plot

The plot counts profiles whose normalized cellulose trait z-score is above zero. That is a different filter from taxon occurrence. There are 35/130 (26.9%) positive-z increased PD fecal profiles, supported by 22 studies, and 67/123 (54.5%) positive-z decreased PD fecal profiles, supported by 32 studies. Do not use 35 and 67 as denominators for the taxon-occurrence table: its numerators were counted across all 130 and 123 profiles.

## Interpretation

- Roseburia faecis has the strongest repeated-study, species-level decreased-only support in this shortlist: 6 profiles from 5 studies. Its decreased reporting rate is 4.9% in PD fecal profiles versus 0.7% in the non-PD fecal background.
- Serratia has increased-only reports in 6 PD profiles from 3 studies. The genus is reported in both directions outside PD (2 increased and 6 decreased profiles), so PD exclusivity should not be generalized across conditions.
- Blautia obeum has increased-only reports in 5 PD profiles, but these come from 2 studies. The broader Blautia genus has mixed-direction reports (2 increased and 15 decreased PD profiles); exact taxonomic resolution matters.
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
| Micrococcaceae | NCBITaxon:1268 |

## Sources and Verification

- [Original PD cellulose taxon counts](pd_cellulose_taxon_drivers.tsv). All nine PD profile and study counts were recomputed from the raw occurrence table and matched this source.
- [PD and other-feces summary](pd_cellulose_vs_other_feces_summary.tsv). The four PD/background direction-specific denominators were recomputed from metadata and matched this source.
- [Cellulose profile support](pd_cellulose_display_edge_summary.tsv). Source of the 35/67 positive-z profile counts and 22/32 supporting-study counts.
- [Profile metadata](../../study_trait_profile_study_metadata.tsv). Source of condition/body-site restrictions and deduplicated cohort study totals.
- [Taxon occurrence records](../../../data/derived/taxa.tsv). Joined to profile metadata by study, experiment, and direction; deduplicated by taxon ID and profile.
- [Taxa-trait counting implementation](../../../src/bsdbng/taxa_trait_matrix.py). See compute_trait_counts for annotation inheritance.
