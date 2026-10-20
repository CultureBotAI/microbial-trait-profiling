# Integrated evidence table (generated)

Version integrated-evidence-v2. Every method keeps its own row, family and units. Nothing here
combines them: there is no pooled p-value, no averaged effect and no false discovery
rate for a set defined as significant under any method.

## Data dictionary

| Field | Meaning |
| --- | --- |
| `method` | which procedure produced the row; never merged with another |
| `method_version` | the procedure's own version string |
| `inference_level` | what the result is about: corpus, body site, condition and site, profile |
| `feature_type` | trait, trait_pair or guild_pair |
| `trait / reference` | the feature, and the reference where a comparison names one |
| `trait_qualifier / reference_qualifier` | the term ids the analysis declared for each side of the hypothesis |
| `trait_columns / reference_columns` | the raw matrix columns those terms resolved to; this is the join key, and it differs from the declared term for every METPO-qualified composite |
| `population / eligibility_rule` | who is included, and the rule that decides |
| `denominator / background` | what the effect is relative to; empty when none applies |
| `unique_studies / paired_experiments / profiles / testable_units / ties` | distinct counts, never interchangeable |
| `studies_total / studies_excluded / studies_partially_observed` | the conditional population: studies the method could have used, studies its eligibility rule removed, and studies contributing only some experiments. Blank where the method defines no such rule |
| `experiments_eligible / experiments_excluded_zero_count` | experiments the method's own rule kept and removed for a zero count; the corpus-wide zero-total pair exclusion is in the provenance record, not here |
| `effect_value / effect_units` | the estimate and what it measures; units differ by method |
| `study_effect_mean / study_effect_median / study_effect_units` | the distribution of per-study values the test ran over, which describes the input rather than the tested estimate, so it survives on an untested row |
| `interval_low / interval_high / interval_type` | the interval and what it is an interval for; 'unavailable' where none was computed |
| `p_value / q_value / correction_method / family_id / family_size` | the test and the family it was corrected in; families are not comparable |
| `q_threshold / effect_threshold / support_threshold` | the thresholds this row used: two for the decision, and the minimum studies a figure cell needs to be drawn. The support threshold is a display rule and is set only on a cell row, which is what identifies the cell inventory |
| `stratum_kind / stratum / condition / body_site / direction` | what the row is about within the corpus, as one label and as separate identifiers; blank on a pooled row |
| `test_status` | tested, untested, insufficient_support, insufficient_variance, all_tied, no_eligible_pairs, zero_annotation_signal, not_sampled, source_not_supplied or tested_profiles_summarized. These distinguish absent signal, absent sampling, absent source and a descriptive summary of profile tests. Insufficient variance means the native t-test was unavailable for a constant nonzero or nearly constant contrast; it is not a tested negative. Native all_tied requires exactly zero study contrasts |
| `decision` | significant, not_significant or undecided; only ever within one family |
| `display_status` | drawn, masked_below_min_studies, not_drawn or not_applicable; never a result. The first two belong to a figure cell, which is what the cell rows carry; a corpus-level test is drawn nowhere, so it is not_applicable |
| `zero_annotation_state` | how zero or absent annotation was treated |
| `selection_note` | whether the feature was selected on this corpus |
| `calibration_status / sensitivity_status / coverage_status / orientation_status / taxon_influence_status` | the state of #182, #183, #184, #188 and #190 for this row's own procedure. 'pending' means no applicable assessment was supplied or mapped; it does not establish that execution has not occurred. 'not_calibrated_here' means it ran and excluded this procedure; native influence completion reports execution, never a scientific pass. A report with no procedure verdict says so. Per-scenario and per-contrast-class detail is in integrated_evidence_assessments.tsv |
| `source` | the input the row was computed from |
| `family_scope / family_source` | complete internal family, source_subset, or profile_test_summary; external families retain their full source checksum |
| `statistical_test / native_pair_p_value / native_pair_q_value` | actual native test; composition pair-level values are secondary to study-level inference and use a separate BH family of the same size |
| `n_nonzero_experiments / n_zero_total_experiments / n_zero_total_excluded` | native composition eligibility counts; historical inclusion is explicit, not silently equated to S1 exclusion |
| `full_screen_bh_q / full_screen_by_q / full_screen_family_size / full_screen_source` | compatible full-trait directional screen values with matching input hashes, study counts and contrast; BY is its primary correction |
| `n_tests / n_significant_inc / n_significant_dec / inc_dec_ratio` | profile-level enrichment tests and discoveries by signature direction; descriptive, not aggregate significance; ratio blank when decreased-hit count is zero |
| `best_profile_q / median_profile_q / median_hit_fold / max_hit_fold` | summaries conditional on enrichment hits at the native q and fold cutoffs; not inferential trait-level q-values |
| `n_depleted_inc / n_depleted_dec / depletion_q_threshold / depletion_max_fold / depletion_family_size / depletion_source` | lower-tail discoveries in a separate global family; missing source leaves unknown counts, not zero |

## What each method contributed

| Method | Records | Families | Tested | Profile summaries | Not tested | Significant | Not significant | Undecided |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| aggregated_profile_support | 6 |  | 0 | 0 | 6 | 0 | 0 | 6 |
| composition_arm_native_transformed_contrast | 6 | native_composition_pooled_study | 6 | 0 | 0 | 6 | 0 | 0 |
| directional_consistency_fraction_sign_test | 6 | table1_six_fractions | 6 | 0 | 0 | 6 | 0 | 0 |
| directional_consistency_log_ratio_sign_test | 15 | table1_fifteen_ratios | 15 | 0 | 0 | 3 | 12 | 0 |
| enrichment_arm_profile_hypergeometric | 6 | native_profile_ORA_enrichment_global | 0 | 6 | 0 | 0 | 0 | 6 |
| oxygen_guild_cell_support | 1872 |  | 0 | 0 | 1872 | 0 | 0 | 1872 |
| oxygen_guild_mcnemar | 1 | guild_mcnemar | 1 | 0 | 0 | 0 | 0 | 1 |

## Figure cells, by what the figure does with them

Cell rows report support and display, never a test. They are listed here as
counts rather than one row each; the whole inventory is in the TSV.

| Display | Test state | Cells |
| --- | --- | --- |
| drawn | untested | 316 |
| drawn | zero_annotation_signal | 7 |
| masked_below_min_studies | insufficient_support | 209 |
| not_drawn | not_sampled | 1340 |

Native ORA rows summarize tested profiles (`tested_profiles_summarized`).
Their hit counts and best/median profile q-values do not define an aggregate
trait or cell test. The corresponding descriptive support rows remain untested.

## Records

| Method | Feature | Level | Effect | Units | p | q (family) | Decision | Test state |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| directional_consistency_fraction_sign_test | facultatively_aerobic | corpus | 0.5814 | share of non-tied studies shifting toward the increased side | 9.44e-07 | 5.66e-06 (table1_six_fractions) | significant | tested |
| directional_consistency_fraction_sign_test | cellulose_catabolic_process | corpus | 0.4302 | share of non-tied studies shifting toward the increased side | 0.000117 | 0.000235 (table1_six_fractions) | significant | tested |
| directional_consistency_fraction_sign_test | coccus_shaped | corpus | 0.545 | share of non-tied studies shifting toward the increased side | 0.00728 | 0.00728 (table1_six_fractions) | significant | tested |
| directional_consistency_fraction_sign_test | sphere_shaped | corpus | 0.554 | share of non-tied studies shifting toward the increased side | 0.00192 | 0.00288 (table1_six_fractions) | significant | tested |
| directional_consistency_fraction_sign_test | arsenate_reduction | corpus | 0.5603 | share of non-tied studies shifting toward the increased side | 0.00319 | 0.00382 (table1_six_fractions) | significant | tested |
| directional_consistency_fraction_sign_test | moderately_halophilic | corpus | 0.5718 | share of non-tied studies shifting toward the increased side | 3.48e-05 | 0.000104 (table1_six_fractions) | significant | tested |
| directional_consistency_log_ratio_sign_test | facultatively_aerobic/cellulose_catabolic_process | corpus | 0.5953 | share of non-tied eligible studies with a higher log2 ratio | 4.43e-05 | 0.0011 (table1_fifteen_ratios) | significant | tested |
| directional_consistency_log_ratio_sign_test | facultatively_aerobic/coccus_shaped | corpus | 0.5081 | share of non-tied eligible studies with a higher log2 ratio | 0.7 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | facultatively_aerobic/sphere_shaped | corpus | 0.5529 | share of non-tied eligible studies with a higher log2 ratio | 0.0124 | 0.124 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | facultatively_aerobic/arsenate_reduction | corpus | 0.5412 | share of non-tied eligible studies with a higher log2 ratio | 0.21 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | facultatively_aerobic/moderately_halophilic | corpus | 0.5318 | share of non-tied eligible studies with a higher log2 ratio | 0.135 | 0.748 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | cellulose_catabolic_process/coccus_shaped | corpus | 0.4615 | share of non-tied eligible studies with a higher log2 ratio | 0.106 | 0.657 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | cellulose_catabolic_process/sphere_shaped | corpus | 0.4075 | share of non-tied eligible studies with a higher log2 ratio | 0.000154 | 0.00255 (table1_fifteen_ratios) | significant | tested |
| directional_consistency_log_ratio_sign_test | cellulose_catabolic_process/arsenate_reduction | corpus | 0.3651 | share of non-tied eligible studies with a higher log2 ratio | 3.38e-05 | 0.0011 (table1_fifteen_ratios) | significant | tested |
| directional_consistency_log_ratio_sign_test | cellulose_catabolic_process/moderately_halophilic | corpus | 0.435 | share of non-tied eligible studies with a higher log2 ratio | 0.0107 | 0.124 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | coccus_shaped/sphere_shaped | corpus | 0.5355 | share of non-tied eligible studies with a higher log2 ratio | 0.0876 | 0.623 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | coccus_shaped/arsenate_reduction | corpus | 0.4784 | share of non-tied eligible studies with a higher log2 ratio | 0.531 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | coccus_shaped/moderately_halophilic | corpus | 0.5479 | share of non-tied eligible studies with a higher log2 ratio | 0.0269 | 0.223 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | sphere_shaped/arsenate_reduction | corpus | 0.4821 | share of non-tied eligible studies with a higher log2 ratio | 0.614 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | sphere_shaped/moderately_halophilic | corpus | 0.501 | share of non-tied eligible studies with a higher log2 ratio | 1 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | arsenate_reduction/moderately_halophilic | corpus | 0.5 | share of non-tied eligible studies with a higher log2 ratio | 1 | 1 (table1_fifteen_ratios) | not_significant | tested |
| oxygen_guild_mcnemar | aerobe/anaerobe | corpus | 1.433 | odds ratio of discordant study signs | 2.24e-06 | NA (guild_mcnemar) | undecided | tested |
| composition_arm_native_transformed_contrast | facultatively_aerobic | study | 0.2027 | study mean of increased-minus-decreased transformed annotation contrast | 6.81e-09 | 2.15e-06 (native_composition_pooled_study) | significant | tested |
| composition_arm_native_transformed_contrast | cellulose_catabolic_process | study | -0.1558 | study mean of increased-minus-decreased transformed annotation contrast | 3.19e-05 | 0.00126 (native_composition_pooled_study) | significant | tested |
| composition_arm_native_transformed_contrast | coccus_shaped | study | 0.1652 | study mean of increased-minus-decreased transformed annotation contrast | 8.97e-06 | 0.000706 (native_composition_pooled_study) | significant | tested |
| composition_arm_native_transformed_contrast | sphere_shaped | study | 0.1237 | study mean of increased-minus-decreased transformed annotation contrast | 0.000517 | 0.0116 (native_composition_pooled_study) | significant | tested |
| composition_arm_native_transformed_contrast | arsenate_reduction | study | 0.1861 | study mean of increased-minus-decreased transformed annotation contrast | 2.23e-06 | 0.000351 (native_composition_pooled_study) | significant | tested |
| composition_arm_native_transformed_contrast | moderately_halophilic | study | 0.1541 | study mean of increased-minus-decreased transformed annotation contrast | 5.99e-06 | 0.000629 (native_composition_pooled_study) | significant | tested |
| enrichment_arm_profile_hypergeometric | facultatively_aerobic | profile | NA | NA | NA | NA (native_profile_ORA_enrichment_global) | undecided | tested_profiles_summarized |
| aggregated_profile_support | facultatively_aerobic | corpus | NA | NA | NA | NA (NA) | undecided | untested |
| enrichment_arm_profile_hypergeometric | cellulose_catabolic_process | profile | NA | NA | NA | NA (native_profile_ORA_enrichment_global) | undecided | tested_profiles_summarized |
| aggregated_profile_support | cellulose_catabolic_process | corpus | NA | NA | NA | NA (NA) | undecided | untested |
| enrichment_arm_profile_hypergeometric | coccus_shaped | profile | NA | NA | NA | NA (native_profile_ORA_enrichment_global) | undecided | tested_profiles_summarized |
| aggregated_profile_support | coccus_shaped | corpus | NA | NA | NA | NA (NA) | undecided | untested |
| enrichment_arm_profile_hypergeometric | sphere_shaped | profile | NA | NA | NA | NA (native_profile_ORA_enrichment_global) | undecided | tested_profiles_summarized |
| aggregated_profile_support | sphere_shaped | corpus | NA | NA | NA | NA (NA) | undecided | untested |
| enrichment_arm_profile_hypergeometric | arsenate_reduction | profile | NA | NA | NA | NA (native_profile_ORA_enrichment_global) | undecided | tested_profiles_summarized |
| aggregated_profile_support | arsenate_reduction | corpus | NA | NA | NA | NA (NA) | undecided | untested |
| enrichment_arm_profile_hypergeometric | moderately_halophilic | profile | NA | NA | NA | NA (native_profile_ORA_enrichment_global) | undecided | tested_profiles_summarized |
| aggregated_profile_support | moderately_halophilic | corpus | NA | NA | NA | NA (NA) | undecided | untested |

## Native composition population

The primary native result uses study-level t-tests, BH over 315 support/variance-eligible traits. It retains 2335 pairs (959 studies), including 94 zero-total pairs; 2241 pairs have positive totals on both sides. Pair-level p/q values are secondary fields. This historical population differs from S1.

## Full-trait directional screen (same S1 population)

The screen's primary correction is BY; BH is also shown. Within-six q-values above are unchanged.

| Trait | Family | BH q | BY q |
| --- | --- | --- | --- |
| facultatively_aerobic | 317 | 0.000189459 | 0.00120073 |
| cellulose_catabolic_process | 317 | 0.0046464 | 0.0294475 |
| coccus_shaped | 317 | 0.0366373 | 0.232196 |
| sphere_shaped | 317 | 0.0141507 | 0.0896829 |
| arsenate_reduction | 317 | 0.0213112 | 0.135064 |
| moderately_halophilic | 317 | 0.00157528 | 0.00998364 |

## Native profile ORA support

Enrichment: global BH q ≤ 0.001 and fold ≥ 1.5. Depletion: separate lower-tail global BH q ≤ 0.001 and fold ≤ 0.5. Counts describe annotation representation in reported signatures, not physiological activity.

| Trait | Tests | Enriched inc | Enriched dec | Inc/dec | Best hit q | Median hit q | Median hit fold | Depleted inc | Depleted dec |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| facultatively_aerobic | 4670 | 187 | 150 | 1.24667 | 3.34734e-22 | 4.73153e-05 | 4.54538 | 0 | 0 |
| cellulose_catabolic_process | 4670 | 181 | 186 | 0.973118 | 6.89036e-15 | 7.03782e-05 | 14.0188 | 0 | 0 |
| coccus_shaped | 4670 | 243 | 200 | 1.215 | 2.73875e-23 | 3.62787e-05 | 5.29681 | 0 | 0 |
| sphere_shaped | 4670 | 235 | 180 | 1.30556 | 3.63785e-23 | 3.71063e-05 | 9.90965 | 0 | 0 |
| arsenate_reduction | 4670 | 130 | 74 | 1.75676 | 6.54511e-15 | 4.73153e-05 | 34.4379 | 0 | 0 |
| moderately_halophilic | 4670 | 222 | 179 | 1.24022 | 1.68932e-32 | 1.79747e-05 | 9.51673 | 0 | 0 |

Reading rules: a blank q is not a negative result; `source_not_supplied` means the
method's output was not given to this run; a decision holds only inside its own family;
and the six Table 1 traits were selected on this corpus, so no row confirms them.
A row with a p and no q was not corrected for multiplicity, so it carries no
threshold decision either.
`integrated_evidence_method_results.tsv` contains all non-cell rows for version control;
the complete generated TSV also includes the descriptive figure-cell inventory.
The family inventory lists primary row families. Auxiliary pair, full-screen and
depletion families retain their own sizes and sources in the relevant columns.
