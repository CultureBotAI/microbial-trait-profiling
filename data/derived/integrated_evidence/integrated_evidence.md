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
| `experiments_eligible / experiments_excluded_zero_count` | experiments the method's own rule kept and removed for a zero count; the zero-total pair exclusion the tests share is in the provenance record, with the methods it applies to, and the figure-cell rows do not apply it |
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
| oxygen_guild_cell_support | 4160 |  | 0 | 0 | 4160 | 0 | 0 | 4160 |
| oxygen_guild_mcnemar | 1 | guild_mcnemar | 1 | 0 | 0 | 0 | 0 | 1 |

## Figure cells, by what the figure does with them

Cell rows report support and display, never a test. They are listed here as
counts rather than one row each; the whole inventory is in the TSV.

The grid is the oxygen guild figure's own: 52 condition rows (51 conditions and the pooled Other) by 20 body-site blocks (19 sites and the pooled Other), each block
4 guild-direction columns wide. It comes from the figure's builder at its
cut: the top 60 conditions and top 20 body sites, the rest pooled into Other,
cells backed by fewer than 3 studies masked, and rows and blocks with nothing drawn left out.
The figure draws every paired human profile, so unlike the tests in this table its
cells include the experiment pairs with a zero retained-annotation total.

| Display | Test state | Cells |
| --- | --- | --- |
| drawn | untested | 340 |
| drawn | zero_annotation_signal | 20 |
| masked_below_min_studies | insufficient_support | 296 |
| not_drawn | not_sampled | 3504 |

Native ORA rows summarize tested profiles (`tested_profiles_summarized`).
Their hit counts and best/median profile q-values do not define an aggregate
trait or cell test. The corresponding descriptive support rows remain untested.

## Records

| Method | Feature | Level | Effect | Units | p | q (family) | Decision | Test state |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| directional_consistency_fraction_sign_test | facultatively_aerobic | corpus | 0.5802 | share of non-tied studies shifting toward the increased side | 1.74e-06 | 1.04e-05 (table1_six_fractions) | significant | tested |
| directional_consistency_fraction_sign_test | cellulose_catabolic_process | corpus | 0.4325 | share of non-tied studies shifting toward the increased side | 0.000217 | 0.000435 (table1_six_fractions) | significant | tested |
| directional_consistency_fraction_sign_test | coccus_shaped | corpus | 0.5446 | share of non-tied studies shifting toward the increased side | 0.00827 | 0.00827 (table1_six_fractions) | significant | tested |
| directional_consistency_fraction_sign_test | sphere_shaped | corpus | 0.5558 | share of non-tied studies shifting toward the increased side | 0.00142 | 0.00213 (table1_six_fractions) | significant | tested |
| directional_consistency_fraction_sign_test | arsenate_reduction | corpus | 0.5627 | share of non-tied studies shifting toward the increased side | 0.00229 | 0.00274 (table1_six_fractions) | significant | tested |
| directional_consistency_fraction_sign_test | moderately_halophilic | corpus | 0.5742 | share of non-tied studies shifting toward the increased side | 2.19e-05 | 6.57e-05 (table1_six_fractions) | significant | tested |
| directional_consistency_log_ratio_sign_test | facultatively_aerobic/cellulose_catabolic_process | corpus | 0.5974 | share of non-tied eligible studies with a higher log2 ratio | 3.28e-05 | 0.00138 (table1_fifteen_ratios) | significant | tested |
| directional_consistency_log_ratio_sign_test | facultatively_aerobic/coccus_shaped | corpus | 0.5045 | share of non-tied eligible studies with a higher log2 ratio | 0.846 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | facultatively_aerobic/sphere_shaped | corpus | 0.5499 | share of non-tied eligible studies with a higher log2 ratio | 0.019 | 0.189 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | facultatively_aerobic/arsenate_reduction | corpus | 0.5339 | share of non-tied eligible studies with a higher log2 ratio | 0.313 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | facultatively_aerobic/moderately_halophilic | corpus | 0.5271 | share of non-tied eligible studies with a higher log2 ratio | 0.209 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | cellulose_catabolic_process/coccus_shaped | corpus | 0.46 | share of non-tied eligible studies with a higher log2 ratio | 0.0942 | 0.586 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | cellulose_catabolic_process/sphere_shaped | corpus | 0.409 | share of non-tied eligible studies with a higher log2 ratio | 0.000212 | 0.00351 (table1_fifteen_ratios) | significant | tested |
| directional_consistency_log_ratio_sign_test | cellulose_catabolic_process/arsenate_reduction | corpus | 0.3682 | share of non-tied eligible studies with a higher log2 ratio | 5.53e-05 | 0.00138 (table1_fifteen_ratios) | significant | tested |
| directional_consistency_log_ratio_sign_test | cellulose_catabolic_process/moderately_halophilic | corpus | 0.4394 | share of non-tied eligible studies with a higher log2 ratio | 0.0181 | 0.189 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | coccus_shaped/sphere_shaped | corpus | 0.5376 | share of non-tied eligible studies with a higher log2 ratio | 0.0721 | 0.513 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | coccus_shaped/arsenate_reduction | corpus | 0.4741 | share of non-tied eligible studies with a higher log2 ratio | 0.449 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | coccus_shaped/moderately_halophilic | corpus | 0.5496 | share of non-tied eligible studies with a higher log2 ratio | 0.023 | 0.191 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | sphere_shaped/arsenate_reduction | corpus | 0.4818 | share of non-tied eligible studies with a higher log2 ratio | 0.611 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | sphere_shaped/moderately_halophilic | corpus | 0.5051 | share of non-tied eligible studies with a higher log2 ratio | 0.857 | 1 (table1_fifteen_ratios) | not_significant | tested |
| directional_consistency_log_ratio_sign_test | arsenate_reduction/moderately_halophilic | corpus | 0.5 | share of non-tied eligible studies with a higher log2 ratio | 1 | 1 (table1_fifteen_ratios) | not_significant | tested |
| oxygen_guild_mcnemar | aerobe/anaerobe | corpus | 1.429 | odds ratio of discordant study signs | 3.61e-06 | NA (guild_mcnemar) | undecided | tested |
| composition_arm_native_transformed_contrast | facultatively_aerobic | study | 0.2048 | study mean of increased-minus-decreased transformed annotation contrast | 1.4e-08 | 4.37e-06 (native_composition_pooled_study) | significant | tested |
| composition_arm_native_transformed_contrast | cellulose_catabolic_process | study | -0.1537 | study mean of increased-minus-decreased transformed annotation contrast | 5.29e-05 | 0.00207 (native_composition_pooled_study) | significant | tested |
| composition_arm_native_transformed_contrast | coccus_shaped | study | 0.1629 | study mean of increased-minus-decreased transformed annotation contrast | 2.01e-05 | 0.000934 (native_composition_pooled_study) | significant | tested |
| composition_arm_native_transformed_contrast | sphere_shaped | study | 0.1279 | study mean of increased-minus-decreased transformed annotation contrast | 0.000434 | 0.0097 (native_composition_pooled_study) | significant | tested |
| composition_arm_native_transformed_contrast | arsenate_reduction | study | 0.1886 | study mean of increased-minus-decreased transformed annotation contrast | 2.42e-06 | 0.000379 (native_composition_pooled_study) | significant | tested |
| composition_arm_native_transformed_contrast | moderately_halophilic | study | 0.1575 | study mean of increased-minus-decreased transformed annotation contrast | 5.46e-06 | 0.000569 (native_composition_pooled_study) | significant | tested |
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

The primary native result uses study-level t-tests, BH over 313 support/variance-eligible traits. It retains 2257 pairs (932 studies), including 55 zero-total pairs; 2202 pairs have positive totals on both sides. Pair-level p/q values are secondary fields. This historical population differs from S1.

## Full-trait directional screen (same S1 population)

The screen's primary correction is BY; BH is also shown. Within-six q-values above are unchanged.

| Trait | Family | BH q | BY q |
| --- | --- | --- | --- |
| facultatively_aerobic | 317 | 0.000424958 | 0.00269326 |
| cellulose_catabolic_process | 317 | 0.00643282 | 0.0407693 |
| coccus_shaped | 317 | 0.0374701 | 0.237474 |
| sphere_shaped | 317 | 0.010718 | 0.0679272 |
| arsenate_reduction | 317 | 0.0151052 | 0.0957324 |
| moderately_halophilic | 317 | 0.000992065 | 0.00628741 |

## Native profile ORA support

Enrichment: global BH q ≤ 0.001 and fold ≥ 1.5. Depletion: separate lower-tail global BH q ≤ 0.001 and fold ≤ 0.5. Counts describe annotation representation in reported signatures, not physiological activity.

| Trait | Tests | Enriched inc | Enriched dec | Inc/dec | Best hit q | Median hit q | Median hit fold | Depleted inc | Depleted dec |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| facultatively_aerobic | 4514 | 190 | 152 | 1.25 | 3.23552e-22 | 5.45221e-05 | 4.53528 | 0 | 0 |
| cellulose_catabolic_process | 4514 | 178 | 183 | 0.972678 | 6.67711e-15 | 6.86411e-05 | 14.0188 | 0 | 0 |
| coccus_shaped | 4514 | 243 | 198 | 1.22727 | 2.64726e-23 | 3.53893e-05 | 5.30552 | 0 | 0 |
| sphere_shaped | 4514 | 235 | 179 | 1.31285 | 3.51633e-23 | 3.63897e-05 | 9.94043 | 0 | 0 |
| arsenate_reduction | 4514 | 130 | 74 | 1.75676 | 6.34264e-15 | 4.61512e-05 | 34.4379 | 0 | 0 |
| moderately_halophilic | 4514 | 222 | 178 | 1.24719 | 1.63289e-32 | 1.75365e-05 | 9.57934 | 0 | 0 |

Reading rules: a blank q is not a negative result; `source_not_supplied` means the
method's output was not given to this run; a decision holds only inside its own family;
and the six Table 1 traits were selected on this corpus, so no row confirms them.
A row with a p and no q was not corrected for multiplicity, so it carries no
threshold decision either.
`integrated_evidence_method_results.tsv` contains all non-cell rows for version control;
the complete generated TSV also includes the descriptive figure-cell inventory.
The family inventory lists primary row families. Auxiliary pair, full-screen and
depletion families retain their own sizes and sources in the relevant columns.
