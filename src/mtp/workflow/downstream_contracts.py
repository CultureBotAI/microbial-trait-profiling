"""Executable downstream contracts for analysis specs, ledgers, bases, and verdicts (#776)."""

from __future__ import annotations

import json
import math
import re
from collections.abc import Iterable, Mapping, Sequence

from mtp.workflow.models import (
    AnalysisBasis,
    AnalysisSpec,
    CalibrationOutcome,
    CalibrationVerdict,
    ContrastPopulation,
    DirectionBasis,
    EligibilityLedger,
    ScientificResultRecord,
)

QUALIFIED_FEATURE_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*:[A-Za-z0-9_.:/|-]+$")

STANDARD_ANALYSIS_SPECS: dict[str, AnalysisSpec] = {
    "noel_trait_direction_screen_sign_test_v1": AnalysisSpec(
        id="noel_trait_direction_screen_sign_test_v1",
        estimand="paired_study_majority_direction_proportion",
        kernel="two_sided_exact_binomial_sign_test",
        transformation="study_majority_vote_of_contrast_fraction_difference",
        sampling_unit="study",
        weighting="equal_study_weight_after_within_study_majority_vote",
        eligibility_rule="paired_human_16S_WMS_nonzero_directional_difference_studies",
        multiple_testing_policy="BH_and_BY_and_eigenvalue_effective_M_across_retained_non_MediaDive_traits",
    ),
    "noel_table1_six_trait_sign_test_v1": AnalysisSpec(
        id="noel_table1_six_trait_sign_test_v1",
        estimand="paired_study_majority_direction_proportion",
        kernel="two_sided_exact_binomial_sign_test",
        transformation="study_majority_vote_of_contrast_fraction_difference",
        sampling_unit="study",
        weighting="equal_study_weight_after_within_study_majority_vote",
        eligibility_rule="paired_human_16S_WMS_nonzero_directional_difference_studies",
        multiple_testing_policy="BH_q_across_6_focal_traits_and_full_screen_reference",
    ),
    "jed_table1_paired_log_ratio_sign_test_v1": AnalysisSpec(
        id="jed_table1_paired_log_ratio_sign_test_v1",
        estimand="paired_study_log_ratio_direction_and_mean_shift",
        kernel="two_sided_exact_binomial_sign_test_and_paired_wilcoxon",
        transformation="within_signature_log_count_ratio_aggregated_to_study",
        sampling_unit="study",
        weighting="equal_study_weight",
        eligibility_rule="paired_human_16S_WMS_positive_numerator_and_denominator_in_both_directions",
        multiple_testing_policy="BH_and_BY_q_across_15_trait_pairs",
    ),
    "jed_composition_arcsine_sqrt_ttest_calibration_v1": AnalysisSpec(
        id="jed_composition_arcsine_sqrt_ttest_calibration_v1",
        estimand="paired_transformed_trait_fraction_shift",
        kernel="paired_student_t_test_and_wilcoxon_signed_rank",
        transformation="arcsine_sqrt_or_clr_or_alr_fraction",
        sampling_unit="study_or_contrast_pair",
        weighting="equal_pair_weight",
        eligibility_rule="paired_human_16S_WMS_retaining_zero_total_as_zero_fraction",
        multiple_testing_policy="BH_q_across_6_focal_traits",
    ),
    "jed_pooled_contrast_fraction_diff_v1": AnalysisSpec(
        id="jed_pooled_contrast_fraction_diff_v1",
        estimand="pooled_contrast_trait_fraction_difference",
        kernel="paired_wilcoxon_and_sign_test_and_paired_t",
        transformation="signature_trait_fraction_difference_inc_minus_dec",
        sampling_unit="paired_experiment_contrast",
        weighting="equal_contrast_weight",
        eligibility_rule="paired_human_16S_WMS_excluding_zero_total_pairs",
        multiple_testing_policy="BH_q_across_retained_non_MediaDive_traits",
    ),
    "ilyas_signature_hypergeometric_ora_v1": AnalysisSpec(
        id="ilyas_signature_hypergeometric_ora_v1",
        estimand="signature_trait_over_and_under_representation",
        kernel="one_sided_hypergeometric_test",
        transformation="binary_taxon_trait_presence_in_signature_vs_universe",
        sampling_unit="directional_signature",
        weighting="unweighted_unique_signature_taxa",
        eligibility_rule="human_16S_WMS_signatures_with_annotated_taxa",
        multiple_testing_policy="global_BH_fdr_across_signature_trait_tests",
    ),
    "historical_native_taxon_influence_experiment_weighted_v1": AnalysisSpec(
        id="historical_native_taxon_influence_experiment_weighted_v1",
        estimand="leave_one_taxon_out_influence_on_table1_signals",
        kernel="experiment_weighted_composition_and_assertion_count_ora",
        transformation="raw_experiment_fraction_difference_and_assertion_counts",
        sampling_unit="experiment_contrast_and_assertion",
        weighting="equal_experiment_weight",
        eligibility_rule="historical_experiment_paired_contrasts",
        multiple_testing_policy="descriptive_influence_delta_audit",
    ),
}


def get_standard_analysis_spec(spec_id: str) -> AnalysisSpec:
    """Return a copy of a declared standard AnalysisSpec or raise KeyError."""
    if spec_id not in STANDARD_ANALYSIS_SPECS:
        raise KeyError(f"Unknown standard AnalysisSpec id: {spec_id}")
    return STANDARD_ANALYSIS_SPECS[spec_id].model_copy(deep=True)


def check_analysis_spec_compatibility(
    production_spec: AnalysisSpec,
    evaluated_spec: AnalysisSpec,
) -> tuple[bool, list[str]]:
    """Compare two AnalysisSpec records across all methodological dimensions."""
    mismatches: list[str] = []
    for field_name in (
        "estimand",
        "kernel",
        "transformation",
        "sampling_unit",
        "weighting",
        "eligibility_rule",
        "multiple_testing_policy",
    ):
        prod_val = getattr(production_spec, field_name)
        eval_val = getattr(evaluated_spec, field_name)
        if prod_val != eval_val:
            mismatches.append(f"{field_name}: production={prod_val!r} vs evaluated={eval_val!r}")
    return (len(mismatches) == 0, mismatches)


def validate_strict_numeric_value(
    raw_value: object,
    *,
    field_name: str,
    row_id: str,
) -> float:
    """Parse a numeric cell strictly without coercing blank, NaN, or bad strings to zero."""
    if isinstance(raw_value, bool):
        raise ValueError(f"Boolean value is not valid numeric in {field_name!r} for {row_id!r}")
    if isinstance(raw_value, (int, float)):
        val = float(raw_value)
    elif isinstance(raw_value, str):
        stripped = raw_value.strip()
        if not stripped:
            raise ValueError(f"Blank numeric value in field {field_name!r} for record {row_id!r}")
        try:
            val = float(stripped)
        except ValueError as exc:
            raise ValueError(
                f"Malformed numeric value {raw_value!r} in field {field_name!r} "
                f"for record {row_id!r}"
            ) from exc
    else:
        raise ValueError(
            f"Unsupported numeric type {type(raw_value).__name__} in {field_name!r} "
            f"for record {row_id!r}"
        )
    if not math.isfinite(val):
        raise ValueError(
            f"Non-finite numeric value {raw_value!r} in field {field_name!r} for record {row_id!r}"
        )
    return val


def validate_unique_record_ids(record_ids: Sequence[str], *, context: str) -> list[str]:
    """Reject duplicate record identifiers instead of silently averaging duplicates."""
    seen: set[str] = set()
    duplicates: list[str] = []
    ordered: list[str] = []
    for raw_id in record_ids:
        rec_id = raw_id.strip()
        if not rec_id:
            raise ValueError(f"Empty record identifier in {context}")
        if rec_id in seen:
            duplicates.append(rec_id)
        else:
            seen.add(rec_id)
            ordered.append(rec_id)
    if duplicates:
        dup_preview = ", ".join(sorted(set(duplicates))[:5])
        raise ValueError(f"Duplicate record identifiers in {context}: {dup_preview}")
    return ordered


def validate_qualified_feature_ids(
    feature_ids: Iterable[str],
    *,
    context: str,
) -> list[str]:
    """Require qualified feature identifiers (prefix:id) rather than bare display labels."""
    validated: list[str] = []
    seen: set[str] = set()
    for raw_feature in feature_ids:
        feature_id = raw_feature.strip()
        if not QUALIFIED_FEATURE_PATTERN.match(feature_id):
            raise ValueError(
                f"Unqualified feature identifier {raw_feature!r} in {context}; "
                "expected qualified CURIE/pair identifier (e.g. GO:0030245, gold:facultative)"
            )
        if feature_id in seen:
            raise ValueError(f"Duplicate qualified feature identifier {feature_id!r} in {context}")
        seen.add(feature_id)
        validated.append(feature_id)
    if not validated:
        raise ValueError(f"No qualified feature identifiers supplied in {context}")
    return validated


def build_eligibility_ledger(
    *,
    ledger_id: str,
    analysis_spec_id: str,
    contrast_population: ContrastPopulation | str,
    included_ids: Iterable[str],
    exclusion_counts: Mapping[str, int],
    common_support_count: int | None = None,
) -> EligibilityLedger:
    """Construct a deterministic EligibilityLedger with reconciled exclusion counts."""
    sorted_included = sorted(
        validate_unique_record_ids(list(included_ids), context=f"EligibilityLedger({ledger_id})")
    )
    normalized_exclusions: dict[str, int] = {}
    total_excluded = 0
    for reason, count in sorted(exclusion_counts.items()):
        if count < 0:
            raise ValueError(
                f"Negative exclusion count for {reason!r} in "
                f"EligibilityLedger({ledger_id}): {count}"
            )
        normalized_exclusions[str(reason)] = int(count)
        total_excluded += int(count)
    return EligibilityLedger(
        id=ledger_id,
        analysis_spec_id=analysis_spec_id,
        contrast_population=ContrastPopulation(contrast_population),
        included_count=len(sorted_included),
        excluded_count=total_excluded,
        included_ids=sorted_included,
        exclusion_counts_json=json.dumps(
            normalized_exclusions, sort_keys=True, separators=(",", ":")
        ),
        common_support_count=common_support_count,
    )


def compute_common_support_ids(
    ledger_a: EligibilityLedger,
    ledger_b: EligibilityLedger,
) -> list[str]:
    """Return sorted intersection of included unit IDs across two EligibilityLedgers."""
    ids_a = set(ledger_a.included_ids or [])
    ids_b = set(ledger_b.included_ids or [])
    return sorted(ids_a & ids_b)


def build_analysis_basis(
    *,
    basis_id: str,
    feature_catalog_sha256: str,
    feature_filter_id: str,
    retained_feature_count: int,
    profile_unit: str,
    aggregation_policy: str,
    rollup_policy: str,
    denominator_policy: str,
    background_inventory_id: str,
    weighting_policy: str,
) -> AnalysisBasis:
    """Construct a schema-validated AnalysisBasis contract."""
    return AnalysisBasis(
        id=basis_id,
        feature_catalog_sha256=feature_catalog_sha256,
        feature_filter_id=feature_filter_id,
        retained_feature_count=retained_feature_count,
        profile_unit=profile_unit,
        aggregation_policy=aggregation_policy,
        rollup_policy=rollup_policy,
        denominator_policy=denominator_policy,
        background_inventory_id=background_inventory_id,
        weighting_policy=weighting_policy,
    )


def build_calibration_verdict(
    *,
    verdict_id: str,
    production_spec: AnalysisSpec,
    evaluated_spec: AnalysisSpec,
    replicate_budget_sufficient: bool,
    null_control_passed: bool | None,
    power_passed: bool | None,
    anti_conservative_detected: bool = False,
) -> CalibrationVerdict:
    """Evaluate calibration status with strict estimator-match enforcement."""
    estimator_match, mismatches = check_analysis_spec_compatibility(production_spec, evaluated_spec)
    if not estimator_match:
        mismatch_str = "; ".join(mismatches)
        return CalibrationVerdict(
            id=verdict_id,
            production_analysis_spec_id=production_spec.id,
            evaluated_analysis_spec_id=evaluated_spec.id,
            estimator_match=False,
            outcome=CalibrationOutcome.different_estimator_uncalibrated,
            null_control_passed=null_control_passed,
            power_passed=power_passed,
            replicate_budget_sufficient=replicate_budget_sufficient,
            summary_message=(
                f"Calibration evaluated {evaluated_spec.id!r} rather than production estimator "
                f"{production_spec.id!r} ({mismatch_str}); production estimator remains "
                "uncalibrated by this harness."
            ),
        )
    if not replicate_budget_sufficient:
        return CalibrationVerdict(
            id=verdict_id,
            production_analysis_spec_id=production_spec.id,
            evaluated_analysis_spec_id=evaluated_spec.id,
            estimator_match=True,
            outcome=CalibrationOutcome.insufficient_budget,
            null_control_passed=null_control_passed,
            power_passed=power_passed,
            replicate_budget_sufficient=False,
            summary_message=(
                f"Calibration harness for {production_spec.id!r} ran below the protocol "
                "minimum replicate budget."
            ),
        )
    if anti_conservative_detected or null_control_passed is False:
        return CalibrationVerdict(
            id=verdict_id,
            production_analysis_spec_id=production_spec.id,
            evaluated_analysis_spec_id=evaluated_spec.id,
            estimator_match=True,
            outcome=CalibrationOutcome.anti_conservative,
            null_control_passed=False,
            power_passed=power_passed,
            replicate_budget_sufficient=True,
            summary_message=(
                f"Empirical type-I error or null FDR for {production_spec.id!r} exceeded "
                "declared tolerance."
            ),
        )
    if null_control_passed is True and power_passed is True:
        return CalibrationVerdict(
            id=verdict_id,
            production_analysis_spec_id=production_spec.id,
            evaluated_analysis_spec_id=evaluated_spec.id,
            estimator_match=True,
            outcome=CalibrationOutcome.calibrated,
            null_control_passed=True,
            power_passed=True,
            replicate_budget_sufficient=True,
            summary_message=(
                f"Production estimator {production_spec.id!r} passed both null control "
                "and power criteria."
            ),
        )
    if null_control_passed is True and power_passed is False:
        return CalibrationVerdict(
            id=verdict_id,
            production_analysis_spec_id=production_spec.id,
            evaluated_analysis_spec_id=evaluated_spec.id,
            estimator_match=True,
            outcome=CalibrationOutcome.null_controlled_low_power,
            null_control_passed=True,
            power_passed=False,
            replicate_budget_sufficient=True,
            summary_message=(
                f"Production estimator {production_spec.id!r} controlled null error/FDR "
                "but fell below target power."
            ),
        )
    return CalibrationVerdict(
        id=verdict_id,
        production_analysis_spec_id=production_spec.id,
        evaluated_analysis_spec_id=evaluated_spec.id,
        estimator_match=True,
        outcome=CalibrationOutcome.partially_calibrated,
        null_control_passed=null_control_passed,
        power_passed=power_passed,
        replicate_budget_sufficient=True,
        summary_message=(
            f"Production estimator {production_spec.id!r} has partial calibration coverage "
            f"(null_control_passed={null_control_passed}, power_passed={power_passed})."
        ),
    )


def build_scientific_result_record(
    *,
    record_id: str,
    contrast_population: ContrastPopulation | str,
    direction_basis: DirectionBasis | str,
    analysis_spec_id: str,
    eligibility_ledger_id: str,
    analysis_basis_id: str,
    calibration_verdict_id: str | None,
    feature_id: str,
    feature_label: str,
    effect_estimate: float | None,
    effect_scale: str,
    p_value: float | None,
    q_value: float | None,
    effective_q_value: float | None,
    n_eligible: int,
    positive_count: int | None = None,
    negative_count: int | None = None,
    tie_or_zero_count: int | None = None,
) -> ScientificResultRecord:
    """Create a validated ScientificResultRecord with qualified feature ID check."""
    validate_qualified_feature_ids([feature_id], context=f"ScientificResultRecord({record_id})")
    return ScientificResultRecord(
        id=record_id,
        contrast_population=ContrastPopulation(contrast_population),
        direction_basis=DirectionBasis(direction_basis),
        analysis_spec_id=analysis_spec_id,
        eligibility_ledger_id=eligibility_ledger_id,
        analysis_basis_id=analysis_basis_id,
        calibration_verdict_id=calibration_verdict_id,
        feature_id=feature_id,
        feature_label=feature_label,
        effect_estimate=effect_estimate,
        effect_scale=effect_scale,
        p_value=p_value,
        q_value=q_value,
        effective_q_value=effective_q_value,
        n_eligible=n_eligible,
        positive_count=positive_count,
        negative_count=negative_count,
        tie_or_zero_count=tie_or_zero_count,
    )
