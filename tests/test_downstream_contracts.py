"""Unit tests for shared downstream workflow contracts B1-B7 (#776)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mtp.datamodel.workflow_pydantic import (
    CalibrationOutcome,
    ContrastPopulation,
    DirectionBasis,
)
from mtp.integrated_evidence import read_status
from mtp.workflow.arms.__main__ import build_parser
from mtp.workflow.downstream_contracts import (
    STANDARD_ANALYSIS_SPECS,
    build_analysis_basis,
    build_calibration_verdict,
    build_eligibility_ledger,
    build_scientific_result_record,
    check_analysis_spec_compatibility,
    compute_common_support_ids,
    validate_qualified_feature_ids,
    validate_strict_numeric_value,
    validate_unique_record_ids,
)


def test_b1_analysis_spec_compatibility_detects_jed_estimator_mismatch() -> None:
    log_ratio_spec = STANDARD_ANALYSIS_SPECS["jed_table1_paired_log_ratio_sign_test_v1"]
    transformed_t_spec = STANDARD_ANALYSIS_SPECS[
        "jed_composition_arcsine_sqrt_ttest_calibration_v1"
    ]
    matched, mismatches_same = check_analysis_spec_compatibility(log_ratio_spec, log_ratio_spec)
    assert matched is True
    assert mismatches_same == []

    matched_diff, mismatches = check_analysis_spec_compatibility(log_ratio_spec, transformed_t_spec)
    assert matched_diff is False
    assert any("estimand:" in item for item in mismatches)
    assert any("kernel:" in item for item in mismatches)
    assert any("transformation:" in item for item in mismatches)


def test_b2_strict_readers_and_eligibility_ledger_and_common_support() -> None:
    assert validate_strict_numeric_value(
        "1.25", field_name="delta", row_id="row1"
    ) == pytest.approx(1.25)
    for bad in ("", "NA", "NaN", "inf", "abc", None):
        with pytest.raises(ValueError, match="delta"):
            validate_strict_numeric_value(bad, field_name="delta", row_id="row1")

    validate_unique_record_ids(["E1", "E2", "E3"], context="experiments")
    with pytest.raises(ValueError, match="Duplicate record identifiers"):
        validate_unique_record_ids(["E1", "E2", "E1"], context="experiments")

    ledger_a = build_eligibility_ledger(
        ledger_id="ledger:a",
        analysis_spec_id="noel_table1_six_trait_sign_test_v1",
        contrast_population=ContrastPopulation.all_paired,
        included_ids=["P3", "P1", "P2"],
        exclusion_counts={"zero_retained_total": 1},
    )
    assert ledger_a.included_count == 3
    assert ledger_a.excluded_count == 1
    assert ledger_a.included_ids == ["P1", "P2", "P3"]

    ledger_b = build_eligibility_ledger(
        ledger_id="ledger:b",
        analysis_spec_id="jed_table1_paired_log_ratio_sign_test_v1",
        contrast_population=ContrastPopulation.all_paired,
        included_ids=["P2", "P3"],
        exclusion_counts={"zero_component": 1, "zero_retained_total": 1},
    )
    assert compute_common_support_ids(ledger_a, ledger_b) == ["P2", "P3"]


def test_b3_analysis_basis_rejects_display_labels_and_duplicates() -> None:
    valid_ids = validate_qualified_feature_ids(
        ["METPO:2000103|GO:0030245", "METPO:1000602"],
        context="test_b3",
    )
    assert valid_ids == ["METPO:2000103|GO:0030245", "METPO:1000602"]
    with pytest.raises(ValueError, match="Unqualified feature identifier"):
        validate_qualified_feature_ids(["Cellulose degradation"], context="test_b3")

    basis = build_analysis_basis(
        basis_id="basis:test",
        feature_catalog_sha256="a" * 64,
        feature_filter_id="sha256:" + ("b" * 64),
        retained_feature_count=317,
        profile_unit="count_profiles",
        aggregation_policy="sum_taxon_trait_assertions_per_signature",
        rollup_policy="ncbi_taxonomy_closure_rollup_off",
        denominator_policy="sum_retained_non_mediadive_trait_assertions",
        background_inventory_id="human_host_retained_profiles",
        weighting_policy="equal_experiment_weight",
    )
    assert basis.retained_feature_count == 317


def test_b4_calibration_verdict_enforces_estimator_match() -> None:
    noel_spec = STANDARD_ANALYSIS_SPECS["noel_table1_six_trait_sign_test_v1"]
    matched = build_calibration_verdict(
        verdict_id="verdict:matched",
        production_spec=noel_spec,
        evaluated_spec=noel_spec,
        replicate_budget_sufficient=True,
        null_control_passed=True,
        power_passed=True,
    )
    assert matched.estimator_match is True
    assert matched.outcome == CalibrationOutcome.calibrated

    jed_lr = STANDARD_ANALYSIS_SPECS["jed_table1_paired_log_ratio_sign_test_v1"]
    jed_t = STANDARD_ANALYSIS_SPECS["jed_composition_arcsine_sqrt_ttest_calibration_v1"]
    mismatched = build_calibration_verdict(
        verdict_id="verdict:mismatched",
        production_spec=jed_lr,
        evaluated_spec=jed_t,
        replicate_budget_sufficient=True,
        null_control_passed=True,
        power_passed=True,
    )
    assert mismatched.estimator_match is False
    assert mismatched.outcome == CalibrationOutcome.different_estimator_uncalibrated
    assert "uncalibrated by this harness" in mismatched.summary_message


def test_b5_read_status_binds_sensitivity_and_coverage_input_hashes(tmp_path: Path) -> None:
    expected = {
        "profile_sha256": "1" * 64,
        "experiments_sha256": "2" * 64,
        "full_dump_sha256": "3" * 64,
    }
    sens_file = tmp_path / "sensitivity_status.json"
    sens_file.write_text(
        json.dumps(
            {
                "sensitivity_status": "assessed_clear #183",
                "inputs": {
                    "profile_sha256": "9" * 64,
                    "experiments_sha256": "2" * 64,
                },
            }
        ),
        encoding="utf-8",
    )
    status_bad = read_status(
        sens_file,
        "sensitivity_status",
        "#183",
        expected_input_sha256=expected,
    )
    assert "unreadable #183" in status_bad.default
    assert status_bad.source is not None
    assert status_bad.source["binding"] == "mismatch_or_missing_identity"

    cov_file = tmp_path / "coverage_status.json"
    cov_file.write_text(
        json.dumps(
            {
                "coverage_status": "assessed_clear #184",
                "inputs": {
                    "profile_sha256": "1" * 64,
                    "experiments_sha256": "2" * 64,
                },
            }
        ),
        encoding="utf-8",
    )
    status_ok = read_status(
        cov_file,
        "coverage_status",
        "#184",
        expected_input_sha256=expected,
    )
    assert status_ok.default == "assessed_clear #184"
    assert status_ok.source is not None
    assert status_ok.source["binding"] == "matched"


def test_b6_arms_cli_exposes_compare_and_verify_comparison() -> None:
    parser = build_parser()
    parsed = parser.parse_args(
        [
            "compare",
            "core/run.json",
            "noel/arm_run.json",
            "jed/arm_run.json",
            "ilyas/arm_run.json",
            "out/comparison",
        ]
    )
    assert parsed.command == "compare"
    assert parsed.core_run == Path("core/run.json")
    assert parsed.noel_run == Path("noel/arm_run.json")
    assert parsed.jed_run == Path("jed/arm_run.json")
    assert parsed.ilyas_run == Path("ilyas/arm_run.json")
    assert parsed.output == Path("out/comparison")


def test_b7_scientific_result_record_rejects_unqualified_feature() -> None:
    direction_val = next(iter(DirectionBasis))
    with pytest.raises(ValueError, match="Unqualified feature identifier"):
        build_scientific_result_record(
            record_id="res:bad",
            contrast_population=ContrastPopulation.all_paired,
            direction_basis=direction_val,
            analysis_spec_id="noel_table1_six_trait_sign_test_v1",
            eligibility_ledger_id="ledger:test",
            analysis_basis_id="basis:test",
            calibration_verdict_id=None,
            feature_id="Cellulose catabolism",
            feature_label="Cellulose catabolism",
            effect_estimate=-0.02,
            effect_scale="paired_fraction_difference",
            p_value=1e-6,
            q_value=1e-5,
            effective_q_value=None,
            n_eligible=500,
        )
    rec = build_scientific_result_record(
        record_id="res:ok",
        contrast_population=ContrastPopulation.all_paired,
        direction_basis=direction_val,
        analysis_spec_id="noel_table1_six_trait_sign_test_v1",
        eligibility_ledger_id="ledger:test",
        analysis_basis_id="basis:test",
        calibration_verdict_id="verdict:noel_table1_six_trait_sign_test_v1",
        feature_id="METPO:2000103|GO:0030245",
        feature_label="Cellulose catabolism",
        effect_estimate=-0.02,
        effect_scale="paired_fraction_difference",
        p_value=1e-6,
        q_value=1e-5,
        effective_q_value=2e-5,
        n_eligible=500,
    )
    assert rec.feature_id == "METPO:2000103|GO:0030245"
