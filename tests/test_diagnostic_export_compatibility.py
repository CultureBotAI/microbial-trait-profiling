"""Historical diagnostic exports remain readable without the new schema label (#1487)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from mtp.integrated_evidence import assessment_rows, read_status

FIXTURES = Path(__file__).parent / "fixtures"


def test_legacy_calibration_export_remains_a_not_certified_assessment() -> None:
    path = FIXTURES / "calibration_status_quick_run.json"
    before = path.read_bytes()
    payload = json.loads(before)
    assert "export_schema_version" not in payload
    assert payload["version"] == "calibration-v2"

    status = read_status(path, "calibration_status", "#182")
    assert status.source["producer_version"] == "calibration-v2"
    assert status.source["sha256"] == hashlib.sha256(before).hexdigest()
    for procedure in ("fraction_sign_test", "log_ratio_sign_test", "guild_mcnemar"):
        assert status.for_procedure(procedure).startswith("not_certified #182 [")
        assert "insufficient_budget" in status.for_procedure(procedure)
    exported = assessment_rows([], {"calibration_status": status})
    assert any(
        row["scope"] == "complete_null_study_flip"
        and row["state_variable"] == "replicates_completed"
        and row["state_value"] == "5"
        for row in exported
    )
    assert path.read_bytes() == before


def test_legacy_sensitivity_export_retains_separate_family_evidence() -> None:
    path = FIXTURES / "sensitivity_status_c9fb6f1.json"
    before = path.read_bytes()
    payload = json.loads(before)
    assert "export_schema_version" not in payload
    assert payload["version"] == "log-ratio-sensitivity-v1"

    status = read_status(path, "sensitivity_status", "#183")
    assert status.source["sha256"] == hashlib.sha256(before).hexdigest()
    assert status.for_procedure("log_ratio_sign_test") == "scenarios reported, no verdict #183"
    assert status.for_procedure("fraction_sign_test") == "pending #183"
    exported = {
        (row["scope"], row["state_variable"], row["state_value"])
        for row in assessment_rows([], {"sensitivity_status": status})
    }
    assert {
        ("S4:six_fractions", "correction", "Benjamini-Hochberg"),
        ("S4:six_fractions", "family_size", "6"),
        ("S4:six_fractions", "q_at_or_below_0.05", "2"),
        ("S4:fifteen_ratios", "correction", "Benjamini-Yekutieli"),
        ("S4:fifteen_ratios", "family_size", "15"),
        ("S4:fifteen_ratios", "q_at_or_below_0.05", "4"),
    } <= exported
    assert path.read_bytes() == before
