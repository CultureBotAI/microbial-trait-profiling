"""Distinguish caller-supplied, checkout-defaulted and absent assessments (#1766)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mtp import integrated_evidence
from mtp.integrated_provenance import verify_integrated_provenance
from mtp.paths import file_sha256
from test_integrated_portable_provenance import _arguments
from test_integrated_portable_provenance import portable_case as portable_case

_DEFAULTS = (
    (
        "--calibration-status",
        "DEFAULT_CALIBRATION_STATUS",
        "calibration_status",
        "calibration.json",
    ),
    (
        "--sensitivity-status",
        "DEFAULT_SENSITIVITY_STATUS",
        "sensitivity_status",
        "sensitivity.json",
    ),
    (
        "--orientation-provenance",
        "DEFAULT_ORIENTATION_STATUS",
        "orientation_status",
        "orientation.json",
    ),
)


@pytest.mark.parametrize("present", [False, True])
def test_omitted_assessments_retain_checkout_default_origin(
    portable_case: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch, present: bool
) -> None:
    root, _ = portable_case
    (root / "sensitivity.json").write_text(json.dumps({"scenario_status": {}}))
    output = root / f"defaulted-{present}"
    arguments = _arguments(root, output)
    candidates: dict[str, Path] = {}
    for option, constant, key, filename in _DEFAULTS:
        position = arguments.index(option)
        del arguments[position : position + 2]
        candidate = root / (filename if present else f"missing-default-{filename}")
        monkeypatch.setattr(integrated_evidence, constant, candidate)
        candidates[key] = candidate

    assert integrated_evidence.main(arguments) == 0
    record_path = output / "integrated_evidence_provenance.json"
    record = json.loads(record_path.read_text())
    for key, candidate in candidates.items():
        source = record["assessment_status"][key]["source"]
        assert source.get("supplied") is False
        assert source.get("input_origin") == "checkout_default"
        if present:
            assert source["path"] == candidate.name
            assert source["sha256"] == file_sha256(candidate)
        else:
            assert source["path"] is None
            assert source["binding"] == "not_supplied"
    result = verify_integrated_provenance(record_path)
    if present:
        verified = result["verified_references"]
        assert isinstance(verified, list)
        assert {f"assessment_status.{key}.source" for key in candidates} <= set(verified)
    else:
        skipped = result["skipped_references"]
        assert isinstance(skipped, list)
        assert {
            (f"assessment_status.{key}.source", "defaulted_but_missing") for key in candidates
        } <= {(item["location"], item["reason"]) for item in skipped}
    absent = record["assessment_status"]["coverage_status"]["source"]
    assert absent.get("input_origin") == "absent"
    assert absent.get("supplied") is False


def test_explicit_missing_assessments_remain_supplied(
    portable_case: tuple[Path, Path],
) -> None:
    root, _ = portable_case
    output = root / "explicit-missing"
    arguments = _arguments(root, output)
    for option, _, _, filename in _DEFAULTS:
        arguments[arguments.index(option) + 1] = str(root / f"explicit-missing-{filename}")
    assert integrated_evidence.main(arguments) == 0
    record_path = output / "integrated_evidence_provenance.json"
    record = json.loads(record_path.read_text())
    for _, _, key, _ in _DEFAULTS:
        source = record["assessment_status"][key]["source"]
        assert source == {
            "path": None,
            "binding": "missing",
            "supplied": True,
            "input_origin": "supplied",
        }
    result = verify_integrated_provenance(record_path)
    skipped = result["skipped_references"]
    assert isinstance(skipped, list)
    assert {
        (f"assessment_status.{key}.source", "supplied_but_missing") for _, _, key, _ in _DEFAULTS
    } <= {(item["location"], item["reason"]) for item in skipped}


def test_unreadable_checkout_default_retains_its_source_bytes(
    portable_case: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    root, _ = portable_case
    candidate = root / "unreadable-orientation.json"
    candidate.write_bytes(b"this is not JSON\n")
    expected_digest = file_sha256(candidate)
    monkeypatch.setattr(integrated_evidence, "DEFAULT_ORIENTATION_STATUS", candidate)
    output = root / "defaulted-unreadable"
    arguments = _arguments(root, output)
    position = arguments.index("--orientation-provenance")
    del arguments[position : position + 2]
    assert integrated_evidence.main(arguments) == 0
    record_path = output / "integrated_evidence_provenance.json"
    record = json.loads(record_path.read_text())
    assessment = record["assessment_status"]["orientation_status"]
    assert assessment["status_for_a_procedure_it_is_silent_about"].startswith("unreadable")
    source = assessment["source"]
    assert source["input_origin"] == "checkout_default"
    assert source["supplied"] is False
    assert source.get("path") == candidate.name
    assert source.get("sha256") == expected_digest
    result = verify_integrated_provenance(record_path)
    verified = result["verified_references"]
    assert isinstance(verified, list)
    assert "assessment_status.orientation_status.source" in verified
    candidate.write_bytes(b"changed source bytes\n")
    with pytest.raises(ValueError, match=r"hash|SHA"):
        verify_integrated_provenance(record_path)
