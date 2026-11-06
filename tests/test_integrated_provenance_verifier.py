"""Whole-record integrated provenance verification, including relocation (#1755)."""

from __future__ import annotations

import importlib
import json
import shutil
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from mtp.paths import file_sha256
from test_integrated_portable_provenance import portable_case as portable_case


@pytest.fixture
def api() -> ModuleType:
    assert importlib.util.find_spec("mtp.integrated_provenance") is not None
    module = importlib.import_module("mtp.integrated_provenance")
    assert module.__file__ is not None
    assert Path(module.__file__).resolve() == (
        Path(__file__).resolve().parents[1] / "src/mtp/integrated_provenance.py"
    )
    return module


@pytest.fixture
def complete_record(tmp_path: Path) -> tuple[Path, dict[str, Any]]:
    root = tmp_path / "original"
    root.mkdir()
    references: dict[str, dict[str, str]] = {}
    for name in ("native", "ora", "assessment", "review", "protocol", "nested"):
        source = root / f"{name}.tsv"
        source.write_text(f"{name} bytes\n", encoding="utf-8")
        references[name] = {"path": source.name, "sha256": file_sha256(source)}
    record = {
        "version": "integrated-evidence-v2",
        "path_basis": {"kind": "root-relative-v1", "root": ".."},
        "native_sources": {
            "arm1": references["native"],
            "absent": {"path": None, "sha256": None, "binding": "not_supplied"},
            "missing": {"path": None, "sha256": None, "binding": "missing", "supplied": True},
        },
        "ora_eligibility_source": references["ora"],
        "assessment_status": {
            "calibration": {
                "source": {
                    **references["assessment"],
                    "protocol": references["protocol"],
                    "review": references["review"],
                    "metadata": json.dumps({"nested_source": references["nested"]}),
                }
            },
            "unsupplied": {"source": {"path": None, "binding": "not_supplied"}},
        },
    }
    output = root / "record/provenance.json"
    output.parent.mkdir()
    output.write_text(json.dumps(record), encoding="utf-8")
    return output, record


def test_whole_record_relocates_and_verifies_every_reference(
    complete_record: tuple[Path, dict[str, Any]], tmp_path: Path, api: ModuleType
) -> None:
    path, _ = complete_record
    root = path.parent.parent
    relocated = tmp_path / "relocated"
    shutil.copytree(root, relocated)
    root.rename(tmp_path / "hidden-original")
    record = relocated / "record/provenance.json"
    before = {p.relative_to(relocated): p.read_bytes() for p in relocated.rglob("*") if p.is_file()}
    result = api.verify_integrated_provenance(record)
    assert result["record_sha256"] == file_sha256(record)
    assert result["verified_reference_count"] == 6
    assert set(result["verified_references"]) == {
        "native_sources.arm1",
        "ora_eligibility_source",
        "assessment_status.calibration.source",
        "assessment_status.calibration.source.review",
        "assessment_status.calibration.source.protocol",
        "assessment_status.calibration.source.metadata.nested_source",
    }
    assert result["skipped_reference_count"] == 3
    assert {item["reason"] for item in result["skipped_references"]} == {
        "not_supplied",
        "supplied_but_missing",
    }
    after = {p.relative_to(relocated): p.read_bytes() for p in relocated.rglob("*") if p.is_file()}
    assert after == before


def test_actual_producer_record_passes_whole_verification(
    portable_case: tuple[Path, Path], api: ModuleType
) -> None:
    _, output = portable_case
    result = api.verify_integrated_provenance(output / "integrated_evidence_provenance.json")
    assert result["verified_reference_count"] == 5
    assert "assessment_status.taxon_influence_status.source" in result["verified_references"]
    assert {
        "location": "assessment_status.calibration_status.source.protocol",
        "reason": "producer_relative_path",
    } in result["skipped_references"]


@pytest.mark.parametrize("source", ["native", "ora", "assessment", "review", "protocol", "nested"])
@pytest.mark.parametrize("change", ["changed", "missing"])
def test_each_record_reference_rejects_changed_or_missing_bytes(
    complete_record: tuple[Path, dict[str, Any]], source: str, change: str, api: ModuleType
) -> None:
    path, _ = complete_record
    target = path.parent.parent / f"{source}.tsv"
    if change == "changed":
        target.write_text("changed bytes\n", encoding="utf-8")
    else:
        target.unlink()
    with pytest.raises(ValueError, match=r"source|reference|hash|SHA|missing|read"):
        api.verify_integrated_provenance(path)


@pytest.mark.parametrize(
    "bad_reference",
    [
        {"path": "native.tsv"},
        {"path": "native.tsv", "sha256": None},
        {"path": None, "sha256": "a" * 64},
        {"path": 5, "sha256": "a" * 64},
        {"sha256": "a" * 64},
        {"producer_relative_path": "methods/a.tsv", "path": "native.tsv", "sha256": "a" * 64},
        {"producer_relative_path": "", "sha256": "a" * 64},
    ],
)
def test_malformed_reference_is_not_silently_skipped(
    complete_record: tuple[Path, dict[str, Any]], bad_reference: dict[str, object], api: ModuleType
) -> None:
    path, record = complete_record
    record["native_sources"]["arm1"] = bad_reference
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match=r"source|reference|path|SHA|digest"):
        api.verify_integrated_provenance(path)


def test_producer_relative_locators_are_explicit_skips(
    complete_record: tuple[Path, dict[str, Any]], api: ModuleType
) -> None:
    path, record = complete_record
    source = record["assessment_status"]["calibration"]["source"]
    for key in ("review", "protocol"):
        source[key] = {"producer_relative_path": f"methods/{key}.tsv", "sha256": "a" * 64}
    path.write_text(json.dumps(record))
    result = api.verify_integrated_provenance(path)
    assert result["verified_reference_count"] == 4
    skipped = {item["location"]: item["reason"] for item in result["skipped_references"]}
    assert skipped["assessment_status.calibration.source.review"] == "producer_relative_path"
    assert skipped["assessment_status.calibration.source.protocol"] == "producer_relative_path"


def test_legacy_absolute_references_and_relative_protocol_remain_readable(
    complete_record: tuple[Path, dict[str, Any]], api: ModuleType
) -> None:
    path, record = complete_record
    record.pop("path_basis")

    def absolute_refs(value: Any) -> Any:
        if isinstance(value, dict):
            result = {key: absolute_refs(item) for key, item in value.items()}
            if isinstance(result.get("path"), str):
                result["path"] = str(path.parent.parent / result["path"])
            return result
        if isinstance(value, str) and value.startswith("{"):
            return json.dumps(absolute_refs(json.loads(value)))
        return value

    record = absolute_refs(record)
    source = record["assessment_status"]["calibration"]["source"]
    source["protocol"] = "methods/protocol.json"
    source["protocol_sha256"] = "a" * 64
    source["review"] = {"path": "methods/review.tsv", "sha256": "a" * 64}
    path.write_text(json.dumps(record))
    result = api.verify_integrated_provenance(path)
    assert result["verified_reference_count"] == 4
    assert {
        "location": "assessment_status.calibration.source.protocol",
        "reason": "legacy_producer_relative_protocol",
    } in result["skipped_references"]
    assert {
        "location": "assessment_status.calibration.source.review",
        "reason": "legacy_producer_relative_review",
    } in result["skipped_references"]


def test_legacy_absolute_protocol_is_hash_checked(
    complete_record: tuple[Path, dict[str, Any]], api: ModuleType
) -> None:
    path, record = complete_record
    record.pop("path_basis")
    # The record lives beside all source files for legacy relative references.
    path = path.parent.parent / "legacy.json"
    source = record["assessment_status"]["calibration"]["source"]
    source["protocol"] = str(path.parent / "protocol.tsv")
    source["protocol_sha256"] = file_sha256(path.parent / "protocol.tsv")
    source["review"]["path"] = str(path.parent / "review.tsv")
    path.write_text(json.dumps(record))
    result = api.verify_integrated_provenance(path)
    assert result["verified_reference_count"] == 6
    assert "assessment_status.calibration.source.protocol" in result["verified_references"]
    assert "assessment_status.calibration.source.review" in result["verified_references"]
    (path.parent / "protocol.tsv").write_text("changed\n")
    with pytest.raises(ValueError, match=r"hash|SHA|reference"):
        api.verify_integrated_provenance(path)


def test_new_record_cannot_leave_protocol_in_ambiguous_string_form(
    complete_record: tuple[Path, dict[str, Any]], api: ModuleType
) -> None:
    path, record = complete_record
    record["assessment_status"]["calibration"]["source"]["protocol"] = "methods/protocol.json"
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match=r"protocol|scope|reference"):
        api.verify_integrated_provenance(path)


@pytest.mark.parametrize("field", ["review", "protocol"])
@pytest.mark.parametrize("value", [False, 5, [], {}, {"unexpected": "metadata"}])
def test_malformed_optional_assessment_locator_is_refused(
    complete_record: tuple[Path, dict[str, Any]], field: str, value: object, api: ModuleType
) -> None:
    """An optional locator may be absent/null, but a declared value needs a scope."""
    path, record = complete_record
    record["assessment_status"]["calibration"]["source"][field] = value
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match=r"review|protocol|reference|locator"):
        api.verify_integrated_provenance(path)


@pytest.mark.parametrize("field", ["review", "protocol"])
def test_absent_optional_assessment_locator_is_allowed(
    complete_record: tuple[Path, dict[str, Any]], field: str, api: ModuleType
) -> None:
    """Current producers serialize null for unused review/protocol metadata."""
    path, record = complete_record
    record["assessment_status"]["calibration"]["source"][field] = None
    path.write_text(json.dumps(record))
    result = api.verify_integrated_provenance(path)
    assert result["verified_reference_count"] == 5


@pytest.mark.parametrize("bad_record", [[], {}, {"native_sources": []}, {"version": "other"}])
def test_malformed_record_cannot_report_success(
    tmp_path: Path, bad_record: object, api: ModuleType
) -> None:
    path = tmp_path / "record.json"
    path.write_text(json.dumps(bad_record))
    with pytest.raises(ValueError, match=r"record|integrated|version|native"):
        api.verify_integrated_provenance(path)


def test_invalid_basis_is_refused_even_with_no_available_references(
    tmp_path: Path, api: ModuleType
) -> None:
    path = tmp_path / "record.json"
    path.write_text(
        json.dumps(
            {
                "version": "integrated-evidence-v2",
                "path_basis": {"kind": "other", "root": "."},
                "native_sources": {},
                "assessment_status": {},
                "ora_eligibility_source": None,
            }
        )
    )
    with pytest.raises(ValueError, match=r"basis|root|declaration"):
        api.verify_integrated_provenance(path)


def test_present_null_basis_cannot_downgrade_a_record_to_legacy_paths(
    complete_record: tuple[Path, dict[str, Any]], api: ModuleType
) -> None:
    """#1755: legacy compatibility requires an absent declaration, not a null one."""
    path, record = complete_record

    def absolute_refs(value: Any) -> Any:
        if isinstance(value, dict):
            result = {key: absolute_refs(item) for key, item in value.items()}
            if isinstance(result.get("path"), str):
                result["path"] = str(path.parent.parent / result["path"])
            return result
        if isinstance(value, str) and value.startswith("{"):
            return json.dumps(absolute_refs(json.loads(value)))
        return value

    record = absolute_refs(record)
    # All source bytes and hashes match; only the malformed declaration must fail.
    record["path_basis"] = None
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match=r"path_basis|basis|declaration"):
        api.verify_integrated_provenance(path)


def test_cli_reports_counts_and_refuses_bad_hash(
    complete_record: tuple[Path, dict[str, Any]],
    api: ModuleType,
    capsys: pytest.CaptureFixture[str],
) -> None:
    path, _ = complete_record
    assert api.main([str(path)]) == 0
    assert json.loads(capsys.readouterr().out)["verified_reference_count"] == 6
    (path.parent.parent / "native.tsv").write_text("changed\n")
    assert api.main([str(path)]) == 2
    captured = capsys.readouterr()
    assert not captured.out
    assert "hash" in captured.err.lower() or "sha" in captured.err.lower()


def test_external_root_requires_and_accepts_explicit_cli_binding(
    complete_record: tuple[Path, dict[str, Any]],
    tmp_path: Path,
    api: ModuleType,
    capsys: pytest.CaptureFixture[str],
) -> None:
    path, record = complete_record
    root = path.parent.parent
    record["path_basis"] = {"kind": "supplied-root-v1"}
    external = tmp_path / "outside.json"
    external.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="--provenance-root"):
        api.verify_integrated_provenance(external)
    assert api.main([str(external), "--provenance-root", str(root)]) == 0
    assert json.loads(capsys.readouterr().out)["verified_reference_count"] == 6


def test_declared_record_scans_non_reference_metadata_for_local_paths(
    complete_record: tuple[Path, dict[str, Any]], api: ModuleType
) -> None:
    """N25: valid file references cannot mask a leak in an unrelated scalar."""
    path, record = complete_record
    record["producer_note"] = str(Path("/") / "home" / "fixture-person" / "work")
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="Provenance contains"):
        api.verify_integrated_provenance(path)


def test_otherwise_valid_record_requires_supported_version(
    complete_record: tuple[Path, dict[str, Any]], api: ModuleType
) -> None:
    """N26: every reference remains valid, leaving only the version check."""
    path, record = complete_record
    record["version"] = "integrated-evidence-v999"
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="Unsupported integrated provenance record version"):
        api.verify_integrated_provenance(path)


def test_legacy_missing_binding_requires_true_supplied_flag(
    complete_record: tuple[Path, dict[str, Any]], api: ModuleType
) -> None:
    """N28: use a legacy availability shape, with no newer origin cross-check."""
    path, record = complete_record
    del record["native_sources"]["missing"]["supplied"]
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="Missing source does not declare supplied=true"):
        api.verify_integrated_provenance(path)


def test_declared_basis_rejects_unscoped_protocol_in_nested_metadata(
    complete_record: tuple[Path, dict[str, Any]], api: ModuleType
) -> None:
    """N29: nested metadata bypasses the direct assessment-field shape check."""
    path, record = complete_record
    record["extra_metadata"] = json.dumps({"upstream": {"protocol": "methods/protocol.json"}})
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="Protocol locator lacks an explicit reference scope"):
        api.verify_integrated_provenance(path)


def test_missing_checkout_default_is_distinguished_from_absence(
    complete_record: tuple[Path, dict[str, Any]], api: ModuleType
) -> None:
    path, record = complete_record
    record["assessment_status"]["defaulted"] = {
        "source": {
            "path": None,
            "binding": "not_supplied",
            "supplied": False,
            "input_origin": "checkout_default",
        }
    }
    path.write_text(json.dumps(record))
    result = api.verify_integrated_provenance(path)
    assert {
        "location": "assessment_status.defaulted.source",
        "reason": "defaulted_but_missing",
    } in result["skipped_references"]


@pytest.mark.parametrize(
    ("origin", "supplied"),
    [
        (None, False),
        ("unknown", False),
        ("supplied", False),
        ("absent", True),
        ("supplied", 1),
        ("checkout_default", 0),
        ("absent", "False"),
        ("absent", None),
    ],
)
def test_invalid_or_contradictory_input_origin_is_refused(
    complete_record: tuple[Path, dict[str, Any]], origin: object, supplied: object, api: ModuleType
) -> None:
    path, record = complete_record
    record["assessment_status"]["unsupplied"]["source"].update(
        input_origin=origin, supplied=supplied
    )
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match=r"origin|selection"):
        api.verify_integrated_provenance(path)


@pytest.mark.parametrize("origin", ["supplied", "checkout_default"])
def test_present_assessment_is_hash_checked_for_each_input_origin(
    complete_record: tuple[Path, dict[str, Any]], origin: str, api: ModuleType
) -> None:
    path, record = complete_record
    record["assessment_status"]["calibration"]["source"].update(
        input_origin=origin, supplied=origin == "supplied"
    )
    path.write_text(json.dumps(record))
    result = api.verify_integrated_provenance(path)
    assert "assessment_status.calibration.source" in result["verified_references"]
    (path.parent.parent / "assessment.tsv").write_text("changed bytes\n")
    with pytest.raises(ValueError, match=r"hash|SHA"):
        api.verify_integrated_provenance(path)
