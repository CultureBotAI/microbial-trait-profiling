"""Portable integrated source identities without regenerating tracked data (#1664)."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

import test_integrated_evidence as fixtures
import test_integrated_native_sources as native_fixtures
from mtp import integrated_evidence
from mtp.paths import file_sha256
from test_no_absolute_local_paths import ABSOLUTE_LOCAL_PATH_PATTERN


@pytest.fixture
def portable_case(tmp_path: Path, request: pytest.FixtureRequest) -> tuple[Path, Path]:
    root = tmp_path / "original"
    root.mkdir()
    profiles, experiments = fixtures._corpus(root)
    dump = root / "dump.csv"
    dump.write_text("unversioned synthetic dump\n", encoding="utf-8")
    review = root / "review.tsv"
    review.write_text("reviewed\n", encoding="utf-8")
    expected = {
        "profile_sha256": file_sha256(profiles),
        "experiments_sha256": file_sha256(experiments),
        "full_dump_sha256": file_sha256(dump),
    }
    native = root / "pooled.tsv"
    native_rows = native_fixtures.composition_rows()
    for row in native_rows:
        row.update(
            profile_sha256=expected["profile_sha256"],
            experiments_sha256=expected["experiments_sha256"],
            pair_count="16",
            profile_count="32",
            study_count="16",
            n_nonzero_experiments="15",
        )
    native_fixtures.write_source(native, native_rows)
    calibration = root / "calibration.json"
    calibration.write_text(
        json.dumps(fixtures._bound_status_payload("calibration_status", expected))
    )
    orientation = root / "orientation.json"
    payload = fixtures._bound_status_payload("orientation_status", expected)
    payload["review"] = {
        "path": str(review),
        "sha256": file_sha256(review),
        "reviewed_experiments": 1,
        "agreed": 1,
        "disagreed": 0,
    }
    orientation.write_text(json.dumps(payload))
    influence = root / "influence.json"
    influence.write_text(
        json.dumps({"native_manifest_sha256": "a" * 64, "tables": {}})
        if getattr(request, "param", None) == "missing_native"
        else "invalid JSON",
        encoding="utf-8",
    )
    output = root / "out"
    arguments = [
        "--profile-input",
        str(profiles),
        "--experiments-input",
        str(experiments),
        "--full-dump-input",
        str(dump),
        "--calibration-status",
        str(calibration),
        "--orientation-provenance",
        str(orientation),
        "--sensitivity-status",
        str(root / "missing-sensitivity.json"),
        "--influence-status",
        str(influence),
        "--output-dir",
        str(output),
        "--arm1-pooled-stats",
        str(native),
    ]
    try:
        result = integrated_evidence.main(arguments)
    except ValueError as exc:
        raise AssertionError(
            "The valid synthetic integration must publish portable outputs"
        ) from exc
    assert result == 0
    assert Path(integrated_evidence.__file__).resolve() == (
        Path(__file__).resolve().parents[1] / "src/mtp/integrated_evidence.py"
    )
    return root, output


def test_saved_json_and_assessments_are_portable(portable_case: tuple[Path, Path]) -> None:
    root, output = portable_case
    provenance_path = output / "integrated_evidence_provenance.json"
    for name in (provenance_path.name, "integrated_evidence_assessments.tsv"):
        text = (output / name).read_text(encoding="utf-8")
        assert not ABSOLUTE_LOCAL_PATH_PATTERN.search(text)
        assert str(root) not in text
    provenance = json.loads(provenance_path.read_text())
    assert provenance["path_basis"] == {"kind": "root-relative-v1", "root": ".."}
    assert provenance["native_sources"]["arm1_pooled_stats"] == {
        "path": "pooled.tsv",
        "sha256": file_sha256(root / "pooled.tsv"),
    }
    statuses = provenance["assessment_status"]
    source = statuses["orientation_status"]["source"]
    assert source["path"] == "orientation.json"
    assert source["review"]["path"] == "review.tsv"
    assert source["review"]["sha256"] == file_sha256(root / "review.tsv")
    assert statuses["calibration_status"]["source"]["binding"] == "matched"
    assert statuses["sensitivity_status"]["source"]["path"] is None
    with (output / "integrated_evidence_assessments.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sources = {row["state_value"] for row in rows if row["state_variable"] == "path"}
    assert {"calibration.json", "orientation.json", "influence.json"} <= sources
    assert any(row["state_variable"] == "sha256" for row in rows)


def test_declared_sources_relocate_and_detect_changed_bytes(
    portable_case: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    from mtp.provenance_paths import verify_reference

    root, _ = portable_case
    relocated = tmp_path / "relocated"
    shutil.copytree(root, relocated)
    root.rename(tmp_path / "original-hidden")
    path = relocated / "out/integrated_evidence_provenance.json"
    record = json.loads(path.read_text())
    source = record["assessment_status"]["orientation_status"]["source"]
    assert verify_reference(path, source, path_basis=record["path_basis"]) == (
        relocated / "orientation.json"
    )
    assert verify_reference(path, source["review"], path_basis=record["path_basis"]) == (
        relocated / "review.tsv"
    )
    (relocated / "review.tsv").write_text("changed!\n")
    with pytest.raises(ValueError, match=r"hash|SHA|digest"):
        verify_reference(path, source["review"], path_basis=record["path_basis"])


@pytest.mark.parametrize("portable_case", ["missing_native"], indirect=True)
def test_unreadable_influence_error_stays_portable_and_permissive(
    portable_case: tuple[Path, Path],
) -> None:
    root, output = portable_case
    text = (output / "integrated_evidence_provenance.json").read_text()
    record = json.loads(text)
    status = record["assessment_status"]["taxon_influence_status"]
    assert status["status_for_a_procedure_it_is_silent_about"].startswith("unreadable #190")
    assert "run_manifest.json" in status["status_for_a_procedure_it_is_silent_about"]
    assert str(root) not in text
    assert not ABSOLUTE_LOCAL_PATH_PATTERN.search(text)


def test_historical_relative_review_path_keeps_its_producer_scope(tmp_path: Path) -> None:
    from mtp.provenance_paths import PortablePaths

    status = integrated_evidence.AssessmentStatus(
        "#188", {}, "pending", source={"review": {"path": "methods/review.tsv", "sha256": None}}
    )
    converted = integrated_evidence._portable_status(
        status, PortablePaths.create(tmp_path / "record.json", [])
    )
    assert converted.source["review"] == {
        "producer_relative_path": "methods/review.tsv",
        "sha256": None,
    }


def test_json_encoded_influence_detail_rebases_only_locations(tmp_path: Path) -> None:
    from mtp.provenance_paths import PortablePaths

    evidence = tmp_path / "evidence.json"
    evidence.write_text("retained bytes\n")
    detail = {"source": {"path": str(evidence), "sha256": file_sha256(evidence)}, "count": 7}
    status = integrated_evidence.AssessmentStatus(
        "#190", {}, "pending", detail=(("execution", "native", "metadata", json.dumps(detail)),)
    )
    converted = integrated_evidence._portable_status(
        status, PortablePaths.create(tmp_path / "record.json", [evidence])
    )
    actual = json.loads(converted.detail[0][3])
    assert actual == {
        "source": {"path": "evidence.json", "sha256": file_sha256(evidence)},
        "count": 7,
    }
