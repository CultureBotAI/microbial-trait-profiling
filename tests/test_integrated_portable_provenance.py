"""Portable integrated source identities without regenerating tracked data (#1664)."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

import test_integrated_evidence as fixtures
import test_integrated_native_sources as native_fixtures
from mtp import integrated_evidence, provenance_paths
from mtp.integrated_provenance import verify_integrated_provenance
from mtp.paths import file_sha256
from test_no_absolute_local_paths import ABSOLUTE_LOCAL_PATH_PATTERN


def _arguments(root: Path, output: Path) -> list[str]:
    return [
        "--profile-input",
        str(root / "study_trait_profiles_wide.tsv"),
        "--experiments-input",
        str(root / "experiments.tsv"),
        "--full-dump-input",
        str(root / "dump.csv"),
        "--calibration-status",
        str(root / "calibration.json"),
        "--orientation-provenance",
        str(root / "orientation.json"),
        "--sensitivity-status",
        str(root / "missing-sensitivity.json"),
        "--influence-status",
        str(root / "influence.json"),
        "--output-dir",
        str(output),
        "--arm1-pooled-stats",
        str(root / "pooled.tsv"),
    ]


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
    arguments = _arguments(root, output)
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
        "supplied": True,
        "binding": "present",
    }
    statuses = provenance["assessment_status"]
    source = statuses["orientation_status"]["source"]
    assert source["path"] == "orientation.json"
    assert source["review"]["path"] == "review.tsv"
    assert source["review"]["sha256"] == file_sha256(root / "review.tsv")
    assert statuses["calibration_status"]["source"]["binding"] == "matched"
    assert statuses["sensitivity_status"]["source"]["path"] is None
    assert statuses["sensitivity_status"]["source"]["binding"] == "missing"
    assert statuses["sensitivity_status"]["source"]["supplied"] is True
    assert provenance["native_sources"]["full_screen"] == {
        "path": None,
        "sha256": None,
        "supplied": False,
        "binding": "not_supplied",
    }
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


@pytest.mark.parametrize("key", ["coverage_status", "taxon_influence_status"])
def test_missing_assessment_is_distinct_from_not_supplied(tmp_path: Path, key: str) -> None:
    absent = integrated_evidence.read_status(None, key, "#1757")
    missing = integrated_evidence.read_status(tmp_path / "absent.json", key, "#1757")
    assert absent.source == {"path": None, "binding": "not_supplied"}
    assert missing.source == {"path": None, "binding": "missing", "supplied": True}
    if key == "taxon_influence_status":
        assert "not supplied" in absent.default
        assert "supplied influence evidence missing" in missing.default


def test_missing_native_records_do_not_add_a_location_to_root_inference(tmp_path: Path) -> None:
    from mtp.provenance_paths import PortablePaths

    helper = getattr(integrated_evidence, "_native_source_record", None)
    assert callable(helper), "Native source availability must be explicit"
    root = tmp_path / "sources"
    root.mkdir()
    paths = PortablePaths.create(root / "record.json", [])
    missing = tmp_path / "outside" / "missing.tsv"
    assert helper(missing, supplied=True, paths=paths) == {
        "path": None,
        "sha256": None,
        "supplied": True,
        "binding": "missing",
    }
    assert helper(missing, supplied=False, paths=paths) == {
        "path": None,
        "sha256": None,
        "supplied": False,
        "binding": "not_supplied",
    }


def test_unreadable_supplied_assessment_retains_source_identity(
    portable_case: tuple[Path, Path],
) -> None:
    root, output = portable_case
    record = json.loads((output / "integrated_evidence_provenance.json").read_text())
    source = record["assessment_status"]["taxon_influence_status"]["source"]
    assert source["path"] == "influence.json"
    assert source["sha256"] == file_sha256(root / "influence.json")
    assert source["binding"] == "supplied"


@pytest.mark.parametrize("absolute", [False, True])
def test_calibration_protocol_locator_has_explicit_scope(tmp_path: Path, absolute: bool) -> None:
    from mtp.provenance_paths import PortablePaths

    protocol = tmp_path / "protocol.json"
    protocol.write_text("{}\n")
    digest = file_sha256(protocol)
    status = integrated_evidence.AssessmentStatus(
        "#182",
        {},
        "pending",
        source={
            "protocol": str(protocol) if absolute else "methods/protocol.json",
            "protocol_sha256": digest,
        },
    )
    converted = integrated_evidence._portable_status(
        status, PortablePaths.create(tmp_path / "record.json", [protocol])
    )
    assert converted.source["protocol"] == {
        "path" if absolute else "producer_relative_path": "protocol.json"
        if absolute
        else "methods/protocol.json",
        "sha256": digest,
    }
    assert converted.source["protocol_sha256"] == digest


def test_declared_provenance_root_is_used_by_integrated_cli(
    portable_case: tuple[Path, Path],
) -> None:
    root, _ = portable_case
    output = root / "declared" / "out"
    try:
        result = integrated_evidence.main(
            [*_arguments(root, output), "--provenance-root", str(root.parent)]
        )
    except SystemExit as error:
        raise AssertionError(
            "The integrated CLI must accept an explicit provenance root"
        ) from error
    assert result == 0
    record = json.loads((output / "integrated_evidence_provenance.json").read_text())
    assert record["path_basis"] == {"kind": "root-relative-v1", "root": "../../.."}
    assert record["native_sources"]["arm1_pooled_stats"]["path"] == "original/pooled.tsv"


def test_final_provenance_backstop_rejects_an_unencoded_field(
    portable_case: tuple[Path, Path],
) -> None:
    root, _ = portable_case
    output = root / "refused"
    leaked = str(Path("/") / "home" / "fixture-user" / "revision")
    with pytest.raises(ValueError, match="Provenance contains"):
        integrated_evidence.main([*_arguments(root, output), "--code-commit", leaked])
    assert not (output / "integrated_evidence_provenance.json").exists()


def test_explicit_source_root_can_keep_report_outside_source_tree(
    portable_case: tuple[Path, Path],
) -> None:
    root, _ = portable_case
    output = root.parent / "external-report"
    try:
        result = integrated_evidence.main(
            [*_arguments(root, output), "--provenance-root", str(root)]
        )
    except SystemExit as error:
        raise AssertionError("The integrated CLI must accept an external source root") from error
    assert result == 0
    text = (output / "integrated_evidence_provenance.json").read_text()
    record = json.loads(text)
    assert record["path_basis"] == {"kind": "supplied-root-v1"}
    assert record["native_sources"]["arm1_pooled_stats"]["path"] == "pooled.tsv"
    assert str(root) not in text


def test_external_output_requires_an_explicit_safe_source_root(
    portable_case: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    root, _ = portable_case
    output = root.parent / "temporary-output" / "report"
    # Model separate checkout and temp trees without creating machine-root files.
    # Their common ancestor is a private boundary and cannot be published.
    monkeypatch.setattr(provenance_paths, "_private_roots", lambda: {root.parent.resolve()})
    record_path = output / "integrated_evidence_provenance.json"
    with pytest.raises(ValueError, match="below a home or temporary root"):
        provenance_paths.PortablePaths.create(record_path, [root / "pooled.tsv"])
    with pytest.raises(ValueError, match="Choose --provenance-root") as error:
        integrated_evidence.main(_arguments(root, output))
    assert str(root) not in str(error.value)
    assert str(output) not in str(error.value)
    assert not output.exists()
    result = integrated_evidence.main([*_arguments(root, output), "--provenance-root", str(root)])
    assert result == 0
    text = record_path.read_text()
    record = json.loads(text)
    assert record["path_basis"] == {"kind": "supplied-root-v1"}
    assert record["native_sources"]["arm1_pooled_stats"]["path"] == "pooled.tsv"
    assert str(root) not in text
    assert str(output) not in text
    with pytest.raises(ValueError, match="requires --provenance-root"):
        verify_integrated_provenance(record_path)
    assert (
        verify_integrated_provenance(record_path, provenance_root=root)["verified_reference_count"]
        == 5
    )
    (root / "review.tsv").write_text("changed reviewed bytes\n")
    with pytest.raises(ValueError, match=r"hash|SHA"):
        verify_integrated_provenance(record_path, provenance_root=root)


@pytest.mark.parametrize("account", ["root", "ci", "app", "sample", "variance"])
def test_real_producer_accepts_account_names_in_nonpath_record_text(
    portable_case: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch, account: str
) -> None:
    root, _ = portable_case
    monkeypatch.setattr(provenance_paths, "_account_names", lambda: {account})
    output = root / "ordinary-record"
    try:
        result = integrated_evidence.main(_arguments(root, output))
    except ValueError as error:
        raise AssertionError("An account name in ordinary record prose must be accepted") from error
    assert result == 0
    record_path = output / "integrated_evidence_provenance.json"
    assert verify_integrated_provenance(record_path)["verified_reference_count"] == 5


@pytest.mark.parametrize("template", ["-Users-{}-proj", "pytest-of-{}"])
def test_real_producer_refuses_short_writer_names_in_source_directories(
    portable_case: tuple[Path, Path],
    monkeypatch: pytest.MonkeyPatch,
    template: str,
) -> None:
    root, _ = portable_case
    directory = root / template.format("ANA")
    directory.mkdir()
    native = directory / "pooled.tsv"
    native.write_bytes((root / "pooled.tsv").read_bytes())
    output = root / "refused-short-account"
    arguments = _arguments(root, output)
    arguments[arguments.index("--arm1-pooled-stats") + 1] = str(native)
    monkeypatch.setattr(provenance_paths, "_account_names", lambda: {"ana"})
    with pytest.raises(ValueError, match="local account name"):
        integrated_evidence.main(arguments)
    assert not output.exists()


@pytest.mark.parametrize("locator", ["review", "protocol"])
def test_declared_upstream_file_locator_keeps_its_format_suffix(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    locator: str,
) -> None:
    monkeypatch.setattr(provenance_paths, "_account_names", lambda: {"json"})
    value = "methods/assessment.json"
    source: dict[str, object] = (
        {"review": {"path": value}} if locator == "review" else {"protocol": value}
    )
    status = integrated_evidence.AssessmentStatus("#1778", {}, "pending", source=source)
    paths = provenance_paths.PortablePaths.create(tmp_path / "record.json", [])
    try:
        converted = integrated_evidence._portable_status(status, paths)
    except ValueError as error:
        raise AssertionError("An upstream file format must not be treated as a login") from error
    reference = converted.source[locator]
    assert isinstance(reference, dict)
    assert reference["producer_relative_path"] == value


@pytest.mark.parametrize("locator", ["review", "protocol"])
@pytest.mark.parametrize("suffix", ["/", "/."])
@pytest.mark.parametrize("account", ["tsv", "json"])
def test_upstream_directory_spelling_cannot_use_a_file_suffix_exception(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    locator: str,
    suffix: str,
    account: str,
) -> None:
    monkeypatch.setattr(provenance_paths, "_account_names", lambda: {account})
    value = "methods/data." + account + suffix
    source: dict[str, object] = (
        {"review": {"path": value}} if locator == "review" else {"protocol": value}
    )
    status = integrated_evidence.AssessmentStatus("#1778", {}, "pending", source=source)
    paths = provenance_paths.PortablePaths.create(tmp_path / "record.json", [])
    with pytest.raises(ValueError, match="local account name"):
        integrated_evidence._portable_status(status, paths)


def test_verification_does_not_apply_the_readers_account_name(
    portable_case: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    _, output = portable_case
    record_path = output / "integrated_evidence_provenance.json"
    original = record_path.read_bytes()
    # The writer's review.tsv is unrelated to a reader whose login is "review".
    monkeypatch.setattr(provenance_paths, "_account_names", lambda: {"review"})
    try:
        result = verify_integrated_provenance(record_path)
    except ValueError as error:
        raise AssertionError(
            "Reader account names must not reject portable source records"
        ) from error
    assert result["verified_reference_count"] == 5
    assert record_path.read_bytes() == original


@pytest.mark.parametrize("locator", ["review", "protocol"])
def test_producer_relative_locators_reject_writer_account_tokens(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, locator: str
) -> None:
    monkeypatch.setattr(provenance_paths, "_account_names", lambda: {"fixtureperson"})
    value = "methods/pytest-of-FIXTUREPERSON/assessment.tsv"
    source: dict[str, object] = (
        {"review": {"path": value}} if locator == "review" else {"protocol": value}
    )
    status = integrated_evidence.AssessmentStatus("#1773", {}, "pending", source=source)
    paths = provenance_paths.PortablePaths.create(tmp_path / "record.json", [])
    with pytest.raises(ValueError, match="local account name"):
        integrated_evidence._portable_status(status, paths)


def test_inference_still_refuses_sources_split_across_a_private_boundary(
    portable_case: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    root, _ = portable_case
    external = root.parent / "other-source-tree"
    external.mkdir()
    review = external / "review.tsv"
    review.write_text("other reviewed bytes\n")
    orientation = root / "orientation.json"
    payload = json.loads(orientation.read_text())
    payload["review"].update(path=str(review), sha256=file_sha256(review))
    orientation.write_text(json.dumps(payload))
    monkeypatch.setattr(provenance_paths, "_private_roots", lambda: {root.parent.resolve()})
    output = root / "refused-cross-tree"
    with pytest.raises(ValueError, match="below a home or temporary root"):
        integrated_evidence.main(_arguments(root, output))
    assert not output.exists()


def test_invalid_supplied_root_is_not_replaced_by_inference(
    portable_case: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    root, _ = portable_case
    monkeypatch.setattr(provenance_paths, "_private_roots", lambda: {root.parent.resolve()})
    output = root / "refused-explicit-root"
    with pytest.raises(ValueError, match="below a home or temporary root"):
        integrated_evidence.main([*_arguments(root, output), "--provenance-root", str(root.parent)])
    assert not output.exists()
