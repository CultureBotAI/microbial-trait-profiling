"""Check retained assessment locators before workflow writes (#1841, #1853)."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

from mtp import integrated_evidence
from mtp.influence_assessment import ASSESSMENT_VERSION
from mtp.paths import file_sha256
from mtp.provenance_paths import PortablePaths
from mtp.workflow.provenance import preflight_integrated_sources


@dataclass(frozen=True)
class Inputs:
    root: Path
    record: Path
    manifest: Path
    counts: Path
    ora: Path
    experiments: Path
    dump: Path

    def argv(self) -> list[str]:
        return [
            "--core-run",
            str(self.manifest),
            "--profile-input",
            str(self.counts),
            "--experiments-input",
            str(self.experiments),
            "--full-dump-input",
            str(self.dump),
            "--output-dir",
            str(self.record.parent),
        ]

    def identities(self, kind: str) -> dict[str, str]:
        if kind == "influence":
            return {
                "study_trait_profiles_wide.tsv": file_sha256(self.counts),
                "experiments.tsv": file_sha256(self.experiments),
            }
        return {
            "profile_sha256": file_sha256(self.counts),
            "experiments_sha256": file_sha256(self.experiments),
            "full_dump_sha256": file_sha256(self.dump),
        }


@pytest.fixture
def inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Inputs:
    """Only path-selection inputs; the caller has already verified the core."""
    root = tmp_path / "source-bundle"
    root.mkdir()
    result = Inputs(
        root,
        root / "output/record.json",
        root / "run.json",
        root / "count.tsv",
        root / "ora.tsv",
        root / "experiments.tsv",
        root / "dump.csv",
    )
    for path in (result.manifest, result.counts, result.ora, result.experiments, result.dump):
        path.write_text(path.name + " synthetic identity\n", encoding="utf-8")
    monkeypatch.setenv("TMPDIR", str(tmp_path))
    return result


def _status_payload(kind: str, inputs: Inputs, locator: Path) -> dict[str, object]:
    reference = {"path": str(locator), "sha256": file_sha256(locator)}
    if kind == "influence":
        return {
            "schema_version": ASSESSMENT_VERSION,
            "input_sha256": inputs.identities(kind),
            "procedure_verdicts": {
                "synthetic_procedure": {
                    "verdict": "pending",
                    "basis": "retained input evidence",
                    "applies_to": "all_rows_of_procedure",
                    "metadata": {"nested": reference, "count": 7},
                }
            },
        }
    producer = {
        "code_commit": "synthetic-producer",
        "analysis_source_sha256": hashlib.sha256(b"synthetic analysis source").hexdigest(),
    }
    if kind == "calibration":
        return {
            "version": "synthetic-calibration",
            **producer,
            "inputs": {
                **inputs.identities(kind),
                "analysis_source_sha256": producer["analysis_source_sha256"],
                "protocol_sha256": reference["sha256"],
            },
            "protocol": str(locator),
            "calibration_status": {},
        }
    assert kind == "orientation"
    return {
        "classifier_version": "synthetic-orientation",
        **producer,
        **inputs.identities(kind),
        "class_counts": [{"contrast_class": "synthetic", "experiments": 1, "studies": 1}],
        "review": {**reference, "reviewed_experiments": 1, "agreed": 1, "disagreed": 0},
    }


def _status(inputs: Inputs, kind: str, payload: dict[str, object]) -> tuple[Path, str, str]:
    path = inputs.root / f"{kind}.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    option, key = {
        "calibration": ("--calibration-status", "calibration_status"),
        "orientation": ("--orientation-provenance", "orientation_status"),
        "influence": ("--influence-status", "taxon_influence_status"),
    }[kind]
    return path, option, key


def _preflight(inputs: Inputs, arguments: list[str], *, planned: tuple[Path, ...] = ()) -> None:
    preflight_integrated_sources(
        inputs.record,
        arguments,
        count_profiles=inputs.counts,
        ora_profiles=inputs.ora,
        experiments=inputs.experiments,
        full_dump=inputs.dump,
        root=inputs.root,
        planned_files=planned,
    )


@pytest.mark.parametrize("kind", ["calibration", "orientation", "influence"])
def test_nested_retained_locator_outside_root_is_rejected_before_output(
    inputs: Inputs,
    tmp_path: Path,
    kind: str,
) -> None:
    outside = tmp_path / "external-evidence.tsv"
    outside.write_text("retained evidence\n", encoding="utf-8")
    status, option, key = _status(inputs, kind, _status_payload(kind, inputs, outside))
    parsed = integrated_evidence.read_status(
        status,
        key,
        "#1841",
        expected_input_sha256=inputs.identities(kind),
        supplied=True,
    )
    retained = integrated_evidence._provenance_sources({}, statuses=(parsed,))
    assert outside in retained
    with pytest.raises(
        ValueError, match=r"^Source is outside the declared provenance root"
    ) as error:
        PortablePaths.create(inputs.record, retained, root=inputs.root)
    actual_message = None
    try:
        _preflight(inputs, [*inputs.argv(), option, str(status)])
    except ValueError as actual:
        actual_message = str(actual)
    assert actual_message == str(error.value), "The known nested source must be rejected"
    assert not inputs.record.parent.exists()


@pytest.mark.parametrize("kind", ["calibration", "orientation", "influence"])
def test_nested_retained_locator_inside_root_keeps_real_source_hashes(
    inputs: Inputs,
    kind: str,
) -> None:
    evidence = inputs.root / "retained.tsv"
    evidence.write_text("retained evidence\n", encoding="utf-8")
    status, option, key = _status(inputs, kind, _status_payload(kind, inputs, evidence))
    before = status.read_bytes()
    parsed = integrated_evidence.read_status(
        status,
        key,
        "#1841",
        expected_input_sha256=inputs.identities(kind),
        supplied=True,
    )
    sources = integrated_evidence._provenance_sources({}, statuses=(parsed,))
    assert evidence in sources
    portable = PortablePaths.create(inputs.record, sources, root=inputs.root)
    converted = integrated_evidence._portable_status(parsed, portable)
    assert converted.source["sha256"] == file_sha256(status)
    assert str(inputs.root) not in json.dumps((converted.source, converted.detail))
    try:
        _preflight(inputs, [*inputs.argv(), option, str(status)])
    except ValueError as error:
        raise AssertionError(
            "The producer's accepted nested-source layout must remain accepted"
        ) from error
    assert status.read_bytes() == before
    assert not inputs.record.parent.exists()


def test_identity_mismatch_does_not_import_unread_influence_detail(
    inputs: Inputs,
    tmp_path: Path,
) -> None:
    outside = tmp_path / "unread-evidence.tsv"
    outside.write_text("irrelevant evidence\n", encoding="utf-8")
    payload = _status_payload("influence", inputs, outside)
    payload["input_sha256"] = {
        **inputs.identities("influence"),
        "study_trait_profiles_wide.tsv": hashlib.sha256(b"another corpus").hexdigest(),
    }
    status, option, key = _status(inputs, "influence", payload)
    parsed = integrated_evidence.read_status(
        status,
        key,
        "#1841",
        expected_input_sha256=inputs.identities("influence"),
        supplied=True,
    )
    assert parsed.default.startswith("unreadable")
    assert outside not in integrated_evidence._provenance_sources({}, statuses=(parsed,))
    try:
        _preflight(inputs, [*inputs.argv(), option, str(status)])
    except ValueError as error:
        raise AssertionError(
            "Preflight must not bypass identity checks to import discarded detail"
        ) from error
    assert not inputs.record.parent.exists()


def test_preflight_does_not_invent_locators_from_discarded_fields_or_plain_prose(
    inputs: Inputs,
    tmp_path: Path,
) -> None:
    outside = tmp_path / "not-a-retained-source.tsv"
    outside.write_text("unrelated\n", encoding="utf-8")
    orientation = inputs.root / "orientation.json"
    orientation.write_text(json.dumps({"ignored_metadata": {"path": str(outside)}}))
    coverage = inputs.root / "coverage.json"
    coverage.write_text(json.dumps({"coverage_status": f"missing_reason: no entry for {outside}"}))
    try:
        _preflight(
            inputs,
            [
                *inputs.argv(),
                "--orientation-provenance",
                str(orientation),
                "--coverage-status",
                str(coverage),
            ],
        )
    except ValueError as error:
        raise AssertionError(
            "Discarded metadata and ordinary prose are not source locators"
        ) from error
    assert not inputs.record.parent.exists()


@pytest.mark.parametrize("kind", ["calibration", "orientation"])
def test_relative_upstream_locator_keeps_producer_scope(inputs: Inputs, kind: str) -> None:
    evidence = inputs.root / "retained.tsv"
    evidence.write_text("retained evidence\n", encoding="utf-8")
    payload = _status_payload(kind, inputs, evidence)
    if kind == "calibration":
        payload["protocol"] = "methods/protocol.json"
    else:
        payload["review"] = {"path": "methods/review.tsv", "sha256": file_sha256(evidence)}
    status, option, key = _status(inputs, kind, payload)
    parsed = integrated_evidence.read_status(
        status,
        key,
        "#1841",
        expected_input_sha256=inputs.identities(kind),
        supplied=True,
    )
    sources = integrated_evidence._provenance_sources({}, statuses=(parsed,))
    assert set(sources) == {status}
    converted = integrated_evidence._portable_status(
        parsed, PortablePaths.create(inputs.record, sources, root=inputs.root)
    )
    locator = converted.source["protocol" if kind == "calibration" else "review"]
    assert isinstance(locator, dict)
    assert "producer_relative_path" in locator
    assert "path" not in locator
    _preflight(inputs, [*inputs.argv(), option, str(status)])
    assert not inputs.record.parent.exists()


def test_future_assessment_file_is_not_opened_for_unknown_metadata(
    inputs: Inputs,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    future = inputs.root / "future/status.json"
    original = Path.open

    def guarded(path: Path, *args: Any, **kwargs: Any) -> Any:
        assert path != future, "A future assessment's contents cannot be inspected before it exists"
        return original(path, *args, **kwargs)

    def unnecessary_hash(_path: Path) -> str:
        raise AssertionError("Absent assessments need no core input fingerprints in preflight")

    monkeypatch.setattr(Path, "open", guarded)
    monkeypatch.setattr(integrated_evidence, "_sha256", unnecessary_hash)
    _preflight(
        inputs,
        [*inputs.argv(), "--coverage-status", str(future)],
        planned=(future,),
    )
    assert not future.exists()
    assert not inputs.record.parent.exists()


def test_workflow_preflight_requires_core_but_standalone_parser_still_accepts_explicit_inputs(
    inputs: Inputs,
) -> None:
    arguments = inputs.argv()[2:]
    standalone = integrated_evidence._build_parser().parse_args(arguments)
    assert standalone.core_run is None
    assert standalone.profile_input == inputs.counts
    refusal = None
    try:
        _preflight(inputs, arguments)
    except ValueError as error:
        refusal = str(error)
    assert refusal is not None, "Workflow provenance preflight requires --core-run"
    assert re.search(r"requires --core-run", refusal)
    assert not inputs.record.parent.exists()
