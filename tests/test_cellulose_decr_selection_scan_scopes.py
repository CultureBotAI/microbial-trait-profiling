"""Relocated v2 scan evidence remains fully verified when scopes are paired (#1080)."""

from __future__ import annotations

import dataclasses
import hashlib
import io
import json
import shutil
from collections.abc import Iterator
from datetime import timedelta
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from mtp.cellulose_decr import selection_scan as scan
from mtp.cellulose_decr import selection_scan_scopes as scopes
from mtp.datamodel.workflow_pydantic import ConsumerRecord, WorkflowRun
from mtp.workflow import consumers
from mtp.workflow.archive import ArchiveResolver, create_archive
from test_cellulose_decr_selection_scan import cores as cores
from test_cellulose_decr_selection_scan import scans as scans


@pytest.fixture(scope="module")
def loaded_scans(scans: dict[str, Path]) -> dict[str, scopes.ScanRun]:
    return {scope: scopes.load_run(directory) for scope, directory in scans.items()}


@pytest.fixture(scope="module")
def retained(
    scans: dict[str, Path], cores: dict[str, Path], tmp_path_factory: pytest.TempPathFactory
) -> ArchiveResolver:
    root = tmp_path_factory.mktemp("relocated-scan-evidence")
    archived = create_archive(
        root / "original-archive-location",
        core_manifests=[cores[f"tier1_2_{scope}"] for scope in ("off", "phenotype")],
        consumer_records=[directory / scan.CONSUMER_RECORD for directory in scans.values()],
        retained_paths=[],
    )
    relocated = root / "restored-at-another-location"
    archived.parent.rename(relocated)
    return ArchiveResolver(relocated / archived.name)


@pytest.fixture
def isolated_retained(
    retained: ArchiveResolver, scans: dict[str, Path], tmp_path: Path
) -> Iterator[ArchiveResolver]:
    """Each byte mutation gets an independent archive; shared source fixtures stay immutable."""
    originals = (
        retained.files
        | {retained.root / "archive.json"}
        | {path for directory in scans.values() for path in directory.rglob("*") if path.is_file()}
    )
    before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in originals}
    restored = tmp_path / "private-archive"
    shutil.copytree(retained.root, restored)
    yield ArchiveResolver(restored / "archive.json")
    assert {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in originals} == before


def forbid_original_reads(monkeypatch: pytest.MonkeyPatch, retained: ArchiveResolver) -> None:
    original_open = Path.open

    def retained_only(path: Path, *args: Any, **kwargs: Any) -> Any:
        if any(path.absolute().is_relative_to(original) for original in retained.mappings):
            pytest.fail(f"A relocated scan read an original path: {path}")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", retained_only)


@pytest.mark.parametrize("scope", ["off", "phenotype"])
def test_load_run_resolves_every_output_and_bound_core_without_live_reads(
    scans: dict[str, Path],
    loaded_scans: dict[str, scopes.ScanRun],
    retained: ArchiveResolver,
    monkeypatch: pytest.MonkeyPatch,
    scope: str,
) -> None:
    before = loaded_scans[scope]
    forbid_original_reads(monkeypatch, retained)
    after = scopes.load_run(scans[scope], resolver=retained.resolve)
    assert after.directory == scans[scope]
    assert after.core_run_id == before.core_run_id
    assert after.summary == before.summary
    assert after.provenance == before.provenance
    assert after.files == before.files
    assert after.protocol == before.protocol
    pd.testing.assert_frame_equal(after.traits, before.traits)


def test_resolver_refusal_does_not_fall_back_to_existing_original_core(
    scans: dict[str, Path], cores: dict[str, Path], retained: ArchiveResolver
) -> None:
    original_core = cores["tier1_2_off"]
    assert original_core.is_file()

    def without_core(path: Path) -> Path:
        if path == original_core:
            raise ValueError("Core manifest deliberately not retained")
        return retained.resolve(path)

    with pytest.raises(ValueError, match="deliberately not retained"):
        scopes.load_run(scans["off"], resolver=without_core)


def one_byte_change(name: str, raw: bytes) -> bytes:
    """Keep valid syntax and byte count so rejection must depend on the output digest."""
    if name in {"summary", "provenance"}:
        offset = raw.index(b"\n ") + 1
        changed = raw[:offset] + b"\t" + raw[offset + 1 :]
        assert json.loads(changed) == json.loads(raw)
    else:
        header, first_row, *_ = raw.splitlines()
        position = header.split(b"\t").index(b"trait_label")
        fields = first_row.split(b"\t")
        label = fields[position]
        character = next(
            index
            for index, value in enumerate(label)
            if chr(value).isascii() and chr(value).isalpha()
        )
        offset = len(header) + 1 + sum(len(field) + 1 for field in fields[:position]) + character
        replacement = bytes([label[character]]).swapcase()
        changed = raw[:offset] + replacement + raw[offset + 1 :]
        before = pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str, keep_default_na=False)
        after = pd.read_csv(io.BytesIO(changed), sep="\t", dtype=str, keep_default_na=False)
        expected = before.copy()
        original_label = str(before.loc[0, "trait_label"])
        expected.loc[0, "trait_label"] = (
            original_label[:character]
            + replacement.decode("ascii")
            + original_label[character + 1 :]
        )
        pd.testing.assert_frame_equal(after, expected)
    assert len(changed) == len(raw)
    assert sum(left != right for left, right in zip(raw, changed, strict=True)) == 1
    return changed


@pytest.mark.parametrize("name", scopes.READ)
def test_loaded_bytes_must_still_match_the_verified_output_inventory(
    scans: dict[str, Path],
    retained: ArchiveResolver,
    isolated_retained: ArchiveResolver,
    monkeypatch: pytest.MonkeyPatch,
    name: str,
) -> None:
    source = scans["off"] / scan.OUTPUTS[name]
    retained_source = retained.resolve(source)
    path = isolated_retained.resolve(source)
    assert path.is_relative_to(isolated_retained.root)
    assert not path.is_relative_to(retained.root)
    original_bytes = path.read_bytes()
    assert source.read_bytes() == retained_source.read_bytes() == original_bytes
    changed = one_byte_change(name, original_bytes)
    original_verify = consumers.verify_consumer_record

    def verify_then_change(*args: Any, **kwargs: Any) -> Any:
        record = original_verify(*args, **kwargs)
        path.write_bytes(changed)
        return record

    monkeypatch.setattr(consumers, "verify_consumer_record", verify_then_change)
    try:
        with pytest.raises(ValueError, match=r"(?i)changed|hash|bound|bytes|digest"):
            scopes.load_run(scans["off"], resolver=isolated_retained.resolve)
    finally:
        assert source.read_bytes() == retained_source.read_bytes() == original_bytes


def test_loaded_consumer_record_must_equal_the_verified_record(
    scans: dict[str, Path],
    isolated_retained: ArchiveResolver,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = isolated_retained.resolve(scans["off"] / scan.CONSUMER_RECORD)
    original_verify = consumers.verify_consumer_record

    def verify_then_change(*args: Any, **kwargs: Any) -> Any:
        record = original_verify(*args, **kwargs)
        changed = record.model_copy(update={"created_at": record.created_at + timedelta(seconds=1)})
        path.write_text(changed.model_dump_json(), encoding="utf-8")
        return record

    monkeypatch.setattr(consumers, "verify_consumer_record", verify_then_change)
    with pytest.raises(ValueError, match="consumer record changed after verification"):
        scopes.load_run(scans["off"], resolver=isolated_retained.resolve)


def test_verified_snapshots_must_still_declare_v2_after_preflight(
    scans: dict[str, Path],
    isolated_retained: ArchiveResolver,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original_verify = consumers.verify_consumer_record
    original = scans["off"]
    consumer_path = isolated_retained.resolve(original / scan.CONSUMER_RECORD)
    paths = {
        name: isolated_retained.resolve(original / scan.OUTPUTS[name])
        for name in ("summary", "provenance")
    }
    protocol = scan.load_protocol(scan.PROTOCOL_V1_PATH)
    rollup = scan.scope_record("off", "recorded", protocol)
    del rollup["scope_source"]

    def change_then_verify(*args: Any, **kwargs: Any) -> Any:
        # The loader has already accepted the original v2 preflight by this call.
        summary = json.loads(paths["summary"].read_bytes())
        provenance = json.loads(paths["provenance"].read_bytes())
        summary.update(protocol=scan.VERSION_V1, **{scan.STRAIN_ROLLUP: rollup})
        provenance.update(
            protocol=scan.protocol_record(scan.VERSION_V1), **{scan.STRAIN_ROLLUP: rollup}
        )
        assert scan.recorded_protocol(summary, provenance) == protocol
        for name, value in (("summary", summary), ("provenance", provenance)):
            paths[name].write_text(json.dumps(value), encoding="utf-8")
        record = ConsumerRecord.model_validate_json(consumer_path.read_bytes())
        outputs = []
        for output in record.outputs:
            path = isolated_retained.resolve(Path(output.path))
            if path in paths.values():
                raw = path.read_bytes()
                output = output.model_copy(
                    update={"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}
                )
            outputs.append(output)
        consumer_path.write_text(
            record.model_copy(update={"outputs": outputs}).model_dump_json(), encoding="utf-8"
        )
        return original_verify(*args, **kwargs)

    monkeypatch.setattr(consumers, "verify_consumer_record", change_then_verify)
    with pytest.raises(ValueError, match=scan.VERSION_V1):
        scopes.load_run(original, resolver=isolated_retained.resolve)


@pytest.mark.parametrize("version_record", ["summary", "provenance", "both"])
def test_v1_is_named_and_refused_before_scientific_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, version_record: str
) -> None:
    versions = {
        name: scan.VERSION_V1 if version_record in {name, "both"} else scan.VERSION
        for name in ("summary", "provenance")
    }
    (tmp_path / scan.OUTPUTS["summary"]).write_text(json.dumps({"protocol": versions["summary"]}))
    (tmp_path / scan.OUTPUTS["provenance"]).write_text(
        json.dumps({"protocol": {"version": versions["provenance"]}})
    )

    def forbid_verification(*_args: Any, **_kwargs: Any) -> Any:
        pytest.fail("Historical v1 must be refused before scientific core verification")

    monkeypatch.setattr(consumers, "verify_consumer_record", forbid_verification)
    with pytest.raises(ValueError, match=scan.VERSION_V1):
        scopes.load_run(tmp_path)


@pytest.mark.parametrize("field", ["environment", "methods", "chunk"])
@pytest.mark.parametrize("missing", [False, True], ids=["different", "missing"])
def test_scopes_require_the_same_numerical_environment_methods_and_chunk(
    loaded_scans: dict[str, scopes.ScanRun], field: str, missing: bool
) -> None:
    off = loaded_scans["off"]
    provenance = json.loads(json.dumps(off.provenance))
    if missing:
        del provenance[field]
    elif field == "environment":
        provenance[field]["numpy"] = "different-version"
    elif field == "methods":
        provenance[field]["sha256"] = "0" * 64
    else:
        provenance[field] += 1
    changed = dataclasses.replace(off, provenance=provenance)
    with pytest.raises(ValueError, match=field):
        scopes.pair_scopes([loaded_scans["phenotype"], changed])


def bound_records(
    cores: dict[str, Path], protocol: scan.Protocol, *, source: str | None
) -> tuple[WorkflowRun, dict[str, Any], dict[str, Any], str]:
    raw = cores["tier1_2_off"].read_bytes()
    run = WorkflowRun.model_validate_json(raw)
    digest = hashlib.sha256(raw).hexdigest()
    rollup: dict[str, Any] = {"scope": "off"}
    if source is not None:
        rollup["scope_source"] = source
    summary = {"core_run": {"run_id": run.run_id}, scan.STRAIN_ROLLUP: rollup}
    provenance = {"core_run": {"run_id": run.run_id, "manifest_sha256": digest}}
    assert scan.core_scope(run, protocol) == ("off", "recorded")
    return run, summary, provenance, digest


def test_v2_requires_scope_source_even_if_presentation_flag_is_false(
    cores: dict[str, Path],
) -> None:
    protocol = dataclasses.replace(scan.load_protocol(), side_by_side=False)
    run, summary, provenance, digest = bound_records(cores, protocol, source=None)
    with pytest.raises(ValueError, match="bound core"):
        scopes.check_bound_core(Path("off"), summary, provenance, run, digest, protocol)


@pytest.mark.parametrize("version", [scan.VERSION, scan.VERSION_V1])
def test_scope_source_validation_is_keyed_by_version_not_presentation_flag(
    cores: dict[str, Path], version: str
) -> None:
    current = version == scan.VERSION
    path = scan.PROTOCOL_PATH if current else scan.PROTOCOL_V1_PATH
    protocol = dataclasses.replace(scan.load_protocol(path), side_by_side=not current)
    run, summary, provenance, digest = bound_records(
        cores, protocol, source="recorded" if current else None
    )
    scopes.check_bound_core(Path("off"), summary, provenance, run, digest, protocol)


def test_open_tier1_2_exposes_named_scope_and_scope_source(cores: dict[str, Path]) -> None:
    opened = scan.open_tier1_2(cores["tier1_2_off"])
    assert opened.bundle.manifest == cores["tier1_2_off"]
    assert opened.scope == "off"
    assert opened.scope_source == "recorded"
    assert tuple(opened) == (opened.bundle, opened.scope, opened.scope_source)


@pytest.mark.parametrize("invalid", ["duplicate", "unknown"])
def test_cli_archive_mappings_are_checked_before_any_input_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, invalid: str
) -> None:
    first, second, unknown = (tmp_path / name for name in ("off", "phenotype", "unknown"))
    archive = tmp_path / "archive.json"
    output = tmp_path / "uncreated" / "report"

    def forbid_read(*_args: Any, **_kwargs: Any) -> Any:
        pytest.fail("Invalid archive mappings must be refused before input reads")

    monkeypatch.setattr(scopes, "ArchiveResolver", forbid_read, raising=False)
    monkeypatch.setattr(scopes, "load_run", forbid_read)
    mapped = [first, first] if invalid == "duplicate" else [unknown]
    argv = ["--scan", str(first), "--scan", str(second), "--output-dir", str(output)]
    for directory in mapped:
        argv.extend(["--archive", str(directory), str(archive)])
    with pytest.raises(ValueError, match=rf"(?i){invalid}|not.*declared|not.*scan"):
        scopes.main([*argv, "--allow-temporary-output"])
    assert not output.parent.exists()


def test_cli_explicit_archives_pair_the_original_scan_directories(
    scans: dict[str, Path],
    loaded_scans: dict[str, scopes.ScanRun],
    retained: ArchiveResolver,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "report"
    argv = [
        "--output-dir",
        str(output),
        "--code-commit",
        "pair-fixture",
        "--allow-temporary-output",
    ]
    for scope in ("off", "phenotype"):
        argv.extend(
            [
                "--scan",
                str(scans[scope]),
                "--archive",
                str(scans[scope]),
                str(retained.root / "archive.json"),
            ]
        )
    ordered = scopes.pair_scopes(list(loaded_scans.values()))
    expected_table, expected_report = scopes.side_by_side(
        ordered, scan.load_protocol(), "pair-fixture"
    )
    forbid_original_reads(monkeypatch, retained)
    assert scopes.main(argv) == 0
    actual = pd.read_csv(
        output / scopes.OUTPUTS["table"], sep="\t", dtype=str, keep_default_na=False
    )
    pd.testing.assert_frame_equal(actual, expected_table)
    report = json.loads((output / scopes.OUTPUTS["summary"]).read_bytes())
    for key in ("protocol", "strain_rollup", "rule", "scopes", "notes"):
        assert report[key] == expected_report[key]
