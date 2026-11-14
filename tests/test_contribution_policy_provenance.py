"""Policy provenance records the parsed metadata bytes without claiming core verification."""

from __future__ import annotations

import hashlib
import json
from dataclasses import FrozenInstanceError, is_dataclass
from pathlib import Path
from typing import Any

import pytest

from mtp import contribution_inputs
from mtp.datamodel.workflow_pydantic import FileInventoryEntry, WorkflowRun
from mtp.trait_key_policy import DEFAULT_TRAIT_KEY_POLICY, TraitKeyPolicy
from mtp.workflow.trait_policy import HISTORICAL_MATRIX_POLICIES, MATRIX_SOURCE


def _resolve(
    inputs: Path, requested: str | None = None, *, core_manifest: Path | None = None
) -> Any:
    resolver = getattr(contribution_inputs, "resolve_inputs_trait_key_policy_record", None)
    assert callable(resolver), "Policy resolution must expose its consumed metadata provenance"
    return resolver(inputs, requested, core_manifest=core_manifest)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _manifest(path: Path, policy: TraitKeyPolicy = "raw_object_v1") -> Path:
    fixture = Path(__file__).parent / "fixtures/trait_key_policy/historical_raw_core/core/run.json"
    run = WorkflowRun.model_validate_json(fixture.read_bytes())
    identity = next(key for key, value in HISTORICAL_MATRIX_POLICIES.items() if value == policy)
    run.code_inventory = [
        FileInventoryEntry(path=MATRIX_SOURCE, sha256=identity[0], size_bytes=identity[1])
    ]
    path.write_text(json.dumps(run.model_dump(mode="json")) + "\n", encoding="utf-8")
    return path


def _columns(path: Path, policy: str | None = "raw_object_v1") -> Path:
    text = (
        f"column_id\ttrait_key_policy\nexample\t{policy}\n"
        if policy is not None
        else "column_id\nexample\n"
    )
    path.write_text(text, encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("manifest_present", "column_kind", "requested", "basis", "policy"),
    [
        (False, "absent", None, "default", DEFAULT_TRAIT_KEY_POLICY),
        (False, "absent", "raw_object_v1", "explicit", "raw_object_v1"),
        (False, "no_policy", None, "default", DEFAULT_TRAIT_KEY_POLICY),
        (False, "policy", None, "column_metadata", "raw_object_v1"),
        (True, "absent", None, "core_manifest", "raw_object_v1"),
        (True, "policy", None, "core_manifest", "raw_object_v1"),
        (True, "no_policy", None, "core_manifest", "raw_object_v1"),
        (True, "policy", "raw_object_v1", "explicit", "raw_object_v1"),
        (False, "policy", "raw_object_v1", "explicit", "raw_object_v1"),
    ],
)
def test_policy_resolution_records_basis_and_every_consumed_source(
    tmp_path: Path,
    manifest_present: bool,
    column_kind: str,
    requested: str | None,
    basis: str,
    policy: str,
) -> None:
    sources = {}
    if manifest_present:
        sources["core_manifest"] = _manifest(tmp_path / "run.json")
    if column_kind != "absent":
        sources["taxa_trait_matrix_columns.tsv"] = _columns(
            tmp_path / "taxa_trait_matrix_columns.tsv",
            None if column_kind == "no_policy" else "raw_object_v1",
        )
    expected_hashes = {key: _sha256(path.read_bytes()) for key, path in sources.items()}

    record = _resolve(tmp_path, requested)

    assert record.policy == policy
    assert record.basis == basis
    assert record.source_sha256 == expected_hashes
    assert record.as_dict() == {
        "policy": policy,
        "basis": basis,
        "source_sha256": expected_hashes,
    }
    assert json.loads(json.dumps(record.as_dict())) == record.as_dict()
    assert str(tmp_path) not in json.dumps(record.as_dict())
    assert contribution_inputs.resolve_inputs_trait_key_policy(tmp_path, requested) == policy


def test_unrelated_manifest_is_parsed_metadata_not_core_verification(tmp_path: Path) -> None:
    inputs = tmp_path / "empty-inputs"
    inputs.mkdir()
    unrelated = _manifest(tmp_path / "unrelated-historical-core.json")
    before = unrelated.read_bytes()

    assert (
        contribution_inputs.resolve_inputs_trait_key_policy(inputs, core_manifest=unrelated)
        == "raw_object_v1"
    )
    record = _resolve(inputs, core_manifest=unrelated)

    assert record.policy == "raw_object_v1"
    assert record.basis == "core_manifest"
    assert record.source_sha256 == {"core_manifest": _sha256(before)}
    assert list(inputs.iterdir()) == []
    assert unrelated.read_bytes() == before


@pytest.mark.parametrize("source_kind", ["manifest", "columns", "both"])
def test_policy_source_hashes_bind_the_same_buffers_that_were_parsed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source_kind: str
) -> None:
    sources = {}
    changed = {}
    if source_kind in {"manifest", "both"}:
        path = _manifest(tmp_path / "run.json")
        sources["core_manifest"] = path
        changed[path] = _manifest(tmp_path / "changed.json", "facultative_union_v1").read_bytes()
    if source_kind in {"columns", "both"}:
        path = _columns(tmp_path / "taxa_trait_matrix_columns.tsv")
        sources["taxa_trait_matrix_columns.tsv"] = path
        changed[path] = b"column_id\ttrait_key_policy\nexample\tfacultative_union_v1\n"
    original = {path: path.read_bytes() for path in sources.values()}
    reads = dict.fromkeys(original, 0)
    original_read_bytes = Path.read_bytes
    original_open = Path.open

    def racing_read_bytes(path: Path) -> bytes:
        if path not in original:
            return original_read_bytes(path)
        reads[path] += 1
        return original[path] if reads[path] == 1 else changed[path]

    def forbid_independent_parse(path: Path, *args: Any, **kwargs: Any) -> Any:
        assert path not in original, "Metadata must be parsed from its single captured bytes buffer"
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_bytes", racing_read_bytes)
    monkeypatch.setattr(Path, "open", forbid_independent_parse)

    record = _resolve(tmp_path, "raw_object_v1")

    assert record.policy == "raw_object_v1"
    assert record.basis == "explicit"
    assert record.source_sha256 == {key: _sha256(original[path]) for key, path in sources.items()}
    assert reads == dict.fromkeys(original, 1)


def test_policy_resolution_record_and_export_are_independent_immutable_values(
    tmp_path: Path,
) -> None:
    _columns(tmp_path / "taxa_trait_matrix_columns.tsv")
    record = _resolve(tmp_path)
    expected = record.as_dict()
    assert is_dataclass(type(record))
    with pytest.raises(FrozenInstanceError):
        record.basis = "default"
    with pytest.raises(TypeError):
        record.source_sha256["forged"] = "not-a-hash"
    exported = record.as_dict()
    exported["source_sha256"]["forged"] = "not-a-hash"
    assert record.as_dict() == expected


@pytest.mark.parametrize("contradiction", ["explicit", "manifest_columns"])
def test_policy_provenance_still_refuses_conflicting_declarations(
    tmp_path: Path, contradiction: str
) -> None:
    _manifest(tmp_path / "run.json")
    _columns(
        tmp_path / "taxa_trait_matrix_columns.tsv",
        "facultative_union_v1" if contradiction == "manifest_columns" else "raw_object_v1",
    )
    requested = "facultative_union_v1" if contradiction == "explicit" else None
    message = ""
    try:
        _resolve(tmp_path, requested)
    except ValueError as error:
        message = str(error)
    assert "contradicts" in message


@pytest.mark.parametrize(
    ("data", "expected"),
    [
        (b"trait_key_policy\ttrait_key_policy\nraw_object_v1\traw_object_v1\n", "duplicate"),
        (b"column_id\ttrait_key_policy\nexample\traw_object_v1\textra\n", "Wrong field count"),
        (b"", "Empty"),
    ],
)
def test_buffered_column_metadata_preserves_strict_table_validation(
    tmp_path: Path, data: bytes, expected: str
) -> None:
    (tmp_path / "taxa_trait_matrix_columns.tsv").write_bytes(data)
    message = ""
    try:
        _resolve(tmp_path)
    except ValueError as error:
        message = str(error)
    assert expected in message
