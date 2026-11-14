"""Standalone core/arm verification binds the parsed manifest bytes (#1216)."""

from __future__ import annotations

import hashlib
import inspect
import re
from collections.abc import Callable
from pathlib import Path

import pytest

from mtp.datamodel.workflow_pydantic import ArmRun, WorkflowRun
from mtp.workflow import runtime, verification
from mtp.workflow.arms import runtime as arm_runtime
from test_comparison_manifest_bytes import _metadata_swap
from test_comparison_run import stage as stage


def _requires_pin(function: Callable[..., object]) -> None:
    """An absent keyword is an API assertion, never a claimed checksum refusal."""
    parameter = inspect.signature(function).parameters.get("expected_sha256")
    assert parameter is not None, "Missing expected_sha256 API (structural failure)"
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is None


def test_standalone_manifest_pin_apis_and_shared_core_digest_are_explicit() -> None:
    for function in (runtime.validate_run, verification.verify_core_run, arm_runtime.verify_arm):
        _requires_pin(function)
    parameter = inspect.signature(arm_runtime._verify_arm_with_core).parameters.get(
        "core_manifest_sha256"
    )
    assert parameter is not None, "Shared-core traversal must carry its verified digest"
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is inspect.Parameter.empty


@pytest.mark.parametrize("entrypoint", ["validate", "verify"])
def test_standalone_core_parses_the_byte_snapshot_not_a_second_text_read(
    stage: dict[str, Path], monkeypatch: pytest.MonkeyPatch, entrypoint: str
) -> None:
    target = stage["core"].resolve()
    original = target.read_bytes()
    changed = _metadata_swap(original)
    expected = WorkflowRun.model_validate_json(original)
    read_text, read_bytes = Path.read_text, Path.read_bytes

    def swapped_text(path: Path, encoding: str | None = None, errors: str | None = None) -> str:
        if path == target:
            return changed.decode("utf-8")
        return read_text(path, encoding=encoding, errors=errors)

    monkeypatch.setattr(Path, "read_text", swapped_text)
    actual = (
        runtime.validate_run(target)
        if entrypoint == "validate"
        else verification.verify_core_run(target)
    )
    assert actual.started_at == expected.started_at
    assert actual.git_dirty == expected.git_dirty
    assert read_bytes(target) == original


def test_standalone_arm_parses_the_byte_snapshot_not_a_second_text_read(
    stage: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    target = stage["noel"].resolve()
    original = target.read_bytes()
    changed = _metadata_swap(original)
    expected = ArmRun.model_validate_json(original)
    read_text, read_bytes = Path.read_text, Path.read_bytes

    def swapped_text(path: Path, encoding: str | None = None, errors: str | None = None) -> str:
        if path == target:
            return changed.decode("utf-8")
        return read_text(path, encoding=encoding, errors=errors)

    monkeypatch.setattr(Path, "read_text", swapped_text)
    actual = arm_runtime.verify_arm(target)
    assert actual.started_at == expected.started_at
    assert actual.git_dirty == expected.git_dirty
    assert read_bytes(target) == original


@pytest.mark.parametrize("role", ["core", "noel"])
def test_standalone_expected_pin_refuses_transient_first_bytes(
    stage: dict[str, Path], monkeypatch: pytest.MonkeyPatch, role: str
) -> None:
    verifier = verification.verify_core_run if role == "core" else arm_runtime.verify_arm
    _requires_pin(verifier)
    target = stage[role].resolve()
    original = target.read_bytes()
    changed = _metadata_swap(original)
    digest = hashlib.sha256(original).hexdigest()
    read_bytes = Path.read_bytes
    calls = 0

    def swapped_first_bytes(path: Path) -> bytes:
        nonlocal calls
        if path == target:
            calls += 1
            return changed if calls == 1 else original
        return read_bytes(path)

    monkeypatch.setattr(Path, "read_bytes", swapped_first_bytes)
    with pytest.raises(ValueError, match=r"(?i)(sha-?256|digest|hash)"):
        verifier(target, expected_sha256=digest)
    assert calls >= 1
    assert read_bytes(target) == original


def test_standalone_arm_checks_core_first_bytes_against_its_stored_pin(
    stage: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    target = stage["core"].resolve()
    original = target.read_bytes()
    changed = _metadata_swap(original)
    arm_bytes = stage["noel"].read_bytes()
    arm = ArmRun.model_validate_json(arm_bytes)
    assert arm.core_manifest_sha256 == hashlib.sha256(original).hexdigest()
    read_bytes = Path.read_bytes
    calls = 0

    def swapped_first_bytes(path: Path) -> bytes:
        nonlocal calls
        if path == target:
            calls += 1
            return changed if calls == 1 else original
        return read_bytes(path)

    monkeypatch.setattr(Path, "read_bytes", swapped_first_bytes)
    with pytest.raises(ValueError, match=r"(?i)(sha-?256|digest|hash)"):
        arm_runtime.verify_arm(stage["noel"])
    assert calls >= 1
    assert read_bytes(target) == original
    assert read_bytes(stage["noel"]) == arm_bytes


def test_shared_core_digest_contradiction_refused_even_with_matching_model(
    stage: dict[str, Path],
) -> None:
    parameter = inspect.signature(arm_runtime._verify_arm_with_core).parameters.get(
        "core_manifest_sha256"
    )
    assert parameter is not None, "Missing shared-core digest API (structural failure)"
    original = stage["core"].read_bytes()
    core = verification.verify_core_run(stage["core"])
    assert core.run_id == ArmRun.model_validate_json(stage["noel"].read_bytes()).core_run_id
    contradictory = hashlib.sha256(_metadata_swap(original)).hexdigest()
    assert contradictory != hashlib.sha256(original).hexdigest()
    with pytest.raises(ValueError, match=r"(?i)(sha-?256|digest|hash)"):
        arm_runtime._verify_arm_with_core(stage["noel"], core, core_manifest_sha256=contradictory)


def test_standalone_exact_pins_preserve_real_core_and_arm_records(
    stage: dict[str, Path],
) -> None:
    for function in (runtime.validate_run, verification.verify_core_run, arm_runtime.verify_arm):
        _requires_pin(function)
    core_raw, arm_raw = stage["core"].read_bytes(), stage["noel"].read_bytes()
    core_digest, arm_digest = (
        hashlib.sha256(core_raw).hexdigest(),
        hashlib.sha256(arm_raw).hexdigest(),
    )
    core = verification.verify_core_run(stage["core"], expected_sha256=core_digest)
    assert runtime.validate_run(stage["core"], expected_sha256=core_digest) == core
    assert core == WorkflowRun.model_validate_json(core_raw)
    arm = arm_runtime.verify_arm(stage["noel"], expected_sha256=arm_digest)
    assert arm == ArmRun.model_validate_json(arm_raw)
    assert re.fullmatch(r"[0-9a-f]{64}", arm.core_manifest_sha256)
    assert stage["core"].read_bytes() == core_raw
    assert stage["noel"].read_bytes() == arm_raw
