"""Comparison metadata must come from the exact manifest bytes it verifies (#1216)."""

from __future__ import annotations

import hashlib
import inspect
import json
import re
from pathlib import Path

import pytest

from mtp.datamodel.workflow_pydantic import ArmRun, ComparisonRun, WorkflowRun
from mtp.multiplicity import sources as multiplicity_sources
from mtp.workflow import comparison
from test_comparison_run import stage as stage


def _metadata_swap(raw: bytes) -> bytes:
    """Change only two top-level metadata fields without changing the byte width."""
    original = json.loads(raw)
    started = original["started_at"]
    dirty = original["git_dirty"]
    assert isinstance(started, str)
    assert isinstance(dirty, bool)
    year = "2000" if not started.startswith("2000") else "1996"
    changed_started = year + started[4:]
    text = raw.decode("utf-8")
    started_pattern = rf'("started_at"\s*:\s*)"{re.escape(started)}"'
    text, count = re.compile(started_pattern).subn(
        lambda match: f'{match.group(1)}"{changed_started}"', text
    )
    assert count == 1
    dirty_pattern = rf'("git_dirty"\s*:\s*){str(dirty).lower()}'
    text, count = re.compile(dirty_pattern).subn(
        lambda match: match.group(1) + str(not dirty).lower(), text
    )
    assert count == 1
    changed = text.encode("utf-8")
    # Core JSON is compact; arm/comparison JSON is indented. Their writers all
    # append LF, which can absorb the one-byte true/false width difference.
    if len(changed) == len(raw) + 1:
        assert changed.endswith(b"\n")
        changed = changed[:-1]
    elif len(changed) == len(raw) - 1:
        changed += b" "
    assert len(changed) == len(raw)
    assert changed != raw
    modified = json.loads(changed)
    assert set(modified) == set(original)
    assert {name for name in original if original[name] != modified[name]} == {
        "started_at",
        "git_dirty",
    }
    return changed


def _captured_parent(
    stage: dict[str, Path], role: str, tmp_path: Path
) -> tuple[Path, bytes, bytes, WorkflowRun | ArmRun, bytes]:
    target = stage[role].resolve()
    original = target.read_bytes()
    changed = _metadata_swap(original)
    model = WorkflowRun if role == "core" else ArmRun
    before = model.model_validate_json(original)
    after = model.model_validate_json(changed)
    assert before.started_at != after.started_at
    assert before.git_dirty != after.git_dirty
    comparison_raw = (stage["comparison"] / "comparison_run.json").read_bytes()
    (tmp_path / f"{role}.A.json").write_bytes(original)
    (tmp_path / f"{role}.B.json").write_bytes(changed)
    (tmp_path / "comparison.A.json").write_bytes(comparison_raw)
    return target, original, changed, before, comparison_raw


@pytest.mark.parametrize("role", ["core", "noel", "jed", "ilyas"])
def test_comparison_parent_metadata_comes_from_hashed_bytes_not_text_reread(
    stage: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str
) -> None:
    """The old text parser sees B while every file/hash read still sees valid A."""
    target, original, changed, expected, comparison_raw = _captured_parent(stage, role, tmp_path)
    manifest = stage["comparison"] / "comparison_run.json"
    read_text = Path.read_text
    read_bytes = Path.read_bytes

    def swapped_text(path: Path, encoding: str | None = None, errors: str | None = None) -> str:
        if path == target:
            return changed.decode("utf-8")
        return read_text(path, encoding=encoding, errors=errors)

    monkeypatch.setattr(Path, "read_text", swapped_text)
    record, core, arms = comparison.verify_comparison_sources(
        manifest, manifest_bytes=comparison_raw
    )
    actual = core if role == "core" else arms[role]
    assert actual.started_at == expected.started_at
    assert actual.git_dirty == expected.git_dirty
    assert record == ComparisonRun.model_validate_json(comparison_raw)
    assert read_bytes(target) == original
    assert read_bytes(manifest) == comparison_raw


@pytest.mark.parametrize("role", ["core", "noel", "jed", "ilyas"])
def test_comparison_parent_second_byte_read_cannot_replace_verified_metadata(
    stage: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str
) -> None:
    """A future byte parser must reuse A or refuse B; a second read cannot supply B."""
    target, original, changed, expected, comparison_raw = _captured_parent(stage, role, tmp_path)
    manifest = stage["comparison"] / "comparison_run.json"
    read_bytes = Path.read_bytes
    target_reads = 0

    def swapped_bytes(path: Path) -> bytes:
        nonlocal target_reads
        if path == target:
            target_reads += 1
            return changed if target_reads == 2 else original
        return read_bytes(path)

    monkeypatch.setattr(Path, "read_bytes", swapped_bytes)
    refusal: str | None = None
    actual: WorkflowRun | ArmRun | None = None
    try:
        _, core, arms = comparison.verify_comparison_sources(
            manifest, manifest_bytes=comparison_raw
        )
        actual = core if role == "core" else arms[role]
    except ValueError as error:
        refusal = str(error)
    if refusal is not None:
        # Only an explicit identity refusal after reaching B is an acceptable defense.
        assert target_reads >= 2
        assert re.search(r"(?i)(sha-?256|digest|hash)", refusal), refusal
    else:
        assert actual is not None
        assert actual.started_at == expected.started_at
        assert actual.git_dirty == expected.git_dirty
    assert read_bytes(target) == original
    assert read_bytes(manifest) == comparison_raw


@pytest.mark.parametrize("role", ["core", "noel", "jed", "ilyas"])
def test_comparison_rejects_transient_parent_bytes_against_recorded_pin(
    stage: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str
) -> None:
    """A parent parsed from B must fail its comparison pin even when disk hashes see A."""
    target, original, changed, _, comparison_raw = _captured_parent(stage, role, tmp_path)
    manifest = stage["comparison"] / "comparison_run.json"
    read_bytes = Path.read_bytes
    target_reads = 0

    def swapped_first_bytes(path: Path) -> bytes:
        nonlocal target_reads
        if path == target:
            target_reads += 1
            return changed if target_reads == 1 else original
        return read_bytes(path)

    monkeypatch.setattr(Path, "read_bytes", swapped_first_bytes)
    with pytest.raises(ValueError, match=r"(?i)(sha-?256|digest|hash)"):
        comparison.verify_comparison_sources(manifest, manifest_bytes=comparison_raw)
    assert target_reads >= 1
    assert read_bytes(target) == original
    assert read_bytes(manifest) == comparison_raw


def _require_comparison_pin_api() -> None:
    """An absent API is a structural failure, not evidence that a hash guard worked."""
    parameters = inspect.signature(comparison.verify_comparison_sources).parameters
    assert "expected_sha256" in parameters, "Comparison verifier has no expected-digest API"
    assert parameters["expected_sha256"].kind is inspect.Parameter.KEYWORD_ONLY
    assert parameters["expected_sha256"].default is None


@pytest.mark.parametrize("buffered", [False, True], ids=["path", "buffer"])
def test_comparison_expected_hash_accepts_the_exact_parsed_bytes(
    stage: dict[str, Path], buffered: bool
) -> None:
    _require_comparison_pin_api()
    manifest = stage["comparison"] / "comparison_run.json"
    raw = manifest.read_bytes()
    expected = ComparisonRun.model_validate_json(raw)
    record, core, arms = comparison.verify_comparison_sources(
        manifest,
        manifest_bytes=raw if buffered else None,
        expected_sha256=hashlib.sha256(raw).hexdigest(),
    )
    assert record == expected
    assert core.run_id == expected.core_run_id
    assert set(arms) == {"noel", "jed", "ilyas"}
    assert manifest.read_bytes() == raw


@pytest.mark.parametrize(
    ("case", "buffered"),
    [
        ("changed-buffer", True),
        ("wrong-pin", False),
        ("wrong-pin", True),
        ("empty-pin", False),
        ("empty-pin", True),
        ("empty-buffer", True),
    ],
)
def test_comparison_expected_hash_refuses_before_upstream_traversal(
    stage: dict[str, Path], monkeypatch: pytest.MonkeyPatch, case: str, buffered: bool
) -> None:
    _require_comparison_pin_api()
    manifest = stage["comparison"] / "comparison_run.json"
    raw = manifest.read_bytes()
    supplied = raw if buffered else None
    digest = hashlib.sha256(raw).hexdigest()
    if case == "changed-buffer":
        supplied = _metadata_swap(raw)
        changed = ComparisonRun.model_validate_json(supplied)
        assert changed.started_at != ComparisonRun.model_validate_json(raw).started_at
    elif case == "wrong-pin":
        digest = hashlib.sha256(raw + b"a different expected snapshot").hexdigest()
    elif case == "empty-pin":
        digest = ""
    else:
        assert case == "empty-buffer"
        supplied = b""
    calls: list[tuple[object, ...]] = []

    def unexpected_core(*args: object, **kwargs: object) -> WorkflowRun:
        calls.append((*args, kwargs))
        raise AssertionError("Mismatched manifest bytes reached core verification")

    monkeypatch.setattr(comparison, "verify_core_run", unexpected_core)
    with pytest.raises(ValueError, match=r"(?i)(sha-?256|digest|hash)"):
        comparison.verify_comparison_sources(
            manifest, manifest_bytes=supplied, expected_sha256=digest
        )
    assert calls == []
    assert manifest.read_bytes() == raw


def test_multiplicity_forwards_both_optional_comparison_bindings_explicitly(
    stage: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """The adapter must forward both keywords even when both are omitted by its caller."""
    parameters = inspect.signature(multiplicity_sources.open_verified_sources).parameters
    assert "expected_sha256" in parameters, "Multiplicity adapter has no expected-digest API"
    assert parameters["expected_sha256"].kind is inspect.Parameter.KEYWORD_ONLY
    assert parameters["expected_sha256"].default is None
    manifest = stage["comparison"] / "comparison_run.json"
    raw = manifest.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    verified = comparison.verify_comparison_sources(manifest)
    calls: list[tuple[Path, bytes | None, str | None]] = []

    def strict_verifier(
        path: Path, **keywords: object
    ) -> tuple[ComparisonRun, WorkflowRun, dict[str, ArmRun]]:
        assert set(keywords) == {"manifest_bytes", "expected_sha256"}, (
            "Multiplicity must explicitly forward both comparison binding keywords"
        )
        manifest_bytes = keywords["manifest_bytes"]
        expected_sha256 = keywords["expected_sha256"]
        assert manifest_bytes is None or isinstance(manifest_bytes, bytes)
        assert expected_sha256 is None or isinstance(expected_sha256, str)
        calls.append((path, manifest_bytes, expected_sha256))
        return verified

    monkeypatch.setattr(multiplicity_sources, "verify_comparison_sources", strict_verifier)
    variants = [(None, None), (raw, None), (None, digest), (raw, digest)]
    for captured, pin in variants:
        if captured is None and pin is None:
            loaded = multiplicity_sources.open_verified_sources(manifest)
        else:
            loaded = multiplicity_sources.open_verified_sources(
                manifest, manifest_bytes=captured, expected_sha256=pin
            )
        assert loaded.comparison == verified[0]
        assert loaded.core == verified[1]
        assert loaded.arms == verified[2]
    assert calls == [(manifest, captured, pin) for captured, pin in variants]
