"""The shipped influence probe must distinguish guard kills from invalid runs (#2078)."""

from __future__ import annotations

import importlib.util
import json
import sys
from importlib.machinery import ModuleSpec
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[1]
PROBE = ROOT / "scripts/adhoc/issue2078_probe_influence_guards.py"
TARGET = "tests/test_synthetic.py::test_guard"


def _load_probe() -> ModuleType:
    name = "issue2078_probe_test_subject"
    spec = importlib.util.spec_from_file_location(name, PROBE)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        del sys.modules[name]
    return module


def _case(probe: ModuleType) -> Any:
    return probe.Case(
        name="synthetic",
        module="mtp.influence_noel",
        source="src/mtp/influence_noel.py",
        old="value = 1\n",
        new="value = 2\n",
        test=TARGET,
    )


def _report(nodeid: str, when: str, outcome: str) -> Any:
    return SimpleNamespace(
        nodeid=nodeid,
        when=when,
        outcome=outcome,
        passed=outcome == "passed",
        failed=outcome == "failed",
    )


def _record_call(results: Any, error: BaseException | None = None, *, nodeid: str = TARGET) -> None:
    results.pytest_runtest_logreport(_report(nodeid, "setup", "passed"))
    report = _report(nodeid, "call", "passed" if error is None else "failed")
    results.pytest_runtest_logreport(report)
    if error is not None:
        results.pytest_exception_interact(
            SimpleNamespace(excinfo=SimpleNamespace(value=error)), report
        )
    results.pytest_runtest_logreport(_report(nodeid, "teardown", "passed"))


def test_each_shipped_mutation_has_one_current_anchor_and_one_test_file() -> None:
    probe = _load_probe()
    assert set(probe.CASES) == {f"M{number:02}" for number in range(1, 15)}
    for name, case in probe.CASES.items():
        assert case.name == name
        original = (ROOT / case.source).read_bytes()
        changed = probe.mutated_source(case, original, control=False)
        assert changed != original
        assert changed == original.replace(case.old.encode(), case.new.encode(), 1)
        assert probe.mutated_source(case, original, control=True) == original
        file, separator, node = case.test.partition("::")
        assert separator == "::"
        assert file in {
            "tests/test_influence_extension_status.py",
            "tests/test_influence_extension_noel_check.py",
            "tests/test_influence_extension_noel_core.py",
        }
        assert f"def {node.partition('[')[0]}(" in (ROOT / file).read_text()


@pytest.mark.parametrize("content", [b"value = 0\n", b"value = 1\nvalue = 1\n"])
@pytest.mark.parametrize("control", [False, True])
def test_missing_or_duplicate_anchor_refuses_even_for_control(
    content: bytes, control: bool
) -> None:
    probe = _load_probe()
    with pytest.raises(ValueError, match="exactly one mutation anchor"):
        probe.mutated_source(_case(probe), content, control=control)


@pytest.mark.parametrize("mismatch", ["file", "spec", "missing_file", "missing_spec"])
def test_module_file_and_loader_origin_must_both_be_canonical(
    tmp_path: Path, mismatch: str
) -> None:
    probe = _load_probe()
    canonical = tmp_path / "canonical.py"
    foreign = tmp_path / "foreign.py"
    module = ModuleType("synthetic")
    module.__file__ = str(canonical)
    module.__spec__ = ModuleSpec("synthetic", None, origin=str(canonical))
    probe.require_module_origin(module, canonical)
    if mismatch == "file":
        module.__file__ = str(foreign)
    elif mismatch == "spec":
        module.__spec__.origin = str(foreign)
    elif mismatch == "missing_file":
        del module.__file__
    else:
        module.__spec__ = None
    with pytest.raises(ValueError, match="module origin"):
        probe.require_module_origin(module, canonical)


@pytest.mark.parametrize(
    "error", [AssertionError("guard failed"), pytest.fail.Exception("DID NOT RAISE ValueError")]
)
def test_only_target_call_assertions_count_as_mutation_kills(error: BaseException) -> None:
    probe = _load_probe()
    results = probe.Results()
    _record_call(results, error)
    assert results.accepted(1, TARGET, control=False)
    assert not results.accepted(0, TARGET, control=False)
    assert not results.accepted(1, TARGET, control=True)


@pytest.mark.parametrize(
    "failure", ["collection", "setup", "teardown", "wrong_test", "raw_exception", "skipped"]
)
def test_unrelated_failures_are_never_credited_as_guard_kills(failure: str) -> None:
    probe = _load_probe()
    results = probe.Results()
    if failure == "wrong_test":
        _record_call(
            results, AssertionError("other test"), nodeid="tests/test_synthetic.py::test_other"
        )
    else:
        _record_call(
            results,
            ValueError("unrelated parser failure")
            if failure == "raw_exception"
            else AssertionError("intended guard"),
        )
        if failure == "collection":
            results.pytest_collectreport(_report("tests/test_missing.py", "collect", "failed"))
        elif failure in {"setup", "teardown"}:
            results.pytest_runtest_logreport(_report(TARGET, failure, "failed"))
        elif failure == "skipped":
            results.pytest_runtest_logreport(_report(TARGET + "[skipped]", "setup", "skipped"))
    assert not results.accepted(1, TARGET, control=False)


@pytest.mark.parametrize("exit_code", [0, 1, 2, 3, 4, 5, 125, 137])
def test_no_test_calls_never_verify(exit_code: int) -> None:
    probe = _load_probe()
    assert not probe.Results().accepted(exit_code, TARGET, control=False)
    assert not probe.Results().accepted(exit_code, TARGET, control=True)


def test_control_requires_all_target_calls_to_pass() -> None:
    probe = _load_probe()
    results = probe.Results()
    _record_call(results, nodeid=TARGET + "[first]")
    _record_call(results, nodeid=TARGET + "[second]")
    assert results.accepted(0, TARGET, control=True)
    assert not results.accepted(1, TARGET, control=True)
    assert not results.accepted(0, TARGET, control=False)


def test_snapshot_refuses_changed_source_even_when_size_is_unchanged(tmp_path: Path) -> None:
    probe = _load_probe()
    source = tmp_path / "source.py"
    source.write_text("value = 1\n")
    before = probe.snapshot(tmp_path, ("source.py",))
    source.write_text("value = 2\n")
    with pytest.raises(ValueError, match="Canonical input changed"):
        probe.require_snapshot(tmp_path, before)


def test_snapshot_refuses_identity_replacement_with_identical_bytes(tmp_path: Path) -> None:
    probe = _load_probe()
    source, replacement = tmp_path / "source.py", tmp_path / "replacement.py"
    source.write_text("value = 1\n")
    before = probe.snapshot(tmp_path, ("source.py",))
    replacement.write_bytes(source.read_bytes())
    replacement.replace(source)
    with pytest.raises(ValueError, match="Canonical input changed"):
        probe.require_snapshot(tmp_path, before)


def test_snapshot_refuses_symlinks_and_oversized_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    probe = _load_probe()
    source = tmp_path / "source.py"
    source.write_text("value = 1\n")
    (tmp_path / "alias.py").symlink_to(source)
    with pytest.raises(ValueError, match="bounded regular"):
        probe.snapshot(tmp_path, ("alias.py",))
    monkeypatch.setattr(probe, "MAX_FILE_BYTES", 4)
    with pytest.raises(ValueError, match="bounded regular"):
        probe.snapshot(tmp_path, ("source.py",))


def test_control_receipt_is_bound_to_case_and_current_canonical_inputs(tmp_path: Path) -> None:
    probe = _load_probe()
    source = tmp_path / "source.py"
    source.write_text("value = 1\n")
    before = probe.snapshot(tmp_path, ("source.py",))
    case = _case(probe)
    record = {
        "case": case.name,
        "mode": "control",
        "verified": True,
        "pytest_exit": 0,
        "target": case.test,
        "canonical_before": before,
        "canonical_after": before,
    }
    receipt = tmp_path / "evidence.json"
    receipt.write_text(json.dumps(record))
    probe.require_control_receipt(receipt, case, before)
    for field, changed in (
        ("case", "other"),
        ("verified", False),
        ("pytest_exit", 1),
        ("canonical_after", {}),
        ("canonical_before", {}),
        ("mode", "mutation"),
    ):
        receipt.write_text(json.dumps({**record, field: changed}))
        with pytest.raises(ValueError, match="matching passing control"):
            probe.require_control_receipt(receipt, case, before)
