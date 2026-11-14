"""Regression tests for isolated ad-hoc mutation probes (#966, #1184)."""

from __future__ import annotations

import ast
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

import mtp
from mtp.paths import file_sha256

REPO_ROOT = Path(__file__).resolve().parents[1]
PROBE_959 = REPO_ROOT / "scripts/adhoc/issue959_probe_coverage_mutants.py"
PROBE_1182 = REPO_ROOT / "scripts/adhoc/issue1182_probe_scan_bindings.py"
PROBE_1170 = REPO_ROOT / "scripts/adhoc/issue1170_quadrant_stool_mutant.py"


def _load_probe_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _forbidden_canonical_writes(monkeypatch: pytest.MonkeyPatch, protected: set[Path]) -> None:
    resolved_protected = {path.resolve() for path in protected}
    orig_write_text = Path.write_text
    orig_write_bytes = Path.write_bytes
    orig_open = Path.open

    def guarded_write_text(self: Path, *args: Any, **kwargs: Any) -> int:
        if self.resolve() in resolved_protected:
            raise PermissionError(f"Forbidden write_text to canonical file: {self}")
        return orig_write_text(self, *args, **kwargs)

    def guarded_write_bytes(self: Path, *args: Any, **kwargs: Any) -> int:
        if self.resolve() in resolved_protected:
            raise PermissionError(f"Forbidden write_bytes to canonical file: {self}")
        return orig_write_bytes(self, *args, **kwargs)

    def guarded_open(self: Path, mode: str = "r", *args: Any, **kwargs: Any) -> Any:
        if (
            any(flag in mode for flag in ("w", "a", "+", "x"))
            and self.resolve() in resolved_protected
        ):
            raise PermissionError(f"Forbidden writable open({mode!r}) on canonical file: {self}")
        return orig_open(self, mode, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", guarded_write_text)
    monkeypatch.setattr(Path, "write_bytes", guarded_write_bytes)
    monkeypatch.setattr(Path, "open", guarded_open)


def test_adhoc_mutation_probes_define_isolation_and_checkout_guards() -> None:
    """#966/#1184: probes must define verify_checkout_identity and never write SOURCE in place."""
    for script_path in (PROBE_959, PROBE_1182, PROBE_1170):
        tree = ast.parse(script_path.read_text(encoding="utf-8"), filename=str(script_path))
        func_names = {node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)}
        assert "verify_checkout_identity" in func_names, (
            f"{script_path.name} must define verify_checkout_identity()"
        )
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr in {"write_text", "write_bytes"}
            ):
                receiver = ast.unparse(node.func.value)
                assert receiver not in {"SOURCE", "(ROOT / SOURCE)", "ROOT / SOURCE"}, (
                    f"{script_path.name} still writes directly to {receiver}"
                )


def test_issue959_probe_rejects_foreign_import_origin_before_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#966: a mismatched mtp import path or REPO_ROOT must fail before any file write."""
    probe: Any = _load_probe_module(PROBE_959, "issue959_probe_foreign")
    out_dir = tmp_path / "should_not_be_created"
    monkeypatch.setattr(sys, "argv", [str(PROBE_959), "--output-dir", str(out_dir)])
    monkeypatch.setattr(probe, "REPO_ROOT", tmp_path / "foreign_checkout")
    with pytest.raises(AssertionError, match="different checkout"):
        probe.main()
    assert not out_dir.exists()

    monkeypatch.undo()
    probe2: Any = _load_probe_module(PROBE_959, "issue959_probe_foreign_mtp")
    monkeypatch.setattr(sys, "argv", [str(PROBE_959), "--output-dir", str(out_dir)])
    monkeypatch.setattr(mtp, "__file__", str(tmp_path / "foreign_checkout/src/mtp/__init__.py"))
    with pytest.raises(AssertionError, match="different checkout"):
        probe2.main()
    assert not out_dir.exists()


def test_issue959_probe_is_verified_case_and_failure_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#966: issue959 probe rejects collection/setup errors and records hashes on failure."""
    probe: Any = _load_probe_module(PROBE_959, "issue959_probe_unit")
    source_before_sha = file_sha256(probe.SOURCE)
    test_before_sha = file_sha256(probe.TEST)

    _forbidden_canonical_writes(monkeypatch, {probe.SOURCE, probe.TEST})

    # 1. Verify collection errors or unrelated/setup failures do not count as mutant rejection
    assert (
        probe.is_verified_case(
            exit_code=1,
            expected="Missing segment callout:",
            log="AssertionError: Missing segment callout:",
            failures=[{"phase": "call", "failure": "AssertionError: Missing segment callout:"}],
            collection_errors=["SyntaxError"],
        )
        is False
    )
    assert (
        probe.is_verified_case(
            exit_code=1,
            expected="Missing segment callout:",
            log="AssertionError: unrelated failure",
            failures=[{"phase": "call", "failure": "AssertionError: unrelated failure"}],
            collection_errors=[],
        )
        is False
    )
    assert (
        probe.is_verified_case(
            exit_code=1,
            expected="Missing segment callout:",
            log="AssertionError: Missing segment callout: Mixed",
            failures=[
                {
                    "phase": "setup",
                    "failure": "AssertionError: Missing segment callout: Mixed",
                }
            ],
            collection_errors=[],
        )
        is False
    )

    # 2. Verify early failure path still records before/after hashes in evidence.json
    fail_out = tmp_path / "early_failure_out"
    monkeypatch.setattr(
        probe,
        "MUTANTS",
        (
            (
                "missing_target_mutant",
                "nonexistent_target_expression",
                "replacement",
                "Expected marker",
            ),
        ),
    )
    monkeypatch.setattr(sys, "argv", [str(PROBE_959), "--output-dir", str(fail_out)])
    with pytest.raises(ValueError, match="Expected exactly one mutation target"):
        probe.main()
    fail_evidence = json.loads((fail_out / "evidence.json").read_text(encoding="utf-8"))
    assert fail_evidence["source_sha256_before"] == source_before_sha
    assert fail_evidence["source_sha256_after"] == source_before_sha
    assert fail_evidence["test_sha256_before"] == test_before_sha
    assert fail_evidence["test_sha256_after"] == test_before_sha
    assert fail_evidence["canonical_unchanged"] is True


def test_issue959_probe_runs_isolated_mutations_and_records_failure_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#966: issue959 probe executes baseline and all 4 mutants in isolated copies."""
    probe_full: Any = _load_probe_module(PROBE_959, "issue959_probe_full")
    source_before_stat = probe_full.SOURCE.stat()
    test_before_stat = probe_full.TEST.stat()
    source_before_sha = file_sha256(probe_full.SOURCE)
    test_before_sha = file_sha256(probe_full.TEST)

    _forbidden_canonical_writes(monkeypatch, {probe_full.SOURCE, probe_full.TEST})
    full_out = tmp_path / "full_out"
    monkeypatch.setattr(sys, "argv", [str(PROBE_959), "--output-dir", str(full_out)])
    assert probe_full.main() == 0

    evidence = json.loads((full_out / "evidence.json").read_text(encoding="utf-8"))
    assert evidence["source_sha256_before"] == source_before_sha
    assert evidence["source_sha256_after"] == source_before_sha
    assert evidence["test_sha256_before"] == test_before_sha
    assert evidence["test_sha256_after"] == test_before_sha
    assert evidence["canonical_unchanged"] is True
    assert [row["case"] for row in evidence["results"]] == [
        "baseline",
        "dropped_label",
        "callout_on_bar",
        "no_leader",
        "scope_hardcoded_at_entrypoint",
    ]
    for row in evidence["results"]:
        assert row["verified"] is True
        origin = Path(row["loaded_module_origin"]).resolve()
        assert origin.is_relative_to(full_out.resolve())
        assert origin != probe_full.SOURCE.resolve()
        if row["case"] == "baseline":
            assert row["isolated_source_sha256"] == source_before_sha
        else:
            assert row["isolated_source_sha256"] != source_before_sha

    assert probe_full.SOURCE.stat().st_mtime_ns == source_before_stat.st_mtime_ns
    assert probe_full.TEST.stat().st_mtime_ns == test_before_stat.st_mtime_ns


def test_issue1182_probe_isolates_mutations_without_touching_canonical_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#1184: issue1182 probe stages isolated module copies and never writes canonical source."""
    probe: Any = _load_probe_module(PROBE_1182, "issue1182_probe_live")
    canonical_source = (probe.ROOT / probe.SOURCE).resolve()
    before_sha = file_sha256(canonical_source)
    before_mtime = canonical_source.stat().st_mtime_ns

    _forbidden_canonical_writes(monkeypatch, {canonical_source, probe.ROOT / probe.TEST_FILE})

    # Verify foreign checkout rejection
    monkeypatch.setattr(mtp, "__file__", str(tmp_path / "foreign/src/mtp/__init__.py"))
    with pytest.raises(AssertionError, match="different checkout"):
        probe.verify_checkout_identity()
    monkeypatch.setattr(mtp, "__file__", str(REPO_ROOT / "src/mtp/__init__.py"))

    # Verify load_isolated_module swaps and restores mtp.cellulose_decr.selection_scan_scopes
    from mtp.cellulose_decr import selection_scan_scopes as canonical_scopes

    isolated_copy = tmp_path / "isolated_selection_scan_scopes.py"
    isolated_copy.write_bytes(canonical_source.read_bytes())
    with probe.load_isolated_module(isolated_copy) as loaded:
        assert Path(loaded.__file__).resolve() == isolated_copy.resolve()
        assert sys.modules["mtp.cellulose_decr.selection_scan_scopes"] is loaded
    assert sys.modules["mtp.cellulose_decr.selection_scan_scopes"] is canonical_scopes
    assert Path(canonical_scopes.__file__).resolve() == canonical_source

    # Verify main() stages isolated copies for baseline + all 4 mutations
    out_dir = tmp_path / "scan_probe_out"

    class _InProcessFakeProcess:
        def __init__(self, target: Any, args: tuple[Any, ...]) -> None:
            self._target = target
            self._args = args
            self.exitcode: int | None = None

        def start(self) -> None:
            directory, isolated_source = self._args
            assert Path(isolated_source).resolve().is_relative_to(out_dir.resolve())
            assert Path(isolated_source).resolve() != canonical_source
            with probe.load_isolated_module(Path(isolated_source)) as mod:
                origin = str(Path(mod.__file__).resolve())
            case = directory.name
            if case == "baseline":
                payload = {
                    "entrypoint": "pytest.main",
                    "arguments": [],
                    "exit_code": 0,
                    "loaded_module_origin": origin,
                    "counts": {"passed": len(probe.TESTS)},
                    "failures": [],
                }
                self.exitcode = 0
            else:
                expected = sorted(probe.EXPECTED_FAILURES[case])
                msg = (
                    "A relocated scan read an original path: run.json"
                    if case == "bypass_bound_core_resolver"
                    else "DID NOT RAISE ValueError"
                )
                payload = {
                    "entrypoint": "pytest.main",
                    "arguments": [],
                    "exit_code": 1,
                    "loaded_module_origin": origin,
                    "counts": {
                        "passed": len(probe.TESTS) - len(expected),
                        "failed": len(expected),
                    },
                    "failures": [
                        {
                            "nodeid": f"{probe.TEST_FILE}::{name}",
                            "phase": "call",
                            "failure": msg,
                        }
                        for name in expected
                    ],
                }
                self.exitcode = 1
            (directory / "pytest.json").write_text(json.dumps(payload), encoding="utf-8")

        def join(self, timeout: float | None = None) -> None:
            _ = timeout
            return None

        def is_alive(self) -> bool:
            return False

    class _FakeContext:
        def Process(self, target: Any, args: tuple[Any, ...]) -> _InProcessFakeProcess:
            return _InProcessFakeProcess(target=target, args=args)

    monkeypatch.setattr(probe.multiprocessing, "get_context", lambda _mode: _FakeContext())
    monkeypatch.setattr(
        sys, "argv", [str(PROBE_1182), "--output", str(out_dir), "--allow-temporary-output"]
    )
    assert probe.main() == 0

    summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["source_restored"] is True
    assert summary["source_unchanged"] is True
    assert summary["input_sha256"] == summary["final_sha256"]
    assert len(summary["cases"]) == 1 + len(probe.MUTATIONS)
    for case_record in summary["cases"]:
        assert case_record["intended_outcome"] is True
        assert Path(case_record["loaded_module_origin"]).resolve().is_relative_to(out_dir.resolve())

    assert file_sha256(canonical_source) == before_sha
    assert canonical_source.stat().st_mtime_ns == before_mtime


def test_issue1170_probe_runs_isolated_mutations_without_touching_canonical_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#1184: issue1170 probe executes against an isolated tables.py copy without touching src/."""
    probe: Any = _load_probe_module(PROBE_1170, "issue1170_probe_live")
    canonical_source = probe.SOURCE.resolve()
    before_sha = file_sha256(canonical_source)
    before_mtime = canonical_source.stat().st_mtime_ns

    _forbidden_canonical_writes(
        monkeypatch, {canonical_source, (probe.ROOT / probe.TEST_FILE).resolve()}
    )

    # Verify foreign checkout rejection
    monkeypatch.setattr(mtp, "__file__", str(tmp_path / "foreign/src/mtp/__init__.py"))
    with pytest.raises(AssertionError, match="different checkout"):
        probe.verify_checkout_identity()
    monkeypatch.setattr(mtp, "__file__", str(REPO_ROOT / "src/mtp/__init__.py"))

    # Run full probe with in-process execution of probe.worker
    class _InProcessWorker:
        def __init__(self, target: Any, args: tuple[Any, ...]) -> None:
            self._target = target
            self._args = args
            self.exitcode: int | None = None

        def start(self) -> None:
            try:
                self._target(*self._args)
            except SystemExit as exc:
                self.exitcode = int(exc.code) if isinstance(exc.code, int) else 1

        def join(self, timeout: float | None = None) -> None:
            _ = timeout
            return None

        def is_alive(self) -> bool:
            return False

    class _InProcessContext:
        def Process(self, target: Any, args: tuple[Any, ...]) -> _InProcessWorker:
            return _InProcessWorker(target=target, args=args)

    out_dir = tmp_path / "stool_probe_out"
    monkeypatch.setattr(probe.multiprocessing, "get_context", lambda _mode: _InProcessContext())
    monkeypatch.setattr(sys, "argv", [str(PROBE_1170), "--output", str(out_dir)])
    assert probe.main() == 0

    summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["source_restored"] is True
    assert summary["source_unchanged"] is True
    assert [row["case"] for row in summary["cases"]] == ["baseline", "contains_quadrant_only"]
    for row in summary["cases"]:
        assert row["intended_outcome"] is True
        origin = Path(row["loaded_module_origin"]).resolve()
        assert origin.is_relative_to(out_dir.resolve())
        assert origin != canonical_source

    assert file_sha256(canonical_source) == before_sha
    assert canonical_source.stat().st_mtime_ns == before_mtime
