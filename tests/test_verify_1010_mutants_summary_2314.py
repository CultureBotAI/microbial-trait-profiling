"""Verify that verify_1010_mutants_agy_20261109.py separates and reports kill classes (#2314)."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
VERIFIER_PATH = REPO_ROOT / "scripts/adhoc/verify_1010_mutants_agy_20261109.py"


def _load_verifier() -> ModuleType:
    spec = importlib.util.spec_from_file_location("verify_1010_mutants_summary_2314", VERIFIER_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _make_report(
    nodeid: str,
    when: str,
    *,
    passed: bool = False,
    failed: bool = False,
    message: str = "",
) -> SimpleNamespace:
    return SimpleNamespace(
        nodeid=nodeid,
        when=when,
        passed=passed,
        failed=failed,
        longrepr=SimpleNamespace(reprcrash=SimpleNamespace(message=message)) if failed else None,
    )


def _setup_fake_repo(tmp_path: Path, verifier: ModuleType) -> Path:
    repo = tmp_path / "repo"
    for rel in (verifier.TARGET_REL, verifier.HISTOGRAM_TEST_REL, verifier.TAXONOMY_TEST_REL):
        target = repo / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("line_a = 1\nline_b = 2\nline_c = 3\n", encoding="utf-8")
    return repo


def test_observer_counts_assertion_error_and_did_not_raise_separately() -> None:
    verifier = _load_verifier()
    plugin = verifier._OutcomeCapturePlugin(
        Path("tests/test_taxa_coverage_taxonomy.py"),
        Path("scripts/create_taxa_coverage_histogram.py"),
    )
    reports = [
        _make_report("tests/t.py::test_pass", "call", passed=True),
        _make_report(
            "tests/t.py::test_assert_exc",
            "call",
            failed=True,
            message="AssertionError: expected resolved_root",
        ),
        _make_report(
            "tests/t.py::test_assert_stmt",
            "call",
            failed=True,
            message="assert 'a' == 'b'",
        ),
        _make_report(
            "tests/t.py::test_embedded_assert",
            "call",
            failed=True,
            message="Failed: AssertionError inside helper",
        ),
        _make_report(
            "tests/t.py::test_dnr_value_error",
            "call",
            failed=True,
            message="Failed: DID NOT RAISE <class 'ValueError'>",
        ),
        _make_report(
            "tests/t.py::test_dnr_assertion_error",
            "call",
            failed=True,
            message="Failed: DID NOT RAISE <class 'AssertionError'>",
        ),
        # Non-call phases and non-assertion call failures must remain in other_failures.
        _make_report(
            "tests/t.py::test_setup_ae",
            "setup",
            failed=True,
            message="AssertionError: broken fixture",
        ),
        _make_report(
            "tests/t.py::test_teardown_dnr",
            "teardown",
            failed=True,
            message="Failed: DID NOT RAISE ValueError",
        ),
        _make_report(
            "tests/t.py::test_collect_err",
            "collect",
            failed=True,
            message="AssertionError: import guard failed",
        ),
        _make_report(
            "tests/t.py::test_call_value_err",
            "call",
            failed=True,
            message="ValueError: unexpected runtime exception",
        ),
    ]
    for report in reports:
        plugin.pytest_runtest_logreport(report)

    assert plugin.passed == 1
    assert plugin.assertion_failures == 3
    assert plugin.did_not_raise_failures == 2
    assert plugin.other_failures == [
        "tests/t.py::test_setup_ae:setup:AssertionError: broken fixture",
        "tests/t.py::test_teardown_dnr:teardown:Failed: DID NOT RAISE ValueError",
        "tests/t.py::test_collect_err:collect:AssertionError: import guard failed",
        "tests/t.py::test_call_value_err:call:ValueError: unexpected runtime exception",
    ]


def test_main_prints_truthful_kill_class_counts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    verifier = _load_verifier()
    fake_repo = _setup_fake_repo(tmp_path, verifier)
    mutants = (
        ("M1_assert_only", "line_a = 1", "line_a = 0"),
        ("M2_dnr_only", "line_b = 2", "line_b = 0"),
        ("M3_mixed", "line_c = 3", "line_c = 0"),
    )
    monkeypatch.setattr(verifier, "REPO_ROOT", fake_repo)
    monkeypatch.setattr(verifier, "MUTANTS", mutants)

    def fake_stage_scratch(scratch: Path) -> None:
        dst = scratch / verifier.TARGET_REL
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(
            (fake_repo / verifier.TARGET_REL).read_text(encoding="utf-8"),
            encoding="utf-8",
        )

    def make_plugin(passed: int, assertion_msgs: list[str], dnr_msgs: list[str]) -> object:
        plugin = verifier._OutcomeCapturePlugin(
            fake_repo / verifier.TAXONOMY_TEST_REL,
            fake_repo / verifier.TARGET_REL,
        )
        plugin.origin_checked = True
        for idx in range(passed):
            plugin.pytest_runtest_logreport(_make_report(f"t::pass_{idx}", "call", passed=True))
        for idx, msg in enumerate(assertion_msgs):
            plugin.pytest_runtest_logreport(
                _make_report(f"t::ae_{idx}", "call", failed=True, message=msg)
            )
        for idx, msg in enumerate(dnr_msgs):
            plugin.pytest_runtest_logreport(
                _make_report(f"t::dnr_{idx}", "call", failed=True, message=msg)
            )
        return plugin

    run_index = 0

    def fake_run_in_scratch(scratch: Path) -> tuple[int, object]:
        nonlocal run_index
        assert (scratch / verifier.TARGET_REL).is_file()
        current = run_index
        run_index += 1
        if current == 0:
            return 0, make_plugin(10, [], [])
        if current == 1:
            return 1, make_plugin(8, ["AssertionError: bad1", "assert 1 == 2"], [])
        if current == 2:
            return 1, make_plugin(9, [], ["Failed: DID NOT RAISE <class 'ValueError'>"])
        return 1, make_plugin(
            8,
            ["AssertionError: bad2"],
            ["Failed: DID NOT RAISE ValueError"],
        )

    monkeypatch.setattr(verifier, "_stage_scratch", fake_stage_scratch)
    monkeypatch.setattr(verifier, "_run_in_scratch", fake_run_in_scratch)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = verifier.main()
    lines = buf.getvalue().splitlines()

    assert rc == 0
    assert (
        "CONTROL: rc=0 origin_checked=True passed=10 assertion_failures=0 "
        "did_not_raise_failures=0 other_failures=0"
    ) in lines
    assert (
        "M1_assert_only: KILLED (rc=1) origin_checked=True passed=8 "
        "assertion_failures=2 did_not_raise_failures=0 other_failures=0"
    ) in lines
    assert (
        "M2_dnr_only: KILLED (rc=1) origin_checked=True passed=9 "
        "assertion_failures=0 did_not_raise_failures=1 other_failures=0"
    ) in lines
    assert (
        "M3_mixed: KILLED (rc=1) origin_checked=True passed=8 "
        "assertion_failures=1 did_not_raise_failures=1 other_failures=0"
    ) in lines
    assert not any("BY ASSERTIONERROR" in line for line in lines)
    assert not any("by AssertionError only" in line for line in lines)
    assert lines[-1] == (
        "ALL 3/3 MUTANTS KILLED (1 matched AssertionError/assert markers only, "
        "1 matched DID NOT RAISE marker only, 1 matched both; "
        "3 AssertionError/assert-marker and 2 DID-NOT-RAISE-marker call-phase failures); "
        "WORKTREE UNTOUCHED."
    )


@pytest.mark.parametrize("phase", ["setup", "teardown", "collect"])
def test_main_rejects_non_call_failures_even_with_kill_messages(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, phase: str
) -> None:
    verifier = _load_verifier()
    fake_repo = _setup_fake_repo(tmp_path, verifier)
    mutants = (("M1_non_call", "line_a = 1", "line_a = 0"),)
    monkeypatch.setattr(verifier, "REPO_ROOT", fake_repo)
    monkeypatch.setattr(verifier, "MUTANTS", mutants)

    def fake_stage_scratch(scratch: Path) -> None:
        dst = scratch / verifier.TARGET_REL
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(
            (fake_repo / verifier.TARGET_REL).read_text(encoding="utf-8"),
            encoding="utf-8",
        )

    run_index = 0

    def fake_run_in_scratch(_scratch: Path) -> tuple[int, object]:
        nonlocal run_index
        current = run_index
        run_index += 1
        plugin = verifier._OutcomeCapturePlugin(
            fake_repo / verifier.TAXONOMY_TEST_REL,
            fake_repo / verifier.TARGET_REL,
        )
        plugin.origin_checked = True
        if current == 0:
            plugin.pytest_runtest_logreport(_make_report("t::ctrl", "call", passed=True))
            return 0, plugin
        plugin.pytest_runtest_logreport(
            _make_report(
                "t::bad",
                phase,
                failed=True,
                message="Failed: DID NOT RAISE ValueError",
            )
        )
        return 1, plugin

    monkeypatch.setattr(verifier, "_stage_scratch", fake_stage_scratch)
    monkeypatch.setattr(verifier, "_run_in_scratch", fake_run_in_scratch)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = verifier.main()
    lines = buf.getvalue().splitlines()

    assert rc == 1
    assert any(line.startswith("M1_non_call: SURVIVED_OR_ERROR") for line in lines)
    assert any(line.startswith(f"OTHER_FAILURE: t::bad:{phase}:") for line in lines)
    assert not any(line.startswith("ALL ") for line in lines)


def test_main_rejects_control_did_not_raise_and_mixed_teardown_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    verifier = _load_verifier()
    fake_repo = _setup_fake_repo(tmp_path, verifier)
    mutants = (("M1_mixed", "line_a = 1", "line_a = 0"),)
    monkeypatch.setattr(verifier, "REPO_ROOT", fake_repo)
    monkeypatch.setattr(verifier, "MUTANTS", mutants)

    def fake_stage_scratch(scratch: Path) -> None:
        dst = scratch / verifier.TARGET_REL
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(
            (fake_repo / verifier.TARGET_REL).read_text(encoding="utf-8"),
            encoding="utf-8",
        )

    def make_plugin(*reports: object) -> object:
        plugin = verifier._OutcomeCapturePlugin(
            fake_repo / verifier.TAXONOMY_TEST_REL,
            fake_repo / verifier.TARGET_REL,
        )
        plugin.origin_checked = True
        for report in reports:
            plugin.pytest_runtest_logreport(report)
        return plugin

    monkeypatch.setattr(verifier, "_stage_scratch", fake_stage_scratch)

    # Case 1: control fails with DID NOT RAISE -> abort immediately.
    monkeypatch.setattr(
        verifier,
        "_run_in_scratch",
        lambda _scratch: (
            1,
            make_plugin(
                _make_report(
                    "t::ctrl",
                    "call",
                    failed=True,
                    message="Failed: DID NOT RAISE ValueError",
                )
            ),
        ),
    )
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc_ctrl = verifier.main()
    ctrl_lines = buf.getvalue().splitlines()
    assert rc_ctrl == 1
    assert any(
        line.startswith("CONTROL:") and "did_not_raise_failures=1" in line for line in ctrl_lines
    )
    assert not any(line.startswith("M1_mixed:") for line in ctrl_lines)
    assert not any(line.startswith("ALL ") for line in ctrl_lines)

    # Case 2: mutant has a call-phase kill AND a teardown failure -> must not count as killed.
    call_idx = 0

    def fake_run_mixed(_scratch: Path) -> tuple[int, object]:
        nonlocal call_idx
        current = call_idx
        call_idx += 1
        if current == 0:
            return 0, make_plugin(_make_report("t::ctrl", "call", passed=True))
        return 1, make_plugin(
            _make_report("t::kill", "call", failed=True, message="AssertionError: killed"),
            _make_report("t::kill", "teardown", failed=True, message="RuntimeError: cleanup"),
        )

    monkeypatch.setattr(verifier, "_run_in_scratch", fake_run_mixed)
    buf_mixed = io.StringIO()
    with contextlib.redirect_stdout(buf_mixed):
        rc_mixed = verifier.main()
    mixed_lines = buf_mixed.getvalue().splitlines()
    assert rc_mixed == 1
    assert any(line.startswith("M1_mixed: SURVIVED_OR_ERROR") for line in mixed_lines)
    assert any(
        "OTHER_FAILURE: t::kill:teardown:RuntimeError: cleanup" in line for line in mixed_lines
    )
    assert not any(line.startswith("ALL ") for line in mixed_lines)


def test_main_reports_non_assertion_error_containing_marker_as_marker_match(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    verifier = _load_verifier()
    fake_repo = _setup_fake_repo(tmp_path, verifier)
    mutants = (("M1_runtime_marker", "line_a = 1", "line_a = 0"),)
    monkeypatch.setattr(verifier, "REPO_ROOT", fake_repo)
    monkeypatch.setattr(verifier, "MUTANTS", mutants)

    def fake_stage_scratch(scratch: Path) -> None:
        dst = scratch / verifier.TARGET_REL
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(
            (fake_repo / verifier.TARGET_REL).read_text(encoding="utf-8"),
            encoding="utf-8",
        )

    call_idx = 0

    def fake_run(_scratch: Path) -> tuple[int, object]:
        nonlocal call_idx
        current = call_idx
        call_idx += 1
        plugin = verifier._OutcomeCapturePlugin(
            fake_repo / verifier.TAXONOMY_TEST_REL,
            fake_repo / verifier.TARGET_REL,
        )
        plugin.origin_checked = True
        if current == 0:
            plugin.pytest_runtest_logreport(_make_report("t::ctrl", "call", passed=True))
            return 0, plugin
        plugin.pytest_runtest_logreport(
            _make_report(
                "t::runtime_marker",
                "call",
                failed=True,
                message="RuntimeError: AssertionError inside helper",
            )
        )
        return 1, plugin

    monkeypatch.setattr(verifier, "_stage_scratch", fake_stage_scratch)
    monkeypatch.setattr(verifier, "_run_in_scratch", fake_run)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = verifier.main()
    lines = buf.getvalue().splitlines()

    assert rc == 0
    assert (
        "M1_runtime_marker: KILLED (rc=1) origin_checked=True passed=0 "
        "assertion_failures=1 did_not_raise_failures=0 other_failures=0"
    ) in lines
    assert lines[-1] == (
        "ALL 1/1 MUTANTS KILLED (1 matched AssertionError/assert markers only, "
        "0 matched DID NOT RAISE marker only, 0 matched both; "
        "1 AssertionError/assert-marker and 0 DID-NOT-RAISE-marker call-phase failures); "
        "WORKTREE UNTOUCHED."
    )
