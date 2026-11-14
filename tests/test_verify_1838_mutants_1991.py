"""Verify that verify_1838_mutants_agy_20261107.py reports the matched assertion marker (#1991)."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
VERIFIER_PATH = REPO_ROOT / "scripts/adhoc/verify_1838_mutants_agy_20261107.py"


def load_verifier(path: Path | None = None) -> ModuleType:
    target = VERIFIER_PATH if path is None else path
    spec = importlib.util.spec_from_file_location("verify_1838_mutants_1991", target)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_main_prints_matched_marker_from_assertion_markers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    verifier = load_verifier()
    co_source = (REPO_ROOT / verifier.CO_REL).read_text(encoding="utf-8")

    def fake_setup_scratch(tmp_dir: Path) -> None:
        target = tmp_dir / verifier.CO_REL
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(co_source, encoding="utf-8")

    call_index = 0
    markers: tuple[str, ...] = verifier.ASSERTION_MARKERS

    def fake_run_case(tmp_dir: Path) -> tuple[int, str]:
        nonlocal call_index
        assert (tmp_dir / verifier.CO_REL).is_file()
        current = call_index
        call_index += 1
        if current == 0:
            return 0, "PYTEST_ORIGIN_OK\n1 passed\n"
        mutant = verifier.MUTANTS[current - 1]
        if mutant.name == "M15_unknown_ratio_trait_accepted_1898":
            return 1, "PYTEST_ORIGIN_OK\nE   Failed: DID NOT RAISE <class 'ValueError'>\n"
        marker = markers[(current - 1) % len(markers)]
        return 1, f"PYTEST_ORIGIN_OK\nE   {marker}: synthetic failure for {mutant.name}\n"

    monkeypatch.setattr(verifier, "_setup_scratch", fake_setup_scratch)
    monkeypatch.setattr(verifier, "_run_case", fake_run_case)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = verifier.main()
    output = buf.getvalue()

    assert rc == 0
    assert (
        "KILLED M15_unknown_ratio_trait_accepted_1898: rc=1 Failed: DID NOT RAISE PYTEST_ORIGIN_OK"
        in output
    )
    assert (
        "KILLED M15_unknown_ratio_trait_accepted_1898: rc=1 AssertionError PYTEST_ORIGIN_OK"
        not in output
    )
    for idx, mutant in enumerate(verifier.MUTANTS):
        expected_marker = (
            "Failed: DID NOT RAISE"
            if mutant.name == "M15_unknown_ratio_trait_accepted_1898"
            else markers[idx % len(markers)]
        )
        assert f"KILLED {mutant.name}: rc=1 {expected_marker} PYTEST_ORIGIN_OK" in output


def test_main_fails_when_no_assertion_marker_matches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    verifier = load_verifier()
    co_source = (REPO_ROOT / verifier.CO_REL).read_text(encoding="utf-8")

    def fake_setup_scratch(tmp_dir: Path) -> None:
        target = tmp_dir / verifier.CO_REL
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(co_source, encoding="utf-8")

    call_index = 0

    def fake_run_case(tmp_dir: Path) -> tuple[int, str]:
        nonlocal call_index
        assert (tmp_dir / verifier.CO_REL).is_file()
        current = call_index
        call_index += 1
        if current == 0:
            return 0, "PYTEST_ORIGIN_OK\n1 passed\n"
        return 1, "PYTEST_ORIGIN_OK\nunrecognized failure text\n"

    monkeypatch.setattr(verifier, "_setup_scratch", fake_setup_scratch)
    monkeypatch.setattr(verifier, "_run_case", fake_run_case)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = verifier.main()
    output = buf.getvalue()

    assert rc == 1
    assert "no assertion marker in output" in output
    assert "KILLED " not in output
