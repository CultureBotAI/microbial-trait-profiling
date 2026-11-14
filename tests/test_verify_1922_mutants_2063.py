"""The #1922 verifier reports the marker it actually matched (#2063)."""

import contextlib
import importlib.util
import io
import sys
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[1]
VERIFIER_PATH = ROOT / "scripts/adhoc/verify_1922_mutants_agy_20261108.py"
CASES: list[tuple[int, str, str | None]] = [
    (1, "AssertionError: mismatch", "AssertionError"),
    (1, "Failed: DID NOT RAISE <class 'ValueError'>", "Failed: DID NOT RAISE"),
    (1, "Failed: DID NOT RAISE\nAssertionError", "AssertionError"),
    (1, "AssertionError\nFailed: DID NOT RAISE", "AssertionError"),
    (1, "unrecognized failure", None),
    (1, "KeyError: broken fixture\nAssertionError", None),
    (1, "ImportError: broken fixture\nFailed: DID NOT RAISE", None),
    (0, "AssertionError", None),
    (2, "AssertionError", None),
]


def load_verifier() -> ModuleType:
    spec = importlib.util.spec_from_file_location("verifier2063", VERIFIER_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(("code", "reason", "marker"), CASES)
def test_main_reports_only_verified_matched_marker(
    monkeypatch: pytest.MonkeyPatch, code: int, reason: str, marker: str | None
) -> None:
    verifier = load_verifier()
    mutant = verifier.Mutant(
        "M8_probe", Path("target.py"), "before", "after", Path("test.py"), "probe"
    )
    monkeypatch.setattr(verifier, "MUTANTS", (mutant,))

    def setup(root: Path) -> None:
        (root / "target.py").write_text("before")

    calls = 0

    def run(root: Path, _test: Path, _filter: str) -> tuple[int, str]:
        nonlocal calls
        calls += 1
        if calls <= 2:
            return 0, "PYTEST_ORIGIN_OK\ncontrol passed"
        assert (root / "target.py").read_text() == "after"
        return code, f"PYTEST_ORIGIN_OK\n{reason}"

    monkeypatch.setattr(verifier, "_setup_scratch", setup)
    monkeypatch.setattr(verifier, "_run_case", run)
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        result = verifier.main()
    assert calls == 3, f"expected 3 runner calls, got {calls}"
    lines = output.getvalue().splitlines()
    killed = [line for line in lines if line.startswith("KILLED ")]
    if marker is None:
        assert result == 1, f"expected verifier rc=1 for {(code, reason)}, got {result}"
        assert killed == [], f"unexpected KILLED line for rejected case: {killed}"
        assert not any(line.startswith("ALL ") for line in lines), (
            f"unexpected ALL summary line for rejected case: {lines}"
        )
    else:
        expected = f"KILLED M8_probe: rc=1 {marker} PYTEST_ORIGIN_OK"
        assert result == 0, f"expected verifier rc=0 for {(code, reason)}, got {result}"
        assert killed == [expected], f"expected {[expected]}, got {killed}"
        assert "ALL 1 #1921/#1922 MUTANTS KILLED CLEANLY." in lines, (
            f"missing ALL summary line in {lines}"
        )
