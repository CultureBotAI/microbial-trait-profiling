"""Verify #2063 marker reporting with scratch copies (2026-11-10).

Inputs: canonical verifier/tests, --scratch-dir, optional --canary.
Outputs: retained copies and verdicts.json; original files remain unchanged.
"""

import argparse
import importlib.util
import json
from pathlib import Path
from tempfile import mkdtemp
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[2]
HELPER = Path("scripts/adhoc/verify_1922_mutants_agy_20261108.py")
TEST = Path("tests/test_verify_1922_mutants_2063.py")
SELECTION = "matched_marker = next(marker for marker in ASSERTION_MARKERS if marker in out)"
PRINT = 'print(f"KILLED {mutant.name}: rc=1 {matched_marker} PYTEST_ORIGIN_OK")'
MUTANTS = (
    (PRINT, 'print(f"KILLED {mutant.name}: rc=1 AssertionError PYTEST_ORIGIN_OK")'),
    (SELECTION, "matched_marker = ASSERTION_MARKERS[0]"),
    (SELECTION, "matched_marker = ASSERTION_MARKERS[-1]"),
    (
        SELECTION,
        "matched_marker = next(marker for marker in reversed(ASSERTION_MARKERS) if marker in out)",
    ),
    ("if rc != 1:", "if rc not in (0, 1):"),
    ("if any(err in out for err in ERROR_MARKERS):", "if False:"),
    (PRINT, 'print(f"KILLED {mutant.name}: rc=1 {matched_marker}")'),
)


def load_tests(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("checks2063", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load tests")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check(root: Path) -> None:
    tests = load_tests(root / TEST)
    for code, reason, marker in tests.CASES:
        with pytest.MonkeyPatch.context() as patch:
            try:
                tests.test_main_reports_only_verified_matched_marker(patch, code, reason, marker)
            except AssertionError as exc:
                raise AssertionError(f"case={(code, reason, marker)}: {exc}") from exc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--canary", action="store_true")
    args = parser.parse_args()
    originals = {path: (ROOT / path).read_bytes() for path in (HELPER, TEST)}
    source = originals[HELPER].decode()
    mutants = MUTANTS[:1] if args.canary else MUTANTS
    for before, _after in mutants:
        if source.count(before) != 1:
            raise RuntimeError(f"Mutation anchor is not unique: {before}")
    args.scratch_dir.mkdir(parents=True, exist_ok=True)
    run = Path(mkdtemp(prefix="verify2063_", dir=args.scratch_dir))
    verdicts: list[dict[str, object]] = []
    try:
        for index in range(len(mutants) + 1):
            root = run / str(index)
            for path, content in originals.items():
                target = root / path
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(content)
            if index:
                before, after = mutants[index - 1]
                (root / HELPER).write_text(source.replace(before, after, 1))
            try:
                check(root)
            except AssertionError as exc:
                if not index:
                    raise
                detail = str(exc).strip()
                if not detail:
                    raise RuntimeError(f"Mutant {index} raised empty AssertionError") from exc
                verdict = "killed by AssertionError"
                verdicts.append({"case": index, "verdict": verdict, "detail": detail})
            else:
                if index:
                    raise RuntimeError(f"Mutant {index} survived")
                verdict = "all case controls passed"
                verdicts.append({"case": index, "verdict": verdict})
            print(index, verdict)
    finally:
        for path, content in originals.items():
            if (ROOT / path).read_bytes() != content:
                raise RuntimeError(f"Canonical file changed: {path}")
    (run / "verdicts.json").write_text(json.dumps(verdicts, indent=2) + "\n")
    print(f"{len(mutants)}/{len(mutants)} killed; canonical files unchanged; evidence: {run}")


if __name__ == "__main__":
    main()
