#!/usr/bin/env python3
"""Adversarial mutation verifier for issue #1722 (2026-11-10).

Purpose:
    Verify that the full 40-character SHA regression tests in
    tests/test_comparison_spec_fixture_extractor.py kill at least 6 distinct
    mutations in isolated scratch copies of scripts/adhoc/extract_comparison_specs_f0cb5b7.py
    and tests/fixtures/comparison_specs/README.md.
Inputs:
    Repository root (defaults to two parent directories above this script).
Outputs:
    Exits 0 when the baseline passes and every mutant is killed by AssertionError.
"""

from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path
from types import ModuleType

import pytest


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module spec from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _run_suite_against_copies(
    test_module: ModuleType,
    helper_path: Path,
    readme_path: Path,
    scratch_dir: Path,
) -> None:
    extractor = _load_module("mutant_extractor_1722", helper_path)
    monkeypatch = pytest.MonkeyPatch()
    try:
        monkeypatch.setattr(test_module, "HELPER", helper_path)
        monkeypatch.setattr(test_module, "README", readme_path)
        test_module.test_source_commit_and_readme_pin_full_40_char_sha(extractor)
        run_dir = scratch_dir / "run_tmp"
        run_dir.mkdir(parents=True, exist_ok=True)
        test_module.test_pinned_git_extraction_reproduces_fixture_bytes_and_refuses_overwrite(
            extractor,
            monkeypatch,
            run_dir,
        )
    finally:
        monkeypatch.undo()


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    helper_src = repo_root / "scripts/adhoc/extract_comparison_specs_f0cb5b7.py"
    readme_src = repo_root / "tests/fixtures/comparison_specs/README.md"
    test_src = repo_root / "tests/test_comparison_spec_fixture_extractor.py"

    test_module = _load_module("test_comparison_spec_fixture_extractor_1722", test_src)
    full_sha = "".join(("f0cb5", "b7a50", "0b113", "20ea8", "4e63b", "faf12", "5fe15", "20df9"))
    helper_text = helper_src.read_text(encoding="utf-8")
    readme_text = readme_src.read_text(encoding="utf-8")

    with tempfile.TemporaryDirectory(prefix="verify_1722_baseline_") as base_tmp:
        _run_suite_against_copies(test_module, helper_src, readme_src, Path(base_tmp))

    mutants: list[tuple[str, str, str]] = [
        (
            "mutant_1_short_sha_in_source_commit",
            helper_text.replace(
                'SOURCE_COMMIT = "".join('
                '("f0cb5", "b7a50", "0b113", "20ea8", "4e63b", "faf12", "5fe15", "20df9"))',
                'SOURCE_COMMIT = "f0cb5b7"',
                1,
            ),
            readme_text,
        ),
        (
            "mutant_2_wrong_40_char_sha_in_source_commit",
            helper_text.replace('"20df9"', '"20df0"', 1),
            readme_text,
        ),
        (
            "mutant_3_short_sha_in_extractor_docstring",
            helper_text.replace(f"commit {full_sha}.", "commit f0cb5b7.", 1),
            readme_text,
        ),
        (
            "mutant_4_short_sha_in_readme_intro",
            helper_text,
            readme_text.replace(f"at commit `{full_sha}`.", "at commit `f0cb5b7`.", 1),
        ),
        (
            "mutant_5_short_sha_in_readme_git_show",
            helper_text,
            readme_text.replace(
                f"`git show {full_sha}:src/mtp/workflow/downstream_contracts.py`",
                "`git show f0cb5b7:src/mtp/workflow/downstream_contracts.py`",
                1,
            ),
        ),
        (
            "mutant_6_short_sha_in_readme_local_commit",
            helper_text,
            readme_text.replace(
                f"and commit `{full_sha}` available locally",
                "and commit `f0cb5b7` available locally",
                1,
            ),
        ),
        (
            "mutant_7_renamed_fixture_in_readme",
            helper_text,
            readme_text.replace("analysis_specs_f0cb5b7.json", "analysis_specs_full.json"),
        ),
    ]

    killed = 0
    for name, mutated_helper, mutated_readme in mutants:
        if mutated_helper == helper_text and mutated_readme == readme_text:
            raise RuntimeError(f"Mutant {name} did not alter any input text")
        with tempfile.TemporaryDirectory(prefix=f"verify_1722_{name}_") as raw_dir:
            scratch = Path(raw_dir)
            helper_copy = scratch / "extract_comparison_specs_f0cb5b7.py"
            readme_copy = scratch / "README.md"
            helper_copy.write_text(mutated_helper, encoding="utf-8")
            readme_copy.write_text(mutated_readme, encoding="utf-8")
            try:
                _run_suite_against_copies(test_module, helper_copy, readme_copy, scratch)
            except AssertionError:
                killed += 1
            else:
                raise RuntimeError(f"Mutant survived: {name}")

    if killed != len(mutants) or killed < 6:
        raise RuntimeError(f"Expected at least 6 killed mutants, got {killed}/{len(mutants)}")
    print(f"All {killed}/{len(mutants)} mutants killed by AssertionError.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
