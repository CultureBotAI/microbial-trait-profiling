"""Probe both main()-level direction-forwarding regressions (#955; 2026-10-23).

Inputs: this checkout's profile heatmap script and synthetic main-level test.
Outputs: per-case pytest logs and source-restoration/expected-assertion evidence.
Use an otherwise idle worktree: the two call sites are mutated separately and
restored byte-for-byte. No production inputs or outputs are used.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import io
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT, file_sha256

TEST_PATH = REPO_ROOT / "tests/test_study_trait_profile_heatmap.py"
TEST_NAME = "test_main_preserves_direction_basis_in_figures_and_drawn_values"
SOURCE = REPO_ROOT / "scripts/create_study_trait_profile_heatmap.py"
FORWARDING = "direction_basis=args.direction_basis,"


def run_case(output: Path, name: str, expected: str | None) -> dict[str, object]:
    """Require a real assertion failure, not a crash or collection error, for each mutant."""
    Path(importlib.util.cache_from_source(str(SOURCE))).unlink(missing_ok=True)
    importlib.invalidate_caches()
    for module_name, module in sys.modules.copy().items():
        filename = getattr(module, "__file__", None)
        if filename and Path(filename).resolve() == TEST_PATH.resolve():
            sys.modules.pop(module_name)
    captured = io.StringIO()
    with redirect_stdout(captured), redirect_stderr(captured):
        code = int(
            pytest.main(
                [
                    "-q",
                    "-p",
                    "no:cacheprovider",
                    "--basetemp",
                    str(output / name),
                    f"{TEST_PATH}::{TEST_NAME}",
                ]
            )
        )
    log = captured.getvalue()
    (output / f"{name}.log").write_text(log, encoding="utf-8")
    return {
        "case": name,
        "exit_code": code,
        "expected_assertion": expected,
        "verified": code == (0 if expected is None else 1)
        and (expected is None or f"AssertionError: {expected}" in log),
    }


def main() -> int:
    """Run a passing baseline and each independently hard-coded call site."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    output = parser.parse_args().output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    original = SOURCE.read_bytes()
    before = file_sha256(SOURCE)
    source = original.decode("utf-8")
    if source.count(FORWARDING) != 2:
        raise ValueError("Expected exactly two main() direction-forwarding call sites")
    results = []
    try:
        results.append(run_case(output, "baseline", None))
        for name, position, expected in (
            ("profile", source.find(FORWARDING), "CLI basis did not reach profile legend"),
            ("paired", source.rfind(FORWARDING), "CLI basis did not reach paired colorbar"),
        ):
            mutated = (
                source[:position]
                + 'direction_basis="group_1_vs_group_0",'
                + source[position + len(FORWARDING) :]
            )
            try:
                SOURCE.write_text(mutated, encoding="utf-8")
                results.append(run_case(output, name, expected))
            finally:
                SOURCE.write_bytes(original)
    finally:
        if SOURCE.read_bytes() != original:
            SOURCE.write_bytes(original)
        Path(importlib.util.cache_from_source(str(SOURCE))).unlink(missing_ok=True)
        evidence = {
            "results": results,
            "source_sha256_before": before,
            "source_sha256_after": file_sha256(SOURCE),
        }
        (output / "evidence.json").write_text(
            json.dumps(evidence, indent=2) + "\n", encoding="utf-8"
        )
    print(json.dumps(evidence, indent=2))
    return int(
        len(results) != 3
        or not all(result["verified"] for result in results)
        or evidence["source_sha256_after"] != before
    )


if __name__ == "__main__":
    raise SystemExit(main())
