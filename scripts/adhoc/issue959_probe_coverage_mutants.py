"""Probe #959 coverage callout and scope regressions (2026-10-23).

Inputs: this idle worktree's coverage producer and synthetic geometry/CLI tests.
Outputs: per-mutant pytest logs and source-restoration/expected-assertion evidence.
Mutations are sequential and restored byte-for-byte; no production data are used.
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

SOURCE = REPO_ROOT / "scripts/create_taxa_coverage_histogram.py"
TEST = REPO_ROOT / "tests/test_taxa_coverage_histogram.py"
SELECTION = "narrow_segments or scope_is_explicitly"
MUTANTS = (
    (
        "dropped_label",
        "        text.remove()\n        callout = axis.text(",
        "        text.remove()\n        continue\n        callout = axis.text(",
        "Missing segment callout:",
    ),
    (
        "callout_on_bar",
        "            row + bar_height / 2 + 0.08,",
        "            row,",
        "Callout is on its bar:",
    ),
    (
        "no_leader",
        "            axis.annotate(\n",
        "            continue\n            axis.annotate(\n",
        "Missing segment leader:",
    ),
    (
        "scope_hardcoded_at_entrypoint",
        "    run_coverage(**vars(parser.parse_args()))",
        '    arguments = vars(parser.parse_args())\n    arguments["strain_rollup"] = "off"\n'
        "    run_coverage(**arguments)",
        "AssertionError: assert 'off' == 'phenotype'",
    ),
)


def run_case(output: Path, name: str, expected: str | None) -> dict[str, object]:
    """Require the intended assertion, not a collection failure or other crash."""
    Path(importlib.util.cache_from_source(str(SOURCE))).unlink(missing_ok=True)
    importlib.invalidate_caches()
    for module_name, module in sys.modules.copy().items():
        filename = getattr(module, "__file__", None)
        if filename and Path(filename).resolve() == TEST.resolve():
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
                    str(TEST),
                    "-k",
                    SELECTION,
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
        and (expected is None or (expected in log and "AssertionError" in log)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    output = parser.parse_args().output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    original = SOURCE.read_bytes()
    before = file_sha256(SOURCE)
    source = original.decode("utf-8")
    for name, old, _new, _expected in MUTANTS:
        if source.count(old) != 1:
            raise ValueError(f"Expected exactly one mutation target for {name}")
    results = []
    try:
        results.append(run_case(output, "baseline", None))
        for name, old, new, expected in MUTANTS:
            try:
                SOURCE.write_text(source.replace(old, new), encoding="utf-8")
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
        len(results) != len(MUTANTS) + 1
        or not all(result["verified"] for result in results)
        or evidence["source_sha256_after"] != before
    )


if __name__ == "__main__":
    raise SystemExit(main())
