"""Probe #967 direction verification using isolated source copies (2026-10-23).

Inputs: a mutation name and this checkout's synthetic parser/replay tests.
Outputs: pytest log and JSON evidence in a fresh requested directory.
The tracked runtime is never rewritten; run each case in a fresh locked process.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from mtp.paths import file_sha256
from mtp.workflow import arms
from mtp.workflow.arms import runtime

ROOT = Path(__file__).resolve().parents[2]
CASES = {
    "known_only": (
        "parser.parse_args(arguments) if registered else parser.parse_known_args(arguments)[0]",
        "parser.parse_known_args(arguments)[0]",
        "arguments_rejected_by_real_figure_parser or rehashed_record",
        "Failed: DID NOT RAISE ValueError",
    ),
    "case_sensitive": (
        "name.casefold() == relative.casefold()",
        "name == relative",
        "case_changed_registered_script or rehashed_record",
        "Failed: DID NOT RAISE ValueError",
    ),
    "hide_recorded_commit": (
        'f"{record.git_commit} and locked environment; retain the original checkout "',
        '"unknown revision and locked environment; retain the original checkout "',
        "missing_required_heatmap",
        "AssertionError: assert",
    ),
    "unregister_profile": (
        '        "scripts/create_study_trait_profile_heatmap.py",\n',
        "",
        "every_direction_aware_jed",
        "AssertionError: assert",
    ),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("baseline", *CASES))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = ROOT / "src/mtp/workflow/arms/runtime.py"
    if Path(runtime.__file__).resolve() != source:
        raise AssertionError("Imported runtime is outside this checkout")
    original = source.read_bytes()
    changed = original.decode("utf-8")
    selection = (
        "direction_aware_jed or real_figure_parser or case_changed "
        "or missing_required_heatmap or rehashed_record"
    )
    marker: str | None = None
    if args.case != "baseline":
        target, replacement, selection, marker = CASES[args.case]
        if changed.count(target) != 1:
            raise ValueError("Expected exactly one mutation target")
        changed = changed.replace(target, replacement)
    copied = output / "runtime.py"
    copied.write_text(changed, encoding="utf-8")
    spec = importlib.util.spec_from_file_location(runtime.__name__, copied)
    if spec is None or spec.loader is None:
        raise AssertionError("Cannot load isolated runtime")
    module = importlib.util.module_from_spec(spec)
    try:
        sys.modules[runtime.__name__] = module
        arms.runtime = module
        spec.loader.exec_module(module)
        with (
            (output / "pytest.log").open("w", encoding="utf-8") as log,
            redirect_stdout(log),
            redirect_stderr(log),
        ):
            code = int(
                pytest.main(
                    [
                        "-q",
                        "-p",
                        "no:cacheprovider",
                        "--basetemp",
                        str(output / "fixtures"),
                        str(ROOT / "tests/test_arm_direction_parser.py"),
                        str(ROOT / "tests/test_workflow_subset_integration.py"),
                        "-k",
                        selection,
                    ]
                )
            )
    finally:
        sys.modules[runtime.__name__] = runtime
        arms.runtime = runtime
    evidence = {
        "case": args.case,
        "pytest_exit": code,
        "expected_assertion": marker,
        "verified": code == (0 if marker is None else 1)
        and (marker is None or marker in (output / "pytest.log").read_text()),
        "source_unchanged": source.read_bytes() == original,
        "source_sha256": file_sha256(source),
        "mutant_sha256": file_sha256(copied),
    }
    (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))
    if not evidence["source_unchanged"]:
        raise AssertionError("Tracked source changed during isolated probe")
    return int(not evidence["verified"])


if __name__ == "__main__":
    raise SystemExit(main())
