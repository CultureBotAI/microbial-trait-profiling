"""Probe the strict profile-loader tests against mutants and old revisions (2026-10-30).

Serves #1295 and #1300 (PR #1292). Each case copies this checkout's ``src`` to
``<output-dir>/<case>/src``, replaces ``mtp/trait_guilds.py`` there with a mutant
or with the file at an older revision, and runs the loader tests of
``tests/test_trait_guilds.py`` in a subprocess whose PYTHONPATH points at the
copy. The tracked source is never rewritten.

Inputs: case names (default: all) and ``--output-dir`` (must not exist yet).
Outputs: per case, ``junit.xml``, ``pytest.log`` and ``evidence.json``; a
``summary.json`` over all cases. Exit status 0 only when every case verified:
the ``fixed`` case passes all selected tests and every other case fails at least
the tests it targets.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TARGET = "src/mtp/trait_guilds.py"
SELECTION = "method_loader"
HEAD_REV = "79a0da103f73321c851b706f49259f0fa1f3d359"
BASE_REV = f"{HEAD_REV}~1"

UNKNOWN_DROPPED = (
    "test_method_loader_rejects_unknown_direction[False-non_human_row]",
    "test_method_loader_rejects_unknown_direction[False-unmatched_row]",
    "test_method_loader_rejects_unknown_direction[True-non_human_row]",
    "test_method_loader_rejects_unknown_direction[True-unmatched_row]",
)
GZ = ("test_method_loader_and_column_resolver_read_tsv_gz_like_tsv",)
METADATA = ("test_method_loader_validates_yet_excludes_metadata_columns_from_total",)
BLANK = ("test_method_loader_rejects_blank_by_default_and_allows_sparse_zero_override",)
MARKER = ("test_method_loader_relabels_only_missing_count_errors_when_the_path_names_the_marker",)
EXACT_LABELS = (
    "test_method_loader_rejects_duplicate_profile_keys[False]",
    "test_method_loader_rejects_duplicate_profile_keys[True]",
    "test_method_loader_rejects_unknown_direction[False-kept_human_row]",
    "test_method_loader_rejects_unknown_direction[True-kept_human_row]",
    *UNKNOWN_DROPPED,
)

DIRECTION_CHECK = """            if inc_dec not in DIRECTIONS:
                msg = (
                    f"Unknown profile direction in {profile_path}:{parsed_row.line_number}: "
                    f"{inc_dec!r}; expected one of {DIRECTIONS}"
                )
                raise ValueError(msg)
"""
HOST_FILTER = """            if human_only and meta["host_species"].strip() != HUMAN_HOST_SPECIES:
                dropped_non_human += 1
                continue
"""

# name -> (list of (target, replacement) edits, tests that must fail)
MUTANTS: dict[str, tuple[tuple[tuple[str, str], ...], tuple[str, ...]]] = {
    # M14 of the review: direction check after the metadata join and host filter.
    "direction_after_host_filter": (
        ((DIRECTION_CHECK, ""), (HOST_FILTER, HOST_FILTER + "\n" + DIRECTION_CHECK)),
        UNKNOWN_DROPPED,
    ),
    # M10: the header reader opens the path as plain text again.
    "header_reader_plain_open": (
        (
            (
                "    with open_input_text(profile_path) as handle:\n"
                "        header = next(csv.reader(handle, delimiter=",
                '    with profile_path.open(encoding="utf-8", newline="") as handle:\n'
                "        header = next(csv.reader(handle, delimiter=",
            ),
        ),
        GZ,
    ),
    # M7: every ValueError from the row loop is relabelled.
    "relabel_every_error": (
        (("        if message.startswith(_READER_MISSING_COUNT_PREFIX):\n", "        if True:\n"),),
        EXACT_LABELS,
    ),
    # #1300: the phrase anywhere in the message triggers the relabel.
    "relabel_on_substring": (
        (
            (
                "        if message.startswith(_READER_MISSING_COUNT_PREFIX):\n",
                "        if _READER_MISSING_COUNT_PREFIX in message:\n",
            ),
        ),
        MARKER,
    ),
    # #1300: the doubled statement that missing is not zero comes back.
    "doubled_wording": (
        (
            (
                'raise ValueError(f"Missing trait count in {detail}") from exc',
                'raise ValueError(f"Missing trait count in {detail}; missing is not zero") '
                "from exc",
            ),
        ),
        BLANK + MARKER,
    ),
}
REVISIONS = {
    # The PR head before this fix: the #1300 label tests must fail there.
    "pr_head_79a0da1": (HEAD_REV, BLANK + MARKER + METADATA),
    # The base before the PR: no gz loader, metadata cells never parsed.
    "base_79a0da1_parent": (BASE_REV, GZ + METADATA),
}
CASES = ("fixed", *MUTANTS, *REVISIONS)


def _source_for(case: str) -> tuple[str, tuple[str, ...]]:
    current = (ROOT / TARGET).read_text(encoding="utf-8")
    if case == "fixed":
        return current, ()
    if case in REVISIONS:
        rev, expected = REVISIONS[case]
        shown = subprocess.run(
            ["git", "-C", str(ROOT), "show", f"{rev}:{TARGET}"],
            check=True,
            capture_output=True,
            text=True,
        )
        return shown.stdout, expected
    edits, expected = MUTANTS[case]
    mutated = current
    for target, replacement in edits:
        if mutated.count(target) != 1:
            msg = f"{case}: expected exactly one mutation target, got {mutated.count(target)}"
            raise ValueError(msg)
        mutated = mutated.replace(target, replacement)
    return mutated, expected


def _junit_results(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Return each test's outcome and, for failures, the first line of the reason."""
    results: dict[str, str] = {}
    reasons: dict[str, str] = {}
    for case in ET.parse(path).getroot().iter("testcase"):
        outcome = "passed"
        problem = case.find("failure")
        if problem is None:
            problem = case.find("error")
        if problem is not None:
            outcome = "failed"
            reasons[case.attrib["name"]] = (problem.attrib.get("message") or "").splitlines()[0]
        elif case.find("skipped") is not None:
            outcome = "skipped"
        results[case.attrib["name"]] = outcome
    return results, reasons


def run_case(case: str, output: Path) -> dict[str, object]:
    case_dir = output / case
    shutil.copytree(
        ROOT / "src",
        case_dir / "src",
        ignore=shutil.ignore_patterns("__pycache__", "*.egg-info"),
    )
    source, expected = _source_for(case)
    (case_dir / TARGET).write_text(source, encoding="utf-8")
    env = {**os.environ, "PYTHONPATH": str(case_dir / "src"), "PYTHONDONTWRITEBYTECODE": "1"}
    imported = subprocess.run(
        [sys.executable, "-c", "import bsdbng.trait_guilds as t; print(t.__file__)"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
        cwd=case_dir,
    ).stdout.strip()
    if Path(imported).resolve() != (case_dir / TARGET).resolve():
        msg = f"{case}: imported {imported}, not the isolated copy"
        raise AssertionError(msg)
    junit = case_dir / "junit.xml"
    with (case_dir / "pytest.log").open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "-q",
                "-p",
                "no:cacheprovider",
                "--rootdir",
                str(ROOT),
                "--basetemp",
                str(case_dir / "pytest-tmp"),
                "--junitxml",
                str(junit),
                str(ROOT / "tests/test_trait_guilds.py"),
                "-k",
                SELECTION,
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            cwd=ROOT,
            check=False,
        )
    results, reasons = _junit_results(junit)
    failed = sorted(name for name, outcome in results.items() if outcome == "failed")
    missing = sorted(name for name in expected if name not in results)
    if case == "fixed":
        verified = completed.returncode == 0 and not failed and bool(results)
    else:
        verified = not missing and all(results[name] == "failed" for name in expected)
    evidence: dict[str, object] = {
        "case": case,
        "imported_module": imported,
        "pytest_exit": completed.returncode,
        "selected": len(results),
        "failed_count": len(failed),
        "failed": failed,
        "expected_to_fail": list(expected),
        "expected_failure_reasons": {name: reasons.get(name, "") for name in expected},
        "expected_missing_from_run": missing,
        "verified": verified,
    }
    (case_dir / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    return evidence


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="*", help=f"default: every case of {CASES}")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    unknown = [case for case in args.cases if case not in CASES]
    if unknown:
        parser.error(f"unknown cases {unknown}; choose from {CASES}")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    tracked = (ROOT / TARGET).read_bytes()
    summary = [run_case(case, output) for case in args.cases or CASES]
    if (ROOT / TARGET).read_bytes() != tracked:
        msg = "Tracked source changed during the isolated probe"
        raise AssertionError(msg)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for row in summary:
        print(
            f"{row['case']}: verified={row['verified']} exit={row['pytest_exit']} "
            f"selected={row['selected']} failed={row['failed_count']}"
        )
    return int(not all(row["verified"] for row in summary))


if __name__ == "__main__":
    raise SystemExit(main())
