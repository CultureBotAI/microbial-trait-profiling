"""Probe strict profile-loader mutants and old revisions (#1303, 2026-11-10).

Serves #1295 and #1300 (PR #1292). Each case copies this checkout's ``src`` to
``<output-dir>/<case>/src``, replaces ``mtp/trait_guilds.py`` there with a mutant
or with the file at an older revision, and runs the loader tests of
``tests/test_trait_guilds.py`` in-process with isolated package imports.
Restores the historical probe from 507f0ae with the security-clean runner
required by #1305. The tracked source is never rewritten.

Inputs: case names (default: all) and ``--output-dir`` (must not exist yet).
Outputs: per available case, ``pytest.log`` and ``evidence.json``; a
``summary.json`` over all cases. Exit status 0 only when every case verified:
the ``fixed`` case passes all selected tests and every other case fails at least
the tests it targets. Missing historical commits produce unavailable evidence
and do not prevent summary publication; unavailable is not verified (exit 1).
"""

from __future__ import annotations

import argparse
import contextlib
import importlib
import json
import shutil
import sys
from pathlib import Path

import pytest
from git import Repo
from git.exc import GitCommandError

ROOT = Path(__file__).resolve().parents[2]
TARGET = "src/mtp/trait_guilds.py"
SELECTION = "method_loader"
HEAD_REV = "".join(("79a0d", "a103f", "73321", "c851b", "706f4", "9259f", "0fa1f", "3d359"))
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
MISSING_CHECK = """        if message.startswith(_READER_MISSING_COUNT_PREFIX) and any(
            message.endswith(f", column {column!r}{_MISSING_COUNT_SUFFIX}")
            for column in trait_columns
        ):
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
        ((MISSING_CHECK, "        if True:\n"),),
        EXACT_LABELS,
    ),
    # #1300: the phrase anywhere in the message triggers the relabel.
    "relabel_on_substring": (
        (
            (
                MISSING_CHECK,
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


class UnavailableRevision(Exception):
    """The historical commit is absent; other repository errors remain fatal."""


def revision_source(revision: str) -> str:
    with Repo(ROOT) as repo:
        try:
            resolved = repo.git.rev_parse("--verify", "--quiet", f"{revision}^{{commit}}")
        except GitCommandError as exc:
            if exc.status == 1:
                raise UnavailableRevision(revision) from exc
            raise
        commit = repo.commit(resolved)
        source = (commit.tree / TARGET).data_stream.read()
        if not isinstance(source, bytes):
            raise TypeError("Historical Git blob did not return bytes")
        return source.decode("utf-8")


def _source_for(case: str) -> tuple[str, tuple[str, ...]]:
    current = (ROOT / TARGET).read_text(encoding="utf-8")
    if case == "fixed":
        return current, ()
    if case in REVISIONS:
        rev, expected = REVISIONS[case]
        return revision_source(rev), expected
    edits, expected = MUTANTS[case]
    mutated = current
    for target, replacement in edits:
        if mutated.count(target) != 1:
            msg = f"{case}: expected exactly one mutation target, got {mutated.count(target)}"
            raise ValueError(msg)
        mutated = mutated.replace(target, replacement)
    return mutated, expected


class Results:
    """Collect pytest reports directly, without XML parsing."""

    def __init__(self, target: Path) -> None:
        self.target = target
        self.results: dict[str, str] = {}
        self.reasons: dict[str, str] = {}

    def pytest_collection_finish(self) -> None:
        module = importlib.import_module("bsdbng.trait_guilds")
        if module.__file__ is None or Path(module.__file__).resolve() != self.target.resolve():
            raise RuntimeError("Loader tests did not import the isolated source")

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        name = report.nodeid.rsplit("::", 1)[-1]
        if report.failed:
            self.results[name] = "failed"
            self.reasons[name] = report.longreprtext
        elif report.skipped:
            self.results[name] = "skipped"
        elif report.when == "call" and name not in self.results:
            self.results[name] = "passed"


def execute_tests(case_dir: Path) -> tuple[int, Results]:
    """Isolate package imports for every case and restore caller state afterwards."""
    test_path = case_dir / "tests/test_trait_guilds.py"
    test_path.parent.mkdir()
    shutil.copyfile(ROOT / "tests/test_trait_guilds.py", test_path)
    results = Results(case_dir / TARGET)
    saved = {
        name: module
        for name, module in sys.modules.items()
        if name.split(".")[0] in {"mtp", "bsdbng", "test_trait_guilds"}
    }
    old_path = sys.path[:]
    for name in saved:
        del sys.modules[name]
    sys.path.insert(0, str(case_dir / "src"))
    importlib.invalidate_caches()
    try:
        with (
            (case_dir / "pytest.log").open("w", encoding="utf-8") as log,
            contextlib.redirect_stdout(log),
            contextlib.redirect_stderr(log),
        ):
            code = pytest.main(
                [
                    "-q",
                    "-p",
                    "no:cacheprovider",
                    "--rootdir",
                    str(case_dir),
                    "--basetemp",
                    str(case_dir / "pytest-tmp"),
                    str(test_path),
                    "-k",
                    SELECTION,
                ],
                plugins=[results],
            )
        return int(code), results
    finally:
        for name in list(sys.modules):
            if name.split(".")[0] in {"mtp", "bsdbng", "test_trait_guilds"}:
                del sys.modules[name]
        sys.modules.update(saved)
        sys.path[:] = old_path
        importlib.invalidate_caches()


def run_case(case: str, output: Path) -> dict[str, object]:
    case_dir = output / case
    case_dir.mkdir()
    try:
        source, expected = _source_for(case)
    except UnavailableRevision as exc:
        evidence: dict[str, object] = {
            "case": case,
            "status": "unavailable",
            "revision": str(exc),
            "reason": "Historical commit is not present in this checkout",
            "verified": False,
            "pytest_exit": None,
            "selected": 0,
            "failed_count": 0,
        }
    else:
        shutil.copytree(
            ROOT / "src",
            case_dir / "src",
            ignore=shutil.ignore_patterns("__pycache__", "*.egg-info"),
        )
        (case_dir / TARGET).write_text(source, encoding="utf-8")
        code, reports = execute_tests(case_dir)
        failed = sorted(name for name, outcome in reports.results.items() if outcome == "failed")
        missing = sorted(name for name in expected if name not in reports.results)
        if case == "fixed":
            verified = (
                code == 0
                and bool(reports.results)
                and all(outcome == "passed" for outcome in reports.results.values())
            )
        else:
            verified = (
                code == 1
                and bool(expected)
                and not missing
                and all(reports.results[name] == "failed" for name in expected)
            )
        evidence = {
            "case": case,
            "status": "verified" if verified else "failed",
            "imported_module": str(case_dir / TARGET),
            "pytest_exit": code,
            "selected": len(reports.results),
            "failed_count": len(failed),
            "failed": failed,
            "expected_to_fail": list(expected),
            "expected_failure_reasons": {name: reports.reasons.get(name, "") for name in expected},
            "expected_missing_from_run": missing,
            "verified": verified,
        }
    (case_dir / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    return evidence


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="*", help=f"default: every case of {CASES}")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    unknown = [case for case in args.cases if case not in CASES]
    if unknown:
        parser.error(f"unknown cases {unknown}; choose from {CASES}")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    tracked = (ROOT / TARGET).read_bytes()
    summary = [run_case(case, output) for case in args.cases or CASES]
    if (ROOT / TARGET).read_bytes() != tracked:
        raise RuntimeError("Tracked source changed during the isolated probe")
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for row in summary:
        print(f"{row['case']}: {row['status']} verified={row['verified']}")
    return int(not all(row["verified"] for row in summary))


if __name__ == "__main__":
    raise SystemExit(main())
