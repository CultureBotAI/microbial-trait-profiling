"""Verify isolated audit_trait_coupling docstring mutants (#2015, 2026-11-10).

Purpose:
  Verify that scripts/audit_trait_coupling.py states
  "Evaluates the named-removal comparison against the Table 1 minimum adjacent gap"
  and omits "formal stability decision rule", "inter-row spacing", and "near-twin",
  and that tests/test_audit_trait_coupling_2015.py fails with a native AssertionError
  when any of those invariants or guard checks are mutated.

Inputs:
  - Bounded repository source files (read-only)
  - --scratch-dir / --output-dir: exclusive directory outside the repository for copies and logs
  - --canary: run only the first declared mutation as a canary check
  - --mutation: optional single mutation name or "all" (default)
  - --pytest: optional path to pytest executable

Outputs:
  - Per-mutation JSON receipts, control/mutant pytest logs, and SHA-256 restoration checks
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import stat
import sys
import time
from collections.abc import Generator
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
TEST = "tests/test_audit_trait_coupling_2015.py"
SCRIPT = "scripts/audit_trait_coupling.py"
SUPPORT = (
    "tests/conftest.py",
    "pyproject.toml",
    "uv.lock",
)
LIMIT = 4 * 1024 * 1024
PLUGIN = "guard_2015_probe_plugin"
PREFIX = "MTP_2015_"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read(path: Path) -> bytes:
    before = path.lstat()
    require(
        stat.S_ISREG(before.st_mode) and before.st_size <= LIMIT,
        f"Expected a bounded regular file: {path}",
    )
    with path.open("rb") as handle:
        raw = handle.read(LIMIT + 1)
    after = path.lstat()
    require(
        len(raw) == before.st_size
        and all(
            getattr(before, name) == getattr(after, name)
            for name in ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns")
        ),
        "Input changed while being read",
    )
    return raw


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def save_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")


def mutations() -> dict[str, tuple[str, str, str, str, str]]:
    """Map mutation names to (target, old, new, test_node, expected_assertion_site)."""
    uses = "test_audit_trait_coupling_docstring_uses_named_removal_comparison"
    rejects = "test_audit_trait_coupling_docstring_guard_rejects_legacy_wording"
    return {
        "cli_revert_to_formal_stability_decision_rule": (
            SCRIPT,
            "- Evaluates the named-removal comparison against the Table 1 minimum adjacent gap.",
            "- Evaluates the formal stability decision rule against Table 1 inter-row spacing.",
            uses,
            "assert EXPECTED_BULLET in doc",
        ),
        "cli_omit_named_removal_bullet": (
            SCRIPT,
            "- Evaluates the named-removal comparison against the Table 1 minimum adjacent gap.\n",
            "",
            uses,
            "assert EXPECTED_BULLET in doc",
        ),
        "cli_use_inter_row_spacing_wording": (
            SCRIPT,
            "the Table 1 minimum adjacent gap",
            "Table 1 inter-row spacing",
            uses,
            "assert EXPECTED_BULLET in doc",
        ),
        "cli_add_formal_stability_comment_in_body": (
            SCRIPT,
            "def parse_args(",
            "# Evaluates the formal stability decision rule\ndef parse_args(",
            uses,
            "assert forbidden not in lower_source",
        ),
        "cli_add_near_twin_in_docstring": (
            SCRIPT,
            "Removes coccus + sphere simultaneously",
            "Removes near-twin coccus + sphere simultaneously",
            uses,
            "assert forbidden not in lower_source",
        ),
        "guard_skip_docstring_presence_check": (
            TEST,
            (
                "    assert doc is not None, "
                '"scripts/audit_trait_coupling.py must define a module docstring"'
            ),
            "    if doc is None:\n        doc = EXPECTED_BULLET",
            f"{rejects}[missing-module-docstring]",
            'raise AssertionError(f"Docstring guard accepted legacy wording for {case_id}")',
        ),
        "guard_skip_expected_bullet_check": (
            TEST,
            "    assert EXPECTED_BULLET in doc, (",
            '    assert "Addresses Issue #260:" in doc, (',
            f"{rejects}[missing-named-removal-bullet]",
            'raise AssertionError(f"Docstring guard accepted legacy wording for {case_id}")',
        ),
        "guard_skip_forbidden_wording_check": (
            TEST,
            "        assert forbidden not in lower_source, (",
            "        assert forbidden == forbidden, (",
            f"{rejects}[legacy-decision-rule-in-body-comment]",
            'raise AssertionError(f"Docstring guard accepted legacy wording for {case_id}")',
        ),
    }


def inventory() -> dict[str, str]:
    """Bound only the CLI script, test file, and declared synthetic test support."""
    names = {TEST, SCRIPT, *SUPPORT}
    names.update(spec[0] for spec in mutations().values())
    names.update(path.relative_to(ROOT).as_posix() for path in (ROOT / "src").rglob("*.py"))
    require(len(names) <= 2048, "Implementation inventory is too large")
    result: dict[str, str] = {}
    total = 0
    for name in sorted(names):
        path = ROOT / name
        require(path.resolve().is_relative_to(ROOT), "Input path escapes the repository")
        raw = read(path)
        total += len(raw)
        require(total <= 16 * 1024 * 1024, "Implementation inventory exceeds 16 MiB")
        result[name] = sha(raw)
    return result


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_call(item: pytest.Item) -> None:
    if PREFIX + "ROOT" not in os.environ:
        return
    root = Path(os.environ[PREFIX + "ROOT"])
    selected_test = os.environ[PREFIX + "TEST"]
    require(selected_test == TEST, "Undeclared selected test file")
    module = getattr(item, "module", None)
    if module is None:
        raise RuntimeError("Selected test has no module")
    require(
        Path(str(module.__file__)).resolve() == root / selected_test,
        "Selected test import escaped the isolated copy",
    )
    code = getattr(getattr(item, "obj", None), "__code__", None)
    require(
        code is not None and Path(code.co_filename).resolve() == root / selected_test,
        "Selected test callable escaped the isolated copy",
    )
    require(Path(module.ROOT).resolve() == root, "Test ROOT escaped the isolated copy")
    require(
        Path(module.SCRIPT).resolve() == root / SCRIPT,
        "Test SCRIPT path escaped the isolated copy",
    )


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(
    item: pytest.Item, call: pytest.CallInfo[Any]
) -> Generator[None, Any, None]:
    outcome = yield
    report = outcome.get_result()
    destination = os.environ.get(PREFIX + "EVENTS")
    if destination is None:
        return
    error = call.excinfo
    selected_test = Path(os.environ[PREFIX + "ROOT"]) / os.environ[PREFIX + "TEST"]
    record = {
        "nodeid": item.nodeid,
        "phase": call.when,
        "outcome": report.outcome,
        "exception": error.typename if error else None,
        "sites": [
            {
                "path": str(entry.path),
                "statement": str(entry.statement)
                if Path(str(entry.path)).resolve() == selected_test
                else "",
            }
            for entry in error.traceback
        ]
        if error
        else [],
    }
    with Path(destination).open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")


def run_pytest(
    executable: Path, repository: Path, output: Path, label: str, target: str
) -> dict[str, Any]:
    test_file, separator, test_node = target.partition("::")
    require(
        separator == "::"
        and bool(test_node)
        and test_file == TEST
        and (repository / test_file).is_file(),
        "Pytest invocation does not name a declared copied test node",
    )
    log, events = output / f"{label}.log", output / f"{label}.events.jsonl"
    temporary = output / label
    temporary.mkdir()
    events.touch(exist_ok=False)
    environment = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith(("PYTEST_", "PYTHONPATH", "GIT_", PREFIX))
    }
    environment.update(
        PYTHONPATH=os.pathsep.join(
            (str(repository / "src"), str(repository / "tests"), str(output))
        ),
        PYTHONDONTWRITEBYTECODE="1",
        PYTEST_DISABLE_PLUGIN_AUTOLOAD="1",
        TMPDIR=str(temporary),
        GIT_CONFIG_GLOBAL=os.devnull,
        GIT_CONFIG_NOSYSTEM="1",
        MTP_2015_ROOT=str(repository),
        MTP_2015_TEST=test_file,
        MTP_2015_EVENTS=str(events),
    )
    command = [
        str(executable),
        "-q",
        "--tb=short",
        "-p",
        PLUGIN,
        "-p",
        "no:cacheprovider",
        f"--basetemp={temporary / 'pytest'}",
        target,
    ]
    with log.open("xb") as handle:
        original_cwd = Path.cwd()
        try:
            os.chdir(repository)
            pid = os.posix_spawn(
                str(executable),
                command,
                environment,
                file_actions=[
                    (os.POSIX_SPAWN_DUP2, handle.fileno(), 1),
                    (os.POSIX_SPAWN_DUP2, handle.fileno(), 2),
                ],
            )
        finally:
            os.chdir(original_cwd)
        deadline = time.monotonic() + 120
        while True:
            waited, status = os.waitpid(pid, os.WNOHANG)
            if waited:
                code = os.waitstatus_to_exitcode(status)
                break
            if time.monotonic() >= deadline:
                os.kill(pid, signal.SIGKILL)
                os.waitpid(pid, 0)
                raise TimeoutError("Isolated test exceeded 120 seconds")
            time.sleep(0.05)
    phases = [json.loads(line) for line in read(events).splitlines()]
    calls = [row for row in phases if row["phase"] == "call"]
    nodes = {row["nodeid"] for row in calls}
    complete = (
        bool(calls)
        and len(phases) == 3 * len(calls)
        and all(
            sum(row["nodeid"] == node and row["phase"] == phase for row in phases) == 1
            for node in nodes
            for phase in ("setup", "call", "teardown")
        )
        and all(row["outcome"] == "passed" for row in phases if row["phase"] != "call")
    )
    result = {
        "command": command,
        "exit": code,
        "phases_valid": complete,
        "calls": calls,
        "log_sha256": sha(read(log)),
        "events_sha256": sha(read(events)),
    }
    save_json(output / f"{label}.json", result)
    return result


def verify_mutation(name: str, executable: Path, output: Path) -> bool:
    specs = mutations()
    require(name in specs, f"Unknown mutation: {name}")
    relative_target, old, new, node, assertion_site = specs[name]
    test_target = node if "::" in node else f"{TEST}::{node}"
    selected_test = test_target.split("::", 1)[0]
    require(selected_test == TEST, "Mutation targets an undeclared test file")
    require(
        executable.is_absolute()
        and executable.is_file()
        and not executable.is_symlink()
        and os.access(executable, os.X_OK),
        "An absolute pytest executable is required",
    )
    before = inventory()
    anchor_counts = {}
    target_hashes = {}
    for candidate, (candidate_target, anchor, replacement, _node, _site) in specs.items():
        canonical_target = read(ROOT / candidate_target)
        require(
            sha(canonical_target) == before[candidate_target],
            "Mutation target changed during preflight",
        )
        count = canonical_target.decode().count(anchor)
        require(
            bool(anchor) and count == 1 and anchor != replacement,
            f"Mutation preflight failed for {candidate}: expected one anchor, found {count}",
        )
        anchor_counts[candidate] = count
        target_hashes[candidate_target] = sha(canonical_target)
    save_json(
        output / "ANCHOR_PREFLIGHT.json",
        {
            "target_sha256": target_hashes,
            "anchor_counts": anchor_counts,
        },
    )
    repository = output / "repository"
    repository.mkdir()
    for relative, expected in before.items():
        raw = read(ROOT / relative)
        require(sha(raw) == expected, "Canonical source changed before copying")
        destination = repository / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as handle:
            handle.write(raw)
    helper = read(Path(__file__))
    with (output / f"{PLUGIN}.py").open("xb") as handle:
        handle.write(helper)
    save_json(output / "INPUTS.json", before)
    save_json(
        output / "TOOLS.json",
        {
            "helper_sha256": sha(helper),
            "pytest": str(executable),
            "pytest_sha256": sha(read(executable)),
        },
    )
    source = repository / relative_target
    original = read(source)
    require(bool(old) and original.decode().count(old) == 1, "Mutation anchor is not unique")
    changed = original.decode().replace(old, new).encode()
    require(changed != original, "Mutation does not change source")
    result = {
        "mutation": name,
        "source": relative_target,
        "target": test_target,
        "killed": False,
        "original_sha256": sha(original),
        "mutated_sha256": sha(changed),
    }
    try:
        control = run_pytest(executable, repository, output, "control", test_target)
        require(
            control["exit"] == 0
            and control["phases_valid"]
            and all(row["outcome"] == "passed" for row in control["calls"]),
            "Positive control failed; mutant was not executed",
        )
        source.write_bytes(changed)
        mutant = run_pytest(executable, repository, output, "mutant", test_target)
        result["killed"] = (
            mutant["exit"] == 1
            and mutant["phases_valid"]
            and {row["nodeid"] for row in mutant["calls"]}
            == {row["nodeid"] for row in control["calls"]}
            and all(
                row["outcome"] == "failed"
                and row["exception"] == "AssertionError"
                and any(
                    Path(site["path"]).resolve() == repository / selected_test
                    and assertion_site in site["statement"]
                    for site in row["sites"]
                )
                for row in mutant["calls"]
            )
        )
    finally:
        source.write_bytes(original)
        canonical = inventory()
        restored = {relative: sha(read(repository / relative)) for relative in before}
        save_json(output / "CANONICAL_AFTER.json", canonical)
        save_json(output / "RESTORED.json", restored)
        save_json(output / "result.json", result)
        require(canonical == restored == before, "Canonical/restored source hashes changed")
        require(read(Path(__file__)) == helper, "Verifier changed during execution")
    print(json.dumps(result, indent=2))
    return bool(result["killed"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scratch-dir",
        "--output-dir",
        dest="scratch_dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--canary",
        action="store_true",
        help="Run only the first declared mutation as a canary check.",
    )
    parser.add_argument("--mutation", default="all")
    parser.add_argument("--pytest", type=Path, default=Path(sys.executable).with_name("pytest"))
    args = parser.parse_args()
    scratch = args.scratch_dir.resolve()
    require(
        not scratch.is_relative_to(ROOT) and not ROOT.is_relative_to(scratch),
        "Evidence output must be outside the canonical repository",
    )
    require(not scratch.exists() and not scratch.is_symlink(), "Evidence output must be exclusive")
    scratch.mkdir(parents=True)
    executable = args.pytest.resolve(strict=True)
    specs = mutations()
    if args.canary and args.mutation == "all":
        first = next(iter(specs))
        return 0 if verify_mutation(first, executable, scratch) else 1
    if args.mutation != "all":
        return 0 if verify_mutation(args.mutation, executable, scratch) else 1
    summary: dict[str, bool] = {}
    for name in specs:
        sub = scratch / name
        sub.mkdir()
        killed = verify_mutation(name, executable, sub)
        summary[name] = killed
        require(killed, f"Mutant survived: {name}")
    save_json(scratch / "summary.json", summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
