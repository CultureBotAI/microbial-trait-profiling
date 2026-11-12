"""Verify isolated scan cross-check mutants (#1183, 2026-11-10).

Inputs: the bounded source tree, declared synthetic tests, and a pytest executable.
Outputs: exclusive scratch copies, phase-level test evidence and restoration hashes.
Canonical sources are read only; no Git operation or scientific run is performed.
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
TEST = "tests/test_cellulose_decr_selection_scan_scopes.py"
SOURCE = "src/mtp/cellulose_decr/selection_scan_scopes.py"
PRODUCERS = {
    "tests/test_cellulose_decr_selection_scan_scopes.py": (
        "mtp.cellulose_decr.selection_scan_scopes",
        SOURCE,
        "load_run",
    ),
}
SUPPORT = (
    "tests/test_cellulose_decr_selection_scan.py",
    "methods/calibration_protocol.json",
    "methods/cellulose_decr_selection_scan_protocol.json",
    "methods/cellulose_decr_selection_scan_protocol_v2.json",
    "methods/cellulose_decr_selection_scan.md",
)
LIMIT = 4 * 1024 * 1024


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


def save(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")


def inventory() -> dict[str, str]:
    """Bound only implementation files and the declared synthetic test support."""
    names = {TEST, SOURCE, "tests/conftest.py", "pyproject.toml", "uv.lock", *SUPPORT}
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


def mutations() -> dict[str, tuple[str, str, str, str, str]]:
    """Map names to target, old/new anchors, complete test node and assertion site."""
    cells = "test_1183_trait_crosscheck_rejects_changed_cells"
    summary = "test_1183_trait_crosscheck_rejects_summary_or_lead_disagreement"
    refused = 'raise AssertionError("Inconsistent summary/trait evidence was accepted")'
    exact_number = "return observed.is_finite() and recorded.is_finite() and observed == recorded"
    lead_flags = (
        '    if any(row["is_lead_trait"] is not '
        '(row["trait_id"] == protocol.lead_trait) for row in records):'
    )
    return {
        "tolerate_changed_number": (
            SOURCE,
            exact_number,
            exact_number.replace(
                "observed == recorded", 'abs(observed - recorded) < Decimal("1e-6")'
            ),
            cells + "[tiny-numeric-change]",
            refused,
        ),
        "accept_nan_as_blank": (
            SOURCE,
            '        return actual == ""',
            '        return actual in ("", "NaN")',
            cells + "[null-nan]",
            refused,
        ),
        "coerce_boolean_spelling": (
            SOURCE,
            "        return actual == str(expected)",
            "        return actual == str(expected) or (\n"
            "            isinstance(expected, bool) and actual == str(int(expected))\n"
            "        )",
            cells + "[boolean]",
            refused,
        ),
        "accept_duplicate_traits": (
            SOURCE,
            "    if len(set(identifiers)) != len(identifiers):",
            "    if False:",
            "test_1183_trait_crosscheck_rejects_changed_inventory[duplicate-trait]",
            refused,
        ),
        "skip_lead_identity": (
            SOURCE,
            '    if summary.get("lead_trait") != protocol.lead_trait:',
            "    if False:",
            summary + "[lead-id]",
            refused,
        ),
        "skip_lead_flags": (
            SOURCE,
            lead_flags,
            "    if False:",
            summary + "[other-flag]",
            refused,
        ),
        "skip_load_crosscheck": (
            SOURCE,
            "    _check_trait_summary(directory, summary, traits, protocol)\n",
            "",
            "test_1183_load_run_rechecks_bound_trait_rows_after_verification",
            'raise AssertionError("Loading accepted bound trait rows that contradict the summary")',
        ),
        "skip_pair_crosscheck": (
            SOURCE,
            "        _check_trait_summary(run.directory, run.summary, run.traits, protocol)\n",
            "",
            "test_1183_pairing_rechecks_preloaded_trait_rows",
            'raise AssertionError("Pairing accepted contradictory preloaded trait rows")',
        ),
        "label_future_as_historical": (
            SOURCE,
            "            if versions == (VERSION_V1, VERSION_V1)\n",
            "            if True\n",
            "test_1183_protocol_refusal_only_calls_matching_v1_historical[future]",
            'assert ("records remain historical evidence" in message) == historical',
        ),
    }


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_call(item: pytest.Item) -> None:
    if "MTP_1183_ROOT" not in os.environ:
        return
    root = Path(os.environ["MTP_1183_ROOT"])
    selected_test = os.environ["MTP_1183_TEST"]
    require(selected_test in PRODUCERS, "Undeclared selected test file")
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
    module_name, producer_source, function_name = PRODUCERS[selected_test]
    producer = sys.modules.get(module_name)
    if producer is None:
        raise RuntimeError("Selected producer module is not imported")
    require(
        Path(str(producer.__file__)).resolve() == root / producer_source,
        "Producer module escaped the isolated copy",
    )
    function = getattr(producer, function_name, None)
    code = getattr(function, "__code__", None)
    require(
        code is not None and Path(code.co_filename).resolve() == root / producer_source,
        "Producer callable escaped the isolated copy",
    )
    paths = sys.modules.get("bsdbng.paths")
    require(paths is not None and root == paths.REPO_ROOT, "Repository root escaped the copy")


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(
    item: pytest.Item, call: pytest.CallInfo[Any]
) -> Generator[None, Any, None]:
    outcome = yield
    report = outcome.get_result()
    destination = os.environ.get("MTP_1183_EVENTS")
    if destination is None:
        return
    error = call.excinfo
    selected_test = Path(os.environ["MTP_1183_ROOT"]) / os.environ["MTP_1183_TEST"]
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
        and test_file in PRODUCERS
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
        if not key.startswith(("PYTEST_", "PYTHONPATH", "GIT_"))
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
        MTP_1183_ROOT=str(repository),
        MTP_1183_TEST=target.split("::", 1)[0],
        MTP_1183_EVENTS=str(events),
    )
    command = [
        str(executable),
        "-q",
        "--tb=short",
        "-p",
        "guard_1183_probe_plugin",
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
                raise TimeoutError("Isolated reader test exceeded 120 seconds")
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
    save(output / f"{label}.json", result)
    return result


def verify_mutation(name: str, executable: Path, output: Path) -> bool:
    specs = mutations()
    require(name in specs, "Unknown mutation")
    relative_target, old, new, node, assertion_site = specs[name]
    test_target = node if "::" in node else f"{TEST}::{node}"
    selected_test = test_target.split("::", 1)[0]
    require(selected_test in PRODUCERS, "Mutation targets an undeclared test file")
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
    save(
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
    with (output / "guard_1183_probe_plugin.py").open("xb") as handle:
        handle.write(helper)
    save(output / "INPUTS.json", before)
    save(
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
        save(output / "CANONICAL_AFTER.json", canonical)
        save(output / "RESTORED.json", restored)
        save(output / "result.json", result)
        require(canonical == restored == before, "Canonical/restored source hashes changed")
        require(read(Path(__file__)) == helper, "Verifier changed during execution")
    print(json.dumps(result, indent=2))
    return bool(result["killed"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mutation", required=True)
    parser.add_argument("--pytest", type=Path, default=Path(sys.executable).with_name("pytest"))
    args = parser.parse_args()
    output = args.output_dir.resolve()
    require(
        not output.resolve().is_relative_to(ROOT) and not ROOT.is_relative_to(output.resolve()),
        "Evidence output must be outside the canonical repository",
    )
    require(not output.exists() and not output.is_symlink(), "Evidence output must be exclusive")
    output.mkdir(parents=True)
    return 0 if verify_mutation(args.mutation, args.pytest.resolve(strict=True), output) else 1


if __name__ == "__main__":
    raise SystemExit(main())
