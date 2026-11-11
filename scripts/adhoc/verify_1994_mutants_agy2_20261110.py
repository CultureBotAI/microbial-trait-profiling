"""Verify isolated route-guard suffix and arity mutants (#1994, 2026-11-10).

Purpose:
  Verify that the AST route guard in tests/test_cellulose_tufte_historical_labels.py
  accepts save() with >=3 parameters and Name-bound suffix sequences in save.body or
  tree.body, while rejecting <3 parameters, unbound/non-sequence suffix bindings, and
  ambiguous duplicate bindings.

Inputs:
  - Bounded repository source files (read-only)
  - --scratch-dir: exclusive directory outside the repository for isolated copies and logs
  - --mutation: optional single mutation name (canary) or "all" (default)
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
TEST = "tests/test_cellulose_tufte_historical_labels.py"
SUPPORT = (
    "tests/conftest.py",
    "pyproject.toml",
    "uv.lock",
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


def save_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")


def mutations() -> dict[str, tuple[str, str, str, str, str]]:
    """Map mutation names to (target, old, new, test_node, expected_assertion_site)."""
    tracks = "test_route_guard_tracks_outputs_not_table_layout_or_assumed_svg"
    refuses = "test_route_guard_refuses_underarity_save_or_invalid_suffix_binding"
    return {
        "require_exact_four_args": (
            TEST,
            "    assert len(args) >= 3",
            "    assert len(args) == 4",
            f"{tracks}[save-optional-dpi]",
            "assert len(args) == 4",
        ),
        "require_at_least_four_args": (
            TEST,
            "    assert len(args) >= 3",
            "    assert len(args) >= 4",
            f"{tracks}[save-three-args]",
            "assert len(args) >= 4",
        ),
        "allow_two_save_args": (
            TEST,
            "    assert len(args) >= 3",
            "    assert len(args) >= 2",
            f"{refuses}[two-args]",
            "Expected AssertionError for",
        ),
        "skip_name_iterator_resolution": (
            TEST,
            "    if isinstance(iterable, ast.Name):",
            "    if False:",
            f"{tracks}[module-suffixes]",
            "assert isinstance(iterable, ast.Tuple | ast.List)",
        ),
        "skip_local_suffix_lookup": (
            TEST,
            "        iterable = _assigned_expression(name, save.body)",
            "        iterable = None",
            f"{tracks}[local-suffixes]",
            "assert isinstance(iterable, ast.Tuple | ast.List)",
        ),
        "skip_module_suffix_lookup": (
            TEST,
            "            iterable = _assigned_expression(name, tree.body)",
            "            iterable = None",
            f"{tracks}[module-suffixes]",
            "assert isinstance(iterable, ast.Tuple | ast.List)",
        ),
        "allow_parameter_shadowed_suffixes": (
            TEST,
            (
                "            assert name not in param_names, "
                'f"Parameter shadows suffix binding: {name}"'
            ),
            "            assert True",
            f"{refuses}[shadowed-parameter-suffixes]",
            'raise AssertionError(f"Route guard accepted invalid save AST: {invalid}")',
        ),
        "allow_ambiguous_suffix_binding": (
            TEST,
            '    assert len(values) <= 1, f"Ambiguous route binding: {name}"',
            "    assert True",
            f"{refuses}[duplicate-suffixes]",
            'raise AssertionError(f"Route guard accepted invalid save AST: {invalid}")',
        ),
    }


def inventory() -> dict[str, str]:
    """Bound only implementation files and the declared synthetic test support."""
    names = {TEST, *SUPPORT}
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
    if "MTP_1994_ROOT" not in os.environ:
        return
    root = Path(os.environ["MTP_1994_ROOT"])
    selected_test = os.environ["MTP_1994_TEST"]
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
    producer = sys.modules.get("mtp.cellulose_decr.tufte_plates")
    if producer is None:
        raise RuntimeError("Selected producer module is not imported")
    require(
        Path(str(producer.__file__)).resolve() == root / "src/mtp/cellulose_decr/tufte_plates.py",
        "Producer module escaped the isolated copy",
    )


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(
    item: pytest.Item, call: pytest.CallInfo[Any]
) -> Generator[None, Any, None]:
    outcome = yield
    report = outcome.get_result()
    destination = os.environ.get("MTP_1994_EVENTS")
    if destination is None:
        return
    error = call.excinfo
    selected_test = Path(os.environ["MTP_1994_ROOT"]) / os.environ["MTP_1994_TEST"]
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
        if not key.startswith(("PYTEST_", "PYTHONPATH", "GIT_"))
    }
    environment.update(
        PYTHONPATH=os.pathsep.join(
            (str(repository / "src"), str(repository / "tests"), str(output))
        ),
        PYTHONDONTWRITEBYTECODE="1",
        PYTEST_DISABLE_PLUGIN_AUTOLOAD="1",
        MPLBACKEND="Agg",
        MPLCONFIGDIR=str(temporary / "mpl"),
        TMPDIR=str(temporary),
        GIT_CONFIG_GLOBAL=os.devnull,
        GIT_CONFIG_NOSYSTEM="1",
        MTP_1994_ROOT=str(repository),
        MTP_1994_TEST=test_file,
        MTP_1994_EVENTS=str(events),
    )
    command = [
        str(executable),
        "-q",
        "--tb=short",
        "-p",
        "guard_1994_probe_plugin",
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
    with (output / "guard_1994_probe_plugin.py").open("xb") as handle:
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
