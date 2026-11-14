"""Prove historical-source guard mutations in isolated copies (#1993/#1995).

Inputs: bounded repository sources, historical fixtures, and a pytest executable.
Outputs: exclusive control/mutant logs, native exception sites, and hash receipts.
The original worktree is read only; no Git operation or figure rendering is used.
Date: 2026-11-10.
"""

from __future__ import annotations

import argparse
import ast
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
TEST = "tests/test_cellulose_tufte_statistics_sources.py"
HISTORICAL = "tests/fixtures/cellulose_tufte_historical_sources"
SUPPORT = (
    "tests/conftest.py",
    "tests/test_cellulose_tufte_historical_labels.py",
    "pyproject.toml",
    "uv.lock",
    "methods/calibration_protocol.json",
    "docs/reviews/2026-10-19-ora-provenance-refresh.md",
)
LIMIT = 4 * 1024 * 1024
TOTAL_LIMIT = 16 * 1024 * 1024
PLUGIN = "historical_guard_mutation_plugin"
PREFIX = "MTP_HISTORICAL_MUTATION_"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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
    path.chmod(0o444)


def mutations() -> dict[str, tuple[str, str, str, str, str]]:
    """Target file, unique old/new anchor, selected node, intended assertion."""
    nested = "test_historical_manifest_accepts_documented_nested_sources"
    hidden = "test_historical_manifest_ignores_hidden_entries"
    inventory_site = "assert set(records) == set(_historical_source_names(HISTORICAL))"
    refusal_site = "assert failure is not None, message"
    hidden_parts = 'any(part.startswith(".") for part in relative.parts)'
    byte_guard = 'assert actual == expected, f"Historical fixture changed: {fixture}"'
    return {
        "ora_reads_live_document": (
            TEST,
            'ORA_AUDIT = HISTORICAL / "ora_provenance_refresh_781e92e4.md"',
            'ORA_AUDIT = ROOT / "docs/reviews/2026-10-19-ora-provenance-refresh.md"',
            "test_plate_five_ora_audit_is_a_documented_historical_fixture",
            "assert ORA_AUDIT.parent == HISTORICAL",
        ),
        "omit_ora_parent_guard": (
            TEST,
            "    assert ORA_AUDIT.parent == HISTORICAL, (\n"
            '        "Plate V\'s ORA numbers must read the historical snapshot, not HEAD"\n'
            "    )\n",
            "",
            "test_historical_source_parent_guard_refuses_live_ora_audit",
            refusal_site,
        ),
        "change_frozen_ora_depletion_count": (
            f"{HISTORICAL}/ora_provenance_refresh_781e92e4.md",
            "| Pooled depletion | 4,670 | 1,480,390 | 23 |",
            "| Pooled depletion | 4,670 | 1,480,390 | 24 |",
            "test_historical_fixture_bytes_equal_the_documented_git_blob"
            "[ora_provenance_refresh_781e92e4.md]",
            byte_guard,
        ),
        "scan_top_level_only": (
            TEST,
            'for path in directory.rglob("*"):',
            'for path in directory.glob("*"):',
            nested + "[nested-source]",
            inventory_site,
        ),
        "discard_relative_directories": (
            TEST,
            "sources.append(relative.as_posix())",
            "sources.append(path.name)",
            nested + "[same-basename]",
            inventory_site,
        ),
        "include_hidden_files": (
            TEST,
            hidden_parts,
            "False",
            hidden + "[.DS_Store]",
            inventory_site,
        ),
        "ignore_only_hidden_leaves": (
            TEST,
            hidden_parts,
            'relative.name.startswith(".")',
            hidden + "[.cache/source.md]",
            inventory_site,
        ),
        "exclude_nested_readme_sources": (
            TEST,
            'relative == Path("README.md")',
            'relative.name == "README.md"',
            nested + "[nested-readme]",
            inventory_site,
        ),
        "accept_unlisted_nested_sources": (
            TEST,
            inventory_site,
            "assert set(records) == set(records)",
            "test_historical_manifest_refuses_undocumented_nested_source",
            refusal_site,
        ),
        "accept_changed_nested_bytes": (
            TEST,
            byte_guard,
            'assert expected == expected, f"Historical fixture changed: {fixture}"',
            "test_historical_byte_guard_refuses_changed_nested_source",
            refusal_site,
        ),
        "byte_decorator_scans_top_level_only": (
            TEST,
            "    _historical_source_names(HISTORICAL),\n",
            "    sorted(\n"
            "        path.name\n"
            "        for path in HISTORICAL.iterdir()\n"
            '        if path.is_file() and path.name != "README.md"\n'
            '        and not path.name.startswith(".")\n'
            "    ),\n",
            "test_historical_byte_parametrization_includes_nested_sources",
            "assert list(parametrizations[0].args[1]) == sorted(sources)",
        ),
    }


def inventory() -> dict[str, str]:
    names = {TEST, *SUPPORT}
    names.update(spec[0] for spec in mutations().values())
    for directory, pattern in ((ROOT / "src", "*.py"), (ROOT / HISTORICAL, "*")):
        for path in directory.rglob(pattern):
            require(not path.is_symlink(), "Source inventory contains a symbolic link")
            if path.is_file():
                names.add(path.relative_to(ROOT).as_posix())
                require(len(names) <= 2048, "Source inventory exceeds 2048 files")
    result = {}
    total = 0
    for name in sorted(names):
        path = ROOT / name
        require(path.resolve().is_relative_to(ROOT), "Input escapes the repository")
        require(
            not any(parent.is_symlink() for parent in path.parents if parent != ROOT),
            "Input has a symbolic-link ancestor",
        )
        raw = read(path)
        total += len(raw)
        require(total <= TOTAL_LIMIT, "Source inventory exceeds 16 MiB")
        result[name] = sha(raw)
    return result


def preflight(before: dict[str, str]) -> dict[str, Any]:
    anchors = {}
    for name, (target, old, new, node, site) in mutations().items():
        require(target in before and bool(old) and old != new, "Invalid mutation declaration")
        raw = read(ROOT / target)
        require(sha(raw) == before[target], "Source changed during anchor preflight")
        count = raw.decode("utf-8").count(old)
        require(count == 1, f"Mutation {name} needs one anchor; found {count}")
        if target.endswith(".py"):
            ast.parse(raw.decode("utf-8").replace(old, new), filename=target)
        require(bool(node) and "::" not in node and bool(site), "Invalid test declaration")
        require(site in read(ROOT / TEST).decode("utf-8"), "Intended assertion site is absent")
        anchors[name] = {"source": target, "matches": count, "sha256": sha(raw)}
    require(len(anchors) >= 6, "At least six declared mutations are required")
    return anchors


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_call(item: pytest.Item) -> None:
    if PREFIX + "ROOT" not in os.environ:
        return
    root = Path(os.environ[PREFIX + "ROOT"])
    target = os.environ[PREFIX + "TARGET"]
    require(item.nodeid == target, "Unexpected selected test node")
    module = getattr(item, "module", None)
    if module is None:
        raise RuntimeError("Selected test has no module")
    require(
        Path(str(module.__file__)).resolve() == root / TEST,
        "Test module import escaped the copied repository",
    )
    function = getattr(item, "obj", None)
    code = getattr(function, "__code__", None)
    require(
        code is not None and Path(code.co_filename).resolve() == root / TEST,
        "Test callable escaped the copied repository",
    )
    require(root == module.ROOT, "Test root escaped the copied repository")
    for name in ("_fixture_manifest", "_historical_source_names"):
        function = getattr(module, name, None)
        code = getattr(function, "__code__", None)
        require(
            code is not None and Path(code.co_filename).resolve() == root / TEST,
            f"Guard callable escaped the copied repository: {name}",
        )
    support = sys.modules.get("test_cellulose_tufte_historical_labels")
    require(
        support is not None
        and Path(str(support.__file__)).resolve()
        == root / "tests/test_cellulose_tufte_historical_labels.py",
        "Supporting test module escaped the copied repository",
    )
    for name, imported in tuple(sys.modules.items()):
        if name in {"mtp", "bsdbng"} or name.startswith(("mtp.", "bsdbng.")):
            location = getattr(imported, "__file__", None)
            require(
                location is not None and Path(location).resolve().is_relative_to(root / "src"),
                f"Producer import escaped the copied repository: {name}",
            )
    paths = sys.modules.get("mtp.paths")
    require(paths is not None and root == paths.REPO_ROOT, "Producer root escaped the copy")


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
    selected_test = Path(os.environ[PREFIX + "ROOT"]) / TEST
    record = {
        "nodeid": item.nodeid,
        "phase": call.when,
        "outcome": report.outcome,
        "exception": error.typename if error else None,
        "native_assertion": error.type is AssertionError if error else False,
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
    executable: Path, repository: Path, output: Path, label: str, node: str
) -> dict[str, Any]:
    target = f"{TEST}::{node}"
    require(bool(node) and "::" not in node, "Expected a declared test function node")
    require((repository / TEST).is_file(), "Copied test source is missing")
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
        MPLCONFIGDIR=str(temporary / "matplotlib"),
        XDG_CACHE_HOME=str(temporary / "cache"),
        GIT_CONFIG_GLOBAL=os.devnull,
        GIT_CONFIG_NOSYSTEM="1",
    )
    environment.update(
        {
            PREFIX + "ROOT": str(repository),
            PREFIX + "TARGET": target,
            PREFIX + "EVENTS": str(events),
        }
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
                raise TimeoutError("Isolated guard test exceeded 120 seconds")
            time.sleep(0.05)
    phases = [json.loads(line) for line in read(events).splitlines()]
    complete = (
        len(phases) == 3
        and all(row["nodeid"] == target for row in phases)
        and {row["phase"] for row in phases} == {"setup", "call", "teardown"}
        and all(row["outcome"] == "passed" for row in phases if row["phase"] != "call")
    )
    result = {
        "command": command,
        "exit": code,
        "phases_valid": complete,
        "calls": [row for row in phases if row["phase"] == "call"],
        "log_sha256": sha(read(log)),
        "events_sha256": sha(read(events)),
    }
    log.chmod(0o444)
    events.chmod(0o444)
    save(output / f"{label}.json", result)
    return result


def verify_mutation(
    name: str, executable: Path, output: Path, before: dict[str, str], helper: bytes
) -> dict[str, Any]:
    relative_target, old, new, node, site = mutations()[name]
    repository = output / "repository"
    repository.mkdir()
    for relative, expected in before.items():
        raw = read(ROOT / relative)
        require(sha(raw) == expected, "Canonical source changed before copying")
        destination = repository / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as handle:
            handle.write(raw)
    with (output / f"{PLUGIN}.py").open("xb") as handle:
        handle.write(helper)
    source = repository / relative_target
    original = read(source)
    require(original.decode("utf-8").count(old) == 1, "Copied mutation anchor changed")
    changed = original.decode("utf-8").replace(old, new).encode("utf-8")
    result = {
        "mutation": name,
        "source": relative_target,
        "target": f"{TEST}::{node}",
        "assertion_site": site,
        "killed": False,
        "original_sha256": sha(original),
        "mutated_sha256": sha(changed),
    }
    try:
        control = run_pytest(executable, repository, output, "control", node)
        require(
            control["exit"] == 0
            and control["phases_valid"]
            and all(row["outcome"] == "passed" for row in control["calls"]),
            "Positive control failed; mutant was not executed",
        )
        source.write_bytes(changed)
        require(read(source) == changed, "Mutation write did not persist")
        mutant = run_pytest(executable, repository, output, "mutant", node)
        result["killed"] = (
            mutant["exit"] == 1
            and mutant["phases_valid"]
            and all(
                row["outcome"] == "failed"
                and row["native_assertion"]
                and any(
                    Path(entry["path"]).resolve() == repository / TEST
                    and site in entry["statement"]
                    for entry in row["sites"]
                )
                for row in mutant["calls"]
            )
        )
    finally:
        source.write_bytes(original)
        restored = {relative: sha(read(repository / relative)) for relative in before}
        save(output / "RESTORED.json", restored)
        save(output / "result.json", result)
        require(restored == before, "Restored source inventory differs from the original")
        require(
            read(output / f"{PLUGIN}.py") == helper,
            "Copied verifier changed during execution",
        )
    return result


def seal_receipts(output: Path) -> None:
    """Bind the small proof files, excluding copied sources and disposable caches."""
    directories = [output, *(path for path in output.iterdir() if path.is_dir())]
    receipt = {}
    for directory in directories:
        for path in directory.iterdir():
            if path.is_file():
                receipt[path.relative_to(output).as_posix()] = sha(read(path))
                path.chmod(0o444)
    save(output / "RECEIPTS.json", receipt)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--mutation", choices=tuple(mutations()))
    parser.add_argument("--pytest", type=Path, default=Path(sys.executable).with_name("pytest"))
    args = parser.parse_args()
    output = args.scratch_dir.resolve()
    require(
        not output.is_relative_to(ROOT) and not ROOT.is_relative_to(output),
        "Evidence output must be outside the canonical repository",
    )
    require(not output.exists() and not output.is_symlink(), "Evidence output must be exclusive")
    executable = args.pytest.resolve(strict=True)
    require(
        executable.is_file() and os.access(executable, os.X_OK) and executable.name == "pytest",
        "An executable named pytest is required",
    )
    before = inventory()
    helper = read(Path(__file__))
    executable_bytes = read(executable)
    output.mkdir(parents=True)
    save(output / "INPUTS.json", before)
    save(
        output / "TOOLS.json",
        {
            "helper_sha256": sha(helper),
            "pytest": str(executable),
            "pytest_sha256": sha(executable_bytes),
            "python_version": sys.version,
        },
    )
    results = []
    try:
        save(output / "ANCHOR_PREFLIGHT.json", preflight(before))
        names = [args.mutation] if args.mutation else list(mutations())
        for name in names:
            attempt = output / name
            attempt.mkdir()
            result = verify_mutation(name, executable, attempt, before, helper)
            results.append(result)
            require(bool(result["killed"]), f"Mutation survived or failed incorrectly: {name}")
    finally:
        canonical = inventory()
        unchanged = (
            canonical == before
            and read(Path(__file__)) == helper
            and read(executable) == executable_bytes
        )
        save(output / "CANONICAL_AFTER.json", canonical)
        save(output / "SUMMARY.json", {"results": results, "canonical_unchanged": unchanged})
        seal_receipts(output)
        require(unchanged, "Canonical source or verification tool changed during execution")
    print(
        json.dumps(
            {
                "mutants_killed": len(results),
                "canonical_unchanged": True,
                "receipts_sha256": sha(read(output / "RECEIPTS.json")),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
