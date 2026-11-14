"""Check #948/#961 co-assertion mutants in a private copy (2026-11-09).

Inputs: tracked source, synthetic audit tests and one selected mutation.
Outputs: exclusive durable copy, control/mutant logs and JSON proof receipts.
The canonical repository is read only. A private Git snapshot supplies genuine
provenance; canonical Git state and working-tree source are never written.
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
from collections.abc import Generator, Sequence
from pathlib import Path
from typing import Any

import pytest
from git import Actor, Repo

ROOT = Path(__file__).resolve().parents[2]
SOURCE = "src/mtp/coassertion_audit.py"
MATRIX = "src/mtp/taxa_trait_matrix.py"
TEST = "tests/test_coassertion_audit.py"
SCOPE_CASE = "test_scope_counts_preserve_sources_and_do_not_propagate_to_ancestors"
GUILD_CASE = "test_guild_pair_accepts_either_member_per_side"
PARENT_CASE = "test_strain_links_require_ncbi_parent"
REQUIRED = {
    SOURCE,
    MATRIX,
    TEST,
    "tests/conftest.py",
    "pyproject.toml",
    "uv.lock",
    "methods/coassertion_audit.md",
    "methods/calibration_protocol.json",
}
LIMIT = 4 * 1024 * 1024
MUTATIONS = {
    "off_strain_evidence": (
        SOURCE,
        "else admitted.get(subject, set())",
        "else graph.strain_traits.get(subject, set())",
        1,
        SCOPE_CASE,
    ),
    "off_strain_attribution": (
        SOURCE,
        "if admitted.get(strain, set()).intersection",
        "if graph.strain_traits.get(strain, set()).intersection",
        2,
        SCOPE_CASE,
    ),
    "off_strain_inventory": (
        SOURCE,
        "for strain, keys in admitted.items()",
        "for strain, keys in graph.strain_traits.items()",
        1,
        SCOPE_CASE,
    ),
    "same_strain_union": (
        SOURCE,
        "(direct_left and direct_right) or same_strains",
        "(direct_left and direct_right)",
        1,
        SCOPE_CASE,
    ),
    "guild_requires_all_left_keys": (
        SOURCE,
        "if keys.intersection(pair.left)}",
        "if set(pair.left).issubset(keys)}",
        1,
        GUILD_CASE,
    ),
    "guild_uses_first_key_per_side": (
        SOURCE,
        "            left = {node for node, keys in scoped.items() "
        "if keys.intersection(pair.left)}\n"
        "            right = {node for node, keys in scoped.items() "
        "if keys.intersection(pair.right)}",
        "            left = {node for node, keys in scoped.items() if pair.left[0] in keys}\n"
        "            right = {node for node, keys in scoped.items() if pair.right[0] in keys}",
        1,
        GUILD_CASE,
    ),
    "strain_parent_accepts_non_ncbi": (
        MATRIX,
        "if predicate == SUBCLASS_OF and _is_taxon(obj):",
        "if predicate == SUBCLASS_OF:",
        1,
        PARENT_CASE,
    ),
}


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def regular_bytes(path: Path) -> bytes:
    before = path.lstat()
    if not stat.S_ISREG(before.st_mode) or before.st_size > LIMIT:
        raise ValueError(f"Expected a bounded regular file, not a symlink: {path}")
    with path.open("rb") as handle:
        raw = handle.read(LIMIT + 1)
    after = path.lstat()
    if len(raw) != before.st_size or any(
        getattr(before, key) != getattr(after, key)
        for key in ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns")
    ):
        raise ValueError(f"File changed during inspection: {path}")
    return raw


def source_path(root: Path, relative: str) -> Path:
    name = Path(relative)
    if name.is_absolute() or ".." in name.parts or not name.parts:
        raise ValueError("Source path must be a bound relative path")
    current = root
    for part in name.parts:
        current /= part
        if current.is_symlink():
            raise ValueError(f"Source symlink is not permitted: {current}")
    return current


def write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def isolate_copy(repository: Path, output: Path) -> Path:
    """Copy bounded tracked inputs once, without shared inodes or Git metadata."""
    repository = repository.resolve(strict=True)
    output = output.absolute()
    resolved = output.resolve()
    if resolved.is_relative_to(repository) or repository.is_relative_to(resolved):
        raise ValueError("Isolation output must be outside the canonical repository")
    if output.exists() or output.is_symlink():
        raise FileExistsError("Isolation output already exists; require an exclusive directory")
    with Repo(repository, search_parent_directories=False) as repo:
        entries = {
            (os.fspath(name), stage): entry for (name, stage), entry in repo.index.entries.items()
        }
    names = {
        name
        for (name, stage) in entries
        if stage == 0
        and (
            (name.startswith("src/") and name.endswith(".py"))
            or name.startswith("schema/")
            or name in REQUIRED
        )
    }
    if not names >= REQUIRED or any(stage for (name, stage) in entries if name in names):
        raise ValueError("Required tracked inputs are absent or unmerged")
    if len(names) > 2048:
        raise ValueError("Source inventory exceeds the bounded copy limit")
    inputs = {}
    total = 0
    for name in sorted(names):
        if entries[name, 0].mode not in (0o100644, 0o100755):
            raise ValueError(f"Tracked input is not a regular file: {name}")
        raw = regular_bytes(source_path(repository, name))
        inputs[name] = digest(raw)
        total += len(raw)
    if total > 16 * 1024 * 1024:
        raise ValueError("Source inventory exceeds 16 MiB")
    output.mkdir(parents=True)
    isolated = output / "repository"
    isolated.mkdir()
    for name, expected in inputs.items():
        raw = regular_bytes(source_path(repository, name))
        if digest(raw) != expected:
            raise ValueError("Canonical source changed before copying")
        target = isolated / name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as handle:
            handle.write(raw)
    write_json(
        output / "ISOLATION.json",
        {
            "canonical_root": str(repository),
            "repository": str(isolated.resolve()),
            "files": inputs,
            "bytes": total,
        },
    )
    return isolated


def isolation(isolated: Path) -> dict[str, Any]:
    receipt = isolated.parent / "ISOLATION.json"
    if isolated.is_symlink() or not receipt.is_file():
        raise ValueError("Mutation requires an isolated repository receipt")
    record: dict[str, Any] = json.loads(regular_bytes(receipt))
    if record["repository"] != str(isolated.resolve()) or record["canonical_root"] == str(
        isolated.resolve()
    ):
        raise ValueError("Mutation path is not the recorded isolated repository")
    return record


def apply_mutation(
    isolated: Path, relative: str, old: str, new: str, expected_count: int
) -> dict[str, str]:
    """Apply a checked textual mutant only to an original, inventoried copy."""
    record = isolation(isolated)
    path = source_path(isolated, relative)
    if relative not in record["files"]:
        raise ValueError("Mutation source is outside the bound inventory")
    raw = regular_bytes(path)
    if digest(raw) != record["files"][relative]:
        raise ValueError("Mutation source no longer matches its isolated original")
    text = raw.decode("utf-8")
    if expected_count < 1 or not old or text.count(old) != expected_count:
        raise ValueError("Mutation target anchor count changed; inspect the implementation")
    changed = text.replace(old, new).encode("utf-8")
    if changed == raw:
        raise ValueError("Mutation target did not change")
    path.write_bytes(changed)
    return {"original_source_sha256": digest(raw), "mutated_source_sha256": digest(changed)}


def pytest_runtest_call(item: pytest.Item) -> None:
    """Refuse an installed/canonical import before the test body runs."""
    if "MTP_COASSERTION_REPOSITORY" not in os.environ:
        return
    root = Path(os.environ["MTP_COASSERTION_REPOSITORY"])
    module = getattr(item, "module", None)
    if module is None or Path(module.__file__).resolve() != root / TEST:
        raise AssertionError("Test origin escaped isolation")
    for name, relative, function in (
        ("mtp.coassertion_audit", SOURCE, "run_audit"),
        ("mtp.taxa_trait_matrix", MATRIX, "read_kg_trait_edges"),
    ):
        implementation = sys.modules.get(name)
        if implementation is None:
            raise AssertionError(f"Isolated implementation is missing: {name}")
        if Path(str(implementation.__file__)).resolve() != root / relative:
            raise AssertionError("Implementation origin escaped isolation")
        spec = implementation.__spec__
        if spec is None or Path(str(spec.origin)).resolve() != root / relative:
            raise AssertionError("Implementation spec escaped isolation")
        if (
            Path(getattr(implementation, function).__code__.co_filename).resolve()
            != root / relative
        ):
            raise AssertionError("Callable origin escaped isolation")
    paths = sys.modules.get("mtp.paths")
    if paths is None or root != paths.REPO_ROOT:
        raise AssertionError("Repository root escaped isolation")


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(
    item: pytest.Item, call: pytest.CallInfo[Any]
) -> Generator[None, Any, None]:
    outcome = yield
    report = outcome.get_result()
    destination = os.environ.get("MTP_COASSERTION_EVENTS")
    if destination is None:
        return
    error = call.excinfo
    record = {
        "nodeid": item.nodeid,
        "when": call.when,
        "outcome": report.outcome,
        "exception": error.typename if error else None,
        "sites": [
            {"path": str(entry.path), "name": entry.name, "statement": str(entry.statement)}
            for entry in error.traceback
        ]
        if error
        else [],
    }
    with Path(destination).open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")


def execute_pytest(
    pytest_path: Path, isolated: Path, output: Path, label: str, case: str
) -> dict[str, Any]:
    """Execute the absolute pytest executable in a fresh process without a shell."""
    log, events = output / f"{label}.log", output / f"{label}.events.jsonl"
    events.touch(exist_ok=False)
    temporary = output / label
    temporary.mkdir()
    env = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith(("PYTEST_", "PYTHONPATH", "GIT_"))
    }
    env.update(
        PYTHONPATH=os.pathsep.join((str(isolated / "src"), str(isolated / "tests"), str(output))),
        PYTHONDONTWRITEBYTECODE="1",
        PYTEST_DISABLE_PLUGIN_AUTOLOAD="1",
        TMPDIR=str(temporary),
        MTP_COASSERTION_REPOSITORY=str(isolated),
        MTP_COASSERTION_EVENTS=str(events),
    )
    command = [
        str(pytest_path),
        "-q",
        "--tb=short",
        "-p",
        "coassertion_probe_plugin",
        "-p",
        "no:cacheprovider",
        f"--basetemp={temporary / 'pytest'}",
        f"{TEST}::{case}",
    ]
    with log.open("xb") as handle:
        before = Path.cwd()
        try:
            os.chdir(isolated)
            pid = os.posix_spawn(
                str(pytest_path),
                command,
                env,
                file_actions=[
                    (os.POSIX_SPAWN_DUP2, handle.fileno(), 1),
                    (os.POSIX_SPAWN_DUP2, handle.fileno(), 2),
                ],
            )
        finally:
            os.chdir(before)
        deadline = time.monotonic() + 120
        while True:
            waited, status = os.waitpid(pid, os.WNOHANG)
            if waited:
                code = os.waitstatus_to_exitcode(status)
                break
            if time.monotonic() >= deadline:
                os.kill(pid, signal.SIGKILL)
                os.waitpid(pid, 0)
                raise TimeoutError("Isolated pytest exceeded 120 seconds")
            time.sleep(0.05)
    records = [json.loads(line) for line in regular_bytes(events).splitlines()]
    calls = [row for row in records if row["when"] == "call"]
    nodeids = {row["nodeid"] for row in calls}
    valid = (
        bool(calls)
        and len(records) == 3 * len(calls)
        and all(
            sum(row["nodeid"] == node and row["when"] == phase for row in records) == 1
            for node in nodeids
            for phase in ("setup", "call", "teardown")
        )
        and all(row["outcome"] == "passed" for row in records if row["when"] != "call")
    )
    result = {
        "command": command,
        "exit": code,
        "calls": calls,
        "phases_valid": valid,
        "log_sha256": digest(regular_bytes(log)),
        "events_sha256": digest(regular_bytes(events)),
    }
    write_json(output / f"{label}.json", result)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mutation", choices=sorted(MUTATIONS), required=True)
    parser.add_argument("--pytest", type=Path, default=Path(sys.executable).with_name("pytest"))
    args = parser.parse_args(argv)
    executable = args.pytest.absolute()
    if not executable.is_file() or executable.is_symlink() or not os.access(executable, os.X_OK):
        raise ValueError("Supply an executable absolute pytest entry point with --pytest")
    output = args.output_dir.absolute()
    isolated = isolate_copy(ROOT, output)
    source, old, new, expected, case = MUTATIONS[args.mutation]
    original = regular_bytes(isolated / source)
    plugin = regular_bytes(Path(__file__))
    (output / "coassertion_probe_plugin.py").write_bytes(plugin)
    write_json(
        output / "TOOLS.json",
        {
            "helper_sha256": digest(plugin),
            "pytest": str(executable),
            "pytest_sha256": digest(regular_bytes(executable)),
        },
    )
    # A new private fixture repository supplies real provenance. Canonical Git
    # metadata is never copied or modified; GitPython's index commit runs no hooks.
    with Repo.init(isolated) as repo:
        repo.index.add(list(isolation(isolated)["files"]))
        identity = Actor("Fixture", "fixture@example.invalid")
        commit = repo.index.commit(
            "Isolated co-assertion fixture", author=identity, committer=identity
        )
        write_json(
            output / "FIXTURE_COMMIT.json", {"sha": commit.hexsha, "canonical_git_written": False}
        )
    record = isolation(isolated)
    result: dict[str, Any] = {
        "mutation": args.mutation,
        "source": source,
        "target": f"{TEST}::{case}",
        "killed": False,
        "canonical_written": False,
    }
    try:
        control = execute_pytest(executable, isolated, output, "control", case)
        if (
            not control["phases_valid"]
            or control["exit"] != 0
            or any(row["outcome"] != "passed" for row in control["calls"])
        ):
            raise ValueError("Positive control failed; no mutant was run")
        result.update(apply_mutation(isolated, source, old, new, expected))
        mutant = execute_pytest(executable, isolated, output, args.mutation, case)
        controls = {row["nodeid"] for row in control["calls"]}
        result["killed"] = (
            mutant["exit"] == 1
            and mutant["phases_valid"]
            and controls == {row["nodeid"] for row in mutant["calls"]}
            and all(
                row["outcome"] == "failed"
                and row["exception"] == "AssertionError"
                and any(
                    Path(site["path"]).resolve() == isolated / TEST and site["name"] == case
                    for site in row["sites"]
                )
                for row in mutant["calls"]
            )
        )
        result["control_calls"] = len(control["calls"])
        result["mutant_calls"] = len(mutant["calls"])
    finally:
        (isolated / source).write_bytes(original)
        canonical_after = {
            name: digest(regular_bytes(source_path(ROOT, name))) for name in record["files"]
        }
        restored = {
            name: digest(regular_bytes(source_path(isolated, name))) for name in record["files"]
        }
        result["restored_source_sha256"] = digest(regular_bytes(isolated / source))
        write_json(output / "CANONICAL_AFTER.json", canonical_after)
        write_json(output / "RESTORED.json", restored)
        write_json(output / "result.json", result)
        if canonical_after != record["files"] or restored != record["files"]:
            raise ValueError("Canonical/copy source bindings changed")
        if regular_bytes(Path(__file__)) != plugin:
            raise ValueError("Mutation helper changed during execution")
    print(json.dumps(result, indent=2))
    return 0 if result["killed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
