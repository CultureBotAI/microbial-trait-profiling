"""Verify #1026 output guards in isolated source copies, 2026-11-09.

Inputs: final source/tests plus one declared mutation. Outputs: an exclusive
durable source copy, fresh-pytest control/mutant logs and hash-bound receipts.
No canonical source, output tree or Git metadata is written. Invoke one explicit
mutation per exclusive durable output directory, under the approved RSS cap.
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
from git import Repo

ROOT = Path(__file__).resolve().parents[2]
PRODUCER = "scripts/create_trait_profile_infographic.py"
CONFTEST = "tests/conftest.py"
LIMIT = 4 * 1024 * 1024
ORIGINAL_MKDIR = Path.mkdir
INFO_TEST = "tests/test_trait_profile_infographic.py"
GUARD_TEST = "tests/test_protected_output_guard.py"
WRITER_CASE = (
    INFO_TEST + "::test_write_stratified_trait_stats_emits_feces_excluded_and_body_site_files"
)
FIXTURE_CASE = INFO_TEST + "::test_isolation_discovers_new_output_constants"
GUARD_CASE = GUARD_TEST + "::test_session_guard_refuses_protected_output_changes"
MANIFEST = "trait_profile_paired_contrast_trait_stats_manifest.tsv"
CANARY_FILES = [
    MANIFEST,
    *[
        f"trait_profile_paired_contrast_trait_stats_{token}.tsv"
        for token in ("feces", "feces_excluded", "oral")
    ],
    "uncapped/" + MANIFEST,
]


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
        "Input changed during read",
    )
    return raw


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def save(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def path_in(root: Path, name: str) -> Path:
    relative = Path(name)
    require(
        not relative.is_absolute() and bool(relative.parts) and ".." not in relative.parts,
        "Unsafe relative source/artifact path",
    )
    path = root
    for part in relative.parts:
        path /= part
        require(not path.is_symlink(), "Symlink inside bound source/artifact path")
    return path


def specifications() -> list[dict[str, Any]]:
    """Eight checked faults; nested-guard expectations are outer call assertions."""
    rows = [
        (
            "omit_visualization_root",
            CONFTEST,
            '_PROTECTED_OUTPUT_DIRS = (GUARDED_DATA_DIR, REPO_ROOT / "trait-visualization")',
            "_PROTECTED_OUTPUT_DIRS = (GUARDED_DATA_DIR,)",
            GUARD_TEST + "::test_default_protected_roots_include_data_and_visualization",
            "assert result.returncode == 0",
            [],
        ),
        (
            "bypass_teardown_comparison",
            CONFTEST,
            "if before[root].get(name) != after[root].get(name)",
            "if False",
            GUARD_CASE + "[modify-trait-visualization]",
            "assert result.returncode ==",
            [],
        ),
        (
            "miss_new_files",
            CONFTEST,
            "before[root].keys() | after[root].keys()",
            "before[root].keys()",
            GUARD_CASE + "[create-trait-visualization]",
            "assert result.returncode ==",
            [],
        ),
        (
            "miss_deleted_files",
            CONFTEST,
            "before[root].keys() | after[root].keys()",
            "after[root].keys()",
            GUARD_CASE + "[delete-trait-visualization]",
            "assert result.returncode ==",
            [],
        ),
        (
            "import_time_default",
            PRODUCER,
            "    output_dir: Path | None = None,",
            "    output_dir: Path = INVESTIGATION_DIR,",
            WRITER_CASE + "[default]",
            "Writer escaped the test output root",
            CANARY_FILES,
        ),
        (
            "ignore_explicit_output",
            PRODUCER,
            "    if output_dir is None:\n        output_dir = INVESTIGATION_DIR\n",
            "    output_dir = INVESTIGATION_DIR\n",
            WRITER_CASE + "[explicit]",
            "assert output.exists()",
            CANARY_FILES,
        ),
        (
            "omit_visualization_rebasing",
            INFO_TEST,
            "        original = producer.VIS_DIR",
            "        original = producer.INVESTIGATION_DIR",
            FIXTURE_CASE,
            "VIS_DIR escaped isolation",
            [],
        ),
        (
            "omit_tuple_rebasing",
            INFO_TEST,
            "            if isinstance(value, Path | tuple):",
            "            if isinstance(value, Path):",
            FIXTURE_CASE,
            "p.is_relative_to(tmp_path) for p in infographic.OUTPUT_PAIR",
            [],
        ),
    ]
    return [
        {
            "name": name,
            "source": source,
            "old": old,
            "new": new,
            "count": 1,
            "target": target,
            "assertion_site": assertion,
            "artifacts": artifacts,
        }
        for name, source, old, new, target, assertion, artifacts in rows
    ]


def copy_inputs(output: Path, specs: list[dict[str, Any]]) -> tuple[Path, dict[str, str]]:
    require(
        not output.resolve().is_relative_to(ROOT) and not ROOT.is_relative_to(output.resolve()),
        "Output must be outside the canonical repository",
    )
    require(not output.exists() and not output.is_symlink(), "Output must be exclusive")
    required = {
        PRODUCER,
        CONFTEST,
        "pyproject.toml",
        "uv.lock",
        "methods/calibration_protocol.json",
    }
    required.update(row["target"].split("::", 1)[0] for row in specs)
    with Repo(ROOT, search_parent_directories=False) as repo:
        entries = {
            (os.fspath(name), stage): entry for (name, stage), entry in repo.index.entries.items()
        }
    names = {
        name
        for name, stage in entries
        if stage == 0
        and (
            name in required
            or (name.startswith("src/") and name.endswith(".py"))
            or name.startswith("schema/")
        )
    }
    require(
        names >= required and len(names) <= 2048,
        "Required tracked source missing or copy too large",
    )
    require(not any(stage for name, stage in entries if name in names), "Unmerged source entry")
    before, total = {}, 0
    for name in sorted(names):
        require(entries[name, 0].mode in (0o100644, 0o100755), "Nonregular tracked source")
        raw = read(path_in(ROOT, name))
        before[name] = sha(raw)
        total += len(raw)
    require(total <= 16 * 1024 * 1024, "Source copy exceeds 16 MiB")
    output.mkdir(parents=True)
    repository = output / "repository"
    repository.mkdir()
    for name, expected in before.items():
        raw = read(path_in(ROOT, name))
        require(sha(raw) == expected, "Canonical source changed while copying")
        target = repository / name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as handle:
            handle.write(raw)
    save(
        output / "INPUTS.json",
        {"canonical": str(ROOT), "repository": str(repository), "files": before, "bytes": total},
    )
    return repository, before


def pytest_runtest_call(item: pytest.Item) -> None:
    destination = os.environ.get("MTP_1026_EVENTS")
    if destination is None:
        return
    root = Path(os.environ["MTP_1026_ROOT"])
    module = getattr(item, "module", None)
    expected = root / item.nodeid.split("::", 1)[0]
    if module is None or Path(module.__file__).resolve() != expected:
        raise AssertionError("Test origin escaped the #1026 isolated copy")
    conftest = sys.modules.get("conftest")
    if conftest is None or Path(str(conftest.__file__)).resolve() != root / CONFTEST:
        raise AssertionError("Conftest origin escaped the #1026 isolated copy")
    for name, imported in tuple(sys.modules.items()):
        if name == "mtp" or name.startswith(("mtp.", "bsdbng.")):
            filename = getattr(imported, "__file__", None)
            if filename and not Path(filename).resolve().is_relative_to(root / "src"):
                raise AssertionError(f"Package origin escaped the #1026 copy: {name}")
    producer = getattr(module, "infographic", None)
    if (
        producer is not None
        and Path(producer.write_stratified_trait_stats.__code__.co_filename).resolve()
        != root / PRODUCER
    ):
        raise AssertionError("Infographic callable origin escaped the #1026 copy")


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(
    item: pytest.Item, call: pytest.CallInfo[Any]
) -> Generator[None, Any, None]:
    outcome = yield
    report = outcome.get_result()
    destination = os.environ.get("MTP_1026_EVENTS")
    if destination is None:
        return
    error = call.excinfo
    record: dict[str, Any] = {
        "nodeid": item.nodeid,
        "phase": call.when,
        "outcome": report.outcome,
        "exception": error.typename if error else None,
        "sites": [
            {"path": str(entry.path), "name": entry.name, "statement": str(entry.statement)}
            for entry in error.traceback
        ]
        if error
        else [],
    }
    if call.when == "call":
        module = getattr(item, "module", None)
        producer = getattr(module, "infographic", None)
        record["origin"] = {
            "test": getattr(module, "__file__", None),
            "producer": getattr(producer, "__file__", None),
            "conftest": getattr(sys.modules.get("conftest"), "__file__", None),
        }
    if call.when == "call" and report.passed and os.environ["MTP_1026_LABEL"] == "control":
        artifacts = json.loads(os.environ["MTP_1026_ARTIFACTS"])
        if artifacts:
            temporary = getattr(item, "funcargs", {})["tmp_path"]
            evidence = Path(destination).parent / "control-artifacts"
            # The writer test deliberately replaces Path.mkdir with a write
            # barrier. This independent evidence channel retains the original.
            ORIGINAL_MKDIR(evidence)
            record["artifacts"] = {}
            for relative in artifacts:
                source = path_in(temporary, relative)
                raw = read(source)
                require(bool(raw), "Canary artifact is empty")
                target = path_in(evidence, relative)
                current = evidence
                for part in Path(relative).parts[:-1]:
                    current /= part
                    if not current.exists():
                        ORIGINAL_MKDIR(current)
                with target.open("xb") as handle:
                    handle.write(raw)
                require(read(target) == raw, "Canary artifact did not persist")
                record["artifacts"][relative] = {
                    "source": str(source),
                    "saved": str(target),
                    "sha256": sha(raw),
                    "size_bytes": len(raw),
                }
            record["captured_before_teardown"] = True
    with Path(destination).open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record) + "\n")


def run(
    executable: Path, repository: Path, output: Path, spec: dict[str, Any], label: str
) -> dict[str, Any]:
    temporary = output / label
    temporary.mkdir()
    events, log = output / f"{label}.events.jsonl", output / f"{label}.log"
    events.touch(exist_ok=False)
    env = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith(("PYTEST_", "PYTHONPATH", "GIT_"))
    }
    env.update(
        PYTHONPATH=os.pathsep.join(
            (str(repository / "src"), str(repository / "tests"), str(output))
        ),
        PYTHONDONTWRITEBYTECODE="1",
        PYTEST_DISABLE_PLUGIN_AUTOLOAD="1",
        TMPDIR=str(temporary),
        MTP_1026_EVENTS=str(events),
        MTP_1026_ROOT=str(repository),
        MTP_1026_LABEL=label,
        MTP_1026_ARTIFACTS=json.dumps(spec["artifacts"]),
    )
    command = [
        str(executable),
        "-q",
        "--tb=short",
        "-p",
        "probe_1026_plugin",
        "-p",
        "no:cacheprovider",
        f"--basetemp={temporary / 'pytest'}",
        spec["target"],
    ]
    with log.open("xb") as handle:
        original_cwd = Path.cwd()
        try:
            os.chdir(repository)
            pid = os.posix_spawn(
                str(executable),
                command,
                env,
                file_actions=[
                    (os.POSIX_SPAWN_DUP2, handle.fileno(), 1),
                    (os.POSIX_SPAWN_DUP2, handle.fileno(), 2),
                ],
            )
        finally:
            os.chdir(original_cwd)
        deadline = time.monotonic() + 180
        while True:
            waited, status = os.waitpid(pid, os.WNOHANG)
            if waited:
                code = os.waitstatus_to_exitcode(status)
                break
            if time.monotonic() >= deadline:
                os.kill(pid, signal.SIGKILL)
                os.waitpid(pid, 0)
                raise TimeoutError("Isolated regression exceeded 180 seconds")
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


def check_canary(output: Path, spec: dict[str, Any], control: dict[str, Any]) -> None:
    """Bind persisted before-teardown artifacts to the actual passing call."""
    if not spec["artifacts"]:
        return
    require(len(control["calls"]) == 1, "Artifact canary requires one selected call")
    call = control["calls"][0]
    require(
        call.get("captured_before_teardown") is True
        and set(call.get("artifacts", {})) == set(spec["artifacts"]),
        "Canary capture is missing",
    )
    for relative, binding in call["artifacts"].items():
        target = path_in(output / "control-artifacts", relative)
        source = Path(binding["source"])
        require(
            source.is_relative_to(output / "control" / "pytest")
            and binding["saved"] == str(target),
            "Canary artifact path escaped its run",
        )
        raw = read(target)
        require(
            raw == read(source)
            and len(raw) == binding["size_bytes"] > 0
            and sha(raw) == binding["sha256"],
            "Canary artifact bytes changed",
        )
    save(
        output / "CANARY.json",
        {
            "status": "CANARY_PASSED",
            "nodeid": call["nodeid"],
            "control_sha256": sha(read(output / "control.json")),
            "artifacts": call["artifacts"],
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mutation", required=True)
    parser.add_argument("--pytest", type=Path, default=Path(sys.executable).with_name("pytest"))
    args = parser.parse_args()
    specs = specifications()
    selected = [row for row in specs if row["name"] == args.mutation]
    require(len(selected) == 1, "Unknown/duplicate mutation")
    spec = selected[0]
    executable = args.pytest.absolute()
    require(
        executable.is_file() and not executable.is_symlink() and os.access(executable, os.X_OK),
        "An executable absolute pytest entry point is required",
    )
    output = args.output_dir.absolute()
    repository, before = copy_inputs(output, specs)
    helper = read(Path(__file__))
    (output / "probe_1026_plugin.py").write_bytes(helper)
    save(
        output / "TOOLS.json",
        {
            "helper_sha256": sha(helper),
            "pytest": str(executable),
            "pytest_sha256": sha(read(executable)),
        },
    )
    require(spec["source"] in before, "Mutation source is not bound")
    source = path_in(repository, spec["source"])
    original = read(source)
    require(
        spec["count"] > 0 and original.decode().count(spec["old"]) == spec["count"],
        "Stale mutation anchor",
    )
    changed = original.decode().replace(spec["old"], spec["new"]).encode()
    require(changed != original, "Mutation has no effect")
    result = {
        "mutation": spec["name"],
        "source": spec["source"],
        "target": spec["target"],
        "original_sha256": sha(original),
        "mutated_sha256": sha(changed),
        "killed": False,
    }
    try:
        control = run(executable, repository, output, spec, "control")
        require(
            control["phases_valid"]
            and control["exit"] == 0
            and all(row["outcome"] == "passed" for row in control["calls"]),
            "Control failed; no mutant executed",
        )
        check_canary(output, spec, control)
        source.write_bytes(changed)
        mutant = run(executable, repository, output, spec, "mutant")
        result["killed"] = (
            mutant["phases_valid"]
            and mutant["exit"] == 1
            and {row["nodeid"] for row in mutant["calls"]}
            == {row["nodeid"] for row in control["calls"]}
            and all(
                row["outcome"] == "failed"
                and row["exception"] == "AssertionError"
                and any(
                    Path(site["path"]).resolve() == repository / spec["target"].split("::", 1)[0]
                    and spec["assertion_site"] in site["statement"]
                    for site in row["sites"]
                )
                for row in mutant["calls"]
            )
        )
    finally:
        source.write_bytes(original)
        canonical = {name: sha(read(path_in(ROOT, name))) for name in before}
        restored = {name: sha(read(path_in(repository, name))) for name in before}
        save(output / "CANONICAL_AFTER.json", canonical)
        save(output / "RESTORED.json", restored)
        save(output / "result.json", result)
        require(canonical == restored == before, "Canonical/restoration hash check failed")
        require(read(Path(__file__)) == helper, "Verifier changed during execution")
    print(json.dumps(result, indent=2))
    return 0 if result["killed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
