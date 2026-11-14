"""Derive the pinned guard reference and verify isolated mutants (#803, 2026-11-10).

Inputs: this checkout's committed guard cases and existing version pins.
Outputs: exclusive evidence directories with source hashes and a reference
snapshot, or isolated mutation receipts. Canonical files are read only.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import pprint
import signal
import stat
import sys
import time
from collections import Counter
from collections.abc import Generator
from dataclasses import asdict
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
TEST = "tests/test_contrast_orientation.py"
SOURCE = "src/mtp/contrast_orientation.py"
PINS = "tests/fixtures/contrast_orientation_guard_digests.json"
REVIEW = "methods/contrast_orientation_review.tsv"
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


def inventory(*, support: bool = False) -> dict[str, str]:
    """Bound only the small Python implementation tree and declared guard inputs."""
    names = {TEST, PINS, REVIEW, "pyproject.toml", "uv.lock"}
    names.update(path.relative_to(ROOT).as_posix() for path in (ROOT / "src").rglob("*.py"))
    if support:
        names.update(
            {
                "tests/conftest.py",
                "methods/contrast_orientation.md",
                "methods/calibration_protocol.json",
            }
        )
        names.update(
            path.relative_to(ROOT).as_posix()
            for path in (ROOT / "schema").rglob("*")
            if path.is_file()
        )
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


def load_guard() -> ModuleType:
    path = ROOT / TEST
    spec = importlib.util.spec_from_file_location("guard_reference_803", path)
    require(spec is not None and spec.loader is not None, "Cannot load the guard module")
    if spec is None or spec.loader is None:
        raise ValueError("Cannot load the guard module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    require(
        Path(module.classify.__code__.co_filename).resolve() == ROOT / SOURCE,
        "Classifier import did not come from the requested checkout",
    )
    require(Path(str(module.__file__)).resolve() == path, "Guard import escaped the checkout")
    return module


def emit_reference(output: Path) -> None:
    """Emit only when both existing legacy hashes match, without repinning them."""
    before = inventory()
    module = load_guard()
    cases = module._guard_cases()
    version = module.CLASSIFIER_VERSION
    pins = json.loads(read(ROOT / PINS))
    found = module._guard_digests(cases)
    require(found == pins.get(version), "Guard cases or decisions differ from the existing pins")
    rows = []
    for record in cases:
        decision = module.classify(record)
        flags = ",".join(sorted(flag for flag in decision.flags.split(",") if flag))
        rows.append(
            {
                "record": asdict(record),
                "decision": [
                    decision.contrast_class,
                    decision.orientation,
                    decision.contrast_type,
                    decision.control_tier,
                    flags,
                ],
            }
        )
    require(inventory() == before, "Inputs changed while deriving the reference")
    save(output / "reference.json", rows)
    save(
        output / "reference_receipt.json",
        {
            "classifier_version": version,
            "case_count": len(rows),
            "legacy_digests": found,
            "reference_sha256": sha(read(output / "reference.json")),
            "inputs": before,
        },
    )
    print(
        json.dumps(
            {
                "classifier_version": version,
                "case_count": len(rows),
                "legacy_digests": found,
                "reference": str(output / "reference.json"),
            }
        )
    )


def format_reference(reference: Path, output: Path) -> None:
    """Compact an already qualified snapshot, without importing evolving test code."""
    raw = read(reference)
    receipt = json.loads(read(reference.with_name("reference_receipt.json")))
    require(sha(raw) == receipt["reference_sha256"], "Reference bytes differ from the receipt")
    rows = json.loads(raw)
    require(isinstance(rows, list) and bool(rows), "Reference must contain case rows")
    fields = set(rows[0]["record"])
    inputs, decisions = [], []
    for row in rows:
        record, decision = row["record"], row["decision"]
        require(
            isinstance(record, dict)
            and set(record) == fields
            and all(
                isinstance(key, str) and isinstance(value, str) for key, value in record.items()
            ),
            "Reference records have inconsistent fields",
        )
        require(
            isinstance(decision, list)
            and len(decision) == 5
            and all(isinstance(value, str) for value in decision),
            "Invalid reference decision",
        )
        case = json.dumps(record, ensure_ascii=False, sort_keys=True)
        inputs.append(case)
        decisions.append("\t".join([case, *decision]))
    digests = {
        "inputs": sha("\n".join(sorted(inputs)).encode()),
        "decisions": sha("\n".join(sorted(decisions)).encode()),
    }
    version = receipt["classifier_version"]
    pins = json.loads(read(ROOT / PINS))
    require(
        digests == receipt["legacy_digests"] == pins.get(version),
        "Reference does not reproduce both existing legacy pins",
    )
    require(len(rows) == receipt["case_count"], "Reference count differs from the receipt")
    common = {}
    for field in sorted(fields):
        counts = Counter(row["record"][field] for row in rows)
        common[field] = min(counts, key=lambda value: (-counts[value], value))
    catalog = tuple(sorted({tuple(row["decision"]) for row in rows}))
    cases = tuple(
        (
            {key: value for key, value in row["record"].items() if value != common[key]},
            catalog.index(tuple(row["decision"])),
        )
        for row in rows
    )
    expanded = [
        {"record": {**common, **changes}, "decision": list(catalog[index])}
        for changes, index in cases
    ]
    require(expanded == rows, "Compact snapshot changed the qualified rows")
    fragment = (
        "# Frozen from the existing committed guard cases only (#803); no new export rows.\n"
        f"_GUARD_REFERENCE_VERSION = {version!r}\n"
        "_GUARD_REFERENCE_COMMON: dict[str, str] = " + pprint.pformat(common, width=88) + "\n\n"
        "_GUARD_REFERENCE_DECISIONS: tuple[tuple[str, str, str, str, str], ...] = "
        + pprint.pformat(catalog, width=88)
        + "\n\n"
        "_GUARD_REFERENCE_CASES: tuple[tuple[dict[str, str], int], ...] = "
        + pprint.pformat(cases, width=88)
        + "\n\n"
        "_GUARD_REFERENCE_SNAPSHOT: list[dict[str, Any]] = [\n"
        "    {\n"
        '        "record": {**_GUARD_REFERENCE_COMMON, **changes},\n'
        '        "decision": list(_GUARD_REFERENCE_DECISIONS[index]),\n'
        "    }\n"
        "    for changes, index in _GUARD_REFERENCE_CASES\n"
        "]\n"
    )
    destination = output / "reference_fragment.txt"
    with destination.open("x", encoding="utf-8") as handle:
        handle.write(fragment)
    save(
        output / "format_receipt.json",
        {
            "classifier_version": version,
            "reference_sha256": sha(raw),
            "fragment_sha256": sha(read(destination)),
            "case_count": len(rows),
            "unique_decisions": len(catalog),
            "legacy_digests": digests,
            "helper_sha256": sha(read(Path(__file__))),
        },
    )
    print(
        json.dumps(
            {
                "classifier_version": version,
                "case_count": len(rows),
                "unique_decisions": len(catalog),
                "fragment": str(destination),
            }
        )
    )


def mutations() -> dict[str, tuple[str, str, str, str]]:
    """Return unique old/new anchors, selected node and intended assertion site."""
    duplicate_test = "test_guard_decision_rows_preserve_duplicate_keys_and_ignore_unused_fields"
    decision_prefix = (
        "def _guard_decision(record: ExperimentRecord) -> _GuardDecision:\n"
        "    found = classify(record)\n"
    )
    flags_anchor = decision_prefix + (
        '    flags = ",".join(sorted(flag for flag in found.flags.split(",") if flag))'
    )
    return {
        "skip_frozen_decisions": (
            "    decision_changes = _guard_changed_reference_decisions(reference)",
            "    decision_changes = []",
            "test_guard_reports_each_changed_classification_field",
            'assert "classification decision changes" in message',
        ),
        "skip_fixture_drift": (
            "    input_changes = _guard_input_differences(cases, reference)",
            "    input_changes = []",
            "test_guard_fixture_drift_does_not_request_a_classifier_version_bump",
            'assert "fixture-input drift" in message',
        ),
        "skip_reference_integrity": (
            "    assert _guard_reference_digests(reference) == pinned, (",
            "    assert True, (",
            "test_guard_rejects_a_corrupt_reference_before_diagnosing_current_rules",
            'assert "frozen reference" in message',
        ),
        "skip_reference_version": (
            "    assert reference_version == CLASSIFIER_VERSION, (",
            "    assert True, (",
            "test_guard_rejects_a_reference_from_a_different_classifier_version",
            'raise AssertionError("The classifier guard did not refuse the changed case")',
        ),
        "discard_duplicate_decisions": (
            "    return {key: tuple(sorted(values)) for key, values in rows.items()}",
            "    return {key: tuple(sorted(set(values))) for key, values in rows.items()}",
            duplicate_test,
            "assert rows == {key: (expected, expected)}",
        ),
        "discard_flags": (
            flags_anchor,
            decision_prefix + '    flags = ""',
            duplicate_test,
            "assert rows == {key: (expected, expected)}",
        ),
        "retain_unsorted_empty_flags": (
            flags_anchor,
            decision_prefix + "    flags = found.flags",
            duplicate_test,
            "assert rows == {key: (expected, expected)}",
        ),
    }


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_call(item: pytest.Item) -> None:
    if "MTP_803_ROOT" not in os.environ:
        return
    root = Path(os.environ["MTP_803_ROOT"])
    module = getattr(item, "module", None)
    require(
        module is not None and Path(str(module.__file__)).resolve() == root / TEST,
        "Selected test import escaped the isolated copy",
    )
    code = getattr(getattr(item, "obj", None), "__code__", None)
    require(
        code is not None and Path(code.co_filename).resolve() == root / TEST,
        "Selected test callable escaped the isolated copy",
    )
    classifier = sys.modules.get("mtp.contrast_orientation")
    require(
        classifier is not None
        and Path(str(classifier.__file__)).resolve() == root / SOURCE
        and Path(classifier.classify.__code__.co_filename).resolve() == root / SOURCE,
        "Classifier module or callable escaped the isolated copy",
    )
    for name in ("_assert_classifier_guard", "_guard_decision_rows"):
        function = getattr(module, name, None)
        code = getattr(function, "__code__", None)
        require(
            code is not None and Path(code.co_filename).resolve() == root / TEST,
            "Guard callable escaped the isolated copy",
        )
    paths = sys.modules.get("mtp.paths")
    require(paths is not None and root == paths.REPO_ROOT, "Repository root escaped the copy")


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(
    item: pytest.Item, call: pytest.CallInfo[Any]
) -> Generator[None, Any, None]:
    outcome = yield
    report = outcome.get_result()
    destination = os.environ.get("MTP_803_EVENTS")
    if destination is None:
        return
    error = call.excinfo
    selected_test = Path(os.environ["MTP_803_ROOT"]) / TEST
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
        MTP_803_ROOT=str(repository),
        MTP_803_EVENTS=str(events),
    )
    command = [
        str(executable),
        "-q",
        "--tb=short",
        "-p",
        "guard_803_probe_plugin",
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
    old, new, node, assertion_site = specs[name]
    require(
        executable.is_file() and not executable.is_symlink() and os.access(executable, os.X_OK),
        "An absolute pytest executable is required",
    )
    before = inventory(support=True)
    canonical_target = read(ROOT / TEST)
    require(sha(canonical_target) == before[TEST], "Mutation target changed during preflight")
    anchor_counts = {}
    for candidate, (anchor, replacement, _node, _site) in specs.items():
        count = canonical_target.decode().count(anchor)
        require(
            bool(anchor) and count == 1 and anchor != replacement,
            f"Mutation preflight failed for {candidate}: expected one anchor, found {count}",
        )
        anchor_counts[candidate] = count
    save(
        output / "ANCHOR_PREFLIGHT.json",
        {
            "target_sha256": sha(canonical_target),
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
    (output / "guard_803_probe_plugin.py").write_bytes(helper)
    save(output / "INPUTS.json", before)
    save(
        output / "TOOLS.json",
        {
            "helper_sha256": sha(helper),
            "pytest": str(executable),
            "pytest_sha256": sha(read(executable)),
        },
    )
    source = repository / TEST
    original = read(source)
    require(bool(old) and original.decode().count(old) == 1, "Mutation anchor is not unique")
    changed = original.decode().replace(old, new).encode()
    require(changed != original, "Mutation does not change source")
    result = {
        "mutation": name,
        "target": f"{TEST}::{node}",
        "killed": False,
        "original_sha256": sha(original),
        "mutated_sha256": sha(changed),
    }
    try:
        control = run_pytest(executable, repository, output, "control", f"{TEST}::{node}")
        require(
            control["exit"] == 0
            and control["phases_valid"]
            and all(row["outcome"] == "passed" for row in control["calls"]),
            "Positive control failed; mutant was not executed",
        )
        source.write_bytes(changed)
        mutant = run_pytest(executable, repository, output, "mutant", f"{TEST}::{node}")
        result["killed"] = (
            mutant["exit"] == 1
            and mutant["phases_valid"]
            and {row["nodeid"] for row in mutant["calls"]}
            == {row["nodeid"] for row in control["calls"]}
            and all(
                row["outcome"] == "failed"
                and row["exception"] == "AssertionError"
                and any(
                    Path(site["path"]).resolve() == source and assertion_site in site["statement"]
                    for site in row["sites"]
                )
                for row in mutant["calls"]
            )
        )
    finally:
        source.write_bytes(original)
        canonical = inventory(support=True)
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
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--emit-reference", action="store_true")
    modes.add_argument("--format-reference", type=Path)
    modes.add_argument("--mutation")
    parser.add_argument("--pytest", type=Path, default=Path(sys.executable).with_name("pytest"))
    args = parser.parse_args()
    output = args.output_dir.resolve()
    require(
        not output.resolve().is_relative_to(ROOT) and not ROOT.is_relative_to(output.resolve()),
        "Evidence output must be outside the canonical repository",
    )
    require(not output.exists() and not output.is_symlink(), "Evidence output must be exclusive")
    output.mkdir(parents=True)
    if args.emit_reference:
        emit_reference(output)
    elif args.format_reference is not None:
        format_reference(args.format_reference, output)
    else:
        return 0 if verify_mutation(args.mutation, args.pytest.resolve(strict=True), output) else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
