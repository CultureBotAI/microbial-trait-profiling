"""Verify isolated committed-result README mutants (#799, 2026-11-10).

Inputs: the bounded source tree, declared synthetic tests, and a pytest executable.
Outputs: exclusive scratch copies, phase-level test evidence and restoration hashes.
Canonical sources are read only; no Git operation or scientific run is performed.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
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
TEST = "tests/test_cellulose_readme_numbers_799.py"
SOURCE = "data/analysis/cellulose_decr/README.md"
DATA = Path(SOURCE).parent
SOURCE_NAMES = (
    "cellulose_decr_summary.json",
    "cellulose_decr_provenance.json",
    "pairwise_jaccard_significance.tsv",
    "cellulose_decreased_taxa_catalog.tsv",
    "cellulose_decreased_taxa_catalog_tier1_2.tsv",
    "compartment_jaccard_and_direction_stats.tsv",
    "significant_results_summary.tsv",
    "condition_definition_sensitivity.tsv",
    "condition_definition_overlap.tsv",
    "compartment_direction_interaction.tsv",
    "condition_definition_interaction.tsv",
    "sensitivity_definitions_jaccard.tsv",
    "trait_convergence_benchmark.tsv",
    "profile_and_study_pair_jaccard.tsv",
    "label_source_corrections.tsv",
)
SUPPORT = tuple((DATA / name).as_posix() for name in SOURCE_NAMES)
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
    """Map names to target path, unique old/new anchors, test node and assertion site."""
    bound = "test_readme_quantitative_blocks_match_committed_sources"
    different = "assert _normal(value) == _normal(binding.text)"
    pair = "| PD vs. HIV | 38, 24 | 21 | 41 | 0.512 | 0.488 |"
    null = (
        "| U6: whole study signatures | 0.527 (0.420\u20130.634) "
        "| 0.19 | 18.3 (12\u201324) | 0.73 |"
    )
    scenario = (
        "| Synonym labels counted under their disease | 38 / 27 / 37 | 20 / 50 | 0.14 | 0.54 |"
    )
    core = "three_way_core\t3-way core (PD ∩ HIV ∩ CRC)\t\t\t38\t24\t17\t47"
    draws = ' "monte_carlo": {\n  "permutations": 100000,'
    return {
        "changed_overlap_count": (
            SOURCE,
            pair,
            pair.replace("| 21 |", "| 22 |"),
            bound,
            different,
        ),
        "changed_percent_tolerance": (
            SOURCE,
            "each set size within max(2, 10% of observed)",
            "each set size within max(2, 11% of observed)",
            bound,
            different,
        ),
        "changed_signature_null_probability": (
            SOURCE,
            null,
            null.replace("| 0.19 |", "| 0.001 |"),
            bound,
            different,
        ),
        "changed_z_threshold": (
            SOURCE,
            "Feces dec, profiles with cellulose z > 0 only",
            "Feces dec, profiles with cellulose z > 1 only",
            bound,
            different,
        ),
        "changed_scenario_core_ratio": (
            SOURCE,
            scenario,
            scenario.replace("20 / 50", "21 / 50"),
            bound,
            different,
        ),
        "wrong_binding_section": (
            SOURCE,
            "## 2. Is the overlap unusual?",
            "## 2. Altered population overlap",
            bound,
            "assert section == binding.section",
        ),
        "new_unbound_z_claim": (
            SOURCE,
            "## 10. Reproduce\n",
            "The unsupported analysis reports z = -3.2.\n\n## 10. Reproduce\n",
            "test_readme_has_no_unbound_quantitative_claims",
            "assert not _unbound_numbers(text)",
        ),
        "changed_source_core_count": (
            (DATA / "pairwise_jaccard_significance.tsv").as_posix(),
            core,
            core.replace("\t17\t47", "\t18\t47"),
            "test_catalog_and_overlap_counts_agree",
            'assert sum(row["is_3way_core"] == "True" for row in catalog)',
        ),
        "changed_source_draw_count": (
            (DATA / "cellulose_decr_summary.json").as_posix(),
            draws,
            draws.replace("100000", "200000"),
            bound,
            different,
        ),
    }


def emit_bindings(output: Path) -> None:
    """Serialize source-derived blocks for review; never rewrite the README."""
    before = inventory()
    helper = read(Path(__file__))
    spec = importlib.util.spec_from_file_location("readme_799_binding_emitter", ROOT / TEST)
    require(spec is not None and spec.loader is not None, "Cannot load the binding module")
    if spec is None or spec.loader is None:
        raise ValueError("Cannot load the binding module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    require(Path(str(module.__file__)).resolve() == ROOT / TEST, "Binding import escaped checkout")
    require(
        module.DATA == ROOT / DATA and module.README == ROOT / SOURCE,
        "Binding data paths escaped checkout",
    )
    require(tuple(module.SOURCES) == SOURCE_NAMES, "Binding source allowlist differs")
    require(
        Path(module._expected_blocks.__code__.co_filename).resolve() == ROOT / TEST,
        "Binding renderer escaped checkout",
    )
    bindings = module._expected_blocks()
    result = {
        name: {"section": value.section, "text": value.text} for name, value in bindings.items()
    }
    save(output / "bindings.json", result)
    require(inventory() == before, "Canonical inputs changed while rendering bindings")
    require(read(Path(__file__)) == helper, "Verifier changed while rendering bindings")
    save(
        output / "binding_receipt.json",
        {
            "input_sha256": before,
            "helper_sha256": sha(helper),
            "bindings_sha256": sha(read(output / "bindings.json")),
            "binding_count": len(result),
        },
    )
    print(json.dumps({"binding_count": len(result), "output": str(output / "bindings.json")}))


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_call(item: pytest.Item) -> None:
    if "MTP_799_ROOT" not in os.environ:
        return
    root = Path(os.environ["MTP_799_ROOT"])
    module = getattr(item, "module", None)
    if module is None:
        raise RuntimeError("Selected test has no module")
    require(
        Path(str(module.__file__)).resolve() == root / TEST,
        "Selected test import escaped the isolated copy",
    )
    code = getattr(getattr(item, "obj", None), "__code__", None)
    require(
        code is not None and Path(code.co_filename).resolve() == root / TEST,
        "Selected test callable escaped the isolated copy",
    )
    require(
        root / DATA == module.DATA and root / SOURCE == module.README,
        "Binding data paths escaped the isolated copy",
    )
    require(tuple(module.SOURCES) == SOURCE_NAMES, "Binding source allowlist differs")
    code = module._expected_blocks.__code__
    require(
        Path(code.co_filename).resolve() == root / TEST,
        "Binding renderer escaped the isolated copy",
    )
    paths = sys.modules.get("bsdbng.paths")
    require(paths is not None and root == paths.REPO_ROOT, "Repository root escaped the copy")


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(
    item: pytest.Item, call: pytest.CallInfo[Any]
) -> Generator[None, Any, None]:
    outcome = yield
    report = outcome.get_result()
    destination = os.environ.get("MTP_799_EVENTS")
    if destination is None:
        return
    error = call.excinfo
    selected_test = Path(os.environ["MTP_799_ROOT"]) / TEST
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
        MTP_799_ROOT=str(repository),
        MTP_799_EVENTS=str(events),
    )
    command = [
        str(executable),
        "-q",
        "--tb=short",
        "-p",
        "guard_799_probe_plugin",
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
    for candidate, (target, anchor, replacement, _node, _site) in specs.items():
        canonical_target = read(ROOT / target)
        require(sha(canonical_target) == before[target], "Mutation target changed during preflight")
        count = canonical_target.decode().count(anchor)
        require(
            bool(anchor) and count == 1 and anchor != replacement,
            f"Mutation preflight failed for {candidate}: expected one anchor, found {count}",
        )
        anchor_counts[candidate] = count
        target_hashes[target] = sha(canonical_target)
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
    with (output / "guard_799_probe_plugin.py").open("xb") as handle:
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
                    Path(site["path"]).resolve() == repository / TEST
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
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--mutation")
    modes.add_argument("--emit-bindings", action="store_true")
    parser.add_argument("--pytest", type=Path, default=Path(sys.executable).with_name("pytest"))
    args = parser.parse_args()
    output = args.output_dir.resolve()
    require(
        not output.resolve().is_relative_to(ROOT) and not ROOT.is_relative_to(output.resolve()),
        "Evidence output must be outside the canonical repository",
    )
    require(not output.exists() and not output.is_symlink(), "Evidence output must be exclusive")
    output.mkdir(parents=True)
    if args.emit_bindings:
        emit_bindings(output)
        return 0
    return 0 if verify_mutation(args.mutation, args.pytest.resolve(strict=True), output) else 1


if __name__ == "__main__":
    raise SystemExit(main())
