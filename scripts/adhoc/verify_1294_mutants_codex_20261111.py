"""Verify #1294 with bounded scratch copies and native pytest failures (2026-11-11).

Inputs: four declared source/test files and an external --scratch-dir.
Outputs: retained copies, input hashes, logs, JUnit, and JSON receipts. No scientific
data are copied. --canary selects one mutant; --all and the default select all.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import signal
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
TEST = Path("tests/test_count_input_fast_path_1294.py")
PRODUCER = Path("src/mtp/count_input.py")
HELPER = Path("scripts/adhoc/verify_1294_mutants_codex_20261111.py")
INPUTS = (Path("src/mtp/__init__.py"), PRODUCER, TEST, HELPER)
MAX_BYTES = 2 * 1024 * 1024
TIMEOUT_SECONDS = 90
PLUGIN = "verify_1294_mutants_codex_20261111"


def validate_scratch_parent(repo_root: Path, scratch_parent: Path) -> Path:
    """Refuse lexical and resolved containment without creating anything."""
    lexical_root = Path(os.path.abspath(repo_root.expanduser()))
    real_root = lexical_root.resolve()
    lexical_parent = Path(os.path.abspath(scratch_parent.expanduser()))
    real_parent = lexical_parent.resolve()
    if (
        lexical_parent.is_relative_to(lexical_root)
        or lexical_parent.is_relative_to(real_root)
        or real_parent.is_relative_to(real_root)
    ):
        raise ValueError("Scratch must be outside the repository, lexically and after resolution")
    return real_parent


def _inventory(root: Path) -> dict[str, bytes]:
    contents: dict[str, bytes] = {}
    total = 0
    for relative in INPUTS:
        path = root / relative
        if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
            raise ValueError(f"Input is not a regular in-repository source: {relative}")
        size = path.stat().st_size
        if size > MAX_BYTES or total + size > MAX_BYTES:
            raise ValueError("The fixed four-file input inventory exceeds its byte cap")
        with path.open("rb") as stream:
            data = stream.read(MAX_BYTES + 1)
        total += len(data)
        if len(data) != size or total > MAX_BYTES:
            raise ValueError("An input changed size or exceeded the byte cap while reading")
        contents[relative.as_posix()] = data
    return contents


def _hashes(contents: dict[str, bytes]) -> dict[str, str]:
    return {name: hashlib.sha256(value).hexdigest() for name, value in contents.items()}


def _write_json(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True) + "\n")


def classify_pytest_result(
    exit_code: int,
    call_reports: list[dict[str, object]],
    non_call_failures: list[dict[str, str]],
) -> bool:
    """Only one bound test's native call assertion is a valid mutation kill."""
    if exit_code != int(pytest.ExitCode.TESTS_FAILED) or non_call_failures:
        return False
    if len(call_reports) != 1:
        return False
    report = call_reports[0]
    return (
        report.get("outcome") == "failed"
        and report.get("native_assertion") is True
        and report.get("test_origin_ok") is True
        and report.get("assertion_in_test") is True
    )


class ProofPlugin:
    """Record actual collection origin, call exceptions and non-call failures."""

    def __init__(self, root: Path, output: Path, expected_node: str) -> None:
        self.root = root
        self.output = output
        self.expected_node = expected_node
        self.origins: dict[str, bool] = {}
        self.calls: list[dict[str, object]] = []
        self.exceptions: dict[str, dict[str, object]] = {}
        self.non_call_failures: list[dict[str, str]] = []

    def pytest_runtest_setup(self, item: pytest.Item) -> None:
        expected = (self.root / self.expected_node.split("::", 1)[0]).resolve()
        module = getattr(item, "module", None)
        function = getattr(item, "obj", None)
        code = getattr(function, "__code__", None)
        if module is None or code is None:
            raise RuntimeError("Selected test has no inspectable module/function origin")
        valid = (
            item.nodeid == self.expected_node
            and Path(module.__file__).resolve() == expected
            and Path(code.co_filename).resolve() == expected
        )
        if self.expected_node.split("::", 1)[0] == TEST.as_posix():
            producer = sys.modules.get("mtp.count_input")
            if producer is None:
                raise RuntimeError("The selected test did not import the copied producer")
            source = (self.root / PRODUCER).resolve()
            valid = (
                valid
                and Path(producer.__file__ or "").resolve() == source
                and Path(producer.CountTableReader.__iter__.__code__.co_filename).resolve()
                == source
            )
        self.origins[item.nodeid] = valid
        if not valid:
            raise RuntimeError("Selected test or producer origin escaped the scratch copy")

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call":
            self.calls.append(
                {
                    "nodeid": report.nodeid,
                    "outcome": report.outcome,
                    "test_origin_ok": self.origins.get(report.nodeid, False),
                }
            )
        elif report.failed:
            self.non_call_failures.append({"phase": report.when, "nodeid": report.nodeid})

    def pytest_exception_interact(
        self, node: object, call: pytest.CallInfo[Any], report: pytest.TestReport
    ) -> None:
        _ = node
        info = call.excinfo
        if info is None:
            return
        if report.when != "call":
            self.non_call_failures.append(
                {"phase": report.when, "exception": info.type.__name__, "detail": str(info.value)}
            )
            return
        expected = (self.root / self.expected_node.split("::", 1)[0]).resolve()
        sites = [
            {"line": entry.lineno + 1, "statement": str(entry.statement)}
            for entry in info.traceback
            if Path(str(entry.path)).resolve() == expected
        ]
        self.exceptions[report.nodeid] = {
            "native_assertion": info.type is AssertionError,
            "assertion_in_test": (
                bool(sites) and Path(str(info.traceback[-1].path)).resolve() == expected
            ),
            "exception": info.type.__name__,
            "exception_detail": str(info.value),
            "sites": sites,
        }

    def pytest_sessionfinish(self, session: pytest.Session, exitstatus: int) -> None:
        _ = session
        calls = [{**row, **self.exceptions.get(str(row["nodeid"]), {})} for row in self.calls]
        _write_json(
            self.output,
            {
                "exit_code": int(exitstatus),
                "call_reports": calls,
                "non_call_failures": self.non_call_failures,
            },
        )


def pytest_configure(config: pytest.Config) -> None:
    """Activated only when this helper is explicitly loaded as the child plugin."""
    root = os.environ.get("COUNT1294_PROOF_ROOT")
    output = os.environ.get("COUNT1294_PROOF_OUTPUT")
    node = os.environ.get("COUNT1294_PROOF_NODE")
    if root and output and node:
        config.pluginmanager.register(ProofPlugin(Path(root), Path(output), node), "proof1294")


def _spawn_pytest(argv: list[str], environment: dict[str, str], log: Path) -> int:
    actions = [
        (os.POSIX_SPAWN_OPEN, 1, str(log), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600),
        (os.POSIX_SPAWN_DUP2, 1, 2),
    ]
    pid = os.posix_spawn(argv[0], argv, environment, file_actions=actions, setpgroup=0)
    deadline = time.monotonic() + TIMEOUT_SECONDS
    while time.monotonic() < deadline:
        finished, status = os.waitpid(pid, os.WNOHANG)
        if finished:
            return os.waitstatus_to_exitcode(status)
        time.sleep(0.05)
    os.killpg(pid, signal.SIGKILL)
    os.waitpid(pid, 0)
    return 124


def run_targeted_pytest(
    scratch_repo: Path,
    nodeid: str,
    junit_xml: Path,
    *,
    repo_root: Path | None = None,
) -> dict[str, Any]:
    """Run the actual pytest executable; retain both valid and invalid receipts."""
    canonical = ROOT if repo_root is None else repo_root
    scratch_repo = validate_scratch_parent(canonical, scratch_repo)
    validate_scratch_parent(canonical, junit_xml.parent)
    if junit_xml.is_symlink():
        raise ValueError("Evidence file must not be a symlink")
    junit_xml = validate_scratch_parent(canonical, junit_xml)
    before = _inventory(canonical)
    log = junit_xml.with_suffix(".log")
    proof = junit_xml.with_suffix(".proof.json")
    receipt_path = junit_xml.with_suffix(".receipt.json")
    plugin_dir = junit_xml.with_suffix(".plugin")
    temporary = junit_xml.with_suffix(".temporary")
    for path in (junit_xml, log, proof, receipt_path, plugin_dir, temporary):
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Refusing to overwrite evidence: {path}")
    junit_xml.parent.mkdir(parents=True, exist_ok=True)
    plugin_dir.mkdir()
    (plugin_dir / HELPER.name).write_bytes(before[HELPER.as_posix()])
    executable = Path(
        os.environ.get("COUNT1294_PYTEST", str(Path(sys.executable).parent / "pytest"))
    )
    argv = [
        str(executable.absolute()),
        "-q",
        "-p",
        PLUGIN,
        "-p",
        "no:cacheprovider",
        "-o",
        "addopts=",
        "--confcutdir",
        str(scratch_repo),
        "--rootdir",
        str(scratch_repo),
        "--basetemp",
        str(temporary),
        f"--junitxml={junit_xml}",
        str(scratch_repo / nodeid.split("::", 1)[0]) + "::" + nodeid.split("::", 1)[1],
    ]
    environment = dict(os.environ)
    environment.update(
        {
            "PYTHONPATH": os.pathsep.join((str(plugin_dir), str(scratch_repo / "src"))),
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1",
            "PYTEST_ADDOPTS": "",
            "COUNT1294_PROOF_ROOT": str(scratch_repo),
            "COUNT1294_PROOF_OUTPUT": str(proof),
            "COUNT1294_PROOF_NODE": nodeid,
        }
    )
    environment.pop("PYTEST_PLUGINS", None)
    calls: list[dict[str, object]] = []
    failures: list[dict[str, str]] = []
    exit_code = -1
    error = ""
    try:
        exit_code = _spawn_pytest(argv, environment, log)
        if proof.is_file():
            recorded = json.loads(proof.read_text(encoding="utf-8"))
            calls = recorded["call_reports"]
            failures = recorded["non_call_failures"]
            if recorded["exit_code"] != exit_code:
                error = "Pytest hook exit and process exit differ"
        else:
            error = "No native pytest proof was written"
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
    finally:
        unchanged = _inventory(canonical) == before
    killed = (
        unchanged
        and not error
        and junit_xml.is_file()
        and junit_xml.stat().st_size > 0
        and classify_pytest_result(exit_code, calls, failures)
    )
    detail = (
        f"exit={exit_code}; calls={calls!r}; non_call_failures={failures!r}; "
        f"canonical_unchanged={unchanged}; error={error or 'none'}"
    )
    receipt: dict[str, Any] = {
        "killed": killed,
        "status": "KILLED" if killed else "SURVIVED_OR_INVALID",
        "exit_code": exit_code,
        "detail": detail,
        "canonical_unchanged": unchanged,
        "call_reports": calls,
        "non_call_failures": failures,
        "input_sha256": _hashes(before),
        "argv": argv,
    }
    _write_json(receipt_path, receipt)
    return receipt


@dataclass(frozen=True)
class Mutation:
    name: str
    old: str
    new: str
    node: str


def mutations() -> tuple[Mutation, ...]:
    """Seven production faults, each selected by an independent test assertion."""
    node = TEST.as_posix() + "::"
    mode_guard = "if self.mode == COMPLETE_COUNT_MODE and all("
    return (
        Mutation(
            "bypass_complete_fast_path",
            mode_guard,
            "if False and all(",
            node + "test_complete_valid_rows_bypass_parser",
        ),
        Mutation(
            "accept_non_ascii_decimal",
            "cell.isascii() and cell.isdecimal() for cell in stripped_counts",
            "cell.isdecimal() for cell in stripped_counts",
            node + r"test_complete_invalid_cells_preserve_exact_parser_error[0-\u0662]",
        ),
        Mutation(
            "omit_count_whitespace_stripping",
            "stripped_counts = [values[position].strip() for position in self._count_positions]",
            "stripped_counts = [values[position] for position in self._count_positions]",
            node + "test_complete_valid_rows_bypass_parser",
        ),
        Mutation(
            "fast_path_in_legacy_mode",
            mode_guard,
            "if all(",
            node + "test_legacy_always_uses_per_cell_parser",
        ),
        Mutation(
            "increment_fast_path_counts",
            "counts = tuple(map(int, stripped_counts))",
            "counts = tuple(int(cell) + 1 for cell in stripped_counts)",
            node + "test_complete_rows_match_per_cell_reference",
        ),
        Mutation(
            "omit_legacy_blank_accounting",
            "self._blank_counts += 1",
            "self._blank_counts += 0",
            node + "test_legacy_always_uses_per_cell_parser",
        ),
        Mutation(
            "reverse_fallback_column_names",
            "self.count_columns, self._count_positions, strict=True",
            "reversed(self.count_columns), self._count_positions, strict=True",
            node + "test_declared_count_order_preserves_first_error",
        ),
    )


def _stage(root: Path, contents: dict[str, bytes]) -> None:
    root.mkdir()
    for relative, data in contents.items():
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as stream:
            stream.write(data)


def verify_mutants(
    repo_root: Path, scratch_parent: Path, *, canary_only: bool = False
) -> tuple[Path, list[dict[str, Any]]]:
    """Run controls and mutations without writing beneath the canonical repository."""
    parent = validate_scratch_parent(repo_root, scratch_parent)
    repo_root = repo_root.resolve()
    before = _inventory(repo_root)
    catalog = mutations()
    if len(catalog) < 6:
        raise RuntimeError("At least six declared mutations are required")
    original = before[PRODUCER.as_posix()].decode("utf-8")
    for mutant in catalog:
        if not mutant.old or original.count(mutant.old) != 1 or mutant.old == mutant.new:
            raise RuntimeError(f"Mutation anchor is not unique and changing: {mutant.name}")
        ast.parse(original.replace(mutant.old, mutant.new, 1))
    parent.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix="verify1294_", dir=parent)).resolve()
    _write_json(run / "inputs.json", _hashes(before))
    verdicts: list[dict[str, Any]] = []
    try:
        for mutant in catalog[:1] if canary_only else catalog:
            directory = run / mutant.name
            directory.mkdir()
            copied = directory / "repo"
            _stage(copied, before)
            control = run_targeted_pytest(
                copied, mutant.node, directory / "control.xml", repo_root=repo_root
            )
            if (
                control["exit_code"] != 0
                or control["canonical_unchanged"] is not True
                or len(control["call_reports"]) != 1
                or control["call_reports"][0].get("outcome") != "passed"
                or control["call_reports"][0].get("test_origin_ok") is not True
                or control["non_call_failures"]
            ):
                raise RuntimeError(f"Positive control failed: {mutant.name}: {control['detail']}")
            target = copied / PRODUCER
            try:
                target.write_text(original.replace(mutant.old, mutant.new, 1), encoding="utf-8")
                result = run_targeted_pytest(
                    copied, mutant.node, directory / "mutant.xml", repo_root=repo_root
                )
            finally:
                target.write_bytes(before[PRODUCER.as_posix()])
            restored = _inventory(copied) == before
            verdicts.append({"mutation": mutant.name, **result, "scratch_restored": restored})
            if result["killed"] is not True or not restored:
                raise RuntimeError(f"Mutant survived or evidence is invalid: {mutant.name}")
    finally:
        unchanged = _inventory(repo_root) == before
        _write_json(
            run / "verdicts.json",
            {"canonical_unchanged": unchanged, "verdicts": verdicts},
        )
        if not unchanged:
            raise RuntimeError("Canonical source bytes changed during verification")
    return run / "verdicts.json", verdicts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch-dir", required=True, type=Path)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--canary", action="store_true")
    selection.add_argument("--all", action="store_true")
    args = parser.parse_args()
    path, verdicts = verify_mutants(ROOT, args.scratch_dir, canary_only=args.canary)
    print(f"Verified {len(verdicts)} native assertion kills; retained evidence: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
