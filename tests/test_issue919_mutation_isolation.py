"""The #919 audit probes must never rewrite the tracked checkout (#1184)."""

from __future__ import annotations

import builtins
import hashlib
import importlib.util
import io
import json
import os
import shutil
import stat
import sys
from collections.abc import Callable
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import pytest
from git import Repo

import mtp
from mtp.paths import REPO_ROOT, file_sha256

PROBE = REPO_ROOT / "scripts/adhoc/issue919_audit_mutations.py"
SOURCE = REPO_ROOT / "src/mtp/label_support_audit.py"
TEST = REPO_ROOT / "tests/test_label_support_audit.py"
MODULE = "mtp.label_support_audit"


def load_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def forbid_checkout_writes(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Record and reject path-based writes, even when nested pytest catches the error.

    Guard both lexical entries and their resolved targets. Relative dir_fd writes
    are conservatively refused because their roots cannot be resolved portably.
    This observes the listed Python APIs, not arbitrary native code or subprocesses.
    """
    original_open, original_io_open, original_os_open = builtins.open, io.open, os.open
    protected = {REPO_ROOT.absolute(), REPO_ROOT.resolve()}
    attempts: list[str] = []

    def check(
        path: Any, writable: bool, *, operation: str = "write", dir_fd: int | None = None
    ) -> None:
        if not writable or isinstance(path, int):
            return
        candidate = Path(os.fsdecode(path))
        descriptor_relative = dir_fd is not None and not candidate.is_absolute()
        locations = {Path(os.path.abspath(candidate)), candidate.resolve()}
        if descriptor_relative or any(
            location.is_relative_to(root) for location in locations for root in protected
        ):
            reason = f"{operation}: {path}"
            if descriptor_relative:
                reason += " (unresolved dir_fd)"
            attempts.append(reason)
            raise PermissionError(f"Forbidden write inside tracked checkout: {reason}")

    def guarded_open(path: Any, mode: str = "r", *args: Any, **kwargs: Any) -> Any:
        check(path, any(flag in mode for flag in "wax+"))
        return original_open(path, mode, *args, **kwargs)

    def guarded_io_open(path: Any, mode: str = "r", *args: Any, **kwargs: Any) -> Any:
        check(path, any(flag in mode for flag in "wax+"))
        return original_io_open(path, mode, *args, **kwargs)

    def guarded_os_open(path: Any, flags: int, *args: Any, **kwargs: Any) -> int:
        check(
            path,
            bool(flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC)),
            dir_fd=kwargs.get("dir_fd"),
        )
        return original_os_open(path, flags, *args, **kwargs)

    def transfer(
        name: str, function: Callable[..., Any], *, mutates_source: bool
    ) -> Callable[..., Any]:
        def guarded(src: Any, dst: Any, *args: Any, **kwargs: Any) -> Any:
            check(src, mutates_source, operation=f"{name} source", dir_fd=kwargs.get("src_dir_fd"))
            check(dst, True, operation=f"{name} destination", dir_fd=kwargs.get("dst_dir_fd"))
            return function(src, dst, *args, **kwargs)

        return guarded

    monkeypatch.setattr(builtins, "open", guarded_open)
    monkeypatch.setattr(io, "open", guarded_io_open)
    monkeypatch.setattr(os, "open", guarded_os_open)
    for name in ("rename", "replace"):
        monkeypatch.setattr(os, name, transfer(name, getattr(os, name), mutates_source=True))
    for name in ("move", "copyfile", "copy", "copy2", "copytree", "copymode", "copystat"):
        monkeypatch.setattr(
            shutil, name, transfer(name, getattr(shutil, name), mutates_source=name == "move")
        )
    return attempts


@pytest.mark.parametrize("operation", ["rename", "replace", "path_rename", "path_replace", "move"])
@pytest.mark.parametrize("direction", ["into", "out_of", "symlink_entry"])
def test_guard_blocks_moves_at_either_checkout_endpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, operation: str, direction: str
) -> None:
    protected = tmp_path / "checkout"
    protected.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    source, destination = outside / "source", protected / "destination"
    if direction == "out_of":
        source, destination = protected / "source", outside / "destination"
    source.write_bytes(b"source")
    if direction == "symlink_entry":
        referent = outside / "referent"
        referent.write_bytes(b"referent")
        destination.symlink_to(referent)
    monkeypatch.setattr(sys.modules[__name__], "REPO_ROOT", protected)
    forbid_checkout_writes(monkeypatch)

    def move() -> None:
        if operation == "rename":
            os.rename(source, destination)
        elif operation == "replace":
            os.replace(source, destination)
        elif operation == "path_rename":
            source.rename(destination)
        elif operation == "path_replace":
            source.replace(destination)
        else:
            shutil.move(source, destination)

    with pytest.raises(PermissionError, match="Forbidden"):
        move()
    assert source.read_bytes() == b"source"
    if direction == "symlink_entry":
        assert destination.is_symlink()
        assert destination.read_bytes() == b"referent"
    else:
        assert not destination.exists()


@pytest.mark.parametrize(
    "operation", ["copyfile", "copy", "copy2", "copytree", "copymode", "copystat"]
)
def test_guard_blocks_copies_before_any_checkout_destination_is_created(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, operation: str
) -> None:
    protected = tmp_path / "checkout"
    protected.mkdir()
    source = tmp_path / "source"
    if operation == "copytree":
        source.mkdir()
        (source / "data").write_bytes(b"source")
    else:
        source.write_bytes(b"source")
    destination = protected / "destination"
    if operation in {"copymode", "copystat"}:
        destination.write_bytes(b"preserved")
    monkeypatch.setattr(sys.modules[__name__], "REPO_ROOT", protected)
    forbid_checkout_writes(monkeypatch)
    with pytest.raises(PermissionError, match="Forbidden"):
        getattr(shutil, operation)(source, destination)
    if operation in {"copymode", "copystat"}:
        assert destination.read_bytes() == b"preserved"
    else:
        assert not destination.exists()


@pytest.mark.parametrize("operation", ["rename", "replace"])
def test_guard_blocks_descriptor_relative_moves(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, operation: str
) -> None:
    protected = tmp_path / "checkout"
    protected.mkdir()
    source = tmp_path / "source"
    source.write_bytes(b"source")
    descriptor = os.open(protected, os.O_RDONLY)
    monkeypatch.setattr(sys.modules[__name__], "REPO_ROOT", protected)
    forbid_checkout_writes(monkeypatch)
    try:
        with pytest.raises(PermissionError, match="Forbidden"):
            getattr(os, operation)(source, "destination", dst_dir_fd=descriptor)
    finally:
        os.close(descriptor)
    assert source.read_bytes() == b"source"
    assert not (protected / "destination").exists()


@pytest.mark.parametrize("operation", ["rename", "copy2"])
def test_guard_blocks_keyword_transfer_destinations(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, operation: str
) -> None:
    protected = tmp_path / "checkout"
    protected.mkdir()
    source = tmp_path / "source"
    source.write_bytes(b"source")
    destination = protected / "destination"
    monkeypatch.setattr(sys.modules[__name__], "REPO_ROOT", protected)
    attempts = forbid_checkout_writes(monkeypatch)
    function = os.rename if operation == "rename" else shutil.copy2
    with pytest.raises(PermissionError, match="Forbidden"):
        function(src=source, dst=destination)
    assert attempts
    assert source.read_bytes() == b"source"
    assert not destination.exists()


def test_guard_allows_checkout_reads_and_scratch_moves(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    protected = tmp_path / "checkout"
    protected.mkdir()
    source = protected / "source"
    source.write_bytes(b"source")
    monkeypatch.setattr(sys.modules[__name__], "REPO_ROOT", protected)
    attempts = forbid_checkout_writes(monkeypatch)
    copied = tmp_path / "copied"
    shutil.copy2(src=source, dst=copied)
    renamed = tmp_path / "renamed"
    os.rename(src=copied, dst=renamed)
    shutil.move(src=renamed, dst=tmp_path / "moved")
    assert (tmp_path / "moved").read_bytes() == source.read_bytes()
    assert attempts == []


@pytest.mark.parametrize("mutation", ["all_input_taxa", "expanded_focal_family"])
def test_probe_never_writes_checkout_and_tests_import_the_mutant(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mutation: str
) -> None:
    monkeypatch.syspath_prepend(str(PROBE.parent))
    probe = load_module(PROBE, "audit_probe_under_test")
    monkeypatch.delitem(sys.modules, MODULE, raising=False)
    monkeypatch.delattr(mtp, "label_support_audit", raising=False)
    stale = ModuleType("stale_audit_tests")
    stale.__file__ = str(TEST)
    monkeypatch.setitem(sys.modules, stale.__name__, stale)
    before = {path: (file_sha256(path), path.stat().st_mtime_ns) for path in (SOURCE, TEST)}
    output = tmp_path / mutation
    monkeypatch.setattr(
        sys, "argv", [str(PROBE), "--mutation", mutation, "--output-dir", str(output)]
    )
    origins = []

    def fake_pytest(arguments: list[str], plugins: list[Any]) -> int:
        assert str(TEST) in arguments[-1]
        assert "no:cacheprovider" in arguments
        assert stale.__name__ not in sys.modules
        audit = sys.modules[MODULE]
        assert vars(mtp)["label_support_audit"] is audit
        origin = Path(str(audit.__file__)).resolve()
        assert origin.is_relative_to(output.resolve())
        assert origin != SOURCE.resolve()
        text = origin.read_text(encoding="utf-8")
        if mutation == "all_input_taxa":
            assert 'bound_taxa = set(occurrences["taxon_id"])' in text
        else:
            assert text.index("FOCAL_TRAITS = tuple(AUDIT_TRAITS)") > text.index("AUDIT_TRAITS =")
            assert len(audit.FOCAL_TRAITS) > 6
        fresh = load_module(TEST, "fresh_audit_tests")
        assert fresh.label_support_audit is audit
        for plugin in plugins:
            plugin.pytest_collection_modifyitems([SimpleNamespace(module=fresh)])
        origins.append(origin)
        print("Synthetic pytest result; actual mutant rejection is validated separately.")
        return 1

    monkeypatch.setattr(pytest, "main", fake_pytest)
    forbid_checkout_writes(monkeypatch)
    assert probe.main() == 0
    assert MODULE not in sys.modules
    assert not hasattr(mtp, "label_support_audit")
    assert sys.modules[stale.__name__] is stale
    assert before == {path: (file_sha256(path), path.stat().st_mtime_ns) for path in (SOURCE, TEST)}
    receipt = json.loads((output / f"{mutation}.json").read_text(encoding="utf-8"))
    assert receipt["loaded_module_origin"] == str(origins[0])
    assert receipt["original_source_sha256"] == receipt["restored_source_sha256"]
    assert receipt["canonical_unchanged"] is True
    assert (output / f"{mutation}.log").stat().st_size > 0


def tracked_tree_snapshot() -> dict[str, tuple[int, int, str]]:
    """Hash every index entry, retaining modes and symlink text without following links."""
    result = {}
    with Repo(REPO_ROOT, search_parent_directories=False) as repository:
        for (name, stage), entry in repository.index.entries.items():
            assert stage == 0, "The isolation test requires an unmerged-entry-free checkout"
            relative = os.fsdecode(name)
            path = REPO_ROOT / relative
            mode = path.lstat().st_mode
            if stat.S_ISLNK(mode):
                digest = hashlib.sha256(os.fsencode(os.readlink(path))).hexdigest()
            else:
                assert stat.S_ISREG(mode), f"Unexpected tracked file type: {relative}"
                digest = file_sha256(path)
            result[relative] = (entry.mode, mode, digest)
    assert SOURCE.relative_to(REPO_ROOT).as_posix() in result
    assert TEST.relative_to(REPO_ROOT).as_posix() in result
    return result


def test_probe_real_nested_pytest_never_writes_tracked_tree(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Observe the real mutant test through teardown, including swallowed write errors."""
    monkeypatch.syspath_prepend(str(PROBE.parent))
    monkeypatch.delitem(sys.modules, MODULE, raising=False)
    monkeypatch.delattr(mtp, "label_support_audit", raising=False)
    monkeypatch.setattr(sys, "dont_write_bytecode", True)
    output = tmp_path / "real-audit-probe"
    monkeypatch.setattr(
        sys, "argv", [str(PROBE), "--mutation", "all_input_taxa", "--output-dir", str(output)]
    )
    before = tracked_tree_snapshot()
    attempts = forbid_checkout_writes(monkeypatch)
    try:
        probe = load_module(PROBE, "audit_probe_real_pytest")
        exit_code = probe.main()
    finally:
        assert tracked_tree_snapshot() == before, "Tracked-tree hashes or modes changed"
    assert not attempts, f"Tracked-tree write attempts: {attempts}"
    assert exit_code == 0
    receipt = json.loads((output / "all_input_taxa.json").read_text(encoding="utf-8"))
    assert receipt["pytest_exit"] == 1
    assert receipt["killed"] is True
    assert Path(receipt["loaded_module_origin"]).is_relative_to(output)
    assert (output / "fixtures").is_dir()
    log = (output / "all_input_taxa.log").read_text(encoding="utf-8")
    assert 'assert cell["source_unknown_taxa"].eq(0).all()' in log
    assert "1 failed" in log
    assert "PermissionError" not in log


@pytest.mark.parametrize("present", [False, True])
@pytest.mark.parametrize("failure", ["none", "load", "body"])
def test_loader_restores_independent_bindings_even_on_exceptions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, present: bool, failure: str
) -> None:
    helper = load_module(PROBE.parent / "mutation_probe_isolation.py", "isolation_under_test")
    cached, attribute = ModuleType(MODULE), ModuleType(MODULE)
    if present:
        monkeypatch.setitem(sys.modules, MODULE, cached)
        monkeypatch.setattr(mtp, "label_support_audit", attribute, raising=False)
    else:
        monkeypatch.delitem(sys.modules, MODULE, raising=False)
        monkeypatch.delattr(mtp, "label_support_audit", raising=False)
    stale = ModuleType("previous_audit_tests")
    stale.__file__ = str(TEST)
    monkeypatch.setitem(sys.modules, stale.__name__, stale)
    bytecode_before = sys.dont_write_bytecode
    payload = b"raise RuntimeError('load failure')\n" if failure == "load" else b"VALUE = 42\n"
    context = helper.isolated_module(REPO_ROOT, MODULE, payload, tmp_path / "copy", [TEST])

    def exercise() -> None:
        with context as module:
            assert sys.modules[MODULE] is vars(mtp)["label_support_audit"] is module
            assert stale.__name__ not in sys.modules
            assert module.VALUE == 42
            fresh = ModuleType("fresh_audit_tests")
            fresh.__file__ = str(TEST)
            sys.modules[fresh.__name__] = fresh
            if failure == "body":
                raise RuntimeError("body failure")

    forbid_checkout_writes(monkeypatch)
    if failure == "none":
        exercise()
    else:
        with pytest.raises(RuntimeError, match=f"{failure} failure"):
            exercise()
    assert sys.dont_write_bytecode is bytecode_before
    assert sys.modules[stale.__name__] is stale
    assert "fresh_audit_tests" not in sys.modules
    if present:
        assert sys.modules[MODULE] is cached
        assert vars(mtp)["label_support_audit"] is attribute
    else:
        assert MODULE not in sys.modules
        assert not hasattr(mtp, "label_support_audit")


@pytest.mark.parametrize(
    "destination", ["inside", "ancestor", "symlink", "dangling_symlink", "existing"]
)
def test_loader_refuses_unsafe_or_reused_output_before_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, destination: str
) -> None:
    helper = load_module(PROBE.parent / "mutation_probe_isolation.py", "isolation_under_test")
    if destination == "inside":
        output = REPO_ROOT / "forbidden-mutation-copy"
    elif destination == "ancestor":
        output = REPO_ROOT.parent
    elif destination == "symlink":
        alias = tmp_path / "checkout-alias"
        alias.symlink_to(REPO_ROOT, target_is_directory=True)
        output = alias / "forbidden-mutation-copy"
    elif destination == "dangling_symlink":
        output = tmp_path / "dangling"
        output.symlink_to(tmp_path / "not-created", target_is_directory=True)
    else:
        output = tmp_path / "existing"
        output.mkdir()
        (output / SOURCE.name).write_bytes(b"preserved existing evidence")
    before = file_sha256(SOURCE), SOURCE.stat().st_mtime_ns
    forbid_checkout_writes(monkeypatch)
    with (
        pytest.raises((ValueError, FileExistsError), match=r"outside|exists"),
        helper.isolated_module(REPO_ROOT, MODULE, b"VALUE = 1\n", output, [TEST]),
    ):
        pytest.fail("Unsafe output was accepted")
    assert (file_sha256(SOURCE), SOURCE.stat().st_mtime_ns) == before
    if destination == "existing":
        assert (output / SOURCE.name).read_bytes() == b"preserved existing evidence"
    elif destination != "ancestor":
        assert not output.exists()


@pytest.mark.parametrize("foreign", ["mtp", "paths"])
def test_probe_rejects_foreign_checkout_before_output_creation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, foreign: str
) -> None:
    from mtp import paths

    monkeypatch.syspath_prepend(str(PROBE.parent))
    probe = load_module(PROBE, "audit_probe_foreign")
    if foreign == "mtp":
        monkeypatch.setattr(mtp, "__file__", str(tmp_path / "src/mtp/__init__.py"))
    else:
        monkeypatch.setattr(paths, "REPO_ROOT", tmp_path)
    output = tmp_path / "not-created"
    monkeypatch.setattr(
        sys, "argv", [str(PROBE), "--mutation", "all_input_taxa", "--output-dir", str(output)]
    )
    with pytest.raises(AssertionError, match="different checkout"):
        probe.main()
    assert not output.exists()


@pytest.mark.parametrize("problem", ["unverified", "wrong_module", "wrong_function", "crash"])
def test_probe_does_not_report_isolation_failures_as_mutant_kills(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, problem: str
) -> None:
    monkeypatch.syspath_prepend(str(PROBE.parent))
    probe = load_module(PROBE, "audit_probe_failure")
    monkeypatch.delitem(sys.modules, MODULE, raising=False)
    monkeypatch.delattr(mtp, "label_support_audit", raising=False)
    output = tmp_path / "failed"
    monkeypatch.setattr(
        sys, "argv", [str(PROBE), "--mutation", "all_input_taxa", "--output-dir", str(output)]
    )

    def fake_pytest(arguments: list[str], plugins: list[Any]) -> int:
        assert str(TEST) in arguments[-1]
        if problem == "crash":
            raise RuntimeError("pytest crashed")
        if problem != "unverified":
            audit = sys.modules[MODULE]
            fresh = SimpleNamespace(
                __file__=str(TEST),
                label_support_audit=object() if problem == "wrong_module" else audit,
                run_label_support_audit=object()
                if problem == "wrong_function"
                else audit.run_label_support_audit,
            )
            plugins[0].pytest_collection_modifyitems([SimpleNamespace(module=fresh)])
        return 1

    monkeypatch.setattr(pytest, "main", fake_pytest)
    forbid_checkout_writes(monkeypatch)
    with pytest.raises((AssertionError, RuntimeError), match=r"isolated|crashed"):
        probe.main()
    assert not (output / "all_input_taxa.json").exists()
    assert MODULE not in sys.modules
    assert not hasattr(mtp, "label_support_audit")
