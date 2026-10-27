"""Run identity, non-destructive output and replay guards for #1044."""

from __future__ import annotations

import ctypes
import errno
import json
import os
import stat
import sys
import tempfile
import uuid
from collections.abc import Iterator
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import Mock

import pandas as pd
import pytest
from git import Actor, Repo

from mtp.datamodel.multiple_testing_pydantic import PostArmMultiplicityRun
from mtp.datamodel.workflow_pydantic import (
    ArmRun,
    ComparisonRun,
    ContrastPopulation,
    DirectionBasis,
    ExecutionStatus,
    WorkflowRun,
)
from mtp.multiplicity.sources import VerifiedSources
from mtp.multiplicity.stage import _verify_post_arm_from_sources, run_post_arm
from mtp.paths import file_sha256
from mtp.three_test_table import stage
from mtp.three_test_table.carriers import carrier_table as carrier_reader
from mtp.three_test_table.native import native_tables as native_reader
from multiplicity_fixture import fixture_sources, write_table


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    """Keep positive-stage artifacts durable independently of pytest's temporary root (#1058)."""
    root = Path(__file__).resolve().parents[1] / ".tmp" / "three-test-stage-tests"
    stage.require_durable(root).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as temporary:
        path = Path(temporary)
        with stage.durable_scratch(path / "default-scratch"):
            yield path


def sources_for_stage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> VerifiedSources:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    sources.core.git_commit = "a" * 40
    sources.comparison.id = "fixture-comparison"
    sources.comparison.git_commit = "b" * 40
    sources.comparison.status = ExecutionStatus.completed
    sources.comparison.started_at = datetime(2026, 10, 25, tzinfo=UTC)
    sources.comparison.git_dirty = False
    sources.comparison.analysis_specs = []
    sources.comparison.eligibility_ledgers = []
    sources.comparison.calibration_verdicts = []
    checked = ComparisonRun.model_validate(sources.comparison.model_dump())
    Path(checked.manifest_path).write_text(checked.model_dump_json())
    for name, arm in sources.arms.items():
        arm.id = f"fixture-{name}"
        arm.git_commit = "c" * 40
    monkeypatch.setattr(stage, "open_verified_sources", lambda _path: sources)
    return sources


@pytest.mark.parametrize("error_number", [errno.ENOTSUP, errno.ENOSYS])
def test_rename_capability_failure_precedes_source_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, error_number: int
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    forbidden = Mock(side_effect=AssertionError("source verification ran before capability check"))
    monkeypatch.setattr(stage, "open_verified_sources", forbidden)
    monkeypatch.setattr(
        stage, "_publish_directory", Mock(side_effect=OSError(error_number, "unsupported rename"))
    )
    with pytest.raises(OSError, match="unsupported rename"):
        stage.run_table(
            Path(sources.comparison.manifest_path),
            tmp_path / "output",
            scratch=tmp_path / "scratch",
        )
    forbidden.assert_not_called()
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".output-probe-*"))


@pytest.mark.parametrize("mutation", ["overwrite", "no_op", "always_exists"])
def test_capability_probe_checks_no_overwrite_and_real_move(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mutation: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    forbidden = Mock(side_effect=AssertionError("source verification ran before capability check"))
    monkeypatch.setattr(stage, "open_verified_sources", forbidden)

    def wrong_rename(source: Path, destination: Path) -> None:
        if mutation == "overwrite":
            source.rename(destination)
        elif mutation == "always_exists" or destination.exists():
            raise FileExistsError(errno.EEXIST, "probe always exists")

    monkeypatch.setattr(stage, "_publish_directory", wrong_rename)
    with pytest.raises((OSError, RuntimeError), match=r"(probe|exclusive rename)"):
        stage.run_table(
            Path(sources.comparison.manifest_path),
            tmp_path / "output",
            scratch=tmp_path / "scratch",
        )
    forbidden.assert_not_called()
    assert not (tmp_path / "output").exists()
    assert not list(tmp_path.glob(".output-probe-*"))


def test_successful_probe_uses_only_private_siblings_and_cleans_before_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    original = stage._publish_directory
    calls: list[tuple[Path, Path]] = []

    def inspect(source: Path, destination: Path) -> None:
        assert source.parent == destination.parent
        assert source.parent.parent == tmp_path
        assert stat.S_IMODE(source.parent.stat().st_mode) == 0o700
        calls.append((source, destination))
        original(source, destination)

    def stop_after_probe(_comparison: Path) -> VerifiedSources:
        assert len(calls) == 2
        assert not list(tmp_path.glob(".output-probe-*"))
        raise ValueError("reached verification after successful probe")

    monkeypatch.setattr(stage, "_publish_directory", inspect)
    monkeypatch.setattr(stage, "open_verified_sources", stop_after_probe)
    with pytest.raises(ValueError, match="reached verification after successful probe"):
        stage.run_table(
            Path(sources.comparison.manifest_path),
            tmp_path / "output",
            scratch=tmp_path / "scratch",
        )
    assert not (tmp_path / "output").exists()


def test_staging_mkdir_failure_does_not_publish_or_remove_unowned_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    output = tmp_path / "output"
    original = Path.mkdir

    def reject_staging(
        path: Path, mode: int = 0o777, parents: bool = False, exist_ok: bool = False
    ) -> None:
        if path.name.startswith(".output-staging-"):
            raise PermissionError("injected staging mkdir failure")
        original(path, mode, parents, exist_ok)

    monkeypatch.setattr(Path, "mkdir", reject_staging)
    with pytest.raises(PermissionError, match="injected staging mkdir failure"):
        stage.run_table(
            Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
        )
    assert not output.exists()
    assert not list(tmp_path.glob(".output-staging-*"))


@pytest.mark.parametrize("collision", ["directory", "symlink", "dangling_symlink"])
def test_staging_mkdir_collision_preserves_existing_directory_or_symlink(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, collision: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    output = tmp_path / "output"
    staging = tmp_path / ".output-staging-collision"
    target = staging if collision == "directory" else tmp_path / "other-directory"
    if collision != "dangling_symlink":
        target.mkdir()
        (target / "keep.txt").write_text("unowned sentinel", encoding="utf-8")
    if collision != "directory":
        staging.symlink_to(target, target_is_directory=True)
    original = staging.lstat()
    monkeypatch.setattr(uuid, "uuid4", lambda: SimpleNamespace(hex="collision"))
    with pytest.raises(FileExistsError):
        stage.run_table(
            Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
        )
    assert not output.exists()
    after = staging.lstat()
    assert (after.st_dev, after.st_ino, after.st_mode) == (
        original.st_dev,
        original.st_ino,
        original.st_mode,
    )
    assert staging.is_symlink() == (collision != "directory")
    if collision == "dangling_symlink":
        assert not target.exists()
    else:
        assert (target / "keep.txt").read_text(encoding="utf-8") == "unowned sentinel"
    assert not list(tmp_path.glob(".output-probe-*"))


def test_successful_publication_does_not_cleanup_recreated_staging_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    output = tmp_path / "output"
    publish = stage._publish_directory
    replacement: list[Path] = []

    def publish_then_replace(staging: Path, destination: Path) -> None:
        publish(staging, destination)
        if destination == output:
            staging.mkdir()
            (staging / "keep.txt").write_text(
                "replacement owned by another process", encoding="utf-8"
            )
            replacement.append(staging)

    monkeypatch.setattr(stage, "_publish_directory", publish_then_replace)
    record = stage.run_table(
        Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
    )
    assert len(replacement) == 1
    assert (replacement[0] / "keep.txt").read_text(
        encoding="utf-8"
    ) == "replacement owned by another process"
    assert {path.name for path in output.iterdir()} == {stage.MANIFEST, *record["outputs"]}
    assert not list(tmp_path.glob(".output-probe-*"))


@pytest.mark.parametrize("location", ["core", "comparison", "noel", "jed", "ilyas"])
def test_output_source_overlap_is_rejected_before_any_probe_or_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, location: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    protected = (
        Path(sources.core.manifest_path).parent
        if location == "core"
        else comparison.parent
        if location == "comparison"
        else Path(sources.arms[location].manifest_path).parent
    )
    before = set(protected.rglob("*"))
    probe = Mock(side_effect=AssertionError("probe ran before source overlap checks"))
    verification = Mock(side_effect=AssertionError("verification ran before source overlap checks"))
    monkeypatch.setattr(stage, "_probe_publication", probe)
    monkeypatch.setattr(stage, "open_verified_sources", verification)
    with pytest.raises(ValueError, match="disjoint"):
        stage.run_table(comparison, protected / "output", scratch=tmp_path / "scratch")
    probe.assert_not_called()
    verification.assert_not_called()
    assert set(protected.rglob("*")) == before


@pytest.mark.parametrize(
    "role", ["comparison", "core", "noel", "jed", "ilyas", "multiplicity", "label_support"]
)
@pytest.mark.parametrize("destination", ["output", "scratch"])
@pytest.mark.parametrize("protected_spelling", ["target", "alias"])
def test_manifest_file_aliases_protect_both_locations_before_any_write(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    role: str,
    destination: str,
    protected_spelling: str,
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    record = json.loads(comparison.read_text())
    options: dict[str, Path] = {}
    if role == "comparison":
        original = comparison
    elif role == "core":
        original = Path(sources.core.manifest_path)
    elif role in sources.arms:
        original = Path(sources.arms[role].manifest_path)
    else:
        original = tmp_path / "optional-evidence" / role / "manifest.json"
        original.parent.mkdir(parents=True)
        original.write_text("{}\n")
    alias = tmp_path / "manifest-aliases" / f"{role}.json"
    alias.parent.mkdir()
    alias.symlink_to(original)
    if role == "comparison":
        comparison = alias
    elif role in {"core", "noel", "jed", "ilyas"}:
        record[f"{role}_manifest_path"] = str(alias)
        comparison.write_text(json.dumps(record))
    else:
        options[role] = alias
    protected = original.parent if protected_spelling == "target" else alias.parent
    before = set(protected.rglob("*"))
    nested = protected / "must-not-be-created" / destination
    output = nested if destination == "output" else tmp_path / "output"
    scratch = nested if destination == "scratch" else tmp_path / "scratch"
    probe = Mock(side_effect=AssertionError("probe ran before alias overlap checks"))
    verification = Mock(side_effect=AssertionError("verification ran before alias overlap checks"))
    monkeypatch.setattr(stage, "_probe_publication", probe)
    monkeypatch.setattr(stage, "open_verified_sources", verification)
    with pytest.raises(ValueError, match="disjoint"):
        stage.run_table(
            comparison,
            output,
            scratch=scratch,
            multiplicity=options.get("multiplicity"),
            label_support=options.get("label_support"),
        )
    probe.assert_not_called()
    verification.assert_not_called()
    assert not nested.parent.exists()
    assert set(protected.rglob("*")) == before


@pytest.mark.parametrize("role", ["multiplicity", "label_support"])
def test_optional_directory_alias_keeps_disjoint_siblings_available(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    target = tmp_path / "optional-evidence" / role
    target.mkdir(parents=True)
    alias = tmp_path / "directory-aliases" / role
    alias.parent.mkdir()
    alias.symlink_to(target, target_is_directory=True)
    options = {role: alias}
    for roots in (
        stage._declared_source_roots(
            comparison,
            multiplicity=options.get("multiplicity"),
            label_support=options.get("label_support"),
        ),
        stage._source_roots(
            sources,
            multiplicity=options.get("multiplicity"),
            label_support=options.get("label_support"),
        ),
    ):
        with pytest.raises(ValueError, match="disjoint"):
            stage._disjoint(target / "nested", roots)
        stage._disjoint(alias.parent / "disjoint-output", roots)


@pytest.mark.parametrize("destination", ["output", "scratch"])
def test_self_declared_comparison_location_is_protected_before_any_write(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, destination: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    original = Path(sources.comparison.manifest_path)
    relocated = tmp_path / "relocated-comparison" / original.name
    relocated.parent.mkdir()
    relocated.write_bytes(original.read_bytes())
    protected = original.parent
    before = set(protected.rglob("*"))
    nested = protected / "must-not-be-created" / destination
    output = nested if destination == "output" else tmp_path / "output"
    scratch = nested if destination == "scratch" else tmp_path / "scratch"
    probe = Mock(side_effect=AssertionError("probe ran before self-declared source preflight"))
    verifier = Mock(side_effect=AssertionError("verification ran before source preflight"))
    monkeypatch.setattr(stage, "_probe_publication", probe)
    monkeypatch.setattr(stage, "open_verified_sources", verifier)
    with pytest.raises(ValueError, match="disjoint"):
        stage.run_table(relocated, output, scratch=scratch)
    probe.assert_not_called()
    verifier.assert_not_called()
    assert not nested.parent.exists()
    assert set(protected.rglob("*")) == before


@pytest.mark.parametrize(
    "role", ["comparison", "core", "noel", "jed", "ilyas", "multiplicity", "label_support"]
)
@pytest.mark.parametrize("protected_spelling", ["target", "alias"])
def test_verified_manifest_file_aliases_protect_both_locations(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str, protected_spelling: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    records: dict[str, ComparisonRun | WorkflowRun | ArmRun] = {
        "comparison": sources.comparison,
        "core": sources.core,
        **sources.arms,
    }
    options: dict[str, Path] = {}
    if role in records:
        original = Path(records[role].manifest_path)
    else:
        original = tmp_path / "optional-evidence" / role / "manifest.json"
        original.parent.mkdir(parents=True)
        original.write_text("{}\n")
    alias = tmp_path / "manifest-aliases" / f"{role}.json"
    alias.parent.mkdir()
    alias.symlink_to(original)
    if role in records:
        records[role].manifest_path = str(alias)
    else:
        options[role] = alias
    roots = stage._source_roots(
        sources,
        multiplicity=options.get("multiplicity"),
        label_support=options.get("label_support"),
    )
    protected = original.parent if protected_spelling == "target" else alias.parent
    with pytest.raises(ValueError, match="disjoint"):
        stage._disjoint(protected / "output", roots)


def test_real_verification_rejects_an_uncompleted_comparison_before_publication(
    tmp_path: Path,
) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    comparison = inputs / "comparison_run.json"
    failed = ComparisonRun(
        id="uncompleted-fixture",
        contrast_population=ContrastPopulation.all_paired,
        direction_basis=DirectionBasis.group_1_vs_group_0,
        status=ExecutionStatus.failed,
        started_at=datetime(2026, 10, 25, tzinfo=UTC),
        manifest_path=str(comparison),
        core_manifest_path=str(inputs / "absent-core.json"),
        core_manifest_sha256="a" * 64,
        core_run_id="absent-core",
        noel_manifest_path=str(inputs / "absent-noel.json"),
        noel_manifest_sha256="b" * 64,
        jed_manifest_path=str(inputs / "absent-jed.json"),
        jed_manifest_sha256="c" * 64,
        ilyas_manifest_path=str(inputs / "absent-ilyas.json"),
        ilyas_manifest_sha256="d" * 64,
        git_commit="e" * 40,
        git_dirty=False,
        analysis_specs=[],
        eligibility_ledgers=[],
        analysis_bases=[],
        calibration_verdicts=[],
        outputs=[],
    )
    comparison.write_text(failed.model_dump_json())
    with pytest.raises(ValueError, match="has status failed, expected completed"):
        stage.run_table(comparison, tmp_path / "output", scratch=tmp_path / "scratch")
    assert not (tmp_path / "output").exists()
    assert not (tmp_path / "scratch").exists()


def test_producer_ignores_untracked_files_but_binds_tracked_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repository_root = tmp_path / "producer"
    with Repo.init(repository_root) as repository:
        module = repository_root / "src/mtp/three_test_table/stage.py"
        module.parent.mkdir(parents=True)
        for name in (
            "src/mtp/three_test_table/stage.py",
            stage.METHODS,
            "uv.lock",
            "pyproject.toml",
            "Justfile",
        ):
            path = repository_root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("fixture\n")
        repository.index.add(
            [
                str(path.relative_to(repository_root))
                for path in repository_root.rglob("*")
                if path.is_file() and ".git" not in path.parts
            ]
        )
        actor = Actor("fixture", "fixture@example.invalid")
        repository.index.commit("Synthetic producer inventory", author=actor, committer=actor)
        monkeypatch.setattr(stage, "__file__", str(module))
        before = stage._producer(None)
        assert before["code_dirty"] is False
        (repository_root / "untracked-note.txt").write_text("unrelated\n")
        assert stage._producer(None) == before
        module.write_text("changed tracked source\n")
        changed = stage._producer(None)
        assert changed["code_dirty"] is True
        assert changed["implementation_sha256"] != before["implementation_sha256"]


@pytest.mark.parametrize("failure", ["write", "late_producer"])
def test_failed_publication_leaves_destination_absent_and_allows_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    output = tmp_path / "output"
    original_writer = pd.DataFrame.to_csv
    original_producer = stage._producer
    calls = 0

    def fail_write(frame: pd.DataFrame, *args: Any, **kwargs: Any) -> Any:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected publication write failure")
        return original_writer(frame, *args, **kwargs)

    def change_late(commit: str | None) -> dict[str, Any]:
        nonlocal calls
        calls += 1
        result = original_producer(commit)
        if calls == 3:
            result["implementation_sha256"] = "changed-during-write"
        return result

    expected: type[Exception]
    with monkeypatch.context() as failure_patch:
        if failure == "write":
            failure_patch.setattr(pd.DataFrame, "to_csv", fail_write)
            expected, message = OSError, "injected publication write failure"
        else:
            failure_patch.setattr(stage, "_producer", change_late)
            expected, message = ValueError, "changed before provenance publication"
        with pytest.raises(expected, match=message):
            stage.run_table(comparison, output, scratch=tmp_path / "scratch")
    assert not output.exists()
    assert not list(tmp_path.glob(".output-staging-*"))
    record = stage.run_table(comparison, output, scratch=tmp_path / "scratch")
    assert stage.verify_table(output, comparison, scratch=tmp_path / "scratch") == record


def test_destination_is_invisible_until_all_artifacts_are_ready(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    output = tmp_path / "output"
    writer = pd.DataFrame.to_csv
    written: list[Path] = []

    def inspect_write(frame: pd.DataFrame, path: Path, *args: Any, **kwargs: Any) -> Any:
        assert not output.exists(), "Partial output must not be visible at the destination"
        assert Path(path).parent.parent == output.parent
        written.append(Path(path))
        return writer(frame, path, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_csv", inspect_write)
    record = stage.run_table(
        Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
    )
    assert len(written) == 3
    assert {path.name for path in output.iterdir()} == {stage.MANIFEST, *record["outputs"]}
    assert not list(tmp_path.glob(".output-staging-*"))


@pytest.mark.parametrize("competitor", ["empty_directory", "directory", "file", "symlink"])
def test_concurrent_destination_is_never_overwritten(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, competitor: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    output = tmp_path / "output"
    publish = stage._publish_directory

    def create_competitor(staging: Path, destination: Path) -> None:
        if destination != output:
            publish(staging, destination)
            return
        if competitor in {"empty_directory", "directory"}:
            destination.mkdir()
            if competitor == "directory":
                (destination / "keep.txt").write_text("preserved")
        elif competitor == "file":
            destination.write_text("preserved")
        else:
            target = tmp_path / "other-output"
            target.mkdir()
            (target / "keep.txt").write_text("preserved")
            destination.symlink_to(target, target_is_directory=True)
        publish(staging, destination)

    monkeypatch.setattr(stage, "_publish_directory", create_competitor)
    with pytest.raises(FileExistsError):
        stage.run_table(
            Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
        )
    if competitor == "empty_directory":
        assert output.is_dir()
        assert not list(output.iterdir())
    elif competitor == "file":
        assert output.read_text() == "preserved"
    else:
        assert (output / "keep.txt").read_text() == "preserved"
        assert output.is_symlink() == (competitor == "symlink")
    assert not list(tmp_path.glob(".output-staging-*"))


@pytest.mark.parametrize("system", ["darwin", "linux"])
@pytest.mark.parametrize("error_number", [errno.EEXIST, errno.EACCES, errno.EINVAL, errno.ENOSYS])
def test_exclusive_publication_preserves_syscall_errors_without_fallback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, system: str, error_number: int
) -> None:
    source, destination = tmp_path / "staging", tmp_path / "output"
    source.mkdir()
    (source / "keep.txt").write_text("staged")

    def fail(*_args: Any) -> int:
        ctypes.set_errno(error_number)
        return -1

    function = Mock(side_effect=fail)
    library = SimpleNamespace(renamex_np=function, renameat2=function)
    monkeypatch.setattr(sys, "platform", system)
    monkeypatch.setattr(ctypes, "CDLL", Mock(return_value=library))
    with pytest.raises(OSError, match=rf"\[Errno {error_number}\]") as caught:
        stage._publish_directory(source, destination)
    assert caught.value.errno == error_number
    assert (source / "keep.txt").read_text() == "staged"
    assert not destination.exists()
    if system == "darwin":
        function.assert_called_once_with(os.fsencode(source), os.fsencode(destination), 4)
    else:
        function.assert_called_once_with(
            -100, os.fsencode(source), -100, os.fsencode(destination), 1
        )


@pytest.mark.parametrize("system", ["unsupported", "darwin", "linux"])
def test_exclusive_publication_fails_closed_when_unavailable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, system: str
) -> None:
    source, destination = tmp_path / "staging", tmp_path / "output"
    source.mkdir()
    monkeypatch.setattr(sys, "platform", system)
    monkeypatch.setattr(ctypes, "CDLL", Mock(return_value=SimpleNamespace()))
    with pytest.raises(NotImplementedError, match=r"exclusive.*rename"):
        stage._publish_directory(source, destination)
    assert source.is_dir()
    assert not destination.exists()


@pytest.mark.parametrize(
    "location", ["core", "comparison", "noel", "jed", "ilyas", "multiplicity", "label"]
)
def test_run_rejects_scratch_inside_a_source_before_creating_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, location: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    options: dict[str, Path] = {}
    if location == "core":
        protected = Path(sources.core.manifest_path).parent
    elif location == "comparison":
        protected = comparison.parent
    elif location in sources.arms:
        protected = Path(sources.arms[location].manifest_path).parent
    else:
        protected = tmp_path / location
        protected.mkdir()
        options["multiplicity" if location == "multiplicity" else "label_support"] = protected
    before = set(protected.rglob("*"))
    scratch = protected / "must-not-be-created"
    probe = Mock(side_effect=AssertionError("probe ran before source overlap checks"))
    monkeypatch.setattr(stage, "_probe_publication", probe, raising=False)
    verification = Mock(side_effect=AssertionError("verification ran before source overlap checks"))
    monkeypatch.setattr(stage, "open_verified_sources", verification)
    with pytest.raises(ValueError, match="disjoint"):
        stage.run_table(
            comparison,
            tmp_path / "output",
            scratch=scratch,
            multiplicity=options.get("multiplicity"),
            label_support=options.get("label_support"),
        )
    assert set(protected.rglob("*")) == before
    assert not (tmp_path / "output").exists()
    probe.assert_not_called()
    verification.assert_not_called()


@pytest.mark.parametrize("traits", [["a", "a"], ["a"], ["a", "b", "extra"]])
def test_component_merge_rejects_duplicate_missing_or_extra_traits(traits: list[str]) -> None:
    component = pd.DataFrame({"trait": traits, "value": range(len(traits))})
    with pytest.raises(ValueError, match="complete catalog exactly once"):
        stage._merge_tables(("a", "b"), [component])


def test_composition_host_selection_strips_whitespace_like_native(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    expected = stage.composition_profiles(sources)
    experiments = sources.product("experiments")
    metadata = pd.read_csv(experiments, sep="\t")
    metadata["host_species"] = "  Homo sapiens  "
    metadata.to_csv(experiments, sep="\t", index=False)
    selected = replace(sources, recorded_hashes={})
    pd.testing.assert_frame_equal(stage.composition_profiles(selected), expected)


@pytest.mark.parametrize("direction", ["increased", "INC"])
def test_composition_rejects_aliases_the_native_population_excludes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, direction: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    path = sources.product("count_profiles")
    profiles = pd.read_csv(path, sep="\t")
    profiles.loc[0, "inc_dec"] = direction
    profiles.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match=r"native.*direction"):
        stage.composition_profiles(replace(sources, recorded_hashes={}))


def test_run_rejects_native_overlap_mismatch_before_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)

    def altered_native(selected: VerifiedSources) -> tuple[pd.DataFrame, tuple[Path, ...]]:
        table, paths = native_reader(selected)
        tested = table.index[table.native_selected_status.eq("tested")]
        assert len(tested) > 0
        table.loc[tested[0], "native_selected_study_t_statistic"] += 1
        return table, paths

    monkeypatch.setattr(stage, "native_tables", altered_native)
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="overlap"):
        stage.run_table(
            Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
        )
    assert not (output / "three_test_run.json").exists()


def test_run_persists_complete_table_and_replays_without_machine_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    output = tmp_path / "evidence"
    before = file_sha256(sources.product("count_profiles"))
    record = stage.run_table(comparison, output, scratch=tmp_path / "scratch")
    table = pd.read_csv(output / "per_trait_evidence.tsv", sep="\t")
    assert list(table.trait) == list(sources.features["count_profiles"])
    assert len(table) == 7
    assert table.ora_calibration_status.eq(
        "uncalibrated (methods/calibration.md excludes native ORA)"
    ).all()
    assert table.carrier_status.eq("not_supplied").all()
    assert record["core_run_id"] == "fixture-core"
    assert record["comparison_run_id"] == "fixture-comparison"
    assert record["source_run_ids"] == {
        "noel": "fixture-noel",
        "jed": "fixture-jed",
        "ilyas": "fixture-ilyas",
    }
    assert len(record["code_commit"]) == 40
    assert record["primary_ranking"] == "undecided"
    assert str(tmp_path) not in (output / "three_test_run.json").read_text()
    assert file_sha256(sources.product("count_profiles")) == before
    assert stage.verify_table(output, comparison, scratch=tmp_path / "scratch") == record


def test_existing_output_and_source_ancestor_are_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    existing = tmp_path / "existing"
    existing.mkdir()
    (existing / "keep.txt").write_text("keep")
    with pytest.raises(FileExistsError):
        stage.run_table(comparison, existing, scratch=tmp_path / "scratch")
    assert (existing / "keep.txt").read_text() == "keep"
    with pytest.raises(ValueError, match="disjoint"):
        stage.run_table(comparison, comparison.parent / "nested", scratch=tmp_path / "scratch")


@pytest.mark.parametrize("path", ["/tmp/mtp-table", "/private/tmp/mtp-table"])
def test_system_temporary_output_rejected_before_source_read(
    path: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(tempfile, "gettempdir", lambda: "/tmp")
    with pytest.raises(ValueError, match="durable"):
        stage.run_table(Path("absent.json"), Path(path), scratch=Path(path) / "scratch")


def test_falsely_requested_code_commit_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="commit"):
        stage.run_table(
            Path(sources.comparison.manifest_path),
            tmp_path / "output",
            scratch=tmp_path / "scratch",
            code_commit="0" * 40,
        )
    assert not (tmp_path / "output" / "three_test_run.json").exists()


def test_verification_rejects_table_tampering_even_if_local_hash_is_updated(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    output = tmp_path / "output"
    stage.run_table(comparison, output, scratch=tmp_path / "scratch")
    table_path = output / "per_trait_evidence.tsv"
    table = pd.read_csv(table_path, sep="\t")
    table.loc[0, "study_mean_delta_z"] = 200.0
    table.to_csv(table_path, sep="\t", index=False)
    manifest = output / "three_test_run.json"
    record = json.loads(manifest.read_text())
    record["outputs"]["per_trait_evidence.tsv"]["sha256"] = file_sha256(table_path)
    record["outputs"]["per_trait_evidence.tsv"]["size_bytes"] = table_path.stat().st_size
    manifest.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="replay"):
        stage.verify_table(output, comparison, scratch=tmp_path / "scratch")


def test_source_tampering_is_not_promoted_to_a_new_snapshot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    table = sources.product("count_profiles")
    table.write_text(table.read_text().replace("\t2\t", "\t9\t", 1))
    with pytest.raises(ValueError, match=r"(recorded|hash|Source)"):
        stage.run_table(
            Path(sources.comparison.manifest_path),
            tmp_path / "output",
            scratch=tmp_path / "scratch",
        )
    assert not (tmp_path / "output" / "three_test_run.json").exists()


def test_validated_population_and_direction_strings_publish_successfully(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    # Assignment validation applies the same use_enum_values model configuration
    # as reading a manifest; model_construct in the base fixture bypasses it.
    sources.comparison.contrast_population = ContrastPopulation.all_paired
    sources.comparison.direction_basis = DirectionBasis.group_1_vs_group_0
    assert type(sources.comparison.contrast_population) is str
    assert type(sources.comparison.direction_basis) is str
    output = tmp_path / "output"
    record = stage.run_table(
        Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
    )
    assert record["contrast_population"] == "all_paired"
    assert record["direction_basis"] == "group_1_vs_group_0"
    table = pd.read_csv(output / "per_trait_evidence.tsv", sep="\t")
    assert table.contrast_population.eq("all_paired").all()
    assert table.direction_basis.eq("group_1_vs_group_0").all()


@pytest.mark.parametrize("prefix", ["", "biolink:has_phenotype|"])
@pytest.mark.parametrize(
    ("identifier", "correct", "incorrect"),
    [
        ("METPO:1000603", "anaerobic", "obligately anaerobic"),
        ("METPO:1000607", "obligately anaerobic", "anaerobic"),
    ],
)
def test_oxygen_identity_pin_applies_to_bare_and_qualified_traits(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    prefix: str,
    identifier: str,
    correct: str,
    incorrect: str,
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    feature = prefix + identifier
    selected = replace(sources, features={**sources.features, "count_profiles": (feature,)})
    metadata = sources.product("feature_metadata")
    write_table(metadata, [{"column_id": feature, "object_name": correct, "column_label": ""}])
    assert stage._labels(selected).to_dict("records") == [
        {"trait": feature, "trait_label": correct}
    ]
    write_table(metadata, [{"column_id": feature, "object_name": incorrect, "column_label": ""}])
    with pytest.raises(ValueError, match=identifier):
        stage._labels(selected)


@pytest.mark.parametrize(
    "location", ["core", "comparison", "noel", "jed", "ilyas", "output", "multiplicity", "label"]
)
def test_verify_rejects_overlapping_scratch_before_creating_any_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, location: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    output = tmp_path / "output"
    stage.run_table(comparison, output, scratch=tmp_path / "initial-scratch")
    options: dict[str, Path] = {}
    if location == "core":
        protected = Path(sources.core.manifest_path).parent
    elif location == "comparison":
        protected = comparison.parent
    elif location in sources.arms:
        protected = Path(sources.arms[location].manifest_path).parent
    elif location == "output":
        protected = output
    else:
        protected = tmp_path / location
        protected.mkdir()
        options["multiplicity" if location == "multiplicity" else "label_support"] = protected
    before = {path.relative_to(protected): path.is_dir() for path in protected.rglob("*")}
    scratch = protected / "must-not-be-created"
    with pytest.raises(ValueError, match="disjoint"):
        stage.verify_table(output, comparison, scratch=scratch, **options)
    after = {path.relative_to(protected): path.is_dir() for path in protected.rglob("*")}
    assert after == before
    assert not scratch.exists()


@pytest.mark.parametrize(
    "role", ["comparison_file", "comparison_directory", "core", "noel", "jed", "ilyas"]
)
@pytest.mark.parametrize("protected_spelling", ["target", "alias"])
def test_verify_preserves_declared_alias_roots_when_verified_records_are_canonical(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str, protected_spelling: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    if role.startswith("comparison_"):
        original = comparison
    elif role == "core":
        original = Path(sources.core.manifest_path)
    else:
        original = Path(sources.arms[role].manifest_path)
    alias = tmp_path / "declared-aliases" / original.name
    alias.parent.mkdir()
    alias.symlink_to(original)
    if role.startswith("comparison_"):
        comparison = alias.parent if role == "comparison_directory" else alias
    else:
        declared = json.loads(comparison.read_text())
        declared[f"{role}_manifest_path"] = str(alias)
        comparison.write_text(json.dumps(declared))
    # The source fixture returns canonical model paths, as upstream verification
    # can do even though the caller or the comparison record used a file alias.
    output = tmp_path / "output"
    stage.run_table(comparison, output, scratch=tmp_path / "initial-scratch")
    verified = Mock(return_value=sources)
    monkeypatch.setattr(stage, "open_verified_sources", verified)
    replay = Mock(side_effect=AssertionError("replay ran before alias scratch preflight"))
    monkeypatch.setattr(stage, "run_table", replay)
    protected = original.parent if protected_spelling == "target" else alias.parent
    before = set(protected.rglob("*"))
    scratch = protected / "must-not-be-created" / "scratch"
    with pytest.raises(ValueError, match="disjoint"):
        stage.verify_table(output, comparison, scratch=scratch)
    verified.assert_not_called()
    replay.assert_not_called()
    assert not scratch.parent.exists()
    assert set(protected.rglob("*")) == before


@pytest.mark.parametrize("protected_spelling", ["target", "alias"])
def test_verify_output_manifest_file_alias_protects_both_bundle_locations(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, protected_spelling: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    original = tmp_path / "output"
    record = stage.run_table(comparison, original, scratch=tmp_path / "initial-scratch")
    alias_root = tmp_path / "output-alias"
    alias_root.mkdir()
    for name in record["outputs"]:
        (alias_root / name).write_bytes((original / name).read_bytes())
    alias = alias_root / stage.MANIFEST
    alias.symlink_to(original / stage.MANIFEST)
    # A manifest file alias with matching sibling tables is an accepted bundle;
    # replay must keep this inventory lookup behavior for a disjoint scratch.
    assert stage.verify_table(alias, comparison, scratch=tmp_path / "valid-replay") == record
    verifier = Mock(side_effect=AssertionError("source verification ran before output preflight"))
    replay = Mock(side_effect=AssertionError("replay ran before output scratch preflight"))
    monkeypatch.setattr(stage, "open_verified_sources", verifier)
    monkeypatch.setattr(stage, "run_table", replay)
    protected = original if protected_spelling == "target" else alias_root
    before = set(protected.rglob("*"))
    scratch = protected / "must-not-be-created" / "scratch"
    with pytest.raises(ValueError, match="disjoint"):
        stage.verify_table(alias, comparison, scratch=scratch)
    verifier.assert_not_called()
    replay.assert_not_called()
    assert not scratch.parent.exists()
    assert set(protected.rglob("*")) == before


@pytest.mark.parametrize("mutation", ["artifact", "manifest", "extra", "remove", "symlink"])
def test_optional_multiplicity_change_after_verification_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mutation: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    sidecar = tmp_path / "multiplicity"
    original = run_post_arm(comparison, sidecar)
    target = sidecar / original.families[0].corrected.identity
    verify_optional = _verify_post_arm_from_sources

    def verified_then_changed(path: Path, bound: VerifiedSources) -> PostArmMultiplicityRun:
        record = verify_optional(path, bound)
        if mutation == "artifact":
            target.write_bytes(target.read_bytes() + b"\n")
        elif mutation == "manifest":
            manifest = sidecar / "multiplicity_run.json"
            manifest.write_bytes(manifest.read_bytes() + b"\n")
        elif mutation == "extra":
            (sidecar / "unexpected.txt").write_text("unexpected")
        elif mutation == "remove":
            target.unlink()
        else:
            (sidecar / "unexpected-link").symlink_to(target)
        return record

    monkeypatch.setattr(stage, "_verify_post_arm_from_sources", verified_then_changed)
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="Multiplicity"):
        stage.run_table(comparison, output, scratch=tmp_path / "scratch", multiplicity=sidecar)
    assert not (output / stage.MANIFEST).exists()


def test_verified_optional_multiplicity_is_bound_and_replays(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    sidecar = tmp_path / "multiplicity"
    multiplicity = run_post_arm(comparison, sidecar)
    output = tmp_path / "output"
    record = stage.run_table(comparison, output, scratch=tmp_path / "scratch", multiplicity=sidecar)
    manifest = sidecar / "multiplicity_run.json"
    identity = f"{manifest.name}@sha256:{file_sha256(manifest)}"
    assert record["inputs"][identity]["sha256"] == file_sha256(manifest)
    for family in multiplicity.families:
        for artifact in (family.declared, family.raw, family.corrected):
            identity = f"{Path(artifact.identity).name}@sha256:{artifact.sha256}"
            assert record["inputs"][identity]["sha256"] == artifact.sha256
    assert (
        stage.verify_table(output, comparison, scratch=tmp_path / "scratch", multiplicity=sidecar)
        == record
    )


def test_optional_inventory_change_during_adaptation_is_rejected_before_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    sidecar = tmp_path / "multiplicity"
    run_post_arm(comparison, sidecar)

    def change_after_initial_binding(
        bound: VerifiedSources, audit: Path | None
    ) -> tuple[pd.DataFrame, tuple[Path, ...]]:
        result = carrier_reader(bound, audit)
        (sidecar / "late-extra.txt").write_text("unexpected")
        return result

    monkeypatch.setattr(stage, "carrier_table", change_after_initial_binding)
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="Multiplicity input inventory"):
        stage.run_table(comparison, output, scratch=tmp_path / "scratch", multiplicity=sidecar)
    assert not (output / stage.MANIFEST).exists()


@pytest.mark.parametrize("source", ["process", "environment"])
def test_declared_temporary_roots_are_rejected_as_written_and_resolved(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source: str
) -> None:
    process = tmp_path / "process-temporary"
    environment = tmp_path / "environment-temporary"
    monkeypatch.setattr(tempfile, "gettempdir", lambda: str(process))
    monkeypatch.setenv("TMPDIR", str(environment))
    temporary = process if source == "process" else environment
    for candidate in (temporary, temporary / "output"):
        with pytest.raises(ValueError, match="durable"):
            stage.require_durable(candidate)
    sibling = tmp_path / f"{temporary.name}-sibling" / "output"
    assert stage.require_durable(sibling) == sibling.resolve()


@pytest.mark.parametrize("direction", ["into_temporary", "out_of_temporary"])
def test_temporary_path_aliases_cannot_hide_written_or_resolved_roots(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, direction: str
) -> None:
    temporary = tmp_path / "declared-temporary"
    durable = tmp_path / "durable"
    temporary.mkdir()
    durable.mkdir()
    monkeypatch.setattr(tempfile, "gettempdir", lambda: str(temporary))
    monkeypatch.setenv("TMPDIR", str(temporary))
    if direction == "into_temporary":
        alias = durable / "link"
        alias.symlink_to(temporary, target_is_directory=True)
    else:
        alias = temporary / "link"
        alias.symlink_to(durable, target_is_directory=True)
    with pytest.raises(ValueError, match="durable"):
        stage.require_durable(alias / "output")


@pytest.mark.skipif(sys.platform != "linux", reason="Documents the Linux-only #1051 gap")
@pytest.mark.parametrize("path", ["/tmp/mtp-table", "/var/tmp/mtp-table"])
def test_linux_root_not_declared_by_process_or_environment_remains_a_known_gap(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, path: str
) -> None:
    custom = tmp_path / "custom-temporary"
    monkeypatch.setattr(tempfile, "gettempdir", lambda: str(custom))
    monkeypatch.setenv("TMPDIR", str(custom))
    # No directory is created: this states the accepted compatibility limitation.
    assert stage.require_durable(Path(path)) == Path(path).resolve()


def test_verify_nested_replay_preserves_inherited_temporary_environment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    output = tmp_path / "output"
    original = stage.run_table(comparison, output, scratch=tmp_path / "scratch")
    inherited = os.environ["TMPDIR"]
    assert stage.verify_table(output, comparison, scratch=tmp_path / "scratch") == original
    assert os.environ["TMPDIR"] == inherited
    assert tempfile.gettempdir() == inherited
