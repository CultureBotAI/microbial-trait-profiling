"""Real producer publication must retain the parent bytes it verified (#1216)."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Sequence
from pathlib import Path

import pytest

from mtp import integrated_evidence
from mtp.workflow.arms.runtime import execute_arm
from mtp.workflow.arms.types import ArmContext, Command
from mtp.workflow.comparison import execute_comparison_run
from test_comparison_run import _consumed_plan
from test_comparison_run import stage as stage


def _changed_metadata(original: bytes) -> bytes:
    """Change only a valid top-level timestamp, retaining every scientific field."""
    payload = json.loads(original)
    started = payload["started_at"]
    assert isinstance(started, str)
    assert len(started) >= 10
    # Leap years keep even a February 29 timestamp valid without changing width.
    year = "1996" if started.startswith("2000") else "2000"
    replacement = year + started[4:]
    pattern = rb'("started_at"\s*:\s*)' + re.escape(json.dumps(started).encode())
    changed, replacements = re.subn(
        pattern,
        lambda match: match.group(1) + json.dumps(replacement).encode(),
        original,
    )
    assert replacements == 1
    assert len(changed) == len(original)
    assert hashlib.sha256(changed).digest() != hashlib.sha256(original).digest()
    assert json.loads(changed) == {**payload, "started_at": replacement}
    return changed


def _safe_root(
    stage_paths: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> Path:
    """Keep the real private-root guard active with a common synthetic parent."""
    root = tmp_path.parent.resolve()
    assert stage_paths["comparison"].parent.parent.resolve() == root
    # The fixture's core, external raw/KG inputs and this test's output are siblings
    # below the pytest base. Declare its parent as TMPDIR, not the provenance root.
    monkeypatch.setenv("TMPDIR", str(root.parent))
    assert all(path.resolve().is_relative_to(root) for path in stage_paths.values())
    core = json.loads(stage_paths["core"].read_bytes())
    for artifact in core["artifacts"]:
        path = Path(artifact["path"])
        if path.is_absolute():
            assert path.resolve().is_relative_to(root)
    return root


def _assert_no_completed_record(path: Path) -> None:
    if path.exists():
        assert json.loads(path.read_bytes())["status"] == "failed"


def test_arm_producer_refuses_core_change_after_bundle_verification(
    stage: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = _safe_root(stage, tmp_path, monkeypatch)
    manifest = stage["core"]
    original = manifest.read_bytes()
    changed = _changed_metadata(original)
    output = tmp_path / "arm-after-core-change"
    real_plan = _consumed_plan("noel")
    injected = False

    def plan(context: ArmContext) -> list[Command]:
        nonlocal injected
        commands = real_plan(context)
        assert commands
        assert context.manifest == manifest.resolve()
        assert manifest.read_bytes() == original
        manifest.write_bytes(changed)
        injected = True
        return commands

    try:
        with pytest.raises(ValueError, match=r"(?i)(manifest|checksum|sha256|digest)"):
            execute_arm(
                "noel",
                manifest,
                output,
                plan,
                replicates=2,
                seed=0,
                figures=False,
                provenance_root=root,
            )
        assert injected, "The real planner must run after core-bundle verification"
        _assert_no_completed_record(output / "arm_run.json")
    finally:
        manifest.write_bytes(original)
        assert manifest.read_bytes() == original


@pytest.mark.parametrize("parent", ["core", "noel", "jed", "ilyas"])
def test_comparison_producer_refuses_parent_change_before_publication(
    stage: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, parent: str
) -> None:
    root = _safe_root(stage, tmp_path, monkeypatch)
    manifest = stage[parent]
    original = manifest.read_bytes()
    changed = _changed_metadata(original)
    output = tmp_path / f"comparison-after-{parent}-change"
    pooled = (
        stage["comparison"]
        / "arm1_pooled_stats"
        / "trait_profile_paired_contrast_trait_stats_pooled.tsv"
    )
    real_main = integrated_evidence.main
    injected = False

    def finish_then_change(argv: Sequence[str] | None = None) -> int:
        nonlocal injected
        result = real_main(argv)
        assert result == 0
        produced = output / "integrated_evidence" / "integrated_evidence_provenance.json"
        assert produced.is_file()
        assert produced.stat().st_size > 0
        assert manifest.read_bytes() == original
        manifest.write_bytes(changed)
        injected = True
        return result

    monkeypatch.setattr(integrated_evidence, "main", finish_then_change)
    try:
        with pytest.raises(ValueError, match=r"(?i)(manifest|checksum|sha256|digest)"):
            execute_comparison_run(
                stage["core"],
                stage["noel"],
                stage["jed"],
                stage["ilyas"],
                output,
                arm1_pooled_stats=pooled,
                provenance_root=root,
            )
        assert injected, "The real integrated-evidence producer must finish before the swap"
        _assert_no_completed_record(output / "comparison_run.json")
    finally:
        manifest.write_bytes(original)
        assert manifest.read_bytes() == original
