"""Noel checks privacy tokens in selected outputs before producing them (#1855)."""

from __future__ import annotations

import re
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import NoReturn

import pytest

from mtp import provenance_paths
from mtp.workflow.arms import runtime
from mtp.workflow.arms.types import ArmContext, Command
from workflow_provenance_fixtures import routing_core


class ReachedAfterPreflight(Exception):
    """A valid plan reached code inventory before any output or command execution."""


@pytest.mark.parametrize("private_token", [True, False], ids=["private", "safe"])
def test_noel_checks_tokens_in_selected_earlier_outputs_before_execution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, private_token: bool
) -> None:
    root = tmp_path / "source-bundle"
    root.mkdir()
    manifest = root / "core.json"
    products = {
        "count_profiles": root / "count.tsv",
        "ora_profiles": root / "ora.tsv",
        "experiments": root / "experiments.tsv",
        "raw_dump": root / "dump.csv",
    }
    for path in (manifest, *products.values()):
        path.write_text("synthetic routing input\n", encoding="utf-8")
    core = routing_core(manifest, products)
    bundle = SimpleNamespace(manifest=manifest, run=core, path=products.__getitem__)
    monkeypatch.setenv("TMPDIR", str(tmp_path))
    monkeypatch.setattr(provenance_paths, "_account_names", lambda: {"researcher"})
    monkeypatch.setattr(runtime, "open_core_bundle", lambda *_args, **_kwargs: bundle)
    monkeypatch.setattr(runtime, "preflight_consumer_outputs", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(runtime, "repository_state", lambda _root: ("a" * 40, False))

    output = root / "noel-output"
    component = "researcher" if private_token else "selected"
    future = output / component / "screen.tsv"
    record = output / "integrated_evidence/integrated_evidence_provenance.json"
    planned: list[Command] = []
    inventoried: list[Path] = []
    executed: list[str] = []

    def planner(context: ArmContext) -> list[Command]:
        assert context.provenance_root == root.resolve()
        assert not future.exists()
        commands = [
            Command(
                "planned_screen",
                (
                    sys.executable,
                    "-m",
                    "mtp.trait_screen",
                    "--core-run",
                    str(context.manifest),
                    "--output-dir",
                    str(future.parent),
                ),
                (context.manifest, context.products["count_profiles"]),
                (future,),
            ),
            Command(
                "integrated_evidence",
                (
                    sys.executable,
                    "-m",
                    "mtp.integrated_evidence",
                    "--core-run",
                    str(context.manifest),
                    "--full-screen",
                    str(future),
                    "--output-dir",
                    str(record.parent),
                    "--provenance-root",
                    str(root),
                ),
                (context.manifest, context.products["count_profiles"], future),
                (record,),
            ),
        ]
        planned.extend(commands)
        return commands

    def after_preflight(repository: Path) -> NoReturn:
        inventoried.append(repository)
        raise ReachedAfterPreflight

    async def forbidden_execution(command: Command, _repository: Path, _log: Path) -> int:
        executed.append(command.id)
        raise AssertionError("The planned-source test must not execute scientific commands")

    monkeypatch.setattr(runtime, "_code_inventory", after_preflight)
    monkeypatch.setattr(runtime, "_run_process", forbidden_execution)

    def invoke() -> None:
        runtime.execute_arm("noel", manifest, output, planner, provenance_root=root)

    if private_token:
        with pytest.raises(
            ValueError, match=r"^Provenance contains a local account name"
        ) as policy:
            provenance_paths.PortablePaths.create(
                record, (products["ora_profiles"], future), root=root
            )
        expected = str(policy.value)
        try:
            with pytest.raises(ValueError, match=f"^{re.escape(expected)}$") as stopped:
                invoke()
        except ReachedAfterPreflight as error:
            raise AssertionError(
                "A selected private future source must be refused before code inventory"
            ) from error
        assert str(stopped.value) == expected
        assert inventoried == []
    else:
        try:
            with pytest.raises(ReachedAfterPreflight):
                invoke()
        except (ValueError, OSError) as error:
            raise AssertionError("A safe selected future source must pass preflight") from error
        assert len(inventoried) == 1

    assert [command.id for command in planned] == ["planned_screen", "integrated_evidence"]
    assert planned[0].outputs == (future,)
    assert future in planned[1].inputs
    assert not future.exists()
    assert not output.exists()
    assert executed == []
