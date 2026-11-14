"""Pin the supplied-root decision at the real execution-record constructors (#1820).

Scientific inputs use routing fakes; the separate real execution tests remain
the authority for serialized records, verification and historical compatibility.
"""

from __future__ import annotations

import importlib.metadata
import json
import sys
from pathlib import Path
from typing import NoReturn

import pytest

from mtp import integrated_evidence
from mtp.workflow import comparison
from mtp.workflow.arms import runtime
from mtp.workflow.arms.types import ArmContext, Command
from test_workflow_provenance_root_preflight import Inputs
from test_workflow_provenance_root_preflight import inputs as inputs


class CapturedRecord(Exception):
    """The actual producer supplied all constructor keyword arguments."""


@pytest.mark.parametrize("supplied", [False, True], ids=["omitted", "supplied"])
def test_arm_constructor_receives_the_actual_root_decision(
    inputs: Inputs, monkeypatch: pytest.MonkeyPatch, supplied: bool
) -> None:
    output = inputs.root / "noel-output"
    root = inputs.root if supplied else None
    captured: list[dict[str, object]] = []
    contexts: list[ArmContext] = []

    def planner(context: ArmContext) -> list[Command]:
        contexts.append(context)
        arguments: tuple[str, ...] = (
            sys.executable,
            "-m",
            "mtp.integrated_evidence",
            "--core-run",
            str(context.manifest),
            "--output-dir",
            str(output / "integrated_evidence"),
        )
        if root is not None:
            arguments += ("--provenance-root", str(root))
        return [
            Command(
                "integrated_evidence",
                arguments,
                (context.manifest, context.products["count_profiles"]),
                (output / "integrated_evidence/record.json",),
            )
        ]

    def capture(**values: object) -> NoReturn:
        captured.append(values)
        raise CapturedRecord

    monkeypatch.setattr(runtime, "ArmRun", capture)
    monkeypatch.setattr(runtime, "_code_inventory", lambda _root: [])
    monkeypatch.setattr(importlib.metadata, "distributions", lambda: ())
    with pytest.raises(CapturedRecord):
        runtime.execute_arm("noel", inputs.manifests["core"], output, planner, provenance_root=root)

    assert len(contexts) == 1
    assert contexts[0].provenance_root == root
    assert len(captured) == 1
    assert captured[0].get("provenance_root_supplied") is supplied
    assert {
        name: value for name, value in captured[0].items() if name.startswith("provenance_root")
    } == {"provenance_root_supplied": supplied}
    assert not output.exists()


@pytest.mark.parametrize("supplied", [False, True], ids=["omitted", "supplied"])
def test_comparison_constructor_receives_the_actual_root_decision(
    inputs: Inputs, monkeypatch: pytest.MonkeyPatch, supplied: bool
) -> None:
    output = inputs.root / "comparison-output"
    root = inputs.root if supplied else None
    captured: list[dict[str, object]] = []
    parsed_roots: list[Path | None] = []

    def integrated(arguments: list[str]) -> int:
        parsed = integrated_evidence._build_parser().parse_args(arguments)
        assert parsed.core_run == inputs.manifests["core"]
        assert parsed.arm1_pooled_stats == inputs.pooled
        assert parsed.provenance_root == root
        parsed_roots.append(parsed.provenance_root)
        directory = Path(parsed.output_dir)
        directory.mkdir()
        (directory / "integrated_evidence_provenance.json").write_text(
            json.dumps({"synthetic_scientific_stage": True}), encoding="utf-8"
        )
        return 0

    def capture(**values: object) -> NoReturn:
        captured.append(values)
        raise CapturedRecord

    async def unexpected_pooled_generation(*_args: object) -> int:
        raise AssertionError("The existing pooled table must avoid a scientific subprocess")

    monkeypatch.setattr(integrated_evidence, "main", integrated)
    monkeypatch.setattr(comparison, "ComparisonRun", capture)
    monkeypatch.setattr(comparison, "reconcile_ledger", lambda *_args: None)
    monkeypatch.setattr(comparison, "_exec_command", unexpected_pooled_generation)
    with pytest.raises(CapturedRecord):
        comparison.execute_comparison_run(
            inputs.manifests["core"],
            inputs.manifests["noel"],
            inputs.manifests["jed"],
            inputs.manifests["ilyas"],
            output,
            arm1_pooled_stats=inputs.pooled,
            provenance_root=root,
        )

    assert parsed_roots == [root]
    assert len(captured) == 1
    assert captured[0].get("provenance_root_supplied") is supplied
    assert {
        name: value for name, value in captured[0].items() if name.startswith("provenance_root")
    } == {"provenance_root_supplied": supplied}
    assert (output / "integrated_evidence/integrated_evidence_provenance.json").is_file()
    assert not (output / "comparison_run.json").exists()
