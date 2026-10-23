"""Noel's rerun cannot silently select checkout-local scientific inputs (#658)."""

from __future__ import annotations

import argparse
import importlib
import sys
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT
from mtp.workflow.arms.noel import build_commands
from mtp.workflow.arms.types import ArmContext, Command


@pytest.fixture
def context(tmp_path: Path) -> ArmContext:
    products = {
        key: tmp_path / "core" / filename
        for key, filename in {
            "count_profiles": "study_trait_profiles_wide.tsv",
            "experiments": "experiments.tsv",
            "raw_dump": "full_dump.csv",
            "feature_metadata": "taxa_trait_matrix_columns.tsv",
            "kg_nodes": "merged-kg_nodes.tsv",
        }.items()
    }
    return ArmContext(
        REPO_ROOT,
        tmp_path / "core" / "run.json",
        tmp_path / "rerun",
        tmp_path / "rerun" / "inputs",
        products,
        "a" * 40,
        None,
        None,
        True,
    )


def _option(command: Command, option: str) -> str:
    return command.argv[command.argv.index(option) + 1]


def test_every_noel_command_has_explicit_inputs_and_outputs(context: ArmContext) -> None:
    commands = build_commands(context)
    assert len({command.id for command in commands}) == len(commands) == 10
    for command in commands:
        assert command.argv[:2] == (sys.executable, "-m")
        assert command.inputs
        assert command.outputs
        assert all(path.is_relative_to(context.output) for path in command.outputs)
        assert any(path.suffix in {".tsv", ".json"} for path in command.outputs)
        assert not any(str(REPO_ROOT / "data") in argument for argument in command.argv)
        if "--profile-input" in command.argv:
            assert _option(command, "--profile-input") == str(context.products["count_profiles"])
            assert _option(command, "--experiments-input") == str(context.products["experiments"])
            assert _option(command, "--full-dump-input") == str(context.products["raw_dump"])
        if "--inputs" in command.argv:
            assert _option(command, "--inputs") == str(context.inputs)
    by_id = {command.id: command for command in commands}
    for identifier in ("directional_screen", "validation_evidence", "integrated_evidence"):
        assert _option(by_id[identifier], "--core-run") == str(context.manifest)
    assert _option(by_id["validation_evidence"], "--nodes-input") == str(
        context.products["kg_nodes"]
    )
    assert _option(by_id["taxon_contributions"], "--kg-dir") == str(context.inputs)
    assert not any("investigate_" in argument for command in commands for argument in command.argv)


def test_integrated_statuses_are_run_local_and_dependencies_are_hashed(
    context: ArmContext,
) -> None:
    command = build_commands(context)[-1]
    for option in (
        "--calibration-status",
        "--orientation-provenance",
        "--sensitivity-status",
        "--coverage-status",
        "--influence-status",
    ):
        assert Path(_option(command, option)).is_relative_to(context.output)
    for option in ("--calibration-status", "--orientation-provenance", "--influence-status"):
        assert Path(_option(command, option)) in command.inputs
    # Absent assessments are explicitly pending; no copied legacy status enters this run.
    for option in ("--sensitivity-status", "--coverage-status"):
        assert Path(_option(command, option)).parent.name == "pending_assessments"


class _Parsed(Exception):
    def __init__(self, namespace: argparse.Namespace) -> None:
        super().__init__()
        self.namespace = namespace


def _parse_command(command: Command, monkeypatch: pytest.MonkeyPatch) -> argparse.Namespace:
    """Use each real parser, stopping immediately before scientific execution."""
    parse_args = argparse.ArgumentParser.parse_args

    def capture(
        parser: argparse.ArgumentParser,
        args: Sequence[str] | None = None,
        namespace: argparse.Namespace | None = None,
    ) -> argparse.Namespace:
        raise _Parsed(parse_args(parser, args, namespace))

    module = importlib.import_module(command.argv[2])
    with monkeypatch.context() as patch:
        patch.setattr(argparse.ArgumentParser, "parse_args", capture)
        patch.setattr(sys, "argv", [command.argv[2], *command.argv[3:]])
        with pytest.raises(_Parsed) as parsed:
            module.main()
    return parsed.value.namespace


def test_real_parsers_accept_all_commands_and_preserve_statistical_defaults(
    context: ArmContext, monkeypatch: pytest.MonkeyPatch
) -> None:
    parsed = {
        command.id: _parse_command(command, monkeypatch) for command in build_commands(context)
    }
    assert parsed["empirical_calibration"].replicates == 1400
    assert parsed["empirical_calibration"].power_replicates is None
    assert parsed["empirical_calibration"].seed == 182
    assert parsed["screen_calibration"].replicates == 100
    assert parsed["screen_calibration"].seed == 182
    assert parsed["shared_flip_assignments"].flip_assignments_only
    assert parsed["native_taxon_influence"].code_commit == "a" * 40


def test_overrides_are_explicit_without_changing_fixed_shared_flip_protocol(
    context: ArmContext, monkeypatch: pytest.MonkeyPatch
) -> None:
    context = replace(context, seed=19, replicates=3, code_commit="a" * 40 + "+dirty")
    commands = {command.id: command for command in build_commands(context)}
    for identifier in ("empirical_calibration", "screen_calibration"):
        parsed = _parse_command(commands[identifier], monkeypatch)
        assert parsed.replicates == 3
        assert parsed.seed == 19
    flips = commands["shared_flip_assignments"]
    assert "--replicates" not in flips.argv
    assert _option(flips, "--seed") == "19"
    native = _parse_command(commands["native_taxon_influence"], monkeypatch)
    assert native.code_commit == "a" * 40
    assert native.code_dirty is True


def test_figures_flag_does_not_invoke_historical_branch_figures(context: ArmContext) -> None:
    assert build_commands(context) == build_commands(replace(context, figures=False))
