"""Noel's rerun cannot silently select checkout-local scientific inputs (#658)."""

from __future__ import annotations

import argparse
import importlib
import sys
from collections.abc import Callable, Sequence
from dataclasses import replace
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT
from mtp.taxa_trait_matrix import StrainRollup
from mtp.workflow.arms import runtime
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


@pytest.mark.parametrize("scope", ["off", "phenotype", "full"])
def test_the_ledger_step_reconciles_under_the_core_scope(
    context: ArmContext, monkeypatch: pytest.MonkeyPatch, scope: StrainRollup
) -> None:
    """Noel's Table 1 ledger rebuilds the matrix under its core's strain rollup scope (#788)."""
    context = replace(context, strain_rollup=scope)
    commands = build_commands(context)
    ledger = next(command for command in commands if command.id == "taxon_contributions")

    assert _option(ledger, "--strain-rollup") == scope
    assert _parse_command(ledger, monkeypatch).strain_rollup == scope
    runtime._validate_plan(commands, context)


def test_a_ledger_step_under_another_scope_than_the_core_is_refused(context: ArmContext) -> None:
    commands = build_commands(replace(context, strain_rollup="phenotype"))

    with pytest.raises(ValueError, match="strain rollup"):
        runtime._validate_plan(commands, replace(context, strain_rollup="off"))


def _ledger_args(commands: list[Command], edit: Callable[[list[str]], list[str]]) -> list[Command]:
    """The plan with the arguments of its Table 1 ledger step edited."""
    return [
        replace(command, argv=(*command.argv[:3], *edit(list(command.argv[3:]))))
        if command.id == "taxon_contributions"
        else command
        for command in commands
    ]


def _without_scope(args: list[str]) -> list[str]:
    index = args.index("--strain-rollup")
    return [*args[:index], *args[index + 2 :]]


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        pytest.param(
            lambda args: [*args, "--strain-rollup", "phenotype"], "scope 'phenotype'", id="repeated"
        ),
        pytest.param(
            lambda args: [*_without_scope(args), "--strain-rollup=phenotype"],
            "scope 'phenotype'",
            id="equals",
        ),
        pytest.param(
            lambda args: [*_without_scope(args), "--strain", "phenotype"],
            "scope 'phenotype'",
            id="abbreviated",
        ),
        pytest.param(_without_scope, "does not state", id="missing"),
        pytest.param(
            lambda args: [*_without_scope(args), "--strain-rollup"],
            "expected one argument",
            id="no_value",
        ),
        pytest.param(
            lambda args: [*_without_scope(args), "--strain-rollup", "strains"],
            "invalid choice",
            id="invalid_value",
        ),
    ],
)
def test_the_plan_check_reads_the_ledger_scope_as_the_ledger_parses_it(
    context: ArmContext, edit: Callable[[list[str]], list[str]], message: str
) -> None:
    """An off core's plan cannot run the ledger under another scope by a later repeat, the
    equals form or an abbreviation; a new plan must state the scope, well formed (#881)."""
    commands = _ledger_args(build_commands(context), edit)

    with pytest.raises(ValueError, match=message):
        runtime._validate_plan(commands, context)


@pytest.mark.parametrize(
    "edit",
    [
        pytest.param(lambda args: [*args, "--strain-rollup", "off"], id="repeated"),
        pytest.param(lambda args: [*_without_scope(args), "--strain-rollup=off"], id="equals"),
    ],
)
def test_a_form_the_ledger_reads_as_the_core_scope_is_accepted(
    context: ArmContext, edit: Callable[[list[str]], list[str]]
) -> None:
    runtime._validate_plan(_ledger_args(build_commands(context), edit), context)


def test_only_a_plan_recorded_before_788_may_omit_the_scope(context: ArmContext) -> None:
    """Its ledger ran under off, on a core that recorded no scope (#788, #881)."""
    commands = _ledger_args(build_commands(context), _without_scope)

    runtime._validate_plan(commands, context, absent_scope_is_off=True)
    with pytest.raises(ValueError, match="does not state"):
        runtime._validate_plan(commands, replace(context, strain_rollup="phenotype"))
