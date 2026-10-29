"""The Jed command plan must route real CLI parsers to verified inputs (#658)."""

from __future__ import annotations

import argparse
import importlib
import runpy
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from mtp.workflow.arms.jed import build_commands
from mtp.workflow.arms.shared import human_count_path
from mtp.workflow.arms.types import ArmContext, Command

REPOSITORY = Path(__file__).resolve().parents[1]


@pytest.fixture
def context(tmp_path: Path) -> ArmContext:
    products = {
        name: tmp_path / "immutable core" / f"{name}.tsv"
        for name in (
            "raw_dump",
            "count_profiles",
            "experiments",
            "feature_metadata",
            "reported_taxa",
            "taxa_matrix",
            "study_table",
        )
    }
    return ArmContext(
        repository=REPOSITORY,
        manifest=tmp_path / "immutable core" / "run.json",
        output=tmp_path / "arm outputs",
        inputs=tmp_path / "input view",
        products=products,
        code_commit="test-revision",
        replicates=None,
        seed=None,
        figures=True,
    )


class _Parsed(Exception):
    """Stop a real command immediately after argument parsing, before data I/O."""


def _parse_command(command: Command, monkeypatch: pytest.MonkeyPatch) -> argparse.Namespace:
    parsed: list[argparse.Namespace] = []
    original = argparse.ArgumentParser.parse_args

    def capture(
        parser: argparse.ArgumentParser,
        args: list[str] | None = None,
        namespace: argparse.Namespace | None = None,
    ) -> argparse.Namespace:
        parsed.append(original(parser, args, namespace))
        raise _Parsed

    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", capture)
    monkeypatch.syspath_prepend(str(REPOSITORY / "scripts"))
    if command.argv[1] == "-m":
        monkeypatch.setattr(sys, "argv", [command.argv[2], *command.argv[3:]])
        with pytest.raises(_Parsed):
            importlib.import_module(command.argv[2]).main()
    else:
        monkeypatch.setattr(sys, "argv", list(command.argv[1:]))
        with pytest.raises(_Parsed):
            runpy.run_path(command.argv[1], run_name="__main__")
    return parsed[0]


@pytest.mark.parametrize("command_index", range(14))
def test_actual_cli_paths_are_explicit_and_declared(
    context: ArmContext, command_index: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An added/defaulted CLI path must not silently read or write the checkout."""
    command = build_commands(context)[command_index]
    parsed = _parse_command(command, monkeypatch)
    for value in vars(parsed).values():
        if not isinstance(value, Path):
            continue
        assert value in command.inputs or value in command.outputs, (command.id, value)
        assert str(value) in command.argv, (command.id, "implicit path", value)
    assert command.argv[0] == sys.executable
    assert all(path.is_relative_to(context.output) for path in command.outputs)
    assert not set(command.inputs) & set(command.outputs)


def test_calibration_retains_existing_defaults_and_honors_overrides(context: ArmContext) -> None:
    from mtp.composition_calibration.composition import build_parser

    baseline = next(c for c in build_commands(context) if c.id == "composition_calibration")
    parsed = build_parser().parse_args(baseline.argv[3:])
    assert (parsed.replicates, parsed.seed) == (100, 182)
    custom = replace(context, replicates=3, seed=77)
    overridden = next(c for c in build_commands(custom) if c.id == "composition_calibration")
    parsed = build_parser().parse_args(overridden.argv[3:])
    assert (parsed.replicates, parsed.seed) == (3, 77)


def test_figures_are_explicitly_optional_and_scoped(context: ArmContext) -> None:
    commands = build_commands(replace(context, figures=False))
    assert {command.id for command in commands} == {
        "human_counts",
        "composition",
        "composition_sensitivity",
        "composition_calibration",
        "density_adjustment",
        "trait_coupling",
    }
    outputs = [command.outputs[0] for command in build_commands(context)]
    assert len(outputs) == len(set(outputs))
    for left in outputs:
        for right in outputs:
            if left != right:
                assert not left.is_relative_to(right)


@pytest.mark.parametrize("direction_basis", ["group_1_vs_group_0", "disease_vs_comparison"])
@pytest.mark.parametrize("command_id", ["taxa_guild_heatmap", "profile_heatmap"])
def test_direction_heatmaps_receive_the_core_direction_basis(
    context: ArmContext, direction_basis: str, command_id: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Both plotting CLIs receive the core's directions (#795, #842)."""
    commands = build_commands(replace(context, direction_basis=direction_basis))
    command = next(command for command in commands if command.id == command_id)
    assert _parse_command(command, monkeypatch).direction_basis == direction_basis


def test_profiles_are_count_products_and_sensitivity_protocol_is_bound(context: ArmContext) -> None:
    commands = {command.id: command for command in build_commands(context)}
    for name in ("composition", "guild_maps", "table1_traits", "human_counts"):
        assert context.products["count_profiles"] in commands[name].inputs
    for name in (
        "composition_sensitivity",
        "composition_calibration",
        "density_adjustment",
        "trait_coupling",
        "taxa_guild_heatmap",
        "profile_heatmap",
        "condition_trait_heatmap",
        "taxa_coverage",
    ):
        assert human_count_path(context) in commands[name].inputs
    assert (REPOSITORY / "methods" / "table1_log_ratio_sensitivity_protocol.json") in commands[
        "composition_sensitivity"
    ].inputs
    for name in ("composition", "table1_traits", "guild_maps"):
        command = commands[name]
        assert command.argv[command.argv.index("--core-run") + 1] == str(context.manifest)
        assert command.outputs[0] / "core_consumption.json" in command.outputs
