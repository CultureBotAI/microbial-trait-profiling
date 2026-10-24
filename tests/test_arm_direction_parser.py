"""Issue967: verifier parsing follows actual planned figure entrypoints."""

from __future__ import annotations

import argparse
import runpy
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT
from mtp.workflow.arms import jed, runtime
from mtp.workflow.arms.types import ArmContext, Command

BASES = ("group_1_vs_group_0", "disease_vs_comparison")


@pytest.fixture(params=BASES)
def figure_commands(tmp_path: Path, request: pytest.FixtureRequest) -> list[Command]:
    products = {
        name: tmp_path / "core" / f"{name}.tsv"
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
    context = ArmContext(
        repository=REPO_ROOT,
        manifest=tmp_path / "core/run.json",
        output=tmp_path / "figures",
        inputs=tmp_path / "inputs",
        products=products,
        code_commit="fixture",
        replicates=None,
        seed=None,
        figures=True,
        direction_basis=request.param,
    )
    commands = [
        command for command in jed.build_commands(context) if "--direction-basis" in command.argv
    ]
    assert len(commands) == 2
    return commands


class _Parsed(Exception):
    """Stop the real main immediately after parsing and before any data access."""


def _actual_arguments(command: Command, monkeypatch: pytest.MonkeyPatch) -> argparse.Namespace:
    parsed = []
    original = argparse.ArgumentParser.parse_args

    def capture(
        parser: argparse.ArgumentParser,
        args: list[str] | None = None,
        namespace: argparse.Namespace | None = None,
    ) -> argparse.Namespace:
        parsed.append(original(parser, args, namespace))
        raise _Parsed

    with monkeypatch.context() as temporary:
        temporary.setattr(argparse.ArgumentParser, "parse_args", capture)
        temporary.setattr(sys, "argv", list(command.argv[1:]))
        with pytest.raises(_Parsed):
            runpy.run_path(command.argv[1], run_name="__main__")
    return parsed[0]


def test_every_direction_aware_jed_figure_has_registered_parser(
    figure_commands: list[Command],
) -> None:
    for command in figure_commands:
        relative = Path(command.argv[1]).relative_to(REPO_ROOT).as_posix()
        assert relative in runtime.DIRECTION_REQUIRED_SCRIPTS


@pytest.mark.parametrize("index", [0, 1])
@pytest.mark.parametrize("form", ["split", "equals", "abbreviated", "repeated"])
def test_real_figure_parser_and_verifier_agree_on_effective_basis(
    figure_commands: list[Command], index: int, form: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    command = figure_commands[index]
    argv = list(command.argv)
    position = argv.index("--direction-basis")
    basis = argv[position + 1]
    opposite = next(value for value in BASES if value != basis)
    options = {
        "split": ["--direction-basis", basis],
        "equals": [f"--direction-basis={basis}"],
        "abbreviated": ["--direction-b", basis],
        "repeated": ["--direction-basis", opposite, "--direction-basis", basis],
    }[form]
    argv[position : position + 2] = options
    changed = replace(command, argv=tuple(argv))
    actual = _actual_arguments(changed, monkeypatch)
    assert actual.direction_basis == basis
    assert runtime._step_direction_basis(changed, REPO_ROOT) == basis


@pytest.mark.parametrize("index", [0, 1])
@pytest.mark.parametrize("suffix", [["--bogus-option"], ["--", "--direction-basis", "other"]])
def test_verifier_rejects_arguments_rejected_by_real_figure_parser(
    figure_commands: list[Command],
    index: int,
    suffix: list[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    command = figure_commands[index]
    changed = replace(command, argv=(*command.argv, *suffix))
    with pytest.raises(SystemExit) as actual:
        _actual_arguments(changed, monkeypatch)
    assert actual.value.code == 2
    with pytest.raises(ValueError, match="unrecognized arguments"):
        runtime._step_direction_basis(changed, REPO_ROOT)


@pytest.mark.parametrize("index", [0, 1])
@pytest.mark.parametrize("tamper", ["missing", "abbreviated_opposite"])
def test_case_changed_registered_script_cannot_bypass_direction_validation(
    figure_commands: list[Command], index: int, tamper: str
) -> None:
    command = figure_commands[index]
    argv = list(command.argv)
    script = Path(argv[1])
    # Works on case-sensitive CI too: registry identity must not depend on the FS.
    argv[1] = str(script.with_name(script.name.replace("create_", "Create_", 1)))
    position = argv.index("--direction-basis")
    basis = argv[position + 1]
    argv[position : position + 2] = (
        [] if tamper == "missing" else ["--direction-b", next(b for b in BASES if b != basis)]
    )
    with pytest.raises(ValueError, match="canonical script path"):
        runtime._step_direction_basis(replace(command, argv=tuple(argv)), REPO_ROOT)
