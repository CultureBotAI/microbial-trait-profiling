"""Issue967: verifier parsing follows actual planned figure entrypoints."""

from __future__ import annotations

import argparse
import runpy
import sys
from dataclasses import replace
from pathlib import Path
from types import ModuleType
from typing import Any

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
            "kg_edges",
            "kg_nodes",
            "source_identity",
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
    assert len(commands) == 3
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

    had_sibling = "create_study_trait_profile_heatmap" in sys.modules
    prior_sibling = sys.modules.get("create_study_trait_profile_heatmap")
    try:
        with monkeypatch.context() as temporary:
            # A real `python path/to/script.py` puts that directory first on sys.path.
            temporary.syspath_prepend(str(Path(command.argv[1]).parent))
            temporary.setattr(argparse.ArgumentParser, "parse_args", capture)
            temporary.setattr(sys, "argv", list(command.argv[1:]))
            with pytest.raises(_Parsed):
                runpy.run_path(command.argv[1], run_name="__main__")
    finally:
        if not had_sibling:
            sys.modules.pop("create_study_trait_profile_heatmap", None)
        elif prior_sibling is not None:
            sys.modules["create_study_trait_profile_heatmap"] = prior_sibling
    return parsed[0]


def test_every_direction_aware_jed_figure_has_registered_parser(
    figure_commands: list[Command],
) -> None:
    for command in figure_commands:
        relative = Path(command.argv[1]).relative_to(REPO_ROOT).as_posix()
        assert relative in runtime.DIRECTION_REQUIRED_SCRIPTS | runtime.FULL_PARSER_SCRIPTS


@pytest.mark.parametrize("index", [0, 1, 2])
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


@pytest.mark.parametrize("index", [0, 1, 2])
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


@pytest.mark.parametrize("index", [0, 1, 2])
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


@pytest.mark.parametrize("index", [0, 1, 2])
def test_verifier_rejects_help_without_printing(
    figure_commands: list[Command], index: int, capsys: pytest.CaptureFixture[str]
) -> None:
    command = figure_commands[index]
    changed = replace(command, argv=(*command.argv, "--help"))
    with pytest.raises(ValueError, match="argument parsing stopped with status 0"):
        runtime._step_direction_basis(changed, REPO_ROOT)
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


@pytest.mark.parametrize("direction", ["all", "inc", "dec"])
def test_condition_module_uses_its_execution_parser(
    direction: str, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    command = Command(
        "condition_module",
        (
            sys.executable,
            str(REPO_ROOT / "scripts/create_condition_trait_module_heatmap.py"),
            "--direction",
            direction,
        ),
        (),
        (),
    )
    scripts_dir = (REPO_ROOT / "scripts").resolve()
    with monkeypatch.context() as isolated:
        isolated.delitem(sys.modules, "create_study_trait_profile_heatmap", raising=False)
        isolated.setattr(
            sys,
            "path",
            [entry for entry in sys.path if Path(entry or ".").resolve() != scripts_dir],
        )
        try:
            assert runtime._step_direction_basis(command, REPO_ROOT) is None
        except ModuleNotFoundError as exc:
            raise AssertionError(
                f"condition_module parser required sibling import on load: {exc}"
            ) from exc
        assert "create_study_trait_profile_heatmap" not in sys.modules
    assert _actual_arguments(command, monkeypatch).direction == direction
    invalid = replace(command, argv=(*command.argv, "--direction-basis", BASES[0]))
    with pytest.raises(SystemExit) as actual:
        _actual_arguments(invalid, monkeypatch)
    assert actual.value.code == 2
    capsys.readouterr()
    with pytest.raises(ValueError, match="unrecognized arguments"):
        runtime._step_direction_basis(invalid, REPO_ROOT)
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_condition_module_parser_loads_without_sibling_import_regardless_of_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#1023: _step_direction_basis loads condition_module parser without sibling import."""
    script_path = REPO_ROOT / "scripts/create_condition_trait_module_heatmap.py"
    scripts_dir = script_path.parent.resolve()
    command = Command(
        "condition_module",
        (
            sys.executable,
            str(script_path),
            "--direction",
            "inc",
        ),
        (),
        (),
    )

    # Simulate prior test execution that cached the sibling module and scripts/ on sys.path.
    sentinel_sibling = ModuleType("create_study_trait_profile_heatmap")
    monkeypatch.setitem(sys.modules, "create_study_trait_profile_heatmap", sentinel_sibling)
    monkeypatch.syspath_prepend(str(scripts_dir))
    assert sys.modules.get("create_study_trait_profile_heatmap") is sentinel_sibling
    assert any(Path(entry or ".").resolve() == scripts_dir for entry in sys.path)

    with monkeypatch.context() as isolated:
        isolated.delitem(sys.modules, "create_study_trait_profile_heatmap", raising=False)
        isolated.setattr(
            sys,
            "path",
            [entry for entry in sys.path if Path(entry or ".").resolve() != scripts_dir],
        )
        assert "create_study_trait_profile_heatmap" not in sys.modules
        assert all(Path(entry or ".").resolve() != scripts_dir for entry in sys.path)
        try:
            basis = runtime._step_direction_basis(command, REPO_ROOT)
        except ModuleNotFoundError as exc:
            raise AssertionError(
                f"condition_module parser required sibling script import on load: {exc}"
            ) from exc
        assert basis is None
        assert "create_study_trait_profile_heatmap" not in sys.modules
        assert all(Path(entry or ".").resolve() != scripts_dir for entry in sys.path)

        # Also verify that --include-media-traits fails at parser.error() in main()
        # before importing create_study_trait_profile_heatmap.
        namespace = runpy.run_path(str(script_path))
        assert "create_study_trait_profile_heatmap" not in sys.modules
        isolated.setattr(sys, "argv", [str(script_path), "--include-media-traits"])
        try:
            with pytest.raises(SystemExit) as media_exit:
                namespace["main"]()
        except ModuleNotFoundError as exc:
            raise AssertionError(
                f"main() imported sibling before --include-media-traits check: {exc}"
            ) from exc
        assert media_exit.value.code == 2
        assert "create_study_trait_profile_heatmap" not in sys.modules

        # Verify that main() itself still defers and invokes load_raw_experiment_metadata
        # after argument parsing when execution proceeds past the --include-media-traits check.
        called_with: list[Any] = []

        class _StopAfterSiblingCall(Exception):
            pass

        fake_sibling: Any = ModuleType("create_study_trait_profile_heatmap")

        def _fake_load_raw_experiment_metadata(path: Path) -> Any:
            called_with.append(path)
            raise _StopAfterSiblingCall

        fake_sibling.load_raw_experiment_metadata = _fake_load_raw_experiment_metadata
        sentinel_dump = tmp_path / "sentinel_dump.csv"
        isolated.setitem(sys.modules, "create_study_trait_profile_heatmap", fake_sibling)
        isolated.setattr(
            sys,
            "argv",
            [str(script_path), "--full-dump-input", str(sentinel_dump)],
        )
        try:
            namespace["main"]()
        except _StopAfterSiblingCall:
            pass
        except Exception as exc:
            raise AssertionError(
                f"main() did not call deferred load_raw_experiment_metadata: {exc}"
            ) from exc
        assert called_with == [sentinel_dump]

    # Confirm outer state was restored after the isolated context exited.
    assert sys.modules.get("create_study_trait_profile_heatmap") is sentinel_sibling
    assert any(Path(entry or ".").resolve() == scripts_dir for entry in sys.path)
