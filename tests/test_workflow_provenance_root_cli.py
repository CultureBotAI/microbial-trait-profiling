"""Route explicit provenance roots through the supported workflows (#1780)."""

from __future__ import annotations

import argparse
import inspect
from collections.abc import Callable, Sequence
from dataclasses import fields, replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from mtp import integrated_evidence
from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow import comparison
from mtp.workflow.arms import __main__ as cli
from mtp.workflow.arms import ilyas, jed, noel, runtime
from mtp.workflow.arms.types import ArmContext, Command

Call = tuple[str, tuple[object, ...], dict[str, object]]


def _arguments(command: str, root: Path) -> list[str]:
    manifests = [str(root / "core run.json")]
    if command == "compare":
        manifests.extend(str(root / f"{arm} run.json") for arm in ("noel", "jed", "ilyas"))
    return [command, *manifests, str(root / "fresh output"), "--allow-temporary-output"]


def _accepted_cli(arguments: list[str]) -> int:
    try:
        return cli.main(arguments)
    except SystemExit as error:
        raise AssertionError("The supported workflow must accept the root option") from error


def _root_keyword(function: Callable[..., object]) -> None:
    parameter = inspect.signature(function).parameters.get("provenance_root")
    assert parameter is not None, "The execution API must expose the optional provenance root"
    assert parameter.kind == inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is None


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> list[Call]:
    observed: list[Call] = []

    def arm(*args: object, **kwargs: object) -> SimpleNamespace:
        observed.append(("arm", args, kwargs))
        return SimpleNamespace(manifest_path="arm_run.json")

    def compare(*args: object, **kwargs: object) -> SimpleNamespace:
        observed.append(("compare", args, kwargs))
        return SimpleNamespace(manifest_path="comparison_run.json")

    monkeypatch.setattr(cli, "execute_arm", arm)
    monkeypatch.setattr(cli, "execute_comparison_run", compare)
    return observed


@pytest.mark.parametrize("command", ["noel", "compare"])
def test_supported_cli_passes_the_complete_root_path_to_execution(
    tmp_path: Path, calls: list[Call], command: str
) -> None:
    root = tmp_path / "source root with spaces"
    assert _accepted_cli([*_arguments(command, tmp_path), "--provenance-root", str(root)]) == 0
    assert len(calls) == 1
    kind, positional, keywords = calls[0]
    assert keywords["provenance_root"] == root
    if command == "noel":
        assert kind == "arm"
        assert positional == (
            "noel",
            tmp_path / "core run.json",
            tmp_path / "fresh output",
            noel.build_commands,
        )
    else:
        assert kind == "compare"
        assert positional == tuple(
            tmp_path / filename
            for filename in (
                "core run.json",
                "noel run.json",
                "jed run.json",
                "ilyas run.json",
                "fresh output",
            )
        )


@pytest.mark.parametrize("command", ["ilyas", "noel", "jed", "compare"])
def test_omitted_option_preserves_existing_execution_defaults(
    tmp_path: Path, calls: list[Call], command: str
) -> None:
    assert _accepted_cli(_arguments(command, tmp_path)) == 0
    assert len(calls) == 1
    kind, positional, keywords = calls[0]
    assert keywords.get("provenance_root") is None
    previous_keywords = {key: value for key, value in keywords.items() if key != "provenance_root"}
    if command == "compare":
        assert kind == "compare"
        assert previous_keywords == {"arm1_pooled_stats": None}
    else:
        assert kind == "arm"
        assert positional[0] == command
        assert (
            positional[3]
            is {
                "ilyas": ilyas.build_commands,
                "noel": noel.build_commands,
                "jed": jed.build_commands,
            }[command]
        )
        assert previous_keywords == {"replicates": None, "seed": None, "figures": True}


@pytest.mark.parametrize("arm", ["jed", "ilyas"])
@pytest.mark.parametrize("spelling", ["separate", "equals", "abbreviated"])
def test_other_arm_clis_refuse_the_option_without_executing(
    tmp_path: Path,
    calls: list[Call],
    capsys: pytest.CaptureFixture[str],
    arm: str,
    spelling: str,
) -> None:
    output = tmp_path / "fresh output"
    root = str(tmp_path / "sources")
    option = {
        "separate": ["--provenance-root", root],
        "equals": ["--provenance-root=" + root],
        "abbreviated": ["--prov", root],
    }[spelling]
    with pytest.raises(SystemExit) as stopped:
        cli.main([*_arguments(arm, tmp_path), *option])
    assert stopped.value.code == 2
    assert "--provenance-root is supported only for the Noel arm" in capsys.readouterr().err
    assert calls == []
    assert not output.exists()


@pytest.mark.parametrize("function", [runtime.execute_arm, comparison.execute_comparison_run])
def test_execution_root_keywords_are_optional_and_keyword_only(
    function: Callable[..., object],
) -> None:
    _root_keyword(function)


@pytest.mark.parametrize("arm", ["jed", "ilyas"])
def test_other_arm_python_apis_refuse_before_core_or_planner_work(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, arm: str
) -> None:
    _root_keyword(runtime.execute_arm)
    touched: list[str] = []

    def forbidden_core(*_args: object, **_kwargs: object) -> None:
        touched.append("core")
        raise AssertionError("Unsupported root must be refused before core verification")

    def forbidden_plan(_context: ArmContext) -> list[Command]:
        touched.append("planner")
        raise AssertionError("Unsupported root must be refused before planning")

    monkeypatch.setattr(runtime, "open_core_bundle", forbidden_core)
    output = tmp_path / "fresh output"
    with pytest.raises(ValueError, match=r"^--provenance-root is supported only for the Noel arm$"):
        runtime.execute_arm(
            arm,
            tmp_path / "missing-core.json",
            output,
            forbidden_plan,
            provenance_root=tmp_path,
        )
    assert touched == []
    assert not output.exists()


@pytest.fixture
def context(tmp_path: Path) -> ArmContext:
    products = {
        name: tmp_path / "core" / filename
        for name, filename in {
            "count_profiles": "study_trait_profiles_wide.tsv",
            "ora_profiles": "study_trait_profiles_ora_wide.tsv",
            "experiments": "experiments.tsv",
            "raw_dump": "full_dump.csv",
            "feature_metadata": "taxa_trait_matrix_columns.tsv",
            "kg_nodes": "merged-kg_nodes.tsv",
        }.items()
    }
    return ArmContext(
        REPO_ROOT,
        tmp_path / "core/run.json",
        tmp_path / "noel",
        tmp_path / "noel/inputs",
        products,
        "a" * 40,
        None,
        None,
        True,
    )


def _require_context_option() -> None:
    options = fields(ArmContext)
    assert options[-1].name == "provenance_root", "Append the option after existing context fields"
    assert options[-1].default is None


def test_default_noel_plan_has_no_root_flag_and_explicit_none_is_identical(
    context: ArmContext,
) -> None:
    _require_context_option()
    baseline = noel.build_commands(context)
    assert noel.build_commands(replace(context, provenance_root=None)) == baseline
    assert all("--provenance-root" not in command.argv for command in baseline)


def test_root_changes_only_the_integrated_child_arguments(context: ArmContext) -> None:
    _require_context_option()
    root = context.output.parent / "source root with spaces"
    baseline = noel.build_commands(context)
    explicit = noel.build_commands(replace(context, provenance_root=root))
    assert [command.id for command in explicit] == [command.id for command in baseline]
    restored: list[Command] = []
    integrated = 0
    for command in explicit:
        if command.id == "integrated_evidence":
            integrated += 1
            assert command.argv.count("--provenance-root") == 1
            index = command.argv.index("--provenance-root")
            assert command.argv[index + 1] == str(root)
            assert root not in command.inputs
            restored.append(replace(command, argv=command.argv[:index] + command.argv[index + 2 :]))
        else:
            assert "--provenance-root" not in command.argv
            restored.append(command)
    assert integrated == 1
    assert restored == baseline


class _Parsed(Exception):
    def __init__(self, namespace: argparse.Namespace) -> None:
        super().__init__()
        self.namespace = namespace


def test_real_integrated_parser_accepts_the_explicit_noel_root(
    context: ArmContext, monkeypatch: pytest.MonkeyPatch
) -> None:
    _require_context_option()
    root = context.output.parent / "source root with spaces"
    command = next(
        item
        for item in noel.build_commands(replace(context, provenance_root=root))
        if item.id == "integrated_evidence"
    )
    parse_args = argparse.ArgumentParser.parse_args

    def capture(
        parser: argparse.ArgumentParser,
        args: Sequence[str] | None = None,
        namespace: argparse.Namespace | None = None,
    ) -> argparse.Namespace:
        raise _Parsed(parse_args(parser, args, namespace))

    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", capture)
    with pytest.raises(_Parsed) as parsed:
        integrated_evidence.main(list(command.argv[3:]))
    assert parsed.value.namespace.provenance_root == root


@pytest.mark.parametrize("explicit", [False, True], ids=["omitted", "explicit"])
def test_comparison_forwards_its_root_to_the_real_integrated_parser(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, explicit: bool
) -> None:
    """Pin argument routing only; verified execution is covered by the real E2E suite."""
    _root_keyword(comparison.execute_comparison_run)
    sources = tmp_path / "source root"
    sources.mkdir()
    paths = {name: sources / f"{name}.json" for name in ("core", "noel", "jed", "ilyas")}
    for path in paths.values():
        path.write_text("fixture bytes for routing only\n", encoding="utf-8")
    pooled = sources / "pooled statistics.tsv"
    pooled.write_text("fixture statistics for routing only\n", encoding="utf-8")
    products = {
        "count_profiles": sources / "study_trait_profiles_wide.tsv",
        "ora_profiles": sources / "study_trait_profiles_ora_wide.tsv",
    }
    for path in products.values():
        path.write_text("fixture profile bytes for routing only\n", encoding="utf-8")
    core = SimpleNamespace(
        run_id="fixture-core",
        contrast_population="all_paired",
        direction_basis="group_1_vs_group_0",
        artifacts=[
            SimpleNamespace(product_id=name, path=str(path.resolve()))
            for name, path in products.items()
        ],
    )
    arms = {
        path: SimpleNamespace(
            arm=name,
            core_run_id=core.run_id,
            core_manifest_sha256=file_sha256(paths["core"]),
            contrast_population=core.contrast_population,
            direction_basis=core.direction_basis,
        )
        for name, path in paths.items()
        if name != "core"
    }
    monkeypatch.setattr(comparison, "verify_core_run", lambda _path: core)
    monkeypatch.setattr(comparison, "verify_arm", lambda path: arms[path])
    monkeypatch.setattr(comparison, "repository_state", lambda _root: ("a" * 40, False))
    monkeypatch.setattr(
        comparison, "recorded_output", lambda arm, step, name: sources / arm.arm / step / name
    )
    for name in (
        "noel_calibration_verdict",
        "jed_calibration_verdict",
        "shared_basis_ledger",
        "core_analysis_bases",
    ):
        monkeypatch.setattr(comparison, name, lambda *_args: None)

    async def forbidden_command(*_args: object) -> None:
        raise AssertionError("Manual pooled statistics must not launch a generator")

    monkeypatch.setattr(comparison, "_exec_command", forbidden_command)
    parse_args = argparse.ArgumentParser.parse_args

    def capture(
        parser: argparse.ArgumentParser,
        args: Sequence[str] | None = None,
        namespace: argparse.Namespace | None = None,
    ) -> argparse.Namespace:
        raise _Parsed(parse_args(parser, args, namespace))

    real_integrated_main = integrated_evidence.main

    def capture_main(argv: Sequence[str] | None = None) -> int:
        with monkeypatch.context() as producer_patch:
            producer_patch.setattr(argparse.ArgumentParser, "parse_args", capture)
            return real_integrated_main(argv)

    monkeypatch.setattr(integrated_evidence, "main", capture_main)
    supplied = sources / "unused" / ".." if explicit else None
    output = tmp_path / "comparison output"
    try:
        with pytest.raises(_Parsed) as parsed:
            comparison.execute_comparison_run(
                paths["core"],
                paths["noel"],
                paths["jed"],
                paths["ilyas"],
                output,
                arm1_pooled_stats=pooled,
                provenance_root=supplied,
            )
    except (ValueError, OSError) as error:
        raise AssertionError("Valid argument routing must reach the integrated parser") from error
    arguments = parsed.value.namespace
    assert arguments.provenance_root == (supplied.resolve() if supplied is not None else None)
    assert arguments.core_run == paths["core"].resolve()
    assert arguments.arm1_pooled_stats == pooled.resolve()
    assert arguments.output_dir == output.resolve() / "integrated_evidence"
    assert arguments.require_complete_sources is True
    assert list(output.iterdir()) == []


class _Planned(Exception):
    """Stop the real runtime at its planner, before execution or record creation."""


@pytest.mark.parametrize("explicit", [False, True], ids=["omitted", "explicit"])
def test_noel_execution_passes_the_normalized_root_to_its_context(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, context: ArmContext, explicit: bool
) -> None:
    _root_keyword(runtime.execute_arm)
    _require_context_option()
    bundle = SimpleNamespace(
        run=SimpleNamespace(
            artifacts=[SimpleNamespace(product_id=name) for name in context.products],
            parameters=[],
            contrast_population="all_paired",
            direction_basis="group_1_vs_group_0",
        ),
        path=lambda product: context.products[product],
    )
    for path in (context.manifest, *context.products.values()):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("synthetic routing input\n", encoding="utf-8")
    monkeypatch.setattr(runtime, "open_core_bundle", lambda *_args, **_kwargs: bundle)
    monkeypatch.setattr(runtime, "preflight_consumer_outputs", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(runtime, "repository_state", lambda _root: ("a" * 40, False))
    supplied = tmp_path / "source root" / ".." if explicit else None
    captured: list[ArmContext] = []

    def plan(context: ArmContext) -> list[Command]:
        captured.append(context)
        raise _Planned

    output = tmp_path / "fresh output"
    try:
        with pytest.raises(_Planned):
            runtime.execute_arm(
                "noel",
                context.manifest,
                output,
                plan,
                provenance_root=supplied,
            )
    except (ValueError, OSError) as error:
        raise AssertionError("Valid Noel root routing must reach its planner") from error
    assert len(captured) == 1
    expected = supplied.resolve() if supplied is not None else None
    assert captured[0].provenance_root == expected
    assert not output.exists()
    command = next(
        item for item in noel.build_commands(captured[0]) if item.id == "integrated_evidence"
    )
    parse_args = argparse.ArgumentParser.parse_args

    def capture(
        parser: argparse.ArgumentParser,
        args: Sequence[str] | None = None,
        namespace: argparse.Namespace | None = None,
    ) -> argparse.Namespace:
        raise _Parsed(parse_args(parser, args, namespace))

    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", capture)
    with pytest.raises(_Parsed) as parsed:
        integrated_evidence.main(list(command.argv[3:]))
    assert parsed.value.namespace.provenance_root == expected
