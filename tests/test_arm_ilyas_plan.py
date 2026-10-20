"""Ilyas reruns preserve scientific CLI settings and never read stale figure inputs (#658)."""

from __future__ import annotations

import runpy
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from mtp.ora_calibration.harness import build_parser
from mtp.paths import REPO_ROOT
from mtp.stratification import parse_args as parse_strata
from mtp.workflow.arms.ilyas import (
    FIGURE_CORE_FILES,
    FIGURE_FILES,
    FIGURE_RESULT_FILES,
    build_commands,
    main,
    materialize_figure_inputs,
)
from mtp.workflow.arms.shared import human_count_path
from mtp.workflow.arms.types import ArmContext


@pytest.fixture
def context(tmp_path: Path) -> ArmContext:
    return ArmContext(
        repository=REPO_ROOT,
        manifest=tmp_path / "core" / "run.json",
        output=tmp_path / "ilyas results",
        inputs=tmp_path / "core inputs",
        products={
            "count_profiles": tmp_path / "core" / "counts.tsv",
            "experiments": tmp_path / "core" / "experiments.tsv",
            "ora_profiles": tmp_path / "core" / "ora.tsv",
            "taxa_matrix": tmp_path / "core" / "matrix.tsv",
        },
        code_commit="abc123",
        replicates=None,
        seed=None,
        figures=True,
    )


def test_ora_analysis_and_rankings_use_existing_parsers(context: ArmContext) -> None:
    enrichment_parser = runpy.run_path(str(REPO_ROOT / "scripts/trait_enrichment.py"))["parse_args"]
    ranking_parser = runpy.run_path(str(REPO_ROOT / "scripts/export_trait_enrichment_rankings.py"))[
        "parse_args"
    ]
    commands = build_commands(context)
    for index, tail, name in ((0, "upper", "enrichment"), (2, "lower", "depletion")):
        command = commands[index]
        args = enrichment_parser(["--tail", tail, *command.argv[2:]])
        assert args.core_run == context.manifest
        assert args.input is None
        assert args.taxa_matrix is None
        assert args.fdr
        assert args.fdr_scope == "global"
        assert args.p_threshold == 0.001
        assert args.min_fold == 1.5 if tail == "upper" else args.max_fold == 0.5
        assert args.output == context.output / f"trait_{name}_report.txt"
        assert args.tsv_output in command.outputs
        assert args.all_results_output in command.outputs
        ranking = commands[index + 1]
        rank_args = ranking_parser(list(ranking.argv[2:]))
        assert rank_args.input == rank_args.results == args.all_results_output
        assert rank_args.unfiltered
        assert rank_args.tail == tail
        assert rank_args.fdr_scope == "global"
        assert ranking.inputs == (args.all_results_output,)
        assert rank_args.overall_output in ranking.outputs
        assert rank_args.pair_output in ranking.outputs


def test_strata_and_calibration_preserve_recipe_defaults(context: ArmContext) -> None:
    commands = {command.id: command for command in build_commands(context)}
    strata = parse_strata(list(commands["strata"].argv[2:]))
    calibration = build_parser().parse_args(commands["calibration"].argv[3:])
    assert strata.strata == ["pooled", "feces", "feces_excluded"]
    assert strata.fdr
    assert strata.fdr_scope == "global"
    assert strata.p_threshold == 0.001
    assert strata.min_fold == 1.5
    assert (
        strata.profiles
        == calibration.profiles
        == context.inputs / ("study_trait_profiles_ora_wide.tsv")
    )
    assert strata.taxa_matrix == calibration.taxa_matrix == context.inputs / "taxa_trait_matrix.tsv"
    assert strata.experiments == context.inputs / "experiments.tsv"
    assert len(commands["strata"].outputs) == 20
    assert calibration.replicates == 10
    assert calibration.seed == 42
    assert calibration.scenario == "all"
    assert calibration.test_type == "enrichment"
    assert calibration.alpha == calibration.fdr_q == 0.05
    assert calibration.signal_strength == 0.6
    assert calibration.n_spike_traits == 3
    assert calibration.min_fold == 1.5
    assert calibration.max_fold == 0.67


def test_calibration_overrides_and_no_figures(context: ArmContext) -> None:
    commands = build_commands(replace(context, replicates=3, seed=0, figures=False))
    assert [command.id for command in commands] == [
        "enrichment",
        "enrichment_rankings",
        "depletion",
        "depletion_rankings",
        "strata",
        "calibration",
    ]
    args = build_parser().parse_args(commands[-1].argv[3:])
    assert args.replicates == 3
    assert args.seed == 0


def test_all_inputs_outputs_are_bound_and_figures_use_fresh_density(context: ArmContext) -> None:
    commands = build_commands(context)
    available = {context.manifest, *context.products.values()}
    available.update(context.inputs / name for name in (*FIGURE_CORE_FILES,))
    available.add(context.inputs / "study_trait_profiles_wide.tsv")
    produced: set[Path] = set()
    for command in commands:
        assert command.argv[0] == sys.executable
        assert all(path in available for path in command.inputs), command.id
        assert all(path.is_relative_to(context.output) for path in command.outputs)
        assert not produced.intersection(command.outputs)
        produced.update(command.outputs)
        available.update(command.outputs)
        assert not any("data/derived" in arg for arg in command.argv)
    by_id = {command.id: command for command in commands}
    density = context.output / "density_adjustment/table1_density_adjusted_effects.tsv"
    assert density in by_id["figure_density_adjustment"].outputs
    assert density in by_id["figure_inputs"].inputs
    assert human_count_path(context) in by_id["figure_density_adjustment"].inputs
    assert by_id["figures"].inputs == by_id["figure_inputs"].outputs
    assert {path.name for path in by_id["figures"].outputs} == set(FIGURE_FILES)


def _figure_sources(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    for name in FIGURE_CORE_FILES:
        (input_dir / name).write_text(f"current core: {name}\n")
    results = tuple(tmp_path / f"result-{index}.tsv" for index in range(3))
    for path in results:
        path.write_text(f"current result: {path.name}\n")
    return input_dir, results[0], results[1], results[2]


def test_figure_view_cli_copies_only_explicit_sources(tmp_path: Path) -> None:
    sources = _figure_sources(tmp_path)
    output = tmp_path / "view"
    # The execution runner creates declared output parents before invoking a command.
    (output / "density_adjustment").mkdir(parents=True)
    main(
        [
            "--input-dir",
            str(sources[0]),
            "--enrichment",
            str(sources[1]),
            "--depletion",
            str(sources[2]),
            "--density-effects",
            str(sources[3]),
            "--output-dir",
            str(output),
        ]
    )
    assert {str(path.relative_to(output)) for path in output.rglob("*") if path.is_file()} == {
        *FIGURE_CORE_FILES,
        *FIGURE_RESULT_FILES,
    }
    for name in FIGURE_CORE_FILES:
        assert (output / name).read_bytes() == (sources[0] / name).read_bytes()
    for source, name in zip(sources[1:], FIGURE_RESULT_FILES, strict=True):
        assert (output / name).read_bytes() == source.read_bytes()


def test_missing_source_rejects_before_any_figure_copy(tmp_path: Path) -> None:
    sources = _figure_sources(tmp_path)
    sources[-1].unlink()
    output = tmp_path / "view"
    with pytest.raises(FileNotFoundError, match="Required figure source"):
        materialize_figure_inputs(*sources, output)
    assert not output.exists()


@pytest.mark.parametrize("symlink", [False, True])
def test_stale_or_aliased_figure_view_is_rejected(tmp_path: Path, symlink: bool) -> None:
    sources = _figure_sources(tmp_path)
    output = tmp_path / "view"
    if symlink:
        output.symlink_to(sources[0], target_is_directory=True)
    else:
        output.mkdir()
        (output / "stale.tsv").write_text("stale data\n")
    with pytest.raises(FileExistsError, match="previous files"):
        materialize_figure_inputs(*sources, output)
