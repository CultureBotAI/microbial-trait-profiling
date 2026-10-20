"""Explicit Ilyas ORA rerun commands bound to one verified core execution (#658)."""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from mtp.workflow.arms.shared import human_count_command, human_count_path
from mtp.workflow.arms.types import ArmContext, Command

FIGURE_CORE_FILES = (
    "study_trait_profiles_ora_wide.tsv",
    "experiments.tsv",
    "taxa_trait_matrix.tsv",
    "taxa_trait_matrix_columns.tsv",
    "bugsigdb_taxa.tsv",
)
FIGURE_RESULT_FILES = (
    "trait_enrichment_results.tsv",
    "trait_depletion_results.tsv",
    "density_adjustment/table1_density_adjusted_effects.tsv",
)
FIGURE_FILES = (
    "figure1_condition_anatomy_ora_heatmap.png",
    "figure2_taxon_calibrated_ora_infographic.png",
    "figure3_lead_findings_ora_evidence.png",
    "figure4_annotation_density_pc1_landscape.png",
    "figure5_pd_crc_taxon_attribution_bubble.png",
)


def _script(context: ArmContext, name: str) -> tuple[str, ...]:
    return sys.executable, str(context.repository / "scripts" / name)


def _ora_commands(context: ArmContext, direction: str) -> list[Command]:
    prefix = f"trait_{direction}"
    report = context.output / f"{prefix}_report.txt"
    results = context.output / f"{prefix}_results.tsv"
    all_results = context.output / f"{prefix}_all_results.tsv"
    overall = context.output / f"{prefix}_traits_overall_ranked.tsv"
    pairs = context.output / f"{prefix}_trait_study_pairs_ranked.tsv"
    tail = "upper" if direction == "enrichment" else "lower"
    fold = ("--min-fold", "1.5") if direction == "enrichment" else ("--max-fold", "0.5")
    binding = report.with_suffix(".core_consumption.json")
    return [
        Command(
            id=direction,
            argv=(
                *_script(context, f"trait_{direction}.py"),
                "--core-run",
                str(context.manifest),
                "--output",
                str(report),
                "--tsv-output",
                str(results),
                "--all-results-output",
                str(all_results),
                "--p-threshold",
                "0.001",
                "--fdr-scope",
                "global",
                *fold,
            ),
            inputs=(
                context.manifest,
                context.products["ora_profiles"],
                context.products["taxa_matrix"],
            ),
            outputs=(report, results, all_results, binding, binding.with_suffix(".coverage.tsv")),
        ),
        Command(
            id=f"{direction}_rankings",
            argv=(
                *_script(context, "export_trait_enrichment_rankings.py"),
                "--input",
                str(all_results),
                "--results",
                str(all_results),
                "--all-results",
                "--fdr-scope",
                "global",
                "--tail",
                tail,
                "--overall-output",
                str(overall),
                "--pair-output",
                str(pairs),
            ),
            inputs=(all_results,),
            outputs=(overall, pairs),
        ),
    ]


def _figure_commands(context: ArmContext) -> list[Command]:
    density = context.output / "density_adjustment"
    density_effects = density / "table1_density_adjusted_effects.tsv"
    view = context.output / "figure_inputs"
    figures = context.output / "deck_figures"
    enrichment = context.output / "trait_enrichment_results.tsv"
    depletion = context.output / "trait_depletion_results.tsv"
    view_files = tuple(view / name for name in (*FIGURE_CORE_FILES, *FIGURE_RESULT_FILES))
    return [
        Command(
            id="figure_density_adjustment",
            argv=(
                sys.executable,
                "-m",
                "mtp.density_adjustment",
                "--profiles",
                str(human_count_path(context)),
                "--columns",
                str(context.inputs / "taxa_trait_matrix_columns.tsv"),
                "--output-dir",
                str(density),
                "--code-commit",
                context.code_commit,
            ),
            inputs=(
                human_count_path(context),
                context.inputs / "taxa_trait_matrix_columns.tsv",
            ),
            outputs=(
                density_effects,
                density / "density_adjustment_summary.tsv",
                density / "density_adjustment_report.txt",
            ),
        ),
        Command(
            id="figure_inputs",
            argv=(
                sys.executable,
                "-m",
                "mtp.workflow.arms.ilyas",
                "--input-dir",
                str(context.inputs),
                "--enrichment",
                str(enrichment),
                "--depletion",
                str(depletion),
                "--density-effects",
                str(density_effects),
                "--output-dir",
                str(view),
            ),
            inputs=(
                *(context.inputs / name for name in FIGURE_CORE_FILES),
                enrichment,
                depletion,
                density_effects,
            ),
            outputs=view_files,
        ),
        Command(
            id="figures",
            argv=(
                *_script(context, "generate_deck_visualizations.py"),
                "--derived-dir",
                str(view),
                "--output-dir",
                str(figures),
                "--min-studies",
                "5",
            ),
            inputs=view_files,
            outputs=tuple(figures / name for name in FIGURE_FILES),
        ),
    ]


def build_commands(context: ArmContext) -> list[Command]:
    """Preserve ORA thresholds, complete rankings, strata, and calibration settings."""
    profiles = context.inputs / "study_trait_profiles_ora_wide.tsv"
    matrix = context.inputs / "taxa_trait_matrix.tsv"
    experiments = context.inputs / "experiments.tsv"
    strata = context.output / "strata"
    calibration = context.output / "ora_calibration"
    commands = [*_ora_commands(context, "enrichment"), *_ora_commands(context, "depletion")]
    strata_outputs = tuple(
        strata / name
        for stratum in ("pooled", "feces", "feces_excluded")
        for name in (
            f"study_trait_profiles_{stratum}.tsv",
            f"trait_enrichment_results_{stratum}.tsv",
            f"trait_enrichment_report_{stratum}.txt",
            f"trait_enrichment_all_results_{stratum}.tsv",
            f"trait_enrichment_traits_overall_ranked_{stratum}.tsv",
            f"trait_enrichment_trait_study_pairs_ranked_{stratum}.tsv",
        )
    )
    commands.extend(
        [
            Command(
                id="strata",
                argv=(
                    *_script(context, "run_stratified_enrichment.py"),
                    "--profiles",
                    str(profiles),
                    "--experiments",
                    str(experiments),
                    "--taxa-matrix",
                    str(matrix),
                    "--output-dir",
                    str(strata),
                    "--p-threshold",
                    "0.001",
                    "--min-fold",
                    "1.5",
                    "--fdr-scope",
                    "global",
                    "--strata",
                    "pooled",
                    "feces",
                    "feces_excluded",
                ),
                inputs=(profiles, experiments, matrix),
                outputs=(
                    *strata_outputs,
                    strata / "stratification_sensitivity_summary.tsv",
                    strata / "stratification_sensitivity_summary.txt",
                ),
            ),
            Command(
                id="calibration",
                argv=(
                    sys.executable,
                    "-m",
                    "mtp.ora_calibration.harness",
                    "--profiles",
                    str(profiles),
                    "--taxa-matrix",
                    str(matrix),
                    "--output-dir",
                    str(calibration),
                    "--replicates",
                    str(context.replicates if context.replicates is not None else 10),
                    "--seed",
                    str(context.seed if context.seed is not None else 42),
                    "--scenario",
                    "all",
                    "--alpha",
                    "0.05",
                    "--fdr-q",
                    "0.05",
                    "--signal-strength",
                    "0.6",
                    "--n-spike-traits",
                    "3",
                    "--test-type",
                    "enrichment",
                    "--min-fold",
                    "1.5",
                    "--max-fold",
                    "0.67",
                ),
                inputs=(profiles, matrix),
                outputs=tuple(
                    calibration / name
                    for name in (
                        "calibration_summary.json",
                        "calibration_replicates.tsv",
                        "calibration_null_pvalues.tsv",
                        "calibration_ecdf.tsv",
                    )
                ),
            ),
        ]
    )
    if context.figures:
        commands.append(human_count_command(context))
        commands.extend(_figure_commands(context))
    return commands


def materialize_figure_inputs(
    input_dir: Path,
    enrichment: Path,
    depletion: Path,
    density_effects: Path,
    output_dir: Path,
) -> None:
    """Build the renderer's combined view exclusively from this execution's inputs."""
    sources = [*(input_dir / name for name in FIGURE_CORE_FILES)]
    sources.extend((enrichment, depletion, density_effects))
    names = (*FIGURE_CORE_FILES, *FIGURE_RESULT_FILES)
    for source in sources:
        if not source.is_file():
            raise FileNotFoundError(f"Required figure source does not exist: {source}")
    if output_dir.is_symlink() or (
        output_dir.exists()
        and any(path.is_file() or path.is_symlink() for path in output_dir.rglob("*"))
    ):
        raise FileExistsError(f"Figure input view must contain no previous files: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    for source, name in zip(sources, names, strict=True):
        destination = output_dir / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)


def main(argv: list[str] | None = None) -> None:
    """Materialize explicitly named figure sources; never consult legacy defaults."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--enrichment", type=Path, required=True)
    parser.add_argument("--depletion", type=Path, required=True)
    parser.add_argument("--density-effects", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    materialize_figure_inputs(
        args.input_dir, args.enrichment, args.depletion, args.density_effects, args.output_dir
    )


if __name__ == "__main__":
    main()
