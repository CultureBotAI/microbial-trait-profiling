"""Explicit core inputs for Jed's current composition and figure commands (#658)."""

from __future__ import annotations

import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

from mtp.workflow.arms.shared import human_count_command, human_count_path
from mtp.workflow.arms.types import ArmContext, Command


def _module(name: str) -> tuple[str, ...]:
    return sys.executable, "-m", f"mtp.{name}"


def _script(context: ArmContext, name: str) -> tuple[str, ...]:
    return sys.executable, str(context.repository / "scripts" / name)


def _arguments(options: Mapping[str, Path | str | int]) -> tuple[str, ...]:
    return tuple(part for option, value in options.items() for part in (option, str(value)))


def _directory_command(
    context: ArmContext,
    identifier: str,
    executable: tuple[str, ...],
    options: Mapping[str, Path | str | int],
    inputs: Sequence[Path],
    required_outputs: Sequence[str],
) -> Command:
    output = context.output / identifier
    return Command(
        identifier,
        (*executable, *_arguments(options), "--output-dir", str(output)),
        tuple(inputs),
        (output, *(output / filename for filename in required_outputs)),
    )


def _heatmap_command(
    context: ArmContext,
    identifier: str,
    script: str,
    inputs: Mapping[str, Path],
    outputs: Mapping[str, str],
) -> Command:
    output = context.output / identifier
    files = {option: output / filename for option, filename in outputs.items()}
    return Command(
        identifier,
        (*_script(context, script), *_arguments(inputs), *_arguments(files)),
        tuple(inputs.values()),
        (output, *files.values()),
    )


def build_commands(context: ArmContext) -> list[Command]:
    """Preserve statistical defaults while binding each numerical input explicitly."""
    products = context.products
    profiles = human_count_path(context)
    experiments = products["experiments"]
    raw_dump = products["raw_dump"]
    columns = products["feature_metadata"]
    taxa = products["reported_taxa"]
    profile_inputs = (profiles, experiments, raw_dump)
    core_profile_inputs = (products["count_profiles"], experiments, raw_dump)
    profile_options: dict[str, Path | str | int] = {
        "--profile-input": profiles,
        "--experiments-input": experiments,
        "--full-dump-input": raw_dump,
        "--code-commit": context.code_commit,
    }
    matrix_options: dict[str, Path | str | int] = {
        "--profiles": profiles,
        "--columns": columns,
    }
    calibration_options = dict(matrix_options)
    if context.replicates is not None:
        calibration_options["--replicates"] = context.replicates
    if context.seed is not None:
        calibration_options["--seed"] = context.seed

    commands = [
        human_count_command(context),
        _directory_command(
            context,
            "composition",
            _module("trait_composition"),
            {"--core-run": context.manifest, "--code-commit": context.code_commit},
            (context.manifest, *core_profile_inputs),
            ("trait_log_ratio_stats.tsv", "trait_log_ratio_studies.tsv", "core_consumption.json"),
        ),
        _directory_command(
            context,
            "composition_sensitivity",
            _module("log_ratio_sensitivity"),
            {**profile_options, "--taxa-input": taxa},
            (
                *profile_inputs,
                taxa,
                context.repository / "methods" / "table1_log_ratio_sensitivity_protocol.json",
            ),
            ("sensitivity_results.tsv", "leave_one_study_out.tsv", "sensitivity_status.json"),
        ),
        _directory_command(
            context,
            "composition_calibration",
            _module("composition_calibration.composition"),
            calibration_options,
            (profiles, columns),
            ("composition_calibration_summary.json", "composition_calibration_replicates.tsv"),
        ),
        _directory_command(
            context,
            "density_adjustment",
            _script(context, "audit_trait_density_correlation.py"),
            {**matrix_options, "--code-commit": context.code_commit},
            (profiles, columns),
            ("density_adjustment_summary.tsv", "ontology_redundancy_summary.tsv"),
        ),
        _directory_command(
            context,
            "trait_coupling",
            _script(context, "audit_trait_coupling.py"),
            {**matrix_options, "--code-commit": context.code_commit},
            (profiles, columns),
            ("trait_coupling_loto_matrix.tsv", "trait_coupling_lco_summary.tsv"),
        ),
    ]
    if not context.figures:
        return commands

    commands.extend(
        [
            _heatmap_command(
                context,
                "genus_heatmap",
                "create_clustered_heatmap.py",
                {"--input": taxa, "--full-dump-input": raw_dump},
                {"--output": "heatmap.png", "--matrix-output": "matrix.tsv"},
            ),
            _heatmap_command(
                context,
                "condition_taxon_heatmap",
                "create_condition_taxon_heatmap.py",
                {"--taxa-input": taxa, "--full-dump-input": raw_dump},
                {"--output": "heatmap.png", "--matrix-output": "matrix.tsv"},
            ),
            _heatmap_command(
                context,
                "taxa_guild_heatmap",
                "create_taxa_guild_heatmap.py",
                {
                    "--bugsigdb-taxa-input": taxa,
                    "--taxa-trait-matrix-input": products["taxa_matrix"],
                    "--taxa-trait-columns-input": columns,
                    "--full-dump-input": raw_dump,
                    "--study-trait-profiles-input": profiles,
                },
                {
                    "--output": "heatmap.png",
                    "--matrix-output": "matrix.tsv",
                    "--taxa-metadata-output": "taxa_metadata.tsv",
                    "--trait-metadata-output": "trait_metadata.tsv",
                },
            ),
            _heatmap_command(
                context,
                "profile_heatmap",
                "create_study_trait_profile_heatmap.py",
                {
                    "--study-trait-profiles-input": profiles,
                    "--trait-columns-input": columns,
                    "--studies-input": products["study_table"],
                    "--full-dump-input": raw_dump,
                },
                {
                    "--output": "heatmap.png",
                    "--matrix-output": "matrix.tsv",
                    "--study-metadata-output": "profile_metadata.tsv",
                    "--trait-metadata-output": "trait_metadata.tsv",
                    "--paired-output": "paired_heatmap.png",
                    "--paired-matrix-output": "paired_matrix.tsv",
                    "--paired-metadata-output": "paired_metadata.tsv",
                    "--paired-required-trait-check-output": "paired_required_traits.tsv",
                },
            ),
            _heatmap_command(
                context,
                "condition_trait_heatmap",
                "create_condition_trait_module_heatmap.py",
                {
                    "--study-trait-profiles-input": profiles,
                    "--trait-columns-input": columns,
                    "--studies-input": products["study_table"],
                    "--full-dump-input": raw_dump,
                },
                {
                    "--output": "heatmap.png",
                    "--matrix-output": "matrix.tsv",
                    "--group-metadata-output": "group_metadata.tsv",
                    "--trait-metadata-output": "trait_metadata.tsv",
                    "--guild-definition-table-output": "guild_definitions.png",
                    "--guild-definition-table-svg-output": "guild_definitions.svg",
                },
            ),
            _directory_command(
                context,
                "guild_maps",
                _script(context, "create_guild_maps.py"),
                {"--core-run": context.manifest, "--code-commit": context.code_commit},
                (context.manifest, *core_profile_inputs),
                (
                    "guild_condition_body_site_heatmap.png",
                    "guild_condition_body_site_matrix.tsv",
                    "core_consumption.json",
                ),
            ),
            _directory_command(
                context,
                "table1_traits",
                _script(context, "create_table1_trait_maps.py"),
                {"--core-run": context.manifest, "--code-commit": context.code_commit},
                (context.manifest, *core_profile_inputs),
                ("table1_traits_direction_shift.tsv", "core_consumption.json"),
            ),
        ]
    )
    return commands
