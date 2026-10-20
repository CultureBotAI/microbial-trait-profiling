"""Explicit, verified-input command plan for Noel's main-branch analyses (#658)."""

from __future__ import annotations

import sys
from pathlib import Path

from mtp.workflow.arms.types import ArmContext, Command


def build_commands(context: ArmContext) -> list[Command]:
    """Keep statistical protocols while isolating every data and status path.

    The influence protocol deliberately still requires its declared KG bytes.
    A different verified KG cannot silently inherit that protocol's eligibility.
    Historical branch scripts and anchored condition follow-ups are separate work.
    """
    output = context.output
    products = context.products
    profile = products["count_profiles"]
    experiments = products["experiments"]
    raw = products["raw_dump"]
    methods = context.repository / "methods"
    data_inputs = (profile, experiments, raw)
    data_args = (
        "--profile-input",
        str(profile),
        "--experiments-input",
        str(experiments),
        "--full-dump-input",
        str(raw),
    )
    commit_args = ("--code-commit", context.code_commit)
    seed_args = () if context.seed is None else ("--seed", str(context.seed))
    replicate_args = () if context.replicates is None else ("--replicates", str(context.replicates))

    def command(
        name: str,
        module: str,
        args: tuple[str, ...],
        inputs: tuple[Path, ...],
        outputs: tuple[Path, ...],
    ) -> Command:
        return Command(name, (sys.executable, "-m", module, *args), inputs, outputs)

    screen = output / "trait_direction_screen"
    flips = output / "shared_flips"
    screen_calibration = output / "screen_calibration"
    calibration = output / "calibration"
    orientation = output / "contrast_orientation"
    influence = output / "table1_taxon_influence"
    validation = output / "validation_evidence"
    integrated = output / "integrated_evidence"
    protocol = methods / "calibration_protocol.json"
    flip_table = flips / "calibration_flip_assignments.tsv"
    source_inputs = tuple(
        context.inputs / filename
        for filename in (
            "study_trait_profiles_wide.tsv",
            "experiments.tsv",
            "bugsigdb_taxa.tsv",
            "taxa_trait_matrix.tsv",
        )
    )
    contribution_protocol = methods / "table1_taxon_influence_protocol.json"
    native_protocol = methods / "table1_native_influence_protocol.json"
    native_commit = context.code_commit.removesuffix("+dirty")
    native_commit_args: tuple[str, ...] = ("--code-commit", native_commit)
    if native_commit != context.code_commit:
        native_commit_args += ("--code-dirty",)
    orientation_module = "mtp.contrast_orientation"
    orientation_extra_args: tuple[str, ...] = ()
    orientation_extra_inputs: tuple[Path, ...] = ()
    if context.contrast_population != "all_paired":
        orientation_module = "mtp.workflow.orientation_assessment"
        orientation_extra_inputs = (
            context.manifest,
            products["contrast_selection"],
            products["studies"],
            products["source_studies"],
        )
        orientation_extra_args = (
            "--core-run",
            str(context.manifest),
            "--population",
            context.contrast_population,
            "--selection-input",
            str(products["contrast_selection"]),
            "--studies-input",
            str(products["studies"]),
            "--source-studies-input",
            str(products["source_studies"]),
        )

    commands = [
        command(
            "directional_screen",
            "mtp.trait_screen",
            (
                "--core-run",
                str(context.manifest),
                *data_args,
                "--column-labels",
                str(products["feature_metadata"]),
                "--output-dir",
                str(screen),
                *commit_args,
            ),
            (context.manifest, *data_inputs, products["feature_metadata"]),
            (
                screen,
                screen / "trait_direction_screen.tsv",
                screen / "trait_screen_provenance.json",
                screen / "core_consumption.json",
            ),
        ),
        command(
            "shared_flip_assignments",
            "mtp.calibration",
            (
                *data_args,
                "--output-dir",
                str(flips),
                "--flip-assignments-only",
                *seed_args,
                *commit_args,
            ),
            (*data_inputs, protocol),
            (flips, flip_table, flips / "calibration_flip_assignments.json"),
        ),
        command(
            "screen_calibration",
            "mtp.screen_calibration",
            (
                *data_args,
                "--output-dir",
                str(screen_calibration),
                "--shared-flips",
                str(flip_table),
                *replicate_args,
                *seed_args,
                *commit_args,
            ),
            (*data_inputs, protocol, flip_table),
            (
                screen_calibration,
                screen_calibration / "screen_null_replicates.tsv",
                screen_calibration / "screen_null_summary.tsv",
                screen_calibration / "screen_null_provenance.json",
            ),
        ),
        command(
            "empirical_calibration",
            "mtp.calibration",
            (
                *data_args,
                "--output-dir",
                str(calibration),
                *replicate_args,
                *seed_args,
                *commit_args,
            ),
            (*data_inputs, protocol),
            (
                calibration,
                calibration / "calibration_results.tsv",
                calibration / "calibration_criteria.tsv",
                calibration / "calibration_status.json",
            ),
        ),
        command(
            "contrast_orientation",
            orientation_module,
            (
                *data_args,
                *orientation_extra_args,
                "--output-dir",
                str(orientation),
                "--review-input",
                str(methods / "contrast_orientation_review.tsv"),
                *commit_args,
            ),
            (*data_inputs, *orientation_extra_inputs, methods / "contrast_orientation_review.tsv"),
            (
                orientation,
                orientation / "contrast_manifest.tsv",
                orientation / "contrast_orientation_provenance.json",
            ),
        ),
        command(
            "taxon_contributions",
            "mtp.taxon_contributions",
            (
                "--kg-dir",
                str(context.inputs),
                "--inputs",
                str(context.inputs),
                "--output",
                str(influence),
                "--protocol",
                str(contribution_protocol),
            ),
            (
                *source_inputs,
                context.inputs / "merged-kg_edges.tsv",
                context.inputs / "merged-kg_nodes.tsv",
                contribution_protocol,
            ),
            (influence, influence / "manifest.json"),
        ),
        command(
            "native_taxon_influence",
            "mtp.native_taxon_influence",
            (
                "--inputs",
                str(context.inputs),
                "--output",
                str(influence),
                "--protocol",
                str(native_protocol),
                *native_commit_args,
            ),
            (*source_inputs, influence, native_protocol),
            (
                influence,
                influence / "native_manifest.json",
                influence / "native_composition_comparison.tsv",
                influence / "native_ora_fold_comparison.tsv",
                influence / "native_summary.md",
            ),
        ),
        command(
            "native_influence_review",
            "mtp.native_influence_publication",
            ("--source", str(influence), "--destination", str(influence / "review")),
            (influence,),
            (influence / "review", influence / "review" / "publication_manifest.json"),
        ),
        command(
            "validation_evidence",
            "mtp.validation_evidence",
            (
                "--core-run",
                str(context.manifest),
                *data_args,
                "--nodes-input",
                str(products["kg_nodes"]),
                "--external-cohorts",
                str(methods / "external_cohorts_sun2024.tsv"),
                "--output-dir",
                str(validation),
                "--require-nodes",
            ),
            (
                context.manifest,
                *data_inputs,
                products["kg_nodes"],
                methods / "external_cohorts_sun2024.tsv",
            ),
            (
                validation,
                validation / "validation_evidence_provenance.json",
                validation / "inferred_function_traits.tsv",
                validation / "core_consumption.json",
            ),
        ),
        command(
            "integrated_evidence",
            "mtp.integrated_evidence",
            (
                *data_args,
                "--output-dir",
                str(integrated),
                "--calibration-status",
                str(calibration / "calibration_status.json"),
                "--orientation-provenance",
                str(orientation / "contrast_orientation_provenance.json"),
                "--sensitivity-status",
                str(output / "pending_assessments" / "sensitivity_status.json"),
                "--coverage-status",
                str(output / "pending_assessments" / "coverage_status.json"),
                "--influence-status",
                str(influence / "native_manifest.json"),
                *commit_args,
            ),
            (
                *data_inputs,
                calibration / "calibration_status.json",
                orientation / "contrast_orientation_provenance.json",
                influence / "native_manifest.json",
            ),
            (
                integrated,
                integrated / "integrated_evidence.tsv",
                integrated / "integrated_evidence_provenance.json",
                integrated / "integrated_evidence_assessments.tsv",
            ),
        ),
    ]
    return commands
