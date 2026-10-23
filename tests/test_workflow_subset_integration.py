"""Raw-to-core subset and downstream provenance integration checks for #658."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest
import yaml

from mtp.datamodel.workflow_pydantic import WorkflowRun
from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow.__main__ import build_parser, main
from mtp.workflow.arms.noel import build_commands
from mtp.workflow.arms.runtime import execute_arm, verify_arm
from mtp.workflow.arms.types import ArmContext, Command
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import artifact_path, content_hash, execute_workflow, validate_run
from mtp.workflow.verification import verify_core_run

POPULATION = "disease_vs_healthy_tier1_2_subject_controls"
FIXTURE = REPO_ROOT / "tests/fixtures/workflow"
INPUTS = {
    "raw_dump": FIXTURE / "contrast_population.csv",
    "kg_edges": FIXTURE / "merged-kg_edges.tsv",
    "kg_nodes": FIXTURE / "merged-kg_nodes.tsv",
}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


@pytest.fixture(scope="module")
def cores(tmp_path_factory: pytest.TempPathFactory) -> dict[str, WorkflowRun]:
    output = tmp_path_factory.mktemp("subset-core")
    return {
        name: execute_workflow(
            REPO_ROOT / "workflows" / spec,
            INPUTS,
            output / name,
            OPERATIONS,
            parameter_overrides={
                "source_mode": "fixture",
                **({"contrast_population": population} if name == "subset" else {}),
            },
            validator=validate_core_product,
            repository_root=REPO_ROOT,
        )
        for name, spec, population in (
            ("default", "bugsigdb_core.yaml", "all_paired"),
            ("subset", "bugsigdb_core_subset.yaml", POPULATION),
        )
    }


def test_subset_retains_source_background_and_exact_oriented_profile_values(
    cores: dict[str, WorkflowRun],
) -> None:
    original, selected = cores["default"], cores["subset"]
    assert verify_core_run(Path(selected.manifest_path)) == selected
    assert original.contrast_population == "all_paired"
    assert original.direction_basis == "group_1_vs_group_0"
    assert selected.contrast_population == POPULATION
    assert selected.direction_basis == "disease_vs_comparison"
    for product in ("taxa_matrix", "feature_metadata", "feature_catalog"):
        assert file_sha256(artifact_path(original, product)) == file_sha256(
            artifact_path(selected, product)
        )
    experiments = _rows(artifact_path(selected, "experiments"))
    assert {row["experiment_id"] for row in experiments} == {"bsdb:1-1", "bsdb:1-2"}
    reversed_metadata = next(row for row in experiments if row["experiment_id"] == "bsdb:1-2")
    assert reversed_metadata["group_1_name"] == "Healthy controls"
    for product in ("count_profiles", "ora_profiles"):
        native = {
            (row["study_id"], row["experiment_id"], row["inc_dec"]): row
            for row in _rows(artifact_path(original, product))
        }
        projected = _rows(artifact_path(selected, product))
        assert len(projected) == 4
        for row in projected:
            direction = row["inc_dec"]
            if row["experiment_id"] == "bsdb:1-2":
                direction = {"inc": "dec", "dec": "inc"}[direction]
            before = native[(row["study_id"], row["experiment_id"], direction)]
            assert row == {**before, "inc_dec": row["inc_dec"]}
    audit = json.loads(artifact_path(selected, "contrast_selection").read_text())
    assert audit["reversed_experiments"] == [{"study_id": "bsdb:1", "experiment_id": "bsdb:1-2"}]
    assert audit["counts"]["source"]["experiments"] > audit["counts"]["selected"]["experiments"]
    assert file_sha256(artifact_path(original, "reported_taxa")) == file_sha256(
        artifact_path(selected, "source_reported_taxa")
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("contrast_population", "all_paired"),
        ("contrast_population", "case_vs_control_subject"),
        ("contrast_population", None),
        ("direction_basis", "group_1_vs_group_0"),
        ("direction_basis", None),
    ],
)
def test_core_record_cannot_relabel_its_executed_population(
    cores: dict[str, WorkflowRun], field: str, value: str | None
) -> None:
    path = Path(cores["subset"].manifest_path)
    original = path.read_bytes()
    changed = json.loads(original)
    changed[field] = value
    try:
        path.write_text(json.dumps(changed))
        with pytest.raises(ValueError, match="contrast population or direction basis"):
            validate_run(path)
    finally:
        path.write_bytes(original)


def test_nondefault_parameter_cannot_label_an_unfiltered_execution(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Undeclared parameter overrides"):
        execute_workflow(
            REPO_ROOT / "workflows/bugsigdb_core.yaml",
            INPUTS,
            tmp_path / "runs",
            OPERATIONS,
            parameter_overrides={"source_mode": "fixture", "contrast_population": POPULATION},
            validator=validate_core_product,
        )
    assert not (tmp_path / "runs").exists()


def test_declared_population_requires_an_actual_selection_stage(tmp_path: Path) -> None:
    definition = yaml.safe_load((REPO_ROOT / "workflows/bugsigdb_core.yaml").read_text())
    definition["parameters"].append({"name": "contrast_population", "default": POPULATION})
    definition["steps"][0]["parameters"].append("contrast_population")
    spec = tmp_path / "mislabelled.yaml"
    spec.write_text(yaml.safe_dump(definition))
    with pytest.raises(ValueError, match="requires the subset selection contract"):
        execute_workflow(spec, INPUTS, tmp_path / "runs", OPERATIONS)
    assert not (tmp_path / "runs").exists()


@pytest.mark.parametrize(
    "mutation",
    [
        "selected_background",
        "matrix_reads_selected_taxa",
        "exports_read_selected_studies",
        "unverified_source_studies",
        "unverified_source_taxa",
        "missing_source_identity",
        "wrong_ingest_operation",
        "incomplete_feature_lineage",
    ],
)
def test_subset_requires_full_source_background_before_execution(
    tmp_path: Path, mutation: str
) -> None:
    definition = yaml.safe_load((REPO_ROOT / "workflows/bugsigdb_core_subset.yaml").read_text())
    steps = {step["id"]: step for step in definition["steps"]}
    products = {product["id"]: product for product in definition["products"]}
    if mutation in {"selected_background", "matrix_reads_selected_taxa"}:
        steps["taxa_matrix"].update(
            operation="taxa_matrix", inputs=["reported_taxa", "kg_edges", "kg_nodes"]
        )
        if mutation == "selected_background":
            removed = set(steps["source_export_tables"]["outputs"])
            definition["steps"].remove(steps["source_export_tables"])
            definition["products"] = [
                product for product in definition["products"] if product["id"] not in removed
            ]
    elif mutation == "exports_read_selected_studies":
        steps["source_export_tables"]["inputs"] = ["studies"]
    elif mutation == "unverified_source_studies":
        products["source_studies"]["kind"] = "structural"
    elif mutation == "unverified_source_taxa":
        products["source_reported_taxa"]["kind"] = "structural"
    elif mutation == "missing_source_identity":
        definition["steps"].remove(steps["identify_sources"])
        definition["products"].remove(products["source_identity"])
        definition["parameters"] = [
            parameter
            for parameter in definition["parameters"]
            if parameter["name"] != "source_mode"
        ]
        steps["ingest"].pop("depends_on")
    elif mutation == "wrong_ingest_operation":
        steps["ingest"]["operation"] = "ingest"
    elif mutation == "incomplete_feature_lineage":
        steps["feature_catalog"]["inputs"] = ["feature_metadata"]
    spec = tmp_path / "invalid-background.yaml"
    spec.write_text(yaml.safe_dump(definition))

    with pytest.raises(ValueError, match="requires the full-source background contract"):
        execute_workflow(
            spec,
            INPUTS,
            tmp_path / "runs",
            OPERATIONS,
            parameter_overrides={"contrast_population": POPULATION},
            validator=validate_core_product,
        )
    assert not (tmp_path / "runs").exists()


def test_serialized_subset_cannot_redeclare_a_selected_taxa_background(
    cores: dict[str, WorkflowRun],
) -> None:
    path = Path(cores["subset"].manifest_path)
    original = path.read_bytes()
    changed = json.loads(original)
    matrix = next(
        step for step in changed["definition"]["steps"] if step["operation"] == "source_taxa_matrix"
    )
    matrix.update(operation="taxa_matrix", inputs=["reported_taxa", "kg_edges", "kg_nodes"])
    changed["definition_sha256"] = content_hash(changed["definition"])
    try:
        path.write_text(json.dumps(changed))
        with pytest.raises(ValueError, match="requires the full-source background contract"):
            validate_run(path)
    finally:
        path.write_bytes(original)


def test_subset_fixture_cli_executes_the_complete_selection_graph(tmp_path: Path) -> None:
    assert (
        main(["fixture", "--output-root", str(tmp_path), "--contrast-population", POPULATION]) == 0
    )
    manifests = list(tmp_path.glob("run-*/run.json"))
    assert len(manifests) == 1
    run = verify_core_run(manifests[0])
    assert run.contrast_population == POPULATION
    assert len(_rows(artifact_path(run, "experiments"))) == 2


@pytest.mark.parametrize("source", [[], ["fixture"]])
def test_rerun_core_accepts_population_without_optional_source_position(source: list[str]) -> None:
    args = build_parser().parse_args(
        [
            "prepare-core",
            "/out",
            "/raw",
            "/edges",
            "/nodes",
            *source,
            "--contrast-population",
            POPULATION,
        ]
    )
    assert args.contrast_population == POPULATION
    assert args.source_mode == ("fixture" if source else "pinned")


@pytest.mark.parametrize("arm", ["jed", "noel", "ilyas"])
def test_each_arm_inherits_and_verifies_subset_identity(
    cores: dict[str, WorkflowRun], tmp_path: Path, arm: str
) -> None:
    core = cores["subset"]
    product = "ora_profiles" if arm == "ilyas" else "count_profiles"

    def plan(context: ArmContext) -> list[Command]:
        assert context.contrast_population == POPULATION
        assert context.direction_basis == "disease_vs_comparison"
        source, output = context.products[product], context.output / "observed.tsv"
        return [
            Command(
                "observe_subset",
                (
                    sys.executable,
                    str(REPO_ROOT / "scripts/adhoc/issue658_arm_test_command.py"),
                    "--input",
                    str(source),
                    "--output",
                    str(output),
                ),
                (source,),
                (output,),
            )
        ]

    result = execute_arm(arm, Path(core.manifest_path), tmp_path / arm, plan, figures=False)
    assert result.contrast_population == POPULATION
    assert result.direction_basis == "disease_vs_comparison"
    assert verify_arm(Path(result.manifest_path)) == result
    changed = json.loads(Path(result.manifest_path).read_text())
    changed["contrast_population"] = "all_paired"
    Path(result.manifest_path).write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="Arm contrast population"):
        verify_arm(Path(result.manifest_path))


def test_noel_plan_consumes_selection_evidence_for_already_oriented_core(
    cores: dict[str, WorkflowRun], tmp_path: Path
) -> None:
    core = cores["subset"]
    context = ArmContext(
        repository=REPO_ROOT,
        manifest=Path(core.manifest_path),
        output=tmp_path,
        inputs=tmp_path / "inputs",
        products={
            item.product_id: artifact_path(core, item.product_id) for item in core.artifacts or []
        },
        code_commit="fixture",
        replicates=1,
        seed=182,
        figures=False,
        contrast_population=POPULATION,
    )
    command = next(item for item in build_commands(context) if item.id == "contrast_orientation")
    assert "mtp.workflow.orientation_assessment" in command.argv
    assert command.argv[command.argv.index("--population") + 1] == POPULATION
    for name in ("contrast_selection", "studies", "source_studies"):
        assert context.products[name] in command.inputs
