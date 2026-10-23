"""Thread the BacDive strain rollup scope through the core workflow and its contracts (#788)."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path
from typing import cast

import pytest

from mtp.study_trait_profiles import collect_profile_taxa
from mtp.taxa_trait_matrix import (
    PRIMARY_STRAIN_ROLLUP,
    STRAIN_ROLLUPS,
    StrainRollup,
    create_taxa_trait_matrix,
)
from mtp.workflow.__main__ import main as workflow_main
from mtp.workflow.comparison import core_analysis_bases
from mtp.workflow.coverage import STRAIN_ATTRIBUTED_EVIDENCE, UNRESOLVED_EVIDENCE
from mtp.workflow.definition import read_definition
from mtp.workflow.models import WorkflowRun
from mtp.workflow.operations import (
    OPERATIONS,
    Operation,
    feature_catalog_operation,
    trace_profile_feature,
    validate_core_product,
)
from mtp.workflow.runtime import WorkflowExecutionError, artifact_path, execute_workflow
from mtp.workflow.verification import verify_core_run

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "workflows" / "bugsigdb_core.yaml"
FIXTURE = ROOT / "tests" / "fixtures" / "workflow"
UNANNOTATED = "NCBITaxon:300"
# Two BacDive strains of the fixture's unannotated taxon. Strain 1 asserts a phenotype,
# a clinical class the phenotype scope excludes and glucose fermentation; strain 2
# asserts that glucose is not fermented, which the full scope drops.
STRAIN_EDGES = [
    ("strain1", "kgmicrobe.strain:1", "biolink:subclass_of", UNANNOTATED),
    ("strain1_shape", "kgmicrobe.strain:1", "biolink:has_phenotype", "METPO:1000678"),
    ("strain1_clinical", "kgmicrobe.strain:1", "biolink:has_phenotype", "METPO:1004004"),
    ("strain1_ferments", "kgmicrobe.strain:1", "METPO:2000011", "CHEBI:17234"),
    ("strain2", "kgmicrobe.strain:2", "biolink:subclass_of", UNANNOTATED),
    ("strain2_absent", "kgmicrobe.strain:2", "METPO:2000037", "CHEBI:17234"),
]
STRAIN_NODES = [
    ("METPO:1000678", "biolink:PhenotypicQuality", "oval shaped"),
    ("METPO:1004004", "biolink:PhenotypicQuality", "human pathogen"),
    ("METPO:2000011", "biolink:OntologyClass", "ferments"),
    ("METPO:2000037", "biolink:OntologyClass", "does not ferment"),
    ("CHEBI:17234", "biolink:SmallMolecule", "glucose"),
]
SHAPE = "METPO:1000678"
TRAIT = "METPO:1000608"
MEDIUM = "METPO:2000510|mediadive.medium:1"
ADDED = {
    "off": set(),
    "phenotype": {SHAPE},
    "full": {SHAPE, "METPO:1004004", "METPO:2000011|CHEBI:17234"},
}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _header(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


@pytest.fixture(scope="module")
def strain_inputs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """The workflow fixture with two strain records appended to its KG."""
    root = tmp_path_factory.mktemp("strain-kg")
    edges, nodes = root / "merged-kg_edges.tsv", root / "merged-kg_nodes.tsv"
    shutil.copyfile(FIXTURE / "merged-kg_edges.tsv", edges)
    shutil.copyfile(FIXTURE / "merged-kg_nodes.tsv", nodes)
    with edges.open("a", encoding="utf-8") as handle:
        handle.writelines("\t".join([*edge, "fixture:strain"]) + "\n" for edge in STRAIN_EDGES)
    with nodes.open("a", encoding="utf-8") as handle:
        handle.writelines("\t".join(node) + "\n" for node in STRAIN_NODES)
    return {"raw_dump": FIXTURE / "full_dump.csv", "kg_edges": edges, "kg_nodes": nodes}


def _run(
    inputs: dict[str, Path], root: Path, scope: StrainRollup, operations: dict[str, Operation]
) -> WorkflowRun:
    return execute_workflow(
        SPEC,
        inputs,
        root,
        operations,
        parameter_overrides={"source_mode": "fixture", "strain_rollup": scope},
        validator=validate_core_product,
        repository_root=ROOT,
    )


@pytest.fixture(scope="module")
def scoped_runs(
    strain_inputs: dict[str, Path], tmp_path_factory: pytest.TempPathFactory
) -> dict[StrainRollup, WorkflowRun]:
    root = tmp_path_factory.mktemp("strain-runs")
    return {
        scope: _run(strain_inputs, root / scope, scope, dict(OPERATIONS))
        for scope in STRAIN_ROLLUPS
    }


def test_both_core_specs_default_to_the_phenotype_scope() -> None:
    for spec in (SPEC, ROOT / "workflows" / "bugsigdb_core_subset.yaml"):
        definition = read_definition(spec)
        defaults = {item.name: item.default for item in definition.parameters or []}
        assert defaults["strain_rollup"] == PRIMARY_STRAIN_ROLLUP == "phenotype"
        steps = {step.id: step for step in definition.steps}
        for step_id in ("taxa_matrix", "feature_catalog"):
            assert steps[step_id].parameters == ["strain_rollup"]


@pytest.mark.parametrize("scope", ["off", "phenotype", "full"])
def test_each_scope_publishes_a_matrix_its_contracts_accept(
    scoped_runs: dict[StrainRollup, WorkflowRun], scope: StrainRollup
) -> None:
    """Every product passed its contract at publication and again on re-verification,
    including the KG source reconciliation of each scope's matrix."""
    run = scoped_runs[scope]
    recorded = {item.name: (item.value, item.source) for item in run.parameters or []}
    assert recorded["strain_rollup"] == (scope, "override")
    steps = {step.step_id: step for step in run.steps or []}
    for step_id in ("taxa_matrix", "feature_catalog"):
        assert [(p.name, p.value) for p in steps[step_id].parameters or []] == [
            ("strain_rollup", scope)
        ]
    verify_core_run(Path(run.manifest_path))

    baseline = set(_header(artifact_path(scoped_runs["off"], "taxa_matrix"))[3:])
    columns = set(_header(artifact_path(run, "taxa_matrix"))[3:])
    assert columns == baseline | ADDED[scope]
    assert "METPO:2000037|CHEBI:17234" not in columns
    profiles = set(_header(artifact_path(run, "count_profiles"))[3:])
    assert profiles == columns


@pytest.mark.parametrize(
    ("scope", "attributed"),
    [("off", set()), ("phenotype", {SHAPE, TRAIT}), ("full", set(ADDED["full"]) | {TRAIT, MEDIUM})],
)
def test_the_catalog_marks_features_whose_counts_may_include_strain_records(
    scoped_runs: dict[StrainRollup, WorkflowRun], scope: StrainRollup, attributed: set[str]
) -> None:
    catalog = {
        row["column_id"]: row for row in _rows(artifact_path(scoped_runs[scope], "feature_catalog"))
    }
    assert {
        feature
        for feature, row in catalog.items()
        if row["annotation_evidence"] == STRAIN_ATTRIBUTED_EVIDENCE
    } == attributed
    assert all(
        row["annotation_evidence"] == UNRESOLVED_EVIDENCE
        for feature, row in catalog.items()
        if feature not in attributed
    )
    assert set(catalog) == set(_header(artifact_path(scoped_runs[scope], "taxa_matrix"))[3:])


def test_a_matrix_built_under_another_scope_is_rejected(
    strain_inputs: dict[str, Path], tmp_path: Path
) -> None:
    """The source contract applies the declared scope, so an off matrix published by a
    phenotype run fails before publication."""

    def off_matrix(
        inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
    ) -> None:
        assert parameters == {"strain_rollup": "phenotype"}
        create_taxa_trait_matrix(
            taxa_path=inputs["reported_taxa"],
            edges_path=inputs["kg_edges"],
            nodes_path=inputs["kg_nodes"],
            output_path=outputs["taxa_matrix"],
            column_metadata_path=outputs["feature_metadata"],
            strain_rollup="off",
        )

    with pytest.raises(WorkflowExecutionError, match="KG"):
        _run(
            strain_inputs, tmp_path / "runs", "phenotype", {**OPERATIONS, "taxa_matrix": off_matrix}
        )


def test_a_catalog_written_under_another_scope_is_rejected(
    strain_inputs: dict[str, Path], tmp_path: Path
) -> None:
    def off_catalog(
        inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
    ) -> None:
        assert parameters == {"strain_rollup": "phenotype"}
        feature_catalog_operation(inputs, outputs, {"strain_rollup": "off"})

    with pytest.raises(WorkflowExecutionError, match="evidence"):
        _run(
            strain_inputs,
            tmp_path / "runs",
            "phenotype",
            {**OPERATIONS, "feature_catalog": off_catalog},
        )


def _strain_profile(run: WorkflowRun) -> tuple[str, str, str]:
    profiles, _ = collect_profile_taxa(artifact_path(run, "studies"))
    return next(key for key, taxa in sorted(profiles.items()) if UNANNOTATED in taxa)


def test_the_trace_reconciles_strain_attributed_counts(
    scoped_runs: dict[StrainRollup, WorkflowRun], strain_inputs: dict[str, Path]
) -> None:
    run = scoped_runs["phenotype"]
    key = _strain_profile(run)
    arguments = (
        artifact_path(run, "studies"),
        artifact_path(run, "taxa_matrix"),
        strain_inputs["kg_edges"],
        strain_inputs["kg_nodes"],
        key,
        SHAPE,
    )

    result = trace_profile_feature(*arguments, strain_rollup="phenotype")

    assert result["count_profile_value"] == 1
    ledger = cast(list[dict[str, object]], result["ledger"])
    assert [row["source_taxon_id"] for row in ledger] == [UNANNOTATED]
    edges = cast(list[dict[str, str]], ledger[0]["evidence_edges"])
    assert [edge["id"] for edge in edges] == ["strain1_shape"]
    with pytest.raises(ValueError, match="does not reconcile"):
        trace_profile_feature(*arguments)


def test_cell_lineage_traces_under_the_recorded_scope(
    scoped_runs: dict[StrainRollup, WorkflowRun], capsys: pytest.CaptureFixture[str]
) -> None:
    run = scoped_runs["phenotype"]
    key = _strain_profile(run)

    workflow_main(
        ["lineage", run.manifest_path, "count_profiles", "--profile", *key, "--feature", SHAPE]
    )

    assert json.loads(capsys.readouterr().out)["count_profile_value"] == 1


def test_the_fixture_and_run_commands_take_a_scope(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The sensitivity and exploratory cores are one option away from the primary (#788)."""
    workflow_main(["fixture", "--output-root", str(tmp_path / "runs"), "--strain-rollup", "off"])
    run = verify_core_run(Path(capsys.readouterr().out.strip()))
    assert {item.name: item.value for item in run.parameters or []}["strain_rollup"] == "off"

    with pytest.raises(SystemExit):
        workflow_main(
            [
                "run",
                "--output-root",
                str(tmp_path / "other"),
                "--parameter",
                "strain_rollup=full",
                "--strain-rollup",
                "off",
            ]
        )
    assert "--strain-rollup conflicts" in capsys.readouterr().err


@pytest.mark.parametrize(("scope", "retained"), [("off", 1), ("phenotype", 2), ("full", 4)])
def test_comparison_bases_record_the_core_scope_and_trait_count(
    scoped_runs: dict[StrainRollup, WorkflowRun], scope: StrainRollup, retained: int
) -> None:
    """A ComparisonRun names the scope its core declared and that core's own m (#788)."""
    bases = core_analysis_bases(scoped_runs[scope])

    assert [basis.profile_unit for basis in bases] == ["count_profiles", "ora_profiles"]
    assert {basis.rollup_policy for basis in bases} == {f"ncbi_taxonomy_closure_rollup_{scope}"}
    assert {basis.retained_feature_count for basis in bases} == {retained}
