"""Thread the BacDive strain rollup scope through the core workflow and its contracts (#788)."""

from __future__ import annotations

import csv
import json
import shutil
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import cast

import pytest

from mtp.paths import file_sha256
from mtp.study_trait_profiles import collect_profile_taxa
from mtp.taxa_trait_matrix import (
    PRIMARY_STRAIN_ROLLUP,
    STRAIN_ROLLUPS,
    StrainRollup,
    create_taxa_trait_matrix,
)
from mtp.workflow import runtime
from mtp.workflow.__main__ import main as workflow_main
from mtp.workflow.arms import runtime as arms_runtime
from mtp.workflow.arms.runtime import execute_arm, verify_arm
from mtp.workflow.arms.types import ArmContext, Command
from mtp.workflow.comparison import core_analysis_bases
from mtp.workflow.coverage import (
    STRAIN_ATTRIBUTED_EVIDENCE,
    UNRESOLVED_EVIDENCE,
    feature_class,
    feature_filter_id,
)
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
            assert steps[step_id].parameters == ["strain_rollup", "trait_key_policy"]


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
            ("strain_rollup", scope),
            ("trait_key_policy", "facultative_union_v1"),
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


@pytest.mark.parametrize(
    ("scope", "evidence", "attributed"),
    [
        ("off", UNRESOLVED_EVIDENCE, 0),
        ("phenotype", STRAIN_ATTRIBUTED_EVIDENCE, 2),
        ("full", STRAIN_ATTRIBUTED_EVIDENCE, 5),
    ],
)
def test_the_recorded_contract_facts_name_the_evidence_the_scope_admits(
    scoped_runs: dict[StrainRollup, WorkflowRun],
    scope: StrainRollup,
    evidence: str,
    attributed: int,
) -> None:
    """A phenotype or full run's validation facts must not describe its counts as only
    direct or descendant assertions, as an off run's do."""
    artifacts = {item.product_id: item for item in scoped_runs[scope].artifacts or []}
    for product_id in ("feature_metadata", "feature_catalog", "taxa_matrix"):
        (facts,) = (
            json.loads(item.details_json or "{}")
            for item in artifacts[product_id].validations or []
            if item.validator.startswith("core:")
        )
        assert (facts["strain_rollup"], facts["annotation_evidence"]) == (scope, evidence)
        if product_id != "taxa_matrix":
            assert facts["strain_attributed_features"] == attributed


def test_a_matrix_built_under_another_scope_is_rejected(
    strain_inputs: dict[str, Path], tmp_path: Path
) -> None:
    """The source contract applies the declared scope, so an off matrix published by a
    phenotype run fails before publication."""

    def off_matrix(
        inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
    ) -> None:
        assert parameters == {
            "strain_rollup": "phenotype",
            "trait_key_policy": "facultative_union_v1",
        }
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
        assert parameters == {
            "strain_rollup": "phenotype",
            "trait_key_policy": "facultative_union_v1",
        }
        feature_catalog_operation(inputs, outputs, {**parameters, "strain_rollup": "off"})

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


NEW_TRAIT = "METPO:1000999"
NEW_MEDIUM = "METPO:2000510|mediadive.medium:2"


def _write(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows([header, *rows])


def _with_features(
    run: WorkflowRun,
    root: Path,
    order: list[str],
    *,
    profile_order: list[str] | None = None,
    relabel: bool = False,
) -> WorkflowRun:
    """A copy of a core run whose feature catalog lists ``order`` and whose count and ORA
    profiles hold ``profile_order`` (by default the same); an added feature counts zero."""
    root.mkdir()
    source = artifact_path(run, "feature_catalog")
    header = _header(source)
    catalog = {row["column_id"]: row for row in _rows(source)}
    template = next(iter(catalog.values()))
    rows = [
        catalog.get(
            feature, {**template, "column_id": feature, "feature_class": feature_class(feature)}
        )
        for feature in order
    ]
    if relabel:
        rows = [{**row, "column_label": f"{row['column_label']} (renamed)"} for row in rows]
    paths = {"feature_catalog": root / source.name}
    _write(paths["feature_catalog"], header, [[row[column] for column in header] for row in rows])
    for product, leading in (("count_profiles", 3), ("ora_profiles", 4)):
        source = artifact_path(run, product)
        columns = [*_header(source)[:leading], *(order if profile_order is None else profile_order)]
        paths[product] = root / source.name
        _write(
            paths[product],
            columns,
            [[row.get(column, "0") for column in columns] for row in _rows(source)],
        )
    artifacts = [
        item.model_copy(update={"path": str(paths[item.product_id])})
        if item.product_id in paths
        else item
        for item in run.artifacts or []
    ]
    return run.model_copy(update={"artifacts": artifacts})


@pytest.mark.parametrize("scope", STRAIN_ROLLUPS)
def test_each_basis_identifies_its_retained_features_in_profile_order(
    scoped_runs: dict[StrainRollup, WorkflowRun], tmp_path: Path, scope: StrainRollup
) -> None:
    """feature_filter_id identifies the retained features a basis analyses, in its profile
    table's column order; the catalog's bytes are only its feature_catalog_sha256 (#860)."""
    run = scoped_runs[scope]
    features = _header(artifact_path(run, "count_profiles"))[3:]
    retained = [feature for feature in features if feature_class(feature) == "trait_annotation"]

    def identities(name: str, order: list[str], *, relabel: bool = False) -> list[str]:
        copy = _with_features(run, tmp_path / name, order, relabel=relabel)
        bases = core_analysis_bases(copy)
        assert {basis.feature_catalog_sha256 for basis in bases} == {
            file_sha256(artifact_path(copy, "feature_catalog"))
        }
        return [basis.feature_filter_id for basis in bases]

    bases = core_analysis_bases(run)
    assert [basis.profile_unit for basis in bases] == ["count_profiles", "ora_profiles"]
    unchanged = [feature_filter_id(retained)] * 2
    assert [basis.feature_filter_id for basis in bases] == unchanged
    # Relabelling the catalog, or dropping or adding a MediaDive medium, keeps the family.
    assert identities("relabelled", features, relabel=True) == unchanged
    assert identities("without_media", retained) == unchanged
    assert identities("more_media", [*features, NEW_MEDIUM]) == unchanged
    # Adding, dropping or moving a retained feature makes another family.
    for name, order, family in (
        ("added", [*features, NEW_TRAIT], [*retained, NEW_TRAIT]),
        ("moved", [NEW_TRAIT, *features], [NEW_TRAIT, *retained]),
        ("dropped", [f for f in features if f != retained[0]], retained[1:]),
    ):
        observed = identities(name, order)
        assert observed == [feature_filter_id(family)] * 2, name
        assert observed != unchanged, name


@pytest.mark.parametrize("change", ["membership", "order"])
def test_a_basis_whose_profiles_disagree_with_the_catalog_is_refused(
    scoped_runs: dict[StrainRollup, WorkflowRun], tmp_path: Path, change: str
) -> None:
    """Each profile table's retained features must be the catalog's, in its order (#860)."""
    run = scoped_runs[PRIMARY_STRAIN_ROLLUP]
    features = _header(artifact_path(run, "count_profiles"))[3:]
    profiles = [*features, NEW_TRAIT] if change == "membership" else features[::-1]
    copy = _with_features(run, tmp_path / change, features, profile_order=profiles)
    with pytest.raises(ValueError, match="feature catalog"):
        core_analysis_bases(copy)


@pytest.mark.parametrize("scope", STRAIN_ROLLUPS)
def test_arms_inherit_their_core_scope(
    scoped_runs: dict[StrainRollup, WorkflowRun], tmp_path: Path, scope: StrainRollup
) -> None:
    """Every arm plans under the scope its core recorded, at execution and on
    verification, so Noel's ledger reconciles each scope's matrix (#788)."""

    def plan(context: ArmContext) -> list[Command]:
        assert context.strain_rollup == scope
        source, output = context.products["count_profiles"], context.output / "observed.tsv"
        helper = ROOT / "scripts" / "adhoc" / "issue658_arm_test_command.py"
        argv = (sys.executable, str(helper), "--input", str(source), "--output", str(output))
        return [Command("observe", argv, (source,), (output,))]

    core = Path(scoped_runs[scope].manifest_path)
    result = execute_arm("noel", core, tmp_path / "noel", plan, figures=False)

    assert verify_arm(Path(result.manifest_path)) == result


def test_verification_reads_an_absent_ledger_scope_as_off_only_before_788(
    scoped_runs: dict[StrainRollup, WorkflowRun], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A recorded plan may omit the scope only if its core recorded none, as every core
    did before #788; a core that records its scope holds its plans to it (#881)."""

    def plan(context: ArmContext) -> list[Command]:
        source, output = context.products["count_profiles"], context.output / "observed.tsv"
        helper = ROOT / "scripts" / "adhoc" / "issue658_arm_test_command.py"
        argv = (sys.executable, str(helper), "--input", str(source), "--output", str(output))
        return [Command("observe", argv, (source,), (output,))]

    core = scoped_runs["off"]
    result = execute_arm("noel", Path(core.manifest_path), tmp_path / "noel", plan, figures=False)
    validate = arms_runtime._validate_plan
    readings: list[bool] = []

    def recorded(
        commands: Sequence[Command], context: ArmContext, *, absent_scope_is_off: bool = False
    ) -> None:
        readings.append(absent_scope_is_off)
        validate(commands, context, absent_scope_is_off=absent_scope_is_off)

    monkeypatch.setattr(arms_runtime, "_validate_plan", recorded)
    verify_arm(Path(result.manifest_path))
    parameters = [item for item in core.parameters or [] if item.name != "strain_rollup"]
    before_788 = core.model_copy(update={"parameters": parameters})

    def verify_core_before_788(manifest: Path) -> WorkflowRun:
        assert manifest == Path(core.manifest_path)
        return before_788

    monkeypatch.setattr(arms_runtime, "verify_core_run", verify_core_before_788)
    verify_arm(Path(result.manifest_path))

    assert readings == [False, True]


def test_a_subset_core_recorded_before_788_is_read_as_off() -> None:
    """Its producers declared no strain_rollup; any other parameter set still fails."""
    subset = read_definition(ROOT / "workflows" / "bugsigdb_core_subset.yaml")
    legacy = subset.model_copy(deep=True)
    legacy.parameters = [item for item in legacy.parameters or [] if item.name != "strain_rollup"]
    for step in legacy.steps:
        if step.id in {"taxa_matrix", "feature_catalog"}:
            step.parameters = []

    runtime._require_full_source_background(legacy)
    runtime._require_full_source_background(subset)
    next(step for step in legacy.steps if step.id == "taxa_matrix").parameters = ["source_mode"]
    with pytest.raises(ValueError, match="invalid producer bindings for taxa_matrix"):
        runtime._require_full_source_background(legacy)
