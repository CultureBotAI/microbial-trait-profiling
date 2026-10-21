"""Offline raw-to-core scientific and provenance integration tests for issue #629."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import cast

import pytest
import yaml

from mtp.count_input import read_profile_coverage
from mtp.ingest import ingest
from mtp.study_trait_profiles import collect_profile_taxa, read_taxa_trait_matrix
from mtp.taxa_trait_matrix import read_kg_edges
from mtp.workflow.consumers import build_consumer_coverage, open_core_bundle
from mtp.workflow.coverage import build_profile_coverage
from mtp.workflow.definition import read_definition, validate_definition
from mtp.workflow.kg_inputs import audit_kg_edges, iter_kg_edges
from mtp.workflow.models import WorkflowRun
from mtp.workflow.operations import (
    OPERATIONS,
    identify_sources,
    ingest_operation,
    ora_profiles_operation,
    trace_profile_feature,
    validate_core_product,
)
from mtp.workflow.runtime import (
    WorkflowExecutionError,
    artifact_path,
    execute_workflow,
    validate_run,
)

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "workflows" / "bugsigdb_core.yaml"
FIXTURE = Path(__file__).parent / "fixtures" / "workflow"
INPUTS = {
    "raw_dump": FIXTURE / "full_dump.csv",
    "kg_edges": FIXTURE / "merged-kg_edges.tsv",
    "kg_nodes": FIXTURE / "merged-kg_nodes.tsv",
}
TRAIT = "METPO:1000608"
MEDIUM = "METPO:2000510|mediadive.medium:1"
OVERLAP_KEY = ("bsdb:1", "bsdb:1-1", "inc")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _profiles(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    return {(r["study_id"], r["experiment_id"], r["inc_dec"]): r for r in _rows(path)}


@pytest.fixture(scope="module")
def completed_core(tmp_path_factory: pytest.TempPathFactory) -> WorkflowRun:
    return execute_workflow(
        SPEC,
        INPUTS,
        tmp_path_factory.mktemp("core-execution"),
        OPERATIONS,
        parameter_overrides={"source_mode": "fixture"},
        validator=validate_core_product,
        repository_root=ROOT,
    )


def test_raw_core_spec_declares_scientific_units_and_dependencies() -> None:
    definition = read_definition(SPEC)
    order = validate_definition(definition)
    assert order.index("ingest") < order.index("export_tables") < order.index("taxa_matrix")
    assert order.index("taxa_matrix") < order.index("count_profiles")
    products = {product.id: product for product in definition.products}
    assert products["count_profiles"].unit == "rolled_up_source_taxon_count"
    assert products["ora_profiles"].unit == "distinct_matched_reported_taxon_count"
    assert products["ora_profiles"].cohort_parameter == "ora_host_species"
    assert {p.id for p in definition.products if p.external} == set(INPUTS)


@pytest.mark.parametrize(
    ("product_id", "field", "value"),
    [
        ("taxa_matrix", "unit", "physiologic_cell_density"),
        ("taxa_matrix", "feature_policy", "exclude_mediadive"),
        ("count_profiles", "coverage_policy", "unknown"),
        ("count_profiles", "unit", "distinct_matched_reported_taxon_count"),
        ("ora_profiles", "schema_ref", "schema/core_products.yaml#CountProfileRow"),
    ],
)
def test_false_scientific_descriptors_fail_before_product_publication(
    tmp_path: Path, product_id: str, field: str, value: str
) -> None:
    definition = yaml.safe_load(SPEC.read_text(encoding="utf-8"))
    product = next(item for item in definition["products"] if item["id"] == product_id)
    product[field] = value
    changed_spec = tmp_path / "changed-spec.yaml"
    changed_spec.write_text(yaml.safe_dump(definition), encoding="utf-8")
    with pytest.raises(WorkflowExecutionError, match=f"declared {field}") as failure:
        execute_workflow(
            changed_spec,
            INPUTS,
            tmp_path / "runs",
            OPERATIONS,
            validator=validate_core_product,
            repository_root=ROOT,
        )
    run = failure.value.run
    assert run.status == "failed"
    assert product_id not in {artifact.product_id for artifact in run.artifacts or []}


@pytest.mark.parametrize("step_id", ["identify_sources", "count_profiles", "ora_profiles"])
def test_declared_but_unimplemented_operation_parameter_fails(tmp_path: Path, step_id: str) -> None:
    definition = yaml.safe_load(SPEC.read_text(encoding="utf-8"))
    definition["parameters"].append({"name": "annotation_rollup", "default": "direct_only"})
    step = next(item for item in definition["steps"] if item["id"] == step_id)
    step["parameters"].append("annotation_rollup")
    changed_spec = tmp_path / "changed-parameter.yaml"
    changed_spec.write_text(yaml.safe_dump(definition), encoding="utf-8")
    with pytest.raises(WorkflowExecutionError, match="requires exactly"):
        execute_workflow(
            changed_spec,
            INPUTS,
            tmp_path / "runs",
            OPERATIONS,
            validator=validate_core_product,
            repository_root=ROOT,
        )


def test_raw_to_core_preserves_count_ora_denominator_and_cohorts(
    completed_core: WorkflowRun,
) -> None:
    assert completed_core.status == "completed"
    validate_run(completed_core.manifest_path)
    counts = _profiles(artifact_path(completed_core, "count_profiles"))
    ora = _profiles(artifact_path(completed_core, "ora_profiles"))
    assert counts[OVERLAP_KEY][TRAIT] == "4"
    assert ora[OVERLAP_KEY][TRAIT] == "2"
    assert ora[OVERLAP_KEY]["taxa_count"] == "2"
    assert len(counts) == 9
    assert len(ora) == 7
    assert ("bsdb:2", "bsdb:2-1", "inc") in counts
    assert ("bsdb:2", "bsdb:2-1", "inc") not in ora
    # Core generation preserves incomplete pairs and culture-medium annotations.
    assert ("bsdb:1", "bsdb:1-3", "inc") in ora
    assert ("bsdb:1", "bsdb:1-3", "dec") not in ora
    assert counts[("bsdb:1", "bsdb:1-2", "inc")][MEDIUM] == "1"
    assert ora[("bsdb:1", "bsdb:1-2", "inc")][MEDIUM] == "1"
    occurrences = _rows(artifact_path(completed_core, "reported_taxa"))
    assert sum(row["taxon_id"] == "NCBITaxon:100" for row in occurrences) == 2
    studies = artifact_path(completed_core, "studies")
    first = yaml.safe_load((studies / "bsdb_1.yaml").read_text(encoding="utf-8"))
    assert first["experiments"][0]["signatures"][0]["taxa"][0]["taxon_name"] == "Examplegenus"


def test_source_identity_and_catalog_do_not_claim_unverified_pin_or_physiology(
    completed_core: WorkflowRun,
) -> None:
    identity = json.loads(artifact_path(completed_core, "source_identity").read_text())
    assert identity["declared_source_mode"] == "fixture"
    assert identity["pinned_verified"] is False
    assert identity["verification_status"] == "origin_unverified"
    assert "verified_release" not in identity
    assert set(identity["consumed_files"]) == set(INPUTS)
    catalog = {
        row["column_id"]: row for row in _rows(artifact_path(completed_core, "feature_catalog"))
    }
    assert catalog[MEDIUM]["feature_class"] == "culture_medium"
    assert catalog[TRAIT]["feature_class"] == "trait_annotation"
    assert catalog[TRAIT]["annotation_evidence"] == "direct_or_descendant_unresolved"


def test_coverage_separates_matrix_match_from_zero_annotations(
    completed_core: WorkflowRun,
) -> None:
    coverage = read_profile_coverage(artifact_path(completed_core, "ora_coverage"))
    assert coverage[OVERLAP_KEY].reported_distinct_taxa == 2
    assert coverage[OVERLAP_KEY].lookup_matched_taxa == 2
    for key in (("bsdb:1", "bsdb:1-1", "dec"), ("bsdb:1", "bsdb:1-2", "dec")):
        # One KG-present/all-zero taxon and one KG-absent taxon both have matrix rows.
        assert coverage[key].coverage_status == "matched_without_annotations"
        assert coverage[key].lookup_matched_taxa == 1
        assert coverage[key].lookup_unmatched_taxa == 0
    media_key = ("bsdb:1", "bsdb:1-2", "inc")
    assert coverage[media_key].coverage_status == "complete_annotation_coverage"
    matrix = read_taxa_trait_matrix(artifact_path(completed_core, "taxa_matrix"))
    profiles, _ = collect_profile_taxa(artifact_path(completed_core, "studies"))
    filtered = build_profile_coverage(profiles, matrix, feature_columns=[TRAIT])
    assert filtered[media_key].coverage_status == "matched_without_annotations"
    assert filtered[media_key].trait_filter_id != coverage[media_key].trait_filter_id


def test_selected_feature_source_ledger_exposes_rank_overlap(completed_core: WorkflowRun) -> None:
    result = trace_profile_feature(
        artifact_path(completed_core, "studies"),
        artifact_path(completed_core, "taxa_matrix"),
        INPUTS["kg_edges"],
        INPUTS["kg_nodes"],
        OVERLAP_KEY,
        TRAIT,
    )
    assert result["count_profile_value"] == 4
    assert result["ora_profile_value"] == result["ora_taxa_count"] == 2
    assert result["unique_annotation_source_taxa"] == 3
    assert result["overlapping_source_taxa"] == ["NCBITaxon:101"]
    assert result["source_taxon_multiplicity"] == {
        "NCBITaxon:101": 2,
        "NCBITaxon:102": 1,
        "NCBITaxon:103": 1,
    }
    ledger = cast(list[dict[str, object]], result["ledger"])
    a_rows = [row for row in ledger if row["source_taxon_id"] == "NCBITaxon:101"]
    assert {row["annotation_origin"] for row in a_rows} == {"direct", "descendant"}
    assert all(len(cast(list[object], row["evidence_edges"])) == 2 for row in a_rows)


@pytest.mark.parametrize(
    ("experiment", "expected"), [("bsdb:1-1", []), ("bsdb:1-2", ["NCBITaxon:999"])]
)
def test_source_ledger_distinguishes_present_zero_from_kg_absent(
    completed_core: WorkflowRun, experiment: str, expected: list[str]
) -> None:
    result = trace_profile_feature(
        artifact_path(completed_core, "studies"),
        artifact_path(completed_core, "taxa_matrix"),
        INPUTS["kg_edges"],
        INPUTS["kg_nodes"],
        ("bsdb:1", experiment, "dec"),
        TRAIT,
    )
    assert result["reported_taxa_absent_from_kg"] == expected
    assert result["matrix_unmatched_taxon_ids"] == []
    assert result["ora_taxa_count"] == 1
    assert result["count_profile_value"] == result["ora_profile_value"] == 0


def test_matrix_unmatched_taxon_remains_explicit_in_coverage(
    completed_core: WorkflowRun,
    tmp_path: Path,
) -> None:
    original = artifact_path(completed_core, "taxa_matrix")
    missing_row_matrix = tmp_path / "matrix.tsv"
    missing_row_matrix.write_text(
        "\n".join(
            line
            for line in original.read_text().splitlines()
            if not line.startswith("NCBITaxon:999\t")
        )
        + "\n"
    )
    outputs = {"ora_profiles": tmp_path / "profiles.tsv", "ora_coverage": tmp_path / "coverage.tsv"}
    ora_profiles_operation(
        {"studies": artifact_path(completed_core, "studies"), "taxa_matrix": missing_row_matrix},
        outputs,
        {"ora_host_species": "Homo sapiens"},
    )
    key = ("bsdb:1", "bsdb:1-2", "dec")
    row = _profiles(outputs["ora_profiles"])[key]
    assert row["taxa_count"] == row[TRAIT] == "0"
    coverage = read_profile_coverage(outputs["ora_coverage"])[key]
    assert coverage.coverage_status == "no_lookup_matches"
    assert coverage.unmatched_taxon_ids == ("NCBITaxon:999",)


def test_host_override_changes_bound_cohort_and_row_membership(tmp_path: Path) -> None:
    run = execute_workflow(
        SPEC,
        INPUTS,
        tmp_path / "runs",
        OPERATIONS,
        parameter_overrides={"count_host_species": "Homo sapiens", "ora_host_species": "all"},
        validator=validate_core_product,
        repository_root=ROOT,
    )
    assert len(_rows(artifact_path(run, "count_profiles"))) == 7
    assert len(_rows(artifact_path(run, "ora_profiles"))) == 9
    artifacts = {artifact.product_id: artifact for artifact in run.artifacts or []}
    assert (
        artifacts["count_profiles"].cohort == artifacts["count_coverage"].cohort == "Homo sapiens"
    )
    assert artifacts["ora_profiles"].cohort == artifacts["ora_coverage"].cohort == "all"


def test_ingest_uses_explicit_snapshot_and_rejects_stale_collection(tmp_path: Path) -> None:
    stale = tmp_path / "studies"
    stale.mkdir()
    (stale / "stale.yaml").write_text("id: bsdb:404\n")
    with pytest.raises(ValueError, match="stale study YAML"):
        ingest_operation(
            {"raw_dump": INPUTS["raw_dump"]},
            {"studies": stale, "ingest_log": tmp_path / "log.json"},
            {},
        )
    # Explicit input must win over an unrelated default full_dump.csv nearby.
    (tmp_path / "full_dump.csv").write_text("not,the,fixture\n")
    written = ingest(
        raw_dir=tmp_path,
        study_dir=tmp_path / "fresh",
        derived_dir=tmp_path / "derived",
        full_dump_path=INPUTS["raw_dump"],
        log_path=tmp_path / "events" / "ingest.json",
    )
    assert len(written) == 3
    assert (tmp_path / "events" / "ingest.json").is_file()


def test_unverified_fixture_cannot_be_declared_pinned(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="does not match verified release bytes"):
        identify_sources(
            INPUTS, {"source_identity": tmp_path / "identity.json"}, {"source_mode": "pinned"}
        )


def test_known_raw_kg_continuations_preserve_native_numeric_inputs(tmp_path: Path) -> None:
    """The production export's split metadata cannot change taxonomy/trait counting."""
    malformed = tmp_path / "raw.tsv"
    malformed.write_bytes(
        b"id\tsubject\tpredicate\tobject\trelation\tagent_type\tagent_type\r"
        b"\tknowledge_level\tknowledge_source\tmeta\tprimary_knowledge_source\n"
        b"e1\tNCBITaxon:101\tbiolink:subclass_of\tNCBITaxon:100\t\tmanual_agent\t"
        b"\tknowledge_assertion\tfixture\t\tfixture:taxonomy\n"
        b"e2\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000608\t\tautomated_agent\t"
        b"manual_validation_of_automated_agent\r\tknowledge_assertion\tfixture\t\tfixture:annotation\n"
    )
    audit = audit_kg_edges(malformed)
    assert audit["records"] == 2
    assert audit["header_continuations_reconstructed"] == 1
    assert audit["agent_field_continuations_reconstructed"] == 1
    assert audit["excluded_records"] == 0
    canonical = list(iter_kg_edges(malformed))
    assert canonical[1]["agent_type"] == "automated_agent"
    assert canonical[1]["agent_type_secondary"] == "manual_validation_of_automated_agent"
    assert canonical[1]["primary_knowledge_source"] == "fixture:annotation"
    normalized = tmp_path / "normalized.tsv"
    with normalized.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["id", "subject", "predicate", "object"])
        writer.writerows(
            [row[name] for name in ("id", "subject", "predicate", "object")] for row in canonical
        )
    assert read_kg_edges(malformed) == read_kg_edges(normalized)


def test_unknown_malformed_raw_kg_record_fails_explicitly(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("id\tsubject\tpredicate\tobject\ne1\tNCBITaxon:101\tMETPO:2000103\n")
    with pytest.raises(ValueError, match="unknown malformed KG edge"):
        audit_kg_edges(path)


def test_only_16s_and_shotgun_experiments_enter_the_profiles(
    completed_core: WorkflowRun, tmp_path: Path
) -> None:
    """An ITS experiment stays in the exported tables but not in either profile product,
    and the profile and coverage contracts rebuild the same cohort (#742)."""
    raw = tmp_path / "full_dump.csv"
    text = INPUTS["raw_dump"].read_text(encoding="utf-8")
    assert text.count(",WMS\n") == 2
    raw.write_text(text.replace(",WMS\n", ",ITS / ITS2\n"), encoding="utf-8")
    run = execute_workflow(
        SPEC,
        {**INPUTS, "raw_dump": raw},
        tmp_path / "runs",
        OPERATIONS,
        validator=validate_core_product,
        repository_root=ROOT,
    )
    experiments = {r["experiment_id"]: r for r in _rows(artifact_path(run, "experiments"))}
    assert experiments["bsdb:3-1"]["sequencing_type"] == "ITS / ITS2"
    for product in ("count_profiles", "ora_profiles", "count_coverage", "ora_coverage"):
        before = _rows(artifact_path(completed_core, product))
        after = _rows(artifact_path(run, product))
        assert any(row["study_id"] == "bsdb:3" for row in before)
        assert after == [row for row in before if row["study_id"] != "bsdb:3"]
    # Consumer coverage and feature lineage apply the same rule (#763)
    coverage = read_profile_coverage(
        build_consumer_coverage(open_core_bundle(Path(run.manifest_path)), tmp_path / "cov.tsv")
    )
    assert coverage
    assert {key[0] for key in coverage} == {
        row["study_id"] for row in _rows(artifact_path(run, "count_profiles"))
    }
    with pytest.raises(ValueError, match="absent from selected study cohort"):
        trace_profile_feature(
            artifact_path(run, "studies"),
            artifact_path(run, "taxa_matrix"),
            INPUTS["kg_edges"],
            INPUTS["kg_nodes"],
            ("bsdb:3", "bsdb:3-1", "inc"),
            TRAIT,
        )
