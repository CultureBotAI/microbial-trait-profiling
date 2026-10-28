"""Fault-injection checks for raw and flattened source boundaries (#651, 2026-10-19)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
import yaml

from mtp.workflow import source_contracts
from mtp.workflow.definition import read_definition
from mtp.workflow.operations import (
    OPERATIONS,
    export_tables_operation,
    ingest_operation,
    validate_core_product,
    validate_execution_product,
)
from mtp.workflow.runtime import WorkflowExecutionError, artifact_path, execute_workflow, inventory
from mtp.workflow.verification import verify_core_run

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "workflows" / "bugsigdb_core.yaml"
FIXTURE = ROOT / "tests" / "fixtures" / "workflow"
PRODUCTS = {product.id: product for product in read_definition(SPEC).products}
EXPORTS = (
    "study_table",
    "experiments",
    "experiment_ontology_terms",
    "signatures",
    "reported_taxa",
    "taxon_profiles_wide",
    "taxon_profiles_long",
)


@pytest.fixture
def source_bundle(tmp_path: Path) -> dict[str, Path]:
    paths = {
        "raw_dump": FIXTURE / "full_dump.csv",
        "studies": tmp_path / "studies",
        "ingest_log": tmp_path / "ingest_log.json",
        **{name: tmp_path / PRODUCTS[name].filename for name in EXPORTS},
    }
    ingest_operation(
        {"raw_dump": paths["raw_dump"]},
        {name: paths[name] for name in ("studies", "ingest_log")},
        {},
    )
    export_tables_operation(
        {"studies": paths["studies"]}, {name: paths[name] for name in EXPORTS}, {}
    )
    return paths


def _validate_ingest(paths: dict[str, Path]) -> None:
    validate_core_product(
        PRODUCTS["studies"],
        paths["studies"],
        paths,
        trait_key_policy="raw_object_v1",
        include_policy_metadata=False,
    )


def _validate_exports(paths: dict[str, Path]) -> None:
    validate_core_product(
        PRODUCTS["study_table"],
        paths["study_table"],
        paths,
        trait_key_policy="raw_object_v1",
        include_policy_metadata=False,
    )


@pytest.mark.parametrize(
    "mutation", ["study", "signature", "taxon", "blank_name", "placeholder_name", "host"]
)
def test_raw_source_rejects_schema_valid_ingest_loss_or_corruption(
    source_bundle: dict[str, Path], mutation: str
) -> None:
    path = source_bundle["studies"] / "bsdb_1.yaml"
    record = yaml.safe_load(path.read_text(encoding="utf-8"))
    signature = record["experiments"][0]["signatures"][0]
    if mutation == "study":
        path.unlink()
    else:
        if mutation == "signature":
            record["experiments"][0]["signatures"].pop(0)
        elif mutation == "taxon":
            signature["taxa"].pop(0)
        elif mutation == "host":
            record["experiments"][0]["host_species"] = "Mus musculus"
        else:
            signature["taxa"][0]["taxon_name"] = "" if mutation == "blank_name" else "taxon_100"
        path.write_text(yaml.safe_dump(record), encoding="utf-8")
    with pytest.raises(ValueError, match=r"source reconciliation|taxon name"):
        _validate_ingest(source_bundle)


@pytest.mark.parametrize(
    ("product", "column"),
    [
        ("study_table", "title"),
        ("experiments", "condition"),
        ("reported_taxa", "taxon_name"),
        ("taxon_profiles_wide", "Examplegenus"),
        ("taxon_profiles_long", "taxon"),
    ],
)
def test_exports_reconcile_values_and_repeated_observations_to_study_source(
    source_bundle: dict[str, Path], product: str, column: str
) -> None:
    path = source_bundle[product]
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    rows[0][column] = "999" if product == "taxon_profiles_wide" else "altered source value"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="source reconciliation"):
        _validate_exports(source_bundle)


def test_unchanged_fixture_source_boundaries_preserve_duplicate_observations(
    source_bundle: dict[str, Path],
) -> None:
    _validate_ingest(source_bundle)
    _validate_exports(source_bundle)


def _drop_study(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    ingest_operation(inputs, outputs, parameters)
    (outputs["studies"] / "bsdb_2.yaml").unlink()


def test_ingest_source_failure_prevents_stage_publication(tmp_path: Path) -> None:
    with pytest.raises(WorkflowExecutionError, match="source reconciliation") as failure:
        execute_workflow(
            SPEC,
            {
                "raw_dump": FIXTURE / "full_dump.csv",
                "kg_edges": FIXTURE / "merged-kg_edges.tsv",
                "kg_nodes": FIXTURE / "merged-kg_nodes.tsv",
            },
            tmp_path / "runs",
            {**OPERATIONS, "ingest": _drop_study},
            validator=validate_execution_product,
            repository_root=ROOT,
        )
    run = failure.value.run
    assert run.status == "failed"
    assert "studies" not in {artifact.product_id for artifact in run.artifacts or []}
    assert not (Path(run.manifest_path).parent / "steps" / "ingest").exists()


def test_fabricated_ingest_skip_cannot_justify_a_dropped_valid_signature(
    source_bundle: dict[str, Path],
) -> None:
    path = source_bundle["ingest_log"]
    events = json.loads(path.read_text(encoding="utf-8"))
    events.append(
        {
            "level": "skip",
            "entity": "signature",
            "id": "bsdb:1-1-1",
            "reason": "no parseable taxa",
        }
    )
    path.write_text(json.dumps(events), encoding="utf-8")
    with pytest.raises(ValueError, match="source reconciliation"):
        _validate_ingest(source_bundle)


@pytest.mark.parametrize("remove_skip", [False, True])
def test_native_normalization_and_explicit_source_empty_exclusions(
    tmp_path: Path, remove_skip: bool
) -> None:
    with (FIXTURE / "full_dump.csv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(line for line in handle if not line.startswith("#"))
        fields = [*(reader.fieldnames or []), "DOI", "URL", "EFO ID", "PMID"]
        rows = list(reader)
    rows[0].update(
        {
            "DOI": " https://doi.org/10.1234/example ",
            "URL": " https://example.org/study ",
            "EFO ID": "EFO:0001, NA, MONDO:0002",
            "PMID": "123",
            "NCBI Taxonomy IDs": "2|100;2|100;101",
            "MetaPhlAn taxon names": (
                "k__Bacteria|g__Examplegenus,k__Bacteria|g__Examplegenus,s__Examplegenus_alpha"
            ),
        }
    )
    rows.append({**rows[0], "BSDB ID": "bsdb:99/1/1", "NCBI Taxonomy IDs": "NA"})
    rows.append({**rows[0], "BSDB ID": "bsdb:99/2/1", "Abundance in Group 1": "NA"})
    raw = tmp_path / "raw.csv"
    with raw.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    paths = {"raw_dump": raw, "studies": tmp_path / "studies", "ingest_log": tmp_path / "log.json"}
    ingest_operation(
        {"raw_dump": raw}, {name: paths[name] for name in ("studies", "ingest_log")}, {}
    )
    if remove_skip:
        events = json.loads(paths["ingest_log"].read_text(encoding="utf-8"))
        omitted = next(index for index, event in enumerate(events) if event["level"] == "skip")
        events.pop(omitted)
        paths["ingest_log"].write_text(json.dumps(events), encoding="utf-8")
        with pytest.raises(ValueError, match="skip events"):
            _validate_ingest(paths)
    else:
        _validate_ingest(paths)


def test_revalidation_checks_raw_source_even_with_updated_artifact_inventory(
    tmp_path: Path,
) -> None:
    run = execute_workflow(
        SPEC,
        {
            "raw_dump": FIXTURE / "full_dump.csv",
            "kg_edges": FIXTURE / "merged-kg_edges.tsv",
            "kg_nodes": FIXTURE / "merged-kg_nodes.tsv",
        },
        tmp_path / "runs",
        OPERATIONS,
        validator=validate_execution_product,
        repository_root=ROOT,
    )
    studies = artifact_path(run, "studies")
    (studies / "bsdb_2.yaml").unlink()
    artifact = next(item for item in run.artifacts or [] if item.product_id == "studies")
    artifact.sha256, artifact.size_bytes, artifact.inventory = inventory(studies, "collection")
    Path(run.manifest_path).write_text(run.model_dump_json(), encoding="utf-8")
    with pytest.raises(ValueError, match="source reconciliation"):
        verify_core_run(Path(run.manifest_path))


def test_a_table_without_sequencing_type_is_named_as_predating_742(tmp_path: Path) -> None:
    """A core export written before #742 names the missing column and the cause (#762)."""
    table = tmp_path / "experiments.tsv"
    table.write_text("study_id\texperiment_id\nbsdb:1\tbsdb:1-1\n", encoding="utf-8")
    with pytest.raises(
        ValueError, match="missing columns: sequencing_type; the core predates #742"
    ):
        source_contracts._compare_rows(
            table, ["study_id", "experiment_id", "sequencing_type"], [["bsdb:1", "bsdb:1-1", "16S"]]
        )
