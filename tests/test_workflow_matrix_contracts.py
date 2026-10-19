"""Reject scientifically corrupted matrix products before publication (#652)."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from mtp.workflow.operations import (
    OPERATIONS,
    feature_catalog_operation,
    taxa_matrix_operation,
    validate_core_product,
)
from mtp.workflow.runtime import WorkflowExecutionError, artifact_path, execute_workflow
from mtp.workflow.verification import verify_core_run

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "workflow"
TRAIT = "METPO:1000608"


def _table(path: Path) -> list[list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def _write(path: Path, rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)


@pytest.mark.parametrize(
    "damage",
    [
        "inflated_count",
        "erased_count",
        "invented_count",
        "dropped_feature",
        "taxon_label",
        "taxon_rank",
        "feature_label",
    ],
)
def test_matrix_source_corruption_fails_before_publication(tmp_path: Path, damage: str) -> None:
    """Internally consistent downstream profiles cannot legitimize a wrong matrix."""

    def faulty_matrix(
        inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
    ) -> None:
        taxa_matrix_operation(inputs, outputs, parameters)
        matrix = _table(outputs["taxa_matrix"])
        metadata = _table(outputs["feature_metadata"])
        feature = matrix[0].index(TRAIT)
        genus = next(row for row in matrix[1:] if row[0] == "NCBITaxon:100")
        unannotated = next(row for row in matrix[1:] if row[0] == "NCBITaxon:300")
        if damage == "inflated_count":
            genus[feature] = "4"
        elif damage == "erased_count":
            genus[feature] = "0"
        elif damage == "invented_count":
            unannotated[feature] = "1"
        elif damage == "dropped_feature":
            for row in matrix:
                del row[feature]
            metadata = [row for row in metadata if row[0] != TRAIT]
        elif damage == "taxon_label":
            genus[1] = "Incorrect biological identity"
        elif damage == "taxon_rank":
            genus[2] = "species"
        else:
            next(row for row in metadata[1:] if row[0] == TRAIT)[1] = "Wrong phenotype label"
        _write(outputs["taxa_matrix"], matrix)
        _write(outputs["feature_metadata"], metadata)

    operations = {**OPERATIONS, "taxa_matrix": faulty_matrix}
    with pytest.raises(WorkflowExecutionError, match=r"(KG|reported taxon)") as failure:
        execute_workflow(
            ROOT / "workflows" / "bugsigdb_core.yaml",
            {
                "raw_dump": FIXTURE / "full_dump.csv",
                "kg_edges": FIXTURE / "merged-kg_edges.tsv",
                "kg_nodes": FIXTURE / "merged-kg_nodes.tsv",
            },
            tmp_path / "runs",
            operations,
            parameter_overrides={"source_mode": "fixture"},
            validator=validate_core_product,
            repository_root=ROOT,
        )
    assert "taxa_matrix" not in {
        artifact.product_id for artifact in failure.value.run.artifacts or []
    }


def test_catalog_cannot_change_source_feature_meaning(tmp_path: Path) -> None:
    def faulty_catalog(
        inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
    ) -> None:
        feature_catalog_operation(inputs, outputs, parameters)
        rows = _table(outputs["feature_catalog"])
        name = rows[0].index("object_name")
        next(row for row in rows[1:] if row[0] == TRAIT)[name] = "Unrelated biological trait"
        _write(outputs["feature_catalog"], rows)

    with pytest.raises(WorkflowExecutionError, match="source feature metadata"):
        execute_workflow(
            ROOT / "workflows" / "bugsigdb_core.yaml",
            {
                "raw_dump": FIXTURE / "full_dump.csv",
                "kg_edges": FIXTURE / "merged-kg_edges.tsv",
                "kg_nodes": FIXTURE / "merged-kg_nodes.tsv",
            },
            tmp_path / "runs",
            {**OPERATIONS, "feature_catalog": faulty_catalog},
            validator=validate_core_product,
            repository_root=ROOT,
        )


def test_direct_higher_rank_evidence_and_multiple_paths_reconcile(tmp_path: Path) -> None:
    """A direct genus annotation adds one source; a second taxonomy path adds none."""
    edges = tmp_path / "edges.tsv"
    edges.write_text(
        (FIXTURE / "merged-kg_edges.tsv").read_text(encoding="utf-8")
        + "direct_genus\tNCBITaxon:100\tbiolink:has_phenotype\tMETPO:1000608\tfixture:direct\n"
        + "second_path\tNCBITaxon:101\tbiolink:subclass_of\tNCBITaxon:900\tfixture:taxonomy\n",
        encoding="utf-8",
    )
    run = execute_workflow(
        ROOT / "workflows" / "bugsigdb_core.yaml",
        {
            "raw_dump": FIXTURE / "full_dump.csv",
            "kg_edges": edges,
            "kg_nodes": FIXTURE / "merged-kg_nodes.tsv",
        },
        tmp_path / "runs",
        OPERATIONS,
        validator=validate_core_product,
        repository_root=ROOT,
    )
    verify_core_run(Path(run.manifest_path))
    matrix = _table(artifact_path(run, "taxa_matrix"))
    feature = matrix[0].index(TRAIT)
    assert {row[0]: row[feature] for row in matrix[1:]} == {
        "NCBITaxon:900": "4",
        "NCBITaxon:100": "4",
        "NCBITaxon:101": "1",
        "NCBITaxon:103": "1",
        "NCBITaxon:200": "0",
        "NCBITaxon:999": "0",
        "NCBITaxon:300": "0",
    }
