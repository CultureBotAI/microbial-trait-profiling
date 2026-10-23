"""Core consumers must not import an unrelated checkout KG identity (#629)."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

from mtp import corpus_identity, table1_trait_maps, trait_composition
from mtp.paths import REPO_ROOT
from mtp.trait_guilds import TABLE1_TRAITS
from mtp.workflow.consumers import verify_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import artifact_path, execute_workflow


def _render_placeholder(*args: object, **_kwargs: object) -> Path:
    """Rendering is unrelated to the CLI provenance contract under test."""
    path = args[1]
    assert isinstance(path, Path)
    path.write_text("fixture rendering\n", encoding="utf-8")
    return path


@pytest.mark.parametrize("consumer", ["composition", "table1_maps"])
def test_core_consumer_identity_comes_only_from_bound_kg(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, consumer: str
) -> None:
    source = REPO_ROOT / "tests" / "fixtures" / "workflow"
    edges = tmp_path / "custom-kg.tsv"
    edges.write_bytes((source / "merged-kg_edges.tsv").read_bytes())
    with edges.open("a", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        for index, identifiers in enumerate(TABLE1_TRAITS.values()):
            term = identifiers[0]
            if term == "METPO:1000608":
                continue
            # Phenotype classes are has_phenotype assertions. The fixture must
            # not depend on collapsing a capable_of relation into a class key.
            predicate = "biolink:has_phenotype" if term.startswith("METPO:") else "METPO:2000103"
            writer.writerow([f"extra-{index}", "NCBITaxon:103", predicate, term, "fixture:custom"])
    run = execute_workflow(
        REPO_ROOT / "workflows" / "bugsigdb_core.yaml",
        {
            "raw_dump": source / "full_dump.csv",
            "kg_edges": edges,
            "kg_nodes": source / "merged-kg_nodes.tsv",
        },
        tmp_path / "runs",
        OPERATIONS,
        parameter_overrides={"source_mode": "fixture"},
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    with artifact_path(run, "taxa_matrix").open(encoding="utf-8", newline="") as stream:
        features = set(next(csv.reader(stream, delimiter="\t")))
    assert {
        "METPO:1000668",
        "METPO:1000683",
        "METPO:1000623",
        "METPO:2000103|GO:0030245",
        "METPO:2000103|pathways:arsenate_reduction",
    } <= features
    legacy_reads: list[Path] = []

    def unrelated_checkout_identity(path: Path) -> tuple[str, list[corpus_identity.KgInput]]:
        legacy_reads.append(path)
        return "checkout-only-release", []

    monkeypatch.setattr(corpus_identity, "read_kg_inputs", unrelated_checkout_identity)
    output = tmp_path / consumer
    arguments = ["--core-run", run.manifest_path, "--output-dir", str(output)]
    if consumer == "composition":
        assert trait_composition.main(arguments) == 0
        provenance = json.loads((output / "trait_log_ratio_provenance.json").read_text())
        identity = provenance["corpus_identity"]
    else:
        monkeypatch.setattr(table1_trait_maps, "write_figure", _render_placeholder)
        monkeypatch.setattr(table1_trait_maps, "write_shift_figure", _render_placeholder)
        monkeypatch.setattr(table1_trait_maps, "FOCUS_PANELS", ())
        monkeypatch.setattr(
            sys,
            "argv",
            ["table1_trait_maps", *arguments, "--min-studies", "1", "--min-site-conditions", "1"],
        )
        assert table1_trait_maps.main() == 0
        identity = json.loads((output / "corpus_identity.json").read_text())
    assert legacy_reads == []
    kg = identity["kg_microbe"]
    assert kg["release_tag"] is None
    assert kg["provenance_path"] == str(artifact_path(run, "source_identity"))
    artifacts = {artifact.product_id: artifact for artifact in run.artifacts or []}
    assert {item["sha256"] for item in kg["inputs"]} == {
        artifacts["kg_edges"].sha256,
        artifacts["kg_nodes"].sha256,
    }
    record = verify_consumer_record(output / "core_consumption.json")
    assert "source_identity" in record.product_ids
