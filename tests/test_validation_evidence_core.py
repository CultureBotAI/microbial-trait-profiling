"""Validation labels follow actual verified core KG bytes, without a release claim (#658)."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT
from mtp.trait_guilds import TABLE1_TRAITS
from mtp.validation_evidence import NODES_CORE_VERIFIED, main
from mtp.workflow.consumers import open_core_bundle, verify_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_execution_product
from mtp.workflow.runtime import execute_workflow


@pytest.fixture
def core_manifest(tmp_path: Path) -> Path:
    source = REPO_ROOT / "tests" / "fixtures" / "workflow"
    inputs = tmp_path / "sources"
    inputs.mkdir()
    filenames = {
        "raw_dump": "full_dump.csv",
        "kg_edges": "merged-kg_edges.tsv",
        "kg_nodes": "merged-kg_nodes.tsv",
    }
    for filename in filenames.values():
        shutil.copyfile(source / filename, inputs / filename)
    raw = inputs / filenames["raw_dump"]
    with raw.open(newline="") as stream:
        comment = stream.readline()
        reader = csv.DictReader(stream)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    metadata_fields = ["PMID", "DOI", "Authors list", "Location of subjects", "Group 1 definition"]
    with raw.open("w", newline="") as stream:
        stream.write(comment)
        writer = csv.DictWriter(stream, fieldnames=[*fieldnames, *metadata_fields])
        writer.writeheader()
        writer.writerows({**row, **dict.fromkeys(metadata_fields, "NA")} for row in rows)
    # This consumer needs all six declared traits; extend the source graph before
    # executing the core so every additional matrix cell remains source-validated.
    with (
        (inputs / filenames["kg_nodes"]).open("a") as nodes,
        (inputs / filenames["kg_edges"]).open("a") as edges,
    ):
        nodes.write("METPO:2000103\tbiolink:OntologyClass\tFixture function\n")
        for index, terms in enumerate(TABLE1_TRAITS.values()):
            if terms[0] == "METPO:1000608":
                continue
            nodes.write(f"{terms[0]}\tbiolink:PhenotypicQuality\tFixture trait {index}\n")
            predicate = (
                "biolink:has_phenotype" if terms[0].startswith("METPO:") else "METPO:2000103"
            )
            edges.write(f"extra{index}\tNCBITaxon:101\t{predicate}\t{terms[0]}\tfixture:test\n")
    run = execute_workflow(
        REPO_ROOT / "workflows" / "bugsigdb_core.yaml",
        {key: inputs / filename for key, filename in filenames.items()},
        tmp_path / "core",
        OPERATIONS,
        validator=validate_execution_product,
        repository_root=REPO_ROOT,
    )
    return Path(run.manifest_path)


def test_core_nodes_without_download_sidecar_produce_bound_inventory(
    core_manifest: Path, tmp_path: Path
) -> None:
    output = tmp_path / "validation"
    assert (
        main(["--core-run", str(core_manifest), "--output-dir", str(output), "--require-nodes"])
        == 0
    )
    provenance = json.loads((output / "validation_evidence_provenance.json").read_text())
    bundle = open_core_bundle(core_manifest)
    assert provenance["nodes_status"] == NODES_CORE_VERIFIED
    assert provenance["nodes_input"] == str(bundle.path("kg_nodes"))
    assert provenance["inferred_function_traits_interpretable"] is True
    record = verify_consumer_record(output / "core_consumption.json")
    assert record.consumer == "mtp.validation_evidence"
    assert "source_identity" in record.product_ids
    assert len(provenance["nodes_sha256"]) == 64
    assert (output / "evidence_availability_by_study.tsv").read_text().count("\n") > 1


def test_core_rejects_nodes_from_another_run_before_writing(
    core_manifest: Path, tmp_path: Path
) -> None:
    unrelated = tmp_path / "unrelated.tsv"
    unrelated.write_text("id\tname\nMETPO:1\tunrelated\n")
    output = tmp_path / "validation"
    with pytest.raises(SystemExit) as error:
        main(
            [
                "--core-run",
                str(core_manifest),
                "--nodes-input",
                str(unrelated),
                "--output-dir",
                str(output),
            ]
        )
    assert error.value.code == 2
    assert not output.exists()


def test_core_rejects_changed_nodes_before_writing(core_manifest: Path, tmp_path: Path) -> None:
    nodes = open_core_bundle(core_manifest).path("kg_nodes")
    nodes.write_text(nodes.read_text().replace("facultatively anaerobic", "fermentation"))
    output = tmp_path / "validation"
    with pytest.raises(ValueError, match="content or inventory changed"):
        main(["--core-run", str(core_manifest), "--output-dir", str(output)])
    assert not output.exists()


def test_core_output_cannot_overwrite_run(core_manifest: Path) -> None:
    with pytest.raises(ValueError, match="core"):
        main(["--core-run", str(core_manifest), "--output-dir", str(core_manifest.parent)])


def test_snapshot_regeneration_cannot_claim_core_execution(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as error:
        main(["--core-run", str(tmp_path / "run.json"), "--snapshot-from", "snapshot.tsv"])
    assert error.value.code == 2
