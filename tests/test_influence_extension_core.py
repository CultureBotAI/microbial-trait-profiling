"""Synthetic full/subset and off/phenotype execution bindings for issue #771."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from mtp.influence_extension import run_extension
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import TABLE1_TRAITS
from mtp.workflow.consumers import open_core_bundle, verify_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow

SUBSET = "disease_vs_healthy_tier1_2_subject_controls"


def source_inputs(root: Path) -> dict[str, Path]:
    source = REPO_ROOT / "tests/fixtures/workflow"
    root.mkdir()
    edges = (source / "merged-kg_edges.tsv").read_text()
    nodes = (source / "merged-kg_nodes.tsv").read_text()
    nodes += "METPO:2000103\tbiolink:OntologyClass\tcapable of\n"
    for index, term in enumerate(sorted({t for terms in TABLE1_TRAITS.values() for t in terms})):
        if term == "METPO:1000608":
            continue
        predicate = "biolink:has_phenotype" if term.startswith("METPO:") else "METPO:2000103"
        edges += f"focal{index}\tNCBITaxon:101\t{predicate}\t{term}\tfixture:annotation\n"
        nodes += f"{term}\tbiolink:PhenotypicQuality\tFixture trait {index}\n"
    nodes += "METPO:1000678\tbiolink:PhenotypicQuality\toval shaped\n"
    edges += "strain\tkgmicrobe.strain:1\tbiolink:subclass_of\tNCBITaxon:300\tfixture:strain\n"
    edges += "shape\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1000678\tfixture:strain\n"
    paths = {
        name: root / filename
        for name, filename in (
            ("raw_dump", "full_dump.csv"),
            ("kg_edges", "edges.tsv"),
            ("kg_nodes", "nodes.tsv"),
        )
    }
    paths["kg_edges"].write_text(edges)
    paths["kg_nodes"].write_text(nodes)
    with (source / "full_dump.csv").open() as handle:
        reader = csv.DictReader(line for line in handle if not line.startswith("#"))
        template = next(reader)
    for name in (
        "PMID",
        "DOI",
        "Study design",
        "EFO ID",
        "UBERON ID",
        "Group 1 definition",
        "Group 0 sample size",
        "Group 1 sample size",
        "Description",
    ):
        template[name] = ""
    with paths["raw_dump"].open("w", newline="") as handle:
        handle.write("# BugSigDB fixture-771, License: CC0 synthetic test data\n")
        writer = csv.DictWriter(handle, fieldnames=list(template))
        writer.writeheader()
        for study in range(1, 9):
            for number, direction in enumerate(("increased", "decreased"), 1):
                reverse = study == 2
                high = (direction == "increased") != reverse
                row = dict(template)
                row.update(
                    {
                        "BSDB ID": f"bsdb:{study}/1/{number}",
                        "Study": f"Study {study}",
                        "Condition": "Parkinson's disease",
                        "Group 0 name": "Parkinson's disease" if reverse else "healthy controls",
                        "Group 1 name": "healthy controls" if reverse else "Parkinson's disease",
                        "Abundance in Group 1": direction,
                        "NCBI Taxonomy IDs": "101;100" if high else "102;300",
                        "MetaPhlAn taxon names": (
                            "s__Examplegenus_alpha,g__Examplegenus"
                            if high
                            else "s__Examplegenus_beta,s__Unannotated_bacterium"
                        ),
                    }
                )
                if study == 3:
                    row["Group 0 name"], row["Group 1 name"] = (
                        "untreated patients",
                        "treated patients",
                    )
                writer.writerow(row)
    return paths


@pytest.mark.parametrize("scope", ["off", "phenotype"])
@pytest.mark.parametrize("population", ["all_paired", SUBSET])
def test_verified_extension_reconciles_scope_subset_and_orientation(
    tmp_path: Path, scope: str, population: str
) -> None:
    inputs = source_inputs(tmp_path / "inputs")
    subset = population != "all_paired"
    parameters = {"source_mode": "fixture", "strain_rollup": scope}
    if subset:
        parameters["contrast_population"] = population
    spec = "bugsigdb_core_subset.yaml" if subset else "bugsigdb_core.yaml"
    core = execute_workflow(
        REPO_ROOT / "workflows" / spec,
        inputs,
        tmp_path / "core",
        OPERATIONS,
        parameter_overrides=parameters,
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    manifest = Path(core.manifest_path)
    bundle = open_core_bundle(manifest)
    originals = {
        name: file_sha256(bundle.path(name))
        for name in ("count_profiles", "reported_taxa", "taxa_matrix", "experiments")
    }
    output = tmp_path / "extension"
    provenance_path = run_extension(manifest, output)
    provenance = json.loads(provenance_path.read_text())
    record = verify_consumer_record(output / "core_consumption.json")
    assert provenance["core_run_id"] == record.run_id
    assert provenance["strain_rollup"] == scope
    assert provenance["contrast_population"] == population
    assert provenance["reconciliation_mismatches"] == 0
    assert ("METPO:1000678" in provenance["retained_features"]) == (scope == "phenotype")
    with (output / "baseline_fractions.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    pooled = [row for row in rows if row["scope"] == "pooled"]
    assert len(pooled) == 6
    assert {int(row["studies"]) for row in pooled} == ({7} if subset else {8})
    assert {row["family_size"] for row in pooled} == {"6"}
    cellulose = next(row for row in pooled if row["trait_key"].endswith("GO:0030245"))
    # The reversed second study changes from a negative to a positive shift;
    # the treatment-only third study leaves the selected population entirely.
    assert (int(cellulose["positive"]), int(cellulose["negative"])) == (
        (7, 0) if subset else (7, 1)
    )
    for name, digest in originals.items():
        assert file_sha256(bundle.path(name)) == digest
    assert "Native baseline was not supplied" in (output / "native_summary.md").read_text()
    with pytest.raises(FileExistsError):
        run_extension(manifest, output)
    altered = bundle.path("taxa_matrix")
    altered.write_text(altered.read_text() + "\n")
    with pytest.raises(ValueError, match="changed"):
        run_extension(manifest, tmp_path / "must-not-write")
    assert not (tmp_path / "must-not-write").exists()
