"""Synthetic full/subset and off/phenotype execution bindings for issue #771."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

from mtp import influence_extension as extension
from mtp.influence_extension import run_extension
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import TABLE1_TRAITS
from mtp.workflow.arms.runtime import execute_arm, verify_arm
from mtp.workflow.arms.types import ArmContext, Command
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
    assert {
        "count_profiles",
        "experiments",
        "reported_taxa",
        "taxa_matrix",
        "kg_edges",
        "kg_nodes",
    } <= set(record.product_ids)
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
    assert (
        "Native baseline was not supplied"
        in (output / "influence_extension_summary.md").read_text()
    )
    assert not (output / "native_summary.md").exists()
    with (output / "lineage_fraction_effects.tsv").open() as handle:
        effects = list(csv.DictReader(handle, delimiter="\t"))
    joint = [
        row
        for row in effects
        if row["selection_kind"] == "declared_joint" and row["scope"] == "pooled"
    ]
    assert len(joint) == 6
    # Neither declared family is in this synthetic graph: distinguish no removed
    # annotation mass from a supported removal that happens to leave an effect unchanged.
    assert {row["removal_status"] for row in joint} == {"no_removed_annotation_mass"}
    assert {int(row["removed_annotation_mass"]) for row in joint} == {0}
    assert {int(row["removed_trait_annotation_mass"]) for row in joint} == {0}
    for row in joint:
        if float(row["before_effect"]) != 0:
            assert float(row["effect_reduction_fraction"]) == pytest.approx(0.0)
            assert float(row["effect_reduction_percent"]) == pytest.approx(0.0)
            assert row["effect_reduction_status"] == "defined"
        else:
            assert row["effect_reduction_fraction"] == ""
            assert row["effect_reduction_percent"] == ""
            assert row["effect_reduction_status"] == "zero_before_effect"
    mass_path = output / "lineage_removed_annotation_mass.tsv"
    with mass_path.open() as handle:
        mass = list(csv.DictReader(handle, delimiter="\t"))
    pooled_mass = [
        row
        for row in mass
        if row["scope"] == "pooled" and row["selection_kind"] == "declared_joint"
    ]
    assert len(pooled_mass) == 2 * len(provenance["retained_features"])
    assert {row["trait_key"] for row in pooled_mass} == set(provenance["retained_features"])
    assert {row["direction"] for row in pooled_mass} == {"inc", "dec"}
    assert {int(row["removed_annotation_mass"]) for row in pooled_mass} == {0}
    assert all(int(row["before_annotation_mass"]) > 0 for row in pooled_mass)
    for direction in ("inc", "dec"):
        directional = [row for row in pooled_mass if row["direction"] == direction]
        assert sum(int(row["before_trait_annotation_mass"]) for row in directional) == int(
            directional[0]["before_annotation_mass"]
        )
    assert provenance["outputs"][mass_path.name] == file_sha256(mass_path)
    assert str(mass_path.resolve()) in {item.path for item in record.outputs}
    with pytest.raises(FileExistsError):
        run_extension(manifest, output)
    altered = bundle.path("taxa_matrix")
    altered.write_text(altered.read_text() + "\n")
    with pytest.raises(ValueError, match="changed"):
        run_extension(manifest, tmp_path / "must-not-write")
    assert not (tmp_path / "must-not-write").exists()


def extension_core(root: Path) -> Path:
    """One real fixture core for protocol and verified-Noel integration checks."""
    inputs = source_inputs(root / "inputs")
    core = execute_workflow(
        REPO_ROOT / "workflows/bugsigdb_core.yaml",
        inputs,
        root / "core",
        OPERATIONS,
        parameter_overrides={"source_mode": "fixture", "strain_rollup": "off"},
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    return Path(core.manifest_path)


def test_protocol_rule_tampering_is_rejected_before_output_creation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    core = extension_core(tmp_path)
    changed = json.loads(extension.PROTOCOL.read_text())
    assert changed["version"] == "771-extension-v1"
    changed["top_per_trait_direction_rank_scope"] = 1
    tampered = tmp_path / "same-version-changed-rule.json"
    tampered.write_text(json.dumps(changed))
    monkeypatch.setattr(extension, "PROTOCOL", tampered)

    def unexpected_analysis(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("Analysis started before the protocol pin was validated")

    monkeypatch.setattr(extension, "load_base", unexpected_analysis)
    output = tmp_path / "must-not-write"
    with pytest.raises(ValueError, match="protocol hash differs from the declared rule"):
        run_extension(core, output)
    assert not output.exists()


def test_cli_includes_a_real_verified_same_core_noel_summary(tmp_path: Path) -> None:
    """Exercise immutable arm binding with a minimal fixture producer, no verifier stub."""
    core = extension_core(tmp_path)

    def plan(context: ArmContext) -> list[Command]:
        source = context.products["count_profiles"]
        report = context.output / "native_summary.md"
        return [
            Command(
                "native_taxon_influence",
                (
                    sys.executable,
                    str(REPO_ROOT / "scripts/adhoc/issue658_arm_test_command.py"),
                    "--input",
                    str(source),
                    "--output",
                    str(report),
                    "--payload",
                    "\n# Synthetic native baseline\nVerified fixture binding marker.\n",
                ),
                (source,),
                (report,),
            )
        ]

    produced = execute_arm("noel", core, tmp_path / "noel", plan, figures=False)
    noel_manifest = Path(produced.manifest_path)
    verified = verify_arm(noel_manifest)
    assert verified.steps
    assert verified.steps[0].id == "native_taxon_influence"
    assert verified.steps[0].outputs
    assert len(verified.steps[0].outputs) == 1
    native_path = Path(verified.steps[0].outputs[0].path)
    assert native_path.name == "native_summary.md"
    native_text, native_hash = native_path.read_text(), file_sha256(native_path)
    output = tmp_path / "with-noel"
    assert (
        extension.main(
            [
                "--core-run",
                str(core),
                "--output-dir",
                str(output),
                "--noel-run",
                str(noel_manifest),
            ]
        )
        == 0
    )
    assert native_text in (output / "influence_extension_summary.md").read_text()
    assert not (output / "native_summary.md").exists()
    provenance = json.loads((output / "extension_provenance.json").read_text())
    assert provenance["input_sha256"]["native_summary"] == native_hash
    assert provenance["input_sha256"]["noel_manifest"] == file_sha256(noel_manifest)
    record = verify_consumer_record(output / "core_consumption.json")
    assert {
        "count_profiles",
        "experiments",
        "reported_taxa",
        "taxa_matrix",
        "kg_edges",
        "kg_nodes",
    } <= set(record.product_ids)
    assert file_sha256(native_path) == native_hash
    assert verify_arm(noel_manifest) == verified
