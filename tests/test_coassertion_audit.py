"""Scope and evidence accounting for the independent KG co-assertion census (#846)."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import pytest

from mtp import coassertion_audit as audit
from mtp import taxa_trait_matrix as matrix_api
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import OXYGEN_GUILDS


@pytest.fixture
def inputs(tmp_path: Path) -> tuple[Path, Path]:
    edges = tmp_path / "edges.tsv"
    nodes = tmp_path / "nodes.tsv"
    rows = [
        # Direct broad oxygen co-assertion; duplicate source records count one node.
        ("NCBITaxon:1", "biolink:has_phenotype", "METPO:1000602", "one|two"),
        ("NCBITaxon:1", "biolink:has_phenotype", "METPO:1000602", "one"),
        ("NCBITaxon:1", "biolink:has_phenotype", "METPO:1000603", "one"),
        # A directly co-asserting taxon also has a supporting strain. The off
        # census must retain the taxon while excluding every strain attribution.
        ("kgmicrobe.strain:s", "biolink:subclass_of", "NCBITaxon:1", "bacdive"),
        ("kgmicrobe.strain:s", "biolink:has_phenotype", "METPO:1000602", "bacdive"),
        ("kgmicrobe.strain:s", "biolink:has_phenotype", "METPO:1000603", "bacdive"),
        # Different strains of taxon2 supply different classes.
        ("kgmicrobe.strain:a", "biolink:subclass_of", "NCBITaxon:2", "bacdive"),
        ("kgmicrobe.strain:b", "biolink:subclass_of", "NCBITaxon:2", "bacdive"),
        ("kgmicrobe.strain:a", "biolink:has_phenotype", "METPO:1000602", "bacdive"),
        ("kgmicrobe.strain:b", "biolink:has_phenotype", "METPO:1000603", "bacdive"),
        # One strain co-asserts both, and maps to two distinct taxon nodes.
        ("kgmicrobe.strain:c", "biolink:subclass_of", "NCBITaxon:3", "bacdive"),
        ("kgmicrobe.strain:c", "biolink:subclass_of", "NCBITaxon:4", "bacdive"),
        ("kgmicrobe.strain:c", "biolink:has_phenotype", "METPO:1000615", "bacdive"),
        ("kgmicrobe.strain:c", "biolink:has_phenotype", "METPO:1000614", "bacdive"),
        ("kgmicrobe.strain:c", "biolink:has_phenotype", "METPO:1000605", "bacdive"),
        ("kgmicrobe.strain:c", "biolink:has_phenotype", "METPO:1000606", "bacdive"),
        # Direct and strain assertions combine, without one source node carrying both.
        ("NCBITaxon:5", "biolink:has_phenotype", "METPO:1000606", "gold"),
        ("kgmicrobe.strain:d", "biolink:subclass_of", "NCBITaxon:5", "bacdive"),
        ("kgmicrobe.strain:d", "biolink:has_phenotype", "METPO:1000603", "bacdive"),
        # An unlinked co-asserting strain contributes no taxon node.
        ("kgmicrobe.strain:u", "biolink:has_phenotype", "METPO:1000602", "bacdive"),
        ("kgmicrobe.strain:u", "biolink:has_phenotype", "METPO:1000603", "bacdive"),
        # Neither a negated object nor an ancestor union is a direct positive node.
        ("NCBITaxon:6", "METPO:2000022", "METPO:1000602", "negative"),
        ("NCBITaxon:6", "biolink:has_phenotype", "METPO:1000603", "one"),
        ("NCBITaxon:1", "biolink:subclass_of", "NCBITaxon:9", "taxonomy"),
        ("NCBITaxon:2", "biolink:subclass_of", "NCBITaxon:9", "taxonomy"),
    ]
    with edges.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "subject", "predicate", "object", "primary_knowledge_source"])
        for index, row in enumerate(rows):
            writer.writerow([f"edge{index}", *row])
    classes = {key for pair in audit.PAIRS for key in (*pair.left, *pair.right)}
    with nodes.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "category", "name"])
        for key in sorted(classes):
            writer.writerow([key, "biolink:PhenotypicQuality", f"fixture label {key}"])
        for number in range(1, 10):
            writer.writerow(
                [f"NCBITaxon:{number}", "biolink:OrganismTaxon", f"Fixture species {number}"]
            )
    return edges, nodes


def table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_scope_counts_preserve_sources_and_do_not_propagate_to_ancestors(
    inputs: tuple[Path, Path], tmp_path: Path
) -> None:
    edges, nodes = inputs
    original = (file_sha256(edges), file_sha256(nodes))
    output = tmp_path / "audit"
    audit.run_audit(edges, nodes, output)
    summary = {
        (row["scope"], row["pair_id"]): row for row in table(output / "coassertion_summary.tsv")
    }
    assert len(summary) == len(audit.PAIRS) * 3
    for scope, expected in [("off", 1), ("phenotype", 2), ("full", 2)]:
        assert (
            int(summary[scope, "aerobic_class_and_anaerobic_class"]["coasserted_taxon_nodes"])
            == expected
        )
    assert summary["phenotype", "mesophilic_and_psychrophilic"]["coasserted_taxon_nodes"] == "2"
    assert summary["phenotype", "obligately_aerobic_and_anaerobic"]["coasserted_taxon_nodes"] == "1"
    assert (
        summary["phenotype", "facultatively_anaerobic_and_obligately_aerobic"][
            "coasserted_taxon_nodes"
        ]
        == "2"
    )
    records = table(output / "coasserted_taxon_nodes.tsv")
    assert not {"NCBITaxon:6", "NCBITaxon:9"} & {row["taxon_id"] for row in records}
    by_node = {(row["scope"], row["pair_id"], row["taxon_id"]): row for row in records}
    direct = by_node["off", "aerobic_class_and_anaerobic_class", "NCBITaxon:1"]
    assert direct["direct_taxon_has_both"] == "True"
    assert direct["same_strain_has_both"] == "False"
    assert direct["requires_union_across_subjects"] == "False"
    for field in ("left_strain_ids", "right_strain_ids", "same_strain_ids"):
        assert json.loads(direct[field]) == []
    supported_direct = by_node["phenotype", "aerobic_class_and_anaerobic_class", "NCBITaxon:1"]
    assert supported_direct["direct_taxon_has_both"] == "True"
    assert supported_direct["same_strain_has_both"] == "True"
    assert supported_direct["requires_union_across_subjects"] == "False"
    mixed = by_node["phenotype", "aerobic_class_and_anaerobic_class", "NCBITaxon:2"]
    assert mixed["requires_union_across_subjects"] == "True"
    assert mixed["same_strain_has_both"] == "False"
    same = by_node["phenotype", "mesophilic_and_psychrophilic", "NCBITaxon:3"]
    assert same["same_strain_has_both"] == "True"
    assert same["direct_taxon_has_both"] == "False"
    assert same["requires_union_across_subjects"] == "False"
    strains = table(output / "coasserted_strain_nodes.tsv")
    assert not [row for row in strains if row["scope"] == "off"]
    unlinked = next(
        row
        for row in strains
        if row["scope"] == "phenotype" and row["strain_id"] == "kgmicrobe.strain:u"
    )
    assert unlinked["mapped_taxon_count"] == "0"
    evidence = table(output / "coassertion_evidence.tsv")
    off_evidence = [row for row in evidence if row["scope"] == "off"]
    assert off_evidence
    assert all(row["origin"] == "direct_taxon" for row in off_evidence)
    assert all(row["assertion_subject"].startswith("NCBITaxon:") for row in off_evidence)
    assert {"one", "two"} <= {row["primary_source"] for row in evidence}
    assert any(row["assertion_subject"] == "kgmicrobe.strain:b" for row in evidence)
    unmapped_evidence = [
        row for row in evidence if row["assertion_subject"] == "kgmicrobe.strain:u"
    ]
    assert unmapped_evidence
    assert all(
        row["taxon_id"] == "" and row["origin"] == "unmapped_strain" for row in unmapped_evidence
    )
    assert all(row["edge_id"] for row in evidence)
    assert original == (file_sha256(edges), file_sha256(nodes))
    manifest = json.loads((output / "provenance.json").read_text())
    assert manifest["inputs"]["kg_edges"]["sha256"] == original[0]
    assert manifest["analysis_kind"] == "kg_annotation_census_not_a_core_or_arm_rerun"
    assert set(manifest["outputs"]) == {path.name for path in output.iterdir()} - {
        "provenance.json"
    }
    for name, digest in manifest["outputs"].items():
        assert file_sha256(output / name) == digest
    report = (output / "README.md").read_text()
    assert "counts of annotations" in report
    assert "not a measurement" in report
    assert "aerobic_class_and_anaerobic_class | phenotype | 2" in report


def test_canonical_guild_definition_and_shared_scope_policy(
    inputs: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    guild = next(pair for pair in audit.PAIRS if pair.id == "aerobe_guild_and_anaerobe_guild")
    assert guild.left == OXYGEN_GUILDS["aerobe"]
    assert guild.right == OXYGEN_GUILDS["anaerobe"]
    original = matrix_api.scoped_taxon_traits
    called = []

    def observed(
        taxa: dict[str, set[str]],
        links: dict[str, set[str]],
        strains: dict[str, set[str]],
        scope: matrix_api.StrainRollup,
    ) -> dict[str, set[str]]:
        called.append(scope)
        return original(taxa, links, strains, scope)

    monkeypatch.setattr(matrix_api, "scoped_taxon_traits", observed)
    audit.run_audit(*inputs, tmp_path / "audit")
    assert called == ["off", "phenotype", "full"]


@pytest.mark.parametrize("scope", matrix_api.STRAIN_ROLLUPS)
def test_target_carriers_match_canonical_edge_reader_before_ancestor_closure(
    inputs: tuple[Path, Path], scope: matrix_api.StrainRollup
) -> None:
    edges, _ = inputs
    with edges.open("a", encoding="utf-8") as handle:
        handle.write(
            "off_target\tNCBITaxon:7\tbiolink:has_phenotype\tMETPO:1000668\tfixture\n"
            "activity\tNCBITaxon:7\tMETPO:2000001\tCHEBI:123\tfixture\n"
            "non_taxon\texample:organism\tbiolink:has_phenotype\tMETPO:1000602\tfixture\n"
            "other_predicate\tNCBITaxon:8\tbiolink:related_to\tMETPO:1000602\tfixture\n"
        )
    graph = audit.read_graph(edges)
    observed = matrix_api.scoped_taxon_traits(
        graph.taxon_traits, graph.strain_taxa, graph.strain_traits, scope
    )
    _, canonical = matrix_api.read_kg_edges(edges, strain_rollup=scope)
    expected = {
        node: target_keys
        for node, keys in canonical.items()
        if (target_keys := keys & audit.TARGET_KEYS)
    }
    assert observed == expected


def test_no_target_labels_yields_explicit_zero_rows(
    inputs: tuple[Path, Path], tmp_path: Path
) -> None:
    edges, nodes = inputs
    edges.write_text(
        "id\tsubject\tpredicate\tobject\tprimary_knowledge_source\n"
        "only\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1000668\tfixture\n"
    )
    output = tmp_path / "empty"
    audit.run_audit(edges, nodes, output)
    rows = table(output / "coassertion_summary.tsv")
    assert len(rows) == 3 * len(audit.PAIRS)
    assert all(row["coasserted_taxon_nodes"] == "0" for row in rows)
    assert not table(output / "coasserted_taxon_nodes.tsv")


def test_rejects_missing_class_metadata_and_output_reuse(
    inputs: tuple[Path, Path], tmp_path: Path
) -> None:
    edges, nodes = inputs
    output = tmp_path / "audit"
    audit.run_audit(edges, nodes, output)
    with pytest.raises(ValueError, match="new or empty"):
        audit.run_audit(edges, nodes, output)
    nodes.write_text("id\tcategory\tname\n")
    with pytest.raises(ValueError, match="class metadata"):
        audit.run_audit(edges, nodes, tmp_path / "missing")
    assert not (tmp_path / "missing").exists()


def test_actual_recipe_reads_known_cr_serialization_and_keeps_source_rows(
    inputs: tuple[Path, Path], tmp_path: Path
) -> None:
    edges, nodes = inputs
    records = table(edges)
    header = (
        "id\tsubject\tpredicate\tobject\trelation\tagent_type\tagent_type\r"
        "\tknowledge_level\tknowledge_source\tmeta\tprimary_knowledge_source\n"
    )
    lines = []
    for index, row in enumerate(records):
        agent = "manual_validation_of_automated_agent\r" if index == 4 else "manual_agent"
        lines.append(
            "\t".join(
                [
                    row["id"],
                    row["subject"],
                    row["predicate"],
                    row["object"],
                    "",
                    "automated_agent",
                    agent,
                    "knowledge_assertion",
                    "",
                    "",
                    row["primary_knowledge_source"],
                ]
            )
            + "\n"
        )
    edges.write_bytes((header + "".join(lines)).encode())
    output = tmp_path / "CLI audit with spaces"
    completed = subprocess.run(
        [
            "just",
            "audit-coassertions",
            "--edges",
            str(edges),
            "--nodes",
            str(nodes),
            "--output-dir",
            str(output),
        ],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    facts = json.loads((output / "provenance.json").read_text())["kg_edge_serialization"]
    assert facts["records"] == len(records)
    assert facts["header_continuations_reconstructed"] == 1
    assert facts["agent_field_continuations_reconstructed"] == 1
    assert facts["excluded_records"] == 0
    assert table(output / "coasserted_taxon_nodes.tsv")


def test_inputs_changed_during_audit_fail_before_publication(
    inputs: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = audit.build_tables

    def changed(
        graph: audit.Graph, names: dict[str, str]
    ) -> dict[str, tuple[tuple[str, ...], list[dict[str, object]]]]:
        result = original(graph, names)
        inputs[0].write_bytes(inputs[0].read_bytes() + b"changed\n")
        return result

    monkeypatch.setattr(audit, "build_tables", changed)
    with pytest.raises(ValueError, match="changed during execution"):
        audit.run_audit(*inputs, tmp_path / "changed")
    assert not (tmp_path / "changed").exists()
