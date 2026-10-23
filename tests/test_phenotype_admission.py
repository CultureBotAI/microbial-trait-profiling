"""Future-KG admission and explicit-negation regressions for issue #848."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from mtp import taxa_trait_matrix as matrix_api
from mtp.contribution_inputs import load_graph
from mtp.study_trait_profiles import collect_profile_taxa
from mtp.trait_terms import unlisted_negations
from mtp.workflow.matrix_contracts import validate_matrix_sources
from mtp.workflow.operations import trace_profile_feature
from test_taxa_trait_matrix import ALL_EDGES, _write_strain_fixture


def build(directory: Path, scope: matrix_api.StrainRollup) -> list[str]:
    output = directory / "matrix.tsv"
    matrix_api.create_taxa_trait_matrix(
        taxa_path=directory / "bugsigdb_taxa.tsv",
        kg_dir=directory,
        output_path=output,
        column_metadata_path=directory / "columns.tsv",
        edges_path=directory / "edges.tsv",
        nodes_path=directory / "nodes.tsv",
        strain_rollup=scope,
    )
    return output.read_text().splitlines()[0].split("\t")[3:]


@pytest.mark.parametrize("subject", ["kgmicrobe.strain:1", "NCBITaxon:101"])
def test_only_has_phenotype_can_admit_a_metpo_class(tmp_path: Path, subject: str) -> None:
    edge = f"future\t{subject}\tMETPO:2000103\tMETPO:1000682"
    _, _, edges = _write_strain_fixture(tmp_path, [*ALL_EDGES, edge])
    raw = matrix_api.read_kg_trait_edges(edges, strains=True)
    keys = (
        raw.strain_traits[subject] if subject.startswith("kgmicrobe") else raw.taxon_traits[subject]
    )
    assert "METPO:2000103|METPO:1000682" in keys
    assert "METPO:1000682" not in keys
    assert not any("METPO:1000682" in key for key in build(tmp_path, "phenotype"))


@pytest.mark.parametrize("scope", matrix_api.STRAIN_ROLLUPS)
def test_handling_classes_are_excluded_from_direct_phenotype_edges(
    tmp_path: Path, scope: matrix_api.StrainRollup
) -> None:
    edges = [
        f"handling-{number}\tNCBITaxon:101\tbiolink:has_phenotype\t{trait}"
        for number, trait in enumerate(matrix_api.PHENOTYPE_EXCLUDED_TRAITS)
    ]
    _write_strain_fixture(tmp_path, [*ALL_EDGES, *edges])
    present = set(build(tmp_path, scope)) & matrix_api.PHENOTYPE_EXCLUDED_TRAITS.keys()
    assert present == (set() if scope == "phenotype" else set(matrix_api.PHENOTYPE_EXCLUDED_TRAITS))


@pytest.mark.parametrize("scope", matrix_api.STRAIN_ROLLUPS)
@pytest.mark.parametrize("field", ["negated", "is_negated", "meta"])
@pytest.mark.parametrize("subject", ["NCBITaxon:101", "kgmicrobe.strain:1"])
def test_positive_predicate_with_negated_assertion_requires_review_in_every_scope(
    tmp_path: Path, scope: matrix_api.StrainRollup, field: str, subject: str
) -> None:
    edge = f"conflict-edge\t{subject}\tbiolink:has_phenotype\tMETPO:1000682"
    _, _, edges = _write_strain_fixture(tmp_path, [*ALL_EDGES, edge])
    lines = edges.read_text().splitlines()
    value = json.dumps({"negated": True}) if field == "meta" else "true"
    edges.write_text(
        lines[0]
        + "\t"
        + field
        + "\n"
        + "\n".join(
            line + "\t" + (value if line.startswith("conflict-edge\t") else "")
            for line in lines[1:]
        )
        + "\n"
    )
    with pytest.raises(ValueError, match=r"conflict-edge.*negat.*review"):
        build(tmp_path, scope)
    assert not (tmp_path / "matrix.tsv").exists()


@pytest.mark.parametrize("field", ["synonym", "description", "negated"])
def test_unlisted_relation_negation_is_found_outside_its_label(tmp_path: Path, field: str) -> None:
    path = tmp_path / "nodes.tsv"
    evidence = {
        "synonym": "does not swim",
        "description": "A relation asserting the absence of swimming capability.",
        "negated": "true",
    }
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["id", "category", "name", field], delimiter="\t"
        )
        writer.writeheader()
        writer.writerow(
            {
                "id": "METPO:2999999",
                "category": "biolink:OntologyClass",
                "name": "novel relation",
                field: evidence[field],
            }
        )
    assert unlisted_negations(path) == ["METPO:2999999"]


@pytest.mark.parametrize("scope", matrix_api.STRAIN_ROLLUPS)
def test_unlisted_negative_predicate_stops_matrix_publication(
    tmp_path: Path, scope: matrix_api.StrainRollup
) -> None:
    edge = "unlisted-edge\tNCBITaxon:101\tMETPO:2999999\tCHEBI:17234"
    _, nodes, _ = _write_strain_fixture(tmp_path, [*ALL_EDGES, edge])
    lines = nodes.read_text().splitlines()
    nodes.write_text(
        lines[0]
        + "\tsynonym\n"
        + "\n".join(line + "\t" for line in lines[1:])
        + "\nMETPO:2999999\tbiolink:OntologyClass\tnovel relation\tdoes not ferment\n"
    )
    with pytest.raises(ValueError, match=r"METPO:2999999.*negat.*review"):
        build(tmp_path, scope)
    assert not (tmp_path / "matrix.tsv").exists()


def test_negative_words_in_phenotype_descriptions_are_not_relation_negation(tmp_path: Path) -> None:
    nodes = tmp_path / "nodes.tsv"
    nodes.write_text(
        "id\tcategory\tname\tdescription\tsynonym\n"
        "METPO:1000699\tbiolink:PhenotypicQuality\tgram negative\t"
        "A negative Gram staining result.\tGram-negative\n"
        "METPO:1999999\tbiolink:PhenotypicQuality\tnonmotile\t"
        "A phenotype describing absence of motility.\tnon-motile\n"
    )
    assert unlisted_negations(nodes) == []


@pytest.mark.parametrize("scope", matrix_api.STRAIN_ROLLUPS)
def test_reviewed_negative_predicate_never_becomes_its_metpo_object(
    tmp_path: Path, scope: matrix_api.StrainRollup
) -> None:
    _write_strain_fixture(
        tmp_path, [*ALL_EDGES, "absent\tNCBITaxon:101\tMETPO:2000037\tMETPO:1000682"]
    )
    columns = build(tmp_path, scope)
    assert "METPO:1000682" not in columns
    assert ("METPO:2000037|METPO:1000682" in columns) == (scope == "off")


@pytest.mark.parametrize("scope", matrix_api.STRAIN_ROLLUPS)
def test_independent_readers_share_polarity_guard_even_for_excluded_strains(
    tmp_path: Path, scope: matrix_api.StrainRollup
) -> None:
    _, nodes, edges = _write_strain_fixture(tmp_path, ALL_EDGES)
    columns = set(build(tmp_path, scope))
    lines = edges.read_text().splitlines()
    edges.write_text(
        lines[0]
        + "\tmeta\n"
        + "\n".join(line + "\t" for line in lines[1:])
        + "\nconflict\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1000603\t"
        + '{"is_negated": true}\n'
    )
    with pytest.raises(ValueError, match=r"conflict.*negat.*review"):
        validate_matrix_sources(
            tmp_path / "matrix.tsv",
            tmp_path / "bugsigdb_taxa.tsv",
            edges,
            nodes,
            tmp_path / "columns.tsv",
            strain_rollup=scope,
        )
    with pytest.raises(ValueError, match=r"conflict.*negat.*review"):
        load_graph(edges, nodes, columns, columns, strain_rollup=scope)
    studies = tmp_path / "studies"
    studies.mkdir()
    (studies / "bsdb_1.yaml").write_text(
        "id: bsdb:1\nexperiments:\n- id: bsdb:1/1/1\n"
        "  host_species: Homo sapiens\n  sequencing_type: 16S\n  signatures:\n"
        "  - direction: increased\n    taxa:\n    - id: NCBITaxon:101\n"
    )
    profiles, _ = collect_profile_taxa(studies)
    with pytest.raises(ValueError, match=r"conflict.*negat.*review"):
        trace_profile_feature(
            studies,
            tmp_path / "matrix.tsv",
            edges,
            nodes,
            next(iter(profiles)),
            "METPO:1000603",
            strain_rollup=scope,
        )


def test_bare_cr_header_and_agent_continuation_preserve_negation_metadata(tmp_path: Path) -> None:
    _, _, edges = _write_strain_fixture(tmp_path, ALL_EDGES)
    edges.write_bytes(
        b"id\tsubject\tpredicate\tobject\trelation\tagent_type\tagent_type\r"
        b"\tknowledge_level\tknowledge_source\tmeta\tprimary_knowledge_source\n"
        b"conflict-cr\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000603\tR\t\t"
        b"manual_validation_of_automated_agent\r\tobservation\tKG\t"
        b'{"negated": true}\tfixture\n'
    )
    with pytest.raises(ValueError, match=r"conflict-cr.*negat.*review"):
        build(tmp_path, "phenotype")


def test_positive_relation_definition_caveat_is_not_assertion_negation() -> None:
    """Pinned grows-in wording qualifies a positive relation (#915)."""
    term = {
        "name": "grows in",
        "description": (
            "A relation between a microbe and a growth medium indicating that the microbe "
            "has been demonstrated to grow in that medium under some set of conditions. "
            "This relation represents a generalization from experimental observations "
            "reported in the literature or databases; it does not imply that the microbe "
            "will grow in the medium under all conditions."
        ),
    }
    assert not matrix_api.negative_relation_metadata(term)
    assert (
        matrix_api.assertion_trait_key(
            {"id": "positive", "predicate": "METPO:2000517", "object": "mediadive.medium:1"},
            {"METPO:2000517": term},
        )
        == "METPO:2000517|mediadive.medium:1"
    )


def test_phenotype_absence_label_is_not_a_negative_predicate(tmp_path: Path) -> None:
    """An affirmed absence phenotype is not assertion negation (#914)."""
    nodes = tmp_path / "nodes.tsv"
    nodes.write_text(
        "id\tcategory\tname\tsynonym\n"
        "METPO:1999998\tbiolink:PhenotypicQuality\tlacks flagella\tdoes not swim\n"
    )
    assert unlisted_negations(nodes) == []
    # If this term actually serves as a predicate, its category must not bypass review.
    edges = tmp_path / "edges.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\n"
        "used-as-relation\tNCBITaxon:1\tMETPO:1999998\tMETPO:1000603\n"
    )
    assert unlisted_negations(nodes, edges) == ["METPO:1999998"]
    terms = matrix_api.read_metpo_term_metadata(nodes)
    with pytest.raises(ValueError, match=r"METPO:1999998.*negat.*review"):
        matrix_api.assertion_trait_key(
            {"id": "used-as-relation", "predicate": "METPO:1999998", "object": "METPO:1000603"},
            terms,
        )


def test_unreviewed_negative_predicate_outside_metpo_is_checked(tmp_path: Path) -> None:
    """A METPO class object is consumed even when its predicate uses another namespace (#916)."""
    _, nodes, _ = _write_strain_fixture(
        tmp_path,
        [*ALL_EDGES, "outside\tNCBITaxon:101\texample:negative_relation\tMETPO:1000603"],
    )
    with nodes.open("a") as handle:
        handle.write("example:negative_relation\tbiolink:OntologyClass\tdoes not exhibit\n")
    with pytest.raises(ValueError, match=r"example:negative_relation.*negat.*review"):
        build(tmp_path, "full")


def test_unrelated_node_metadata_is_not_parsed_as_assertion_polarity(tmp_path: Path) -> None:
    nodes, edges = tmp_path / "nodes.tsv", tmp_path / "edges.tsv"
    nodes.write_text(
        "id\tname\tmeta\nNCBITaxon:1\tSpecies one\tfree-form organism metadata\n"
        "biolink:has_phenotype\thas phenotype\t{}\n"
    )
    edges.write_text(
        "id\tsubject\tpredicate\tobject\n"
        "positive\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1000603\n"
    )
    terms = matrix_api.read_assertion_term_metadata(nodes, edges)
    assert set(terms) == {"biolink:has_phenotype"}
    assert matrix_api.read_kg_trait_edges(edges, nodes_path=nodes).taxon_traits == {
        "NCBITaxon:1": {"METPO:1000603"}
    }
