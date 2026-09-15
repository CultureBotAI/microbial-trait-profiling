"""The pinned KG's known multiline records must not silently lose traits."""

from pathlib import Path

import pytest

from bsdbng.contribution_inputs import Graph, count, load_graph, rows


def test_pinned_continuations_preserve_distinct_assertions(tmp_path: Path) -> None:
    edges, nodes = tmp_path / "edges.tsv", tmp_path / "nodes.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\trelation\tagent_type\tagent_type\n"
        "\tknowledge_level\tknowledge_source\tmeta\tprimary_knowledge_source\n"
        "e1\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1000608\tR\tmanual_agent\t\t"
        "observation\tGraph\t\tpublic:1\n"
        "e2\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1000608\tR\t\t"
        "manual_validation_of_automated_agent\n\tlogical_entailment\tGraph\t\tpublic:2\n",
        encoding="utf-8",
    )
    nodes.write_text("id\tname\nNCBITaxon:1\tSpecies one\n", encoding="utf-8")
    graph = load_graph(edges, nodes, {"METPO:1000608"}, {"METPO:1000608"})
    assert graph.traits["NCBITaxon:1"] == {"METPO:1000608"}
    assert len(graph.assertions[("NCBITaxon:1", "METPO:1000608")]) == 2
    assert graph.header_repaired
    assert graph.continued_records == 1


@pytest.mark.parametrize("invalid", ["", "-1", "1.5", "NaN", " 1", "\u0661"])
def test_invalid_counts_are_errors(invalid: str) -> None:
    with pytest.raises(ValueError, match="nonnegative integer"):
        count(invalid)


def test_duplicate_headers_and_cycles_are_errors(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.tsv"
    path.write_text("trait\ttrait\n1\t2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate header"):
        list(rows(path))
    graph = Graph({"a": {"b"}, "b": {"a"}}, {}, {}, {}, False)
    with pytest.raises(ValueError, match="cycle"):
        graph.ancestors("a")


def test_unknown_malformed_record_is_not_skipped(tmp_path: Path) -> None:
    edges, nodes = tmp_path / "edges.tsv", tmp_path / "nodes.tsv"
    edges.write_text("id\tsubject\tpredicate\tobject\ne1\tNCBITaxon:1\n", encoding="utf-8")
    nodes.write_text("id\tname\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Malformed"):
        load_graph(edges, nodes, set(), set())
