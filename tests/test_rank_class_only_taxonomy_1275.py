"""Rank vocabulary cannot establish organism ancestry (#1275)."""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/create_taxa_coverage_histogram.py"
FIXTURE = ROOT / "tests/fixtures/coverage/taxonomy_rank_class.tsv"
SPEC = importlib.util.spec_from_file_location("coverage_rank_class_1275", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load coverage taxonomy script")
coverage = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = coverage
SPEC.loader.exec_module(coverage)


@pytest.mark.parametrize("case", ["species", "multiple", "fixture"])
def test_rank_class_only_is_rejected(tmp_path: Path, case: str) -> None:
    path = tmp_path / "rank_only.tsv"
    if case == "fixture":
        path.write_bytes(FIXTURE.read_bytes())
    else:
        ranks = ["species"] if case == "species" else ["species", "genus", "biotype"]
        path.write_text(
            "id\tsubject\tpredicate\tobject\n"
            + "".join(
                f"r{index}\tNCBITaxon:{rank}\tbiolink:subclass_of\tOBO:NCBITaxon#_taxonomic_rank\n"
                for index, rank in enumerate(ranks)
            )
        )
    error = None
    try:
        coverage.CoverageTaxonomy(path)
    except ValueError as exc:
        error = str(exc)
    assert error == "Taxonomy input has no NCBITaxon subclass_of edges"


def test_rank_classes_are_counted_and_excluded_with_organism_edge(tmp_path: Path) -> None:
    path = tmp_path / "mixed.tsv"
    # Keep the recorded fixture's split header and full edge encoding intact.
    organism = ["organism", "NCBITaxon:100", "biolink:subclass_of", "NCBITaxon:2"]
    path.write_bytes(FIXTURE.read_bytes() + ("\t".join(organism + [""] * 7) + "\n").encode())
    taxonomy = None
    error = None
    try:
        taxonomy = coverage.CoverageTaxonomy(path)
    except ValueError as exc:
        error = str(exc)
    assert error is None
    assert taxonomy is not None
    assert taxonomy.audit["taxonomy_rank_class_edge_records_skipped"] == 1
    assert taxonomy.audit["taxonomy_edge_records"] == 1
    assert taxonomy.audit["unique_taxonomy_edges"] == 1
    assert taxonomy.audit["taxonomy_node_count"] == 2
    assert taxonomy.audit["excluded_records"] == 0
    assert taxonomy.nodes == {"NCBITaxon:100", "NCBITaxon:2"}
    assert taxonomy.classify("NCBITaxon:100") == ("bacteria", "resolved_root", ("NCBITaxon:2",))
    assert taxonomy.classify("NCBITaxon:biotype").status == "non_ncbi_identifier"
