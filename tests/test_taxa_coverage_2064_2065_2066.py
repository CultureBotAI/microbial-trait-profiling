"""Regression tests for coverage rank-root handling and audit wording (#2064, #2065, #2066)."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "create_taxa_coverage_histogram.py"
WORKFLOW_DOC_PATH = REPO_ROOT / "docs" / "workflow.md"

_SPEC = importlib.util.spec_from_file_location("taxa_coverage_histogram_2064", SCRIPT_PATH)
assert _SPEC is not None
assert _SPEC.loader is not None
coverage = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(coverage)


def _write_minimal_coverage_inputs(base_dir: Path) -> dict[str, Path]:
    experiments = base_dir / "experiments.tsv"
    profiles = base_dir / "profiles.tsv"
    matrix = base_dir / "matrix.tsv"
    taxa = base_dir / "taxa.tsv"

    pd.DataFrame(
        [("s1", "e1", "Homo sapiens", "16S")],
        columns=["study_id", "experiment_id", "host_species", "sequencing_type"],
    ).to_csv(experiments, sep="\t", index=False)
    pd.DataFrame(
        [("s1", "e1", "inc", 2, 1)],
        columns=["study_id", "experiment_id", "inc_dec", "taxa_count", "METPO:1"],
    ).to_csv(profiles, sep="\t", index=False)
    pd.DataFrame(
        [
            ("NCBITaxon:100", "Escherichia coli", "species", 1),
            ("NCBITaxon:101", "uncultured bacterium", "species", 0),
        ],
        columns=["taxon_id", "taxon_name", "taxonomic_rank", "METPO:1"],
    ).to_csv(matrix, sep="\t", index=False)
    pd.DataFrame(
        [
            ("s1", "e1", "sig1", "increased", "NCBITaxon:100", "Escherichia coli", "species"),
            ("s1", "e1", "sig1", "increased", "NCBITaxon:101", "uncultured bacterium", "species"),
        ],
        columns=[
            "study_id",
            "experiment_id",
            "signature_id",
            "direction",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
        ],
    ).to_csv(taxa, sep="\t", index=False)
    return {
        "experiments": experiments,
        "profiles": profiles,
        "matrix": matrix,
        "taxa": taxa,
    }


def test_obo_form_rank_classes_under_taxonomic_rank_skip_without_counting_as_ncbitaxon_rank(
    tmp_path: Path,
) -> None:
    """Issue #2064: OBO-form rank-class subjects take the non-NCBITaxon skip."""
    edges = tmp_path / "taxonomy.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\n"
        "r0\tNCBITaxon:species\tbiolink:subclass_of\tOBO:NCBITaxon#_taxonomic_rank\n"
        "r1\tOBO:NCBITaxon#_species_group\tbiolink:subclass_of\tOBO:NCBITaxon#_taxonomic_rank\n"
        "r2\tOBO:NCBITaxon#_species_subgroup\tbiolink:subclass_of\tOBO:NCBITaxon#_taxonomic_rank\n"
        "e0\tNCBITaxon:91626\tbiolink:subclass_of\tNCBITaxon:2759\n",
        encoding="utf-8",
    )
    try:
        taxonomy = coverage.CoverageTaxonomy(edges)
    except ValueError as exc:
        raise AssertionError(
            f"Expected OBO-form rank-class subjects to skip without ValueError, got: {exc}"
        ) from exc

    assert taxonomy.audit["excluded_records"] == 0
    assert taxonomy.audit["taxonomy_rank_class_edge_records_skipped"] == 1
    assert taxonomy.audit["taxonomy_edge_records"] == 1
    assert taxonomy.audit["unique_taxonomy_edges"] == 1
    assert taxonomy.audit["taxonomy_node_count"] == 2
    assert taxonomy.nodes == {"NCBITaxon:91626", "NCBITaxon:2759"}


def test_excluded_records_note_pins_always_zero_and_ncbitaxon_name_rank_scope(
    tmp_path: Path,
) -> None:
    """Issues #2064 and #2066: note pins (always 0) and NCBITaxon:<name> vs OBO-form scope."""
    edges = tmp_path / "taxonomy.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\n"
        "r0\tNCBITaxon:species\tbiolink:subclass_of\tOBO:NCBITaxon#_taxonomic_rank\n"
        "e0\tNCBITaxon:91626\tbiolink:subclass_of\tNCBITaxon:2759\n",
        encoding="utf-8",
    )
    taxonomy = coverage.CoverageTaxonomy(edges)
    note = str(taxonomy.audit["excluded_records_note"])
    assert (
        "excluded_records counts malformed KG edge records dropped by the strict "
        "edge parser (always 0)"
    ) in note
    assert (
        "taxonomy_rank_class_edge_records_skipped counts well-formed NCBITaxon:<name> "
        "ontology rank-class subclass_of records skipped during organism-ancestry "
        "graph projection, while OBO-form rank-class subjects take the non-NCBITaxon "
        "skip (#1008, #1274)."
    ) in note


def test_non_rank_ncbitaxon_subject_under_taxonomic_rank_raises_distinct_error(
    tmp_path: Path,
) -> None:
    """Issue #2065: non-rank NCBITaxon subject under taxonomic_rank reports distinct error."""
    for bad_subject in ("NCBITaxon:10", "NCBITaxon:0"):
        edges = tmp_path / f"bad_{bad_subject.replace(':', '_')}.tsv"
        edges.write_text(
            "id\tsubject\tpredicate\tobject\n"
            "e0\tNCBITaxon:91626\tbiolink:subclass_of\tNCBITaxon:2759\n"
            f"e1\t{bad_subject}\tbiolink:subclass_of\tOBO:NCBITaxon#_taxonomic_rank\n",
            encoding="utf-8",
        )
        raised_message: str | None = None
        try:
            coverage.CoverageTaxonomy(edges)
        except ValueError as exc:
            raised_message = str(exc)

        assert raised_message == (
            f"Rank-class edge e1 has a non-rank NCBITaxon subject: {bad_subject}"
        )


def test_readme_pins_always_zero_strict_edge_parser_wording(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Issue #2066: generated README explains excluded_records is always 0."""
    if not (coverage.REPO_ROOT / ".git").exists():
        import mtp

        monkeypatch.setattr(coverage, "REPO_ROOT", Path(mtp.__file__).resolve().parents[2])
    inputs = _write_minimal_coverage_inputs(tmp_path)
    taxonomy_edges = tmp_path / "taxonomy.tsv"
    taxonomy_edges.write_text(
        "id\tsubject\tpredicate\tobject\n"
        "r0\tNCBITaxon:species\tbiolink:subclass_of\tOBO:NCBITaxon#_taxonomic_rank\n"
        "e0\tNCBITaxon:100\tbiolink:subclass_of\tNCBITaxon:2\n"
        "e1\tNCBITaxon:101\tbiolink:subclass_of\tNCBITaxon:2\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "out"
    coverage.run_coverage(
        output_dir=output_dir,
        tables_only=True,
        strain_rollup="off",
        experiments=inputs["experiments"],
        profiles=inputs["profiles"],
        matrix=inputs["matrix"],
        taxa=inputs["taxa"],
        taxonomy_edges=taxonomy_edges,
    )
    provenance = json.loads((output_dir / "provenance.json").read_text(encoding="utf-8"))
    audit = provenance["taxonomy_diagnosis"]["parser_audit"]
    assert audit["excluded_records"] == 0
    assert audit["taxonomy_rank_class_edge_records_skipped"] == 1

    readme = (output_dir / "README.md").read_text(encoding="utf-8")
    assert (
        "excluded_records is always 0: the strict KG edge parser raises on a "
        "malformed record instead of dropping it)."
    ) in readme
    assert "excluded_records: 0 confirms the strict KG edge parser dropped no record" not in readme


def test_workflow_doc_853_pins_non_rank_ncbitaxon_subjects_under_taxonomic_rank() -> None:
    """Issue #2065: docs/workflow.md documents non-rank NCBITaxon subjects under taxonomic_rank."""
    workflow_text = WORKFLOW_DOC_PATH.read_text(encoding="utf-8")
    assert (
        "Malformed NCBITaxon identifiers in organism ancestry edges and non-rank "
        "NCBITaxon subjects under OBO:NCBITaxon#_taxonomic_rank still fail."
    ) in workflow_text
