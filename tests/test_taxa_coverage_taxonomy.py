"""Lineage-backed coverage diagnosis and input binding regressions (#957)."""

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from test_taxa_coverage_histogram import coverage


def write_edges(path: Path, pairs: list[tuple[str, str]]) -> Path:
    path.write_text(
        "id\tsubject\tpredicate\tobject\n"
        + "".join(
            f"e{i}\tNCBITaxon:{child}\tbiolink:subclass_of\tNCBITaxon:{parent}\n"
            for i, (child, parent) in enumerate(pairs)
        )
    )
    return path


@pytest.fixture
def taxonomic_inputs(tmp_path: Path) -> dict[str, Path]:
    """Four reported zero-annotation taxa; every occurrence is retained."""
    names = {
        "NCBITaxon:91626": "Mucor ambiguus",
        "NCBITaxon:4896": "Schizosaccharomyces pombe",
        "NCBITaxon:155900": "uncultured organism",
        "NCBITaxon:32644": "unidentified",
    }
    frames = {
        "matrix": pd.DataFrame(
            [(taxon, name, "species", 0) for taxon, name in names.items()],
            columns=["taxon_id", "taxon_name", "taxonomic_rank", "METPO:1"],
        ),
        "taxa": pd.DataFrame(
            [
                ("s1", "e1", "sig1", "increased", taxon, name, "species")
                for taxon, name in names.items()
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
        ),
        "profiles": pd.DataFrame(
            [("s1", "e1", "inc", 4, 0)],
            columns=["study_id", "experiment_id", "inc_dec", "taxa_count", "METPO:1"],
        ),
        "experiments": pd.DataFrame(
            [("s1", "e1", "Homo sapiens", "16S")],
            columns=["study_id", "experiment_id", "host_species", "sequencing_type"],
        ),
    }
    paths = {}
    for key, frame in frames.items():
        paths[key] = tmp_path / f"{key}.tsv"
        frame.to_csv(paths[key], sep="\t", index=False)
    return paths


@pytest.mark.parametrize(
    ("taxon", "expected"),
    [
        ("NCBITaxon:91626", "Non-Bacterial / Fungal Eukaryote"),
        ("NCBITaxon:4896", "Non-Bacterial / Fungal Eukaryote"),
        ("NCBITaxon:155900", "Unresolved Taxonomic Label"),
        ("NCBITaxon:32644", "Unresolved Taxonomic Label"),
    ],
)
def test_four_reviewed_taxa_use_lineage_or_unresolved_label(
    taxonomic_inputs: dict[str, Path],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    taxon: str,
    expected: str,
) -> None:
    edges = write_edges(
        tmp_path / "taxonomy.tsv",
        [
            ("91626", "4751"),
            ("4896", "4751"),
            ("4751", "2759"),
            ("2759", "1"),
            ("155900", "1"),
            ("32644", "1"),
        ],
    )
    # The old API accepts these inputs but misdiagnoses all four rows. The new
    # default is redirected to a fixture; neither version opens a production KG.
    monkeypatch.setattr(coverage, "TAXONOMY_EDGES_PATH", edges, raising=False)
    _, _, _, _, _, zero = coverage.load_and_process_data(**taxonomic_inputs)
    row = zero.set_index("taxon_id").loc[taxon]
    assert row.diagnosis_category == expected


@pytest.fixture
def taxonomy(tmp_path: Path) -> Any:
    return coverage.CoverageTaxonomy(
        write_edges(
            tmp_path / "domains.tsv",
            [
                ("10", "2"),
                ("11", "2157"),
                ("12", "2759"),
                ("13", "10239"),
                ("2", "1"),
                ("2157", "1"),
                ("2759", "1"),
                ("10239", "1"),
                ("14", "15"),
                ("16", "10"),
                ("16", "2"),
                ("17", "1"),
            ],
        )
    )


@pytest.mark.parametrize(
    ("taxon", "group", "status", "roots"),
    [
        ("NCBITaxon:10", "bacteria", "resolved_root", ("NCBITaxon:2",)),
        ("NCBITaxon:11", "archaea", "resolved_root", ("NCBITaxon:2157",)),
        ("NCBITaxon:12", "eukaryota", "resolved_root", ("NCBITaxon:2759",)),
        ("NCBITaxon:13", "viruses", "resolved_root", ("NCBITaxon:10239",)),
        ("NCBITaxon:2", "bacteria", "resolved_root", ("NCBITaxon:2",)),
        ("NCBITaxon:16", "bacteria", "resolved_root", ("NCBITaxon:2",)),
        ("NCBITaxon:14", "unresolved", "no_recognized_root_path", ()),
        ("NCBITaxon:15", "unresolved", "no_recognized_root_path", ()),
        ("NCBITaxon:999", "unresolved", "not_in_taxonomy_edges", ()),
        ("local:10", "unresolved", "non_ncbi_identifier", ()),
    ],
)
def test_taxonomic_state_is_evidence_not_a_name_guess(
    taxonomy: Any, taxon: str, group: str, status: str, roots: tuple[str, ...]
) -> None:
    assert taxonomy.classify(taxon) == coverage.TaxonomyState(group, status, roots)


@pytest.mark.parametrize(
    ("taxon", "name", "category"),
    [
        ("12", "Previously unnamed fungus", "Non-Bacterial / Fungal Eukaryote"),
        ("12", "Candidatus unfamiliar lineage", "Non-Bacterial / Fungal Eukaryote"),
        ("12", "Viridiplantae", "Non-Bacterial / Fungal Eukaryote"),
        ("13", "Candidatus virus", "Non-Bacterial / Viral Taxon"),
        ("12", "uncultured organism", "Unresolved Taxonomic Label"),
        ("12", "unidentified", "Unresolved Taxonomic Label"),
        ("13", "unclassified", "Unresolved Taxonomic Label"),
        ("13", "unknown", "Unresolved Taxonomic Label"),
        ("12", "Fungal alpha / fungal beta", "Ambiguous Taxon / Bin Slash"),
        ("13", "Viral alpha / viral beta", "Ambiguous Taxon / Bin Slash"),
        ("10", "Candida", "Rare / Novel Isolate (Lacks KG Phenotype)"),
        ("11", "Archaeal isolate", "Rare / Novel Isolate (Lacks KG Phenotype)"),
        ("10", "Candidatus example", "Candidatus Lineage"),
        ("10", "uncultured bacterium", "Uncultivated / Metagenomic MAG"),
        ("10", "Unknown / composite", "Ambiguous Taxon / Bin Slash"),
        ("17", "  Unclassified  ", "Unresolved Taxonomic Label"),
        ("17", "UNKNOWN", "Unresolved Taxonomic Label"),
        ("999", "Fungus with missing lineage", "Taxonomy Unresolved"),
        ("999", "Candidatus unknown placement", "Taxonomy Unresolved"),
        ("541000", "Alias candidate", "NCBI TaxID Desynchronization"),
    ],
)
def test_lineage_and_label_priority(taxonomy: Any, taxon: str, name: str, category: str) -> None:
    identifier = f"NCBITaxon:{taxon}"
    state = taxonomy.classify(identifier)
    diagnosis, remediation = coverage.diagnose_zero_trait_taxon(identifier, name, "species", state)
    assert diagnosis == category
    if taxon == "541000":
        assert state.status == "not_in_taxonomy_edges"
        assert state.root_ids == ()
        assert remediation.startswith("Review possible alias before transferring traits")
    if category == "Unresolved Taxonomic Label":
        assert "LCA" not in remediation


@pytest.mark.parametrize(
    ("pairs", "message"),
    [
        ([("91626", "2"), ("91626", "2759")], "Conflicting recognized taxonomy roots"),
        ([("91626", "4751"), ("4751", "91626")], "Taxonomy cycle"),
        ([("91626", "91626")], "Taxonomy cycle"),
    ],
)
def test_bad_reachable_taxonomy_rejects_before_output(
    taxonomic_inputs: dict[str, Path],
    tmp_path: Path,
    pairs: list[tuple[str, str]],
    message: str,
) -> None:
    output = tmp_path / "unpublished"
    edges = write_edges(tmp_path / "invalid.tsv", pairs)
    with pytest.raises(ValueError, match=message):
        coverage.run_coverage(
            output_dir=output,
            strain_rollup="off",
            tables_only=True,
            taxonomy_edges=edges,
            **taxonomic_inputs,
        )
    assert not output.exists()


def test_absent_root_is_not_resolved_from_graph_self_inclusion(tmp_path: Path) -> None:
    taxonomy = coverage.CoverageTaxonomy(write_edges(tmp_path / "only-bacteria.tsv", [("10", "2")]))
    assert taxonomy.classify("NCBITaxon:2759").status == "not_in_taxonomy_edges"


def test_alias_candidate_does_not_inherit_suggested_target_ancestry(tmp_path: Path) -> None:
    taxonomy = coverage.CoverageTaxonomy(
        write_edges(tmp_path / "alias-target.tsv", [("216572", "2")])
    )
    state = taxonomy.classify("NCBITaxon:541000")
    assert taxonomy.classify("NCBITaxon:216572").group == "bacteria"
    assert state == coverage.TaxonomyState("unresolved", "not_in_taxonomy_edges", ())
    category, remediation = coverage.diagnose_zero_trait_taxon(
        "NCBITaxon:541000", "Alias candidate", "family", state
    )
    assert category == "NCBI TaxID Desynchronization"
    assert "Review possible alias" in remediation


def test_mapped_taxa_also_undergo_taxonomy_validation(
    taxonomic_inputs: dict[str, Path], tmp_path: Path
) -> None:
    matrix = pd.read_csv(taxonomic_inputs["matrix"], sep="\t")
    matrix.loc[matrix.taxon_id.eq("NCBITaxon:91626"), "METPO:1"] = 1
    matrix.to_csv(taxonomic_inputs["matrix"], sep="\t", index=False)
    profiles = pd.read_csv(taxonomic_inputs["profiles"], sep="\t")
    profiles["METPO:1"] = 1
    profiles.to_csv(taxonomic_inputs["profiles"], sep="\t", index=False)
    edges = write_edges(tmp_path / "conflicting-mapped.tsv", [("91626", "2"), ("91626", "2759")])
    output = tmp_path / "mapped-conflict"
    with pytest.raises(
        ValueError, match="Conflicting recognized taxonomy roots for NCBITaxon:91626"
    ):
        coverage.run_coverage(
            output_dir=output,
            strain_rollup="off",
            tables_only=True,
            taxonomy_edges=edges,
            **taxonomic_inputs,
        )
    assert not output.exists()


def test_missing_taxonomy_file_rejects_before_output(
    taxonomic_inputs: dict[str, Path], tmp_path: Path
) -> None:
    output = tmp_path / "missing-taxonomy"
    with pytest.raises(FileNotFoundError):
        coverage.run_coverage(
            output_dir=output,
            strain_rollup="off",
            tables_only=True,
            taxonomy_edges=tmp_path / "absent.tsv",
            **taxonomic_inputs,
        )
    assert not output.exists()


def test_strict_parser_preserves_known_export_repairs_and_deduplicates_edges(
    tmp_path: Path,
) -> None:
    edges = tmp_path / "split.tsv"
    edges.write_bytes(
        b"id\tsubject\tpredicate\tobject\trelation\tagent_type\tagent_type\r"
        b"\tknowledge_level\tknowledge_source\tmeta\tprimary_knowledge_source\n"
        b"e1\tNCBITaxon:10\tbiolink:subclass_of\tNCBITaxon:2759\t\t\t"
        b"manual_validation_of_automated_agent\n\t\t\t\t\n"
        b"e2\tNCBITaxon:10\tbiolink:subclass_of\tNCBITaxon:2759\t\t\t\t\t\t\t\n"
        b"e3\tMETPO:10\tbiolink:subclass_of\tMETPO:1\t\t\t\t\t\t\t\n"
    )
    taxonomy = coverage.CoverageTaxonomy(edges)
    assert taxonomy.classify("NCBITaxon:10").group == "eukaryota"
    assert taxonomy.audit["header_continuations_reconstructed"] == 1
    assert taxonomy.audit["agent_field_continuations_reconstructed"] == 1
    assert taxonomy.audit["records"] == 3
    assert taxonomy.audit["taxonomy_edge_records"] == 2
    assert taxonomy.audit["unique_taxonomy_edges"] == 1


@pytest.mark.parametrize("rank_records", [1, 2])
def test_rank_class_export_records_are_counted_without_becoming_taxa(
    tmp_path: Path, rank_records: int
) -> None:
    # The record line (including empty trailing fields) is copied verbatim from
    # KG-Microbe 20260120, edges MD5 555ecf5ba5d89288cdc44e368375374f (#1008).
    # The native header also retains its duplicate agent field and embedded CR.
    fixture = Path(__file__).parent / "fixtures/coverage/taxonomy_rank_class.tsv"
    header, rank_line = fixture.read_bytes().split(b"\n", 1)
    edges = tmp_path / "rank-and-taxon.tsv"
    edges.write_bytes(
        header
        + b"\n"
        + rank_line * rank_records
        + b"taxon\tNCBITaxon:91626\tbiolink:subclass_of\tNCBITaxon:2759"
        b"\t\t\t\t\t\t\t\n"
    )

    taxonomy = coverage.CoverageTaxonomy(edges)

    assert taxonomy.classify("NCBITaxon:91626").group == "eukaryota"
    assert taxonomy.nodes == {"NCBITaxon:91626", "NCBITaxon:2759"}
    assert taxonomy.audit["records"] == rank_records + 1
    assert taxonomy.audit["header_continuations_reconstructed"] == 1
    assert taxonomy.audit["duplicate_agent_fields_preserved"] == 1
    assert taxonomy.audit["taxonomy_rank_class_edge_records_skipped"] == rank_records
    assert taxonomy.audit["taxonomy_edge_records"] == 1
    assert taxonomy.audit["unique_taxonomy_edges"] == 1


def test_rank_class_skip_count_is_emitted_in_coverage_provenance(
    taxonomic_inputs: dict[str, Path], tmp_path: Path
) -> None:
    fixture = Path(__file__).parent / "fixtures/coverage/taxonomy_rank_class.tsv"
    edges = tmp_path / "taxonomy.tsv"
    edges.write_bytes(
        fixture.read_bytes() + b"taxon\tNCBITaxon:91626\tbiolink:subclass_of\tNCBITaxon:2759"
        b"\t\t\t\t\t\t\t\n"
    )
    output = tmp_path / "coverage"
    coverage.run_coverage(
        output_dir=output,
        strain_rollup="off",
        tables_only=True,
        taxonomy_edges=edges,
        **taxonomic_inputs,
    )

    provenance = json.loads((output / "provenance.json").read_text())
    audit = provenance["taxonomy_diagnosis"]["parser_audit"]
    assert audit["taxonomy_rank_class_edge_records_skipped"] == 1
    assert audit["records"] == 2
    assert audit["taxonomy_edge_records"] == 1


@pytest.mark.parametrize(
    "content",
    [
        "subject\tpredicate\tobject\nNCBITaxon:10\tbiolink:subclass_of\tNCBITaxon:2\n",
        "id\tsubject\tpredicate\tobject\ne1\tNCBITaxon:nope\tbiolink:subclass_of\tNCBITaxon:2\n",
        "id\tsubject\tpredicate\tobject\ne1\tNCBITaxon:biotype\tbiolink:subclass_of\tOBO:NCBITaxon#_taxonomic_rank_extra\n",
        "id\tsubject\tpredicate\tobject\ne1\tNCBITaxon:10\tbiolink:subclass_of\tNCBITaxon:biotype\n",
        "id\tsubject\tpredicate\tobject\ne1\tNCBITaxon:10\tbiolink:subclass_of\n",
        "id\tsubject\tpredicate\tobject\ne1\tNCBITaxon:10\tbiolink:capable_of\tNCBITaxon:2\n",
    ],
)
def test_malformed_or_unusable_taxonomy_fails_closed(tmp_path: Path, content: str) -> None:
    edges = tmp_path / "bad.tsv"
    edges.write_text(content)
    with pytest.raises(
        ValueError, match=r"expected leading|Invalid taxonomy|malformed|no NCBITaxon"
    ):
        coverage.CoverageTaxonomy(edges)


def test_taxonomy_changes_diagnosis_but_preserves_all_coverage_counts(
    taxonomic_inputs: dict[str, Path], tmp_path: Path
) -> None:
    bacterial = write_edges(tmp_path / "bacterial.tsv", [("91626", "2"), ("4896", "2")])
    fungal = write_edges(tmp_path / "fungal.tsv", [("91626", "2759"), ("4896", "2759")])
    frames = []
    for edges in (bacterial, fungal):
        frames.append(coverage.load_and_process_data(taxonomy_edges=edges, **taxonomic_inputs))
    for before, after in zip(frames[0][:-1], frames[1][:-1], strict=True):
        pd.testing.assert_frame_equal(before, after)
    assert frames[0][-1].diagnosis_category.tolist() != frames[1][-1].diagnosis_category.tolist()
    assert frames[1][-1].taxonomy_group.tolist().count("eukaryota") == 2
    assert frames[1][-1].taxonomy_root_ids.tolist().count('["NCBITaxon:2759"]') == 2


def test_taxonomy_provenance_and_final_input_mutation_guard(
    taxonomic_inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    edges = write_edges(tmp_path / "taxonomy.tsv", [("91626", "2759"), ("4896", "2759")])
    output = tmp_path / "valid"
    kwargs = dict(strain_rollup="off", tables_only=True, taxonomy_edges=edges, **taxonomic_inputs)
    coverage.run_coverage(output_dir=output, **kwargs)
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["inputs"]["taxonomy_edges"]["sha256"] == coverage.file_sha256(edges)
    assert provenance["taxonomy_diagnosis"]["matrix_kg_identity"] == "not_verified"
    assert provenance["taxonomy_diagnosis"]["parser_audit"]["taxonomy_edge_records"] == 2
    for field, path in [
        ("taxonomy_parser_sha256", "src/mtp/workflow/kg_inputs.py"),
        ("taxonomy_graph_sha256", "src/mtp/contribution_inputs.py"),
    ]:
        assert provenance["producer"][field] == coverage.file_sha256(coverage.REPO_ROOT / path)
    original = coverage.write_readme

    def change_taxonomy(*args: Any, **kwargs: Any) -> None:
        original(*args, **kwargs)
        edges.write_text(edges.read_text() + "e3\tNCBITaxon:10\tbiolink:subclass_of\tNCBITaxon:2\n")

    monkeypatch.setattr(coverage, "write_readme", change_taxonomy)
    with pytest.raises(
        ValueError, match="Input changed during coverage generation: taxonomy_edges"
    ):
        coverage.run_coverage(output_dir=tmp_path / "changed", **kwargs)
    assert not (tmp_path / "changed/provenance.json").exists()
