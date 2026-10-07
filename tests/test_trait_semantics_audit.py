"""Tests for Table 1 trait semantics, ontology overlap, and source evidence audit."""

from __future__ import annotations

from pathlib import Path

import pytest

from bsdbng.trait_semantics import (
    DEFAULT_OUTPUT_TSV,
    PINNED_KG_PROVENANCE,
    TABLE_1_TRAIT_KEYS,
    TRAIT_SEMANTIC_SPECS,
    CoccusSphereComparison,
    build_semantic_audit_records,
    read_semantic_audit_tsv,
    verify_provenance_checksums,
    write_semantic_audit_tsv,
)


def test_provenance_verification() -> None:
    """Verify that KG-Microbe provenance checksums match the pinned build."""
    assert verify_provenance_checksums() is True
    assert "merged-kg_edges.tsv" in PINNED_KG_PROVENANCE
    assert "merged-kg_nodes.tsv" in PINNED_KG_PROVENANCE
    assert PINNED_KG_PROVENANCE["merged-kg_edges.tsv"]["md5"] == "555ecf5ba5d89288cdc44e368375374f"
    assert PINNED_KG_PROVENANCE["merged-kg_nodes.tsv"]["md5"] == "dd76b4656f9341db25a7c4afa96cd58e"


def test_semantic_audit_records_coverage() -> None:
    """Verify that all six canonical Table 1 traits are covered with required attributes."""
    records = build_semantic_audit_records()
    assert len(records) == 6

    record_keys = [r.trait_key for r in records]
    assert record_keys == TABLE_1_TRAIT_KEYS

    for r in records:
        assert r.feature_name
        assert r.ontology_id
        assert r.ontology_label
        assert r.ontology_definition
        assert r.kg_edge_predicate
        assert r.source_knowledge_source
        assert r.original_source_term
        assert r.evidence_type
        assert r.taxonomic_resolution
        assert r.direct_kg_taxa_count > 0
        assert r.bugsigdb_matrix_taxa_count > 0
        assert r.profile_support_count > 0
        assert 0.0 < r.profile_support_pct <= 100.0
        assert r.evidence_status
        assert r.decision in {"keep_separate", "unresolved", "merge"}
        assert r.approved_wording
        assert r.unsupported_interpretations


def test_coccus_sphere_comparison_metrics() -> None:
    """Verify quantitative overlap and semantic distinction between coccus and sphere."""
    comp = CoccusSphereComparison()
    assert comp.ontology_relationship == "Sibling classes under METPO:1000666 (cell shape)"
    assert comp.coccus_direct_taxa == 5614
    assert comp.sphere_direct_taxa == 381
    assert comp.direct_taxa_intersection == 266
    assert pytest.approx(comp.direct_taxa_jaccard, rel=1e-3) == 0.0464

    assert comp.coccus_bugsigdb_taxa == 675
    assert comp.sphere_bugsigdb_taxa == 283
    assert comp.bugsigdb_taxa_intersection == 262
    assert pytest.approx(comp.bugsigdb_taxa_jaccard, rel=1e-3) == 0.3764

    assert comp.coccus_profile_nonzero_count == 4409
    assert comp.sphere_profile_nonzero_count == 3541
    assert pytest.approx(comp.pearson_correlation, rel=1e-3) == 0.9110
    assert pytest.approx(comp.spearman_correlation, rel=1e-3) == 0.8468
    assert comp.decision == "keep_separate"


def test_arsenate_reduction_semantic_constraints() -> None:
    """Verify that arsenate reduction cannot be classified as arsC or arrAB respiration."""
    spec = TRAIT_SEMANTIC_SPECS["METPO:2000103|pathways:arsenate_reduction"]
    assert spec["decision"] == "unresolved"
    assert "unresolved mechanism" in spec["evidence_status"]
    assert "0/51 KG edges resolve gene or pathway mechanism" in spec["evidence_status"]
    assert "Do not classify as arsC detoxification" in spec["unsupported_interpretations"]
    assert (
        "Do not classify as arrAB dissimilatory respiration" in spec["unsupported_interpretations"]
    )
    assert "Do not infer gene presence, plasmid carriage" in spec["unsupported_interpretations"]


def test_moderately_halophilic_semantic_constraints() -> None:
    """Verify that moderately halophilic requires salt and cannot be read as halotolerant."""
    spec = TRAIT_SEMANTIC_SPECS["METPO:1000623"]
    assert spec["decision"] == "keep_separate"
    assert "requires high levels of sodium chloride" in spec["ontology_definition"]
    assert "Do not read as 'halotolerant'" in spec["unsupported_interpretations"]
    assert "unresolved" in spec["evidence_status"]


def test_write_and_read_semantic_audit_tsv(tmp_path: Path) -> None:
    """Verify round-trip export and read of semantic audit TSV."""
    records = build_semantic_audit_records()
    tmp_tsv = tmp_path / "table1_semantic_audit_test.tsv"

    exported_path = write_semantic_audit_tsv(records, tmp_tsv)
    assert exported_path.exists()

    loaded = read_semantic_audit_tsv(tmp_tsv)
    assert len(loaded) == len(records)
    for orig, roundtripped in zip(records, loaded, strict=True):
        assert orig == roundtripped


def test_derived_audit_tsv_exists_and_matches() -> None:
    """Verify that the committed derived TSV matches the authoritative audit specification."""
    assert DEFAULT_OUTPUT_TSV.exists(), f"Missing required derived TSV: {DEFAULT_OUTPUT_TSV}"
    records = read_semantic_audit_tsv(DEFAULT_OUTPUT_TSV)
    expected = build_semantic_audit_records()
    assert len(records) == 6
    assert records == expected
