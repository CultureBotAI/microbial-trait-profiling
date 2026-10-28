"""Carrier imports preserve bound taxa, complete evidence and unknown record counts (#1044)."""

from __future__ import annotations

import csv
import io
import json
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from mtp.datamodel.workflow_pydantic import ArtifactInstance, ComparisonRun, WorkflowRun
from mtp.label_support_audit import AUDIT_TRAITS, run_label_support_audit
from mtp.multiplicity.sources import VerifiedSources
from mtp.paths import file_sha256
from mtp.three_test_table import carriers
from mtp.three_test_table.carriers import carrier_table
from multiplicity_fixture import fixture_sources, write_table

CELL = "METPO:2000103|GO:0030245"
ANAEROBIC = "METPO:1000603"
OUTSIDE = "METPO:1000429"
LIGNIN = "METPO:2000103|GO:0046274"
HIGH_TEMP = "METPO:1000447"


def audit_fixture(root: Path) -> tuple[VerifiedSources, Path]:
    """Generate a real audit with nested carriers and out-of-profile occurrences."""
    root.mkdir()
    products = {
        name: root / f"{name}.tsv"
        for name in ("count_profiles", "taxa_matrix", "reported_taxa", "kg_edges", "kg_nodes")
    }
    matrix = [
        {
            "taxon_id": f"NCBITaxon:{taxon}",
            "taxon_name": f"Fixture {taxon}",
            "taxonomic_rank": "genus" if taxon in (10, 20) else "species",
            **dict.fromkeys(AUDIT_TRAITS, 0),
            CELL: cell,
            OUTSIDE: other,
        }
        for taxon, cell, other in ((10, 2, 1), (11, 1, 1), (20, 1, 0))
    ]
    write_table(products["taxa_matrix"], matrix)
    profiles = [
        {
            "study_id": "S",
            "experiment_id": "E",
            "inc_dec": direction,
            **{feature: row[feature] for feature in (*AUDIT_TRAITS, OUTSIDE)},
        }
        for direction, row in zip(("inc", "dec"), matrix[:2], strict=True)
    ]
    write_table(products["count_profiles"], profiles)
    occurrences = [
        {
            "study_id": study,
            "experiment_id": "E",
            "direction": direction,
            "taxon_id": f"NCBITaxon:{taxon}",
        }
        for study, direction, taxon in (
            ("S", "increased", 10),
            ("S", "increased", 10),
            ("S", "decreased", 11),
            ("out-of-scope", "increased", 20),
        )
    ]
    write_table(products["reported_taxa"], occurrences)
    edges = [
        {
            "id": f"edge-{index}",
            "subject": subject,
            "predicate": predicate,
            "object": obj,
            "primary_knowledge_source": source,
        }
        for index, (subject, predicate, obj, source) in enumerate(
            [
                ("NCBITaxon:11", "biolink:subclass_of", "NCBITaxon:10", "infores:ncbi"),
                ("NCBITaxon:12", "biolink:subclass_of", "NCBITaxon:10", "infores:ncbi"),
                ("NCBITaxon:21", "biolink:subclass_of", "NCBITaxon:20", "infores:ncbi"),
                (
                    "NCBITaxon:11",
                    "METPO:2000103",
                    "GO:0030245",
                    "infores:madin_etal|infores:bactotraits",
                ),
                ("NCBITaxon:12", "METPO:2000103", "GO:0030245", "infores:madin_etal"),
                ("NCBITaxon:21", "METPO:2000103", "GO:0030245", "infores:madin_etal"),
                ("NCBITaxon:11", "biolink:has_phenotype", OUTSIDE, "infores:bactotraits"),
            ]
        )
    ]
    write_table(products["kg_edges"], edges)
    write_table(
        products["kg_nodes"],
        [
            {"id": f"NCBITaxon:{taxon}", "name": f"Fixture {taxon}"}
            for taxon in (10, 11, 12, 20, 21)
        ],
    )
    audit = root / "audit"
    run_label_support_audit(
        taxa_path=products["reported_taxa"],
        matrix_path=products["taxa_matrix"],
        profiles_path=products["count_profiles"],
        edges_path=products["kg_edges"],
        nodes_path=products["kg_nodes"],
        output_dir=audit,
        strain_rollup="off",
    )
    core_manifest = root / "core.json"
    core_manifest.write_text('{"synthetic": true}\n', encoding="utf-8")
    core = WorkflowRun.model_construct(
        manifest_path=str(core_manifest),
        artifacts=[
            ArtifactInstance.model_construct(product_id=name, path=str(path))
            for name, path in products.items()
        ],
    )
    sources = VerifiedSources(
        ComparisonRun.model_construct(),
        core,
        {},
        {},
        {"count_profiles": (*AUDIT_TRAITS, OUTSIDE)},
        "off",
        {path.resolve(): file_sha256(path) for path in products.values()},
    )
    return sources, audit


def rewrite_output(audit: Path, name: str, frame: pd.DataFrame) -> None:
    """Refresh self-declared hashes to exercise semantic validation beyond checksums."""
    path = audit / name
    frame.to_csv(path, sep="\t", index=False)
    manifest = audit / "provenance.json"
    record = json.loads(manifest.read_text(encoding="utf-8"))
    record["outputs"][name].update(sha256=file_sha256(path), size_bytes=path.stat().st_size)
    if name == "taxon_trait_label_support.tsv":
        record["support_cells"] = len(frame)
    manifest.write_text(json.dumps(record), encoding="utf-8")


def test_carrier_positive_audit_counts_unique_carriers_not_source_records(tmp_path: Path) -> None:
    sources, audit = audit_fixture(tmp_path / "inputs")
    table, paths = carrier_table(sources, audit)
    result = table.set_index("trait")
    cell = result.loc[CELL]
    assert cell.carrier_status == "audit_hash_bound_focal_coverage"
    assert cell.carrier_direct_reported_taxa == 1
    assert cell.carrier_inherited_reported_taxa == 1
    assert cell.carrier_unique_taxa == 2
    assert pd.isna(cell.carrier_source_record_count)
    assert cell.carrier_population == "reported_taxa_observed_in_supplied_count_profiles"
    assert cell.carrier_source_record_status == "unresolved_original_record_identifiers_unavailable"
    assert result.loc[ANAEROBIC, "carrier_unique_taxa"] == 0
    assert result.loc[OUTSIDE, "carrier_status"] == "not_in_focal_audit"
    assert pd.isna(result.loc[OUTSIDE, "carrier_unique_taxa"])
    assert audit / "provenance.json" in paths
    assert audit / "carrier_source_evidence.tsv" in paths


@pytest.mark.parametrize("supplied", [False, True])
def test_carrier_nullable_counts_preserve_unresolved_and_observed_values(
    tmp_path: Path, supplied: bool
) -> None:
    sources, audit = audit_fixture(tmp_path / "inputs")
    table, _ = carrier_table(sources, audit if supplied else None)
    columns = (
        "carrier_direct_reported_taxa",
        "carrier_inherited_reported_taxa",
        "carrier_unique_taxa",
        "carrier_source_record_count",
    )
    assert all(str(table[column].dtype) == "Int64" for column in columns), (
        "Carrier counts must use nullable integers"
    )
    exported = {
        row["trait"]: row
        for row in csv.DictReader(io.StringIO(table.to_csv(sep="\t", index=False)), delimiter="\t")
    }
    assert all(row["carrier_source_record_count"] == "" for row in exported.values())
    if supplied:
        assert [exported[CELL][column] for column in columns[:3]] == ["1", "1", "2"]
        assert [exported[ANAEROBIC][column] for column in columns[:3]] == ["0", "0", "0"]
        assert all(exported[OUTSIDE][column] == "" for column in columns)
    else:
        assert all(row[column] == "" for row in exported.values() for column in columns)


@pytest.mark.parametrize(
    ("damage", "expected_error"),
    [
        ("missing_support", "^Label-support audit is missing positive core matrix cells$"),
        (
            "missing_evidence",
            "^Label-support carrier evidence disagrees with support counts$",
        ),
        ("duplicate_evidence", "^Duplicate label-support carrier evidence$"),
        ("orphan_evidence", "^Orphan or malformed label-support carrier evidence$"),
        ("wrong_carrier_total", "^Label-support carrier counts differ from the core matrix$"),
        ("wrong_inherited_count", "^Label-support carrier counts differ from the core matrix$"),
        ("false_direct", "^Label-support carrier counts differ from the core matrix$"),
        (
            "scope_occurrence",
            "^Label-support signature occurrence count differs from bound scope$",
        ),
        (
            "negative_occurrence",
            r"^Expected a finite nonnegative integer count in .*taxon_trait_label_support\.tsv:0, "
            r"column 'signature_occurrences'; got '-1'$",
        ),
        (
            "fractional_count",
            r"^Expected a finite nonnegative integer count in .*taxon_trait_label_support\.tsv:0, "
            r"column 'inherited_carrier_count'; got '0\.5'$",
        ),
        ("blank_carrier", "^Orphan or malformed label-support carrier evidence$"),
        ("wrong_scope_facts", "^Label-support audit scope facts differ from core inputs$"),
        (
            "missing_carrier_column",
            r"^Label-support carrier_source_evidence\.tsv lacks required columns$",
        ),
    ],
)
def test_carrier_rejects_self_rehashed_incomplete_or_false_evidence(
    tmp_path: Path, damage: str, expected_error: str
) -> None:
    sources, audit = audit_fixture(tmp_path / "inputs")
    support_name = "taxon_trait_label_support.tsv"
    evidence_name = "carrier_source_evidence.tsv"
    support = pd.read_csv(audit / support_name, sep="\t", dtype=str, keep_default_na=False)
    evidence = pd.read_csv(audit / evidence_name, sep="\t", dtype=str, keep_default_na=False)
    parent = support.index[support.taxon_id == "NCBITaxon:10"][0]
    if damage == "missing_support":
        rewrite_output(audit, support_name, support.drop(index=parent))
    elif damage == "missing_evidence":
        rewrite_output(audit, evidence_name, evidence.iloc[1:])
    elif damage == "duplicate_evidence":
        rewrite_output(
            audit, evidence_name, pd.concat([evidence, evidence.iloc[:1]], ignore_index=True)
        )
    elif damage == "orphan_evidence":
        evidence.loc[0, "reported_taxon_id"] = "NCBITaxon:999"
        rewrite_output(audit, evidence_name, evidence)
    elif damage == "blank_carrier":
        evidence.loc[0, "carrier_taxon_id"] = ""
        rewrite_output(audit, evidence_name, evidence)
    elif damage == "missing_carrier_column":
        rewrite_output(audit, evidence_name, evidence.drop(columns="carrier_taxon_id"))
    elif damage == "wrong_scope_facts":
        manifest = audit / "provenance.json"
        record = json.loads(manifest.read_text())
        record["scope"]["bound_profiles"] = 999
        manifest.write_text(json.dumps(record))
    else:
        column, value = {
            "wrong_carrier_total": ("carrier_count", "99"),
            "wrong_inherited_count": ("inherited_carrier_count", "0"),
            "false_direct": ("carrier_at_reported_taxon", "True"),
            "scope_occurrence": ("signature_occurrences", "0"),
            "negative_occurrence": ("signature_occurrences", "-1"),
            "fractional_count": ("inherited_carrier_count", "0.5"),
        }[damage]
        support.loc[parent, column] = value
        rewrite_output(audit, support_name, support)
    with pytest.raises(ValueError, match=expected_error):
        carrier_table(sources, audit)


@pytest.mark.parametrize(
    ("taxon", "direct", "inherited"),
    [("NCBITaxon:10", "True", "1"), ("NCBITaxon:11", "False", "1")],
)
def test_carrier_rejects_count_preserving_direct_carrier_conflict(
    tmp_path: Path, taxon: str, direct: str, inherited: str
) -> None:
    """Evidence must reconcile direct identity even when all carrier totals agree (#1077)."""
    sources, audit = audit_fixture(tmp_path / "inputs")
    name = "taxon_trait_label_support.tsv"
    support = pd.read_csv(audit / name, sep="\t", dtype=str, keep_default_na=False)
    index = support.index[(support.taxon_id == taxon) & (support.trait_key == CELL)][0]
    assert support.loc[index, "carrier_at_reported_taxon"] != direct
    support.loc[index, "carrier_at_reported_taxon"] = direct
    support.loc[index, "inherited_carrier_count"] = inherited
    assert int(support.loc[index, "carrier_count"]) == int(inherited) + int(direct == "True")
    rewrite_output(audit, name, support)
    with pytest.raises(
        ValueError, match=r"^Label-support carrier evidence disagrees with support counts$"
    ):
        carrier_table(sources, audit)


def test_carrier_rejects_provenance_changed_during_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A valid inventory does not authorize a later replacement of its manifest (#1077)."""
    sources, audit = audit_fixture(tmp_path / "inputs")
    manifest = audit / "provenance.json"
    evidence = audit / "carrier_source_evidence.tsv"
    changed = False

    def hash_then_change_manifest(path: Path) -> str:
        nonlocal changed
        digest = file_sha256(path)
        if path == evidence and not changed:
            manifest.write_bytes(manifest.read_bytes() + b"\n")
            changed = True
        return digest

    monkeypatch.setattr(carriers, "file_sha256", hash_then_change_manifest)
    with pytest.raises(
        ValueError, match=r"^Label-support provenance changed while reading the audit$"
    ):
        carrier_table(sources, audit)
    assert changed


@pytest.mark.parametrize("role", ["profiles", "matrix", "taxa", "kg_edges", "kg_nodes"])
def test_carrier_rejects_audit_from_other_core_input(tmp_path: Path, role: str) -> None:
    sources, audit = audit_fixture(tmp_path / "inputs")
    manifest = audit / "provenance.json"
    record = json.loads(manifest.read_text())
    record["inputs"][role]["sha256"] = "0" * 64
    manifest.write_text(json.dumps(record))
    with pytest.raises(ValueError, match=f"verified core: {role}"):
        carrier_table(sources, audit)


def test_carrier_rejects_rollup_scope_and_uninventoried_changes(tmp_path: Path) -> None:
    sources, audit = audit_fixture(tmp_path / "inputs")
    with pytest.raises(ValueError, match="scope"):
        carrier_table(replace(sources, scope="phenotype"), audit)
    with (audit / "carrier_source_evidence.tsv").open("a") as handle:
        handle.write("tampered\n")
    with pytest.raises(ValueError, match="inventory"):
        carrier_table(sources, audit)


def test_missing_audit_is_unavailable_not_a_zero_count(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    table, paths = carrier_table(sources, None)
    assert list(table.trait) == list(sources.features["count_profiles"])
    assert table.carrier_status.eq("not_supplied").all()
    assert table.carrier_direct_reported_taxa.isna().all()
    assert table.carrier_source_record_count.isna().all()
    assert paths == ()


def test_incomplete_audit_fails_instead_of_supplying_zero_carriers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    audit = tmp_path / "audit"
    audit.mkdir()
    (audit / "provenance.json").write_text('{"analysis": "label_support_sidecar"}')
    with pytest.raises(ValueError, match=r"(scope|provenance|audit)"):
        carrier_table(sources, audit)


def _declared_list_audit_fixture(root: Path) -> tuple[VerifiedSources, Path, Path]:
    """Generate an audit with declared traits from taxonomic inheritance and strain rollup."""
    sources, _ = audit_fixture(root)
    products = {
        pid: sources.product(pid)
        for pid in ("taxa_matrix", "count_profiles", "reported_taxa", "kg_edges", "kg_nodes")
    }
    matrix_df = pd.read_csv(products["taxa_matrix"], sep="\t")
    matrix_df[LIGNIN] = [1 if tid == "NCBITaxon:10" else 0 for tid in matrix_df["taxon_id"]]
    matrix_df[HIGH_TEMP] = [
        1 if tid in ("NCBITaxon:10", "NCBITaxon:11") else 0 for tid in matrix_df["taxon_id"]
    ]
    matrix_df.to_csv(products["taxa_matrix"], sep="\t", index=False)

    profiles_df = pd.read_csv(products["count_profiles"], sep="\t")
    profiles_df[LIGNIN] = [1 if d == "inc" else 0 for d in profiles_df["inc_dec"]]
    profiles_df[HIGH_TEMP] = [1 for _ in profiles_df["inc_dec"]]
    profiles_df.to_csv(products["count_profiles"], sep="\t", index=False)

    with products["kg_edges"].open("a", encoding="utf-8") as handle:
        handle.write("edge-lig-12\tNCBITaxon:12\tMETPO:2000103\tGO:0046274\tinfores:madin_etal\n")
        handle.write(
            "edge-strain-sub\tkgmicrobe.strain:thermo\tbiolink:subclass_of\t"
            "NCBITaxon:11\tinfores:bacdive\n"
        )
        handle.write(
            f"edge-strain-temp\tkgmicrobe.strain:thermo\tbiolink:has_phenotype\t"
            f"{HIGH_TEMP}\tinfores:bacdive\n"
        )
    with products["kg_nodes"].open("a", encoding="utf-8") as handle:
        handle.write("GO:0046274\tlignin catabolic process\n")
        handle.write(f"{HIGH_TEMP}\ttemperature optimum high\n")

    trait_list_path = root / "headline_ten_phenotype.txt"
    trait_list_path.write_text(f"{CELL}\n{LIGNIN}\n{HIGH_TEMP}\n", encoding="utf-8")
    audit_declared = root / "audit_declared"
    run_label_support_audit(
        taxa_path=products["reported_taxa"],
        matrix_path=products["taxa_matrix"],
        profiles_path=products["count_profiles"],
        edges_path=products["kg_edges"],
        nodes_path=products["kg_nodes"],
        output_dir=audit_declared,
        strain_rollup="phenotype",
        trait_list=trait_list_path,
    )
    updated_sources = VerifiedSources(
        sources.comparison,
        sources.core,
        {},
        {},
        {"count_profiles": (*AUDIT_TRAITS, OUTSIDE, LIGNIN, HIGH_TEMP)},
        "phenotype",
        {path.resolve(): file_sha256(path) for path in products.values()},
    )
    return updated_sources, audit_declared, trait_list_path


def test_carrier_table_imports_declared_trait_list_taxonomic_inheritance_and_strain_rollup(
    tmp_path: Path,
) -> None:
    sources, audit, trait_list_path = _declared_list_audit_fixture(tmp_path / "inputs")
    table, _ = carrier_table(sources, audit)
    result = table.set_index("trait")
    expected_coverage = f"focal_and_selected_oxygen_traits_plus_{trait_list_path.name}"
    assert result["carrier_catalog_coverage"].eq(expected_coverage).all()

    lignin = result.loc[LIGNIN]
    assert lignin.carrier_status == "audit_hash_bound_focal_coverage"
    assert lignin.carrier_direct_reported_taxa == 0
    assert lignin.carrier_inherited_reported_taxa == 1
    assert lignin.carrier_unique_taxa == 1
    assert pd.isna(lignin.carrier_source_record_count)
    assert lignin.carrier_madin_status == "not_supplied_source_support_unknown"

    high_temp = result.loc[HIGH_TEMP]
    assert high_temp.carrier_status == "audit_hash_bound_focal_coverage"
    assert high_temp.carrier_direct_reported_taxa == 1
    assert high_temp.carrier_inherited_reported_taxa == 1
    assert high_temp.carrier_unique_taxa == 1
    assert pd.isna(high_temp.carrier_source_record_count)
    assert high_temp.carrier_madin_status == "not_supplied_source_support_unknown"

    outside = result.loc[OUTSIDE]
    assert outside.carrier_status == "not_in_focal_audit"
    assert pd.isna(outside.carrier_unique_taxa)


def test_carrier_table_rejects_tampered_trait_list_provenance(tmp_path: Path) -> None:
    sources, audit, _ = _declared_list_audit_fixture(tmp_path / "inputs")
    manifest = audit / "provenance.json"
    record = json.loads(manifest.read_text(encoding="utf-8"))
    record["trait_list"]["sha256"] = "0" * 64
    manifest.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(
        ValueError, match=r"^Malformed or mismatched label-support trait_list provenance$"
    ):
        carrier_table(sources, audit)


@pytest.mark.parametrize("missing_side", ["record", "inputs"])
def test_carrier_table_rejects_incomplete_trait_list_provenance(
    tmp_path: Path, missing_side: str
) -> None:
    sources, audit, _ = _declared_list_audit_fixture(tmp_path / "inputs")
    manifest = audit / "provenance.json"
    record = json.loads(manifest.read_text(encoding="utf-8"))
    if missing_side == "record":
        del record["trait_list"]
    else:
        del record["inputs"]["trait_list"]
    manifest.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(
        ValueError, match=r"^Label-support audit trait_list provenance is incomplete$"
    ):
        carrier_table(sources, audit)


def test_carrier_table_rejects_trait_list_path_mismatch_against_inputs(tmp_path: Path) -> None:
    sources, audit, trait_list_path = _declared_list_audit_fixture(tmp_path / "inputs")
    manifest = audit / "provenance.json"
    record = json.loads(manifest.read_text(encoding="utf-8"))
    record["trait_list"]["path"] = f"other_dir/{trait_list_path.name}"
    manifest.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(
        ValueError, match=r"^Malformed or mismatched label-support trait_list provenance$"
    ):
        carrier_table(sources, audit)


def test_carrier_table_rejects_trait_list_name_mismatch_against_path(tmp_path: Path) -> None:
    sources, audit, _ = _declared_list_audit_fixture(tmp_path / "inputs")
    manifest = audit / "provenance.json"
    record = json.loads(manifest.read_text(encoding="utf-8"))
    record["trait_list"]["name"] = "wrong_name.txt"
    manifest.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(
        ValueError, match=r"^Malformed or mismatched label-support trait_list provenance$"
    ):
        carrier_table(sources, audit)


@pytest.mark.parametrize("mutation", ["wrong_count", "unstripped_id"])
def test_carrier_table_rejects_trait_list_ids_or_count_mismatch(
    tmp_path: Path, mutation: str
) -> None:
    sources, audit, _ = _declared_list_audit_fixture(tmp_path / "inputs")
    manifest = audit / "provenance.json"
    record = json.loads(manifest.read_text(encoding="utf-8"))
    if mutation == "wrong_count":
        record["trait_list"]["trait_count"] = 99
    else:
        record["trait_list"]["trait_ids"] = [f" {CELL} ", LIGNIN, HIGH_TEMP]
    manifest.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(
        ValueError, match=r"^Label-support trait_list identifiers or count mismatch$"
    ):
        carrier_table(sources, audit)


def test_carrier_table_rejects_trait_list_outside_core_profile_features(tmp_path: Path) -> None:
    sources, audit, _ = _declared_list_audit_fixture(tmp_path / "inputs")
    truncated_sources = replace(
        sources,
        features={"count_profiles": (*AUDIT_TRAITS, OUTSIDE, LIGNIN)},
    )
    with pytest.raises(
        ValueError, match=r"^Label-support trait_list contains traits outside core profiles$"
    ):
        carrier_table(truncated_sources, audit)
