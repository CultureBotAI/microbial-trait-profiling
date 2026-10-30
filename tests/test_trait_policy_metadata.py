"""Versioned feature semantics preserve historical matrices and expose unions (#1129)."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Literal

import pytest

from mtp.taxa_trait_matrix import (
    NodeMetadata,
    column_metadata_rows,
    create_taxa_trait_matrix,
    make_trait_key,
    read_kg_edges,
)
from mtp.workflow.contracts import read_table, validate_product
from mtp.workflow.coverage import write_feature_catalog
from mtp.workflow.matrix_contracts import validate_matrix_sources

Policy = Literal["raw_object_v1", "facultative_union_v1"]
RAW: Policy = "raw_object_v1"
UNION: Policy = "facultative_union_v1"
MEMBERS = ["METPO:1000605", "METPO:1000608", "METPO:1000612"]
CANONICAL = MEMBERS[1]
POLICY_COLUMNS = (
    "trait_key_policy",
    "trait_member_ids_json",
    "trait_source_vocabularies_json",
)
OLD_COLUMNS = (
    "column_id",
    "column_label",
    "predicate_id",
    "predicate_name",
    "object_id",
    "object_name",
    "object_category",
)


def _write(path: Path, header: list[str], rows: list[list[str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


@pytest.fixture
def sources(tmp_path: Path) -> dict[str, Path]:
    edges = [
        ["t1", "NCBITaxon:1", "biolink:subclass_of", "NCBITaxon:10", "fixture:taxonomy"],
        ["t2", "NCBITaxon:2", "biolink:subclass_of", "NCBITaxon:10", "fixture:taxonomy"],
        ["t3", "NCBITaxon:3", "biolink:subclass_of", "NCBITaxon:10", "fixture:taxonomy"],
        ["a1", "NCBITaxon:1", "biolink:has_phenotype", MEMBERS[0], "infores:bacdive"],
        ["a2", "NCBITaxon:1", "biolink:has_phenotype", MEMBERS[1], "infores:madin_etal"],
        ["a3", "NCBITaxon:1", "biolink:has_phenotype", MEMBERS[2], "infores:bactotraits"],
        ["a4", "NCBITaxon:1", "biolink:has_phenotype", MEMBERS[0], "fixture:duplicate"],
        ["a5", "NCBITaxon:2", "biolink:has_phenotype", MEMBERS[2], "infores:bactotraits"],
        ["s2", "kgmicrobe.strain:two", "biolink:subclass_of", "NCBITaxon:2", "infores:bacdive"],
        ["s3", "kgmicrobe.strain:three", "biolink:subclass_of", "NCBITaxon:3", "infores:bacdive"],
        ["s4", "kgmicrobe.strain:two", "biolink:has_phenotype", MEMBERS[0], "infores:bacdive"],
        ["s5", "kgmicrobe.strain:three", "biolink:has_phenotype", MEMBERS[0], "infores:bacdive"],
    ]
    taxa = [
        ["increased", f"NCBITaxon:{i}", f"Named taxon {i}", "genus" if i == 10 else "species"]
        for i in (10, 1, 2, 3)
    ]
    nodes = [
        [member, "biolink:PhenotypicQuality", name]
        for member, name in zip(
            MEMBERS,
            ("facultatively anaerobic", "facultatively aerobic", "facultative oxygen preference"),
            strict=True,
        )
    ]
    return {
        "edges": _write(
            tmp_path / "edges.tsv",
            ["id", "subject", "predicate", "object", "primary_knowledge_source"],
            edges,
        ),
        "taxa": _write(
            tmp_path / "taxa.tsv", ["direction", "taxon_id", "taxon_name", "taxonomic_rank"], taxa
        ),
        "nodes": _write(tmp_path / "nodes.tsv", ["id", "category", "name"], nodes),
    }


def test_policy_is_explicit_and_cannot_leak_between_calls(sources: dict[str, Path]) -> None:
    for policy, expected in ((RAW, MEMBERS[0]), (UNION, CANONICAL), (RAW, MEMBERS[0])):
        assert (
            make_trait_key("biolink:has_phenotype", MEMBERS[0], trait_key_policy=policy) == expected
        )
        _, traits = read_kg_edges(sources["edges"], trait_key_policy=policy)
        assert traits["NCBITaxon:1"] == (set(MEMBERS) if policy == RAW else {CANONICAL})
    negated = "METPO:2000021"
    assert make_trait_key(negated, MEMBERS[0], trait_key_policy=RAW) == f"{negated}|{MEMBERS[0]}"
    assert make_trait_key(negated, MEMBERS[0], trait_key_policy=UNION) == f"{negated}|{CANONICAL}"
    assert (
        make_trait_key("METPO:2000001", "GO:0030245", trait_key_policy=UNION)
        == "METPO:2000001|GO:0030245"
    )


@pytest.mark.parametrize("policy", [RAW, UNION])
@pytest.mark.parametrize("scope", ["off", "phenotype"])
@pytest.mark.parametrize("modern", [False, True])
def test_matrix_and_independent_source_check_use_declared_policy(
    sources: dict[str, Path],
    tmp_path: Path,
    policy: Policy,
    scope: Literal["off", "phenotype"],
    modern: bool,
) -> None:
    matrix, metadata = tmp_path / "matrix.tsv", tmp_path / "metadata.tsv"
    create_taxa_trait_matrix(
        taxa_path=sources["taxa"],
        edges_path=sources["edges"],
        nodes_path=sources["nodes"],
        output_path=matrix,
        column_metadata_path=metadata,
        strain_rollup=scope,
        trait_key_policy=policy,
        include_policy_metadata=modern,
    )
    _, rows = read_table(matrix)
    genus = next(row for row in rows if row["taxon_id"] == "NCBITaxon:10")
    if policy == RAW:
        assert [int(genus[member]) for member in MEMBERS] == (
            [1, 1, 2] if scope == "off" else [3, 1, 2]
        )
    else:
        assert int(genus[CANONICAL]) == (2 if scope == "off" else 3)
        assert MEMBERS[0] not in genus
        assert MEMBERS[2] not in genus
    header, _ = read_table(metadata)
    assert header == ((*OLD_COLUMNS, *POLICY_COLUMNS) if modern else OLD_COLUMNS)
    validate_matrix_sources(
        matrix,
        sources["taxa"],
        sources["edges"],
        sources["nodes"],
        metadata,
        strain_rollup=scope,
        trait_key_policy=policy,
        include_policy_metadata=modern,
    )
    other = UNION if policy == RAW else RAW
    with pytest.raises(ValueError, match=r"(vocabulary|metadata|policy|count)"):
        validate_matrix_sources(
            matrix,
            sources["taxa"],
            sources["edges"],
            sources["nodes"],
            metadata,
            strain_rollup=scope,
            trait_key_policy=other,
            include_policy_metadata=modern,
        )


def test_union_metadata_discloses_declared_members_and_vocabularies_and_labels_union() -> None:
    from mtp.trait_key_policy import display_trait_label

    nodes = {CANONICAL: NodeMetadata("facultatively aerobic", "biolink:PhenotypicQuality")}
    header, rows = column_metadata_rows({"taxon": Counter({CANONICAL: 1})}, nodes)
    row = dict(zip(header, rows[0], strict=True))
    assert row["object_name"] == "facultatively aerobic"
    assert row["column_label"] == f"facultatively aerobic {CANONICAL}"
    assert row["trait_key_policy"] == UNION
    assert json.loads(row["trait_member_ids_json"]) == MEMBERS
    assert (
        display_trait_label(CANONICAL, row["object_name"], row["trait_member_ids_json"])
        == "facultative (aerobe or anaerobe)"
    )
    assert (
        display_trait_label(CANONICAL, row["column_label"], row["trait_member_ids_json"])
        == f"facultative (aerobe or anaerobe) {CANONICAL}"
    )
    mappings = json.loads(row["trait_source_vocabularies_json"])
    assert {
        tuple(mapping[key] for key in ("object_id", "source_id", "source_label"))
        for mapping in mappings
    } == {
        (MEMBERS[0], "infores:madin_etal", "facultative"),
        (MEMBERS[0], "infores:bacdive", "facultative anaerobe"),
        (MEMBERS[1], "infores:madin_etal", "facultative"),
        (MEMBERS[1], "infores:bacdive", "facultative aerobe"),
        (MEMBERS[2], "infores:bactotraits", "Ox_facultative_aerobe_anaerobe"),
    }
    raw_header, raw_rows = column_metadata_rows(
        {"taxon": Counter({CANONICAL: 1})},
        nodes,
        trait_key_policy=RAW,
    )
    raw = dict(zip(raw_header, raw_rows[0], strict=True))
    assert raw["object_name"] == "facultatively aerobic"
    assert raw["column_label"] == f"facultatively aerobic {CANONICAL}"
    assert json.loads(raw["trait_member_ids_json"]) == [CANONICAL]
    assert (
        display_trait_label(CANONICAL, raw["object_name"], raw["trait_member_ids_json"])
        == "facultatively aerobic"
    )
    assert display_trait_label(CANONICAL, raw["object_name"], None) == "facultatively aerobic"
    assert len(json.loads(raw["trait_source_vocabularies_json"])) == 2


def test_union_relabel_changes_only_labels_not_numerical_outputs(
    sources: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify display_trait_label relabels only union-member runs (#1138, #1232, #1235)."""
    import importlib.util

    import pandas as pd

    from mtp.density_adjustment import TABLE_1_TRAITS as DENSITY_TABLE1_TRAITS
    from mtp.density_adjustment import run_density_audit
    from mtp.label_support_audit import (
        AUDIT_TRAITS,
        build_label_support_tables,
        load_kg_provenance_graph,
    )
    from mtp.ontology_audit import run_ontology_redundancy_audit
    from mtp.table1_trait_maps import compute_direction_shifts, write_shift_table
    from mtp.trait_guilds import TABLE1_TRAIT_LABELS, Profile
    from mtp.trait_lists import trait_labels
    from mtp.trait_screen import (
        _report,
        load_column_labels,
        run_screen,
    )

    monkeypatch.setenv("GIT_CONFIG_GLOBAL", "/dev/null")
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")

    raw_label = "facultatively aerobic"
    raw_col_label = f"{raw_label} {CANONICAL}"
    union_label = "facultative (aerobe or anaerobe)"
    union_col_label = f"{union_label} {CANONICAL}"

    # Static defaults keep raw KG label (#1232).
    assert TABLE1_TRAIT_LABELS["facultatively_aerobic"] == f"{raw_label} ({CANONICAL})"
    assert DENSITY_TABLE1_TRAITS[CANONICAL] == raw_label
    assert AUDIT_TRAITS[CANONICAL] == raw_label

    union_matrix_path = tmp_path / "union_matrix.tsv"
    union_metadata_path = tmp_path / "union_metadata.tsv"
    create_taxa_trait_matrix(
        taxa_path=sources["taxa"],
        edges_path=sources["edges"],
        nodes_path=sources["nodes"],
        output_path=union_matrix_path,
        column_metadata_path=union_metadata_path,
        strain_rollup="phenotype",
        trait_key_policy=UNION,
    )
    raw_matrix_path = tmp_path / "raw_matrix.tsv"
    raw_metadata_path = tmp_path / "raw_metadata.tsv"
    create_taxa_trait_matrix(
        taxa_path=sources["taxa"],
        edges_path=sources["edges"],
        nodes_path=sources["nodes"],
        output_path=raw_matrix_path,
        column_metadata_path=raw_metadata_path,
        strain_rollup="phenotype",
        trait_key_policy=RAW,
    )

    union_cols_df = pd.read_csv(union_metadata_path, sep="\t")
    raw_cols_df = pd.read_csv(raw_metadata_path, sep="\t")
    assert trait_labels(union_cols_df, [CANONICAL]) == {CANONICAL: union_label}
    assert trait_labels(raw_cols_df, [CANONICAL]) == {CANONICAL: raw_label}

    loaded_union = load_column_labels(union_metadata_path)
    loaded_raw = load_column_labels(raw_metadata_path)
    assert loaded_union[CANONICAL].object_name == union_label
    assert loaded_union[CANONICAL].label == union_col_label
    assert loaded_raw[CANONICAL].object_name == raw_label
    assert loaded_raw[CANONICAL].label == raw_col_label

    # Compare trait_screen numerical outputs and both Table 1 + candidate sections (#1232).
    screen_profiles = [
        Profile("S1", "E1", "inc", "cond", "feces", {CANONICAL: 4}, 5),
        Profile("S1", "E1", "dec", "cond", "feces", {CANONICAL: 1}, 5),
        Profile("S2", "E1", "inc", "cond", "feces", {CANONICAL: 3}, 4),
        Profile("S2", "E1", "dec", "cond", "feces", {CANONICAL: 1}, 4),
    ]
    res_union = run_screen(screen_profiles, [CANONICAL], labels=loaded_union)
    res_raw = run_screen(screen_profiles, [CANONICAL], labels=loaded_raw)
    row_union, row_raw = res_union.rows[0], res_raw.rows[0]
    assert row_union.trait_label == union_col_label
    assert row_union.object_name == union_label
    assert row_raw.trait_label == raw_col_label
    assert row_raw.object_name == raw_label
    for field in (
        "studies_up",
        "studies_down",
        "studies_tied",
        "studies_scored",
        "share_up",
        "ci_low",
        "ci_high",
        "p_value",
        "by_q_nominal",
        "by_q_effective",
        "test_status",
    ):
        assert getattr(row_union, field) == getattr(row_raw, field)
    md_union = _report(res_union, vintage="fixture")
    md_raw = _report(res_raw, vintage="fixture")
    assert f"| {union_label} ({CANONICAL}) |" in md_union
    assert f"| {union_col_label} | true |" in md_union
    assert raw_label not in md_union
    assert f"| {raw_label} ({CANONICAL}) |" in md_raw
    assert f"| {raw_col_label} | true |" in md_raw
    assert union_label not in md_raw

    # Check table1_trait_maps shift table label and numerical values for union vs raw metadata.
    guild_profiles = [
        Profile("S1", "E1", "inc", "cond", "feces", {"facultatively_aerobic": 4}, 5),
        Profile("S1", "E1", "dec", "cond", "feces", {"facultatively_aerobic": 1}, 5),
        Profile("S2", "E1", "inc", "cond", "feces", {"facultatively_aerobic": 3}, 4),
        Profile("S2", "E1", "dec", "cond", "feces", {"facultatively_aerobic": 1}, 4),
    ]
    shifts = compute_direction_shifts(guild_profiles, trait_order=("facultatively_aerobic",))
    shift_union_tsv, _ = write_shift_table(
        shifts,
        tmp_path / "shifts_union",
        vintage="fixture",
        column_metadata_path=union_metadata_path,
    )
    shift_raw_tsv, _ = write_shift_table(
        shifts,
        tmp_path / "shifts_raw",
        vintage="fixture",
        column_metadata_path=raw_metadata_path,
    )
    _, shift_union_rows = read_table(shift_union_tsv)
    _, shift_raw_rows = read_table(shift_raw_tsv)
    assert shift_union_rows[0]["trait_label"] == f"{union_label} ({CANONICAL})"
    assert shift_raw_rows[0]["trait_label"] == f"{raw_label} ({CANONICAL})"
    assert shift_union_rows[0]["studies_up"] == shift_raw_rows[0]["studies_up"] == "2"

    # Check label_support_audit uses union label only when columns_path records union members.
    graph = load_kg_provenance_graph(sources["edges"], sources["nodes"], strain_rollup="phenotype")
    occ_df = pd.DataFrame(
        [
            {"study_id": "S1", "experiment_id": "E1", "inc_dec": "inc", "taxon_id": "NCBITaxon:1"},
            {"study_id": "S1", "experiment_id": "E1", "inc_dec": "dec", "taxon_id": "NCBITaxon:2"},
        ]
    )
    union_matrix_df = pd.read_csv(union_matrix_path, sep="\t")
    support_union, _, _ = build_label_support_tables(
        graph, union_matrix_df, occ_df, {}, column_metadata_path=union_metadata_path
    )
    support_raw, _, _ = build_label_support_tables(
        graph, union_matrix_df, occ_df, {}, column_metadata_path=raw_metadata_path
    )
    assert (
        next(r for r in support_union if r["trait_key"] == CANONICAL)["trait_label"] == union_label
    )
    assert next(r for r in support_raw if r["trait_key"] == CANONICAL)["trait_label"] == raw_label

    # Check density_adjustment and scripts/audit_trait_density_correlation.py (mutant 1).
    wide_profiles_path = tmp_path / "wide_profiles.tsv"
    _write(
        wide_profiles_path,
        [
            "study_id",
            "experiment_id",
            "inc_dec",
            "taxa_count",
            CANONICAL,
            "METPO:2000103|GO:0030245",
            "METPO:1000668",
            "METPO:1000683",
            "METPO:2000103|pathways:arsenate_reduction",
            "METPO:1000623",
        ],
        [
            ["S1", "E1", "inc", "10", "5", "1", "4", "4", "2", "3"],
            ["S1", "E1", "dec", "10", "1", "4", "2", "2", "1", "1"],
            ["S2", "E1", "inc", "10", "6", "1", "5", "4", "3", "4"],
            ["S2", "E1", "dec", "10", "2", "5", "1", "1", "1", "1"],
            ["S3", "E1", "inc", "10", "4", "2", "3", "3", "2", "2"],
            ["S3", "E1", "dec", "10", "1", "3", "1", "2", "0", "1"],
        ],
    )
    other_rows = [
        [
            "METPO:2000103|GO:0030245",
            "cellulose catabolic process",
            "cellulose catabolic process",
            "biolink:BiologicalProcess",
            '["GO:0030245"]',
        ],
        [
            "METPO:1000668",
            "coccus shaped METPO:1000668",
            "coccus shaped",
            "biolink:PhenotypicQuality",
            '["METPO:1000668"]',
        ],
        [
            "METPO:1000683",
            "sphere shaped METPO:1000683",
            "sphere shaped",
            "biolink:PhenotypicQuality",
            '["METPO:1000683"]',
        ],
        [
            "METPO:2000103|pathways:arsenate_reduction",
            "arsenate reduction",
            "arsenate_reduction",
            "biolink:BiologicalProcess",
            '["pathways:arsenate_reduction"]',
        ],
        [
            "METPO:1000623",
            "moderately halophilic METPO:1000623",
            "moderately halophilic",
            "biolink:PhenotypicQuality",
            '["METPO:1000623"]',
        ],
    ]
    union_full_cols = tmp_path / "union_full_columns.tsv"
    _write(
        union_full_cols,
        ["column_id", "column_label", "object_name", "object_category", "trait_member_ids_json"],
        [
            [
                CANONICAL,
                raw_col_label,
                raw_label,
                "biolink:PhenotypicQuality",
                json.dumps(MEMBERS),
            ],
            *other_rows,
        ],
    )
    raw_full_cols = tmp_path / "raw_full_columns.tsv"
    _write(
        raw_full_cols,
        ["column_id", "column_label", "object_name", "object_category", "trait_member_ids_json"],
        [
            [
                CANONICAL,
                raw_col_label,
                raw_label,
                "biolink:PhenotypicQuality",
                json.dumps([CANONICAL]),
            ],
            *other_rows,
        ],
    )
    audit_script_path = (
        Path(__file__).resolve().parents[1] / "scripts/audit_trait_density_correlation.py"
    )
    audit_spec = importlib.util.spec_from_file_location(
        "audit_trait_density_correlation", audit_script_path
    )
    assert audit_spec is not None
    assert audit_spec.loader is not None
    audit_script = importlib.util.module_from_spec(audit_spec)
    audit_spec.loader.exec_module(audit_script)

    for mode_name, cols_file, expected_lbl in (
        ("union", union_full_cols, union_label),
        ("raw", raw_full_cols, raw_label),
    ):
        out_dir = tmp_path / f"density_{mode_name}"
        density_res = run_density_audit(
            profiles_path=wide_profiles_path,
            columns_path=cols_file,
            output_dir=out_dir,
        )
        _, t1_rows = read_table(density_res["table1_tsv"])
        fac_row = next(row for row in t1_rows if row["trait_id"] == CANONICAL)
        assert fac_row["trait_label"] == expected_lbl
        report_text = density_res["report_txt"].read_text(encoding="utf-8")
        assert f"Test Case Pair: {expected_lbl} x moderately halophilic" in report_text

        ontology_res = run_ontology_redundancy_audit(
            profiles_path=wide_profiles_path,
            columns_path=cols_file,
            output_dir=out_dir,
        )
        consolidated_path = audit_script.write_consolidated_report(
            out_dir, density_res, ontology_res
        )
        consolidated_text = consolidated_path.read_text(encoding="utf-8")
        assert f"Test Case Pair: {expected_lbl} x moderately halophilic" in consolidated_text
        assert (
            f"Regressing out PC1 collapses the spurious {expected_lbl} x moderately halophilic"
            in consolidated_text
        )


@pytest.mark.parametrize("modern", [False, True])
def test_catalog_preserves_matching_metadata_format(
    sources: dict[str, Path], tmp_path: Path, modern: bool
) -> None:
    metadata, catalog = tmp_path / "metadata.tsv", tmp_path / "catalog.tsv"
    create_taxa_trait_matrix(
        taxa_path=sources["taxa"],
        edges_path=sources["edges"],
        nodes_path=sources["nodes"],
        output_path=tmp_path / "matrix.tsv",
        column_metadata_path=metadata,
        trait_key_policy=RAW,
        include_policy_metadata=modern,
    )
    write_feature_catalog(metadata, catalog, trait_key_policy=RAW, include_policy_metadata=modern)
    header, rows = read_table(catalog)
    expected = (*OLD_COLUMNS, *POLICY_COLUMNS) if modern else OLD_COLUMNS
    assert header == (*expected, "feature_class", "annotation_evidence")
    facts = validate_product(
        catalog,
        "feature_catalog",
        context={"feature_metadata": metadata, "kg_nodes": sources["nodes"]},
        trait_key_policy=RAW,
        include_policy_metadata=modern,
    )
    assert facts["trait_key_policy"] == RAW
    assert rows[0]["object_name"]


@pytest.mark.parametrize(
    "damage", ["policy", "members", "vocabulary", "partial", "legacy_downgrade"]
)
def test_new_metadata_cannot_drop_or_forge_policy_disclosure(
    sources: dict[str, Path],
    tmp_path: Path,
    damage: str,
) -> None:
    metadata = tmp_path / "metadata.tsv"
    create_taxa_trait_matrix(
        taxa_path=sources["taxa"],
        edges_path=sources["edges"],
        nodes_path=sources["nodes"],
        output_path=tmp_path / "matrix.tsv",
        column_metadata_path=metadata,
    )
    header, rows = read_table(metadata)
    edited = list(header)
    if damage == "policy":
        rows[0]["trait_key_policy"] = RAW
    elif damage == "members":
        rows[0]["trait_member_ids_json"] = json.dumps([CANONICAL])
    elif damage == "vocabulary":
        rows[0]["trait_source_vocabularies_json"] = "[]"
    elif damage == "partial":
        edited.remove("trait_member_ids_json")
    else:
        edited = list(OLD_COLUMNS)
    _write(metadata, edited, [[row[column] for column in edited] for row in rows])
    with pytest.raises(ValueError, match=r"(metadata|policy|columns|vocabular|member)"):
        validate_product(
            metadata,
            "feature_metadata",
            trait_key_policy=UNION,
            include_policy_metadata=True,
        )


def test_unknown_policy_is_rejected() -> None:
    from mtp.trait_key_policy import parse_trait_key_policy

    for value in ("", "facultative_union_v2", None, 4):
        with pytest.raises(ValueError, match="trait_key_policy"):
            parse_trait_key_policy(value)
