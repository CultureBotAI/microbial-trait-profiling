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
    nodes = {CANONICAL: NodeMetadata("facultatively aerobic", "biolink:PhenotypicQuality")}
    header, rows = column_metadata_rows({"taxon": Counter({CANONICAL: 1})}, nodes)
    row = dict(zip(header, rows[0], strict=True))
    assert row["object_name"] == "facultative (aerobe or anaerobe)"
    assert row["column_label"] == f"facultative (aerobe or anaerobe) {CANONICAL}"
    assert row["trait_key_policy"] == UNION
    assert json.loads(row["trait_member_ids_json"]) == MEMBERS
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
    assert len(json.loads(raw["trait_source_vocabularies_json"])) == 2


def test_union_relabel_changes_only_labels_not_numerical_outputs(
    sources: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify on a fixture that relabelling METPO:1000608 changes labels only (#1138)."""
    from types import SimpleNamespace
    from typing import cast

    import pandas as pd

    monkeypatch.setenv("GIT_CONFIG_GLOBAL", "/dev/null")
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")

    from mtp.density_adjustment import TABLE_1_TRAITS as DENSITY_TABLE1_TRAITS
    from mtp.density_adjustment import run_density_audit
    from mtp.label_support_audit import AUDIT_TRAITS
    from mtp.ontology_audit import RedundancyCluster, build_redundant_pairs_table
    from mtp.publication_figures.inputs import _labels
    from mtp.table1_trait_maps import compute_direction_shifts, write_shift_table
    from mtp.three_test_table.reader import CompletedTable
    from mtp.trait_guilds import TABLE1_TRAIT_LABELS, Profile
    from mtp.trait_lists import trait_labels
    from mtp.trait_screen import (
        ColumnLabel,
        _report,
        load_column_labels,
        run_screen,
    )

    union_label = "facultative (aerobe or anaerobe)"
    union_col_label = f"{union_label} {CANONICAL}"
    assert TABLE1_TRAIT_LABELS["facultatively_aerobic"] == f"{union_label} ({CANONICAL})"
    assert DENSITY_TABLE1_TRAITS[CANONICAL] == union_label
    assert AUDIT_TRAITS[CANONICAL] == union_label

    matrix_path = tmp_path / "matrix.tsv"
    metadata_path = tmp_path / "metadata.tsv"
    create_taxa_trait_matrix(
        taxa_path=sources["taxa"],
        edges_path=sources["edges"],
        nodes_path=sources["nodes"],
        output_path=matrix_path,
        column_metadata_path=metadata_path,
        strain_rollup="phenotype",
        trait_key_policy=UNION,
    )
    cols_df = pd.read_csv(metadata_path, sep="\t")
    assert trait_labels(cols_df, [CANONICAL]) == {CANONICAL: union_label}
    loaded_labels = load_column_labels(metadata_path)
    assert loaded_labels[CANONICAL].object_name == union_label
    assert loaded_labels[CANONICAL].label == union_col_label

    # Compare trait_screen numerical outputs between old and new labels on the same fixture.
    screen_profiles = [
        Profile("S1", "E1", "inc", "cond", "feces", {CANONICAL: 4}, 5),
        Profile("S1", "E1", "dec", "cond", "feces", {CANONICAL: 1}, 5),
        Profile("S2", "E1", "inc", "cond", "feces", {CANONICAL: 3}, 4),
        Profile("S2", "E1", "dec", "cond", "feces", {CANONICAL: 1}, 4),
    ]
    legacy_labels = {
        CANONICAL: ColumnLabel(
            label=f"facultatively aerobic {CANONICAL}",
            object_name="facultatively aerobic",
            object_category="biolink:PhenotypicQuality",
        )
    }
    res_new = run_screen(screen_profiles, [CANONICAL], labels=loaded_labels)
    res_old = run_screen(screen_profiles, [CANONICAL], labels=legacy_labels)
    row_new, row_old = res_new.rows[0], res_old.rows[0]
    assert row_new.trait_label == union_col_label
    assert row_new.object_name == union_label
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
        assert getattr(row_new, field) == getattr(row_old, field)
    md = _report(res_new, vintage="fixture")
    assert f"| {union_label} ({CANONICAL}) |" in md

    # Check table1_trait_maps shift table label and numerical values.
    guild_profiles = [
        Profile("S1", "E1", "inc", "cond", "feces", {"facultatively_aerobic": 4}, 5),
        Profile("S1", "E1", "dec", "cond", "feces", {"facultatively_aerobic": 1}, 5),
        Profile("S2", "E1", "inc", "cond", "feces", {"facultatively_aerobic": 3}, 4),
        Profile("S2", "E1", "dec", "cond", "feces", {"facultatively_aerobic": 1}, 4),
    ]
    shifts = compute_direction_shifts(guild_profiles, trait_order=("facultatively_aerobic",))
    shift_out_dir = tmp_path / "shifts_out"
    shift_tsv, _ = write_shift_table(shifts, shift_out_dir, vintage="fixture")
    _, shift_rows = read_table(shift_tsv)
    assert shift_rows[0]["trait_label"] == f"{union_label} ({CANONICAL})"
    assert shift_rows[0]["studies_up"] == "2"
    assert shift_rows[0]["studies_down"] == "0"

    # Check ontology_audit redundant pairs table uses updated column_label without changing r.
    profiles_df = pd.DataFrame(
        {
            CANONICAL: [4, 1, 3, 1],
            "METPO:1000623": [2, 1, 2, 0],
        }
    )
    extra_cols_path = tmp_path / "extra_columns.tsv"
    _write(
        extra_cols_path,
        ["column_id", "column_label", "object_name"],
        [
            [CANONICAL, union_col_label, union_label],
            ["METPO:1000623", "moderately halophilic METPO:1000623", "moderately halophilic"],
        ],
    )
    clusters = [
        RedundancyCluster(
            cluster_id="C1",
            cluster_label="test",
            tier="lexical_synonym",
            surviving_trait=CANONICAL,
            surviving_label=union_label,
            collapsed_traits=["METPO:1000623"],
            selection_rule="max_support",
            rationale="fixture",
        )
    ]
    pairs_df = build_redundant_pairs_table(profiles_df, extra_cols_path, clusters)
    assert pairs_df["label_1"].tolist() == [union_col_label]
    assert pairs_df["support_trait_1"].tolist() == [9]

    # Check publication_figures.inputs._labels reads the updated object_name from feature_metadata.
    parent = cast(
        CompletedTable,
        SimpleNamespace(
            record={"retained_traits": [CANONICAL]},
            evidence=pd.DataFrame({"trait": [CANONICAL], "trait_label": [union_label]}),
        ),
    )
    assert _labels(cols_df, parent) == {CANONICAL: union_label}

    # Check density_adjustment report and Table 1 adjusted effects table.
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
    full_cols_path = tmp_path / "full_columns.tsv"
    _write(
        full_cols_path,
        ["column_id", "column_label", "object_name", "object_category"],
        [
            [CANONICAL, union_col_label, union_label, "biolink:PhenotypicQuality"],
            [
                "METPO:2000103|GO:0030245",
                "cellulose catabolic process",
                "cellulose catabolic process",
                "biolink:BiologicalProcess",
            ],
            [
                "METPO:1000668",
                "coccus shaped METPO:1000668",
                "coccus shaped",
                "biolink:PhenotypicQuality",
            ],
            [
                "METPO:1000683",
                "sphere shaped METPO:1000683",
                "sphere shaped",
                "biolink:PhenotypicQuality",
            ],
            [
                "METPO:2000103|pathways:arsenate_reduction",
                "arsenate reduction",
                "arsenate_reduction",
                "biolink:BiologicalProcess",
            ],
            [
                "METPO:1000623",
                "moderately halophilic METPO:1000623",
                "moderately halophilic",
                "biolink:PhenotypicQuality",
            ],
        ],
    )
    density_out = tmp_path / "density_out"
    density_res = run_density_audit(
        profiles_path=wide_profiles_path,
        columns_path=full_cols_path,
        output_dir=density_out,
    )
    _, t1_rows = read_table(density_res["table1_tsv"])
    fac_row = next(row for row in t1_rows if row["trait_id"] == CANONICAL)
    assert fac_row["trait_label"] == union_label
    report_text = density_res["report_txt"].read_text(encoding="utf-8")
    assert "Test Case Pair: facultative (aerobe or anaerobe) x moderately halophilic" in report_text


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
