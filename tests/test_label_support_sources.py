"""Refuse ambiguous label evidence and preserve single-source behavior (#1362)."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from mtp.label_support_audit import build_label_support_tables, run_label_support_audit
from test_label_support_audit import inputs as inputs

FAC = "METPO:1000608"
COCCUS = "METPO:1000668"
DECLARED = "METPO:1000447"
TAXON = "NCBITaxon:816"
UNION_MEMBERS = ("METPO:1000605", FAC, "METPO:1000612")
RAW_LABEL = "facultatively aerobic"
UNION_LABEL = "facultative (aerobe or anaerobe)"


@pytest.fixture
def tables() -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """One Bacteroides node appears in both support and focal output tables."""
    keys = (FAC, COCCUS, DECLARED)
    graph = {
        "traits": {TAXON: set(keys)},
        "sources": {key: {TAXON: {"infores:fixture"}} for key in keys},
        "children": {},
        "names": {
            TAXON: "Bacteroides",
            FAC: RAW_LABEL,
            COCCUS: "coccus shaped",
            DECLARED: "temperature optimum high",
        },
    }
    matrix = pd.DataFrame(
        [
            {
                "taxon_id": TAXON,
                "taxon_name": "Bacteroides",
                "taxonomic_rank": "genus",
                FAC: 1,
                COCCUS: 1,
                DECLARED: 1,
            }
        ]
    )
    occurrences = pd.DataFrame(
        [
            {
                "study_id": "fixture-study",
                "experiment_id": "fixture-experiment",
                "inc_dec": "inc",
                "taxon_id": TAXON,
            }
        ]
    )
    return graph, matrix, occurrences


def _metadata(path: Path, members: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["column_id", "object_id", "trait_member_ids_json"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            {
                "column_id": key,
                "object_id": key,
                "trait_member_ids_json": json.dumps(values),
            }
            for key, values in (
                (FAC, members),
                (COCCUS, (COCCUS,)),
                (DECLARED, (DECLARED,)),
            )
        )


@pytest.mark.parametrize("path_exists", [True, False], ids=["existing-path", "missing-path"])
@pytest.mark.parametrize(
    "members",
    [{}, {FAC: UNION_MEMBERS}, {FAC: (FAC,)}],
    ids=["empty-map", "consistent-map", "conflicting-map"],
)
def test_builder_refuses_two_sources_even_when_the_mapping_is_empty(
    tmp_path: Path,
    tables: tuple[dict[str, Any], pd.DataFrame, pd.DataFrame],
    path_exists: bool,
    members: dict[str, tuple[str, ...]],
) -> None:
    """A mapping must not silently override a supplied path, even an unreadable one."""
    columns = tmp_path / "columns.tsv"
    if path_exists:
        _metadata(columns, UNION_MEMBERS)
    with pytest.raises(
        ValueError,
        match=r"column_metadata_path.*member_ids_by_trait|member_ids_by_trait.*column_metadata_path",
    ):
        build_label_support_tables(
            *tables,
            {},
            trait_list=(DECLARED,),
            column_metadata_path=columns,
            member_ids_by_trait=members,
        )


@pytest.mark.parametrize(
    ("members", "expected_label"),
    [((FAC,), RAW_LABEL), (UNION_MEMBERS, UNION_LABEL)],
    ids=["raw-members", "union-members"],
)
def test_path_and_mapping_sources_produce_identical_tables(
    tmp_path: Path,
    tables: tuple[dict[str, Any], pd.DataFrame, pd.DataFrame],
    members: tuple[str, ...],
    expected_label: str,
) -> None:
    """Equivalent supplied evidence changes only the supported display label."""
    columns = tmp_path / "columns.tsv"
    _metadata(columns, members)
    # This mapping is independently declared, not obtained by invoking the path reader.
    mapping = {FAC: members, COCCUS: (COCCUS,), DECLARED: (DECLARED,)}
    from_path = build_label_support_tables(
        *tables, {}, trait_list=(DECLARED,), column_metadata_path=columns
    )
    from_mapping = build_label_support_tables(
        *tables, {}, trait_list=(DECLARED,), member_ids_by_trait=mapping
    )
    assert from_path == from_mapping
    support, focal, evidence = from_path
    assert {row["trait_key"] for row in support} == {FAC, COCCUS, DECLARED}
    for rows in (support, focal):
        row = next(row for row in rows if row["taxon_id"] == TAXON and row["trait_key"] == FAC)
        assert row["trait_label"] == expected_label
        assert row["carrier_count"] == 1
        assert row["signature_occurrences"] == 1
    assert len(evidence) == 3
    assert {row["carrier_taxon_id"] for row in evidence} == {TAXON}


def test_empty_mapping_alone_preserves_the_no_metadata_fallback(
    tables: tuple[dict[str, Any], pd.DataFrame, pd.DataFrame],
) -> None:
    without_source = build_label_support_tables(*tables, {}, trait_list=(DECLARED,))
    empty_mapping = build_label_support_tables(
        *tables, {}, trait_list=(DECLARED,), member_ids_by_trait={}
    )
    assert empty_mapping == without_source
    labels = {
        row["trait_label"]
        for rows in empty_mapping[:2]
        for row in rows
        if row["taxon_id"] == TAXON and row["trait_key"] == FAC
    }
    assert labels == {RAW_LABEL}


@pytest.mark.parametrize("source", ["missing", "raw", "union"])
def test_public_provenance_names_only_audited_traits_whose_label_changed(
    inputs: dict[str, Path], tmp_path: Path, source: str
) -> None:
    """Available metadata rows are not the count of labels changed by those rows."""
    columns = inputs["matrix"].with_name("taxa_trait_matrix_columns.tsv")
    absent_facultative = "METPO:2000103|METPO:1000608"
    if source != "missing":
        _metadata(columns, UNION_MEMBERS if source == "union" else (FAC,))
        # This absent composite would acquire a union label if it were audited.
        # The other metadata rows include a present coccus and an absent non-focal
        # trait. None may be counted as an actual label change in these outputs.
        with columns.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([absent_facultative, absent_facultative, json.dumps(UNION_MEMBERS)])
    output = tmp_path / "source-provenance"
    record = run_label_support_audit(
        taxa_path=inputs["taxa"],
        matrix_path=inputs["matrix"],
        profiles_path=inputs["profiles"],
        edges_path=inputs["edges"],
        nodes_path=inputs["nodes"],
        output_dir=output,
        strain_rollup="off",
    )
    expected = [FAC] if source == "union" else []
    metadata = record["inputs"]["column_metadata"]
    assert "label_changed_trait_ids" in metadata
    assert metadata["label_changed_trait_ids"] == expected
    assert metadata["member_id_rows"] == (0 if source == "missing" else 4)
    saved = json.loads((output / "provenance.json").read_text(encoding="utf-8"))
    assert saved["inputs"]["column_metadata"]["label_changed_trait_ids"] == expected
    support = pd.read_csv(output / "taxon_trait_label_support.tsv", sep="\t")
    assert set(support.loc[support["trait_key"] == FAC, "trait_label"]) == {
        UNION_LABEL if source == "union" else RAW_LABEL
    }
    assert set(support.loc[support["trait_key"] == COCCUS, "trait_label"]) == {"coccus shaped"}
    assert absent_facultative not in set(support["trait_key"])
