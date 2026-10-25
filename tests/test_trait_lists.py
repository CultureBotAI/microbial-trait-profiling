"""Declared trait lists for the six-trait analyses (#1054)."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pytest

from mtp.density_adjustment import TABLE_1_TRAITS
from mtp.trait_lists import load_trait_list, parse_trait_list, trait_labels


def test_trait_list_keeps_declared_order_and_hashes_the_file(tmp_path: Path) -> None:
    text = (
        "# Ten headline traits\n"
        "METPO:1000605\n"
        "\n"
        "  METPO:2000103|pathways:arsenate_reduction  \n"
        "METPO:2000103|GO:0030245\n"
    )
    path = tmp_path / "traits.txt"
    path.write_text(text, encoding="utf-8")

    declared = load_trait_list(path)

    assert declared.trait_ids == (
        "METPO:1000605",
        "METPO:2000103|pathways:arsenate_reduction",
        "METPO:2000103|GO:0030245",
    )
    assert declared.sha256 == hashlib.sha256(text.encode("utf-8")).hexdigest()
    assert declared.provenance() == {
        "trait_list_file": "traits.txt",
        "trait_list_sha256": declared.sha256,
        "trait_count": 3,
        "trait_ids": list(declared.trait_ids),
    }


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("# only a comment\n\n", "declares no traits"),
        ("METPO:1000605\nMETPO:1000605\n", "repeats"),
        ("METPO:1000605\nfacultatively anaerobic\n", "line 2 is not a trait identifier"),
        ("METPO:2000103|\n", "line 1 is not a trait identifier"),
    ],
)
def test_trait_list_rejects_malformed_declarations(text: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        parse_trait_list(text)


def test_trait_labels_come_from_object_name_and_never_from_identifiers() -> None:
    columns = pd.DataFrame(
        {
            "column_id": ["METPO:1000605", "METPO:1000607", "METPO:1000603"],
            "object_name": ["facultatively anaerobic", "  ", "anaerobic"],
        }
    )
    assert trait_labels(columns, ["METPO:1000603", "METPO:1000605"]) == {
        "METPO:1000603": "anaerobic",
        "METPO:1000605": "facultatively anaerobic",
    }
    with pytest.raises(ValueError, match="without an object_name label"):
        trait_labels(columns, ["METPO:1000607"])
    with pytest.raises(ValueError, match=r"absent from the trait columns table: \['METPO:9'\]"):
        trait_labels(columns, ["METPO:1000605", "METPO:9"])
    with pytest.raises(ValueError, match="needs column_id and object_name"):
        trait_labels(columns.drop(columns="object_name"), ["METPO:1000605"])


def test_committed_trait_lists_parse_and_nest() -> None:
    """The #1045 lists: ten headline traits, and those plus Table 1 and METPO:1000607."""
    lists = Path(__file__).resolve().parents[1] / "methods" / "trait_lists"
    headline = load_trait_list(lists / "headline_ten_phenotype.txt").trait_ids
    density = load_trait_list(lists / "density_headline_table1_phenotype.txt").trait_ids

    assert len(headline) == 10
    assert density[: len(headline)] == headline
    assert set(density) == set(headline) | set(TABLE_1_TRAITS) | {"METPO:1000607"}
