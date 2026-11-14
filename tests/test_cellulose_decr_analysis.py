"""Tests for the unified literal stool predicate across analyses (#1157)."""

from __future__ import annotations

import pandas as pd
import pytest

import mtp.cellulose_decr.analysis as cellulose_analysis
import mtp.delta_z_derived.tables as delta_z_tables
import mtp.density_adjustment as density_adjustment
import mtp.matrix_overlap_coverage as matrix_overlap_coverage
from mtp.cellulose_decr.analysis import (
    FECES,
    NON_FECAL_GI,
    ORAL,
    OTHER,
    compartment_of,
    is_stool_body_site,
    validate_stool_body_site,
)


def test_shared_stool_predicate_reexported_identically_across_modules() -> None:
    """All four analysis modules export the same literal stool predicate (#1157)."""
    assert delta_z_tables.STOOL_BODY_SITE == "Feces"
    assert cellulose_analysis.STOOL_BODY_SITE == "Feces"
    assert matrix_overlap_coverage.STOOL_BODY_SITE == "Feces"
    assert density_adjustment.STOOL_BODY_SITE == "Feces"
    assert delta_z_tables.is_stool_body_site is is_stool_body_site
    assert cellulose_analysis.is_stool_body_site is is_stool_body_site
    assert matrix_overlap_coverage.is_stool_body_site is is_stool_body_site
    assert density_adjustment.is_stool_body_site is is_stool_body_site
    assert delta_z_tables.validate_stool_body_site is validate_stool_body_site
    assert cellulose_analysis.validate_stool_body_site is validate_stool_body_site
    assert matrix_overlap_coverage.validate_stool_body_site is validate_stool_body_site
    assert density_adjustment.validate_stool_body_site is validate_stool_body_site


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("Feces", True),
        ("feces", False),
        ("FECES", False),
        (" Feces ", False),
        ("Feces,Blood", False),
        ("Rectal,Feces", False),
        ("Oral", False),
        ("", False),
        (None, False),
        (123, False),
    ],
)
def test_is_stool_body_site_accepts_only_exact_literal_feces(value: object, expected: bool) -> None:
    """is_stool_body_site returns True only for the exact literal 'Feces' (#1157)."""
    assert is_stool_body_site(value) is expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("Feces", True),
        ("Feces,Blood", False),
        ("Rectal,Feces", False),
        ("Oral", False),
        ("Colon", False),
        ("", False),
        (None, False),
        (123, False),
    ],
)
def test_validate_stool_body_site_accepts_literal_and_non_stool_sites(
    value: object, expected: bool
) -> None:
    """validate_stool_body_site returns True for 'Feces' and False for non-stool sites (#1157)."""
    assert validate_stool_body_site(value) is expected


@pytest.mark.parametrize(
    "variant",
    [
        "feces",
        "FECES",
        " Feces",
        "Feces ",
        " Feces ",
        "\tFeces\n",
        " feces ",
    ],
)
def test_validate_stool_body_site_refuses_case_and_whitespace_variants(variant: str) -> None:
    """validate_stool_body_site refuses case and whitespace variants of 'Feces' (#1157)."""
    message = ""
    try:
        validate_stool_body_site(variant)
    except ValueError as exc:
        message = str(exc)
    assert "Feces" in message, f"Expected ValueError for {variant!r}"


@pytest.mark.parametrize(
    ("site", "compartment"),
    [
        ("Feces", FECES),
        ("Colon", NON_FECAL_GI),
        ("Rectal", NON_FECAL_GI),
        ("Sigmoid colon", NON_FECAL_GI),
        ("Terminal ileum", NON_FECAL_GI),
        ("Tongue", ORAL),
        ("Supragingival dental plaque", ORAL),
        ("Oropharynx", ORAL),
        ("Blood", OTHER),
        ("Vagina", OTHER),
        ("Colon,Lung", OTHER),
        ("Rectal,Feces", OTHER),
        ("Feces,Blood", OTHER),
        ("Skin,Nasal cavity", OTHER),
    ],
)
def test_compartment_of_accepts_literal_feces_and_non_stool_sites(
    site: str, compartment: str
) -> None:
    """compartment_of assigns 'Feces' to FECES and compound/non-stool sites without folding."""
    assert compartment_of(site) == compartment


@pytest.mark.parametrize(
    "variant",
    [
        "feces",
        "FECES",
        " Feces",
        "Feces ",
        " Feces ",
        " feces ",
    ],
)
def test_compartment_of_refuses_case_and_whitespace_feces_variants(variant: str) -> None:
    """compartment_of refuses lowercase or whitespace-padded 'Feces' variants (#1157)."""
    message = ""
    try:
        compartment_of(variant)
    except ValueError as exc:
        message = str(exc)
    assert "Feces" in message, f"Expected ValueError for {variant!r}"


def test_delta_z_tables_use_shared_is_stool_body_site() -> None:
    """delta_z_derived tables classify only exact literal 'Feces' as stool (#1157)."""
    series = pd.Series(["Feces", "feces", "Feces,Blood", "Oral", ""])
    mapped = series.map(delta_z_tables.is_stool_body_site).tolist()
    assert mapped == [True, False, False, False, False]
