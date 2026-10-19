"""Reading the Table 1 sign tests that the deck and the plates quote (#579)."""

from __future__ import annotations

import itertools

import pytest

from mtp.table1_directions import ALL_STUDIES, CELLULOSE, CODES, ORIENTED, populations


def _rows(population: str, studies: int, passing: set[frozenset[str]]) -> list[dict[str, str]]:
    """Six fractions and 15 log-ratios; the ``passing`` pairs get BY q = 0.01."""
    rows = []
    for trait in CODES:
        rows.append(
            {
                "population": population,
                "comparison": "fraction",
                "numerator": trait,
                "reference": "",
                "population_studies": str(studies),
                "studies_scored": "100",
                "studies_up": "40" if trait == CELLULOSE else "60",
                "studies_down": "60" if trait == CELLULOSE else "40",
                "p_value": "0.001",
                "q_value": "0.002",
            }
        )
    for numerator, reference in itertools.combinations(CODES, 2):
        rows.append(
            {
                "population": population,
                "comparison": "log_ratio",
                "numerator": numerator,
                "reference": reference,
                "population_studies": str(studies),
                "studies_scored": "50",
                "studies_up": "20",
                "studies_down": "30",
                "p_value": "0.2",
                "q_value": "0.01" if frozenset((numerator, reference)) in passing else "0.3",
            }
        )
    return rows


CEL_HAL = frozenset((CELLULOSE, "moderately_halophilic"))
FA_CEL = frozenset(("facultatively_aerobic", CELLULOSE))


def test_populations_read_fractions_and_ratios() -> None:
    table = populations(_rows(ALL_STUDIES, 959, {FA_CEL}) + _rows(ORIENTED, 592, {FA_CEL, CEL_HAL}))
    assert table[ALL_STUDIES].studies == 959
    oriented = table[ORIENTED]
    assert oriented.studies == 592
    cellulose = oriented.fractions[CELLULOSE]
    assert (cellulose.label, cellulose.share_down) == ("CEL", 0.6)
    # A pair is found in either order and labelled in the table's order.
    assert oriented.ratio("moderately_halophilic", CELLULOSE).label == "CEL/HAL"
    assert [test.label for test in oriented.ratio_survivors()] == ["FA/CEL", "CEL/HAL"]
    assert [test.label for test in table[ALL_STUDIES].ratio_survivors()] == ["FA/CEL"]


def test_a_partial_population_is_refused() -> None:
    rows = _rows(ALL_STUDIES, 959, set()) + _rows(ORIENTED, 592, set())[:-1]
    with pytest.raises(ValueError, match="disease_vs_healthy_oriented: expected the six"):
        populations(rows)


def test_both_quoted_populations_are_needed() -> None:
    with pytest.raises(ValueError, match="no disease_vs_healthy_oriented rows"):
        populations(_rows(ALL_STUDIES, 959, set()))


def test_an_untested_row_is_named() -> None:
    rows = _rows(ALL_STUDIES, 959, set()) + _rows(ORIENTED, 592, set())
    rows[-1]["q_value"] = ""
    with pytest.raises(ValueError, match="arsenate_reduction/moderately_halophilic: untested"):
        populations(rows)


def test_populations_must_agree_on_their_size() -> None:
    rows = _rows(ALL_STUDIES, 959, set()) + _rows(ORIENTED, 592, set())
    rows[-1]["population_studies"] = "591"
    with pytest.raises(ValueError, match="rows disagree on the number of studies"):
        populations(rows)
