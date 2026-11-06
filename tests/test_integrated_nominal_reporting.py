"""Integrated Records tables label nominal p and count displayed tests (#1430, R3).

Small synthetic records exercise the real Markdown and TSV writers without
loading any pipeline inputs or changing correction families and decisions.
"""

from __future__ import annotations

import csv
import re
from dataclasses import fields, replace
from pathlib import Path

import pytest

from mtp.integrated_evidence import (
    DRAWN,
    NOT_SIGNIFICANT,
    PROFILE_SUMMARY,
    TESTED,
    UNDECIDED,
    Record,
    placeholder_records,
    summarize,
    write_table,
    write_view,
)


def _record(trait: str, p_value: float | None) -> Record:
    """A displayed method row whose source family is larger than its exported subset."""
    base = placeholder_records("nominal_fixture", "study", "bounded report fixture", [trait])[0]
    return replace(
        base,
        method_version="fixture-v1",
        effect_value=-0.123456,
        effect_units="fraction difference",
        p_value=p_value,
        q_value=0.05004,
        correction_method="BH",
        family_id="fixture-source-family",
        family_size=1000,
        family_scope="source_subset",
        family_source="fixture.tsv",
        q_threshold=0.05,
        test_status=TESTED,
        decision=NOT_SIGNIFICANT,
        testable_units=73,
        paired_experiments=97,
        native_pair_p_value=0.0002,
        n_tests=701,
    )


def _cell() -> Record:
    return replace(
        _record("hidden-cell", 0.03125),
        method="oxygen_guild_cell_support",
        support_threshold=3,
        display_status=DRAWN,
        condition="fixture condition",
        body_site="fixture site",
    )


def _ora_summary() -> Record:
    return replace(
        _record("ora-summary", None),
        method="enrichment_arm_profile_hypergeometric",
        test_status=PROFILE_SUMMARY,
        decision=UNDECIDED,
        p_value=None,
        q_value=None,
        effect_value=None,
        n_tests=701,
        n_significant_inc=11,
        n_significant_dec=7,
        inc_dec_ratio=11 / 7,
        best_profile_q=0.0001,
        median_profile_q=0.0004,
        median_hit_fold=2.0,
        n_depleted_inc=2,
        n_depleted_dec=3,
    )


def _records_section(view: str) -> str:
    return view.split("## Records\n", 1)[1].split("\n## ", 1)[0]


def _table_rows(view: str) -> list[list[str]]:
    lines = [line for line in _records_section(view).splitlines() if line.startswith("| ")]
    return [[cell.strip() for cell in line.strip("|").split("|")] for line in lines]


def _assert_count(view: str, expected: int) -> None:
    statement = (
        f"Finite nominal tests shown in the Records table: {expected}. "
        "This displayed-row count is not a correction family size."
    )
    assert statement in _records_section(view)
    assert view.count(statement) == 1


def test_records_table_counts_only_finite_primary_p_in_displayed_rows(tmp_path: Path) -> None:
    records = [
        _record("finite", 0.01234567),
        replace(_record("zero", 0.0), effect_value=None, q_value=None, decision=UNDECIDED),
        _record("missing", None),
        _record("nan", float("nan")),
        _record("positive-infinity", float("inf")),
        _record("negative-infinity", -float("inf")),
        replace(
            _record("native-subset", 0.5),
            method="composition_arm_native_transformed_contrast",
            unique_studies=43,
            n_zero_total_experiments=2,
            n_nonzero_experiments=95,
        ),
        _cell(),
        _ora_summary(),
    ]
    table = write_table(records, tmp_path / "before.tsv")
    before = table.read_bytes()
    assert before
    path = write_view(records, summarize(records), tmp_path / "evidence.md")
    view = path.read_text(encoding="utf-8")
    _assert_count(view, 3)
    rows = _table_rows(view)
    assert rows[0] == [
        "Method",
        "Feature",
        "Level",
        "Effect",
        "Units",
        "nominal p",
        "q (family)",
        "Decision",
        "Test state",
    ]
    data = rows[2:]
    assert [row[1] for row in data] == [
        "finite",
        "zero",
        "missing",
        "nan",
        "positive-infinity",
        "negative-infinity",
        "native-subset",
        "ora-summary",
    ]
    assert [row[5] for row in data] == ["0.0123", "0", "NA", "nan", "inf", "-inf", "0.5", "NA"]
    assert data[0][3] == "-0.1235"
    assert data[0][6] == "0.05 (fixture-source-family)"
    assert data[0][7:] == [NOT_SIGNIFICANT, TESTED]
    assert data[1][3] == "NA"
    assert data[1][6] == "NA (fixture-source-family)"
    assert "hidden-cell" not in _records_section(view)
    assert "BH over 1000 support/variance-eligible traits" in view
    assert "| ora-summary | 701 | 11 | 7 |" in view
    # Presentation must not alter source values, selection, decisions or TSV schema.
    assert table.read_bytes() == before
    assert write_table(records, tmp_path / "after.tsv").read_bytes() == before
    with table.open(encoding="utf-8", newline="") as handle:
        exported = list(csv.DictReader(handle, delimiter="\t"))
    assert list(exported[0]) == [field.name for field in fields(Record)]
    assert [row["trait"] for row in exported] == [record.trait for record in records]
    assert exported[0]["p_value"] == "0.0123457"
    assert exported[0]["q_value"] == "0.05004"
    assert exported[0]["decision"] == NOT_SIGNIFICANT
    assert exported[0]["q_threshold"] == "0.05"
    assert exported[0]["correction_method"] == "BH"
    assert exported[0]["family_size"] == "1000"
    assert exported[0]["native_pair_p_value"] == "0.0002"


@pytest.mark.parametrize("case", ["empty", "unavailable", "cells_only", "ora_only"])
def test_records_table_reports_zero_when_no_finite_primary_tests_are_displayed(
    tmp_path: Path, case: str
) -> None:
    records: list[Record]
    if case == "empty":
        records = []
    elif case == "unavailable":
        records = [
            _record("missing", None),
            _record("nan", float("nan")),
            _record("positive-infinity", float("inf")),
            _record("negative-infinity", -float("inf")),
        ]
    elif case == "cells_only":
        records = [_cell()]
    else:
        records = [_ora_summary()]
    path = write_view(records, summarize(records), tmp_path / "evidence.md")
    view = path.read_text(encoding="utf-8")
    _assert_count(view, 0)
    assert _table_rows(view)[0][5] == "nominal p"


def test_records_table_counts_zero_and_one_without_requiring_q_values(tmp_path: Path) -> None:
    records = [
        replace(_record(f"trait-{index}", p), q_value=None, decision=UNDECIDED)
        for index, p in enumerate([0.0, 1.0, 0.05, 0.99])
    ]
    path = write_view(records, summarize(records), tmp_path / "evidence.md")
    view = path.read_text(encoding="utf-8")
    _assert_count(view, 4)
    assert [row[5] for row in _table_rows(view)[2:]] == ["0", "1", "0.05", "0.99"]


def test_rendered_dictionary_distinguishes_nominal_p_from_corrected_q(tmp_path: Path) -> None:
    path = write_view([], [], tmp_path / "evidence.md")
    view = path.read_text(encoding="utf-8")
    dictionary = view.split("## Data dictionary\n", 1)[1].split("\n## ", 1)[0]
    line = next(
        line
        for line in dictionary.splitlines()
        if "`p_value / q_value / correction_method / family_id / family_size`" in line
    )
    meaning = line.strip("|").split("|")[1]
    assert re.search(r"\bp_value\b[^;.]*\bnominal\b", meaning)
    assert re.search(r"\bq_value\b[^;.]*\bcorrected\b", meaning)
