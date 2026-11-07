"""Declare both exported p columns without changing their TSV bytes (#1662)."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path

import mtp.integrated_evidence as evidence
from mtp.integrated_provenance import verify_integrated_provenance

FULL_TABLE = "integrated_evidence.tsv"
METHOD_TABLE = "integrated_evidence_method_results.tsv"
VALUE_POLICY = {
    "finite": "nominal",
    "missing_or_nonfinite": "no_finite_value",
    "correction": "none",
}


def _record(trait: str, primary: float | None, secondary: float | None) -> evidence.Record:
    base = evidence.placeholder_records("p_declaration_fixture", "study", "synthetic", [trait])[0]
    return replace(
        base,
        p_value=primary,
        native_pair_p_value=secondary,
        # These columns deliberately have different availability from both p columns.
        q_value=None if trait == "A" else 0.5,
        native_pair_q_value=0.7 if trait == "A" else None,
        test_status=evidence.TESTED,
        decision=evidence.UNDECIDED,
    )


def _declare(
    records: Sequence[evidence.Record], method_records: Sequence[evidence.Record]
) -> dict[str, object]:
    helper = getattr(evidence, "nominal_p_columns", None)
    assert callable(helper), "The producer must expose its nominal-p column declaration"
    result = helper(records, method_records)
    assert isinstance(result, dict)
    return result


def _expected_table(
    row_count: int, primary: tuple[int, int], secondary: tuple[int, int]
) -> dict[str, object]:
    """Only package literal expected counts; never classify fixture records here."""
    return {
        "row_count": row_count,
        "columns": {
            "p_value": {
                "finite_nominal_rows": primary[0],
                "missing_or_nonfinite_rows": primary[1],
            },
            "native_pair_p_value": {
                "finite_nominal_rows": secondary[0],
                "missing_or_nonfinite_rows": secondary[1],
            },
        },
    }


def test_nominal_p_columns_have_independent_value_and_export_scopes(tmp_path: Path) -> None:
    ordinary = [
        _record("A", 0.0, None),
        _record("B", 1.0, 0.0),
        _record("C", 0.25, None),
        replace(_record("D", None, 0.5), test_status=evidence.INSUFFICIENT),
        _record("E", float("nan"), None),
        _record("F", float("inf"), float("nan")),
        _record("G", -float("inf"), float("inf")),
        _record("H", 0.75, 1.0),
        _record("I", None, -float("inf")),
    ]
    cell = replace(
        _record("cell", 0.9, None),
        method="oxygen_guild_cell_support",
        support_threshold=3,
        display_status=evidence.DRAWN,
    )
    ora = replace(
        _record("ORA-summary", None, None),
        method="enrichment_arm_profile_hypergeometric",
        test_status=evidence.PROFILE_SUMMARY,
        n_tests=701,
    )
    records = [*ordinary, cell, ora]
    method_records = [*ordinary, ora]
    full = evidence.write_table(records, tmp_path / FULL_TABLE)
    methods = evidence.write_table(method_records, tmp_path / METHOD_TABLE)
    before_full, before_methods = full.read_bytes(), methods.read_bytes()

    assert _declare(records, method_records) == {
        "version": "nominal-p-columns-v1",
        "value_policy": VALUE_POLICY,
        "tables": {
            FULL_TABLE: _expected_table(11, (5, 6), (3, 8)),
            METHOD_TABLE: _expected_table(10, (4, 6), (3, 7)),
        },
    }
    # Zero is finite even without q; primary test state cannot classify native pair p.
    assert records[0].p_value == 0.0
    assert records[0].q_value is None
    assert records[3].p_value is None
    assert records[3].native_pair_p_value == 0.5
    assert full.read_bytes() == before_full
    assert methods.read_bytes() == before_methods
    assert evidence.write_table(records, tmp_path / "full-after.tsv").read_bytes() == before_full
    assert (
        evidence.write_table(method_records, tmp_path / "methods-after.tsv").read_bytes()
        == before_methods
    )


def test_nominal_p_columns_use_final_untested_values_without_erasing_secondary_p() -> None:
    original = _record("replaced", 0.25, 0.5)
    final = evidence.untested(replace(original, test_status=evidence.INSUFFICIENT))
    assert original.p_value == 0.25
    assert final.p_value is None
    assert final.native_pair_p_value == 0.5
    assert _declare([final], [final]) == {
        "version": "nominal-p-columns-v1",
        "value_policy": VALUE_POLICY,
        "tables": {
            FULL_TABLE: _expected_table(1, (0, 1), (1, 0)),
            METHOD_TABLE: _expected_table(1, (0, 1), (1, 0)),
        },
    }


def test_empty_exports_keep_both_p_columns_and_their_policy() -> None:
    assert _declare([], []) == {
        "version": "nominal-p-columns-v1",
        "value_policy": VALUE_POLICY,
        "tables": {
            FULL_TABLE: _expected_table(0, (0, 0), (0, 0)),
            METHOD_TABLE: _expected_table(0, (0, 0), (0, 0)),
        },
    }


def test_dictionary_explains_both_nominal_columns_and_the_missing_value_policy() -> None:
    meanings = dict(evidence.DICTIONARY)
    primary = meanings["p_value / q_value / correction_method / family_id / family_size"]
    secondary = meanings["statistical_test / native_pair_p_value / native_pair_q_value"]
    assert re.search(r"\bp_value\b[^;.]*\bnominal\b", primary)
    assert re.search(r"\bq_value\b[^;.]*\bcorrected\b", primary)
    assert re.search(r"\bnative_pair_p_value\b[^;.]*\bnominal\b", secondary)
    assert re.search(r"\bnative_pair_q_value\b[^;.]*\bcorrected\b", secondary)
    policy = meanings.get("nominal_p_columns (provenance)")
    assert policy is not None, "The provenance declaration needs its own dictionary entry"
    assert re.search(r"\bno_finite_value\b[^.]*\bno finite\b", policy)
    assert "export" in policy
    assert re.search(r"\bnot\b[^.]*\b(?:why|reason)\b", policy)


def test_nominal_p_metadata_does_not_become_a_provenance_source(tmp_path: Path) -> None:
    rows = [_record("primary", 0.0, None), _record("secondary", None, 0.5)]
    declaration = _declare(rows, rows)
    pending: list[object] = [declaration]
    while pending:
        value = pending.pop()
        if isinstance(value, Mapping):
            assert {"path", "sha256", "producer_relative_path"}.isdisjoint(value), (
                "Nominal-p metadata must not use keys that identify provenance sources"
            )
            pending.extend(value.values())
        elif isinstance(value, list):
            pending.extend(value)

    source = tmp_path / "source.tsv"
    source_bytes = b"fixture\tvalue\nkept\t1\n"
    source.write_bytes(source_bytes)
    base_record = {
        "version": "integrated-evidence-v2",
        "native_sources": {
            "fixture": {
                "path": source.name,
                "sha256": hashlib.sha256(source_bytes).hexdigest(),
            }
        },
        "ora_eligibility_source": None,
        "assessment_status": {},
    }
    path = tmp_path / "provenance.json"
    expected_inventory = {
        "verified_reference_count": 1,
        "verified_references": ["native_sources.fixture"],
        "skipped_reference_count": 0,
        "skipped_references": [],
    }
    for record in (base_record, {**base_record, "nominal_p_columns": declaration}):
        path.write_text(json.dumps(record), encoding="utf-8")
        try:
            verified = verify_integrated_provenance(path)
        except ValueError as exc:
            raise AssertionError("Count-only metadata must preserve source verification") from exc
        assert {key: verified[key] for key in expected_inventory} == expected_inventory
