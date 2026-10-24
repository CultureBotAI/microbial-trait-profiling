"""Removing unused guild matrices preserves same-seed scientific outputs (#784)."""

from __future__ import annotations

import importlib.util
import inspect
import json
from copy import deepcopy
from types import ModuleType
from typing import Any

import pytest

from mtp.guild_maps import compute_correlation
from mtp.integrated_evidence import guild_records
from mtp.paths import REPO_ROOT


def _fixture_helper() -> ModuleType:
    script = REPO_ROOT / "scripts/adhoc/issue784_calibration_fixture.py"
    spec = importlib.util.spec_from_file_location("issue784_fixture", script)
    assert spec is not None
    assert spec.loader is not None
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


def _baseline() -> dict[str, Any]:
    snapshot = json.loads(
        (REPO_ROOT / "tests/fixtures/calibration/issue784_before_matrix_removal.json").read_text()
    )
    assert isinstance(snapshot, dict)
    return snapshot


def _scientific_results(snapshot: dict[str, Any]) -> dict[str, Any]:
    """Keep the ledger, outcomes and numeric/decision criteria, without prose."""
    return {
        **snapshot,
        "criteria": [
            {key: value for key, value in row.items() if key not in {"detail", "requirement"}}
            for row in snapshot["criteria"]
        ],
    }


def test_calibration_matches_before_removal_fixture_without_matrix() -> None:
    assert _scientific_results(_fixture_helper().capture()) == _scientific_results(_baseline())


def test_fixture_capture_has_no_obsolete_legacy_interface() -> None:
    assert not inspect.signature(_fixture_helper().capture).parameters


def test_fixture_comparison_ignores_only_explanatory_prose() -> None:
    baseline = _baseline()
    reworded = deepcopy(baseline)
    for row in reworded["criteria"]:
        row["detail"] = "Reworded explanation."
        row["requirement"] = "Reworded requirement."
    assert _scientific_results(reworded) == _scientific_results(baseline)


@pytest.mark.parametrize(
    ("table", "field", "value"),
    [
        ("ledger", "tests_evaluated", 999),
        ("guild_outcomes", "p_value", 0.5),
        ("criteria", "estimate", 0.5),
        ("criteria", "tolerance", 0.5),
        ("criteria", "verdict", "pass"),
        ("criteria", "procedure", "different_procedure"),
    ],
)
def test_fixture_comparison_retains_scientific_results(
    table: str, field: str, value: object
) -> None:
    baseline = _baseline()
    changed = deepcopy(baseline)
    changed[table][0][field] = value
    assert _scientific_results(changed) != _scientific_results(baseline)


def test_guild_statistics_have_no_ignored_matrix_argument() -> None:
    assert list(inspect.signature(compute_correlation).parameters) == ["profiles"]
    assert "matrix" not in inspect.signature(guild_records).parameters
