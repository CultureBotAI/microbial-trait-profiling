"""Removing unused guild matrices preserves same-seed scientific outputs (#784)."""

from __future__ import annotations

import importlib.util
import inspect
import json

from mtp.guild_maps import compute_correlation
from mtp.integrated_evidence import guild_records
from mtp.paths import REPO_ROOT


def test_calibration_matches_before_removal_fixture_without_matrix() -> None:
    script = REPO_ROOT / "scripts/adhoc/issue784_calibration_fixture.py"
    spec = importlib.util.spec_from_file_location("issue784_fixture", script)
    assert spec is not None
    assert spec.loader is not None
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    expected = json.loads(
        (REPO_ROOT / "tests/fixtures/calibration/issue784_before_matrix_removal.json").read_text()
    )
    observed = helper.capture()
    assert observed == expected


def test_guild_statistics_have_no_ignored_matrix_argument() -> None:
    assert list(inspect.signature(compute_correlation).parameters) == ["profiles"]
    assert "matrix" not in inspect.signature(guild_records).parameters
