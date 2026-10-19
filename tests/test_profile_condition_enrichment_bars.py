"""Standalone bar charts stay out of the committed profile-condition-bars folder (#690)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_profile_condition_enrichment_bars.py"
)
SPEC = importlib.util.spec_from_file_location("profile_condition_enrichment_bars", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
bars = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = bars
SPEC.loader.exec_module(bars)


def test_standalone_bars_are_refused_in_the_committed_folder(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(sys, "argv", ["bars", "--standalone-bars"])
    with pytest.raises(SystemExit):
        bars.parse_args()


def test_standalone_bars_are_allowed_elsewhere(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(sys, "argv", ["bars", "--standalone-bars", "--output-dir", str(tmp_path)])
    args = bars.parse_args()
    assert args.standalone_bars
    assert args.output_dir == tmp_path
