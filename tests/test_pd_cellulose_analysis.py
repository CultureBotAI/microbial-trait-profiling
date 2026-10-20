"""The PD cellulose analysis refuses stale six-shared-trait Sankey edge tables (#687)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "create_pd_cellulose_analysis.py"
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("pd_cellulose_analysis", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
cellulose = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = cellulose
SPEC.loader.exec_module(cellulose)


def _display(increased: tuple[int, int], decreased: tuple[int, int]) -> pd.DataFrame:
    rows = [
        (direction, cellulose.CELLULOSE_LABEL, cellulose.FECES, profiles, studies)
        for direction, (profiles, studies) in [("increased", increased), ("decreased", decreased)]
    ]
    rows.append(("decreased", "coccus shaped", cellulose.FECES, 99, 99))
    return pd.DataFrame(
        rows,
        columns=["direction", "trait", "body_site", "focus_enriched_profiles", "focus_studies"],
    )


def _recounted(increased: tuple[int, int], decreased: tuple[int, int]) -> pd.DataFrame:
    return pd.DataFrame(
        [("increased", *increased), ("decreased", *decreased)],
        columns=["direction", "focus_enriched_profiles", "focus_studies"],
    )


def test_matching_edges_pass() -> None:
    cellulose.check_display_edges_current(
        _display((34, 22), (69, 34)), _recounted((34, 22), (69, 34))
    )


@pytest.mark.parametrize(
    "recounted",
    [_recounted((34, 22), (70, 34)), _recounted((34, 21), (69, 34))],
)
def test_stale_edges_are_refused_with_the_recipe_to_run(recounted: pd.DataFrame) -> None:
    with pytest.raises(ValueError, match=r"just six-shared-trait-sankeys"):
        cellulose.check_display_edges_current(_display((34, 22), (69, 34)), recounted)


def test_a_missing_edge_is_refused() -> None:
    display = _display((34, 22), (69, 34))
    display = display[~display["direction"].eq("increased")]
    with pytest.raises(ValueError, match="stale"):
        cellulose.check_display_edges_current(display, _recounted((34, 22), (69, 34)))
