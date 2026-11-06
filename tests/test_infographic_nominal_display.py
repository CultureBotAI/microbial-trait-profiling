"""Pin nominal labels and the exact displayed statistic inventory for #1430."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_trait_profile_infographic.py"
)
SPEC = importlib.util.spec_from_file_location("infographic_nominal_display", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
infographic = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(infographic)


def _stats(study: bool) -> pd.DataFrame:
    """Supplied, displayed, unavailable, and alternate-unit counts differ deliberately."""
    frame = pd.DataFrame(
        {
            "trait_label": [f"Trait {index}" for index in range(6)],
            "p_value": [0.1, 0.2, 0.3, np.nan, 0.4, 0.5],
            "minus_log10_p": -np.log10([0.1, 0.2, 0.3, np.nan, 0.4, 0.5]),
            "mean_increased_minus_decreased": [0.9, -0.6, 0.2, -0.1, 0.4, -0.3],
            "ci95_lower": [0.8, -0.7, 0.1, -0.2, 0.3, -0.4],
            "ci95_upper": [1.0, -0.5, 0.3, 0.0, 0.5, -0.2],
            "q_value": [0.2, 0.4, 0.6, np.nan, 0.8, 1.0],
        }
    )
    if study:
        frame["study_p_value"] = [0.0001, 0.001, np.nan, 0.01, np.inf, 0.1]
        frame["study_minus_log10_p"] = [4.0, 3.0, np.nan, 2.0, np.nan, 1.0]
        frame["study_mean_increased_minus_decreased"] = [0.5, -0.3, 0.2, -0.1, 0.4, -0.3]
        frame["study_ci95_lower"] = [0.4, -0.4, 0.1, -0.2, 0.3, -0.4]
        frame["study_ci95_upper"] = [0.6, -0.2, 0.3, 0.0, 0.5, -0.2]
        frame["study_q_value"] = [0.0006, 0.003, np.nan, 0.03, np.nan, 0.2]
    return frame


def _panel(row: int, column: int) -> tuple[Figure, Axes]:
    figure = Figure(figsize=infographic.FIGURE_SIZE, dpi=80)
    FigureCanvasAgg(figure)
    grid = figure.add_gridspec(2, 2, **infographic.PANEL_GRID)
    return figure, figure.add_subplot(grid[row, column])


@pytest.mark.parametrize("study", [True, False], ids=["study", "contrast-fallback"])
def test_top_trait_axis_labels_nominal_p_and_counts_the_active_source(
    study: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(infographic, "TOP_N_TRAITS", 3)
    frame = _stats(study)
    original = frame.copy(deep=True)
    _, axis = _panel(0, 0)
    infographic.plot_top_traits(axis, frame)
    unit = "study-level" if study else "contrast-level"
    assert f"-log10({unit} p; nominal)" in axis.get_xlabel()
    shown, supplied = (2, 4) if study else (3, 5)
    assert f"{shown} finite tests shown / {supplied} finite tests supplied" in axis.get_xlabel()
    score = "study_minus_log10_p" if study else "minus_log10_p"
    expected = frame.head(3).iloc[::-1]
    widths = []
    for bar in axis.patches:
        assert isinstance(bar, Rectangle)
        widths.append(bar.get_width())
    np.testing.assert_allclose(widths, expected[score], equal_nan=True)
    assert [tick.get_text() for tick in axis.get_yticklabels()] == expected["trait_label"].tolist()
    pd.testing.assert_frame_equal(frame, original)


@pytest.mark.parametrize("study", [True, False], ids=["study", "contrast-fallback"])
def test_score_only_input_does_not_invent_raw_p_test_counts(
    study: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(infographic, "TOP_N_TRAITS", 3)
    frame = _stats(study).drop(columns="study_p_value" if study else "p_value")
    _, axis = _panel(0, 0)
    infographic.plot_top_traits(axis, frame)
    shown, supplied = (2, 4) if study else (3, 5)
    label = axis.get_xlabel()
    assert f"{shown} finite scores shown / {supplied} finite scores supplied" in label
    assert "Raw p-values unavailable; counts use -log10(p) scores" in label
    assert "finite tests" not in label


@pytest.mark.parametrize("study", [True, False], ids=["study", "contrast-fallback"])
def test_story_distinguishes_nominal_p_bh_q_and_actual_contrast_unit(study: bool) -> None:
    frame = _stats(study)
    _, axis = _panel(1, 1)
    infographic.plot_story_box(
        axis,
        frame,
        pd.DataFrame({"condition": ["Condition A"]}),
        pd.DataFrame({"body_site_group": ["Feces"]}),
        pair_count=999,
        study_count=111,
        trait_count=len(frame),
    )
    text = " ".join(item.get_text().replace("-\n", "-").replace("\n", " ") for item in axis.texts)
    basis = "study means of paired contrasts" if study else "paired experiment contrasts"
    assert f"Nominal p-values in A describe {basis}" in text
    assert f"BH-adjusted q-values here describe {basis}" in text
    prefix = "study" if study else "contrast"
    expected_q = ("6.0e-04", "0.003") if study else ("0.200", "0.400")
    assert all(f"{prefix} q={value}" in text for value in expected_q)
    if not study:
        assert "study q=" not in text
        assert "study-averaged" not in text


@pytest.mark.parametrize("study", [True, False], ids=["study", "contrast-fallback"])
@pytest.mark.parametrize("score_only", [False, True], ids=["raw-p", "scores-only"])
def test_added_nominal_text_stays_inside_the_figure_and_story_box(
    study: bool, score_only: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(infographic, "TOP_N_TRAITS", 3)
    frame = _stats(study)
    if score_only:
        frame = frame.drop(columns="study_p_value" if study else "p_value")
    figure, axis = _panel(0, 0)
    infographic.plot_top_traits(axis, frame)
    figure.canvas.draw()
    label_box = axis.xaxis.label.get_window_extent()
    assert figure.bbox.contains(label_box.x0, label_box.y0)
    assert figure.bbox.contains(label_box.x1, label_box.y1)
    assert label_box.x0 >= axis.bbox.x0
    assert label_box.x1 <= axis.bbox.x1
    story_figure, story_axis = _panel(1, 1)
    infographic.plot_story_box(
        story_axis,
        frame,
        pd.DataFrame({"condition": ["Condition A"]}),
        pd.DataFrame({"body_site_group": ["Feces"]}),
        pair_count=999,
        study_count=111,
        trait_count=len(frame),
    )
    story_figure.canvas.draw()
    caveat = next(item for item in story_axis.texts if item.get_text().startswith("- Caveat:"))
    extent = caveat.get_window_extent()
    box = story_axis.patches[0].get_window_extent()
    assert box.contains(extent.x0, extent.y0)
    assert box.contains(extent.x1, extent.y1)
