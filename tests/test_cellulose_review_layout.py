"""Actual Agg layout regressions for the current cellulose figures (#1866, #1868)."""

from __future__ import annotations

import itertools
import math
from typing import Any

import matplotlib
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.text import Text

from mtp.cellulose_decr import figures
from mtp.cellulose_decr.analysis import ANALYSIS_DIR
from test_cellulose_nominal_labels import synthetic_results


def _results(source: str) -> figures.Results:
    if source == "committed":
        return figures.load(ANALYSIS_DIR)
    assert source == "scientific_and_floor"
    results = synthetic_results()
    # These are distinct display cases, not new statistical calculations. In
    # particular, the Monte Carlo bound is 1/(100,000 + 1), not a short 0.040.
    for null in results.summary["null_models"].values():
        null["iterations"] = 100_000
    values = {
        "u1_hypergeom_p": 1.1e-14,
        "u2_hypergeom_p": 1.8e-11,
        "u3_hypergeom_p": math.nan,
        "u1_mc_p": 0.0,
        "u2_mc_p": 1.8e-11,
        "u3_mc_p": math.nan,
        "u4_shared_trait_p": 0.0,
        "u5_emp_weighted_p": 0.00042,
        "u6_signature_perm_p": 0.13,
        "u6_size_matched_p": math.nan,
    }
    for column, value in values.items():
        results.pairwise[column] = value
    return results


def _draw(figure: Figure) -> Any:
    canvas: Any = FigureCanvasAgg(figure)
    canvas.draw()
    return canvas.get_renderer()


@pytest.mark.parametrize("source", ["committed", "scientific_and_floor"])
def test_detailed_null_labels_fit_real_magnitudes(source: str) -> None:
    """Pin both dimensions, including adjacent scientific and floor labels."""
    with matplotlib.rc_context(figures.STYLE):
        figure = figures.detailed_nulls(_results(source))
        try:
            renderer = _draw(figure)
            table = figure.axes[0]
            labels = [text for text in table.texts if "nominal p" in text.get_text()]
            assert len(labels) == 28
            if source == "scientific_and_floor":
                texts = [label.get_text() for label in labels]
                assert any("1.1\N{MULTIPLICATION SIGN}10⁻¹⁴" in text for text in texts)
                assert any("1.8\N{MULTIPLICATION SIGN}10⁻¹¹" in text for text in texts)
                assert any(
                    "≤" in text and "1.0\N{MULTIPLICATION SIGN}10⁻⁵" in text for text in texts
                )
                assert any("unavailable" in text for text in texts)

            rows: dict[float, list[Text]] = {}
            columns: dict[float, list[Text]] = {}
            for label in labels:
                x, y = label.get_position()
                rows.setdefault(y, []).append(label)
                columns.setdefault(x, []).append(label)
                assert label.get_visible()
                assert float(label.get_fontsize()) >= 7.4
                bounds = label.get_window_extent(renderer)
                # Eight equal-width table columns include the observed column.
                # Check the allocated cell, not only whether neighbours happen
                # to be shorter or were omitted by a future renderer.
                cell_left = table.transAxes.transform((x - 1 / 16, y))[0]
                cell_right = table.transAxes.transform((x + 1 / 16, y))[0]
                assert bounds.x0 > cell_left, label.get_text()
                assert bounds.x1 < cell_right, label.get_text()
                expected = [
                    text
                    for text in table.texts
                    if text.get_position()[0] == x
                    and text.get_position()[1] > y
                    and "nominal p" not in text.get_text()
                ]
                assert expected
                nearest = min(expected, key=lambda text: text.get_position()[1] - y)
                assert bounds.y1 < nearest.get_window_extent(renderer).y0, label.get_text()

            assert len(rows) == 4
            assert len(columns) == 7
            for row in rows.values():
                assert len(row) == 7
                ordered = sorted(row, key=lambda text: text.get_position()[0])
                for left, right in itertools.pairwise(ordered):
                    assert (
                        left.get_window_extent(renderer).x1 < right.get_window_extent(renderer).x0
                    ), (left.get_text(), right.get_text())
            for column in columns.values():
                assert len(column) == 4
                ordered = sorted(column, key=lambda text: text.get_position()[1])
                for lower, upper in itertools.pairwise(ordered):
                    assert (
                        lower.get_window_extent(renderer).y1 < upper.get_window_extent(renderer).y0
                    ), (lower.get_text(), upper.get_text())
        finally:
            figure.clear()


@pytest.mark.parametrize("source", ["committed", "scientific_and_floor"])
def test_detailed_compartment_footnote_clears_axis_label(source: str) -> None:
    """The five-line footer needs a visible gap, without smaller type."""
    with matplotlib.rc_context(figures.STYLE):
        figure = figures.detailed_compartments(_results(source))
        try:
            renderer = _draw(figure)
            footnotes = [
                text
                for text in figure.texts
                if text.get_text().startswith("Experiments within a study are not independent;")
            ]
            xlabels = [
                axis.xaxis.label
                for axis in figure.axes
                if axis.get_xlabel().startswith("odds ratio of the ↓ share")
            ]
            assert len(footnotes) == 1
            assert len(xlabels) == 1
            footnote, xlabel = footnotes[0], xlabels[0]
            assert "Nominal p-values:" in footnote.get_text()
            assert float(footnote.get_fontsize()) >= 7.5
            assert float(xlabel.get_fontsize()) >= 8
            assert footnote.get_visible()
            assert xlabel.get_visible()
            footer_bounds = footnote.get_window_extent(renderer)
            xlabel_bounds = xlabel.get_window_extent(renderer)
            # The annotation spans the x-label horizontally, so vertical
            # separation must supply the gap. Two points survive export DPI.
            assert footer_bounds.x0 < xlabel_bounds.x1
            assert xlabel_bounds.x0 < footer_bounds.x1
            assert footer_bounds.y1 + renderer.points_to_pixels(2) <= xlabel_bounds.y0
            for bounds in (footer_bounds, xlabel_bounds):
                assert bounds.x0 >= figure.bbox.x0
                assert bounds.x1 <= figure.bbox.x1
                assert bounds.y0 >= figure.bbox.y0
                assert bounds.y1 <= figure.bbox.y1
        finally:
            figure.clear()
