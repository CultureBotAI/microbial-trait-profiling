"""Publication renderers preserve the verified estimands and complete populations (#1046)."""

from __future__ import annotations

from collections.abc import Callable
from itertools import pairwise
from types import SimpleNamespace
from typing import TYPE_CHECKING, cast

import matplotlib
import numpy as np
import pandas as pd
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg, RendererAgg
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.figure import Figure

from mtp.publication_figures.plots import (
    forest_pages,
    quadrant_pages,
    table_pages,
    validate_authored_wording,
)

if TYPE_CHECKING:
    from mtp.publication_figures.inputs import FigureInputs


def evidence(count: int = 14) -> pd.DataFrame:
    effects = [9, 8, 7, 6, 5, 5, 1, -1, -2, -3, -4, -5, -5, 0]
    rows = []
    for index in range(count):
        effect = float(effects[index % len(effects)])
        mean = -effect / 3 + index / 100
        rows.append(
            {
                "trait": f"METPO:{index:07d}",
                "trait_label": "untrusted stale table label",
                "experiment_mean_delta_z": effect,
                "study_mean_delta_z": mean,
                "study_mean_delta_z_ci95_low": mean - 0.3,
                "study_mean_delta_z_ci95_high": mean + 0.7,
                "study_t_status": "tested",
                "study_t_reason": "",
                "study_t_q_by": 0.04,
                "composition_studies": 17,
                "composition_pairs": 29,
                "noel_status": "tested",
                "noel_studies_up": 2,
                "noel_studies_down": 7,
                "noel_studies_tied": 3,
                "noel_q_by_nominal": 0.2,
                "ora_direction_status": "tested",
                "ora_direction_studies_up": 4,
                "ora_direction_studies_down": 1,
                "ora_direction_studies_tied": 2,
                "ora_direction_studies_eligible": 7,
                "ora_direction_vote_preference": "increase",
            }
        )
    result = pd.DataFrame(rows)
    for field in ("experiment_mean_delta_z", "study_mean_delta_z"):
        result[f"{field}_rank_ascending"] = result[field].rank(method="min")
        result[f"{field}_rank_descending"] = result[field].rank(method="min", ascending=False)
    return result


def bundle(frame: pd.DataFrame | None = None, pairs: pd.DataFrame | None = None) -> FigureInputs:
    frame = evidence() if frame is None else frame
    labels = {trait: f"Verified annotated label {trait}" for trait in frame.trait}
    labels.update(
        {"METPO:1000608": "verified union label", "METPO:1000607": "verified obligate label"}
    )
    tables = {} if pairs is None else {"quadrant_pairs.tsv": pairs}
    if pairs is not None:
        rows = []
        for population, subset in (("all", pairs), ("stool", pairs.loc[pairs.is_stool])):
            counts = subset.quadrant.value_counts()
            rows.append(
                {
                    "population": population,
                    "fa_trait": "METPO:1000608",
                    "oa_trait": "METPO:1000607",
                    "fa_label": labels["METPO:1000608"],
                    "oa_label": labels["METPO:1000607"],
                    "n_comparisons": len(subset),
                    "n_studies": subset.study_id.nunique(),
                    **{
                        name: int(counts.get(name, 0))
                        for name in (
                            "fa_up_oa_down",
                            "fa_down_oa_up",
                            "both_up",
                            "both_down",
                            "axis",
                        )
                    },
                    "sign_test_status": "tested",
                    "sign_test_p": 0.125,
                    "pearson_status": "tested",
                    "pearson_r": -0.75,
                    "p_adjustment": "none_nominal",
                }
            )
        tables["quadrant_summary.tsv"] = pd.DataFrame(rows)
    return cast(
        "FigureInputs",
        SimpleNamespace(
            parent=SimpleNamespace(evidence=frame, record={"run_id": "parent-1"}),
            record={"strain_rollup": "phenotype", "contrast_population": "test-population"},
            tables=tables,
            labels=labels,
        ),
    )


def text(figure: Figure) -> str:
    return "\n".join(item.get_text() for item in figure.findobj(matplotlib.text.Text))


def points(figure: Figure, gid: str) -> PathCollection:
    candidates = [item for item in figure.findobj(PathCollection) if item.get_gid() == gid]
    assert len(candidates) == 1
    return candidates[0]


def test_forest_retains_boundary_ties_and_separates_estimates() -> None:
    source = bundle()
    pages = list(forest_pages(source))
    assert len(pages) == 1
    page = pages[0]
    expected = [f"METPO:{number:07d}" for number in [0, 1, 2, 3, 4, 5, 11, 12, 10, 9, 8]]
    assert page.row_ids == expected
    assert page.details["ranking"] == "comparison_mean"
    assert page.details["positive_count"] == 6
    assert page.details["negative_count"] == 5
    selected = source.parent.evidence.set_index("trait").loc[expected]
    np.testing.assert_array_equal(
        np.asarray(points(page.figure, "forest-study").get_offsets())[:, 0],
        selected.study_mean_delta_z,
    )
    np.testing.assert_array_equal(
        np.asarray(points(page.figure, "forest-comparison").get_offsets())[:, 0],
        selected.experiment_mean_delta_z,
    )
    intervals = [item for item in page.figure.axes[0].collections if item.get_gid() == "forest-ci"]
    assert len(intervals) == 1
    assert isinstance(intervals[0], LineCollection)
    np.testing.assert_allclose(
        np.asarray(intervals[0].get_segments())[:, :, 0],
        selected[["study_mean_delta_z_ci95_low", "study_mean_delta_z_ci95_high"]],
    )
    assert "untrusted stale" not in text(page.figure)
    assert "Verified annotated label" in text(page.figure)
    assert "Comparison mean" in text(page.figure)
    assert "Study mean" in text(page.figure)


def test_study_ranking_is_explicit_and_can_select_other_traits() -> None:
    source = bundle()
    first = next(forest_pages(source, ranking="comparison_mean", top_n=1))
    second = next(forest_pages(source, ranking="study_mean", top_n=1))
    assert first.row_ids == ["METPO:0000000", "METPO:0000011", "METPO:0000012"]
    assert second.row_ids == ["METPO:0000012", "METPO:0000000"]
    assert second.details["ranking"] == "study_mean"


@pytest.mark.parametrize("count", [10, 22])
def test_forest_footer_legend_and_axis_label_have_separate_strips(count: int) -> None:
    frame = evidence(count)
    effects = [5, 4, 3, 2, 1, -1, -2, -3, -4, -5] if count == 10 else [1] * 11 + [-1] * 11
    frame["experiment_mean_delta_z"] = effects
    frame["experiment_mean_delta_z_rank_ascending"] = frame.experiment_mean_delta_z.rank(
        method="min"
    )
    frame["experiment_mean_delta_z_rank_descending"] = frame.experiment_mean_delta_z.rank(
        method="min", ascending=False
    )
    source = bundle(frame)
    source.labels.update(
        {trait: f"Verified annotation {index}" for index, trait in enumerate(frame.trait)}
    )
    page = next(forest_pages(source))
    assert len(page.row_ids) == count
    assert page.mark_count == 2 * count
    canvas = cast(FigureCanvasAgg, page.figure.canvas)
    cast(Callable[[], None], canvas.draw)()
    renderer = cast(Callable[[], RendererAgg], canvas.get_renderer)()
    footer = next(item for item in page.figure.texts if item.get_text().startswith("Full-family"))
    bounds = [
        footer.get_window_extent(renderer),
        page.figure.legends[0].get_window_extent(renderer),
        page.figure.axes[0].xaxis.label.get_window_extent(renderer),
    ]
    for lower, upper in pairwise(bounds):
        assert lower.y1 + 0.04 * page.figure.dpi < upper.y0
        assert not lower.overlaps(upper)
    assert bounds[0].y0 >= 0
    for item in bounds:
        assert item.x0 >= 0
        assert item.x1 <= page.figure.bbox.width


@pytest.mark.parametrize("top_n", [0, -1, True, 1.5])
def test_invalid_top_n_refused(top_n: object) -> None:
    with pytest.raises(ValueError, match="positive integer"):
        list(forest_pages(bundle(), top_n=cast(int, top_n)))


def test_forest_untestable_interval_stays_missing() -> None:
    frame = evidence()
    frame.loc[
        0, ["study_mean_delta_z_ci95_low", "study_mean_delta_z_ci95_high", "study_t_q_by"]
    ] = np.nan
    frame.loc[0, "study_t_status"] = "untestable"
    frame.loc[0, "study_t_reason"] = "insufficient_variance"
    page = next(forest_pages(bundle(frame), top_n=1))
    assert "METPO:0000000" in page.row_ids
    intervals = next(
        item for item in page.figure.axes[0].collections if item.get_gid() == "forest-ci"
    )
    assert isinstance(intervals, LineCollection)
    assert len(intervals.get_segments()) == 2
    assert "untestable" in text(page.figure).lower()


def quadrant_rows() -> pd.DataFrame:
    # More than common plotting caps; repeated coordinates and exact axis points all survive.
    count = 1207
    frame = pd.DataFrame(
        {
            "study_id": [f"s{index // 3}" for index in range(count)],
            "experiment_id": [f"e{index}" for index in range(count)],
            "body_site": ["Feces" if index % 3 else "Feces,Blood" for index in range(count)],
            "fa_trait": "METPO:1000608",
            "oa_trait": "METPO:1000607",
            "fa_delta_z": 1.0,
            "oa_delta_z": -2.0,
            "quadrant": "fa_up_oa_down",
        }
    )
    frame.loc[1, "fa_delta_z"] = 0.0
    frame.loc[1, "quadrant"] = "axis"
    frame.loc[2, "fa_delta_z"] = -0.4
    frame.loc[2, "oa_delta_z"] = 3.0
    frame.loc[2, "quadrant"] = "fa_down_oa_up"
    frame["is_stool"] = frame.body_site.eq("Feces")
    return frame


def test_quadrant_keeps_every_pair_axis_and_duplicate_without_jitter() -> None:
    rows = quadrant_rows()
    page = next(quadrant_pages(bundle(pairs=rows)))
    all_points = np.asarray(points(page.figure, "quadrant-all").get_offsets())
    stool_points = np.asarray(points(page.figure, "quadrant-stool").get_offsets())
    np.testing.assert_array_equal(all_points, rows[["fa_delta_z", "oa_delta_z"]])
    np.testing.assert_array_equal(
        stool_points, rows.loc[rows.is_stool, ["fa_delta_z", "oa_delta_z"]]
    )
    assert len(all_points) == 1207
    assert len(stool_points) == 804
    assert page.mark_count == 2011
    assert len(page.row_ids) == len(set(page.row_ids)) == 2011
    assert "stool:s0/e0" not in page.row_ids
    assert "stool:s0/e1" in page.row_ids
    assert page.figure.axes[0].get_xlim() == page.figure.axes[1].get_xlim()
    assert page.figure.axes[0].get_ylim() == page.figure.axes[1].get_ylim()
    assert "verified union label" in text(page.figure)
    assert "facultatively aerobic" not in text(page.figure)
    assert "Stool" in text(page.figure)


def test_table_paginates_complete_catalog_with_no_significance_selection() -> None:
    frame = evidence(87).iloc[::-1].reset_index(drop=True)
    frame.loc[0, "noel_q_by_nominal"] = 1.0
    pages = list(table_pages(bundle(frame)))
    ids = [trait for page in pages for trait in page.row_ids]
    assert len(pages) > 1
    assert ids == sorted(frame.trait)
    assert len(ids) == len(set(ids)) == 87
    assert sum(page.mark_count for page in pages) == 87
    for page in pages:
        content = text(page.figure)
        assert "Noel" in content
        assert "Study mean" in content
        assert "ORA" in content
        assert "uncalibrated" in content.lower()
        assert "untrusted stale" not in content
        assert page.figure.get_size_inches()[0] == pytest.approx(7.2)
        assert all(
            float(item.get_fontsize()) >= 8 for item in page.figure.findobj(matplotlib.text.Text)
        )


def test_table_missing_is_not_zero_and_verified_words_are_not_censored() -> None:
    frame = evidence(1)
    frame.loc[
        0, ["study_t_q_by", "study_mean_delta_z_ci95_low", "study_mean_delta_z_ci95_high"]
    ] = np.nan
    frame.loc[0, "study_t_status"] = "untestable"
    frame.loc[0, "noel_q_by_nominal"] = 0.0
    source = bundle(frame)
    source.labels["METPO:0000000"] = "verified redox annotation label"
    page = next(table_pages(source))
    content = text(page.figure)
    assert "verified redox" in content
    assert "BY q = 0" in content
    assert "95% CI: unavailable" in content
    assert "untestable" in content
    # Authored interpretation never adds prohibited mechanism claims.
    assert "coupled" not in content.lower()
    assert "physiological" not in content.lower()


def test_table_text_stays_inside_page_and_outside_footer() -> None:
    source = bundle(evidence(21))
    source.labels["METPO:0000003"] = "A verified label whose full wording must remain visible " * 3
    source.parent.evidence.loc[0, "noel_studies_up"] = 12345
    pages = list(table_pages(source))
    assert "12345 / 7 / 3" in text(pages[0].figure)
    for page in pages:
        canvas = cast(FigureCanvasAgg, page.figure.canvas)
        # Matplotlib's concrete Agg methods are unannotated; bind their known interfaces.
        cast(Callable[[], None], canvas.draw)()
        renderer = cast(Callable[[], RendererAgg], canvas.get_renderer)()
        for item in page.figure.findobj(matplotlib.text.Text):
            if not item.get_text():
                continue
            bounds = item.get_window_extent(renderer)
            assert bounds.x0 >= 0
            assert bounds.x1 <= page.figure.bbox.width
            assert bounds.y0 >= 0
            assert bounds.y1 <= page.figure.bbox.height
            if (item.get_gid() or "").startswith("table-cell:"):
                assert bounds.y0 > 0.7 * page.figure.dpi


def test_quadrant_empty_population_is_explicit_and_still_draws_two_panels() -> None:
    page = next(quadrant_pages(bundle(pairs=quadrant_rows().head(0))))
    assert page.mark_count == 0
    assert page.row_ids == []
    assert len(page.figure.axes) == 2
    assert "0 comparisons; 0 studies" in text(page.figure)


@pytest.mark.parametrize("authored", ["Coupled oxygen response", "Redox shift", "Causes disease"])
def test_authored_wording_guard_rejects_mechanism_claims(authored: str) -> None:
    with pytest.raises(ValueError, match="authored"):
        validate_authored_wording(authored)
