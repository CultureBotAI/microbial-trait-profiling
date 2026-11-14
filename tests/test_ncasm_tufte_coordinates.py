"""Check plate 21's rendered artists against independent run-table reads (#2381).

The oracle never calls the production loader. Mutants change the live artists at the
save boundary, so a loader-only check cannot satisfy these regressions.
"""

from __future__ import annotations

import csv
from collections.abc import Callable, Iterator
from itertools import pairwise
from pathlib import Path
from typing import NamedTuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.axes import Axes
from matplotlib.collections import PathCollection
from matplotlib.figure import Figure
from matplotlib.patches import Polygon
from numpy.typing import NDArray

from mtp.ncasm_tufte_plates import story_from_tables, story_plots

FIXTURE_RUN = Path(__file__).parent / "fixtures/ncasm_tufte_plates/run_tier12_off_511c891"
ARM1 = "comparison/arm1_pooled_stats/trait_profile_paired_contrast_trait_stats_pooled.tsv"
SCREEN = "arms/noel/trait_direction_screen/trait_direction_screen.tsv"
Coordinates = NDArray[np.float64]


class Expected(NamedTuple):
    points: Coordinates
    rings: Coordinates
    focal: Coordinates


@pytest.fixture(scope="module")
def expected() -> Expected:
    def rows(relative: str) -> list[dict[str, str]]:
        with (FIXTURE_RUN / relative).open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    arm1 = rows(ARM1)
    votes = {row["trait"]: row for row in rows(SCREEN)}
    points = np.array(
        [
            (
                float(row["mean_increased_minus_decreased"]),
                100 * float(votes[row["trait"]]["share_up"]),
            )
            for row in arm1
        ],
        dtype=float,
    )
    rings = np.array([float(votes[row["trait"]]["by_q_nominal"]) <= 0.05 for row in arm1])
    focal = np.array([votes[row["trait"]]["in_table1"] == "true" for row in arm1])
    assert points.shape == (311, 2)
    assert len(votes) == 317
    assert (int(rings.sum()), int(focal.sum())) == (4, 6)
    return Expected(points, points[rings], points[focal])


@pytest.fixture(autouse=True)
def _small_clean_figures(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    monkeypatch.setattr(story_plots, "SLIDE_DPI", 40)
    with mpl.style.context("default"):
        try:
            yield
        finally:
            plt.close("all")


def _main_collections(fig: Figure) -> list[PathCollection]:
    axes = [ax for ax in fig.axes if ax.get_ylabel() == "Arm 3 Share of Studies Voting Up (%)"]
    assert len(axes) == 1
    return [artist for artist in axes[0].collections if isinstance(artist, PathCollection)]


def _offsets(collection: PathCollection) -> Coordinates:
    offsets = collection.get_offsets()
    assert not np.any(np.ma.getmask(offsets))
    return np.asarray(offsets, dtype=float)


def _same(actual: Coordinates, wanted: Coordinates, message: str) -> None:
    np.testing.assert_allclose(actual, wanted, rtol=0, atol=1e-12, err_msg=message)


def _sorted_pairs(values: Coordinates) -> Coordinates:
    return values[np.lexsort((values[:, 1], values[:, 0]))]


def _assert_marginal(ax: Axes, values: Coordinates, *, horizontal: bool) -> None:
    """Read the saved polygon's upper step at each bin midpoint, then count TSV values.

    The observed bin edges are layout; expected counts come directly from the fixture,
    without calling the renderer's histogram implementation or using its input arrays.
    """
    polygons = [patch for patch in ax.patches if isinstance(patch, Polygon)]
    assert len(polygons) == 1
    vertices = np.asarray(polygons[0].get_xy(), dtype=float)
    if horizontal:
        vertices = vertices[:, ::-1]
    edges = np.unique(vertices[:, 0])
    assert len(edges) == (22 if horizontal else 25)
    assert edges[0] <= values.min() <= values.max() <= edges[-1]
    heights, counts = [], []
    for index, (left, right) in enumerate(pairwise(edges)):
        midpoint = (left + right) / 2
        heights.append(
            max(
                start[1]
                for start, stop in pairwise(vertices)
                if start[1] == stop[1]
                and min(start[0], stop[0]) < midpoint < max(start[0], stop[0])
            )
        )
        counts.append(
            sum(
                left <= value and (value < right or (index == len(edges) - 2 and value <= right))
                for value in values
            )
        )
    assert sum(counts) == 311
    np.testing.assert_array_equal(heights, counts, err_msg="marginal histogram counts")


def _assert_rendered_coordinates(fig: Figure, expected: Expected) -> None:
    collections = _main_collections(fig)
    # Gray cloud, hollow BY rings, two rugs, and six individually coloured focal dots.
    assert len(collections) == 10
    _same(_offsets(collections[0]), expected.points, "cloud offsets")
    assert len(collections[1].get_facecolor()) == 0, "BY markers must be hollow rings"
    _same(_offsets(collections[1]), expected.rings, "BY ring offsets")
    _same(_offsets(collections[2])[:, 0], expected.points[:, 0], "x rug offsets")
    _same(_offsets(collections[3])[:, 1], expected.points[:, 1], "y rug offsets")
    focal = np.concatenate([_offsets(artist) for artist in collections[4:]])
    _same(_sorted_pairs(focal), _sorted_pairs(expected.focal), "focal dot offsets")

    histogram_axes = [ax for ax in fig.axes if any(isinstance(p, Polygon) for p in ax.patches)]
    assert len(histogram_axes) == 2
    top, right = histogram_axes
    _assert_marginal(top, expected.points[:, 0], horizontal=False)
    _assert_marginal(right, expected.points[:, 1], horizontal=True)
    top_focal = [artist for artist in top.collections if isinstance(artist, PathCollection)]
    right_focal = [artist for artist in right.collections if isinstance(artist, PathCollection)]
    assert len(top_focal) == len(right_focal) == 1
    _same(np.sort(_offsets(top_focal[0])[:, 0]), np.sort(expected.focal[:, 0]), "top focal rug")
    _same(np.sort(_offsets(right_focal[0])[:, 1]), np.sort(expected.focal[:, 1]), "right focal rug")


def _render(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    expected: Expected,
    mutate: Callable[[Figure], None] | None = None,
) -> None:
    save = story_plots.save_and_verify_story_plot
    checked: list[Figure] = []

    def checked_save(
        fig: Figure, stem: str, out_dir: Path, *, random_marks: bool
    ) -> story_plots.StoryPlotInfo:
        checked.append(fig)
        fig.canvas.draw()
        if mutate is not None:
            mutate(fig)
        _assert_rendered_coordinates(fig, expected)
        return save(fig, stem, out_dir, random_marks=random_marks)

    monkeypatch.setattr(story_from_tables, "save_and_verify_story_plot", checked_save)
    drawn = story_from_tables.plot_01_trait_cloud(tmp_path / "out", FIXTURE_RUN)
    assert len(checked) == 1
    for path in (drawn.info["png_path"], drawn.info["svg_path"], drawn.info["pdf_path"]):
        assert Path(path).stat().st_size > 0


def test_plate_21_saves_every_point_ring_and_marginal_at_the_tabulated_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, expected: Expected
) -> None:
    _render(tmp_path, monkeypatch, expected)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("reverse-cloud-x", "cloud offsets"),
        ("reverse-ring-x", "BY ring offsets"),
        ("wrong-ring-membership", "BY ring offsets"),
    ],
)
def test_rendered_coordinate_oracle_rejects_mutated_artists(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, expected: Expected, mutation: str, message: str
) -> None:
    changed: list[str] = []

    def mutate(fig: Figure) -> None:
        collections = _main_collections(fig)
        target = collections[0 if mutation == "reverse-cloud-x" else 1]
        original = _offsets(target)
        wrong = original.copy()
        if mutation == "wrong-ring-membership":
            wrong[0] = next(
                point
                for point in _offsets(collections[0])
                if not np.any(np.all(original == point, axis=1))
            )
        else:
            wrong[:, 0] = wrong[::-1, 0]
        assert not np.array_equal(wrong, original), "mutation must alter plotted coordinates"
        target.set_offsets(wrong)
        changed.append(mutation)

    with pytest.raises(AssertionError, match=message):
        _render(tmp_path, monkeypatch, expected, mutate)
    assert changed == [mutation]
