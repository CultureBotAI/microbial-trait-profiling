"""Data-story plots drawn from run tables instead of a random generator (#2264, #775).

Each builder here takes the output folder and a run root, reads the tables through
``mtp.ncasm_tufte_plates.run_tables`` and draws every point, count and label from what
it read: nothing printed on the plate is typed here except layout. The speaker note is
written from the same values (#776).

Plot 01 keeps the design of the original joint cloud. Its random 233-point background,
typed focal positions and typed summary numbers are replaced by the run's per-trait
values, so its title, counts and callouts differ from the 2026-10-22 render.
"""

from __future__ import annotations

import itertools
import math
from pathlib import Path
from typing import NamedTuple

import matplotlib as mpl
import numpy as np
from matplotlib.ticker import FuncFormatter, MaxNLocator

from mtp.ncasm_tufte_plates.data_status import InputFile
from mtp.ncasm_tufte_plates.run_tables import RunIdentity, TraitCloud, TraitPoint, load_trait_cloud
from mtp.ncasm_tufte_plates.story_common import (
    ACCENT,
    FAINT,
    INK,
    MUTED,
    PAPER,
    RULE,
    SECONDARY,
    SLIDE_H,
    SLIDE_SIZES,
    TERTIARY,
    add_panel,
    despine,
)
from mtp.ncasm_tufte_plates.story_plots import (
    StoryPlotInfo,
    draw_story_callout,
    init_story_plot,
    save_and_verify_story_plot,
)

__all__ = ["TableStoryPlot", "plot_01_trait_cloud", "trait_cloud_notes"]

MINUS = "\u2212"
# Plot 01 layout: the scatter axes in figure fractions, and the label bands in points.
AXES_RECT = (0.08, 0.14, 0.55, 0.45)
LABEL_STEP_PT = 15.5
LABEL_PAD_PT = 9.0
RUG_PT = 9.0
COUNT_WORDS = ("No", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten")


class TableStoryPlot(NamedTuple):
    """A plot drawn from run tables, with what it read and the note written from it."""

    info: StoryPlotInfo
    inputs: tuple[InputFile, ...]
    run: RunIdentity
    speaker_notes: str


def _signed(value: float, digits: int = 2) -> str:
    return f"{value:+.{digits}f}".replace("-", MINUS)


def _count_word(count: int) -> str:
    return COUNT_WORDS[count] if count < len(COUNT_WORDS) else str(count)


def _q_value(value: float) -> str:
    return f"{value:.3f}" if value >= 0.001 else f"{value:.1e}"


def _focal_label(point: TraitPoint) -> str:
    return (
        f"{point.name}  {_signed(point.mean_shift)} \u00b7 {point.percent_up:.1f}% up "
        f"({point.studies_up}:{point.studies_down}) \u00b7 BY q = {_q_value(point.by_q)}"
    )


def _leader_start(
    dot: tuple[float, float],
    label_left: float,
    label_right: float,
    label_y: float,
    gap: float,
    half_height: float,
) -> tuple[float, float]:
    """Where a leader line leaves its label: the side of the label nearest the dot."""
    if dot[0] > label_right:
        return label_right + gap, label_y
    if dot[0] < label_left:
        return label_left - gap, label_y
    return dot[0], label_y + (half_height if dot[1] > label_y else -half_height)


def _crosses(
    a: tuple[float, float], b: tuple[float, float], c: tuple[float, float], d: tuple[float, float]
) -> bool:
    """Whether segment a-b properly crosses segment c-d."""

    def side(p: tuple[float, float], q: tuple[float, float], r: tuple[float, float]) -> float:
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])

    return side(a, b, c) * side(a, b, d) < 0 and side(c, d, a) * side(c, d, b) < 0


def _hits_box(
    a: tuple[float, float], b: tuple[float, float], box: tuple[float, float, float, float]
) -> bool:
    """Whether segment a-b touches the rectangle (left, bottom, right, top)."""
    left, bottom, right, top = box
    if any(left <= p[0] <= right and bottom <= p[1] <= top for p in (a, b)):
        return True
    corners = [(left, bottom), (right, bottom), (right, top), (left, top)]
    return any(_crosses(a, b, corners[i], corners[(i + 1) % 4]) for i in range(4))


def _rows_without_crossings(
    dots: list[tuple[float, float]],
    label_left: float,
    label_ends: list[float],
    rows: list[float],
    gap: float,
    half_height: float,
) -> list[int]:
    """The row for each label that keeps leader lines off other labels and each other.

    Every assignment is tried for up to seven labels. The fewest lines over another
    label win, then the fewest crossing lines, then the shortest total length, then the
    first assignment found. More labels keep their given order."""
    count = len(dots)
    if count > 7:
        return list(range(count))
    best: tuple[int, int, float] | None = None
    best_order = list(range(count))
    for order in itertools.permutations(range(count)):
        starts = [(label_ends[i] + gap, rows[row]) for i, row in enumerate(order)]
        boxes = [
            (label_left, rows[row] - half_height, label_ends[i], rows[row] + half_height)
            for i, row in enumerate(order)
        ]
        over_labels = sum(
            _hits_box(starts[i], dots[i], boxes[j])
            for i in range(count)
            for j in range(count)
            if i != j
        )
        crossings = sum(
            _crosses(starts[i], dots[i], starts[j], dots[j])
            for i in range(count)
            for j in range(i + 1, count)
        )
        length = sum(math.dist(starts[i], dots[i]) for i in range(count))
        if best is None or (over_labels, crossings, length) < best:
            best, best_order = (over_labels, crossings, length), list(order)
    return best_order


def trait_cloud_notes(cloud: TraitCloud) -> str:
    """The speaker note for plot 01, written from the same values the plate shows."""
    table1 = cloud.table1
    passing = [point.name for point in table1 if point.passes_by]
    return (
        f"Data Story Plot 01 (joint trait cloud): each of the {len(cloud.points)} traits "
        f"tested by both arms is placed by its mean paired shift across "
        f"{cloud.comparisons:,} comparisons ({cloud.arm1_studies} studies) and by the share "
        f"of studies in which it rises in the study sign screen ({cloud.screen_traits} "
        f"traits, up to {cloud.screen_studies_max} studies). The two measures agree "
        f"(Pearson r = {cloud.correlation:.2f}). {len(cloud.screen_by_hits)} of "
        f"{cloud.screen_traits} traits pass Benjamini-Yekutieli q <= 0.05 in the sign "
        f"screen: {', '.join(cloud.screen_by_hits)}. Of the {len(table1)} Table 1 traits, "
        f"{len(passing)} pass ({', '.join(passing) if passing else 'none'}). These are "
        "alignments between reported taxa and trait annotations, not measured physiology."
    )


def plot_01_trait_cloud(out_dir: Path, run_root: Path) -> TableStoryPlot:
    """Plot 01: every trait's pooled shift against its share of studies voting up."""
    cloud = load_trait_cloud(run_root)
    points = cloud.points
    table1 = cloud.table1
    passing_table1 = [point for point in table1 if point.passes_by]
    x = np.array([point.mean_shift for point in points])
    y = np.array([point.percent_up for point in points])
    by_hit = np.array([point.passes_by for point in points])

    fig = init_story_plot(
        1,
        "Tufte Joint Scatter + Marginal Histograms (cf. tufte_023, tufte_101)",
        f"Where the {_count_word(len(table1))} Table 1 Traits Sit Among {len(points)} Traits "
        "Tested by Both Arms",
        f"One dot per trait: mean paired shift across {cloud.comparisons:,} comparisons "
        "(\u0394z) against the share of studies in which the trait rises.",
        f"Source: Arm 1 pooled shift ({cloud.comparisons:,} comparisons, "
        f"{cloud.arm1_studies} studies) and Arm 3 study sign screen ({cloud.screen_traits} "
        f"traits); run {cloud.identity.git_commit[:7]}.",
    )

    ax = fig.add_axes(AXES_RECT, facecolor=PAPER)
    ax_top = fig.add_axes((0.08, 0.595, 0.55, 0.085), facecolor=PAPER)
    ax_right = fig.add_axes((0.635, 0.14, 0.075, 0.45), facecolor=PAPER)
    fig.text(
        0.08,
        0.735,
        f"A. {len(points)} Traits: Mean Paired Shift (\u0394z) vs. Share of Studies Voting Up",
        fontsize=SLIDE_SIZES["panel_title"],
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.08,
        0.703,
        f"Gray dots = {len(points) - len(table1)} other traits; dark rings = the "
        f"{int(by_hit.sum())} traits passing BY q \u2264 0.05 in the sign screen",
        fontsize=12.0,
        color=MUTED,
        ha="left",
        va="top",
    )

    # Limits from the data. The x range is symmetric. In y, the labels of the rising
    # Table 1 traits are stacked in a band above the highest point and the labels of the
    # falling ones in a band below the lowest, so no label can cover a dot whatever the
    # data are. Band heights are fixed in points and converted with the axes height.
    rising = sorted((p for p in table1 if p.mean_shift >= 0), key=lambda p: -p.mean_shift)
    falling = sorted((p for p in table1 if p.mean_shift < 0), key=lambda p: p.mean_shift)
    axes_height_pt = AXES_RECT[3] * SLIDE_H * 72.0
    band_top = (len(rising) * LABEL_STEP_PT + LABEL_PAD_PT) / axes_height_pt
    band_bottom = (len(falling) * LABEL_STEP_PT + LABEL_PAD_PT + RUG_PT) / axes_height_pt
    data_low, data_high = float(np.floor(y.min())) - 1.0, float(np.ceil(y.max())) + 1.0
    y_span = (data_high - data_low) / (1.0 - band_top - band_bottom)
    y_low = data_low - band_bottom * y_span
    y_high = data_high + band_top * y_span
    step = LABEL_STEP_PT / axes_height_pt * y_span
    x_limit = 1.30 * float(np.abs(x).max())
    x_span = 2.0 * x_limit

    ax.axvline(0.0, color=RULE, lw=1.1, ls="--", zorder=1)
    ax.axhline(50.0, color=RULE, lw=1.1, ls="--", zorder=1)
    ax.scatter(x, y, color=FAINT, edgecolor="#8E887E", linewidth=0.4, s=46, alpha=0.72, zorder=2)
    ax.scatter(
        x[by_hit], y[by_hit], facecolor="none", edgecolor=INK, linewidth=1.4, s=150, zorder=4
    )
    # Rug ticks along the bottom and left frames
    ax.scatter(x, np.full_like(x, y_low + 0.02 * y_span), marker="|", color=MUTED, s=28, alpha=0.45)
    ax.scatter(
        np.full_like(y, -x_limit + 0.015 * x_span), y, marker="_", color=MUTED, s=28, alpha=0.45
    )

    ax.set_xlim(-x_limit, x_limit)
    ax.set_ylim(y_low, y_high)
    label_x = -x_limit + 0.035 * x_span
    gap = 0.008 * x_span
    placements: list[tuple[TraitPoint, float, str]] = []
    for group, rows, color in (
        (rising, [y_high - (row + 0.6) * step for row in range(len(rising))], ACCENT),
        (falling, [data_low - (row + 0.7) * step for row in range(len(falling))], SECONDARY),
    ):
        # A label's width is known only once it has been drawn, and the row order that
        # keeps leader lines from crossing depends on those widths.
        drawn = [
            ax.text(
                label_x,
                rows[0] if rows else 0.0,
                _focal_label(point),
                fontsize=11.2,
                fontweight="bold",
                color=color,
                ha="left",
                va="center",
                zorder=6,
            )
            for point in group
        ]
        fig.canvas.draw()
        to_data = ax.transData.inverted()
        ends = [float(to_data.transform((label.get_window_extent().x1, 0.0))[0]) for label in drawn]
        order = _rows_without_crossings(
            [(point.mean_shift, point.percent_up) for point in group],
            label_x,
            ends,
            rows,
            gap,
            0.45 * step,
        )
        for index, row in enumerate(order):
            point, label, text_y = group[index], drawn[index], rows[row]
            label.set_y(text_y)
            placements.append((point, text_y, color))
            ax.scatter(
                [point.mean_shift],
                [point.percent_up],
                color=color,
                s=125,
                edgecolor=PAPER,
                lw=1.3,
            ).set_zorder(5)
            start_x, start_y = _leader_start(
                (point.mean_shift, point.percent_up), label_x, ends[index], text_y, gap, 0.45 * step
            )
            ax.plot(
                [start_x, point.mean_shift],
                [start_y, point.percent_up],
                color=color,
                lw=1.1,
                zorder=4,
            )

    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, symmetric=True))
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda value, _position: "0.00" if value == 0 else _signed(value))
    )
    ax.set_yticks([tick for tick in range(0, 101, 5) if data_low <= tick <= data_high])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _position: f"{value:.0f}%"))
    ax.set_xlabel(
        f"Arm 1 Mean Paired Shift \u0394z, Increased {MINUS} Decreased "
        f"({cloud.comparisons:,} Comparisons)",
        fontsize=SLIDE_SIZES["axis_label"],
    )
    ax.set_ylabel("Arm 3 Share of Studies Voting Up (%)", fontsize=SLIDE_SIZES["axis_label"])
    despine(ax, top=True, right=True, left=False, bottom=False)

    focal_x = [point.mean_shift for point, _text_y, _color in placements]
    focal_y = [point.percent_up for point, _text_y, _color in placements]
    focal_colors = [color for _point, _text_y, color in placements]

    counts_x, _edges, _patches = ax_top.hist(
        x,
        bins=np.linspace(-x_limit, x_limit, 25).tolist(),
        color=FAINT,
        edgecolor=INK,
        lw=0.8,
        histtype="stepfilled",
        alpha=0.65,
    )
    ax_top.axvline(0.0, color=RULE, lw=1.0, ls="--")
    ax_top.scatter(
        focal_x, np.full(len(focal_x), 0.92 * float(np.max(counts_x))), color=focal_colors, s=45
    ).set_zorder(5)
    ax_top.set_xlim(-x_limit, x_limit)
    ax_top.axis("off")

    counts_y, _edges, _patches = ax_right.hist(
        y,
        bins=np.linspace(y_low, y_high, 22).tolist(),
        orientation="horizontal",
        color=FAINT,
        edgecolor=INK,
        lw=0.8,
        histtype="stepfilled",
        alpha=0.65,
    )
    ax_right.axhline(50.0, color=RULE, lw=1.0, ls="--")
    ax_right.scatter(
        np.full(len(focal_y), 0.92 * float(np.max(counts_y))), focal_y, color=focal_colors, s=45
    ).set_zorder(5)
    ax_right.set_ylim(y_low, y_high)
    ax_right.axis("off")

    ax_info = add_panel(fig, (0.74, 0.15, 0.22, 0.58), "B. Summary")
    ax_info.set_xlim(0, 1)
    ax_info.set_ylim(0, 1)
    ax_info.axis("off")
    draw_story_callout(
        ax_info,
        0.0,
        0.70,
        1.0,
        0.25,
        f"{len(points)} / {cloud.screen_traits}",
        "Traits tested by both arms",
        f"{cloud.screen_traits - len(points)} screen traits lack an Arm 1 row",
        SECONDARY,
    )
    draw_story_callout(
        ax_info,
        0.0,
        0.37,
        1.0,
        0.25,
        f"r = {cloud.correlation:.2f}",
        "Agreement between arms",
        "Mean shift vs. share voting up",
        ACCENT,
    )
    draw_story_callout(
        ax_info,
        0.0,
        0.04,
        1.0,
        0.25,
        f"{len(cloud.screen_by_hits)} / {cloud.screen_traits}",
        "Pass BY q \u2264 0.05 (sign test)",
        f"{len(passing_table1)} of the {len(table1)} Table 1 traits pass",
        TERTIARY,
    )

    # Text stays text in the SVG, so the plate's numbers can be read back and checked.
    with mpl.rc_context({"svg.fonttype": "none"}):
        info = save_and_verify_story_plot(
            fig, "plot_01_joint_effect_consistency_cloud", out_dir, random_marks=False
        )
    return TableStoryPlot(
        info=info,
        inputs=cloud.inputs,
        run=cloud.identity,
        speaker_notes=trait_cloud_notes(cloud),
    )
