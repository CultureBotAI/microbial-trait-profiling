"""Layout, table, callout and geometry-check helpers shared by NCASM plates I-XX.

Ported from ``tufte_plates_common.py`` in the slide-wizard deck folder
``slides/NCASM-disease-microbiome-traits/`` (sha256 ``d03019d9ee1e``, 2026-10-19; #2264).
The drawing calls are unchanged. The differences are the imports (the style helpers now
live in ``mtp.ncasm_tufte_plates.tufte_style``), type annotations, and the figure width.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import TypedDict

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from mtp.ncasm_tufte_plates import tufte_style
from mtp.ncasm_tufte_plates.tufte_style import (
    ACCENT_CRIMSON,
    ACCENT_OCHRE,
    ACCENT_SLATE,
    BG_CARD,
    BG_PAPER,
    INK_BODY,
    INK_MAIN,
    RULE_LIGHT,
    TINT_SLATE,
    add_panel_header,
    make_tufte_plate,
    relax_vertical_positions,
    save_plate,
)

# The builders ask for a 13.333 x 7.5 in canvas. The 2026-10-22 renders were drawn with
# matplotlib's macOS backend, which resizes a figure to whole screen pixels: 1,333 px at
# 100 dpi, so 13.33 in. That size is fixed here, so the headless Agg backend the Justfile
# exports draws the same pixels as those renders (#2264, #719).
PLATE_WIDTH_IN = 13.33
PLATE_HEIGHT_IN = 7.5
PLATE_DPI = 300


class SlopeSeries(TypedDict, total=False):
    """One slopegraph trajectory, as the plate code writes it."""

    name: str
    values: list[float]
    color: str
    lw: float
    ls: str
    marker_size: float
    zorder: int
    left_label: str
    right_label: str
    label_weight: str


class PlateInfo(TypedDict):
    """What one plate build wrote, with the layout check's findings."""

    plate_num: int
    stem: str
    png_path: str
    svg_path: str
    pdf_path: str
    out_of_bounds: list[str]
    overlaps: list[str]


# Color aliases used across Part 1 & Part 2
TUFTE_PAPER = BG_PAPER
TUFTE_SAND = BG_CARD
TUFTE_INK = INK_MAIN
TUFTE_MUTED = INK_BODY
TUFTE_LIGHT = RULE_LIGHT
TUFTE_FAINT = "#E5E0D5"
TUFTE_HAIRLINE = "#EAE5DC"
TUFTE_GHOST = TINT_SLATE
TUFTE_CRIMSON = ACCENT_CRIMSON
TUFTE_SLATE = ACCENT_SLATE
TUFTE_AMBER = ACCENT_OCHRE
TUFTE_FOREST = "#2D5A3C"

ROMAN = [
    "I",
    "II",
    "III",
    "IV",
    "V",
    "VI",
    "VII",
    "VIII",
    "IX",
    "X",
    "XI",
    "XII",
    "XIII",
    "XIV",
    "XV",
    "XVI",
    "XVII",
    "XVIII",
    "XIX",
    "XX",
    "XXI",
    "XXII",
    "XXIII",
    "XXIV",
    "XXV",
    "XXVI",
    "XXVII",
    "XXVIII",
    "XXIX",
    "XXX",
]


def sanitize_text(s: str) -> str:
    """Replace Unicode symbols missing from macOS Palatino with clean typographic equivalents."""
    if not s:
        return s
    replacements = {
        "\u2191": "up",
        "\u2193": "down",
        "\u2192": "->",
        "\u2194": "vs.",
        "\u2229": "&",
        "\u2208": "in",
        "\u1d62": "_i",
        "\u2713": "(ok)",
        "\u2514": "-",
    }
    for k, v in replacements.items():
        if k in s:
            s = s.replace(k, v)
    return s


def new_plate(
    plate_num: int,
    module_tag: str,
    role_tag: str,
    title: str,
    subtitle: str,
    footer: str,
) -> Figure:
    """Create a 16:9 widescreen Tufte plate with ``tufte_style.make_tufte_plate``."""
    kicker = sanitize_text(
        f"PLATE {ROMAN[plate_num - 1]}   ·   {module_tag.upper()}   ·   {role_tag.upper()}"
    )
    plate_id = f"NCASM · Plate {ROMAN[plate_num - 1]} ({plate_num}/20)"
    clean_footer = sanitize_text(footer.removeprefix("Source: ").strip())
    fig = make_tufte_plate(
        kicker=kicker,
        title=sanitize_text(title),
        subtitle=sanitize_text(subtitle),
        source_note=clean_footer,
        plate_id=plate_id,
        figsize=(PLATE_WIDTH_IN, PLATE_HEIGHT_IN),
    )
    return fig


def two_panels(
    fig: Figure,
    left_rect: tuple[float, float, float, float] = (0.055, 0.145, 0.415, 0.610),
    right_rect: tuple[float, float, float, float] = (0.535, 0.145, 0.415, 0.610),
) -> tuple[Axes, Axes]:
    """Add two generously spaced side-by-side panels on the archival canvas."""
    ax_a = fig.add_axes(left_rect, facecolor=TUFTE_PAPER)
    ax_b = fig.add_axes(right_rect, facecolor=TUFTE_PAPER)
    return ax_a, ax_b


def title_panel(ax: Axes, title: str, subtitle: str = "") -> None:
    """Apply the standard panel title and subtitle with ``tufte_style.add_panel_header``."""
    add_panel_header(ax, title=sanitize_text(title), subcaption=sanitize_text(subtitle), gap_pt=5.0)


def style_tufte_axis(
    ax: Axes,
    *,
    keep_bottom: bool = True,
    keep_left: bool = True,
) -> None:
    """Configure clean minimalist Tufte spines."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(keep_bottom)
    ax.spines["left"].set_visible(keep_left)
    if not keep_left:
        ax.tick_params(axis="y", length=0)


def apply_range_frame(
    ax: Axes,
    x_bounds: tuple[float, float] | None = None,
    y_bounds: tuple[float, float] | None = None,
    outward_pt: float = 5.0,
    keep_labels: bool = False,
) -> None:
    """Wrap ``tufte_style.apply_range_frame``; y-category labels go in the panel's left column."""
    if y_bounds is None and not keep_labels:
        yticks = list(ax.get_yticks())
        raw_labels = [sanitize_text(t.get_text().strip()) for t in ax.get_yticklabels()]
        if yticks and any(raw_labels):
            xmin_orig, xmax_orig = ax.get_xlim()
            span_orig = xmax_orig - xmin_orig
            max_chars = max(
                (len(line) for lbl in raw_labels for line in lbl.split("\n")),
                default=18,
            )
            pos = ax.get_position()
            panel_w_pt = pos.width * 13.333 * 72.0
            lbl_w_pt = max_chars * 4.15 + 14.0
            f_left = min(0.43, max(0.31, lbl_w_pt / max(panel_w_pt, 200.0)))
            xmax_new = xmin_orig + span_orig * 1.22
            right_span = xmax_new - xmin_orig
            xmin_new = xmin_orig - right_span * (f_left / (1.0 - f_left))
            ax.set_xlim(xmin_new, xmax_new)
            ax.set_yticks([])
            for y_val, lbl_str in zip(yticks, raw_labels, strict=False):
                if not lbl_str:
                    continue
                ax.text(
                    xmin_new,
                    y_val,
                    lbl_str,
                    ha="left",
                    va="center",
                    fontsize=8.2,
                    color=TUFTE_INK,
                    linespacing=1.10,
                )

    tufte_style.apply_range_frame(
        ax,
        x_bounds=x_bounds,
        y_bounds=y_bounds,
        outward_pt=outward_pt,
        keep_labels=keep_labels,
    )


def setup_canvas_panel(ax: Axes) -> None:
    """Configure an axis as a clean 0..1 typographic/diagram canvas without ticks or spines."""
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_labeled_hline(
    ax: Axes,
    y: float,
    label: str,
    *,
    x_text: float = 0.5,
    dy: float = 1.0,
    color: str = TUFTE_MUTED,
    size: float = 8.8,
) -> None:
    """Draw a dashed horizontal reference line with a direct text label."""
    ax.axhline(y, color=color, lw=0.9, ls="--", zorder=1)
    ax.text(
        x_text, y + dy, sanitize_text(label), fontsize=size, color=color, ha="center", va="bottom"
    )


def draw_labeled_vline(
    ax: Axes,
    x: float,
    label: str,
    *,
    y_text: float = 0.5,
    dx: float = 0.05,
    color: str = TUFTE_MUTED,
    size: float = 8.8,
) -> None:
    """Draw a dashed vertical reference line with a direct text label."""
    ax.axvline(x, color=color, lw=0.9, ls="--", zorder=1)
    ax.text(
        x + dx, y_text, sanitize_text(label), fontsize=size, color=color, ha="left", va="center"
    )


def draw_slopegraph_series(
    ax: Axes,
    x_cols: Sequence[float],
    series: Sequence[SlopeSeries],
    *,
    min_dy: float = 4.0,
    left_pad: float = 0.05,
    annotation_size: float = 8.4,
) -> None:
    """Draw a multi-series slopegraph: compact left/mid numbers, full right series labels."""
    left_targets = [float(s["values"][0]) for s in series]
    right_targets = [float(s["values"][-1]) for s in series]
    left_relaxed = relax_vertical_positions(left_targets, min_gap=min_dy * 0.75)
    right_relaxed = relax_vertical_positions(right_targets, min_gap=min_dy)
    mid_relaxed = (
        relax_vertical_positions([float(s["values"][1]) for s in series], min_gap=min_dy * 0.78)
        if len(x_cols) >= 3
        else []
    )

    for idx, (s, yl, yr) in enumerate(zip(series, left_relaxed, right_relaxed, strict=True)):
        vals = [float(v) for v in s["values"]]
        col = str(s.get("color", TUFTE_INK))
        lw = float(s.get("lw", 2.0))
        ls = str(s.get("ls", "-"))
        ms = float(s.get("marker_size", 6.0))
        fw = str(s.get("label_weight", "normal"))
        is_focal = fw == "bold"

        name_str = sanitize_text(str(s.get("name", "")))
        right_raw = sanitize_text(str(s.get("right_label", name_str)))
        if name_str and name_str.lower() not in right_raw.lower():
            full_right = f"{name_str}: {right_raw}"
        else:
            full_right = right_raw

        tufte_style.draw_slopegraph_series(
            ax,
            x_positions=x_cols,
            y_values=vals,
            color=col,
            right_label=full_right,
            right_y_relaxed=yr,
            linewidth=lw,
            linestyle=ls,
            markersize=ms,
            is_focal=is_focal,
            label_pad_pt=7.0,
        )
        ax.text(
            x_cols[0] - left_pad,
            yl,
            f"{vals[0]:.1f}%",
            ha="right",
            va="center",
            fontsize=annotation_size,
            fontweight=fw,
            color=col,
        )
        if len(x_cols) >= 3 and mid_relaxed:
            ax.text(
                x_cols[1] - 0.04,
                mid_relaxed[idx],
                f"{vals[1]:.1f}%",
                ha="right",
                va="center",
                fontsize=annotation_size - 0.3,
                color=col,
            )


def draw_tufte_quartile_rule(
    ax: Axes,
    values: Iterable[float],
    y: float = 0.33,
) -> None:
    """Wrapper around ``tufte_style.draw_tufte_quartile_rule``."""
    tufte_style.draw_tufte_quartile_rule(
        ax,
        values=values,
        y_rule=y,
        unit_suffix="%",
        tick_pt=5.0,
        label_pad_pt=5.0,
        label_width_pt=38.0,
    )


def draw_tufte_table(
    ax: Axes,
    headers: Sequence[tuple[float, str, str]],
    rows: Sequence[Sequence[str]],
    *,
    y_top: float = 0.92,
    row_step: float = 0.125,
    font_size: float = 9.0,
    header_size: float = 9.0,
    row_accents: dict[int, str] | None = None,
    col_colors: dict[int, str] | None = None,
    bold_cols: set[int] | None = None,
) -> float:
    """Draw a clean horizontal-hairline Tufte table inside a 0..1 canvas axis."""
    row_accents = row_accents or {}
    col_colors = col_colors or {}
    bold_cols = bold_cols or {0}

    ax.plot([0.0, 1.0], [y_top + 0.045, y_top + 0.045], color=TUFTE_INK, lw=1.0, clip_on=False)
    for x_pos, h_label, ha in headers:
        ax.text(
            x_pos,
            y_top + 0.012,
            sanitize_text(h_label),
            ha=ha,
            va="bottom",
            fontsize=header_size,
            fontweight="bold",
            color=TUFTE_MUTED,
        )
    ax.plot([0.0, 1.0], [y_top - 0.008, y_top - 0.008], color=TUFTE_LIGHT, lw=0.8, clip_on=False)

    y_curr = y_top - 0.035
    for r_idx, row in enumerate(rows):
        accent = row_accents.get(r_idx)
        max_lines = max(cell.count("\n") + 1 for cell in row)
        step = row_step * (1.0 + 0.36 * (max_lines - 1))
        if accent:
            ax.add_patch(
                Rectangle(
                    (0.0, y_curr - step + 0.018),
                    1.0,
                    step - 0.010,
                    facecolor=accent,
                    edgecolor="none",
                    alpha=0.08,
                    zorder=0,
                )
            )
        for c_idx, ((x_pos, _, ha), cell) in enumerate(zip(headers, row, strict=False)):
            c_color = col_colors.get(c_idx, TUFTE_INK)
            if r_idx in row_accents and c_idx == 0:
                c_color = row_accents[r_idx]
            ax.text(
                x_pos,
                y_curr - 0.008,
                sanitize_text(cell),
                ha=ha,
                va="top",
                fontsize=font_size,
                fontweight="bold" if c_idx in bold_cols else "normal",
                color=c_color,
                linespacing=1.20,
            )
        y_curr -= step
        ax.plot(
            [0.0, 1.0],
            [y_curr + 0.014, y_curr + 0.014],
            color=TUFTE_HAIRLINE,
            lw=0.6,
            clip_on=False,
        )

    ax.plot([0.0, 1.0], [y_curr + 0.014, y_curr + 0.014], color=TUFTE_INK, lw=0.9, clip_on=False)
    return y_curr


def draw_marginal_note(
    ax: Axes,
    x: float,
    y_top: float,
    kicker: str,
    body_lines: Sequence[str],
    *,
    accent_color: str = TUFTE_SLATE,
    line_step: float = 0.046,
    font_size: float = 9.0,
) -> float:
    """Draw a borderless Tufte marginal callout with a left vertical accent bar."""
    total_h = 0.042 + len(body_lines) * line_step
    y_bot = y_top - total_h
    ax.plot([x, x], [y_bot + 0.008, y_top], color=accent_color, lw=2.4, solid_capstyle="butt")
    ax.text(
        x + 0.022,
        y_top - 0.002,
        sanitize_text(kicker),
        ha="left",
        va="top",
        fontsize=font_size,
        fontweight="bold",
        color=accent_color,
    )
    y_line = y_top - 0.042
    for line in body_lines:
        ax.text(
            x + 0.022,
            y_line,
            sanitize_text(line),
            ha="left",
            va="top",
            fontsize=font_size - 0.2,
            color=TUFTE_INK,
        )
        y_line -= line_step
    return y_bot


def save_and_verify_plate(
    fig: Figure,
    plate_num: int,
    slug: str,
    out_dir: Path,
) -> PlateInfo:
    """Sanitize all text, check its bounding boxes for out-of-bounds or overlaps, and save.

    Writes PNG, SVG and PDF."""
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"plate_{plate_num:02d}_{slug}"
    stem_path = out_dir / stem

    all_texts: list[Text] = list(fig.texts)
    for ax in fig.axes:
        all_texts.extend(ax.texts)
        all_texts.extend(ax.get_xticklabels())
        all_texts.extend(ax.get_yticklabels())

    for t in all_texts:
        orig = t.get_text()
        clean = sanitize_text(orig)
        if clean != orig:
            t.set_text(clean)

    fig.canvas.draw()
    fig_bbox = fig.bbox
    w_px, h_px = fig_bbox.width, fig_bbox.height

    out_of_bounds: list[str] = []
    texts_info: list[tuple[str, tuple[float, float, float, float]]] = []

    check_texts: list[Text] = list(fig.texts)
    for ax in fig.axes:
        check_texts.extend(ax.texts)

    for t in check_texts:
        if not t.get_visible():
            continue
        s = t.get_text().strip()
        if not s:
            continue
        bb = t.get_window_extent()
        x0, y0, x1, y1 = bb.x0 / w_px, bb.y0 / h_px, bb.x1 / w_px, bb.y1 / h_px
        if x0 < 0.005 or y0 < 0.005 or x1 > 0.995 or y1 > 0.995:
            out_of_bounds.append(f"{s[:30]!r} at ({x0:.3f},{y0:.3f})-({x1:.3f},{y1:.3f})")
        texts_info.append((s, (x0, y0, x1, y1)))

    overlaps: list[str] = []
    for i in range(len(texts_info)):
        s1, (ax0, ay0, ax1, ay1) = texts_info[i]
        for j in range(i + 1, len(texts_info)):
            s2, (bx0, by0, bx1, by1) = texts_info[j]
            ix0, it0 = max(ax0, bx0), max(ay0, by0)
            ix1, it1 = min(ax1, bx1), min(ay1, by1)
            if ix1 - ix0 > 0.004 and it1 - it0 > 0.004:
                overlaps.append(f"{s1[:22]!r} <-> {s2[:22]!r}")

    save_plate(
        fig,
        stem=stem_path,
        formats=("png", "svg", "pdf"),
        dpi=PLATE_DPI,
        svg_fonttype="none",
        salt=stem,
        facecolor=TUFTE_PAPER,
    )
    plt.close(fig)

    return {
        "plate_num": plate_num,
        "stem": stem,
        "png_path": str(out_dir / f"{stem}.png"),
        "svg_path": str(out_dir / f"{stem}.svg"),
        "pdf_path": str(out_dir / f"{stem}.pdf"),
        "out_of_bounds": out_of_bounds,
        "overlaps": overlaps,
    }
