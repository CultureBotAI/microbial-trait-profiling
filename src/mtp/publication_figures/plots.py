"""Publication pages drawn exclusively from verified evidence (#1046).

Rendering never recalculates inference, chooses a significance-selected family,
clips data, samples comparisons, or rewrites verified feature labels.
"""

from __future__ import annotations

import math
import re
import textwrap
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.typing import RcKeyType

if TYPE_CHECKING:
    from mtp.publication_figures.inputs import FigureInputs

WIDTH = 7.2
FONT = 8
INK = "#263746"
ACCENT = "#996B25"
GRAY = "#777777"
STYLE: dict[RcKeyType, Any] = {
    "font.family": "DejaVu Sans",
    "font.size": FONT,
    "axes.labelsize": FONT,
    "axes.titlesize": 9,
    "axes.linewidth": 0.6,
    "xtick.labelsize": FONT,
    "ytick.labelsize": FONT,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "lines.linewidth": 0.7,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "svg.hashsalt": "mtp-publication-figures-v1",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
}


@dataclass(frozen=True)
class Page:
    """One complete page, with the identities and quantities actually drawn."""

    figure: Figure
    row_ids: list[str]
    mark_count: int
    details: dict[str, Any]


def _figure(height: float) -> Figure:
    figure = Figure(figsize=(WIDTH, height))
    FigureCanvasAgg(figure)
    return figure


def _wrap(value: str, width: int) -> str:
    # Wrap long qualified IDs too; never shorten them or discard label characters.
    return "\n".join(
        textwrap.wrap(value, width=width, break_long_words=True, break_on_hyphens=False)
    )


def _number(value: object) -> float | None:
    if value is None or value is pd.NA or value == "":
        return None
    try:
        number = float(str(value))
    except ValueError as error:
        raise ValueError(f"Non-numeric verified figure value: {value!r}") from error
    return number if math.isfinite(number) else None


def _format(value: object) -> str:
    number = _number(value)
    return "unavailable" if number is None else f"{number:.3g}"


def validate_authored_wording(value: str) -> str:
    """Guard these annotation figures, not literal labels or source-study language.

    The rendered contrasts use increase/decrease and stool. They do not measure
    abundance or physiology, and enrichment/depletion belongs to a different
    quantity. Verified trait labels are kept verbatim outside this guard.
    """
    if re.search(
        r"\b(redox|coupled|dysbiosis|causes|physiological|enriched|depleted|"
        r"abundances?|fecal|faecal)\b",
        value,
        re.IGNORECASE,
    ):
        raise ValueError("Unsupported claim in authored figure wording")
    return value


def _count(value: object) -> str:
    number = _number(value)
    if number is None:
        return "unavailable"
    if number < 0 or not number.is_integer():
        raise ValueError("Figure support count must be a nonnegative integer")
    return str(int(number))


def _label(bundle: FigureInputs, trait: str, width: int = 30) -> str:
    return _wrap(bundle.labels[trait], width) + "\n" + _wrap(trait, width + 4)


def _context(bundle: FigureInputs) -> str:
    return (
        f"Rollup: {bundle.record.get('strain_rollup', 'unspecified')}; "
        f"population: {bundle.record.get('contrast_population', 'unspecified')}"
    )


def _records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    return [{str(key): value for key, value in row.items()} for row in frame.to_dict("records")]


def _ranked(frame: pd.DataFrame, field: str, top_n: int) -> tuple[list[dict[str, Any]], int]:
    selected = []
    positive_count = 0
    for positive, suffix in ((True, "descending"), (False, "ascending")):
        rank_name = f"{field}_rank_{suffix}"
        group = []
        for row in _records(frame):
            effect, rank = _number(row[field]), _number(row[rank_name])
            if effect is None or rank is None or rank > top_n:
                continue
            if (positive and effect > 0) or (not positive and effect < 0):
                group.append(row)
        group.sort(key=lambda row: (float(row[rank_name]), str(row["trait"])))
        selected.extend(group)
        if positive:
            positive_count = len(group)
    return selected, positive_count


def forest_pages(
    bundle: FigureInputs, *, ranking: str = "comparison_mean", top_n: int = 5
) -> Iterator[Page]:
    """Select native full-family effect ranks; display both distinct estimands."""
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n <= 0:
        raise ValueError("top_n must be a positive integer")
    fields = {"comparison_mean": "experiment_mean_delta_z", "study_mean": "study_mean_delta_z"}
    if ranking not in fields:
        raise ValueError("ranking must be comparison_mean or study_mean")
    rows, n_positive = _ranked(bundle.parent.evidence, fields[ranking], top_n)
    row_ids = [str(row["trait"]) for row in rows]
    labels = [
        _label(bundle, str(row["trait"]))
        + (
            "\n"
            + _wrap(
                validate_authored_wording(f"Study CI unavailable ({row['study_t_status']})"), 30
            )
            if row["study_t_status"] != "tested"
            else ""
        )
        for row in rows
    ]
    label_lines = [max(2, len(label.splitlines())) for label in labels]
    row_heights = [max(0.42, lines * 0.135) for lines in label_lines]
    # Footer, legend, and axis-label strips use inches, independent of row count.
    height = max(3.5, 1.95 + sum(row_heights))
    with matplotlib.rc_context(STYLE):
        figure = _figure(height)
        axis = figure.add_axes((0.46, 1.05 / height, 0.50, (height - 1.95) / height))
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.axvline(0, color=GRAY, linewidth=0.6, zorder=0)
        positions = np.cumsum(row_heights) - np.asarray(row_heights) / 2
        axis.set_ylim(sum(row_heights) + 0.1, -0.1)
        axis.set_yticks(positions, labels)
        axis.tick_params(axis="y", length=0, pad=8)
        study_x, study_y, comparison_x, comparison_y = [], [], [], []
        lows, highs, interval_y = [], [], []
        for row, position in zip(rows, positions, strict=True):
            study, comparison = (
                _number(row["study_mean_delta_z"]),
                _number(row["experiment_mean_delta_z"]),
            )
            if study is not None:
                study_x.append(study)
                study_y.append(float(position) - 0.04)
            if comparison is not None:
                comparison_x.append(comparison)
                comparison_y.append(float(position) + 0.04)
            lower = _number(row["study_mean_delta_z_ci95_low"])
            upper = _number(row["study_mean_delta_z_ci95_high"])
            if row["study_t_status"] == "tested" and lower is not None and upper is not None:
                lows.append(lower)
                highs.append(upper)
                interval_y.append(float(position) - 0.04)
        axis.hlines(interval_y, lows, highs, color=INK, linewidth=0.7, gid="forest-ci")
        axis.scatter(study_x, study_y, marker="o", s=14, color=INK, gid="forest-study", zorder=3)
        axis.scatter(
            comparison_x,
            comparison_y,
            marker="D",
            s=15,
            facecolors="none",
            edgecolors=ACCENT,
            linewidths=0.7,
            gid="forest-comparison",
            zorder=3,
        )
        axis.set_xlabel(validate_authored_wording("Mean Δz (increase minus decrease)"))
        axis.margins(x=0.08)
        if not rows:
            axis.text(
                0.5,
                0.5,
                validate_authored_wording("No traits meet the declared rank and direction rule"),
                ha="center",
                transform=axis.transAxes,
            )
        title_rank = "comparison" if ranking == "comparison_mean" else "study"
        figure.text(
            0.04,
            1 - 0.22 / height,
            validate_authored_wording(
                "Largest signed shifts in transformed annotation composition"
            ),
            fontsize=10,
            va="top",
        )
        figure.text(
            0.04,
            1 - 0.52 / height,
            _wrap(
                validate_authored_wording(f"Top {top_n} ranks per direction. " + _context(bundle)),
                115,
            ),
            va="top",
        )
        handles = [
            Line2D(
                [],
                [],
                marker="o",
                color=INK,
                markersize=4,
                label=validate_authored_wording("Study mean and 95% CI"),
            ),
            Line2D(
                [],
                [],
                marker="D",
                linestyle="none",
                markerfacecolor="none",
                markeredgecolor=ACCENT,
                markersize=4,
                label=validate_authored_wording("Comparison mean (no interval)"),
            ),
        ]
        figure.legend(
            handles=handles,
            loc="lower left",
            bbox_to_anchor=(0.04, 0.28 / height),
            borderaxespad=0,
            borderpad=0.2,
            frameon=False,
            ncol=2,
        )
        figure.text(
            0.04,
            0.10 / height,
            validate_authored_wording(
                f"Full-family {title_rank}-mean ranks; ties retained. "
                "Δz: standardized log-transformed annotation fractions."
            ),
            va="bottom",
        )
        details = {
            "ranking": ranking,
            "top_n": top_n,
            "selection_rule": (
                "native rank(min) <= top_n, separately positive and negative; zero excluded"
            ),
            "positive_count": n_positive,
            "negative_count": len(rows) - n_positive,
            "retained_family_count": len(bundle.parent.evidence),
            "interval_count": len(lows),
            "interval_estimand": "study_mean_delta_z",
            "units": "standardized_log1p_1000_retained_fraction_difference",
            "point_estimands": ["study_mean_delta_z", "experiment_mean_delta_z"],
        }
        page = Page(figure, row_ids, len(study_x) + len(comparison_x), details)
    yield page


def _extent(values: Sequence[float]) -> tuple[float, float]:
    lower, upper = min([0.0, *values]), max([0.0, *values])
    padding = (upper - lower) * 0.08 if upper != lower else 1.0
    return lower - padding, upper + padding


def quadrant_pages(bundle: FigureInputs) -> Iterator[Page]:
    """Draw every comparison in each declared population, preserving exact coordinates."""
    pairs = bundle.tables["quadrant_pairs.tsv"]
    summary = {
        str(row["population"]): row for row in _records(bundle.tables["quadrant_summary.tsv"])
    }
    stool = pairs.loc[pairs.body_site.eq("Feces")]
    xlim = _extent(pairs.fa_delta_z.astype(float).tolist())
    ylim = _extent(pairs.oa_delta_z.astype(float).tolist())
    all_ids: list[str] = []
    populations: dict[str, Any] = {}
    with matplotlib.rc_context(STYLE):
        # Keep count explanations in a separate inch-based strip, never over observations.
        height = 5.8
        figure = _figure(height)
        for index, (population, frame) in enumerate((("all", pairs), ("stool", stool))):
            row = summary[population]
            axis = figure.add_axes((0.13 + index * 0.49, 2.35 / height, 0.32, 2.3 / height))
            axis.spines[["top", "right"]].set_visible(False)
            axis.axvline(0, color=GRAY, linewidth=0.6, zorder=0)
            axis.axhline(0, color=GRAY, linewidth=0.6, zorder=0)
            axis.scatter(
                frame.fa_delta_z,
                frame.oa_delta_z,
                s=6,
                alpha=0.4,
                color=INK,
                edgecolors="none",
                gid=f"quadrant-{population}",
            )
            axis.set_xlim(*xlim)
            axis.set_ylim(*ylim)
            axis.set_xlabel(
                validate_authored_wording("Δz: ") + _label(bundle, str(row["fa_trait"]), 29)
            )
            axis.set_ylabel(
                validate_authored_wording("Δz: ") + _label(bundle, str(row["oa_trait"]), 27)
            )
            title = "All comparisons" if population == "all" else "Stool (exact Feces)"
            axis.set_title(
                validate_authored_wording(
                    f"{title}\n{int(row['n_comparisons'])} comparisons; "
                    f"{int(row['n_studies'])} studies"
                )
            )
            figure.text(
                0.13 + index * 0.49,
                1.24 / height,
                validate_authored_wording(
                    f"x increase / y increase: {int(row['both_up'])}\n"
                    f"x decrease / y increase: {int(row['fa_down_oa_up'])}\n"
                    f"x decrease / y decrease: {int(row['both_down'])}\n"
                    f"x increase / y decrease: {int(row['fa_up_oa_down'])}\n"
                    f"Axis comparisons (separate): {int(row['axis'])}"
                ),
                va="top",
                gid=f"quadrant-count-{population}",
            )
            all_ids.extend(
                f"{population}:{study}/{experiment}"
                for study, experiment in frame[["study_id", "experiment_id"]].itertuples(
                    index=False, name=None
                )
            )
            populations[population] = {
                "n_comparisons": len(frame),
                "n_studies": int(row["n_studies"]),
                "axis": int(row["axis"]),
                "sign_test_status": row["sign_test_status"],
                "sign_test_p": _number(row["sign_test_p"]),
                "pearson_status": row["pearson_status"],
                "pearson_r": _number(row["pearson_r"]),
            }
        figure.text(
            0.04,
            1.46 / height,
            validate_authored_wording(
                "Comparisons by quadrant (axis excluded); x/y: horizontal/vertical Δz"
            ),
            weight="bold",
            va="top",
            gid="quadrant-count-key",
        )
        figure.text(
            0.04,
            1 - 0.22 / height,
            validate_authored_wording(
                "Oxygen-label contrasts in transformed annotation composition"
            ),
            fontsize=10,
            va="top",
        )
        figure.text(
            0.04,
            1 - 0.55 / height,
            _wrap(validate_authored_wording(_context(bundle)), 115),
            va="top",
        )
        figure.text(
            0.04,
            0.27 / height,
            validate_authored_wording(
                "Δz: standardized log-transformed annotation fractions, increase minus decrease."
            ),
        )
        figure.text(
            0.04,
            0.08 / height,
            validate_authored_wording(
                "Shared denominators constrain correlation; every comparison is shown. "
                "Studies can repeat."
            ),
        )
        page = Page(
            figure,
            all_ids,
            len(all_ids),
            {
                "populations": populations,
                "stool_rule": "body_site == Feces",
                "sampling": "none",
                "jitter": "none",
                "axis_sign_threshold": 0.0,
                "count_placement": "separate_strip_outside_data_axes",
                "count_unit": "comparisons_by_quadrant_axis_excluded_and_reported_separately",
                "units": "standardized_log1p_1000_retained_fraction_difference",
                "x_trait": str(summary["all"]["fa_trait"]),
                "y_trait": str(summary["all"]["oa_trait"]),
                "p_adjustment": "none_nominal",
                "test_unit": "comparison",
            },
        )
    yield page


def _votes(row: dict[str, Any], prefix: str) -> str:
    return " / ".join(
        _count(row[f"{prefix}_studies_{suffix}"]) for suffix in ("up", "down", "tied")
    )


def _cells(bundle: FigureInputs, row: dict[str, Any]) -> list[str]:
    noel = (
        f"Increase / decrease / tied:\n{_votes(row, 'noel')}\n"
        f"BY q = {_format(row['noel_q_by_nominal'])}\n{row['noel_status']}"
    )
    low = _number(row["study_mean_delta_z_ci95_low"])
    high = _number(row["study_mean_delta_z_ci95_high"])
    interval = (
        f"95% CI [{_format(low)}, {_format(high)}]"
        if low is not None and high is not None
        else "95% CI: unavailable"
    )
    study = (
        f"Δz = {_format(row['study_mean_delta_z'])}\n{interval}\n"
        f"BY q = {_format(row['study_t_q_by'])}\n{row['study_t_status']}"
    )
    ora = (
        f"Increase / decrease / tied:\n{_votes(row, 'ora_direction')}\n"
        f"Eligible studies: {_count(row['ora_direction_studies_eligible'])}\n"
        f"{row['ora_direction_status']}"
    )
    return [
        _label(bundle, str(row["trait"]), 29),
        *[
            "\n".join(_wrap(line, width) for line in validate_authored_wording(cell).splitlines())
            for cell, width in ((noel, 24), (study, 25), (ora, 22))
        ],
    ]


def table_pages(bundle: FigureInputs) -> Iterator[Page]:
    """Paginate the complete retained trait family without significance selection."""
    rows = sorted(_records(bundle.parent.evidence), key=lambda row: str(row["trait"]))
    groups: list[list[tuple[dict[str, Any], list[str], float]]] = [[]]
    used = 0.0
    for row in rows:
        cells = _cells(bundle, row)
        height = max(len(cell.splitlines()) for cell in cells) * 0.14 + 0.15
        if used + height > 6.6 and groups[-1]:
            groups.append([])
            used = 0.0
        groups[-1].append((row, cells, height))
        used += height
    for number, group in enumerate(groups, start=1):
        with matplotlib.rc_context(STYLE):
            figure = _figure(9.0)
            figure.text(
                0.04,
                0.97,
                validate_authored_wording(
                    f"Per-trait annotation evidence · page {number} of {len(groups)}"
                ),
                fontsize=10,
                va="top",
            )
            figure.text(
                0.04, 0.935, _wrap(validate_authored_wording(_context(bundle)), 115), va="top"
            )
            columns = [0.04, 0.39, 0.59, 0.81]
            headers = [
                "Verified label / qualified ID",
                "Noel study votes\nFull-family BY",
                "Study mean Δz\n95% CI; full-family BY",
                "ORA direction counts\nUncalibrated",
            ]
            for x, header in zip(columns, headers, strict=True):
                figure.text(x, 0.88, validate_authored_wording(header), va="top", weight="bold")
            y = 0.825
            for row, cells, height in group:
                for column, (x, cell) in enumerate(zip(columns, cells, strict=True)):
                    figure.text(x, y, cell, va="top", gid=f"table-cell:{row['trait']}:{column}")
                y -= height / 9.0
                figure.add_artist(
                    Line2D(
                        [0.04, 0.98],
                        [y + 0.012, y + 0.012],
                        transform=figure.transFigure,
                        color="#dddddd",
                        linewidth=0.5,
                    )
                )
            if not group:
                figure.text(0.04, 0.78, validate_authored_wording("No retained traits"))
            figure.text(
                0.04,
                0.055,
                validate_authored_wording(
                    "Δz: standardized log-transformed annotation fractions, "
                    "increase minus decrease."
                ),
            )
            figure.text(
                0.04,
                0.03,
                validate_authored_wording(
                    "Study means have equal weight. ORA is directional only; "
                    "unavailable values are not zero."
                ),
            )
            ids = [str(row["trait"]) for row, _cells_text, _height in group]
            page = Page(
                figure,
                ids,
                len(ids),
                {
                    "page": number,
                    "page_count": len(groups),
                    "page_row_count": len(ids),
                    "retained_family_count": len(rows),
                    "selection": "complete_retained_family",
                    "ordering": "qualified_trait_id_ascending",
                    "ora_role": "uncalibrated_directional_summary",
                    "units": "standardized_log1p_1000_retained_fraction_difference",
                    "correction": "full_family_BY_in_each_native_family",
                },
            )
        yield page
