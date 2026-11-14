"""Paginated slide views of the three detailed cellulose reports (#1885).

Only saved results are displayed. Pagination preserves rows and quantitative
artists, rather than shrinking text or selecting a smaller scientific result.
The public entry point and export provenance live in ``figures.py``.
"""

from __future__ import annotations

import math
from collections.abc import Iterator, Sequence
from typing import Any
from urllib.parse import quote

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from mtp.cellulose_decr import figures as f
from mtp.cellulose_decr.analysis import ALL_CONTRASTS, ORIENTED, POOLED_STRATIFIED

SIZE = (13.33, 7.5)
NULL_SOURCE = "fig1_pairwise_and_3way_jaccard_nulls"
SITE_SOURCE = "fig2_direction_by_body_site_and_definition"
CORE_SOURCE = "fig3_core_taxa_and_trait_benchmark"
ANNOTATIONS = (
    "Reported signature taxa linked to cellulose annotations, including taxonomic rollup; "
    "not measured cellulose degradation."
)


def _new(title: str, subtitle: str, footer: str) -> Figure:
    figure = Figure(figsize=SIZE)
    figure.text(
        0.035, 0.965, f.wrap(title, 88), fontsize=20, fontweight="bold", color=f.INK, va="top"
    )
    figure.text(0.035, 0.855, f.wrap(subtitle, 142), fontsize=12, color=f.INK2, va="top")
    figure.text(0.035, 0.075, f.wrap(footer, 142), fontsize=12, color=f.INK2, va="bottom")
    return figure


def _text(
    figure: Figure,
    x: float,
    y: float,
    text: str,
    *,
    gid: str | None = None,
    ha: str = "left",
    bold: bool = False,
) -> None:
    figure.text(
        x,
        y,
        text,
        fontsize=12,
        color=f.INK2,
        va="center",
        ha=ha,
        fontweight="bold" if bold else "normal",
        gid=gid,
    )


def _axes(figure: Figure, rectangle: tuple[float, float, float, float]) -> Axes:
    ax = figure.add_axes(rectangle)
    ax.tick_params(labelsize=12, width=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    f.recessive_grid(ax)
    return ax


def _page(
    results: f.Results,
    figure: Figure,
    source: str,
    section: str,
    number: int,
    title: str,
    takeaway: str,
    caveat: str,
    evidence: str,
) -> f.SlidePage:
    population = {
        "compartment_shares": "all contrasts + oriented disease-vs-healthy subject controls",
        "sensitivity": "definition-specific populations shown in the row labels",
        "interactions": "all contrasts; stool vs other sites",
    }.get(section, "all contrasts; stool")
    figure.text(
        0.035,
        0.018,
        f"BugSigDB {results.summary['bugsigdb_vintage']} · human profiles · {population}",
        fontsize=12,
        color=f.INK2,
        va="bottom",
    )
    # Gids keep numerical artists traceable in SVG without invalid whitespace IDs.
    for artist in figure.findobj():
        if artist.get_gid() is not None:
            artist.set_gid(quote(str(artist.get_gid()), safe=":._-"))
    return f.SlidePage(
        f"{source}_{section}_{number:02d}",
        source,
        section,
        figure,
        f.Caption(
            title,
            takeaway,
            f"{caveat} Population: {population}. {f.vintage_note(results)}",
            evidence,
        ),
    )


def _groups(rows: Sequence[dict[str, Any]], size: int) -> Iterator[list[dict[str, Any]]]:
    for start in range(0, len(rows), size):
        yield list(rows[start : start + size])


def _null_table(results: f.Results) -> f.SlidePage:
    title = "Pairwise overlap against each saved null model"
    comparisons = [f"{a} vs. {b}" for a, b in f.PAIRS]
    rows = [
        f.only(results.pairwise[results.pairwise["comparison"].eq(name)]) for name in comparisons
    ]
    rows.append(
        f.only(results.pairwise[results.pairwise["endpoint"].eq("mean_pairwise_jaccard_distance")])
    )
    endpoints: dict[tuple[str, ...], Any] = {}
    footer = (
        "Cells show expected Jaccard distance and nominal p. U1 to U3 pairs use exact "
        "hypergeometric tails; other p-values use Monte Carlo draws. " + ANNOTATIONS
    )
    figure = _new(
        title,
        "Lower Jaccard distance means more overlap; observed values are in the column headings.",
        footer,
    )
    footer_artist = figure.texts[-1]
    xs = (0.40, 0.56, 0.72, 0.88)
    for x, row in zip(xs, rows, strict=True):
        _text(
            figure,
            x,
            0.755,
            f"{row['comparison']}\nobserved {row['jaccard_distance']:.3f}",
            ha="center",
            bold=True,
        )
    labels = dict(f.NULL_ROWS)
    for index, (column, key, _) in enumerate(f.TABLE_NULLS):
        y = 0.665 - index * 0.073
        _text(figure, 0.035, y, f.wrap(f.null_label(results, key, labels[key]), 38))
        for x, row in zip(xs, rows, strict=True):
            comparison = str(row["comparison"])
            exact = column in ("u1", "u2", "u3") and not comparison.startswith("Mean")
            suffix = "hypergeom_p" if exact else "mc_p" if column in ("u1", "u2", "u3") else "p"
            value = row[f"{column}_{suffix}"]
            p_text = (
                f.nominal_p(value)
                if exact
                else f.mc_p(value, results.summary["null_models"][key]["iterations"])
            )
            _text(
                figure,
                x,
                y + 0.016,
                f"{row[f'{column}_expected_jaccard_distance']:.3f}",
                gid=f"null_expected:{comparison}:{column}",
                ha="center",
                bold=True,
            )
            _text(figure, x, y - 0.018, p_text, ha="center")
            endpoints[(key, comparison, "exact" if exact else "monte_carlo")] = value
    scope = f.nominal_scope(endpoints)
    footer_artist.set_text(f.wrap(f"{scope} {footer}", 142))
    return _page(
        results,
        figure,
        NULL_SOURCE,
        "null_table",
        1,
        title,
        "All four comparisons and all seven saved null columns are shown.",
        f"{footer} {scope}",
        "pairwise_jaccard_significance.tsv",
    )


def _null_core(results: f.Results) -> f.SlidePage:
    title = "Taxa shared by all three conditions against the saved nulls"
    observed = int(results.summary["observed"]["three_way_core"])
    rows = [(key, results.summary["null_models"][key]["three_way_core"]) for key, _ in f.NULL_ROWS]
    scope = f.nominal_scope({("null", key, "three_way_core"): row["p"] for key, row in rows})
    figure = _new(
        title,
        f"Observed core: {observed} taxa. Dots: null mean; bars: central 95% of saved draws.",
        f"{scope} {ANNOTATIONS}",
    )
    ax = _axes(figure, (0.45, 0.23, 0.33, 0.53))
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_yticks([])
    ax.set_xlim(0, f.nice_limit(max(observed, *(float(row["p97_5"]) for _, row in rows)) * 1.15, 5))
    ax.set_xlabel("taxa in all three conditions", fontsize=12)
    ax.axvline(observed, color=f.INK, linestyle="--", linewidth=1, gid="null_observed")
    labels = dict(f.NULL_ROWS)
    for index, (key, row) in enumerate(rows):
        y = len(rows) - 1 - index
        ax.plot(
            [row["p2_5"], row["p97_5"]],
            [y, y],
            color=f.MUTED,
            linewidth=2,
            gid=f"null_interval:{key}",
        )
        ax.plot(row["expected"], y, "o", color=f.INK, markersize=6, gid=f"null_mean:{key}")
        position = 0.23 + 0.53 * (y + 0.6) / (len(rows) + 0.2)
        _text(figure, 0.035, position, f.wrap(f.null_label(results, key, labels[key]), 48))
        _text(
            figure,
            0.97,
            position,
            f.mc_p(row["p"], results.summary["null_models"][key]["iterations"]),
            ha="right",
        )
    return _page(
        results,
        figure,
        NULL_SOURCE,
        "null_core",
        1,
        title,
        f"The observed core contains {observed} taxa; every saved null is displayed.",
        f"{scope} {ANNOTATIONS}",
        "cellulose_decr_summary.json",
    )


def _region_pages(results: f.Results) -> Iterator[f.SlidePage]:
    number = 0
    for label, _, names in f.regions(results.catalog):
        chunks = [names[start : start + 14] for start in range(0, len(names), 14)] or [[]]
        for part, chunk in enumerate(chunks, 1):
            number += 1
            title = f"Cellulose-annotated taxa: {label}"
            figure = _new(
                title,
                f"{len(names)} taxa in this set region; page {part} of {len(chunks)}.",
                "Membership means reported in at least one decreased stool signature. "
                + ANNOTATIONS,
            )
            if not chunk:
                _text(figure, 0.035, 0.73, "No taxa in this set region.")
            for index, name in enumerate(chunk):
                _text(
                    figure,
                    0.035,
                    0.755 - index * 0.039,
                    f.wrap(name, 125),
                    gid=f"taxon_name:{name}",
                )
            yield _page(
                results,
                figure,
                NULL_SOURCE,
                "taxon_regions",
                number,
                title,
                f"All {len(names)} taxa in this region are retained across its pages.",
                ANNOTATIONS,
                "cellulose_decreased_taxa_catalog.tsv",
            )


def _share_pages(results: f.Results) -> Iterator[f.SlidePage]:
    series = (
        (ALL_CONTRASTS, "experiment", "o", f.INK, False, "All contrasts, experiments"),
        (ALL_CONTRASTS, "study", "s", f.INK2, False, "All contrasts, studies"),
        (ORIENTED, "experiment", "D", f.INK, True, "Disease vs healthy, experiments"),
        (ORIENTED, "study", "^", f.INK2, True, "Disease vs healthy, studies"),
    )
    for number, compartment in enumerate(f.COMPARTMENT_LABEL, 1):
        title = f"Cellulose annotation direction: {f.COMPARTMENT_LABEL[compartment]}"
        footer = (
            "Ties excluded; exact 95% binomial intervals. "
            "Experiments within studies are not independent. "
            "Right counts: all-contrast experiments. "
            "All contrasts retain reported directions; disease-vs-healthy rows "
            "are oriented subject controls. " + ANNOTATIONS
        )
        figure = _new(
            title,
            "Share of untied pairs with the larger cellulose annotation fraction "
            "in the decreased signature.",
            footer,
        )
        ax = _axes(figure, (0.16, 0.37, 0.61, 0.39))
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, 3.5)
        ax.set_yticks([3, 2, 1, 0], [*f.SHORT, "Pooled"], fontsize=12)
        ax.set_xlabel("share with larger annotation fraction in decreased signatures", fontsize=12)
        ax.axvline(0.5, color=f.INK2, linestyle="--", linewidth=0.8)
        for i, short in enumerate((*f.SHORT, "Pooled")):
            for j, (population, unit, marker, color, hollow, _) in enumerate(series):
                row = f.compartment_row(results, population, compartment, short)
                if row is None:
                    continue
                down, up = f.untied(row, unit)
                share, low, high = f.share_interval(down, up)
                y = 3 - i + 0.24 - j * 0.16
                identity = f"{population}:{compartment}:{short}:{unit}"
                if math.isfinite(share):
                    ax.plot(
                        [low, high],
                        [y, y],
                        color=color,
                        linewidth=1.5,
                        gid=f"share_interval:{identity}",
                    )
                    ax.plot(
                        share,
                        y,
                        marker,
                        color=color,
                        markerfacecolor="white" if hollow else color,
                        markersize=6,
                        gid=f"share_mean:{identity}",
                    )
                if j == 0:
                    ax.text(
                        1.04,
                        y,
                        f"{down} ↓ / {up} ↑" if down + up else "no untied pairs",
                        transform=ax.get_yaxis_transform(),
                        va="center",
                        fontsize=12,
                        color=f.INK2,
                    )
        for index, (_, _, marker, color, hollow, label) in enumerate(series):
            x, y = 0.05 + (index % 2) * 0.48, 0.255 - (index // 2) * 0.035
            figure.add_artist(
                Line2D(
                    [x],
                    [y],
                    transform=figure.transFigure,
                    marker=marker,
                    color=color,
                    markerfacecolor="white" if hollow else color,
                    markersize=6,
                )
            )
            _text(figure, x + 0.02, y, label)
        yield _page(
            results,
            figure,
            SITE_SOURCE,
            "compartment_shares",
            number,
            title,
            "All available experiment and study shares for both populations are retained.",
            footer,
            "compartment_jaccard_and_direction_stats.tsv",
        )


def _sensitivity(results: f.Results) -> f.SlidePage:
    title = "Overlap under every saved set definition"
    footer = (
        "Jaccard distances and shared / union counts are descriptive sensitivity summaries, "
        "not tests. 'Reference' includes every stool-decreased taxon regardless of trait. "
        + ANNOTATIONS
    )
    figure = _new(
        title,
        "All three pair distances, their mean, and the three-condition core / union are shown.",
        footer,
    )
    for x, text in zip(
        (0.50, 0.60, 0.70, 0.80, 0.93),
        ("PD / HIV", "PD / CRC", "HIV / CRC", "Mean", "Core / union"),
        strict=True,
    ):
        _text(figure, x, 0.76, text, ha="center", bold=True)
    for index, (slug, label) in enumerate(f.SENSITIVITY_LABELS.items()):
        y = 0.70 - index * 0.047
        _text(figure, 0.035, y, f.wrap(label, 47))
        part = results.sensitivity[results.sensitivity["definition_slug"].eq(slug)]
        for x, comparison in zip(
            (0.50, 0.60, 0.70, 0.80),
            [f"{a} vs. {b}" for a, b in f.PAIRS] + ["Mean of 3 pairs"],
            strict=True,
        ):
            row = f.only(
                part[part["comparison"].str.startswith("Mean")]
                if comparison.startswith("Mean")
                else part[part["comparison"].eq(comparison)]
            )
            saved_comparison = str(row["comparison"])
            _text(
                figure,
                x,
                y,
                f"{row['jaccard_distance']:.3f}",
                gid=f"sensitivity_value:{slug}:{saved_comparison}",
                ha="center",
            )
            if comparison.startswith("Mean"):
                _text(
                    figure, 0.93, y, f"{int(row['intersection'])}/{int(row['union'])}", ha="center"
                )
    return _page(
        results,
        figure,
        SITE_SOURCE,
        "sensitivity",
        1,
        title,
        "Every saved set definition is displayed without recomputing its distances.",
        footer,
        "sensitivity_definitions_jaccard.tsv",
    )


def _interaction_pages(results: f.Results) -> Iterator[f.SlidePage]:
    selected = results.interaction[
        results.interaction["contrast_population"].eq(ALL_CONTRASTS)
        & results.interaction["short_condition"].isin(["CRC", POOLED_STRATIFIED])
    ]
    rows_per_page = 4
    for number, rows in enumerate(_groups(f.records(selected), rows_per_page), 1):
        title = "Cellulose annotation direction: stool vs other sites"
        scope = f.nominal_scope(dict(f.interaction_endpoint(row) for row in rows))
        footer = (
            "Exact conditional 95% intervals. Arrows mark limits outside 0.01 to 1,000. "
            "Pooled rows estimate a common odds ratio stratified by condition. " + scope
        )
        figure = _new(
            title,
            "Odds ratio of the decreased-direction share; "
            f"saved rows {(number - 1) * rows_per_page + 1} "
            f"to {(number - 1) * rows_per_page + len(rows)}.",
            footer,
        )
        ax = _axes(figure, (0.38, 0.24, 0.29, 0.51))
        ax.set_xscale("log")
        ax.set_xlim(0.01, 1000)
        ax.set_xticks(
            [0.01, 0.1, 1, 10, 100, 1000], ["0.01", "0.1", "1", "10", "100", "1000"], fontsize=12
        )
        ax.minorticks_off()
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_yticks([])
        ax.set_xlabel("odds ratio (log scale)", fontsize=12)
        ax.axvline(1, color=f.INK2, linestyle="--", linewidth=0.8)
        saved_values: list[str] = []
        for index, row in enumerate(rows):
            y = len(rows) - index - 1
            position = 0.24 + 0.51 * (y + 0.6) / (len(rows) + 0.2)
            condition = (
                "All three, stratified" if row["short_condition"] == POOLED_STRATIFIED else "CRC"
            )
            _text(
                figure,
                0.035,
                position,
                f"{condition}\n{f.COMPARTMENT_LABEL[row['compartment_b']]} · {row['unit']}",
            )
            identity = f"{row['short_condition']}:{row['compartment_b']}:{row['unit']}"
            estimate = float(row["conditional_or"])
            note = (
                f"{int(row['a_down'])}↓/{int(row['a_up'])}↑ vs "
                f"{int(row['b_down'])}↓/{int(row['b_up'])}↑\n"
                f"{f.nominal_p(row['fisher_p_two_sided'])}"
            )
            if math.isnan(estimate):
                note += "\n" + f.wrap(f.not_estimable(row), 34)
            else:
                low, high = (
                    float(row["conditional_or_ci95_low"]),
                    float(row["conditional_or_ci95_high"]),
                )
                saved_values.append(
                    f"{condition} / {row['compartment_b']} / {row['unit']}: "
                    f"Saved OR {estimate}; exact 95% interval [{low}, {high}]."
                )
                # Clip both ends independently, including infinite bounds. This
                # preserves a visible segment without reversing wholly off-axis CIs.
                ax.plot(
                    [min(max(low, 0.01), 1000), min(max(high, 0.01), 1000)],
                    [y, y],
                    color=f.MUTED,
                    linewidth=2,
                    gid=f"interaction_interval:{identity}",
                )
                for part, value in (("low", low), ("high", high), ("estimate", estimate)):
                    if value < 0.01 or value > 1000:
                        ax.plot(
                            0.01 if value < 0.01 else 1000,
                            y,
                            "<" if value < 0.01 else ">",
                            color=f.INK if part == "estimate" else f.MUTED,
                            markersize=6,
                            clip_on=False,
                            gid=f"interaction_{part}_limit:{identity}",
                        )
                if estimate == 0 or math.isinf(estimate):
                    note += "\n" + f.wrap(f.zero_cell(row), 34)
                elif 0.01 <= estimate <= 1000:
                    ax.plot(
                        estimate,
                        y,
                        "o",
                        color=f.INK,
                        markersize=6,
                        gid=f"interaction_mean:{identity}",
                    )
            _text(figure, 0.70, position, note)
        yield _page(
            results,
            figure,
            SITE_SOURCE,
            "interactions",
            number,
            title,
            "All saved CRC and condition-stratified interaction rows are retained.",
            f"{footer} {ANNOTATIONS} {' '.join(saved_values)}",
            "compartment_direction_interaction.tsv",
        )


def _shared_pages(results: f.Results) -> Iterator[f.SlidePage]:
    selected = results.catalog[results.catalog["n_conditions_feces_dec"].ge(2)].copy()
    selected["total_dec"] = sum(selected[f"{short}_dec_studies"] for short in f.SHORT)
    selected["block"] = ~selected["is_3way_core"].astype(bool)
    selected = selected.sort_values(
        ["block", "total_dec", "taxon_name"], ascending=[True, False, True]
    )
    for number, rows in enumerate(_groups(f.records(selected), 9), 1):
        title = "Cellulose-annotated taxa reported decreased in at least two conditions"
        footer = (
            "Cells: decreased / increased stool studies. "
            "Flags: core = all three sets; direct = own KG edge; "
            "net↓ = more decreased studies in all three. Rank is shown explicitly. " + ANNOTATIONS
        )
        figure = _new(
            title,
            f"All {len(selected)} shared taxa retained; "
            f"page {number} of {math.ceil(len(selected) / 9)}.",
            footer,
        )
        for x, label in (
            (0.035, "Taxon"),
            (0.42, "Rank"),
            (0.57, "PD"),
            (0.65, "HIV"),
            (0.73, "CRC"),
            (0.81, "Flags"),
        ):
            _text(figure, x, 0.77, label, bold=True)
        for index, row in enumerate(rows):
            y = 0.70 - index * 0.058
            name = str(row["taxon_name"])
            _text(figure, 0.035, y, f.wrap(name, 48), gid=f"shared_name:{name}")
            _text(figure, 0.42, y, str(row["taxonomic_rank"]))
            for x, short in zip((0.57, 0.65, 0.73), f.SHORT, strict=True):
                _text(
                    figure,
                    x,
                    y,
                    f"{int(row[f'{short}_dec_studies'])}/{int(row[f'{short}_inc_studies'])}",
                    gid=f"shared_count:{name}:{short}",
                )
            flags = [
                label
                for column, label in (
                    ("is_3way_core", "core"),
                    ("has_direct_trait_edge", "direct"),
                    ("is_strictly_net_decreased_all3", "net↓"),
                )
                if row[column]
            ]
            _text(figure, 0.81, y, f.wrap(" · ".join(flags) or "—", 20))
        yield _page(
            results,
            figure,
            CORE_SOURCE,
            "shared_taxa",
            number,
            title,
            "Every shared taxon and all three decreased/increased study-count cells are retained.",
            footer,
            "cellulose_decreased_taxa_catalog.tsv",
        )


def _benchmark_pages(results: f.Results) -> Iterator[f.SlidePage]:
    rows = f.records(
        results.benchmark[~results.benchmark["is_growth_medium"].astype(bool)].head(15)
    )
    for number, chunk in enumerate(_groups(rows, 8), 1):
        title = "Fifteen lowest saved non-growth-medium trait distances"
        footer = (
            "Descriptive ranking, not a test: traits share taxa. "
            "Names retain their annotation predicates; "
            "MediaDive growth-medium columns are excluded. " + ANNOTATIONS
        )
        figure = _new(
            title,
            f"Rows {(number - 1) * 8 + 1} to {(number - 1) * 8 + len(chunk)} of {len(rows)}; "
            "smaller distance means more overlap.",
            footer,
        )
        ax = _axes(figure, (0.57, 0.23, 0.32, 0.53))
        ax.set_ylim(-0.6, len(chunk) - 0.4)
        ax.set_xlim(
            0,
            f.nice_limit(
                max(float(row["mean_pairwise_jaccard_distance"]) for row in rows) + 0.1, 0.1
            ),
        )
        ax.set_yticks([])
        ax.set_xlabel("mean pairwise Jaccard distance", fontsize=12)
        for index, row in enumerate(chunk):
            y = len(chunk) - 1 - index
            value = float(row["mean_pairwise_jaccard_distance"])
            color = f.INK if row["trait_id"] == results.summary["trait_id"] else f.INK2
            ax.plot([0, value], [y, y], color=color, linewidth=2)
            ax.plot(
                value, y, "o", color=color, markersize=6, gid=f"benchmark_mean:{row['trait_id']}"
            )
            position = 0.23 + 0.53 * (y + 0.6) / (len(chunk) + 0.2)
            _text(figure, 0.035, position, f.wrap(f.trait_display(row), 64))
            _text(figure, 0.97, position, f"{value:.3f}", ha="right")
        yield _page(
            results,
            figure,
            CORE_SOURCE,
            "trait_benchmark",
            number,
            title,
            "The first fifteen saved non-growth-medium rows are preserved across pages.",
            footer,
            "trait_convergence_benchmark.tsv",
        )


def pages(results: f.Results) -> Iterator[f.SlidePage]:
    """Yield report sections lazily; callers save and release each page in turn."""
    yield _null_table(results)
    yield _null_core(results)
    yield from _region_pages(results)
    yield from _share_pages(results)
    yield _sensitivity(results)
    yield from _interaction_pages(results)
    yield from _shared_pages(results)
    yield from _benchmark_pages(results)
