"""Create static Cytoscape-inspired PNG graph views for the six shared traits."""

from __future__ import annotations

from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, PathPatch, RegularPolygon
from matplotlib.path import Path as MplPath

from create_six_shared_trait_graphs import (
    GRAPH_VIEWS,
    OUTPUT_DIR,
    REPO_ROOT,
    color_from_signal,
    joined_unique_count,
    read_display_edges,
    value_for_color,
)

MANIFEST_FILENAME = "six_shared_trait_static_graph_manifest.tsv"

SLATE_50 = "#f8fafc"
SLATE_200 = "#e2e8f0"
SLATE_400 = "#94a3b8"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
NODE_FILL = "#ffffff"
CONDITION_TRAIT_WIDTH_RANGE = (1.8, 9.8)
TRAIT_BODY_WIDTH_RANGE = (1.6, 10.6)
CONDITION_TRAIT_ALPHA = 0.66
TRAIT_BODY_ALPHA = 0.82
BODY_SITE_CONNECTOR_ALPHA = 0.90


def wrapped_plot_text(value: str, width: int = 150) -> str:
    """Wrap long figure annotations so they do not widen the saved canvas."""
    return "\n".join(wrap(value, width=width))


def wrap_node_label(value: object, width: int = 22) -> str:
    """Wrap a graph-node label with a compact width."""
    text = str(value or "").strip() or "unknown"
    if len(text) > 42:
        text = f"{text[:39].rstrip()}..."
    return "\n".join(wrap(text, width=max(10, width))) or "unknown"


def weighted_signal(group: pd.DataFrame) -> float:
    """Return count-weighted color intensity for an aggregated graph edge."""
    values = group["color_value"].astype(float)
    weights = group["focus_enriched_profile_count"].astype(float)
    if float(weights.sum()) <= 0:
        return float(values.mean()) if not values.empty else 0.0
    return float(np.average(values, weights=weights))


def ordered_node_labels(stats: pd.DataFrame) -> list[str]:
    """Return node labels ordered by trait breadth and positive-z support."""
    return (
        stats.sort_values(
            ["traits", "n", "links", "rank_score"],
            ascending=[False, False, False, False],
        )
        .index.astype(str)
        .tolist()
    )


def evenly_spaced_positions(
    labels: list[str],
    *,
    x: float,
    top: float = 0.82,
    bottom: float = 0.16,
) -> dict[str, tuple[float, float]]:
    """Return node positions along a vertical column."""
    if not labels:
        return {}
    if len(labels) == 1:
        return {labels[0]: (x, 0.50)}
    return {
        label: (x, y)
        for label, y in zip(labels, np.linspace(top, bottom, len(labels)), strict=True)
    }


def scale_values(
    values: pd.Series,
    *,
    output_min: float,
    output_max: float,
) -> pd.Series:
    """Scale positive values with a square-root transform."""
    numeric = values.astype(float)
    positive = numeric[numeric > 0]
    if positive.empty:
        return pd.Series(output_min, index=values.index, dtype=float)
    min_count = float(positive.min())
    max_count = float(positive.max())
    if max_count <= min_count:
        return pd.Series((output_min + output_max) / 2, index=values.index, dtype=float)
    scaled = (np.sqrt(numeric.clip(lower=min_count)) - np.sqrt(min_count)) / (
        np.sqrt(max_count) - np.sqrt(min_count)
    )
    return output_min + scaled.clip(0.0, 1.0) * (output_max - output_min)


def parallel_offsets(
    count: int,
    *,
    spacing: float = 0.009,
    max_offset: float = 0.024,
) -> list[float]:
    """Return top-to-bottom offsets for body-site trait bundle nodes."""
    if count <= 1:
        return [0.0]
    center = (count - 1) / 2
    return [
        min(max_offset, max(-max_offset, (center - position) * spacing))
        for position in range(count)
    ]


def node_stats(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    """Summarize support for condition or body-site nodes."""
    records: list[dict[str, object]] = []
    for label, group in frame.groupby(column, sort=False):
        records.append(
            {
                column: label,
                "links": len(group),
                "traits": group["trait_id"].nunique(),
                "n": int(group["focus_enriched_profile_count"].sum()),
                "s": joined_unique_count(group["focus_enriched_study_ids"]),
                "rank_score": float(group["rank_score"].max()),
            }
        )
    return pd.DataFrame.from_records(records).set_index(column)


def trait_stats(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarize support for the six central trait nodes."""
    records: list[dict[str, object]] = []
    for trait_id, group in frame.groupby("trait_id", sort=False):
        records.append(
            {
                "trait_id": trait_id,
                "trait_label": str(group["trait_label"].iloc[0]),
                "trait_slug": str(group["trait_slug"].iloc[0]),
                "trait_color": str(group["trait_color"].iloc[0]),
                "trait_rank": int(group["trait_rank"].iloc[0]),
                "condition_count": group["condition"].nunique(),
                "body_site_count": group["body_site_group"].nunique(),
                "links": len(group),
                "n": int(group["focus_enriched_profile_count"].sum()),
                "s": joined_unique_count(group["focus_enriched_study_ids"]),
                "rank_score": float(group["rank_score"].max()),
            }
        )
    return (
        pd.DataFrame.from_records(records)
        .sort_values("trait_rank")
        .set_index("trait_id", drop=False)
    )


def body_site_bundle_stats(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarize one accumulated bundle per body-site and trait."""
    records: list[dict[str, object]] = []
    body_site_groups = frame.groupby(
        ["body_site_group", "trait_id", "trait_label", "trait_color", "trait_rank"],
        sort=False,
    )
    for (body_site, trait_id, trait_label, trait_color, trait_rank), group in body_site_groups:
        records.append(
            {
                "body_site_group": body_site,
                "trait_id": trait_id,
                "trait_label": trait_label,
                "trait_color": trait_color,
                "trait_rank": int(trait_rank),
                "bundle_label": f"{trait_label} -> {body_site}",
                "n": int(group["focus_enriched_profile_count"].sum()),
                "bg": int(group["background_no_enriched_profile_count"].sum()),
                "s": joined_unique_count(group["focus_enriched_study_ids"]),
                "links": len(group),
                "signal": weighted_signal(group),
                "rank_score": float(group["rank_score"].max()),
            }
        )
    bundles = pd.DataFrame.from_records(records)
    bundles["bundle_rank"] = 1
    bundles["bundle_offset"] = 0.0
    for _body_site, group in bundles.groupby("body_site_group", sort=False):
        ordered = group.sort_values("trait_rank").index.astype(int).tolist()
        for rank, (index, offset) in enumerate(
            zip(ordered, parallel_offsets(len(ordered)), strict=True),
            start=1,
        ):
            bundles.loc[index, "bundle_rank"] = rank
            bundles.loc[index, "bundle_offset"] = offset
    return bundles


def aggregate_graph_edges(frame: pd.DataFrame, bundles: pd.DataFrame) -> pd.DataFrame:
    """Aggregate displayed Sankey edges into context and body-site bundle edges."""
    records: list[dict[str, object]] = []
    condition_groups = frame.groupby(
        ["condition", "trait_id", "trait_label", "trait_color", "trait_rank"],
        sort=False,
    )
    for (condition, trait_id, trait_label, trait_color, trait_rank), group in condition_groups:
        records.append(
            {
                "edge_type": "condition_trait",
                "source_label": condition,
                "target_label": trait_label,
                "body_site_group": "",
                "trait_id": trait_id,
                "trait_label": trait_label,
                "trait_color": trait_color,
                "trait_rank": int(trait_rank),
                "n": int(group["focus_enriched_profile_count"].sum()),
                "bg": int(group["background_no_enriched_profile_count"].sum()),
                "s": joined_unique_count(group["focus_enriched_study_ids"]),
                "links": len(group),
                "signal": weighted_signal(group),
                "rank_score": float(group["rank_score"].max()),
            }
        )
    for row in bundles.itertuples(index=False):
        records.append(
            {
                "body_site_group": row.body_site_group,
                "trait_id": row.trait_id,
                "trait_label": row.trait_label,
                "trait_color": row.trait_color,
                "trait_rank": int(row.trait_rank),
                "n": int(row.n),
                "bg": int(row.bg),
                "s": int(row.s),
                "links": int(row.links),
                "signal": float(row.signal),
                "rank_score": float(row.rank_score),
                "edge_type": "trait_body_bundle",
                "source_label": row.trait_label,
                "target_label": row.bundle_label,
            }
        )
        records.append(
            {
                "body_site_group": row.body_site_group,
                "trait_id": row.trait_id,
                "trait_label": row.trait_label,
                "trait_color": row.trait_color,
                "trait_rank": int(row.trait_rank),
                "n": int(row.n),
                "bg": int(row.bg),
                "s": int(row.s),
                "links": int(row.links),
                "signal": float(row.signal),
                "rank_score": float(row.rank_score),
                "edge_type": "bundle_body_site",
                "source_label": row.bundle_label,
                "target_label": row.body_site_group,
            }
        )
    return pd.DataFrame.from_records(records)


def assign_edge_widths(graph_edges: pd.DataFrame) -> pd.DataFrame:
    """Scale source links and accumulated body-site edges separately."""
    graph_edges = graph_edges.copy()
    graph_edges["width"] = 0.0
    width_ranges = {
        "condition_trait": CONDITION_TRAIT_WIDTH_RANGE,
        "trait_body_bundle": TRAIT_BODY_WIDTH_RANGE,
        "bundle_body_site": TRAIT_BODY_WIDTH_RANGE,
    }
    for edge_type, (output_min, output_max) in width_ranges.items():
        mask = graph_edges["edge_type"] == edge_type
        if bool(mask.any()):
            graph_edges.loc[mask, "width"] = scale_values(
                graph_edges.loc[mask, "n"],
                output_min=output_min,
                output_max=output_max,
            )
    return graph_edges


def draw_ordered_edges(graph_edges: pd.DataFrame) -> pd.DataFrame:
    """Return graph edges in back-to-front drawing order."""
    edge_order = {
        "condition_trait": 0,
        "trait_body_bundle": 1,
        "bundle_body_site": 2,
    }
    return (
        graph_edges.assign(_edge_order=graph_edges["edge_type"].map(edge_order).fillna(99))
        .sort_values(["_edge_order", "n"], ascending=[True, True])
        .drop(columns="_edge_order")
    )


def body_site_bundle_positions(
    bundles: pd.DataFrame,
    body_site_positions: dict[str, tuple[float, float]],
    *,
    x: float = 0.84,
) -> dict[tuple[str, str], tuple[float, float]]:
    """Return positions for body-site trait bundles just before body-site nodes."""
    return {
        (str(row.body_site_group), str(row.trait_id)): (
            x,
            body_site_positions[str(row.body_site_group)][1] + float(row.bundle_offset),
        )
        for row in bundles.itertuples(index=False)
    }


def bundle_row_key(row: object) -> tuple[str, str]:
    """Return the position lookup key for a bundle-like row."""
    return (str(row.body_site_group), str(row.trait_id))


def trait_positions(traits: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """Return stable vertically aligned positions for central trait hub nodes."""
    ordered_traits = traits.sort_values("trait_rank")
    y_values = [0.50] if len(ordered_traits) == 1 else np.linspace(0.82, 0.18, len(ordered_traits))
    return {
        row.trait_label: (0.51, float(y))
        for row, y in zip(ordered_traits.itertuples(index=False), y_values, strict=True)
    }


def text_color_for_fill(color: str) -> str:
    """Return readable text color for a filled trait node."""
    r, g, b = to_rgb(color)
    luminance = (0.2126 * r) + (0.7152 * g) + (0.0722 * b)
    return "white" if luminance < 0.47 else SLATE_900


def draw_curve(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str,
    linewidth: float,
    alpha: float,
    zorder: float,
    bend: float,
) -> None:
    """Draw one Cytoscape-like curved graph edge."""
    x0, y0 = start
    x1, y1 = end
    dx = x1 - x0
    control_offset = bend * abs(dx)
    path = MplPath(
        [
            (x0, y0),
            (x0 + dx * 0.45, y0 + control_offset),
            (x1 - dx * 0.45, y1 + control_offset),
            (x1, y1),
        ],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    ax.add_patch(
        PathPatch(
            path,
            facecolor="none",
            edgecolor=color,
            linewidth=linewidth,
            alpha=alpha,
            capstyle="round",
            joinstyle="round",
            zorder=zorder,
        )
    )


def draw_circle_node(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    radius: float,
    facecolor: str,
    edgecolor: str,
    linewidth: float = 1.2,
    zorder: float = 5,
) -> None:
    """Draw a circular graph node."""
    ax.add_patch(
        Circle(
            (x, y),
            radius,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
        )
    )


def draw_diamond_node(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    radius: float,
    facecolor: str,
    edgecolor: str,
    linewidth: float = 0.8,
    zorder: float = 5,
) -> None:
    """Draw a compact diamond node for accumulated trait/body-site bundles."""
    ax.add_patch(
        RegularPolygon(
            (x, y),
            numVertices=4,
            radius=radius,
            orientation=np.pi / 4,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
        )
    )


def render_static_graph(
    *,
    view_key: str,
    view_label: str,
    condition_mode: str,
    focus_direction: str,
    display_edges: pd.DataFrame,
    output_path: Path,
    edge_table_path: Path,
    node_table_path: Path,
    shared_signal_max: float,
) -> dict[str, object]:
    """Render one static graph PNG and write its audit node/edge tables."""
    frame = display_edges.copy()
    frame["color_value"] = frame.apply(value_for_color, axis=1)
    condition_frame = node_stats(frame, "condition")
    body_site_frame = node_stats(frame, "body_site_group")
    trait_frame = trait_stats(frame)
    bundles = body_site_bundle_stats(frame)
    bundles["color"] = [
        color_from_signal(row.trait_color, float(row.signal), shared_signal_max)
        for row in bundles.itertuples(index=False)
    ]
    graph_edges = aggregate_graph_edges(frame, bundles)
    graph_edges["color"] = [
        color_from_signal(row.trait_color, float(row.signal), shared_signal_max)
        for row in graph_edges.itertuples(index=False)
    ]
    graph_edges = assign_edge_widths(graph_edges)

    condition_order = ordered_node_labels(condition_frame)
    body_site_order = ordered_node_labels(body_site_frame)
    condition_positions = evenly_spaced_positions(condition_order, x=0.15)
    body_site_positions = evenly_spaced_positions(body_site_order, x=0.88)
    center_positions = trait_positions(trait_frame)
    bundle_positions = body_site_bundle_positions(bundles, body_site_positions)

    node_records: list[dict[str, object]] = []
    for label, row in condition_frame.iterrows():
        x, y = condition_positions[str(label)]
        node_records.append(
            {
                "node_type": "condition",
                "label": label,
                "x": x,
                "y": y,
                "n": int(row.n),
                "s": int(row.s),
                "links": int(row.links),
                "traits": int(row.traits),
            }
        )
    for row in trait_frame.itertuples(index=False):
        x, y = center_positions[row.trait_label]
        node_records.append(
            {
                "node_type": "trait",
                "label": row.trait_label,
                "trait_id": row.trait_id,
                "x": x,
                "y": y,
                "n": int(row.n),
                "s": int(row.s),
                "links": int(row.links),
                "condition_count": int(row.condition_count),
                "body_site_count": int(row.body_site_count),
            }
        )
    for row in bundles.itertuples(index=False):
        x, y = bundle_positions[bundle_row_key(row)]
        node_records.append(
            {
                "node_type": "trait_body_bundle",
                "label": row.bundle_label,
                "body_site_group": row.body_site_group,
                "trait_id": row.trait_id,
                "trait_label": row.trait_label,
                "x": x,
                "y": y,
                "n": int(row.n),
                "s": int(row.s),
                "links": int(row.links),
                "bundle_rank": int(row.bundle_rank),
            }
        )
    for label, row in body_site_frame.iterrows():
        x, y = body_site_positions[str(label)]
        node_records.append(
            {
                "node_type": "body_site",
                "label": label,
                "x": x,
                "y": y,
                "n": int(row.n),
                "s": int(row.s),
                "links": int(row.links),
                "traits": int(row.traits),
            }
        )
    node_columns = [
        "node_type",
        "label",
        "body_site_group",
        "trait_id",
        "trait_label",
        "condition_count",
        "body_site_count",
        "bundle_rank",
        "traits",
        "x",
        "y",
        "n",
        "s",
        "links",
    ]
    pd.DataFrame.from_records(node_records).to_csv(
        node_table_path,
        sep="\t",
        index=False,
        columns=node_columns,
    )
    graph_edges.to_csv(edge_table_path, sep="\t", index=False)

    fig_height = max(11.0, 4.0 + 0.36 * max(len(condition_order), len(body_site_order)))
    fig, ax = plt.subplots(figsize=(18.0, fig_height))
    fig.subplots_adjust(left=0.035, right=0.965, top=0.90, bottom=0.12)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")

    for edge in draw_ordered_edges(graph_edges).itertuples(index=False):
        if edge.edge_type == "condition_trait":
            start = condition_positions[str(edge.source_label)]
            end = center_positions[str(edge.target_label)]
            alpha = CONDITION_TRAIT_ALPHA
            zorder = 1.0
        elif edge.edge_type == "trait_body_bundle":
            start = center_positions[str(edge.source_label)]
            end = bundle_positions[bundle_row_key(edge)]
            alpha = TRAIT_BODY_ALPHA
            zorder = 2.4
        elif edge.edge_type == "bundle_body_site":
            start = bundle_positions[bundle_row_key(edge)]
            end = body_site_positions[str(edge.target_label)]
            alpha = BODY_SITE_CONNECTOR_ALPHA
            zorder = 2.8
        else:
            continue
        bend = (
            0.0 if edge.edge_type == "bundle_body_site" else (float(edge.trait_rank) - 2.5) * 0.010
        )
        draw_curve(
            ax,
            start,
            end,
            color=str(edge.color),
            linewidth=float(edge.width),
            alpha=alpha,
            zorder=zorder,
            bend=bend,
        )

    condition_sizes = scale_values(condition_frame["n"], output_min=0.010, output_max=0.019)
    body_site_sizes = scale_values(body_site_frame["n"], output_min=0.011, output_max=0.021)
    trait_sizes = scale_values(trait_frame["n"], output_min=0.031, output_max=0.043)
    bundle_sizes = scale_values(bundles["n"], output_min=0.0055, output_max=0.0115)

    for label in condition_order:
        x, y = condition_positions[label]
        row = condition_frame.loc[label]
        draw_circle_node(
            ax,
            x=x,
            y=y,
            radius=float(condition_sizes.loc[label]),
            facecolor=NODE_FILL,
            edgecolor=SLATE_500,
            zorder=4.0,
        )
        ax.text(
            x - 0.022,
            y,
            (
                f"{wrap_node_label(label)}\n"
                f"links={int(row.links)}, traits={int(row.traits)}, n={int(row.n):,}"
            ),
            ha="right",
            va="center",
            fontsize=7.2,
            color=SLATE_900,
            linespacing=1.03,
        )

    for index, row in bundles.iterrows():
        x, y = bundle_positions[(str(row.body_site_group), str(row.trait_id))]
        draw_diamond_node(
            ax,
            x=x,
            y=y,
            radius=float(bundle_sizes.loc[index]),
            facecolor=str(row.color),
            edgecolor="none",
            linewidth=0.0,
            zorder=5.4,
        )

    for label in body_site_order:
        x, y = body_site_positions[label]
        row = body_site_frame.loc[label]
        draw_circle_node(
            ax,
            x=x,
            y=y,
            radius=float(body_site_sizes.loc[label]),
            facecolor=NODE_FILL,
            edgecolor=SLATE_500,
            zorder=4.0,
        )
        ax.text(
            x + 0.022,
            y,
            (
                f"{wrap_node_label(label)}\n"
                f"links={int(row.links)}, traits={int(row.traits)}, n={int(row.n):,}"
            ),
            ha="left",
            va="center",
            fontsize=7.2,
            color=SLATE_900,
            linespacing=1.03,
        )

    for row in trait_frame.itertuples(index=False):
        x, y = center_positions[row.trait_label]
        radius = float(trait_sizes.loc[row.trait_id])
        draw_circle_node(
            ax,
            x=x,
            y=y,
            radius=radius,
            facecolor=row.trait_color,
            edgecolor=SLATE_900,
            linewidth=1.4,
            zorder=6.0,
        )
        ax.text(
            x,
            y,
            wrap_node_label(row.trait_label, width=15),
            ha="center",
            va="center",
            fontsize=7.6,
            color=text_color_for_fill(row.trait_color),
            fontweight="bold",
            linespacing=1.02,
            zorder=7.0,
        )
        ax.text(
            x,
            y - radius - 0.018,
            (
                f"n={int(row.n):,}\n"
                f"conditions={int(row.condition_count)}, sites={int(row.body_site_count)}"
            ),
            ha="center",
            va="top",
            fontsize=6.8,
            color=SLATE_700,
            linespacing=1.05,
            zorder=7.0,
        )

    fig.text(
        0.035,
        0.965,
        f"Six shared traits Cytoscape-inspired graph: {focus_direction} profiles",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.035,
        0.925,
        wrapped_plot_text(
            f"Graph view from displayed Sankey edges ({view_label}). Left nodes are "
            f"conditions using {condition_mode}; center nodes are the six q <= 0.001 "
            "trait hubs. Same-trait links accumulate into small colored bundle nodes "
            "adjacent to each body site, then connect directly to the corresponding "
            "body-site node. Width uses separate scales for condition-to-trait and "
            "trait-to-body-site links. n = positive-z profile memberships; "
            "s = unique studies; bg = z <= 0 "
            "background profiles in the audit TSV.",
            width=156,
        ),
        ha="left",
        va="top",
        fontsize=9.5,
        color=SLATE_700,
        linespacing=1.15,
    )

    ax.text(0.15, 0.93, "Condition nodes", ha="center", va="center", fontsize=11.5, color=SLATE_700)
    ax.text(0.51, 0.93, "Six trait hubs", ha="center", va="center", fontsize=11.5, color=SLATE_700)
    ax.text(
        0.90,
        0.93,
        "Body-site nodes",
        ha="center",
        va="center",
        fontsize=11.5,
        color=SLATE_700,
    )

    source_edges = graph_edges[graph_edges["edge_type"] == "condition_trait"]
    bundle_edges = graph_edges[graph_edges["edge_type"] == "trait_body_bundle"]
    source_width_examples = sorted(
        {
            int(source_edges["n"].min()),
            int(source_edges["n"].median()),
            int(source_edges["n"].max()),
        }
    )
    bundle_width_examples = sorted(
        {
            int(bundle_edges["n"].min()),
            int(bundle_edges["n"].median()),
            int(bundle_edges["n"].max()),
        }
    )
    source_width_lookup = scale_values(
        pd.Series(source_width_examples, index=source_width_examples),
        output_min=CONDITION_TRAIT_WIDTH_RANGE[0],
        output_max=CONDITION_TRAIT_WIDTH_RANGE[1],
    )
    bundle_width_lookup = scale_values(
        pd.Series(bundle_width_examples, index=bundle_width_examples),
        output_min=TRAIT_BODY_WIDTH_RANGE[0],
        output_max=TRAIT_BODY_WIDTH_RANGE[1],
    )
    legend_handles = [
        Line2D(
            [0],
            [0],
            color=SLATE_700,
            linewidth=float(source_width_lookup.loc[count]),
            alpha=0.72,
            solid_capstyle="round",
            label=f"condition-trait n={count}",
        )
        for count in source_width_examples
    ] + [
        Line2D(
            [0],
            [0],
            color=SLATE_700,
            linewidth=float(bundle_width_lookup.loc[count]),
            alpha=0.80,
            solid_capstyle="round",
            label=f"trait-body-site n={count}",
        )
        for count in bundle_width_examples
    ]
    legend = ax.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.50, -0.085),
        ncol=3,
        frameon=False,
        fontsize=8,
        title="Separate width scales: condition-to-trait vs. trait-to-body-site",
        title_fontsize=9,
    )
    ax.add_artist(legend)

    fig.text(
        0.035,
        0.035,
        wrapped_plot_text(
            f"Static graph PNG generated from {view_key}. Edge TSV aggregates the "
            "displayed Sankey links into condition-trait context edges, accumulated "
            "trait-to-bundle links, and short bundle-to-body-site connectors.",
            width=155,
        ),
        ha="left",
        va="bottom",
        fontsize=8.0,
        color=SLATE_500,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_output_path = output_path.with_name(f".{output_path.stem}.tmp{output_path.suffix}")
    try:
        fig.savefig(temp_output_path, dpi=240, bbox_inches="tight", facecolor="white")
        temp_output_path.replace(output_path)
    finally:
        if temp_output_path.exists():
            temp_output_path.unlink()
    plt.close(fig)

    return {
        "view_key": view_key,
        "view_label": view_label,
        "focus_direction": focus_direction,
        "condition_mode": condition_mode,
        "condition_node_count": len(condition_order),
        "trait_node_count": len(trait_frame),
        "trait_body_bundle_count": len(bundles),
        "body_site_node_count": len(body_site_order),
        "graph_edge_count": len(graph_edges),
        "condition_trait_edge_count": int((graph_edges["edge_type"] == "condition_trait").sum()),
        "trait_body_bundle_edge_count": int(
            (graph_edges["edge_type"] == "trait_body_bundle").sum()
        ),
        "bundle_body_site_edge_count": int((graph_edges["edge_type"] == "bundle_body_site").sum()),
        "png_path": output_path.relative_to(REPO_ROOT).as_posix(),
        "node_table_path": node_table_path.relative_to(REPO_ROOT).as_posix(),
        "edge_table_path": edge_table_path.relative_to(REPO_ROOT).as_posix(),
    }


def append_readme_section(output_dir: Path, manifest: pd.DataFrame) -> None:
    """Append or replace the static graph PNG section in the README."""
    readme_path = output_dir / "README.md"
    if not readme_path.exists():
        return
    marker = "## Static Cytoscape-Inspired PNG Graphs"
    content = readme_path.read_text(encoding="utf-8").split(marker, maxsplit=1)[0].rstrip()
    lines = [
        "",
        "",
        marker,
        "",
        f"- Manifest: `{MANIFEST_FILENAME}`",
        "",
        "These PNGs are graph versions of the six shared trait Sankey views. "
        "They use six central trait hubs rather than condition-to-body-site flow "
        "bands. Trait-to-body-site links accumulate into small colored bundle nodes "
        "adjacent to each body site, then connect directly to the corresponding "
        "body-site node. Condition-to-trait and trait-to-body-site widths use "
        "separate scales.",
        "",
        (
            "| view | PNG | condition nodes | trait hubs | body-site nodes | "
            "trait/body-site bundles | condition-trait edges | trait-body links | "
            "body-site connectors |"
        ),
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in manifest.itertuples(index=False):
        lines.append(
            f"| {row.view_label} | `{row.png_path}` | {row.condition_node_count} | "
            f"{row.trait_node_count} | {row.body_site_node_count} | "
            f"{row.trait_body_bundle_count} | {row.condition_trait_edge_count} | "
            f"{row.trait_body_bundle_edge_count} | {row.bundle_body_site_edge_count} |"
        )
    lines.append("")
    readme_path.write_text(content + "\n".join(lines), encoding="utf-8")


def write_six_shared_trait_static_graphs(
    output_dir: Path = OUTPUT_DIR,
) -> pd.DataFrame:
    """Write Cytoscape-inspired static graph PNGs and audit tables."""
    frames = [read_display_edges(output_dir / view.display_edges_file) for view in GRAPH_VIEWS]
    signals = pd.concat(
        [frame.apply(value_for_color, axis=1) for frame in frames],
        ignore_index=True,
    )
    shared_signal_max = max(0.4, float(signals.max() if not signals.empty else 0.0))

    records: list[dict[str, object]] = []
    for view, frame in zip(GRAPH_VIEWS, frames, strict=True):
        stem = f"six_shared_traits_{view.key}_cytoscape_hub_graph"
        records.append(
            render_static_graph(
                view_key=view.key,
                view_label=view.label,
                condition_mode=view.condition_mode,
                focus_direction=view.focus_direction,
                display_edges=frame,
                output_path=output_dir / f"{stem}.png",
                edge_table_path=output_dir / f"{stem}_edges.tsv",
                node_table_path=output_dir / f"{stem}_nodes.tsv",
                shared_signal_max=shared_signal_max,
            )
        )
    manifest = pd.DataFrame.from_records(records)
    manifest.to_csv(output_dir / MANIFEST_FILENAME, sep="\t", index=False)
    append_readme_section(output_dir, manifest)
    return manifest


def main() -> None:
    """Write static Cytoscape-inspired graph outputs."""
    manifest = write_six_shared_trait_static_graphs()
    print(f"Wrote static graph PNGs: {len(manifest)}")
    print(f"Wrote manifest: {OUTPUT_DIR / MANIFEST_FILENAME}")


if __name__ == "__main__":
    main()
