"""Create all-BugSigDB condition/body-site maps by profile direction."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from textwrap import shorten, wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path as MplPath

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"
OUTPUT_DIR = INVESTIGATION_DIR / "all-bugsigdb-condition-body-site-maps"
METADATA_PATH = VIS_DIR / "study_trait_profile_study_metadata.tsv"

MIN_FOCUS_PROFILES = 5
MIN_FOCUS_STUDIES = 2
MAX_DISPLAY_EDGES = 24
MAX_DISPLAY_CONDITIONS = 16
MAX_LABEL_WIDTH = 38
MAX_OTHER_SITES_PER_CONDITION = 2

DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}
DISPLAY_EXCLUDED_CONDITION_LABELS = frozenset(
    {
        "Body fat distribution",
        "Clinical treatment",
        "Exposure",
        "Health study participation",
        "Health trait",
        "Microbiome",
        "Physical activity",
        "Population",
        "Timepoint",
    }
)

SLATE_50 = "#f8fafc"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
EDGE_NEUTRAL = "#475569"


@dataclass(frozen=True)
class DirectionMapConfig:
    """Configuration for one all-data direction map."""

    key: str
    focus_direction: str
    comparator_direction: str
    title: str
    output_stem: str

    @property
    def focus_label(self) -> str:
        return DIRECTION_LABELS[self.focus_direction]

    @property
    def comparator_label(self) -> str:
        return DIRECTION_LABELS[self.comparator_direction]

    @property
    def edge_tsv_path(self) -> Path:
        return OUTPUT_DIR / f"{self.output_stem}_edges.tsv"

    @property
    def edge_v2_tsv_path(self) -> Path:
        return OUTPUT_DIR / f"{self.output_stem}_v2_edges.tsv"

    @property
    def png_path(self) -> Path:
        return OUTPUT_DIR / f"{self.output_stem}.png"

    @property
    def png_v2_path(self) -> Path:
        return OUTPUT_DIR / f"{self.output_stem}_v2.png"


CONFIGS = (
    DirectionMapConfig(
        key="increased",
        focus_direction="inc",
        comparator_direction="dec",
        title="All BugSigDB Condition-to-Body-Site Map: Increased Profiles",
        output_stem="all_bugsigdb_increased_profile_condition_body_site_map",
    ),
    DirectionMapConfig(
        key="decreased",
        focus_direction="dec",
        comparator_direction="inc",
        title="All BugSigDB Condition-to-Body-Site Map: Decreased Profiles",
        output_stem="all_bugsigdb_decreased_profile_condition_body_site_map",
    ),
)


def clean_text(value: object) -> str:
    """Normalize missing display values."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def compact_label(value: object, width: int = MAX_LABEL_WIDTH) -> str:
    """Return a compact one-line label."""
    return shorten(clean_text(value), width=width, placeholder="...")


def wrap_label(value: object, width: int = MAX_LABEL_WIDTH) -> str:
    """Wrap a node label for display."""
    label = compact_label(value, width=width)
    return "\n".join(wrap(label, width=max(12, width // 2))) or "unknown"


def join_unique(values: pd.Series) -> str:
    """Join unique non-empty values in first-seen order."""
    seen: list[str] = []
    for value in values:
        text = clean_text(value)
        if text != "unknown" and text not in seen:
            seen.append(text)
    return "; ".join(seen)


def split_joined_values(value: object) -> list[str]:
    """Split a semicolon-delimited value emitted by join_unique."""
    text = clean_text(value)
    if text == "unknown":
        return []
    return [part.strip() for part in text.split(";") if part.strip()]


def is_display_body_site(value: object) -> bool:
    """Return whether a body-site group is useful as a displayed anatomical node."""
    return clean_text(value).lower() not in {"unknown", "other"}


def is_display_condition(value: object) -> bool:
    """Return whether a condition label is specific enough for the displayed map."""
    label = clean_text(value)
    lower_label = label.lower()
    broad_terms = ("exposure", "measurement", "population")
    return (
        label != "unknown"
        and label not in DISPLAY_EXCLUDED_CONDITION_LABELS
        and not any(term in lower_label for term in broad_terms)
    )


def load_profile_metadata(path: Path = METADATA_PATH) -> pd.DataFrame:
    """Load all study/experiment/direction profiles used by the trait-profile figures."""
    metadata = pd.read_csv(path, sep="\t", index_col=0)
    metadata = metadata.copy()
    metadata["profile_id"] = metadata.index
    for column in ["study_id", "experiment_id", "inc_dec", "condition", "body_site_group"]:
        metadata[column] = metadata[column].map(clean_text)
    metadata["direction_label"] = (
        metadata["inc_dec"].map(DIRECTION_LABELS).fillna(metadata["inc_dec"])
    )
    if "body_site" in metadata:
        metadata["body_site"] = metadata["body_site"].map(clean_text)
    else:
        metadata["body_site"] = "unknown"
    return metadata[metadata["inc_dec"].isin(DIRECTION_LABELS)].copy()


def direction_count(group: pd.DataFrame, direction: str, column: str) -> int:
    """Count unique values within one profile direction."""
    return int(group.loc[group["inc_dec"] == direction, column].nunique())


def direction_ids(group: pd.DataFrame, direction: str, column: str) -> str:
    """Join unique values within one profile direction."""
    return join_unique(group.loc[group["inc_dec"] == direction, column])


def summarize_condition_body_sites(metadata: pd.DataFrame) -> pd.DataFrame:
    """Summarize all BugSigDB profile support by condition, body site, and direction."""
    records: list[dict[str, object]] = []
    for (condition, body_site_group), group in metadata.groupby(
        ["condition", "body_site_group"],
        sort=False,
        dropna=False,
    ):
        increased = group[group["inc_dec"] == "inc"]
        decreased = group[group["inc_dec"] == "dec"]
        profile_count = len(group)
        increased_profile_count = len(increased)
        decreased_profile_count = len(decreased)
        records.append(
            {
                "condition": clean_text(condition),
                "body_site_group": clean_text(body_site_group),
                "profile_count": profile_count,
                "study_count": group["study_id"].nunique(),
                "study_ids": join_unique(group["study_id"]),
                "experiment_count": group["experiment_id"].nunique(),
                "body_sites": join_unique(group["body_site"]),
                "increased_profile_count": increased_profile_count,
                "increased_study_count": direction_count(group, "inc", "study_id"),
                "increased_study_ids": direction_ids(group, "inc", "study_id"),
                "decreased_profile_count": decreased_profile_count,
                "decreased_study_count": direction_count(group, "dec", "study_id"),
                "decreased_study_ids": direction_ids(group, "dec", "study_id"),
                "increased_profile_fraction": increased_profile_count / profile_count,
                "decreased_profile_fraction": decreased_profile_count / profile_count,
                "direction_balance_increased_minus_decreased": (
                    increased_profile_count - decreased_profile_count
                )
                / profile_count,
                "display_eligible_condition": is_display_condition(condition),
                "display_eligible_body_site": is_display_body_site(body_site_group),
            }
        )

    summary = pd.DataFrame.from_records(records)
    summary["display_eligible"] = (
        summary["display_eligible_condition"] & summary["display_eligible_body_site"]
    )
    return summary.sort_values(
        ["profile_count", "study_count", "condition", "body_site_group"],
        ascending=[False, False, True, True],
    )


def build_direction_edges(
    summary: pd.DataFrame,
    config: DirectionMapConfig,
) -> pd.DataFrame:
    """Build one display-candidate edge table for a profile direction."""
    edges = summary[summary["display_eligible"]].copy()
    focus_label = config.focus_label
    comparator_label = config.comparator_label
    edges["focus_direction"] = focus_label
    edges["comparator_direction"] = comparator_label
    edges["focus_profile_count"] = edges[f"{focus_label}_profile_count"]
    edges["focus_study_count"] = edges[f"{focus_label}_study_count"]
    edges["focus_study_ids"] = edges[f"{focus_label}_study_ids"]
    edges["comparator_profile_count"] = edges[f"{comparator_label}_profile_count"]
    edges["comparator_study_count"] = edges[f"{comparator_label}_study_count"]
    edges["comparator_study_ids"] = edges[f"{comparator_label}_study_ids"]
    edges = edges[edges["focus_profile_count"] > 0].copy()
    edges["focus_profile_fraction"] = edges["focus_profile_count"] / edges["profile_count"]
    direction_balance = (edges["focus_profile_count"] - edges["comparator_profile_count"]) / edges[
        "profile_count"
    ]
    edges["focus_minus_comparator_profile_balance"] = direction_balance
    edges["rank_score"] = (
        edges["focus_profile_count"]
        * np.log1p(edges["focus_study_count"])
        * (1.0 + np.maximum(direction_balance, 0.0))
    )
    edges["displayed_in_map"] = False
    sort_columns = ["rank_score", "focus_profile_count", "focus_study_count", "condition"]
    return edges.sort_values(sort_columns, ascending=[False, False, False, True])


def select_display_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Choose a readable subset while preserving body-site diversity."""
    candidates = edges[
        (edges["focus_profile_count"] >= MIN_FOCUS_PROFILES)
        & (edges["focus_study_count"] >= MIN_FOCUS_STUDIES)
    ].copy()
    if candidates.empty:
        candidates = edges.head(MAX_DISPLAY_EDGES).copy()

    candidates = candidates.sort_values(
        ["rank_score", "focus_profile_count", "focus_study_count", "condition"],
        ascending=[False, False, False, True],
    )
    per_body_site = candidates.groupby("body_site_group", sort=False).head(2)
    selected_indices: list[int] = []
    selected_conditions: set[str] = set()

    def add_candidate(index: int, row: pd.Series) -> None:
        if index in selected_indices or len(selected_indices) >= MAX_DISPLAY_EDGES:
            return
        condition = str(row["condition"])
        if (
            condition not in selected_conditions
            and len(selected_conditions) >= MAX_DISPLAY_CONDITIONS
        ):
            return
        selected_indices.append(index)
        selected_conditions.add(condition)

    for index, row in per_body_site.iterrows():
        add_candidate(index, row)
    for index, row in candidates.iterrows():
        add_candidate(index, row)

    return candidates.loc[selected_indices].sort_values(
        ["rank_score", "focus_profile_count", "focus_study_count", "condition"],
        ascending=[False, False, False, True],
    )


def select_uncapped_display_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Return every supported edge for the uncapped v2 PNG."""
    candidates = edges[
        (edges["focus_profile_count"] >= MIN_FOCUS_PROFILES)
        & (edges["focus_study_count"] >= MIN_FOCUS_STUDIES)
    ].copy()
    return candidates.sort_values(
        ["rank_score", "focus_profile_count", "focus_study_count", "condition"],
        ascending=[False, False, False, True],
    )


def scaled_line_width(values: pd.Series) -> pd.Series:
    """Scale profile counts to readable line widths."""
    counts = values.astype(float)
    if counts.max() == counts.min():
        return pd.Series(4.0, index=values.index)
    scaled = (np.sqrt(counts) - np.sqrt(counts.min())) / (
        np.sqrt(counts.max()) - np.sqrt(counts.min())
    )
    return 1.1 + scaled * 8.0


def node_positions(
    labels: list[str],
    *,
    top: float = 0.88,
    bottom: float = 0.14,
) -> dict[str, float]:
    """Return evenly spaced vertical positions for node labels."""
    if len(labels) == 1:
        return {labels[0]: 0.52}
    positions = np.linspace(top, bottom, len(labels))
    return dict(zip(labels, positions, strict=True))


def count_joined_ids(group: pd.DataFrame, column: str, fallback_column: str) -> int:
    """Count unique identifiers from joined-id columns."""
    values: set[str] = set()
    if column in group:
        for value in group[column]:
            values.update(split_joined_values(value))
    if values:
        return len(values)
    if fallback_column in group and not group[fallback_column].empty:
        return int(group[fallback_column].max())
    return 0


def weighted_node_stats(edges: pd.DataFrame, node_column: str) -> pd.DataFrame:
    """Summarize node support over displayed edges."""
    rows: list[dict[str, object]] = []
    for label, group in edges.groupby(node_column, sort=False):
        focus_count = int(group["focus_profile_count"].sum())
        total_count = int(group["profile_count"].sum())
        balance = float(
            np.average(
                group["focus_minus_comparator_profile_balance"].to_numpy(dtype=float),
                weights=group["profile_count"].to_numpy(dtype=float),
            )
        )
        rows.append(
            {
                node_column: label,
                "focus_profile_count": focus_count,
                "profile_count": total_count,
                "focus_study_count": count_joined_ids(
                    group,
                    "focus_study_ids",
                    "focus_study_count",
                ),
                "focus_profile_fraction": focus_count / total_count if total_count else 0.0,
                "focus_minus_comparator_profile_balance": balance,
            }
        )
    return pd.DataFrame.from_records(rows).set_index(node_column)


def draw_curve(
    ax: plt.Axes,
    *,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    linewidth: float,
) -> None:
    """Draw one left-to-right Bezier edge."""
    x0, y0 = start
    x1, y1 = end
    control_dx = (x1 - x0) * 0.46
    path = MplPath(
        [
            (x0, y0),
            (x0 + control_dx, y0),
            (x1 - control_dx, y1),
            (x1, y1),
        ],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    patch = PathPatch(
        path,
        facecolor="none",
        edgecolor=color,
        linewidth=linewidth,
        alpha=0.72,
        capstyle="round",
        joinstyle="round",
    )
    ax.add_patch(patch)


def build_sampled_site_notes(display_edges: pd.DataFrame, summary: pd.DataFrame) -> dict[str, str]:
    """List sampled body sites for displayed conditions that are not drawn."""
    notes: dict[str, str] = {}
    displayed_by_condition = {
        condition: set(group["body_site_group"])
        for condition, group in display_edges.groupby("condition", sort=False)
    }
    for condition, displayed_sites in displayed_by_condition.items():
        condition_sites = summary[
            (summary["condition"] == condition)
            & summary["display_eligible_body_site"]
            & ~summary["body_site_group"].isin(displayed_sites)
        ].copy()
        if condition_sites.empty:
            continue
        condition_sites = condition_sites.sort_values(
            ["profile_count", "study_count", "body_site_group"],
            ascending=[False, False, True],
        )
        entries = [
            (
                f"{compact_label(row.body_site_group, width=18)} "
                f"n={int(row.profile_count)} s={int(row.study_count)}"
            )
            for row in condition_sites.head(MAX_OTHER_SITES_PER_CONDITION).itertuples()
        ]
        remaining = len(condition_sites) - len(entries)
        if remaining:
            entries.append(f"+{remaining} more")
        notes[condition] = "other sampled sites: " + "; ".join(entries)
    return notes


def draw_node(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    label: str,
    stats_row: pd.Series,
    align: str,
    note: str = "",
) -> None:
    """Draw a node marker and its text label."""
    ax.add_patch(
        Circle(
            (x, y),
            0.012,
            facecolor="white",
            edgecolor=SLATE_700,
            linewidth=1.0,
            zorder=5,
        )
    )
    text_x = x - 0.025 if align == "right" else x + 0.025
    label_lines = [
        wrap_label(label),
        (
            f"n={int(stats_row.focus_profile_count)}, "
            f"s={int(stats_row.focus_study_count)}, "
            f"focus={stats_row.focus_profile_fraction:.0%}"
        ),
    ]
    if note:
        label_lines.extend(wrap(note, width=42))
    ax.text(
        text_x,
        y,
        "\n".join(label_lines),
        ha=align,
        va="center",
        fontsize=7.3 if note else 8.2,
        color=SLATE_900,
        linespacing=1.08,
    )


def plot_direction_map(
    display_edges: pd.DataFrame,
    summary: pd.DataFrame,
    config: DirectionMapConfig,
    output_path: Path,
    *,
    title_suffix: str = "",
    selection_note: str | None = None,
    dense_layout: bool = False,
) -> None:
    """Render an all-BugSigDB direction map."""
    if display_edges.empty:
        msg = f"No display edges available for {config.key}."
        raise ValueError(msg)

    condition_order = (
        display_edges.groupby("condition", sort=False)["rank_score"]
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    body_site_order = (
        display_edges.groupby("body_site_group", sort=False)["focus_profile_count"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )
    node_top = 0.98 if dense_layout else 0.88
    node_bottom = 0.05 if dense_layout else 0.14
    header_y = 0.985 if dense_layout else 0.93
    condition_y = node_positions(condition_order, top=node_top, bottom=node_bottom)
    body_site_y = node_positions(body_site_order, top=node_top, bottom=node_bottom)
    condition_stats = weighted_node_stats(display_edges, "condition")
    body_site_stats = weighted_node_stats(display_edges, "body_site_group")
    condition_notes = build_sampled_site_notes(display_edges, summary)

    height = max(
        9.0,
        1.4 + 0.62 * max(len(condition_order), len(body_site_order)) + 0.10 * len(condition_notes),
    )
    fig, ax = plt.subplots(figsize=(14.6, height))
    if dense_layout:
        fig.subplots_adjust(left=0.05, right=0.90, top=0.91, bottom=0.07)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    cmap = plt.get_cmap("RdBu_r")
    widths = scaled_line_width(display_edges["focus_profile_count"])

    left_x = 0.34
    right_x = 0.68
    for edge, linewidth in zip(display_edges.itertuples(), widths, strict=True):
        draw_curve(
            ax,
            start=(left_x, condition_y[edge.condition]),
            end=(right_x, body_site_y[edge.body_site_group]),
            color=cmap(norm(edge.focus_minus_comparator_profile_balance)),
            linewidth=float(linewidth),
        )

    for condition in condition_order:
        draw_node(
            ax,
            x=left_x,
            y=condition_y[condition],
            label=condition,
            stats_row=condition_stats.loc[condition],
            align="right",
            note=condition_notes.get(condition, ""),
        )
    for body_site in body_site_order:
        draw_node(
            ax,
            x=right_x,
            y=body_site_y[body_site],
            label=body_site,
            stats_row=body_site_stats.loc[body_site],
            align="left",
        )

    fig.text(
        0.05,
        0.965,
        f"{config.title}{title_suffix}",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
    )
    fig.text(
        0.05,
        0.927,
        (
            "Edges summarize all BugSigDB study/experiment/direction profiles after label "
            f"filtering. Width = {config.focus_label}-profile support; color = "
            f"({config.focus_label} - {config.comparator_label}) / all profiles; "
            "s = unique studies."
        ),
        ha="left",
        va="top",
        fontsize=10,
        color=SLATE_700,
    )
    ax.text(0.21, header_y, "Conditions", ha="center", va="center", fontsize=12, color=SLATE_700)
    ax.text(
        0.79,
        header_y,
        "Body-site groups",
        ha="center",
        va="center",
        fontsize=12,
        color=SLATE_700,
    )

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.030, pad=0.02, shrink=0.55)
    cbar.set_label(
        f"{config.focus_label} profile balance",
        fontsize=9,
    )
    cbar.ax.tick_params(labelsize=8)

    representative_counts = sorted(
        {
            int(display_edges["focus_profile_count"].min()),
            int(display_edges["focus_profile_count"].median()),
            int(display_edges["focus_profile_count"].max()),
        }
    )
    width_lookup = scaled_line_width(pd.Series(representative_counts, index=representative_counts))
    legend_lines = [
        Line2D(
            [0],
            [0],
            color=EDGE_NEUTRAL,
            linewidth=float(width_lookup.loc[count]),
            alpha=0.72,
            solid_capstyle="round",
            label=f"n={count}",
        )
        for count in representative_counts
    ]
    legend = ax.legend(
        handles=legend_lines,
        title=f"{config.focus_label} profiles",
        loc="lower center",
        bbox_to_anchor=(0.50, -0.015),
        ncol=len(legend_lines),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    ax.add_artist(legend)

    if selection_note is None:
        selection_note = (
            f"Displayed edges require at least {MIN_FOCUS_PROFILES} {config.focus_label} "
            f"profiles and {MIN_FOCUS_STUDIES} focus-direction studies, then are ranked by "
            "profile support, study support, and direction balance. Other/unknown body sites "
            "and broad metadata-like condition labels remain in the summary TSV but are not "
            "drawn in the PNG."
        )

    fig.text(
        0.05,
        0.035,
        selection_note,
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=SLATE_500,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_readme(summary: pd.DataFrame, edge_tables: dict[str, pd.DataFrame]) -> None:
    """Write a short methods README for the all-data direction maps."""
    metadata_path = METADATA_PATH.relative_to(REPO_ROOT).as_posix()
    metadata_profile_count = int(summary["profile_count"].sum())
    study_count = len(
        {study_id for value in summary["study_ids"] for study_id in split_joined_values(value)}
    )
    lines = [
        "# All BugSigDB Condition/Body-Site Direction Maps",
        "",
        (
            "These maps summarize all BugSigDB study/experiment/direction profiles "
            f"available in `{metadata_path}`. They do not use a selected trait set or "
            "trait-enrichment score; they are direction-support maps for increased "
            "and decreased profiles."
        ),
        "",
        f"Profiles summarized: {metadata_profile_count:,}",
        f"Studies summarized: {study_count:,}",
        "",
        "## Files",
        "",
        (
            "- `all_bugsigdb_condition_body_site_direction_summary.tsv`: all "
            "condition/body-site pairs, including broad labels and `Other`/`unknown` "
            "body-site groups."
        ),
    ]
    for config in CONFIGS:
        edges = edge_tables[config.key]
        displayed = int(edges["displayed_in_map"].sum())
        lines.extend(
            [
                f"- `{config.output_stem}.png`: displayed {config.focus_label}-profile map.",
                f"- `{config.output_stem}_edges.tsv`: display-candidate edge table ",
                f"  for the {config.focus_label}-profile map; {displayed} rows are drawn.",
                (f"- `{config.output_stem}_v2.png`: uncapped v2 {config.focus_label}-profile map."),
                (
                    f"- `{config.output_stem}_v2_edges.tsv`: uncapped v2 edge table "
                    f"for the {config.focus_label}-profile map."
                ),
            ]
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            (
                "Each row of the source profile metadata is one BugSigDB "
                "study/experiment/direction profile. Edges connect condition labels to "
                "grouped body sites. Edge width is the number of profiles in the focus "
                "direction. Edge color is the direction balance: `(focus profiles - "
                "comparator profiles) / all profiles` for that condition/body-site pair. "
                "Node labels use `n=` for focus-profile count, `s=` for unique "
                "focus-direction study count, and `focus=` for the focus-direction "
                "fraction among profiles represented by displayed edges."
            ),
            "",
            (
                "The TSV preserves all sampled condition/body-site pairs. The PNG is "
                "intentionally filtered for readability: broad metadata-like condition "
                "labels and `Other`/`unknown` body-site groups are kept in the TSV but "
                "not drawn."
            ),
            "",
            (
                "The `_v2` PNGs remove the readability caps on displayed edges, "
                "conditions, body-site groups, and per-body-site edge selection. They "
                "still use the same support threshold and label filters as the displayed "
                "maps."
            ),
        ]
    )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Write all-data direction maps and companion TSVs."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    metadata = load_profile_metadata()
    summary = summarize_condition_body_sites(metadata)
    edge_tables: dict[str, pd.DataFrame] = {}

    for config in CONFIGS:
        edges = build_direction_edges(summary, config)
        display_edges = select_display_edges(edges)
        display_edges_v2 = select_uncapped_display_edges(edges)
        display_keys = set(
            zip(display_edges["condition"], display_edges["body_site_group"], strict=True)
        )
        display_v2_keys = set(
            zip(display_edges_v2["condition"], display_edges_v2["body_site_group"], strict=True)
        )
        edge_output = edges.copy()
        edge_output["displayed_in_map"] = [
            (condition, body_site) in display_keys
            for condition, body_site in zip(
                edge_output["condition"], edge_output["body_site_group"], strict=True
            )
        ]
        edge_v2_output = edges.copy()
        edge_v2_output["displayed_in_map"] = [
            (condition, body_site) in display_v2_keys
            for condition, body_site in zip(
                edge_v2_output["condition"], edge_v2_output["body_site_group"], strict=True
            )
        ]
        edge_tables[config.key] = edge_output
        edge_output.to_csv(config.edge_tsv_path, sep="\t", index=False)
        edge_v2_output.to_csv(config.edge_v2_tsv_path, sep="\t", index=False)
        plot_direction_map(display_edges, summary, config, config.png_path)
        plot_direction_map(
            display_edges_v2,
            summary,
            config,
            config.png_v2_path,
            title_suffix=" (v2 uncapped)",
            selection_note=(
                f"V2 uncapped display: edges require at least {MIN_FOCUS_PROFILES} "
                f"{config.focus_label} profiles and {MIN_FOCUS_STUDIES} focus-direction "
                "studies. Every edge meeting those criteria is drawn, with no edge, "
                "condition, body-site, or per-body-site display cap. Other/unknown body "
                "sites and broad metadata-like condition labels remain in the summary TSV "
                "but are not drawn in the PNG."
            ),
            dense_layout=True,
        )
        print(f"Wrote: {config.edge_tsv_path}")
        print(f"Wrote: {config.png_path}")
        print(f"Wrote: {config.edge_v2_tsv_path}")
        print(f"Wrote: {config.png_v2_path}")

        summary[f"displayed_in_{config.key}_map"] = [
            (condition, body_site) in display_keys
            for condition, body_site in zip(
                summary["condition"], summary["body_site_group"], strict=True
            )
        ]
        summary[f"displayed_in_{config.key}_v2_map"] = [
            (condition, body_site) in display_v2_keys
            for condition, body_site in zip(
                summary["condition"], summary["body_site_group"], strict=True
            )
        ]

    summary.to_csv(
        OUTPUT_DIR / "all_bugsigdb_condition_body_site_direction_summary.tsv",
        sep="\t",
        index=False,
    )
    write_readme(summary, edge_tables)
    print(f"Wrote: {OUTPUT_DIR / 'all_bugsigdb_condition_body_site_direction_summary.tsv'}")
    print(f"Wrote: {OUTPUT_DIR / 'README.md'}")


if __name__ == "__main__":
    main()
