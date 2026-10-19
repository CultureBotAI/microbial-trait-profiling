"""Create condition-level bar plots from profile trait enrichment signals."""

from __future__ import annotations

import argparse
from pathlib import Path
from textwrap import shorten

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"
DEFAULT_MATRIX_PATH = VIS_DIR / "study_trait_profile_clustered_heatmap_matrix.tsv"
DEFAULT_METADATA_PATH = VIS_DIR / "study_trait_profile_study_metadata.tsv"
DEFAULT_OUTPUT_DIR = INVESTIGATION_DIR / "profile-condition-bars"

PROFILE_DIRECTION_LABELS = {
    "increased": "inc",
    "decreased": "dec",
}
DIRECTION_LABELS = {
    "dec": "decreased profiles",
    "inc": "increased profiles",
}
DIRECTION_COLORS = {
    "dec": "#e74c3c",
    "inc": "#3498db",
}
SIGNAL_COLORS = {
    "dec": "#7c3aed",
    "inc": "#d97706",
}
LOW_COVERAGE_STUDY_THRESHOLD = 2
FOCUS_OUTPUT_LABELS = {
    "dec": "depleted",
    "inc": "increased",
}
EXCLUDED_CONDITION_LABELS = {
    "Age",
    "Aging",
    "Alcohol drinking",
    "Antimicrobial agent",
    "Body fat distribution",
    "Clinical treatment",
    "Delivery method",
    "Diet",
    "Environmental factor",
    "Ethnic group",
    "Exposure",
    "Health study participation",
    "Health trait",
    "Microbiome",
    "Population",
    "Physical activity",
    "Response to diet",
    "Timepoint",
}
EXCLUDED_CONDITION_TERMS = (
    "diet",
    "exposure",
    "measurement",
    "placebo",
    "pollution",
    "population",
    "protocol",
    "sampling",
    "smoking behavior",
)


def clean_text(value: object) -> str:
    """Normalize display strings."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return "unknown"
    return text


def compact_label(value: object, max_length: int = 42) -> str:
    """Return a compact label for horizontal bar plots."""
    return shorten(clean_text(value).replace("\n", " "), width=max_length, placeholder="...")


def profile_id_from_matrix_label(value: object) -> str:
    """Convert matrix display row labels into stable profile IDs."""
    parts = [part.strip() for part in clean_text(value).split(" | ")]
    if len(parts) != 3:
        msg = f"Could not parse profile matrix row label: {value!r}"
        raise ValueError(msg)
    direction = PROFILE_DIRECTION_LABELS.get(parts[2].lower())
    if direction is None:
        msg = f"Unknown profile direction in matrix row label: {value!r}"
        raise ValueError(msg)
    return f"{parts[0]}|{parts[1]}|{direction}"


def direction_prefix(direction: str) -> str:
    """Return the summary column prefix for a profile direction."""
    if direction == "dec":
        return "decreased"
    if direction == "inc":
        return "increased"
    msg = f"Unsupported profile direction: {direction!r}"
    raise ValueError(msg)


def comparator_direction(direction: str) -> str:
    """Return the opposite profile direction."""
    if direction == "dec":
        return "inc"
    if direction == "inc":
        return "dec"
    msg = f"Unsupported profile direction: {direction!r}"
    raise ValueError(msg)


def focus_output_label(direction: str) -> str:
    """Return the filename/title label for a focus direction."""
    return FOCUS_OUTPUT_LABELS[direction]


def is_excluded_condition_label(value: object) -> bool:
    """Return whether a condition label is broad metadata rather than a condition."""
    label = clean_text(value)
    lower_label = label.lower()
    return (
        label in EXCLUDED_CONDITION_LABELS
        or label == "unknown"
        or any(term in lower_label for term in EXCLUDED_CONDITION_TERMS)
    )


def load_profile_matrix(matrix_path: Path) -> pd.DataFrame:
    """Load the selected profile z-score matrix with stable profile IDs."""
    matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
    matrix.index = [profile_id_from_matrix_label(index_value) for index_value in matrix.index]
    return matrix


def load_profile_scores(matrix_path: Path, metadata_path: Path) -> pd.DataFrame:
    """Load profile z-scores and metadata, then count positive trait signals."""
    matrix = load_profile_matrix(matrix_path)
    positive_trait_count = matrix.gt(0).sum(axis=1).rename("positive_selected_trait_count")

    metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
    scored = metadata.join(positive_trait_count, how="inner").copy()
    scored["condition"] = scored["condition"].map(clean_text)
    scored["inc_dec"] = scored["inc_dec"].map(clean_text)
    scored["direction_label"] = scored["inc_dec"].map(DIRECTION_LABELS).fillna(scored["inc_dec"])
    scored["has_positive_selected_trait"] = scored["positive_selected_trait_count"] > 0
    scored["raw_trait_count"] = pd.to_numeric(scored["trait_count"], errors="coerce").fillna(0)
    scored["total_trait_signal"] = pd.to_numeric(
        scored["total_trait_signal"],
        errors="coerce",
    ).fillna(0)
    return scored


def filter_scored_profiles(scored: pd.DataFrame, body_site: str | None) -> pd.DataFrame:
    """Optionally restrict profile rows to a body-site group or exact body-site label."""
    if body_site is None:
        return scored

    requested = clean_text(body_site).lower()
    body_site_columns = [
        column for column in ("body_site_group", "body_site") if column in scored.columns
    ]
    if not body_site_columns:
        msg = "Metadata does not contain body_site_group or body_site columns."
        raise ValueError(msg)

    mask = pd.Series(False, index=scored.index)
    for column in body_site_columns:
        normalized = scored[column].map(lambda value: clean_text(value).lower())
        mask = mask | normalized.eq(requested)

    filtered = scored[mask].copy()
    if filtered.empty:
        available_labels: set[str] = set()
        for column in body_site_columns:
            available_labels.update(scored[column].map(clean_text).dropna().unique())
        available = ", ".join(sorted(available_labels)[:30])
        msg = f"No profiles matched body site {body_site!r}. Available examples: {available}"
        raise ValueError(msg)
    return filtered


def summarize_conditions(
    scored: pd.DataFrame,
    *,
    min_profiles: int,
    min_studies: int,
    include_broad_labels: bool,
) -> pd.DataFrame:
    """Summarize positive trait signals by condition and profile direction."""
    selected = scored.copy()
    if not include_broad_labels:
        selected = selected[~selected["condition"].map(is_excluded_condition_label)].copy()

    records: list[dict[str, object]] = []
    for condition, group in selected.groupby("condition", sort=True, dropna=False):
        profile_count = len(group)
        study_count = group["study_id"].nunique()
        if profile_count < min_profiles or study_count < min_studies:
            continue

        record: dict[str, object] = {
            "condition": condition,
            "profile_count": profile_count,
            "study_count": study_count,
            "positive_profile_count": int(group["has_positive_selected_trait"].sum()),
            "positive_profile_fraction": float(group["has_positive_selected_trait"].mean()),
            "mean_positive_selected_trait_count": group["positive_selected_trait_count"].mean(),
            "median_positive_selected_trait_count": group["positive_selected_trait_count"].median(),
            "mean_raw_trait_count": group["raw_trait_count"].mean(),
            "mean_total_trait_signal": group["total_trait_signal"].mean(),
        }

        total_direction_mean = 0.0
        for direction in ["dec", "inc"]:
            direction_group = group[group["inc_dec"] == direction]
            prefix = direction_prefix(direction)
            direction_profile_count = len(direction_group)
            direction_positive_count = int(direction_group["has_positive_selected_trait"].sum())
            direction_mean = (
                float(direction_group["positive_selected_trait_count"].mean())
                if direction_profile_count
                else 0.0
            )
            total_direction_mean += direction_mean
            record[f"{prefix}_profile_count"] = direction_profile_count
            record[f"{prefix}_positive_profile_count"] = direction_positive_count
            record[f"{prefix}_positive_profile_fraction"] = (
                direction_positive_count / direction_profile_count
                if direction_profile_count
                else 0.0
            )
            record[f"{prefix}_mean_positive_selected_trait_count"] = direction_mean
            record[f"{prefix}_median_positive_selected_trait_count"] = (
                float(direction_group["positive_selected_trait_count"].median())
                if direction_profile_count
                else 0.0
            )

        record["total_direction_mean_positive_selected_trait_count"] = total_direction_mean
        records.append(record)

    if not records:
        return pd.DataFrame()
    summary = pd.DataFrame(records)
    for focus_direction in ["dec", "inc"]:
        focus_prefix = direction_prefix(focus_direction)
        comparator_prefix = direction_prefix(comparator_direction(focus_direction))
        summary[f"{focus_prefix}_minus_{comparator_prefix}_mean_positive_selected_trait_count"] = (
            summary[f"{focus_prefix}_mean_positive_selected_trait_count"]
            - summary[f"{comparator_prefix}_mean_positive_selected_trait_count"]
        )
        summary[f"{focus_prefix}_minus_{comparator_prefix}_positive_profile_count"] = (
            summary[f"{focus_prefix}_positive_profile_count"]
            - summary[f"{comparator_prefix}_positive_profile_count"]
        )
        summary[f"{focus_prefix}_minus_{comparator_prefix}_positive_profile_fraction"] = (
            summary[f"{focus_prefix}_positive_profile_fraction"]
            - summary[f"{comparator_prefix}_positive_profile_fraction"]
        )
        positive_signal_profiles = (
            summary[f"{focus_prefix}_positive_profile_count"]
            + summary[f"{comparator_prefix}_positive_profile_count"]
        )
        summary[focus_coverage_weighted_signal_column(focus_direction)] = (
            summary[f"{focus_prefix}_minus_{comparator_prefix}_mean_positive_selected_trait_count"]
            * positive_signal_profiles.pow(0.5)
            * summary["study_count"].pow(0.5)
        )

    return summary.sort_values(
        [
            "total_direction_mean_positive_selected_trait_count",
            "positive_profile_count",
            "condition",
        ],
        ascending=[False, False, True],
    )


def save_stacked_condition_bar(
    summary: pd.DataFrame,
    *,
    output_path: Path,
    value_columns: tuple[str, str],
    total_column: str,
    title: str,
    xlabel: str,
    top_n: int,
    label_format: str,
) -> None:
    """Save one horizontal stacked condition bar plot."""
    if summary.empty:
        return

    plot_data = summary.sort_values(total_column, ascending=False).head(top_n).iloc[::-1].copy()
    labels = [compact_label(condition, 44) for condition in plot_data["condition"]]
    decreased = plot_data[value_columns[0]]
    increased = plot_data[value_columns[1]]

    sns.set_theme(style="whitegrid")
    fig_height = max(5.8, 0.44 * len(plot_data) + 1.8)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    ax.barh(
        labels,
        decreased,
        color=DIRECTION_COLORS["dec"],
        alpha=0.84,
        label=DIRECTION_LABELS["dec"],
    )
    ax.barh(
        labels,
        increased,
        left=decreased,
        color=DIRECTION_COLORS["inc"],
        alpha=0.84,
        label=DIRECTION_LABELS["inc"],
    )
    ax.set_title(title, fontsize=15, fontweight="bold", pad=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("")
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8)
    ax.grid(axis="y", visible=False)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=9)
    ax.tick_params(axis="x", labelsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#e2e8f0")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.03), frameon=False, ncol=2)

    max_total = max(float(plot_data[total_column].max()), 1.0)
    for y_index, row in enumerate(plot_data.itertuples(index=False)):
        total = float(getattr(row, total_column))
        ax.text(
            total + max_total * 0.015,
            y_index,
            label_format.format(total),
            va="center",
            ha="left",
            fontsize=8.4,
            color="#334155",
            fontweight="bold",
        )
    ax.set_xlim(right=max_total * 1.18)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def focus_delta_column(focus_direction: str) -> str:
    """Return the mean positive-trait delta column for a focus direction."""
    focus_prefix = direction_prefix(focus_direction)
    comparator_prefix = direction_prefix(comparator_direction(focus_direction))
    return f"{focus_prefix}_minus_{comparator_prefix}_mean_positive_selected_trait_count"


def focus_coverage_weighted_signal_column(focus_direction: str) -> str:
    """Return the coverage-weighted focus signal column."""
    return f"{direction_prefix(focus_direction)}_coverage_weighted_signal"


def focus_delta_axis_label(focus_direction: str) -> str:
    """Return a readable axis label for focus-direction mean differences."""
    focus_label = DIRECTION_LABELS[focus_direction]
    comparator_label = DIRECTION_LABELS[comparator_direction(focus_direction)]
    return f"More selected trait signals in {focus_label} than {comparator_label}"


def focus_delta_sort_label(focus_direction: str) -> str:
    """Return a readable description of the default focus-delta ranking."""
    focus_label = DIRECTION_LABELS[focus_direction]
    comparator_label = DIRECTION_LABELS[comparator_direction(focus_direction)]
    return f"how much stronger selected traits are in {focus_label} than {comparator_label}"


def focus_signal_panel_title(focus_direction: str) -> str:
    """Return a plain-language title for focus-delta signal panels."""
    return f"Trait signal stronger in {DIRECTION_LABELS[focus_direction]}"


def focus_sort_label(focus_direction: str, sort_by: str) -> str:
    """Return a readable sort description for figure subtitles."""
    if sort_by in {"depletion-bias", "direction-delta"}:
        return focus_delta_sort_label(focus_direction)
    if sort_by == "coverage-weighted-signal":
        return "coverage-weighted signal: difference x sqrt(profiles) x sqrt(studies)"
    return sort_by.replace("-", " ")


def low_coverage_label() -> str:
    """Return the label used for rows with limited study support."""
    return f"limited study support ({LOW_COVERAGE_STUDY_THRESHOLD} studies)"


def style_limited_support_bars(
    bars: object,
    low_support: pd.Series,
    *,
    alpha: float,
) -> None:
    """Apply hatching and lower opacity to bars with limited study support."""
    for bar, is_limited in zip(bars, low_support, strict=True):
        if not bool(is_limited):
            continue
        bar.set_alpha(alpha)
        bar.set_hatch("//")
        bar.set_edgecolor("#475569")
        bar.set_linewidth(0.6)


def signal_legend_handles(focus_direction: str, has_low_support: bool) -> list[Patch]:
    """Return legend handles for the focus signal panel."""
    handles = [
        Patch(
            facecolor=SIGNAL_COLORS[focus_direction],
            edgecolor="none",
            alpha=0.9,
            label="trait-signal difference",
        )
    ]
    if has_low_support:
        handles.append(
            Patch(
                facecolor=SIGNAL_COLORS[focus_direction],
                edgecolor="#475569",
                hatch="//",
                alpha=0.55,
                label=low_coverage_label(),
            )
        )
    return handles


def resolve_focus_sort_column(focus_direction: str, sort_by: str) -> str:
    """Map a sort option to a summary table column."""
    focus_prefix = direction_prefix(focus_direction)
    aliases = {
        "depletion-bias": focus_delta_column(focus_direction),
        "direction-delta": focus_delta_column(focus_direction),
        "decreased-count": f"{focus_prefix}_positive_profile_count",
        "focus-count": f"{focus_prefix}_positive_profile_count",
        "decreased-mean": f"{focus_prefix}_mean_positive_selected_trait_count",
        "focus-mean": f"{focus_prefix}_mean_positive_selected_trait_count",
        "coverage": "positive_profile_count",
        "coverage-weighted-signal": focus_coverage_weighted_signal_column(focus_direction),
    }
    if sort_by not in aliases:
        msg = f"Unsupported focus sort option: {sort_by!r}"
        raise ValueError(msg)
    return aliases[sort_by]


def select_focus_conditions(
    summary: pd.DataFrame,
    *,
    focus_direction: str,
    sort_by: str,
    top_n: int | None,
    min_decreased_profiles: int,
    require_positive_delta: bool,
) -> pd.DataFrame:
    """Select conditions for the focus-direction signal and coverage diagrams."""
    if summary.empty:
        return summary

    selected = summary.copy()
    focus_prefix = direction_prefix(focus_direction)
    comparator_prefix = direction_prefix(comparator_direction(focus_direction))
    delta_column = focus_delta_column(focus_direction)
    focus_mean_column = f"{focus_prefix}_mean_positive_selected_trait_count"
    focus_profile_column = f"{focus_prefix}_profile_count"
    sort_column = resolve_focus_sort_column(focus_direction, sort_by)

    selected = selected[selected[focus_profile_column] >= min_decreased_profiles].copy()
    if require_positive_delta:
        selected = selected[(selected[delta_column] > 0) & (selected[focus_mean_column] > 0)].copy()

    if selected.empty:
        return selected

    selected = selected.sort_values(
        [
            sort_column,
            f"{focus_prefix}_positive_profile_count",
            f"{comparator_prefix}_positive_profile_count",
            "study_count",
            "condition",
        ],
        ascending=[False, False, False, False, True],
    )
    if top_n is None:
        return selected.copy()
    return selected.head(top_n).copy()


def set_signal_axis_limits(ax: plt.Axes, values: pd.Series) -> None:
    """Set padded x limits for a signal axis that may include negative values."""
    min_value = min(float(values.min()), 0.0)
    max_value = max(float(values.max()), 0.0)
    span = max(max_value - min_value, 1.0)
    ax.set_xlim(left=min_value - span * 0.08, right=max_value + span * 0.2)


def save_focus_signal_bar(
    focus_conditions: pd.DataFrame,
    *,
    output_path: Path,
    focus_direction: str,
    sort_by: str,
) -> None:
    """Save the primary focus-direction signal bar chart."""
    if focus_conditions.empty:
        return

    delta_column = focus_delta_column(focus_direction)
    plot_data = focus_conditions.iloc[::-1].copy()
    labels = [compact_label(condition, 42) for condition in plot_data["condition"]]
    values = plot_data[delta_column].astype(float)
    low_support = plot_data["study_count"] <= LOW_COVERAGE_STUDY_THRESHOLD

    sns.set_theme(style="whitegrid")
    fig_height = max(5.8, 0.48 * len(plot_data) + 1.8)
    fig, ax = plt.subplots(figsize=(11.5, fig_height))
    bars = ax.barh(labels, values, color=SIGNAL_COLORS[focus_direction], alpha=0.9)
    style_limited_support_bars(bars, low_support, alpha=0.55)
    ax.axvline(0, color="#64748b", linewidth=0.9)
    ax.set_title(
        focus_signal_panel_title(focus_direction),
        fontsize=15,
        fontweight="bold",
        pad=10,
    )
    ax.set_xlabel(focus_delta_axis_label(focus_direction))
    ax.set_ylabel("")
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8)
    ax.grid(axis="y", visible=False)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=9)
    ax.tick_params(axis="x", labelsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#e2e8f0")
    ax.legend(
        handles=signal_legend_handles(focus_direction, bool(low_support.any())),
        loc="lower right",
        bbox_to_anchor=(1.0, 1.03),
        frameon=False,
        ncol=2,
    )

    set_signal_axis_limits(ax, values)
    x_min, x_max = ax.get_xlim()
    x_pad = (x_max - x_min) * 0.015
    for y_index, value in enumerate(values):
        if value >= 0:
            x_position = value + x_pad
            ha = "left"
        else:
            x_position = value - x_pad
            ha = "right"
        ax.text(
            x_position,
            y_index,
            f"{value:.1f}",
            va="center",
            ha=ha,
            fontsize=8.4,
            color="#334155",
            fontweight="bold",
        )

    fig.text(
        0.01,
        0.01,
        f"Rows ranked by {focus_sort_label(focus_direction, sort_by)}.",
        fontsize=8.5,
        color="#475569",
    )
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_focus_coverage_bar(
    focus_conditions: pd.DataFrame,
    *,
    output_path: Path,
    focus_direction: str,
) -> None:
    """Save the focus-condition profile coverage companion chart."""
    if focus_conditions.empty:
        return

    focus_prefix = direction_prefix(focus_direction)
    comparator = comparator_direction(focus_direction)
    comparator_prefix = direction_prefix(comparator)
    focus_count_column = f"{focus_prefix}_positive_profile_count"
    comparator_count_column = f"{comparator_prefix}_positive_profile_count"
    plot_data = focus_conditions.iloc[::-1].copy()
    labels = [compact_label(condition, 42) for condition in plot_data["condition"]]
    focus_counts = plot_data[focus_count_column].astype(float)
    comparator_counts = plot_data[comparator_count_column].astype(float)
    totals = focus_counts + comparator_counts
    low_support = plot_data["study_count"] <= LOW_COVERAGE_STUDY_THRESHOLD

    sns.set_theme(style="whitegrid")
    fig_height = max(5.8, 0.48 * len(plot_data) + 1.8)
    fig, ax = plt.subplots(figsize=(11.5, fig_height))
    focus_bars = ax.barh(
        labels,
        focus_counts,
        color=DIRECTION_COLORS[focus_direction],
        alpha=0.86,
        label=DIRECTION_LABELS[focus_direction],
    )
    comparator_bars = ax.barh(
        labels,
        comparator_counts,
        left=focus_counts,
        color=DIRECTION_COLORS[comparator],
        alpha=0.86,
        label=DIRECTION_LABELS[comparator],
    )
    style_limited_support_bars(focus_bars, low_support, alpha=0.55)
    style_limited_support_bars(comparator_bars, low_support, alpha=0.55)
    ax.set_title(
        f"Profile Coverage for {focus_output_label(focus_direction).title()}-Signal Conditions",
        fontsize=15,
        fontweight="bold",
        pad=10,
    )
    ax.set_xlabel("Profiles with at least one positive selected trait signal")
    ax.set_ylabel("")
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8)
    ax.grid(axis="y", visible=False)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=9)
    ax.tick_params(axis="x", labelsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#e2e8f0")
    handles, labels_for_legend = ax.get_legend_handles_labels()
    if bool(low_support.any()):
        handles.append(Patch(facecolor="#f8fafc", edgecolor="#475569", hatch="//"))
        labels_for_legend.append(low_coverage_label())
    ax.legend(
        handles=handles,
        labels=labels_for_legend,
        loc="lower right",
        bbox_to_anchor=(1.0, 1.03),
        frameon=False,
        ncol=3,
    )

    max_total = max(float(totals.max()), 1.0)
    for y_index, row in enumerate(plot_data.itertuples(index=False)):
        total = float(getattr(row, focus_count_column)) + float(
            getattr(row, comparator_count_column)
        )
        ax.text(
            total + max_total * 0.015,
            y_index,
            f"{total:.0f}; {int(row.study_count)} studies",
            va="center",
            ha="left",
            fontsize=8.3,
            color="#334155",
        )
    ax.set_xlim(right=max_total * 1.28)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_focus_signal_diagram(
    focus_conditions: pd.DataFrame,
    *,
    output_path: Path,
    focus_direction: str,
    sort_by: str,
    min_decreased_profiles: int,
    min_studies: int,
    body_site: str | None,
) -> None:
    """Save the two-panel focus signal plus coverage diagram."""
    if focus_conditions.empty:
        return

    focus_prefix = direction_prefix(focus_direction)
    comparator = comparator_direction(focus_direction)
    comparator_prefix = direction_prefix(comparator)
    delta_column = focus_delta_column(focus_direction)
    focus_count_column = f"{focus_prefix}_positive_profile_count"
    comparator_count_column = f"{comparator_prefix}_positive_profile_count"

    plot_data = focus_conditions.iloc[::-1].copy()
    y_positions = list(range(len(plot_data)))
    labels = [compact_label(condition, 38) for condition in plot_data["condition"]]
    signal_values = plot_data[delta_column].astype(float)
    focus_counts = plot_data[focus_count_column].astype(float)
    comparator_counts = plot_data[comparator_count_column].astype(float)
    totals = focus_counts + comparator_counts
    low_support = plot_data["study_count"] <= LOW_COVERAGE_STUDY_THRESHOLD

    sns.set_theme(style="whitegrid")
    fig_height = max(6.2, 0.48 * len(plot_data) + 2.2)
    fig, (signal_ax, coverage_ax) = plt.subplots(
        ncols=2,
        sharey=True,
        figsize=(15.5, fig_height),
        gridspec_kw={"width_ratios": [1.05, 1.25]},
    )
    signal_bars = signal_ax.barh(
        y_positions,
        signal_values,
        color=SIGNAL_COLORS[focus_direction],
        alpha=0.9,
    )
    style_limited_support_bars(signal_bars, low_support, alpha=0.55)
    signal_ax.axvline(0, color="#64748b", linewidth=0.9)
    signal_ax.set_yticks(y_positions, labels)
    signal_ax.set_xlabel(focus_delta_axis_label(focus_direction))
    signal_ax.set_title(
        f"A. {focus_signal_panel_title(focus_direction)}",
        fontsize=11,
        fontweight="bold",
    )
    signal_ax.grid(axis="x", color="#e2e8f0", linewidth=0.8)
    signal_ax.grid(axis="y", visible=False)
    signal_ax.spines[["top", "right", "left"]].set_visible(False)
    signal_ax.spines["bottom"].set_color("#e2e8f0")
    signal_ax.tick_params(axis="y", labelsize=8.8)
    signal_ax.tick_params(axis="x", labelsize=8.8)
    signal_ax.legend(
        handles=signal_legend_handles(focus_direction, bool(low_support.any())),
        loc="lower right",
        bbox_to_anchor=(1.0, 1.03),
        frameon=False,
        ncol=1,
    )
    set_signal_axis_limits(signal_ax, signal_values)

    x_min, x_max = signal_ax.get_xlim()
    x_pad = (x_max - x_min) * 0.015
    for y_index, value in enumerate(signal_values):
        if value >= 0:
            x_position = value + x_pad
            ha = "left"
        else:
            x_position = value - x_pad
            ha = "right"
        signal_ax.text(
            x_position,
            y_index,
            f"{value:.1f}",
            va="center",
            ha=ha,
            fontsize=8.0,
            color="#334155",
            fontweight="bold",
        )

    focus_bars = coverage_ax.barh(
        y_positions,
        focus_counts,
        color=DIRECTION_COLORS[focus_direction],
        alpha=0.86,
        label=DIRECTION_LABELS[focus_direction],
    )
    comparator_bars = coverage_ax.barh(
        y_positions,
        comparator_counts,
        left=focus_counts,
        color=DIRECTION_COLORS[comparator],
        alpha=0.86,
        label=DIRECTION_LABELS[comparator],
    )
    style_limited_support_bars(focus_bars, low_support, alpha=0.55)
    style_limited_support_bars(comparator_bars, low_support, alpha=0.55)
    coverage_ax.set_xlabel("Positive-signal profiles")
    coverage_ax.set_title("B. Profile coverage", fontsize=11, fontweight="bold")
    coverage_ax.grid(axis="x", color="#e2e8f0", linewidth=0.8)
    coverage_ax.grid(axis="y", visible=False)
    coverage_ax.spines[["top", "right", "left"]].set_visible(False)
    coverage_ax.spines["bottom"].set_color("#e2e8f0")
    coverage_ax.tick_params(axis="x", labelsize=8.8)
    coverage_ax.tick_params(axis="y", labelleft=False, left=False)
    coverage_handles, coverage_labels = coverage_ax.get_legend_handles_labels()
    if bool(low_support.any()):
        coverage_handles.append(Patch(facecolor="#f8fafc", edgecolor="#475569", hatch="//"))
        coverage_labels.append(low_coverage_label())
    coverage_ax.legend(
        handles=coverage_handles,
        labels=coverage_labels,
        loc="lower right",
        bbox_to_anchor=(1.0, 1.03),
        frameon=False,
        ncol=3,
    )

    max_total = max(float(totals.max()), 1.0)
    for y_index, row in enumerate(plot_data.itertuples(index=False)):
        total = float(getattr(row, focus_count_column)) + float(
            getattr(row, comparator_count_column)
        )
        coverage_ax.text(
            total + max_total * 0.018,
            y_index,
            f"{total:.0f}; {int(row.study_count)} studies",
            va="center",
            ha="left",
            fontsize=7.8,
            color="#334155",
        )
    coverage_ax.set_xlim(right=max_total * 1.35)

    site_text = f"; body site: {body_site}" if body_site else ""
    fig.suptitle(
        f"{focus_output_label(focus_direction).title()}-Profile Conditions: "
        "Trait Signal and Coverage",
        fontsize=16,
        fontweight="bold",
        y=0.985,
    )
    fig.text(
        0.5,
        0.94,
        f"Ranked by {focus_sort_label(focus_direction, sort_by)}; filters: >= "
        f"{min_decreased_profiles} {focus_prefix} profiles, >= {min_studies} studies"
        f"{site_text}.",
        ha="center",
        fontsize=9.2,
        color="#475569",
    )
    fig.text(
        0.01,
        0.026,
        "Panel A: signal bars show the difference between the "
        f"{DIRECTION_LABELS[focus_direction]} and {DIRECTION_LABELS[comparator]} averages.",
        fontsize=8.2,
        color="#475569",
    )
    fig.text(
        0.01,
        0.01,
        "Coverage counts profiles with at least one positive selected trait z-score.",
        fontsize=8.2,
        color="#475569",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.925])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_signal_coverage_table(
    focus_conditions: pd.DataFrame,
    *,
    focus_direction: str,
) -> pd.DataFrame:
    """Build the uncapped table used by the signal-vs-coverage scatterplot."""
    if focus_conditions.empty:
        return focus_conditions

    focus_prefix = direction_prefix(focus_direction)
    comparator_prefix = direction_prefix(comparator_direction(focus_direction))
    delta_column = focus_delta_column(focus_direction)
    focus_positive_column = f"{focus_prefix}_positive_profile_count"
    comparator_positive_column = f"{comparator_prefix}_positive_profile_count"
    focus_mean_column = f"{focus_prefix}_mean_positive_selected_trait_count"
    comparator_mean_column = f"{comparator_prefix}_mean_positive_selected_trait_count"
    focus_profile_column = f"{focus_prefix}_profile_count"
    comparator_profile_column = f"{comparator_prefix}_profile_count"
    weighted_signal_column = focus_coverage_weighted_signal_column(focus_direction)

    table = focus_conditions.copy()
    table["signal_difference"] = table[delta_column].astype(float)
    table["focus_positive_signal_profiles"] = table[focus_positive_column].astype(int)
    table["comparator_positive_signal_profiles"] = table[comparator_positive_column].astype(int)
    table["positive_signal_profiles"] = (
        table["focus_positive_signal_profiles"] + table["comparator_positive_signal_profiles"]
    )
    table["focus_mean_positive_selected_trait_count"] = table[focus_mean_column].astype(float)
    table["comparator_mean_positive_selected_trait_count"] = table[comparator_mean_column].astype(
        float
    )
    table["focus_profile_count"] = table[focus_profile_column].astype(int)
    table["comparator_profile_count"] = table[comparator_profile_column].astype(int)
    table["limited_study_support"] = table["study_count"] <= LOW_COVERAGE_STUDY_THRESHOLD
    table["coverage_weighted_signal"] = table[weighted_signal_column].astype(float)

    ordered_columns = [
        "condition",
        "signal_difference",
        "positive_signal_profiles",
        "focus_positive_signal_profiles",
        "comparator_positive_signal_profiles",
        "focus_mean_positive_selected_trait_count",
        "comparator_mean_positive_selected_trait_count",
        "focus_profile_count",
        "comparator_profile_count",
        "profile_count",
        "study_count",
        "positive_profile_fraction",
        "limited_study_support",
        "coverage_weighted_signal",
    ]
    return table[ordered_columns].sort_values(
        ["signal_difference", "positive_signal_profiles", "study_count", "condition"],
        ascending=[False, False, False, True],
    )


def save_signal_vs_coverage_scatter(
    scatter_table: pd.DataFrame,
    *,
    output_path: Path,
    focus_direction: str,
    min_decreased_profiles: int,
    min_studies: int,
    body_site: str | None,
) -> None:
    """Save an uncapped scatterplot comparing signal strength and coverage."""
    if scatter_table.empty:
        return

    plot_data = scatter_table.copy()
    sizes = 40 + plot_data["positive_signal_profiles"].pow(0.5) * 12
    low_support = plot_data["limited_study_support"].astype(bool)

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    points = ax.scatter(
        plot_data["signal_difference"],
        plot_data["positive_signal_profiles"],
        c=plot_data["study_count"],
        s=sizes,
        cmap="viridis",
        alpha=0.82,
        edgecolors="#334155",
        linewidths=0.35,
    )
    if bool(low_support.any()):
        ax.scatter(
            plot_data.loc[low_support, "signal_difference"],
            plot_data.loc[low_support, "positive_signal_profiles"],
            s=sizes[low_support] + 36,
            facecolors="none",
            edgecolors="#111827",
            marker="s",
            linewidths=1.0,
            label=low_coverage_label(),
        )
        ax.legend(loc="upper right", frameon=False)

    ax.axvline(0, color="#64748b", linewidth=0.9)
    ax.set_title(
        f"{focus_output_label(focus_direction).title()}-Profile Signal vs Coverage",
        fontsize=15,
        fontweight="bold",
        pad=10,
    )
    ax.set_xlabel(focus_delta_axis_label(focus_direction))
    ax.set_ylabel("Positive-signal profiles")
    ax.grid(axis="both", color="#e2e8f0", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_color("#e2e8f0")
    ax.spines["bottom"].set_color("#e2e8f0")
    ax.set_ylim(bottom=0)
    colorbar = fig.colorbar(points, ax=ax, pad=0.02)
    colorbar.set_label("Studies")

    site_text = f"; body site: {body_site}" if body_site else ""
    fig.text(
        0.5,
        0.01,
        f"All {len(plot_data)} focus-eligible conditions after filters: >= "
        f"{min_decreased_profiles} focus profiles, >= {min_studies} studies"
        f"{site_text}. Point size follows positive-signal profile count.",
        ha="center",
        fontsize=8.4,
        color="#475569",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_condition_trait_signal_heatmap(
    matrix: pd.DataFrame,
    scored: pd.DataFrame,
    focus_conditions: pd.DataFrame,
    *,
    focus_direction: str,
) -> pd.DataFrame:
    """Build an uncapped condition x trait signal-difference matrix."""
    if focus_conditions.empty:
        return pd.DataFrame()

    comparator = comparator_direction(focus_direction)
    positive_matrix = matrix.gt(0).astype(float)
    scored_subset = scored[scored["condition"].isin(focus_conditions["condition"])].copy()
    zero_trait_signal = pd.Series(0.0, index=positive_matrix.columns)

    records: list[pd.Series] = []
    for row in focus_conditions.itertuples(index=False):
        condition_profiles = scored_subset[scored_subset["condition"] == row.condition]
        focus_ids = condition_profiles[
            condition_profiles["inc_dec"] == focus_direction
        ].index.intersection(positive_matrix.index)
        comparator_ids = condition_profiles[
            condition_profiles["inc_dec"] == comparator
        ].index.intersection(positive_matrix.index)

        focus_mean = (
            positive_matrix.loc[focus_ids].mean(axis=0) if len(focus_ids) else zero_trait_signal
        )
        comparator_mean = (
            positive_matrix.loc[comparator_ids].mean(axis=0)
            if len(comparator_ids)
            else zero_trait_signal
        )
        records.append((focus_mean - comparator_mean).rename(row.condition))

    return pd.DataFrame(records)


def save_condition_trait_signal_heatmap(
    heatmap_matrix: pd.DataFrame,
    *,
    output_path: Path,
    focus_direction: str,
) -> None:
    """Save an uncapped condition x selected-trait signal heatmap."""
    if heatmap_matrix.empty:
        return

    comparator = comparator_direction(focus_direction)
    row_labels = [compact_label(condition, 48) for condition in heatmap_matrix.index]
    column_labels = [compact_label(trait, 34) for trait in heatmap_matrix.columns]
    row_count, trait_count = heatmap_matrix.shape

    sns.set_theme(style="white")
    fig_width = max(14.0, 0.34 * trait_count + 5.0)
    fig_height = max(8.5, 0.23 * row_count + 3.2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        heatmap_matrix,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.12,
        linecolor="#f1f5f9",
        xticklabels=column_labels,
        yticklabels=row_labels,
        cbar_kws={
            "label": (
                "Positive-signal fraction difference "
                f"({DIRECTION_LABELS[focus_direction]} - {DIRECTION_LABELS[comparator]})"
            )
        },
    )
    ax.set_title(
        f"{focus_output_label(focus_direction).title()}-Profile Condition x Trait "
        "Signal Differences",
        fontsize=15,
        fontweight="bold",
        pad=30,
    )
    ax.set_xlabel("Selected traits")
    ax.set_ylabel("Conditions")
    ax.tick_params(axis="x", labelrotation=45, labelsize=6.4)
    ax.tick_params(axis="y", labelsize=5.8 if row_count > 80 else 6.6)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")

    fig.text(
        0.5,
        0.955,
        f"Rows include all {row_count} focus-eligible conditions; columns include all "
        f"{trait_count} selected traits.",
        ha="center",
        fontsize=8.8,
        color="#475569",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_output_path = output_path.with_name(f"{output_path.stem}.tmp{output_path.suffix}")
    if tmp_output_path.exists():
        tmp_output_path.unlink()
    fig.savefig(tmp_output_path, dpi=300, bbox_inches="tight")
    try:
        tmp_output_path.replace(output_path)
    except PermissionError:
        if not output_path.exists():
            raise
        tmp_output_path.unlink(missing_ok=True)
        print(f"Skipped locked heatmap PNG: {output_path}")
    plt.close(fig)


def write_readme(
    output_path: Path,
    summary: pd.DataFrame,
    focus_conditions: pd.DataFrame,
    focus_conditions_all: pd.DataFrame,
    *,
    matrix_path: Path,
    metadata_path: Path,
    min_profiles: int,
    min_studies: int,
    top_n: int,
    include_broad_labels: bool,
    focus_direction: str,
    sort_by: str,
    min_decreased_profiles: int,
    body_site: str | None,
    trait_count: int,
    standalone_bars: bool,
) -> None:
    """Write a short methods/readout note for the generated bar plots."""
    top_count = summary.sort_values("positive_profile_count", ascending=False).head(5)
    top_mean = summary.sort_values(
        "total_direction_mean_positive_selected_trait_count",
        ascending=False,
    ).head(5)
    focus_prefix = direction_prefix(focus_direction)
    comparator_prefix = direction_prefix(comparator_direction(focus_direction))
    focus_delta = focus_delta_column(focus_direction)
    focus_label = focus_output_label(focus_direction)
    focus_condition_count = len(focus_conditions_all)
    lines = [
        "# Profile-Condition Enrichment Bars",
        "",
        "These plots are profile-data companions to Ilyas's ORA enrichment infographic "
        "from `data/derived/enrichment_infographic.png` on the "
        "`enrichment-analysis-heatmap-cluster` branch.",
        "",
        "They do not rerun the hypergeometric ORA. Instead, they use the selected "
        "profile heatmap z-score matrix. A selected trait is counted when its z-score "
        "is positive, meaning that trait's fraction is above the dataset-average "
        "fraction for the selected trait column.",
        "",
        "## Inputs",
        "",
        f"- Matrix: `{matrix_path.relative_to(REPO_ROOT)}`",
        f"- Metadata: `{metadata_path.relative_to(REPO_ROOT)}`",
        f"- Minimum condition profiles: {min_profiles}",
        f"- Minimum condition studies: {min_studies}",
        f"- Broad metadata labels included: {include_broad_labels}",
        f"- Body-site filter: {body_site if body_site else 'none'}",
        f"- Focus direction: {focus_prefix} profiles",
        f"- Focus-condition sort: {sort_by}",
        f"- Minimum focus-direction profiles: {min_decreased_profiles}",
        "",
        "## Outputs",
        "",
    ]
    if standalone_bars:
        lines.extend(
            [
                "- `condition_positive_profile_count_bar.png`: condition x number of profiles "
                "with at least one positive selected trait signal, stacked by "
                "increased/decreased BugSigDB profile direction.",
                "- `condition_average_positive_traits_bar.png`: condition x average number of "
                "positive selected traits, stacked as decreased-profile mean plus "
                "increased-profile mean in the same style as the ORA infographic condition "
                "panel.",
                f"- `condition_{focus_label}_trait_signal_bar.png`: focus-condition bars "
                f"ranked by the selected sort option; signal is `{focus_delta}`.",
                f"- `condition_{focus_label}_profile_coverage_bar.png`: profile coverage for "
                "the same focus conditions, stacked by profile direction.",
            ]
        )
    lines += [
        f"- `condition_{focus_label}_profile_signal_diagram.png`: two-panel "
        "diagram combining direction-specific trait signal and coverage.",
        f"- `condition_{focus_label}_signal_vs_coverage_scatter.png`: uncapped "
        "scatterplot of direction-specific signal difference versus positive-signal "
        "profile coverage.",
        f"- `condition_{focus_label}_signal_vs_coverage_scatter.tsv`: full data table "
        "used by the scatterplot.",
        f"- `condition_{focus_label}_trait_signal_heatmap.png`: uncapped condition x "
        "selected-trait heatmap of each trait's positive-signal fraction difference "
        "between directions.",
        f"- `condition_{focus_label}_trait_signal_heatmap_matrix.tsv`: full condition x "
        "trait matrix used by the heatmap.",
        "- `condition_profile_trait_summary.tsv`: profile-data summary table used by the plots.",
        "",
        "## Top Conditions By Profile Count",
        "",
        "| condition | positive profiles | increased | decreased | studies |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in top_count.itertuples(index=False):
        lines.append(
            f"| {row.condition} | {int(row.positive_profile_count)} | "
            f"{int(row.increased_positive_profile_count)} | "
            f"{int(row.decreased_positive_profile_count)} | {int(row.study_count)} |"
        )

    lines.extend(
        [
            "",
            "## Top Conditions By Average Positive Trait Depth",
            "",
            "| condition | decreased mean | increased mean | direction-mean total | "
            "profiles | studies |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in top_mean.itertuples(index=False):
        lines.append(
            f"| {row.condition} | "
            f"{row.decreased_mean_positive_selected_trait_count:.2f} | "
            f"{row.increased_mean_positive_selected_trait_count:.2f} | "
            f"{row.total_direction_mean_positive_selected_trait_count:.2f} | "
            f"{int(row.profile_count)} | {int(row.study_count)} |"
        )

    if not focus_conditions.empty:
        lines.extend(
            [
                "",
                f"## Top {focus_label.title()}-Profile Conditions",
                "",
                "| condition | direction delta | focus mean | comparator mean | "
                "focus profiles | comparator profiles | studies |",
                "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            ]
        )
        focus_mean_column = f"{focus_prefix}_mean_positive_selected_trait_count"
        comparator_mean_column = f"{comparator_prefix}_mean_positive_selected_trait_count"
        focus_count_column = f"{focus_prefix}_positive_profile_count"
        comparator_count_column = f"{comparator_prefix}_positive_profile_count"
        for row in focus_conditions.head(5).itertuples(index=False):
            lines.append(
                f"| {row.condition} | {getattr(row, focus_delta):.2f} | "
                f"{getattr(row, focus_mean_column):.2f} | "
                f"{getattr(row, comparator_mean_column):.2f} | "
                f"{int(getattr(row, focus_count_column))} | "
                f"{int(getattr(row, comparator_count_column))} | "
                f"{int(row.study_count)} |"
            )

    lines.extend(
        [
            "",
            "## Uncapped Follow-Up Figures",
            "",
            f"The scatterplot and condition x trait heatmap include all {focus_condition_count} "
            "focus-eligible condition rows after the shared filters, minimum "
            "focus-direction profile filter, and positive direction-delta filter. "
            f"The heatmap includes all {trait_count} selected trait columns from the input "
            "matrix.",
            "",
            "The scatterplot is intended to separate large but sparsely supported "
            "direction deltas from conditions with broader profile/study support. "
            "`coverage_weighted_signal` in the scatter TSV is a descriptive product of "
            "the direction delta, square-root positive-signal profile count, and "
            "square-root study count; it is not used as a filter.",
            "",
            "The heatmap decomposes the focus-direction signal by trait. Each cell is "
            "the focus-direction fraction of profiles with a positive selected-trait "
            "z-score minus the comparator-direction fraction for the same condition "
            "and trait.",
        ]
    )

    coverage_views = (
        "The count plot shows where profile coverage is deepest; the average-depth plot shows "
        "where profiles have the broadest above-average selected trait signal."
        if standalone_bars
        else "The two condition tables above show where profile coverage is deepest and where "
        "profiles have the broadest above-average selected trait signal."
    )
    capped_views = (
        "The diagram and bar plots display" if standalone_bars else "The diagram displays"
    )
    lines.extend(
        [
            "",
            "## Interpretation Note",
            "",
            f"Use these as a profile-level analogue to the ORA condition panel. {coverage_views} "
            "Because the metric is z-score based, it should be described as "
            "`positive selected trait signal` unless it is later replaced with formal "
            "ORA q-value enrichment counts.",
            "",
            f"The {focus_label}-condition diagram ranks conditions by "
            f"{focus_sort_label(focus_direction, sort_by)} while keeping the "
            "comparator visible as a guardrail. A positive direction delta means "
            "the focus direction has a higher mean positive selected-trait count "
            "than the comparator direction. Hatched rows have limited study support.",
            "",
            f"{capped_views} the top {top_n} conditions after filtering. "
            "The scatterplot and heatmap are not top-N capped.",
            "",
        ]
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Create condition bar plots from profile-level trait signals.",
    )
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX_PATH)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--top-n", type=int, default=15)
    parser.add_argument("--min-profiles", type=int, default=10)
    parser.add_argument("--min-studies", type=int, default=2)
    parser.add_argument(
        "--focus-direction",
        choices=["dec", "inc"],
        default="dec",
        help="Profile direction for the focused signal and coverage diagrams.",
    )
    parser.add_argument(
        "--sort-by",
        choices=[
            "depletion-bias",
            "direction-delta",
            "decreased-count",
            "focus-count",
            "decreased-mean",
            "focus-mean",
            "coverage",
            "coverage-weighted-signal",
        ],
        default="depletion-bias",
        help="Ranking used for the focused condition diagrams.",
    )
    parser.add_argument(
        "--min-decreased-profiles",
        type=int,
        default=3,
        help="Minimum number of profiles in the focus direction for focused diagrams.",
    )
    parser.add_argument(
        "--body-site",
        help="Optional body_site_group or body_site value to filter before condition summaries.",
    )
    parser.add_argument(
        "--include-negative-delta",
        action="store_true",
        help="Keep focus conditions even when the focus direction is not above the comparator.",
    )
    parser.add_argument(
        "--include-broad-labels",
        action="store_true",
        help="Keep broad study-design/metadata labels in condition rankings.",
    )
    parser.add_argument(
        "--standalone-bars",
        action="store_true",
        help=(
            "Also draw the four standalone bar charts. The committed folder keeps only the "
            "combined diagram, scatterplot and heatmap for each direction."
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run profile-condition summary and plotting."""
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matrix = load_profile_matrix(args.matrix)
    scored = load_profile_scores(args.matrix, args.metadata)
    scored = filter_scored_profiles(scored, args.body_site)
    summary = summarize_conditions(
        scored,
        min_profiles=args.min_profiles,
        min_studies=args.min_studies,
        include_broad_labels=args.include_broad_labels,
    )
    if summary.empty:
        msg = "No condition groups passed the requested profile/study filters."
        raise ValueError(msg)

    focus_conditions_all = select_focus_conditions(
        summary,
        focus_direction=args.focus_direction,
        sort_by=args.sort_by,
        top_n=None,
        min_decreased_profiles=args.min_decreased_profiles,
        require_positive_delta=not args.include_negative_delta,
    )
    if focus_conditions_all.empty:
        msg = "No focus-direction condition groups passed the requested filters."
        raise ValueError(msg)
    focus_conditions = focus_conditions_all.head(args.top_n).copy()

    summary.to_csv(args.output_dir / "condition_profile_trait_summary.tsv", sep="\t", index=False)
    focus_label = focus_output_label(args.focus_direction)
    scatter_table = build_signal_coverage_table(
        focus_conditions_all,
        focus_direction=args.focus_direction,
    )
    scatter_table.to_csv(
        args.output_dir / f"condition_{focus_label}_signal_vs_coverage_scatter.tsv",
        sep="\t",
        index=False,
    )
    heatmap_matrix = build_condition_trait_signal_heatmap(
        matrix,
        scored,
        focus_conditions_all,
        focus_direction=args.focus_direction,
    )
    heatmap_matrix_for_output = heatmap_matrix.rename(
        columns=lambda value: clean_text(value).replace("\n", " ")
    )
    heatmap_matrix_for_output.to_csv(
        args.output_dir / f"condition_{focus_label}_trait_signal_heatmap_matrix.tsv",
        sep="\t",
    )
    if args.standalone_bars:
        save_stacked_condition_bar(
            summary,
            output_path=args.output_dir / "condition_positive_profile_count_bar.png",
            value_columns=(
                "decreased_positive_profile_count",
                "increased_positive_profile_count",
            ),
            total_column="positive_profile_count",
            title="Profile Data: Conditions by Positive Trait-Signal Coverage",
            xlabel="Profiles with at least one positive selected trait signal",
            top_n=args.top_n,
            label_format="{:.0f}",
        )
        save_stacked_condition_bar(
            summary,
            output_path=args.output_dir / "condition_average_positive_traits_bar.png",
            value_columns=(
                "decreased_mean_positive_selected_trait_count",
                "increased_mean_positive_selected_trait_count",
            ),
            total_column="total_direction_mean_positive_selected_trait_count",
            title="Profile Data: Avg Positive Trait Signals per Condition",
            xlabel=(
                "Mean positive selected traits per profile (decreased + increased direction means)"
            ),
            top_n=args.top_n,
            label_format="{:.1f}",
        )
        save_focus_signal_bar(
            focus_conditions,
            output_path=args.output_dir / f"condition_{focus_label}_trait_signal_bar.png",
            focus_direction=args.focus_direction,
            sort_by=args.sort_by,
        )
        save_focus_coverage_bar(
            focus_conditions,
            output_path=args.output_dir / f"condition_{focus_label}_profile_coverage_bar.png",
            focus_direction=args.focus_direction,
        )
    save_focus_signal_diagram(
        focus_conditions,
        output_path=args.output_dir / f"condition_{focus_label}_profile_signal_diagram.png",
        focus_direction=args.focus_direction,
        sort_by=args.sort_by,
        min_decreased_profiles=args.min_decreased_profiles,
        min_studies=args.min_studies,
        body_site=args.body_site,
    )
    save_signal_vs_coverage_scatter(
        scatter_table,
        output_path=args.output_dir / f"condition_{focus_label}_signal_vs_coverage_scatter.png",
        focus_direction=args.focus_direction,
        min_decreased_profiles=args.min_decreased_profiles,
        min_studies=args.min_studies,
        body_site=args.body_site,
    )
    save_condition_trait_signal_heatmap(
        heatmap_matrix,
        output_path=args.output_dir / f"condition_{focus_label}_trait_signal_heatmap.png",
        focus_direction=args.focus_direction,
    )
    write_readme(
        args.output_dir / "README.md",
        summary,
        focus_conditions,
        focus_conditions_all,
        matrix_path=args.matrix,
        metadata_path=args.metadata,
        min_profiles=args.min_profiles,
        min_studies=args.min_studies,
        top_n=args.top_n,
        include_broad_labels=args.include_broad_labels,
        focus_direction=args.focus_direction,
        sort_by=args.sort_by,
        min_decreased_profiles=args.min_decreased_profiles,
        body_site=args.body_site,
        trait_count=matrix.shape[1],
        standalone_bars=args.standalone_bars,
    )

    print(f"Wrote profile-condition bar plots to: {args.output_dir}")
    print(f"Profiles analyzed: {len(scored):,}")
    print(f"Condition rows: {len(summary):,}")
    print(f"Focus-condition rows: {len(focus_conditions_all):,}")
    print(
        "Top average-depth condition: "
        f"{summary.iloc[0]['condition']} "
        f"({summary.iloc[0]['total_direction_mean_positive_selected_trait_count']:.1f})"
    )
    print(
        f"Top {focus_label}-profile condition: "
        f"{focus_conditions.iloc[0]['condition']} "
        f"({focus_conditions.iloc[0][focus_delta_column(args.focus_direction)]:.1f})"
    )


if __name__ == "__main__":
    main()
