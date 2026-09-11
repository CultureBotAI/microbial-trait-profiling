"""Create combined increased/decreased Sankey maps for six shared traits."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, replace
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

from create_condition_body_site_guild_maps import (
    FACULTATIVE_CONFIG,
    MIN_STUDIES,
    SLATE_500,
    SLATE_700,
    SLATE_900,
    V3_BACKGROUND_EDGE_Y_OFFSET,
    V3_GROUPED_MODE,
    V3_MIN_FOCUS_ENRICHED_PROFILES,
    add_v3_display_rollups,
    clean_text,
    draw_curve,
    draw_node_v3,
    join_unique,
    nan_mean,
    node_positions,
    scaled_line_width_from_range,
    select_v3_display_edges,
    summarize_v3_condition_body_site_edges_for_mode,
    weighted_v3_node_stats,
    wrap_label,
)
from create_six_shared_trait_cytoscape_heatmaps import (
    write_six_shared_trait_cytoscape_heatmaps,
)
from create_six_shared_trait_graphs import HTML_FILENAME, write_six_shared_trait_graphs
from create_six_shared_trait_static_graphs import write_six_shared_trait_static_graphs

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"
OUTPUT_DIR = INVESTIGATION_DIR / "six-shared-traits"

PROFILE_MATRIX_PATH = REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv"
PROFILE_METADATA_PATH = VIS_DIR / "study_trait_profile_study_metadata.tsv"
TRAIT_STATS_PATH = INVESTIGATION_DIR / "trait_profile_paired_contrast_trait_stats.tsv"

PROFILE_METADATA_COLUMNS = ["study_id", "experiment_id", "inc_dec"]
DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}
COMPARATOR_DIRECTIONS = {"inc": "dec", "dec": "inc"}

EDGE_BACKGROUND = "#cbd5e1"
NODE_EDGE = "#475569"
NODE_FILL = "white"
LIGHTNESS_FRACTION = 0.18
SHARED_PATTERN_MIN_TRAITS = 3
SHARED_PATTERN_HIGH_SUPPORT_MIN_PROFILES = 30
SHARED_PATTERN_MAX_PAIRS = 12
PARALLEL_EDGE_SPACING = 0.004
MAX_PARALLEL_EDGE_OFFSET = 0.010
BODY_SITE_BUNDLE_X = 0.585
BODY_SITE_BUNDLE_EDGE_SPACING = 0.007
MAX_BODY_SITE_BUNDLE_EDGE_OFFSET = 0.018
BIDIRECTIONAL_ANALYSIS = "six_shared_traits_combined_shared_high_support"
BIDIRECTIONAL_STEM = f"{BIDIRECTIONAL_ANALYSIS}_profile_sankey"
BIDIRECTIONAL_README_HEADING = "## Increased and Decreased High-Support Overlay"


def wrapped_plot_text(value: str, width: int = 150) -> str:
    """Wrap long figure annotations so they do not widen the saved canvas."""
    return "\n".join(wrap(value, width=width))


@dataclass(frozen=True)
class SharedTrait:
    """One q-significant trait included in the combined Sankey pair."""

    trait_id: str
    label: str
    slug: str
    base_color: str


SHARED_TRAITS = (
    SharedTrait(
        trait_id="METPO:1000608",
        label="facultatively aerobic",
        slug="facultatively_aerobic",
        base_color="#0072b2",
    ),
    SharedTrait(
        trait_id="METPO:2000103|GO:0030245",
        label="cellulose catabolic process",
        slug="cellulose_catabolic_process",
        base_color="#d55e00",
    ),
    SharedTrait(
        trait_id="METPO:1000668",
        label="coccus shaped",
        slug="coccus_shaped",
        base_color="#009e73",
    ),
    SharedTrait(
        trait_id="METPO:1000683",
        label="sphere shaped",
        slug="sphere_shaped",
        base_color="#cc79a7",
    ),
    SharedTrait(
        trait_id="METPO:2000103|pathways:arsenate_reduction",
        label="arsenate reduction",
        slug="arsenate_reduction",
        base_color="#e69f00",
    ),
    SharedTrait(
        trait_id="METPO:1000623",
        label="moderately halophilic",
        slug="moderately_halophilic",
        base_color="#6a3d9a",
    ),
)

TRAIT_RANK_BY_ID = {trait.trait_id: index for index, trait in enumerate(SHARED_TRAITS)}
TRAIT_SORT_COLUMNS = ["focus_enriched_profile_count", "rank_score", "trait_rank"]


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    """Read a TSV file with local default options."""
    return pd.read_csv(path, sep="\t", **kwargs)


def reset_output_dir() -> None:
    """Clear the dedicated generated-output folder before regeneration."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    resolved_output = OUTPUT_DIR.resolve()
    resolved_parent = INVESTIGATION_DIR.resolve()
    if resolved_output.name != "six-shared-traits" or resolved_parent not in (
        resolved_output,
        *resolved_output.parents,
    ):
        msg = f"Refusing to clean unexpected output path: {resolved_output}"
        raise RuntimeError(msg)

    for path in OUTPUT_DIR.iterdir():
        if path.is_file() or path.is_symlink():
            path.unlink()
        elif path.is_dir():
            msg = f"Refusing to remove nested directory in output path: {path}"
            raise RuntimeError(msg)


def load_profile_metadata() -> pd.DataFrame:
    """Load the human-host profile metadata used by the current Sankey maps."""
    metadata = read_tsv(PROFILE_METADATA_PATH, index_col=0)
    metadata = metadata.copy()
    metadata["profile_id"] = metadata.index
    for column in ["study_id", "experiment_id", "inc_dec", "condition", "body_site_group"]:
        metadata[column] = metadata[column].map(clean_text)
    if "body_site" in metadata:
        metadata["body_site"] = metadata["body_site"].map(clean_text)
    else:
        metadata["body_site"] = "unknown"
    metadata["total_trait_signal"] = pd.to_numeric(
        metadata["total_trait_signal"],
        errors="coerce",
    ).fillna(0.0)
    return metadata[metadata["inc_dec"].isin(DIRECTION_LABELS)].copy()


def load_raw_trait_counts(traits: tuple[SharedTrait, ...], metadata: pd.DataFrame) -> pd.DataFrame:
    """Load raw count columns for the requested traits and align to profile metadata."""
    trait_ids = [trait.trait_id for trait in traits]
    profile_counts = read_tsv(PROFILE_MATRIX_PATH, usecols=PROFILE_METADATA_COLUMNS + trait_ids)
    profile_ids = [
        f"{study_id}|{experiment_id}|{inc_dec}"
        for study_id, experiment_id, inc_dec in profile_counts[PROFILE_METADATA_COLUMNS].itertuples(
            index=False, name=None
        )
    ]
    values = profile_counts[trait_ids].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    values.index = profile_ids
    values = values.groupby(level=0, sort=False).sum()

    missing_profiles = metadata.index.difference(values.index)
    if not missing_profiles.empty:
        examples = ", ".join(missing_profiles[:5].tolist())
        msg = (
            f"{len(missing_profiles):,} metadata profiles were not found in "
            f"{PROFILE_MATRIX_PATH}: {examples}"
        )
        raise ValueError(msg)
    return values.loc[metadata.index, trait_ids].copy()


def build_zscore_matrix(raw_trait_counts: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Create z-scores from log-scaled trait fractions using all-trait row totals."""
    row_sums = metadata["total_trait_signal"].replace(0, pd.NA)
    composition = raw_trait_counts.div(row_sums, axis=0).fillna(0.0)
    scaled = (composition * 1000).apply(lambda column: column.map(math.log1p))
    zscores = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0, ddof=0).replace(0, pd.NA)
    return zscores.fillna(0.0)


def load_q_significance() -> pd.DataFrame:
    """Load pooled paired-contrast q statistics for the six traits."""
    stats = read_tsv(TRAIT_STATS_PATH)
    stats = stats[stats["trait"].isin([trait.trait_id for trait in SHARED_TRAITS])].copy()
    numeric_columns = [
        "q_value",
        "p_value",
        "mean_increased_minus_decreased",
        "positive_pair_fraction",
        "negative_pair_fraction",
        "pair_count",
        "profile_count",
        "n_profiles_where_annotated",
    ]
    for column in numeric_columns:
        stats[column] = pd.to_numeric(stats[column], errors="coerce")
    return stats.set_index("trait", drop=False)


def score_profiles(
    metadata: pd.DataFrame,
    zscores: pd.DataFrame,
    trait: SharedTrait,
) -> pd.DataFrame:
    """Attach one trait's z-score fields using the guild-map column names."""
    scored = metadata.copy()
    scored["trait_id"] = trait.trait_id
    scored["trait_label"] = trait.label
    scored["trait_z"] = zscores[trait.trait_id]
    scored["guild_z_mean"] = scored["trait_z"]
    scored["guild_positive_trait_count"] = (scored["trait_z"] > 0.0).astype(int)
    scored["has_guild_enriched_trait"] = scored["guild_positive_trait_count"] > 0
    return scored


def single_trait_map_config(trait: SharedTrait, focus_direction: str):
    """Return a MapConfig compatible with the existing grouped Sankey utilities."""
    focus_label = DIRECTION_LABELS[focus_direction]
    comparator_direction = COMPARATOR_DIRECTIONS[focus_direction]
    comparator_label = DIRECTION_LABELS[comparator_direction]
    title = f"{trait.label}: condition-to-body-site Sankey for {focus_label} profiles"
    return replace(
        FACULTATIVE_CONFIG,
        title=title,
        focus_label=focus_label,
        comparator_label=comparator_label,
        focus_direction=focus_direction,
        comparator_direction=comparator_direction,
    )


def output_stem(trait: SharedTrait, focus_direction: str) -> str:
    """Return the base filename for one individual trait/direction map."""
    return f"{trait.slug}_{DIRECTION_LABELS[focus_direction]}_profile_sankey"


def add_trait_fields(frame: pd.DataFrame, trait: SharedTrait, analysis_key: str) -> pd.DataFrame:
    """Attach trait metadata to one edge or rollup table."""
    output = frame.copy()
    output["analysis"] = analysis_key
    output.insert(1, "trait_id", trait.trait_id)
    output.insert(2, "trait_label", trait.label)
    output.insert(3, "trait_slug", trait.slug)
    output.insert(4, "trait_color", trait.base_color)
    output.insert(5, "trait_rank", TRAIT_RANK_BY_ID[trait.trait_id])
    return output


def add_displayed_flag_by_trait(
    edges: pd.DataFrame,
    display_edges: pd.DataFrame,
) -> pd.DataFrame:
    """Return an edge table with a trait-aware displayed flag."""
    display_keys = set(
        zip(
            display_edges["trait_id"],
            display_edges["condition"],
            display_edges["body_site_group"],
            strict=True,
        )
    )
    output = edges.copy()
    output["displayed_in_map"] = [
        (trait_id, condition, body_site) in display_keys
        for trait_id, condition, body_site in zip(
            output["trait_id"],
            output["condition"],
            output["body_site_group"],
            strict=True,
        )
    ]
    return output


def summarize_rollups_from_source(
    source: pd.DataFrame,
    config,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize condition/body-site rollups from a prepared display source table."""
    condition_rollups = (
        source.groupby(
            ["raw_condition", "condition", "condition_rollup_rule"],
            dropna=False,
            sort=False,
        )
        .agg(
            profile_count=("profile_id", "count"),
            study_count=("study_id", "nunique"),
            focus_profile_count=(
                "inc_dec",
                lambda values: int((values == config.focus_direction).sum()),
            ),
            selected_guild_enriched_profile_count=(
                "has_guild_enriched_trait",
                lambda values: int(values.sum()),
            ),
        )
        .reset_index()
        .sort_values(
            ["profile_count", "study_count", "raw_condition"],
            ascending=[False, False, True],
        )
    )
    body_site_rollups = (
        source.groupby(
            [
                "raw_body_site",
                "raw_body_site_group",
                "body_site_group",
                "body_site_rollup_rule",
            ],
            dropna=False,
            sort=False,
        )
        .agg(
            profile_count=("profile_id", "count"),
            study_count=("study_id", "nunique"),
            focus_profile_count=(
                "inc_dec",
                lambda values: int((values == config.focus_direction).sum()),
            ),
            selected_guild_enriched_profile_count=(
                "has_guild_enriched_trait",
                lambda values: int(values.sum()),
            ),
        )
        .reset_index()
        .sort_values(
            ["profile_count", "study_count", "raw_body_site"],
            ascending=[False, False, True],
        )
    )
    return condition_rollups, body_site_rollups


def summarize_edges_from_source(source: pd.DataFrame, config) -> pd.DataFrame:
    """Build condition/body-site edges from a prepared display source table."""
    records: list[dict[str, object]] = []
    for (condition, body_site_group), group in source.groupby(
        ["condition", "body_site_group"],
        sort=False,
        dropna=False,
    ):
        focus = group[group["inc_dec"] == config.focus_direction]
        comparator = group[group["inc_dec"] == config.comparator_direction]
        focus_enriched = focus[focus["has_guild_enriched_trait"]]
        background_no_enriched = group[~group["has_guild_enriched_trait"]]
        focus_no_enriched = focus[~focus["has_guild_enriched_trait"]]
        comparator_enriched = comparator[comparator["has_guild_enriched_trait"]]

        focus_mean = nan_mean(focus["guild_z_mean"])
        focus_enriched_mean = nan_mean(focus_enriched["guild_z_mean"])
        comparator_mean = nan_mean(comparator["guild_z_mean"])
        signal_focus_mean = focus_enriched_mean if pd.notna(focus_enriched_mean) else focus_mean
        contrast = (
            signal_focus_mean - comparator_mean
            if pd.notna(signal_focus_mean) and pd.notna(comparator_mean)
            else signal_focus_mean
        )
        records.append(
            {
                "analysis": config.key,
                "condition": clean_text(condition),
                "body_site_group": clean_text(body_site_group),
                "raw_conditions": join_unique(group["raw_condition"]),
                "raw_body_sites": join_unique(group["raw_body_site"]),
                "raw_body_site_groups": join_unique(group["raw_body_site_group"]),
                "profile_count": len(group),
                "study_count": group["study_id"].nunique(),
                "study_ids": join_unique(group["study_id"]),
                "focus_direction": config.focus_label,
                "comparator_direction": config.comparator_label,
                "focus_profile_count": len(focus),
                "focus_study_count": focus["study_id"].nunique(),
                "focus_study_ids": join_unique(focus["study_id"]),
                "focus_enriched_profile_count": len(focus_enriched),
                "focus_enriched_study_count": focus_enriched["study_id"].nunique(),
                "focus_enriched_study_ids": join_unique(focus_enriched["study_id"]),
                "focus_no_enriched_trait_profile_count": len(focus_no_enriched),
                "comparator_profile_count": len(comparator),
                "comparator_study_count": comparator["study_id"].nunique(),
                "comparator_study_ids": join_unique(comparator["study_id"]),
                "comparator_enriched_profile_count": len(comparator_enriched),
                "background_no_enriched_profile_count": len(background_no_enriched),
                "background_no_enriched_study_count": background_no_enriched["study_id"].nunique(),
                "background_no_enriched_study_ids": join_unique(background_no_enriched["study_id"]),
                "increased_profile_count": int((group["inc_dec"] == "inc").sum()),
                "decreased_profile_count": int((group["inc_dec"] == "dec").sum()),
                "selected_guild_enriched_profile_count": int(
                    group["has_guild_enriched_trait"].sum()
                ),
                "mean_guild_z": nan_mean(group["guild_z_mean"]),
                "focus_mean_guild_z": focus_mean,
                "focus_enriched_mean_guild_z": focus_enriched_mean,
                "comparator_mean_guild_z": comparator_mean,
                "focus_enriched_minus_comparator_guild_z": contrast,
            }
        )

    if not records:
        return pd.DataFrame()

    edges = pd.DataFrame.from_records(records)
    edges = edges[edges["focus_enriched_profile_count"] > 0].copy()
    if edges.empty:
        return edges
    edges["rank_signal"] = edges["focus_enriched_minus_comparator_guild_z"].fillna(
        edges["focus_enriched_mean_guild_z"]
    )
    edges["rank_score"] = (
        edges["focus_enriched_profile_count"]
        * np.log1p(edges["focus_enriched_study_count"].clip(lower=1))
        * (1.0 + np.maximum(edges["rank_signal"], 0.0))
    )
    edges["displayed_in_map"] = False
    return edges.sort_values(
        [
            "focus_enriched_profile_count",
            "rank_score",
            "rank_signal",
            "study_count",
        ],
        ascending=[False, False, False, False],
    )


def individual_condition_grouped_body_source(
    scored: pd.DataFrame,
    config,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Keep raw conditions while preserving the grouped v3 body-site rollups."""
    source, _, _ = add_v3_display_rollups(
        scored,
        config,
        display_mode=V3_GROUPED_MODE,
    )
    source = source.copy()
    source["condition"] = source["raw_condition"].map(clean_text)
    source["condition_rollup_rule"] = (
        "kept raw individual condition label; body sites use grouped v3 rollups"
    )
    source = source[source["condition"] != "unknown"].copy()
    condition_rollups, body_site_rollups = summarize_rollups_from_source(source, config)
    return source, condition_rollups, body_site_rollups


def build_trait_direction_tables(
    metadata: pd.DataFrame,
    zscores: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build all trait-direction edge, display, and rollup tables."""
    all_edges: list[pd.DataFrame] = []
    all_display_edges: list[pd.DataFrame] = []
    all_condition_rollups: list[pd.DataFrame] = []
    all_body_site_rollups: list[pd.DataFrame] = []

    for trait in SHARED_TRAITS:
        scored = score_profiles(metadata, zscores, trait)
        for focus_direction in ["inc", "dec"]:
            focus_label = DIRECTION_LABELS[focus_direction]
            analysis_key = f"{trait.slug}_{focus_label}"
            config = single_trait_map_config(trait, focus_direction)

            _source, condition_rollups, body_site_rollups = add_v3_display_rollups(
                scored,
                config,
                display_mode=V3_GROUPED_MODE,
            )
            edges = summarize_v3_condition_body_site_edges_for_mode(
                scored,
                config,
                display_mode=V3_GROUPED_MODE,
            )
            edges = add_trait_fields(edges, trait, analysis_key)
            display_edges = select_v3_display_edges(edges)
            edge_output = add_displayed_flag_by_trait(edges, display_edges)

            condition_rollups.insert(0, "analysis", analysis_key)
            condition_rollups.insert(1, "trait_id", trait.trait_id)
            condition_rollups.insert(2, "trait_label", trait.label)
            condition_rollups.insert(3, "trait_slug", trait.slug)
            condition_rollups.insert(4, "trait_color", trait.base_color)
            condition_rollups.insert(5, "focus_direction", focus_label)
            body_site_rollups.insert(0, "analysis", analysis_key)
            body_site_rollups.insert(1, "trait_id", trait.trait_id)
            body_site_rollups.insert(2, "trait_label", trait.label)
            body_site_rollups.insert(3, "trait_slug", trait.slug)
            body_site_rollups.insert(4, "trait_color", trait.base_color)
            body_site_rollups.insert(5, "focus_direction", focus_label)

            all_edges.append(edge_output)
            all_display_edges.append(display_edges)
            all_condition_rollups.append(condition_rollups)
            all_body_site_rollups.append(body_site_rollups)

    return (
        pd.concat(all_edges, ignore_index=True),
        pd.concat(all_display_edges, ignore_index=True),
        pd.concat(all_condition_rollups, ignore_index=True),
        pd.concat(all_body_site_rollups, ignore_index=True),
    )


def build_individual_condition_direction_tables(
    metadata: pd.DataFrame,
    zscores: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build six-trait tables with individual conditions and grouped body sites."""
    all_edges: list[pd.DataFrame] = []
    all_display_edges: list[pd.DataFrame] = []
    all_condition_rollups: list[pd.DataFrame] = []
    all_body_site_rollups: list[pd.DataFrame] = []

    for trait in SHARED_TRAITS:
        scored = score_profiles(metadata, zscores, trait)
        for focus_direction in ["inc", "dec"]:
            focus_label = DIRECTION_LABELS[focus_direction]
            analysis_key = f"{trait.slug}_{focus_label}_individual_conditions"
            config = single_trait_map_config(trait, focus_direction)

            source, condition_rollups, body_site_rollups = individual_condition_grouped_body_source(
                scored, config
            )
            edges = summarize_edges_from_source(source, config)
            edges = add_trait_fields(edges, trait, analysis_key)
            display_edges = select_v3_display_edges(edges)
            edge_output = add_displayed_flag_by_trait(edges, display_edges)

            condition_rollups.insert(0, "analysis", analysis_key)
            condition_rollups.insert(1, "trait_id", trait.trait_id)
            condition_rollups.insert(2, "trait_label", trait.label)
            condition_rollups.insert(3, "trait_slug", trait.slug)
            condition_rollups.insert(4, "trait_color", trait.base_color)
            condition_rollups.insert(5, "focus_direction", focus_label)
            body_site_rollups.insert(0, "analysis", analysis_key)
            body_site_rollups.insert(1, "trait_id", trait.trait_id)
            body_site_rollups.insert(2, "trait_label", trait.label)
            body_site_rollups.insert(3, "trait_slug", trait.slug)
            body_site_rollups.insert(4, "trait_color", trait.base_color)
            body_site_rollups.insert(5, "focus_direction", focus_label)

            all_edges.append(edge_output)
            all_display_edges.append(display_edges)
            all_condition_rollups.append(condition_rollups)
            all_body_site_rollups.append(body_site_rollups)

    return (
        pd.concat(all_edges, ignore_index=True),
        pd.concat(all_display_edges, ignore_index=True),
        pd.concat(all_condition_rollups, ignore_index=True),
        pd.concat(all_body_site_rollups, ignore_index=True),
    )


def lighten_color(color: str, fraction: float = LIGHTNESS_FRACTION) -> str:
    """Return a light tint that still belongs to the same hue family."""
    rgb = to_rgb(color)
    mixed = tuple(1.0 - ((1.0 - channel) * fraction) for channel in rgb)
    return to_hex(mixed)


def value_for_color(edge: pd.Series) -> float:
    """Return the nonnegative signal value used for edge color intensity."""
    contrast = edge.get("focus_enriched_minus_comparator_guild_z")
    mean = edge.get("focus_enriched_mean_guild_z")
    value = contrast if pd.notna(contrast) else mean
    if pd.isna(value):
        return 0.0
    return max(0.0, float(value))


def parallel_offset_sequence(
    count: int,
    *,
    spacing: float = PARALLEL_EDGE_SPACING,
    max_offset: float = MAX_PARALLEL_EDGE_OFFSET,
) -> list[float]:
    """Return top-to-bottom offsets for parallel links, ordered thick to thin."""
    if count <= 1:
        return [0.0]
    center = (count - 1) / 2
    return [
        min(
            max_offset,
            max(-max_offset, (center - position) * spacing),
        )
        for position in range(count)
    ]


def add_parallel_offsets(
    display_edges: pd.DataFrame, *, separate_directions: bool = False
) -> pd.DataFrame:
    """Attach per-pair offsets ordered by edge width from largest to smallest."""
    output = display_edges.copy()
    output["parallel_rank"] = 1
    output["parallel_offset"] = 0.0
    for _key, group in output.groupby(["condition", "body_site_group"], sort=False):
        ordered_indices = group.sort_values(
            TRAIT_SORT_COLUMNS,
            ascending=[False, False, True],
        ).index.tolist()
        offsets = parallel_offset_sequence(len(ordered_indices))
        if separate_directions:
            ordered_indices = (
                group.assign(
                    direction_rank=group["focus_direction"].map({"increased": 0, "decreased": 1})
                )
                .sort_values(["trait_rank", "direction_rank"])
                .index.tolist()
            )
            offsets = (
                np.linspace(0.070, -0.070, len(ordered_indices)).tolist()
                if len(ordered_indices) > 1
                else [0.0]
            )
        for rank, (index, offset) in enumerate(
            zip(ordered_indices, offsets, strict=True),
            start=1,
        ):
            output.loc[index, "parallel_rank"] = rank
            output.loc[index, "parallel_offset"] = offset
    return output


def weighted_color_value(group: pd.DataFrame) -> float:
    """Return count-weighted color intensity for a combined bundle."""
    values = group["color_value"].astype(float)
    weights = group["focus_enriched_profile_count"].astype(float)
    if float(weights.sum()) <= 0:
        return float(values.mean()) if not values.empty else 0.0
    return float(np.average(values, weights=weights))


def body_site_trait_bundles(display_edges: pd.DataFrame) -> pd.DataFrame:
    """Aggregate displayed links into one body-site-side bundle per trait color."""
    records: list[dict[str, object]] = []
    for (body_site_group, trait_id), group in display_edges.groupby(
        ["body_site_group", "trait_id"],
        sort=False,
    ):
        records.append(
            {
                "body_site_group": body_site_group,
                "trait_id": trait_id,
                "trait_label": group["trait_label"].iloc[0],
                "trait_color": group["trait_color"].iloc[0],
                "trait_rank": int(group["trait_rank"].iloc[0]),
                "focus_enriched_profile_count": int(group["focus_enriched_profile_count"].sum()),
                "background_no_enriched_profile_count": int(
                    group["background_no_enriched_profile_count"].sum()
                ),
                "rank_score": float(group["rank_score"].max()),
                "color_value": weighted_color_value(group),
            }
        )
    bundles = pd.DataFrame.from_records(records)
    if bundles.empty:
        return bundles
    bundles["body_site_bundle_rank"] = 1
    bundles["body_site_bundle_offset"] = 0.0
    for _body_site, group in bundles.groupby("body_site_group", sort=False):
        ordered_indices = group.sort_values(
            TRAIT_SORT_COLUMNS,
            ascending=[False, False, True],
        ).index.tolist()
        offsets = parallel_offset_sequence(
            len(ordered_indices),
            spacing=BODY_SITE_BUNDLE_EDGE_SPACING,
            max_offset=MAX_BODY_SITE_BUNDLE_EDGE_OFFSET,
        )
        for rank, (index, offset) in enumerate(
            zip(ordered_indices, offsets, strict=True),
            start=1,
        ):
            bundles.loc[index, "body_site_bundle_rank"] = rank
            bundles.loc[index, "body_site_bundle_offset"] = offset
    return bundles


def offset_node_y(y: float, offset: float = 0.0) -> float:
    """Offset parallel trait edges while keeping them inside the plot area."""
    return min(0.925, max(0.035, y + offset))


def count_joined_unique(values: pd.Series) -> int:
    """Count unique semicolon-delimited IDs from already summarized edge rows."""
    unique_values: set[str] = set()
    for value in values.dropna():
        for token in str(value).split(";"):
            token = token.strip()
            if token:
                unique_values.add(token)
    return len(unique_values)


def combined_node_stats(display_edges: pd.DataFrame, node_column: str) -> pd.DataFrame:
    """Summarize combined Sankey node support over displayed trait-specific edges."""
    records: list[dict[str, object]] = []
    for label, group in display_edges.groupby(node_column, sort=False):
        records.append(
            {
                node_column: label,
                "display_trait_edge_count": len(group),
                "display_trait_count": group["trait_id"].nunique(),
                "focus_enriched_profile_count": int(group["focus_enriched_profile_count"].sum()),
                "focus_enriched_study_count": count_joined_unique(
                    group["focus_enriched_study_ids"]
                ),
            }
        )
    return pd.DataFrame.from_records(records).set_index(node_column)


def draw_combined_node(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    label: str,
    stats_row: pd.Series,
    align: str,
) -> None:
    """Draw one combined Sankey node with trait-edge membership support."""
    ax.add_patch(
        Circle(
            (x, y),
            0.012,
            facecolor=NODE_FILL,
            edgecolor=NODE_EDGE,
            linewidth=1.0,
            zorder=5,
        )
    )
    text_x = x - 0.025 if align == "right" else x + 0.025
    label_lines = [
        wrap_label(label),
        (
            f"links={int(stats_row.display_trait_edge_count)}, "
            f"traits={int(stats_row.display_trait_count)}, "
            f"n={int(stats_row.focus_enriched_profile_count):,}"
        ),
    ]
    ax.text(
        text_x,
        y,
        "\n".join(label_lines),
        ha=align,
        va="center",
        fontsize=7.6,
        color=SLATE_900,
        linespacing=1.08,
    )


def trait_color_scales(
    display_edges: pd.DataFrame,
    *,
    shared_vmax: float | None = None,
) -> dict[str, tuple[LinearSegmentedColormap, Normalize]]:
    """Build one light-to-dark monochrome scale per trait."""
    scales: dict[str, tuple[LinearSegmentedColormap, Normalize]] = {}
    shared_vmax = max(0.4, float(shared_vmax)) if shared_vmax is not None else None
    for trait in SHARED_TRAITS:
        trait_edges = display_edges[display_edges["trait_id"] == trait.trait_id]
        values = trait_edges.apply(value_for_color, axis=1) if not trait_edges.empty else []
        vmax = shared_vmax or max(0.4, float(max(values) if len(values) else 0.0))
        cmap = LinearSegmentedColormap.from_list(
            f"{trait.slug}_combined_scale",
            [lighten_color(trait.base_color), trait.base_color],
        )
        scales[trait.trait_id] = (cmap, Normalize(vmin=0.0, vmax=vmax))
    return scales


def shared_color_vmax(display_edges: pd.DataFrame) -> float:
    """Return one color-scale maximum for the supplied displayed edges."""
    if display_edges.empty:
        return 0.4
    values = display_edges.apply(value_for_color, axis=1)
    return max(0.4, float(values.max() if len(values) else 0.0))


def add_color_scale_legend(
    fig: plt.Figure,
    scales: dict[str, tuple[LinearSegmentedColormap, Normalize]],
    *,
    shared_range: bool = False,
) -> None:
    """Add six compact per-trait color scales to the right side of the figure."""
    fig.text(
        0.835,
        0.865,
        "Trait color scales",
        ha="left",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color=SLATE_700,
    )
    description = (
        "Shared numeric range for all traits.\nShade = mean positive z-score."
        if shared_range
        else "Shade = mean positive\nfocus-profile z-score."
    )
    fig.text(
        0.835,
        0.845,
        description,
        ha="left",
        va="top",
        fontsize=7.4,
        color=SLATE_500,
        linespacing=1.10,
    )
    for index, trait in enumerate(SHARED_TRAITS):
        y = 0.790 - index * 0.066
        fig.text(
            0.835,
            y + 0.022,
            trait.label,
            ha="left",
            va="bottom",
            fontsize=7.8,
            color=SLATE_700,
        )
        cax = fig.add_axes([0.835, y, 0.125, 0.010])
        cmap, norm = scales[trait.trait_id]
        cbar = fig.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cax,
            orientation="horizontal",
        )
        vmax = norm.vmax
        cbar.set_ticks([0.0, vmax])
        cbar.set_ticklabels(["0", f"{vmax:.1f}"])
        cbar.ax.tick_params(labelsize=6.5, length=2, pad=1)
        cbar.outline.set_linewidth(0.3)


def plot_single_trait_sankey(
    display_edges: pd.DataFrame,
    trait: SharedTrait,
    focus_direction: str,
    output_path: Path,
) -> None:
    """Render one grouped condition-to-body-site Sankey with a single trait scale."""
    if display_edges.empty:
        msg = f"No display edges available for {trait.label} {focus_direction}."
        raise ValueError(msg)

    focus_label = DIRECTION_LABELS[focus_direction]
    comparator_label = DIRECTION_LABELS[COMPARATOR_DIRECTIONS[focus_direction]]
    display_edges = display_edges.copy()
    display_edges["color_value"] = display_edges.apply(value_for_color, axis=1)

    condition_order = (
        display_edges.groupby("condition", sort=False)
        .agg(
            focus_enriched_profile_count=("focus_enriched_profile_count", "sum"),
            rank_score=("rank_score", "max"),
        )
        .sort_values(
            ["focus_enriched_profile_count", "rank_score"],
            ascending=[False, False],
        )
        .index.tolist()
    )
    body_site_order = (
        display_edges.groupby("body_site_group", sort=False)
        .agg(
            focus_enriched_profile_count=("focus_enriched_profile_count", "sum"),
            background_no_enriched_profile_count=(
                "background_no_enriched_profile_count",
                "sum",
            ),
        )
        .sort_values(
            ["focus_enriched_profile_count", "background_no_enriched_profile_count"],
            ascending=[False, False],
        )
        .index.tolist()
    )

    condition_y = node_positions(condition_order, top=0.90, bottom=0.06)
    body_site_y = node_positions(body_site_order, top=0.90, bottom=0.06)
    condition_stats = weighted_v3_node_stats(display_edges, "condition")
    body_site_stats = weighted_v3_node_stats(display_edges, "body_site_group")

    height = max(9.0, 1.5 + 0.66 * max(len(condition_order), len(body_site_order)))
    fig, ax = plt.subplots(figsize=(15.0, height))
    fig.subplots_adjust(left=0.05, right=0.90, top=0.90, bottom=0.08)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    values = display_edges["color_value"]
    vmax = max(0.4, float(values.max()))
    cmap = LinearSegmentedColormap.from_list(
        f"{trait.slug}_{focus_direction}_monochrome",
        [lighten_color(trait.base_color), trait.base_color],
    )
    norm = Normalize(vmin=0.0, vmax=vmax)

    all_counts = pd.concat(
        [
            display_edges["focus_enriched_profile_count"],
            display_edges["background_no_enriched_profile_count"],
        ],
        ignore_index=True,
    )
    positive_counts = all_counts[all_counts > 0]
    min_count = float(positive_counts.min())
    max_count = float(positive_counts.max())
    enriched_widths = scaled_line_width_from_range(
        display_edges["focus_enriched_profile_count"],
        min_count=min_count,
        max_count=max_count,
        min_width=1.1,
        max_width=9.5,
    )
    background_widths = scaled_line_width_from_range(
        display_edges["background_no_enriched_profile_count"],
        min_count=min_count,
        max_count=max_count,
        min_width=1.6,
        max_width=12.0,
    )

    left_x = 0.34
    right_x = 0.68
    for edge, linewidth in zip(display_edges.itertuples(), background_widths, strict=True):
        if linewidth <= 0:
            continue
        draw_curve(
            ax,
            start=(left_x, condition_y[edge.condition] + V3_BACKGROUND_EDGE_Y_OFFSET),
            end=(right_x, body_site_y[edge.body_site_group] + V3_BACKGROUND_EDGE_Y_OFFSET),
            color=EDGE_BACKGROUND,
            linewidth=float(linewidth),
            alpha=0.42,
            zorder=1.0,
        )

    for edge, linewidth in zip(display_edges.itertuples(), enriched_widths, strict=True):
        signal = max(0.0, float(edge.color_value))
        draw_curve(
            ax,
            start=(left_x, condition_y[edge.condition]),
            end=(right_x, body_site_y[edge.body_site_group]),
            color=cmap(norm(signal)),
            linewidth=float(linewidth),
            alpha=0.78,
            zorder=1.5,
        )

    for condition in condition_order:
        draw_node_v3(
            ax,
            x=left_x,
            y=condition_y[condition],
            label=condition,
            stats_row=condition_stats.loc[condition],
            align="right",
        )
    for body_site in body_site_order:
        draw_node_v3(
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
        f"{trait.label}: {focus_label} profile Sankey",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.05,
        0.925,
        wrapped_plot_text(
            f"Colored edges = {focus_label} profiles with {trait.label} z-score "
            "> 0. Parallel grey edges = same condition/body-site group profiles "
            f"with {trait.label} z-score <= 0. Color = mean positive "
            f"focus-profile z-score; TSV includes {focus_label} - "
            f"{comparator_label} contrast. Node key: n = positive-z profiles, "
            "s = unique studies for n, bg = z <= 0 background profiles."
        ),
        ha="left",
        va="top",
        fontsize=10,
        color=SLATE_700,
    )
    ax.text(
        0.21,
        0.965,
        "Condition display groups",
        ha="center",
        va="center",
        fontsize=12,
        color=SLATE_700,
    )
    ax.text(
        0.79,
        0.965,
        "Body-site display groups",
        ha="center",
        va="center",
        fontsize=12,
        color=SLATE_700,
    )

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.030, pad=0.02, shrink=0.55)
    cbar.set_label(f"Mean positive {focus_label} z-score", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    representative_counts = sorted(
        {
            int(display_edges["focus_enriched_profile_count"].min()),
            int(display_edges["focus_enriched_profile_count"].median()),
            int(display_edges["focus_enriched_profile_count"].max()),
        }
    )
    width_lookup = scaled_line_width_from_range(
        pd.Series(representative_counts, index=representative_counts),
        min_count=min_count,
        max_count=max_count,
        min_width=1.1,
        max_width=9.5,
    )
    legend_lines = [
        Line2D(
            [0],
            [0],
            color=trait.base_color,
            linewidth=float(width_lookup.loc[count]),
            alpha=0.78,
            solid_capstyle="round",
            label=f"n={count} profiles",
        )
        for count in representative_counts
    ]
    legend_lines.append(
        Line2D(
            [0],
            [0],
            color=EDGE_BACKGROUND,
            linewidth=6.5,
            alpha=0.42,
            solid_capstyle="round",
            label="bg profiles: z <= 0",
        )
    )
    legend = ax.legend(
        handles=legend_lines,
        title=f"Edge width: {focus_label} positive-z profiles (n)",
        loc="lower center",
        bbox_to_anchor=(0.50, -0.02),
        ncol=min(len(legend_lines), 4),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    ax.add_artist(legend)

    fig.text(
        0.05,
        0.035,
        (
            "Grouped display assigns every condition and body site to v3 rollup labels. "
            f"Displayed colored edges require at least {V3_MIN_FOCUS_ENRICHED_PROFILES} "
            f"{focus_label} positive-z profiles and {MIN_STUDIES} studies where possible."
        ),
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=SLATE_500,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def draw_lane_node(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    label: str,
    stats_row: pd.Series,
    align: str,
) -> None:
    """Draw a compact node for a trait-lane panel."""
    ax.add_patch(
        Circle(
            (x, y),
            0.010,
            facecolor=NODE_FILL,
            edgecolor=NODE_EDGE,
            linewidth=0.9,
            zorder=5,
        )
    )
    text_x = x - 0.020 if align == "right" else x + 0.020
    text_lines = [
        wrap_label(label, width=28),
        (
            f"n={int(stats_row.focus_enriched_profile_count)}, "
            f"s={int(stats_row.focus_enriched_study_count)}, "
            f"bg={int(stats_row.background_no_enriched_profile_count)}"
        ),
    ]
    ax.text(
        text_x,
        y,
        "\n".join(text_lines),
        ha=align,
        va="center",
        fontsize=6.3,
        color=SLATE_900,
        linespacing=1.03,
    )


def plot_trait_lane_direction_sankey(
    display_edges: pd.DataFrame,
    focus_direction: str,
    output_path: Path,
) -> None:
    """Render one six-lane Sankey figure that preserves all displayed trait links."""
    if display_edges.empty:
        msg = f"No display edges available for {DIRECTION_LABELS[focus_direction]} lanes."
        raise ValueError(msg)

    focus_label = DIRECTION_LABELS[focus_direction]
    display_edges = display_edges.copy()
    display_edges["color_value"] = display_edges.apply(value_for_color, axis=1)

    condition_order = (
        display_edges.groupby("condition", sort=False)
        .agg(
            display_trait_count=("trait_id", "nunique"),
            focus_enriched_profile_count=("focus_enriched_profile_count", "sum"),
            rank_score=("rank_score", "max"),
        )
        .sort_values(
            ["display_trait_count", "focus_enriched_profile_count", "rank_score"],
            ascending=[False, False, False],
        )
        .index.tolist()
    )
    body_site_order = (
        display_edges.groupby("body_site_group", sort=False)
        .agg(
            display_trait_count=("trait_id", "nunique"),
            focus_enriched_profile_count=("focus_enriched_profile_count", "sum"),
            background_no_enriched_profile_count=(
                "background_no_enriched_profile_count",
                "sum",
            ),
        )
        .sort_values(
            [
                "display_trait_count",
                "focus_enriched_profile_count",
                "background_no_enriched_profile_count",
            ],
            ascending=[False, False, False],
        )
        .index.tolist()
    )

    counts = pd.concat(
        [
            display_edges["focus_enriched_profile_count"],
            display_edges["background_no_enriched_profile_count"],
        ],
        ignore_index=True,
    )
    positive_counts = counts[counts > 0]
    min_count = float(positive_counts.min())
    max_count = float(positive_counts.max())
    scales = trait_color_scales(display_edges)

    fig, axes = plt.subplots(
        nrows=len(SHARED_TRAITS),
        ncols=1,
        figsize=(17.5, 32.0),
    )
    fig.subplots_adjust(left=0.06, right=0.80, top=0.925, bottom=0.060, hspace=0.44)
    fig.patch.set_facecolor("white")

    fig.text(
        0.04,
        0.985,
        f"Six shared traits: {focus_label} profile Sankey, trait-lane version",
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.04,
        0.965,
        wrapped_plot_text(
            "Each horizontal lane is one trait. This preserves the same displayed "
            "trait-specific links as the full overlay while separating crossings "
            "by trait. Node key: n = positive-z profiles, s = unique studies "
            "for n, bg = z <= 0 background profiles."
        ),
        ha="left",
        va="top",
        fontsize=9.6,
        color=SLATE_700,
    )

    left_x = 0.35
    right_x = 0.66
    for ax, trait in zip(axes, SHARED_TRAITS, strict=True):
        trait_edges = display_edges[display_edges["trait_id"] == trait.trait_id].copy()
        trait_edges = trait_edges.sort_values(
            ["focus_enriched_profile_count", "rank_score"],
            ascending=[False, False],
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.set_facecolor("white")

        ax.text(
            0.02,
            0.965,
            trait.label,
            ha="left",
            va="top",
            fontsize=10.8,
            fontweight="bold",
            color=trait.base_color,
        )
        ax.text(
            0.02,
            0.885,
            f"{len(trait_edges)} displayed links",
            ha="left",
            va="top",
            fontsize=7.3,
            color=SLATE_500,
        )
        ax.text(
            0.21,
            0.955,
            "Conditions",
            ha="center",
            va="center",
            fontsize=8.6,
            color=SLATE_700,
        )
        ax.text(
            0.79,
            0.955,
            "Body sites",
            ha="center",
            va="center",
            fontsize=8.6,
            color=SLATE_700,
        )

        if trait_edges.empty:
            continue

        lane_conditions = [
            condition for condition in condition_order if condition in set(trait_edges["condition"])
        ]
        lane_body_sites = [
            body_site
            for body_site in body_site_order
            if body_site in set(trait_edges["body_site_group"])
        ]
        condition_y = node_positions(lane_conditions, top=0.86, bottom=0.06)
        body_site_y = node_positions(lane_body_sites, top=0.86, bottom=0.06)
        condition_stats = weighted_v3_node_stats(trait_edges, "condition")
        body_site_stats = weighted_v3_node_stats(trait_edges, "body_site_group")
        enriched_widths = scaled_line_width_from_range(
            trait_edges["focus_enriched_profile_count"],
            min_count=min_count,
            max_count=max_count,
            min_width=0.65,
            max_width=5.5,
        )
        background_widths = scaled_line_width_from_range(
            trait_edges["background_no_enriched_profile_count"],
            min_count=min_count,
            max_count=max_count,
            min_width=0.8,
            max_width=6.5,
        )
        cmap, norm = scales[trait.trait_id]

        for edge, linewidth in zip(trait_edges.itertuples(), background_widths, strict=True):
            if linewidth <= 0:
                continue
            draw_curve(
                ax,
                start=(
                    left_x,
                    condition_y[edge.condition] + V3_BACKGROUND_EDGE_Y_OFFSET,
                ),
                end=(
                    right_x,
                    body_site_y[edge.body_site_group] + V3_BACKGROUND_EDGE_Y_OFFSET,
                ),
                color=EDGE_BACKGROUND,
                linewidth=float(linewidth),
                alpha=0.34,
                zorder=1.0,
            )

        for edge, linewidth in zip(trait_edges.itertuples(), enriched_widths, strict=True):
            signal = max(0.0, float(edge.color_value))
            draw_curve(
                ax,
                start=(left_x, condition_y[edge.condition]),
                end=(right_x, body_site_y[edge.body_site_group]),
                color=cmap(norm(signal)),
                linewidth=float(linewidth),
                alpha=0.78,
                zorder=1.8,
            )

        for condition in lane_conditions:
            draw_lane_node(
                ax,
                x=left_x,
                y=condition_y[condition],
                label=condition,
                stats_row=condition_stats.loc[condition],
                align="right",
            )
        for body_site in lane_body_sites:
            draw_lane_node(
                ax,
                x=right_x,
                y=body_site_y[body_site],
                label=body_site,
                stats_row=body_site_stats.loc[body_site],
                align="left",
            )

        position = ax.get_position()
        cax = fig.add_axes([0.835, position.y0 + position.height * 0.41, 0.115, 0.007])
        cbar = fig.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cax,
            orientation="horizontal",
        )
        cbar.set_ticks([0.0, norm.vmax])
        cbar.set_ticklabels(["0", f"{norm.vmax:.1f}"])
        cbar.ax.tick_params(labelsize=6.2, length=2, pad=1)
        cbar.outline.set_linewidth(0.3)
        fig.text(
            0.835,
            position.y0 + position.height * 0.50,
            "mean positive z",
            ha="left",
            va="bottom",
            fontsize=6.5,
            color=SLATE_500,
        )

    representative_counts = sorted(
        {
            int(display_edges["focus_enriched_profile_count"].min()),
            int(display_edges["focus_enriched_profile_count"].median()),
            int(display_edges["focus_enriched_profile_count"].max()),
        }
    )
    width_lookup = scaled_line_width_from_range(
        pd.Series(representative_counts, index=representative_counts),
        min_count=min_count,
        max_count=max_count,
        min_width=0.65,
        max_width=5.5,
    )
    legend_lines = [
        Line2D(
            [0],
            [0],
            color=SLATE_700,
            linewidth=float(width_lookup.loc[count]),
            alpha=0.78,
            solid_capstyle="round",
            label=f"n={count} profiles",
        )
        for count in representative_counts
    ]
    legend_lines.append(
        Line2D(
            [0],
            [0],
            color=EDGE_BACKGROUND,
            linewidth=5.8,
            alpha=0.34,
            solid_capstyle="round",
            label="bg profiles: z <= 0",
        )
    )
    fig.legend(
        handles=legend_lines,
        title=f"Edge width: {focus_label} positive-z profiles per trait link (n)",
        loc="lower center",
        bbox_to_anchor=(0.49, 0.020),
        ncol=min(len(legend_lines), 4),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def select_shared_high_support_display_edges(
    display_edges: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select a smaller overlay focused on shared/high-support condition-body pairs."""
    pair_stats = (
        display_edges.groupby(["condition", "body_site_group"], sort=False)
        .agg(
            shared_pattern_trait_count=("trait_id", "nunique"),
            shared_pattern_link_count=("trait_id", "count"),
            shared_pattern_focus_profile_count=("focus_enriched_profile_count", "sum"),
            shared_pattern_max_focus_profile_count=(
                "focus_enriched_profile_count",
                "max",
            ),
            shared_pattern_max_rank_score=("rank_score", "max"),
        )
        .reset_index()
    )
    candidates = pair_stats[
        (pair_stats["shared_pattern_trait_count"] >= SHARED_PATTERN_MIN_TRAITS)
        | (
            pair_stats["shared_pattern_max_focus_profile_count"]
            >= SHARED_PATTERN_HIGH_SUPPORT_MIN_PROFILES
        )
    ].copy()
    candidates = candidates.sort_values(
        [
            "shared_pattern_trait_count",
            "shared_pattern_focus_profile_count",
            "shared_pattern_max_focus_profile_count",
            "shared_pattern_max_rank_score",
        ],
        ascending=[False, False, False, False],
    ).head(SHARED_PATTERN_MAX_PAIRS)
    candidates.insert(0, "shared_pattern_rank", range(1, len(candidates) + 1))

    keys = ["condition", "body_site_group"]
    filtered_edges = display_edges.merge(candidates, on=keys, how="inner")
    filtered_edges = filtered_edges.sort_values(
        ["shared_pattern_rank", "trait_rank"],
        ascending=[True, True],
    )
    return filtered_edges, candidates


def plot_combined_direction_sankey(
    display_edges: pd.DataFrame,
    focus_direction: str,
    output_path: Path,
    *,
    title_prefix: str = "Six shared traits",
    subtitle_extra: str = "",
    footer_extra: str = "",
    condition_axis_label: str = "Condition display groups",
    body_site_axis_label: str = "Body-site display groups",
    display_filter_note: str | None = None,
    color_scale_vmax: float | None = None,
    subplot_top: float = 0.89,
    bundle_by_trait_at_body_site: bool = False,
) -> None:
    """Render one combined condition-to-body-site Sankey for all six traits."""
    both_directions = focus_direction == "both"
    if both_directions and bundle_by_trait_at_body_site:
        raise ValueError("The bidirectional overlay requires separate direction links.")
    if display_edges.empty:
        msg = f"No display edges available for {focus_direction}."
        raise ValueError(msg)

    focus_label = "increased + decreased" if both_directions else DIRECTION_LABELS[focus_direction]
    comparator_label = (
        "" if both_directions else DIRECTION_LABELS[COMPARATOR_DIRECTIONS[focus_direction]]
    )
    display_edges = display_edges.copy()
    display_edges["color_value"] = (
        display_edges["focus_enriched_mean_guild_z"].fillna(0.0).clip(lower=0.0)
        if both_directions
        else display_edges.apply(value_for_color, axis=1)
    )
    display_edges = add_parallel_offsets(display_edges, separate_directions=both_directions)
    bundles = pd.DataFrame()
    if bundle_by_trait_at_body_site:
        bundles = body_site_trait_bundles(display_edges)
        display_edges = display_edges.merge(
            bundles[
                [
                    "body_site_group",
                    "trait_id",
                    "body_site_bundle_rank",
                    "body_site_bundle_offset",
                ]
            ],
            on=["body_site_group", "trait_id"],
            how="left",
        )

    condition_order_frame = (
        display_edges.groupby("condition", sort=False)
        .agg(
            focus_enriched_profile_count=("focus_enriched_profile_count", "sum"),
            display_trait_count=("trait_id", "nunique"),
            display_trait_edge_count=("trait_id", "count"),
            rank_score=("rank_score", "max"),
        )
        .sort_values(
            [
                "display_trait_count",
                "focus_enriched_profile_count",
                "display_trait_edge_count",
                "rank_score",
            ],
            ascending=[False, False, False, False],
        )
    )
    body_site_order_frame = (
        display_edges.groupby("body_site_group", sort=False)
        .agg(
            focus_enriched_profile_count=("focus_enriched_profile_count", "sum"),
            display_trait_count=("trait_id", "nunique"),
            display_trait_edge_count=("trait_id", "count"),
            background_no_enriched_profile_count=(
                "background_no_enriched_profile_count",
                "sum",
            ),
        )
        .sort_values(
            [
                "display_trait_count",
                "focus_enriched_profile_count",
                "display_trait_edge_count",
                "background_no_enriched_profile_count",
            ],
            ascending=[False, False, False, False],
        )
    )
    condition_order = condition_order_frame.index.tolist()
    body_site_order = body_site_order_frame.index.tolist()

    condition_y = node_positions(condition_order, top=0.89, bottom=0.06)
    body_site_y = node_positions(body_site_order, top=0.89, bottom=0.06)
    condition_stats = combined_node_stats(display_edges, "condition")
    body_site_stats = combined_node_stats(display_edges, "body_site_group")

    height = max(11.0, 2.2 + 0.54 * max(len(condition_order), len(body_site_order)))
    if both_directions:
        height = max(13.0, height)
    fig, ax = plt.subplots(figsize=(17.5, height))
    fig.subplots_adjust(left=0.04, right=0.80, top=subplot_top, bottom=0.105)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    source_counts = pd.concat(
        [
            display_edges["focus_enriched_profile_count"],
            display_edges["background_no_enriched_profile_count"],
        ],
        ignore_index=True,
    )
    source_positive_counts = source_counts[source_counts > 0]
    if both_directions:
        source_positive_counts = display_edges["focus_enriched_profile_count"]
    source_min_count = float(source_positive_counts.min())
    source_max_count = float(source_positive_counts.max())
    source_max_enriched_width = 4.6 if bundle_by_trait_at_body_site else 5.7
    source_max_background_width = 5.1 if bundle_by_trait_at_body_site else 6.6
    enriched_widths = scaled_line_width_from_range(
        display_edges["focus_enriched_profile_count"],
        min_count=source_min_count,
        max_count=source_max_count,
        min_width=0.75,
        max_width=source_max_enriched_width,
    )
    background_widths = scaled_line_width_from_range(
        display_edges["background_no_enriched_profile_count"],
        min_count=source_min_count,
        max_count=source_max_count,
        min_width=0.9,
        max_width=source_max_background_width,
    )
    if both_directions:
        background_widths[:] = 0.0
    bundle_enriched_widths = pd.Series(dtype=float)
    bundle_background_widths = pd.Series(dtype=float)
    if bundle_by_trait_at_body_site:
        bundle_counts = pd.concat(
            [
                bundles["focus_enriched_profile_count"],
                bundles["background_no_enriched_profile_count"],
            ],
            ignore_index=True,
        )
        bundle_positive_counts = bundle_counts[bundle_counts > 0]
        bundle_min_count = float(bundle_positive_counts.min())
        bundle_max_count = float(bundle_positive_counts.max())
        bundle_enriched_widths = scaled_line_width_from_range(
            bundles["focus_enriched_profile_count"],
            min_count=bundle_min_count,
            max_count=bundle_max_count,
            min_width=1.6,
            max_width=10.0,
        )
        bundle_background_widths = scaled_line_width_from_range(
            bundles["background_no_enriched_profile_count"],
            min_count=bundle_min_count,
            max_count=bundle_max_count,
            min_width=2.0,
            max_width=11.0,
        )
    scales = trait_color_scales(display_edges, shared_vmax=color_scale_vmax)

    left_x = 0.33
    right_x = 0.64
    source_end_x = BODY_SITE_BUNDLE_X if bundle_by_trait_at_body_site else right_x
    ordered = display_edges.sort_values(
        [
            "condition",
            "body_site_group",
            "parallel_rank",
            "focus_enriched_profile_count",
            "rank_score",
            "trait_rank",
        ],
        ascending=[True, True, True, False, False, True],
    )
    enriched_widths = enriched_widths.loc[ordered.index]
    background_widths = background_widths.loc[ordered.index]

    for edge, linewidth in zip(ordered.itertuples(), background_widths, strict=True):
        if linewidth <= 0:
            continue
        y0 = offset_node_y(condition_y[edge.condition], edge.parallel_offset)
        body_offset = (
            edge.body_site_bundle_offset if bundle_by_trait_at_body_site else edge.parallel_offset
        )
        y1 = offset_node_y(body_site_y[edge.body_site_group], body_offset)
        draw_curve(
            ax,
            start=(left_x, y0),
            end=(source_end_x, y1),
            color=EDGE_BACKGROUND,
            linewidth=float(linewidth),
            alpha=0.12 if bundle_by_trait_at_body_site else 0.18,
            zorder=1.0,
        )

    if bundle_by_trait_at_body_site:
        ordered_bundles = bundles.sort_values(
            [
                "body_site_group",
                "body_site_bundle_rank",
                "focus_enriched_profile_count",
                "rank_score",
                "trait_rank",
            ],
            ascending=[True, True, False, False, True],
        )
        bundle_background_widths = bundle_background_widths.loc[ordered_bundles.index]
        for bundle, linewidth in zip(
            ordered_bundles.itertuples(),
            bundle_background_widths,
            strict=True,
        ):
            if linewidth <= 0:
                continue
            y = offset_node_y(
                body_site_y[bundle.body_site_group],
                bundle.body_site_bundle_offset,
            )
            draw_curve(
                ax,
                start=(BODY_SITE_BUNDLE_X, y),
                end=(right_x, y),
                color=EDGE_BACKGROUND,
                linewidth=float(linewidth),
                alpha=0.18,
                zorder=1.1,
            )

    for edge, linewidth in zip(ordered.itertuples(), enriched_widths, strict=True):
        signal = max(0.0, float(edge.color_value))
        cmap, norm = scales[edge.trait_id]
        y0 = offset_node_y(condition_y[edge.condition], edge.parallel_offset)
        body_offset = (
            edge.body_site_bundle_offset if bundle_by_trait_at_body_site else edge.parallel_offset
        )
        y1 = offset_node_y(body_site_y[edge.body_site_group], body_offset)
        if both_directions:
            y0 = condition_y[edge.condition]
            y1 = body_site_y[edge.body_site_group]
        draw_curve(
            ax,
            start=(left_x, y0),
            end=(source_end_x, y1),
            color=cmap(norm(signal)),
            linewidth=float(linewidth),
            alpha=0.92 if both_directions else (0.46 if bundle_by_trait_at_body_site else 0.72),
            zorder=2.0 + ((len(SHARED_TRAITS) - int(edge.parallel_rank)) * 0.03),
            linestyle="--" if both_directions and edge.focus_direction == "decreased" else "-",
            capstyle="butt" if both_directions else "round",
            control_y_offset=edge.parallel_offset if both_directions else 0.0,
        )

    if bundle_by_trait_at_body_site:
        bundle_enriched_widths = bundle_enriched_widths.loc[ordered_bundles.index]
        for bundle, linewidth in zip(
            ordered_bundles.itertuples(),
            bundle_enriched_widths,
            strict=True,
        ):
            signal = max(0.0, float(bundle.color_value))
            cmap, norm = scales[bundle.trait_id]
            y = offset_node_y(
                body_site_y[bundle.body_site_group],
                bundle.body_site_bundle_offset,
            )
            draw_curve(
                ax,
                start=(BODY_SITE_BUNDLE_X, y),
                end=(right_x, y),
                color=cmap(norm(signal)),
                linewidth=float(linewidth),
                alpha=0.88,
                zorder=3.0 + ((len(SHARED_TRAITS) - int(bundle.body_site_bundle_rank)) * 0.04),
            )

    for condition in condition_order:
        draw_combined_node(
            ax,
            x=left_x,
            y=condition_y[condition],
            label=condition,
            stats_row=condition_stats.loc[condition],
            align="right",
        )
    for body_site in body_site_order:
        draw_combined_node(
            ax,
            x=right_x,
            y=body_site_y[body_site],
            label=body_site,
            stats_row=body_site_stats.loc[body_site],
            align="left",
        )

    fig.text(
        0.04,
        0.965,
        f"{title_prefix}: {focus_label} profile Sankey",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.04,
        0.925,
        wrapped_plot_text(
            "Solid = increased profiles; dashed = decreased profiles. "
            "Color identifies the trait; width = positive-z profile count; "
            "shade = mean positive z-score. "
            "Each direction retains its own top-12 pair selection. "
            "Node n sums profile-trait memberships across both directions."
            if both_directions
            else (
                f"Colored curves are {focus_label} profiles with trait z-score > 0; "
                "parallel trait colors distinguish the six q <= 0.001 traits. "
                "Width = positive-z profile count; shade = mean positive z-score "
                "on a shared numeric scale for all six traits. Within each parallel "
                "condition/body-site stack, curves are ordered thick-to-thin. "
                "Light grey underlays show z <= 0 background for the same "
                f"condition/body-site/trait; TSVs retain the {focus_label} - "
                f"{comparator_label} contrast. Node key: links = plotted "
                "trait-specific links, traits = distinct traits, n = positive-z "
                "profile memberships summed across displayed trait links."
                + (
                    " Body-site trunks sum same-color links."
                    if bundle_by_trait_at_body_site
                    else ""
                )
                + (f" {subtitle_extra}" if subtitle_extra else "")
            )
        ),
        ha="left",
        va="top",
        fontsize=9.6,
        color=SLATE_700,
    )
    ax.text(
        0.18,
        0.965,
        condition_axis_label,
        ha="center",
        va="center",
        fontsize=11.5,
        color=SLATE_700,
    )
    ax.text(
        0.77,
        0.965,
        body_site_axis_label,
        ha="center",
        va="center",
        fontsize=11.5,
        color=SLATE_700,
    )

    representative_counts = sorted(
        {
            int(display_edges["focus_enriched_profile_count"].min()),
            int(display_edges["focus_enriched_profile_count"].median()),
            int(display_edges["focus_enriched_profile_count"].max()),
        }
    )
    width_lookup = scaled_line_width_from_range(
        pd.Series(representative_counts, index=representative_counts),
        min_count=source_min_count,
        max_count=source_max_count,
        min_width=0.75,
        max_width=source_max_enriched_width,
    )
    legend_lines = [
        Line2D(
            [0],
            [0],
            color=SLATE_700,
            linewidth=float(width_lookup.loc[count]),
            alpha=0.72,
            solid_capstyle="round",
            label=f"n={count} profiles",
        )
        for count in representative_counts
    ]
    if bundle_by_trait_at_body_site:
        legend_lines.append(
            Line2D(
                [0],
                [0],
                color=SLATE_700,
                linewidth=8.0,
                alpha=0.88,
                solid_capstyle="round",
                label="body-site trunk: summed n",
            )
        )
    background_legend = Line2D(
        [0],
        [0],
        color=EDGE_BACKGROUND,
        linewidth=5.8,
        alpha=0.28,
        solid_capstyle="round",
        label="bg profiles: z <= 0",
    )
    if not both_directions:
        legend_lines.append(background_legend)
    else:
        fig.legend(
            handles=[
                Line2D([0], [0], color=SLATE_700, linewidth=2.5, label="Increased (solid)"),
                Line2D(
                    [0],
                    [0],
                    color=SLATE_700,
                    linewidth=2.5,
                    linestyle="--",
                    dash_capstyle="butt",
                    label="Decreased (dashed)",
                ),
            ],
            title="Profile direction",
            loc="upper left",
            bbox_to_anchor=(0.83, 0.39),
            frameon=False,
            fontsize=9,
            title_fontsize=10,
            handlelength=3.5,
        )
    legend = ax.legend(
        handles=legend_lines,
        title=(
            f"Edge width: {focus_label} positive-z profiles per trait link (n)"
            if not bundle_by_trait_at_body_site
            else "Edge width: source n; body-site trunks sum n per trait"
        ),
        loc="lower center",
        bbox_to_anchor=(0.49, -0.035),
        ncol=min(len(legend_lines), 4),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    ax.add_artist(legend)
    add_color_scale_legend(fig, scales, shared_range=color_scale_vmax is not None)

    if display_filter_note is None:
        display_filter_note = (
            "Display filtering uses the existing grouped v3 Sankey rules: at "
            f"least {V3_MIN_FOCUS_ENRICHED_PROFILES} positive-z focus profiles "
            f"and {MIN_STUDIES} studies where possible."
        )
    fig.text(
        0.04,
        0.035,
        wrapped_plot_text(
            "Every line is trait-specific. Node counts are profile-trait "
            "memberships across displayed links, not deduplicated profile totals. "
            f"{display_filter_note}"
            + (
                " Body-site trunks sum all displayed condition links by trait "
                "color; each body site can show up to six trunks."
                if bundle_by_trait_at_body_site
                else ""
            )
            + (f" {footer_extra}" if footer_extra else "")
        ),
        ha="left",
        va="bottom",
        fontsize=8.1,
        color=SLATE_500,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_trait_summary(
    metadata: pd.DataFrame,
    raw_trait_counts: pd.DataFrame,
    zscores: pd.DataFrame,
    q_stats: pd.DataFrame,
) -> pd.DataFrame:
    """Build a compact audit table for the six selected traits."""
    records: list[dict[str, object]] = []
    for trait in SHARED_TRAITS:
        counts = raw_trait_counts[trait.trait_id]
        z = zscores[trait.trait_id]
        q_row = q_stats.loc[trait.trait_id]
        record = {
            "trait_id": trait.trait_id,
            "trait_label": trait.label,
            "slug": trait.slug,
            "base_color": trait.base_color,
            "paired_q_value": q_row["q_value"],
            "paired_mean_increased_minus_decreased": q_row["mean_increased_minus_decreased"],
            "paired_direction": q_row["direction"],
            "profile_count": len(metadata),
            "increased_profile_count": int((metadata["inc_dec"] == "inc").sum()),
            "decreased_profile_count": int((metadata["inc_dec"] == "dec").sum()),
            "annotated_profile_count": int((counts > 0).sum()),
            "increased_annotated_profile_count": int(
                ((metadata["inc_dec"] == "inc") & (counts > 0)).sum()
            ),
            "decreased_annotated_profile_count": int(
                ((metadata["inc_dec"] == "dec") & (counts > 0)).sum()
            ),
            "positive_z_profile_count": int((z > 0).sum()),
            "increased_positive_z_profile_count": int(
                ((metadata["inc_dec"] == "inc") & (z > 0)).sum()
            ),
            "decreased_positive_z_profile_count": int(
                ((metadata["inc_dec"] == "dec") & (z > 0)).sum()
            ),
        }
        records.append(record)
    return pd.DataFrame.from_records(records)


def build_trait_selection(q_stats: pd.DataFrame) -> pd.DataFrame:
    """Write a minimal trait-selection audit table."""
    records: list[dict[str, object]] = []
    for trait in SHARED_TRAITS:
        row = q_stats.loc[trait.trait_id]
        records.append(
            {
                "trait_id": trait.trait_id,
                "trait_name": trait.label,
                "slug": trait.slug,
                "base_color": trait.base_color,
                "paired_q_value": row["q_value"],
                "paired_p_value": row["p_value"],
                "paired_mean_increased_minus_decreased": row["mean_increased_minus_decreased"],
                "paired_direction": row["direction"],
            }
        )
    return pd.DataFrame.from_records(records)


def write_direction_outputs(
    all_edges: pd.DataFrame,
    all_display_edges: pd.DataFrame,
    condition_rollups: pd.DataFrame,
    body_site_rollups: pd.DataFrame,
) -> pd.DataFrame:
    """Write the two combined Sankeys and their direction-specific TSVs."""
    output_records: list[dict[str, object]] = []
    color_vmax = shared_color_vmax(all_display_edges)
    for focus_direction, focus_label in DIRECTION_LABELS.items():
        analysis_key = f"six_shared_traits_{focus_label}"
        stem = f"{analysis_key}_profile_sankey"
        edges = all_edges[all_edges["focus_direction"] == focus_label].copy()
        display_edges = all_display_edges[
            all_display_edges["focus_direction"] == focus_label
        ].copy()
        condition_subset = condition_rollups[
            condition_rollups["focus_direction"] == focus_label
        ].copy()
        body_site_subset = body_site_rollups[
            body_site_rollups["focus_direction"] == focus_label
        ].copy()

        edges_path = OUTPUT_DIR / f"{stem}_edges.tsv"
        display_edges_path = OUTPUT_DIR / f"{stem}_display_edges.tsv"
        condition_rollups_path = OUTPUT_DIR / f"{stem}_condition_rollups.tsv"
        body_site_rollups_path = OUTPUT_DIR / f"{stem}_body_site_rollups.tsv"
        png_path = OUTPUT_DIR / f"{stem}.png"

        edges.to_csv(edges_path, sep="\t", index=False)
        display_edges.to_csv(display_edges_path, sep="\t", index=False)
        condition_subset.to_csv(condition_rollups_path, sep="\t", index=False)
        body_site_subset.to_csv(body_site_rollups_path, sep="\t", index=False)
        plot_combined_direction_sankey(
            display_edges,
            focus_direction,
            png_path,
            color_scale_vmax=color_vmax,
            subplot_top=0.855,
            bundle_by_trait_at_body_site=True,
        )

        output_records.append(
            {
                "analysis": analysis_key,
                "figure_type": "combined_direction",
                "trait_id": "all",
                "trait_label": "all six shared traits",
                "trait_slug": "all",
                "focus_direction": focus_label,
                "trait_count": display_edges["trait_id"].nunique(),
                "edge_count": len(edges),
                "display_edge_count": len(display_edges),
                "png_path": png_path.relative_to(REPO_ROOT).as_posix(),
                "edges_path": edges_path.relative_to(REPO_ROOT).as_posix(),
                "display_edges_path": display_edges_path.relative_to(REPO_ROOT).as_posix(),
                "condition_rollups_path": condition_rollups_path.relative_to(REPO_ROOT).as_posix(),
                "body_site_rollups_path": body_site_rollups_path.relative_to(REPO_ROOT).as_posix(),
            }
        )
    return pd.DataFrame.from_records(output_records)


def write_trait_lane_outputs(
    all_edges: pd.DataFrame,
    all_display_edges: pd.DataFrame,
) -> pd.DataFrame:
    """Write two combined six-lane Sankey figures."""
    output_records: list[dict[str, object]] = []
    for focus_direction, focus_label in DIRECTION_LABELS.items():
        analysis_key = f"six_shared_traits_{focus_label}_trait_lanes"
        stem = f"{analysis_key}_profile_sankey"
        edges = all_edges[all_edges["focus_direction"] == focus_label].copy()
        display_edges = all_display_edges[
            all_display_edges["focus_direction"] == focus_label
        ].copy()

        display_edges_path = OUTPUT_DIR / f"{stem}_display_edges.tsv"
        png_path = OUTPUT_DIR / f"{stem}.png"

        display_edges.to_csv(display_edges_path, sep="\t", index=False)
        plot_trait_lane_direction_sankey(display_edges, focus_direction, png_path)

        output_records.append(
            {
                "analysis": analysis_key,
                "figure_type": "combined_trait_lanes",
                "trait_id": "all",
                "trait_label": "all six shared traits",
                "trait_slug": "all",
                "focus_direction": focus_label,
                "trait_count": display_edges["trait_id"].nunique(),
                "edge_count": len(edges),
                "display_edge_count": len(display_edges),
                "png_path": png_path.relative_to(REPO_ROOT).as_posix(),
                "edges_path": (
                    OUTPUT_DIR / f"six_shared_traits_{focus_label}_profile_sankey_edges.tsv"
                )
                .relative_to(REPO_ROOT)
                .as_posix(),
                "display_edges_path": display_edges_path.relative_to(REPO_ROOT).as_posix(),
                "condition_rollups_path": (
                    OUTPUT_DIR
                    / f"six_shared_traits_{focus_label}_profile_sankey_condition_rollups.tsv"
                )
                .relative_to(REPO_ROOT)
                .as_posix(),
                "body_site_rollups_path": (
                    OUTPUT_DIR
                    / f"six_shared_traits_{focus_label}_profile_sankey_body_site_rollups.tsv"
                )
                .relative_to(REPO_ROOT)
                .as_posix(),
            }
        )
    return pd.DataFrame.from_records(output_records)


def write_shared_high_support_outputs(
    all_edges: pd.DataFrame,
    all_display_edges: pd.DataFrame,
) -> pd.DataFrame:
    """Write two reduced overlays focused on shared/high-support links."""
    output_records: list[dict[str, object]] = []
    color_vmax = shared_color_vmax(all_display_edges)
    for focus_direction, focus_label in DIRECTION_LABELS.items():
        analysis_key = f"six_shared_traits_{focus_label}_shared_high_support"
        stem = f"{analysis_key}_profile_sankey"
        direction_edges = all_edges[all_edges["focus_direction"] == focus_label].copy()
        direction_display_edges = all_display_edges[
            all_display_edges["focus_direction"] == focus_label
        ].copy()
        display_edges, pair_selection = select_shared_high_support_display_edges(
            direction_display_edges
        )

        selected_keys = pair_selection[
            [
                "condition",
                "body_site_group",
                "shared_pattern_rank",
                "shared_pattern_trait_count",
                "shared_pattern_focus_profile_count",
                "shared_pattern_max_focus_profile_count",
            ]
        ]
        edges = direction_edges.merge(
            selected_keys,
            on=["condition", "body_site_group"],
            how="inner",
        )
        display_keys = set(
            zip(
                display_edges["trait_id"],
                display_edges["condition"],
                display_edges["body_site_group"],
                strict=True,
            )
        )
        edges["displayed_in_shared_high_support_map"] = [
            (trait_id, condition, body_site) in display_keys
            for trait_id, condition, body_site in zip(
                edges["trait_id"],
                edges["condition"],
                edges["body_site_group"],
                strict=True,
            )
        ]

        edges_path = OUTPUT_DIR / f"{stem}_edges.tsv"
        display_edges_path = OUTPUT_DIR / f"{stem}_display_edges.tsv"
        pair_selection_path = OUTPUT_DIR / f"{stem}_selected_pairs.tsv"
        png_path = OUTPUT_DIR / f"{stem}.png"

        edges.to_csv(edges_path, sep="\t", index=False)
        display_edges.to_csv(display_edges_path, sep="\t", index=False)
        pair_selection.to_csv(pair_selection_path, sep="\t", index=False)
        plot_combined_direction_sankey(
            display_edges,
            focus_direction,
            png_path,
            title_prefix="Six shared traits shared/high-support",
            subtitle_extra=(
                f"Reduced overlay: top {len(pair_selection)} condition/body-site "
                "pairs ranked by trait recurrence and profile support."
            ),
            footer_extra=(
                f"Pair selection starts from links shared by at least "
                f"{SHARED_PATTERN_MIN_TRAITS} traits or with at least "
                f"{SHARED_PATTERN_HIGH_SUPPORT_MIN_PROFILES} positive-z profiles, "
                f"then keeps the top {SHARED_PATTERN_MAX_PAIRS} pairs."
            ),
            color_scale_vmax=color_vmax,
        )

        output_records.append(
            {
                "analysis": analysis_key,
                "figure_type": "combined_shared_high_support",
                "trait_id": "all",
                "trait_label": "all six shared traits",
                "trait_slug": "all",
                "focus_direction": focus_label,
                "trait_count": display_edges["trait_id"].nunique(),
                "edge_count": len(edges),
                "display_edge_count": len(display_edges),
                "selected_pair_count": len(pair_selection),
                "png_path": png_path.relative_to(REPO_ROOT).as_posix(),
                "edges_path": edges_path.relative_to(REPO_ROOT).as_posix(),
                "display_edges_path": display_edges_path.relative_to(REPO_ROOT).as_posix(),
                "selected_pairs_path": pair_selection_path.relative_to(REPO_ROOT).as_posix(),
                "condition_rollups_path": "",
                "body_site_rollups_path": "",
            }
        )
    return pd.DataFrame.from_records(output_records)


def write_bidirectional_high_support_outputs(
    all_edges: pd.DataFrame,
    all_display_edges: pd.DataFrame,
) -> pd.DataFrame:
    """Overlay the existing per-direction selections without merging their counts."""
    selected_frames = []
    pair_frames = []
    for direction in DIRECTION_LABELS.values():
        selected, pairs = select_shared_high_support_display_edges(
            all_display_edges[all_display_edges["focus_direction"] == direction]
        )
        selected_frames.append(selected)
        pair_frames.append(pairs.assign(focus_direction=direction))
    display_edges = pd.concat(selected_frames, ignore_index=True)
    pair_selection = pd.concat(pair_frames, ignore_index=True)
    direction_keys = ["condition", "body_site_group", "focus_direction"]
    edge_keys = [*direction_keys, "trait_id"]
    edges = all_edges.merge(pair_selection, on=direction_keys, how="inner", validate="many_to_one")
    edges = edges.merge(
        display_edges[edge_keys].assign(displayed_in_shared_high_support_map=True),
        on=edge_keys,
        how="left",
        validate="one_to_one",
    )
    edges["displayed_in_shared_high_support_map"] = edges[
        "displayed_in_shared_high_support_map"
    ].eq(True)
    for frame in [edges, display_edges]:
        frame["line_style"] = frame["focus_direction"].map(
            {"increased": "solid", "decreased": "dashed"}
        )

    png_path = OUTPUT_DIR / f"{BIDIRECTIONAL_STEM}.png"
    edges_path = OUTPUT_DIR / f"{BIDIRECTIONAL_STEM}_edges.tsv"
    display_path = OUTPUT_DIR / f"{BIDIRECTIONAL_STEM}_display_edges.tsv"
    pairs_path = OUTPUT_DIR / f"{BIDIRECTIONAL_STEM}_selected_pairs.tsv"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    edges.to_csv(edges_path, sep="\t", index=False)
    display_edges.to_csv(display_path, sep="\t", index=False)
    pair_selection.to_csv(pairs_path, sep="\t", index=False)
    plot_combined_direction_sankey(
        display_edges,
        "both",
        png_path,
        title_prefix="Six shared traits, high support",
        subplot_top=0.865,
        color_scale_vmax=max(0.4, float(all_display_edges["focus_enriched_mean_guild_z"].max())),
        footer_extra=(
            f"Top {SHARED_PATTERN_MAX_PAIRS} pairs per direction, among pairs shared by "
            f"at least {SHARED_PATTERN_MIN_TRAITS} traits or with at least "
            f"{SHARED_PATTERN_HIGH_SUPPORT_MIN_PROFILES} positive-z profiles in one trait. "
            "Background counts remain in the TSVs; grey underlays are omitted."
        ),
    )
    return pd.DataFrame.from_records(
        [
            {
                "analysis": BIDIRECTIONAL_ANALYSIS,
                "figure_type": "bidirectional_shared_high_support",
                "trait_id": "all",
                "trait_label": "all six shared traits",
                "trait_slug": "all",
                "focus_direction": "both",
                "trait_count": display_edges["trait_id"].nunique(),
                "edge_count": len(edges),
                "display_edge_count": len(display_edges),
                "selected_pair_count": len(pair_selection[direction_keys[:2]].drop_duplicates()),
                "png_path": png_path.relative_to(REPO_ROOT).as_posix(),
                "edges_path": edges_path.relative_to(REPO_ROOT).as_posix(),
                "display_edges_path": display_path.relative_to(REPO_ROOT).as_posix(),
                "selected_pairs_path": pairs_path.relative_to(REPO_ROOT).as_posix(),
                "condition_rollups_path": "",
                "body_site_rollups_path": "",
            }
        ]
    )


def write_bidirectional_readme_note() -> None:
    """Update only this variant's note, preserving the other generated analyses."""
    path = OUTPUT_DIR / "README.md"
    text = path.read_text(encoding="utf-8") if path.exists() else "# Six Shared Traits Sankeys\n"
    section = (
        f"{BIDIRECTIONAL_README_HEADING}\n\n"
        f"![Combined high-support Sankey]({BIDIRECTIONAL_STEM}.png)\n\n"
        "Solid colored curves show increased profiles and dashed colored curves show "
        "decreased profiles. Both run from condition to body site; the line style "
        "describes the profile direction, not a negative trait z-score. A colored "
        "curve represents positive-z support for one trait in one direction. "
        "Hue identifies the trait, width encodes its positive-z profile count, and "
        "shade encodes its mean positive z-score on a shared numeric scale.\n\n"
        f"The top {SHARED_PATTERN_MAX_PAIRS} pairs are selected independently for each "
        "direction using the same rules as the separate high-support plots. Their "
        "union is displayed without combining increased and decreased counts. "
        "A missing direction means that link was not selected for display, not "
        "necessarily that its support is zero. Node n sums profile-trait memberships "
        "over both directions. Grey background underlays are omitted for readability; "
        "background counts remain in the tables.\n\n"
        f"- Displayed links and direction: `{BIDIRECTIONAL_STEM}_display_edges.tsv`\n"
        f"- All edges for selected direction/pair keys: `{BIDIRECTIONAL_STEM}_edges.tsv`\n"
        f"- Pair ranks within each direction: `{BIDIRECTIONAL_STEM}_selected_pairs.tsv`\n\n"
        "Regenerate this overlay from the existing grouped Sankey edge tables with "
        "`just six-shared-trait-combined-sankey`.\n"
    )
    before, found, after = text.partition(BIDIRECTIONAL_README_HEADING)
    if found:
        _, next_heading, remaining = after.partition("\n## ")
        text = before.rstrip() + "\n\n" + section
        if next_heading:
            text += "\n## " + remaining
    else:
        text = text.rstrip() + "\n\n" + section
    path.write_text(text, encoding="utf-8")


def regenerate_bidirectional_high_support() -> None:
    """Regenerate this overlay from saved tables without resetting the output folder."""
    frames = {}
    for table in ["edges", "display_edges"]:
        frames[table] = pd.concat(
            [
                read_tsv(OUTPUT_DIR / f"six_shared_traits_{direction}_profile_sankey_{table}.tsv")
                for direction in DIRECTION_LABELS.values()
            ],
            ignore_index=True,
        )
    new_manifest = write_bidirectional_high_support_outputs(
        frames["edges"], frames["display_edges"]
    )
    manifest_path = OUTPUT_DIR / "six_shared_trait_sankey_manifest.tsv"
    manifest = read_tsv(manifest_path)
    manifest = manifest[manifest["analysis"] != BIDIRECTIONAL_ANALYSIS]
    pd.concat([manifest, new_manifest], ignore_index=True).to_csv(
        manifest_path, sep="\t", index=False
    )
    write_bidirectional_readme_note()
    print(f"Wrote {OUTPUT_DIR / (BIDIRECTIONAL_STEM + '.png')}")
    print(new_manifest[["display_edge_count", "selected_pair_count"]].to_string(index=False))


def write_individual_condition_outputs(
    all_edges: pd.DataFrame,
    all_display_edges: pd.DataFrame,
    condition_rollups: pd.DataFrame,
    body_site_rollups: pd.DataFrame,
) -> pd.DataFrame:
    """Write two combined Sankeys with individual conditions and grouped body sites."""
    output_records: list[dict[str, object]] = []
    color_vmax = shared_color_vmax(all_display_edges)
    for focus_direction, focus_label in DIRECTION_LABELS.items():
        analysis_key = f"six_shared_traits_{focus_label}_individual_conditions"
        stem = f"{analysis_key}_profile_sankey"
        edges = all_edges[all_edges["focus_direction"] == focus_label].copy()
        display_edges = all_display_edges[
            all_display_edges["focus_direction"] == focus_label
        ].copy()
        condition_subset = condition_rollups[
            condition_rollups["focus_direction"] == focus_label
        ].copy()
        body_site_subset = body_site_rollups[
            body_site_rollups["focus_direction"] == focus_label
        ].copy()

        edges_path = OUTPUT_DIR / f"{stem}_edges.tsv"
        display_edges_path = OUTPUT_DIR / f"{stem}_display_edges.tsv"
        condition_rollups_path = OUTPUT_DIR / f"{stem}_condition_rollups.tsv"
        body_site_rollups_path = OUTPUT_DIR / f"{stem}_body_site_rollups.tsv"
        png_path = OUTPUT_DIR / f"{stem}.png"

        edges.to_csv(edges_path, sep="\t", index=False)
        display_edges.to_csv(display_edges_path, sep="\t", index=False)
        condition_subset.to_csv(condition_rollups_path, sep="\t", index=False)
        body_site_subset.to_csv(body_site_rollups_path, sep="\t", index=False)
        plot_combined_direction_sankey(
            display_edges,
            focus_direction,
            png_path,
            title_prefix="Six shared traits, individual conditions",
            subtitle_extra=(
                "Conditions use raw individual labels; body sites keep grouped v3 labels."
            ),
            footer_extra=(
                "Hybrid display mode keeps individual condition labels on the "
                "left and grouped v3 body-site labels on the right. Same-trait "
                "profiles from one condition into the same body-site group are "
                "summed into one colored link."
            ),
            condition_axis_label="Individual conditions",
            body_site_axis_label="Body-site display groups",
            display_filter_note=(
                "Display filtering uses the existing v3 thresholds after "
                "aggregating by raw condition, grouped body-site label, and trait: "
                f"at least {V3_MIN_FOCUS_ENRICHED_PROFILES} positive-z focus "
                f"profiles and {MIN_STUDIES} studies where possible."
            ),
            color_scale_vmax=color_vmax,
            subplot_top=0.855,
            bundle_by_trait_at_body_site=True,
        )

        output_records.append(
            {
                "analysis": analysis_key,
                "figure_type": "combined_individual_conditions",
                "trait_id": "all",
                "trait_label": "all six shared traits",
                "trait_slug": "all",
                "focus_direction": focus_label,
                "trait_count": display_edges["trait_id"].nunique(),
                "edge_count": len(edges),
                "display_edge_count": len(display_edges),
                "png_path": png_path.relative_to(REPO_ROOT).as_posix(),
                "edges_path": edges_path.relative_to(REPO_ROOT).as_posix(),
                "display_edges_path": display_edges_path.relative_to(REPO_ROOT).as_posix(),
                "condition_rollups_path": condition_rollups_path.relative_to(REPO_ROOT).as_posix(),
                "body_site_rollups_path": body_site_rollups_path.relative_to(REPO_ROOT).as_posix(),
            }
        )
    return pd.DataFrame.from_records(output_records)


def write_individual_outputs(
    all_edges: pd.DataFrame,
    all_display_edges: pd.DataFrame,
    condition_rollups: pd.DataFrame,
    body_site_rollups: pd.DataFrame,
    trait_selection: pd.DataFrame,
) -> pd.DataFrame:
    """Write the twelve individual trait-direction Sankeys and TSVs."""
    output_records: list[dict[str, object]] = []
    for trait in SHARED_TRAITS:
        trait_selection[trait_selection["trait_id"] == trait.trait_id].to_csv(
            OUTPUT_DIR / f"{trait.slug}_trait_selection.tsv", sep="\t", index=False
        )

        for focus_direction, focus_label in DIRECTION_LABELS.items():
            analysis_key = f"{trait.slug}_{focus_label}"
            stem = output_stem(trait, focus_direction)
            edges = all_edges[all_edges["analysis"] == analysis_key].copy()
            display_edges = all_display_edges[all_display_edges["analysis"] == analysis_key].copy()
            condition_subset = condition_rollups[
                condition_rollups["analysis"] == analysis_key
            ].copy()
            body_site_subset = body_site_rollups[
                body_site_rollups["analysis"] == analysis_key
            ].copy()

            edges_path = OUTPUT_DIR / f"{stem}_edges.tsv"
            display_edges_path = OUTPUT_DIR / f"{stem}_display_edges.tsv"
            condition_rollups_path = OUTPUT_DIR / f"{stem}_condition_rollups.tsv"
            body_site_rollups_path = OUTPUT_DIR / f"{stem}_body_site_rollups.tsv"
            png_path = OUTPUT_DIR / f"{stem}.png"

            edges.to_csv(edges_path, sep="\t", index=False)
            display_edges.to_csv(display_edges_path, sep="\t", index=False)
            condition_subset.to_csv(condition_rollups_path, sep="\t", index=False)
            body_site_subset.to_csv(body_site_rollups_path, sep="\t", index=False)
            plot_single_trait_sankey(display_edges, trait, focus_direction, png_path)

            output_records.append(
                {
                    "analysis": analysis_key,
                    "figure_type": "single_trait_direction",
                    "trait_id": trait.trait_id,
                    "trait_label": trait.label,
                    "trait_slug": trait.slug,
                    "focus_direction": focus_label,
                    "trait_count": 1,
                    "edge_count": len(edges),
                    "display_edge_count": len(display_edges),
                    "png_path": png_path.relative_to(REPO_ROOT).as_posix(),
                    "edges_path": edges_path.relative_to(REPO_ROOT).as_posix(),
                    "display_edges_path": display_edges_path.relative_to(REPO_ROOT).as_posix(),
                    "condition_rollups_path": condition_rollups_path.relative_to(
                        REPO_ROOT
                    ).as_posix(),
                    "body_site_rollups_path": body_site_rollups_path.relative_to(
                        REPO_ROOT
                    ).as_posix(),
                }
            )
    return pd.DataFrame.from_records(output_records)


def write_readme(summary: pd.DataFrame, output_manifest: pd.DataFrame) -> None:
    """Write a concise methods note and file index for the generated outputs."""
    lines = [
        "# Six Shared Traits Sankeys",
        "",
        "These outputs render the six pooled paired-contrast traits with "
        "`q_value <= 0.001` as two combined condition-to-body-site Sankey maps "
        "plus readability variants and twelve individual trait-direction Sankey "
        "maps.",
        "",
        "## Inputs",
        "",
        f"- Profile counts: `{PROFILE_MATRIX_PATH.relative_to(REPO_ROOT).as_posix()}`",
        f"- Human profile metadata: `{PROFILE_METADATA_PATH.relative_to(REPO_ROOT).as_posix()}`",
        f"- Paired trait statistics: `{TRAIT_STATS_PATH.relative_to(REPO_ROOT).as_posix()}`",
        "",
        "## Scoring",
        "",
        "For each trait, raw KG-Microbe trait counts were divided by the profile's "
        "total trait signal, scaled per 1,000 counts, log1p-transformed, and "
        "standardized to z-scores across the human-host profile universe used by "
        "the existing study trait profile maps. Each trait is then treated as a "
        "one-trait guild.",
        "",
        "## Combined Readability Variants",
        "",
        "- `combined_direction`: all-edge overlay for one profile direction, "
        "with same-color links bundled into one summed trunk per body site.",
        "- `combined_trait_lanes`: option 1; one figure per direction with six stacked "
        "trait lanes, preserving all displayed trait-specific links while separating "
        "line crossings by trait.",
        "- `combined_shared_high_support`: option 2; one reduced overlay per direction, "
        f"limited to the top {SHARED_PATTERN_MAX_PAIRS} condition/body-site pairs "
        f"selected from pairs shared by at least {SHARED_PATTERN_MIN_TRAITS} traits "
        f"or with at least {SHARED_PATTERN_HIGH_SUPPORT_MIN_PROFILES} positive-z "
        "focus profiles.",
        "- `combined_individual_conditions`: one figure per direction with raw "
        "individual condition labels on the left and grouped v3 body-site labels "
        "on the right. Same-trait records from one condition into one body-site "
        "group are summed into a single colored link, then bundled into one "
        "body-site trunk per trait color.",
        "- `single_trait_direction`: individual one-trait, one-direction Sankeys.",
        "",
        "The combined maps retain the individual trait-direction edge records. "
        "Every colored curve in the combined maps is a trait-specific "
        "condition/body-site link. The individual maps show one trait and one "
        "profile direction per figure. Curve width is the focus-direction "
        "positive-z profile count for that trait and shade is the mean positive "
        "focus-profile z-score. The overlay maps use the same numeric color "
        "range for all six trait-specific monochrome scales. Light grey underlays "
        "show the z <= 0 background for the same condition/body-site/trait link. "
        "Parallel curves sharing a condition/body-site pair are ordered by edge "
        "width from thickest to thinnest. In `combined_direction` and "
        "`combined_individual_conditions`, links approaching each body-site node "
        "merge into at most six colored trunks, one per trait color scale.",
        "",
        "Image node keys spell out the compact labels used in the PNGs: `n=` is "
        "positive-z focus-profile support, `s=` is the unique study count for "
        "`n`, `bg=` is z <= 0 background-profile support, `links=` is plotted "
        "trait-specific link count, and `traits=` is the number of distinct "
        "traits represented at a combined-map node.",
        "",
        "Node counts are profile-trait memberships across displayed links, not "
        "deduplicated profile totals. Use the edge TSVs for exact trait-specific "
        "counts, study IDs, raw rolled-up labels, and increased-minus-decreased "
        "contrasts.",
        "",
        "## Traits",
        "",
        "| trait | ID | color | q-value | paired direction | profile count | positive-z profiles |",
        "| --- | --- | --- | ---: | --- | ---: | ---: |",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.trait_label} | `{row.trait_id}` | `{row.base_color}` | "
            f"{row.paired_q_value:.3g} | {row.paired_direction} | "
            f"{row.profile_count} | {row.positive_z_profile_count} |"
        )

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "| analysis | type | PNG | displayed trait links | selected pairs | all edge rows |",
            "| --- | --- | --- | ---: | ---: | ---: |",
        ]
    )
    for row in output_manifest.itertuples(index=False):
        selected_pairs = getattr(row, "selected_pair_count", pd.NA)
        selected_pairs_text = "" if pd.isna(selected_pairs) else str(int(float(selected_pairs)))
        lines.append(
            f"| {row.analysis} | {row.figure_type} | `{row.png_path}` | "
            f"{row.display_edge_count} | {selected_pairs_text} | {row.edge_count} |"
        )
    lines.append("")
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    reset_output_dir()
    metadata = load_profile_metadata()
    raw_trait_counts = load_raw_trait_counts(SHARED_TRAITS, metadata)
    zscores = build_zscore_matrix(raw_trait_counts, metadata)
    q_stats = load_q_significance()

    summary = build_trait_summary(metadata, raw_trait_counts, zscores, q_stats)
    trait_selection = build_trait_selection(q_stats)
    summary.to_csv(OUTPUT_DIR / "six_shared_trait_summary.tsv", sep="\t", index=False)
    trait_selection.to_csv(
        OUTPUT_DIR / "six_shared_trait_selection.tsv",
        sep="\t",
        index=False,
    )

    all_edges, display_edges, condition_rollups, body_site_rollups = build_trait_direction_tables(
        metadata, zscores
    )
    (
        individual_condition_edges,
        individual_condition_display_edges,
        individual_condition_rollups,
        individual_condition_body_site_rollups,
    ) = build_individual_condition_direction_tables(metadata, zscores)
    combined_manifest = write_direction_outputs(
        all_edges,
        display_edges,
        condition_rollups,
        body_site_rollups,
    )
    trait_lane_manifest = write_trait_lane_outputs(all_edges, display_edges)
    shared_high_support_manifest = write_shared_high_support_outputs(all_edges, display_edges)
    bidirectional_manifest = write_bidirectional_high_support_outputs(all_edges, display_edges)
    individual_condition_manifest = write_individual_condition_outputs(
        individual_condition_edges,
        individual_condition_display_edges,
        individual_condition_rollups,
        individual_condition_body_site_rollups,
    )
    individual_manifest = write_individual_outputs(
        all_edges,
        display_edges,
        condition_rollups,
        body_site_rollups,
        trait_selection,
    )
    output_manifest = pd.concat(
        [
            combined_manifest,
            trait_lane_manifest,
            shared_high_support_manifest,
            bidirectional_manifest,
            individual_condition_manifest,
            individual_manifest,
        ],
        ignore_index=True,
    )
    output_manifest.to_csv(
        OUTPUT_DIR / "six_shared_trait_sankey_manifest.tsv",
        sep="\t",
        index=False,
    )
    write_readme(summary, output_manifest)
    write_bidirectional_readme_note()
    graph_manifest = write_six_shared_trait_graphs(OUTPUT_DIR, REPO_ROOT)
    static_graph_manifest = write_six_shared_trait_static_graphs(OUTPUT_DIR)
    cytoscape_heatmap_manifest = write_six_shared_trait_cytoscape_heatmaps(
        OUTPUT_DIR,
        REPO_ROOT,
    )

    print(f"Wrote output directory: {OUTPUT_DIR}")
    print(f"Wrote PNGs: {len(output_manifest)}")
    print(f"Wrote manifest: {OUTPUT_DIR / 'six_shared_trait_sankey_manifest.tsv'}")
    print(f"Wrote graph explorer: {OUTPUT_DIR / HTML_FILENAME}")
    print(f"Wrote graph views: {len(graph_manifest)}")
    print(f"Wrote static graph PNGs: {len(static_graph_manifest)}")
    print(f"Wrote Cytoscape heatmaps: {len(cytoscape_heatmap_manifest)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--combined-high-support-only",
        action="store_true",
        help="Combine directions from saved grouped tables without regenerating other figures.",
    )
    args = parser.parse_args()
    if args.combined_high_support_only:
        regenerate_bidirectional_high_support()
    else:
        main()
