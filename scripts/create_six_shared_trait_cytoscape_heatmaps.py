"""Create inferential heatmaps from the six-shared-trait Cytoscape inputs."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import fisher_exact

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization" / "trait-profile-investigation" / "six-shared-traits"

GRAPH_MANIFEST_FILENAME = "six_shared_trait_graph_manifest.tsv"
TRAIT_SELECTION_FILENAME = "six_shared_trait_selection.tsv"
MANIFEST_FILENAME = "six_shared_trait_cytoscape_heatmap_manifest.tsv"
PATTERN_TABLE_FILENAME = "six_shared_trait_condition_body_site_trait_pattern_matrix.tsv"
CONTRAST_TABLE_FILENAME = "six_shared_trait_condition_body_site_direction_contrasts.tsv"
SELECTED_CELLS_FILENAME = "six_shared_trait_condition_body_site_bh_selected_cells.tsv"
CONTRAST_PAGES_FILENAME = "six_shared_trait_signed_direction_heatmap_pages.tsv"
GROUPED_CONTRAST_CLUSTER_ROWS_FILENAME = (
    "six_shared_traits_grouped_signed_direction_cluster_rows.tsv"
)
GROUPED_CONTRAST_CLUSTER_SUMMARY_FILENAME = (
    "six_shared_traits_grouped_signed_direction_cluster_summary.tsv"
)
CONTRAST_ROWS_PER_PAGE = 24
GROUPED_CONTRAST_CLUSTER_COUNT = 6

MIN_POSITIVE_PROFILES = 3
MIN_STUDIES = 2
P_VALUE_CUTOFF = 0.01
Q_VALUE_CUTOFF = 0.05
MAX_MINUS_LOG10_P = 10.0
MAX_PATTERN_ROWS = 24

SLATE_50 = "#f8fafc"
SLATE_100 = "#f1f5f9"
SLATE_200 = "#e2e8f0"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"

CLUSTER_LABEL_REPLACEMENTS = {
    "Oral cavity / dental": "Oral/dental",
    "Sampling site / time protocol": "Sampling/time",
    "Vaginal / reproductive": "Vaginal",
    "Cancer / neoplasia": "Cancer",
    "Smoking / alcohol / drug use": "Smoking/drug",
    "Renal / urinary disease": "Renal/urinary",
    "Urinary tract": "Urinary",
    "Colorectal neoplasia": "Colorectal",
}

LEGACY_CORRELATION_OUTPUTS = (
    "six_shared_trait_condition_body_site_trait_correlations.tsv",
    "six_shared_trait_condition_body_site_trait_vectors.tsv",
    "six_shared_trait_correlation_heatmap_manifest.tsv",
    "six_shared_traits_increased_grouped_condition_body_site_trait_correlation_heatmap.png",
    "six_shared_traits_decreased_grouped_condition_body_site_trait_correlation_heatmap.png",
    (
        "six_shared_traits_increased_individual_conditions_condition_body_site_trait_"
        "correlation_heatmap.png"
    ),
    (
        "six_shared_traits_decreased_individual_conditions_condition_body_site_trait_"
        "correlation_heatmap.png"
    ),
)


@dataclass(frozen=True)
class ContrastMode:
    """One condition-label mode for signed increased-vs-decreased contrasts."""

    key: str
    label: str
    condition_mode: str
    increased_view_key: str
    decreased_view_key: str


CONTRAST_MODES = (
    ContrastMode(
        key="grouped",
        label="Grouped conditions",
        condition_mode="grouped condition labels",
        increased_view_key="increased_grouped",
        decreased_view_key="decreased_grouped",
    ),
    ContrastMode(
        key="individual_conditions",
        label="Individual conditions",
        condition_mode="raw individual condition labels",
        increased_view_key="increased_individual_conditions",
        decreased_view_key="decreased_individual_conditions",
    ),
)


def read_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV with local defaults."""
    return pd.read_csv(path, sep="\t")


def write_tsv(table: pd.DataFrame, path: Path) -> None:
    """Write TSVs with explicit missing values instead of trailing empty fields."""
    table.replace("", pd.NA).to_csv(path, sep="\t", index=False, na_rep="NA")


def relative_path(path: Path, repo_root: Path = REPO_ROOT) -> str:
    """Return a stable repo-relative path."""
    return path.relative_to(repo_root).as_posix()


def resolve_repo_path(path_text: str, repo_root: Path = REPO_ROOT) -> Path:
    """Resolve a manifest path that may be repo-relative or absolute."""
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def full_edges_path(display_path: str, repo_root: Path = REPO_ROOT) -> Path:
    """Resolve the complete edge table paired with a graph display table."""
    path = resolve_repo_path(display_path, repo_root)
    suffix = "_display_edges.tsv"
    if not path.name.endswith(suffix):
        raise ValueError(f"Expected a {suffix} source, got {path.name}")
    return path.with_name(path.name.removesuffix(suffix) + "_edges.tsv")


def wrap_label(value: object, width: int = 24) -> str:
    """Wrap labels for compact heatmap axes."""
    text = str(value)
    return "\n".join(wrap(text, width=width, break_long_words=False)) or text


def compact_cluster_label(value: object, *, width: int = 18) -> str:
    """Return a short anatomy/condition label for dendrogram cluster spans."""
    text = str(value)
    if ":" in text and text.split(":", 1)[0].strip().startswith("C"):
        text = text.split(":", 1)[1].strip()
    parts = [part.strip() for part in text.split(";") if part.strip()]
    compact_parts = [
        CLUSTER_LABEL_REPLACEMENTS.get(part, part).replace(" / ", "/") for part in parts[:2]
    ]
    return "\n".join(wrap_label(part, width=width) for part in compact_parts) or wrap_label(
        text,
        width=width,
    )


def format_p_value(value: float) -> str:
    """Format p-values for cell annotations."""
    if not math.isfinite(value):
        return "NA"
    if value < 0.001:
        return "<0.001"
    if value < P_VALUE_CUTOFF <= round(value, 3):
        return f"<{P_VALUE_CUTOFF:g}"
    return f"{value:.3f}"


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Calculate Benjamini-Hochberg adjusted q-values."""
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    adjusted = np.ones_like(values, dtype=float)
    valid = np.isfinite(values)
    valid_indices = np.flatnonzero(valid)
    if len(valid_indices) == 0:
        return pd.Series(adjusted, index=p_values.index)

    sorted_indices = valid_indices[np.argsort(values[valid])]
    sorted_p = values[sorted_indices]
    ranks = np.arange(1, len(sorted_p) + 1, dtype=float)
    sorted_adjusted = np.minimum.accumulate((sorted_p * len(sorted_p) / ranks)[::-1])[::-1]
    adjusted[sorted_indices] = np.clip(sorted_adjusted, 0.0, 1.0)
    return pd.Series(adjusted, index=p_values.index)


def minus_log10_p_values(p_values: pd.Series) -> pd.Series:
    """Convert finite p-values to positive evidence scores for heatmap coloring."""
    values = pd.to_numeric(p_values, errors="coerce")
    finite = np.isfinite(values)
    output = pd.Series(np.nan, index=p_values.index, dtype=float)
    clipped = values[finite].clip(lower=np.finfo(float).tiny, upper=1.0)
    output.loc[finite] = -np.log10(clipped)
    return output


def finalize_inferential_cells(cells: pd.DataFrame) -> pd.DataFrame:
    """Add p/q-value decisions and heatmap color values to tested cells."""
    output = cells.copy()
    p_values = pd.to_numeric(output["fisher_p_value"], errors="coerce")
    finite_p = np.isfinite(p_values)
    support = output["passes_support"].astype(bool)

    output["fisher_q_value"] = benjamini_hochberg(p_values.fillna(1.0))
    q_values = pd.to_numeric(output["fisher_q_value"], errors="coerce")
    finite_q = np.isfinite(q_values)
    output["minus_log10_fisher_p_value"] = minus_log10_p_values(p_values)
    output["plot_minus_log10_fisher_p_value"] = output["minus_log10_fisher_p_value"].clip(
        upper=MAX_MINUS_LOG10_P
    )
    output["passes_p_value_cutoff"] = finite_p & (p_values < P_VALUE_CUTOFF)
    output["displayed_in_heatmap"] = support & output["passes_p_value_cutoff"]
    output["significant_in_heatmap"] = output["displayed_in_heatmap"]
    output["selected_after_bh"] = support & finite_q & (q_values <= Q_VALUE_CUTOFF)

    output["blank_reason"] = ""
    output.loc[~support, "blank_reason"] = "low positive-z profile or study support"
    output.loc[support & ~finite_p, "blank_reason"] = "Fisher p=NA"
    output.loc[
        support & finite_p & ~output["passes_p_value_cutoff"],
        "blank_reason",
    ] = f"Fisher p>={P_VALUE_CUTOFF:g}"

    if "increased_minus_decreased_positive_fraction" in output:
        delta = pd.to_numeric(
            output["increased_minus_decreased_positive_fraction"],
            errors="coerce",
        ).fillna(0.0)
        direction = np.sign(delta)
        output["signed_minus_log10_fisher_p_value"] = (
            direction * output["minus_log10_fisher_p_value"]
        )
        output["plot_signed_minus_log10_fisher_p_value"] = (
            direction * output["plot_minus_log10_fisher_p_value"]
        )
        output["displayed_in_heatmap"] = support & finite_p
        output.loc[output["displayed_in_heatmap"], "blank_reason"] = ""

    return output


def clean_legacy_outputs(output_dir: Path = OUTPUT_DIR) -> None:
    """Remove the superseded correlation heatmap outputs from the prior workflow."""
    for filename in LEGACY_CORRELATION_OUTPUTS:
        path = output_dir / filename
        if path.exists():
            path.unlink()


def read_trait_order(output_dir: Path = OUTPUT_DIR) -> pd.DataFrame:
    """Read the six shared traits in display order."""
    traits = read_tsv(output_dir / TRAIT_SELECTION_FILENAME)
    traits = traits.rename(columns={"trait_name": "trait_label"})
    traits["trait_rank"] = np.arange(len(traits), dtype=int)
    return traits[["trait_id", "trait_label", "slug", "base_color", "trait_rank"]].copy()


def coerce_display_edge_numbers(edges: pd.DataFrame) -> pd.DataFrame:
    """Coerce display-edge numeric columns used by both heatmap types."""
    output = edges.copy()
    numeric_columns = [
        "trait_rank",
        "profile_count",
        "study_count",
        "focus_profile_count",
        "focus_study_count",
        "focus_enriched_profile_count",
        "focus_enriched_study_count",
        "focus_no_enriched_trait_profile_count",
        "comparator_profile_count",
        "comparator_study_count",
        "comparator_enriched_profile_count",
        "background_no_enriched_profile_count",
        "focus_enriched_mean_guild_z",
        "focus_enriched_minus_comparator_guild_z",
        "rank_score",
    ]
    for column in numeric_columns:
        if column in output:
            output[column] = pd.to_numeric(output[column], errors="coerce").fillna(0.0)
    return output


def signal_value(row: pd.Series) -> float:
    """Return the nonnegative signal used by the Cytoscape color scale."""
    contrast = row.get("focus_enriched_minus_comparator_guild_z", math.nan)
    mean = row.get("focus_enriched_mean_guild_z", math.nan)
    value = contrast if pd.notna(contrast) else mean
    if pd.isna(value):
        return 0.0
    return max(0.0, float(value))


def focus_vs_comparator_p_value(row: pd.Series) -> float:
    """Test whether focus positive-z support exceeds comparator positive-z support."""
    focus_positive = int(row["focus_enriched_profile_count"])
    focus_negative = int(row["focus_no_enriched_trait_profile_count"])
    comparator_positive = int(row["comparator_enriched_profile_count"])
    comparator_negative = max(0, int(row["comparator_profile_count"]) - comparator_positive)
    if focus_positive + focus_negative == 0 or comparator_positive + comparator_negative == 0:
        return math.nan
    result = fisher_exact(
        [[focus_positive, focus_negative], [comparator_positive, comparator_negative]],
        alternative="greater",
    )
    return float(result.pvalue)


def signed_direction_p_value(
    increased_positive: int,
    increased_total: int,
    decreased_positive: int,
    decreased_total: int,
) -> float:
    """Test whether positive-z support differs between increased and decreased profiles."""
    if increased_total <= 0 or decreased_total <= 0:
        return math.nan
    result = fisher_exact(
        [
            [increased_positive, max(0, increased_total - increased_positive)],
            [decreased_positive, max(0, decreased_total - decreased_positive)],
        ],
        alternative="two-sided",
    )
    return float(result.pvalue)


def pair_label(condition: object, body_site_group: object) -> str:
    """Return the row label used by condition/body-site heatmaps."""
    return f"{condition} | {body_site_group}"


def prepare_pattern_cells(
    display_edges: pd.DataFrame,
    view: pd.Series,
    trait_order: pd.DataFrame,
) -> pd.DataFrame:
    """Build long-form cells for one condition/body-site x trait pattern heatmap."""
    trait_rank = dict(zip(trait_order["trait_id"], trait_order["trait_rank"], strict=True))
    edges = coerce_display_edge_numbers(display_edges)
    records: list[dict[str, object]] = []
    for row in edges.itertuples(index=False):
        row_series = pd.Series(row._asdict())
        p_value = focus_vs_comparator_p_value(row_series)
        signal = signal_value(row_series)
        strength = math.log1p(float(row.focus_enriched_profile_count)) * signal
        supported = (
            int(row.focus_enriched_profile_count) >= MIN_POSITIVE_PROFILES
            and int(row.focus_enriched_study_count) >= MIN_STUDIES
        )
        records.append(
            {
                "view_key": view.view_key,
                "view_label": view.view_label,
                "focus_direction": view.focus_direction,
                "condition_mode": view.condition_mode,
                "condition": row.condition,
                "body_site_group": row.body_site_group,
                "condition_body_site_pair": pair_label(row.condition, row.body_site_group),
                "trait_id": row.trait_id,
                "trait_label": row.trait_label,
                "trait_rank": trait_rank.get(row.trait_id, 999),
                "focus_positive_profile_count": int(row.focus_enriched_profile_count),
                "focus_positive_study_count": int(row.focus_enriched_study_count),
                "focus_profile_count": int(row.focus_profile_count),
                "comparator_positive_profile_count": int(row.comparator_enriched_profile_count),
                "comparator_profile_count": int(row.comparator_profile_count),
                "signal": signal,
                "strength": strength,
                "fisher_p_value": p_value,
                "passes_support": supported,
            }
        )
    return finalize_inferential_cells(pd.DataFrame.from_records(records))


def select_pattern_rows(cells: pd.DataFrame, *, max_rows: int = MAX_PATTERN_ROWS) -> list[str]:
    """Select condition/body-site rows with the strongest significant trait recurrence."""
    cells = cells.copy()
    cells["display_evidence"] = cells["plot_minus_log10_fisher_p_value"].where(
        cells["displayed_in_heatmap"],
        0.0,
    )
    cells["display_positive_profiles"] = cells["focus_positive_profile_count"].where(
        cells["displayed_in_heatmap"],
        0,
    )
    cells["display_strength"] = cells["strength"].where(cells["displayed_in_heatmap"], 0.0)
    stats = (
        cells.groupby("condition_body_site_pair", sort=False)
        .agg(
            significant_trait_count=("displayed_in_heatmap", "sum"),
            bh_selected_trait_count=("selected_after_bh", "sum"),
            total_positive_profiles=("display_positive_profiles", "sum"),
            max_evidence=("display_evidence", "max"),
            total_evidence=("display_evidence", "sum"),
            total_strength=("display_strength", "sum"),
        )
        .sort_values(
            [
                "significant_trait_count",
                "bh_selected_trait_count",
                "total_evidence",
                "total_positive_profiles",
                "total_strength",
                "max_evidence",
            ],
            ascending=[False, False, False, False, False, False],
        )
    )
    stats = stats[stats["significant_trait_count"] > 0]
    if stats.empty:
        stats = (
            cells.groupby("condition_body_site_pair", sort=False)
            .agg(
                significant_trait_count=("displayed_in_heatmap", "sum"),
                total_positive_profiles=("focus_positive_profile_count", "sum"),
                max_evidence=("minus_log10_fisher_p_value", "max"),
                total_evidence=("minus_log10_fisher_p_value", "sum"),
                max_strength=("strength", "max"),
                total_strength=("strength", "sum"),
            )
            .sort_values(
                ["total_evidence", "total_positive_profiles", "total_strength", "max_strength"],
                ascending=[False, False, False, False],
            )
        )
    return stats.head(max_rows).index.astype(str).tolist()


def prepare_contrast_cells(
    increased_edges: pd.DataFrame,
    decreased_edges: pd.DataFrame,
    mode: ContrastMode,
    trait_order: pd.DataFrame,
) -> pd.DataFrame:
    """Build signed increased-vs-decreased cells for one condition-label mode."""
    trait_rank = dict(zip(trait_order["trait_id"], trait_order["trait_rank"], strict=True))
    records: list[dict[str, object]] = []
    for edges in (increased_edges, decreased_edges):
        edges = coerce_display_edge_numbers(edges)
        for row in edges.itertuples(index=False):
            if row.focus_direction == "increased":
                increased_positive = int(row.focus_enriched_profile_count)
                increased_total = int(row.focus_profile_count)
                decreased_positive = int(row.comparator_enriched_profile_count)
                decreased_total = int(row.comparator_profile_count)
            else:
                decreased_positive = int(row.focus_enriched_profile_count)
                decreased_total = int(row.focus_profile_count)
                increased_positive = int(row.comparator_enriched_profile_count)
                increased_total = int(row.comparator_profile_count)
            records.append(
                {
                    "condition": row.condition,
                    "body_site_group": row.body_site_group,
                    "trait_id": row.trait_id,
                    "trait_label": row.trait_label,
                    "trait_rank": trait_rank.get(row.trait_id, 999),
                    "study_count": int(row.study_count),
                    "increased_positive_profile_count": increased_positive,
                    "increased_profile_count": increased_total,
                    "decreased_positive_profile_count": decreased_positive,
                    "decreased_profile_count": decreased_total,
                }
            )
    combined = (
        pd.DataFrame.from_records(records)
        .groupby(["condition", "body_site_group", "trait_id", "trait_label"], as_index=False)
        .agg(
            trait_rank=("trait_rank", "min"),
            study_count=("study_count", "max"),
            increased_positive_profile_count=("increased_positive_profile_count", "max"),
            increased_profile_count=("increased_profile_count", "max"),
            decreased_positive_profile_count=("decreased_positive_profile_count", "max"),
            decreased_profile_count=("decreased_profile_count", "max"),
        )
    )

    output_records: list[dict[str, object]] = []
    for row in combined.itertuples(index=False):
        increased_fraction = (
            row.increased_positive_profile_count / row.increased_profile_count
            if row.increased_profile_count
            else math.nan
        )
        decreased_fraction = (
            row.decreased_positive_profile_count / row.decreased_profile_count
            if row.decreased_profile_count
            else math.nan
        )
        delta = increased_fraction - decreased_fraction
        p_value = signed_direction_p_value(
            int(row.increased_positive_profile_count),
            int(row.increased_profile_count),
            int(row.decreased_positive_profile_count),
            int(row.decreased_profile_count),
        )
        supported = (
            max(
                int(row.increased_positive_profile_count),
                int(row.decreased_positive_profile_count),
            )
            >= MIN_POSITIVE_PROFILES
            and int(row.study_count) >= MIN_STUDIES
        )
        output_records.append(
            {
                "view_key": f"{mode.key}_direction_contrast",
                "view_label": f"{mode.label}: signed increased-vs-decreased contrast",
                "focus_direction": "increased_vs_decreased",
                "condition_mode": mode.condition_mode,
                "condition": row.condition,
                "body_site_group": row.body_site_group,
                "condition_body_site_pair": pair_label(row.condition, row.body_site_group),
                "trait_id": row.trait_id,
                "trait_label": row.trait_label,
                "trait_rank": int(row.trait_rank),
                "study_count": int(row.study_count),
                "increased_positive_profile_count": int(row.increased_positive_profile_count),
                "increased_profile_count": int(row.increased_profile_count),
                "decreased_positive_profile_count": int(row.decreased_positive_profile_count),
                "decreased_profile_count": int(row.decreased_profile_count),
                "increased_positive_fraction": increased_fraction,
                "decreased_positive_fraction": decreased_fraction,
                "increased_minus_decreased_positive_fraction": delta,
                "fisher_p_value": p_value,
                "passes_support": supported,
            }
        )
    return finalize_inferential_cells(pd.DataFrame.from_records(output_records))


def select_contrast_rows(cells: pd.DataFrame) -> list[str]:
    """Group all supported rows by body site, then rank by evidence and coverage."""
    cells = cells.copy()
    cells["abs_delta"] = cells["increased_minus_decreased_positive_fraction"].abs()
    cells["abs_display_evidence"] = (
        cells["plot_signed_minus_log10_fisher_p_value"]
        .abs()
        .where(
            cells["displayed_in_heatmap"],
            0.0,
        )
    )
    cells["positive_support"] = (
        cells["increased_positive_profile_count"] + cells["decreased_positive_profile_count"]
    )
    cells["display_positive_support"] = cells["positive_support"].where(
        cells["displayed_in_heatmap"],
        0,
    )
    cells["display_abs_delta"] = cells["abs_delta"].where(cells["displayed_in_heatmap"], 0.0)
    stats = (
        cells.groupby("condition_body_site_pair", sort=False)
        .agg(
            body_site_group=("body_site_group", "first"),
            displayed_trait_count=("displayed_in_heatmap", "sum"),
            significant_trait_count=("significant_in_heatmap", "sum"),
            bh_selected_trait_count=("selected_after_bh", "sum"),
            positive_support=("display_positive_support", "sum"),
            max_abs_delta=("display_abs_delta", "max"),
            total_abs_delta=("display_abs_delta", "sum"),
            max_evidence=("abs_display_evidence", "max"),
            total_evidence=("abs_display_evidence", "sum"),
        )
        .sort_values(
            [
                "body_site_group",
                "significant_trait_count",
                "displayed_trait_count",
                "bh_selected_trait_count",
                "total_evidence",
                "positive_support",
                "total_abs_delta",
                "max_abs_delta",
            ],
            ascending=[True, False, False, False, False, False, False, False],
        )
    )
    stats = stats[stats["displayed_trait_count"] > 0]
    return stats.index.astype(str).tolist()


def top_values(values: pd.Series, *, max_items: int = 3) -> str:
    """Return a semicolon-delimited mode summary for cluster labels."""
    counts = values.astype(str).value_counts()
    if counts.empty:
        return "NA"
    return "; ".join(f"{label} ({count})" for label, count in counts.head(max_items).items())


def spearman_row_linkage(matrix: pd.DataFrame) -> np.ndarray | None:
    """Cluster rows with average linkage on Spearman-correlation distance."""
    if len(matrix.index) < 2:
        return None
    values = matrix.fillna(0.0).astype(float)
    correlations = values.T.corr(method="spearman").fillna(0.0).clip(-1.0, 1.0)
    correlation_values = correlations.to_numpy(dtype=float, copy=True)
    np.fill_diagonal(correlation_values, 1.0)
    distances = np.clip(1.0 - correlation_values, 0.0, 2.0)
    condensed = squareform(distances, checks=False)
    return linkage(condensed, method="average")


def contrast_row_stats(cells: pd.DataFrame) -> pd.DataFrame:
    """Summarize displayed contrast rows for ordering and cluster audits."""
    frame = cells.copy()
    frame["positive_support"] = (
        frame["increased_positive_profile_count"] + frame["decreased_positive_profile_count"]
    ).where(frame["displayed_in_heatmap"], 0)
    frame["abs_delta"] = (
        frame["increased_minus_decreased_positive_fraction"]
        .abs()
        .where(
            frame["displayed_in_heatmap"],
            0.0,
        )
    )
    frame["abs_evidence"] = (
        frame["plot_signed_minus_log10_fisher_p_value"]
        .abs()
        .where(
            frame["displayed_in_heatmap"],
            0.0,
        )
    )
    return frame.groupby("condition_body_site_pair", sort=False).agg(
        condition=("condition", "first"),
        body_site_group=("body_site_group", "first"),
        displayed_trait_count=("displayed_in_heatmap", "sum"),
        significant_trait_count=("significant_in_heatmap", "sum"),
        bh_selected_trait_count=("selected_after_bh", "sum"),
        positive_support=("positive_support", "sum"),
        total_abs_delta=("abs_delta", "sum"),
        max_abs_delta=("abs_delta", "max"),
        total_evidence=("abs_evidence", "sum"),
        max_evidence=("abs_evidence", "max"),
    )


def clustered_contrast_rows(
    matrix: pd.DataFrame,
    cells: pd.DataFrame,
    *,
    max_clusters: int = GROUPED_CONTRAST_CLUSTER_COUNT,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray | None]:
    """Return clustered row metadata and per-cluster labels for the grouped contrast."""
    linkage_matrix = spearman_row_linkage(matrix)
    if linkage_matrix is None:
        row_order = list(matrix.index)
        raw_clusters = np.ones(len(row_order), dtype=int)
    else:
        dendrogram_data = dendrogram(linkage_matrix, labels=list(matrix.index), no_plot=True)
        row_order = [str(label) for label in dendrogram_data["ivl"]]
        raw_cluster_by_position = fcluster(
            linkage_matrix,
            t=min(max_clusters, len(matrix.index)),
            criterion="maxclust",
        )
        raw_clusters = np.array(
            [raw_cluster_by_position[list(matrix.index).index(label)] for label in row_order],
            dtype=int,
        )

    cluster_id_by_raw: dict[int, int] = {}
    cluster_ids: list[int] = []
    for raw_cluster in raw_clusters:
        raw = int(raw_cluster)
        if raw not in cluster_id_by_raw:
            cluster_id_by_raw[raw] = len(cluster_id_by_raw) + 1
        cluster_ids.append(cluster_id_by_raw[raw])

    stats = contrast_row_stats(cells).reindex(row_order)
    row_records = stats.reset_index().rename(columns={"condition_body_site_pair": "row_label"})
    row_records.insert(0, "row_order", np.arange(1, len(row_records) + 1))
    row_records.insert(1, "spearman_cluster", cluster_ids)

    summary_records: list[dict[str, object]] = []
    for cluster_id, group in row_records.groupby("spearman_cluster", sort=True):
        dominant_body_site = group["body_site_group"].astype(str).value_counts().index[0]
        dominant_condition = group["condition"].astype(str).value_counts().index[0]
        label = f"{dominant_body_site}; {dominant_condition}"
        summary_records.append(
            {
                "spearman_cluster": int(cluster_id),
                "cluster_label": label,
                "row_count": len(group),
                "displayed_cell_count": int(group["displayed_trait_count"].sum()),
                "significant_cell_count": int(group["significant_trait_count"].sum()),
                "bh_selected_cell_count": int(group["bh_selected_trait_count"].sum()),
                "dominant_body_site_group": dominant_body_site,
                "dominant_grouped_condition": dominant_condition,
                "top_body_site_groups": top_values(group["body_site_group"]),
                "top_grouped_conditions": top_values(group["condition"]),
            }
        )
    summary = pd.DataFrame.from_records(summary_records)
    row_records = row_records.merge(
        summary[["spearman_cluster", "cluster_label"]],
        on="spearman_cluster",
        how="left",
    )
    return row_records, summary, linkage_matrix


def cell_matrix(
    cells: pd.DataFrame,
    row_labels: list[str],
    trait_order: pd.DataFrame,
    *,
    value_column: str,
) -> pd.DataFrame:
    """Pivot significant heatmap cells into a row x trait matrix."""
    trait_ids = trait_order["trait_id"].tolist()
    significant = cells[cells["displayed_in_heatmap"]].copy()
    matrix = significant.pivot_table(
        index="condition_body_site_pair",
        columns="trait_id",
        values=value_column,
        aggfunc="max",
    )
    return matrix.reindex(index=row_labels, columns=trait_ids)


def flag_matrix(
    cells: pd.DataFrame,
    row_labels: list[str],
    trait_order: pd.DataFrame,
    *,
    flag_column: str,
) -> pd.DataFrame:
    """Pivot a boolean cell flag into row x trait layout."""
    trait_ids = trait_order["trait_id"].tolist()
    flags = cells[cells[flag_column]].copy()
    matrix = flags.pivot_table(
        index="condition_body_site_pair",
        columns="trait_id",
        values=flag_column,
        aggfunc="max",
    )
    return matrix.reindex(index=row_labels, columns=trait_ids).fillna(False).astype(bool)


def annotation_matrix(
    cells: pd.DataFrame,
    row_labels: list[str],
    trait_order: pd.DataFrame,
    *,
    mode: str,
) -> pd.DataFrame:
    """Return compact annotations for significant cells."""
    trait_ids = trait_order["trait_id"].tolist()
    annotations = pd.DataFrame("", index=row_labels, columns=trait_ids)
    significant = cells[cells["displayed_in_heatmap"]].copy()
    for row in significant.itertuples(index=False):
        if row.condition_body_site_pair not in annotations.index:
            continue
        if row.trait_id not in annotations.columns:
            continue
        p_text = format_p_value(row.fisher_p_value)
        p_label = f"p{p_text}" if p_text.startswith("<") else f"p={p_text}"
        if mode == "pattern":
            annotations.loc[row.condition_body_site_pair, row.trait_id] = (
                f"n={row.focus_positive_profile_count}\n{p_label}"
            )
        else:
            delta = row.increased_minus_decreased_positive_fraction
            annotations.loc[row.condition_body_site_pair, row.trait_id] = f"{delta:+.2f}"
    return annotations


def plot_matrix_heatmap(
    matrix: pd.DataFrame,
    annotations: pd.DataFrame,
    trait_order: pd.DataFrame,
    output_path: Path,
    *,
    title: str,
    subtitle: str,
    colorbar_label: str,
    cmap_name: str,
    vmin: float | None = None,
    vmax: float | None = None,
    compact: bool = False,
    border_mask: pd.DataFrame | None = None,
) -> None:
    """Render one row x six-trait heatmap with blank nonsignificant cells."""
    row_count = max(1, len(matrix.index))
    fig_height = max(6.2, 2.8 + 0.36 * row_count)
    fig, ax = plt.subplots(figsize=(13.8, fig_height))
    fig.subplots_adjust(left=0.34, right=0.90, top=0.82, bottom=0.18)
    if compact:
        fig.subplots_adjust(top=1 - 1.15 / fig_height, bottom=1.2 / fig_height)
    fig.patch.set_facecolor("white")
    ax.set_facecolor(SLATE_100)

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(SLATE_100)
    values = (
        matrix.to_numpy(dtype=float)
        if len(matrix.index)
        else np.full((1, len(trait_order)), np.nan)
    )
    image = ax.imshow(values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    if matrix.empty:
        ax.text(0.5, 0.5, "No qualifying contrasts", transform=ax.transAxes, ha="center")

    ax.set_xticks(
        np.arange(len(trait_order)),
        [wrap_label(label, width=16) for label in trait_order["trait_label"]],
    )
    ax.set_yticks(
        np.arange(len(matrix.index)),
        [wrap_label(label, width=42) for label in matrix.index],
    )
    ax.tick_params(axis="x", labelrotation=35, labelsize=8.8)
    ax.tick_params(axis="y", labelsize=8.6)
    ax.set_xlabel("Six shared traits", color=SLATE_700)
    ax.set_ylabel("Condition | body-site pair", color=SLATE_700)

    ax.set_xticks(np.arange(-0.5, len(trait_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(matrix.index), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.3)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for y_index, row_label in enumerate(matrix.index):
        for x_index, trait_id in enumerate(trait_order["trait_id"]):
            if border_mask is not None and bool(border_mask.loc[row_label, trait_id]):
                ax.add_patch(
                    Rectangle(
                        (x_index - 0.5, y_index - 0.5),
                        1,
                        1,
                        fill=False,
                        edgecolor=SLATE_900,
                        linewidth=1.8,
                    )
                )
            text = annotations.loc[row_label, trait_id]
            if not text:
                continue
            value = matrix.loc[row_label, trait_id]
            text_color = SLATE_900
            if pd.notna(value):
                red, green, blue, _alpha = image.cmap(image.norm(float(value)))
                luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
                text_color = "white" if luminance < 0.56 else SLATE_900
            ax.text(
                x_index,
                y_index,
                text,
                ha="center",
                va="center",
                fontsize=7.2,
                color=text_color,
            )

    fig.text(
        0.04,
        1 - 0.12 / fig_height if compact else 0.965,
        title,
        ha="left",
        va="top",
        fontsize=15.5,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.04,
        1 - 0.45 / fig_height if compact else 0.925,
        "\n".join(wrap(subtitle, width=130)),
        ha="left",
        va="top",
        fontsize=9.2,
        color=SLATE_500,
    )
    cbar = fig.colorbar(image, ax=ax, shrink=0.78, pad=0.02)
    cbar.set_label(colorbar_label, color=SLATE_700)
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_edgecolor(SLATE_200)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def cluster_color(cluster_id: int) -> tuple[float, float, float, float]:
    """Return a stable color for a row cluster."""
    palette = plt.get_cmap("tab20")
    return palette((cluster_id - 1) % 20)


def plot_grouped_clustered_contrast_heatmap(
    matrix: pd.DataFrame,
    annotations: pd.DataFrame,
    border_mask: pd.DataFrame,
    trait_order: pd.DataFrame,
    row_clusters: pd.DataFrame,
    linkage_matrix: np.ndarray | None,
    output_path: Path,
    *,
    preview: bool = False,
) -> None:
    """Render the grouped signed direction heatmap with a left Spearman dendrogram."""
    row_count = max(1, len(matrix.index))
    if preview:
        fig_height = 8.5
        fig_width = 11.0
        left = 0.08
        right = 0.92
        top = 0.86
        bottom = 0.12
        label_width = 0.01
        show_row_labels = False
        show_cell_text = False
    else:
        fig_height = max(9.5, 3.2 + 0.34 * row_count)
        fig_width = 20.5
        left = 0.04
        right = 0.965
        top = 1 - 1.55 / fig_height
        bottom = 0.85 / fig_height
        label_width = 3.2
        show_row_labels = True
        show_cell_text = True

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor="white")
    grid = fig.add_gridspec(
        1,
        5,
        width_ratios=[1.65, 0.14, 6.0, label_width, 0.32],
        left=left,
        right=right,
        top=top,
        bottom=bottom,
        wspace=0.06,
    )
    dendro_ax = fig.add_subplot(grid[0, 0])
    cluster_ax = fig.add_subplot(grid[0, 1])
    ax = fig.add_subplot(grid[0, 2])
    row_label_spacer_ax = fig.add_subplot(grid[0, 3])
    cbar_ax = fig.add_subplot(grid[0, 4])

    y_extent = row_count * 10
    if linkage_matrix is not None:
        dendrogram(
            linkage_matrix,
            orientation="left",
            no_labels=True,
            color_threshold=0,
            above_threshold_color=SLATE_500,
            ax=dendro_ax,
        )
        dendro_ax.set_ylim(y_extent, 0)
    dendro_ax.axis("off")

    for cluster_id, group in row_clusters.groupby("spearman_cluster", sort=False):
        positions = [
            matrix.index.get_loc(label) for label in group["row_label"] if label in matrix.index
        ]
        if not positions:
            continue
        start = min(positions) * 10
        height = (max(positions) - min(positions) + 1) * 10
        dendro_ax.axhspan(
            start,
            start + height,
            facecolor=cluster_color(int(cluster_id)),
            alpha=0.08,
            lw=0,
            zorder=0,
        )
        cluster_ax.add_patch(
            Rectangle((0, start), 1, height, facecolor=cluster_color(int(cluster_id)), lw=0)
        )
        y_mid = start + height / 2
        label_value = group["cluster_label"].dropna().iloc[0]
        label_text = compact_cluster_label(
            label_value,
            width=18 if not preview else 15,
        )
        font_size = 6.8 if not preview else 5.9
        if height <= 14:
            label_parts = [part for part in label_text.split("\n") if part]
            if label_parts[:2] == ["Urinary", "Renal/urinary"]:
                label_text = "Urinary/renal"
            else:
                label_text = " | ".join(label_parts)
            font_size = 5.8 if not preview else 5.2
        dendro_ax.text(
            0.96,
            y_mid,
            label_text,
            transform=dendro_ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=font_size,
            color=SLATE_900,
            fontweight="bold",
            bbox={
                "boxstyle": "round,pad=0.18",
                "facecolor": cluster_color(int(cluster_id)),
                "edgecolor": "white",
                "linewidth": 0.6,
                "alpha": 0.9,
            },
            clip_on=False,
            zorder=5,
        )
    cluster_ax.set_xlim(0, 1)
    cluster_ax.set_ylim(y_extent, 0)
    cluster_ax.set_xticks([])
    cluster_ax.set_yticks([])
    for spine in cluster_ax.spines.values():
        spine.set_visible(False)
    row_label_spacer_ax.patch.set_alpha(0)
    row_label_spacer_ax.set_zorder(-1)
    row_label_spacer_ax.axis("off")

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(SLATE_100)
    values = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
    image = ax.imshow(
        values,
        origin="lower",
        aspect="auto",
        cmap=cmap,
        vmin=-1.0,
        vmax=1.0,
        extent=[-0.5, len(trait_order) - 0.5, 0, y_extent],
    )
    ax.set_ylim(y_extent, 0)
    ax.set_xticks(
        np.arange(len(trait_order)),
        [wrap_label(label, width=16) for label in trait_order["trait_label"]],
    )
    ax.tick_params(axis="x", labelrotation=35, labelsize=8.5 if not preview else 7.2)
    if show_row_labels:
        ax.set_yticks(
            np.arange(len(matrix.index)) * 10 + 5,
            [wrap_label(label, width=48) for label in matrix.index],
        )
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.tick_params(
            axis="y",
            which="major",
            left=False,
            right=False,
            labelleft=False,
            labelright=True,
            labelsize=7.3,
            pad=7,
        )
    else:
        ax.set_yticks([])
    ax.set_xlabel("Six shared traits", color=SLATE_700)
    ax.set_ylabel(
        "" if preview else "Spearman-clustered grouped condition | anatomy rows",
        color=SLATE_700,
        labelpad=18,
    )
    ax.set_xticks(np.arange(-0.5, len(trait_order), 1), minor=True)
    ax.set_yticks(np.arange(0, y_extent + 1, 10), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False, right=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if show_cell_text:
        for y_index, row_label in enumerate(matrix.index):
            y_value = y_index * 10 + 5
            for x_index, trait_id in enumerate(trait_order["trait_id"]):
                if bool(border_mask.loc[row_label, trait_id]):
                    ax.add_patch(
                        Rectangle(
                            (x_index - 0.5, y_index * 10),
                            1,
                            10,
                            fill=False,
                            edgecolor=SLATE_900,
                            linewidth=1.9,
                        )
                    )
                text = annotations.loc[row_label, trait_id]
                if not text:
                    continue
                value = matrix.loc[row_label, trait_id]
                text_color = SLATE_900
                if pd.notna(value):
                    red, green, blue, _alpha = image.cmap(image.norm(float(value)))
                    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
                    text_color = "white" if luminance < 0.56 else SLATE_900
                ax.text(
                    x_index,
                    y_value,
                    text,
                    ha="center",
                    va="center",
                    fontsize=6.7,
                    color=text_color,
                )
    else:
        for y_index, row_label in enumerate(matrix.index):
            for x_index, trait_id in enumerate(trait_order["trait_id"]):
                if bool(border_mask.loc[row_label, trait_id]):
                    ax.add_patch(
                        Rectangle(
                            (x_index - 0.5, y_index * 10),
                            1,
                            10,
                            fill=False,
                            edgecolor=SLATE_900,
                            linewidth=0.85,
                        )
                    )

    title = "Grouped signed direction contrast: Spearman-clustered condition/anatomy rows"
    subtitle = (
        "Rows are grouped condition/body-site pairs retained when at least one six-trait cell "
        "has finite increased-vs-decreased support. Row dendrogram uses average linkage on "
        "Spearman-correlation distance across signed trait deltas, with missing cells set to "
        "neutral zero for clustering. Dark borders mark Fisher p<0.01."
    )
    fig.text(
        0.04,
        0.985 if preview else 1 - 0.12 / fig_height,
        title if not preview else "Grouped signed direction contrast preview",
        ha="left",
        va="top",
        fontsize=15.5 if not preview else 13.0,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.04,
        0.94 if preview else 1 - 0.52 / fig_height,
        "\n".join(wrap(subtitle, width=145 if not preview else 105)),
        ha="left",
        va="top",
        fontsize=8.7 if not preview else 7.6,
        color=SLATE_500,
    )
    cbar = fig.colorbar(image, cax=cbar_ax)
    cbar.set_label("Increased - decreased positive-z fraction", color=SLATE_700)
    cbar.ax.tick_params(labelsize=7.5)
    cbar.outline.set_edgecolor(SLATE_200)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path, dpi=220 if not preview else 180, bbox_inches="tight", facecolor="white"
    )
    plt.close(fig)


def write_pattern_heatmap(
    cells: pd.DataFrame,
    view: pd.Series,
    trait_order: pd.DataFrame,
    output_dir: Path,
) -> tuple[Path, int, int, int]:
    """Write one approach-2 condition/body-site x trait pattern heatmap."""
    row_labels = select_pattern_rows(cells)
    matrix = cell_matrix(
        cells,
        row_labels,
        trait_order,
        value_column="plot_minus_log10_fisher_p_value",
    )
    annotations = annotation_matrix(cells, row_labels, trait_order, mode="pattern")
    output_path = output_dir / (
        f"six_shared_traits_{view.view_key}_condition_body_site_trait_pattern_heatmap.png"
    )
    displayed_cells = int(cells["displayed_in_heatmap"].sum())
    subtitle = (
        "Approach 2: rows are condition/body-site pairs and columns are the six "
        "traits. Cell color is Fisher evidence, -log10(p), so darker cells have "
        "smaller p-values. Cells are shown for exploratory review only when Fisher "
        f"p<{P_VALUE_CUTOFF:g}, n>=3, and s>=2; BH q-values are written to the TSV."
    )
    plot_matrix_heatmap(
        matrix,
        annotations,
        trait_order,
        output_path,
        title=f"Condition/body-site trait pattern: {view.view_label}",
        subtitle=subtitle,
        colorbar_label=f"Fisher evidence (-log10 p, capped at {MAX_MINUS_LOG10_P:g})",
        cmap_name="YlGnBu",
        vmin=0.0,
        vmax=MAX_MINUS_LOG10_P,
    )
    return output_path, len(row_labels), displayed_cells, len(cells)


def write_contrast_heatmap(
    cells: pd.DataFrame,
    mode: ContrastMode,
    trait_order: pd.DataFrame,
    output_dir: Path,
) -> tuple[Path, int, int, int, dict[str, object]]:
    """Write one approach-4 signed direction contrast heatmap."""
    row_labels = select_contrast_rows(cells)
    matrix = cell_matrix(
        cells,
        row_labels,
        trait_order,
        value_column="increased_minus_decreased_positive_fraction",
    )
    annotations = annotation_matrix(cells, row_labels, trait_order, mode="contrast")
    border_mask = flag_matrix(
        cells,
        row_labels,
        trait_order,
        flag_column="significant_in_heatmap",
    )
    output_path = (
        output_dir
        / f"six_shared_traits_{mode.key}_condition_body_site_signed_direction_contrast_heatmap.png"
    )
    displayed_cells = int(cells["displayed_in_heatmap"].sum())
    if mode.key == "grouped":
        row_clusters, cluster_summary, row_linkage = clustered_contrast_rows(matrix, cells)
        ordered_labels = row_clusters["row_label"].astype(str).tolist()
        matrix = matrix.reindex(ordered_labels)
        annotations = annotations.reindex(ordered_labels)
        border_mask = border_mask.reindex(ordered_labels).fillna(False).astype(bool)

        cluster_by_row = row_clusters.set_index("row_label")["spearman_cluster"].to_dict()
        cluster_label_by_row = row_clusters.set_index("row_label")["cluster_label"].to_dict()
        in_rows = cells["condition_body_site_pair"].map(cluster_by_row)
        cells["spearman_cluster"] = in_rows.astype("Int64")
        cells["spearman_cluster_label"] = cells["condition_body_site_pair"].map(
            cluster_label_by_row
        )
        order_by_row = row_clusters.set_index("row_label")["row_order"].to_dict()
        cells["spearman_cluster_row_order"] = (
            cells["condition_body_site_pair"].map(order_by_row).astype("Int64")
        )

        cluster_rows_path = output_dir / GROUPED_CONTRAST_CLUSTER_ROWS_FILENAME
        cluster_summary_path = output_dir / GROUPED_CONTRAST_CLUSTER_SUMMARY_FILENAME
        write_tsv(row_clusters, cluster_rows_path)
        write_tsv(cluster_summary, cluster_summary_path)

        preview_path = output_path.with_stem(f"{output_path.stem}_preview")
        plot_grouped_clustered_contrast_heatmap(
            matrix,
            annotations,
            border_mask,
            trait_order,
            row_clusters,
            row_linkage,
            output_path,
        )
        plot_grouped_clustered_contrast_heatmap(
            matrix,
            annotations,
            border_mask,
            trait_order,
            row_clusters,
            row_linkage,
            preview_path,
            preview=True,
        )
        for page_path in output_dir.glob(f"{output_path.stem}_page_*.png"):
            page_path.unlink()
        extra = {
            "display_rule": "all supported testable cells; p cutoff marks dark borders only",
            "page_count": 0,
            "preview_png_path": preview_path,
            "cluster_rows_path": cluster_rows_path,
            "cluster_summary_path": cluster_summary_path,
            "clustering_method": "row average linkage on Spearman-correlation distance",
        }
        return output_path, len(row_labels), displayed_cells, len(cells), extra

    subtitle = (
        "Color and numbers show increased minus decreased positive-z fractions. "
        "Red favors increased profiles; blue favors decreased profiles. "
        f"Dark borders mark Fisher p<{P_VALUE_CUTOFF:g}. "
        "All testable cells with >=3 positive-z profiles in either direction and >=2 total "
        "studies are shown. Gray = unavailable or insufficient support. "
        "Raw p and BH q are in the TSV."
    )
    plot_matrix_heatmap(
        matrix,
        annotations,
        trait_order,
        output_path,
        title=f"Signed direction contrast: {mode.label}",
        subtitle=subtitle,
        colorbar_label="Increased - decreased positive-z fraction",
        cmap_name="RdBu_r",
        vmin=-1.0,
        vmax=1.0,
        compact=True,
        border_mask=border_mask,
    )
    pages = contrast_row_pages(row_labels)
    for page, labels in enumerate(pages, start=1):
        plot_matrix_heatmap(
            matrix.loc[labels],
            annotations.loc[labels],
            trait_order,
            contrast_page_path(output_path, page),
            title=f"Signed direction contrast: {mode.label} ({page}/{len(pages)})",
            subtitle=subtitle,
            colorbar_label="Increased - decreased positive-z fraction",
            cmap_name="RdBu_r",
            vmin=-1.0,
            vmax=1.0,
            compact=True,
            border_mask=border_mask.loc[labels],
        )
    extra = {
        "display_rule": "all supported testable cells; p cutoff marks dark borders only",
        "page_count": len(pages),
        "preview_png_path": "",
        "cluster_rows_path": "",
        "cluster_summary_path": "",
        "clustering_method": "",
    }
    return output_path, len(row_labels), displayed_cells, len(cells), extra


def contrast_row_pages(row_labels: list[str]) -> list[list[str]]:
    """Partition the complete row order without dropping or duplicating rows."""
    return [
        row_labels[start : start + CONTRAST_ROWS_PER_PAGE]
        for start in range(0, len(row_labels), CONTRAST_ROWS_PER_PAGE)
    ]


def contrast_page_path(output_path: Path, page: int) -> Path:
    """Return a stable path for a readable page of the full contrast overview."""
    return output_path.with_stem(f"{output_path.stem}_page_{page:02d}")


def append_readme_section(output_dir: Path, manifest: pd.DataFrame) -> None:
    """Append or replace the Cytoscape heatmap README section."""
    readme_path = output_dir / "README.md"
    if not readme_path.exists():
        return
    new_marker = "## Condition-Body-Site Cytoscape Heatmaps"
    old_marker = "## Condition-Body-Site Trait Correlation Heatmaps"
    content = readme_path.read_text(encoding="utf-8")
    split_marker = old_marker if old_marker in content else new_marker
    prefix = content.split(split_marker, maxsplit=1)[0].rstrip()
    lines = [
        "",
        "",
        new_marker,
        "",
        f"- Manifest: `{MANIFEST_FILENAME}`",
        f"- Pattern matrix table: `{PATTERN_TABLE_FILENAME}`",
        f"- Signed contrast table: `{CONTRAST_TABLE_FILENAME}`",
        f"- Grouped signed-direction cluster rows: `{GROUPED_CONTRAST_CLUSTER_ROWS_FILENAME}`",
        "- Grouped signed-direction cluster summary: "
        f"`{GROUPED_CONTRAST_CLUSTER_SUMMARY_FILENAME}`",
        f"- BH-selected cell table: `{SELECTED_CELLS_FILENAME}`",
        "",
        "Trait-pattern heatmaps use the six-shared-trait display-edge TSVs that feed "
        "the Cytoscape graph views. Signed contrasts use the full increased and decreased "
        "edge TSVs before graph display selection. They preserve condition/body-site/trait records "
        "instead of correlating already-collapsed condition and body-site nodes.",
        "",
        "Signed contrasts display all supported, testable effects, including nonsignificant "
        "results. Gray cells are unavailable or under-supported. Trait-pattern heatmaps "
        "continue to mask cells that fail the raw p-value threshold.",
        "",
        "- `trait_pattern`: Approach 2. Rows are condition/body-site pairs and "
        "columns are the six traits. Cell color is Fisher evidence, `-log10(p)`, "
        "rendered with a yellow-to-blue sequential palette so darker cells have "
        f"smaller p-values. Cells require Fisher p<{P_VALUE_CUTOFF:g}, n>=3, "
        "and s>=2 for exploratory display.",
        "- `signed_direction_contrast`: Approach 4. Rows are condition/body-site "
        "pairs and columns are traits. Cell color and numbers are "
        "`increased_fraction - decreased_fraction` on a fixed -1 to +1 scale. Red "
        "favors increased profiles and blue favors decreased profiles. "
        f"A dark border marks Fisher p<{P_VALUE_CUTOFF:g}; unbordered cells are descriptive "
        "effects, not significant findings at that cutoff.",
        "- Condition/body-site rows are chosen from the full increased and decreased edge TSVs, "
        "not from the graph display subset. A row is retained when at least one trait cell has "
        "a finite increased-vs-decreased Fisher test, at least 3 positive-z profiles in either "
        "direction, and at least 2 total studies.",
        "- The grouped signed-direction heatmap is a single clustered diagram. Rows use average "
        "linkage on Spearman-correlation distance across the six signed trait deltas; missing "
        "or unsupported cells are neutral zero for clustering. The left dendrogram clusters "
        "condition/anatomy rows, and the cluster audit tables summarize dominant anatomies and "
        "grouped condition labels.",
        "- Individual-condition signed-direction heatmaps still receive readable pages containing "
        f"at most {CONTRAST_ROWS_PER_PAGE} rows each. `{CONTRAST_PAGES_FILENAME}` records only "
        "those individual-condition page paths and counts.",
        f"- BH q-values are computed within each heatmap view. "
        f"`selected_after_bh=True` marks supported cells with q<={Q_VALUE_CUTOFF:g}; "
        "those cells are also collected in the BH-selected cell table.",
        "",
    ]
    grouped = manifest[manifest["view_key"].eq("grouped_direction_contrast")]
    if not grouped.empty and str(grouped.iloc[0].get("preview_png_path", "")):
        preview_path = Path(str(grouped.iloc[0]["preview_png_path"])).name
        full_path = Path(str(grouped.iloc[0]["png_path"])).name
        lines.extend(
            [
                "### Grouped Signed-Direction Cluster Preview",
                "",
                f"[![Grouped signed-direction cluster preview]({preview_path})]({full_path})",
                "",
            ]
        )
    lines.extend(
        [
            "| heatmap | type | PNG | displayed cells | BH-selected cells | tested cells | rows |",
            "| --- | --- | --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in manifest.itertuples(index=False):
        lines.append(
            f"| {row.view_label} | {row.heatmap_type} | `{row.png_path}` | "
            f"{row.displayed_cell_count} | {row.bh_selected_cell_count} | "
            f"{row.tested_cell_count} | {row.row_count} |"
        )
    lines.append("")
    pages_path = output_dir / CONTRAST_PAGES_FILENAME
    if pages_path.exists():
        lines.extend(["### Signed Direction Pages", ""])
        for page in read_tsv(pages_path).itertuples(index=False):
            lines.append(
                f"- [{page.view_label}, page {page.page}]({Path(page.png_path).name}) "
                f"({page.row_count} rows, {page.displayed_cell_count} filled cells)"
            )
        lines.append("")
    readme_path.write_text(prefix + "\n".join(lines), encoding="utf-8")


def write_six_shared_trait_cytoscape_heatmaps(
    output_dir: Path = OUTPUT_DIR,
    repo_root: Path = REPO_ROOT,
) -> pd.DataFrame:
    """Write approach-2 and approach-4 Cytoscape-derived heatmaps."""
    clean_legacy_outputs(output_dir)
    trait_order = read_trait_order(output_dir)
    graph_manifest = read_tsv(output_dir / GRAPH_MANIFEST_FILENAME)
    views_by_key = {row.view_key: row for row in graph_manifest.itertuples(index=False)}

    pattern_tables: list[pd.DataFrame] = []
    contrast_tables: list[pd.DataFrame] = []
    manifest_records: list[dict[str, object]] = []
    page_records: list[dict[str, object]] = []

    for view in graph_manifest.itertuples(index=False):
        display_edges = read_tsv(resolve_repo_path(view.source_display_edges_path, repo_root))
        cells = prepare_pattern_cells(display_edges, view, trait_order)
        pattern_tables.append(cells)
        output_path, row_count, displayed_count, tested_count = write_pattern_heatmap(
            cells,
            view,
            trait_order,
            output_dir,
        )
        manifest_records.append(
            {
                "view_key": view.view_key,
                "view_label": view.view_label,
                "heatmap_type": "trait_pattern",
                "focus_direction": view.focus_direction,
                "condition_mode": view.condition_mode,
                "row_count": row_count,
                "displayed_cell_count": displayed_count,
                "bh_selected_cell_count": int(cells["selected_after_bh"].sum()),
                "tested_cell_count": tested_count,
                "min_positive_profiles": MIN_POSITIVE_PROFILES,
                "min_studies": MIN_STUDIES,
                "p_value_cutoff": P_VALUE_CUTOFF,
                "q_value_cutoff": Q_VALUE_CUTOFF,
                "color_value": "-log10(fisher_p_value)",
                "png_path": relative_path(output_path, repo_root),
                "source_display_edges_path": view.source_display_edges_path,
            }
        )

    for mode in CONTRAST_MODES:
        increased_view = views_by_key[mode.increased_view_key]
        decreased_view = views_by_key[mode.decreased_view_key]
        increased_path = full_edges_path(increased_view.source_display_edges_path, repo_root)
        decreased_path = full_edges_path(decreased_view.source_display_edges_path, repo_root)
        cells = prepare_contrast_cells(
            read_tsv(increased_path),
            read_tsv(decreased_path),
            mode,
            trait_order,
        )
        output_path, row_count, displayed_count, tested_count, extra = write_contrast_heatmap(
            cells,
            mode,
            trait_order,
            output_dir,
        )
        cells["heatmap_page"] = pd.Series(pd.NA, index=cells.index, dtype="Int64")
        if mode.key != "grouped":
            for page, labels in enumerate(
                contrast_row_pages(select_contrast_rows(cells)),
                start=1,
            ):
                in_page = cells["condition_body_site_pair"].isin(labels)
                cells.loc[in_page, "heatmap_page"] = page
                page_records.append(
                    {
                        "view_key": f"{mode.key}_direction_contrast",
                        "view_label": mode.label,
                        "page": page,
                        "row_count": len(labels),
                        "displayed_cell_count": int(
                            cells.loc[in_page, "displayed_in_heatmap"].sum()
                        ),
                        "significant_cell_count": int(
                            cells.loc[in_page, "significant_in_heatmap"].sum()
                        ),
                        "png_path": relative_path(contrast_page_path(output_path, page), repo_root),
                    }
                )
        contrast_tables.append(cells)
        preview_path = extra.get("preview_png_path", "")
        cluster_rows_path = extra.get("cluster_rows_path", "")
        cluster_summary_path = extra.get("cluster_summary_path", "")
        manifest_records.append(
            {
                "view_key": f"{mode.key}_direction_contrast",
                "view_label": f"{mode.label}: signed direction contrast",
                "heatmap_type": "signed_direction_contrast",
                "focus_direction": "increased_vs_decreased",
                "condition_mode": mode.condition_mode,
                "row_count": row_count,
                "displayed_cell_count": displayed_count,
                "bh_selected_cell_count": int(cells["selected_after_bh"].sum()),
                "tested_cell_count": tested_count,
                "min_positive_profiles": MIN_POSITIVE_PROFILES,
                "min_studies": MIN_STUDIES,
                "p_value_cutoff": P_VALUE_CUTOFF,
                "q_value_cutoff": Q_VALUE_CUTOFF,
                "color_value": "increased_minus_decreased_positive_fraction",
                "display_rule": extra["display_rule"],
                "significant_cell_count": int(cells["significant_in_heatmap"].sum()),
                "filled_fraction": displayed_count / (row_count * len(trait_order))
                if row_count
                else 0,
                "page_count": extra["page_count"],
                "png_path": relative_path(output_path, repo_root),
                "preview_png_path": relative_path(preview_path, repo_root)
                if isinstance(preview_path, Path)
                else "",
                "cluster_rows_path": relative_path(cluster_rows_path, repo_root)
                if isinstance(cluster_rows_path, Path)
                else "",
                "cluster_summary_path": relative_path(cluster_summary_path, repo_root)
                if isinstance(cluster_summary_path, Path)
                else "",
                "clustering_method": extra["clustering_method"],
                "source_display_edges_path": "",
                "source_increased_edges_path": relative_path(increased_path, repo_root),
                "source_decreased_edges_path": relative_path(decreased_path, repo_root),
            }
        )

    pattern_table = pd.concat(pattern_tables, ignore_index=True)
    contrast_table = pd.concat(contrast_tables, ignore_index=True)
    manifest = pd.DataFrame.from_records(manifest_records)
    selected_cells = pd.concat(
        [
            pattern_table.assign(heatmap_type="trait_pattern"),
            contrast_table.assign(heatmap_type="signed_direction_contrast"),
        ],
        ignore_index=True,
        sort=False,
    )
    selected_cells = selected_cells[selected_cells["selected_after_bh"]].sort_values(
        ["view_key", "fisher_q_value", "fisher_p_value", "condition_body_site_pair", "trait_rank"],
        ascending=[True, True, True, True, True],
    )

    write_tsv(pattern_table, output_dir / PATTERN_TABLE_FILENAME)
    write_tsv(contrast_table, output_dir / CONTRAST_TABLE_FILENAME)
    write_tsv(selected_cells, output_dir / SELECTED_CELLS_FILENAME)
    write_tsv(manifest, output_dir / MANIFEST_FILENAME)
    write_tsv(
        pd.DataFrame.from_records(
            page_records,
            columns=[
                "view_key",
                "view_label",
                "page",
                "row_count",
                "displayed_cell_count",
                "significant_cell_count",
                "png_path",
            ],
        ),
        output_dir / CONTRAST_PAGES_FILENAME,
    )
    append_readme_section(output_dir, manifest)
    return manifest


def main() -> None:
    """Write the Cytoscape-derived heatmap outputs."""
    manifest = write_six_shared_trait_cytoscape_heatmaps()
    print(f"Wrote Cytoscape heatmaps: {len(manifest)}")
    print(f"Wrote manifest: {OUTPUT_DIR / MANIFEST_FILENAME}")


if __name__ == "__main__":
    main()
