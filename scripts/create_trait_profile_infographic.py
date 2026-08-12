"""Create an infographic summary of paired trait-profile contrasts."""

from __future__ import annotations

import math
import re
from pathlib import Path
from textwrap import fill, shorten

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle, FancyBboxPatch, PathPatch
from matplotlib.path import Path as MplPath
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"

PAIRED_MATRIX_PATH = VIS_DIR / "study_trait_profile_paired_enrichment_matrix.tsv"
PAIRED_METADATA_PATH = VIS_DIR / "study_trait_profile_paired_enrichment_metadata.tsv"
PROFILE_MATRIX_PATH = REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv"
TRAIT_COLUMNS_PATH = REPO_ROOT / "data/derived/taxa_trait_matrix_columns.tsv"
SELECTED_TRAIT_METADATA_PATH = VIS_DIR / "study_trait_profile_trait_metadata.tsv"

OUTPUT_PNG = INVESTIGATION_DIR / "trait_profile_paired_contrast_infographic.png"
TRAIT_STATS_TSV = INVESTIGATION_DIR / "trait_profile_paired_contrast_trait_stats.tsv"
POOLED_TRAIT_STATS_TSV = INVESTIGATION_DIR / "trait_profile_paired_contrast_trait_stats_pooled.tsv"
STRATIFIED_TRAIT_STATS_MANIFEST_TSV = (
    INVESTIGATION_DIR / "trait_profile_paired_contrast_trait_stats_manifest.tsv"
)
CONDITION_SUMMARY_TSV = INVESTIGATION_DIR / "trait_profile_paired_contrast_condition_summary.tsv"
BODY_SITE_SUMMARY_TSV = INVESTIGATION_DIR / "trait_profile_paired_contrast_body_site_summary.tsv"

TOP_N_TRAITS = 10
TOP_N_CONDITIONS = 10
TOP_N_BODY_SITES = 9
STRATIFIED_TRAIT_COUNT: int | None = None
MIN_STRATIFIED_PAIRS = 50
MIN_STRATIFIED_PROFILE_SUPPORT = 6
MAX_STRATIFIED_PROFILE_SUPPORT_FRACTION = 0.95
MIN_CONDITION_PAIRS = 10
MIN_BODY_SITE_PAIRS = 3
PROFILE_METADATA_COLUMNS = ["study_id", "experiment_id", "inc_dec"]

EXCLUDED_CONDITION_LABELS = {
    "Body fat distribution",
    "Clinical treatment",
    "Exposure",
    "Health trait",
    "Microbiome",
    "Population",
    "Physical activity",
}

SLATE_50 = "#f8fafc"
SLATE_200 = "#e2e8f0"
SLATE_400 = "#94a3b8"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
INCREASED_BLUE = "#3b82f6"
DECREASED_RED = "#ef4444"


def clean_text(value: object) -> str:
    """Normalize empty metadata and display values."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def compact_label(value: object, max_length: int = 42) -> str:
    """Return a compact single-line label for axes."""
    return shorten(clean_text(value).replace("\n", " "), width=max_length, placeholder="...")


def trait_label(value: object, max_length: int = 46) -> str:
    """Return the readable part of a trait display label."""
    label = clean_text(value).replace("\n", " ")
    return shorten(label, width=max_length, placeholder="...")


def full_trait_label(value: object) -> str:
    """Return a non-wrapped trait label for TSV outputs."""
    return clean_text(value).replace("\n", " ")


def safe_filename_token(value: object) -> str:
    """Convert a stratum label into a stable lowercase filename token."""
    token = re.sub(r"[^a-z0-9]+", "_", clean_text(value).lower()).strip("_")
    return token or "unknown"


def format_trait_label(column_id: str, row: pd.Series) -> str:
    """Build a compact display label from trait column metadata."""
    name = clean_text(row.get("object_name", "")).replace("_", " ")
    if name == "unknown":
        name = clean_text(row.get("column_label", "")).replace("_", " ")
    short_id = column_id.split("|")[-1]
    if short_id.startswith("pathways:"):
        short_id = short_id.replace("pathways:", "")
    return full_trait_label(f"{name} {short_id}")


def load_trait_labels(path: Path = TRAIT_COLUMNS_PATH) -> pd.DataFrame:
    """Load trait labels for every raw study-trait profile column."""
    labels = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    labels = labels.set_index("column_id", drop=False)
    labels["display_label"] = [
        format_trait_label(column_id, row) for column_id, row in labels.iterrows()
    ]
    labels["trait_label"] = labels["display_label"].map(full_trait_label)
    labels["display_category"] = labels["object_category"].str.removeprefix("biolink:")
    labels.loc[labels["display_category"] == "", "display_category"] = "unknown"
    return labels


def load_selected_trait_metadata(path: Path = SELECTED_TRAIT_METADATA_PATH) -> pd.DataFrame:
    """Load metadata for traits in the existing paired heatmap matrix."""
    metadata = pd.read_csv(path, sep="\t", index_col=0, dtype=str).fillna("")
    if "trait_id" in metadata.columns:
        metadata = metadata.set_index("trait_id", drop=False)
    else:
        metadata["trait_id"] = metadata.index
    if "trait_label" not in metadata.columns:
        source = metadata.get("display_label", metadata["trait_id"])
        metadata["trait_label"] = source.map(full_trait_label)
    for column in ["profile_support", "total_signal", "composition_variance", "selection_score"]:
        if column in metadata.columns:
            metadata[column] = pd.to_numeric(metadata[column], errors="coerce")
    metadata["selected_trait_rank"] = np.arange(1, len(metadata) + 1)
    return metadata


def metadata_from_trait_ids(trait_ids: pd.Index) -> pd.DataFrame:
    """Build minimal trait metadata when no richer metadata table is available."""
    return pd.DataFrame(
        {
            "trait_id": trait_ids,
            "trait_label": [full_trait_label(trait_id) for trait_id in trait_ids],
            "selected_trait_rank": np.arange(1, len(trait_ids) + 1),
        },
        index=trait_ids,
    )


def load_profile_matrix(path: Path = PROFILE_MATRIX_PATH) -> pd.DataFrame:
    """Load the full raw study/experiment/direction trait-count matrix."""
    profiles = pd.read_csv(path, sep="\t")
    trait_columns = [
        column for column in profiles.columns if column not in PROFILE_METADATA_COLUMNS
    ]
    values = profiles[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    profile_ids = [
        f"{study_id}|{experiment_id}|{inc_dec}"
        for study_id, experiment_id, inc_dec in profiles[PROFILE_METADATA_COLUMNS].itertuples(
            index=False,
            name=None,
        )
    ]
    values.index = profile_ids
    return values.groupby(level=0, sort=False).sum()


def split_profile_ids(value: object) -> list[str]:
    """Split the semicolon-delimited profile IDs stored in paired metadata."""
    return [profile_id for profile_id in clean_text(value).split(";") if profile_id != "unknown"]


def profile_ids_for_pairs(pair_metadata: pd.DataFrame) -> list[str]:
    """Return unique increased/decreased profile IDs referenced by matched pairs."""
    profile_ids: list[str] = []
    for row in pair_metadata.itertuples(index=False):
        profile_ids.extend(split_profile_ids(row.increased_profile_ids))
        profile_ids.extend(split_profile_ids(row.decreased_profile_ids))
    return list(dict.fromkeys(profile_ids))


def select_stratum_traits(
    profile_matrix: pd.DataFrame,
    trait_labels: pd.DataFrame,
    *,
    max_traits: int | None = STRATIFIED_TRAIT_COUNT,
    min_profile_support: int = MIN_STRATIFIED_PROFILE_SUPPORT,
    max_profile_support_fraction: float = MAX_STRATIFIED_PROFILE_SUPPORT_FRACTION,
) -> pd.DataFrame:
    """Select supported variable traits within a stratum's profile-count matrix."""
    row_sums = profile_matrix.sum(axis=1).replace(0, pd.NA)
    composition = profile_matrix.div(row_sums, axis=0).fillna(0.0)
    profile_support = profile_matrix.gt(0).sum(axis=0).astype(int)
    total_signal = profile_matrix.sum(axis=0).astype(float)
    composition_variance = composition.var(axis=0)

    trait_metadata = pd.DataFrame(
        {
            "trait_id": composition.columns,
            "profile_support": profile_support,
            "total_signal": total_signal,
            "composition_variance": composition_variance,
        },
        index=composition.columns,
    ).join(
        trait_labels[
            [
                "trait_label",
                "display_label",
                "column_label",
                "object_name",
                "object_category",
                "display_category",
            ]
        ],
        how="left",
    )
    trait_metadata["trait_label"] = trait_metadata["trait_label"].fillna(
        pd.Series(trait_metadata.index, index=trait_metadata.index)
    )
    trait_metadata["selection_score"] = [
        variance * math.log1p(signal)
        for variance, signal in zip(
            trait_metadata["composition_variance"],
            trait_metadata["total_signal"],
            strict=True,
        )
    ]

    max_support = max(2, int(len(profile_matrix) * max_profile_support_fraction))
    eligible = trait_metadata[
        (trait_metadata["profile_support"] >= min_profile_support)
        & (trait_metadata["profile_support"] <= max_support)
        & (trait_metadata["composition_variance"] > 0)
    ].copy()
    eligible = eligible[
        ~eligible.index.to_series().str.contains("mediadive.medium", regex=False)
        & ~eligible["trait_label"].str.contains("mediadive.medium", regex=False)
    ].copy()
    eligible = eligible.sort_values(
        ["selection_score", "composition_variance", "profile_support", "trait_id"],
        ascending=[False, False, False, True],
    )
    if max_traits is not None:
        eligible = eligible.head(max_traits)
    eligible["selected_trait_rank"] = np.arange(1, len(eligible) + 1)
    if len(eligible) < 2:
        msg = "Fewer than two traits passed the stratum selection filters."
        raise ValueError(msg)
    return eligible


def build_zscore_matrix(profile_matrix: pd.DataFrame, selected_traits: pd.Index) -> pd.DataFrame:
    """Create column z-scores from log-scaled trait fractions."""
    row_sums = profile_matrix.sum(axis=1).replace(0, pd.NA)
    composition = profile_matrix[selected_traits].div(row_sums, axis=0).fillna(0.0)
    scaled = (composition * 1000).apply(lambda column: column.map(math.log1p))
    zscores = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0, ddof=0).replace(0, pd.NA)
    return zscores.fillna(0.0)


def annotation_counts_for_traits(profile_matrix: pd.DataFrame, traits: pd.Index) -> pd.Series:
    """Count profiles with nonzero raw annotation signal for each trait."""
    counts = pd.Series(0, index=traits, dtype=int)
    available_traits = traits.intersection(profile_matrix.columns)
    if len(available_traits) > 0:
        counts.loc[available_traits] = profile_matrix[available_traits].gt(0).sum(axis=0)
    return counts


def build_pair_matrix_from_profiles(
    profile_zscores: pd.DataFrame,
    pair_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Build matched increased-minus-decreased rows from profile-level z-scores."""
    rows: list[pd.Series] = []
    for pair_id, row in pair_metadata.iterrows():
        inc_ids = pd.Index(split_profile_ids(row["increased_profile_ids"])).intersection(
            profile_zscores.index
        )
        dec_ids = pd.Index(split_profile_ids(row["decreased_profile_ids"])).intersection(
            profile_zscores.index
        )
        if len(inc_ids) == 0 or len(dec_ids) == 0:
            continue
        pair_values = profile_zscores.loc[inc_ids].mean(axis=0) - profile_zscores.loc[dec_ids].mean(
            axis=0
        )
        pair_values.name = pair_id
        rows.append(pair_values)
    if not rows:
        msg = "No matched increased/decreased pairs were available for this stratum."
        raise ValueError(msg)
    return pd.DataFrame(rows)


def is_excluded_condition_label(value: object) -> bool:
    """Return whether a condition label is too broad for condition-panel ranking."""
    label = clean_text(value)
    lower_label = label.lower()
    broad_terms = ("exposure", "measurement", "population")
    return (
        label in EXCLUDED_CONDITION_LABELS
        or label == "unknown"
        or any(term in lower_label for term in broad_terms)
    )


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Calculate Benjamini-Hochberg adjusted q-values."""
    values = p_values.to_numpy(dtype=float)
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


def load_paired_inputs(
    matrix_path: Path = PAIRED_MATRIX_PATH,
    metadata_path: Path = PAIRED_METADATA_PATH,
    selected_trait_metadata_path: Path = SELECTED_TRAIT_METADATA_PATH,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load paired contrast matrix and metadata, aligned by row order."""
    matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
    if len(matrix) != len(metadata):
        msg = (
            f"Paired matrix row count ({len(matrix)}) does not match metadata "
            f"row count ({len(metadata)})."
        )
        raise ValueError(msg)

    if selected_trait_metadata_path.exists():
        trait_metadata = load_selected_trait_metadata(selected_trait_metadata_path)
        if len(trait_metadata) == len(matrix.columns):
            matrix.columns = trait_metadata.index
        else:
            trait_metadata = metadata_from_trait_ids(pd.Index(matrix.columns))
    else:
        trait_metadata = metadata_from_trait_ids(pd.Index(matrix.columns))

    # The displayed paired matrix uses readable labels as its index, while the
    # metadata keeps stable pair IDs. The files are written in the same order.
    matrix.index = metadata.index
    return matrix, metadata, trait_metadata.loc[matrix.columns]


def compute_trait_stats(
    matrix: pd.DataFrame,
    trait_metadata: pd.DataFrame | None = None,
    annotation_counts: pd.Series | None = None,
    *,
    stratum_kind: str = "pooled",
    stratum: str = "pooled",
    selected_profile_count: int | None = None,
) -> pd.DataFrame:
    """Calculate paired one-sample t-test summaries for each trait column."""
    if trait_metadata is None:
        trait_metadata = metadata_from_trait_ids(pd.Index(matrix.columns))
    trait_metadata = trait_metadata.reindex(matrix.columns)
    if annotation_counts is None:
        annotation_counts = pd.Series(pd.NA, index=matrix.columns)
    annotation_counts = annotation_counts.reindex(matrix.columns)

    records: list[dict[str, object]] = []
    for column in matrix.columns:
        values = matrix[column].to_numpy(dtype=float)
        if len(values) < 3 or np.isclose(values.std(ddof=1), 0.0):
            t_stat = float("nan")
            p_value = float("nan")
        else:
            result = stats.ttest_1samp(values, 0.0)
            t_stat = float(result.statistic)
            p_value = float(result.pvalue)

        mean_contrast = float(values.mean())
        metadata_row = trait_metadata.loc[column]
        annotated_profile_count = annotation_counts.loc[column]
        annotated_profile_fraction = (
            float(annotated_profile_count) / selected_profile_count
            if selected_profile_count and pd.notna(annotated_profile_count)
            else pd.NA
        )
        records.append(
            {
                "trait": str(column),
                "trait_label": full_trait_label(metadata_row.get("trait_label", column)),
                "stratum_kind": stratum_kind,
                "stratum": stratum,
                "selected_trait_rank": metadata_row.get("selected_trait_rank", pd.NA),
                "pair_count": len(values),
                "profile_count": selected_profile_count
                if selected_profile_count is not None
                else pd.NA,
                "n_profiles_where_annotated": annotated_profile_count,
                "annotated_profile_fraction": annotated_profile_fraction,
                "profile_support": metadata_row.get("profile_support", pd.NA),
                "total_signal": metadata_row.get("total_signal", pd.NA),
                "composition_variance": metadata_row.get("composition_variance", pd.NA),
                "selection_score": metadata_row.get("selection_score", pd.NA),
                "mean_increased_minus_decreased": mean_contrast,
                "median_increased_minus_decreased": float(np.median(values)),
                "std_increased_minus_decreased": float(values.std(ddof=1)),
                "positive_pair_fraction": float(np.mean(values > 0)),
                "negative_pair_fraction": float(np.mean(values < 0)),
                "nonzero_pair_fraction": float(np.mean(values != 0)),
                "t_statistic": t_stat,
                "p_value": p_value,
                "direction": (
                    "trait enriched in increased profiles"
                    if mean_contrast >= 0
                    else "trait enriched in decreased profiles"
                ),
            }
        )

    stats_df = pd.DataFrame.from_records(records)
    stats_df["q_value"] = benjamini_hochberg(stats_df["p_value"].fillna(1.0))
    clipped_p = stats_df["p_value"].clip(lower=np.finfo(float).tiny, upper=1.0)
    stats_df["minus_log10_p"] = -np.log10(clipped_p)
    stats_df["abs_t_statistic"] = stats_df["t_statistic"].abs()
    ordered_columns = [
        "trait",
        "trait_label",
        "stratum_kind",
        "stratum",
        "selected_trait_rank",
        "pair_count",
        "profile_count",
        "n_profiles_where_annotated",
        "annotated_profile_fraction",
        "profile_support",
        "total_signal",
        "composition_variance",
        "selection_score",
        "mean_increased_minus_decreased",
        "median_increased_minus_decreased",
        "std_increased_minus_decreased",
        "positive_pair_fraction",
        "negative_pair_fraction",
        "nonzero_pair_fraction",
        "t_statistic",
        "p_value",
        "q_value",
        "minus_log10_p",
        "abs_t_statistic",
        "direction",
    ]
    return stats_df[ordered_columns].sort_values(
        ["p_value", "abs_t_statistic", "trait_label"],
        ascending=[True, False, True],
    )


def add_pair_scores(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Attach per-pair directional trait contrast loads to metadata."""
    scored = metadata.copy()
    scored["increased_profile_higher_load"] = matrix.clip(lower=0).mean(axis=1)
    scored["decreased_profile_higher_load"] = (-matrix.clip(upper=0)).mean(axis=1)
    scored["absolute_contrast_load"] = matrix.abs().mean(axis=1)
    scored["net_increased_minus_decreased"] = matrix.mean(axis=1)
    return scored


def summarize_groups(
    scored: pd.DataFrame,
    group_column: str,
    *,
    min_pairs: int,
) -> pd.DataFrame:
    """Summarize paired profile contrast load by condition or body-site group."""
    source = scored.copy()
    if group_column == "condition":
        source = source[~source[group_column].map(is_excluded_condition_label)]
    else:
        labels = source[group_column].map(clean_text).str.lower()
        source = source[(labels != "unknown") & (labels != "other")]

    summary = (
        source.groupby(group_column, dropna=False)
        .agg(
            pair_count=("study_id", "size"),
            study_count=("study_id", "nunique"),
            increased_profile_higher_load=("increased_profile_higher_load", "mean"),
            decreased_profile_higher_load=("decreased_profile_higher_load", "mean"),
            absolute_contrast_load=("absolute_contrast_load", "mean"),
            net_increased_minus_decreased=("net_increased_minus_decreased", "mean"),
        )
        .reset_index()
    )
    summary = summary[summary["pair_count"] >= min_pairs].copy()
    summary["dominant_direction"] = np.where(
        summary["net_increased_minus_decreased"] >= 0,
        "trait enriched in increased profiles",
        "trait enriched in decreased profiles",
    )
    return summary.sort_values(
        ["absolute_contrast_load", "pair_count", group_column],
        ascending=[False, False, True],
    )


def trait_stats_path(output_dir: Path, token: str) -> Path:
    """Return the conventional trait-stats filename for a stratum token."""
    return output_dir / f"trait_profile_paired_contrast_trait_stats_{token}.tsv"


def manifest_path(output_dir: Path) -> Path:
    """Return the conventional stratified trait-stats manifest path."""
    return output_dir / "trait_profile_paired_contrast_trait_stats_manifest.tsv"


def display_path(path: Path) -> str:
    """Return a repo-relative path when possible, otherwise the absolute path."""
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def compute_stratum_trait_stats(
    profile_matrix: pd.DataFrame,
    pair_metadata: pd.DataFrame,
    trait_labels: pd.DataFrame,
    *,
    stratum_kind: str,
    stratum: str,
    max_traits: int | None = STRATIFIED_TRAIT_COUNT,
    min_profile_support: int = MIN_STRATIFIED_PROFILE_SUPPORT,
    max_profile_support_fraction: float = MAX_STRATIFIED_PROFILE_SUPPORT_FRACTION,
) -> pd.DataFrame:
    """Select traits and calculate paired contrast statistics within one stratum."""
    profile_ids = [
        profile_id
        for profile_id in profile_ids_for_pairs(pair_metadata)
        if profile_id in profile_matrix.index
    ]
    selected_profile_matrix = profile_matrix.loc[profile_ids]
    trait_metadata = select_stratum_traits(
        selected_profile_matrix,
        trait_labels,
        max_traits=max_traits,
        min_profile_support=min_profile_support,
        max_profile_support_fraction=max_profile_support_fraction,
    )
    profile_zscores = build_zscore_matrix(selected_profile_matrix, trait_metadata.index)
    pair_matrix = build_pair_matrix_from_profiles(profile_zscores, pair_metadata)
    annotation_counts = annotation_counts_for_traits(selected_profile_matrix, trait_metadata.index)
    return compute_trait_stats(
        pair_matrix,
        trait_metadata,
        annotation_counts,
        stratum_kind=stratum_kind,
        stratum=stratum,
        selected_profile_count=len(selected_profile_matrix),
    )


def body_site_strata(
    pair_metadata: pd.DataFrame, min_pairs: int
) -> list[tuple[str, str, pd.DataFrame]]:
    """Return body-site-group strata with enough matched pairs."""
    strata: list[tuple[str, str, pd.DataFrame]] = []
    for body_site_group, group in pair_metadata.groupby("body_site_group", sort=True, dropna=False):
        label = clean_text(body_site_group)
        if label.lower() == "unknown" or len(group) < min_pairs:
            continue
        strata.append((safe_filename_token(label), label, group))
    return strata


def write_stratified_trait_stats(
    profile_matrix: pd.DataFrame,
    pair_metadata: pd.DataFrame,
    trait_labels: pd.DataFrame,
    *,
    output_dir: Path = INVESTIGATION_DIR,
    min_pairs: int = MIN_STRATIFIED_PAIRS,
    max_traits: int | None = STRATIFIED_TRAIT_COUNT,
    min_profile_support: int = MIN_STRATIFIED_PROFILE_SUPPORT,
    max_profile_support_fraction: float = MAX_STRATIFIED_PROFILE_SUPPORT_FRACTION,
) -> pd.DataFrame:
    """Write per-stratum trait-stat TSVs and return a manifest of emitted files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []
    emitted_tokens: set[str] = set()

    explicit_strata = [
        (
            "feces",
            "body_site_group",
            "Feces",
            pair_metadata[pair_metadata["body_site_group"].map(clean_text) == "Feces"],
            "primary-analysis",
        ),
        (
            "feces_excluded",
            "body_site_group_excluded",
            "Feces excluded",
            pair_metadata[pair_metadata["body_site_group"].map(clean_text) != "Feces"],
            "robustness-feces-excluded",
        ),
    ]
    for token, stratum_kind, stratum, group, analysis_role in explicit_strata:
        if len(group) < min_pairs:
            continue
        stats_df = compute_stratum_trait_stats(
            profile_matrix,
            group,
            trait_labels,
            stratum_kind=stratum_kind,
            stratum=stratum,
            max_traits=max_traits,
            min_profile_support=min_profile_support,
            max_profile_support_fraction=max_profile_support_fraction,
        )
        path = trait_stats_path(output_dir, token)
        stats_df.to_csv(path, sep="\t", index=False)
        emitted_tokens.add(token)
        manifest_rows.append(
            {
                "stratum_kind": stratum_kind,
                "stratum": stratum,
                "filename_token": token,
                "analysis_role": analysis_role,
                "pair_count": len(group),
                "study_count": group["study_id"].nunique(),
                "selected_trait_count": len(stats_df),
                "path": display_path(path),
            }
        )

    for token, stratum, group in body_site_strata(pair_metadata, min_pairs):
        if token in emitted_tokens:
            continue
        stats_df = compute_stratum_trait_stats(
            profile_matrix,
            group,
            trait_labels,
            stratum_kind="body_site_group",
            stratum=stratum,
            max_traits=max_traits,
            min_profile_support=min_profile_support,
            max_profile_support_fraction=max_profile_support_fraction,
        )
        path = trait_stats_path(output_dir, token)
        stats_df.to_csv(path, sep="\t", index=False)
        manifest_rows.append(
            {
                "stratum_kind": "body_site_group",
                "stratum": stratum,
                "filename_token": token,
                "analysis_role": "body-site",
                "pair_count": len(group),
                "study_count": group["study_id"].nunique(),
                "selected_trait_count": len(stats_df),
                "path": display_path(path),
            }
        )

    manifest = pd.DataFrame.from_records(manifest_rows)
    manifest.to_csv(manifest_path(output_dir), sep="\t", index=False)
    return manifest


def anatomical_site_base(value: object) -> tuple[float, float]:
    """Map body-site group labels to approximate body-diagram coordinates."""
    label = clean_text(value).lower()
    if "nasopharynx" in label or "airway" in label:
        return 0.30, 0.86
    if "subgingival" in label or "dental" in label:
        return 0.275, 0.825
    if "saliva" in label:
        return 0.325, 0.825
    if "mouth" in label or "oral" in label:
        return 0.30, 0.825
    if "colon" in label:
        return 0.285, 0.47
    if "caecum" in label or "cecum" in label:
        return 0.33, 0.46
    if "feces" in label or "faeces" in label or "stool" in label:
        return 0.31, 0.37
    if "vagina" in label or "urogenital" in label:
        return 0.31, 0.25
    return 0.16, 0.47


def add_anatomical_outline(ax: plt.Axes) -> None:
    """Draw an original front-facing body-template outline for site callouts."""
    face = "#ffffff"
    edge = SLATE_400

    def add_path(
        start: tuple[float, float],
        segments: list[tuple[tuple[float, float], ...]],
    ) -> None:
        vertices = [start]
        codes = [MplPath.MOVETO]
        for segment in segments:
            if len(segment) == 1:
                vertices.append(segment[0])
                codes.append(MplPath.LINETO)
            else:
                vertices.extend(segment)
                codes.extend([MplPath.CURVE4] * 3)
        vertices.append(start)
        codes.append(MplPath.CLOSEPOLY)
        ax.add_patch(
            PathPatch(
                MplPath(vertices, codes),
                facecolor=face,
                edgecolor=edge,
                lw=1.2,
                joinstyle="round",
                zorder=1,
            )
        )

    def mirror_point(point: tuple[float, float]) -> tuple[float, float]:
        return 0.60 - point[0], point[1]

    def mirror_segments(
        segments: list[tuple[tuple[float, float], ...]],
    ) -> list[tuple[tuple[float, float], ...]]:
        return [tuple(mirror_point(point) for point in segment) for segment in segments]

    ax.add_patch(Circle((0.30, 0.87), 0.045, facecolor=face, edgecolor=edge, lw=1.2, zorder=1))
    torso_segments = [
        ((0.224, 0.800), (0.207, 0.775), (0.196, 0.735)),
        ((0.178, 0.660), (0.178, 0.565), (0.195, 0.480)),
        ((0.205, 0.420), (0.220, 0.355), (0.240, 0.320)),
        ((0.265, 0.295), (0.335, 0.295), (0.360, 0.320)),
        ((0.380, 0.355), (0.395, 0.420), (0.405, 0.480)),
        ((0.422, 0.565), (0.422, 0.660), (0.404, 0.735)),
        ((0.393, 0.775), (0.376, 0.800), (0.355, 0.805)),
        ((0.335, 0.812), (0.265, 0.812), (0.245, 0.805)),
    ]
    left_arm_segments = [
        ((0.175, 0.665), (0.160, 0.585), (0.142, 0.505)),
        ((0.127, 0.445), (0.112, 0.395), (0.105, 0.365)),
        ((0.113, 0.350), (0.135, 0.346), (0.150, 0.358)),
        ((0.162, 0.395), (0.178, 0.455), (0.198, 0.535)),
        ((0.207, 0.600), (0.213, 0.665), (0.224, 0.715)),
    ]
    left_leg_segments = [
        ((0.258, 0.245), (0.250, 0.170), (0.238, 0.090)),
        ((0.233, 0.060), (0.250, 0.045), (0.278, 0.050)),
        ((0.292, 0.100), (0.295, 0.225), (0.295, 0.300)),
    ]
    add_path((0.245, 0.805), torso_segments)
    add_path((0.205, 0.725), left_arm_segments)
    add_path(mirror_point((0.205, 0.725)), mirror_segments(left_arm_segments))
    add_path((0.265, 0.300), left_leg_segments)
    add_path(mirror_point((0.265, 0.300)), mirror_segments(left_leg_segments))
    ax.plot([0.280, 0.280], [0.805, 0.770], color=edge, lw=1.0, zorder=2)
    ax.plot([0.320, 0.320], [0.805, 0.770], color=edge, lw=1.0, zorder=2)


def plot_top_traits(ax: plt.Axes, trait_stats: pd.DataFrame) -> None:
    """Plot panel A: traits with strongest paired contrast evidence."""
    top_traits = trait_stats.head(TOP_N_TRAITS).iloc[::-1].copy()
    colors = [
        INCREASED_BLUE if value >= 0 else DECREASED_RED
        for value in top_traits["mean_increased_minus_decreased"]
    ]
    labels = [trait_label(label, 38) for label in top_traits["trait_label"]]
    bars = ax.barh(labels, top_traits["minus_log10_p"], color=colors, edgecolor="white")
    ax.set_title("A. Top Profile-Differentiating Traits", fontsize=14, fontweight="bold", pad=8)
    ax.set_xlabel("-log10(p) from paired increased-minus-decreased profile contrast")
    ax.grid(axis="x", color=SLATE_200, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=9)
    ax.tick_params(axis="x", labelsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(SLATE_200)
    ax.legend(
        handles=[
            plt.Rectangle(
                (0, 0), 1, 1, color=INCREASED_BLUE, label="Enriched in increased profiles"
            ),
            plt.Rectangle(
                (0, 0), 1, 1, color=DECREASED_RED, label="Enriched in decreased profiles"
            ),
        ],
        loc="lower right",
        frameon=False,
        fontsize=9,
    )

    x_offset = max(top_traits["minus_log10_p"]) * 0.015
    for bar, row in zip(bars, top_traits.itertuples(index=False), strict=True):
        label = f"mean {row.mean_increased_minus_decreased:+.2f}"
        ax.text(
            bar.get_width() + x_offset,
            bar.get_y() + bar.get_height() / 2,
            label,
            va="center",
            ha="left",
            fontsize=8.2,
            color=SLATE_700,
        )


def plot_condition_summary(ax: plt.Axes, condition_summary: pd.DataFrame) -> None:
    """Plot panel B: condition groups with strongest directional contrast load."""
    plot_data = condition_summary.head(TOP_N_CONDITIONS).iloc[::-1].copy()
    labels = [compact_label(label, 40) for label in plot_data["condition"]]
    decreased = plot_data["decreased_profile_higher_load"]
    increased = plot_data["increased_profile_higher_load"]
    ax.barh(
        labels,
        decreased,
        color=DECREASED_RED,
        alpha=0.82,
        label="Traits enriched in decreased profiles",
    )
    ax.barh(
        labels,
        increased,
        left=decreased,
        color=INCREASED_BLUE,
        alpha=0.82,
        label="Traits enriched in increased profiles",
    )
    ax.set_title(
        "B. Avg Trait-Profile Contrast per Condition", fontsize=14, fontweight="bold", pad=8
    )
    ax.set_xlabel("Mean absolute paired trait-profile contrast")
    ax.grid(axis="x", color=SLATE_200, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=8.8)
    ax.tick_params(axis="x", labelsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(SLATE_200)
    ax.legend(
        loc="lower right",
        bbox_to_anchor=(1.0, 1.03),
        frameon=False,
        fontsize=9,
        ncol=2,
    )

    max_total = float(plot_data["absolute_contrast_load"].max())
    for y_index, row in enumerate(plot_data.itertuples(index=False)):
        ax.text(
            row.absolute_contrast_load + max_total * 0.015,
            y_index,
            f"{row.absolute_contrast_load:.2f}",
            va="center",
            ha="left",
            fontsize=8.2,
            color=SLATE_700,
            fontweight="bold",
        )


def plot_body_site_summary(ax: plt.Axes, body_site_summary: pd.DataFrame) -> None:
    """Plot panel C: anatomical body-site contrast hotspots."""
    plot_sites = body_site_summary.head(TOP_N_BODY_SITES).copy()
    plot_sites["anchor"] = [anatomical_site_base(site) for site in plot_sites["body_site_group"]]
    plot_sites["anchor_x"] = [point[0] for point in plot_sites["anchor"]]
    plot_sites["anchor_y"] = [point[1] for point in plot_sites["anchor"]]
    plot_sites = plot_sites.sort_values(["anchor_y", "absolute_contrast_load"], ascending=False)

    ax.set_title("C. Anatomical Profile-Contrast Hotspots", fontsize=14, fontweight="bold", pad=8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    add_anatomical_outline(ax)

    max_count = max(float(plot_sites["pair_count"].max()), 1.0)
    max_abs_net = max(float(plot_sites["net_increased_minus_decreased"].abs().max()), 0.05)
    norm = TwoSlopeNorm(vmin=-max_abs_net, vcenter=0, vmax=max_abs_net)
    cmap = plt.get_cmap("RdBu")
    label_ys = np.linspace(0.84, 0.18, len(plot_sites)) if len(plot_sites) > 1 else [0.52]

    for label_y, row in zip(label_ys, plot_sites.itertuples(index=False), strict=False):
        signal = float(row.net_increased_minus_decreased)
        color = cmap(norm(signal))
        size = 130 + 950 * math.sqrt(float(row.pair_count) / max_count)
        ax.scatter(
            [row.anchor_x],
            [row.anchor_y],
            s=size,
            c=[color],
            edgecolors=SLATE_900,
            linewidths=0.8,
            zorder=4,
        )
        ax.plot(
            [row.anchor_x, 0.55, 0.60],
            [row.anchor_y, label_y, label_y],
            color=SLATE_400,
            lw=0.9,
            zorder=2,
        )
        ax.text(
            0.61,
            label_y + 0.018,
            compact_label(row.body_site_group, 28),
            ha="left",
            va="bottom",
            fontsize=9.2,
            fontweight="bold",
            color=SLATE_900,
        )
        ax.text(
            0.61,
            label_y - 0.018,
            (
                f"{int(row.pair_count):,} pairs; load "
                f"{float(row.absolute_contrast_load):.2f}; net {signal:+.2f}"
            ),
            ha="left",
            va="top",
            fontsize=8.0,
            color=SLATE_700,
        )

    ax.text(
        0.02,
        0.04,
        "Bubble size = matched pairs; color = net increased-minus-decreased trait contrast.",
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=SLATE_500,
    )


def format_q(value: float) -> str:
    """Format q-values compactly for annotation text."""
    if not np.isfinite(value):
        return "NA"
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}"


def plot_story_box(
    ax: plt.Axes,
    trait_stats: pd.DataFrame,
    condition_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    *,
    pair_count: int,
    study_count: int,
    trait_count: int,
) -> None:
    """Plot panel D: interpretation text connecting this to enrichment analysis."""
    ax.set_title("D. How it Fits the Story", fontsize=14, fontweight="bold", pad=8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box = FancyBboxPatch(
        (0.02, 0.02),
        0.96,
        0.90,
        boxstyle="round,pad=0.025,rounding_size=0.025",
        linewidth=1.0,
        edgecolor=SLATE_400,
        facecolor=SLATE_50,
    )
    ax.add_patch(box)

    top_increased = trait_stats[trait_stats["mean_increased_minus_decreased"] > 0].iloc[0]
    top_decreased = trait_stats[trait_stats["mean_increased_minus_decreased"] < 0].iloc[0]
    top_condition = condition_summary.iloc[0]
    top_body_site = body_site_summary.iloc[0]
    bullets = [
        (
            "Scope",
            (
                f"{trait_count:,} selected profile traits across {pair_count:,} matched "
                f"study/experiment increased-minus-decreased pairs from {study_count:,} studies."
            ),
        ),
        (
            "Strongest split",
            (
                f"{trait_label(top_increased.trait_label, 42)} is enriched in inc "
                f"profiles (mean {top_increased.mean_increased_minus_decreased:+.2f}, "
                f"q={format_q(float(top_increased.q_value))}); "
                f"{trait_label(top_decreased.trait_label, 42)} is enriched in dec "
                f"profiles (mean {top_decreased.mean_increased_minus_decreased:+.2f}, "
                f"q={format_q(float(top_decreased.q_value))})."
            ),
        ),
        (
            "Where it concentrates",
            (
                f"Highest condition burden: {compact_label(top_condition.condition, 42)}. "
                "Highest explicit body-site load: "
                f"{compact_label(top_body_site.body_site_group, 32)}."
            ),
        ),
        (
            "Relationship to Ilyas",
            (
                "His enrichment infographic ranks common enriched traits across profiles; "
                "this ranks matched directional profile contrasts. Differences are expected "
                "and help separate prevalence from direction-specific structure."
            ),
        ),
        (
            "Caveat",
            (
                "P-values and q-values are descriptive screens over selected profile traits, "
                "not new differential abundance tests."
            ),
        ),
    ]

    ax.text(
        0.06,
        0.86,
        "PROFILE CONTRAST SUMMARY",
        fontsize=11.5,
        fontweight="bold",
        color=SLATE_900,
        va="top",
    )
    y = 0.78
    for heading, body in bullets:
        wrapped = fill(f"{heading}: {body}", width=76)
        ax.text(
            0.06,
            y,
            f"- {wrapped}",
            fontsize=9.3,
            color=SLATE_700,
            va="top",
            linespacing=1.18,
        )
        y -= 0.042 + 0.038 * (wrapped.count("\n") + 1)


def create_infographic(
    trait_stats: pd.DataFrame,
    condition_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    *,
    pair_count: int,
    study_count: int,
    trait_count: int,
) -> None:
    """Render and save the four-panel paired trait-profile contrast infographic."""
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 14,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(18.5, 12.5), dpi=180, facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        left=0.07,
        right=0.97,
        top=0.84,
        bottom=0.08,
        hspace=0.42,
        wspace=0.38,
    )
    ax_traits = fig.add_subplot(grid[0, 0])
    ax_conditions = fig.add_subplot(grid[0, 1])
    ax_body = fig.add_subplot(grid[1, 0])
    ax_story = fig.add_subplot(grid[1, 1])

    fig.suptitle(
        "CultureBotAI: Trait-Profile Contrast Infographic",
        y=0.965,
        fontsize=24,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.5,
        0.925,
        (
            "Analog to enrichment infographic using matched increased-minus-decreased "
            "profile trait z-score contrasts"
        ),
        ha="center",
        va="center",
        fontsize=11.5,
        color=SLATE_500,
    )

    plot_top_traits(ax_traits, trait_stats)
    plot_condition_summary(ax_conditions, condition_summary)
    plot_body_site_summary(ax_body, body_site_summary)
    plot_story_box(
        ax_story,
        trait_stats,
        condition_summary,
        body_site_summary,
        pair_count=pair_count,
        study_count=study_count,
        trait_count=trait_count,
    )
    fig.text(
        0.07,
        0.035,
        (
            "Source: study_trait_profiles_wide.tsv and paired enrichment metadata. Rows are "
            "matched BugSigDB study/experiment increased-minus-decreased profile pairs."
        ),
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=SLATE_500,
    )
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    _paired_matrix, metadata, _paired_trait_metadata = load_paired_inputs()
    profile_matrix = load_profile_matrix()
    trait_labels = load_trait_labels()
    pooled_profile_ids = [
        profile_id
        for profile_id in profile_ids_for_pairs(metadata)
        if profile_id in profile_matrix.index
    ]
    pooled_profile_matrix = profile_matrix.loc[pooled_profile_ids]
    pooled_trait_metadata = select_stratum_traits(pooled_profile_matrix, trait_labels)
    pooled_profile_zscores = build_zscore_matrix(pooled_profile_matrix, pooled_trait_metadata.index)
    matrix = build_pair_matrix_from_profiles(pooled_profile_zscores, metadata)
    pooled_annotation_counts = annotation_counts_for_traits(
        pooled_profile_matrix,
        pooled_trait_metadata.index,
    )
    trait_stats = compute_trait_stats(
        matrix,
        pooled_trait_metadata,
        pooled_annotation_counts,
        stratum_kind="pooled",
        stratum="pooled",
        selected_profile_count=len(pooled_profile_matrix),
    )
    scored = add_pair_scores(matrix, metadata)
    condition_summary = summarize_groups(
        scored,
        "condition",
        min_pairs=MIN_CONDITION_PAIRS,
    )
    body_site_summary = summarize_groups(
        scored,
        "body_site_group",
        min_pairs=MIN_BODY_SITE_PAIRS,
    )

    INVESTIGATION_DIR.mkdir(parents=True, exist_ok=True)
    trait_stats.to_csv(TRAIT_STATS_TSV, sep="\t", index=False)
    trait_stats.to_csv(POOLED_TRAIT_STATS_TSV, sep="\t", index=False)
    condition_summary.to_csv(CONDITION_SUMMARY_TSV, sep="\t", index=False)
    body_site_summary.to_csv(BODY_SITE_SUMMARY_TSV, sep="\t", index=False)
    stratified_manifest = write_stratified_trait_stats(
        profile_matrix,
        metadata,
        trait_labels,
    )

    create_infographic(
        trait_stats,
        condition_summary,
        body_site_summary,
        pair_count=len(matrix),
        study_count=metadata["study_id"].nunique(),
        trait_count=matrix.shape[1],
    )
    print(f"Wrote: {OUTPUT_PNG}")
    print(f"Wrote: {TRAIT_STATS_TSV}")
    print(f"Wrote: {POOLED_TRAIT_STATS_TSV}")
    print(f"Wrote: {CONDITION_SUMMARY_TSV}")
    print(f"Wrote: {BODY_SITE_SUMMARY_TSV}")
    if not stratified_manifest.empty:
        for path in stratified_manifest["path"]:
            print(f"Wrote: {path}")
    print(f"Wrote: {STRATIFIED_TRAIT_STATS_MANIFEST_TSV}")


if __name__ == "__main__":
    main()
