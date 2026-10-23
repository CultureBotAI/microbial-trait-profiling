"""Create an infographic summary of paired trait-profile contrasts."""

from __future__ import annotations

import argparse
import json
import math
import re
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path
from textwrap import fill, shorten

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle, FancyBboxPatch, PathPatch
from matplotlib.path import Path as MplPath
from scipy import stats

from mtp.body_site_diagram import (
    anatomical_site_base,
    bubble_radius,
    leader_line,
    leader_order_key,
)
from mtp.composition_denominators import (
    DENOMINATOR_DEFINITION,
    retained_trait_counts,
    trait_fractions,
)
from mtp.paths import file_sha256
from mtp.trait_guilds import load_profiles

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
# Page and 2 x 2 panel grid of the infographic; panel C is grid[1, 0].
FIGURE_SIZE = (18.5, 12.5)
FIGURE_DPI = 180
PANEL_GRID = {
    "left": 0.07,
    "right": 0.97,
    "top": 0.84,
    "bottom": 0.08,
    "hspace": 0.42,
    "wspace": 0.38,
}
# Panel C's leader lines bend at LEADER_ELBOW_X and end at LEADER_END_X.
LEADER_ELBOW_X = 0.55
LEADER_END_X = 0.60

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
    return retained_trait_counts(values.groupby(level=0, sort=False).sum())


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


def core_pair_metadata(profiles_path: Path, experiments_path: Path) -> pd.DataFrame:
    """Derive the historical human paired cohort from actual core input rows.

    Zero-total profiles remain included, matching the native pooled estimand.
    Neither a heatmap selection nor an old paired-metadata file defines this
    cohort. Repeated direction rows refer to the same summed raw profile.
    """
    profiles, corpus = load_profiles(
        profiles_path, experiments_path, {}, reject_missing_counts=True
    )
    if not profiles:
        raise ValueError("No human paired profiles in the verified core inputs")
    groups: dict[tuple[str, str], dict[str, object]] = {}
    for profile in profiles:
        key = (profile.study_id, profile.experiment_id)
        if key not in groups:
            groups[key] = {
                "pair_id": f"{profile.study_id}|{profile.experiment_id}|paired_inc_minus_dec",
                "study_id": profile.study_id,
                "experiment_id": profile.experiment_id,
                "condition": profile.condition,
                "body_site": profile.body_site,
                "host_species": "Homo sapiens",
                "zero_total_policy": "included_in_historical_native_composition",
            }
        direction = "increased" if profile.inc_dec == "inc" else "decreased"
        groups[key][f"{direction}_profile_ids"] = (
            f"{profile.study_id}|{profile.experiment_id}|{profile.inc_dec}"
        )
        groups[key][f"{direction}_profile_count"] = 1
    metadata = pd.DataFrame.from_records(list(groups.values())).set_index("pair_id")
    for name, value in (
        ("input_profile_rows", corpus.profile_rows),
        ("nonhuman_profile_rows_excluded", corpus.dropped_non_human),
        ("unpaired_profile_rows_excluded", corpus.dropped_unpaired),
        ("missing_metadata_profile_rows_excluded", corpus.unmatched_experiments),
    ):
        metadata[name] = value
    return metadata


def pair_denominator_summary(
    profile_matrix: pd.DataFrame, pair_metadata: pd.DataFrame
) -> dict[str, object]:
    """Report historical inclusion and nonzero-pair eligibility from raw retained counts.

    A pair is nonzero only if every referenced profile has a positive retained
    annotation total, matching ``drop_zero_total_pairs``. Historical native
    composition includes zero-total profiles in standardization and inference;
    this audit does not change that estimand or interpret missing signal as
    measured absence. Totals use every retained trait, before trait selection.
    """
    retained = retained_trait_counts(profile_matrix)
    profile_ids = profile_ids_for_pairs(pair_metadata)
    missing = pd.Index(profile_ids).difference(retained.index)
    if len(missing):
        raise ValueError(f"Paired metadata references missing raw profiles: {list(missing[:3])}")
    if not retained.index.is_unique or not pair_metadata.index.is_unique:
        raise ValueError("Profile and pair identifiers must be unique for denominator reporting")
    totals = retained.loc[profile_ids].sum(axis=1)
    nonzero_pairs: list[object] = []
    for pair_id, row in pair_metadata.iterrows():
        inc_ids = split_profile_ids(row["increased_profile_ids"])
        dec_ids = split_profile_ids(row["decreased_profile_ids"])
        if not inc_ids or not dec_ids:
            raise ValueError(f"Pair {pair_id} lacks an increased or decreased profile")
        if totals.loc[inc_ids + dec_ids].gt(0).all():
            nonzero_pairs.append(pair_id)
    study_ids = (
        pair_metadata["study_id"]
        if "study_id" in pair_metadata.columns
        else pd.Series(pair_metadata.index, index=pair_metadata.index)
    )
    n_nonzero_studies = int(study_ids.loc[nonzero_pairs].nunique())
    return {
        "n_experiments": len(pair_metadata),
        "n_nonzero_experiments": len(nonzero_pairs),
        "n_zero_total_experiments": len(pair_metadata) - len(nonzero_pairs),
        "n_zero_total_excluded": 0,
        "n_nonzero_studies": n_nonzero_studies,
        "n_studies_without_nonzero_experiments": int(study_ids.nunique()) - n_nonzero_studies,
        "n_zero_total_profiles": int(totals.le(0).sum()),
        "retained_trait_count": len(retained.columns),
        "zero_total_policy": "included_in_historical_native_composition",
        "denominator_definition": DENOMINATOR_DEFINITION,
    }


def select_stratum_traits(
    profile_matrix: pd.DataFrame,
    trait_labels: pd.DataFrame,
    *,
    max_traits: int | None = STRATIFIED_TRAIT_COUNT,
    min_profile_support: int = MIN_STRATIFIED_PROFILE_SUPPORT,
    max_profile_support_fraction: float = MAX_STRATIFIED_PROFILE_SUPPORT_FRACTION,
) -> pd.DataFrame:
    """Select supported variable traits within a stratum's profile-count matrix."""
    profile_matrix = retained_trait_counts(profile_matrix)
    composition = trait_fractions(profile_matrix)
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
    composition = trait_fractions(profile_matrix).loc[:, selected_traits]
    scaled = (composition * 1000).apply(lambda column: column.map(math.log1p))
    zscores = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0, ddof=0).replace(0, float("nan"))
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
    study_ids: pd.Series | Sequence[object] | None = None,
    stratum_kind: str = "pooled",
    stratum: str = "pooled",
    selected_profile_count: int | None = None,
    denominator_summary: dict[str, object] | None = None,
) -> pd.DataFrame:
    """Calculate paired and study-averaged one-sample t-test summaries for each trait column."""
    if trait_metadata is None:
        trait_metadata = metadata_from_trait_ids(pd.Index(matrix.columns))
    trait_metadata = trait_metadata.reindex(matrix.columns)
    if annotation_counts is None:
        annotation_counts = pd.Series(pd.NA, index=matrix.columns)
    annotation_counts = annotation_counts.reindex(matrix.columns)

    if study_ids is None:
        study_series = pd.Series(matrix.index, index=matrix.index)
    elif isinstance(study_ids, pd.Series):
        if study_ids.index.equals(matrix.index):
            study_series = study_ids.copy()
        elif matrix.index.isin(study_ids.index).all():
            study_series = study_ids.loc[matrix.index].copy()
        else:
            study_series = study_ids.reindex(matrix.index).copy()
    else:
        study_series = pd.Series(list(study_ids), index=matrix.index)

    study_series = study_series.fillna(pd.Series(matrix.index, index=matrix.index))
    study_matrix = matrix.groupby(study_series, sort=False).mean()
    study_count = len(study_matrix)
    denominator_columns = [
        "n_experiments",
        "n_nonzero_experiments",
        "n_zero_total_experiments",
        "n_zero_total_excluded",
        "n_nonzero_studies",
        "n_studies_without_nonzero_experiments",
        "n_zero_total_profiles",
        "retained_trait_count",
        "zero_total_policy",
        "denominator_definition",
    ]
    denominator_fields = {
        column: (denominator_summary or {}).get(column, pd.NA) for column in denominator_columns
    }
    if denominator_summary is not None and denominator_summary["n_experiments"] != len(matrix):
        raise ValueError("Denominator audit and analyzed pair matrix must describe the same pairs")

    records: list[dict[str, object]] = []
    for column in matrix.columns:
        values = matrix[column].to_numpy(dtype=float)
        pair_n = len(values)
        mean_contrast = float(values.mean())
        median_contrast = float(np.median(values))
        std_contrast = float(values.std(ddof=1)) if pair_n >= 2 else 0.0

        if pair_n < 2 or np.isclose(std_contrast, 0.0):
            pair_test_status = (
                "insufficient_support"
                if pair_n < 2
                else (
                    "all_tied"
                    if std_contrast == 0 and mean_contrast == 0
                    else "insufficient_variance"
                )
            )
            t_stat = float("nan")
            p_value = float("nan")
            ci95_lower = float("nan")
            ci95_upper = float("nan")
        else:
            pair_test_status = "tested"
            result = stats.ttest_1samp(values, 0.0)
            t_stat = float(result.statistic)
            p_value = float(result.pvalue)
            se = std_contrast / np.sqrt(pair_n)
            t_crit = float(stats.t.ppf(0.975, pair_n - 1))
            ci95_lower = mean_contrast - t_crit * se
            ci95_upper = mean_contrast + t_crit * se

        study_values = study_matrix[column].to_numpy(dtype=float)
        study_n = len(study_values)
        study_mean = float(study_values.mean())
        study_median = float(np.median(study_values))
        study_std = float(study_values.std(ddof=1)) if study_n >= 2 else 0.0

        if study_n < 2 or np.isclose(study_std, 0.0):
            study_test_status = (
                "insufficient_support"
                if study_n < 2
                else ("all_tied" if study_std == 0 and study_mean == 0 else "insufficient_variance")
            )
            study_t_stat = float("nan")
            study_p_value = float("nan")
            study_ci95_lower = float("nan")
            study_ci95_upper = float("nan")
        else:
            study_test_status = "tested"
            study_result = stats.ttest_1samp(study_values, 0.0)
            study_t_stat = float(study_result.statistic)
            study_p_value = float(study_result.pvalue)
            study_se = study_std / np.sqrt(study_n)
            study_t_crit = float(stats.t.ppf(0.975, study_n - 1))
            study_ci95_lower = study_mean - study_t_crit * study_se
            study_ci95_upper = study_mean + study_t_crit * study_se

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
                "pair_count": pair_n,
                **denominator_fields,
                "test_method": "two_sided_one_sample_t_test",
                "pair_test_status": pair_test_status,
                "study_test_status": study_test_status,
                "correction_method": "BH",
                "correction_family_size": len(matrix.columns),
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
                "median_increased_minus_decreased": median_contrast,
                "std_increased_minus_decreased": std_contrast,
                "ci95_lower": ci95_lower,
                "ci95_upper": ci95_upper,
                "positive_pair_fraction": float(np.mean(values > 0)),
                "negative_pair_fraction": float(np.mean(values < 0)),
                "nonzero_pair_fraction": float(np.mean(values != 0)),
                "t_statistic": t_stat,
                "p_value": p_value,
                "study_count": study_count,
                "study_mean_increased_minus_decreased": study_mean,
                "study_median_increased_minus_decreased": study_median,
                "study_std_increased_minus_decreased": study_std,
                "study_ci95_lower": study_ci95_lower,
                "study_ci95_upper": study_ci95_upper,
                "study_t_statistic": study_t_stat,
                "study_p_value": study_p_value,
                "direction": (
                    "trait enriched in increased profiles"
                    if study_mean >= 0
                    else "trait enriched in decreased profiles"
                ),
            }
        )

    stats_df = pd.DataFrame.from_records(records)
    stats_df["q_value"] = benjamini_hochberg(stats_df["p_value"].fillna(1.0))
    clipped_p = stats_df["p_value"].clip(lower=np.finfo(float).tiny, upper=1.0)
    stats_df["minus_log10_p"] = -np.log10(clipped_p)
    stats_df["abs_t_statistic"] = stats_df["t_statistic"].abs()

    stats_df["study_q_value"] = benjamini_hochberg(stats_df["study_p_value"].fillna(1.0))
    clipped_study_p = stats_df["study_p_value"].clip(lower=np.finfo(float).tiny, upper=1.0)
    stats_df["study_minus_log10_p"] = -np.log10(clipped_study_p)
    stats_df["study_abs_t_statistic"] = stats_df["study_t_statistic"].abs()

    ordered_columns = [
        "trait",
        "trait_label",
        "stratum_kind",
        "stratum",
        "selected_trait_rank",
        "pair_count",
        *denominator_columns,
        "test_method",
        "pair_test_status",
        "study_test_status",
        "correction_method",
        "correction_family_size",
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
        "ci95_lower",
        "ci95_upper",
        "positive_pair_fraction",
        "negative_pair_fraction",
        "nonzero_pair_fraction",
        "t_statistic",
        "p_value",
        "q_value",
        "minus_log10_p",
        "abs_t_statistic",
        "study_count",
        "study_mean_increased_minus_decreased",
        "study_median_increased_minus_decreased",
        "study_std_increased_minus_decreased",
        "study_ci95_lower",
        "study_ci95_upper",
        "study_t_statistic",
        "study_p_value",
        "study_q_value",
        "study_minus_log10_p",
        "study_abs_t_statistic",
        "direction",
    ]
    return stats_df[ordered_columns].sort_values(
        ["study_p_value", "study_abs_t_statistic", "p_value", "abs_t_statistic", "trait_label"],
        ascending=[True, False, True, False, True],
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
    """Return a repo-relative path or an external basename, without private prefixes."""
    try:
        return str(path.resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return path.name


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
    denominator_summary = pair_denominator_summary(profile_matrix, pair_metadata)
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
    study_ids = (
        pair_metadata.loc[pair_matrix.index, "study_id"]
        if "study_id" in pair_metadata.columns
        else None
    )
    return compute_trait_stats(
        pair_matrix,
        trait_metadata,
        annotation_counts,
        stratum_kind=stratum_kind,
        stratum=stratum,
        selected_profile_count=len(selected_profile_matrix),
        study_ids=study_ids,
        denominator_summary=denominator_summary,
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

    manifest = pd.DataFrame.from_records(
        manifest_rows,
        columns=[
            "stratum_kind",
            "stratum",
            "filename_token",
            "analysis_role",
            "pair_count",
            "study_count",
            "selected_trait_count",
            "path",
        ],
    )
    manifest.to_csv(manifest_path(output_dir), sep="\t", index=False)
    return manifest


def write_legacy_stratified_stats(
    profiles_path: Path,
    pair_metadata_path: Path,
    trait_labels_path: Path,
    output_dir: Path,
    *,
    min_pairs: int = MIN_STRATIFIED_PAIRS,
) -> Path:
    """Refresh pooled and site-group tables with explicit file-level provenance (#824).

    Pair metadata must come from the heatmap producer on these exact counts and
    labels. This is a legacy input binding, not a verified core execution record.
    """
    paths = {
        "profiles": profiles_path,
        "pair_metadata": pair_metadata_path,
        "trait_labels": trait_labels_path,
    }
    hashes = {name: file_sha256(path) for name, path in paths.items()}
    metadata = pd.read_csv(pair_metadata_path, sep="\t", index_col=0)
    for source, column in (("profiles", "profile_sha256"), ("trait_labels", "trait_labels_sha256")):
        if column not in metadata or not metadata[column].eq(hashes[source]).all():
            raise ValueError(
                f"Rebuild pair metadata with just trait-heatmaps: missing or mismatched {column}"
            )
    if metadata.empty:
        raise ValueError("Rebuild pair metadata: no matched pairs")
    if metadata.duplicated(["study_id", "experiment_id"]).any():
        raise ValueError("Rebuild pair metadata: duplicate study/experiment pairs")
    heatmap_selection = {}
    for column in metadata.columns:
        if column.startswith("heatmap_parameter_"):
            values = metadata[column].drop_duplicates().tolist()
            if len(values) != 1 or pd.isna(values[0]):
                raise ValueError(f"Rebuild pair metadata: inconsistent {column}")
            heatmap_selection[column.removeprefix("heatmap_parameter_")] = values[0]
    for row in metadata.itertuples():
        for direction, column in (
            ("inc", "increased_profile_ids"),
            ("dec", "decreased_profile_ids"),
        ):
            expected = f"{row.study_id}|{row.experiment_id}|{direction}"
            if split_profile_ids(getattr(row, column)) != [expected]:
                raise ValueError(
                    f"Pair metadata contains mismatched profile identities: {row.Index}"
                )
    profiles = load_profile_matrix(profiles_path)
    labels = load_trait_labels(trait_labels_path)
    audit = pair_denominator_summary(profiles, metadata)
    previous_tables = set(output_dir.glob("trait_profile_paired_contrast_trait_stats*.tsv"))
    provenance = output_dir / "stratified_stats_provenance.json"
    for candidate in [*previous_tables, provenance]:
        if candidate.exists() and any(candidate.samefile(path) for path in paths.values()):
            raise ValueError("Site-statistics output overlaps an input file")
    # An interrupted/failed refresh must not leave a prior success record in place.
    provenance.unlink(missing_ok=True)
    pooled = compute_stratum_trait_stats(
        profiles, metadata, labels, stratum_kind="pooled", stratum="pooled"
    )
    for source, path in paths.items():
        prefix = "profile" if source == "profiles" else source
        pooled[f"{prefix}_source"] = display_path(path)
        pooled[f"{prefix}_sha256"] = hashes[source]
    output_dir.mkdir(parents=True, exist_ok=True)
    pooled_paths = [
        trait_stats_path(output_dir, "pooled"),
        output_dir / "trait_profile_paired_contrast_trait_stats.tsv",
    ]
    for path in pooled_paths:
        pooled.to_csv(path, sep="\t", index=False)
    manifest = write_stratified_trait_stats(
        profiles,
        metadata,
        labels,
        output_dir=output_dir,
        min_pairs=min_pairs,
    )
    outputs = [*pooled_paths, manifest_path(output_dir)] + [
        trait_stats_path(output_dir, token) for token in manifest["filename_token"]
    ]
    if hashes != {name: file_sha256(path) for name, path in paths.items()}:
        raise ValueError("Inputs changed during stratified export; rerun before using outputs")

    def inventory(path: Path) -> dict[str, object]:
        return {
            "path": display_path(path),
            "sha256": file_sha256(path),
            "size_bytes": path.stat().st_size,
        }

    record = {
        "schema_version": "1.1",
        "analysis": "legacy_native_pooled_and_body_site_group_composition",
        "verification_status": "legacy_inputs_hashed_not_core_verified",
        "created_at": datetime.now(UTC).isoformat(),
        "producer": inventory(Path(__file__)),
        "inputs": {name: inventory(path) for name, path in paths.items()},
        "parameters": {
            "min_pairs": min_pairs,
            "max_traits": STRATIFIED_TRAIT_COUNT,
            "min_profile_support": MIN_STRATIFIED_PROFILE_SUPPORT,
            "max_profile_support_fraction": MAX_STRATIFIED_PROFILE_SUPPORT_FRACTION,
            "site_grouping": "body_site_group from source-bound heatmap pair metadata",
            "pooled_population": "all matched pairs in source-bound heatmap metadata",
            "heatmap_selection": heatmap_selection,
            "heatmap_selection_status": "recorded" if heatmap_selection else "unavailable",
            "strain_rollup_status": "unverified; awaiting upstream scope contract (#847)",
            "multiple_testing": "BH separately within each selected stratum family",
        },
        "pair_denominator_audit": audit,
        "outputs": [inventory(path) for path in outputs],
        "not_regenerated": [
            {**inventory(path), "reason": "outside this run's emitted pooled/site tables"}
            for path in sorted(previous_tables - set(outputs))
        ],
    }
    provenance.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return provenance


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
    score_col = (
        "study_minus_log10_p" if "study_minus_log10_p" in top_traits.columns else "minus_log10_p"
    )
    mean_col = (
        "study_mean_increased_minus_decreased"
        if "study_mean_increased_minus_decreased" in top_traits.columns
        else "mean_increased_minus_decreased"
    )
    ci_low_col = "study_ci95_lower" if "study_ci95_lower" in top_traits.columns else "ci95_lower"
    ci_high_col = "study_ci95_upper" if "study_ci95_upper" in top_traits.columns else "ci95_upper"

    colors = [INCREASED_BLUE if value >= 0 else DECREASED_RED for value in top_traits[mean_col]]
    labels = [trait_label(label, 38) for label in top_traits["trait_label"]]
    bars = ax.barh(labels, top_traits[score_col], color=colors, edgecolor="white")
    ax.set_title("A. Top Profile-Differentiating Traits", fontsize=14, fontweight="bold", pad=22)
    if score_col == "study_minus_log10_p":
        ax.set_xlabel(
            "-log10(study-level p) from paired increased-minus-decreased profile contrast"
        )
    else:
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

    max_p = float(top_traits[score_col].max()) if not top_traits.empty else 1.0
    x_offset = max_p * 0.015
    for bar, row in zip(bars, top_traits.itertuples(index=False), strict=True):
        mean_val = getattr(row, mean_col)
        ci_low = getattr(row, ci_low_col, float("nan"))
        ci_high = getattr(row, ci_high_col, float("nan"))
        if pd.notna(ci_low) and pd.notna(ci_high):
            label = f"mean {mean_val:+.2f} [{ci_low:+.2f}, {ci_high:+.2f}]"
        else:
            label = f"mean {mean_val:+.2f}"
        ax.text(
            bar.get_width() + x_offset,
            bar.get_y() + bar.get_height() / 2,
            label,
            va="center",
            ha="left",
            fontsize=8.2,
            color=SLATE_700,
            zorder=4,
            bbox={
                "boxstyle": "round,pad=0.15",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.9,
            },
        )
    ax.set_xlim(0, max_p * 1.38)


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
        "B. Avg Trait-Profile Contrast per Condition", fontsize=14, fontweight="bold", pad=22
    )
    ax.set_xlabel("Mean absolute paired trait-profile contrast")
    ax.grid(axis="x", color=SLATE_200, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=8.8)
    ax.tick_params(axis="x", labelsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(SLATE_200)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        frameon=False,
        fontsize=8.8,
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
    plot_sites["label_order"] = [
        leader_order_key(point, LEADER_ELBOW_X) for point in plot_sites["anchor"]
    ]
    plot_sites = plot_sites.sort_values(["label_order", "absolute_contrast_load"], ascending=False)

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
    # The panel is a quarter of the page, so the labels use its full height and sit
    # close to their leader lines to stay apart (#620).
    label_ys = np.linspace(0.90, 0.10, len(plot_sites)) if len(plot_sites) > 1 else [0.52]

    for label_y, row in zip(label_ys, plot_sites.itertuples(index=False), strict=False):
        signal = float(row.net_increased_minus_decreased)
        color = cmap(norm(signal))
        radius = bubble_radius(float(row.pair_count), max_count)
        ax.add_patch(
            Circle(
                (row.anchor_x, row.anchor_y),
                radius,
                facecolor=color,
                edgecolor=SLATE_900,
                linewidth=0.8,
                zorder=4,
            )
        )
        ax.plot(
            *leader_line(
                (row.anchor_x, row.anchor_y), radius, label_y, LEADER_ELBOW_X, LEADER_END_X
            ),
            color=SLATE_400,
            lw=0.9,
            zorder=5,
        )
        ax.text(
            0.61,
            label_y + 0.008,
            compact_label(row.body_site_group, 28),
            ha="left",
            va="bottom",
            fontsize=9.2,
            fontweight="bold",
            color=SLATE_900,
        )
        ax.text(
            0.61,
            label_y - 0.008,
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
        0.0,
        "Bubble size = matched pairs; color = net increased-minus-decreased trait contrast.",
        ha="left",
        va="top",
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
        (0.03, 0.03),
        0.94,
        0.89,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=1.0,
        edgecolor=SLATE_400,
        facecolor=SLATE_50,
        clip_on=False,
    )
    ax.add_patch(box)

    mean_col = (
        "study_mean_increased_minus_decreased"
        if "study_mean_increased_minus_decreased" in trait_stats.columns
        else "mean_increased_minus_decreased"
    )
    q_col = "study_q_value" if "study_q_value" in trait_stats.columns else "q_value"
    top_increased = trait_stats[trait_stats[mean_col] > 0].iloc[0]
    top_decreased = trait_stats[trait_stats[mean_col] < 0].iloc[0]
    top_condition = condition_summary.iloc[0]
    top_body_site = body_site_summary.iloc[0]
    inc_mean = getattr(top_increased, mean_col)
    inc_q = getattr(top_increased, q_col)
    dec_mean = getattr(top_decreased, mean_col)
    dec_q = getattr(top_decreased, q_col)
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
                f"profiles (mean {inc_mean:+.2f}, "
                f"study q={format_q(float(inc_q))}); "
                f"{trait_label(top_decreased.trait_label, 42)} is enriched in dec "
                f"profiles (mean {dec_mean:+.2f}, "
                f"study q={format_q(float(dec_q))})."
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
                "P-values and q-values reflect study-averaged paired profile contrasts "
                "over selected traits, serving as descriptive screens rather than new "
                "differential abundance tests."
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
    fig = plt.figure(figsize=FIGURE_SIZE, dpi=FIGURE_DPI, facecolor="white")
    grid = fig.add_gridspec(2, 2, **PANEL_GRID)
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


def main(argv: Sequence[str] | None = None) -> None:
    from mtp.workflow.consumers import (
        bind_count_arguments,
        preflight_consumer_outputs,
        write_consumer_record,
    )

    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--stats-only", action="store_true", help="Write only the pooled stats TSV")
    mode.add_argument(
        "--stratified-stats-only",
        action="store_true",
        help="Refresh legacy pooled and site-group stats, manifest and provenance without figures",
    )
    parser.add_argument(
        "--min-stratified-pairs",
        type=int,
        help=(
            f"Minimum site-group pairs for --stratified-stats-only (default {MIN_STRATIFIED_PAIRS})"
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=INVESTIGATION_DIR)
    parser.add_argument("--core-run", type=Path, help="Verified core run for --stats-only")
    parser.add_argument("--profiles", dest="profile_input", type=Path)
    parser.add_argument("--pair-metadata", type=Path)
    parser.add_argument("--trait-labels", dest="column_labels", type=Path)
    parser.set_defaults(experiments_input=None, full_dump_input=None)
    args = parser.parse_args(argv)
    if (
        not args.stats_only
        and not args.stratified_stats_only
        and args.output_dir != INVESTIGATION_DIR
    ):
        parser.error("--output-dir requires --stats-only or --stratified-stats-only")
    if args.core_run is not None and not args.stats_only:
        parser.error("--core-run requires --stats-only")
    if args.core_run is not None and args.pair_metadata is not None:
        parser.error("--core-run derives pair metadata; omit --pair-metadata")
    if args.min_stratified_pairs is not None:
        if not args.stratified_stats_only:
            parser.error("--min-stratified-pairs requires --stratified-stats-only")
        if args.min_stratified_pairs < 1:
            parser.error("--min-stratified-pairs must be positive")
    else:
        args.min_stratified_pairs = MIN_STRATIFIED_PAIRS
    if args.stratified_stats_only:
        provenance = write_legacy_stratified_stats(
            args.profile_input or PROFILE_MATRIX_PATH,
            args.pair_metadata or PAIRED_METADATA_PATH,
            args.column_labels or TRAIT_COLUMNS_PATH,
            args.output_dir,
            min_pairs=args.min_stratified_pairs,
        )
        print(f"Wrote pooled/site-group statistics and provenance: {provenance}")
        return
    bundle = bind_count_arguments(args)
    output_path = trait_stats_path(args.output_dir, "pooled")
    record_path = args.output_dir / "core_consumption.json"
    if bundle is not None:
        args.pair_metadata = args.output_dir / "trait_profile_paired_contrast_pair_metadata.tsv"
        preflight_consumer_outputs(
            bundle,
            record_path,
            [output_path, args.pair_metadata],
            output_directory=args.output_dir,
        )
        metadata = core_pair_metadata(args.profile_input, args.experiments_input)
        for source_name, source_path in (
            ("profile", args.profile_input),
            ("experiments", args.experiments_input),
            ("trait_labels", args.column_labels),
        ):
            metadata[f"{source_name}_source"] = display_path(source_path)
            metadata[f"{source_name}_sha256"] = file_sha256(source_path)
    else:
        args.profile_input = args.profile_input or PROFILE_MATRIX_PATH
        args.pair_metadata = args.pair_metadata or PAIRED_METADATA_PATH
        args.column_labels = args.column_labels or TRAIT_COLUMNS_PATH
        metadata = pd.read_csv(args.pair_metadata, sep="\t", index_col=0)
    profile_matrix = load_profile_matrix(args.profile_input)
    trait_labels = load_trait_labels(args.column_labels)
    denominator_summary = pair_denominator_summary(profile_matrix, metadata)
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
    study_ids = metadata.loc[matrix.index, "study_id"] if "study_id" in metadata.columns else None
    trait_stats = compute_trait_stats(
        matrix,
        pooled_trait_metadata,
        pooled_annotation_counts,
        stratum_kind="pooled",
        stratum="pooled",
        selected_profile_count=len(pooled_profile_matrix),
        study_ids=study_ids,
        denominator_summary=denominator_summary,
    )
    if bundle is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        metadata.to_csv(args.pair_metadata, sep="\t")
        trait_stats["experiments_source"] = display_path(args.experiments_input)
        trait_stats["experiments_sha256"] = file_sha256(args.experiments_input)
    for source_name, source_path in [
        ("profile", args.profile_input),
        ("pair_metadata", args.pair_metadata),
        ("trait_labels", args.column_labels),
    ]:
        trait_stats[f"{source_name}_source"] = display_path(source_path)
        trait_stats[f"{source_name}_sha256"] = file_sha256(source_path)
    if args.stats_only:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "stratified_stats_provenance.json").unlink(missing_ok=True)
        trait_stats.to_csv(output_path, sep="\t", index=False)
        if bundle is not None:
            write_consumer_record(
                bundle,
                record_path,
                consumer="native-pooled-composition-stats",
                outputs=[output_path, args.pair_metadata],
            )
        print(f"Wrote: {output_path}")
        print(f"Pair denominator audit: {denominator_summary}")
        return
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
    (INVESTIGATION_DIR / "stratified_stats_provenance.json").unlink(missing_ok=True)
    trait_stats.to_csv(TRAIT_STATS_TSV, sep="\t", index=False)
    trait_stats.to_csv(POOLED_TRAIT_STATS_TSV, sep="\t", index=False)
    condition_summary.to_csv(CONDITION_SUMMARY_TSV, sep="\t", index=False)
    body_site_summary.to_csv(BODY_SITE_SUMMARY_TSV, sep="\t", index=False)
    stratified_manifest = write_stratified_trait_stats(
        profile_matrix,
        metadata,
        trait_labels,
        output_dir=INVESTIGATION_DIR,
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
