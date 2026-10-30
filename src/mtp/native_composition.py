"""Shared native composition selection, transformation and support (#976).

These functions reproduce the historical producer's quantities. They do not
evaluate test probabilities or change its all-retained-annotation denominator.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from mtp.composition_denominators import retained_trait_counts, trait_fractions
from mtp.trait_guilds import load_profiles

PROFILE_METADATA_COLUMNS = ["study_id", "experiment_id", "inc_dec"]


def clean_text(value: object) -> str:
    """Preserve the producer's metadata normalization and missing-label spelling."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def full_trait_label(value: object) -> str:
    return clean_text(value).replace("\n", " ")


def format_trait_label(column_id: str, row: pd.Series) -> str:
    name = clean_text(row.get("object_name", "")).replace("_", " ")
    if name == "unknown":
        name = clean_text(row.get("column_label", "")).replace("_", " ")
    short_id = column_id.split("|")[-1]
    if short_id.startswith("pathways:"):
        short_id = short_id.replace("pathways:", "")
    return full_trait_label(f"{name} {short_id}")


def load_trait_labels(path: Path) -> pd.DataFrame:
    labels = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    labels = labels.set_index("column_id", drop=False)
    labels["display_label"] = [
        format_trait_label(str(column_id), row) for column_id, row in labels.iterrows()
    ]
    labels["trait_label"] = labels["display_label"].map(full_trait_label)
    labels["display_category"] = labels["object_category"].str.removeprefix("biolink:")
    labels.loc[labels["display_category"] == "", "display_category"] = "unknown"
    return labels


def load_profile_matrix(path: Path) -> pd.DataFrame:
    profiles = pd.read_csv(path, sep="\t")
    columns = [column for column in profiles.columns if column not in PROFILE_METADATA_COLUMNS]
    values = profiles[columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    values.index = [
        f"{study}|{experiment}|{direction}"
        for study, experiment, direction in profiles[PROFILE_METADATA_COLUMNS].itertuples(
            index=False, name=None
        )
    ]
    return retained_trait_counts(values.groupby(level=0, sort=False).sum())


def split_profile_ids(value: object) -> list[str]:
    return [profile_id for profile_id in clean_text(value).split(";") if profile_id != "unknown"]


def profile_ids_for_pairs(pair_metadata: pd.DataFrame) -> list[str]:
    profile_ids: list[str] = []
    for row in pair_metadata.itertuples(index=False):
        profile_ids.extend(split_profile_ids(row.increased_profile_ids))
        profile_ids.extend(split_profile_ids(row.decreased_profile_ids))
    return list(dict.fromkeys(profile_ids))


def core_pair_metadata(profiles_path: Path, experiments_path: Path) -> pd.DataFrame:
    """The producer's human paired cohort, including historical zero-total pairs."""
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


def select_stratum_traits(
    profile_matrix: pd.DataFrame,
    trait_labels: pd.DataFrame,
    *,
    max_traits: int | None = None,
    min_profile_support: int = 6,
    max_profile_support_fraction: float = 0.95,
) -> pd.DataFrame:
    """The native support/variance selection, with all retained columns in totals."""
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
            trait_metadata["composition_variance"], trait_metadata["total_signal"], strict=True
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
        raise ValueError("Fewer than two traits passed the stratum selection filters.")
    return eligible


def build_zscore_matrix(profile_matrix: pd.DataFrame, selected_traits: pd.Index) -> pd.DataFrame:
    composition = trait_fractions(profile_matrix).loc[:, selected_traits]
    scaled = (composition * 1000).apply(lambda column: column.map(math.log1p))
    zscores: pd.DataFrame = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0, ddof=0).replace(
        0, float("nan")
    )
    return zscores.fillna(0.0)


def build_pair_matrix_from_profiles(
    profile_zscores: pd.DataFrame, pair_metadata: pd.DataFrame
) -> pd.DataFrame:
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
        raise ValueError("No matched increased/decreased pairs were available for this stratum.")
    return pd.DataFrame(rows)


def study_contrasts(
    matrix: pd.DataFrame, study_ids: pd.Series | Sequence[object] | None
) -> pd.DataFrame:
    """Collapse the producer's pair rows to equal-study observations, in native order."""
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
    return matrix.groupby(study_series, sort=False).mean()


@dataclass(frozen=True)
class StudyMoments:
    count: int
    mean: float
    median: float
    std: float
    status: str


def study_moments(values: NDArray[np.float64]) -> StudyMoments:
    """Native ddof=1 dispersion and zero-variance reason, without a t test."""
    count = len(values)
    mean = float(values.mean())
    median = float(np.median(values))
    std = float(values.std(ddof=1)) if count >= 2 else 0.0
    if count < 2:
        status = "insufficient_support"
    elif np.isclose(std, 0.0):
        status = "all_tied" if std == 0 and mean == 0 else "insufficient_variance"
    else:
        status = "tested"
    return StudyMoments(count, mean, median, std, status)
