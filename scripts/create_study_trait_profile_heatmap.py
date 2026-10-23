"""Create a clustered heatmap of BugSigDB study/experiment/direction trait profiles."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list, linkage

from mtp.composition_denominators import (
    DENOMINATOR_DEFINITION,
    retained_trait_counts,
    trait_fractions,
)
from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES
from mtp.paths import file_sha256

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization"
PROFILE_METADATA_COLUMNS = ["study_id", "experiment_id", "inc_dec"]
DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}
# These name the already oriented input; plotting never swaps the data (#842).
DIRECTION_BASIS_LABELS = {
    "group_1_vs_group_0": "Group 1 vs Group 0",
    "disease_vs_comparison": "disease vs comparison",
}
PAIR_DIRECTION_LABEL = "increased - decreased"
HUMAN_HOST_SPECIES = "Homo sapiens"
EXPERIMENT_LABEL_COLUMNS = {
    "condition": "Condition",
    "body_site": "Body site",
    "location": "Location of subjects",
    "host_species": "Host species",
    # Read for the corpus's 16S and shotgun rule (#742)
    "sequencing_type": "Sequencing type",
}
REQUIRED_PAIRED_TRAITS = {
    "METPO:1000615": "mesophilic",
    "METPO:1000698": "gram positive",
    "METPO:1000703": "non motile",
}


def clean_value(value: object) -> str:
    """Normalize empty BugSigDB metadata values."""
    text = str(value or "").strip()
    if not text or text.upper() == "NA" or text.lower() == "nan":
        return "unknown"
    return text


def compact_label(value: str, max_length: int = 42) -> str:
    """Shorten long labels for heatmap axes and legends."""
    value = clean_value(value)
    if len(value) <= max_length:
        return value
    return f"{value[: max_length - 3]}..."


def study_id_from_bsdb_id(bsdb_id: str) -> str:
    """Convert a full BugSigDB signature ID into the derived study ID format."""
    compact_id = bsdb_id.strip().removeprefix("bsdb:").split("/", 1)[0]
    return f"bsdb:{compact_id}"


def experiment_id_from_bsdb_id(bsdb_id: str) -> str:
    """Convert a full BugSigDB signature ID into the derived experiment ID format."""
    compact_id = bsdb_id.strip().removeprefix("bsdb:")
    parts = compact_id.split("/")
    if len(parts) < 2:
        return "unknown"
    return f"bsdb:{parts[0]}-{parts[1]}"


def most_common_label(values: list[str]) -> str:
    """Return the most common non-empty label, with deterministic ties."""
    labels = [clean_value(value) for value in values if clean_value(value) != "unknown"]
    if not labels:
        return "unknown"
    counts = Counter(labels)
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def load_raw_study_metadata(full_dump_path: Path) -> pd.DataFrame:
    """Extract each study's majority condition, body-site and location labels.

    These summarize the whole study and are kept only as `study_majority_*` columns.
    A profile's own labels come from its experiment (`load_raw_experiment_metadata`).
    Joining the majority labels to every experiment gave 1,241 of 5,829 profiles
    (803 experiments, 265 studies) another experiment's condition, body site or
    location; 1,112 of them had the wrong condition or body site (#555).
    """
    values_by_study: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: {
            "condition": [],
            "body_site": [],
            "location": [],
        }
    )
    with full_dump_path.open(encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(rows)
        for row in reader:
            bsdb_id = clean_value(row.get("BSDB ID"))
            if bsdb_id == "unknown":
                continue
            study_id = study_id_from_bsdb_id(bsdb_id)
            values_by_study[study_id]["condition"].append(clean_value(row.get("Condition")))
            values_by_study[study_id]["body_site"].append(clean_value(row.get("Body site")))
            values_by_study[study_id]["location"].append(
                clean_value(row.get("Location of subjects"))
            )

    records = []
    for study_id, metadata in values_by_study.items():
        records.append(
            {
                "study_id": study_id,
                "study_majority_condition": most_common_label(metadata["condition"]),
                "study_majority_body_site": most_common_label(metadata["body_site"]),
                "study_majority_location": most_common_label(metadata["location"]),
                "study_condition_count": len(
                    {value for value in metadata["condition"] if clean_value(value) != "unknown"}
                ),
                "study_body_site_count": len(
                    {value for value in metadata["body_site"] if clean_value(value) != "unknown"}
                ),
            }
        )
    return pd.DataFrame.from_records(records).set_index("study_id")


def load_raw_experiment_metadata(full_dump_path: Path) -> pd.DataFrame:
    """Extract condition, body-site, location and host labels for each experiment.

    BugSigDB curates these labels per experiment, so every signature of an experiment
    carries the same values. An experiment with two different values for one label is
    rejected instead of being resolved by majority.
    """
    values_by_experiment: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(
        lambda: {column: set() for column in EXPERIMENT_LABEL_COLUMNS}
    )
    with full_dump_path.open(encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(rows)
        for row in reader:
            bsdb_id = clean_value(row.get("BSDB ID"))
            if bsdb_id == "unknown":
                continue
            experiment_id = experiment_id_from_bsdb_id(bsdb_id)
            if experiment_id == "unknown":
                continue
            study_id = study_id_from_bsdb_id(bsdb_id)
            labels = values_by_experiment[(study_id, experiment_id)]
            for column, source_column in EXPERIMENT_LABEL_COLUMNS.items():
                value = clean_value(row.get(source_column))
                if value != "unknown":
                    labels[column].add(value)

    records = []
    for (study_id, experiment_id), labels in values_by_experiment.items():
        record = {"study_id": study_id, "experiment_id": experiment_id}
        for column, values in labels.items():
            if len(values) > 1:
                msg = (
                    f"Experiment {experiment_id} has more than one {column} label in "
                    f"{full_dump_path}: {sorted(values)}"
                )
                raise ValueError(msg)
            record[column] = next(iter(values), "unknown")
        records.append(record)
    return pd.DataFrame.from_records(records).set_index(["study_id", "experiment_id"])


def load_study_titles(studies_path: Path) -> pd.DataFrame:
    """Load compact title metadata from derived studies.tsv."""
    studies = pd.read_csv(
        studies_path,
        sep="\t",
        usecols=["study_id", "title", "publication_year", "taxon_count", "signature_count"],
        dtype=str,
    )
    for column in ["title", "publication_year"]:
        studies[column] = studies[column].map(clean_value)
    for column in ["taxon_count", "signature_count"]:
        studies[column] = pd.to_numeric(studies[column], errors="coerce").fillna(0).astype(int)
    return studies.set_index("study_id")


def load_trait_labels(columns_path: Path) -> pd.DataFrame:
    """Load display labels and categories for trait columns."""
    columns = pd.read_csv(columns_path, sep="\t", dtype=str).fillna("")
    columns = columns.set_index("column_id", drop=False)
    columns["display_label"] = [
        format_trait_label(column_id, row) for column_id, row in columns.iterrows()
    ]
    columns["display_category"] = columns["object_category"].str.removeprefix("biolink:")
    columns.loc[columns["display_category"] == "", "display_category"] = "unknown"
    return columns


def format_trait_label(column_id: str, row: pd.Series) -> str:
    """Build a compact two-line trait label."""
    name = row.get("object_name", "") or row.get("column_label", "") or column_id
    name = compact_label(str(name).replace("_", " "), 34)
    short_id = column_id.split("|")[-1].replace("pathways:", "").replace("carbon_substrates:", "")
    short_id = compact_label(short_id, 22)
    return f"{name}\n{short_id}"


def group_top_categories(series: pd.Series, max_categories: int) -> pd.Series:
    """Keep the most common categories and group the rest as Other."""
    cleaned = series.map(clean_value)
    counts = cleaned[cleaned != "unknown"].value_counts()
    keep = set(counts.head(max_categories).index)
    grouped = cleaned.where(cleaned.isin(keep), "Other")
    grouped = grouped.where(cleaned != "unknown", "unknown")
    return grouped


def profile_id_from_parts(study_id: str, experiment_id: str, inc_dec: str) -> str:
    """Build a stable row ID for one study/experiment/direction profile."""
    return f"{study_id}|{experiment_id}|{inc_dec}"


def load_profile_matrix(profile_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Load one row per study, experiment, and increased/decreased trait profile."""
    profiles = pd.read_csv(profile_path, sep="\t")
    trait_columns = [
        column for column in profiles.columns if column not in PROFILE_METADATA_COLUMNS
    ]
    values = profiles[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    profile_ids = [
        profile_id_from_parts(str(study_id), str(experiment_id), str(inc_dec))
        for study_id, experiment_id, inc_dec in profiles[PROFILE_METADATA_COLUMNS].itertuples(
            index=False,
            name=None,
        )
    ]
    values.index = profile_ids
    profile_matrix = retained_trait_counts(values.groupby(level=0, sort=False).sum())
    trait_columns = profile_matrix.columns.tolist()

    metadata = profiles[PROFILE_METADATA_COLUMNS].copy()
    metadata["profile_id"] = profile_ids
    metadata = metadata.drop_duplicates(subset=["profile_id"], keep="first").set_index("profile_id")
    metadata["direction_label"] = (
        metadata["inc_dec"].map(DIRECTION_LABELS).fillna(metadata["inc_dec"])
    )
    metadata["total_trait_signal"] = profile_matrix.sum(axis=1)
    metadata["trait_count"] = profile_matrix.gt(0).sum(axis=1)
    metadata["denominator_definition"] = DENOMINATOR_DEFINITION
    return profile_matrix, metadata, trait_columns


def join_profile_labels(
    profile_metadata: pd.DataFrame,
    studies_path: Path,
    full_dump_path: Path,
) -> pd.DataFrame:
    """Attach study titles, each profile's experiment labels and study-majority labels."""
    profile_metadata = profile_metadata.join(
        load_study_titles(studies_path),
        on="study_id",
        how="left",
    )
    profile_metadata = profile_metadata.join(
        load_raw_experiment_metadata(full_dump_path),
        on=["study_id", "experiment_id"],
        how="left",
    )
    profile_metadata = profile_metadata.join(
        load_raw_study_metadata(full_dump_path),
        on="study_id",
        how="left",
    )
    for column in [
        "title",
        "publication_year",
        *EXPERIMENT_LABEL_COLUMNS,
        "study_majority_condition",
        "study_majority_body_site",
        "study_majority_location",
    ]:
        profile_metadata[column] = profile_metadata[column].map(clean_value)
    return profile_metadata


def filter_human_host_profiles(
    profile_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep only human-host 16S and shotgun profiles before profile or trait selection.

    The canonical table is already restricted; a non-default input is restricted here
    too, so its figures match their titles (#742, #754).
    """
    for column in ("host_species", "sequencing_type"):
        if column not in profile_metadata.columns:
            msg = f"Cannot apply human-host corpus filter: missing {column} metadata."
            raise ValueError(msg)

    host_species = profile_metadata["host_species"].map(clean_value)
    sequencing_type = profile_metadata["sequencing_type"].map(clean_value)
    human_profile_ids = profile_metadata.index[
        (host_species == HUMAN_HOST_SPECIES) & sequencing_type.isin(CORPUS_SEQUENCING_TYPES)
    ]
    aligned_profile_ids = profile_matrix.index[profile_matrix.index.isin(human_profile_ids)]
    filtered_matrix = profile_matrix.loc[aligned_profile_ids].copy()
    filtered_metadata = profile_metadata.loc[aligned_profile_ids].copy()
    filtered_metadata["host_species"] = host_species.loc[aligned_profile_ids]
    if filtered_metadata.empty:
        msg = (
            "No Homo sapiens 16S or WMS profiles were available after applying the "
            "human-host filter."
        )
        raise ValueError(msg)
    return filtered_matrix, filtered_metadata


def select_profiles(
    profile_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
    max_profiles: int,
    min_traits: int,
    min_total_signal: int,
) -> pd.Index:
    """Select profile rows for the heatmap; by default this keeps every profile."""
    eligible = profile_metadata[
        (profile_metadata["trait_count"] >= min_traits)
        & (profile_metadata["total_trait_signal"] >= min_total_signal)
    ].copy()
    if len(eligible) < 2:
        msg = "Fewer than two profiles passed the selection filters."
        raise ValueError(msg)
    eligible = eligible.loc[profile_matrix.index.intersection(eligible.index)]
    if max_profiles <= 0:
        return eligible.index

    eligible = eligible.sort_values(
        ["total_trait_signal", "trait_count"],
        ascending=False,
    )
    return eligible.head(max_profiles).index


def select_traits(
    selected_profile_matrix: pd.DataFrame,
    trait_labels: pd.DataFrame,
    max_traits: int,
    min_profile_support: int,
    max_profile_support_fraction: float,
    include_media_traits: bool,
) -> tuple[pd.Index, pd.DataFrame]:
    """Select traits that vary across the selected study profiles."""
    if include_media_traits:
        msg = "MediaDive traits cannot be included in the harmonized composition arm."
        raise ValueError(msg)
    selected_profile_matrix = retained_trait_counts(selected_profile_matrix)
    composition = trait_fractions(selected_profile_matrix)
    support = selected_profile_matrix.gt(0).sum(axis=0).astype(int)
    total_signal = selected_profile_matrix.sum(axis=0).astype(float)
    variance = composition.var(axis=0)

    trait_metadata = pd.DataFrame(
        {
            "trait_id": composition.columns,
            "profile_support": support,
            "total_signal": total_signal,
            "composition_variance": variance,
        },
        index=composition.columns,
    ).join(
        trait_labels[
            [
                "display_label",
                "column_label",
                "object_name",
                "object_category",
                "display_category",
            ]
        ],
        how="left",
    )
    trait_metadata[["display_label", "display_category"]] = trait_metadata[
        ["display_label", "display_category"]
    ].fillna("")
    trait_metadata["selection_score"] = [
        variance_value * math.log1p(signal)
        for variance_value, signal in zip(
            trait_metadata["composition_variance"],
            trait_metadata["total_signal"],
            strict=True,
        )
    ]

    max_support = max(2, int(len(selected_profile_matrix) * max_profile_support_fraction))
    is_required = trait_metadata.index.isin(REQUIRED_PAIRED_TRAITS)
    eligible = trait_metadata[
        (trait_metadata["profile_support"] >= min_profile_support)
        & ((trait_metadata["profile_support"] <= max_support) | is_required)
        & (trait_metadata["composition_variance"] > 0)
    ].copy()
    if not include_media_traits:
        eligible = eligible[
            ~eligible.index.to_series().str.contains("mediadive.medium", regex=False)
            & ~eligible["display_label"].str.contains("mediadive.medium", regex=False)
        ].copy()
    eligible = eligible.sort_values(
        ["selection_score", "composition_variance", "profile_support"],
        ascending=False,
    )
    missing_required = [
        f"{trait_name} ({trait_id})"
        for trait_id, trait_name in REQUIRED_PAIRED_TRAITS.items()
        if trait_id not in eligible.index
    ]
    if missing_required:
        msg = "Required paired-heatmap trait(s) did not pass selection filters: " + ", ".join(
            missing_required
        )
        raise ValueError(msg)

    required_ids = set(REQUIRED_PAIRED_TRAITS)
    if max_traits < len(required_ids):
        msg = (
            f"--max-traits must be at least {len(required_ids)} to include the "
            "required paired-heatmap trait checks."
        )
        raise ValueError(msg)

    selected_set = set(eligible.head(max_traits).index) | required_ids
    optional_keep = [
        trait_id for trait_id in eligible.index if trait_id in selected_set - required_ids
    ][: max_traits - len(required_ids)]
    selected_ids = [
        trait_id
        for trait_id in eligible.index
        if trait_id in required_ids or trait_id in optional_keep
    ]

    eligible = eligible.loc[selected_ids]
    if len(eligible) < 2:
        msg = "Fewer than two traits passed the selection filters."
        raise ValueError(msg)
    return eligible.index, trait_metadata.loc[eligible.index]


def build_plot_matrix(
    selected_profile_matrix: pd.DataFrame,
    selected_traits: pd.Index,
) -> pd.DataFrame:
    """Create column z-scores from log-scaled trait fractions."""
    composition = trait_fractions(selected_profile_matrix).loc[:, selected_traits]
    scaled = (composition * 1000).apply(lambda column: column.map(math.log1p))
    zscores = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0, ddof=0).replace(0, pd.NA)
    return zscores.fillna(0)


def cluster_index(indexed_matrix: pd.DataFrame) -> pd.Index:
    """Return a clustered row order for a block, falling back to current order."""
    if len(indexed_matrix) < 2:
        return indexed_matrix.index
    values = indexed_matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        return indexed_matrix.index
    try:
        row_linkage = linkage(values, method="average", metric="euclidean")
    except ValueError:
        return indexed_matrix.index
    return indexed_matrix.index[leaves_list(row_linkage)]


def order_profiles_by_direction(
    plot_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
) -> pd.Index:
    """Cluster rows within increased/decreased sections and keep sections separate."""
    ordered: list[str] = []
    for direction in ["inc", "dec"]:
        direction_index = profile_metadata.index[profile_metadata["inc_dec"] == direction]
        direction_index = direction_index.intersection(plot_matrix.index)
        ordered.extend(cluster_index(plot_matrix.loc[direction_index]).tolist())

    remaining = plot_matrix.index.difference(pd.Index(ordered), sort=False)
    ordered.extend(cluster_index(plot_matrix.loc[remaining]).tolist())
    return pd.Index(ordered)


def most_common_series_value(series: pd.Series) -> str:
    """Return the most common category value in a display bin."""
    if series.empty:
        return "unknown"
    counts = series.map(clean_value).value_counts()
    if counts.empty:
        return "unknown"
    return str(counts.sort_values(ascending=False).index[0])


def allocate_display_bins(profile_metadata: pd.DataFrame, max_display_rows: int) -> dict[str, int]:
    """Allocate rendered row bins across direction sections."""
    if max_display_rows <= 0 or len(profile_metadata) <= max_display_rows:
        return {
            direction: int((profile_metadata["inc_dec"] == direction).sum())
            for direction in profile_metadata["inc_dec"].dropna().unique()
        }

    directions = [
        direction
        for direction in ["inc", "dec"]
        if int((profile_metadata["inc_dec"] == direction).sum()) > 0
    ]
    directions.extend(
        direction
        for direction in sorted(set(profile_metadata["inc_dec"]) - set(directions))
        if int((profile_metadata["inc_dec"] == direction).sum()) > 0
    )

    remaining_rows = max_display_rows
    allocations: dict[str, int] = {}
    for index, direction in enumerate(directions):
        direction_count = int((profile_metadata["inc_dec"] == direction).sum())
        directions_left = len(directions) - index
        if directions_left == 1:
            allocation = remaining_rows
        else:
            proportional = round(max_display_rows * direction_count / len(profile_metadata))
            allocation = max(1, min(proportional, remaining_rows - directions_left + 1))
        allocation = min(allocation, direction_count)
        allocations[direction] = allocation
        remaining_rows -= allocation

    return allocations


def build_rendered_display_matrix(
    plot_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    max_display_rows: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate adjacent clustered profiles into display bins for a readable PNG."""
    if max_display_rows <= 0 or len(plot_matrix) <= max_display_rows:
        display_matrix = build_display_matrix(plot_matrix, profile_metadata, trait_metadata)
        display_metadata = profile_metadata.copy()
        display_metadata.index = display_matrix.index
        display_metadata["profile_count"] = 1
        display_metadata["study_count_in_bin"] = 1
        return display_matrix, display_metadata

    rows: list[pd.Series] = []
    metadata_rows: list[dict[str, object]] = []
    allocations = allocate_display_bins(profile_metadata, max_display_rows)

    for direction, bin_count in allocations.items():
        direction_index = profile_metadata.index[profile_metadata["inc_dec"] == direction]
        direction_index = direction_index.intersection(plot_matrix.index)
        if len(direction_index) == 0:
            continue

        for bin_number, bin_profile_ids in enumerate(
            np.array_split(direction_index.to_numpy(), bin_count),
            start=1,
        ):
            if len(bin_profile_ids) == 0:
                continue

            bin_index = pd.Index(bin_profile_ids)
            bin_values = plot_matrix.loc[bin_index].mean(axis=0)
            bin_metadata = profile_metadata.loc[bin_index]
            direction_label = DIRECTION_LABELS.get(direction, direction)
            label = (
                f"{direction_label} profile bin {bin_number:03d} | "
                f"n={len(bin_index)} | studies={bin_metadata['study_id'].nunique()}"
            )
            bin_values.name = label
            rows.append(bin_values)
            metadata_rows.append(
                {
                    "display_label": label,
                    "inc_dec": direction,
                    "direction_label": direction_label,
                    "profile_count": len(bin_index),
                    "study_count_in_bin": bin_metadata["study_id"].nunique(),
                    "body_site_group": most_common_series_value(bin_metadata["body_site_group"]),
                    "condition_group": most_common_series_value(bin_metadata["condition_group"]),
                    "first_profile_id": bin_index[0],
                    "last_profile_id": bin_index[-1],
                }
            )

    display_matrix = pd.DataFrame(rows)
    display_matrix.columns = [
        trait_metadata.loc[trait_id, "display_label"] or trait_id
        for trait_id in display_matrix.columns
    ]
    display_metadata = pd.DataFrame.from_records(metadata_rows).set_index("display_label")
    return display_matrix, display_metadata


def build_display_matrix(
    plot_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Apply readable labels for the plotted matrix."""
    display = plot_matrix.copy()
    display.index = [
        " | ".join(
            [
                profile_metadata.loc[profile_id, "study_id"],
                profile_metadata.loc[profile_id, "experiment_id"],
                profile_metadata.loc[profile_id, "direction_label"],
            ]
        )
        for profile_id in display.index
    ]
    display.columns = [
        trait_metadata.loc[trait_id, "display_label"] or trait_id for trait_id in display.columns
    ]
    return display


def build_paired_enrichment_matrix(
    plot_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build matched increased-minus-decreased rows for each study/experiment pair."""
    rows: list[pd.Series] = []
    metadata_rows: list[dict[str, object]] = []
    group_columns = ["study_id", "experiment_id"]
    for (study_id, experiment_id), group in profile_metadata.groupby(group_columns, sort=False):
        inc_ids = group.index[group["inc_dec"] == "inc"].intersection(plot_matrix.index)
        dec_ids = group.index[group["inc_dec"] == "dec"].intersection(plot_matrix.index)
        if len(inc_ids) == 0 or len(dec_ids) == 0:
            continue

        inc_values = plot_matrix.loc[inc_ids].mean(axis=0)
        dec_values = plot_matrix.loc[dec_ids].mean(axis=0)
        pair_id = f"{study_id}|{experiment_id}|paired_inc_minus_dec"
        pair_values = inc_values - dec_values
        pair_values.name = pair_id
        rows.append(pair_values)

        pair_metadata = group.iloc[0]
        metadata_rows.append(
            {
                "pair_id": pair_id,
                "study_id": study_id,
                "experiment_id": experiment_id,
                "direction_label": PAIR_DIRECTION_LABEL,
                "increased_profile_count": len(inc_ids),
                "decreased_profile_count": len(dec_ids),
                "increased_profile_ids": ";".join(inc_ids),
                "decreased_profile_ids": ";".join(dec_ids),
                "mean_paired_trait_enrichment": pair_values.mean(),
                "max_abs_paired_trait_enrichment": pair_values.abs().max(),
                "body_site": pair_metadata["body_site"],
                "body_site_group": pair_metadata["body_site_group"],
                "condition": pair_metadata["condition"],
                "condition_group": pair_metadata["condition_group"],
                "location": pair_metadata["location"],
                "host_species": pair_metadata["host_species"],
                "title": pair_metadata["title"],
                "publication_year": pair_metadata["publication_year"],
            }
        )

    if not rows:
        msg = "No matched increased/decreased study-experiment pairs were available."
        raise ValueError(msg)

    pair_matrix = pd.DataFrame(rows)
    pair_metadata = pd.DataFrame.from_records(metadata_rows).set_index("pair_id")
    return pair_matrix, pair_metadata


def summarize_required_paired_traits(
    pair_matrix: pd.DataFrame,
    trait_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize requested trait checks in the paired enrichment matrix."""
    rows: list[dict[str, object]] = []
    for trait_id, requested_trait_name in REQUIRED_PAIRED_TRAITS.items():
        metadata_row = (
            trait_metadata.loc[trait_id]
            if trait_id in trait_metadata.index
            else pd.Series(dtype=object)
        )
        if trait_id in pair_matrix.columns:
            values = pd.to_numeric(pair_matrix[trait_id], errors="coerce").fillna(0.0)
            mean_contrast = float(values.mean())
            median_contrast = float(values.median())
            std_contrast = float(values.std(ddof=1))
            positive_fraction = float((values > 0).mean())
            negative_fraction = float((values < 0).mean())
            pair_count = len(values)
            direction = (
                "trait enriched in increased profiles"
                if mean_contrast >= 0
                else "trait enriched in decreased profiles"
            )
        else:
            mean_contrast = math.nan
            median_contrast = math.nan
            std_contrast = math.nan
            positive_fraction = math.nan
            negative_fraction = math.nan
            pair_count = 0
            direction = "not present"

        rows.append(
            {
                "trait_id": trait_id,
                "requested_trait_name": requested_trait_name,
                "heatmap_trait_name": metadata_row.get("object_name", ""),
                "display_label": clean_value(metadata_row.get("display_label", "")).replace(
                    "\n", " "
                ),
                "profile_support": metadata_row.get("profile_support", pd.NA),
                "total_signal": metadata_row.get("total_signal", pd.NA),
                "composition_variance": metadata_row.get("composition_variance", pd.NA),
                "selection_score": metadata_row.get("selection_score", pd.NA),
                "included_in_paired_heatmap": trait_id in pair_matrix.columns,
                "pair_count": pair_count,
                "mean_increased_minus_decreased": mean_contrast,
                "median_increased_minus_decreased": median_contrast,
                "std_increased_minus_decreased": std_contrast,
                "positive_pair_fraction": positive_fraction,
                "negative_pair_fraction": negative_fraction,
                "direction": direction,
            }
        )

    return pd.DataFrame.from_records(rows)


def build_paired_display_matrix(
    pair_matrix: pd.DataFrame,
    pair_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Apply readable labels for the plotted paired-enrichment matrix."""
    display = pair_matrix.copy()
    display.index = [
        " | ".join(
            [
                pair_metadata.loc[pair_id, "study_id"],
                pair_metadata.loc[pair_id, "experiment_id"],
                "increased - decreased",
            ]
        )
        for pair_id in display.index
    ]
    display.columns = [
        trait_metadata.loc[trait_id, "display_label"] or trait_id for trait_id in display.columns
    ]
    return display


def build_rendered_paired_enrichment_matrix(
    pair_matrix: pd.DataFrame,
    pair_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    max_display_rows: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate adjacent clustered pairs into display bins for a readable PNG."""
    ordered_pairs = cluster_index(pair_matrix)
    pair_matrix = pair_matrix.loc[ordered_pairs]
    pair_metadata = pair_metadata.loc[ordered_pairs]

    if max_display_rows <= 0 or len(pair_matrix) <= max_display_rows:
        display_matrix = build_paired_display_matrix(pair_matrix, pair_metadata, trait_metadata)
        display_metadata = pair_metadata.copy()
        display_metadata.index = display_matrix.index
        display_metadata["pair_count"] = 1
        display_metadata["study_count_in_bin"] = 1
        display_metadata["first_pair_id"] = ordered_pairs
        display_metadata["last_pair_id"] = ordered_pairs
        return display_matrix, display_metadata

    rows: list[pd.Series] = []
    metadata_rows: list[dict[str, object]] = []
    for bin_number, bin_pair_ids in enumerate(
        np.array_split(ordered_pairs.to_numpy(), max_display_rows),
        start=1,
    ):
        if len(bin_pair_ids) == 0:
            continue

        bin_index = pd.Index(bin_pair_ids)
        bin_values = pair_matrix.loc[bin_index].mean(axis=0)
        bin_metadata = pair_metadata.loc[bin_index]
        label = (
            f"paired enrichment bin {bin_number:03d} | "
            f"pairs={len(bin_index)} | studies={bin_metadata['study_id'].nunique()}"
        )
        bin_values.name = label
        rows.append(bin_values)
        metadata_rows.append(
            {
                "display_label": label,
                "direction_label": PAIR_DIRECTION_LABEL,
                "pair_count": len(bin_index),
                "study_count_in_bin": bin_metadata["study_id"].nunique(),
                "body_site_group": most_common_series_value(bin_metadata["body_site_group"]),
                "condition_group": most_common_series_value(bin_metadata["condition_group"]),
                "mean_paired_trait_enrichment": bin_values.mean(),
                "max_abs_paired_trait_enrichment": bin_values.abs().max(),
                "first_pair_id": bin_index[0],
                "last_pair_id": bin_index[-1],
            }
        )

    display_matrix = pd.DataFrame(rows)
    display_matrix.columns = [
        trait_metadata.loc[trait_id, "display_label"] or trait_id
        for trait_id in display_matrix.columns
    ]
    display_metadata = pd.DataFrame.from_records(metadata_rows).set_index("display_label")
    return display_matrix, display_metadata


def palette_for_categories(categories: pd.Series, palette_name: str) -> dict[str, object]:
    """Build a stable palette for row or column annotation categories."""
    ordered = [
        category for category in sorted(categories.unique()) if category not in {"Other", "unknown"}
    ]
    colors: dict[str, object] = {}
    if ordered:
        colors = dict(
            zip(
                ordered,
                sns.color_palette(palette_name, n_colors=len(ordered)),
                strict=True,
            )
        )
    colors["Other"] = "#b0b0b0"
    colors["unknown"] = "#e5e7eb"
    return colors


def legend_handles(palette: dict[str, object], present_categories: pd.Series) -> list[Patch]:
    """Build legend handles only for categories present in the plot."""
    present = set(present_categories)
    return [Patch(color=color, label=label) for label, color in palette.items() if label in present]


def save_heatmap(
    display_matrix: pd.DataFrame,
    display_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
    full_profile_count: int,
    full_study_count: int,
    direction_basis: str,
) -> None:
    """Render and save the clustered study profile heatmap."""
    sns.set_theme(style="white")
    display_to_trait_id = dict(zip(display_matrix.columns, trait_metadata.index, strict=True))

    direction_palette = {"inc": "#c0392b", "dec": "#2874a6"}
    body_palette = palette_for_categories(display_metadata["body_site_group"], "tab20")
    condition_palette = palette_for_categories(display_metadata["condition_group"], "Set3")
    trait_palette = palette_for_categories(trait_metadata["display_category"], "Set2")

    row_colors = pd.DataFrame(
        {
            "direction": display_metadata.loc[display_matrix.index, "inc_dec"]
            .map(direction_palette)
            .fillna("#7f8c8d")
            .to_numpy(),
            "body_site": display_metadata.loc[display_matrix.index, "body_site_group"]
            .map(body_palette)
            .to_numpy(),
            "condition": display_metadata.loc[display_matrix.index, "condition_group"]
            .map(condition_palette)
            .to_numpy(),
        },
        index=display_matrix.index,
    )
    col_trait_ids = pd.Index([display_to_trait_id[label] for label in display_matrix.columns])
    col_colors = pd.DataFrame(
        {
            "trait_category": trait_metadata.loc[col_trait_ids, "display_category"]
            .map(trait_palette)
            .to_numpy(),
        },
        index=display_matrix.columns,
    )

    cluster = sns.clustermap(
        display_matrix,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        method="average",
        metric="correlation",
        row_cluster=False,
        col_cluster=True,
        row_colors=row_colors,
        col_colors=col_colors,
        linewidths=0.05,
        linecolor="#e5e7eb",
        xticklabels=True,
        yticklabels=len(display_matrix) <= 120,
        figsize=(24, min(40, max(18, len(display_matrix) / 8))),
        cbar_pos=(0.02, 0.80, 0.02, 0.14),
        dendrogram_ratio=(0.12, 0.10),
        colors_ratio=(0.03, 0.03),
        cbar_kws={
            "label": (
                "Column z-score of log1p non-MediaDive trait fraction per 1,000 retained counts; "
                "display capped at +/-1"
            )
        },
    )
    cluster.fig.suptitle(
        "Study Trait Profile Clustering\n"
        f"{full_profile_count} BugSigDB study/experiment/direction profiles from "
        f"{full_study_count} studies, aggregated into {len(display_matrix)} "
        "clustered display bins; "
        "increased and decreased sections are separated",
        y=1.08,
        fontsize=16,
        fontweight="bold",
    )
    cluster.ax_heatmap.set_xlabel("Trait")
    cluster.ax_heatmap.set_ylabel("Clustered profile display bin")
    increased_count = int((display_metadata["inc_dec"] == "inc").sum())
    if 0 < increased_count < len(display_matrix):
        cluster.ax_heatmap.axhline(increased_count, color="#111827", linewidth=0.8)
    for label in cluster.ax_heatmap.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(7)
    for label in cluster.ax_heatmap.get_yticklabels():
        label.set_fontsize(7)

    cluster.ax_row_dendrogram.legend(
        handles=legend_handles(direction_palette, display_metadata["inc_dec"]),
        title=f"BugSigDB direction\n({DIRECTION_BASIS_LABELS[direction_basis]})",
        loc="center",
        bbox_to_anchor=(-0.25, 0.65),
        frameon=False,
    )
    cluster.fig.legend(
        handles=legend_handles(trait_palette, trait_metadata["display_category"]),
        title="Trait category",
        loc="upper center",
        bbox_to_anchor=(0.5, 1.16),
        ncol=4,
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )
    cluster.fig.legend(
        handles=legend_handles(body_palette, display_metadata["body_site_group"]),
        title="Body site",
        loc="lower left",
        bbox_to_anchor=(0.02, 0.18),
        ncol=1,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    cluster.fig.legend(
        handles=legend_handles(condition_palette, display_metadata["condition_group"]),
        title="Condition",
        loc="lower left",
        bbox_to_anchor=(0.02, 0.02),
        ncol=1,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster.fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(cluster.fig)


def save_paired_enrichment_heatmap(
    display_matrix: pd.DataFrame,
    display_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
    full_pair_count: int,
    full_study_count: int,
    direction_basis: str,
) -> None:
    """Render and save the paired increased-minus-decreased enrichment heatmap."""
    sns.set_theme(style="white")
    display_to_trait_id = dict(zip(display_matrix.columns, trait_metadata.index, strict=True))

    body_palette = palette_for_categories(display_metadata["body_site_group"], "tab20")
    condition_palette = palette_for_categories(display_metadata["condition_group"], "Set3")
    trait_palette = palette_for_categories(trait_metadata["display_category"], "Set2")

    row_colors = pd.DataFrame(
        {
            "body_site": display_metadata.loc[display_matrix.index, "body_site_group"]
            .map(body_palette)
            .to_numpy(),
            "condition": display_metadata.loc[display_matrix.index, "condition_group"]
            .map(condition_palette)
            .to_numpy(),
        },
        index=display_matrix.index,
    )
    col_trait_ids = pd.Index([display_to_trait_id[label] for label in display_matrix.columns])
    col_colors = pd.DataFrame(
        {
            "trait_category": trait_metadata.loc[col_trait_ids, "display_category"]
            .map(trait_palette)
            .to_numpy(),
        },
        index=display_matrix.columns,
    )

    cluster = sns.clustermap(
        display_matrix,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        method="average",
        metric="correlation",
        row_cluster=False,
        col_cluster=True,
        row_colors=row_colors,
        col_colors=col_colors,
        linewidths=0.05,
        linecolor="#e5e7eb",
        xticklabels=True,
        yticklabels=len(display_matrix) <= 120,
        figsize=(24, min(40, max(18, len(display_matrix) / 8))),
        cbar_pos=(0.02, 0.80, 0.02, 0.14),
        dendrogram_ratio=(0.12, 0.10),
        colors_ratio=(0.03, 0.03),
        cbar_kws={
            "label": (
                "Paired enrichment statistic: increased minus decreased\n"
                "non-MediaDive profile z-score; display capped at +/-1\n"
                f"BugSigDB direction: {DIRECTION_BASIS_LABELS[direction_basis]}"
            )
        },
    )
    cluster.fig.suptitle(
        "Paired Study Trait Enrichment\n"
        f"{full_pair_count} matched BugSigDB study/experiment pairs from "
        f"{full_study_count} studies, aggregated into {len(display_matrix)} "
        "clustered display bins",
        y=1.08,
        fontsize=16,
        fontweight="bold",
    )
    cluster.ax_heatmap.set_xlabel("Trait")
    cluster.ax_heatmap.set_ylabel("Clustered matched-pair display bin")
    for label in cluster.ax_heatmap.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(7)
    for label in cluster.ax_heatmap.get_yticklabels():
        label.set_fontsize(7)

    cluster.fig.legend(
        handles=legend_handles(trait_palette, trait_metadata["display_category"]),
        title="Trait category",
        loc="upper center",
        bbox_to_anchor=(0.5, 1.16),
        ncol=4,
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )
    cluster.fig.legend(
        handles=legend_handles(body_palette, display_metadata["body_site_group"]),
        title="Body site",
        loc="lower left",
        bbox_to_anchor=(0.02, 0.18),
        ncol=1,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    cluster.fig.legend(
        handles=legend_handles(condition_palette, display_metadata["condition_group"]),
        title="Condition",
        loc="lower left",
        bbox_to_anchor=(0.02, 0.02),
        ncol=1,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster.fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(cluster.fig)


def write_outputs(
    plot_matrix: pd.DataFrame,
    display_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    rendered_rows: int,
    output_path: Path,
    matrix_output_path: Path,
    profile_metadata_output_path: Path,
    trait_metadata_output_path: Path,
) -> None:
    """Write the plotted matrix and source metadata."""
    matrix_output_path.parent.mkdir(parents=True, exist_ok=True)
    display_matrix.to_csv(matrix_output_path, sep="\t")
    profile_ids = plot_matrix.index
    trait_ids = plot_matrix.columns
    profile_metadata.loc[profile_ids].to_csv(profile_metadata_output_path, sep="\t")
    trait_metadata.loc[trait_ids].to_csv(trait_metadata_output_path, sep="\t")
    print(f"Wrote heatmap: {output_path}")
    print(f"Wrote matrix: {matrix_output_path}")
    print(f"Wrote profile metadata: {profile_metadata_output_path}")
    print(f"Wrote trait metadata: {trait_metadata_output_path}")
    print("Measurement: column z-score of log1p trait fraction per 1,000 profile counts")
    print(f"Profiles analyzed and written: {len(profile_ids)}")
    print(f"Rendered heatmap display rows: {rendered_rows}")
    print(f"Studies represented: {profile_metadata.loc[profile_ids, 'study_id'].nunique()}")
    print(
        "Increased/decreased profiles: "
        f"{int((profile_metadata.loc[profile_ids, 'inc_dec'] == 'inc').sum())}/"
        f"{int((profile_metadata.loc[profile_ids, 'inc_dec'] == 'dec').sum())}"
    )
    print(f"Traits plotted: {len(trait_ids)}")


def write_paired_outputs(
    pair_matrix: pd.DataFrame,
    pair_display_matrix: pd.DataFrame,
    pair_metadata: pd.DataFrame,
    required_trait_check: pd.DataFrame,
    rendered_rows: int,
    output_path: Path,
    matrix_output_path: Path,
    metadata_output_path: Path,
    required_trait_check_output_path: Path,
    *,
    source_hashes: dict[str, str] | None = None,
    selection_parameters: dict[str, int | float | str] | None = None,
) -> None:
    """Write the paired enrichment matrix and matched-pair metadata."""
    matrix_output_path.parent.mkdir(parents=True, exist_ok=True)
    pair_display_matrix.to_csv(matrix_output_path, sep="\t")
    metadata = pair_metadata.loc[pair_matrix.index].copy()
    for source, digest in (source_hashes or {}).items():
        metadata[f"{source}_sha256"] = digest
    for parameter, value in (selection_parameters or {}).items():
        metadata[f"heatmap_parameter_{parameter}"] = value
    metadata.to_csv(metadata_output_path, sep="\t")
    required_trait_check.to_csv(required_trait_check_output_path, sep="\t", index=False)
    print(f"Wrote paired enrichment heatmap: {output_path}")
    print(f"Wrote paired enrichment matrix: {matrix_output_path}")
    print(f"Wrote paired enrichment metadata: {metadata_output_path}")
    print(f"Wrote paired requested-trait check: {required_trait_check_output_path}")
    print("Measurement: matched increased-minus-decreased profile z-score per trait")
    print(f"Matched study/experiment pairs analyzed and written: {len(pair_matrix)}")
    print(f"Rendered paired heatmap display rows: {rendered_rows}")
    print(f"Studies represented by matched pairs: {pair_metadata['study_id'].nunique()}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Cluster studies by aggregate trait profiles to inspect whether similar "
            "conditions or body sites converge on similar microbial traits."
        )
    )
    parser.add_argument(
        "--direction-basis",
        choices=sorted(DIRECTION_BASIS_LABELS),
        required=True,
        help=(
            "Direction basis of the input profiles, as recorded by the core: "
            "group_1_vs_group_0 for all paired contrasts or disease_vs_comparison "
            "for an oriented subset. Both heatmap legends name this basis; the "
            "script does not reorient input values."
        ),
    )
    parser.add_argument(
        "--study-trait-profiles-input",
        type=Path,
        default=REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv",
        help="Input wide study trait profiles TSV path.",
    )
    parser.add_argument(
        "--trait-columns-input",
        type=Path,
        default=REPO_ROOT / "data/derived/taxa_trait_matrix_columns.tsv",
        help="Input trait column metadata TSV path.",
    )
    parser.add_argument(
        "--studies-input",
        type=Path,
        default=REPO_ROOT / "data/derived/studies.tsv",
        help="Input derived studies TSV path.",
    )
    parser.add_argument(
        "--full-dump-input",
        type=Path,
        default=REPO_ROOT / "data/raw/full_dump.csv",
        help="Input BugSigDB full_dump.csv path.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_clustered_heatmap.png",
        help="Output heatmap image path.",
    )
    parser.add_argument(
        "--matrix-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_clustered_heatmap_matrix.tsv",
        help="Output plotted z-score matrix TSV path.",
    )
    parser.add_argument(
        "--study-metadata-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_study_metadata.tsv",
        help="Output selected profile metadata TSV path.",
    )
    parser.add_argument(
        "--trait-metadata-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_trait_metadata.tsv",
        help="Output selected trait metadata TSV path.",
    )
    parser.add_argument(
        "--paired-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_paired_enrichment_heatmap.png",
        help="Output paired increased-minus-decreased enrichment heatmap image path.",
    )
    parser.add_argument(
        "--paired-matrix-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_paired_enrichment_matrix.tsv",
        help="Output paired increased-minus-decreased enrichment matrix TSV path.",
    )
    parser.add_argument(
        "--paired-metadata-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_paired_enrichment_metadata.tsv",
        help="Output matched-pair metadata TSV path.",
    )
    parser.add_argument(
        "--paired-required-trait-check-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_paired_required_trait_check.tsv",
        help=(
            "Output requested-trait audit TSV for required paired heatmap traits "
            "(mesophilic, gram positive, non motile)."
        ),
    )
    parser.add_argument(
        "--max-profiles",
        type=int,
        default=0,
        help="Maximum number of profile rows to plot. Use 0 to include every profile.",
    )
    parser.add_argument(
        "--max-studies",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--max-traits", type=int, default=60)
    parser.add_argument(
        "--display-max-rows",
        type=int,
        default=320,
        help=(
            "Maximum rows to render in the PNG. When fewer than all profiles, adjacent "
            "clustered profiles are averaged into direction-specific display bins."
        ),
    )
    parser.add_argument(
        "--paired-display-max-rows",
        type=int,
        default=320,
        help=(
            "Maximum rows to render in the paired enrichment PNG. When fewer than all "
            "matched pairs, adjacent clustered pairs are averaged into display bins."
        ),
    )
    parser.add_argument("--min-traits", type=int, default=0)
    parser.add_argument("--min-total-signal", type=int, default=0)
    parser.add_argument(
        "--min-profile-support",
        "--min-study-support",
        dest="min_profile_support",
        type=int,
        default=6,
    )
    parser.add_argument(
        "--max-profile-support-fraction",
        "--max-study-support-fraction",
        dest="max_profile_support_fraction",
        type=float,
        default=0.95,
    )
    parser.add_argument("--max-body-site-categories", type=int, default=8)
    parser.add_argument("--max-condition-categories", type=int, default=8)
    parser.add_argument(
        "--include-media-traits",
        action="store_true",
        help="Rejected: MediaDive growth-medium traits are excluded before normalization.",
    )
    args = parser.parse_args()
    if args.include_media_traits:
        parser.error("MediaDive is excluded before composition normalization (#174).")

    pair_source_inputs = {
        "profile": args.study_trait_profiles_input,
        "trait_labels": args.trait_columns_input,
        "raw_dump": args.full_dump_input,
        "studies": args.studies_input,
    }
    pair_source_hashes = {name: file_sha256(path) for name, path in pair_source_inputs.items()}
    profile_matrix, profile_metadata, _trait_columns = load_profile_matrix(
        args.study_trait_profiles_input
    )
    profile_metadata = join_profile_labels(
        profile_metadata,
        args.studies_input,
        args.full_dump_input,
    )

    unfiltered_profile_count = len(profile_metadata)
    unfiltered_study_count = profile_metadata["study_id"].nunique()
    profile_matrix, profile_metadata = filter_human_host_profiles(
        profile_matrix,
        profile_metadata,
    )
    print(
        "Human-host corpus filter retained "
        f"{len(profile_metadata):,}/{unfiltered_profile_count:,} profiles from "
        f"{profile_metadata['study_id'].nunique():,}/{unfiltered_study_count:,} studies."
    )

    max_profiles = args.max_profiles
    if args.max_studies is not None and max_profiles == 0:
        max_profiles = args.max_studies
    selected_profiles = select_profiles(
        profile_matrix,
        profile_metadata,
        max_profiles,
        args.min_traits,
        args.min_total_signal,
    )
    selected_profile_matrix = profile_matrix.loc[selected_profiles]
    profile_metadata = profile_metadata.loc[selected_profiles].copy()
    profile_metadata["body_site_group"] = group_top_categories(
        profile_metadata["body_site"],
        args.max_body_site_categories,
    )
    profile_metadata["condition_group"] = group_top_categories(
        profile_metadata["condition"],
        args.max_condition_categories,
    )

    trait_labels = load_trait_labels(args.trait_columns_input)
    selected_traits, trait_metadata = select_traits(
        selected_profile_matrix,
        trait_labels,
        args.max_traits,
        args.min_profile_support,
        args.max_profile_support_fraction,
        args.include_media_traits,
    )
    plot_matrix = build_plot_matrix(selected_profile_matrix, selected_traits)
    ordered_profiles = order_profiles_by_direction(plot_matrix, profile_metadata)
    plot_matrix = plot_matrix.loc[ordered_profiles]
    profile_metadata = profile_metadata.loc[ordered_profiles]
    display_matrix = build_display_matrix(plot_matrix, profile_metadata, trait_metadata)
    rendered_matrix, rendered_metadata = build_rendered_display_matrix(
        plot_matrix,
        profile_metadata,
        trait_metadata,
        args.display_max_rows,
    )
    save_heatmap(
        rendered_matrix,
        rendered_metadata,
        trait_metadata,
        args.output,
        full_profile_count=len(plot_matrix),
        full_study_count=profile_metadata["study_id"].nunique(),
        direction_basis=args.direction_basis,
    )
    write_outputs(
        plot_matrix,
        display_matrix,
        profile_metadata,
        trait_metadata,
        len(rendered_matrix),
        args.output,
        args.matrix_output,
        args.study_metadata_output,
        args.trait_metadata_output,
    )

    pair_matrix, pair_metadata = build_paired_enrichment_matrix(plot_matrix, profile_metadata)
    ordered_pairs = cluster_index(pair_matrix)
    pair_matrix = pair_matrix.loc[ordered_pairs]
    pair_metadata = pair_metadata.loc[ordered_pairs]
    pair_display_matrix = build_paired_display_matrix(
        pair_matrix,
        pair_metadata,
        trait_metadata,
    )
    rendered_pair_matrix, rendered_pair_metadata = build_rendered_paired_enrichment_matrix(
        pair_matrix,
        pair_metadata,
        trait_metadata,
        args.paired_display_max_rows,
    )
    save_paired_enrichment_heatmap(
        rendered_pair_matrix,
        rendered_pair_metadata,
        trait_metadata,
        args.paired_output,
        full_pair_count=len(pair_matrix),
        full_study_count=pair_metadata["study_id"].nunique(),
        direction_basis=args.direction_basis,
    )
    required_trait_check = summarize_required_paired_traits(pair_matrix, trait_metadata)
    if pair_source_hashes != {name: file_sha256(path) for name, path in pair_source_inputs.items()}:
        raise ValueError("Heatmap source inputs changed during execution; rebuild before export")
    write_paired_outputs(
        pair_matrix,
        pair_display_matrix,
        pair_metadata,
        required_trait_check,
        len(rendered_pair_matrix),
        args.paired_output,
        args.paired_matrix_output,
        args.paired_metadata_output,
        args.paired_required_trait_check_output,
        source_hashes=pair_source_hashes,
        selection_parameters={
            "max_profiles": max_profiles,
            "max_body_site_categories": args.max_body_site_categories,
            "max_condition_categories": args.max_condition_categories,
            "min_traits": args.min_traits,
            "min_total_signal": args.min_total_signal,
            "max_traits": args.max_traits,
            "min_profile_support": args.min_profile_support,
            "max_profile_support_fraction": args.max_profile_support_fraction,
            "display_max_rows": args.display_max_rows,
            "paired_display_max_rows": args.paired_display_max_rows,
            "host": HUMAN_HOST_SPECIES,
            "sequencing_types": ",".join(sorted(CORPUS_SEQUENCING_TYPES)),
        },
    )


if __name__ == "__main__":
    main()
