"""Follow up the enrichment-analysis infographic outlier conditions."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path
from textwrap import shorten

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = (
    REPO_ROOT
    / "trait-visualization"
    / "trait-profile-investigation"
    / ("enrichment-infographic-followup")
)
PROFILE_METADATA_COLUMNS = ["study_id", "experiment_id", "inc_dec"]
TARGET_CONDITIONS = (
    "Thalassemia",
    "Sleep Disorder",
    "Head and neck squamous cell carcinoma",
)
INFOGRAPHIC_AVG_ENRICHMENT_DEPTH = {
    "Thalassemia": 225.00,
    "Sleep Disorder": 216.00,
    "Head and neck squamous cell carcinoma": 210.33,
}
INFOGRAPHIC_TOP_TRAITS = {
    "METPO:1000603": "anaerobic",
    "METPO:1000615": "mesophilic",
    "METPO:1000607": "obligately anaerobic",
    "METPO:1000698": "gram positive",
    "METPO:1000699": "gram negative",
    "METPO:1000703": "non motile",
    "METPO:1000432": "GC high",
    "METPO:1000872": "non-spore forming",
    "METPO:1002005": "Fermentation",
    "METPO:1000604": "microaerophilic",
}
GUILD_TRAITS = {
    "aerobe_v2": {
        "METPO:1000602": "aerobic",
        "METPO:1000604": "microaerophilic",
        "METPO:1000608": "facultatively aerobic",
    },
    "anaerobe_v2_depletion": {
        "METPO:1000603": "anaerobic",
        "METPO:1000607": "obligately anaerobic",
    },
}
AUTISM_LABELS = ("Autism", "Autism spectrum disorder")
DIRECTION_TO_INC_DEC = {"increased": "inc", "decreased": "dec"}
INC_DEC_TO_DIRECTION = {"inc": "increased", "dec": "decreased"}
MATRIX_METADATA_COLUMNS = ["taxon_id", "taxon_name", "taxonomic_rank"]
PREFERRED_HEATMAP_RANKS = {"genus", "species", "strain"}
METAPHLAN_RANK_MAP = {
    "k__": "superkingdom",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
    "t__": "strain",
}


def clean_text(value: object) -> str:
    """Normalize empty BugSigDB metadata values."""
    text = str(value or "").strip()
    if not text or text.upper() == "NA" or text.lower() == "nan":
        return "unknown"
    return text


def clean_optional_text(value: object) -> str:
    """Normalize optional fields where a blank string is more useful than unknown."""
    text = clean_text(value)
    return "" if text == "unknown" else text


def compact_label(value: object, max_length: int = 44) -> str:
    """Create a compact one-line display label."""
    return shorten(clean_text(value).replace("\n", " "), width=max_length, placeholder="...")


def study_id_from_raw_study(value: object) -> str:
    """Convert a raw BugSigDB Study field to the local study_id format."""
    return f"bsdb:{clean_text(value)}"


def experiment_number_from_name(value: object) -> str:
    """Extract the integer suffix from a BugSigDB experiment page name."""
    match = re.search(r"(\d+)", clean_text(value))
    return match.group(1) if match else ""


def experiment_id_from_raw_parts(study: object, experiment_name: object) -> str:
    """Build a local experiment_id from raw BugSigDB experiment metadata."""
    experiment_number = experiment_number_from_name(experiment_name)
    if not experiment_number:
        return f"{study_id_from_raw_study(study)}-unknown"
    return f"{study_id_from_raw_study(study)}-{experiment_number}"


def study_id_from_bsdb_id(bsdb_id: object) -> str:
    """Convert a full raw signature ID to a local study_id."""
    compact_id = clean_text(bsdb_id).removeprefix("bsdb:").split("/", 1)[0]
    return f"bsdb:{compact_id}"


def experiment_id_from_bsdb_id(bsdb_id: object) -> str:
    """Convert a full raw signature ID to a local experiment_id."""
    parts = clean_text(bsdb_id).removeprefix("bsdb:").split("/")
    if len(parts) >= 2:
        return f"bsdb:{parts[0]}-{parts[1]}"
    return f"bsdb:{parts[0]}-unknown"


def profile_id_from_parts(study_id: object, experiment_id: object, inc_dec: object) -> str:
    """Build the stable profile ID used by study_trait_profiles_wide.tsv."""
    return f"{study_id}|{experiment_id}|{inc_dec}"


def profile_id_from_raw_row(row: dict[str, str]) -> str:
    """Build a profile ID from a raw full_dump.csv row."""
    direction = clean_text(row.get("Abundance in Group 1")).lower()
    inc_dec = DIRECTION_TO_INC_DEC.get(direction, direction)
    bsdb_id = clean_text(row.get("BSDB ID"))
    return profile_id_from_parts(
        study_id_from_bsdb_id(bsdb_id),
        experiment_id_from_bsdb_id(bsdb_id),
        inc_dec,
    )


def experiment_label(row: pd.Series) -> str:
    """Create a readable experiment comparison label."""
    name = clean_text(row.get("experiment_name"))
    group_0 = clean_text(row.get("group_0_name"))
    group_1 = clean_text(row.get("group_1_name"))
    if group_0 == "unknown" and group_1 == "unknown":
        return name
    return f"{name}: {group_0} vs {group_1}"


def most_common(values: pd.Series) -> str:
    """Return a deterministic most common non-empty value."""
    cleaned = [clean_text(value) for value in values if clean_text(value) != "unknown"]
    if not cleaned:
        return "unknown"
    counts = Counter(cleaned)
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def join_unique(values: pd.Series, max_items: int = 12) -> str:
    """Join unique non-empty values with deterministic ordering."""
    unique_values = sorted(
        {clean_text(value) for value in values if clean_text(value) != "unknown"}
    )
    if len(unique_values) > max_items:
        shown = unique_values[:max_items]
        return "; ".join([*shown, f"+{len(unique_values) - max_items} more"])
    return "; ".join(unique_values)


def load_raw_experiments(path: Path) -> pd.DataFrame:
    """Load raw BugSigDB experiment metadata with local IDs."""
    raw = pd.read_csv(path, dtype=str, comment="#").fillna("")
    raw["study_id"] = raw["Study"].map(study_id_from_raw_study)
    raw["experiment_id"] = [
        experiment_id_from_raw_parts(study, experiment_name)
        for study, experiment_name in raw[["Study", "Experiment page name"]].itertuples(
            index=False,
            name=None,
        )
    ]
    raw = raw.rename(
        columns={
            "Experiment page name": "experiment_name",
            "Location of subjects": "location",
            "Host species": "host_species",
            "Body site": "body_site",
            "UBERON ID": "uberon_id",
            "Condition": "condition",
            "EFO ID": "efo_id",
            "Group 0 name": "group_0_name",
            "Group 1 name": "group_1_name",
            "Group 1 definition": "group_1_definition",
            "Group 0 sample size": "group_0_sample_size",
            "Group 1 sample size": "group_1_sample_size",
            "Sequencing type": "sequencing_type",
        }
    )
    columns = [
        "study_id",
        "experiment_id",
        "experiment_name",
        "location",
        "host_species",
        "body_site",
        "uberon_id",
        "condition",
        "efo_id",
        "group_0_name",
        "group_1_name",
        "group_1_definition",
        "group_0_sample_size",
        "group_1_sample_size",
        "sequencing_type",
    ]
    for column in columns:
        raw[column] = raw[column].map(clean_text)
    raw["experiment_label"] = raw.apply(experiment_label, axis=1)
    return raw[[*columns, "experiment_label"]].drop_duplicates(
        subset=["study_id", "experiment_id"],
        keep="first",
    )


def load_studies(path: Path) -> pd.DataFrame:
    """Load study title metadata."""
    studies = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    columns = [
        "study_id",
        "pmid",
        "title",
        "publication_year",
        "doi",
        "url",
        "experiment_count",
        "signature_count",
        "taxon_count",
    ]
    studies = studies[columns].copy()
    for column in columns:
        studies[column] = studies[column].map(clean_text)
    return studies.drop_duplicates(subset=["study_id"], keep="first")


def load_trait_labels(path: Path) -> pd.DataFrame:
    """Load readable labels and categories for trait IDs."""
    labels = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    labels = labels.set_index("column_id", drop=False)
    labels["display_name"] = [
        clean_text(row.get("object_name"))
        if clean_text(row.get("object_name")) != "unknown"
        else clean_text(row.get("column_label"))
        for _, row in labels.iterrows()
    ]
    labels["display_category"] = labels["object_category"].str.removeprefix("biolink:")
    labels.loc[labels["display_category"] == "", "display_category"] = "unknown"
    labels["heatmap_label"] = [
        f"{compact_label(row['display_name'], 31)}\n{compact_label(column_id, 24)}"
        for column_id, row in labels.iterrows()
    ]
    return labels


def load_profiles(
    profile_path: Path,
    experiments: pd.DataFrame,
    studies: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Load the full study-trait profile matrix and profile metadata."""
    profiles = pd.read_csv(profile_path, sep="\t")
    trait_columns = [
        column for column in profiles.columns if column not in PROFILE_METADATA_COLUMNS
    ]
    values = (
        profiles[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0).astype("float32")
    )
    profile_ids = [
        profile_id_from_parts(study_id, experiment_id, inc_dec)
        for study_id, experiment_id, inc_dec in profiles[PROFILE_METADATA_COLUMNS].itertuples(
            index=False,
            name=None,
        )
    ]
    values.index = profile_ids

    metadata = profiles[PROFILE_METADATA_COLUMNS].copy()
    metadata["profile_id"] = profile_ids
    metadata["direction_label"] = (
        metadata["inc_dec"].map(INC_DEC_TO_DIRECTION).fillna(metadata["inc_dec"])
    )
    metadata["profile_trait_count"] = values.gt(0).sum(axis=1).to_numpy(dtype=int)
    metadata["total_trait_signal"] = values.sum(axis=1).to_numpy(dtype=float)
    metadata = metadata.merge(experiments, on=["study_id", "experiment_id"], how="left")
    metadata = metadata.merge(studies, on="study_id", how="left", suffixes=("", "_study"))
    for column in [
        "condition",
        "body_site",
        "location",
        "host_species",
        "title",
        "publication_year",
        "efo_id",
        "experiment_label",
    ]:
        if column not in metadata:
            metadata[column] = "unknown"
        metadata[column] = metadata[column].map(clean_text)
    metadata = metadata.set_index("profile_id", drop=False)
    return values, metadata, trait_columns


def summarize_condition_depth(metadata: pd.DataFrame, min_profiles: int) -> pd.DataFrame:
    """Summarize trait-depth distributions by exact condition label."""
    source = metadata[metadata["condition"] != "unknown"].copy()
    summary = (
        source.groupby("condition", dropna=False)
        .agg(
            profile_count=("profile_id", "size"),
            nonzero_profile_count=("profile_trait_count", lambda values: int((values > 0).sum())),
            study_count=("study_id", "nunique"),
            experiment_count=("experiment_id", "nunique"),
            mean_profile_trait_count=("profile_trait_count", "mean"),
            median_profile_trait_count=("profile_trait_count", "median"),
            max_profile_trait_count=("profile_trait_count", "max"),
            mean_total_trait_signal=("total_trait_signal", "mean"),
        )
        .reset_index()
    )
    ranked = summary[summary["profile_count"] >= min_profiles].copy()
    ranked["depth_rank"] = ranked["mean_profile_trait_count"].rank(
        ascending=False,
        method="min",
    )
    ranked["depth_percentile"] = ranked["mean_profile_trait_count"].rank(pct=True)
    median = float(ranked["mean_profile_trait_count"].median())
    mad = float((ranked["mean_profile_trait_count"] - median).abs().median())
    if math.isclose(mad, 0.0):
        ranked["robust_z"] = 0.0
    else:
        ranked["robust_z"] = 0.6745 * (ranked["mean_profile_trait_count"] - median) / mad
    ranked["current_profile_depth_outlier"] = (ranked["depth_percentile"] >= 0.95) | (
        ranked["robust_z"] >= 3.0
    )
    ranked["infographic_avg_enrichment_depth"] = ranked["condition"].map(
        INFOGRAPHIC_AVG_ENRICHMENT_DEPTH
    )
    ranked["infographic_top_three_condition"] = ranked["condition"].isin(TARGET_CONDITIONS)
    return ranked.sort_values(
        ["infographic_top_three_condition", "infographic_avg_enrichment_depth", "depth_rank"],
        ascending=[False, False, True],
    )


def write_target_title_table(
    metadata: pd.DataFrame,
    experiments: pd.DataFrame,
    studies: pd.DataFrame,
    output_path: Path,
) -> pd.DataFrame:
    """Write the BugSigDB study and experiment labels behind target conditions."""
    profiles = metadata[metadata["condition"].isin(TARGET_CONDITIONS)].copy()
    profile_summary = (
        profiles.groupby(["study_id", "experiment_id"], dropna=False)
        .agg(
            profile_count=("profile_id", "size"),
            nonzero_profile_count=(
                "profile_trait_count",
                lambda values: int((values > 0).sum()),
            ),
            directions=("direction_label", join_unique),
            mean_profile_trait_count=("profile_trait_count", "mean"),
            max_profile_trait_count=("profile_trait_count", "max"),
            total_trait_signal=("total_trait_signal", "sum"),
        )
        .reset_index()
    )
    table = experiments[experiments["condition"].isin(TARGET_CONDITIONS)].merge(
        studies,
        on="study_id",
        how="left",
    )
    table = table.merge(profile_summary, on=["study_id", "experiment_id"], how="left")
    for column in [
        "profile_count",
        "nonzero_profile_count",
        "mean_profile_trait_count",
        "max_profile_trait_count",
        "total_trait_signal",
    ]:
        table[column] = pd.to_numeric(table[column], errors="coerce").fillna(0)
    ordered_columns = [
        "condition",
        "efo_id",
        "study_id",
        "experiment_id",
        "title",
        "publication_year",
        "doi",
        "url",
        "experiment_label",
        "group_1_definition",
        "body_site",
        "location",
        "host_species",
        "sequencing_type",
        "profile_count",
        "nonzero_profile_count",
        "directions",
        "mean_profile_trait_count",
        "max_profile_trait_count",
        "total_trait_signal",
    ]
    table = table[ordered_columns].sort_values(
        ["condition", "study_id", "experiment_id"],
        kind="stable",
    )
    table.to_csv(output_path, sep="\t", index=False)
    return table


def summarize_autism_labels(metadata: pd.DataFrame, experiments: pd.DataFrame) -> pd.DataFrame:
    """Summarize exact BugSigDB autism and ASD condition labels."""
    profiles = metadata[metadata["condition"].isin(AUTISM_LABELS)].copy()
    experiments = experiments[experiments["condition"].isin(AUTISM_LABELS)].copy()
    rows = []
    for condition in AUTISM_LABELS:
        condition_profiles = profiles[profiles["condition"] == condition]
        condition_experiments = experiments[experiments["condition"] == condition]
        rows.append(
            {
                "exact_condition_label": condition,
                "efo_ids": join_unique(condition_experiments["efo_id"]),
                "profile_count": len(condition_profiles),
                "nonzero_profile_count": int((condition_profiles["profile_trait_count"] > 0).sum()),
                "study_count": condition_profiles["study_id"].nunique(),
                "experiment_count": condition_profiles["experiment_id"].nunique(),
                "body_sites": join_unique(condition_profiles["body_site"]),
                "mean_profile_trait_count": condition_profiles["profile_trait_count"].mean(),
                "interpretation": (
                    "Distinct exact BugSigDB label; analyze separately unless an explicit "
                    "autism/ASD grouping is declared."
                ),
            }
        )
    output = pd.DataFrame(rows)
    output["recommended_group_label"] = "Autism / autism spectrum disorder"
    return output


def l2_normalize(matrix: np.ndarray) -> np.ndarray:
    """L2-normalize rows of a numeric matrix."""
    norms = np.linalg.norm(matrix, axis=1)
    normalized = matrix.copy()
    nonzero = norms > 0
    normalized[nonzero] = normalized[nonzero] / norms[nonzero, None]
    normalized[~nonzero] = 0.0
    return normalized


def compute_condition_centroids(
    normalized_matrix: np.ndarray,
    metadata: pd.DataFrame,
    profile_ids: list[str],
    min_profiles: int,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Compute unit centroids for condition-level profile vectors."""
    active_metadata = metadata.loc[profile_ids]
    centroids: dict[str, np.ndarray] = {}
    records: list[dict[str, object]] = []
    for condition, group in active_metadata.groupby("condition", sort=False):
        if condition == "unknown" or len(group) < min_profiles:
            continue
        positions = [profile_ids.index(profile_id) for profile_id in group.index]
        centroid = normalized_matrix[positions].mean(axis=0)
        norm = np.linalg.norm(centroid)
        if norm <= 0:
            continue
        centroid = centroid / norm
        centroids[condition] = centroid.astype("float32")
        records.append(
            {
                "condition": condition,
                "profile_count": len(group),
                "study_count": group["study_id"].nunique(),
                "experiment_count": group["experiment_id"].nunique(),
            }
        )
    return pd.DataFrame(records), centroids


def compute_profile_neighbors(
    values: pd.DataFrame,
    metadata: pd.DataFrame,
    top_n: int,
    min_condition_profiles: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Find the most similar individual profiles to each target condition centroid."""
    active_profile_ids = metadata.index[metadata["profile_trait_count"] > 0].tolist()
    active_values = np.log1p(values.loc[active_profile_ids].to_numpy(dtype="float32"))
    normalized = l2_normalize(active_values)
    condition_summary, centroids = compute_condition_centroids(
        normalized,
        metadata,
        active_profile_ids,
        min_condition_profiles,
    )

    neighbor_records: list[dict[str, object]] = []
    for query_condition in TARGET_CONDITIONS:
        if query_condition not in centroids:
            continue
        scores = normalized @ centroids[query_condition]
        candidates = metadata.loc[active_profile_ids].copy()
        candidates["cosine_similarity"] = scores
        candidates = candidates[candidates["condition"] != query_condition].copy()
        candidates = candidates.sort_values(
            ["cosine_similarity", "profile_trait_count", "study_id", "experiment_id", "inc_dec"],
            ascending=[False, False, True, True, True],
        ).head(top_n)
        for rank, (_, row) in enumerate(candidates.iterrows(), start=1):
            neighbor_records.append(
                {
                    "query_condition": query_condition,
                    "rank": rank,
                    "neighbor_profile_id": row["profile_id"],
                    "neighbor_condition": row["condition"],
                    "neighbor_is_other_target_condition": row["condition"]
                    in set(TARGET_CONDITIONS) - {query_condition},
                    "cosine_similarity": row["cosine_similarity"],
                    "study_id": row["study_id"],
                    "experiment_id": row["experiment_id"],
                    "direction": row["direction_label"],
                    "study_title": row["title"],
                    "experiment_label": row["experiment_label"],
                    "body_site": row["body_site"],
                    "location": row["location"],
                    "profile_trait_count": row["profile_trait_count"],
                    "total_trait_signal": row["total_trait_signal"],
                }
            )
    neighbors = pd.DataFrame(neighbor_records)

    crosscheck_records: list[dict[str, object]] = []
    for query_condition in TARGET_CONDITIONS:
        hits = neighbors[
            (neighbors["query_condition"] == query_condition)
            & (neighbors["neighbor_is_other_target_condition"])
        ].copy()
        crosscheck_records.append(
            {
                "query_condition": query_condition,
                "other_target_condition_found_in_top5_profiles": bool(len(hits)),
                "matched_other_target_conditions": join_unique(hits["neighbor_condition"])
                if len(hits)
                else "",
                "matched_neighbor_profile_ids": join_unique(hits["neighbor_profile_id"])
                if len(hits)
                else "",
            }
        )
    crosscheck = pd.DataFrame(crosscheck_records)

    centroid_neighbor_records: list[dict[str, object]] = []
    conditions = list(centroids)
    centroid_matrix = np.vstack([centroids[condition] for condition in conditions])
    for query_condition in TARGET_CONDITIONS:
        if query_condition not in centroids:
            continue
        scores = centroid_matrix @ centroids[query_condition]
        order = np.argsort(scores)[::-1]
        rank = 1
        for position in order:
            neighbor_condition = conditions[position]
            if neighbor_condition == query_condition:
                continue
            summary_row = condition_summary[
                condition_summary["condition"] == neighbor_condition
            ].iloc[0]
            centroid_neighbor_records.append(
                {
                    "query_condition": query_condition,
                    "rank": rank,
                    "neighbor_condition": neighbor_condition,
                    "neighbor_is_other_target_condition": neighbor_condition
                    in set(TARGET_CONDITIONS) - {query_condition},
                    "cosine_similarity": float(scores[position]),
                    "neighbor_profile_count": int(summary_row["profile_count"]),
                    "neighbor_study_count": int(summary_row["study_count"]),
                    "neighbor_experiment_count": int(summary_row["experiment_count"]),
                }
            )
            rank += 1
            if rank > top_n:
                break
    centroid_neighbors = pd.DataFrame(centroid_neighbor_records)
    return neighbors, crosscheck, centroid_neighbors


def parse_metaphlan_name(raw: str) -> tuple[str, str | None]:
    """Parse a MetaPhlAn leaf name such as s__Gardnerella vaginalis."""
    for prefix, rank in METAPHLAN_RANK_MAP.items():
        if raw.startswith(prefix):
            return raw[len(prefix) :], rank
    return raw, None


def infer_taxon_rank(name: str) -> str:
    """Infer a coarse rank when MetaPhlAn rank prefixes are unavailable."""
    parts = name.strip().split()
    if len(parts) >= 3:
        return "strain"
    if len(parts) == 2:
        return "species"
    return "genus"


def parse_taxa_from_raw_row(row: dict[str, str]) -> list[dict[str, str]]:
    """Parse most-specific taxon IDs and names from one BugSigDB full-dump row."""
    tax_ids = clean_optional_text(row.get("NCBI Taxonomy IDs"))
    if not tax_ids:
        return []

    id_chains = [chain.strip() for chain in tax_ids.split(";") if chain.strip()]
    names = clean_optional_text(row.get("MetaPhlAn taxon names")) or clean_optional_text(
        row.get("Taxon")
    )
    name_chains = [chain.strip() for chain in names.split(",") if chain.strip()] if names else []

    taxa: list[dict[str, str]] = []
    for index, id_chain in enumerate(id_chains):
        lineage_ids = [value.strip() for value in id_chain.split("|") if value.strip()]
        if not lineage_ids:
            continue

        leaf_id = lineage_ids[-1]
        if leaf_id.isdigit():
            taxon_id = f"NCBITaxon:{leaf_id}"
        elif leaf_id.startswith("NCBITaxon:"):
            taxon_id = leaf_id
        else:
            continue

        taxon_name = f"taxon_{leaf_id.removeprefix('NCBITaxon:')}"
        rank: str | None = None
        if index < len(name_chains):
            name_parts = [part.strip() for part in name_chains[index].split("|") if part.strip()]
            if name_parts:
                taxon_name, rank = parse_metaphlan_name(name_parts[-1])

        taxa.append(
            {
                "taxon_id": taxon_id,
                "taxon_name": taxon_name,
                "taxonomic_rank": rank or infer_taxon_rank(taxon_name),
            }
        )
    return taxa


def load_target_taxa_occurrences(full_dump_path: Path) -> pd.DataFrame:
    """Load target-condition taxon/profile occurrences from raw full_dump.csv."""
    records: list[dict[str, object]] = []
    with full_dump_path.open(encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(rows)
        for row in reader:
            condition = clean_text(row.get("Condition"))
            direction = clean_text(row.get("Abundance in Group 1")).lower()
            if condition not in TARGET_CONDITIONS or direction not in DIRECTION_TO_INC_DEC:
                continue
            profile_id = profile_id_from_raw_row(row)
            for taxon in parse_taxa_from_raw_row(row):
                records.append(
                    {
                        "condition": condition,
                        "profile_id": profile_id,
                        "study_id": study_id_from_bsdb_id(row.get("BSDB ID")),
                        "experiment_id": experiment_id_from_bsdb_id(row.get("BSDB ID")),
                        "direction": direction,
                        **taxon,
                    }
                )
    if not records:
        msg = "No target-condition taxon occurrences were found in full_dump.csv."
        raise ValueError(msg)
    return pd.DataFrame.from_records(records).drop_duplicates()


def select_heatmap_traits(
    values: pd.DataFrame,
    metadata: pd.DataFrame,
    trait_labels: pd.DataFrame,
    max_traits: int,
) -> list[str]:
    """Select infographic and target-condition traits for the companion heatmap."""
    target_profile_ids = metadata.index[metadata["condition"].isin(TARGET_CONDITIONS)]
    target_values = values.loc[target_profile_ids]
    support = target_values.gt(0).sum(axis=0).astype(int)
    total_signal = target_values.sum(axis=0).astype(float)
    condition_support: dict[str, int] = {}
    for trait_id in target_values.columns:
        condition_support[trait_id] = int(
            sum(
                target_values.loc[
                    metadata.loc[target_values.index, "condition"] == condition,
                    trait_id,
                ].sum()
                > 0
                for condition in TARGET_CONDITIONS
            )
        )
    trait_summary = pd.DataFrame(
        {
            "trait_id": target_values.columns,
            "target_profile_support": support,
            "target_total_signal": total_signal,
            "target_condition_support": pd.Series(condition_support),
        },
        index=target_values.columns,
    ).join(trait_labels[["display_name", "display_category"]], how="left")
    trait_summary = trait_summary[
        ~trait_summary.index.to_series().str.contains("mediadive.medium", regex=False)
    ].copy()
    trait_summary["selection_score"] = (
        trait_summary["target_condition_support"] * np.log1p(trait_summary["target_total_signal"])
        + trait_summary["target_profile_support"]
    )
    ordered = list(INFOGRAPHIC_TOP_TRAITS)
    ordered.extend(
        trait_summary.sort_values(
            [
                "selection_score",
                "target_condition_support",
                "target_profile_support",
                "target_total_signal",
            ],
            ascending=False,
        )
        .head(max_traits * 2)
        .index.tolist()
    )
    selected: list[str] = []
    for trait_id in ordered:
        if trait_id in trait_summary.index and trait_id not in selected:
            selected.append(trait_id)
        if len(selected) >= max_traits:
            break
    return selected


def summarize_target_taxa(occurrences: pd.DataFrame) -> pd.DataFrame:
    """Summarize target-condition organism occurrences."""
    grouped = (
        occurrences.groupby(["taxon_id", "taxon_name", "taxonomic_rank"], dropna=False)
        .agg(
            target_profile_count=("profile_id", "nunique"),
            target_study_count=("study_id", "nunique"),
            target_experiment_count=("experiment_id", "nunique"),
            target_condition_count=("condition", "nunique"),
            target_conditions=("condition", join_unique),
            target_directions=("direction", join_unique),
        )
        .reset_index()
        .set_index("taxon_id")
    )
    direction_counts = (
        occurrences[["taxon_id", "profile_id", "direction"]]
        .drop_duplicates()
        .pivot_table(
            index="taxon_id",
            columns="direction",
            values="profile_id",
            aggfunc="nunique",
            fill_value=0,
        )
    )
    for direction in ["increased", "decreased"]:
        if direction not in direction_counts:
            direction_counts[direction] = 0
    grouped = grouped.join(direction_counts[["increased", "decreased"]], how="left").fillna(0)
    grouped["dominant_direction"] = np.select(
        [
            grouped["increased"] > grouped["decreased"],
            grouped["decreased"] > grouped["increased"],
        ],
        ["increased", "decreased"],
        default="mixed",
    )
    condition_counts = (
        occurrences[["taxon_id", "profile_id", "condition"]]
        .drop_duplicates()
        .pivot_table(
            index="taxon_id",
            columns="condition",
            values="profile_id",
            aggfunc="nunique",
            fill_value=0,
        )
    )
    grouped["dominant_condition"] = [
        condition_counts.loc[taxon_id].sort_values(ascending=False).index[0]
        if taxon_id in condition_counts.index
        else "unknown"
        for taxon_id in grouped.index
    ]
    return grouped


def build_traits_vs_organisms_heatmap(
    values: pd.DataFrame,
    metadata: pd.DataFrame,
    trait_labels: pd.DataFrame,
    full_dump_path: Path,
    taxa_trait_matrix_path: Path,
    max_traits: int,
    max_taxa: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build the companion target-condition traits x organisms matrix."""
    selected_traits = select_heatmap_traits(values, metadata, trait_labels, max_traits)
    occurrences = load_target_taxa_occurrences(full_dump_path)
    preferred_occurrences = occurrences[
        occurrences["taxonomic_rank"].isin(PREFERRED_HEATMAP_RANKS)
    ].copy()
    if preferred_occurrences["taxon_id"].nunique() >= max(10, max_taxa // 2):
        occurrences = preferred_occurrences
    taxa_summary = summarize_target_taxa(occurrences)
    usecols = set(MATRIX_METADATA_COLUMNS + selected_traits)
    taxa_trait_matrix = pd.read_csv(
        taxa_trait_matrix_path,
        sep="\t",
        usecols=lambda column: column in usecols,
    )
    taxa_trait_matrix = taxa_trait_matrix[
        taxa_trait_matrix["taxon_id"].isin(taxa_summary.index)
    ].copy()
    taxa_trait_matrix = taxa_trait_matrix.drop_duplicates(subset=["taxon_id"], keep="first")
    taxa_trait_matrix = taxa_trait_matrix.set_index("taxon_id", drop=False)
    trait_values = (
        taxa_trait_matrix[selected_traits].apply(pd.to_numeric, errors="coerce").fillna(0)
    )
    binary = trait_values.gt(0).astype(int)
    taxa_metadata = taxa_summary.join(
        pd.DataFrame({"selected_trait_count": binary.sum(axis=1)}, index=binary.index),
        how="inner",
    )
    taxa_metadata = taxa_metadata[taxa_metadata["selected_trait_count"] > 0].copy()
    taxa_metadata = taxa_metadata.sort_values(
        [
            "target_condition_count",
            "target_profile_count",
            "selected_trait_count",
            "target_study_count",
            "taxon_name",
        ],
        ascending=[False, False, False, False, True],
    ).head(max_taxa)
    selected_matrix = binary.loc[taxa_metadata.index].copy()
    selected_traits_with_support = selected_matrix.columns[selected_matrix.sum(axis=0) > 0].tolist()
    selected_matrix = selected_matrix[selected_traits_with_support]
    selected_trait_metadata = trait_labels.loc[selected_traits_with_support].copy()
    return selected_matrix, taxa_metadata, selected_trait_metadata


def display_taxon_label(taxon_id: str, taxa_metadata: pd.DataFrame) -> str:
    """Create a compact organism heatmap row label."""
    row = taxa_metadata.loc[taxon_id]
    name = clean_text(row["taxon_name"])
    if name.startswith("taxon_"):
        name = taxon_id
    label = f"{name} ({int(row['target_profile_count'])} profiles)"
    return compact_label(label, 64)


def save_traits_vs_organisms_heatmap(
    matrix: pd.DataFrame,
    taxa_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
) -> pd.DataFrame:
    """Render and save the companion traits x organisms heatmap."""
    sns.set_theme(style="white")
    display_matrix = matrix.copy()
    display_matrix.index = [
        display_taxon_label(taxon_id, taxa_metadata) for taxon_id in display_matrix.index
    ]
    display_matrix.columns = [
        trait_metadata.loc[trait_id, "heatmap_label"] for trait_id in display_matrix.columns
    ]

    condition_palette = {
        "Thalassemia": "#b91c1c",
        "Sleep Disorder": "#2563eb",
        "Head and neck squamous cell carcinoma": "#047857",
        "unknown": "#64748b",
    }
    row_taxon_ids = matrix.index.tolist()
    row_colors = taxa_metadata.loc[row_taxon_ids, "dominant_condition"].map(
        lambda condition: condition_palette.get(condition, "#64748b")
    )
    row_colors.index = display_matrix.index
    row_colors.name = ""

    cluster = sns.clustermap(
        display_matrix,
        cmap=sns.color_palette(["#f8fafc", "#155e75"], as_cmap=True),
        method="average",
        metric="jaccard",
        row_colors=row_colors,
        linewidths=0.2,
        linecolor="#e2e8f0",
        xticklabels=True,
        yticklabels=True,
        figsize=(18, 16),
        cbar_pos=(0.02, 0.82, 0.025, 0.12),
        dendrogram_ratio=(0.12, 0.10),
        colors_ratio=(0.025, 0.025),
    )
    cluster.ax_heatmap.set_xlabel("Trait")
    cluster.ax_heatmap.set_ylabel("Organism observed in target-condition profiles")
    cluster.ax_heatmap.set_title("")
    cluster.fig.subplots_adjust(top=0.88)
    cluster.fig.suptitle(
        "Companion Heatmap: Traits vs Organisms in Infographic Outlier Conditions\n"
        "Rows are organisms observed in target BugSigDB profiles; columns are binary "
        "KG-Microbe trait annotations",
        y=0.98,
        fontsize=12,
    )
    cluster.ax_cbar.set_ylabel("Trait present")
    cluster.ax_row_colors.set_ylabel("")
    for label in cluster.ax_heatmap.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(7)
    for label in cluster.ax_heatmap.get_yticklabels():
        label.set_fontsize(7)

    handles = [
        Patch(color=color, label=condition)
        for condition, color in condition_palette.items()
        if condition != "unknown"
    ]
    cluster.ax_row_dendrogram.legend(
        handles=handles,
        title="Dominant target condition",
        loc="center",
        bbox_to_anchor=(-0.08, 0.5),
        frameon=False,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster.fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(cluster.fig)
    return display_matrix


def summarize_guild_connections(values: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Summarize how predefined guild trait sets score target-condition profiles."""
    records: list[dict[str, object]] = []
    target_metadata = metadata[metadata["condition"].isin(TARGET_CONDITIONS)].copy()
    for guild_name, traits in GUILD_TRAITS.items():
        available_traits = [trait_id for trait_id in traits if trait_id in values.columns]
        if not available_traits:
            continue
        guild_values = values.loc[target_metadata.index, available_traits].copy()
        row_totals = values.loc[target_metadata.index].sum(axis=1).replace(0, np.nan)
        relative_values = guild_values.div(row_totals, axis=0).fillna(0)
        for condition, group in target_metadata.groupby("condition", sort=False):
            condition_values = guild_values.loc[group.index]
            condition_relative = relative_values.loc[group.index]
            increased = group[group["inc_dec"] == "inc"]
            decreased = group[group["inc_dec"] == "dec"]
            records.append(
                {
                    "guild": guild_name,
                    "guild_traits": "; ".join(traits.values()),
                    "condition": condition,
                    "profile_count": len(group),
                    "study_count": group["study_id"].nunique(),
                    "experiment_count": group["experiment_id"].nunique(),
                    "mean_guild_trait_signal": condition_values.mean(axis=1).mean(),
                    "mean_relative_guild_signal": condition_relative.mean(axis=1).mean(),
                    "increased_mean_guild_trait_signal": condition_values.loc[increased.index]
                    .mean(axis=1)
                    .mean()
                    if len(increased)
                    else np.nan,
                    "decreased_mean_guild_trait_signal": condition_values.loc[decreased.index]
                    .mean(axis=1)
                    .mean()
                    if len(decreased)
                    else np.nan,
                }
            )
    return pd.DataFrame(records)


def write_report(
    output_path: Path,
    condition_depth: pd.DataFrame,
    titles: pd.DataFrame,
    autism: pd.DataFrame,
    profile_neighbors: pd.DataFrame,
    centroid_neighbors: pd.DataFrame,
    crosscheck: pd.DataFrame,
    guild_summary: pd.DataFrame,
    heatmap_matrix: pd.DataFrame,
) -> None:
    """Write a compact Markdown report for the follow-up artifacts."""
    target_depth = condition_depth[condition_depth["condition"].isin(TARGET_CONDITIONS)].copy()
    target_depth = target_depth.sort_values(
        "infographic_avg_enrichment_depth",
        ascending=False,
    )

    lines: list[str] = [
        "# Enrichment Infographic Follow-Up",
        "",
        "Source figure: "
        "https://github.com/CultureBotAI/bsdbng/blob/"
        "enrichment-analysis-heatmap-cluster/data/derived/enrichment_infographic.png",
        "",
        "The referenced infographic ranks the following three conditions as the top "
        "condition-level enrichment-depth outliers. The local follow-up uses the full "
        "`data/derived/study_trait_profiles_wide.tsv` matrix and raw BugSigDB experiment "
        "metadata to list titles, find similar profiles, and connect enriched traits to "
        "organisms.",
        "",
        "## Target Conditions",
        "",
        "| condition | infographic avg enrichment depth | local profiles | studies | "
        "experiments | local mean nonzero traits/profile | local percentile | local outlier? |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for _, row in target_depth.iterrows():
        lines.append(
            "| {condition} | {infographic:.2f} | {profiles:d} | {studies:d} | "
            "{experiments:d} | {mean_traits:.1f} | {percentile:.1%} | {outlier} |".format(
                condition=row["condition"],
                infographic=float(row["infographic_avg_enrichment_depth"]),
                profiles=int(row["profile_count"]),
                studies=int(row["study_count"]),
                experiments=int(row["experiment_count"]),
                mean_traits=float(row["mean_profile_trait_count"]),
                percentile=float(row["depth_percentile"]),
                outlier="yes" if bool(row["current_profile_depth_outlier"]) else "no",
            )
        )

    lines.extend(
        [
            "",
            "Interpretation: the infographic itself identifies these as high average "
            "enrichment-depth conditions. In this local profile-trait-count proxy, "
            "condition support matters: Thalassemia is concentrated in one study and two "
            "experiments, Sleep Disorder spans two studies and five experiments, and head "
            "and neck squamous cell carcinoma spans more experiments but has a mix of "
            "saliva/oral-cavity/fecal comparisons.",
            "",
            "## Study And Experiment Titles",
            "",
            "| condition | study | year | experiment | profiles |",
            "| --- | --- | ---: | --- | ---: |",
        ]
    )
    for _, row in titles.iterrows():
        lines.append(
            f"| {row['condition']} | {row['title']} ({row['study_id']}) | "
            f"{row['publication_year']} | {row['experiment_label']} | "
            f"{int(row['profile_count'])} |"
        )

    lines.extend(
        [
            "",
            "## Autism Label Check",
            "",
            "| exact label | EFO IDs | profiles | studies | experiments | interpretation |",
            "| --- | --- | ---: | ---: | ---: | --- |",
        ]
    )
    for _, row in autism.iterrows():
        lines.append(
            f"| {row['exact_condition_label']} | {row['efo_ids']} | "
            f"{int(row['profile_count'])} | {int(row['study_count'])} | "
            f"{int(row['experiment_count'])} | {row['interpretation']} |"
        )
    lines.extend(
        [
            "",
            "Autism and autism spectrum disorder should not be silently treated as the "
            "same exact condition label. For a clinical/domain rollup, they can be grouped "
            "as `Autism / autism spectrum disorder` if the merge is declared, but exact "
            "label and EFO-level analyses should keep them separate.",
            "",
            "## Similar Profiles",
            "",
            "Profile similarity uses cosine similarity over log1p-transformed full trait "
            "vectors. For each target condition, the query vector is the centroid of its "
            "active profiles, and the top five individual profiles exclude the same exact "
            "condition label.",
            "",
            "| query condition | rank | neighbor condition | similarity | study title | "
            "experiment |",
            "| --- | ---: | --- | ---: | --- | --- |",
        ]
    )
    for _, row in profile_neighbors.iterrows():
        lines.append(
            f"| {row['query_condition']} | {int(row['rank'])} | "
            f"{row['neighbor_condition']} | {float(row['cosine_similarity']):.3f} | "
            f"{row['study_title']} | {row['experiment_label']} |"
        )

    lines.extend(
        [
            "",
            "### Cross-Target Check",
            "",
            "| query condition | any other target in top 5 profiles? | matched target conditions |",
            "| --- | --- | --- |",
        ]
    )
    for _, row in crosscheck.iterrows():
        lines.append(
            f"| {row['query_condition']} | "
            f"{'yes' if bool(row['other_target_condition_found_in_top5_profiles']) else 'no'} | "
            f"{row['matched_other_target_conditions']} |"
        )

    lines.extend(
        [
            "",
            "Top five condition-centroid neighbors are also written to "
            "`target_condition_centroid_neighbors.tsv`; this checks whether the three "
            "conditions find each other at condition rather than individual-profile level.",
            "",
            "| query condition | rank | centroid neighbor | similarity | other target? |",
            "| --- | ---: | --- | ---: | --- |",
        ]
    )
    for _, row in centroid_neighbors.iterrows():
        lines.append(
            f"| {row['query_condition']} | {int(row['rank'])} | "
            f"{row['neighbor_condition']} | {float(row['cosine_similarity']):.3f} | "
            f"{'yes' if bool(row['neighbor_is_other_target_condition']) else 'no'} |"
        )

    lines.extend(
        [
            "",
            "## How Guilds Connect To Profiles",
            "",
            "A profile is one BugSigDB study/experiment/direction row. Its trait vector "
            "counts KG-Microbe trait annotations among organisms reported as increased "
            "or decreased in that experiment. A guild is a named set of biologically "
            "related traits, such as `aerobe_v2` or `anaerobe_v2_depletion`. A profile "
            "connects to a guild when the guild's trait columns carry signal in that "
            "profile; organism support is then checked by intersecting organisms observed "
            "in the same profile with organism-trait annotations from "
            "`data/derived/taxa_trait_matrix.tsv`.",
            "",
            "| guild | condition | profiles | mean guild signal | mean relative guild signal |",
            "| --- | --- | ---: | ---: | ---: |",
        ]
    )
    for _, row in guild_summary.iterrows():
        lines.append(
            f"| {row['guild']} | {row['condition']} | {int(row['profile_count'])} | "
            f"{float(row['mean_guild_trait_signal']):.3f} | "
            f"{float(row['mean_relative_guild_signal']):.5f} |"
        )

    lines.extend(
        [
            "",
            "## Companion Traits Vs Organisms Heatmap",
            "",
            f"`traits_vs_organisms_heatmap.png` plots {len(heatmap_matrix)} organisms by "
            f"{len(heatmap_matrix.columns)} traits. Rows are organisms observed in the "
            "target-condition BugSigDB profiles; columns are binary KG-Microbe trait "
            "annotations selected from the infographic top traits plus high-signal "
            "target-condition traits.",
            "",
            "## Files",
            "",
            "- `outlier_condition_study_experiment_titles.tsv`",
            "- `outlier_condition_depth_summary.tsv`",
            "- `autism_condition_label_check.tsv`",
            "- `outlier_condition_profile_neighbors.tsv`",
            "- `outlier_condition_neighbor_crosscheck.tsv`",
            "- `target_condition_centroid_neighbors.tsv`",
            "- `guild_profile_connection_summary.tsv`",
            "- `traits_vs_organisms_heatmap_matrix.tsv`",
            "- `traits_vs_organisms_taxa_metadata.tsv`",
            "- `traits_vs_organisms_trait_metadata.tsv`",
            "- `traits_vs_organisms_heatmap.png`",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Analyze the three condition outliers from the enrichment infographic and "
            "create a companion traits-vs-organisms heatmap."
        )
    )
    parser.add_argument(
        "--profiles-input",
        type=Path,
        default=REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv",
        help="Input wide study trait profiles TSV path.",
    )
    parser.add_argument(
        "--raw-experiments-input",
        type=Path,
        default=REPO_ROOT / "data/raw/experiments.csv",
        help="Input raw BugSigDB experiments CSV path.",
    )
    parser.add_argument(
        "--full-dump-input",
        type=Path,
        default=REPO_ROOT / "data/raw/full_dump.csv",
        help="Input raw BugSigDB full dump CSV path.",
    )
    parser.add_argument(
        "--studies-input",
        type=Path,
        default=REPO_ROOT / "data/derived/studies.tsv",
        help="Input derived studies TSV path.",
    )
    parser.add_argument(
        "--taxa-trait-matrix-input",
        type=Path,
        default=REPO_ROOT / "data/derived/taxa_trait_matrix.tsv",
        help="Input taxa trait matrix TSV path.",
    )
    parser.add_argument(
        "--taxa-trait-columns-input",
        type=Path,
        default=REPO_ROOT / "data/derived/taxa_trait_matrix_columns.tsv",
        help="Input taxa trait column metadata TSV path.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Output directory for report artifacts.",
    )
    parser.add_argument(
        "--top-neighbors",
        type=int,
        default=5,
        help="Number of nearest profiles and condition centroids to report.",
    )
    parser.add_argument(
        "--max-heatmap-traits",
        type=int,
        default=34,
        help="Maximum number of trait columns in the companion heatmap.",
    )
    parser.add_argument(
        "--max-heatmap-taxa",
        type=int,
        default=45,
        help="Maximum number of organism rows in the companion heatmap.",
    )
    parser.add_argument(
        "--min-condition-profiles",
        type=int,
        default=2,
        help="Minimum active profiles required for condition centroid comparisons.",
    )
    return parser.parse_args()


def main() -> None:
    """Run the enrichment infographic follow-up analysis."""
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    experiments = load_raw_experiments(args.raw_experiments_input)
    studies = load_studies(args.studies_input)
    trait_labels = load_trait_labels(args.taxa_trait_columns_input)
    values, metadata, _ = load_profiles(args.profiles_input, experiments, studies)

    condition_depth = summarize_condition_depth(metadata, args.min_condition_profiles)
    titles = write_target_title_table(
        metadata,
        experiments,
        studies,
        args.output_dir / "outlier_condition_study_experiment_titles.tsv",
    )
    autism = summarize_autism_labels(metadata, experiments)
    profile_neighbors, crosscheck, centroid_neighbors = compute_profile_neighbors(
        values,
        metadata,
        args.top_neighbors,
        args.min_condition_profiles,
    )
    guild_summary = summarize_guild_connections(values, metadata)
    organism_matrix, taxa_metadata, trait_metadata = build_traits_vs_organisms_heatmap(
        values,
        metadata,
        trait_labels,
        args.full_dump_input,
        args.taxa_trait_matrix_input,
        args.max_heatmap_traits,
        args.max_heatmap_taxa,
    )
    display_matrix = save_traits_vs_organisms_heatmap(
        organism_matrix,
        taxa_metadata,
        trait_metadata,
        args.output_dir / "traits_vs_organisms_heatmap.png",
    )

    condition_depth.to_csv(
        args.output_dir / "outlier_condition_depth_summary.tsv",
        sep="\t",
        index=False,
    )
    autism.to_csv(args.output_dir / "autism_condition_label_check.tsv", sep="\t", index=False)
    profile_neighbors.to_csv(
        args.output_dir / "outlier_condition_profile_neighbors.tsv",
        sep="\t",
        index=False,
    )
    crosscheck.to_csv(
        args.output_dir / "outlier_condition_neighbor_crosscheck.tsv",
        sep="\t",
        index=False,
    )
    centroid_neighbors.to_csv(
        args.output_dir / "target_condition_centroid_neighbors.tsv",
        sep="\t",
        index=False,
    )
    guild_summary.to_csv(
        args.output_dir / "guild_profile_connection_summary.tsv",
        sep="\t",
        index=False,
    )
    display_matrix.to_csv(
        args.output_dir / "traits_vs_organisms_heatmap_matrix.tsv",
        sep="\t",
    )
    taxa_metadata.to_csv(
        args.output_dir / "traits_vs_organisms_taxa_metadata.tsv",
        sep="\t",
    )
    trait_metadata.to_csv(
        args.output_dir / "traits_vs_organisms_trait_metadata.tsv",
        sep="\t",
    )
    write_report(
        args.output_dir / "README.md",
        condition_depth,
        titles,
        autism,
        profile_neighbors,
        centroid_neighbors,
        crosscheck,
        guild_summary,
        organism_matrix,
    )

    print(f"Wrote report: {args.output_dir / 'README.md'}")
    print(f"Wrote heatmap: {args.output_dir / 'traits_vs_organisms_heatmap.png'}")
    print(f"Target title rows: {len(titles)}")
    print(f"Profile neighbor rows: {len(profile_neighbors)}")
    print(f"Heatmap shape: {organism_matrix.shape[0]} taxa x {organism_matrix.shape[1]} traits")


if __name__ == "__main__":
    main()
