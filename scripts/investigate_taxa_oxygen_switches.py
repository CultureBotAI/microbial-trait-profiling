"""Investigate disease-state taxon direction switches by oxygen tolerance."""

from __future__ import annotations

from pathlib import Path
from textwrap import shorten

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from create_condition_body_site_guild_maps import (
    FACULTATIVE_CONFIG,
    v3_body_site_display_label,
    v3_fully_grouped_condition_display_label,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
PROFILE_METADATA_PATH = REPO_ROOT / "trait-visualization/study_trait_profile_study_metadata.tsv"
TAXA_PATH = REPO_ROOT / "data/derived/taxa.tsv"
TAXA_TRAIT_MATRIX_PATH = REPO_ROOT / "data/derived/taxa_trait_matrix.tsv"
OUTPUT_DIR = REPO_ROOT / "trait-visualization/trait-profile-investigation/taxa-oxygen-switches"

AEROBIC_TRAITS = {
    "METPO:1000602": "aerobic",
    "METPO:1000604": "microaerophilic",
    "METPO:1000606": "obligately aerobic",
    "METPO:1000608": "facultatively aerobic",
}
ANAEROBIC_TRAITS = {
    "METPO:1000603": "anaerobic",
    "METPO:1000607": "obligately anaerobic",
    "METPO:1000611": "strictly anaerobic",
}

DISPLAY_RANKS = {"genus", "species"}
MIN_TAXON_CONDITION_PROFILES = 2
MIN_SWITCH_SIDE_PROFILES = 3
MIN_SWITCH_SIDE_STUDIES = 2
MIN_SWITCH_SIDE_CONDITIONS = 1
MIN_CONDITION_CLASS_OCCURRENCES = 12
MIN_CONDITION_STUDIES = 2
TOP_ROWS = 15

NON_DISEASE_EXACT_LABELS = {
    "Age",
    "Aging",
    "Antimicrobial agent",
    "Breastfeeding duration",
    "Cesarean section",
    "Clinical treatment",
    "Delivery method",
    "Diet",
    "Ethnic group",
    "Exposure",
    "Glucose",
    "Health study participation",
    "Health trait",
    "High altitude adaptation",
    "Microbiome",
    "Physical activity",
    "Population",
    "Pregnancy",
    "Species design",
    "Time",
    "Timepoint",
    "Treatment",
    "Tumor stage",
}
NON_DISEASE_TERMS = (
    "measurement",
    "response to ",
    "smoking behavior",
    "smoking behaviour",
    "smoking status",
    "waist circumference",
    "body fat distribution",
    "library preparation",
    "growth condition",
)


def clean_text(value: object) -> str:
    """Normalize empty display values."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def join_unique(values: pd.Series) -> str:
    """Join unique values in first-seen order."""
    seen: list[str] = []
    for value in values:
        text = clean_text(value)
        if text != "unknown" and text not in seen:
            seen.append(text)
    return "; ".join(seen)


def is_disease_state_condition(value: object) -> bool:
    """Return whether a BugSigDB condition label is disease-like enough for this report."""
    label = clean_text(value)
    lower = label.lower()
    if label == "unknown" or label == "Healthy":
        return False
    if label in NON_DISEASE_EXACT_LABELS:
        return False
    return not any(term in lower for term in NON_DISEASE_TERMS)


def profile_id_from_taxa_row(study_id: object, experiment_id: object, direction: object) -> str:
    """Build the profile ID used by the profile metadata table."""
    direction_text = clean_text(direction)
    inc_dec = "inc" if direction_text == "increased" else "dec"
    return f"{study_id}|{experiment_id}|{inc_dec}"


def load_profile_metadata() -> pd.DataFrame:
    """Load human-host profile metadata from the existing trait visualization artifacts."""
    metadata = pd.read_csv(PROFILE_METADATA_PATH, sep="\t", index_col=0)
    metadata = metadata.copy()
    metadata["profile_id"] = metadata.index
    for column in [
        "study_id",
        "experiment_id",
        "inc_dec",
        "direction_label",
        "condition",
        "body_site",
        "body_site_group",
        "host_species",
        "title",
    ]:
        metadata[column] = metadata[column].map(clean_text)
    metadata = metadata[
        (metadata["host_species"] == "Homo sapiens")
        & metadata["inc_dec"].isin({"inc", "dec"})
        & metadata["condition"].map(is_disease_state_condition)
    ].copy()

    condition_pairs = metadata["condition"].map(v3_fully_grouped_condition_display_label)
    metadata["condition_display"] = [label for label, _rule in condition_pairs]
    metadata["condition_rollup_rule"] = [rule for _label, rule in condition_pairs]

    body_pairs = [
        v3_body_site_display_label(row.body_site, row.body_site_group, FACULTATIVE_CONFIG)
        for row in metadata[["body_site", "body_site_group"]].itertuples(index=False)
    ]
    metadata["body_site_display"] = [label for label, _rule in body_pairs]
    metadata["body_site_rollup_rule"] = [rule for _label, rule in body_pairs]
    return metadata


def classify_oxygen(row: pd.Series) -> str:
    """Classify a taxon by direct oxygen-tolerance trait annotations."""
    aerobic_signal = int(row["aerobic_trait_signal"])
    anaerobic_signal = int(row["anaerobic_trait_signal"])
    if aerobic_signal > 0 and anaerobic_signal == 0:
        return "aerobic_only"
    if anaerobic_signal > 0 and aerobic_signal == 0:
        return "anaerobic_only"
    if aerobic_signal > 0 and anaerobic_signal > 0:
        return "mixed_oxygen_annotation"
    return "unannotated"


def trait_names_for_row(row: pd.Series, traits: dict[str, str]) -> str:
    """List oxygen traits present for one taxon."""
    names = [name for trait_id, name in traits.items() if int(row[trait_id]) > 0]
    return "; ".join(names)


def load_oxygen_annotations() -> pd.DataFrame:
    """Load taxon oxygen-tolerance annotations from the taxa trait matrix."""
    trait_ids = [*AEROBIC_TRAITS, *ANAEROBIC_TRAITS]
    annotations = pd.read_csv(
        TAXA_TRAIT_MATRIX_PATH,
        sep="\t",
        usecols=["taxon_id", "taxon_name", "taxonomic_rank", *trait_ids],
    )
    for trait_id in trait_ids:
        annotations[trait_id] = (
            pd.to_numeric(
                annotations[trait_id],
                errors="coerce",
            )
            .fillna(0)
            .astype(int)
        )
    annotations["aerobic_trait_signal"] = annotations[list(AEROBIC_TRAITS)].sum(axis=1)
    annotations["anaerobic_trait_signal"] = annotations[list(ANAEROBIC_TRAITS)].sum(axis=1)
    annotations["oxygen_class"] = annotations.apply(classify_oxygen, axis=1)
    annotations["aerobic_traits"] = annotations.apply(
        lambda row: trait_names_for_row(row, AEROBIC_TRAITS),
        axis=1,
    )
    annotations["anaerobic_traits"] = annotations.apply(
        lambda row: trait_names_for_row(row, ANAEROBIC_TRAITS),
        axis=1,
    )
    return annotations[
        [
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
            "oxygen_class",
            "aerobic_trait_signal",
            "anaerobic_trait_signal",
            "aerobic_traits",
            "anaerobic_traits",
        ]
    ]


def load_taxa_occurrences(metadata: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    """Load unique taxon/profile occurrences and attach metadata plus oxygen class."""
    taxa = pd.read_csv(
        TAXA_PATH,
        sep="\t",
        usecols=[
            "study_id",
            "experiment_id",
            "signature_id",
            "direction",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
        ],
    )
    taxa["profile_id"] = [
        profile_id_from_taxa_row(study_id, experiment_id, direction)
        for study_id, experiment_id, direction in taxa[
            ["study_id", "experiment_id", "direction"]
        ].itertuples(index=False, name=None)
    ]
    taxa = taxa[taxa["taxonomic_rank"].isin(DISPLAY_RANKS)].copy()
    taxa = taxa.drop_duplicates(["profile_id", "taxon_id"])

    metadata_columns = [
        "profile_id",
        "study_id",
        "experiment_id",
        "inc_dec",
        "direction_label",
        "condition",
        "condition_display",
        "condition_rollup_rule",
        "body_site",
        "body_site_group",
        "body_site_display",
        "title",
        "publication_year",
        "host_species",
    ]
    occurrences = taxa.merge(
        metadata[metadata_columns],
        on="profile_id",
        how="inner",
        suffixes=("_taxa", ""),
    )
    occurrences = occurrences.merge(
        annotations.drop(columns=["taxon_name", "taxonomic_rank"]),
        on="taxon_id",
        how="left",
    )
    occurrences["oxygen_class"] = occurrences["oxygen_class"].fillna("unannotated")
    for column in ["aerobic_traits", "anaerobic_traits"]:
        occurrences[column] = occurrences[column].fillna("")
    for column in ["aerobic_trait_signal", "anaerobic_trait_signal"]:
        occurrences[column] = (
            pd.to_numeric(
                occurrences[column],
                errors="coerce",
            )
            .fillna(0)
            .astype(int)
        )
    return occurrences


def direction_counts(group: pd.DataFrame, direction: str) -> int:
    """Count unique taxon-profile occurrences in one direction."""
    return int(group.loc[group["direction"] == direction, ["profile_id", "taxon_id"]].shape[0])


def direction_study_count(group: pd.DataFrame, direction: str) -> int:
    """Count studies in one direction."""
    return int(group.loc[group["direction"] == direction, "study_id"].nunique())


def summarize_condition_oxygen(occurrences: pd.DataFrame) -> pd.DataFrame:
    """Summarize oxygen-class direction balances by disease condition and body site."""
    annotated = occurrences[
        occurrences["oxygen_class"].isin({"aerobic_only", "anaerobic_only"})
    ].copy()
    grouped = (
        annotated.groupby(
            [
                "condition_display",
                "condition",
                "body_site_display",
                "oxygen_class",
            ],
            dropna=False,
            sort=False,
        )
        .agg(
            taxon_profile_count=("profile_id", "count"),
            study_count=("study_id", "nunique"),
            taxon_count=("taxon_id", "nunique"),
            increased_taxon_profile_count=(
                "direction",
                lambda values: int((values == "increased").sum()),
            ),
            decreased_taxon_profile_count=(
                "direction",
                lambda values: int((values == "decreased").sum()),
            ),
            increased_study_count=(
                "direction",
                lambda values: int(
                    annotated.loc[values.index][
                        annotated.loc[values.index, "direction"] == "increased"
                    ]["study_id"].nunique()
                ),
            ),
            decreased_study_count=(
                "direction",
                lambda values: int(
                    annotated.loc[values.index][
                        annotated.loc[values.index, "direction"] == "decreased"
                    ]["study_id"].nunique()
                ),
            ),
            raw_body_sites=("body_site", join_unique),
        )
        .reset_index()
    )
    grouped["direction_balance"] = (
        grouped["increased_taxon_profile_count"] - grouped["decreased_taxon_profile_count"]
    ) / grouped["taxon_profile_count"].replace(0, np.nan)
    return grouped.sort_values(
        ["condition_display", "condition", "body_site_display", "oxygen_class"],
    )


def build_condition_switch_scores(condition_oxygen: pd.DataFrame) -> pd.DataFrame:
    """Build condition/body-site scores for reciprocal aerobic/anaerobic direction shifts."""
    support = condition_oxygen[
        (condition_oxygen["taxon_profile_count"] >= MIN_CONDITION_CLASS_OCCURRENCES)
        & (condition_oxygen["study_count"] >= MIN_CONDITION_STUDIES)
        & condition_oxygen["oxygen_class"].isin({"aerobic_only", "anaerobic_only"})
    ].copy()
    pivot = support.pivot_table(
        index=["condition_display", "condition", "body_site_display"],
        columns="oxygen_class",
        values=[
            "taxon_profile_count",
            "study_count",
            "taxon_count",
            "increased_taxon_profile_count",
            "decreased_taxon_profile_count",
            "direction_balance",
        ],
        aggfunc="first",
    )
    if pivot.empty:
        return pd.DataFrame()

    pivot.columns = [f"{oxygen_class}_{metric}" for metric, oxygen_class in pivot.columns]
    output = pivot.reset_index()
    for oxygen_class in ["aerobic_only", "anaerobic_only"]:
        for metric in [
            "taxon_profile_count",
            "study_count",
            "taxon_count",
            "increased_taxon_profile_count",
            "decreased_taxon_profile_count",
            "direction_balance",
        ]:
            column = f"{oxygen_class}_{metric}"
            if column not in output:
                output[column] = np.nan
    output = output.dropna(
        subset=[
            "aerobic_only_direction_balance",
            "anaerobic_only_direction_balance",
        ]
    ).copy()
    output["aerobe_up_anaerobe_down_score"] = (
        output["aerobic_only_direction_balance"] - output["anaerobic_only_direction_balance"]
    )
    output["anaerobe_up_aerobe_down_score"] = -output["aerobe_up_anaerobe_down_score"]
    output["support_score"] = np.log1p(
        output["aerobic_only_taxon_profile_count"] + output["anaerobic_only_taxon_profile_count"]
    ) * np.log1p(output["aerobic_only_study_count"] + output["anaerobic_only_study_count"])
    output["rank_score"] = output["aerobe_up_anaerobe_down_score"].abs() * output["support_score"]
    return output.sort_values(
        ["rank_score", "support_score", "condition_display", "condition"],
        ascending=[False, False, True, True],
    )


def summarize_taxon_condition_direction(occurrences: pd.DataFrame) -> pd.DataFrame:
    """Summarize per-taxon direction support in each disease condition."""
    rows: list[dict[str, object]] = []
    for (
        condition_display,
        condition,
        taxon_id,
        taxon_name,
        taxonomic_rank,
        oxygen_class,
    ), group in occurrences.groupby(
        [
            "condition_display",
            "condition",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
            "oxygen_class",
        ],
        dropna=False,
        sort=False,
    ):
        increased = group[group["direction"] == "increased"]
        decreased = group[group["direction"] == "decreased"]
        profile_count = len(group)
        if profile_count < MIN_TAXON_CONDITION_PROFILES:
            continue
        increased_count = len(increased)
        decreased_count = len(decreased)
        if increased_count > decreased_count:
            dominant_direction = "increased"
        elif decreased_count > increased_count:
            dominant_direction = "decreased"
        else:
            dominant_direction = "tie"
        rows.append(
            {
                "condition_display": condition_display,
                "condition": condition,
                "taxon_id": taxon_id,
                "taxon_name": taxon_name,
                "taxonomic_rank": taxonomic_rank,
                "oxygen_class": oxygen_class,
                "profile_count": profile_count,
                "study_count": group["study_id"].nunique(),
                "increased_profile_count": increased_count,
                "decreased_profile_count": decreased_count,
                "increased_study_count": increased["study_id"].nunique(),
                "decreased_study_count": decreased["study_id"].nunique(),
                "dominant_direction": dominant_direction,
                "direction_balance": (increased_count - decreased_count) / profile_count,
                "body_site_displays": join_unique(group["body_site_display"]),
                "raw_body_sites": join_unique(group["body_site"]),
                "example_titles": join_unique(group["title"].head(3)),
                "aerobic_traits": join_unique(group["aerobic_traits"]),
                "anaerobic_traits": join_unique(group["anaerobic_traits"]),
            }
        )
    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary
    return summary.sort_values(
        ["profile_count", "study_count", "taxon_name", "condition"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)


def top_condition_text(
    rows: pd.DataFrame,
    direction: str,
    count_column: str,
    max_items: int = 4,
) -> str:
    """Describe top condition evidence for one switch side."""
    subset = rows[rows["dominant_direction"] == direction].copy()
    subset = subset[
        (subset[count_column] >= MIN_SWITCH_SIDE_PROFILES)
        & (subset["study_count"] >= MIN_SWITCH_SIDE_STUDIES)
    ]
    subset = subset.sort_values(
        [count_column, "study_count", "condition"],
        ascending=[False, False, True],
    )
    labels = [
        (
            f"{row['condition']} ({int(row[count_column])} profiles, "
            f"{int(row['study_count'])} studies)"
        )
        for _, row in subset.head(max_items).iterrows()
    ]
    return "; ".join(labels)


def summarize_taxon_switches(taxon_condition: pd.DataFrame) -> pd.DataFrame:
    """Find taxa whose dominant direction flips across disease conditions."""
    if taxon_condition.empty:
        return taxon_condition
    annotated = taxon_condition[
        taxon_condition["oxygen_class"].isin(
            {"aerobic_only", "anaerobic_only", "mixed_oxygen_annotation"}
        )
    ].copy()
    strong = annotated[
        (annotated["dominant_direction"].isin({"increased", "decreased"}))
        & (annotated["study_count"] >= MIN_SWITCH_SIDE_STUDIES)
        & (
            (
                (annotated["dominant_direction"] == "increased")
                & (annotated["increased_profile_count"] >= MIN_SWITCH_SIDE_PROFILES)
            )
            | (
                (annotated["dominant_direction"] == "decreased")
                & (annotated["decreased_profile_count"] >= MIN_SWITCH_SIDE_PROFILES)
            )
        )
    ].copy()

    records: list[dict[str, object]] = []
    for (
        taxon_id,
        taxon_name,
        taxonomic_rank,
        oxygen_class,
    ), group in strong.groupby(
        ["taxon_id", "taxon_name", "taxonomic_rank", "oxygen_class"],
        dropna=False,
        sort=False,
    ):
        increased_rows = group[group["dominant_direction"] == "increased"]
        decreased_rows = group[group["dominant_direction"] == "decreased"]
        if (
            len(increased_rows) < MIN_SWITCH_SIDE_CONDITIONS
            or len(decreased_rows) < MIN_SWITCH_SIDE_CONDITIONS
        ):
            continue
        total_profiles = int(group["profile_count"].sum())
        increased_profiles = int(increased_rows["increased_profile_count"].sum())
        decreased_profiles = int(decreased_rows["decreased_profile_count"].sum())
        switch_score = min(increased_profiles, decreased_profiles) * np.log1p(
            len(increased_rows) + len(decreased_rows)
        )
        records.append(
            {
                "taxon_id": taxon_id,
                "taxon_name": taxon_name,
                "taxonomic_rank": taxonomic_rank,
                "oxygen_class": oxygen_class,
                "switch_score": switch_score,
                "total_supported_profiles": total_profiles,
                "increased_dominant_condition_count": len(increased_rows),
                "decreased_dominant_condition_count": len(decreased_rows),
                "increased_side_profiles": increased_profiles,
                "decreased_side_profiles": decreased_profiles,
                "all_supported_condition_count": group["condition"].nunique(),
                "top_increased_conditions": top_condition_text(
                    group,
                    "increased",
                    "increased_profile_count",
                ),
                "top_decreased_conditions": top_condition_text(
                    group,
                    "decreased",
                    "decreased_profile_count",
                ),
                "aerobic_traits": join_unique(group["aerobic_traits"]),
                "anaerobic_traits": join_unique(group["anaerobic_traits"]),
                "body_site_displays": join_unique(group["body_site_displays"]),
            }
        )
    summary = pd.DataFrame(records)
    if summary.empty:
        return summary
    return summary.sort_values(
        [
            "switch_score",
            "total_supported_profiles",
            "increased_dominant_condition_count",
            "decreased_dominant_condition_count",
            "taxon_name",
        ],
        ascending=[False, False, False, False, True],
    ).reset_index(drop=True)


def summarize_within_condition_conflicts(taxon_condition: pd.DataFrame) -> pd.DataFrame:
    """Find condition/taxon rows with support in both increased and decreased profiles."""
    conflicts = taxon_condition[
        (taxon_condition["increased_profile_count"] >= MIN_SWITCH_SIDE_PROFILES)
        & (taxon_condition["decreased_profile_count"] >= MIN_SWITCH_SIDE_PROFILES)
    ].copy()
    if conflicts.empty:
        return conflicts
    conflicts["conflict_support"] = np.minimum(
        conflicts["increased_profile_count"],
        conflicts["decreased_profile_count"],
    )
    return conflicts.sort_values(
        ["conflict_support", "profile_count", "study_count", "taxon_name"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)


def summarize_global_oxygen(occurrences: pd.DataFrame) -> pd.DataFrame:
    """Summarize global direction support by oxygen class."""
    summary = (
        occurrences.groupby("oxygen_class", dropna=False, sort=False)
        .agg(
            taxon_profile_count=("profile_id", "count"),
            study_count=("study_id", "nunique"),
            taxon_count=("taxon_id", "nunique"),
            increased_taxon_profile_count=(
                "direction",
                lambda values: int((values == "increased").sum()),
            ),
            decreased_taxon_profile_count=(
                "direction",
                lambda values: int((values == "decreased").sum()),
            ),
        )
        .reset_index()
    )
    summary["direction_balance"] = (
        summary["increased_taxon_profile_count"] - summary["decreased_taxon_profile_count"]
    ) / summary["taxon_profile_count"].replace(0, np.nan)
    return summary.sort_values("taxon_profile_count", ascending=False)


def write_top_switch_plot(condition_scores: pd.DataFrame) -> None:
    """Write a compact plot of the strongest condition-level oxygen switches."""
    if condition_scores.empty:
        return
    positive = condition_scores.nlargest(TOP_ROWS, "aerobe_up_anaerobe_down_score")
    negative = condition_scores.nsmallest(TOP_ROWS, "aerobe_up_anaerobe_down_score")
    plot_data = pd.concat([positive, negative], ignore_index=True).drop_duplicates(
        ["condition", "body_site_display"]
    )
    plot_data = plot_data.sort_values("aerobe_up_anaerobe_down_score")
    labels = [
        shorten(f"{row.condition} | {row.body_site_display}", width=58, placeholder="...")
        for row in plot_data.itertuples(index=False)
    ]
    colors = np.where(plot_data["aerobe_up_anaerobe_down_score"] >= 0, "#dc2626", "#2563eb")

    fig_height = max(6, 0.32 * len(plot_data) + 1.8)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    ax.barh(labels, plot_data["aerobe_up_anaerobe_down_score"], color=colors)
    ax.axvline(0, color="#334155", linewidth=0.9)
    ax.set_xlabel("Aerobic-only direction balance minus anaerobic-only direction balance")
    ax.set_ylabel("")
    ax.set_title("Disease-State Oxygen Direction Switches")
    ax.text(
        0,
        1.01,
        "Positive: aerobic-only taxa skew increased while anaerobic-only taxa skew decreased. "
        "Negative: reciprocal pattern.",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color="#475569",
    )
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "condition_oxygen_switch_scores.png", dpi=180)
    plt.close(fig)


def markdown_table(frame: pd.DataFrame, columns: list[str], max_rows: int = 10) -> list[str]:
    """Render a small markdown table."""
    if frame.empty:
        return ["No rows passed the report filters."]
    subset = frame.loc[:, columns].head(max_rows).copy()
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in subset.itertuples(index=False):
        values = [clean_text(value) for value in row]
        lines.append("| " + " | ".join(values) + " |")
    return lines


def write_report(
    metadata: pd.DataFrame,
    occurrences: pd.DataFrame,
    global_oxygen: pd.DataFrame,
    condition_scores: pd.DataFrame,
    switch_summary: pd.DataFrame,
    conflicts: pd.DataFrame,
) -> None:
    """Write the investigation README."""
    top_aerobe_up = condition_scores.nlargest(TOP_ROWS, "aerobe_up_anaerobe_down_score").copy()
    top_anaerobe_up = condition_scores.nsmallest(TOP_ROWS, "aerobe_up_anaerobe_down_score").copy()
    for frame in [top_aerobe_up, top_anaerobe_up]:
        if not frame.empty:
            frame["aerobe_up_anaerobe_down_score"] = frame["aerobe_up_anaerobe_down_score"].map(
                lambda value: f"{value:.3f}"
            )
            frame["aerobic_only_direction_balance"] = frame["aerobic_only_direction_balance"].map(
                lambda value: f"{value:.3f}"
            )
            frame["anaerobic_only_direction_balance"] = frame[
                "anaerobic_only_direction_balance"
            ].map(lambda value: f"{value:.3f}")

    switch_display = switch_summary.copy()
    if not switch_display.empty:
        switch_display["switch_score"] = switch_display["switch_score"].map(
            lambda value: f"{value:.2f}"
        )

    lines = [
        "# Taxa Oxygen Switch Investigation",
        "",
        "## Scope",
        "",
        (
            f"- Human-host disease-like profile universe: {len(metadata):,} profiles, "
            f"{metadata['study_id'].nunique():,} studies."
        ),
        (f"- Genus/species taxon-profile occurrences after metadata join: {len(occurrences):,}."),
        (
            "- Oxygen classes are based on direct KG-Microbe trait columns: "
            f"aerobic side = {', '.join(AEROBIC_TRAITS.values())}; "
            f"anaerobic side = {', '.join(ANAEROBIC_TRAITS.values())}."
        ),
        (
            "- Disease-like filtering excludes Healthy/reference rows, broad cohort/design "
            "labels, measurements, exposures, diet/intervention response labels, and other "
            "non-disease metadata labels."
        ),
        "",
        "## Global Direction Balance by Oxygen Class",
        "",
        *markdown_table(
            global_oxygen.assign(
                direction_balance=global_oxygen["direction_balance"].map(
                    lambda value: f"{value:.3f}"
                )
            ),
            [
                "oxygen_class",
                "taxon_profile_count",
                "taxon_count",
                "study_count",
                "increased_taxon_profile_count",
                "decreased_taxon_profile_count",
                "direction_balance",
            ],
            max_rows=10,
        ),
        "",
        "## Strongest Condition-Level Oxygen Switches",
        "",
        (
            "Positive scores mean aerobic-only taxa skew increased while anaerobic-only taxa "
            "skew decreased."
        ),
        "",
        *markdown_table(
            top_aerobe_up,
            [
                "condition",
                "body_site_display",
                "aerobe_up_anaerobe_down_score",
                "aerobic_only_direction_balance",
                "anaerobic_only_direction_balance",
                "aerobic_only_taxon_profile_count",
                "anaerobic_only_taxon_profile_count",
            ],
            max_rows=8,
        ),
        "",
        (
            "Negative scores mean the reciprocal pattern: anaerobic-only taxa skew increased "
            "while aerobic-only taxa skew decreased."
        ),
        "",
        *markdown_table(
            top_anaerobe_up,
            [
                "condition",
                "body_site_display",
                "aerobe_up_anaerobe_down_score",
                "aerobic_only_direction_balance",
                "anaerobic_only_direction_balance",
                "aerobic_only_taxon_profile_count",
                "anaerobic_only_taxon_profile_count",
            ],
            max_rows=8,
        ),
        "",
        "## Top Cross-Disease Taxon Switches",
        "",
        (
            "These taxa have at least one disease condition where they are predominantly "
            "reported increased and at least one where they are predominantly reported decreased."
        ),
        "",
        *markdown_table(
            switch_display,
            [
                "taxon_name",
                "taxonomic_rank",
                "oxygen_class",
                "switch_score",
                "increased_dominant_condition_count",
                "decreased_dominant_condition_count",
                "top_increased_conditions",
                "top_decreased_conditions",
            ],
            max_rows=10,
        ),
        "",
        "## Within-Condition Direction Conflicts",
        "",
        (
            f"{len(conflicts):,} condition/taxon rows have at least "
            f"{MIN_SWITCH_SIDE_PROFILES} increased and {MIN_SWITCH_SIDE_PROFILES} decreased "
            "profile occurrences. These are not necessarily contradictions; many combine "
            "different experiments, cohorts, or body sites under the same condition label."
        ),
        "",
        "## Files",
        "",
        "- `global_oxygen_direction_summary.tsv`",
        "- `condition_oxygen_direction_summary.tsv`",
        "- `condition_oxygen_switch_scores.tsv`",
        "- `taxon_condition_direction_summary.tsv`",
        "- `taxa_cross_disease_switch_summary.tsv`",
        "- `taxa_within_condition_direction_conflicts.tsv`",
        "- `condition_oxygen_switch_scores.png`",
        "",
    ]
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run the investigation."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    metadata = load_profile_metadata()
    annotations = load_oxygen_annotations()
    occurrences = load_taxa_occurrences(metadata, annotations)

    global_oxygen = summarize_global_oxygen(occurrences)
    condition_oxygen = summarize_condition_oxygen(occurrences)
    condition_scores = build_condition_switch_scores(condition_oxygen)
    taxon_condition = summarize_taxon_condition_direction(occurrences)
    switch_summary = summarize_taxon_switches(taxon_condition)
    conflicts = summarize_within_condition_conflicts(taxon_condition)

    global_oxygen.to_csv(OUTPUT_DIR / "global_oxygen_direction_summary.tsv", sep="\t", index=False)
    condition_oxygen.to_csv(
        OUTPUT_DIR / "condition_oxygen_direction_summary.tsv",
        sep="\t",
        index=False,
    )
    condition_scores.to_csv(
        OUTPUT_DIR / "condition_oxygen_switch_scores.tsv",
        sep="\t",
        index=False,
    )
    taxon_condition.to_csv(
        OUTPUT_DIR / "taxon_condition_direction_summary.tsv",
        sep="\t",
        index=False,
    )
    switch_summary.to_csv(
        OUTPUT_DIR / "taxa_cross_disease_switch_summary.tsv",
        sep="\t",
        index=False,
    )
    conflicts.to_csv(
        OUTPUT_DIR / "taxa_within_condition_direction_conflicts.tsv",
        sep="\t",
        index=False,
    )
    write_top_switch_plot(condition_scores)
    write_report(metadata, occurrences, global_oxygen, condition_scores, switch_summary, conflicts)

    print(f"Wrote investigation: {OUTPUT_DIR}")
    print(f"Profiles: {len(metadata):,}")
    print(f"Taxon-profile occurrences: {len(occurrences):,}")
    print(f"Condition oxygen switch rows: {len(condition_scores):,}")
    print(f"Cross-disease taxon switches: {len(switch_summary):,}")
    print(f"Within-condition conflicts: {len(conflicts):,}")


if __name__ == "__main__":
    main()
