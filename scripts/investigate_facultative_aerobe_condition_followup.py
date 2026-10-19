"""Follow up aerobe v2 trait enrichment for autism and autoimmune labels."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle

from investigate_facultative_aerobe_expansion import (
    GUILD_TRAITS,
    REPO_ROOT,
    TARGET_TAXA,
    add_guild_scores,
    build_kg_microbe_taxon_universe,
    clean_text,
    load_study_heatmap,
    summarize_inferred_guild_evidence,
)
from investigate_facultative_aerobe_expansion import (
    OUTPUT_DIR as FACULTATIVE_OUTPUT_DIR,
)

OUTPUT_DIR = FACULTATIVE_OUTPUT_DIR
CONDITION_GROUPS = {
    "Autism": ("Autism",),
    "Autism spectrum disorder": ("Autism spectrum disorder",),
    "Autoimmune disease": ("Autoimmune disease",),
}
CONDITION_ORDER = [
    "Autism",
    "Autism spectrum disorder",
    "Autoimmune disease",
]
DIRECTION_ORDER = {"increased": 0, "decreased": 1}
HEATMAP_RANKS = {"genus", "species"}
HEATMAP_TOP_TAXA = 45
PATTERN_LABEL_TOP_N = 15
NON_GUILD_SPECIFIC_TAXA_HEATMAP_NOTE = (
    "Guild-specific: no - taxon direction-count heatmap; selected traits are annotation context."
)
PROFILE_SIMILARITY_PATTERN_NOTE = (
    "Colored x labels/strip: red = exact Autism and Autoimmune increased-skew; "
    "blue = broader ASD decreased-skew contrast; purple = both."
)
PATTERN_COLOR_EXACT_RED = "#b91c1c"
PATTERN_COLOR_BLUE = "#1d4ed8"
PATTERN_COLOR_BOTH = "#7c3aed"
PATTERN_COLOR_LABEL = "#111827"


def profile_id_from_parts(study_id: object, experiment_id: object, direction: object) -> str:
    """Build the profile ID used by the study trait profile matrix."""
    inc_dec = "inc" if clean_text(direction) == "increased" else "dec"
    return f"{study_id}|{experiment_id}|{inc_dec}"


def condition_group_lookup() -> dict[str, str]:
    """Map exact BugSigDB condition labels to follow-up groups."""
    lookup: dict[str, str] = {}
    for group, labels in CONDITION_GROUPS.items():
        lookup.update(dict.fromkeys(labels, group))
    return lookup


def load_scored_condition_profiles() -> pd.DataFrame:
    """Load scored profiles for the follow-up condition labels."""
    study_matrix, study_metadata, _ = load_study_heatmap()
    scored_profiles = add_guild_scores(study_matrix, study_metadata)
    lookup = condition_group_lookup()
    profiles = scored_profiles[scored_profiles["condition"].isin(lookup)].copy()
    profiles["condition_group_followup"] = profiles["condition"].map(lookup)
    profiles["profile_id"] = profiles.index
    return profiles


def summarize_profiles(profiles: pd.DataFrame) -> pd.DataFrame:
    """Summarize profile counts and trait-enrichment scores for each exact condition label."""
    rows: list[dict[str, object]] = []
    for condition, group in profiles.groupby("condition", sort=False):
        increased = group[group["direction_label"] == "increased"]
        decreased = group[group["direction_label"] == "decreased"]
        rows.append(
            {
                "condition_group": group["condition_group_followup"].iloc[0],
                "condition": condition,
                "profile_count": len(group),
                "study_count": group["study_id"].nunique(),
                "increased_profile_count": len(increased),
                "decreased_profile_count": len(decreased),
                "increased_study_count": increased["study_id"].nunique(),
                "decreased_study_count": decreased["study_id"].nunique(),
                "mean_guild_z": group["guild_z_mean"].mean(),
                "increased_mean_guild_z": increased["guild_z_mean"].mean(),
                "decreased_mean_guild_z": decreased["guild_z_mean"].mean(),
                "increased_minus_decreased_guild_z": increased["guild_z_mean"].mean()
                - decreased["guild_z_mean"].mean(),
                "body_site_groups": "; ".join(
                    sorted(group["body_site_group"].map(clean_text).unique())
                ),
            }
        )
    summary = pd.DataFrame(rows)
    summary["condition_order"] = summary["condition"].map(
        {condition: index for index, condition in enumerate(CONDITION_ORDER)}
    )
    return summary.sort_values("condition_order").drop(columns="condition_order")


def load_taxa_occurrences(profiles: pd.DataFrame) -> pd.DataFrame:
    """Load taxon/profile occurrences for focused profiles."""
    profile_columns = [
        "profile_id",
        "study_id",
        "experiment_id",
        "condition_group_followup",
        "condition",
        "direction_label",
        "body_site",
        "body_site_group",
        "title",
        "publication_year",
        "guild_z_mean",
        "oxygen_tolerance_z_mean",
    ]
    taxa = pd.read_csv(
        REPO_ROOT / "data/derived/bugsigdb_taxa.tsv",
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
        profile_id_from_parts(study_id, experiment_id, direction)
        for study_id, experiment_id, direction in taxa[
            ["study_id", "experiment_id", "direction"]
        ].itertuples(index=False, name=None)
    ]
    occurrences = taxa.merge(
        profiles[profile_columns],
        on="profile_id",
        how="inner",
        suffixes=("_taxa", ""),
    )
    occurrences["condition_group"] = occurrences["condition_group_followup"]
    return occurrences.drop(columns="condition_group_followup")


def load_guild_trait_annotations(taxon_ids: set[str]) -> pd.DataFrame:
    """Load aerobe v2 guild trait annotations for observed taxa."""
    matrix = pd.read_csv(
        REPO_ROOT / "data/derived/taxa_trait_matrix.tsv",
        sep="\t",
        usecols=lambda column: column in {"taxon_id", *GUILD_TRAITS.keys()},
    )
    matrix = matrix[matrix["taxon_id"].isin(taxon_ids)].copy()
    for trait_id, trait_name in GUILD_TRAITS.items():
        matrix[f"{trait_name}_count"] = pd.to_numeric(matrix[trait_id], errors="coerce").fillna(0)
        matrix[f"{trait_name}_present"] = matrix[f"{trait_name}_count"] > 0

    present_columns = [f"{trait_name}_present" for trait_name in GUILD_TRAITS.values()]
    matrix["guild_trait_presence_count"] = matrix[present_columns].sum(axis=1).astype(int)
    matrix["guild_trait_names"] = [
        "; ".join(
            trait_name for trait_name in GUILD_TRAITS.values() if bool(row[f"{trait_name}_present"])
        )
        for _, row in matrix.iterrows()
    ]
    matrix["is_focal_target_taxon"] = matrix["taxon_id"].isin(TARGET_TAXA)
    matrix["focal_target_taxon_name"] = matrix["taxon_id"].map(TARGET_TAXA).fillna("")
    return matrix[
        [
            "taxon_id",
            "guild_trait_presence_count",
            "guild_trait_names",
            "is_focal_target_taxon",
            "focal_target_taxon_name",
        ]
    ]


def add_taxa_annotations(summary: pd.DataFrame, taxon_ids: set[str]) -> pd.DataFrame:
    """Attach guild-trait annotations to a taxon summary table."""
    annotations = load_guild_trait_annotations(taxon_ids)
    output = summary.merge(annotations, on="taxon_id", how="left")
    output["guild_trait_presence_count"] = (
        pd.to_numeric(output["guild_trait_presence_count"], errors="coerce").fillna(0).astype(int)
    )
    for column in ["guild_trait_names", "focal_target_taxon_name"]:
        output[column] = output[column].fillna("")
    output["is_focal_target_taxon"] = output["is_focal_target_taxon"].fillna(False).astype(bool)
    return output


def summarize_taxa_by_direction(occurrences: pd.DataFrame) -> pd.DataFrame:
    """Summarize exact condition/direction taxon lists."""
    grouped = (
        occurrences.groupby(
            [
                "condition_group",
                "condition",
                "direction",
                "taxon_id",
                "taxon_name",
                "taxonomic_rank",
            ],
            dropna=False,
        )
        .agg(
            profile_count=("profile_id", "nunique"),
            study_count=("study_id", "nunique"),
            experiment_count=("experiment_id", "nunique"),
            signature_count=("signature_id", "nunique"),
            mean_profile_guild_z=("guild_z_mean", "mean"),
        )
        .reset_index()
    )
    output = add_taxa_annotations(grouped, set(occurrences["taxon_id"]))
    output["condition_order"] = output["condition"].map(
        {condition: index for index, condition in enumerate(CONDITION_ORDER)}
    )
    output["direction_order"] = output["direction"].map(DIRECTION_ORDER)
    return (
        output.sort_values(
            [
                "condition_order",
                "direction_order",
                "profile_count",
                "study_count",
                "guild_trait_presence_count",
                "taxon_name",
            ],
            ascending=[True, True, False, False, False, True],
        )
        .drop(columns=["condition_order", "direction_order"])
        .reset_index(drop=True)
    )


def summarize_taxa_by_group_direction(occurrences: pd.DataFrame) -> pd.DataFrame:
    """Summarize grouped condition/direction taxon lists."""
    grouped = (
        occurrences.groupby(
            [
                "condition_group",
                "direction",
                "taxon_id",
                "taxon_name",
                "taxonomic_rank",
            ],
            dropna=False,
        )
        .agg(
            profile_count=("profile_id", "nunique"),
            study_count=("study_id", "nunique"),
            experiment_count=("experiment_id", "nunique"),
            signature_count=("signature_id", "nunique"),
            exact_conditions=(
                "condition",
                lambda values: "; ".join(sorted(map(clean_text, set(values)))),
            ),
            body_site_groups=(
                "body_site_group",
                lambda values: "; ".join(sorted(map(clean_text, set(values)))),
            ),
            mean_profile_guild_z=("guild_z_mean", "mean"),
        )
        .reset_index()
    )
    output = add_taxa_annotations(grouped, set(occurrences["taxon_id"]))
    output["condition_group_order"] = output["condition_group"].map(
        {group: index for index, group in enumerate(CONDITION_GROUPS)}
    )
    output["direction_order"] = output["direction"].map(DIRECTION_ORDER)
    return (
        output.sort_values(
            [
                "condition_group_order",
                "direction_order",
                "profile_count",
                "study_count",
                "guild_trait_presence_count",
                "taxon_name",
            ],
            ascending=[True, True, False, False, False, True],
        )
        .drop(columns=["condition_group_order", "direction_order"])
        .reset_index(drop=True)
    )


def build_signed_taxa_balance(occurrences: pd.DataFrame) -> pd.DataFrame:
    """Build condition/taxon signed profile counts for increased minus decreased."""
    unique_taxon_profiles = occurrences.drop_duplicates(
        subset=["condition", "direction", "taxon_id", "profile_id"]
    ).copy()
    counts = (
        unique_taxon_profiles.groupby(
            [
                "condition_group",
                "condition",
                "taxon_id",
                "taxon_name",
                "taxonomic_rank",
                "direction",
            ],
            dropna=False,
        )["profile_id"]
        .nunique()
        .rename("profile_count")
        .reset_index()
    )
    pivot = counts.pivot_table(
        index=["condition_group", "condition", "taxon_id", "taxon_name", "taxonomic_rank"],
        columns="direction",
        values="profile_count",
        aggfunc="sum",
        fill_value=0,
    )
    for direction in ["increased", "decreased"]:
        if direction not in pivot.columns:
            pivot[direction] = 0
    balance = pivot.reset_index().rename(
        columns={
            "increased": "increased_profile_count",
            "decreased": "decreased_profile_count",
        }
    )
    balance["signed_profile_count"] = (
        balance["increased_profile_count"] - balance["decreased_profile_count"]
    )
    balance["total_profile_count"] = (
        balance["increased_profile_count"] + balance["decreased_profile_count"]
    )
    return add_taxa_annotations(balance, set(occurrences["taxon_id"]))


def compact_taxon_label(row: pd.Series, max_length: int = 34) -> str:
    """Create a compact, unique-enough heatmap label."""
    name = clean_text(row["taxon_name"]).replace("_", " ")
    if name.startswith("taxon "):
        name = row["taxon_id"].replace("NCBITaxon:", "NCBI:")
    if len(name) > max_length:
        name = f"{name[: max_length - 3]}..."
    return f"{name}\n{clean_text(row['taxonomic_rank'])}"


def combine_text_values(values: pd.Series) -> str:
    """Combine semicolon-delimited annotation strings into one stable value."""
    items: set[str] = set()
    for value in values:
        for item in clean_text(value).split(";"):
            item = item.strip()
            if item and item != "unknown":
                items.add(item)
    return "; ".join(sorted(items))


def build_taxon_label_frame(balance: pd.DataFrame) -> pd.DataFrame:
    """Build heatmap axis labels and export-friendly labels for taxa."""
    label_frame = balance[["taxon_id", "taxon_name", "taxonomic_rank"]].drop_duplicates("taxon_id")
    label_frame = label_frame.copy()
    label_frame["heatmap_axis_label"] = label_frame.apply(compact_taxon_label, axis=1)
    label_frame["heatmap_label"] = label_frame["heatmap_axis_label"].str.replace(
        "\n",
        " ",
        regex=False,
    )
    return label_frame


def select_heatmap_taxa(balance: pd.DataFrame, limit: int = HEATMAP_TOP_TAXA) -> pd.Index:
    """Select a readable top genus/species set for the heatmap."""
    candidates = balance[balance["taxonomic_rank"].isin(HEATMAP_RANKS)].copy()
    if candidates.empty:
        candidates = balance.copy()
    taxa_rank = (
        candidates.assign(abs_signed=candidates["signed_profile_count"].abs())
        .groupby(["taxon_id", "taxon_name", "taxonomic_rank"], dropna=False)
        .agg(
            max_abs_signed=("abs_signed", "max"),
            total_profile_count=("total_profile_count", "sum"),
            guild_trait_presence_count=("guild_trait_presence_count", "max"),
        )
        .reset_index()
        .sort_values(
            [
                "max_abs_signed",
                "total_profile_count",
                "guild_trait_presence_count",
                "taxon_name",
            ],
            ascending=[False, False, False, True],
        )
    )
    return pd.Index(taxa_rank.head(limit)["taxon_id"])


def build_heatmap_matrix(balance: pd.DataFrame, selected_taxa: pd.Index) -> pd.DataFrame:
    """Build the signed condition x taxon heatmap matrix."""
    plot_balance = balance[balance["taxon_id"].isin(selected_taxa)].copy()
    label_frame = build_taxon_label_frame(plot_balance)
    label_lookup = label_frame.set_index("taxon_id")["heatmap_axis_label"].to_dict()
    plot_balance["taxon_label"] = plot_balance["taxon_id"].map(label_lookup)
    matrix = plot_balance.pivot_table(
        index="condition",
        columns="taxon_label",
        values="signed_profile_count",
        aggfunc="sum",
        fill_value=0,
    )
    row_order = [condition for condition in CONDITION_ORDER if condition in matrix.index]
    matrix = matrix.loc[row_order]
    column_order = [
        label_lookup[taxon_id]
        for taxon_id in selected_taxa
        if taxon_id in label_lookup and label_lookup[taxon_id] in matrix.columns
    ]
    return matrix.loc[:, column_order]


def append_pattern_rows(
    rows: list[pd.DataFrame],
    source: pd.DataFrame,
    mask: pd.Series,
    pattern_name: str,
    pattern_color: str,
    profile_similarity_role: str,
    pattern_score: pd.Series,
    pattern_order: int,
) -> None:
    """Append one signed-pattern subset to the pattern-label output."""
    selected = source.loc[mask].copy()
    if selected.empty:
        return
    selected["pattern_name"] = pattern_name
    selected["pattern_color"] = pattern_color
    selected["profile_similarity_role"] = profile_similarity_role
    selected["pattern_score"] = pattern_score.loc[selected.index].round(3)
    selected["pattern_order"] = pattern_order
    rows.append(selected)


def build_signed_pattern_labels(
    balance: pd.DataFrame, heatmap_matrix: pd.DataFrame
) -> pd.DataFrame:
    """Extract taxa labels supporting red/blue Autism-Autoimmune profile patterns."""
    taxon_columns = ["taxon_id", "taxon_name", "taxonomic_rank"]
    signed = balance.pivot_table(
        index=taxon_columns,
        columns="condition",
        values="signed_profile_count",
        aggfunc="sum",
        fill_value=0,
    ).reset_index()
    for condition in CONDITION_ORDER:
        if condition not in signed.columns:
            signed[condition] = 0
    signed = signed.rename(
        columns={
            "Autism": "autism_signed_profile_count",
            "Autism spectrum disorder": "asd_signed_profile_count",
            "Autoimmune disease": "autoimmune_signed_profile_count",
        }
    )
    signed_count_columns = [
        "autism_signed_profile_count",
        "asd_signed_profile_count",
        "autoimmune_signed_profile_count",
    ]
    for column in signed_count_columns:
        signed[column] = pd.to_numeric(signed[column], errors="coerce").fillna(0).astype(int)

    annotations = (
        balance.groupby(taxon_columns, dropna=False)
        .agg(
            total_profile_count=("total_profile_count", "sum"),
            guild_trait_presence_count=("guild_trait_presence_count", "max"),
            guild_trait_names=("guild_trait_names", combine_text_values),
            is_focal_target_taxon=("is_focal_target_taxon", "max"),
            focal_target_taxon_name=("focal_target_taxon_name", combine_text_values),
        )
        .reset_index()
    )
    output = signed.merge(annotations, on=taxon_columns, how="left")
    output = output.merge(build_taxon_label_frame(balance), on=taxon_columns, how="left")
    output["in_heatmap"] = output["heatmap_axis_label"].isin(set(heatmap_matrix.columns))
    output["guild_trait_presence_count"] = (
        pd.to_numeric(output["guild_trait_presence_count"], errors="coerce").fillna(0).astype(int)
    )
    output["total_profile_count"] = (
        pd.to_numeric(output["total_profile_count"], errors="coerce").fillna(0).astype(int)
    )
    for column in ["guild_trait_names", "focal_target_taxon_name"]:
        output[column] = output[column].fillna("")
    output["is_focal_target_taxon"] = output["is_focal_target_taxon"].fillna(False).astype(bool)

    autism = output["autism_signed_profile_count"]
    asd = output["asd_signed_profile_count"]
    autoimmune = output["autoimmune_signed_profile_count"]
    guild_bonus = output["guild_trait_presence_count"] * 0.01
    rows: list[pd.DataFrame] = []
    append_pattern_rows(
        rows,
        output,
        (autism > 0) & (autoimmune > 0),
        "shared_red_exact_autism_autoimmune",
        "red",
        "Exact Autism and Autoimmune disease are both increased-skewed for this taxon.",
        autism.clip(lower=0) + autoimmune.clip(lower=0) + guild_bonus,
        0,
    )
    append_pattern_rows(
        rows,
        output,
        (autism < 0) & (autoimmune < 0),
        "shared_blue_exact_autism_autoimmune",
        "blue",
        "Exact Autism and Autoimmune disease are both decreased-skewed for this taxon.",
        -autism.clip(upper=0) - autoimmune.clip(upper=0) + guild_bonus,
        1,
    )
    append_pattern_rows(
        rows,
        output,
        (asd < 0) & ((autism >= 0) | (autoimmune >= 0)),
        "asd_blue_contrast_to_exact_labels",
        "blue",
        (
            "Broader Autism spectrum disorder is decreased-skewed while an exact "
            "comparison label is neutral or increased-skewed."
        ),
        -asd.clip(upper=0) + autism.clip(lower=0) * 0.1 + autoimmune.clip(lower=0) * 0.1,
        2,
    )
    if not rows:
        return pd.DataFrame(
            columns=[
                "pattern_name",
                "pattern_color",
                "profile_similarity_role",
                "taxon_id",
                "taxon_name",
                "taxonomic_rank",
                "heatmap_label",
                *signed_count_columns,
                "total_profile_count",
                "guild_trait_presence_count",
                "guild_trait_names",
                "is_focal_target_taxon",
                "focal_target_taxon_name",
                "in_heatmap",
                "pattern_score",
            ]
        )

    pattern_labels = pd.concat(rows, ignore_index=True)
    pattern_labels = pattern_labels.sort_values(
        [
            "pattern_order",
            "in_heatmap",
            "pattern_score",
            "total_profile_count",
            "guild_trait_presence_count",
            "taxon_name",
        ],
        ascending=[True, False, False, False, False, True],
    )
    return pattern_labels[
        [
            "pattern_name",
            "pattern_color",
            "profile_similarity_role",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
            "heatmap_label",
            *signed_count_columns,
            "total_profile_count",
            "guild_trait_presence_count",
            "guild_trait_names",
            "is_focal_target_taxon",
            "focal_target_taxon_name",
            "in_heatmap",
            "pattern_score",
        ]
    ].reset_index(drop=True)


def heatmap_pattern_style_lookup(pattern_labels: pd.DataFrame) -> dict[str, tuple[str, str]]:
    """Map export-friendly heatmap labels to visible plot highlight styles."""
    if pattern_labels.empty:
        return {}

    lookup: dict[str, tuple[str, str]] = {}
    plotted = pattern_labels[pattern_labels["in_heatmap"]].copy()
    for heatmap_label, group in plotted.groupby("heatmap_label", sort=False):
        pattern_names = set(group["pattern_name"])
        has_shared_red = "shared_red_exact_autism_autoimmune" in pattern_names
        has_blue = bool(
            pattern_names
            & {
                "shared_blue_exact_autism_autoimmune",
                "asd_blue_contrast_to_exact_labels",
            }
        )
        if has_shared_red and has_blue:
            lookup[heatmap_label] = (
                PATTERN_COLOR_BOTH,
                "shared red plus ASD blue contrast",
            )
        elif has_shared_red:
            lookup[heatmap_label] = (
                PATTERN_COLOR_EXACT_RED,
                "shared red: exact Autism and Autoimmune disease",
            )
        elif has_blue:
            lookup[heatmap_label] = (PATTERN_COLOR_BLUE, "blue contrast label")
    return lookup


def save_heatmap(
    matrix: pd.DataFrame,
    profile_summary: pd.DataFrame,
    pattern_labels: pd.DataFrame,
    output_path: Path,
) -> None:
    """Save a signed profile-count heatmap."""
    if matrix.empty:
        return

    row_counts = profile_summary.set_index("condition")["profile_count"].to_dict()
    display_matrix = matrix.copy()
    display_matrix.index = [
        f"{condition} (n={int(row_counts.get(condition, 0))})" for condition in matrix.index
    ]

    max_abs = float(display_matrix.abs().to_numpy().max())
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs) if max_abs else None
    fig_width = max(12.0, len(display_matrix.columns) * 0.34)
    fig_height = max(5.0, len(display_matrix.index) * 0.75 + 2.6)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(display_matrix.to_numpy(), aspect="auto", cmap="RdBu_r", norm=norm)
    pattern_styles = heatmap_pattern_style_lookup(pattern_labels)

    ax.set_xticks(range(len(display_matrix.columns)))
    ax.set_xticklabels(display_matrix.columns, rotation=90, ha="center", fontsize=7)
    for tick_label, column_label in zip(
        ax.get_xticklabels(),
        display_matrix.columns,
        strict=False,
    ):
        plot_label = column_label.replace("\n", " ")
        if plot_label not in pattern_styles:
            continue
        tick_color, _ = pattern_styles[plot_label]
        tick_label.set_color(tick_color)
        tick_label.set_fontweight("bold")
    ax.set_yticks(range(len(display_matrix.index)))
    ax.set_yticklabels(display_matrix.index, fontsize=9)
    ax.set_xlabel("Taxon")
    ax.set_ylabel("BugSigDB condition label")
    ax.set_title(
        "Autism and Autoimmune Disease Follow-Up: Aerobe v2 Taxa Direction Signal\n"
        "Signed unique taxon-profile count: increased profiles minus decreased profiles",
        pad=30,
    )
    ax.text(
        0,
        1.02,
        NON_GUILD_SPECIFIC_TAXA_HEATMAP_NOTE,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#475569",
    )
    for column_index, column_label in enumerate(display_matrix.columns):
        plot_label = column_label.replace("\n", " ")
        if plot_label not in pattern_styles:
            continue
        fill_color, _ = pattern_styles[plot_label]
        ax.add_patch(
            Rectangle(
                (column_index - 0.5, -0.76),
                1.0,
                0.18,
                facecolor=fill_color,
                edgecolor="none",
                clip_on=False,
            )
        )
    ax.text(
        -0.75,
        -0.67,
        "pattern",
        ha="right",
        va="center",
        fontsize=7.5,
        color="#475569",
        clip_on=False,
    )
    ax.set_ylim(len(display_matrix.index) - 0.5, -0.9)

    ax.set_xticks([tick - 0.5 for tick in range(1, len(display_matrix.columns))], minor=True)
    ax.set_yticks([tick - 0.5 for tick in range(1, len(display_matrix.index))], minor=True)
    ax.grid(which="minor", color="#e5e7eb", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    values = display_matrix.to_numpy()
    threshold = max_abs * 0.45 if max_abs else 0
    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            value = int(values[row_index, column_index])
            if value == 0:
                continue
            color = "white" if abs(value) >= threshold else "#111827"
            ax.text(
                column_index,
                row_index,
                str(value),
                ha="center",
                va="center",
                fontsize=6.5,
                color=color,
            )

    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.01)
    colorbar.set_label("increased minus decreased profile count")
    legend_handles = [
        Patch(
            facecolor=PATTERN_COLOR_EXACT_RED,
            label="shared red exact Autism + Autoimmune",
        ),
        Patch(facecolor=PATTERN_COLOR_BLUE, label="ASD blue contrast"),
        Patch(facecolor=PATTERN_COLOR_BOTH, label="both"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.18),
        ncol=3,
        frameon=False,
        fontsize=8,
        handlelength=1.1,
        columnspacing=1.1,
    )
    fig.text(
        0.01,
        0.015,
        PROFILE_SIMILARITY_PATTERN_NOTE,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#334155",
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def append_pattern_label_table(
    lines: list[str],
    pattern_labels: pd.DataFrame,
    pattern_name: str,
    title: str,
) -> None:
    """Append a compact table of red/blue heatmap labels to the Markdown report."""
    selected = pattern_labels[pattern_labels["pattern_name"] == pattern_name].head(
        PATTERN_LABEL_TOP_N
    )
    lines.extend(
        [
            f"### {title}",
            "",
            "| heatmap label | Autism | ASD | Autoimmune disease | selected traits | in heatmap |",
            "| --- | ---: | ---: | ---: | --- | --- |",
        ]
    )
    if selected.empty:
        lines.append("| none |  |  |  |  |  |")
        lines.append("")
        return

    for row in selected.itertuples(index=False):
        guild_traits = clean_text(row.guild_trait_names)
        if guild_traits == "unknown":
            guild_traits = ""
        in_heatmap = "yes" if bool(row.in_heatmap) else "no"
        lines.append(
            f"| {row.heatmap_label} | {int(row.autism_signed_profile_count)} | "
            f"{int(row.asd_signed_profile_count)} | "
            f"{int(row.autoimmune_signed_profile_count)} | {guild_traits} | "
            f"{in_heatmap} |"
        )
    lines.append("")


def trait_z_shift_phrase(label: str, row: pd.Series) -> str:
    """Describe one label's increased-minus-decreased trait z shift with its sign, size and n."""
    shift = float(row["increased_minus_decreased_guild_z"])
    counts = (
        f"{int(row['increased_profile_count'])} increased and "
        f"{int(row['decreased_profile_count'])} decreased profiles from "
        f"{int(row['study_count'])} studies"
    )
    if pd.isna(shift):
        # The shift is the increased mean minus the decreased mean, so it needs both.
        return (
            f"`{label}` has no computable increased-minus-decreased trait z shift, which "
            f"needs both increased and decreased profiles ({counts})"
        )
    sign = "positive" if shift > 0 else "negative" if shift < 0 else "zero"
    return (
        f"`{label}` has a {sign} increased-minus-decreased trait z shift ({shift:+.3f}; {counts})"
    )


def autism_label_contrast_sentence(profile_summary: pd.DataFrame) -> str:
    """Explain why exact Autism and ASD are kept apart, from their computed shifts (#128)."""
    labels = ("Autism", "Autism spectrum disorder")
    rows = profile_summary.set_index("condition")
    lead = (
        "Exact `Autism` and `Autism spectrum disorder` are reported separately rather than pooled"
    )
    if not all(label in rows.index for label in labels):
        return f"{lead} (#128)."
    autism, asd = (rows.loc[label] for label in labels)
    autism_shift = float(autism["increased_minus_decreased_guild_z"])
    asd_shift = float(asd["increased_minus_decreased_guild_z"])
    phrases = [trait_z_shift_phrase(label, rows.loc[label]) for label in labels]
    # An undefined (NaN) shift makes this comparison false, so it never claims divergence.
    if autism_shift * asd_shift < 0:
        return (
            f"{lead} because their directional skews diverge: {phrases[0]}, whereas "
            f"{phrases[1]} (#128)."
        )
    return f"{lead} (#128): {phrases[0]}, and {phrases[1]}."


def write_readme(
    profile_summary: pd.DataFrame,
    grouped_summary: pd.DataFrame,
    balance: pd.DataFrame,
    heatmap_matrix: pd.DataFrame,
    pattern_labels: pd.DataFrame,
    inferred_group_summary: pd.DataFrame,
    kg_microbe_taxon_universe: pd.DataFrame,
) -> None:
    """Write a compact human-readable follow-up summary."""
    lines = [
        "# Autism and Autoimmune Disease Follow-Up: Aerobe v2",
        "",
        "This follow-up starts from aerobe v2 trait enrichment in increased profiles "
        "and inspects the BugSigDB taxa reported in the autism/autism spectrum disorder "
        "and autoimmune disease condition labels.",
        "",
        "## Plot Guild Specificity",
        "",
        "| plot | guild-specific? | basis |",
        "| --- | --- | --- |",
        "| `autism_autoimmune_taxa_signed_heatmap.png` | no | "
        "taxon direction-count heatmap; selected traits are annotation context |",
        "",
        "## Included Condition Labels",
        "",
        "| condition group | exact condition label | profiles | studies | "
        "increased / decreased profiles | increased - decreased trait z | body-site groups |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in profile_summary.itertuples(index=False):
        lines.append(
            f"| {row.condition_group} | {row.condition} | {int(row.profile_count)} | "
            f"{int(row.study_count)} | {int(row.increased_profile_count)} / "
            f"{int(row.decreased_profile_count)} | "
            f"{row.increased_minus_decreased_guild_z:.3f} | {row.body_site_groups} |"
        )

    lines.extend(
        [
            "",
            "## Top Taxa",
            "",
            "The tables below count unique taxon/profile occurrences after deduplicating "
            "the same taxon within the same study/experiment/direction profile. "
            + autism_label_contrast_sentence(profile_summary),
            "",
        ]
    )
    for condition_group in CONDITION_GROUPS:
        lines.extend([f"### {condition_group}", ""])
        for direction in ["increased", "decreased"]:
            selected = grouped_summary[
                (grouped_summary["condition_group"] == condition_group)
                & (grouped_summary["direction"] == direction)
            ].head(12)
            lines.extend(
                [
                    f"Top {direction} taxa:",
                    "",
                    "| taxon | rank | profiles | studies | selected traits |",
                    "| --- | --- | ---: | ---: | --- |",
                ]
            )
            if selected.empty:
                lines.append("| none |  | 0 | 0 |  |")
            for row in selected.itertuples(index=False):
                guild_traits = clean_text(row.guild_trait_names)
                if guild_traits == "unknown":
                    guild_traits = ""
                lines.append(
                    f"| {row.taxon_name} | {row.taxonomic_rank} | "
                    f"{int(row.profile_count)} | {int(row.study_count)} | {guild_traits} |"
                )
            lines.append("")

    lines.extend(
        [
            "## Red/Blue Pattern Labels",
            "",
            "These lists make the profile-similarity interpretation auditable. "
            "Shared red labels have positive signed counts in both exact `Autism` "
            "and `Autoimmune disease`, supporting the increased-profile similarity "
            "readout. Shared blue labels would have negative signed counts in both "
            "exact labels. The ASD blue-contrast list marks taxa where the broader "
            "`Autism spectrum disorder` label is decreased-skewed while at least one "
            "exact comparison label is neutral or increased-skewed.",
            "",
        ]
    )
    append_pattern_label_table(
        lines,
        pattern_labels,
        "shared_red_exact_autism_autoimmune",
        "Shared red labels: exact Autism and Autoimmune disease",
    )
    append_pattern_label_table(
        lines,
        pattern_labels,
        "shared_blue_exact_autism_autoimmune",
        "Shared blue labels: exact Autism and Autoimmune disease",
    )
    append_pattern_label_table(
        lines,
        pattern_labels,
        "asd_blue_contrast_to_exact_labels",
        "ASD blue-contrast labels",
    )

    if not inferred_group_summary.empty:
        lines.extend(
            [
                "## Inferred Guild Evidence",
                "",
                "This section uses a data-derived trait-enrichment view: for each exact condition "
                "label, it identifies v2 traits enriched in that condition's increased "
                "profiles, then lists taxa observed in those same profiles that also "
                "belong to the KG-Microbe-wide candidate universe and carry the enriched "
                "trait annotations. The universe contains "
                f"{len(kg_microbe_taxon_universe):,} taxa with at least one core "
                "oxygen/respiration criterion annotation, not the fixed focal taxa list.",
                "",
                "| condition | increased profiles | universe taxa / observed taxa | "
                "enriched traits | supporting taxa |",
                "| --- | ---: | ---: | --- | --- |",
            ]
        )
        for row in inferred_group_summary.itertuples(index=False):
            lines.append(
                f"| {row.condition} | {int(row.focus_profile_count)} | "
                f"{int(row.kg_microbe_universe_taxon_count)} / "
                f"{int(row.total_observed_taxon_count)} | "
                f"{clean_text(row.top_enriched_traits)} | "
                f"{clean_text(row.top_supporting_taxa)} |"
            )
        lines.append("")

    heatmap_taxa = len(heatmap_matrix.columns)
    full_taxa = balance["taxon_id"].nunique()
    guild_trait_text = ", ".join(GUILD_TRAITS.values())
    positive_conditions = profile_summary[profile_summary["increased_minus_decreased_guild_z"] > 0]
    if positive_conditions.empty:
        positive_condition_text = "No included condition label has a positive v2 guild shift."
    else:
        condition_bits = [
            f"`{row.condition}` ({row.increased_minus_decreased_guild_z:.3f})"
            for row in positive_conditions.itertuples(index=False)
        ]
        positive_condition_text = (
            "Included condition labels with positive increased-minus-decreased v2 trait-enrichment "
            "shifts are " + ", ".join(condition_bits) + "."
        )
    lines.extend(
        [
            "## Interpretation",
            "",
            "This v2 follow-up uses only the narrowed aerobe v2 selected traits: "
            f"{guild_trait_text}. It supports a label-dependent signal check rather than "
            "a broad causal claim.",
            "",
            positive_condition_text,
            "",
            "At the taxa-list level, overlap should still be interpreted cautiously. "
            "The signed heatmap counts direction balance by taxon/profile occurrence, "
            "while the guild annotation columns only indicate whether those taxa carry "
            "one or more of the v2 oxygen-tolerance traits.",
            "",
            "`autism_autoimmune_taxa_signed_pattern_labels.tsv` extracts the heatmap "
            "labels that directly support the red/blue readout above so the visual "
            "pattern can be audited without manually reading the plot.",
            "",
            "The most defensible interpretation is that the exact autism and "
            "autoimmune disease labels both show taxa annotated as facultative or "
            "aerotolerant enriched among increased profiles. This "
            "analysis does not show that autism and autoimmune disease are causally "
            "linked, driven by the same taxa, or even consistently connected across "
            "broader ASD metadata. The autoimmune subset is small, so the taxa-level "
            "results are candidates for follow-up.",
            "",
            "## Heatmap",
            "",
            f"`autism_autoimmune_taxa_signed_heatmap.png` shows the top {heatmap_taxa} "
            f"genus/species taxa selected from {full_taxa} observed taxa. Cell values are "
            "`increased unique taxon-profile count - decreased unique taxon-profile count`; "
            "red is increased-skewed and blue is decreased-skewed. Colored x-axis labels "
            "and the marker strip above the heatmap identify taxa included in the "
            "red/blue pattern-label extraction.",
            "",
            "## Files",
            "",
            "- `autism_autoimmune_profile_summary.tsv`",
            "- `autism_autoimmune_taxa_direction_summary.tsv`",
            "- `autism_autoimmune_taxa_direction_summary_grouped.tsv`",
            "- `autism_autoimmune_taxa_signed_balance.tsv`",
            "- `autism_autoimmune_taxa_signed_heatmap_matrix.tsv`",
            "- `autism_autoimmune_taxa_signed_pattern_labels.tsv`",
            "- `autism_autoimmune_taxa_signed_heatmap.png`",
            "- `autism_autoimmune_inferred_guild_summary.tsv`",
            "- `autism_autoimmune_inferred_guild_trait_evidence.tsv`",
            "- `autism_autoimmune_inferred_guild_taxa_evidence.tsv`",
            "- `kg_microbe_candidate_taxon_universe.tsv`",
            "",
        ]
    )
    (OUTPUT_DIR / "autism_autoimmune_followup.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    study_matrix, _, study_trait_metadata = load_study_heatmap()
    profiles = load_scored_condition_profiles()
    profiles["condition_group"] = profiles["condition_group_followup"]
    if profiles.empty:
        msg = "No matching profiles found for autism or autoimmune disease condition labels."
        raise ValueError(msg)

    profile_summary = summarize_profiles(profiles)
    occurrences = load_taxa_occurrences(profiles)
    taxa_direction_summary = summarize_taxa_by_direction(occurrences)
    grouped_taxa_summary = summarize_taxa_by_group_direction(occurrences)
    signed_balance = build_signed_taxa_balance(occurrences)
    selected_taxa = select_heatmap_taxa(signed_balance)
    heatmap_matrix = build_heatmap_matrix(signed_balance, selected_taxa)
    pattern_labels = build_signed_pattern_labels(signed_balance, heatmap_matrix)
    kg_microbe_taxon_universe = build_kg_microbe_taxon_universe()
    inferred_group_summary, inferred_trait_evidence, inferred_taxa_evidence = (
        summarize_inferred_guild_evidence(
            study_matrix,
            profiles,
            study_trait_metadata,
            occurrences,
            profile_summary[["condition_group", "condition"]],
            ["condition_group", "condition"],
            focus_direction="increased",
            taxon_universe_ids=set(kg_microbe_taxon_universe["taxon_id"]),
        )
    )

    profile_summary.to_csv(
        OUTPUT_DIR / "autism_autoimmune_profile_summary.tsv",
        sep="\t",
        index=False,
    )
    taxa_direction_summary.to_csv(
        OUTPUT_DIR / "autism_autoimmune_taxa_direction_summary.tsv",
        sep="\t",
        index=False,
    )
    grouped_taxa_summary.to_csv(
        OUTPUT_DIR / "autism_autoimmune_taxa_direction_summary_grouped.tsv",
        sep="\t",
        index=False,
    )
    signed_balance.to_csv(
        OUTPUT_DIR / "autism_autoimmune_taxa_signed_balance.tsv",
        sep="\t",
        index=False,
    )
    heatmap_matrix.to_csv(
        OUTPUT_DIR / "autism_autoimmune_taxa_signed_heatmap_matrix.tsv",
        sep="\t",
    )
    pattern_labels.to_csv(
        OUTPUT_DIR / "autism_autoimmune_taxa_signed_pattern_labels.tsv",
        sep="\t",
        index=False,
    )
    kg_microbe_taxon_universe.to_csv(
        OUTPUT_DIR / "kg_microbe_candidate_taxon_universe.tsv",
        sep="\t",
        index=False,
    )
    inferred_group_summary.to_csv(
        OUTPUT_DIR / "autism_autoimmune_inferred_guild_summary.tsv",
        sep="\t",
        index=False,
    )
    inferred_trait_evidence.to_csv(
        OUTPUT_DIR / "autism_autoimmune_inferred_guild_trait_evidence.tsv",
        sep="\t",
        index=False,
    )
    inferred_taxa_evidence.to_csv(
        OUTPUT_DIR / "autism_autoimmune_inferred_guild_taxa_evidence.tsv",
        sep="\t",
        index=False,
    )
    save_heatmap(
        heatmap_matrix,
        profile_summary,
        pattern_labels,
        OUTPUT_DIR / "autism_autoimmune_taxa_signed_heatmap.png",
    )
    write_readme(
        profile_summary,
        grouped_taxa_summary,
        signed_balance,
        heatmap_matrix,
        pattern_labels,
        inferred_group_summary,
        kg_microbe_taxon_universe,
    )

    print(f"Wrote follow-up outputs to: {OUTPUT_DIR}")
    print(f"Profiles analyzed: {len(profiles)}")
    print(f"Taxon/profile occurrence rows: {len(occurrences)}")
    print(f"Unique taxa observed: {occurrences['taxon_id'].nunique()}")
    print(f"KG-Microbe candidate taxon universe rows: {len(kg_microbe_taxon_universe)}")
    print(f"Inferred guild rows: {len(inferred_group_summary)}")
    print(f"Heatmap taxa plotted: {len(heatmap_matrix.columns)}")
    print(f"Pattern label rows: {len(pattern_labels)}")


if __name__ == "__main__":
    main()
