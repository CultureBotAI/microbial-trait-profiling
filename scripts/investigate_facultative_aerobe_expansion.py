"""Investigate facultative aerobe and inflammatory expansion signals."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
OUTPUT_DIR = VIS_DIR / "trait-profile-investigation" / "facultative-aerobe-inflammatory-expansion"

GUILD_TRAITS = {
    "METPO:1000608": "facultatively aerobic",
    "METPO:1000602": "aerobic",
    "METPO:1000699": "gram negative",
    "METPO:1000702": "motile",
    "METPO:1000604": "microaerophilic",
    "METPO:2000103|pathways:nitrate_reduction": "nitrate_reduction",
}

TARGET_TAXA = {
    "NCBITaxon:33033": "Parvimonas micra",
    "NCBITaxon:573": "Klebsiella pneumoniae",
    "NCBITaxon:1328": "Streptococcus anginosus",
    "NCBITaxon:562": "Escherichia coli",
    "NCBITaxon:851": "Fusobacterium nucleatum",
    "NCBITaxon:208479": "Enterocloster bolteae",
    "NCBITaxon:84112": "Eggerthella lenta",
    "NCBITaxon:28131": "Prevotella intermedia",
    "NCBITaxon:837": "Porphyromonas gingivalis",
    "NCBITaxon:39778": "Veillonella dispar",
}


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def clean_text(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return math.nan
    return numerator / denominator


def profile_id_from_taxa_row(row: pd.Series) -> str:
    direction = "inc" if row["direction"] == "increased" else "dec"
    return f"{row['study_id']}|{row['experiment_id']}|{direction}"


def load_study_heatmap() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix = read_tsv(VIS_DIR / "study_trait_profile_clustered_heatmap_matrix.tsv", index_col=0)
    metadata = read_tsv(VIS_DIR / "study_trait_profile_study_metadata.tsv", index_col=0)
    trait_metadata = read_tsv(VIS_DIR / "study_trait_profile_trait_metadata.tsv", index_col=0)

    matrix.index = metadata.index
    column_labels = {
        row["display_label"]: row["object_name"] for _, row in trait_metadata.iterrows()
    }
    matrix = matrix.rename(columns=column_labels)
    return matrix, metadata, trait_metadata


def add_guild_scores(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    available_traits = [trait for trait in GUILD_TRAITS.values() if trait in matrix.columns]
    scored = metadata.copy()
    scored["guild_z_mean"] = matrix[available_traits].mean(axis=1)
    scored["oxygen_tolerance_z_mean"] = matrix[
        ["facultatively aerobic", "aerobic", "microaerophilic", "nitrate_reduction"]
    ].mean(axis=1)
    scored["gram_negative_motility_z_mean"] = matrix[["gram negative", "motile"]].mean(axis=1)
    for trait_name in available_traits:
        scored[f"z_{trait_name}"] = matrix[trait_name]
    return scored


def summarize_direction(metadata: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for trait_name in GUILD_TRAITS.values():
        column = f"z_{trait_name}"
        inc_mean = metadata.loc[metadata["inc_dec"] == "inc", column].mean()
        dec_mean = metadata.loc[metadata["inc_dec"] == "dec", column].mean()
        rows.append(
            {
                "trait": trait_name,
                "increased_profile_mean_z": inc_mean,
                "decreased_profile_mean_z": dec_mean,
                "increased_minus_decreased_z": inc_mean - dec_mean,
            }
        )
    for score_name in [
        "guild_z_mean",
        "oxygen_tolerance_z_mean",
        "gram_negative_motility_z_mean",
    ]:
        inc_mean = metadata.loc[metadata["inc_dec"] == "inc", score_name].mean()
        dec_mean = metadata.loc[metadata["inc_dec"] == "dec", score_name].mean()
        rows.append(
            {
                "trait": score_name,
                "increased_profile_mean_z": inc_mean,
                "decreased_profile_mean_z": dec_mean,
                "increased_minus_decreased_z": inc_mean - dec_mean,
            }
        )
    return pd.DataFrame(rows)


def summarize_group(metadata: pd.DataFrame, group_column: str, min_profiles: int) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for group_value, group in metadata.groupby(group_column, dropna=False):
        if len(group) < min_profiles:
            continue
        inc = group[group["inc_dec"] == "inc"]
        dec = group[group["inc_dec"] == "dec"]
        records.append(
            {
                group_column: clean_text(group_value),
                "profile_count": len(group),
                "study_count": group["study_id"].nunique(),
                "increased_profile_count": len(inc),
                "decreased_profile_count": len(dec),
                "mean_guild_z": group["guild_z_mean"].mean(),
                "increased_mean_guild_z": inc["guild_z_mean"].mean(),
                "decreased_mean_guild_z": dec["guild_z_mean"].mean(),
                "increased_minus_decreased_guild_z": inc["guild_z_mean"].mean()
                - dec["guild_z_mean"].mean(),
                "increased_mean_oxygen_tolerance_z": inc["oxygen_tolerance_z_mean"].mean(),
                "increased_mean_gram_negative_motility_z": inc[
                    "gram_negative_motility_z_mean"
                ].mean(),
            }
        )
    return pd.DataFrame(records).sort_values(
        ["increased_mean_guild_z", "profile_count"],
        ascending=[False, False],
    )


def load_taxa_direction_summary() -> pd.DataFrame:
    taxa_metadata = read_tsv(VIS_DIR / "taxa_guild_discovery_taxa_metadata.tsv", index_col=0)
    summary = taxa_metadata.loc[taxa_metadata.index.intersection(TARGET_TAXA.keys())].copy()
    summary["target_taxon_name"] = summary.index.map(TARGET_TAXA)
    summary["increased_fraction"] = [
        safe_ratio(float(row["increased"]), float(row["bugsigdb_profile_count"]))
        for _, row in summary.iterrows()
    ]
    summary["expansion_bias"] = [
        safe_ratio(
            float(row["increased"]) - float(row["decreased"]),
            float(row["bugsigdb_profile_count"]),
        )
        for _, row in summary.iterrows()
    ]
    summary["increased_to_decreased_ratio"] = [
        safe_ratio(float(row["increased"]), float(row["decreased"]))
        for _, row in summary.iterrows()
    ]
    return summary[
        [
            "target_taxon_name",
            "taxon_name",
            "taxonomic_rank",
            "bugsigdb_study_count",
            "bugsigdb_profile_count",
            "decreased",
            "increased",
            "dominant_direction",
            "increased_fraction",
            "expansion_bias",
            "increased_to_decreased_ratio",
            "trait_count",
        ]
    ].sort_values("expansion_bias", ascending=False)


def load_target_taxa_trait_counts() -> pd.DataFrame:
    usecols = ["taxon_id", "taxon_name", "taxonomic_rank", *GUILD_TRAITS.keys()]
    matrix = read_tsv(REPO_ROOT / "data/derived/taxa_trait_matrix.tsv", usecols=usecols)
    matrix = matrix[matrix["taxon_id"].isin(TARGET_TAXA)].copy()
    taxa_metadata = read_tsv(VIS_DIR / "taxa_guild_discovery_taxa_metadata.tsv", index_col=0)
    matrix["target_taxon_name"] = matrix["taxon_id"].map(TARGET_TAXA)
    matrix["heatmap_taxon_name"] = matrix["taxon_id"].map(taxa_metadata["taxon_name"])
    matrix["heatmap_taxonomic_rank"] = matrix["taxon_id"].map(taxa_metadata["taxonomic_rank"])

    for trait_id, trait_name in GUILD_TRAITS.items():
        matrix[f"{trait_name}_count"] = pd.to_numeric(matrix[trait_id], errors="coerce").fillna(0)
        matrix[f"{trait_name}_present"] = matrix[f"{trait_name}_count"] > 0

    present_columns = [f"{trait_name}_present" for trait_name in GUILD_TRAITS.values()]
    matrix["guild_trait_presence_count"] = matrix[present_columns].sum(axis=1)
    output_columns = [
        "taxon_id",
        "target_taxon_name",
        "heatmap_taxonomic_rank",
        "guild_trait_presence_count",
    ]
    for trait_name in GUILD_TRAITS.values():
        output_columns.extend([f"{trait_name}_present", f"{trait_name}_count"])
    output = matrix[output_columns].rename(columns={"heatmap_taxonomic_rank": "taxonomic_rank"})
    return output.sort_values(
        ["guild_trait_presence_count", "target_taxon_name"],
        ascending=[False, True],
    )


def load_target_taxa_occurrences(scored_profiles: pd.DataFrame) -> pd.DataFrame:
    taxa = read_tsv(
        REPO_ROOT / "data/derived/taxa.tsv",
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
    taxa = taxa[taxa["taxon_id"].isin(TARGET_TAXA)].copy()
    taxa["profile_id"] = taxa.apply(profile_id_from_taxa_row, axis=1)
    taxa = taxa.drop_duplicates(subset=["taxon_id", "profile_id"])
    taxa["target_taxon_name"] = taxa["taxon_id"].map(TARGET_TAXA)
    profile_columns = [
        "study_id",
        "experiment_id",
        "inc_dec",
        "direction_label",
        "title",
        "publication_year",
        "condition",
        "body_site",
        "body_site_group",
        "condition_group",
        "location",
        "host_species",
        "guild_z_mean",
        "oxygen_tolerance_z_mean",
        "gram_negative_motility_z_mean",
    ]
    profiles = scored_profiles[profile_columns].copy()
    return taxa.merge(
        profiles,
        left_on="profile_id",
        right_index=True,
        how="left",
        suffixes=("_taxa", ""),
    )


def summarize_taxa_conditions(occurrences: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        occurrences.groupby(["target_taxon_name", "condition", "body_site_group", "direction"])[
            ["profile_id", "study_id"]
        ]
        .nunique()
        .rename(columns={"profile_id": "profile_count", "study_id": "study_count"})
        .reset_index()
    )
    if grouped.empty:
        return grouped

    profile_pivot = grouped.pivot_table(
        index=["target_taxon_name", "condition", "body_site_group"],
        columns="direction",
        values="profile_count",
        aggfunc="sum",
        fill_value=0,
    )
    study_pivot = grouped.pivot_table(
        index=["target_taxon_name", "condition", "body_site_group"],
        columns="direction",
        values="study_count",
        aggfunc="sum",
        fill_value=0,
    )
    for direction in ["decreased", "increased"]:
        if direction not in profile_pivot.columns:
            profile_pivot[direction] = 0
        if direction not in study_pivot.columns:
            study_pivot[direction] = 0

    output = profile_pivot.rename(
        columns={
            "decreased": "decreased_profile_count",
            "increased": "increased_profile_count",
        }
    ).join(
        study_pivot.rename(
            columns={
                "decreased": "decreased_study_count",
                "increased": "increased_study_count",
            }
        )
    )
    output = output.reset_index()
    output["total_profile_count"] = (
        output["decreased_profile_count"] + output["increased_profile_count"]
    )
    output["expansion_bias"] = (
        output["increased_profile_count"] - output["decreased_profile_count"]
    ) / output["total_profile_count"].replace(0, pd.NA)
    return output.sort_values(
        ["expansion_bias", "increased_study_count", "total_profile_count"],
        ascending=[False, False, False],
    )


def make_plots(
    taxa_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    direction_summary: pd.DataFrame,
) -> None:
    plt.rcParams.update({"figure.dpi": 160, "savefig.dpi": 200})

    fig, ax = plt.subplots(figsize=(9, 5.4))
    plot_taxa = taxa_summary.sort_values("expansion_bias")
    ax.barh(plot_taxa["target_taxon_name"], plot_taxa["expansion_bias"], color="#b91c1c")
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("(increased profiles - decreased profiles) / total profiles")
    ax.set_title("Facultative Aerobe Candidate Taxa Are Expansion-Biased")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "candidate_taxa_expansion_bias.png")
    plt.close(fig)

    top_conditions = (
        condition_summary[
            (condition_summary["increased_profile_count"] >= 3)
            & (condition_summary["study_count"] >= 2)
            & (condition_summary["increased_minus_decreased_guild_z"] > 0)
            & (condition_summary["increased_mean_guild_z"] > 0)
        ]
        .sort_values(
            ["increased_minus_decreased_guild_z", "increased_mean_guild_z"],
            ascending=[False, False],
        )
        .head(12)
    )
    if top_conditions.empty:
        top_conditions = condition_summary[condition_summary["increased_profile_count"] >= 3].head(
            12
        )
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.barh(
        top_conditions["condition"][::-1],
        top_conditions["increased_mean_guild_z"][::-1],
        color="#dc2626",
    )
    ax.set_xlabel("Mean facultative-aerobe guild z-score in increased profiles")
    ax.set_title("Conditions With High Increased-Profile Guild Signal")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "condition_increased_guild_signal.png")
    plt.close(fig)

    plot_direction = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("increased_minus_decreased_z")
    fig, ax = plt.subplots(figsize=(8, 4.4))
    ax.barh(
        plot_direction["trait"],
        plot_direction["increased_minus_decreased_z"],
        color="#be123c",
    )
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("increased profile mean z - decreased profile mean z")
    ax.set_title("Guild Traits Are Higher in Increased Profiles")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "guild_trait_direction_shift.png")
    plt.close(fig)


def write_summary_report(
    direction_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    taxa_summary: pd.DataFrame,
    condition_taxa_summary: pd.DataFrame,
    top_expansion_profiles: pd.DataFrame,
) -> None:
    direction_rows = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("increased_minus_decreased_z", ascending=False)
    top_conditions = (
        condition_summary[
            (condition_summary["increased_profile_count"] >= 3)
            & (condition_summary["study_count"] >= 2)
            & (condition_summary["increased_minus_decreased_guild_z"] > 0)
            & (condition_summary["increased_mean_guild_z"] > 0)
        ]
        .sort_values(
            ["increased_minus_decreased_guild_z", "increased_mean_guild_z"],
            ascending=[False, False],
        )
        .head(8)
    )
    if top_conditions.empty:
        top_conditions = condition_summary[condition_summary["increased_profile_count"] >= 3].head(
            8
        )
    top_taxa = taxa_summary.head(10)
    top_taxon_conditions = condition_taxa_summary[
        (condition_taxa_summary["increased_profile_count"] >= 3)
        & (condition_taxa_summary["increased_study_count"] >= 2)
    ].head(12)

    lines = [
        "# Facultative Aerobe / Inflammatory Expansion Investigation",
        "",
        "This investigation uses the existing study trait profile clustered heatmap matrix, "
        "study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.",
        "",
        "## Scope",
        "",
        "- Guild traits: "
        + ", ".join(f"{name} ({trait_id})" for trait_id, name in GUILD_TRAITS.items()),
        "- Target taxa: " + ", ".join(TARGET_TAXA.values()),
        "- Context: all body sites, with body-site and condition summaries to separate "
        "gut, oral, airway, and other patterns.",
        "",
        "## Main Signal",
        "",
        "The selected guild traits are higher in increased profiles than in decreased "
        "profiles, supporting a facultative-aerobe/pathobiont expansion signal.",
        "",
        "| trait | increased - decreased z |",
        "| --- | ---: |",
    ]
    for row in direction_rows.itertuples(index=False):
        lines.append(f"| {row.trait} | {row.increased_minus_decreased_z:.3f} |")

    lines.extend(
        [
            "",
            "## Strongest Expansion-Biased Taxa",
            "",
            "| taxon | increased / decreased profiles | expansion bias |",
            "| --- | ---: | ---: |",
        ]
    )
    for row in top_taxa.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {int(row.increased)} / {int(row.decreased)} | "
            f"{row.expansion_bias:.3f} |"
        )

    lines.extend(
        [
            "",
            "## Recurrent Target-Taxon Expansions By Condition",
            "",
            "| taxon | condition | body site | "
            "increased / decreased profiles | increased studies |",
            "| --- | --- | --- | ---: | ---: |",
        ]
    )
    for row in top_taxon_conditions.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {clean_text(row.condition)} | "
            f"{clean_text(row.body_site_group)} | "
            f"{int(row.increased_profile_count)} / {int(row.decreased_profile_count)} | "
            f"{int(row.increased_study_count)} |"
        )

    lines.extend(
        [
            "",
            "## Contexts With High Increased-Profile Guild Signal",
            "",
            "| condition | profiles | studies | increased mean guild z | increased - decreased z |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in top_conditions.itertuples(index=False):
        diff = row.increased_minus_decreased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        lines.append(
            f"| {row.condition} | {int(row.profile_count)} | {int(row.study_count)} | "
            f"{row.increased_mean_guild_z:.3f} | {diff_text} |"
        )

    lines.extend(
        [
            "",
            "## Body-Site Check",
            "",
            "| body_site_group | profiles | increased mean guild z | increased - decreased z |",
            "| --- | ---: | ---: | ---: |",
        ]
    )
    for row in body_site_summary.head(10).itertuples(index=False):
        diff = row.increased_minus_decreased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        lines.append(
            f"| {row.body_site_group} | {int(row.profile_count)} | "
            f"{row.increased_mean_guild_z:.3f} | {diff_text} |"
        )

    lines.extend(
        [
            "",
            "## Top Increased Profiles To Inspect",
            "",
            "| profile | condition | body site | guild z | title |",
            "| --- | --- | --- | ---: | --- |",
        ]
    )
    for row in top_expansion_profiles.head(10).itertuples(index=False):
        title = clean_text(row.title).replace("|", "/")
        if len(title) > 90:
            title = f"{title[:87]}..."
        profile_id = str(row.profile_id).replace("|", "\\|")
        lines.append(
            f"| {profile_id} | {clean_text(row.condition)} | {clean_text(row.body_site)} | "
            f"{row.guild_z_mean:.3f} | {title} |"
        )

    lines.extend(
        [
            "",
            "## Suggested Follow-Up",
            "",
            "1. Split the signal by body site before treating it as inflammatory expansion; "
            "oral and airway profiles can dominate aerobic/microaerophilic traits.",
            "2. In gut-only analyses, test the target taxa against inflammation markers, "
            "nitrate respiration, oxygen tolerance, virulence genes, and resistome features.",
            "3. Compare expansion-biased taxa against the anaerobe-depletion guild to see "
            "whether the same conditions show reciprocal loss and expansion patterns.",
            "",
            "## Files",
            "",
            "- `guild_trait_direction_shift.tsv`",
            "- `body_site_guild_summary.tsv`",
            "- `condition_guild_summary.tsv`",
            "- `candidate_taxa_direction_summary.tsv`",
            "- `candidate_taxa_trait_presence.tsv`",
            "- `candidate_taxa_condition_direction.tsv`",
            "- `candidate_taxa_profile_occurrences.tsv`",
            "- `top_increased_guild_profiles.tsv`",
            "- `candidate_taxa_expansion_bias.png`",
            "- `condition_increased_guild_signal.png`",
            "- `guild_trait_direction_shift.png`",
            "",
        ]
    )
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    study_matrix, study_metadata, study_trait_metadata = load_study_heatmap()
    scored_profiles = add_guild_scores(study_matrix, study_metadata)

    direction_summary = summarize_direction(scored_profiles)
    body_site_summary = summarize_group(scored_profiles, "body_site_group", min_profiles=10)
    condition_summary = summarize_group(scored_profiles, "condition", min_profiles=5)
    taxa_summary = load_taxa_direction_summary()
    trait_presence = load_target_taxa_trait_counts()
    occurrences = load_target_taxa_occurrences(scored_profiles)
    condition_taxa_summary = summarize_taxa_conditions(occurrences)

    trait_definitions = study_trait_metadata.loc[
        study_trait_metadata.index.intersection(GUILD_TRAITS.keys())
    ].copy()
    trait_definitions["investigation_trait_name"] = trait_definitions.index.map(GUILD_TRAITS)

    top_expansion_profiles = scored_profiles[scored_profiles["inc_dec"] == "inc"].copy()
    top_expansion_profiles = top_expansion_profiles.sort_values(
        ["guild_z_mean", "total_trait_signal"],
        ascending=[False, False],
    ).head(100)
    top_expansion_profiles = top_expansion_profiles.reset_index(names="profile_id")

    direction_summary.to_csv(OUTPUT_DIR / "guild_trait_direction_shift.tsv", sep="\t", index=False)
    body_site_summary.to_csv(OUTPUT_DIR / "body_site_guild_summary.tsv", sep="\t", index=False)
    condition_summary.to_csv(OUTPUT_DIR / "condition_guild_summary.tsv", sep="\t", index=False)
    taxa_summary.to_csv(OUTPUT_DIR / "candidate_taxa_direction_summary.tsv", sep="\t")
    trait_presence.to_csv(OUTPUT_DIR / "candidate_taxa_trait_presence.tsv", sep="\t", index=False)
    condition_taxa_summary.to_csv(
        OUTPUT_DIR / "candidate_taxa_condition_direction.tsv", sep="\t", index=False
    )
    occurrences.sort_values(["target_taxon_name", "profile_id"]).to_csv(
        OUTPUT_DIR / "candidate_taxa_profile_occurrences.tsv",
        sep="\t",
        index=False,
    )
    top_expansion_profiles.to_csv(
        OUTPUT_DIR / "top_increased_guild_profiles.tsv",
        sep="\t",
        index=False,
    )
    trait_definitions.to_csv(OUTPUT_DIR / "guild_trait_definitions.tsv", sep="\t")

    make_plots(taxa_summary, condition_summary, direction_summary)
    write_summary_report(
        direction_summary,
        body_site_summary,
        condition_summary,
        taxa_summary,
        condition_taxa_summary,
        top_expansion_profiles,
    )

    print(f"Wrote investigation outputs to: {OUTPUT_DIR}")
    print(f"Study profiles analyzed: {len(scored_profiles)}")
    print(f"Target taxa occurrence rows: {len(occurrences)}")
    print(f"Condition/taxon rows: {len(condition_taxa_summary)}")


if __name__ == "__main__":
    main()
