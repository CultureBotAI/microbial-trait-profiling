"""Investigate gut anaerobe and fermentation depletion signals."""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
OUTPUT_DIR = VIS_DIR / "trait-profile-investigation" / "gut-anaerobe-fermentation-depletion"

GUILD_TRAITS = {
    "METPO:1000603": "anaerobic",
    "METPO:1000607": "obligately anaerobic",
    "METPO:1002005": "Fermentation",
    "METPO:2000103|GO:0030245": "cellulose catabolic process",
    "METPO:2000103|GO:0045493": "xylan catabolic process",
}

TARGET_TAXA = {
    "NCBITaxon:853": "Faecalibacterium prausnitzii",
    "NCBITaxon:39488": "Anaerobutyricum hallii",
    "NCBITaxon:39485": "Lachnospira eligens",
    "NCBITaxon:39491": "Agathobacter rectalis",
    "NCBITaxon:166486": "Roseburia intestinalis",
    "NCBITaxon:40518": "Ruminococcus bromii",
    "NCBITaxon:360807": "Roseburia inulinivorans",
}

GUT_BODY_SITE_RE = re.compile(
    r"feces|faeces|stool|colon|caecum|cecum|gut|intestin|digestive",
    re.IGNORECASE,
)


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def clean_text(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def is_gut_body_site(value: object) -> bool:
    return bool(GUT_BODY_SITE_RE.search(clean_text(value)))


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
    scored["gut_context"] = scored["body_site"].map(is_gut_body_site) | scored[
        "body_site_group"
    ].map(is_gut_body_site)
    scored["guild_z_mean"] = matrix[available_traits].mean(axis=1)
    scored["anaerobe_fermentation_z_mean"] = matrix[
        ["anaerobic", "obligately anaerobic", "Fermentation"]
    ].mean(axis=1)
    scored["fiber_degradation_z_mean"] = matrix[
        ["cellulose catabolic process", "xylan catabolic process"]
    ].mean(axis=1)
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
                "decreased_minus_increased_z": dec_mean - inc_mean,
            }
        )
    for score_name in [
        "guild_z_mean",
        "anaerobe_fermentation_z_mean",
        "fiber_degradation_z_mean",
    ]:
        inc_mean = metadata.loc[metadata["inc_dec"] == "inc", score_name].mean()
        dec_mean = metadata.loc[metadata["inc_dec"] == "dec", score_name].mean()
        rows.append(
            {
                "trait": score_name,
                "increased_profile_mean_z": inc_mean,
                "decreased_profile_mean_z": dec_mean,
                "decreased_minus_increased_z": dec_mean - inc_mean,
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
                "decreased_minus_increased_guild_z": dec["guild_z_mean"].mean()
                - inc["guild_z_mean"].mean(),
                "decreased_mean_fiber_z": dec["fiber_degradation_z_mean"].mean(),
                "decreased_mean_anaerobe_fermentation_z": dec[
                    "anaerobe_fermentation_z_mean"
                ].mean(),
            }
        )
    return pd.DataFrame(records).sort_values(
        ["decreased_mean_guild_z", "profile_count"],
        ascending=[False, False],
    )


def load_taxa_direction_summary() -> pd.DataFrame:
    taxa_metadata = read_tsv(VIS_DIR / "taxa_guild_discovery_taxa_metadata.tsv", index_col=0)
    summary = taxa_metadata.loc[taxa_metadata.index.intersection(TARGET_TAXA.keys())].copy()
    summary["target_taxon_name"] = summary.index.map(TARGET_TAXA)
    summary["decreased_fraction"] = [
        safe_ratio(float(row["decreased"]), float(row["bugsigdb_profile_count"]))
        for _, row in summary.iterrows()
    ]
    summary["depletion_bias"] = [
        safe_ratio(
            float(row["decreased"]) - float(row["increased"]),
            float(row["bugsigdb_profile_count"]),
        )
        for _, row in summary.iterrows()
    ]
    summary["decreased_to_increased_ratio"] = [
        safe_ratio(float(row["decreased"]), float(row["increased"]))
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
            "decreased_fraction",
            "depletion_bias",
            "decreased_to_increased_ratio",
            "trait_count",
        ]
    ].sort_values("depletion_bias", ascending=False)


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
        "gut_context",
        "guild_z_mean",
        "anaerobe_fermentation_z_mean",
        "fiber_degradation_z_mean",
    ]
    profiles = scored_profiles[profile_columns].copy()
    merged = taxa.merge(
        profiles,
        left_on="profile_id",
        right_index=True,
        how="left",
        suffixes=("_taxa", ""),
    )
    return merged


def summarize_taxa_conditions(occurrences: pd.DataFrame) -> pd.DataFrame:
    gut = occurrences[occurrences["gut_context"].fillna(False)].copy()
    grouped = (
        gut.groupby(["target_taxon_name", "condition", "body_site_group", "direction"])[
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
    output["depletion_bias"] = (
        output["decreased_profile_count"] - output["increased_profile_count"]
    ) / output["total_profile_count"].replace(0, pd.NA)
    return output.sort_values(
        ["depletion_bias", "decreased_study_count", "total_profile_count"],
        ascending=[False, False, False],
    )


def make_plots(
    taxa_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    direction_summary: pd.DataFrame,
) -> None:
    plt.rcParams.update({"figure.dpi": 160, "savefig.dpi": 200})

    fig, ax = plt.subplots(figsize=(9, 4.8))
    plot_taxa = taxa_summary.sort_values("depletion_bias")
    ax.barh(plot_taxa["target_taxon_name"], plot_taxa["depletion_bias"], color="#8c2d04")
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("(decreased profiles - increased profiles) / total profiles")
    ax.set_title("Gut Anaerobe Candidate Taxa Are Depletion-Biased")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "candidate_taxa_depletion_bias.png")
    plt.close(fig)

    top_conditions = condition_summary[
        (condition_summary["decreased_profile_count"] >= 3)
        & (condition_summary["study_count"] >= 2)
    ].head(12)
    if top_conditions.empty:
        top_conditions = condition_summary[condition_summary["decreased_profile_count"] >= 3].head(
            12
        )
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.barh(
        top_conditions["condition"][::-1],
        top_conditions["decreased_mean_guild_z"][::-1],
        color="#2563eb",
    )
    ax.set_xlabel("Mean gut anaerobe/fermentation guild z-score in decreased profiles")
    ax.set_title("Conditions With High Decreased-Profile Guild Signal")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "condition_decreased_guild_signal.png")
    plt.close(fig)

    plot_direction = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("decreased_minus_increased_z")
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.barh(
        plot_direction["trait"],
        plot_direction["decreased_minus_increased_z"],
        color="#047857",
    )
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("decreased profile mean z - increased profile mean z")
    ax.set_title("Guild Traits Are Higher in Decreased Profiles")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "guild_trait_direction_shift.png")
    plt.close(fig)


def write_summary_report(
    direction_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    taxa_summary: pd.DataFrame,
    condition_taxa_summary: pd.DataFrame,
    top_depletion_profiles: pd.DataFrame,
) -> None:
    direction_rows = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("decreased_minus_increased_z", ascending=False)
    top_conditions = condition_summary[
        (condition_summary["decreased_profile_count"] >= 3)
        & (condition_summary["study_count"] >= 2)
    ].head(8)
    if top_conditions.empty:
        top_conditions = condition_summary[condition_summary["decreased_profile_count"] >= 3].head(
            8
        )
    top_taxa = taxa_summary.head(7)
    top_taxon_conditions = condition_taxa_summary[
        (condition_taxa_summary["decreased_profile_count"] >= 3)
        & (condition_taxa_summary["decreased_study_count"] >= 2)
    ].head(10)

    lines = [
        "# Gut Anaerobe / Fermentation Depletion Investigation",
        "",
        "This investigation uses the existing study trait profile clustered heatmap matrix, "
        "study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.",
        "",
        "## Scope",
        "",
        "- Guild traits: "
        + ", ".join(f"{name} ({trait_id})" for trait_id, name in GUILD_TRAITS.items()),
        "- Target taxa: " + ", ".join(TARGET_TAXA.values()),
        "- Gut context: body-site labels matching feces/faeces/stool/colon/caecum/cecum/gut/"
        "intestin/digestive, plus matching grouped body-site labels.",
        "",
        "## Main Signal",
        "",
        "The selected guild traits are consistently higher in decreased profiles than in "
        "increased profiles, which supports the gut anaerobe/fermentation depletion hypothesis.",
        "",
        "| trait | decreased - increased z |",
        "| --- | ---: |",
    ]
    for row in direction_rows.itertuples(index=False):
        lines.append(f"| {row.trait} | {row.decreased_minus_increased_z:.3f} |")

    lines.extend(
        [
            "",
            "## Strongest Depletion-Biased Taxa",
            "",
            "| taxon | decreased / increased profiles | depletion bias |",
            "| --- | ---: | ---: |",
        ]
    )
    for row in top_taxa.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {int(row.decreased)} / {int(row.increased)} | "
            f"{row.depletion_bias:.3f} |"
        )

    lines.extend(
        [
            "",
            "## Recurrent Target-Taxon Depletions By Condition",
            "",
            "| taxon | condition | body site | "
            "decreased / increased profiles | decreased studies |",
            "| --- | --- | --- | ---: | ---: |",
        ]
    )
    for row in top_taxon_conditions.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {clean_text(row.condition)} | "
            f"{clean_text(row.body_site_group)} | "
            f"{int(row.decreased_profile_count)} / {int(row.increased_profile_count)} | "
            f"{int(row.decreased_study_count)} |"
        )

    lines.extend(
        [
            "",
            "## Gut Contexts With High Decreased-Profile Guild Signal",
            "",
            "| condition | profiles | studies | decreased mean guild z | decreased - increased z |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in top_conditions.itertuples(index=False):
        diff = row.decreased_minus_increased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        lines.append(
            f"| {row.condition} | {int(row.profile_count)} | {int(row.study_count)} | "
            f"{row.decreased_mean_guild_z:.3f} | {diff_text} |"
        )

    lines.extend(
        [
            "",
            "## Body-Site Check",
            "",
            "| body_site_group | profiles | decreased mean guild z | decreased - increased z |",
            "| --- | ---: | ---: | ---: |",
        ]
    )
    for row in body_site_summary.head(10).itertuples(index=False):
        diff = row.decreased_minus_increased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        lines.append(
            f"| {row.body_site_group} | {int(row.profile_count)} | "
            f"{row.decreased_mean_guild_z:.3f} | {diff_text} |"
        )

    lines.extend(
        [
            "",
            "## Top Decreased Profiles To Inspect",
            "",
            "| profile | condition | body site | guild z | title |",
            "| --- | --- | --- | ---: | --- |",
        ]
    )
    for row in top_depletion_profiles.head(10).itertuples(index=False):
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
            "1. Treat the target taxa as a depletion guild and test enrichment by "
            "condition/body site "
            "rather than analyzing each taxon in isolation.",
            "2. Prioritize conditions where decreased profiles have high guild scores "
            "and target taxa "
            "are repeatedly decreased across multiple studies.",
            "3. Re-run the same summaries after adding explicit SCFA/butyrate "
            "production annotations, "
            "because the current selected heatmap traits capture fermentation and "
            "fiber degradation "
            "but do not fully encode butyrate production.",
            "",
            "## Files",
            "",
            "- `guild_trait_direction_shift.tsv`",
            "- `gut_body_site_guild_summary.tsv`",
            "- `gut_condition_guild_summary.tsv`",
            "- `candidate_taxa_direction_summary.tsv`",
            "- `candidate_taxa_trait_presence.tsv`",
            "- `candidate_taxa_condition_direction.tsv`",
            "- `candidate_taxa_profile_occurrences.tsv`",
            "- `top_decreased_guild_profiles.tsv`",
            "- `candidate_taxa_depletion_bias.png`",
            "- `condition_decreased_guild_signal.png`",
            "- `guild_trait_direction_shift.png`",
            "",
        ]
    )
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    study_matrix, study_metadata, study_trait_metadata = load_study_heatmap()
    scored_profiles = add_guild_scores(study_matrix, study_metadata)
    gut_profiles = scored_profiles[scored_profiles["gut_context"]].copy()

    direction_summary = summarize_direction(scored_profiles)
    body_site_summary = summarize_group(gut_profiles, "body_site_group", min_profiles=10)
    condition_summary = summarize_group(gut_profiles, "condition", min_profiles=5)
    taxa_summary = load_taxa_direction_summary()
    trait_presence = load_target_taxa_trait_counts()
    occurrences = load_target_taxa_occurrences(scored_profiles)
    condition_taxa_summary = summarize_taxa_conditions(occurrences)

    trait_definitions = study_trait_metadata.loc[
        study_trait_metadata.index.intersection(GUILD_TRAITS.keys())
    ].copy()
    trait_definitions["investigation_trait_name"] = trait_definitions.index.map(GUILD_TRAITS)

    top_depletion_profiles = gut_profiles[gut_profiles["inc_dec"] == "dec"].copy()
    top_depletion_profiles = top_depletion_profiles.sort_values(
        ["guild_z_mean", "total_trait_signal"],
        ascending=[False, False],
    ).head(100)
    top_depletion_profiles = top_depletion_profiles.reset_index(names="profile_id")

    direction_summary.to_csv(OUTPUT_DIR / "guild_trait_direction_shift.tsv", sep="\t", index=False)
    body_site_summary.to_csv(OUTPUT_DIR / "gut_body_site_guild_summary.tsv", sep="\t", index=False)
    condition_summary.to_csv(OUTPUT_DIR / "gut_condition_guild_summary.tsv", sep="\t", index=False)
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
    top_depletion_profiles.to_csv(
        OUTPUT_DIR / "top_decreased_guild_profiles.tsv",
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
        top_depletion_profiles,
    )

    print(f"Wrote investigation outputs to: {OUTPUT_DIR}")
    print(f"Study profiles analyzed: {len(scored_profiles)}")
    print(f"Gut-context profiles analyzed: {len(gut_profiles)}")
    print(f"Target taxa occurrence rows: {len(occurrences)}")
    print(f"Condition/taxon rows: {len(condition_taxa_summary)}")


if __name__ == "__main__":
    main()
