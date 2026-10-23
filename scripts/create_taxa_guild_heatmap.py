"""Create a leaf-rank taxa guild discovery heatmap from BugSigDB trait data."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization"
MATRIX_METADATA_COLUMNS = ["taxon_id", "taxon_name", "taxonomic_rank"]
HUMAN_HOST_SPECIES = "Homo sapiens"
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
# The direction legend's title for each run direction_basis (#795). An all-paired run
# keeps the curated Abundance in Group 1 direction; a contrast-population subset flips
# its reversed experiments so that increased means higher in the disease group.
DIRECTION_LEGEND_TITLES = {
    "group_1_vs_group_0": "BugSigDB direction\n(Group 1 vs Group 0)",
    "disease_vs_comparison": "BugSigDB direction\n(disease vs comparison)",
}


def clean_value(value: object) -> str:
    """Normalize empty BugSigDB metadata values."""
    text = str(value or "").strip()
    if not text or text.upper() == "NA" or text.lower() == "nan":
        return ""
    return text


def signature_id_from_bsdb_id(bsdb_id: str) -> str:
    """Convert a full BugSigDB signature ID into the derived signature ID format."""
    return bsdb_id.strip().replace("/", "-")


def parse_metaphlan_name(raw: str) -> tuple[str, str | None]:
    """Parse a MetaPhlAn leaf name such as ``s__Gardnerella vaginalis``."""
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
    tax_ids = clean_value(row.get("NCBI Taxonomy IDs"))
    if not tax_ids:
        return []

    id_chains = [chain.strip() for chain in tax_ids.split(";") if chain.strip()]
    names = clean_value(row.get("MetaPhlAn taxon names")) or clean_value(row.get("Taxon"))
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


def load_corpus_taxon_names(full_dump_path: Path) -> dict[str, dict[str, tuple[str, str]]]:
    """Map each corpus signature in full_dump.csv to its taxa's raw names and ranks.

    Only Homo sapiens 16S and WMS experiments are kept, the analysis corpus's rule; a
    blank sequencing type is left out (#742).
    """
    if not full_dump_path.exists():
        # Without the dump no corpus rule applies; never fall back to every host (#747)
        msg = (
            f"{full_dump_path} is missing; it decides which signatures belong to "
            f"{HUMAN_HOST_SPECIES} 16S or WMS experiments"
        )
        raise FileNotFoundError(msg)

    corpus_taxa: dict[str, dict[str, tuple[str, str]]] = {}
    with full_dump_path.open(encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(rows)
        for row in reader:
            bsdb_id = clean_value(row.get("BSDB ID"))
            host_species = clean_value(row.get("Host species"))
            sequencing_type = clean_value(row.get("Sequencing type"))
            if (
                not bsdb_id
                or host_species != HUMAN_HOST_SPECIES
                or sequencing_type not in CORPUS_SEQUENCING_TYPES
            ):
                continue
            corpus_taxa.setdefault(signature_id_from_bsdb_id(bsdb_id), {}).update(
                {
                    taxon["taxon_id"]: (taxon["taxon_name"], taxon["taxonomic_rank"])
                    for taxon in parse_taxa_from_raw_row(row)
                }
            )
    if not corpus_taxa:
        # An export that yields no corpus row must not silently drop both filters (#753)
        msg = (
            f"{full_dump_path} has no {HUMAN_HOST_SPECIES} 16S or WMS taxon rows; check its "
            "Host species and Sequencing type columns"
        )
        raise ValueError(msg)
    return corpus_taxa


def load_bugsigdb_taxa_records(
    bugsigdb_taxa_path: Path,
    full_dump_path: Path,
) -> pd.DataFrame:
    """Load the run's own increased and decreased taxon rows inside the analysis corpus.

    The rows and their directions come from the run's taxa table, so a contrast
    population summarizes only its own signatures, on its own direction basis (#768).
    The full dump keeps the table's corpus signatures and supplies the raw names and
    ranks, joined on signature and taxon ID.
    """
    corpus_taxa = load_corpus_taxon_names(full_dump_path)
    taxa = pd.read_csv(
        bugsigdb_taxa_path,
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
        dtype=str,
    )
    taxa = taxa[
        taxa["signature_id"].isin(corpus_taxa.keys())
        & taxa["direction"].isin(["increased", "decreased"])
    ].reset_index(drop=True)
    if taxa.empty:
        msg = (
            f"{bugsigdb_taxa_path} has no increased or decreased taxon rows in the "
            f"{HUMAN_HOST_SPECIES} 16S or WMS signatures of {full_dump_path}"
        )
        raise ValueError(msg)

    raw_names = [
        corpus_taxa[signature_id].get(taxon_id, (taxon_name, taxonomic_rank))
        for signature_id, taxon_id, taxon_name, taxonomic_rank in zip(
            taxa["signature_id"],
            taxa["taxon_id"],
            taxa["taxon_name"],
            taxa["taxonomic_rank"],
            strict=True,
        )
    ]
    taxa["taxon_name"] = [taxon_name for taxon_name, _ in raw_names]
    taxa["taxonomic_rank"] = [taxonomic_rank for _, taxonomic_rank in raw_names]
    return taxa


def choose_taxonomic_rank(
    taxa: pd.DataFrame,
    requested_rank: str,
    fallback_ranks: list[str],
) -> str:
    """Select the requested rank, or the first fallback rank present in the data."""
    available_ranks = set(taxa["taxonomic_rank"].dropna().astype(str))
    for candidate in [requested_rank, *fallback_ranks]:
        if candidate in available_ranks:
            return candidate

    available_display = ", ".join(sorted(available_ranks)) or "none"
    requested_display = ", ".join([requested_rank, *fallback_ranks])
    raise ValueError(
        "None of the requested taxonomic ranks are present. "
        f"Requested: {requested_display}. Available: {available_display}."
    )


def load_bugsigdb_taxa_summary(
    bugsigdb_taxa_path: Path,
    rank: str,
    fallback_ranks: list[str],
    full_dump_path: Path,
) -> tuple[pd.DataFrame, str, int, int]:
    """Summarize BugSigDB occurrence and direction counts by taxon."""
    taxa = load_bugsigdb_taxa_records(bugsigdb_taxa_path, full_dump_path)
    studies_analyzed = taxa["study_id"].nunique()
    profiles_analyzed = taxa[["study_id", "experiment_id", "direction"]].drop_duplicates().shape[0]
    selected_rank = choose_taxonomic_rank(taxa, rank, fallback_ranks)
    taxa = taxa[taxa["taxonomic_rank"] == selected_rank].copy()
    taxa["profile_id"] = (
        taxa["study_id"].fillna("")
        + "|"
        + taxa["experiment_id"].fillna("")
        + "|"
        + taxa["direction"].fillna("")
    )

    summary = taxa.groupby("taxon_id", sort=False).agg(
        taxon_name=("taxon_name", "first"),
        taxonomic_rank=("taxonomic_rank", "first"),
        bugsigdb_study_count=("study_id", "nunique"),
        bugsigdb_profile_count=("profile_id", "nunique"),
    )

    direction_counts = (
        taxa[["taxon_id", "profile_id", "direction"]]
        .drop_duplicates()
        .pivot_table(
            index="taxon_id",
            columns="direction",
            values="profile_id",
            aggfunc="nunique",
            fill_value=0,
        )
    )
    summary = summary.join(direction_counts, how="left").fillna(0)
    for column in ["increased", "decreased"]:
        if column not in summary.columns:
            summary[column] = 0
        summary[column] = summary[column].astype(int)
    summary["dominant_direction"] = summary.apply(classify_direction, axis=1)
    return summary, selected_rank, studies_analyzed, profiles_analyzed


def classify_direction(row: pd.Series) -> str:
    """Classify whether a taxon is mostly increased, decreased, or mixed."""
    increased = int(row.get("increased", 0))
    decreased = int(row.get("decreased", 0))
    if increased == 0 and decreased == 0:
        return "other"
    if increased >= max(2, decreased * 1.25):
        return "increased"
    if decreased >= max(2, increased * 1.25):
        return "decreased"
    return "mixed"


def load_profile_support(study_trait_profiles_path: Path, trait_columns: list[str]) -> pd.Series:
    """Count how many study profiles have each trait."""
    available_columns = pd.read_csv(study_trait_profiles_path, sep="\t", nrows=0).columns
    selected_columns = [column for column in trait_columns if column in available_columns]
    support = pd.Series(0, index=selected_columns, dtype="int64")
    if not selected_columns:
        return support

    for chunk in pd.read_csv(
        study_trait_profiles_path,
        sep="\t",
        usecols=selected_columns,
        chunksize=100,
    ):
        support = support.add(chunk.gt(0).sum(axis=0).astype("int64"), fill_value=0).astype("int64")
    return support


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
    """Build a compact two-line heatmap trait label."""
    name = row.get("object_name", "") or row.get("column_label", "") or column_id
    name = str(name).replace("_", " ")
    if len(name) > 36:
        name = f"{name[:33]}..."
    short_id = column_id.split("|")[-1].replace("pathways:", "").replace("carbon_substrates:", "")
    if len(short_id) > 24:
        short_id = f"{short_id[:21]}..."
    return f"{name}\n{short_id}"


def select_taxa(
    trait_matrix: pd.DataFrame,
    taxa_summary: pd.DataFrame,
    trait_columns: list[str],
    max_taxa: int,
    min_profiles: int,
    min_traits: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select BugSigDB taxa with enough profile and trait support."""
    matrix = trait_matrix[trait_matrix["taxon_id"].isin(taxa_summary.index)].copy()
    matrix = matrix.drop_duplicates(subset=["taxon_id"], keep="first")
    matrix = matrix.set_index("taxon_id", drop=False)
    values = matrix[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    binary = values.gt(0).astype(int)

    metadata = taxa_summary.join(
        pd.DataFrame({"trait_count": binary.sum(axis=1)}, index=binary.index),
        how="inner",
    )
    metadata = metadata[
        (metadata["bugsigdb_profile_count"] >= min_profiles)
        & (metadata["trait_count"] >= min_traits)
    ].copy()
    metadata = metadata.sort_values(
        ["bugsigdb_profile_count", "bugsigdb_study_count", "trait_count"],
        ascending=False,
    ).head(max_taxa)
    if len(metadata) < 2:
        msg = "Fewer than two taxa passed the selection filters."
        raise ValueError(msg)
    return binary.loc[metadata.index], metadata


def select_traits(
    binary_taxa_traits: pd.DataFrame,
    trait_labels: pd.DataFrame,
    study_profile_support: pd.Series,
    max_traits: int,
    min_taxa_support: int,
    min_profile_support: int,
    include_media_traits: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select informative traits supported by both taxa and study-profile data."""
    taxa_support = binary_taxa_traits.sum(axis=0)
    trait_summary = pd.DataFrame(
        {
            "trait_id": taxa_support.index,
            "bugsigdb_taxa_support": taxa_support.astype(int),
            "study_profile_support": study_profile_support.reindex(taxa_support.index)
            .fillna(0)
            .astype(int),
        },
        index=taxa_support.index,
    )
    trait_summary = trait_summary.join(
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
    trait_summary[["display_label", "display_category"]] = trait_summary[
        ["display_label", "display_category"]
    ].fillna("")
    trait_summary["selection_score"] = [
        int(taxa_count) * math.log1p(int(profile_count))
        for taxa_count, profile_count in zip(
            trait_summary["bugsigdb_taxa_support"],
            trait_summary["study_profile_support"],
            strict=True,
        )
    ]

    eligible = trait_summary[
        (trait_summary["bugsigdb_taxa_support"] >= min_taxa_support)
        & (trait_summary["study_profile_support"] >= min_profile_support)
    ].copy()
    if not include_media_traits:
        eligible = eligible[
            ~eligible.index.to_series().str.contains("mediadive.medium", regex=False)
            & ~eligible["display_label"].str.contains("mediadive.medium", regex=False)
        ].copy()
    eligible = eligible.sort_values(
        ["selection_score", "bugsigdb_taxa_support", "study_profile_support"],
        ascending=False,
    ).head(max_traits)
    if len(eligible) < 2:
        msg = "Fewer than two traits passed the selection filters."
        raise ValueError(msg)

    selected_matrix = binary_taxa_traits[eligible.index].copy()
    selected_matrix = selected_matrix.loc[selected_matrix.sum(axis=1) > 0]
    return selected_matrix, eligible


def build_display_matrix(
    selected_matrix: pd.DataFrame,
    taxa_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Apply readable row and column labels to the selected matrix."""
    display_matrix = selected_matrix.copy()
    display_matrix.index = [
        f"{taxa_metadata.loc[taxon_id, 'taxon_name']} ({taxon_id})"
        for taxon_id in display_matrix.index
    ]
    display_matrix.columns = [
        trait_metadata.loc[trait_id, "display_label"] or trait_id
        for trait_id in display_matrix.columns
    ]
    return display_matrix


def save_heatmap(
    display_matrix: pd.DataFrame,
    taxa_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
    selected_rank: str,
    requested_rank: str,
    studies_analyzed: int,
    profiles_analyzed: int,
    direction_basis: str,
) -> None:
    """Render and save the clustered heatmap, naming the taxa table's direction basis."""
    sns.set_theme(style="white")

    direction_palette = {
        "increased": "#c0392b",
        "decreased": "#2874a6",
        "mixed": "#7f8c8d",
        "other": "#8e44ad",
    }
    trait_categories = sorted(trait_metadata["display_category"].fillna("unknown").unique())
    category_palette = dict(
        zip(
            trait_categories,
            sns.color_palette("Set2", n_colors=len(trait_categories)),
            strict=True,
        )
    )
    row_taxon_ids = display_matrix.index.map(extract_taxon_id)
    row_colors = taxa_metadata.loc[row_taxon_ids, "dominant_direction"].map(direction_palette)
    row_colors.index = display_matrix.index
    trait_category_by_label = trait_metadata.set_index("display_label")["display_category"]
    col_colors = trait_category_by_label.loc[display_matrix.columns].map(category_palette)
    col_colors.index = display_matrix.columns

    cluster = sns.clustermap(
        display_matrix,
        cmap=sns.color_palette(["#f8f9fa", "#0f766e"], as_cmap=True),
        method="average",
        metric="jaccard",
        row_colors=row_colors,
        col_colors=col_colors,
        linewidths=0.15,
        linecolor="#e5e7eb",
        xticklabels=True,
        yticklabels=True,
        figsize=(24, 18),
        cbar_pos=(0.02, 0.82, 0.02, 0.12),
        dendrogram_ratio=(0.12, 0.10),
        colors_ratio=(0.025, 0.025),
    )

    cluster.ax_heatmap.set_xlabel("Trait")
    cluster.ax_heatmap.set_ylabel(f"BugSigDB {selected_rank} taxon")
    rank_note = (
        f"requested {requested_rank}; using {selected_rank}"
        if selected_rank != requested_rank
        else f"{selected_rank} rank"
    )
    cluster.ax_heatmap.set_title(
        "Taxa Guild Discovery: BugSigDB Taxa Clustered by Trait Vectors\n"
        f"Rows are high-support human-host 16S and shotgun taxa ({rank_note}); columns are "
        "binary "
        "KG-Microbe traits; "
        f"{studies_analyzed} studies and {profiles_analyzed} study/experiment/direction "
        "profiles analyzed",
        pad=18,
    )
    cluster.ax_cbar.set_ylabel("Trait present")
    for label in cluster.ax_heatmap.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(7)
    for label in cluster.ax_heatmap.get_yticklabels():
        label.set_fontsize(7)

    direction_handles = [
        Patch(color=color, label=label)
        for label, color in direction_palette.items()
        if label in set(taxa_metadata["dominant_direction"])
    ]
    category_handles = [
        Patch(color=color, label=label) for label, color in category_palette.items()
    ]
    # Left of the row dendrogram, whose root line the two-line title crossed (#896)
    cluster.ax_row_dendrogram.legend(
        handles=direction_handles,
        title=DIRECTION_LEGEND_TITLES[direction_basis],
        loc="center right",
        bbox_to_anchor=(-0.02, 0.5),
        frameon=False,
    )
    cluster.ax_col_dendrogram.legend(
        handles=category_handles,
        title="Trait category",
        loc="center",
        bbox_to_anchor=(0.5, 1.25),
        ncol=3,
        frameon=False,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster.fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(cluster.fig)


def extract_taxon_id(display_label: str) -> str:
    """Extract the taxon ID from a display label."""
    return display_label.rsplit("(", 1)[-1].rstrip(")")


def write_outputs(
    selected_matrix: pd.DataFrame,
    display_matrix: pd.DataFrame,
    taxa_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
    matrix_output_path: Path,
    taxa_metadata_output_path: Path,
    trait_metadata_output_path: Path,
    selected_rank: str,
    requested_rank: str,
    studies_analyzed: int,
    profiles_analyzed: int,
) -> None:
    """Write the selected matrix and metadata tables."""
    matrix_output_path.parent.mkdir(parents=True, exist_ok=True)
    display_matrix.to_csv(matrix_output_path, sep="\t")
    taxa_metadata.loc[selected_matrix.index].to_csv(taxa_metadata_output_path, sep="\t")
    trait_metadata.loc[selected_matrix.columns].to_csv(trait_metadata_output_path, sep="\t")
    print(f"Wrote heatmap: {output_path}")
    print(f"Wrote matrix: {matrix_output_path}")
    print(f"Wrote taxa metadata: {taxa_metadata_output_path}")
    print(f"Wrote trait metadata: {trait_metadata_output_path}")
    print(f"Requested rank: {requested_rank}")
    print(f"Taxonomic rank plotted: {selected_rank}")
    print(f"BugSigDB studies analyzed: {studies_analyzed}")
    print(f"BugSigDB profiles analyzed: {profiles_analyzed}")
    print(f"Taxa plotted: {len(selected_matrix)}")
    print(f"Traits plotted: {len(selected_matrix.columns)}")


def build_parser(
    parser_class: type[argparse.ArgumentParser] = argparse.ArgumentParser,
) -> argparse.ArgumentParser:
    """Expose the execution parser to verified arm planning and replay."""
    parser = parser_class(
        description=(
            "Create a taxa guild discovery clustered heatmap by clustering BugSigDB taxa "
            "from the taxa-trait matrix and filtering traits with study-profile support."
        )
    )
    parser.add_argument(
        "--bugsigdb-taxa-input",
        type=Path,
        default=REPO_ROOT / "data/derived/bugsigdb_taxa.tsv",
        help=(
            "Input BugSigDB taxa TSV path: the run's own taxon rows and directions, so a "
            "contrast-population subset plots only its own signatures."
        ),
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
        help="Input taxa trait matrix column metadata TSV path.",
    )
    parser.add_argument(
        "--full-dump-input",
        type=Path,
        default=REPO_ROOT / "data/raw/full_dump.csv",
        help=(
            "Input BugSigDB full_dump.csv path (required). It keeps the taxa TSV's Homo "
            "sapiens 16S and WMS signatures, and its parsed taxon names/ranks replace the "
            "taxa TSV's."
        ),
    )
    parser.add_argument(
        "--study-trait-profiles-input",
        type=Path,
        default=REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv",
        help="Input wide study trait profiles TSV path.",
    )
    parser.add_argument(
        "--direction-basis",
        choices=sorted(DIRECTION_LEGEND_TITLES),
        required=True,
        help=(
            "Direction basis of the taxa TSV, the core run's direction_basis: "
            "group_1_vs_group_0 for all paired contrasts, disease_vs_comparison for a "
            "contrast-population subset. The direction legend names it. There is no "
            "default, so leaving it out cannot label a subset's directions Group 1 vs "
            "Group 0 (#840)."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIR / "taxa_guild_discovery_clustered_heatmap.png",
        help="Output heatmap image path.",
    )
    parser.add_argument(
        "--matrix-output",
        type=Path,
        default=OUTPUT_DIR / "taxa_guild_discovery_clustered_heatmap_matrix.tsv",
        help="Output selected taxa x trait matrix path.",
    )
    parser.add_argument(
        "--taxa-metadata-output",
        type=Path,
        default=OUTPUT_DIR / "taxa_guild_discovery_taxa_metadata.tsv",
        help="Output selected taxa metadata path.",
    )
    parser.add_argument(
        "--trait-metadata-output",
        type=Path,
        default=OUTPUT_DIR / "taxa_guild_discovery_trait_metadata.tsv",
        help="Output selected trait metadata path.",
    )
    parser.add_argument(
        "--rank",
        default="strain",
        help="Taxonomic rank to cluster.",
    )
    parser.add_argument(
        "--fallback-ranks",
        nargs="*",
        default=["species"],
        help="Fallback ranks to use when --rank is absent from the BugSigDB export.",
    )
    parser.add_argument(
        "--max-taxa",
        type=int,
        default=80,
        help="Maximum number of BugSigDB taxa to plot.",
    )
    parser.add_argument(
        "--max-traits",
        type=int,
        default=60,
        help="Maximum number of traits to plot.",
    )
    parser.add_argument(
        "--min-profiles",
        type=int,
        default=4,
        help="Minimum BugSigDB profile count required for a taxon.",
    )
    parser.add_argument(
        "--min-traits",
        type=int,
        default=5,
        help="Minimum trait count required for a taxon.",
    )
    parser.add_argument(
        "--min-taxa-support",
        type=int,
        default=4,
        help="Minimum selected-taxa support required for a trait.",
    )
    parser.add_argument(
        "--min-profile-support",
        type=int,
        default=4,
        help="Minimum study-profile support required for a trait.",
    )
    parser.add_argument(
        "--include-media-traits",
        action="store_true",
        help="Include MediaDive growth-medium traits in the selected traits.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    taxa_summary, selected_rank, studies_analyzed, profiles_analyzed = load_bugsigdb_taxa_summary(
        args.bugsigdb_taxa_input,
        args.rank,
        args.fallback_ranks,
        args.full_dump_input,
    )
    trait_matrix = pd.read_csv(args.taxa_trait_matrix_input, sep="\t")
    trait_columns = [
        column for column in trait_matrix.columns if column not in MATRIX_METADATA_COLUMNS
    ]
    trait_labels = load_trait_labels(args.taxa_trait_columns_input)
    study_profile_support = load_profile_support(args.study_trait_profiles_input, trait_columns)

    selected_taxa_traits, taxa_metadata = select_taxa(
        trait_matrix,
        taxa_summary,
        trait_columns,
        args.max_taxa,
        args.min_profiles,
        args.min_traits,
    )
    selected_matrix, trait_metadata = select_traits(
        selected_taxa_traits,
        trait_labels,
        study_profile_support,
        args.max_traits,
        args.min_taxa_support,
        args.min_profile_support,
        args.include_media_traits,
    )
    display_matrix = build_display_matrix(selected_matrix, taxa_metadata, trait_metadata)
    save_heatmap(
        display_matrix,
        taxa_metadata,
        trait_metadata,
        args.output,
        selected_rank,
        args.rank,
        studies_analyzed,
        profiles_analyzed,
        args.direction_basis,
    )
    write_outputs(
        selected_matrix,
        display_matrix,
        taxa_metadata,
        trait_metadata,
        args.output,
        args.matrix_output,
        args.taxa_metadata_output,
        args.trait_metadata_output,
        selected_rank,
        args.rank,
        studies_analyzed,
        profiles_analyzed,
    )


if __name__ == "__main__":
    main()
