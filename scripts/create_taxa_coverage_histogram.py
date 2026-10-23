# ruff: noqa: E501
"""Generate human 16S/WMS taxa coverage figures, tables and provenance (#820/#821).

Inputs are the rebuilt profile table, taxon matrix, taxon occurrences and experiment
metadata. Both profiles and occurrences are filtered explicitly under #742; no
contrast-tier or pairing selection is applied. All counts, plot annotations and the
README are computed from the selected data. Coverage means recorded annotations,
not measured physiology. Zero-trait profiles remain in this coverage audit; only
the explicitly named mapped-profile panels omit them.
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES, HUMAN_HOST_SPECIES
from mtp.paths import file_sha256

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"
OUTPUT_DIR = INVESTIGATION_DIR / "taxa-coverage"

TAXA_PATH = REPO_ROOT / "data" / "derived" / "bugsigdb_taxa.tsv"
TAXA_TRAIT_MATRIX_PATH = REPO_ROOT / "data" / "derived" / "taxa_trait_matrix.tsv"
PROFILES_WIDE_PATH = REPO_ROOT / "data" / "derived" / "study_trait_profiles_wide.tsv"
EXPERIMENTS_PATH = REPO_ROOT / "data" / "derived" / "experiments.tsv"

# Standard taxonomic hierarchy (ascending order)
RANK_ORDER = ["species", "genus", "family", "order", "class", "phylum", "superkingdom", "other"]
RANK_DISPLAY = {
    "species": "Species",
    "genus": "Genus",
    "family": "Family",
    "order": "Order",
    "class": "Class",
    "phylum": "Phylum",
    "superkingdom": "Superkingdom",
    "other": "Other / unrecorded",
}

# Colorblind-safe palette (Okabe-Ito & Viridis inspired)
COLOR_RICH = "#0072b2"  # Deep Blue: >= 3 traits / All-mapped
COLOR_SPARSE = "#56b4e9"  # Sky Blue: 1-2 traits / Mixed mapped+zero
COLOR_ZERO = "#d55e00"  # Vermilion: 0 traits / Zero-only
COLOR_ACCENT = "#009e73"  # Bluish Green for highlights
COLOR_PURPLE = "#cc79a7"  # Reddish Purple
COLOR_AMBER = "#e69f00"  # Orange/Amber
COLOR_SLATE_DARK = "#0f172a"
COLOR_SLATE_MID = "#475569"
COLOR_SLATE_LIGHT = "#94a3b8"
COLOR_BG_LIGHT = "#f8fafc"

RANK_COLORS = {
    "species": "#0072b2",  # Deep Blue
    "genus": "#009e73",  # Bluish Green
    "family": "#56b4e9",  # Sky Blue
    "order": "#e69f00",  # Amber
    "class": "#cc79a7",  # Purple
    "phylum": "#f0e442",  # Yellow
    "superkingdom": "#64748b",  # Slate
    "other": "#94a3b8",
}

plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["svg.hashsalt"] = "taxa-coverage-histogram-v2"


def read_tsv(path: Path, **kwargs: Any) -> pd.DataFrame:
    """Read a local TSV file."""
    return pd.read_csv(path, sep="\t", **kwargs)


def write_tsv(frame: pd.DataFrame, output_path: Path) -> None:
    """Write a DataFrame to a TSV file with parent directories created."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, sep="\t", index=False)


def save_figure(fig: plt.Figure, output_path: Path) -> None:
    """Save a figure as both 300 DPI PNG and clean SVG."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output_path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def diagnose_zero_trait_taxon(taxon_id: str, name: str, rank: str) -> tuple[str, str]:
    """Diagnose the informatics root cause and remediation for a zero-trait BugSigDB taxon."""
    lower_name = name.lower()

    if taxon_id in {"NCBITaxon:541000", "NCBITaxon:81850", "NCBITaxon:508215"}:
        canonical_map = {
            "NCBITaxon:541000": "NCBITaxon:216572 (Oscillospiraceae)",
            "NCBITaxon:81850": "NCBITaxon:33958 (Lactobacillaceae)",
            "NCBITaxon:508215": "NCBITaxon:32207 (Rothia)",
        }
        twin = canonical_map.get(taxon_id, "canonical TaxID")
        return (
            "NCBI TaxID Desynchronization",
            f"Review possible alias against the pinned taxonomy: {twin}; do not remap automatically",
        )

    if "/" in name or re.search(r"\b(unidentified|uncultured organism|unknown)\b", lower_name):
        return (
            "Ambiguous Taxon / Bin Slash",
            "Resolve to lowest common ancestor (LCA) or treat as composite multi-taxon marker",
        )

    if "candidatus" in lower_name or lower_name.startswith("ca. "):
        return (
            "Candidatus Lineage",
            "Infer phenotypic traits from MAG-derived metabolic models or phylogenetic neighbors",
        )

    if re.search(
        r"\b(candida|malassezia|ascomycota|saccharomyces|fungus|fungi|basidiomycota|streptophyta|viridiplantae)\b",
        lower_name,
    ) or (rank in {"kingdom", "superkingdom", "phylum"} and "fung" in lower_name):
        return (
            "Non-Bacterial / Fungal Eukaryote",
            "Integrate specialized microeukaryotic trait ontologies (e.g. FungiTraits) into KG-Microbe",
        )

    if (
        "uncultured" in lower_name
        or re.search(r"\b(cag:\d+|uba\d+|ucg-\d+|rf39|gut metagenome)\b", lower_name)
        or re.search(r"\s+bacterium\b", lower_name)
        or re.search(r"\s+sp\.\s*$", lower_name)
    ):
        return (
            "Uncultivated / Metagenomic MAG",
            "Seek direct evidence; higher-rank annotations do not establish leaf-taxon traits",
        )

    return (
        "Rare / Novel Isolate (Lacks KG Phenotype)",
        "Target for literature curation in BacDive / KG-Microbe phenotypic assertions",
    )


def percent(numerator: float, denominator: float) -> float:
    """A display percentage, with zero for an empty denominator."""
    return 100 * numerator / denominator if denominator else 0.0


def profile_keys(frame: pd.DataFrame) -> pd.Series:
    """Use the same study/experiment/direction unit as the canonical builder."""
    return frame[["study_id", "experiment_id", "inc_dec"]].astype(str).agg("|".join, axis=1)


def select_corpus(
    taxa: pd.DataFrame, matrix: pd.DataFrame, profiles: pd.DataFrame, experiments: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Filter both granularities to human 16S/WMS and reject unreconciled inputs.

    There is no paired-contrast or disease/control restriction. Matrix columns are
    restricted to the profile vocabulary, including MediaDive if present. Missing
    metadata, taxa, profiles, or disagreeing aggregate counts are errors.
    """
    keys = ["study_id", "experiment_id"]
    if experiments.duplicated(keys).any():
        raise ValueError("Duplicate experiment metadata")
    metadata = experiments[[*keys, "host_species", "sequencing_type"]].copy()
    for column in ("host_species", "sequencing_type"):
        metadata[column] = metadata[column].fillna("").str.strip()
    selected = []
    for frame in (taxa, profiles):
        joined = frame.merge(metadata, on=keys, how="left", validate="many_to_one", indicator=True)
        if joined["_merge"].eq("left_only").any():
            raise ValueError("Input has an experiment without host/sequencing metadata")
        keep = joined.host_species.eq(HUMAN_HOST_SPECIES) & joined.sequencing_type.isin(
            CORPUS_SEQUENCING_TYPES
        )
        selected.append(joined.loc[keep, frame.columns].copy())
    taxa, profiles = selected
    if profiles.empty:
        raise ValueError("No profiles in the human 16S/WMS population")
    if taxa.empty:
        raise ValueError("No reported taxa in the human 16S/WMS population")
    taxa["inc_dec"] = (
        taxa.direction.str.strip()
        .str.lower()
        .map({"increased": "inc", "inc": "inc", "decreased": "dec", "dec": "dec"})
    )
    if taxa.inc_dec.isna().any() or not profiles.inc_dec.isin(["inc", "dec"]).all():
        raise ValueError("Unrecognized profile direction")
    profile_ids = profile_keys(profiles)
    if profile_ids.duplicated().any():
        raise ValueError("Duplicate profile keys")
    taxa_ids = profile_keys(taxa)
    if set(taxa_ids) - set(profile_ids):
        raise ValueError("Canonical taxon occurrences are missing from the profile table")
    if matrix.taxon_id.duplicated().any():
        raise ValueError("Duplicate matrix taxon IDs")
    if set(taxa.taxon_id) - set(matrix.taxon_id):
        raise ValueError("Reported taxa are missing from the trait matrix")
    trait_cols = [c for c in profiles if c not in {*keys, "inc_dec", "taxa_count"}]
    if set(trait_cols) - set(matrix.columns):
        raise ValueError("Profile trait columns are missing from the matrix")
    matrix = matrix.loc[
        matrix.taxon_id.isin(taxa.taxon_id),
        ["taxon_id", "taxon_name", "taxonomic_rank", *trait_cols],
    ].copy()
    for frame in (matrix, profiles):
        values = frame[trait_cols].to_numpy(dtype=float)
        if (
            not np.isfinite(values).all()
            or (values < 0).any()
            or (values != np.floor(values)).any()
        ):
            raise ValueError("Trait counts must be finite nonnegative integers")
    distinct = taxa.assign(profile_key=taxa_ids)[["profile_key", "taxon_id"]].drop_duplicates()
    contributions = distinct.merge(
        matrix[["taxon_id", *trait_cols]], on="taxon_id", validate="many_to_one"
    )
    expected = (
        contributions.groupby("profile_key")[trait_cols].sum().reindex(profile_ids, fill_value=0)
    )
    if not np.array_equal(expected.to_numpy(), profiles[trait_cols].to_numpy()):
        raise ValueError("Profile counts do not reconcile to the reported taxa and matrix")
    expected_taxa = distinct.groupby("profile_key").size().reindex(profile_ids, fill_value=0)
    if "taxa_count" in profiles and not np.array_equal(
        expected_taxa.to_numpy(), profiles.taxa_count.to_numpy()
    ):
        raise ValueError("Profile taxa_count does not reconcile to distinct reported taxa")
    # The matrix supplies one authoritative rank for each reported taxon.
    taxa = taxa.drop(columns="taxonomic_rank").merge(
        matrix[["taxon_id", "taxonomic_rank"]], on="taxon_id", validate="many_to_one"
    )
    return taxa, matrix, profiles


def load_and_process_data(
    *,
    taxa: Path = TAXA_PATH,
    matrix: Path = TAXA_TRAIT_MATRIX_PATH,
    profiles: Path = PROFILES_WIDE_PATH,
    experiments: Path = EXPERIMENTS_PATH,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load taxa.tsv, taxa_trait_matrix.tsv, and study_trait_profiles_wide.tsv and compute statistics."""
    print(f"Reading taxa occurrences from {taxa}...")
    taxa_df = read_tsv(taxa)

    print(f"Reading taxa trait matrix from {matrix}...")
    matrix_df = read_tsv(matrix)

    print(f"Reading study trait profiles from {profiles}...")
    profiles_df = read_tsv(profiles)
    input_profile_count = len(profiles_df)
    input_occurrence_count = len(taxa_df)
    taxa_df, matrix_df, profiles_df = select_corpus(
        taxa_df, matrix_df, profiles_df, read_tsv(experiments)
    )

    meta_cols = ["taxon_id", "taxon_name", "taxonomic_rank"]
    trait_cols = [c for c in matrix_df.columns if c not in meta_cols]
    non_mediadive_cols = [c for c in trait_cols if not c.startswith("mediadive.medium:")]
    print(
        f"  Loaded {len(matrix_df)} unique taxa and {len(trait_cols)} trait columns ({len(non_mediadive_cols)} core non-MediaDive)."
    )

    # Vectorized trait counts per taxon
    trait_counts = (matrix_df[trait_cols] > 0).sum(axis=1).to_numpy()
    core_trait_counts = (matrix_df[non_mediadive_cols] > 0).sum(axis=1).to_numpy()
    taxa_master = matrix_df[meta_cols].copy()
    taxa_master["trait_count"] = trait_counts
    taxa_master["core_trait_count"] = core_trait_counts

    # Standardize rank strings and profile keys
    taxa_master["taxonomic_rank"] = taxa_master["taxonomic_rank"].str.lower().str.strip()
    taxa_df["taxonomic_rank"] = taxa_df["taxonomic_rank"].str.lower().str.strip()
    for frame in (taxa_master, taxa_df):
        frame["taxonomic_rank"] = frame.taxonomic_rank.where(
            frame.taxonomic_rank.isin(RANK_ORDER), "other"
        )
    taxa_df["inc_dec"] = (
        taxa_df["direction"]
        .str.lower()
        .str.strip()
        .map({"increased": "inc", "inc": "inc", "decreased": "dec", "dec": "dec"})
    )
    taxa_df["profile_key"] = (
        taxa_df["study_id"].astype(str)
        + "|"
        + taxa_df["experiment_id"].astype(str)
        + "|"
        + taxa_df["inc_dec"].astype(str)
    )

    # Calculate occurrence frequency, distinct profile counts, and distinct study counts in BugSigDB
    occ_stats = (
        taxa_df.groupby("taxon_id")
        .agg(
            signature_occurrences=("signature_id", "count"),
            profile_count=("profile_key", "nunique"),
            study_count=("study_id", "nunique"),
        )
        .reset_index()
    )

    taxa_master = taxa_master.merge(occ_stats, on="taxon_id", how="left")
    for col in ["signature_occurrences", "profile_count", "study_count"]:
        taxa_master[col] = taxa_master[col].fillna(0).astype(int)

    # Categorize trait coverage
    taxa_master["has_traits"] = taxa_master["trait_count"] > 0
    conditions = [
        taxa_master["trait_count"] == 0,
        (taxa_master["trait_count"] >= 1) & (taxa_master["trait_count"] <= 2),
        taxa_master["trait_count"] >= 3,
    ]
    choices = ["zero_traits", "sparse_1_2", "rich_gte_3"]
    taxa_master["coverage_category"] = np.select(conditions, choices, default="unknown")

    # Process profile-level data
    prof_meta_cols = [
        c
        for c in ["study_id", "experiment_id", "inc_dec", "taxa_count"]
        if c in profiles_df.columns
    ]
    prof_trait_cols = [c for c in profiles_df.columns if c not in prof_meta_cols]
    prof_core_cols = [c for c in prof_trait_cols if not c.startswith("mediadive.medium:")]

    prof_summary = profiles_df[prof_meta_cols].copy()
    prof_summary["profile_key"] = (
        prof_summary["study_id"].astype(str)
        + "|"
        + prof_summary["experiment_id"].astype(str)
        + "|"
        + prof_summary["inc_dec"].astype(str)
    )
    prof_summary["mapped_trait_count"] = (profiles_df[prof_trait_cols] > 0).sum(axis=1).to_numpy()
    prof_summary["mapped_core_trait_count"] = (
        (profiles_df[prof_core_cols] > 0).sum(axis=1).to_numpy()
    )
    prof_summary["total_trait_annotations"] = profiles_df[prof_trait_cols].sum(axis=1).to_numpy()

    # Join taxa_df with taxa_master and deduplicate per (profile_key, taxon_id)
    merged_occ = taxa_df.merge(
        taxa_master[
            ["taxon_id", "trait_count", "core_trait_count", "has_traits", "coverage_category"]
        ],
        on="taxon_id",
        how="left",
    )
    prof_taxa = merged_occ.groupby(
        ["profile_key", "taxon_id", "taxonomic_rank"], as_index=False
    ).agg(
        trait_count=("trait_count", "first"),
        core_trait_count=("core_trait_count", "first"),
        has_traits=("has_traits", "first"),
        coverage_category=("coverage_category", "first"),
        sig_mentions=("signature_id", "count"),
    )

    # Overall taxa counts per profile
    taxa_by_prof = (
        prof_taxa.groupby("profile_key")
        .agg(
            total_unique_taxa=("taxon_id", "count"),
            mapped_trait_taxa=("has_traits", lambda s: int(s.sum())),
            zero_trait_taxa=("has_traits", lambda s: int((~s).sum())),
            sum_taxon_traits=("trait_count", "sum"),
        )
        .reset_index()
    )
    prof_summary = prof_summary.merge(taxa_by_prof, on="profile_key", how="left")
    for col in ["total_unique_taxa", "mapped_trait_taxa", "zero_trait_taxa", "sum_taxon_traits"]:
        prof_summary[col] = prof_summary[col].fillna(0).astype(int)

    # Also compute per-rank mapped & zero taxa counts per profile
    for rank in RANK_ORDER:
        sub_r = prof_taxa[prof_taxa["taxonomic_rank"] == rank]
        r_agg = (
            sub_r.groupby("profile_key")
            .agg(
                **{
                    f"taxa_{rank}": ("taxon_id", "count"),
                    f"mapped_taxa_{rank}": ("has_traits", lambda s: int(s.sum())),
                    f"zero_taxa_{rank}": ("has_traits", lambda s: int((~s).sum())),
                }
            )
            .reset_index()
        )
        prof_summary = prof_summary.merge(r_agg, on="profile_key", how="left")
        for col in [f"taxa_{rank}", f"mapped_taxa_{rank}", f"zero_taxa_{rank}"]:
            prof_summary[col] = prof_summary[col].fillna(0).astype(int)

    prof_summary["has_mapped_traits"] = prof_summary["mapped_trait_count"] > 0

    # Build rank-stratified summary (both taxa-level and profile-level)
    rank_rows: list[dict[str, Any]] = []
    prof_rank_rows: list[dict[str, Any]] = []
    total_profiles = len(prof_summary)
    total_mapped_profiles = int(prof_summary["has_mapped_traits"].sum())

    for rank in RANK_ORDER:
        sub_taxa = taxa_master[taxa_master["taxonomic_rank"] == rank]
        sub_occ = merged_occ[merged_occ["taxonomic_rank"] == rank]

        total_taxa = len(sub_taxa)
        zero_taxa = int((sub_taxa["coverage_category"] == "zero_traits").sum())
        sparse_taxa = int((sub_taxa["coverage_category"] == "sparse_1_2").sum())
        rich_taxa = int((sub_taxa["coverage_category"] == "rich_gte_3").sum())
        trait_taxa = sparse_taxa + rich_taxa

        total_occ = len(sub_occ)
        zero_occ = int((sub_occ["coverage_category"] == "zero_traits").sum())
        sparse_occ = int((sub_occ["coverage_category"] == "sparse_1_2").sum())
        rich_occ = int((sub_occ["coverage_category"] == "rich_gte_3").sum())
        trait_occ = sparse_occ + rich_occ

        # Profile-level statistics for this rank
        profs_with_rank = prof_summary[prof_summary[f"taxa_{rank}"] > 0]
        n_profs_rank = len(profs_with_rank)
        n_profs_all_mapped = int((profs_with_rank[f"zero_taxa_{rank}"] == 0).sum())
        n_profs_mixed = int(
            (
                (profs_with_rank[f"mapped_taxa_{rank}"] > 0)
                & (profs_with_rank[f"zero_taxa_{rank}"] > 0)
            ).sum()
        )
        n_profs_zero_only = int((profs_with_rank[f"mapped_taxa_{rank}"] == 0).sum())
        n_profs_any_mapped = n_profs_all_mapped + n_profs_mixed

        # Among mapped profiles, retain those with >=1 mapped taxon of this rank.
        mapped_profs_with_rank = prof_summary[
            prof_summary["has_mapped_traits"] & (prof_summary[f"mapped_taxa_{rank}"] > 0)
        ]

        rank_rows.append(
            {
                "taxonomic_rank": rank,
                "display_rank": RANK_DISPLAY.get(rank, rank.capitalize()),
                "total_taxa": total_taxa,
                "zero_trait_taxa": zero_taxa,
                "zero_trait_pct": round(zero_taxa / total_taxa * 100, 2) if total_taxa > 0 else 0.0,
                "sparse_1_2_taxa": sparse_taxa,
                "sparse_1_2_pct": round(sparse_taxa / total_taxa * 100, 2)
                if total_taxa > 0
                else 0.0,
                "rich_gte_3_taxa": rich_taxa,
                "rich_gte_3_pct": round(rich_taxa / total_taxa * 100, 2) if total_taxa > 0 else 0.0,
                "trait_bearing_taxa": trait_taxa,
                "trait_bearing_pct": round(trait_taxa / total_taxa * 100, 2)
                if total_taxa > 0
                else 0.0,
                "total_occurrences": total_occ,
                "zero_trait_occurrences": zero_occ,
                "zero_trait_occurrences_pct": round(zero_occ / total_occ * 100, 2)
                if total_occ > 0
                else 0.0,
                "sparse_occurrences": sparse_occ,
                "rich_occurrences": rich_occ,
                "trait_bearing_occurrences": trait_occ,
                "trait_bearing_occurrences_pct": round(trait_occ / total_occ * 100, 2)
                if total_occ > 0
                else 0.0,
                "profiles_with_rank": n_profs_rank,
                "profiles_with_rank_pct_of_all": round(n_profs_rank / total_profiles * 100, 2),
                "profiles_all_mapped": n_profs_all_mapped,
                "profiles_all_mapped_pct": round(n_profs_all_mapped / n_profs_rank * 100, 2)
                if n_profs_rank > 0
                else 0.0,
                "profiles_mixed_mapped_zero": n_profs_mixed,
                "profiles_mixed_pct": round(n_profs_mixed / n_profs_rank * 100, 2)
                if n_profs_rank > 0
                else 0.0,
                "profiles_zero_only": n_profs_zero_only,
                "profiles_zero_only_pct": round(n_profs_zero_only / n_profs_rank * 100, 2)
                if n_profs_rank > 0
                else 0.0,
                "profiles_any_mapped": n_profs_any_mapped,
                "profiles_any_mapped_pct": round(n_profs_any_mapped / n_profs_rank * 100, 2)
                if n_profs_rank > 0
                else 0.0,
            }
        )

        prof_rank_rows.append(
            {
                "taxonomic_rank": rank,
                "display_rank": RANK_DISPLAY.get(rank, rank.capitalize()),
                "total_profiles_in_cohort": total_profiles,
                "profiles_containing_rank": n_profs_rank,
                "pct_of_all_profiles": round(n_profs_rank / total_profiles * 100, 2),
                "profiles_all_taxa_mapped": n_profs_all_mapped,
                "pct_all_taxa_mapped_within_rank": round(n_profs_all_mapped / n_profs_rank * 100, 2)
                if n_profs_rank > 0
                else 0.0,
                "profiles_mixed_mapped_and_zero": n_profs_mixed,
                "pct_mixed_within_rank": round(n_profs_mixed / n_profs_rank * 100, 2)
                if n_profs_rank > 0
                else 0.0,
                "profiles_zero_traits_only_for_rank": n_profs_zero_only,
                "pct_zero_only_within_rank": round(n_profs_zero_only / n_profs_rank * 100, 2)
                if n_profs_rank > 0
                else 0.0,
                "profiles_with_at_least_1_mapped_taxon": n_profs_any_mapped,
                "pct_any_mapped_within_rank": round(n_profs_any_mapped / n_profs_rank * 100, 2)
                if n_profs_rank > 0
                else 0.0,
                "pct_of_mapped_profiles": round(
                    percent(len(mapped_profs_with_rank), total_mapped_profiles), 2
                ),
                "median_mapped_taxa_when_present": float(
                    mapped_profs_with_rank[f"mapped_taxa_{rank}"].median()
                )
                if len(mapped_profs_with_rank) > 0
                else 0.0,
                "mean_mapped_taxa_when_present": round(
                    float(mapped_profs_with_rank[f"mapped_taxa_{rank}"].mean()), 2
                )
                if len(mapped_profs_with_rank) > 0
                else 0.0,
                "median_profile_mapped_traits_when_rank_mapped": float(
                    mapped_profs_with_rank["mapped_trait_count"].median()
                )
                if len(mapped_profs_with_rank) > 0
                else 0.0,
            }
        )

    rank_summary_df = pd.DataFrame(rank_rows)
    prof_rank_df = pd.DataFrame(prof_rank_rows)

    # Top zero-trait taxa diagnosis
    zero_taxa_df = (
        taxa_master[taxa_master["coverage_category"] == "zero_traits"]
        .sort_values(
            by=["signature_occurrences", "profile_count", "study_count"],
            ascending=[False, False, False],
        )
        .reset_index(drop=True)
    )

    top_zero_records: list[dict[str, Any]] = []
    for rank_idx, row in zero_taxa_df.head(100).iterrows():
        cat, rec = diagnose_zero_trait_taxon(
            str(row["taxon_id"]), str(row["taxon_name"]), str(row["taxonomic_rank"])
        )
        top_zero_records.append(
            {
                "zero_trait_rank": rank_idx + 1,
                "taxon_id": row["taxon_id"],
                "taxon_name": row["taxon_name"],
                "taxonomic_rank": row["taxonomic_rank"],
                "signature_occurrences": row["signature_occurrences"],
                "profile_count": row["profile_count"],
                "study_count": row["study_count"],
                "diagnosis_category": cat,
                "suggested_remediation": rec,
            }
        )
    top_zero_df = pd.DataFrame(top_zero_records)
    if top_zero_df.empty:
        top_zero_df = pd.DataFrame(
            columns=[
                "zero_trait_rank",
                "taxon_id",
                "taxon_name",
                "taxonomic_rank",
                "signature_occurrences",
                "profile_count",
                "study_count",
                "diagnosis_category",
                "suggested_remediation",
            ]
        )
    prof_summary.attrs["population_counts"] = {
        "input_profiles": input_profile_count,
        "selected_profiles": len(prof_summary),
        "excluded_profiles": input_profile_count - len(prof_summary),
        "input_occurrences": input_occurrence_count,
        "selected_occurrences": len(taxa_df),
        "excluded_occurrences": input_occurrence_count - len(taxa_df),
    }

    return taxa_master, rank_summary_df, prof_rank_df, prof_summary, prof_taxa, top_zero_df


def plot_rank_stratified_histogram(
    rank_df: pd.DataFrame, output_path: Path, prof_summary: pd.DataFrame
) -> None:
    """Generate Figure 1: 3-Panel Rank-Stratified Coverage (Unique Taxa Counts, Taxa %, & Profile Distribution by Rank)."""
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22, 7.2), gridspec_kw={"wspace": 0.25})

    ranks = [RANK_DISPLAY[r] for r in RANK_ORDER]
    x = np.arange(len(ranks))
    width = 0.62

    # --- Panel A: Absolute Unique Taxa Counts & Percentages ---
    zeros = rank_df["zero_trait_taxa"].to_numpy()
    sparse = rank_df["sparse_1_2_taxa"].to_numpy()
    rich = rank_df["rich_gte_3_taxa"].to_numpy()
    totals = rank_df["total_taxa"].to_numpy()
    trait_pcts = rank_df["trait_bearing_pct"].to_numpy()

    ax1.bar(
        x,
        zeros,
        width,
        label="Mapped, 0 traits (No KG traits)",
        color=COLOR_ZERO,
        edgecolor="white",
        linewidth=0.8,
    )
    ax1.bar(
        x,
        sparse,
        width,
        bottom=zeros,
        label="Mapped, 1-2 traits (Sparse)",
        color=COLOR_SPARSE,
        edgecolor="white",
        linewidth=0.8,
    )
    ax1.bar(
        x,
        rich,
        width,
        bottom=zeros + sparse,
        label="Mapped, ≥ 3 traits (Rich)",
        color=COLOR_RICH,
        edgecolor="white",
        linewidth=0.8,
    )

    for i, total in enumerate(totals):
        ax1.text(
            x[i],
            total + max(totals) * 0.015,
            f"{total:,}\n({trait_pcts[i]:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=8.0,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
        )

    ax1.set_xticks(x)
    ax1.set_xticklabels(ranks, fontsize=10, rotation=20, ha="right")
    ax1.set_ylabel("Number of Unique Taxa in BugSigDB", fontsize=11, fontweight="bold", labelpad=8)
    ax1.set_title(
        "A. Unique Taxa Counts & Mapped % by Rank",
        fontsize=12,
        fontweight="bold",
        pad=12,
        loc="left",
    )
    ax1.set_ylim(0, max(totals) * 1.18)
    ax1.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax1.set_axisbelow(True)

    ax1.text(
        0.32,
        0.95,
        f"Total catalog: N = {int(rank_df.total_taxa.sum()):,} unique taxa\n"
        f"• Trait-bearing (≥1 trait): {int(rank_df.trait_bearing_taxa.sum()):,} ({percent(rank_df.trait_bearing_taxa.sum(), rank_df.total_taxa.sum()):.1f}%)\n"
        f"• Zero traits: {int(rank_df.zero_trait_taxa.sum()):,} ({percent(rank_df.zero_trait_taxa.sum(), rank_df.total_taxa.sum()):.1f}%)\n(Top label shows total N & % mapped)",
        transform=ax1.transAxes,
        fontsize=8.8,
        verticalalignment="top",
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": COLOR_BG_LIGHT,
            "edgecolor": COLOR_SLATE_LIGHT,
            "alpha": 0.92,
        },
    )

    # --- Panel B: 100% Normalized Proportions of Unique Taxa ---
    zero_pcts = rank_df["zero_trait_pct"].to_numpy()
    sparse_pcts = rank_df["sparse_1_2_pct"].to_numpy()
    rich_pcts = rank_df["rich_gte_3_pct"].to_numpy()

    ax2.bar(x, zero_pcts, width, color=COLOR_ZERO, edgecolor="white", linewidth=0.8)
    ax2.bar(
        x,
        sparse_pcts,
        width,
        bottom=zero_pcts,
        color=COLOR_SPARSE,
        edgecolor="white",
        linewidth=0.8,
    )
    ax2.bar(
        x,
        rich_pcts,
        width,
        bottom=zero_pcts + sparse_pcts,
        color=COLOR_RICH,
        edgecolor="white",
        linewidth=0.8,
    )

    for i in range(len(ranks)):
        if zero_pcts[i] >= 12.0:
            ax2.text(
                x[i],
                zero_pcts[i] / 2.0,
                f"{zero_pcts[i]:.1f}%\n({zeros[i]:,})",
                ha="center",
                va="center",
                color="white",
                fontsize=8.2,
                fontweight="bold",
            )
        if rich_pcts[i] >= 14.0:
            ax2.text(
                x[i],
                zero_pcts[i] + sparse_pcts[i] + rich_pcts[i] / 2.0,
                f"{rich_pcts[i]:.1f}%\n({rich[i]:,})",
                ha="center",
                va="center",
                color="white",
                fontsize=8.2,
                fontweight="bold",
            )

    ax2.axhline(50, color=COLOR_SLATE_MID, linestyle=":", linewidth=1.2, alpha=0.7)
    ax2.text(
        len(ranks) - 0.45,
        51.5,
        "50% parity threshold",
        fontsize=8.5,
        color=COLOR_SLATE_MID,
        ha="right",
    )

    ax2.set_xticks(x)
    ax2.set_xticklabels(ranks, fontsize=10, rotation=20, ha="right")
    ax2.set_ylabel("Percentage of Unique Taxa (%)", fontsize=11, fontweight="bold", labelpad=8)
    ax2.set_title(
        "B. Normalized Taxa Coverage Proportion (%)",
        fontsize=12,
        fontweight="bold",
        pad=12,
        loc="left",
    )
    ax2.set_ylim(0, 108)
    ax2.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax2.set_axisbelow(True)

    # --- Panel C: Distribution of Study Profiles Across Ranks by Trait Mapping ---
    prof_all_mapped = rank_df["profiles_all_mapped"].to_numpy()
    prof_mixed = rank_df["profiles_mixed_mapped_zero"].to_numpy()
    prof_zero_only = rank_df["profiles_zero_only"].to_numpy()
    prof_totals = rank_df["profiles_with_rank"].to_numpy()
    prof_any_mapped_pct = rank_df["profiles_any_mapped_pct"].to_numpy()

    ax3.bar(
        x,
        prof_zero_only,
        width,
        color=COLOR_ZERO,
        edgecolor="white",
        linewidth=0.8,
        label="Profile rank has 0-trait taxa only",
    )
    ax3.bar(
        x,
        prof_mixed,
        width,
        bottom=prof_zero_only,
        color=COLOR_SPARSE,
        edgecolor="white",
        linewidth=0.8,
        label="Profile rank has Mixed (Mapped + 0-trait taxa)",
    )
    ax3.bar(
        x,
        prof_all_mapped,
        width,
        bottom=prof_zero_only + prof_mixed,
        color=COLOR_RICH,
        edgecolor="white",
        linewidth=0.8,
        label="Profile rank has All Taxa Mapped (≥1 trait)",
    )

    for i, p_tot in enumerate(prof_totals):
        ax3.text(
            x[i],
            p_tot + max(prof_totals) * 0.015,
            f"{p_tot:,}\n({prof_any_mapped_pct[i]:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=8.0,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
        )
        # Add inline label for large bars
        if p_tot >= 1000:
            ax3.text(
                x[i],
                prof_zero_only[i] + prof_mixed[i] + prof_all_mapped[i] / 2.0,
                f"{prof_all_mapped[i]:,}\n({prof_all_mapped[i] / p_tot * 100:.0f}%)",
                ha="center",
                va="center",
                color="white",
                fontsize=8,
                fontweight="bold",
            )

    ax3.set_xticks(x)
    ax3.set_xticklabels(ranks, fontsize=10, rotation=20, ha="right")
    ax3.set_ylabel(
        "Number of Study Profiles Containing Rank", fontsize=11, fontweight="bold", labelpad=8
    )
    ax3.set_title(
        "C. Profile Distribution by Taxa Rank & Mapping",
        fontsize=12,
        fontweight="bold",
        pad=12,
        loc="left",
    )
    ax3.set_ylim(0, max(prof_totals) * 1.20)
    ax3.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax3.set_axisbelow(True)

    ax3.text(
        0.37,
        0.96,
        f"Across N = {len(prof_summary):,} Study Profiles:\n"
        f"• {int(prof_summary.has_mapped_traits.sum()):,} ({100 * prof_summary.has_mapped_traits.mean():.1f}%) map to ≥1 trait\n"
        f"• {int((~prof_summary.has_mapped_traits).sum()):,} ({100 * (~prof_summary.has_mapped_traits).mean():.1f}%) have 0 traits\n"
        + "\n".join(
            f"• {row.display_rank}: {row.profiles_with_rank:,} profiles ({row.profiles_any_mapped_pct:.1f}% mapped)"
            for row in rank_df.loc[rank_df.taxonomic_rank.isin(["genus", "species"])].itertuples()
        ),
        transform=ax3.transAxes,
        fontsize=8.5,
        verticalalignment="top",
        bbox={
            "boxstyle": "round,pad=0.38",
            "facecolor": "#ecfdf5",
            "edgecolor": COLOR_ACCENT,
            "alpha": 0.94,
        },
    )

    handles = [
        Patch(
            facecolor=COLOR_RICH,
            edgecolor="white",
            label="Mapped Rich (≥3 traits) / All Rank Taxa Mapped in Profile",
        ),
        Patch(
            facecolor=COLOR_SPARSE,
            edgecolor="white",
            label="Mapped Sparse (1-2 traits) / Mixed (Mapped + 0-Trait Taxa in Profile)",
        ),
        Patch(
            facecolor=COLOR_ZERO,
            edgecolor="white",
            label="0 Traits (No KG Traits) / Only 0-Trait Taxa for Rank in Profile",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.03),
        ncol=3,
        fontsize=10,
        frameon=True,
        facecolor=COLOR_BG_LIGHT,
        edgecolor=COLOR_SLATE_LIGHT,
    )

    save_figure(fig, output_path)
    print(f"  Saved Figure 1: {output_path} and .svg")


def plot_occurrences_vs_unique(
    rank_df: pd.DataFrame, taxa_master: pd.DataFrame, prof_summary: pd.DataFrame, output_path: Path
) -> None:
    """Generate Figure 2: Unique Taxa vs. Literature Signature Occurrences vs. Study Profiles."""
    fig, ax = plt.subplots(figsize=(12, 6.2))

    # 1. Unique Taxa counts
    total_taxa = len(taxa_master)
    zero_taxa = int((taxa_master["coverage_category"] == "zero_traits").sum())
    sparse_taxa = int((taxa_master["coverage_category"] == "sparse_1_2").sum())
    rich_taxa = int((taxa_master["coverage_category"] == "rich_gte_3").sum())

    # 2. Signature Occurrences counts
    total_occ = int(rank_df["total_occurrences"].sum())
    zero_occ = int(rank_df["zero_trait_occurrences"].sum())
    sparse_occ = int(rank_df["sparse_occurrences"].sum())
    rich_occ = int(rank_df["rich_occurrences"].sum())

    # 3. Study Trait Profiles counts
    total_prof = len(prof_summary)
    zero_prof = int((prof_summary["mapped_trait_count"] == 0).sum())
    mixed_prof = int(
        ((prof_summary["mapped_trait_count"] > 0) & (prof_summary["zero_trait_taxa"] > 0)).sum()
    )
    all_mapped_prof = int(
        ((prof_summary["mapped_trait_count"] > 0) & (prof_summary["zero_trait_taxa"] == 0)).sum()
    )

    categories = [
        f"Unique Taxa Catalog\n(N = {total_taxa:,} unique taxa)",
        f"Literature Signature Occurrences\n(N = {total_occ:,} taxon mentions)",
        f"Study Trait Profiles (inc/dec)\n(N = {total_prof:,} profiles)",
    ]
    y_pos = np.arange(len(categories))
    bar_height = 0.52

    zero_pcts = [
        zero_taxa / total_taxa * 100,
        zero_occ / total_occ * 100,
        zero_prof / total_prof * 100,
    ]
    mid_pcts = [
        sparse_taxa / total_taxa * 100,
        sparse_occ / total_occ * 100,
        mixed_prof / total_prof * 100,
    ]
    top_pcts = [
        rich_taxa / total_taxa * 100,
        rich_occ / total_occ * 100,
        all_mapped_prof / total_prof * 100,
    ]

    ax.barh(
        y_pos,
        zero_pcts,
        height=bar_height,
        color=COLOR_ZERO,
        edgecolor="white",
        label="0 Traits (Zero-Trait Taxa / Zero-Trait Profiles)",
    )
    ax.barh(
        y_pos,
        mid_pcts,
        height=bar_height,
        left=zero_pcts,
        color=COLOR_SPARSE,
        edgecolor="white",
        label="Sparse (1-2 Traits) / Mixed Profiles (Mapped + 0-Trait Taxa)",
    )
    ax.barh(
        y_pos,
        top_pcts,
        height=bar_height,
        left=np.array(zero_pcts) + np.array(mid_pcts),
        color=COLOR_RICH,
        edgecolor="white",
        label="Rich (≥3 Traits) / Fully Mapped Profiles (100% Taxa Mapped)",
    )

    # Row 0: Unique Taxa
    ax.text(
        zero_pcts[0] / 2,
        y_pos[0],
        f"0 Traits: {zero_pcts[0]:.1f}%\n({zero_taxa:,} taxa)",
        ha="center",
        va="center",
        color="white",
        fontweight="bold",
        fontsize=9.5,
    )
    ax.text(
        zero_pcts[0] + mid_pcts[0] + top_pcts[0] / 2,
        y_pos[0],
        f"Trait-Bearing Taxa: {mid_pcts[0] + top_pcts[0]:.1f}% ({sparse_taxa + rich_taxa:,} / {total_taxa:,})",
        ha="center",
        va="center",
        color="white",
        fontweight="bold",
        fontsize=10,
    )

    # Row 1: Occurrences
    ax.text(
        zero_pcts[1] / 2,
        y_pos[1],
        f"{zero_pcts[1]:.1f}%\n({zero_occ:,})",
        ha="center",
        va="center",
        color="white",
        fontweight="bold",
        fontsize=9,
    )
    ax.text(
        zero_pcts[1] + mid_pcts[1] + top_pcts[1] / 2,
        y_pos[1],
        f"Trait-Bearing Occurrences: {mid_pcts[1] + top_pcts[1]:.2f}% ({sparse_occ + rich_occ:,} / {total_occ:,})",
        ha="center",
        va="center",
        color="white",
        fontweight="bold",
        fontsize=10,
    )

    # Row 2: Study Profiles
    ax.text(
        max(zero_pcts[2] + 0.8, 4.0),
        y_pos[2] - 0.34,
        f"0-Trait Profiles: {zero_pcts[2]:.1f}% ({zero_prof:,})",
        ha="left",
        va="center",
        color=COLOR_ZERO,
        fontweight="bold",
        fontsize=9,
    )
    ax.text(
        zero_pcts[2] + mid_pcts[2] / 2,
        y_pos[2],
        f"Mixed:\n{mid_pcts[2]:.1f}%\n({mixed_prof:,})",
        ha="center",
        va="center",
        color=COLOR_SLATE_DARK,
        fontweight="bold",
        fontsize=8.5,
    )
    ax.text(
        zero_pcts[2] + mid_pcts[2] + top_pcts[2] / 2,
        y_pos[2],
        f"All Taxa Mapped: {top_pcts[2]:.1f}% ({all_mapped_prof:,})\nTotal Trait-Mapped Profiles: {mid_pcts[2] + top_pcts[2]:.1f}% ({mixed_prof + all_mapped_prof:,} / {total_prof:,})",
        ha="center",
        va="center",
        color="white",
        fontweight="bold",
        fontsize=9.8,
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(categories, fontsize=10.5, fontweight="bold")
    ax.set_xlabel("Percentage of Cohort (%)", fontsize=11, fontweight="bold", labelpad=8)
    ax.set_xlim(0, 105)
    ax.set_ylim(-0.45, 2.55)
    ax.set_title(
        "Coverage Across Granularities: Unique Taxa vs. Literature Occurrences vs. Study Profiles",
        fontsize=12.5,
        fontweight="bold",
        pad=15,
        loc="left",
    )
    ax.grid(axis="x", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax.set_axisbelow(True)

    handles = [
        Patch(
            facecolor=COLOR_RICH,
            edgecolor="white",
            label="Rich (≥3 traits) / All Profile Taxa Mapped",
        ),
        Patch(
            facecolor=COLOR_SPARSE,
            edgecolor="white",
            label="Sparse (1-2 traits) / Mixed Profile (Mapped + 0-Trait Taxa)",
        ),
        Patch(
            facecolor=COLOR_ZERO,
            edgecolor="white",
            label=f"0 Traits (No KG Traits / {zero_prof:,} Zero-Trait Profiles)",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.15),
        ncol=3,
        fontsize=9.2,
        frameon=True,
        facecolor=COLOR_BG_LIGHT,
    )

    save_figure(fig, output_path)
    print(f"  Saved Figure 2: {output_path} and .svg")


def plot_profile_distribution_by_rank(
    prof_rank_df: pd.DataFrame, prof_summary: pd.DataFrame, output_path: Path
) -> None:
    """Generate Figure 3: Distribution of Profiles for Each Type of Taxa (Rank) Based on Trait Mapping."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.5, 6.8), gridspec_kw={"wspace": 0.26})

    ranks = [RANK_DISPLAY[r] for r in RANK_ORDER]
    x = np.arange(len(ranks))
    width = 0.62

    # Panel A: 100% Stacked Bar of Profile Trait Mapping Status for Each Taxonomic Rank (with counts & %)
    all_map_pct = prof_rank_df["pct_all_taxa_mapped_within_rank"].to_numpy()
    mixed_pct = prof_rank_df["pct_mixed_within_rank"].to_numpy()
    zero_pct = prof_rank_df["pct_zero_only_within_rank"].to_numpy()

    all_map_n = prof_rank_df["profiles_all_taxa_mapped"].to_numpy()
    mixed_n = prof_rank_df["profiles_mixed_mapped_and_zero"].to_numpy()
    zero_n = prof_rank_df["profiles_zero_traits_only_for_rank"].to_numpy()
    tot_n = prof_rank_df["profiles_containing_rank"].to_numpy()
    any_map_pct = prof_rank_df["pct_any_mapped_within_rank"].to_numpy()

    ax1.bar(
        x,
        zero_pct,
        width,
        color=COLOR_ZERO,
        edgecolor="white",
        linewidth=0.8,
        label="Zero-Only for Rank (0 mapped taxa of this rank)",
    )
    ax1.bar(
        x,
        mixed_pct,
        width,
        bottom=zero_pct,
        color=COLOR_SPARSE,
        edgecolor="white",
        linewidth=0.8,
        label="Mixed for Rank (≥1 mapped + ≥1 zero-trait taxon)",
    )
    ax1.bar(
        x,
        all_map_pct,
        width,
        bottom=zero_pct + mixed_pct,
        color=COLOR_RICH,
        edgecolor="white",
        linewidth=0.8,
        label="All Mapped for Rank (100% of rank taxa have ≥1 trait)",
    )

    for i in range(len(ranks)):
        # Concise top header above bar: total profiles containing this rank & % mapped
        ax1.text(
            x[i],
            101.2,
            f"N={tot_n[i]:,}\n({any_map_pct[i]:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=8.0,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
        )
        if zero_pct[i] >= 8.0:
            ax1.text(
                x[i],
                zero_pct[i] / 2.0,
                f"{zero_pct[i]:.1f}%\n({zero_n[i]:,})",
                ha="center",
                va="center",
                color="white",
                fontsize=8.0,
                fontweight="bold",
            )
        if mixed_pct[i] >= 8.0:
            ax1.text(
                x[i],
                zero_pct[i] + mixed_pct[i] / 2.0,
                f"{mixed_pct[i]:.1f}%\n({mixed_n[i]:,})",
                ha="center",
                va="center",
                color=COLOR_SLATE_DARK,
                fontsize=8.0,
                fontweight="bold",
            )
        if all_map_pct[i] >= 12.0:
            ax1.text(
                x[i],
                zero_pct[i] + mixed_pct[i] + all_map_pct[i] / 2.0,
                f"{all_map_pct[i]:.1f}%\n({all_map_n[i]:,})",
                ha="center",
                va="center",
                color="white",
                fontsize=8.2,
                fontweight="bold",
            )

    ax1.set_xticks(x)
    ax1.set_xticklabels(ranks, fontsize=10.5, rotation=20, ha="right")
    ax1.set_ylabel(
        "Percentage of Profiles Containing Rank (%)", fontsize=11, fontweight="bold", labelpad=8
    )
    ax1.set_title(
        "A. Profile Mapping Composition Within Each Taxa Rank",
        fontsize=12,
        fontweight="bold",
        pad=16,
        loc="left",
    )
    ax1.set_ylim(0, 115)
    ax1.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax1.set_axisbelow(True)

    # Panel B: Trait counts of mapped profiles containing each mapped rank.
    mapped_only = prof_summary[prof_summary["has_mapped_traits"]]
    box_data = []
    box_labels = []
    rank_medians = []
    for r in RANK_ORDER:
        sub_vals = mapped_only[mapped_only[f"mapped_taxa_{r}"] > 0]["mapped_trait_count"].to_numpy()
        if len(sub_vals) > 0:
            box_data.append(sub_vals)
            n_p = len(sub_vals)
            pct_p = n_p / len(mapped_only) * 100
            box_labels.append(f"{RANK_DISPLAY[r]}\n(n={n_p:,}, {pct_p:.1f}%)")
            rank_medians.append(float(np.median(sub_vals)))

    if not box_data:
        ax2.text(
            0.5,
            0.5,
            "No mapped profiles; richness is not estimable",
            ha="center",
            transform=ax2.transAxes,
        )
        ax2.set_axis_off()
        save_figure(fig, output_path)
        return
    bp = ax2.boxplot(
        box_data,
        tick_labels=box_labels,
        patch_artist=True,
        showmeans=True,
        meanprops={
            "marker": "o",
            "markeredgecolor": COLOR_SLATE_DARK,
            "markerfacecolor": "white",
            "markersize": 5.5,
        },
        medianprops={"color": COLOR_SLATE_DARK, "linewidth": 1.8},
        whiskerprops={"color": COLOR_SLATE_MID, "linewidth": 1.2},
        capprops={"color": COLOR_SLATE_MID, "linewidth": 1.2},
        flierprops={
            "marker": ".",
            "markerfacecolor": COLOR_SLATE_LIGHT,
            "markersize": 3,
            "alpha": 0.35,
        },
    )
    present_ranks = [r for r in RANK_ORDER if (mapped_only[f"mapped_taxa_{r}"] > 0).any()]
    for patch, r in zip(bp["boxes"], present_ranks, strict=True):
        patch.set_facecolor(RANK_COLORS[r])
        patch.set_alpha(0.78)
        patch.set_edgecolor(COLOR_SLATE_DARK)

    for i, med in enumerate(rank_medians):
        ax2.text(
            i + 1,
            med * 1.08,
            f"Med: {round(med)}",
            ha="center",
            va="bottom",
            fontsize=8.2,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
            bbox={
                "boxstyle": "round,pad=0.18",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.8,
            },
        )

    ax2.set_yscale("log")
    ax2.set_xticks(np.arange(1, len(box_labels) + 1))
    ax2.set_xticklabels(box_labels, fontsize=9.2, rotation=20, ha="right")
    ax2.set_ylabel(
        "Mapped Traits per Profile (Log10 Scale)", fontsize=11, fontweight="bold", labelpad=8
    )
    ax2.set_title(
        f"B. Mapped Trait Richness of Profiles by Taxa Rank (N = {int(prof_summary.has_mapped_traits.sum()):,})",
        fontsize=12,
        fontweight="bold",
        pad=16,
        loc="left",
    )
    ax2.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax2.set_axisbelow(True)

    handles = [
        Patch(
            facecolor=COLOR_RICH,
            edgecolor="white",
            label="All Rank Taxa Mapped in Profile (100% have ≥1 trait)",
        ),
        Patch(
            facecolor=COLOR_SPARSE,
            edgecolor="white",
            label="Mixed Rank Taxa in Profile (Both Mapped & 0-Trait Taxa)",
        ),
        Patch(
            facecolor=COLOR_ZERO,
            edgecolor="white",
            label="Zero-Only for Rank in Profile (0 Mapped Taxa of that Rank)",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.03),
        ncol=3,
        fontsize=9.8,
        frameon=True,
        facecolor=COLOR_BG_LIGHT,
        edgecolor=COLOR_SLATE_LIGHT,
    )

    save_figure(fig, output_path)
    print(f"  Saved Figure 3 (Profile Distribution by Rank): {output_path} and .svg")


def plot_histogram1_mapped_profiles(prof_summary: pd.DataFrame, output_path: Path) -> None:
    """Generate Histogram 1: mapped traits and taxa per mapped profile."""
    mapped_df = prof_summary[prof_summary["has_mapped_traits"]].copy()
    n_mapped = len(mapped_df)
    n_total = len(prof_summary)
    n_excluded = n_total - n_mapped

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.8, 7.0), gridspec_kw={"wspace": 0.25})

    # --- Panel A: Frequency of mapped traits across mapped profiles ---
    trait_bins = [
        ("1-19", (mapped_df["mapped_trait_count"] >= 1) & (mapped_df["mapped_trait_count"] <= 19)),
        (
            "20-49",
            (mapped_df["mapped_trait_count"] >= 20) & (mapped_df["mapped_trait_count"] <= 49),
        ),
        (
            "50-99",
            (mapped_df["mapped_trait_count"] >= 50) & (mapped_df["mapped_trait_count"] <= 99),
        ),
        (
            "100-199",
            (mapped_df["mapped_trait_count"] >= 100) & (mapped_df["mapped_trait_count"] <= 199),
        ),
        (
            "200-399",
            (mapped_df["mapped_trait_count"] >= 200) & (mapped_df["mapped_trait_count"] <= 399),
        ),
        (
            "400-799",
            (mapped_df["mapped_trait_count"] >= 400) & (mapped_df["mapped_trait_count"] <= 799),
        ),
        (
            "800-1,199",
            (mapped_df["mapped_trait_count"] >= 800) & (mapped_df["mapped_trait_count"] <= 1199),
        ),
        ("≥1,200", mapped_df["mapped_trait_count"] >= 1200),
    ]
    bin_labels_a = [b[0] for b in trait_bins]
    bin_counts_a = [int(b[1].sum()) for b in trait_bins]
    bin_pcts_a = [percent(c, n_mapped) for c in bin_counts_a]

    x_a = np.arange(len(bin_labels_a))
    bars_a = ax1.bar(
        x_a, bin_counts_a, width=0.66, color=COLOR_RICH, edgecolor="white", linewidth=0.9
    )
    bars_a[0].set_facecolor(COLOR_SPARSE)

    for i, (cnt, pct) in enumerate(zip(bin_counts_a, bin_pcts_a, strict=False)):
        ax1.text(
            x_a[i],
            cnt + max(bin_counts_a) * 0.015,
            f"{cnt:,}\n({pct:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=8.6,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
        )

    med_traits, q25_traits, q75_traits, mean_traits = distribution_stats(
        mapped_df["mapped_trait_count"]
    )

    median_trait_bin = int(
        np.searchsorted([20, 50, 100, 200, 400, 800, 1200], med_traits, side="right")
    )
    if n_mapped:
        ax1.axvline(median_trait_bin, color=COLOR_ACCENT, linestyle="--", linewidth=1.8)
    ax1.annotate(
        f"Median: {med_traits:g} traits\n(IQR: {q25_traits:g}-{q75_traits:g})\nMean: {mean_traits:.1f} traits"
        if n_mapped
        else "No mapped profiles",
        xy=(median_trait_bin, max(bin_counts_a) * 0.82),
        xytext=(4.35, max(bin_counts_a) * 0.88),
        arrowprops={"arrowstyle": "->", "color": COLOR_ACCENT, "lw": 1.6},
        fontsize=9.0,
        fontweight="bold",
        color=COLOR_SLATE_DARK,
        bbox={
            "boxstyle": "round,pad=0.38",
            "facecolor": "#ecfdf5",
            "edgecolor": COLOR_ACCENT,
            "alpha": 0.95,
        },
    )

    ax1.set_xticks(x_a)
    ax1.set_xticklabels(bin_labels_a, fontsize=10, rotation=20, ha="right")
    ax1.set_xlabel(
        "Number of Unique Mapped KG-Microbe Traits per Profile",
        fontsize=11,
        fontweight="bold",
        labelpad=8,
    )
    ax1.set_ylabel(
        f"Number of Mapped Study Profiles (N = {n_mapped:,})",
        fontsize=11,
        fontweight="bold",
        labelpad=8,
    )
    ax1.set_title(
        "A. Frequency of Mapped Traits Across Profiles (≥1 Trait)",
        fontsize=12,
        fontweight="bold",
        pad=12,
        loc="left",
    )
    ax1.set_ylim(0, max(1, max(bin_counts_a)) * 1.25)
    ax1.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax1.set_axisbelow(True)

    # Zero-trait exclusion confirmation box in Panel A (placed cleanly over left bins 0-2)
    ax1.text(
        0.02,
        0.96,
        f"This histogram shows mapped profiles:\n• Included (≥1 trait): {n_mapped:,} / {n_total:,} ({percent(n_mapped, n_total):.1f}%)\n• Not shown (0 traits): {n_excluded:,} / {n_total:,} ({percent(n_excluded, n_total):.1f}%)",
        transform=ax1.transAxes,
        fontsize=8.6,
        verticalalignment="top",
        bbox={
            "boxstyle": "round,pad=0.38",
            "facecolor": "#fff7ed",
            "edgecolor": COLOR_ZERO,
            "alpha": 0.94,
        },
    )

    # --- Panel B: Mapped taxa per mapped profile, stacked by rank ---
    taxa_bins = [
        ("1", mapped_df["mapped_trait_taxa"] == 1),
        ("2", mapped_df["mapped_trait_taxa"] == 2),
        ("3-4", (mapped_df["mapped_trait_taxa"] >= 3) & (mapped_df["mapped_trait_taxa"] <= 4)),
        ("5-7", (mapped_df["mapped_trait_taxa"] >= 5) & (mapped_df["mapped_trait_taxa"] <= 7)),
        ("8-11", (mapped_df["mapped_trait_taxa"] >= 8) & (mapped_df["mapped_trait_taxa"] <= 11)),
        ("12-19", (mapped_df["mapped_trait_taxa"] >= 12) & (mapped_df["mapped_trait_taxa"] <= 19)),
        ("20-34", (mapped_df["mapped_trait_taxa"] >= 20) & (mapped_df["mapped_trait_taxa"] <= 34)),
        ("≥35", mapped_df["mapped_trait_taxa"] >= 35),
    ]
    bin_labels_b = [b[0] for b in taxa_bins]
    bin_counts_b = [int(b[1].sum()) for b in taxa_bins]
    bin_pcts_b = [percent(c, n_mapped) for c in bin_counts_b]

    # Compute fractional rank contribution to profiles in each bin so stack heights sum exactly to profile count!
    rank_stack = {r: np.zeros(len(taxa_bins)) for r in RANK_ORDER}
    for b_idx, (_, mask) in enumerate(taxa_bins):
        sub_p = mapped_df[mask]
        if len(sub_p) == 0:
            continue
        denom = sub_p["mapped_trait_taxa"].to_numpy(dtype=float)
        for r in RANK_ORDER:
            r_frac = sub_p[f"mapped_taxa_{r}"].to_numpy(dtype=float) / denom
            rank_stack[r][b_idx] = float(r_frac.sum())

    x_b = np.arange(len(bin_labels_b))
    bottom_b = np.zeros(len(bin_labels_b))
    for r in RANK_ORDER:
        ax2.bar(
            x_b,
            rank_stack[r],
            width=0.66,
            bottom=bottom_b,
            color=RANK_COLORS[r],
            edgecolor="white",
            linewidth=0.7,
            label=RANK_DISPLAY[r],
        )
        bottom_b += rank_stack[r]

    for i, (cnt, pct) in enumerate(zip(bin_counts_b, bin_pcts_b, strict=False)):
        ax2.text(
            x_b[i],
            cnt + max(bin_counts_b) * 0.015,
            f"{cnt:,}\n({pct:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=8.6,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
        )

    med_taxa, q25_taxa, q75_taxa, mean_taxa = distribution_stats(mapped_df["mapped_trait_taxa"])

    median_taxon_bin = int(np.searchsorted([2, 3, 5, 8, 12, 20, 35], med_taxa, side="right"))
    if n_mapped:
        ax2.axvline(median_taxon_bin, color=COLOR_SLATE_DARK, linestyle="--", linewidth=1.8)
    ax2.annotate(
        f"Mapped Taxa / Profile:\n• Median: {med_taxa:g} taxa (IQR: {q25_taxa:g}-{q75_taxa:g})\n• Mean: {mean_taxa:.1f} mapped taxa\n• Stacked by rank composition"
        if n_mapped
        else "No mapped profiles",
        xy=(median_taxon_bin, max(bin_counts_b) * 0.82),
        xytext=(3.85, max(bin_counts_b) * 0.78),
        arrowprops={"arrowstyle": "->", "color": COLOR_SLATE_DARK, "lw": 1.5},
        fontsize=8.2,
        fontweight="bold",
        color=COLOR_SLATE_DARK,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": COLOR_BG_LIGHT,
            "edgecolor": COLOR_SLATE_LIGHT,
            "alpha": 0.95,
        },
    )

    ax2.set_xticks(x_b)
    ax2.set_xticklabels(bin_labels_b, fontsize=10)
    ax2.set_xlabel(
        "Number of Trait-Mapped Taxa per Profile", fontsize=11, fontweight="bold", labelpad=8
    )
    ax2.set_ylabel(
        f"Number of Mapped Study Profiles (N = {n_mapped:,})",
        fontsize=11,
        fontweight="bold",
        labelpad=8,
    )
    ax2.set_title(
        "B. Mapped Taxa per Profile (Stacked by Taxonomic Rank)",
        fontsize=12,
        fontweight="bold",
        pad=12,
        loc="left",
    )
    ax2.set_ylim(0, max(1, max(bin_counts_b)) * 1.34)
    ax2.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax2.set_axisbelow(True)

    handles_b = [
        Patch(facecolor=RANK_COLORS[r], edgecolor="white", label=RANK_DISPLAY[r])
        for r in RANK_ORDER
    ]
    ax2.legend(
        handles=handles_b,
        title="Taxonomic Rank of Mapped Taxa",
        title_fontsize=8.5,
        loc="upper right",
        ncol=2,
        fontsize=8.0,
        frameon=True,
        facecolor="white",
        edgecolor=COLOR_SLATE_LIGHT,
    )

    fig.suptitle(
        f"Histogram 1: Mapped Traits & Mapped Taxa Across Profiles (≥1 Trait, N = {n_mapped:,})",
        fontsize=13.5,
        fontweight="bold",
        y=1.02,
    )

    save_figure(fig, output_path)
    print(f"  Saved Histogram 1: {output_path} and .svg")


def plot_histogram2_taxa_traits_and_zero_frequency(
    taxa_master: pd.DataFrame,
    prof_summary: pd.DataFrame,
    output_path: Path,
    legacy_path: Path | None = None,
) -> None:
    """Generate Histogram 2: Number of Traits from Taxa Profiles & Frequency of Zero-Trait Taxa Across Profiles."""
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22.5, 7.2), gridspec_kw={"wspace": 0.25})

    n_all_taxa = len(taxa_master)
    taxa_master["trait_count"].to_numpy()
    zero_taxa_df = taxa_master[taxa_master["trait_count"] == 0].copy()
    n_zero_taxa = len(zero_taxa_df)

    # --- Panel A: Number of Traits per Taxon (Stacked by Rank, 0-Bin Highlighted, Counts & %) ---
    deg_bins = [
        ("0\n(No traits)", taxa_master["trait_count"] == 0),
        ("1-4", (taxa_master["trait_count"] >= 1) & (taxa_master["trait_count"] <= 4)),
        ("5-9", (taxa_master["trait_count"] >= 5) & (taxa_master["trait_count"] <= 9)),
        ("10-19", (taxa_master["trait_count"] >= 10) & (taxa_master["trait_count"] <= 19)),
        ("20-49", (taxa_master["trait_count"] >= 20) & (taxa_master["trait_count"] <= 49)),
        ("50-99", (taxa_master["trait_count"] >= 50) & (taxa_master["trait_count"] <= 99)),
        ("100-199", (taxa_master["trait_count"] >= 100) & (taxa_master["trait_count"] <= 199)),
        ("200-499", (taxa_master["trait_count"] >= 200) & (taxa_master["trait_count"] <= 499)),
        ("≥500", taxa_master["trait_count"] >= 500),
    ]
    labels_a = [b[0] for b in deg_bins]
    totals_a = [int(b[1].sum()) for b in deg_bins]
    pcts_a = [c / n_all_taxa * 100 for c in totals_a]

    x_a = np.arange(len(labels_a))
    bottom_a = np.zeros(len(labels_a))
    for r in RANK_ORDER:
        r_counts = np.array(
            [int((mask & (taxa_master["taxonomic_rank"] == r)).sum()) for _, mask in deg_bins]
        )
        ax1.bar(
            x_a,
            r_counts,
            width=0.66,
            bottom=bottom_a,
            color=RANK_COLORS[r],
            edgecolor="white",
            linewidth=0.6,
            label=RANK_DISPLAY[r],
        )
        bottom_a += r_counts

    # Outline the 0-trait bar in bold Vermilion to visually distinguish it
    ax1.bar(x_a[0], totals_a[0], width=0.66, facecolor="none", edgecolor=COLOR_ZERO, linewidth=2.5)

    for i, (cnt, pct) in enumerate(zip(totals_a, pcts_a, strict=False)):
        ax1.text(
            x_a[i],
            cnt + max(totals_a) * 0.015,
            f"{cnt:,}\n({pct:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=8.0,
            fontweight="bold",
            color=COLOR_ZERO if i == 0 else COLOR_SLATE_DARK,
        )

    ax1.annotate(
        f"Zero-Trait Taxa:\n{n_zero_taxa:,} taxa ({percent(n_zero_taxa, n_all_taxa):.1f}%)\n"
        + "\n".join(
            f"• {RANK_DISPLAY[r]}: {int(zero_taxa_df.taxonomic_rank.eq(r).sum()):,} ({percent(zero_taxa_df.taxonomic_rank.eq(r).sum(), n_zero_taxa):.1f}%)"
            for r in ("species", "genus")
        )
        + f"\n• Other ranks: {int((~zero_taxa_df.taxonomic_rank.isin(['species', 'genus'])).sum()):,}",
        xy=(x_a[0], totals_a[0] * 0.85),
        xytext=(x_a[0] + 1.35, totals_a[0] * 0.76),
        arrowprops={"arrowstyle": "->", "color": COLOR_ZERO, "lw": 1.6},
        fontsize=8.8,
        fontweight="bold",
        color=COLOR_ZERO,
        bbox={
            "boxstyle": "round,pad=0.38",
            "facecolor": "#fff7ed",
            "edgecolor": COLOR_ZERO,
            "alpha": 0.95,
        },
    )

    ax1.set_xticks(x_a)
    ax1.set_xticklabels(labels_a, fontsize=9, rotation=20, ha="right")
    ax1.set_xlabel(
        "Number of Annotated KG-Microbe Traits per Taxon",
        fontsize=11,
        fontweight="bold",
        labelpad=8,
    )
    ax1.set_ylabel(
        f"Number of Unique Taxa (N = {n_all_taxa:,})", fontsize=11, fontweight="bold", labelpad=8
    )
    ax1.set_title(
        "A. Number of Traits per Taxon (Ranked & 0-Bin Highlighted)",
        fontsize=12,
        fontweight="bold",
        pad=12,
        loc="left",
    )
    ax1.set_ylim(0, max(totals_a) * 1.18)
    ax1.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax1.set_axisbelow(True)

    # --- Panel B: Frequency of taxa with no traits across all selected profiles ---
    n_all_prof = len(prof_summary)
    z_prof_bins = [
        ("0", prof_summary["zero_trait_taxa"] == 0),
        ("1", prof_summary["zero_trait_taxa"] == 1),
        ("2", prof_summary["zero_trait_taxa"] == 2),
        ("3-4", (prof_summary["zero_trait_taxa"] >= 3) & (prof_summary["zero_trait_taxa"] <= 4)),
        ("5-9", (prof_summary["zero_trait_taxa"] >= 5) & (prof_summary["zero_trait_taxa"] <= 9)),
        ("≥10", prof_summary["zero_trait_taxa"] >= 10),
    ]
    labels_b = [b[0] for b in z_prof_bins]
    totals_b = [int(b[1].sum()) for b in z_prof_bins]
    pcts_b = [c / n_all_prof * 100 for c in totals_b]

    # Split each bar into:
    # 1) Profiles that map to >=1 trait.
    # 2) Profiles that have zero mapped traits; these remain in this audit.
    retained_counts = np.array(
        [int((mask & prof_summary["has_mapped_traits"]).sum()) for _, mask in z_prof_bins]
    )
    excluded_counts = np.array(
        [int((mask & (~prof_summary["has_mapped_traits"])).sum()) for _, mask in z_prof_bins]
    )

    x_b = np.arange(len(labels_b))
    # For bin 0 ("0 zero-trait taxa"), color in Rich Blue; for bins >=1, retained is Sky Blue (Mixed) and excluded is Vermilion
    bar_colors_ret = [COLOR_RICH] + [COLOR_SPARSE] * (len(labels_b) - 1)
    ax2.bar(
        x_b,
        retained_counts,
        width=0.64,
        color=bar_colors_ret,
        edgecolor="white",
        linewidth=0.8,
        label=f"Mapped Profile (≥1 Trait; N = {int(prof_summary.has_mapped_traits.sum()):,})",
    )
    ax2.bar(
        x_b,
        excluded_counts,
        width=0.64,
        bottom=retained_counts,
        color=COLOR_ZERO,
        edgecolor="white",
        linewidth=0.8,
        label=f"Zero-Trait Profile (N = {int((~prof_summary.has_mapped_traits).sum()):,})",
    )

    for i, (cnt, pct) in enumerate(zip(totals_b, pcts_b, strict=False)):
        extra = f"\n(0 traits: {excluded_counts[i]})" if excluded_counts[i] > 0 else ""
        ax2.text(
            x_b[i],
            cnt + max(totals_b) * 0.015,
            f"{cnt:,}\n({pct:.1f}%){extra}",
            ha="center",
            va="bottom",
            fontsize=8.1,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
        )

    ax2.text(
        0.22,
        0.95,
        f"Zero-Trait Taxa Across {n_all_prof:,} Profiles:\n"
        f"• {totals_b[0]:,} ({pcts_b[0]:.1f}%) have 0 zero-trait taxa\n"
        f"• {sum(totals_b[1:]):,} ({sum(pcts_b[1:]):.1f}%) have ≥1 zero-trait taxon\n"
        f"• {int(retained_counts[1:].sum()):,} are mixed (also map to traits)\n"
        f"• {int(excluded_counts.sum()):,} have 0 mapped traits",
        transform=ax2.transAxes,
        fontsize=8.6,
        verticalalignment="top",
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": COLOR_BG_LIGHT,
            "edgecolor": COLOR_SLATE_LIGHT,
            "alpha": 0.95,
        },
    )

    ax2.set_xticks(x_b)
    ax2.set_xticklabels(labels_b, fontsize=10)
    ax2.set_xlabel(
        "Number of Zero-Trait Taxa per Study Profile", fontsize=11, fontweight="bold", labelpad=8
    )
    ax2.set_ylabel(
        f"Number of Study Profiles (N = {n_all_prof:,})", fontsize=11, fontweight="bold", labelpad=8
    )
    ax2.set_title(
        f"B. Zero-Trait Taxa Frequency per Profile (All {n_all_prof:,} Profiles)",
        fontsize=12,
        fontweight="bold",
        pad=12,
        loc="left",
    )
    ax2.set_ylim(0, max(totals_b) * 1.18)
    ax2.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax2.set_axisbelow(True)

    # --- Panel C: Profile recurrence of zero-trait taxa, stacked by rank ---
    occ_bins = [
        ("1 profile\n(Singleton)", zero_taxa_df["profile_count"] == 1),
        ("2 profiles", zero_taxa_df["profile_count"] == 2),
        (
            "3-4 profiles",
            (zero_taxa_df["profile_count"] >= 3) & (zero_taxa_df["profile_count"] <= 4),
        ),
        (
            "5-9 profiles",
            (zero_taxa_df["profile_count"] >= 5) & (zero_taxa_df["profile_count"] <= 9),
        ),
        (
            "10-19 profiles",
            (zero_taxa_df["profile_count"] >= 10) & (zero_taxa_df["profile_count"] <= 19),
        ),
        ("≥20 profiles", zero_taxa_df["profile_count"] >= 20),
    ]
    labels_c = [b[0] for b in occ_bins]
    totals_c = [int(b[1].sum()) for b in occ_bins]
    pcts_c = [percent(c, n_zero_taxa) for c in totals_c]

    x_c = np.arange(len(labels_c))
    bottom_c = np.zeros(len(labels_c))
    for r in RANK_ORDER:
        r_counts = np.array(
            [int((mask & (zero_taxa_df["taxonomic_rank"] == r)).sum()) for _, mask in occ_bins]
        )
        ax3.bar(
            x_c,
            r_counts,
            width=0.64,
            bottom=bottom_c,
            color=RANK_COLORS[r],
            edgecolor="white",
            linewidth=0.6,
            label=RANK_DISPLAY[r],
        )
        bottom_c += r_counts

    for i, (cnt, pct) in enumerate(zip(totals_c, pcts_c, strict=False)):
        ax3.text(
            x_c[i],
            cnt + max(totals_c) * 0.015,
            f"{cnt:,}\n({pct:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=8.3,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
        )

    ax3.text(
        0.24,
        0.95,
        f"Recurrence of Zero-Trait Taxa (N = {n_zero_taxa:,}):\n• {totals_c[0]:,} ({pcts_c[0]:.1f}%) appear in 1 profile\n• {totals_c[0] + totals_c[1]:,} ({pcts_c[0] + pcts_c[1]:.1f}%) appear in ≤2 profiles\n• Most recurring zero-trait taxa:\n"
        + "\n".join(
            f"  - {row.taxon_name} ({row.profile_count:,} profiles)"
            for row in zero_taxa_df.sort_values(
                ["profile_count", "taxon_id"], ascending=[False, True]
            )
            .head(2)
            .itertuples()
        ),
        transform=ax3.transAxes,
        fontsize=8.6,
        verticalalignment="top",
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": "#fff7ed",
            "edgecolor": COLOR_ZERO,
            "alpha": 0.95,
        },
    )

    ax3.set_xticks(x_c)
    ax3.set_xticklabels(labels_c, fontsize=9.5, rotation=15, ha="right")
    ax3.set_xlabel(
        "Number of Study Profiles in Which Zero-Trait Taxon Appears",
        fontsize=11,
        fontweight="bold",
        labelpad=8,
    )
    ax3.set_ylabel(
        f"Number of Zero-Trait Taxa (N = {n_zero_taxa:,})",
        fontsize=11,
        fontweight="bold",
        labelpad=8,
    )
    ax3.set_title(
        "C. Profile Recurrence of Zero-Trait Taxa (Stacked by Rank)",
        fontsize=12,
        fontweight="bold",
        pad=12,
        loc="left",
    )
    ax3.set_ylim(0, max(1, max(totals_c)) * 1.18)
    ax3.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax3.set_axisbelow(True)

    handles_rank = [
        Patch(facecolor=RANK_COLORS[r], edgecolor="white", label=RANK_DISPLAY[r])
        for r in RANK_ORDER
    ]
    handles_prof = [
        Patch(facecolor=COLOR_SPARSE, edgecolor="white", label="Mixed Profile (≥1 Trait)"),
        Patch(facecolor=COLOR_ZERO, edgecolor="white", label="0-Trait Profile"),
    ]
    fig.legend(
        handles=handles_rank + handles_prof,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.03),
        ncol=9,
        fontsize=8.8,
        frameon=True,
        facecolor=COLOR_BG_LIGHT,
        edgecolor=COLOR_SLATE_LIGHT,
    )

    fig.suptitle(
        "Histogram 2: Number of Traits from Taxa Profiles & Frequency of Zero-Trait Taxa Across Profiles",
        fontsize=13.5,
        fontweight="bold",
        y=1.08,
    )

    save_figure(fig, output_path)
    print(f"  Saved Histogram 2: {output_path} and .svg")

    if legacy_path is not None:
        # Also update the legacy filename so existing references see the enhanced diagram
        plot_legacy_trait_count_distribution(taxa_master, legacy_path)


def plot_legacy_trait_count_distribution(taxa_master: pd.DataFrame, output_path: Path) -> None:
    """Update legacy Figure `taxa_trait_count_distribution.png` with both counts & percentages."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.2), gridspec_kw={"wspace": 0.25})

    trait_counts = taxa_master["trait_count"].to_numpy()
    n_total = len(taxa_master)
    zero_count = int((trait_counts == 0).sum())

    bin_labels = ["0", "1-4", "5-9", "10-19", "20-49", "50-99", "100-199", "200-499", "≥500"]
    counts_per_bin = [
        int(zero_count),
        int(((trait_counts >= 1) & (trait_counts < 5)).sum()),
        int(((trait_counts >= 5) & (trait_counts < 10)).sum()),
        int(((trait_counts >= 10) & (trait_counts < 20)).sum()),
        int(((trait_counts >= 20) & (trait_counts < 50)).sum()),
        int(((trait_counts >= 50) & (trait_counts < 100)).sum()),
        int(((trait_counts >= 100) & (trait_counts < 200)).sum()),
        int(((trait_counts >= 200) & (trait_counts < 500)).sum()),
        int((trait_counts >= 500).sum()),
    ]
    pcts_per_bin = [c / n_total * 100 for c in counts_per_bin]

    x_b = np.arange(len(bin_labels))
    bar_colors = [COLOR_ZERO] + [
        COLOR_SPARSE if i == 1 else COLOR_RICH for i in range(1, len(bin_labels))
    ]

    ax1.bar(x_b, counts_per_bin, color=bar_colors, edgecolor="white", linewidth=0.8)
    for i, (c, p) in enumerate(zip(counts_per_bin, pcts_per_bin, strict=False)):
        ax1.text(
            x_b[i],
            c + max(counts_per_bin) * 0.015,
            f"{c:,}\n({p:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=8.0,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
        )

    ax1.set_xticks(x_b)
    ax1.set_xticklabels(bin_labels, fontsize=9.5, rotation=25)
    ax1.set_ylabel(
        f"Number of Unique Taxa (N = {n_total:,})", fontsize=11, fontweight="bold", labelpad=8
    )
    ax1.set_xlabel(
        "Number of Annotated KG-Microbe Traits per Taxon",
        fontsize=11,
        fontweight="bold",
        labelpad=8,
    )
    ax1.set_title(
        "A. Trait Degree Distribution Across All Taxa (Counts & %)",
        fontsize=12,
        fontweight="bold",
        pad=12,
        loc="left",
    )
    ax1.set_ylim(0, max(counts_per_bin) * 1.16)
    ax1.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax1.set_axisbelow(True)

    ax1.annotate(
        f"Zero-Trait Taxa:\n{zero_count:,} taxa ({percent(zero_count, n_total):.1f}%)\nhave no recorded KG traits",
        xy=(x_b[0], zero_count),
        xytext=(x_b[0] + 1.2, zero_count * 0.85),
        arrowprops={"arrowstyle": "->", "color": COLOR_ZERO, "lw": 1.5},
        fontsize=9.5,
        fontweight="bold",
        color=COLOR_ZERO,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "#fff7ed",
            "edgecolor": COLOR_ZERO,
            "alpha": 0.9,
        },
    )

    rank_pos_data = []
    plot_ranks = []
    for rank in RANK_ORDER:
        sub = taxa_master[
            (taxa_master["taxonomic_rank"] == rank) & (taxa_master["trait_count"] > 0)
        ]["trait_count"].to_numpy()
        if len(sub) > 0:
            rank_pos_data.append(sub)
            plot_ranks.append(f"{RANK_DISPLAY[rank]}\n(n={len(sub):,})")

    if not rank_pos_data:
        ax2.text(
            0.5,
            0.5,
            "No trait-bearing taxa; richness is not estimable",
            ha="center",
            transform=ax2.transAxes,
        )
        ax2.set_axis_off()
        save_figure(fig, output_path)
        return
    ax2.boxplot(
        rank_pos_data,
        tick_labels=plot_ranks,
        patch_artist=True,
        showmeans=True,
        meanprops={
            "marker": "o",
            "markeredgecolor": COLOR_SLATE_DARK,
            "markerfacecolor": "white",
            "markersize": 5,
        },
        medianprops={"color": COLOR_SLATE_DARK, "linewidth": 1.5},
        boxprops={"facecolor": "#bae6fd", "color": COLOR_RICH, "linewidth": 1.2},
        whiskerprops={"color": COLOR_SLATE_MID, "linewidth": 1.2},
        capprops={"color": COLOR_SLATE_MID, "linewidth": 1.2},
        flierprops={
            "marker": ".",
            "markerfacecolor": COLOR_SLATE_LIGHT,
            "markersize": 3,
            "alpha": 0.4,
        },
    )

    ax2.set_yscale("log")
    ax2.set_xticklabels(plot_ranks, fontsize=9, rotation=15, ha="right")
    ax2.set_ylabel(
        "Annotated Traits per Taxon (Log10 Scale)", fontsize=11, fontweight="bold", labelpad=8
    )
    ax2.set_title(
        "B. Trait Richness Rollup Across Ranks (Trait-Bearing Taxa)",
        fontsize=12,
        fontweight="bold",
        pad=12,
        loc="left",
    )
    ax2.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax2.set_axisbelow(True)

    ax2.text(
        0.03,
        0.96,
        "Trait Rollup in KG-Microbe:\nHigher ranks aggregate annotations\nfrom their classified descendants.\nMedian traits among trait-bearing taxa:\n"
        + "; ".join(
            f"{RANK_DISPLAY[r]}={taxa_master.loc[(taxa_master.taxonomic_rank == r) & (taxa_master.trait_count > 0), 'trait_count'].median():g}"
            for r in ("species", "genus", "phylum")
            if ((taxa_master.taxonomic_rank == r) & (taxa_master.trait_count > 0)).any()
        ),
        transform=ax2.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": COLOR_BG_LIGHT,
            "edgecolor": COLOR_SLATE_LIGHT,
            "alpha": 0.9,
        },
    )

    save_figure(fig, output_path)
    print(f"  Saved updated legacy Figure: {output_path} and .svg")


def distribution_stats(values: pd.Series) -> tuple[float, float, float, float]:
    """Return median, quartiles and mean; empty populations have no estimate."""
    if values.empty:
        return float("nan"), float("nan"), float("nan"), float("nan")
    return (
        float(values.median()),
        float(values.quantile(0.25)),
        float(values.quantile(0.75)),
        float(values.mean()),
    )


def write_readme(
    taxa: pd.DataFrame, ranks: pd.DataFrame, profiles: pd.DataFrame, output: Path
) -> None:
    """Write the coverage narrative from the same frames as the backing tables."""
    mapped = profiles.loc[profiles.has_mapped_traits]
    zero = taxa.loc[taxa.trait_count.eq(0)]
    metrics = {
        "Total profiles": len(profiles),
        "Profiles with mapped traits": len(mapped),
        "Zero-trait profiles": len(profiles) - len(mapped),
        "Unique taxa": len(taxa),
        "Trait-bearing taxa": int(taxa.has_traits.sum()),
        "Zero-trait taxa": len(zero),
        "Taxon signature mentions": int(taxa.signature_occurrences.sum()),
        "Zero-trait taxa recurring in ≥20 profiles": int(zero.profile_count.ge(20).sum()),
    }
    lines = [
        "# Taxa coverage and profile annotation mapping",
        "",
        f"Population: {HUMAN_HOST_SPECIES}; sequencing types {', '.join(CORPUS_SEQUENCING_TYPES)}.",
        "No contrast-class, control-tier, pairing or direction-orientation restriction is applied.",
        "Profiles are study / experiment / curated direction; taxa and signature mentions",
        "are restricted to those same profiles. The matrix is restricted to the profile vocabulary,",
        "including MediaDive columns if present. Source hashes and exclusions are in `provenance.json`.",
        "",
        "| Metric | Count |",
        "| --- | ---: |",
        *(f"| {label} | {value} |" for label, value in metrics.items()),
        "",
        "Zero-trait profiles remain in the coverage tables. Only the mapped-profile histogram",
        "and mapped-profile richness boxplot omit them; this audit does not declare eligibility",
        "for downstream statistical tests. No recorded annotation is not biological absence.",
        "",
        "| Rank | Unique taxa | Zero-trait taxa | Trait-bearing (%) | Signature mentions | Profiles with rank |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in ranks.itertuples():
        lines.append(
            f"| {row.display_rank} | {row.total_taxa} | {row.zero_trait_taxa} | {row.trait_bearing_pct:.2f} | {row.total_occurrences} | {row.profiles_with_rank} |"
        )
    lines.extend(
        ["", "| Distribution | Median | Q1 | Q3 | Mean |", "| --- | ---: | ---: | ---: | ---: |"]
    )
    for label, values in [
        ("Traits per mapped profile", mapped.mapped_trait_count),
        ("Mapped taxa per mapped profile", mapped.mapped_trait_taxa),
        ("Traits per reported taxon", taxa.trait_count),
        ("Traits per trait-bearing taxon", taxa.loc[taxa.has_traits, "trait_count"]),
    ]:
        display = [
            f"{v:g}" if np.isfinite(v) else "not estimable" for v in distribution_stats(values)
        ]
        lines.append(f"| {label} | {' | '.join(display)} |")
    lines.extend(
        [
            "",
            "Annotation counts can include traits aggregated from descendants; they do not measure",
            "physiology, patient-level abundance or activity. Suggested zero-trait diagnoses are",
            "heuristics for review, not verified causes or permission to propagate annotations downward.",
            "",
            "Regenerate after rebuilding `data/derived`: `just taxa-coverage`.",
            "The five TSVs, this README, and generated figures are inventoried by SHA-256 in `provenance.json`.",
            "`profile_coverage_by_rank.tsv` uses `pct_of_mapped_profiles`, whose denominator is the",
            "observed count of mapped profiles. Rank totals include an other/unrecorded category.",
            "The columns consumed by the slide and plate builders remain available in the rank and taxon detail TSVs.",
        ]
    )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_coverage(
    *,
    output_dir: Path,
    taxa: Path = TAXA_PATH,
    matrix: Path = TAXA_TRAIT_MATRIX_PATH,
    profiles: Path = PROFILES_WIDE_PATH,
    experiments: Path = EXPERIMENTS_PATH,
    tables_only: bool = False,
) -> None:
    """Validate first, then write data-derived outputs and a consumed-byte manifest."""
    inputs = {"taxa": taxa, "matrix": matrix, "profiles": profiles, "experiments": experiments}
    input_records = {
        name: {"file": path.name, "sha256": file_sha256(path), "size_bytes": path.stat().st_size}
        for name, path in inputs.items()
    }

    # 1. Load data and compute statistics
    taxa_master, rank_df, prof_rank_df, prof_summary, _prof_taxa, top_zero_df = (
        load_and_process_data(**inputs)
    )
    if (
        tables_only
        and output_dir.exists()
        and any(p.suffix in {".png", ".svg"} for p in output_dir.iterdir())
    ):
        raise ValueError(
            "Use an empty output directory for --tables-only; existing figures would become stale"
        )
    out_dir = output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "provenance.json").unlink(missing_ok=True)

    # 2. Save TSV Audit Tables
    rank_tsv_path = out_dir / "taxa_coverage_by_rank.tsv"
    write_tsv(rank_df, rank_tsv_path)
    print(f"  Written rank coverage summary: {rank_tsv_path}")

    prof_rank_tsv_path = out_dir / "profile_coverage_by_rank.tsv"
    write_tsv(prof_rank_df, prof_rank_tsv_path)
    print(f"  Written profile coverage by rank summary: {prof_rank_tsv_path}")

    prof_summary_tsv_path = out_dir / "profile_trait_mapping_summary.tsv"
    write_tsv(prof_summary, prof_summary_tsv_path)
    print(
        f"  Written profile-level trait mapping summary ({len(prof_summary)} profiles): {prof_summary_tsv_path}"
    )

    details_tsv_path = out_dir / "taxa_kgmicrobe_coverage_details.tsv"
    write_tsv(taxa_master, details_tsv_path)
    print(f"  Written complete taxa audit table ({len(taxa_master)} taxa): {details_tsv_path}")

    top_zero_tsv_path = out_dir / "top_zero_trait_bugsigdb_taxa.tsv"
    write_tsv(top_zero_df, top_zero_tsv_path)
    print(f"  Written top 100 zero-trait diagnosis table: {top_zero_tsv_path}")

    write_readme(taxa_master, rank_df, prof_summary, out_dir / "README.md")

    # 3. Generate visualizations unless explicitly testing/exporting tables only.
    if not tables_only:
        render_figures(taxa_master, rank_df, prof_rank_df, prof_summary, out_dir)
    for name, path in inputs.items():
        if file_sha256(path) != input_records[name]["sha256"]:
            raise ValueError(f"Input changed during coverage generation: {name}")
    output_names = [
        "taxa_coverage_by_rank.tsv",
        "profile_coverage_by_rank.tsv",
        "profile_trait_mapping_summary.tsv",
        "taxa_kgmicrobe_coverage_details.tsv",
        "top_zero_trait_bugsigdb_taxa.tsv",
        "README.md",
    ]
    if not tables_only:
        output_names.extend(
            f"{stem}.{suffix}" for stem in FIGURE_STEMS for suffix in ("png", "svg")
        )
    provenance = {
        "schema_version": "1.0",
        "recipe": "just taxa-coverage",
        "generated_at": datetime.now(UTC).isoformat(),
        "tables_only": tables_only,
        "population": {
            "host_species": HUMAN_HOST_SPECIES,
            "sequencing_types": list(CORPUS_SEQUENCING_TYPES),
            "contrast_filter": "none",
            "direction_basis": "curated",
            "zero_trait_profiles": "included",
        },
        "counts": prof_summary.attrs["population_counts"],
        "inputs": input_records,
        "producer": {
            "script_sha256": file_sha256(Path(__file__)),
            "corpus_policy_sha256": file_sha256(REPO_ROOT / "src/mtp/human_host_corpus.py"),
            "uv_lock_sha256": file_sha256(REPO_ROOT / "uv.lock"),
        },
        "outputs": {
            name: {
                "sha256": file_sha256(out_dir / name),
                "size_bytes": (out_dir / name).stat().st_size,
            }
            for name in output_names
        },
    }
    (out_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )


FIGURE_STEMS = (
    "taxa_coverage_histogram_rank_stratified",
    "taxa_coverage_occurrences_vs_unique",
    "taxa_profile_distribution_by_rank",
    "histogram1_mapped_traits_and_taxa_per_profile",
    "histogram2_taxa_trait_counts_and_zero_trait_frequency",
    "taxa_trait_count_distribution",
)


def render_figures(
    taxa_master: pd.DataFrame,
    rank_df: pd.DataFrame,
    prof_rank_df: pd.DataFrame,
    prof_summary: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Render the existing six figure filenames from one selected corpus."""
    fig1_path = out_dir / "taxa_coverage_histogram_rank_stratified.png"
    plot_rank_stratified_histogram(rank_df, fig1_path, prof_summary)

    fig2_path = out_dir / "taxa_coverage_occurrences_vs_unique.png"
    plot_occurrences_vs_unique(rank_df, taxa_master, prof_summary, fig2_path)

    fig3_prof_rank_path = out_dir / "taxa_profile_distribution_by_rank.png"
    plot_profile_distribution_by_rank(prof_rank_df, prof_summary, fig3_prof_rank_path)

    hist1_path = out_dir / "histogram1_mapped_traits_and_taxa_per_profile.png"
    plot_histogram1_mapped_profiles(prof_summary, hist1_path)

    hist2_path = out_dir / "histogram2_taxa_trait_counts_and_zero_trait_frequency.png"
    legacy_fig3_path = out_dir / "taxa_trait_count_distribution.png"
    plot_histogram2_taxa_traits_and_zero_frequency(
        taxa_master, prof_summary, hist2_path, legacy_path=legacy_fig3_path
    )


def main() -> None:
    """Execute the explicit-I/O coverage recipe."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--taxa", type=Path, default=TAXA_PATH)
    parser.add_argument("--matrix", type=Path, default=TAXA_TRAIT_MATRIX_PATH)
    parser.add_argument("--profiles", type=Path, default=PROFILES_WIDE_PATH)
    parser.add_argument("--experiments", type=Path, default=EXPERIMENTS_PATH)
    parser.add_argument(
        "--tables-only",
        action="store_true",
        help="Write tables, README and provenance without rendering figures (use a separate directory)",
    )
    run_coverage(**vars(parser.parse_args()))


if __name__ == "__main__":
    main()
