"""Matrix trait overlap, Tier 1+2 taxon coverage cascade, and Parkinson's study overlap (#1053).

Consumes a verified core workflow run (``run.json``) and emits five TSV tables plus
``summary.json``, ``core_consumption.json``, and ``provenance.json``:

1. ``matrix_trait_overlap.tsv`` — pairwise Jaccard and taxon counts across a declared
   list of retained (non-MediaDive) traits in ``taxa_trait_matrix.tsv``.
2. ``tier12_coverage_by_rank.tsv`` — Tier 1+2 reported-taxon coverage by taxonomic rank,
   genus-mention coverage across the active signatures, and single-study frequency splits
   for labeled versus unlabeled taxa under the retained non-MediaDive trait set.
3. ``filter_cascade.tsv`` — trait-axis and study-axis filter cascade counts from the
   KG-Microbe trait catalog and full BugSigDB export down to the Tier 1+2 corpus.
4. ``parkinsons_stool_study_overlap_summary.tsv`` — mean and median pairwise study taxon
   Jaccard across Parkinson's disease stool studies in Tier 1+2 and the full human export.
5. ``parkinsons_stool_study_pairs.tsv`` — pairwise study taxon overlap ledger backing
   the Parkinson's disease stool summary.
"""

from __future__ import annotations

import argparse
import json
import statistics
from collections.abc import Sequence
from itertools import combinations
from pathlib import Path
from typing import Any

import pandas as pd

from mtp.human_host_corpus import (
    CORPUS_SEQUENCING_TYPES,
    HUMAN_HOST_SPECIES,
    experiment_in_profile_population,
)
from mtp.matrix_overlap_inputs import (
    EXPERIMENT_KEYS,
    PROFILE_KEYS,
    population_inputs,
    select_experiments,
    validated_matrix,
)
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_lists import DeclaredTraitList, load_trait_list, trait_labels
from mtp.workflow.coverage import feature_class
from mtp.workflow.output_roots import (
    add_temporary_output_option,
    require_durable_output,
    require_durable_output_path,
)

CARBON_SUBSTRATE_PREFIX = "METPO:2000006|"
SCRIPT_RELATIVE_PATH = "src/mtp/matrix_overlap_coverage.py"
TRAIT_LISTS_RELATIVE_PATH = "src/mtp/trait_lists.py"
METHODS_NOTE_RELATIVE_PATH = "methods/matrix_coverage_and_pd_overlap.md"


def validate_durable_output_dir(
    output_dir: Path,
    *,
    allow_temporary_output: bool = False,
) -> Path:
    """Reject temporary output directories unless explicitly opted in (#1049, #1051)."""
    return require_durable_output_path(output_dir.expanduser(), allowed=allow_temporary_output)


def _is_mediadive_column(column_id: str) -> bool:
    return feature_class(column_id) == "culture_medium"


def _classify_retained_trait(column_id: str) -> str:
    if (
        column_id.startswith(CARBON_SUBSTRATE_PREFIX)
        or "CHEBI:" in column_id
        or "carbon_substrates:" in column_id
    ):
        return "chebi_carbon_source_traits"
    if "pathways:" in column_id:
        return "pathway_traits"
    if "GO:" in column_id:
        return "go_process_traits"
    return "metpo_phenotype_traits"


def load_declared_traits(trait_list_path: Path) -> DeclaredTraitList:
    """Load a declared trait list via ``mtp.trait_lists`` and require at least two traits."""
    declared = load_trait_list(trait_list_path)
    if len(declared.trait_ids) < 2:
        raise ValueError("At least two trait IDs are required to compute pairwise overlap")
    mediadive = [trait_id for trait_id in declared.trait_ids if _is_mediadive_column(trait_id)]
    if mediadive:
        raise ValueError(
            f"MediaDive growth-medium columns are excluded from trait overlap: {mediadive}"
        )
    return declared


def compute_matrix_trait_overlap(
    matrix: pd.DataFrame,
    columns: pd.DataFrame,
    trait_ids: Sequence[str],
) -> list[dict[str, Any]]:
    """Compute pairwise Jaccard overlap and taxon counts for declared matrix traits."""
    if "taxon_id" not in matrix.columns:
        raise ValueError("taxa_trait_matrix.tsv must contain a taxon_id column")
    if len(trait_ids) < 2:
        raise ValueError("At least two trait IDs are required for pairwise matrix overlap")
    matrix = validated_matrix(matrix, columns)

    matrix_cols = set(matrix.columns.astype(str))
    for trait_id in trait_ids:
        if _is_mediadive_column(trait_id):
            raise ValueError(
                f"Declared trait {trait_id!r} is a MediaDive growth-medium column; "
                "only retained non-MediaDive traits are valid"
            )
        if trait_id not in matrix_cols:
            raise ValueError(
                f"Declared trait {trait_id!r} is absent from taxa_trait_matrix.tsv columns"
            )
    label_by_id = trait_labels(columns, trait_ids)

    masks: dict[str, pd.Series] = {}
    for trait_id in trait_ids:
        masks[trait_id] = matrix[trait_id] > 0.0

    matrix_taxa_total = len(matrix)
    rows: list[dict[str, Any]] = []
    for trait_a, trait_b in combinations(trait_ids, 2):
        mask_a = masks[trait_a]
        mask_b = masks[trait_b]
        trait_a_taxa = int(mask_a.sum())
        trait_b_taxa = int(mask_b.sum())
        intersection_taxa = int((mask_a & mask_b).sum())
        union_taxa = int((mask_a | mask_b).sum())
        jaccard_fraction = float(intersection_taxa / union_taxa) if union_taxa > 0 else None
        rows.append(
            {
                "trait_a_id": trait_a,
                "trait_a_label": label_by_id[trait_a],
                "trait_b_id": trait_b,
                "trait_b_label": label_by_id[trait_b],
                "matrix_taxa_total": matrix_taxa_total,
                "trait_a_taxa": trait_a_taxa,
                "trait_b_taxa": trait_b_taxa,
                "intersection_taxa": intersection_taxa,
                "union_taxa": union_taxa,
                "jaccard_fraction": jaccard_fraction,
                "jaccard_pct": _percentage(jaccard_fraction),
            }
        )
    return rows


def _percentage(value: float | None) -> float | None:
    return value * 100.0 if value is not None else None


def _pct(numerator: int, denominator: int) -> float | None:
    return float(numerator * 100.0 / denominator) if denominator > 0 else None


def _median_or_none(values: Sequence[int | float]) -> float | None:
    return float(statistics.median(values)) if values else None


def _host_case_variant_counts(
    source_experiments: pd.DataFrame, source_taxa: pd.DataFrame
) -> dict[str, int]:
    """Count source host case variants excluded by the unchanged human-host rule.

    These source-level diagnostics precede sequencing, disease and site selection;
    an experiment without taxon occurrences contributes only to the experiment count.
    """
    hosts = source_experiments["host_species"].fillna("").astype(str).str.strip()
    variant_mask = hosts.str.casefold().eq(HUMAN_HOST_SPECIES.casefold()) & ~hosts.eq(
        HUMAN_HOST_SPECIES
    )
    variants = source_experiments.loc[variant_mask]
    profile_keys = select_experiments(source_taxa[PROFILE_KEYS], variants).drop_duplicates()
    return {
        "source_host_case_variant_experiments": len(variants),
        "source_host_case_variant_profile_keys": len(profile_keys),
    }


def compute_tier12_coverage_and_cascade(
    *,
    matrix: pd.DataFrame,
    columns: pd.DataFrame,
    tier12_experiments: pd.DataFrame,
    tier12_taxa: pd.DataFrame,
    tier12_profile_keys: pd.DataFrame,
    source_experiments: pd.DataFrame,
    source_taxa: pd.DataFrame,
    kg_catalog_traits_total: int | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Compute Tier 1+2 rank coverage, single-study frequencies, and the filter cascade."""
    matrix = validated_matrix(matrix, columns)
    population = population_inputs(
        tier12_experiments, tier12_taxa, tier12_profile_keys, source_experiments, source_taxa
    )
    all_column_ids = list(columns["column_id"])
    retained_trait_cols = [col for col in all_column_ids if not _is_mediadive_column(col)]
    mediadive_matrix_cols = [col for col in all_column_ids if _is_mediadive_column(col)]

    if retained_trait_cols:
        numeric_retained = matrix[retained_trait_cols]
        has_retained_trait = (numeric_retained > 0.0).any(axis=1)
    else:
        has_retained_trait = pd.Series(False, index=matrix.index)

    labeled_taxon_ids: set[str] = {
        str(taxon_id)
        for taxon_id, is_labeled in zip(
            matrix["taxon_id"].astype(str), has_retained_trait, strict=True
        )
        if bool(is_labeled)
    }
    matrix_labeled_taxa_count = len(labeled_taxon_ids)

    active_keys = population.active_keys
    taxa_df = population.active_taxa.copy()
    taxa_df["taxonomic_rank"] = (
        taxa_df["taxonomic_rank"]
        .fillna("unknown")
        .astype(str)
        .str.strip()
        .str.lower()
        .replace({"": "unknown"})
    )
    raw_occurrence_rows = len(taxa_df)

    # Deduplicate repeated within-signature occurrences to distinct (profile, taxon_id) mentions
    mentions_df = taxa_df.drop_duplicates(
        subset=["study_id", "experiment_id", "inc_dec", "taxon_id"]
    ).copy()
    mentions_df["is_labeled"] = mentions_df["taxon_id"].isin(labeled_taxon_ids)

    taxon_summary_rows: list[dict[str, Any]] = []
    for taxon_id, group in mentions_df.groupby("taxon_id", sort=True):
        rank_val = str(group["taxonomic_rank"].iloc[0])
        n_studies = int(group["study_id"].nunique())
        n_mentions = len(group)
        taxon_summary_rows.append(
            {
                "taxon_id": str(taxon_id),
                "taxonomic_rank": rank_val,
                "is_labeled": str(taxon_id) in labeled_taxon_ids,
                "n_studies": n_studies,
                "n_mentions": n_mentions,
            }
        )
    taxon_summary_df = pd.DataFrame(
        taxon_summary_rows,
        columns=["taxon_id", "taxonomic_rank", "is_labeled", "n_studies", "n_mentions"],
    )

    def _build_rank_row(
        rank_label: str, sub_taxa: pd.DataFrame, sub_mentions: pd.DataFrame
    ) -> dict[str, Any]:
        total_taxa = len(sub_taxa)
        if total_taxa > 0:
            labeled_sub = sub_taxa[sub_taxa["is_labeled"]]
            unlabeled_sub = sub_taxa[~sub_taxa["is_labeled"]]
        else:
            labeled_sub = sub_taxa
            unlabeled_sub = sub_taxa

        labeled_taxa = len(labeled_sub)
        unlabeled_taxa = len(unlabeled_sub)

        total_mentions = len(sub_mentions)
        labeled_mentions = int(sub_mentions["is_labeled"].sum()) if total_mentions > 0 else 0
        unlabeled_mentions = total_mentions - labeled_mentions

        lab_studies = (
            [int(v) for v in labeled_sub["n_studies"].tolist()] if labeled_taxa > 0 else []
        )
        unlab_studies = (
            [int(v) for v in unlabeled_sub["n_studies"].tolist()] if unlabeled_taxa > 0 else []
        )

        lab_1 = sum(1 for v in lab_studies if v == 1)
        lab_le2 = sum(1 for v in lab_studies if v <= 2)
        unlab_1 = sum(1 for v in unlab_studies if v == 1)
        unlab_le2 = sum(1 for v in unlab_studies if v <= 2)

        return {
            "taxonomic_rank": rank_label,
            "total_taxa": total_taxa,
            "labeled_taxa": labeled_taxa,
            "unlabeled_taxa": unlabeled_taxa,
            "labeled_taxa_pct": _pct(labeled_taxa, total_taxa),
            "unlabeled_taxa_pct": _pct(unlabeled_taxa, total_taxa),
            "total_mentions": total_mentions,
            "labeled_mentions": labeled_mentions,
            "unlabeled_mentions": unlabeled_mentions,
            "labeled_mentions_pct": _pct(labeled_mentions, total_mentions),
            "unlabeled_mentions_pct": _pct(unlabeled_mentions, total_mentions),
            "labeled_1_study_taxa": lab_1,
            "labeled_1_study_pct": _pct(lab_1, labeled_taxa),
            "labeled_le2_study_taxa": lab_le2,
            "labeled_le2_study_pct": _pct(lab_le2, labeled_taxa),
            "labeled_median_studies": _median_or_none(lab_studies),
            "unlabeled_1_study_taxa": unlab_1,
            "unlabeled_1_study_pct": _pct(unlab_1, unlabeled_taxa),
            "unlabeled_le2_study_taxa": unlab_le2,
            "unlabeled_le2_study_pct": _pct(unlab_le2, unlabeled_taxa),
            "unlabeled_median_studies": _median_or_none(unlab_studies),
        }

    coverage_rows: list[dict[str, Any]] = [_build_rank_row("all", taxon_summary_df, mentions_df)]
    observed_ranks = (
        taxon_summary_df.groupby("taxonomic_rank")["taxon_id"]
        .count()
        .sort_values(ascending=False)
        .index.tolist()
        if not taxon_summary_df.empty
        else []
    )
    priority_ranks = ["genus", "species"]
    ordered_ranks = [r for r in priority_ranks if r in observed_ranks] + [
        r for r in observed_ranks if r not in priority_ranks
    ]
    for rank_name in ordered_ranks:
        coverage_rows.append(
            _build_rank_row(
                str(rank_name),
                taxon_summary_df[taxon_summary_df["taxonomic_rank"] == rank_name],
                mentions_df[mentions_df["taxonomic_rank"] == rank_name],
            )
        )

    # Trait cascade counts
    category_counts = {
        "metpo_phenotype_traits": 0,
        "chebi_carbon_source_traits": 0,
        "pathway_traits": 0,
        "go_process_traits": 0,
    }
    for col in retained_trait_cols:
        category_counts[_classify_retained_trait(col)] += 1

    catalog_total = (
        int(kg_catalog_traits_total) if kg_catalog_traits_total is not None else len(all_column_ids)
    )
    retained_traits_count = len(retained_trait_cols)
    catalog_dropped = catalog_total - retained_traits_count
    matrix_dropped = len(mediadive_matrix_cols)

    # Study cascade counts
    src_exp, src_taxa = population.source_experiments, population.source_taxa
    src_sig_keys = src_taxa[PROFILE_KEYS].drop_duplicates()

    all_hosts_studies = int(src_exp["study_id"].nunique())
    all_hosts_experiments = len(src_exp)
    all_hosts_signatures = len(src_sig_keys)

    if "sequencing_type" not in src_exp:
        raise ValueError("Source experiment metadata lacks sequencing_type")
    human_mask = src_exp.apply(
        lambda row: experiment_in_profile_population(
            row["host_species"],
            row["sequencing_type"],
            host_species="Homo sapiens",
            sequencing_types=CORPUS_SEQUENCING_TYPES,
        ),
        axis=1,
    )
    human_sig_keys = select_experiments(src_sig_keys, src_exp.loc[human_mask])
    full_human_studies = int(human_sig_keys["study_id"].nunique())
    full_human_experiments = len(human_sig_keys[EXPERIMENT_KEYS].drop_duplicates())
    full_human_signatures = len(human_sig_keys)

    dirs_per_exp = human_sig_keys.groupby(EXPERIMENT_KEYS)["inc_dec"].nunique()
    paired_experiments = dirs_per_exp[dirs_per_exp == 2].reset_index()[EXPERIMENT_KEYS]
    paired_sig_keys = select_experiments(human_sig_keys, paired_experiments)
    paired_human_studies = int(paired_sig_keys["study_id"].nunique())
    paired_human_experiments = len(paired_experiments)
    paired_human_signatures = len(paired_sig_keys)

    t12_exp = population.active_experiments

    tier12_studies = int(active_keys["study_id"].nunique())
    tier12_experiments_count = len(t12_exp)
    tier12_signatures = len(active_keys)
    tier12_conditions = int(
        t12_exp["condition"]
        .fillna("")
        .astype(str)
        .str.strip()
        .replace({"": pd.NA})
        .dropna()
        .nunique()
    )
    stool_mask = t12_exp["body_site"].eq("Feces")
    tier12_stool_studies = int(t12_exp.loc[stool_mask, "study_id"].nunique())
    tier12_stool_experiments = int(stool_mask.sum())

    cascade_rows: list[dict[str, Any]] = [
        {
            "dimension": "traits",
            "step_order": 1,
            "step_id": "kg_trait_catalog",
            "step_label": "KG-Microbe trait catalog (feature_catalog / matrix columns)",
            "retained_count": catalog_total,
            "dropped_from_previous": 0,
            "matrix_columns_count": len(all_column_ids),
            "matrix_dropped_from_previous": 0,
            "experiments_count": pd.NA,
            "signatures_count": pd.NA,
            "dropped_signatures_from_previous": pd.NA,
            "conditions_count": pd.NA,
            "stool_studies_count": pd.NA,
            "stool_experiments_count": pd.NA,
            "metpo_phenotype_traits": category_counts["metpo_phenotype_traits"],
            "chebi_carbon_source_traits": category_counts["chebi_carbon_source_traits"],
            "pathway_traits": category_counts["pathway_traits"],
            "go_process_traits": category_counts["go_process_traits"],
            "mediadive_traits": catalog_dropped,
        },
        {
            "dimension": "traits",
            "step_order": 2,
            "step_id": "retained_non_mediadive_traits",
            "step_label": "Retained non-MediaDive traits",
            "retained_count": retained_traits_count,
            "dropped_from_previous": catalog_dropped,
            "matrix_columns_count": retained_traits_count,
            "matrix_dropped_from_previous": matrix_dropped,
            "experiments_count": pd.NA,
            "signatures_count": pd.NA,
            "dropped_signatures_from_previous": pd.NA,
            "conditions_count": pd.NA,
            "stool_studies_count": pd.NA,
            "stool_experiments_count": pd.NA,
            "metpo_phenotype_traits": category_counts["metpo_phenotype_traits"],
            "chebi_carbon_source_traits": category_counts["chebi_carbon_source_traits"],
            "pathway_traits": category_counts["pathway_traits"],
            "go_process_traits": category_counts["go_process_traits"],
            "mediadive_traits": 0,
        },
        {
            "dimension": "studies",
            "step_order": 1,
            "step_id": "all_hosts_bugsigdb_export",
            "step_label": "All-hosts BugSigDB export table (experiments.tsv)",
            "retained_count": all_hosts_studies,
            "dropped_from_previous": 0,
            "matrix_columns_count": pd.NA,
            "matrix_dropped_from_previous": pd.NA,
            "experiments_count": all_hosts_experiments,
            "signatures_count": all_hosts_signatures,
            "dropped_signatures_from_previous": 0,
            "conditions_count": pd.NA,
            "stool_studies_count": pd.NA,
            "stool_experiments_count": pd.NA,
            "metpo_phenotype_traits": pd.NA,
            "chebi_carbon_source_traits": pd.NA,
            "pathway_traits": pd.NA,
            "go_process_traits": pd.NA,
            "mediadive_traits": pd.NA,
        },
        {
            "dimension": "studies",
            "step_order": 2,
            "step_id": "full_human_16s_wms",
            "step_label": "Canonical human 16S/WMS profile corpus",
            "retained_count": full_human_studies,
            "dropped_from_previous": all_hosts_studies - full_human_studies,
            "matrix_columns_count": pd.NA,
            "matrix_dropped_from_previous": pd.NA,
            "experiments_count": full_human_experiments,
            "signatures_count": full_human_signatures,
            "dropped_signatures_from_previous": all_hosts_signatures - full_human_signatures,
            "conditions_count": pd.NA,
            "stool_studies_count": pd.NA,
            "stool_experiments_count": pd.NA,
            "metpo_phenotype_traits": pd.NA,
            "chebi_carbon_source_traits": pd.NA,
            "pathway_traits": pd.NA,
            "go_process_traits": pd.NA,
            "mediadive_traits": pd.NA,
        },
        {
            "dimension": "studies",
            "step_order": 3,
            "step_id": "paired_human_16s_wms",
            "step_label": "Paired human 16S/WMS corpus (both increased and decreased lists)",
            "retained_count": paired_human_studies,
            "dropped_from_previous": full_human_studies - paired_human_studies,
            "matrix_columns_count": pd.NA,
            "matrix_dropped_from_previous": pd.NA,
            "experiments_count": paired_human_experiments,
            "signatures_count": paired_human_signatures,
            "dropped_signatures_from_previous": full_human_signatures - paired_human_signatures,
            "conditions_count": pd.NA,
            "stool_studies_count": pd.NA,
            "stool_experiments_count": pd.NA,
            "metpo_phenotype_traits": pd.NA,
            "chebi_carbon_source_traits": pd.NA,
            "pathway_traits": pd.NA,
            "go_process_traits": pd.NA,
            "mediadive_traits": pd.NA,
        },
        {
            "dimension": "studies",
            "step_order": 4,
            "step_id": "tier12_subject_controls",
            "step_label": "Tier 1+2 disease vs healthy subject controls",
            "retained_count": tier12_studies,
            "dropped_from_previous": paired_human_studies - tier12_studies,
            "matrix_columns_count": pd.NA,
            "matrix_dropped_from_previous": pd.NA,
            "experiments_count": tier12_experiments_count,
            "signatures_count": tier12_signatures,
            "dropped_signatures_from_previous": paired_human_signatures - tier12_signatures,
            "conditions_count": tier12_conditions,
            "stool_studies_count": tier12_stool_studies,
            "stool_experiments_count": tier12_stool_experiments,
            "metpo_phenotype_traits": pd.NA,
            "chebi_carbon_source_traits": pd.NA,
            "pathway_traits": pd.NA,
            "go_process_traits": pd.NA,
            "mediadive_traits": pd.NA,
        },
    ]

    all_row = coverage_rows[0]
    summary = {
        "matrix_taxa_total": len(matrix),
        "matrix_labeled_non_mediadive_taxa": matrix_labeled_taxa_count,
        "kg_catalog_traits_total": catalog_total,
        "retained_non_mediadive_traits": retained_traits_count,
        "mediadive_traits_dropped": catalog_dropped,
        "trait_category_counts": category_counts,
        "tier12_unique_taxa": int(all_row["total_taxa"]),
        "tier12_labeled_taxa": int(all_row["labeled_taxa"]),
        "tier12_unlabeled_taxa": int(all_row["unlabeled_taxa"]),
        "tier12_labeled_taxa_pct": all_row["labeled_taxa_pct"],
        "tier12_unlabeled_taxa_pct": all_row["unlabeled_taxa_pct"],
        "tier12_distinct_profile_taxon_mentions": int(all_row["total_mentions"]),
        "tier12_raw_occurrence_rows": raw_occurrence_rows,
        "tier12_input_occurrence_rows": population.input_occurrence_rows,
        "tier12_excluded_occurrence_rows": population.excluded_occurrence_rows,
        "tier12_excluded_profile_keys": population.excluded_profile_keys,
        "tier12_input_experiments": population.input_experiments,
        "tier12_excluded_experiments_without_active_profiles": population.excluded_experiments,
        "tier12_studies": tier12_studies,
        "tier12_experiments": tier12_experiments_count,
        "tier12_signatures": tier12_signatures,
        "tier12_conditions": tier12_conditions,
        "tier12_stool_studies": tier12_stool_studies,
        "tier12_stool_experiments": tier12_stool_experiments,
        **_host_case_variant_counts(src_exp, src_taxa),
    }
    return coverage_rows, cascade_rows, summary


def _is_parkinsons_stool_row(
    condition: object,
    body_site: object,
    *,
    exact_condition: bool = False,
) -> bool:
    if exact_condition:
        has_pd = str(condition).strip().lower() == "parkinson's disease"
    else:
        cond_tokens = [part.strip().lower() for part in str(condition).split(",") if part.strip()]
        has_pd = any(tok == "parkinson's disease" for tok in cond_tokens)
    is_feces = isinstance(body_site, str) and body_site == "Feces"
    return has_pd and is_feces


def _compute_scope_study_overlap(
    *,
    scope: str,
    scope_label: str,
    experiments_df: pd.DataFrame,
    taxa_df: pd.DataFrame,
    exact_condition: bool = False,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if experiments_df.empty or taxa_df.empty:
        return (
            {
                "scope": scope,
                "scope_label": scope_label,
                "condition_filter": "Parkinson's disease",
                "body_site_filter": "Feces",
                "n_studies": 0,
                "n_experiments": 0,
                "n_study_pairs": 0,
                "mean_jaccard": None,
                "mean_jaccard_pct": None,
                "median_jaccard": None,
                "median_jaccard_pct": None,
                "min_jaccard": None,
                "max_jaccard": None,
                "nonzero_pairs": 0,
                "zero_pairs": 0,
            },
            [],
        )

    exp_sub = experiments_df[
        experiments_df.apply(
            lambda r: _is_parkinsons_stool_row(
                r.get("condition", ""),
                r.get("body_site", ""),
                exact_condition=exact_condition,
            ),
            axis=1,
        )
    ].copy()
    taxa_sub = select_experiments(taxa_df, exp_sub)

    study_taxa: dict[str, set[str]] = {}
    for study_id, group in taxa_sub.groupby(taxa_sub["study_id"].astype(str), sort=True):
        taxa_set = set(group["taxon_id"])
        if taxa_set:
            study_taxa[str(study_id)] = taxa_set

    def _study_sort_key(study_id: str) -> tuple[int, int | str]:
        suffix = study_id.removeprefix("bsdb:")
        if suffix.isdigit():
            return (0, int(suffix))
        return (1, study_id)

    ordered_studies = sorted(study_taxa, key=_study_sort_key)
    active_exp_count = len(taxa_sub[EXPERIMENT_KEYS].drop_duplicates())

    pair_rows: list[dict[str, Any]] = []
    jaccards: list[float] = []
    for study_a, study_b in combinations(ordered_studies, 2):
        set_a = study_taxa[study_a]
        set_b = study_taxa[study_b]
        inter = len(set_a & set_b)
        union = len(set_a | set_b)
        jaccard = float(inter / union)  # both study sets are nonempty
        jaccards.append(jaccard)
        pair_rows.append(
            {
                "scope": scope,
                "study_a_id": study_a,
                "study_b_id": study_b,
                "study_a_taxa": len(set_a),
                "study_b_taxa": len(set_b),
                "intersection_taxa": inter,
                "union_taxa": union,
                "jaccard_fraction": jaccard,
                "jaccard_pct": jaccard * 100.0,
            }
        )

    mean_j = float(statistics.fmean(jaccards)) if jaccards else None
    median_j = float(statistics.median(jaccards)) if jaccards else None
    min_j = float(min(jaccards)) if jaccards else None
    max_j = float(max(jaccards)) if jaccards else None
    nonzero_pairs = sum(1 for v in jaccards if v > 0.0)

    summary_row = {
        "scope": scope,
        "scope_label": scope_label,
        "condition_filter": "Parkinson's disease",
        "body_site_filter": "Feces",
        "n_studies": len(ordered_studies),
        "n_experiments": active_exp_count,
        "n_study_pairs": len(pair_rows),
        "mean_jaccard": mean_j,
        "mean_jaccard_pct": _percentage(mean_j),
        "median_jaccard": median_j,
        "median_jaccard_pct": _percentage(median_j),
        "min_jaccard": min_j,
        "max_jaccard": max_j,
        "nonzero_pairs": nonzero_pairs,
        "zero_pairs": len(pair_rows) - nonzero_pairs,
    }
    return summary_row, pair_rows


def compute_parkinsons_stool_overlap(
    *,
    tier12_experiments: pd.DataFrame,
    tier12_taxa: pd.DataFrame,
    tier12_profile_keys: pd.DataFrame,
    source_experiments: pd.DataFrame,
    source_taxa: pd.DataFrame,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Compute pairwise study taxon Jaccard across Parkinson's stool studies."""
    population = population_inputs(
        tier12_experiments, tier12_taxa, tier12_profile_keys, source_experiments, source_taxa
    )
    t12_exp, t12_taxa = population.active_experiments, population.active_taxa
    source_experiments, source_taxa = population.source_experiments, population.source_taxa

    tier12_summary, tier12_pairs = _compute_scope_study_overlap(
        scope="tier12_paired",
        scope_label="Tier 1+2 paired disease vs healthy subject-control studies",
        experiments_df=t12_exp,
        taxa_df=t12_taxa,
    )
    tier12_summary.update(
        {
            "input_occurrence_rows": population.input_occurrence_rows,
            "excluded_outside_active_profile_rows": population.excluded_occurrence_rows,
            "excluded_outside_active_profile_keys": population.excluded_profile_keys,
            "input_experiments": population.input_experiments,
            "excluded_experiments_without_active_profiles": population.excluded_experiments,
        }
    )

    full_summary, full_pairs = _compute_scope_study_overlap(
        scope="full_export",
        scope_label="Full BugSigDB export studies (source/derived/experiments.tsv)",
        experiments_df=source_experiments.copy(),
        taxa_df=source_taxa,
    )

    human_mask = source_experiments.apply(
        lambda row: experiment_in_profile_population(
            row["host_species"],
            row.get("sequencing_type"),
            host_species="Homo sapiens",
            sequencing_types=None,
        ),
        axis=1,
    )
    human_source_exp = source_experiments.loc[human_mask].copy()
    human_summary, human_pairs = _compute_scope_study_overlap(
        scope="export_human_hosts",
        scope_label="Human-host BugSigDB export studies (host_species == Homo sapiens)",
        experiments_df=human_source_exp,
        taxa_df=source_taxa,
    )
    human_summary.update(_host_case_variant_counts(source_experiments, source_taxa))

    exact_summary, exact_pairs = _compute_scope_study_overlap(
        scope="full_export_exact_condition",
        scope_label="Full BugSigDB export studies (exact condition == Parkinson's disease)",
        experiments_df=source_experiments.copy(),
        taxa_df=source_taxa,
        exact_condition=True,
    )
    return (
        [tier12_summary, full_summary, human_summary, exact_summary],
        tier12_pairs + full_pairs + human_pairs + exact_pairs,
    )


def run_matrix_overlap_coverage(
    core_run: Path,
    output_dir: Path,
    *,
    trait_list_path: Path,
    code_commit: str | None = None,
    allow_temporary_output: bool = False,
) -> dict[str, Any]:
    """Run the #1053 matrix overlap, Tier 1+2 coverage cascade, and PD overlap consumer."""
    from mtp.workflow.consumers import (
        open_core_bundle,
        preflight_consumer_outputs,
        write_consumer_record,
    )
    from mtp.workflow.runtime import repository_state

    output_dir = validate_durable_output_dir(
        output_dir, allow_temporary_output=allow_temporary_output
    )
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise ValueError(f"--output-dir must be a new or empty directory: {output_dir}")
    declared = load_declared_traits(trait_list_path)
    manifest_path = core_run / "run.json" if core_run.is_dir() else core_run
    bundle = open_core_bundle(manifest_path)

    params = bundle.run.parameters or []
    params_by_name = {param.name: param.value for param in params}
    if len(params_by_name) != len(params):
        raise ValueError("Duplicate core run parameter names")
    strain_rollup = params_by_name.get("strain_rollup", "")
    raw_contrast_pop = bundle.run.contrast_population
    default_contrast_pop = (
        str(getattr(raw_contrast_pop, "value", raw_contrast_pop))
        if raw_contrast_pop is not None
        else ""
    )
    contrast_population = default_contrast_pop
    if (
        "contrast_population" in params_by_name
        and params_by_name["contrast_population"] != contrast_population
    ):
        raise ValueError("Core contrast_population parameter contradicts the run field")
    if contrast_population != "disease_vs_healthy_tier1_2_subject_controls":
        raise ValueError(
            "matrix-overlap-coverage requires a Tier 1+2 disease vs healthy subject-control "
            f"core run (contrast_population='disease_vs_healthy_tier1_2_subject_controls'); "
            f"got {contrast_population!r}"
        )
    raw_direction_basis = bundle.run.direction_basis
    default_direction_basis = (
        str(getattr(raw_direction_basis, "value", raw_direction_basis))
        if raw_direction_basis is not None
        else ""
    )
    direction_basis = default_direction_basis
    if "direction_basis" in params_by_name and params_by_name["direction_basis"] != direction_basis:
        raise ValueError("Core direction_basis parameter contradicts the run field")
    if direction_basis != "disease_vs_comparison":
        raise ValueError(
            "matrix-overlap-coverage requires direction_basis='disease_vs_comparison'; "
            f"got {direction_basis!r}"
        )

    overlap_tsv = output_dir / "matrix_trait_overlap.tsv"
    coverage_tsv = output_dir / "tier12_coverage_by_rank.tsv"
    cascade_tsv = output_dir / "filter_cascade.tsv"
    pd_summary_tsv = output_dir / "parkinsons_stool_study_overlap_summary.tsv"
    pd_pairs_tsv = output_dir / "parkinsons_stool_study_pairs.tsv"
    summary_json = output_dir / "summary.json"
    provenance_path = output_dir / "provenance.json"
    primary_outputs = (
        overlap_tsv,
        coverage_tsv,
        cascade_tsv,
        pd_summary_tsv,
        pd_pairs_tsv,
        summary_json,
    )
    preflight_consumer_outputs(
        bundle,
        output_dir / "core_consumption.json",
        outputs=(*primary_outputs, provenance_path),
        output_directory=output_dir,
    )

    available_products = {product.id for product in bundle.run.definition.products}
    required_sources = {"source_experiments", "source_reported_taxa"}
    if not required_sources.issubset(available_products):
        raise ValueError(
            "Full-export comparison requires source_experiments and source_reported_taxa; "
            f"missing {sorted(required_sources - available_products)}"
        )
    matrix_path = bundle.path("taxa_matrix")
    columns_path = bundle.path("feature_metadata")
    count_profiles_path = bundle.path("count_profiles")
    experiments_path = bundle.path("experiments")
    reported_taxa_path = bundle.path("reported_taxa")
    source_experiments_path = bundle.path("source_experiments")
    source_taxa_path = bundle.path("source_reported_taxa")
    feature_catalog_path = (
        bundle.path("feature_catalog") if "feature_catalog" in available_products else None
    )

    matrix_df = pd.read_csv(matrix_path, sep="\t", dtype=str, keep_default_na=False)
    columns_df = pd.read_csv(columns_path, sep="\t", dtype=str, keep_default_na=False)
    profile_keys_df = pd.read_csv(
        count_profiles_path,
        sep="\t",
        usecols=["study_id", "experiment_id", "inc_dec"],
        dtype=str,
        keep_default_na=False,
    )
    tier12_experiments_df = pd.read_csv(
        experiments_path, sep="\t", dtype=str, keep_default_na=False
    )
    tier12_taxa_df = pd.read_csv(reported_taxa_path, sep="\t", dtype=str, keep_default_na=False)
    source_experiments_df = pd.read_csv(
        source_experiments_path, sep="\t", dtype=str, keep_default_na=False
    )
    source_taxa_df = pd.read_csv(source_taxa_path, sep="\t", dtype=str, keep_default_na=False)

    kg_catalog_traits_total: int | None = None
    if feature_catalog_path is not None and feature_catalog_path.is_file():
        catalog_df = pd.read_csv(
            feature_catalog_path,
            sep="\t",
            usecols=["column_id"],
            low_memory=False,
        )
        kg_catalog_traits_total = int(catalog_df["column_id"].nunique())

    overlap_rows = compute_matrix_trait_overlap(matrix_df, columns_df, declared.trait_ids)
    coverage_rows, cascade_rows, coverage_summary = compute_tier12_coverage_and_cascade(
        matrix=matrix_df,
        columns=columns_df,
        tier12_experiments=tier12_experiments_df,
        tier12_taxa=tier12_taxa_df,
        tier12_profile_keys=profile_keys_df,
        source_experiments=source_experiments_df,
        source_taxa=source_taxa_df,
        kg_catalog_traits_total=kg_catalog_traits_total,
    )
    pd_summary_rows, pd_pair_rows = compute_parkinsons_stool_overlap(
        tier12_experiments=tier12_experiments_df,
        tier12_taxa=tier12_taxa_df,
        tier12_profile_keys=profile_keys_df,
        source_experiments=source_experiments_df,
        source_taxa=source_taxa_df,
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(overlap_rows).to_csv(overlap_tsv, sep="\t", index=False, float_format="%.12g")
    pd.DataFrame(coverage_rows).to_csv(coverage_tsv, sep="\t", index=False, float_format="%.12g")
    pd.DataFrame(cascade_rows).to_csv(cascade_tsv, sep="\t", index=False, float_format="%.12g")
    pd.DataFrame(pd_summary_rows).to_csv(
        pd_summary_tsv, sep="\t", index=False, float_format="%.12g"
    )
    pd.DataFrame(
        pd_pair_rows,
        columns=[
            "scope",
            "study_a_id",
            "study_b_id",
            "study_a_taxa",
            "study_b_taxa",
            "intersection_taxa",
            "union_taxa",
            "jaccard_fraction",
            "jaccard_pct",
        ],
    ).to_csv(pd_pairs_tsv, sep="\t", index=False, float_format="%.12g")

    declared_provenance = {
        **declared.provenance(),
        "source_name": declared.source_name,
        "sha256": declared.sha256,
    }
    summary_payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "consumer": "matrix-overlap-coverage",
        "core_run_id": bundle.run.run_id,
        "workflow_id": bundle.run.definition.id,
        "strain_rollup": strain_rollup,
        "contrast_population": contrast_population,
        "declared_traits": {
            **declared_provenance,
            "pair_count": len(overlap_rows),
        },
        "coverage_and_cascade": coverage_summary,
        "parkinsons_stool_overlap": pd_summary_rows,
    }
    summary_json.write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    extra_products = tuple(
        pid
        for pid in (
            "reported_taxa",
            "source_experiments",
            "source_reported_taxa",
        )
        if pid in available_products
    )
    consumed_ids = sorted(
        {
            "count_profiles",
            "experiments",
            "taxa_matrix",
            "feature_metadata",
            "feature_catalog",
            *extra_products,
        }
        & available_products
    )

    if code_commit is not None and code_commit.strip():
        resolved_code_commit = code_commit.strip()
    else:
        producer_commit, producer_dirty = repository_state(REPO_ROOT)
        resolved_code_commit = producer_commit + ("+dirty" if producer_dirty else "")
    script_path = REPO_ROOT / SCRIPT_RELATIVE_PATH
    methods_path = REPO_ROOT / METHODS_NOTE_RELATIVE_PATH
    script_sha256 = file_sha256(script_path)
    methods_sha256 = file_sha256(methods_path)
    run_root = bundle.manifest.parent.resolve()

    provenance_payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "consumer": "matrix-overlap-coverage",
        "core_run_id": bundle.run.run_id,
        "workflow_id": bundle.run.definition.id,
        "strain_rollup": strain_rollup,
        "contrast_population": contrast_population,
        "code_commit": resolved_code_commit,
        "core_git_commit_sha": bundle.run.git_commit,
        "script_sha256": script_sha256,
        "methods_sha256": methods_sha256,
        "scripts": {
            name: {"sha256": file_sha256(REPO_ROOT / name)}
            for name in (
                SCRIPT_RELATIVE_PATH,
                TRAIT_LISTS_RELATIVE_PATH,
                "src/mtp/matrix_overlap_inputs.py",
                "src/mtp/human_host_corpus.py",
                "src/mtp/workflow/coverage.py",
                "src/mtp/workflow/output_roots.py",
            )
        },
        "methods": {
            METHODS_NOTE_RELATIVE_PATH: {"sha256": methods_sha256},
        },
        "declared_traits": declared_provenance,
        "consumed_products": {
            pid: {
                "relative_path": bundle.path(pid).resolve().relative_to(run_root).as_posix(),
                "sha256": file_sha256(bundle.path(pid)),
            }
            for pid in consumed_ids
        },
        "outputs": {
            path.name: {
                "relative_path": path.name,
                "sha256": file_sha256(path),
            }
            for path in primary_outputs
        },
    }
    provenance_path.write_text(
        json.dumps(provenance_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    write_consumer_record(
        bundle,
        output_dir / "core_consumption.json",
        consumer="matrix-overlap-coverage",
        outputs=(*primary_outputs, provenance_path),
        feature_policy="exclude_mediadive",
        additional_products=extra_products,
    )
    return summary_payload


def build_parser() -> argparse.ArgumentParser:
    """Create CLI argument parser for ``mtp.matrix_overlap_coverage``."""
    parser = argparse.ArgumentParser(
        prog="python -m mtp.matrix_overlap_coverage",
        description=(
            "Compute taxon-by-trait matrix overlap, Tier 1+2 coverage by rank and filter "
            "cascade, and Parkinson's disease stool study overlap from a verified core run."
        ),
    )
    parser.add_argument(
        "--core-run",
        type=Path,
        required=True,
        help="Path to a verified core workflow run directory or run.json manifest.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Durable directory where TSV and JSON outputs will be written.",
    )
    parser.add_argument(
        "--trait-list",
        type=Path,
        required=True,
        help="Path to a newline-delimited file of trait column IDs (# comments allowed).",
    )
    parser.add_argument(
        "--code-commit",
        type=str,
        default=None,
        help=(
            "Git commit SHA (with optional +dirty suffix) of the sidecar code. "
            "Computed from the local repository when omitted."
        ),
    )
    add_temporary_output_option(parser)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point for ``mtp.matrix_overlap_coverage``."""
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    args.output_dir = args.output_dir.expanduser()
    require_durable_output(parser, args.output_dir, allowed=bool(args.allow_temporary_output))
    run_matrix_overlap_coverage(
        args.core_run,
        args.output_dir,
        trait_list_path=args.trait_list,
        code_commit=args.code_commit,
        allow_temporary_output=bool(args.allow_temporary_output),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
