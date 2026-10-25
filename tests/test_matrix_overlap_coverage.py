"""Tests for matrix trait overlap, Tier 1+2 coverage, and PD overlap (#1053)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd
import pytest

from mtp.matrix_overlap_coverage import (
    DEFAULT_MATRIX_OVERLAP_TRAITS,
    compute_matrix_trait_overlap,
    compute_parkinsons_stool_overlap,
    compute_tier12_coverage_and_cascade,
    load_declared_traits,
    main,
    run_matrix_overlap_coverage,
    validate_durable_output_dir,
)
from mtp.paths import REPO_ROOT, file_sha256

SUBSET_POPULATION = "disease_vs_healthy_tier1_2_subject_controls"
FAC_AER = "METPO:1000608"
FAC_ANA = "METPO:1000605"
ANA = "METPO:1000603"
CELL = "METPO:2000103|GO:0030245"
XYLAN = "METPO:2000103|GO:0045493"
GRAM_POS = "METPO:1000447"
GLUCOSE = "METPO:2000006|CHEBI:17234"
ARSENATE = "METPO:2000103|pathways:arsenate_reduction"
MEDIADIVE = "METPO:2000517|mediadive.medium:1"


def _synthetic_matrix_and_columns() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build a 6-taxon matrix where one taxon has only MediaDive and one has no traits."""
    matrix = pd.DataFrame(
        [
            {
                "taxon_id": "NCBITaxon:10",
                "taxon_name": "GenusA",
                "taxonomic_rank": "genus",
                FAC_AER: 2,
                FAC_ANA: 1,
                ANA: 0,
                CELL: 1,
                XYLAN: 1,
                GRAM_POS: 1,
                GLUCOSE: 1,
                ARSENATE: 0,
                MEDIADIVE: 3,
            },
            {
                "taxon_id": "NCBITaxon:11",
                "taxon_name": "GenusA sp1",
                "taxonomic_rank": "species",
                FAC_AER: 1,
                FAC_ANA: 1,
                ANA: 0,
                CELL: 0,
                XYLAN: 0,
                GRAM_POS: 1,
                GLUCOSE: 0,
                ARSENATE: 1,
                MEDIADIVE: 0,
            },
            {
                "taxon_id": "NCBITaxon:12",
                "taxon_name": "GenusB sp2",
                "taxonomic_rank": "species",
                FAC_AER: 1,
                FAC_ANA: 0,
                ANA: 1,
                CELL: 1,
                XYLAN: 1,
                GRAM_POS: 0,
                GLUCOSE: 0,
                ARSENATE: 0,
                MEDIADIVE: 1,
            },
            {
                "taxon_id": "NCBITaxon:13",
                "taxon_name": "FamilyC",
                "taxonomic_rank": "family",
                FAC_AER: 0,
                FAC_ANA: 1,
                ANA: 1,
                CELL: 0,
                XYLAN: 1,
                GRAM_POS: 1,
                GLUCOSE: 0,
                ARSENATE: 0,
                MEDIADIVE: 0,
            },
            {
                "taxon_id": "NCBITaxon:14",
                "taxon_name": "GenusMediaOnly",
                "taxonomic_rank": "genus",
                FAC_AER: 0,
                FAC_ANA: 0,
                ANA: 0,
                CELL: 0,
                XYLAN: 0,
                GRAM_POS: 0,
                GLUCOSE: 0,
                ARSENATE: 0,
                MEDIADIVE: 5,
            },
            {
                "taxon_id": "NCBITaxon:15",
                "taxon_name": "SpeciesUnlabeled",
                "taxonomic_rank": "species",
                FAC_AER: 0,
                FAC_ANA: 0,
                ANA: 0,
                CELL: 0,
                XYLAN: 0,
                GRAM_POS: 0,
                GLUCOSE: 0,
                ARSENATE: 0,
                MEDIADIVE: 0,
            },
        ]
    )
    columns = pd.DataFrame(
        [
            {"column_id": FAC_AER, "object_name": "facultative"},
            {"column_id": FAC_ANA, "object_name": "facultatively anaerobic"},
            {"column_id": ANA, "object_name": "anaerobic"},
            {"column_id": CELL, "object_name": "cellulose catabolic process"},
            {"column_id": XYLAN, "object_name": "xylan catabolic process"},
            {"column_id": GRAM_POS, "object_name": "gram positive"},
            {"column_id": GLUCOSE, "object_name": "D-glucose"},
            {"column_id": ARSENATE, "object_name": "arsenate_reduction"},
            {"column_id": MEDIADIVE, "object_name": "MediaDive medium 1"},
        ]
    )
    return matrix, columns


def test_matrix_trait_overlap_computes_exact_pairwise_jaccard_and_counts(tmp_path: Path) -> None:
    matrix, columns = _synthetic_matrix_and_columns()
    traits = (FAC_AER, FAC_ANA, CELL, XYLAN)
    rows = compute_matrix_trait_overlap(matrix, columns, traits)
    assert len(rows) == 6
    by_pair = {(row["trait_a_id"], row["trait_b_id"]): row for row in rows}

    # FAC_AER is positive on {10, 11, 12} (3 taxa); FAC_ANA on {10, 11, 13} (3 taxa).
    # Intersection = {10, 11} (2); Union = {10, 11, 12, 13} (4) -> Jaccard = 0.5 (50.0%).
    fac_pair = by_pair[(FAC_AER, FAC_ANA)]
    assert fac_pair["matrix_taxa_total"] == 6
    assert fac_pair["trait_a_taxa"] == 3
    assert fac_pair["trait_b_taxa"] == 3
    assert fac_pair["intersection_taxa"] == 2
    assert fac_pair["union_taxa"] == 4
    assert fac_pair["jaccard_fraction"] == pytest.approx(0.5)
    assert fac_pair["jaccard_pct"] == pytest.approx(50.0)
    assert fac_pair["trait_a_label"] == "facultative"
    assert fac_pair["trait_b_label"] == "facultatively anaerobic"

    # CELL is positive on {10, 12} (2); XYLAN on {10, 12, 13} (3).
    # Intersection = {10, 12} (2); Union = {10, 12, 13} (3) -> Jaccard = 2/3.
    poly_pair = by_pair[(CELL, XYLAN)]
    assert poly_pair["intersection_taxa"] == 2
    assert poly_pair["union_taxa"] == 3
    assert poly_pair["jaccard_fraction"] == pytest.approx(2.0 / 3.0)

    # Custom trait list file parsing and validation
    list_file = tmp_path / "traits.txt"
    list_file.write_text(f"# Declared traits\n{FAC_AER}\n{FAC_ANA}\n", encoding="utf-8")
    loaded = load_declared_traits(trait_list_path=list_file)
    assert loaded.trait_ids == (FAC_AER, FAC_ANA)
    assert loaded.source_name == "traits.txt"

    with pytest.raises(ValueError, match="MediaDive"):
        compute_matrix_trait_overlap(matrix, columns, (FAC_AER, MEDIADIVE))
    with pytest.raises(ValueError, match="absent"):
        compute_matrix_trait_overlap(matrix, columns, (FAC_AER, "METPO:9999999"))


def test_tier12_coverage_by_rank_and_filter_cascade() -> None:
    matrix, columns = _synthetic_matrix_and_columns()
    # Full source has 4 experiments across 3 human studies + 1 mouse study:
    # - s1/e1: PD Feces, paired (inc + dec), Tier 1+2
    # - s2/e1: PD Feces, paired (inc + dec), Tier 1+2
    # - s2/e2: CRC Oral, paired (inc + dec), dropped at Tier 1+2 selection
    # - s3/e1: PD Feces, unpaired (inc only), dropped at paired step
    # - s4/e1: Mouse 16S, excluded from human 16S/WMS step 1
    source_experiments = pd.DataFrame(
        [
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1-1",
                "host_species": "Homo sapiens",
                "sequencing_type": "16S",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
            },
            {
                "study_id": "bsdb:2",
                "experiment_id": "bsdb:2-1",
                "host_species": "Homo sapiens",
                "sequencing_type": "WMS",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
            },
            {
                "study_id": "bsdb:2",
                "experiment_id": "bsdb:2-2",
                "host_species": "Homo sapiens",
                "sequencing_type": "16S",
                "condition": "Colorectal cancer",
                "body_site": "Saliva",
            },
            {
                "study_id": "bsdb:3",
                "experiment_id": "bsdb:3-1",
                "host_species": "Homo sapiens",
                "sequencing_type": "16S",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
            },
            {
                "study_id": "bsdb:4",
                "experiment_id": "bsdb:4-1",
                "host_species": "Mus musculus",
                "sequencing_type": "16S",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
            },
        ]
    )
    source_taxa = pd.DataFrame(
        [
            # bsdb:1-1 inc: GenusA (labeled), GenusMediaOnly (unlabeled under non-MediaDive)
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1-1",
                "direction": "increased",
                "taxon_id": "NCBITaxon:10",
                "taxon_name": "GenusA",
                "taxonomic_rank": "genus",
            },
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1-1",
                "direction": "increased",
                "taxon_id": "NCBITaxon:14",
                "taxon_name": "GenusMediaOnly",
                "taxonomic_rank": "genus",
            },
            # bsdb:1-1 dec: GenusA sp1 (labeled, species), SpeciesUnlabeled (unlabeled, species)
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1-1",
                "direction": "decreased",
                "taxon_id": "NCBITaxon:11",
                "taxon_name": "GenusA sp1",
                "taxonomic_rank": "species",
            },
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1-1",
                "direction": "decreased",
                "taxon_id": "NCBITaxon:15",
                "taxon_name": "SpeciesUnlabeled",
                "taxonomic_rank": "species",
            },
            # bsdb:2-1 inc: GenusA (labeled, genus - now in 2 studies!)
            {
                "study_id": "bsdb:2",
                "experiment_id": "bsdb:2-1",
                "direction": "increased",
                "taxon_id": "NCBITaxon:10",
                "taxon_name": "GenusA",
                "taxonomic_rank": "genus",
            },
            # bsdb:2-1 dec: GenusB sp2 (labeled, species), FamilyC (labeled, family)
            {
                "study_id": "bsdb:2",
                "experiment_id": "bsdb:2-1",
                "direction": "decreased",
                "taxon_id": "NCBITaxon:12",
                "taxon_name": "GenusB sp2",
                "taxonomic_rank": "species",
            },
            {
                "study_id": "bsdb:2",
                "experiment_id": "bsdb:2-1",
                "direction": "decreased",
                "taxon_id": "NCBITaxon:13",
                "taxon_name": "FamilyC",
                "taxonomic_rank": "family",
            },
            # bsdb:2-2 (paired, but non-Tier12): inc + dec
            {
                "study_id": "bsdb:2",
                "experiment_id": "bsdb:2-2",
                "direction": "increased",
                "taxon_id": "NCBITaxon:10",
                "taxon_name": "GenusA",
                "taxonomic_rank": "genus",
            },
            {
                "study_id": "bsdb:2",
                "experiment_id": "bsdb:2-2",
                "direction": "decreased",
                "taxon_id": "NCBITaxon:11",
                "taxon_name": "GenusA sp1",
                "taxonomic_rank": "species",
            },
            # bsdb:3-1 (unpaired: inc only)
            {
                "study_id": "bsdb:3",
                "experiment_id": "bsdb:3-1",
                "direction": "increased",
                "taxon_id": "NCBITaxon:12",
                "taxon_name": "GenusB sp2",
                "taxonomic_rank": "species",
            },
            # bsdb:4-1 (mouse: inc + dec)
            {
                "study_id": "bsdb:4",
                "experiment_id": "bsdb:4-1",
                "direction": "increased",
                "taxon_id": "NCBITaxon:10",
                "taxon_name": "GenusA",
                "taxonomic_rank": "genus",
            },
            {
                "study_id": "bsdb:4",
                "experiment_id": "bsdb:4-1",
                "direction": "decreased",
                "taxon_id": "NCBITaxon:11",
                "taxon_name": "GenusA sp1",
                "taxonomic_rank": "species",
            },
        ]
    )
    tier12_experiments = source_experiments[
        source_experiments["experiment_id"].isin(["bsdb:1-1", "bsdb:2-1"])
    ].copy()
    tier12_taxa = source_taxa[source_taxa["experiment_id"].isin(["bsdb:1-1", "bsdb:2-1"])].copy()
    tier12_profile_keys = pd.DataFrame(
        [
            {"study_id": "bsdb:1", "experiment_id": "bsdb:1-1", "inc_dec": "inc"},
            {"study_id": "bsdb:1", "experiment_id": "bsdb:1-1", "inc_dec": "dec"},
            {"study_id": "bsdb:2", "experiment_id": "bsdb:2-1", "inc_dec": "inc"},
            {"study_id": "bsdb:2", "experiment_id": "bsdb:2-1", "inc_dec": "dec"},
        ]
    )

    coverage_rows, cascade_rows, summary = compute_tier12_coverage_and_cascade(
        matrix=matrix,
        columns=columns,
        tier12_experiments=tier12_experiments,
        tier12_taxa=tier12_taxa,
        tier12_profile_keys=tier12_profile_keys,
        source_experiments=source_experiments,
        source_taxa=source_taxa,
        kg_catalog_traits_total=12,
    )

    by_rank = {row["taxonomic_rank"]: row for row in coverage_rows}
    # Across Tier 1+2 (bsdb:1-1 and bsdb:2-1):
    # Unique taxa = 6 (10, 11, 12, 13 labeled; 14, 15 unlabeled because 14 is MediaDive-only!)
    assert by_rank["all"]["total_taxa"] == 6
    assert by_rank["all"]["labeled_taxa"] == 4
    assert by_rank["all"]["unlabeled_taxa"] == 2
    assert by_rank["all"]["total_mentions"] == 7
    assert by_rank["all"]["labeled_mentions"] == 5
    assert by_rank["all"]["unlabeled_mentions"] == 2
    # Labeled taxa study counts: 10 in 2 studies; 11, 12, 13 in 1 study -> 3/4 (75%) in 1 study
    assert by_rank["all"]["labeled_1_study_taxa"] == 3
    assert by_rank["all"]["labeled_1_study_pct"] == pytest.approx(75.0)
    # Unlabeled taxa study counts: 14, 15 are both in 1 study -> 2/2 (100%) in 1 study
    assert by_rank["all"]["unlabeled_1_study_taxa"] == 2
    assert by_rank["all"]["unlabeled_1_study_pct"] == pytest.approx(100.0)

    # Genus rank: taxa 10 (labeled, 2 mentions) and 14 (unlabeled, 1 mention) -> 2/3 mentions
    assert by_rank["genus"]["total_taxa"] == 2
    assert by_rank["genus"]["labeled_taxa"] == 1
    assert by_rank["genus"]["unlabeled_taxa"] == 1
    assert by_rank["genus"]["total_mentions"] == 3
    assert by_rank["genus"]["labeled_mentions"] == 2
    assert by_rank["genus"]["labeled_mentions_pct"] == pytest.approx(200.0 / 3.0)

    # Species rank: taxa 11, 12 (labeled) and 15 (unlabeled) -> 1 of 3 unlabeled
    assert by_rank["species"]["total_taxa"] == 3
    assert by_rank["species"]["labeled_taxa"] == 2
    assert by_rank["species"]["unlabeled_taxa"] == 1

    # Check filter cascade
    cascade_by_step = {(row["dimension"], row["step_id"]): row for row in cascade_rows}
    trait_ret = cascade_by_step[("traits", "retained_non_mediadive_traits")]
    assert trait_ret["retained_count"] == 8
    assert trait_ret["dropped_from_previous"] == 4  # 12 catalog - 8 retained
    assert trait_ret["matrix_dropped_from_previous"] == 1  # 1 MediaDive column in matrix
    assert trait_ret["metpo_phenotype_traits"] == 4
    assert trait_ret["chebi_carbon_source_traits"] == 1
    assert trait_ret["pathway_traits"] == 1
    assert trait_ret["go_process_traits"] == 2

    step_full = cascade_by_step[("studies", "full_human_16s_wms")]
    assert step_full["retained_count"] == 3  # bsdb:1, bsdb:2, bsdb:3
    assert step_full["signatures_count"] == 7  # 3 exps * 2 + 1 unpaired

    step_paired = cascade_by_step[("studies", "paired_human_16s_wms")]
    assert step_paired["retained_count"] == 2  # bsdb:1, bsdb:2
    assert step_paired["dropped_from_previous"] == 1  # bsdb:3 dropped
    assert step_paired["experiments_count"] == 3  # bsdb:1-1, bsdb:2-1, bsdb:2-2
    assert step_paired["signatures_count"] == 6
    assert step_paired["dropped_signatures_from_previous"] == 1

    step_tier12 = cascade_by_step[("studies", "tier12_subject_controls")]
    assert step_tier12["retained_count"] == 2  # bsdb:1, bsdb:2
    assert step_tier12["experiments_count"] == 2  # bsdb:1-1, bsdb:2-1
    assert step_tier12["signatures_count"] == 4
    assert step_tier12["dropped_signatures_from_previous"] == 2
    assert step_tier12["conditions_count"] == 1
    assert step_tier12["stool_studies_count"] == 2
    assert summary["tier12_unique_taxa"] == 6


def test_parkinsons_stool_study_overlap_tier12_and_full_export() -> None:
    # 3 PD Feces studies in full export (bsdb:1, bsdb:2, bsdb:3); 2 in Tier 1+2 (bsdb:1, bsdb:2)
    source_experiments = pd.DataFrame(
        [
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1-1",
                "host_species": "Homo sapiens",
                "sequencing_type": "16S",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
            },
            {
                "study_id": "bsdb:2",
                "experiment_id": "bsdb:2-1",
                "host_species": "Homo sapiens",
                "sequencing_type": "WMS",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
            },
            {
                "study_id": "bsdb:3",
                "experiment_id": "bsdb:3-1",
                "host_species": "Homo sapiens",
                "sequencing_type": "16S",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
            },
        ]
    )
    source_taxa = pd.DataFrame(
        [
            # bsdb:1 has {10, 11}
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1-1",
                "direction": "increased",
                "taxon_id": "NCBITaxon:10",
                "taxon_name": "GenusA",
                "taxonomic_rank": "genus",
            },
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1-1",
                "direction": "decreased",
                "taxon_id": "NCBITaxon:11",
                "taxon_name": "GenusA sp1",
                "taxonomic_rank": "species",
            },
            # bsdb:2 has {10, 12} -> Jaccard(1, 2) = 1/3
            {
                "study_id": "bsdb:2",
                "experiment_id": "bsdb:2-1",
                "direction": "increased",
                "taxon_id": "NCBITaxon:10",
                "taxon_name": "GenusA",
                "taxonomic_rank": "genus",
            },
            {
                "study_id": "bsdb:2",
                "experiment_id": "bsdb:2-1",
                "direction": "decreased",
                "taxon_id": "NCBITaxon:12",
                "taxon_name": "GenusB sp2",
                "taxonomic_rank": "species",
            },
            # bsdb:3 has {12} -> Jaccard(1, 3) = 0/3 = 0; Jaccard(2, 3) = 1/2 = 0.5
            {
                "study_id": "bsdb:3",
                "experiment_id": "bsdb:3-1",
                "direction": "increased",
                "taxon_id": "NCBITaxon:12",
                "taxon_name": "GenusB sp2",
                "taxonomic_rank": "species",
            },
        ]
    )
    tier12_experiments = source_experiments[
        source_experiments["experiment_id"].isin(["bsdb:1-1", "bsdb:2-1"])
    ].copy()
    tier12_taxa = source_taxa[source_taxa["experiment_id"].isin(["bsdb:1-1", "bsdb:2-1"])].copy()
    tier12_profile_keys = pd.DataFrame(
        [
            {"study_id": "bsdb:1", "experiment_id": "bsdb:1-1", "inc_dec": "inc"},
            {"study_id": "bsdb:1", "experiment_id": "bsdb:1-1", "inc_dec": "dec"},
            {"study_id": "bsdb:2", "experiment_id": "bsdb:2-1", "inc_dec": "inc"},
            {"study_id": "bsdb:2", "experiment_id": "bsdb:2-1", "inc_dec": "dec"},
        ]
    )

    summary_rows, pair_rows = compute_parkinsons_stool_overlap(
        tier12_experiments=tier12_experiments,
        tier12_taxa=tier12_taxa,
        tier12_profile_keys=tier12_profile_keys,
        source_experiments=source_experiments,
        source_taxa=source_taxa,
    )
    by_scope = {row["scope"]: row for row in summary_rows}
    assert by_scope["tier12_paired"]["n_studies"] == 2
    assert by_scope["tier12_paired"]["n_study_pairs"] == 1
    assert by_scope["tier12_paired"]["mean_jaccard"] == pytest.approx(1.0 / 3.0)
    assert by_scope["tier12_paired"]["median_jaccard"] == pytest.approx(1.0 / 3.0)

    assert by_scope["full_export"]["n_studies"] == 3
    assert by_scope["full_export"]["n_study_pairs"] == 3
    # Jaccards are 1/3, 0, 1/2 -> mean = 5/18, median = 1/3
    assert by_scope["full_export"]["mean_jaccard"] == pytest.approx(5.0 / 18.0)
    assert by_scope["full_export"]["median_jaccard"] == pytest.approx(1.0 / 3.0)
    assert len(pair_rows) == 4


def _build_subset_fixture_inputs(root: Path) -> dict[str, Path]:
    source = REPO_ROOT / "tests" / "fixtures" / "workflow"
    root.mkdir(parents=True, exist_ok=True)
    edges = (source / "merged-kg_edges.tsv").read_text(encoding="utf-8")
    nodes = (source / "merged-kg_nodes.tsv").read_text(encoding="utf-8")
    nodes += "METPO:2000103\tbiolink:OntologyClass\tcapable of\n"
    for idx, term in enumerate(DEFAULT_MATRIX_OVERLAP_TRAITS):
        if term == "METPO:1000608":
            continue
        if "|" in term:
            predicate, obj = term.split("|", 1)
        else:
            predicate, obj = "biolink:has_phenotype", term
        edges += f"extra{idx}\tNCBITaxon:101\t{predicate}\t{obj}\tfixture:annotation\n"
        nodes += f"{obj}\tbiolink:PhenotypicQuality\tFixture trait {idx}\n"
    paths = {
        "raw_dump": root / "full_dump.csv",
        "kg_edges": root / "merged-kg_edges.tsv",
        "kg_nodes": root / "merged-kg_nodes.tsv",
    }
    paths["kg_edges"].write_text(edges, encoding="utf-8")
    paths["kg_nodes"].write_text(nodes, encoding="utf-8")

    with (source / "full_dump.csv").open(encoding="utf-8") as handle:
        reader = csv.DictReader(line for line in handle if not line.startswith("#"))
        template = next(reader)
    for field in (
        "PMID",
        "DOI",
        "Study design",
        "Location of subjects",
        "EFO ID",
        "UBERON ID",
        "Group 1 definition",
        "Group 0 sample size",
        "Group 1 sample size",
        "Description",
    ):
        template[field] = ""
    with paths["raw_dump"].open("w", encoding="utf-8", newline="") as handle:
        handle.write("# BugSigDB fixture-1053, License: CC0 synthetic test data\n")
        writer = csv.DictWriter(handle, fieldnames=list(template))
        writer.writeheader()
        for study in (1, 2, 3):
            for exp_idx, direction in enumerate(("increased", "decreased"), start=1):
                if study == 3 and direction == "decreased":
                    # Study 3 is unpaired in the full export -> dropped at paired step
                    continue
                row = dict(template)
                row.update(
                    {
                        "BSDB ID": f"bsdb:{study}/1/{exp_idx}",
                        "Study": f"Study {study}",
                        "Condition": "Parkinson's disease",
                        "Host species": "Homo sapiens",
                        "Body site": "Feces",
                        "Sequencing type": "16S",
                        "Group 0 name": "healthy controls",
                        "Group 1 name": "Parkinson's disease",
                        "Abundance in Group 1": direction,
                        "NCBI Taxonomy IDs": "101;100" if direction == "increased" else "102",
                        "MetaPhlAn taxon names": (
                            "s__Examplegenus_alpha,g__Examplegenus"
                            if direction == "increased"
                            else "s__Examplegenus_beta"
                        ),
                    }
                )
                writer.writerow(row)
    return paths


def test_end_to_end_core_subset_consumer_provenance_and_temp_guard(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", "/dev/null")
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    from mtp.workflow.consumers import verify_consumer_record
    from mtp.workflow.operations import OPERATIONS, validate_core_product
    from mtp.workflow.runtime import execute_workflow

    inputs = _build_subset_fixture_inputs(tmp_path / "inputs")
    core = execute_workflow(
        REPO_ROOT / "workflows" / "bugsigdb_core_subset.yaml",
        inputs,
        tmp_path / "core",
        OPERATIONS,
        parameter_overrides={
            "source_mode": "fixture",
            "strain_rollup": "phenotype",
            "contrast_population": SUBSET_POPULATION,
        },
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    manifest = Path(core.manifest_path)
    output_dir = tmp_path / "matrix_overlap_coverage"

    # Without --allow-temporary-output, writing under tmp_path must be refused (#1049)
    with pytest.raises(ValueError, match="temporary directory"):
        run_matrix_overlap_coverage(manifest, output_dir)

    assert (
        main(
            [
                "--core-run",
                str(manifest),
                "--output-dir",
                str(output_dir),
                "--allow-temporary-output",
            ]
        )
        == 0
    )

    binding = verify_consumer_record(output_dir / "core_consumption.json")
    assert {
        "count_profiles",
        "experiments",
        "reported_taxa",
        "taxa_matrix",
        "feature_metadata",
        "feature_catalog",
        "source_experiments",
        "source_reported_taxa",
        "contrast_selection",
    } <= set(binding.product_ids)

    provenance_text = (output_dir / "provenance.json").read_text(encoding="utf-8")
    assert "/Users/" not in provenance_text
    assert "/private/" not in provenance_text
    provenance = json.loads(provenance_text)
    assert provenance["core_run_id"] == core.run_id
    assert provenance["contrast_population"] == SUBSET_POPULATION
    assert provenance["strain_rollup"] == "phenotype"
    for name, entry in provenance["outputs"].items():
        assert entry["sha256"] == file_sha256(output_dir / name)

    methods_note = REPO_ROOT / "methods" / "matrix_coverage_and_pd_overlap.md"
    assert methods_note.is_file()
    note_text = methods_note.read_text(encoding="utf-8")
    assert "taxa-coverage" in note_text
    assert "create_taxa_coverage_histogram.py" in note_text


def test_declared_trait_list_validation_and_file_loading(tmp_path: Path) -> None:
    matrix, columns = _synthetic_matrix_and_columns()
    trait_file = tmp_path / "traits.txt"
    trait_file.write_text(f"# comment\n\n{FAC_AER}\n{FAC_ANA}\n", encoding="utf-8")
    loaded = load_declared_traits(trait_list_path=trait_file)
    assert loaded.trait_ids == (FAC_AER, FAC_ANA)
    assert loaded.source_name == "traits.txt"

    with pytest.raises(ValueError, match="Specify at most one of"):
        load_declared_traits(traits=[FAC_AER, FAC_ANA], trait_list_path=trait_file)

    with pytest.raises(ValueError, match="At least two trait IDs"):
        load_declared_traits(traits=[FAC_AER])

    with pytest.raises(ValueError, match="Duplicate"):
        load_declared_traits(traits=[FAC_AER, FAC_AER])

    with pytest.raises(ValueError, match="MediaDive"):
        compute_matrix_trait_overlap(matrix, columns, [FAC_AER, MEDIADIVE])

    with pytest.raises(ValueError, match=r"absent from taxa_trait_matrix\.tsv columns"):
        compute_matrix_trait_overlap(matrix, columns, [FAC_AER, "METPO:9999999"])

    bad_matrix = matrix.copy()
    bad_matrix[FAC_AER] = [-1, 0, 0, 0, 0, 0]
    with pytest.raises(ValueError, match="negative or non-numeric"):
        compute_matrix_trait_overlap(bad_matrix, columns, [FAC_AER, FAC_ANA])


def test_validate_durable_output_dir_refuses_temporary_paths(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="temporary directory"):
        validate_durable_output_dir(tmp_path / "out")

    resolved = validate_durable_output_dir(tmp_path / "out", allow_temporary_output=True)
    assert resolved.name == "out"
    durable = validate_durable_output_dir(REPO_ROOT / "derived" / "matrix_overlap_coverage")
    assert durable.name == "matrix_overlap_coverage"
