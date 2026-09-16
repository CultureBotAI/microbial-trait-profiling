import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts/create_pd_six_trait_taxon_heatmap.py"
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("pd_six_trait_taxon_heatmap", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
heatmap = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(heatmap)


def test_leave_one_species_out_effect_can_differ_by_trait() -> None:
    rows = []
    metadata_rows = []
    for study in ["s1", "s2", "s3"]:
        for direction, counts in [("inc", [50, 5, 45]), ("dec", [5, 50, 45])]:
            identifier = f"{study}|e1|{direction}"
            rows.append([identifier, *counts])
            metadata_rows.append(
                [
                    identifier,
                    study,
                    "e1",
                    direction,
                    heatmap.PD_CONDITION,
                    heatmap.FECES,
                    "Homo sapiens",
                ]
            )
    profiles = pd.DataFrame(rows, columns=["profile_id", "t1", "t2", "other"]).set_index(
        "profile_id"
    )
    metadata = pd.DataFrame(
        metadata_rows,
        columns=[
            "profile_id",
            "study_id",
            "experiment_id",
            "inc_dec",
            "condition",
            "body_site_group",
            "host_species",
        ],
    ).set_index("profile_id")
    baseline, totals, _, standard_deviations = heatmap.fixed_reference(
        profiles, ["t1", "t2"], metadata.index
    )
    pairs = heatmap.complete_pd_pairs(metadata, profiles.index)
    taxon_vector = pd.Series({"t1": 5, "t2": 5, "other": 10})
    influence = heatmap.taxon_influence(
        profiles,
        totals,
        taxon_vector,
        profiles.index,
        ["t1", "t2"],
        standard_deviations,
    )

    trait_one = heatmap.paired_study_contrasts(influence["t1"], pairs)
    trait_two = heatmap.paired_study_contrasts(influence["t2"], pairs)
    one_mean, one_low, one_high = heatmap.mean_interval(trait_one)
    two_mean, two_low, two_high = heatmap.mean_interval(trait_two)

    assert baseline.shape == (6, 2)
    assert one_mean < 0 < two_mean
    assert heatmap.classify_effect(one_mean, one_low, one_high, 3) == "dec"
    assert heatmap.classify_effect(two_mean, two_low, two_high, 3) == "inc"


def test_classification_preserves_uncertainty_and_low_support() -> None:
    assert heatmap.classify_effect(0.4, 0.1, 0.7, 3) == "inc"
    assert heatmap.classify_effect(-0.4, -0.7, -0.1, 3) == "dec"
    assert heatmap.classify_effect(0.1, -0.2, 0.4, 3) == "mixed"
    assert heatmap.classify_effect(0.4, 0.1, 0.7, 2) == "low support"


def test_unannotated_trait_retains_denominator_mediated_influence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    profile_rows = []
    metadata_rows = []
    evidence_rows = []
    for study in ["s1", "s2", "s3"]:
        for direction, counts in [("inc", [50, 5, 45]), ("dec", [5, 50, 45])]:
            identifier = f"{study}|e1|{direction}"
            profile_rows.append([identifier, *counts])
            metadata_rows.append(
                [
                    identifier,
                    study,
                    "e1",
                    direction,
                    heatmap.PD_CONDITION,
                    heatmap.FECES,
                    "Homo sapiens",
                ]
            )
            evidence_rows.append(
                [
                    "t1",
                    "taxon-a",
                    "species",
                    identifier,
                    study,
                    {"inc": "increased", "dec": "decreased"}[direction],
                ]
            )

    profiles = pd.DataFrame(profile_rows, columns=["profile_id", "t1", "t2", "other"]).set_index(
        "profile_id"
    )
    metadata = pd.DataFrame(
        metadata_rows,
        columns=[
            "profile_id",
            "study_id",
            "experiment_id",
            "inc_dec",
            "condition",
            "body_site_group",
            "host_species",
        ],
    ).set_index("profile_id")
    taxon_vectors = pd.DataFrame(
        [["taxon-a", "Species A", "species", 5, 0, 10]],
        columns=["taxon_id", "taxon_name", "taxonomic_rank", "t1", "t2", "other"],
    ).set_index("taxon_id")
    evidence = pd.DataFrame(
        evidence_rows,
        columns=[
            "trait_id",
            "taxon_id",
            "taxonomic_rank",
            "profile_id",
            "study_id",
            "direction",
        ],
    )
    trait_summary = pd.DataFrame(
        [["t1", "Trait one", 1.0], ["t2", "Trait two", -1.0]],
        columns=["trait_id", "trait", "increased_minus_decreased_mean_z"],
    )
    monkeypatch.setattr(
        heatmap,
        "load_taxon_specificity",
        lambda taxon_ids, _metadata: pd.DataFrame({"taxon_id": taxon_ids}),
    )

    cells, _, _ = heatmap.build_influence_cells(
        profiles, metadata, taxon_vectors, evidence, trait_summary
    )
    by_trait = cells.set_index("trait_id")

    assert by_trait.loc["t1", "influence_basis"] == "direct + denominator"
    assert by_trait.loc["t2", "influence_basis"] == "denominator only"
    assert not bool(by_trait.loc["t2", "has_trait_annotation"])
    assert pd.notna(by_trait.loc["t2", "mean_influence_contrast_z"])
    assert pd.notna(by_trait.loc["t2", "ci95_low"])
    assert by_trait.loc["t2", "status"] != "not annotated"


def test_taxon_specificity_compares_pd_feces_to_other_feces() -> None:
    metadata = pd.DataFrame(
        [
            ["s1|e1|inc", "s1", "e1", "inc", heatmap.PD_CONDITION, heatmap.FECES],
            ["s1|e1|dec", "s1", "e1", "dec", heatmap.PD_CONDITION, heatmap.FECES],
            ["s2|e1|inc", "s2", "e1", "inc", "Colorectal cancer", heatmap.FECES],
            ["s3|e1|inc", "s3", "e1", "inc", heatmap.PD_CONDITION, "Oral cavity"],
        ],
        columns=[
            "profile_id",
            "study_id",
            "experiment_id",
            "inc_dec",
            "condition",
            "body_site_group",
        ],
    ).set_index("profile_id", drop=False)
    taxa = pd.DataFrame(
        [
            ["s1", "e1", "increased", "taxon-a"],
            ["s1", "e1", "decreased", "taxon-a"],
            ["s2", "e1", "increased", "taxon-a"],
            ["s3", "e1", "increased", "taxon-a"],
            ["s2", "e1", "increased", "taxon-b"],
        ],
        columns=["study_id", "experiment_id", "direction", "taxon_id"],
    )

    specificity = heatmap.taxon_specificity_from_occurrences(
        taxa,
        metadata,
        pd.Index(["taxon-a", "taxon-b"]),
    ).set_index("taxon_id")

    assert specificity.loc["taxon-a", "pd_other_feces_counts"] == "2/1"
    assert specificity.loc["taxon-a", "pd_feces_studies"] == 1
    assert specificity.loc["taxon-a", "other_feces_studies"] == 1
    assert specificity.loc["taxon-a", "pd_specificity_fraction"] == 2 / 3
    assert specificity.loc["taxon-b", "pd_other_feces_counts"] == "0/1"
