import importlib.util
import sys
from pathlib import Path

import pandas as pd

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


def test_literature_evidence_keeps_citations_separate_and_preserves_conflict() -> None:
    sources = pd.DataFrame(
        [
            ["A", "Earlier study", 2020, "cohort", "10.1000/earlier", "https://doi.org/a"],
            ["B", "Later study", 2022, "meta-analysis", "10.1000/later", "https://doi.org/b"],
        ],
        columns=["source_key", "citation", "year", "study_design", "doi", "url"],
    )
    evidence = pd.DataFrame(
        [
            [
                "NCBITaxon:1",
                "Current name",
                "species",
                "Older name",
                "dec",
                "A",
                "cohort",
                "true",
                "decreased",
            ],
            [
                "NCBITaxon:1",
                "Current name",
                "species",
                "Current name",
                "inc",
                "B",
                "meta-analysis",
                "true",
                "increased",
            ],
            [
                "NCBITaxon:2",
                "Context genus",
                "genus",
                "Context genus",
                "inc",
                "B",
                "meta-analysis",
                "false",
                "not projected to species",
            ],
        ],
        columns=[
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
            "source_reported_name",
            "reported_pd_direction",
            "source_key",
            "evidence_scope",
            "heatmap_eligible",
            "note",
        ],
    )

    result = heatmap.aggregate_literature_evidence(evidence, sources)

    assert result.index.tolist() == ["NCBITaxon:1"]
    assert result.loc["NCBITaxon:1", "literature_direction"] == "mixed"
    assert result.loc["NCBITaxon:1", "literature_sources"] == "A,B"
    assert result.loc["NCBITaxon:1", "literature_dois"] == ("10.1000/earlier;10.1000/later")
    assert result.loc["NCBITaxon:1", "literature_reported_names"] == ("Older name | Current name")
