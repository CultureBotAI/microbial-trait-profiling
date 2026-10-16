import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "create_pd_significant_condition_taxon_trait_network.py"
)
SPEC = importlib.util.spec_from_file_location("pd_significant_network", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
network = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = network
SPEC.loader.exec_module(network)


def test_benjamini_hochberg_adjusts_complete_family() -> None:
    adjusted = network.benjamini_hochberg(pd.Series([0.01, 0.04, 0.03, float("nan")]))

    assert adjusted.iloc[:3].tolist() == pytest.approx([0.03, 0.04, 0.04])
    assert pd.isna(adjusted.iloc[3])


def test_exact_mcnemar_uses_only_discordant_pairs() -> None:
    assert network.exact_mcnemar_p_value(8, 0) == pytest.approx(0.0078125)
    assert network.exact_mcnemar_p_value(4, 4) == pytest.approx(1.0)
    assert network.exact_mcnemar_p_value(0, 0) == pytest.approx(1.0)


def test_evidence_tiers_preserve_the_significance_boundary() -> None:
    assert network.evidence_tier(0.05) == "FDR-significant"
    assert network.evidence_tier(0.05001) == "exploratory"
    assert network.evidence_tier(0.20) == "exploratory"
    assert network.evidence_tier(0.20001) == "other tested"
    assert network.evidence_tier(float("nan")) == "ineligible"


def test_condition_taxon_tests_keep_ranks_separate_and_apply_bh(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    metadata_rows = []
    taxon_rows = []
    for index in range(10):
        study_id = f"study-{index}"
        for short_direction, direction in [("inc", "increased"), ("dec", "decreased")]:
            metadata_rows.append(
                {
                    "study_id": study_id,
                    "experiment_id": f"experiment-{index}-{short_direction}",
                    "direction": direction,
                    "inc_dec": short_direction,
                    "condition": "Condition A",
                }
            )
        if index < 8:
            taxon_rows.append(
                {
                    "study_id": study_id,
                    "experiment_id": f"experiment-{index}-inc",
                    "direction": "increased",
                    "taxon_id": "genus-1",
                }
            )
        if index < 3:
            taxon_rows.extend(
                [
                    {
                        "study_id": study_id,
                        "experiment_id": f"experiment-{index}-inc",
                        "direction": "increased",
                        "taxon_id": "species-1",
                    },
                    {
                        "study_id": study_id,
                        "experiment_id": f"experiment-{index}-dec",
                        "direction": "decreased",
                        "taxon_id": "species-1",
                    },
                ]
            )
    feces_metadata = pd.DataFrame(metadata_rows)
    eligibility = pd.DataFrame(
        {
            "condition": ["Condition A"],
            "paired_study_count": [10],
            "passes_min_paired_studies": [True],
        }
    )
    annotated_taxa = pd.DataFrame(
        {
            "taxon_id": ["genus-1", "species-1"],
            "taxon_name": ["Genus one", "Species one"],
            "taxonomic_rank": ["genus", "species"],
        }
    )
    monkeypatch.setattr(network, "read_tsv", lambda *_args, **_kwargs: pd.DataFrame(taxon_rows))

    tests = network.build_condition_taxon_tests(
        feces_metadata,
        eligibility,
        annotated_taxa,
    ).set_index("taxon_id")

    assert tests.loc["genus-1", "taxonomic_rank"] == "genus"
    assert tests.loc["species-1", "taxonomic_rank"] == "species"
    assert tests.loc["genus-1", "increased_only_studies"] == 8
    assert bool(tests.loc["genus-1", "selected_after_bh"])
    assert tests.loc["species-1", "both_directions_studies"] == 3
    assert not bool(tests.loc["species-1", "selected_after_bh"])


def test_condition_trait_tests_select_consistent_paired_contrasts() -> None:
    metadata_rows = []
    zscore_rows = []
    indices = []
    for index in range(8):
        for direction, value in [("inc", 1.0), ("dec", 0.0)]:
            profile_id = f"study-{index}|experiment-{index}|{direction}"
            indices.append(profile_id)
            metadata_rows.append(
                {
                    "study_id": f"study-{index}",
                    "inc_dec": direction,
                    "condition": "Condition A",
                    "body_site_group": network.FECES,
                    "host_species": network.HUMAN_HOST,
                }
            )
            zscore_rows.append({"trait-1": value})
    metadata = pd.DataFrame(metadata_rows, index=indices)
    zscores = pd.DataFrame(zscore_rows, index=indices)
    traits = pd.DataFrame({"trait_id": ["trait-1"], "trait_label": ["trait one"]})
    eligibility = pd.DataFrame(
        {
            "condition": ["Condition A"],
            "paired_study_count": [8],
            "passes_min_paired_studies": [True],
        }
    )

    tests, contrasts = network.build_condition_trait_tests(
        metadata,
        zscores,
        traits,
        eligibility,
    )

    assert len(contrasts) == 8
    assert tests.iloc[0]["mean_contrast"] == pytest.approx(1.0)
    assert tests.iloc[0]["q_value"] == pytest.approx(0.0078125)
    assert bool(tests.iloc[0]["selected_after_bh"])


def test_annotation_edges_are_explicitly_noninferential() -> None:
    taxon_tests = pd.DataFrame(
        {
            "condition": ["Condition A"],
            "taxon_id": ["taxon-1"],
            "taxon_name": ["Taxon one"],
            "taxonomic_rank": ["species"],
            "direction": ["increased"],
            "directional_difference_fraction": [0.8],
            "p_value": [0.01],
            "q_value": [0.02],
            "reported_study_count": [8],
            "test": ["paired taxon test"],
            "selected_after_bh": [True],
        }
    )
    trait_tests = pd.DataFrame(
        {
            "condition": ["Condition A"],
            "trait_id": ["trait-1"],
            "trait": ["trait one"],
            "direction": ["increased"],
            "mean_contrast": [0.7],
            "p_value": [0.01],
            "q_value": [0.02],
            "paired_study_count": [8],
            "test": ["paired trait test"],
            "selected_after_bh": [True],
        }
    )
    annotated_taxa = pd.DataFrame(
        {
            "taxon_id": ["taxon-1"],
            "taxon_name": ["Taxon one"],
            "taxonomic_rank": ["species"],
            "trait-1": [3],
        }
    )
    traits = pd.DataFrame(
        {"trait_id": ["trait-1"], "trait_label": ["trait one"], "trait_order": [0]}
    )

    nodes, edges, selected_taxa = network.build_network_tables(
        taxon_tests,
        trait_tests,
        annotated_taxa,
        traits,
    )

    annotation = edges[edges["edge_class"].eq("KG-Microbe annotation provenance")].iloc[0]
    assert annotation["line_style"] == "dashed"
    assert not bool(annotation["statistically_significant"])
    assert pd.isna(annotation["q_value"])
    assert "species" in set(nodes["layer"])
    assert selected_taxa["taxon_id"].tolist() == ["taxon-1"]


def test_build_network_tables_handles_empty_fdr_selection() -> None:
    taxon_tests = pd.DataFrame(
        {
            "condition": ["Condition A"],
            "taxon_id": ["taxon-1"],
            "taxon_name": ["Taxon one"],
            "taxonomic_rank": ["species"],
            "direction": ["increased"],
            "directional_difference_fraction": [0.4],
            "p_value": [0.02],
            "q_value": [0.06],
            "reported_study_count": [8],
            "test": ["paired taxon test"],
            "selected_after_bh": [False],
        }
    )
    trait_tests = pd.DataFrame(
        {
            "condition": ["Condition A"],
            "trait_id": ["trait-1"],
            "trait": ["trait one"],
            "direction": ["increased"],
            "mean_contrast": [0.5],
            "p_value": [0.00017],
            "q_value": [0.056],
            "paired_study_count": [45],
            "test": ["paired trait test"],
            "selected_after_bh": [False],
        }
    )
    annotated_taxa = pd.DataFrame(
        {
            "taxon_id": ["taxon-1"],
            "taxon_name": ["Taxon one"],
            "taxonomic_rank": ["species"],
            "trait-1": [3],
        }
    )
    traits = pd.DataFrame(
        {"trait_id": ["trait-1"], "trait_label": ["trait one"], "trait_order": [0]}
    )

    nodes, edges, selected_taxa = network.build_network_tables(
        taxon_tests,
        trait_tests,
        annotated_taxa,
        traits,
    )

    assert nodes.empty
    assert list(nodes.columns) == network.NODE_COLUMNS
    assert edges.empty
    assert list(edges.columns) == network.EDGE_COLUMNS
    assert selected_taxa.empty
