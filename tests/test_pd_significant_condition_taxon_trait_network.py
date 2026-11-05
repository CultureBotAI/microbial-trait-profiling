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
    assert tests.loc["genus-1", "q_value"] == pytest.approx(0.015625)
    assert tests.loc["species-1", "both_directions_studies"] == 3
    assert not bool(tests.loc["species-1", "selected_after_bh"])
    assert pd.isna(tests.loc["species-1", "p_value"])
    assert pd.isna(tests.loc["species-1", "q_value"])
    assert tests.loc["genus-1", "p_value_type"] == "nominal"
    assert tests.loc["species-1", "p_value_type"] == "not_tested"
    assert tests["n_tests"].tolist() == [1, 1]
    assert tests["bh_family_size"].tolist() == [2, 2]


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


def test_network_tables_figure_and_methods_label_nominal_p_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Taxon/trait tables, Panel C labels, and methods summary label p as nominal (#1335)."""
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mplconfig"))
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
                    "experiment_id": f"experiment-{index}-{direction}",
                    "direction": "increased" if direction == "inc" else "decreased",
                    "inc_dec": direction,
                    "condition": network.PD_CONDITION,
                    "body_site_group": network.FECES,
                    "host_species": network.HUMAN_HOST,
                }
            )
            zscore_rows.append({"trait-1": value})
    for index in range(2):
        for direction, value in [("inc", 0.5), ("dec", 0.5)]:
            profile_id = f"other-study-{index}|experiment-{index}|{direction}"
            indices.append(profile_id)
            metadata_rows.append(
                {
                    "study_id": f"other-study-{index}",
                    "experiment_id": f"experiment-{index}-{direction}",
                    "direction": "increased" if direction == "inc" else "decreased",
                    "inc_dec": direction,
                    "condition": "Condition B",
                    "body_site_group": network.FECES,
                    "host_species": network.HUMAN_HOST,
                }
            )
            zscore_rows.append({"trait-1": value})
    for index in range(3):
        for direction, value in [("inc", 0.5), ("dec", 0.5)]:
            profile_id = f"tied-study-{index}|experiment-{index}|{direction}"
            indices.append(profile_id)
            metadata_rows.append(
                {
                    "study_id": f"tied-study-{index}",
                    "experiment_id": f"experiment-{index}-{direction}",
                    "direction": "increased" if direction == "inc" else "decreased",
                    "inc_dec": direction,
                    "condition": "Condition C",
                    "body_site_group": network.FECES,
                    "host_species": network.HUMAN_HOST,
                }
            )
            zscore_rows.append({"trait-1": value})
    metadata = pd.DataFrame(metadata_rows, index=indices)
    zscores = pd.DataFrame(zscore_rows, index=indices)
    traits = pd.DataFrame(
        {
            "trait_id": ["trait-1"],
            "trait_label": ["cellulose catabolic process"],
            "trait_order": [0],
        }
    )
    eligibility = pd.DataFrame(
        {
            "condition": [network.PD_CONDITION, "Condition B", "Condition C"],
            "paired_study_count": [8, 2, 3],
            "passes_min_paired_studies": [True, True, True],
            "displayed_in_network": [True, False, False],
        }
    )
    annotated_taxa = pd.DataFrame(
        {
            "taxon_id": ["genus-1", "species-low", "species-tie"],
            "taxon_name": ["Roseburia", "Low support species", "Tied species"],
            "taxonomic_rank": ["genus", "species", "species"],
            "trait-1": [1, 1, 1],
        }
    )
    taxon_rows = (
        [
            {
                "study_id": f"study-{index}",
                "experiment_id": f"experiment-{index}-inc",
                "direction": "increased",
                "taxon_id": "genus-1",
            }
            for index in range(8)
        ]
        + [
            {
                "study_id": f"study-{index}",
                "experiment_id": f"experiment-{index}-inc",
                "direction": "increased",
                "taxon_id": "species-low",
            }
            for index in range(2)
        ]
        + [
            {
                "study_id": f"study-{index}",
                "experiment_id": f"experiment-{index}-{short_dir}",
                "direction": full_dir,
                "taxon_id": "species-tie",
            }
            for index in range(3)
            for short_dir, full_dir in [("inc", "increased"), ("dec", "decreased")]
        ]
    )
    monkeypatch.setattr(network, "read_tsv", lambda *_args, **_kwargs: pd.DataFrame(taxon_rows))

    taxon_tests = network.build_condition_taxon_tests(metadata, eligibility, annotated_taxa)
    trait_tests, _ = network.build_condition_trait_tests(metadata, zscores, traits, eligibility)

    taxon_by_id = taxon_tests.set_index("taxon_id")
    assert taxon_by_id.loc["genus-1", "p_value_type"] == "nominal"
    assert taxon_by_id.loc["genus-1", "q_value"] == pytest.approx(0.015625)
    assert taxon_by_id.loc["species-low", "p_value_type"] == "nominal"
    assert taxon_by_id.loc["species-tie", "p_value_type"] == "not_tested"
    assert pd.isna(taxon_by_id.loc["species-tie", "p_value"])
    assert pd.isna(taxon_by_id.loc["species-tie", "q_value"])
    assert len(taxon_tests) == 3
    assert taxon_tests["n_tests"].tolist() == [2, 2, 2]
    assert taxon_tests.loc[taxon_tests["passes_min_taxon_studies"], "bh_family_size"].tolist() == [
        2,
        2,
    ]
    assert taxon_tests.loc[~taxon_tests["passes_min_taxon_studies"], "bh_family_size"].isna().all()

    trait_by_cond = trait_tests.set_index("condition")
    assert trait_by_cond.loc[network.PD_CONDITION, "p_value_type"] == "nominal"
    assert trait_by_cond.loc[network.PD_CONDITION, "q_value"] == pytest.approx(0.015625)
    assert trait_by_cond.loc["Condition B", "p_value_type"] == "not_tested"
    assert pd.isna(trait_by_cond.loc["Condition B", "p_value"])
    assert pd.isna(trait_by_cond.loc["Condition B", "q_value"])
    assert trait_by_cond.loc["Condition C", "p_value_type"] == "not_tested"
    assert pd.isna(trait_by_cond.loc["Condition C", "p_value"])
    assert pd.isna(trait_by_cond.loc["Condition C", "q_value"])
    assert len(trait_tests) == 3
    assert trait_tests["n_tests"].tolist() == [1, 1, 1]
    assert trait_tests.loc[trait_tests["passes_min_paired_studies"], "bh_family_size"].tolist() == [
        2,
        2,
    ]
    assert trait_tests.loc[~trait_tests["passes_min_paired_studies"], "bh_family_size"].isna().all()

    nodes, edges, _ = network.build_network_tables(taxon_tests, trait_tests, annotated_taxa, traits)
    import matplotlib.figure

    captured_texts: list[str] = []

    def capture_savefig(
        self: matplotlib.figure.Figure, path: Path, *_args: object, **_kwargs: object
    ) -> None:
        Path(path).write_text("", encoding="utf-8")
        for ax in self.axes:
            captured_texts.extend(t.get_text() for t in ax.texts)
            captured_texts.append(ax.get_title(loc="left"))

    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", capture_savefig)
    fig_path = tmp_path / "landscape.png"
    network.plot_fdr_evidence_landscape(taxon_tests, trait_tests, traits, fig_path)
    joined_texts = "\n".join(captured_texts)
    assert "nominal p=" in joined_texts
    assert "nominal p is unadjusted (1 PD trait tests shown)" in joined_texts
    assert "q is BH-adjusted across 2 condition-trait tests (1 with finite p)." in joined_texts

    methods_path = tmp_path / "methods.md"
    network.write_methods_summary(eligibility, taxon_tests, trait_tests, nodes, edges, methods_path)
    methods_text = methods_path.read_text(encoding="utf-8")
    assert (
        "Reported `p_value` columns and Panel C `nominal p` labels are nominal "
        "(uncorrected) across the `1` tested (`2` BH-family) condition-taxon tests and "
        "`1` tested (`2` BH-family) condition-trait tests"
    ) in methods_text
    assert (
        "`2` eligible condition-trait comparisons (`1` with finite nominal p; "
        "`1` PD trait comparisons shown in Panel C)."
    ) in methods_text
