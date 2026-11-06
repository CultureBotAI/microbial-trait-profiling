import csv
import importlib.util
import math
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

LEGACY_EDGE_COLUMNS = [
    "source",
    "target",
    "source_layer",
    "target_layer",
    "edge_class",
    "evidence_type",
    "direction",
    "effect",
    "p_value",
    "q_value",
    "study_support",
    "annotation_count",
    "statistically_significant",
    "line_style",
]
EDGE_METADATA_COLUMNS = ["p_value_type", "n_tests", "bh_family_size", "displayed_p_value_count"]


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
            "p_value_type": ["nominal"],
            "n_tests": [1],
            "bh_family_size": [1],
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
            "p_value_type": ["nominal"],
            "n_tests": [1],
            "bh_family_size": [1],
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
    assert annotation["p_value_type"] == "not_inferential"
    assert annotation[["p_value", "q_value", *EDGE_METADATA_COLUMNS[1:]]].isna().all()
    assert pd.isna(annotation["q_value"])
    assert "species" in set(nodes["layer"])
    assert selected_taxa["taxon_id"].tolist() == ["taxon-1"]


def test_build_network_tables_handles_empty_fdr_selection(tmp_path: Path) -> None:
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
            "p_value_type": ["nominal"],
            "n_tests": [1],
            "bh_family_size": [1],
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
            "p_value_type": ["nominal"],
            "n_tests": [1],
            "bh_family_size": [1],
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
    assert list(edges.columns) == LEGACY_EDGE_COLUMNS + EDGE_METADATA_COLUMNS
    assert selected_taxa.empty
    output = tmp_path / "empty_edges.tsv"
    network.write_tsv(edges, output)
    restored = pd.read_csv(output, sep="\t")
    assert restored.empty
    assert list(restored.columns) == LEGACY_EDGE_COLUMNS + EDGE_METADATA_COLUMNS


def _edge_export_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Keep source finite counts, padded BH sizes, and displayed counts distinct (#1430)."""
    taxon_tests = pd.DataFrame(
        {
            "condition": [f"Condition {letter}" for letter in "ABCDEF"],
            "taxon_id": ["taxon-1", "taxon-1", "taxon-2", "taxon-3", "taxon-4", "taxon-5"],
            "taxon_name": [
                "Taxon one",
                "Taxon one",
                "Taxon two",
                "Taxon three",
                "Taxon four",
                "Taxon five",
            ],
            "taxonomic_rank": ["species"] * 6,
            "direction": ["increased"] * 6,
            "directional_difference_fraction": [0.8, 0.7, 0.2, 0.0, 0.0, 0.1],
            "p_value": [0.001, 0.002, 0.2, float("nan"), float("nan"), 0.8],
            "q_value": [0.005, 0.005, 1 / 3, float("nan"), float("nan"), float("nan")],
            "reported_study_count": [8, 8, 8, 8, 8, 2],
            "test": ["paired taxon test"] * 6,
            "selected_after_bh": [True, True, False, False, False, False],
            "passes_min_taxon_studies": [True, True, True, True, True, False],
            "p_value_type": [
                "nominal",
                "nominal",
                "nominal",
                "not_tested",
                "not_tested",
                "nominal",
            ],
            "n_tests": [4] * 6,
            "bh_family_size": pd.array([5, 5, 5, 5, 5, pd.NA], dtype="Int64"),
        }
    )
    trait_tests = pd.DataFrame(
        {
            "condition": [f"Condition {letter}" for letter in "ABCDE"],
            "trait_id": ["trait-1"] * 5,
            "trait": ["Trait one"] * 5,
            "direction": ["increased"] * 5,
            "mean_contrast": [0.7, 0.2, 0.0, 0.0, 0.1],
            "p_value": [0.001, 0.2, float("nan"), float("nan"), 0.8],
            "q_value": [0.004, 0.4, float("nan"), float("nan"), float("nan")],
            "paired_study_count": [8, 8, 8, 8, 2],
            "test": ["paired trait test"] * 5,
            "selected_after_bh": [True, False, False, False, False],
            "passes_min_paired_studies": [True, True, True, True, False],
            "p_value_type": ["nominal", "nominal", "not_tested", "not_tested", "nominal"],
            "n_tests": [3] * 5,
            "bh_family_size": pd.array([4, 4, 4, 4, pd.NA], dtype="Int64"),
        }
    )
    annotations = pd.DataFrame(
        {
            "taxon_id": ["taxon-1"],
            "taxon_name": ["Taxon one"],
            "taxonomic_rank": ["species"],
            "trait-1": [3],
            "trait-2": [2],
        }
    )
    traits = pd.DataFrame(
        {"trait_id": ["trait-1", "trait-2"], "trait_label": ["Trait one", "Trait two"]}
    )
    return taxon_tests, trait_tests, annotations, traits


@pytest.mark.parametrize("kind", ["taxon", "trait"])
@pytest.mark.parametrize("selected_untestable", [False, True])
def test_network_edges_preserve_source_statistic_and_separate_test_counts(
    kind: str, selected_untestable: bool
) -> None:
    taxon_tests, trait_tests, annotations, traits = _edge_export_inputs()
    if selected_untestable:
        # The exporter preserves caller selection, even when an upstream test is unavailable.
        # A displayed-test count must count finite p entries, not every exported row.
        taxon_tests.loc[3, "selected_after_bh"] = True
        trait_tests.loc[2, "selected_after_bh"] = True
    original_taxa, original_traits = taxon_tests.copy(deep=True), trait_tests.copy(deep=True)
    _, edges, selected_taxa = network.build_network_tables(
        taxon_tests, trait_tests, annotations, traits
    )
    source = taxon_tests if kind == "taxon" else trait_tests
    selected = source.loc[source["selected_after_bh"]]
    edge_class = f"FDR-significant condition-{kind} " + (
        "direction" if kind == "taxon" else "contrast"
    )
    actual = edges.loc[edges["edge_class"].eq(edge_class)].set_index("source")
    expected_displayed = sum(math.isfinite(float(value)) for value in selected["p_value"])
    expected_source_tests = sum(math.isfinite(float(value)) for value in source["p_value"])
    support = "passes_min_taxon_studies" if kind == "taxon" else "passes_min_paired_studies"
    expected_bh_size = int(source[support].sum())
    assert len({expected_displayed, expected_source_tests, expected_bh_size}) == 3
    assert actual["displayed_p_value_count"].tolist() == [expected_displayed] * len(selected)
    assert actual["n_tests"].tolist() == [expected_source_tests] * len(selected)
    assert actual["bh_family_size"].tolist() == [expected_bh_size] * len(selected)
    assert list(edges.columns) == LEGACY_EDGE_COLUMNS + EDGE_METADATA_COLUMNS
    assert actual.index.tolist() == [f"condition::{value}" for value in selected["condition"]]
    effect_column = "directional_difference_fraction" if kind == "taxon" else "mean_contrast"
    support_column = "reported_study_count" if kind == "taxon" else "paired_study_count"
    for _, row in selected.iterrows():
        edge = actual.loc[f"condition::{row['condition']}"]
        assert edge["target"] == f"{kind}::{row[f'{kind}_id']}"
        assert edge["p_value_type"] == row["p_value_type"]
        assert edge["n_tests"] == row["n_tests"]
        assert edge["bh_family_size"] == row["bh_family_size"]
        for field in ("p_value", "q_value"):
            assert pd.isna(edge[field]) if pd.isna(row[field]) else edge[field] == row[field]
        assert edge["effect"] == row[effect_column]
        assert edge["study_support"] == row[support_column]
        assert edge["evidence_type"] == row["test"]
        assert edge["direction"] == row["direction"]
        assert edge["line_style"] == "solid"
        assert bool(edge["statistically_significant"])
    pd.testing.assert_frame_equal(selected_taxa, original_taxa[original_taxa["selected_after_bh"]])
    pd.testing.assert_frame_equal(taxon_tests, original_taxa)
    pd.testing.assert_frame_equal(trait_tests, original_traits)


def test_network_edge_metadata_roundtrip_keeps_annotations_noninferential(tmp_path: Path) -> None:
    taxon_tests, trait_tests, annotations, traits = _edge_export_inputs()
    _, edges, _ = network.build_network_tables(taxon_tests, trait_tests, annotations, traits)
    path = tmp_path / "edges.tsv"
    network.write_tsv(edges, path)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == LEGACY_EDGE_COLUMNS + EDGE_METADATA_COLUMNS
        rows = list(reader)
    annotation_rows = [
        row for row in rows if row["edge_class"] == "KG-Microbe annotation provenance"
    ]
    assert len(annotation_rows) == 2
    for row in annotation_rows:
        assert row["p_value_type"] == "not_inferential"
        assert all(
            row[field] == "NA" for field in ("p_value", "q_value", *EDGE_METADATA_COLUMNS[1:])
        )
        assert row["statistically_significant"] == "False"
        assert row["line_style"] == "dashed"
    inferential = [row for row in rows if row not in annotation_rows]
    assert len(inferential) == 3
    assert all(row["p_value_type"] == "nominal" for row in inferential)
    for row in inferential:
        kind = "taxon" if "condition-taxon" in row["edge_class"] else "trait"
        expected = {"taxon": (4, 5, 2), "trait": (3, 4, 1)}[kind]
        assert tuple(float(row[field]) for field in EDGE_METADATA_COLUMNS[1:]) == expected


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
    for kind, source in (("taxon", taxon_tests), ("trait", trait_tests)):
        inferential = edges[edges["edge_class"].str.startswith(f"FDR-significant condition-{kind}")]
        selected = source[source["selected_after_bh"]]
        assert inferential["p_value_type"].tolist() == selected["p_value_type"].tolist()
        assert inferential["n_tests"].tolist() == selected["n_tests"].tolist()
        assert inferential["bh_family_size"].tolist() == selected["bh_family_size"].tolist()
        assert inferential["displayed_p_value_count"].tolist() == [1]
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
    for fragment in (
        "`displayed_p_value_count` counts finite p-values in exported edges separately by "
        "inferential edge class",
        "`n_tests` preserves the source table count (including finite tests outside support "
        "eligibility)",
        "`bh_family_size` preserves its support-eligible padded BH family",
        "Annotation edges have no inferential counts.",
    ):
        assert fragment in methods_text
    assert (
        "Reported `p_value` columns and Panel C `nominal p` labels are nominal "
        "(uncorrected) across the `1` tested (`2` BH-family) condition-taxon tests and "
        "`1` tested (`2` BH-family) condition-trait tests"
    ) in methods_text
    assert (
        "`2` eligible condition-trait comparisons (`1` with finite nominal p; "
        "`1` PD trait comparisons shown in Panel C)."
    ) in methods_text
