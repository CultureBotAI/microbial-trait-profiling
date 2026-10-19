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


def pd_pairs_fixture() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Three studies whose t1 is higher in increased profiles and t2 in decreased ones."""
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
    return profiles, metadata, taxon_vectors, evidence


def test_unannotated_trait_retains_denominator_mediated_influence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    profiles, metadata, taxon_vectors, evidence = pd_pairs_fixture()
    trait_summary = pd.DataFrame(
        [["t1", "Trait one", 1.0], ["t2", "Trait two", -1.0]],
        columns=["trait_id", "trait", "increased_minus_decreased_mean_z"],
    )
    monkeypatch.setattr(
        heatmap,
        "load_taxon_specificity",
        lambda taxon_ids, _metadata: pd.DataFrame({"taxon_id": taxon_ids}),
    )

    cells, _, _, removed_low_support_taxa = heatmap.build_influence_cells(
        profiles, metadata, taxon_vectors, evidence, trait_summary
    )
    by_trait = cells.set_index("trait_id")

    assert removed_low_support_taxa == 0
    assert by_trait.loc["t1", "influence_basis"] == "direct + denominator"
    assert by_trait.loc["t2", "influence_basis"] == "denominator only"
    assert not bool(by_trait.loc["t2", "has_trait_annotation"])
    assert pd.notna(by_trait.loc["t2", "mean_influence_contrast_z"])
    assert pd.notna(by_trait.loc["t2", "ci95_low"])
    assert by_trait.loc["t2", "status"] != "not annotated"


def test_traits_carry_the_lollipop_direction_as_data(monkeypatch: pytest.MonkeyPatch) -> None:
    """Direction comes from the lollipop summary even where this figure's baseline disagrees."""
    profiles, metadata, taxon_vectors, evidence = pd_pairs_fixture()
    trait_summary = pd.DataFrame(
        [["t1", "Trait one", -0.2], ["t2", "Trait two", 0.3]],
        columns=["trait_id", "trait", "increased_minus_decreased_mean_z"],
    )
    monkeypatch.setattr(
        heatmap,
        "load_taxon_specificity",
        lambda taxon_ids, _metadata: pd.DataFrame({"taxon_id": taxon_ids}),
    )

    cells, _, traits, _ = heatmap.build_influence_cells(
        profiles, metadata, taxon_vectors, evidence, trait_summary
    )

    assert traits["trait_id"].tolist() == ["t2", "t1"]
    by_trait = traits.set_index("trait_id")
    assert by_trait.loc["t1", "trait_direction"] == "decreased"
    assert by_trait.loc["t1", "increased_minus_decreased_mean_z"] == -0.2
    assert by_trait.loc["t1", "overall_direction"] == "inc"
    assert by_trait.loc["t2", "trait_direction"] == "increased"
    assert cells.groupby("trait_id")["trait_direction"].first().to_dict() == {
        "t1": "decreased",
        "t2": "increased",
    }


def test_column_label_shows_lollipop_direction_and_both_contrasts() -> None:
    label = heatmap.trait_column_label("coccus shaped", "decreased", -0.0813, 0.0119)

    assert label.startswith("↓ coccus shaped")
    assert "lollipop -0.08" in label
    assert label.endswith("baseline Δz\n+0.01")
    assert heatmap.trait_column_label("halo", "increased", 0.49, 0.43).startswith("↑ halo")


def test_filter_all_low_support_taxa_removes_only_fully_low_support_rows() -> None:
    cells = pd.DataFrame(
        [
            ["taxon-a", "t1", "low support"],
            ["taxon-a", "t2", "low support"],
            ["taxon-b", "t1", "low support"],
            ["taxon-b", "t2", "mixed"],
        ],
        columns=["taxon_id", "trait_id", "status"],
    )
    rows = pd.DataFrame(
        [["taxon-a", "Species A"], ["taxon-b", "Species B"]],
        columns=["taxon_id", "taxon_name"],
    )

    filtered_cells, filtered_rows, removed = heatmap.filter_all_low_support_taxa(cells, rows)

    assert removed == 1
    assert filtered_rows["taxon_id"].tolist() == ["taxon-b"]
    assert filtered_cells["taxon_id"].unique().tolist() == ["taxon-b"]


def test_rows_sort_by_bugsigdb_inc_dec_value() -> None:
    rows = pd.DataFrame(
        [
            ["faecalibacterium", "Faecalibacterium prausnitzii", 0, 13, -float("inf")],
            ["roseburia", "Roseburia faecis", 0, 6, -float("inf")],
            ["balanced", "Ruminococcus bromii", 2, 2, 0.0],
            ["moderate-inc", "Mediterraneibacter torques", 3, 1, 1.585],
            ["inc-three", "Methanobrevibacter smithii", 3, 0, float("inf")],
            ["inc-five", "Bifidobacterium longum", 5, 0, float("inf")],
            ["inc-four", "Streptococcus salivarius", 4, 0, float("inf")],
            ["moderate-dec", "Escherichia coli", 1, 2, -1.0],
        ],
        columns=[
            "taxon_id",
            "taxon_name",
            "increased_profiles",
            "decreased_profiles",
            "log2_inc_dec_ratio",
        ],
    )

    ordered = heatmap.sort_rows_by_bugsigdb_inc_dec(rows)

    assert ordered["taxon_id"].tolist() == [
        "inc-five",
        "inc-four",
        "inc-three",
        "moderate-inc",
        "balanced",
        "moderate-dec",
        "roseburia",
        "faecalibacterium",
    ]


def test_literature_comparison_label_makes_disagreement_explicit() -> None:
    assert heatmap.direction_hamming_distance("increased", "increased") == 0
    assert heatmap.direction_hamming_distance("increased", "decreased") == 1
    assert pd.isna(heatmap.direction_hamming_distance("increased", "balanced"))
    assert (
        heatmap.literature_comparison_label("decreased", "differs", 1) == "lit dec\nH=1 different"
    )
    assert heatmap.literature_comparison_label("increased", "agrees", 0) == "lit inc\nH=0 similar"
    assert (
        heatmap.literature_comparison_label("decreased", "balanced in BugSigDB")
        == "lit dec\nbalanced"
    )


def test_influence_display_filter_excludes_taxon_and_small_effects() -> None:
    rows = pd.DataFrame(
        [
            ["NCBITaxon:33039", "Mediterraneibacter torques"],
            ["taxon-small", "Small effect"],
            ["taxon-kept", "Kept effect"],
        ],
        columns=["taxon_id", "taxon_name"],
    )
    cells = pd.DataFrame(
        [
            ["NCBITaxon:33039", "t1", 0.5],
            ["taxon-small", "t1", 0.01],
            ["taxon-small", "t2", -0.01],
            ["taxon-kept", "t1", 0.0101],
            ["taxon-kept", "t2", 0.002],
        ],
        columns=["taxon_id", "trait_id", "mean_influence_contrast_z"],
    )

    filtered_cells, filtered_rows, excluded, below_threshold = heatmap.filter_influence_display(
        cells, rows
    )

    assert filtered_rows["taxon_id"].tolist() == ["taxon-kept"]
    assert filtered_rows["displayed_trait_cells"].tolist() == [1]
    assert filtered_cells["displayed_in_heatmap"].tolist() == [True, False]
    assert excluded == 1
    assert below_threshold == 1


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


def test_pd_taxa_reference_table_adds_literature_and_interpretation(tmp_path: Path) -> None:
    rows = pd.DataFrame(
        [
            [
                "taxon-e",
                "Escherichia coli",
                "species",
                1,
                1,
                2,
                "1/2",
                0.5,
                -1.0,
                3,
                0.8,
                1,
                3,
                "1/3",
                0.25,
                -1.0,
            ],
            [
                "taxon-f",
                "Faecalibacterium prausnitzii",
                "species",
                2,
                0,
                4,
                "0/4",
                0.0,
                -float("inf"),
                4,
                0.6,
                4,
                1,
                "4/1",
                0.8,
                1.0,
            ],
        ],
        columns=[
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
            "annotated_trait_count",
            "increased_profiles",
            "decreased_profiles",
            "inc_dec_counts",
            "inc_dec_ratio",
            "log2_inc_dec_ratio",
            "paired_support_studies",
            "max_absolute_influence_z",
            "pd_feces_profiles",
            "other_feces_profiles",
            "pd_other_feces_counts",
            "pd_specificity_fraction",
            "pd_vs_other_feces_log2_ratio",
        ],
    )
    cells = pd.DataFrame(
        [
            ["taxon-e", "Trait one", "t1", True, 0.2, "mixed"],
            ["taxon-e", "Trait two", "t2", False, 0.8, "inc"],
            ["taxon-f", "Trait one", "t1", True, -0.6, "dec"],
            ["taxon-f", "Trait two", "t2", True, -0.1, "mixed"],
        ],
        columns=[
            "taxon_id",
            "trait",
            "trait_id",
            "has_trait_annotation",
            "mean_influence_contrast_z",
            "status",
        ],
    )

    table = heatmap.build_pd_taxa_reference_table(
        rows,
        cells,
        heatmap.literature_reference_table(),
    ).set_index("taxon_id")

    assert table.loc["taxon-e", "literature_taxon"] == "Escherichia-Shigella"
    assert table.loc["taxon-e", "bugsigdb_pd_direction"] == "decreased"
    assert table.loc["taxon-e", "direction_hamming_distance"] == 1
    assert table.loc["taxon-e", "hamming_interpretation"] == "different"
    assert table.loc["taxon-e", "direction_agreement"] == "differs"
    assert "review cohort" in table.loc["taxon-e", "interpretation"]
    assert table.loc["taxon-f", "literature_direction"] == "decreased"
    assert table.loc["taxon-f", "direction_hamming_distance"] == 0
    assert table.loc["taxon-f", "direction_agreement"] == "agrees"
    assert "10.1002/mds.28119" in table.loc["taxon-f", "doi_citations"]
    assert table.loc["taxon-f", "kg_microbe_direct_trait_count"] == 2

    output_path = tmp_path / "reference.md"
    heatmap.write_reference_interpretation(table.reset_index(), output_path)
    output = output_path.read_text(encoding="utf-8")
    assert "## Agreement Summary" in output
    assert "## Methods: Literature Review and Direction Comparison" in output
    assert "targeted narrative review" in output
    assert "increased-profile count is greater than decreased-profile count" in output
    assert "[10.1002/mds.28119](https://doi.org/10.1002/mds.28119)" in output
