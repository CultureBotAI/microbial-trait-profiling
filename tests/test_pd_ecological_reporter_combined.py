import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_pd_ecological_reporter_combined.py"
)
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("pd_ecological_reporter_combined", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
reporter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = reporter
SPEC.loader.exec_module(reporter)


def _condition_rows(
    condition: str,
    *,
    significant_trait: str | None = None,
    paired_study_count: int = 3,
    significant_composite: bool = False,
) -> list[dict[str, object]]:
    rows = [
        {
            "condition": condition,
            "trait": trait,
            "paired_study_count": paired_study_count,
            "ci_excludes_zero": trait == significant_trait,
        }
        for trait in reporter.SIX_TRAITS
    ]
    rows.append(
        {
            "condition": condition,
            "trait": reporter.REPORTER_COMPOSITE,
            "paired_study_count": paired_study_count,
            "ci_excludes_zero": significant_composite,
        }
    )
    return rows


def test_qualifying_conditions_requires_one_of_six_individual_traits() -> None:
    summary = pd.DataFrame(
        [
            *_condition_rows(
                reporter.PD_CONDITION,
                significant_trait="arsenate reduction",
            ),
            *_condition_rows("Composite only", significant_composite=True),
            *_condition_rows(
                "Too little support",
                significant_trait="coccus shaped",
                paired_study_count=2,
            ),
            *_condition_rows("No change"),
        ]
    )

    selected = reporter.qualifying_conditions(summary)

    assert selected["condition"].tolist() == [reporter.PD_CONDITION]
    assert selected.iloc[0]["qualifying_significant_trait_count"] == 1
    assert selected.iloc[0]["qualifying_significant_traits"] == "arsenate reduction"


def test_panel_b_composite_uses_the_three_supported_positive_pd_traits() -> None:
    assert reporter.INCREASED_SIDE_TRAITS == [
        "moderately halophilic",
        "facultatively aerobic",
        "arsenate reduction",
    ]


def _traits() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "trait_id": ["fac", "cell", "coccus", "sphere", "ars", "halo"],
            "trait_label": [
                "facultatively aerobic",
                "cellulose catabolic process",
                "coccus shaped",
                "sphere shaped",
                "arsenate reduction",
                "moderately halophilic",
            ],
        }
    )


def test_condition_summary_uses_paired_study_contrasts_and_oriented_composite() -> None:
    metadata = pd.DataFrame(
        [
            ["s1", "e1", "inc", reporter.PD_CONDITION, "Neurologic", "Feces"],
            ["s1", "e1", "dec", reporter.PD_CONDITION, "Neurologic", "Feces"],
            ["s2", "e2", "inc", reporter.PD_CONDITION, "Neurologic", "Feces"],
            ["s2", "e2", "dec", reporter.PD_CONDITION, "Neurologic", "Feces"],
            ["s3", "e3", "inc", reporter.PD_CONDITION, "Neurologic", "Oral"],
        ],
        columns=[
            "study_id",
            "experiment_id",
            "inc_dec",
            "condition",
            "condition_group",
            "body_site_group",
        ],
        index=["p1", "p2", "p3", "p4", "p5"],
    )
    zscores = pd.DataFrame(
        [
            [2, 0, 0, 0, 1, 3],
            [0, 2, 0, 0, 0, 1],
            [1, 1, 0, 0, 3, 2],
            [0, 2, 0, 0, 1, 0],
            [100, 100, 100, 100, 100, 100],
        ],
        columns=["fac", "cell", "coccus", "sphere", "ars", "halo"],
        index=metadata.index,
    )
    summary, study = reporter.build_condition_trait_summaries(metadata, zscores, _traits())

    pd_rows = study[study["condition"].eq(reporter.PD_CONDITION)]
    assert pd_rows["study_id"].nunique() == 2
    s1 = pd_rows[pd_rows["study_id"].eq("s1")].set_index("trait")
    assert s1.loc["cellulose catabolic process", "contrast"] == -2
    assert s1.loc[reporter.THREE_TRAIT_COMPOSITE, "contrast"] == pytest.approx(5 / 3)
    assert s1.loc[reporter.REPORTER_COMPOSITE, "contrast"] == pytest.approx(11 / 3)
    composite = summary[summary["trait"].eq(reporter.REPORTER_COMPOSITE)].iloc[0]
    assert composite["paired_study_count"] == 2
    assert composite["increased_profile_count"] == 2
    assert composite["decreased_profile_count"] == 2
    cellulose = summary[summary["trait"].eq("cellulose catabolic process")].iloc[0]
    assert cellulose["study_weighted_mean_contrast"] == pytest.approx(-1.5)
    assert cellulose["contrast_definition"] == reporter.CONTRAST_DEFINITION


def test_confidence_interval_is_a_two_sided_95_percent_t_interval() -> None:
    low, high = reporter._confidence_interval(pd.Series([1.0, 2.0, 3.0]))
    assert low == pytest.approx(-0.484, abs=1e-3)
    assert high == pytest.approx(4.484, abs=1e-3)
    single = reporter._confidence_interval(pd.Series([1.0]))
    assert all(pd.isna(value) for value in single)


def test_condition_summary_interval_calls_and_study_counts() -> None:
    # Three studies, each with one increased and one decreased fecal profile, give
    # these increased-minus-decreased contrasts per trait.
    contrasts = {
        "halo": [1.0, 2.0, 3.0],
        "cell": [-1.0, -1.2, -0.8],
        "fac": [0.0, 0.0, 0.0],
        "ars": [1.0, -1.0, 2.0],
        "coccus": [0.0, 0.0, 0.0],
        "sphere": [0.0, 0.0, 0.0],
    }
    rows, metadata_rows = [], []
    for index, study in enumerate(["s1", "s2", "s3"]):
        for direction in ["inc", "dec"]:
            rows.append(
                {
                    trait: values[index] if direction == "inc" else 0.0
                    for trait, values in contrasts.items()
                }
            )
            metadata_rows.append([study, f"{study}-1", direction, "C", "C", "Feces"])
    metadata = pd.DataFrame(
        metadata_rows,
        columns=[
            "study_id",
            "experiment_id",
            "inc_dec",
            "condition",
            "condition_group",
            "body_site_group",
        ],
    )
    zscores = pd.DataFrame(rows, index=metadata.index)

    summary, _ = reporter.build_condition_trait_summaries(metadata, zscores, _traits())
    row = summary.set_index("trait")

    halophilic = row.loc["moderately halophilic"]
    assert halophilic["ci_low"] == pytest.approx(-0.484, abs=1e-3)
    assert halophilic["ci_high"] == pytest.approx(4.484, abs=1e-3)
    assert not halophilic["ci_excludes_zero"]
    assert halophilic["study_weighted_median_contrast"] == pytest.approx(2.0)
    cellulose = row.loc["cellulose catabolic process"]
    assert cellulose["ci_high"] < 0
    assert cellulose["ci_excludes_zero"]
    assert (cellulose["positive_study_count"], cellulose["negative_study_count"]) == (0, 3)
    assert cellulose["majority_direction_fraction"] == pytest.approx(1.0)
    arsenate = row.loc["arsenate reduction"]
    assert (arsenate["positive_study_count"], arsenate["negative_study_count"]) == (2, 1)
    assert arsenate["zero_study_count"] == 0
    assert arsenate["majority_direction_fraction"] == pytest.approx(2 / 3)
    assert arsenate["study_weighted_median_contrast"] == pytest.approx(1.0)
    aerobic = row.loc["facultatively aerobic"]
    assert aerobic["zero_study_count"] == 3
    assert aerobic["majority_direction_fraction"] == pytest.approx(0.0)
    assert not aerobic["ci_excludes_zero"]
    assert (halophilic["paired_study_count"], halophilic["increased_profile_count"]) == (3, 3)
    assert halophilic["decreased_profile_count"] == 3


def _pd_profiles() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = [
        ["s1|e1|inc", "s1", "e1", "inc", reporter.PD_CONDITION, "Feces"],
        ["s1|e1|dec", "s1", "e1", "dec", reporter.PD_CONDITION, "Feces"],
        ["s2|e2|inc", "s2", "e2", "inc", reporter.PD_CONDITION, "Feces"],
        ["s2|e2|dec", "s2", "e2", "dec", reporter.PD_CONDITION, "Feces"],
        ["s3|e3|inc", "s3", "e3", "inc", reporter.PD_CONDITION, "Saliva"],
        ["s4|e4|inc", "s4", "e4", "inc", "Diet", "Feces"],
    ]
    columns = ["profile_id", "study_id", "experiment_id", "inc_dec", "condition"]
    metadata = pd.DataFrame(rows, columns=[*columns, "body_site_group"]).set_index(
        "profile_id", drop=False
    )
    zscores = pd.DataFrame(0.0, index=metadata.index, columns=_traits()["trait_id"])
    # Halophilic is higher in increased PD fecal profiles, cellulose in decreased ones.
    zscores.loc[["s1|e1|inc", "s2|e2|inc"], "halo"] = [1.0, 2.0]
    zscores.loc["s1|e1|dec", "cell"] = 1.5
    return metadata, zscores


def _reported_taxa(rows: list[list[str]]) -> pd.DataFrame:
    keys = ["study_id", "experiment_id", "direction"]
    return pd.DataFrame(rows, columns=[*keys, "taxon_id", "taxon_name", "taxonomic_rank"])


def _annotations() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ["g1", "Genus one", "genus", 0, 3, 0, 0, 0, 2],
            ["sp1", "Species one", "species", 0, 0, 0, 0, 0, 1],
            ["f1", "Family one", "family", 1, 1, 1, 1, 1, 1],
        ],
        columns=["taxon_id", "taxon_name", "taxonomic_rank", *_traits()["trait_id"]],
    )


def test_pd_taxon_carriers_count_each_taxon_once_per_pd_fecal_profile() -> None:
    metadata, zscores = _pd_profiles()
    taxa = _reported_taxa(
        [
            ["s1", "e1", "increased", "g1", "Genus one", "genus"],
            ["s1", "e1", "increased", "g1", "Genus one", "genus"],
            ["s2", "e2", "increased", "g1", "Genus one", "genus"],
            ["s1", "e1", "decreased", "g1", "Genus one", "genus"],
            ["s2", "e2", "decreased", "sp1", "Species one", "species"],
            ["s3", "e3", "increased", "sp1", "Species one", "species"],
            ["s4", "e4", "increased", "f1", "Family one", "family"],
        ]
    )

    carriers = reporter.build_pd_taxon_carriers(
        metadata, zscores, _traits(), taxa, _annotations()
    ).set_index(["trait", "taxon_id"])

    genus_halo = carriers.loc[("moderately halophilic", "g1")]
    assert genus_halo["annotation_count"] == 2
    assert genus_halo["increased_profiles"] == 2
    assert genus_halo["increased_studies"] == 2
    assert genus_halo["increased_positive_z_profiles"] == 2
    assert genus_halo["decreased_profiles"] == 1
    assert genus_halo["decreased_positive_z_profiles"] == 0
    assert bool(genus_halo["positive_side_in_lollipop"])
    genus_cellulose = carriers.loc[("cellulose catabolic process", "g1")]
    assert genus_cellulose["decreased_positive_z_profiles"] == 1
    assert not bool(genus_cellulose["positive_side_in_lollipop"])
    # The oral PD profile and the Diet profile contribute nothing.
    species_halo = carriers.loc[("moderately halophilic", "sp1")]
    assert species_halo["increased_profiles"] == 0
    assert species_halo["decreased_profiles"] == 1
    assert "f1" not in carriers.index.get_level_values("taxon_id")


def test_pd_taxon_carriers_reject_taxa_missing_from_the_matrix() -> None:
    metadata, zscores = _pd_profiles()
    taxa = _reported_taxa([["s1", "e1", "increased", "unknown", "Unknown", "genus"]])
    with pytest.raises(ValueError, match="missing from trait matrix"):
        reporter.build_pd_taxon_carriers(metadata, zscores, _traits(), taxa, _annotations())


def test_taxon_panel_keeps_genus_and_species_and_marks_unannotated_cells() -> None:
    carriers = pd.DataFrame(
        [
            ["moderately halophilic", "g1", "Genus one", "genus", 2, 4, 1],
            ["cellulose catabolic process", "g1", "Genus one", "genus", 3, 4, 1],
            ["moderately halophilic", "sp1", "Species one", "species", 1, 0, 2],
            ["moderately halophilic", "f1", "Family one", "family", 9, 9, 9],
        ],
        columns=[
            "trait",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
            "annotation_count",
            "increased_studies",
            "decreased_studies",
        ],
    )

    cells, ranking = reporter.prepare_pd_taxon_panel(carriers, reporter.SIX_TRAITS)

    assert ranking["taxon_id"].tolist() == ["g1", "sp1"]
    assert len(cells) == 2 * len(reporter.SIX_TRAITS)
    cell = cells.set_index(["taxon_id", "trait"])
    assert cell.loc[("g1", "cellulose catabolic process"), "inc_dec_study_counts"] == "4/1"
    assert cell.loc[("g1", "moderately halophilic"), "direction_log2_ratio"] == pytest.approx(
        1.5849625
    )
    unannotated = cell.loc[("sp1", "cellulose catabolic process")]
    assert not bool(unannotated["trait_annotated"])
    assert unannotated["blank_reason"] == "no annotation in local taxa-trait matrix"
