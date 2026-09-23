import importlib.util
import math
import sys
from pathlib import Path

import pandas as pd
import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts/create_pd_trait_contributors.py"
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("pd_trait_contributors", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
contributors = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(contributors)


def test_specificity_normalizes_unequal_profile_pools(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    selected = pd.DataFrame({"taxon": ["taxon-a"]})
    specificity = pd.DataFrame(
        {
            "taxon_id": ["taxon-a"],
            "pd_feces_profiles": [1],
            "other_feces_profiles": [2],
            "pd_other_feces_counts": ["1/2"],
            "pd_feces_studies": [1],
            "other_feces_studies": [2],
            "pd_specificity_fraction": [1 / 3],
            "pd_vs_other_feces_log2_ratio": [math.log2(1.5 / 2.5)],
        }
    )
    monkeypatch.setattr(
        contributors,
        "load_taxon_specificity",
        lambda _taxa, _metadata: specificity,
    )
    metadata = pd.DataFrame(
        [
            [f"pd-{index}", contributors.PD_CONDITION, contributors.FECES, "inc"]
            for index in range(2)
        ]
        + [[f"other-{index}", "Other condition", contributors.FECES, "dec"] for index in range(20)],
        columns=["profile_id", "condition", "body_site_group", "inc_dec"],
    )

    result = contributors.augment_selected_taxa(selected, metadata).iloc[0]

    assert result["pd_feces_total_profiles"] == 2
    assert result["other_feces_total_profiles"] == 20
    assert result["pd_feces_prevalence"] == 0.5
    assert result["other_feces_prevalence"] == 0.1
    assert result["pd_vs_other_feces_log2_ratio"] < 0
    assert result["pd_vs_other_feces_log2_prevalence_ratio"] > 0


def test_direction_differences_omit_equal_share_pairs() -> None:
    cells = pd.DataFrame(
        [
            ["taxon-a", "Taxon A", "inc", 20, 0.2, 3],
            ["taxon-a", "Taxon A", "dec", 10, 0.1, 2],
            ["taxon-b", "Taxon B", "inc", 10, 0.1, 2],
            ["taxon-b", "Taxon B", "dec", 20, 0.1, 2],
        ],
        columns=["taxon", "display_name", "direction", "count", "share", "studies"],
    )
    cells = cells.assign(
        scope="pd",
        rank="species",
        name=cells["display_name"],
        trait_key="trait-1",
        trait_slug="trait_one",
        trait_label="trait one",
    )
    coverage = pd.DataFrame({"trait_slug": ["trait_one"], "trait_label": ["trait one"]})

    differences = contributors.build_direction_differences(cells, coverage)

    assert differences["taxon"].tolist() == ["taxon-a"]
    row = differences.iloc[0]
    assert row["increased_minus_decreased_share"] == pytest.approx(0.1)
    assert row["direction_pattern"] == "larger increased share"
    assert not bool(row["direction_exclusive"])


def test_background_circle_grey_area_matches_unresolved_share() -> None:
    coverage = pd.DataFrame(
        [
            ["trait_one", "trait one", "inc", 100, 80, 0.8, 0.2],
            ["trait_one", "trait one", "dec", 50, 25, 0.5, 0.5],
            ["trait_two", "trait two", "inc", 10, 9, 0.9, 0.1],
            ["trait_two", "trait two", "dec", 20, 15, 0.75, 0.25],
        ],
        columns=[
            "trait_slug",
            "trait_label",
            "direction",
            "total_count",
            "resolved_count",
            "resolved_share",
            "unresolved_or_ambiguous_share",
        ],
    )
    columns = contributors.trait_direction_columns(coverage)

    circles = contributors.direction_background_circles(coverage, columns)

    for _, pair in circles.groupby("trait_slug"):
        assert pair["outer_area"].max() == contributors.BACKGROUND_CIRCLE_AREA
    assert (circles["resolved_area"] / circles["outer_area"]).tolist() == pytest.approx(
        circles["resolved_share"].tolist()
    )


def test_cocarriage_display_keeps_only_direction_changes() -> None:
    selected = pd.DataFrame(
        [
            ["species", "taxon-a", "Taxon A", 10],
            ["family", "taxon-b", "Taxon B", 5],
        ],
        columns=["rank", "taxon", "display_name", "total_count"],
    )
    rows = []
    presence = {
        ("taxon-a", "trait_one", "inc"),
        ("taxon-a", "trait_one", "dec"),
        ("taxon-a", "trait_two", "inc"),
        ("taxon-a", "trait_two", "dec"),
        ("taxon-a", "trait_three", "inc"),
        ("taxon-a", "trait_three", "dec"),
        ("taxon-b", "trait_one", "inc"),
        ("taxon-b", "trait_two", "inc"),
    }
    ranks = {"taxon-a": "species", "taxon-b": "family"}
    for taxon in ranks:
        for trait in ["trait_one", "trait_two", "trait_three"]:
            for direction in contributors.DIRECTION_ORDER:
                rows.append(
                    [
                        ranks[taxon],
                        taxon,
                        trait,
                        direction,
                        int((taxon, trait, direction) in presence),
                    ]
                )
    cells = pd.DataFrame(
        rows,
        columns=["rank", "taxon", "trait_slug", "direction", "count"],
    )
    coverage = pd.DataFrame(
        {
            "trait_slug": ["trait_one", "trait_two", "trait_three"],
            "trait_label": ["trait one", "trait two", "trait three"],
        }
    )

    summary, all_edges, changed_edges = contributors.build_trait_cocarriage_tables(
        cells, selected, coverage
    )

    assert len(summary) == 2
    assert len(all_edges) == 3
    assert len(changed_edges) == 1
    changed = changed_edges.iloc[0]
    assert (changed["trait_a_slug"], changed["trait_b_slug"]) == (
        "trait_one",
        "trait_two",
    )
    assert changed["increased_shared_taxa"] == 2
    assert changed["decreased_shared_taxa"] == 1
    assert changed["increased_only_taxa"] == "Taxon B"
    assert changed["decreased_only_taxa"] == ""
