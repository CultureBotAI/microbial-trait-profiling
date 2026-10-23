import importlib.util
import math
import sys
from pathlib import Path

import pandas as pd
import pytest

from mtp.contribution_inputs import Base, Graph
from mtp.strain_rollup_record import write_strain_rollup_record
from mtp.taxon_contributions import Ledger

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
    assert result["pd_minus_other_feces_prevalence_difference"] == pytest.approx(0.4)
    assert result["pd_other_feces_display"] == "1/2"


def test_zero_other_count_is_displayed_as_undefined_not_as_a_ratio(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    selected = pd.DataFrame({"taxon": ["taxon-a"]})
    specificity = pd.DataFrame(
        {
            "taxon_id": ["taxon-a"],
            "pd_feces_profiles": [1],
            "other_feces_profiles": [0],
            "pd_other_feces_counts": ["1/0"],
            "pd_feces_studies": [1],
            "other_feces_studies": [0],
            "pd_specificity_fraction": [1.0],
            "pd_vs_other_feces_log2_ratio": [math.log2(3.0)],
        }
    )
    monkeypatch.setattr(
        contributors,
        "load_taxon_specificity",
        lambda _taxa, _metadata: specificity,
    )
    metadata = pd.DataFrame(
        [
            ["pd-1", contributors.PD_CONDITION, contributors.FECES, "inc"],
            ["pd-2", contributors.PD_CONDITION, contributors.FECES, "dec"],
            ["other-1", "Other condition", contributors.FECES, "inc"],
            ["other-2", "Other condition", contributors.FECES, "dec"],
        ],
        columns=["profile_id", "condition", "body_site_group", "inc_dec"],
    )

    result = contributors.augment_selected_taxa(selected, metadata).iloc[0]

    assert result["pd_other_feces_counts"] == "1/0"
    assert result["pd_other_feces_display"] == "---"
    assert result["pd_minus_other_feces_prevalence_difference"] == pytest.approx(0.5)


def test_large_difference_filter_removes_small_dot_rows() -> None:
    selected = pd.DataFrame(
        {
            "rank": ["species", "species", "family"],
            "taxon": ["taxon-a", "taxon-b", "taxon-c"],
            "display_name": ["Taxon A", "Taxon B", "Taxon C"],
        }
    )
    differences = pd.DataFrame(
        {
            "rank": ["species", "species", "species", "family"],
            "taxon": ["taxon-a", "taxon-a", "taxon-b", "taxon-c"],
            "increased_minus_decreased_share": [0.02, -0.08, 0.049, -0.05],
        }
    )

    result = contributors.select_large_difference_taxa(selected, differences)

    assert set(result["taxon"]) == {"taxon-a", "taxon-c"}
    maxima = result.set_index("taxon")["max_absolute_direction_share_difference"]
    assert maxima["taxon-a"] == pytest.approx(0.08)
    assert maxima["taxon-c"] == pytest.approx(0.05)


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
            "trait_direction": ["increased", "decreased", "decreased"],
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
    assert (changed["trait_a_direction"], changed["trait_b_direction"]) == (
        "increased",
        "decreased",
    )


F, G, H, S1, S2 = "NCBITaxon:100", "NCBITaxon:10", "NCBITaxon:20", "NCBITaxon:1", "NCBITaxon:2"
STRAIN, LONE, SHARED = "NCBITaxon:1001", "NCBITaxon:2001", "NCBITaxon:3"
INCREASED = ("st1", "e1", "inc")
DECREASED = ("st1", "e1", "dec")


def ledger_fixture() -> tuple[Base, Graph, Ledger]:
    """Family F holds genus G (species S1, with a strain, and S2) and genus H, with no species.

    SHARED sits under both S1 and S2, and LONE sits directly under H. The increased profile
    reports G and H; the decreased one reports S1 and, nested above it, G. All but S1 itself
    carry trait T1. The KG and the matrix name S1 differently.
    """
    parents = {
        STRAIN: {S1},
        S1: {G},
        S2: {G},
        G: {F},
        H: {F},
        LONE: {H},
        SHARED: {S1, S2},
    }
    trait_bearing = {STRAIN: {"T1"}, S2: {"T1"}, LONE: {"T1"}, SHARED: {"T1"}}
    graph = Graph(parents, trait_bearing, {}, {S1: "Species one", S2: "Species two"}, False)
    ranks = {F: "family", G: "genus", H: "genus", S1: "species", S2: "species"}
    sources = {
        target: {s for s in trait_bearing if target in graph.ancestors(s)} for target in (G, H, S1)
    }
    base = Base(
        columns=["T1"],
        traits={"facultatively_aerobic": "T1"},
        profiles={INCREASED: {"T1": 4}, DECREASED: {"T1": 5}},
        metadata={},
        reported={
            INCREASED: {G: {"sig-1"}, H: {"sig-1"}},
            DECREASED: {S1: {"sig-2"}, G: {"sig-2"}},
        },
        taxa={},
        names={F: "Family F", S1: "Matrix name one"},
        ranks=ranks,
        eligibility=[],
    )
    return base, graph, Ledger(sources, {}, [])


def test_the_ledger_reconciles_data_derived_under_the_recorded_scope(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Only a BacDive strain of S1 asserts the trait, so inputs built under the primary
    phenotype scope reconcile under the scope the manifest records, and an off ledger is
    refused with its scope named (#788)."""
    trait = "METPO:1000608"
    (tmp_path / "merged-kg_edges.tsv").write_text(
        "id\tsubject\tpredicate\tobject\n"
        f"s1\tkgmicrobe.strain:1\tbiolink:subclass_of\t{S1}\n"
        f"s2\tkgmicrobe.strain:1\tbiolink:has_phenotype\t{trait}\n",
        encoding="utf-8",
    )
    (tmp_path / "merged-kg_nodes.tsv").write_text("id\tname\n", encoding="utf-8")
    protocol = tmp_path / "protocol.json"
    protocol.write_text('{"kg_md5": {}}', encoding="utf-8")
    base = Base(
        columns=[trait],
        traits={"rod_shaped": trait},
        profiles={INCREASED: {trait: 1}},
        metadata={INCREASED[:2]: {"condition": "condition", "body_site": "site"}},
        reported={INCREASED: {S1: {"sig-1"}}},
        taxa={S1: {trait: 1}},
        names={},
        ranks={S1: "species"},
        eligibility=[],
    )

    def load_base(*paths: Path) -> Base:
        assert paths == tuple(tmp_path / name for name in contributors.LEDGER_INPUTS)
        return base

    monkeypatch.setattr(contributors, "load_base", load_base)

    assert contributors.STRAIN_ROLLUP == "phenotype"
    _, _, ledger = contributors.load_contribution_ledger(
        tmp_path, tmp_path, protocol, contributors.STRAIN_ROLLUP
    )
    assert ledger.mismatches == []
    with pytest.raises(ValueError, match="under the off strain rollup scope"):
        contributors.load_contribution_ledger(tmp_path, tmp_path, protocol, "off")


def test_the_ledger_refuses_data_derived_built_under_another_scope(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Inputs that record another scope are named before the KG is read (#847)."""
    for name in contributors.LEDGER_INPUTS[0], contributors.LEDGER_INPUTS[-1]:
        (tmp_path / name).write_text("study_id\texperiment_id\tinc_dec\n", encoding="utf-8")
        write_strain_rollup_record(tmp_path / name, "off")
    protocol = tmp_path / "protocol.json"
    protocol.write_text('{"kg_md5": {}}', encoding="utf-8")

    def unread(*paths: Path) -> Base:
        raise AssertionError(f"read {paths} after a scope mismatch")

    monkeypatch.setattr(contributors, "load_base", unread)
    with pytest.raises(ValueError, match="off strain rollup scope, not the declared phenotype"):
        contributors.load_contribution_ledger(tmp_path, tmp_path, protocol, "phenotype")


def test_contributions_resolve_each_source_to_one_species_and_family() -> None:
    base, graph, ledger = ledger_fixture()

    counts, coverage = contributors.count_contributions(base, graph, ledger, [INCREASED, DECREASED])

    by_cell = counts.set_index(["rank", "taxon", "direction"])["count"].to_dict()
    # Reported nested taxa each count their sources, as the profile trait sums do.
    assert by_cell == {
        ("family", F, "dec"): 5,
        ("family", F, "inc"): 4,
        ("species", S1, "dec"): 2,
        ("species", S1, "inc"): 1,
        ("species", S2, "dec"): 1,
        ("species", S2, "inc"): 1,
    }
    statuses = coverage[coverage["rank"].eq("species") & coverage["direction"].eq("inc")]
    assert statuses.set_index("status")["count"].to_dict() == {
        "ambiguous_rank": 1,
        contributors.RESOLVED_STATUS: 2,
        "unresolved_rank": 1,
    }
    # Every source of the profile's trait count is attributed somewhere at each rank.
    for key in (INCREASED, DECREASED):
        in_profile = coverage[coverage["direction"].eq(key[2])]
        assert int(in_profile["count"].sum()) == 2 * base.profiles[key]["T1"]


def test_contributor_tables_carry_shares_status_and_direction() -> None:
    base, graph, ledger = ledger_fixture()
    counts, coverage_counts = contributors.count_contributions(
        base, graph, ledger, [INCREASED, DECREASED]
    )
    directions = pd.DataFrame({"trait_id": ["T1"], "trait_direction": ["increased"]})
    traits = contributors.trait_frame(base, directions)

    cells, coverage, selected = contributors.build_contributor_tables(
        counts, coverage_counts, traits, contributors.contributor_names(base, graph)
    )

    species_inc = coverage.set_index(["rank", "direction"]).loc[("species", "inc")]
    assert (species_inc["total_count"], species_inc["resolved_count"]) == (4, 2)
    assert (species_inc["unresolved_rank_count"], species_inc["ambiguous_rank_count"]) == (1, 1)
    assert species_inc["unresolved_or_ambiguous_share"] == pytest.approx(0.5)
    species_dec = coverage.set_index(["rank", "direction"]).loc[("species", "dec")]
    assert (species_dec["unresolved_rank_count"], species_dec["ambiguous_rank_count"]) == (0, 2)
    assert selected["taxon"].tolist() == [S1, S2, F]
    assert selected["display_name"].tolist() == [
        "Species one (1)",
        "Species two (2)",
        "Family F (100)",
    ]
    cell = cells.set_index(["taxon", "direction"])
    assert cell.loc[(S1, "inc"), "share"] == pytest.approx(0.25)
    assert cell.loc[(S1, "dec"), "share"] == pytest.approx(0.4)
    assert cell.loc[(S2, "dec"), "status"] == "resolved"
    assert (cell.loc[(F, "inc"), "count"], cell.loc[(F, "inc"), "studies"]) == (4, 1)
    assert set(cells["trait_direction"]) == {"increased"}
    assert set(cells["trait_label"]) == {"facultatively aerobic (METPO:1000608)"}


def test_cells_without_contributions_are_absent() -> None:
    counts = pd.DataFrame(
        [("species", "NCBITaxon:7", "inc", "T1", 3, 2)],
        columns=["rank", "taxon", "direction", "trait_key", "count", "studies"],
    )
    coverage_counts = pd.DataFrame(
        [
            (rank, direction, "T1", contributors.RESOLVED_STATUS, total)
            for rank in contributors.CONTRIBUTOR_RANKS
            for direction, total in (("inc", 4), ("dec", 2))
        ],
        columns=["rank", "direction", "trait_key", "status", "count"],
    )
    traits = pd.DataFrame(
        [
            {
                "trait_key": "T1",
                "trait_slug": "facultatively_aerobic",
                "trait_label": "facultatively aerobic (METPO:1000608)",
                "trait_direction": "increased",
            }
        ]
    )

    cells, _, _ = contributors.build_contributor_tables(
        counts, coverage_counts, traits, {"NCBITaxon:7": "Species seven"}
    )

    by_direction = cells.set_index("direction")
    assert (by_direction.loc["inc", "status"], by_direction.loc["inc", "share"]) == (
        "resolved",
        0.75,
    )
    assert by_direction.loc["dec", ["status", "count", "studies"]].tolist() == ["absent", 0, 0]


def test_traits_without_a_lollipop_direction_are_refused() -> None:
    base, _, _ = ledger_fixture()
    directions = pd.DataFrame({"trait_id": ["other"], "trait_direction": ["increased"]})

    with pytest.raises(ValueError, match="has no direction for T1"):
        contributors.trait_frame(base, directions)


def test_only_paired_human_pd_fecal_profiles_are_counted() -> None:
    base, _, _ = ledger_fixture()
    metadata = pd.DataFrame(
        {
            "profile_id": ["st1|e1|inc", "st1|e1|dec", "st2|e2|inc", "st3|e3|inc"],
            "study_id": ["st1", "st1", "st2", "st3"],
            "experiment_id": ["e1", "e1", "e2", "e3"],
            "inc_dec": ["inc", "dec", "inc", "inc"],
            "condition": ["Parkinson's disease"] * 3 + ["Colorectal cancer"],
            "body_site_group": ["Feces"] * 4,
        }
    )

    assert contributors.pd_contribution_profiles(base, metadata) == [DECREASED, INCREASED]


def test_selection_keeps_top_species_per_cell_and_top_families(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = [("species", f"sp{index}", "inc", "T1", 10 - index, 1) for index in range(7)]
    rows.append(("species", "sp6", "dec", "T1", 50, 1))
    rows += [("family", f"fam{index}", "inc", "T1", 5 + index, 1) for index in range(3)]
    counts = pd.DataFrame(
        rows, columns=["rank", "taxon", "direction", "trait_key", "count", "studies"]
    )
    monkeypatch.setattr(contributors, "TOP_FAMILIES", 2)

    names = {str(taxon): f"Name of {taxon}" for taxon in counts["taxon"]}

    selected = contributors.select_contributor_taxa(counts, names)

    species = selected[selected["rank"].eq("species")]["taxon"].tolist()
    assert species == ["sp6", "sp0", "sp1", "sp2", "sp3", "sp4"]
    assert selected[selected["rank"].eq("family")]["taxon"].tolist() == ["fam2", "fam1"]

    monkeypatch.setattr(contributors, "MAX_SELECTED_SPECIES", 2)
    capped = contributors.select_contributor_taxa(counts, names)
    assert capped[capped["rank"].eq("species")]["taxon"].tolist() == ["sp6", "sp0"]


def test_unnamed_selected_taxa_are_refused() -> None:
    """A CURIE is a placeholder, not a name (#696)."""
    counts = pd.DataFrame(
        [
            ("species", "NCBITaxon:9", "inc", "T1", 5, 1),
            ("species", "NCBITaxon:10", "inc", "T1", 4, 1),
        ],
        columns=["rank", "taxon", "direction", "trait_key", "count", "studies"],
    )

    with pytest.raises(ValueError, match=r"No KG or taxa-trait matrix name .*NCBITaxon:10$"):
        contributors.select_contributor_taxa(counts, {"NCBITaxon:9": "Named"})


@pytest.mark.parametrize("kg_name", ["", "  "])
def test_an_empty_kg_name_does_not_hide_the_matrix_name(kg_name: str) -> None:
    """Some KG taxa have an empty name; the matrix name stands for it (#703)."""
    base, graph, _ = ledger_fixture()
    graph.names[S1] = kg_name

    names = contributors.contributor_names(base, graph)

    assert (names[S1], names[S2], names[F]) == ("Matrix name one", "Species two", "Family F")


def trait_intervals(values: dict[str, float]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "trait_id": list(values),
            "trait": [f"{trait_id} trait" for trait_id in values],
            "increased_minus_decreased_mean_z": list(values.values()),
            "trait_direction": [
                "increased" if value > 0 else "decreased" for value in values.values()
            ],
            "overall_ci95_low": 0.1,
            "overall_ci95_high": 0.2,
        }
    )


LOLLIPOP = {"T1": 0.49, "T2": -0.08}


def lollipop_directions(values: dict[str, float] = LOLLIPOP) -> pd.DataFrame:
    frame = trait_intervals(values)
    return frame[["trait_id", "trait", "increased_minus_decreased_mean_z", "trait_direction"]]


def test_heatmap_trait_table_from_the_current_lollipop_passes() -> None:
    contributors.check_trait_intervals(trait_intervals(LOLLIPOP), lollipop_directions())


@pytest.mark.parametrize(
    ("values", "problem"),
    [
        ({"T1": 0.49}, "missing traits: T2"),
        ({**LOLLIPOP, "T3": 0.1}, "extra traits: T3"),
        ({"T1": 0.49, "T2": 0.08}, "other directions: T2"),
        ({"T1": 0.50, "T2": -0.08}, "other mean-z differences: T1"),
    ],
)
def test_heatmap_trait_table_from_another_lollipop_is_refused(
    values: dict[str, float], problem: str
) -> None:
    """The contributors note quotes this table, so it must match the lollipop (#694)."""
    with pytest.raises(ValueError, match=rf"heatmap_traits\.tsv does not match .*{problem}"):
        contributors.check_trait_intervals(trait_intervals(values), lollipop_directions())


def test_heatmap_trait_table_without_directions_is_refused() -> None:
    old_format = trait_intervals(LOLLIPOP).drop(columns="trait_direction")

    with pytest.raises(ValueError, match="lacks trait_direction"):
        contributors.check_trait_intervals(old_format, lollipop_directions())


def test_selected_taxa_with_equal_totals_are_listed_by_name() -> None:
    counts = pd.DataFrame(
        [
            ("species", "NCBITaxon:9", "inc", "T1", 5, 1),
            ("species", "NCBITaxon:10", "inc", "T1", 5, 1),
        ],
        columns=["rank", "taxon", "direction", "trait_key", "count", "studies"],
    )
    names = {"NCBITaxon:9": "Zeta", "NCBITaxon:10": "Alpha"}

    selected = contributors.select_contributor_taxa(counts, names)

    assert selected["name"].tolist() == ["Alpha", "Zeta"]


def test_trait_labels_lead_with_the_lollipop_direction() -> None:
    coverage = pd.DataFrame(
        {
            "trait_slug": ["a", "b"],
            "trait_label": ["coccus shaped (METPO:1000668)", "moderately halophilic (x)"],
            "trait_direction": ["decreased", "increased"],
        }
    )

    assert contributors.trait_display_labels(coverage) == {
        "a": "↓ coccus shaped",
        "b": "↑ moderately halophilic",
    }
