"""Source-count conservation, lineage attribution and matched-population inference."""

import csv
import gzip
from collections import Counter
from pathlib import Path

import pytest

from mtp.contribution_inputs import Base, Graph, count, load_graph, rows
from mtp.taxon_contributions import export_ledger, reconcile
from mtp.taxon_influence import (
    adjust,
    comparison,
    describe,
    pair_values,
    removal,
    run_influence,
    scenarios,
)

A = "METPO:1000608"
B = "METPO:2000103|GO:0030245"
OTHER = "METPO:1000615"
SPECIES1 = "NCBITaxon:1"
SPECIES2 = "NCBITaxon:2"
GENUS = "NCBITaxon:10"
FAMILY = "NCBITaxon:100"


@pytest.fixture
def corpus() -> tuple[Base, Graph]:
    graph = Graph(
        {SPECIES1: {GENUS}, SPECIES2: {GENUS}, GENUS: {FAMILY}},
        {SPECIES1: {A, OTHER}, SPECIES2: {B, OTHER}, GENUS: {A}, FAMILY: {OTHER}},
        {},
        {},
        False,
    )
    base = Base(
        [A, B, OTHER],
        {"a": A, "b": B},
        {
            ("study", "experiment", "inc"): {A: 5, B: 2, OTHER: 6},
            ("study", "experiment", "dec"): {A: 0, B: 1, OTHER: 1},
        },
        {("study", "experiment"): {"condition": "condition", "body_site": "site"}},
        {
            ("study", "experiment", "inc"): {SPECIES1: {"sig1"}, GENUS: {"sig1"}, FAMILY: {"sig1"}},
            ("study", "experiment", "dec"): {SPECIES2: {"sig2"}},
        },
        {
            SPECIES1: {A: 1, OTHER: 1},
            SPECIES2: {B: 1, OTHER: 1},
            GENUS: {A: 2, B: 1, OTHER: 2},
            FAMILY: {A: 2, B: 1, OTHER: 3},
        },
        {},
        {SPECIES1: "species", SPECIES2: "species", GENUS: "genus", FAMILY: "family"},
        [],
    )
    return base, graph


def test_rollup_reconciles_and_exposes_parent_child_repetition(
    corpus: tuple[Base, Graph],
    tmp_path: Path,
) -> None:
    base, graph = corpus
    ledger = reconcile(base, graph)
    assert ledger.mismatches == []
    assert ledger.reconstructed[FAMILY] == Counter({A: 2, B: 1, OTHER: 3})
    counts = export_ledger(base, graph, ledger, tmp_path)
    assert counts["profile_counts"] == 4
    with gzip.open(tmp_path / "profile_counts.tsv.gz", "rt", encoding="utf-8") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    first = next(r for r in records if r["direction"] == "inc" and r["trait_key"] == A)
    assert first["profile_count"] == "5"
    assert first["distinct_source_taxa"] == "2"
    assert first["repeated_propagated_counts"] == "3"
    assert first["original_assertion_count"] == ""


def test_exported_source_contributions_sum_to_profile_counts(
    corpus: tuple[Base, Graph],
    tmp_path: Path,
) -> None:
    base, graph = corpus
    export_ledger(base, graph, reconcile(base, graph), tmp_path)

    source_counts: Counter[tuple[str, str]] = Counter()
    with gzip.open(tmp_path / "source_contributions.tsv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            source_counts[(row["reported_taxon"], row["trait_key"])] += count(
                row["count_contribution"]
            )

    reconstructed: Counter[tuple[str, str, str, str]] = Counter()
    with gzip.open(tmp_path / "profile_reported_taxa.tsv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for trait in base.traits.values():
                key = (row["study_id"], row["experiment_id"], row["direction"], trait)
                reconstructed[key] += source_counts[(row["reported_taxon"], trait)]

    expected: dict[tuple[str, str, str, str], int] = {}
    with gzip.open(tmp_path / "profile_counts.tsv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (
                row["study_id"],
                row["experiment_id"],
                row["direction"],
                row["trait_key"],
            )
            expected[key] = count(row["profile_count"])

    assert set(reconstructed).issubset(expected)
    assert all(reconstructed[key] == profile_count for key, profile_count in expected.items())


def test_lineage_removal_reaches_higher_rank_records_and_full_denominator(
    corpus: tuple[Base, Graph],
) -> None:
    base, graph = corpus
    remaining, removed, ambiguous = removal(base, graph, reconcile(base, graph), GENUS)
    assert removed[FAMILY] == Counter({A: 2, B: 1, OTHER: 2})
    assert remaining[("study", "experiment", "inc")] == {A: 0, B: 0, OTHER: 1}
    assert remaining[("study", "experiment", "dec")] == {A: 0, B: 0, OTHER: 0}
    assert ambiguous == {FAMILY: {FAMILY}}
    keys = list(base.profiles)
    assert pair_values(remaining, keys, A) == {}
    assert pair_values(remaining, keys, A, native_zero=True) == {("study", "experiment"): 0}


def test_reconciliation_does_not_silently_replace_profile_counts(
    corpus: tuple[Base, Graph],
) -> None:
    base, graph = corpus
    base.profiles[("study", "experiment", "inc")][A] += 1
    mismatch = reconcile(base, graph).mismatches
    assert len(mismatch) == 1
    assert mismatch[0]["level"] == "profile"
    assert mismatch[0]["expected"] == 6
    assert mismatch[0]["reconstructed"] == 5


def test_matching_reaverages_experiments_before_equal_study_weighting() -> None:
    before = {("s1", "e1"): 1.0, ("s1", "e2"): -1.0, ("s2", "e1"): 0.2}
    after = {("s1", "e1"): 0.5, ("s2", "e1"): 0.1}
    result = comparison(before, after)
    assert result["before_effect"] == pytest.approx(0.1)
    assert result["before_matched_effect"] == pytest.approx(0.6)
    assert result["matched_change_effect"] == pytest.approx(-0.3)
    assert result["eligibility_change"] == pytest.approx(0.5)
    assert result["lost_pairs"] == 1


def test_log_ratio_orientation_and_four_positive_rule() -> None:
    profiles = {
        ("s", "e", "inc"): {A: 8, B: 2, OTHER: 90},
        ("s", "e", "dec"): {A: 2, B: 2, OTHER: 0},
    }
    assert pair_values(profiles, list(profiles), A, B) == {("s", "e"): 2.0}
    profiles[("s", "e", "dec")][B] = 0
    assert pair_values(profiles, list(profiles), A, B) == {}


def test_intervals_target_different_parameters_and_keep_untestable() -> None:
    result = describe({("s1", "e"): 0.01, ("s2", "e"): 0.03})
    assert result["effect"] == pytest.approx(0.02)
    assert result["mean_ci_low"] < 0.01
    assert result["share_ci_high"] == 1
    tied = describe({("s1", "e"): 0.0, ("s2", "e"): 0.0})
    assert tied["status"] == "all_tied"
    assert tied["share_ci_low"] == ""
    assert describe({})["effect"] == ""
    assert adjust([0.01, None, 0.5]) == [0.03, None, 0.75]
    assert adjust([0.01, None, 0.5], by=True)[0] == pytest.approx(0.055)


def test_rank_resolution_and_scenario_support_are_not_name_heuristics(
    corpus: tuple[Base, Graph],
) -> None:
    base, graph = corpus
    assert graph.ranked(SPECIES1, "genus", base.ranks)[0] == GENUS
    assert graph.ranked(GENUS, "species", base.ranks) == ("", "unresolved_rank")
    assert all(
        r["status"] == "insufficient_support"
        for r in scenarios(base, graph, reconcile(base, graph), 3)
    )
    graph.parents[SPECIES1].add("NCBITaxon:11")
    graph.ancestor_cache.clear()
    base.ranks["NCBITaxon:11"] = "genus"
    assert graph.ranked(SPECIES1, "genus", base.ranks) == ("", "ambiguous_rank")


def test_pinned_header_and_agent_continuations_preserve_edge_identity(tmp_path: Path) -> None:
    edges, nodes = tmp_path / "edges.tsv", tmp_path / "nodes.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\trelation\tagent_type\tagent_type\n"
        "\tknowledge_level\tknowledge_source\tmeta\tprimary_knowledge_source\n"
        "e1\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1000608\tR\tmanual_agent\t\t"
        "observation\tGraph\t\tpublic:1\n"
        "e2\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1000608\tR\t\t"
        "manual_validation_of_automated_agent\n\tlogical_entailment\tGraph\t\tpublic:2\n",
        encoding="utf-8",
    )
    nodes.write_text("id\tname\nNCBITaxon:1\tSpecies one\n", encoding="utf-8")
    graph = load_graph(edges, nodes, {A}, {A})
    assert graph.traits[SPECIES1] == {A}
    assert len(graph.assertions[(SPECIES1, A)]) == 2
    assert graph.assertions[(SPECIES1, A)][1].source_row == 4
    assert graph.header_repaired
    assert graph.continued_records == 1


@pytest.mark.parametrize("invalid", ["", "-1", "1.5", "NaN", " 1", "\u0661"])
def test_invalid_counts_are_errors(invalid: str) -> None:
    with pytest.raises(ValueError, match="nonnegative integer"):
        count(invalid)


def test_duplicate_headers_and_cycles_are_errors(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.tsv"
    path.write_text("trait\ttrait\n1\t2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate header"):
        list(rows(path))
    graph = Graph({"a": {"b"}, "b": {"a"}}, {}, {}, {}, False)
    with pytest.raises(ValueError, match="cycle"):
        graph.ancestors("a")


def test_influence_exports_counts_backgrounds_missingness_and_summary(
    corpus: tuple[Base, Graph],
    tmp_path: Path,
) -> None:
    base, graph = corpus
    for number in range(4):
        trait = f"METPO:extra{number}"
        base.traits[f"extra{number}"] = trait
        base.columns.append(trait)
        for vector in base.profiles.values():
            vector[trait] = 0
    protocol = {
        "minimum_lineage_studies": 1,
        "top_lineages_per_rank": 1,
        "predicted_genera": [],
        "scopes": [{"id": "pooled", "condition": "*", "body_site": "*"}],
    }
    result = run_influence(base, graph, reconcile(base, graph), protocol, tmp_path)
    assert result["run_scenarios"] == 2
    assert result["effect_rows"] == 81
    background = list(rows(tmp_path / "exclude_genus_NCBITaxon_10" / "ora_background_counts.tsv"))
    original = next(
        r for r in background if r["trait_key"] == OTHER and r["treatment"] == "fixed_baseline"
    )
    remaining = next(
        r
        for r in background
        if r["trait_key"] == OTHER and r["treatment"] == "recomputed_remaining"
    )
    assert original["count"] == "7"
    assert remaining["count"] == "1"
    effects = list(rows(tmp_path / "pooled_fraction_effects.tsv"))
    removed = next(
        r for r in effects if r["scenario"] == f"exclude_genus_{GENUS}" and r["trait_key"] == A
    )
    assert removed["lost_pairs"] == "1"
    assert removed["after_effect"] == ""
    assert removed["after_status"] == "no_eligible_pairs"
    assert "not_rerun" in result["native_ora"]
    assert (tmp_path / "influence_summary.md").is_file()
