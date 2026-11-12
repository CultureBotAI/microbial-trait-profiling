"""Lineage absence is distinct from zero annotation mass in the selected scope (#1030)."""

from __future__ import annotations

import json
from fractions import Fraction

import pytest

from mtp import influence_extension as extension
from mtp.contribution_inputs import Base, Graph
from mtp.taxon_contributions import Record, reconcile
from test_influence_extension_joint import (
    ANNOTATIONS,
    BACKGROUND,
    LACHNO,
    OSCILLO,
    POOLED_JOINT,
    SPECIES_L,
    TRAITS,
)
from test_influence_extension_joint import joint_corpus as joint_corpus

ABSENT = "NCBITaxon:99998"
ABSENT_OTHER = "NCBITaxon:99997"
EMPTY = "NCBITaxon:99999"
CASES = {
    "all_present": ([LACHNO, OSCILLO], "removed_annotation_mass", (5, 24)),
    "partly_absent": ([LACHNO, ABSENT], "partially_absent_lineage", (2, 12)),
    "all_absent": ([ABSENT, ABSENT_OTHER], "absent_lineages", (0, 0)),
    "present_without_contribution": ([EMPTY], "no_removed_annotation_mass", (0, 0)),
}


def _scenario(name: str, lineages: list[str]) -> Record:
    return {
        **POOLED_JOINT,
        "scenario": f"status_{name}",
        "lineage": "+".join(lineages),
        "lineages": json.dumps(lineages),
    }


def _add_parent_only_empty_lineage(graph: Graph) -> None:
    graph.parents["NCBITaxon:77770"] = {EMPTY}
    for lineage in (LACHNO, OSCILLO, EMPTY):
        assert lineage not in graph.names
        assert lineage not in graph.parents
        assert lineage not in graph.traits
        assert any(lineage in parents for parents in graph.parents.values())
    for lineage in (ABSENT, ABSENT_OTHER):
        assert lineage not in graph.names
        assert lineage not in graph.parents
        assert lineage not in graph.traits
        assert all(lineage not in parents for parents in graph.parents.values())


def _rows(base: Base, graph: Graph, scenarios: list[Record]) -> tuple[list[Record], list[Record]]:
    ledger = reconcile(base, graph)
    assert ledger.mismatches == []
    directional: list[Record] = []
    effects, ambiguous = extension.lineage_effects(
        base,
        graph,
        ledger,
        scenarios,
        extension.load_protocol(),
        removed_mass_rows=directional,
    )
    assert ambiguous == []
    return effects, directional


@pytest.mark.parametrize("case", list(CASES))
def test_lineage_presence_status_preserves_effect_and_directional_mass(
    joint_corpus: tuple[Base, Graph], case: str
) -> None:
    base, graph = joint_corpus
    _add_parent_only_empty_lineage(graph)
    lineages, expected_status, removed_by_direction = CASES[case]
    effects, directional = _rows(base, graph, [_scenario(case, lineages)])
    assert len(effects) == 6
    assert len(directional) == 2 * len(base.columns)
    assert {row["removal_status"] for row in effects} == {expected_status}
    assert {row["removal_status"] for row in directional} == {expected_status}
    expected_mass = sum(removed_by_direction)
    trait_masses = {
        "all_present": {TRAITS[0]: (2, 10), TRAITS[1]: (1, 4), BACKGROUND: (2, 10)},
        "partly_absent": {TRAITS[0]: (1, 6), BACKGROUND: (1, 6)},
        "all_absent": {},
        "present_without_contribution": {},
    }[case]
    first_trait = next(row for row in effects if row["trait_key"] == TRAITS[0])
    assert first_trait["before_effect"] == pytest.approx(
        float(Fraction(-1327, 6720)), rel=0, abs=1e-12
    )
    after_effect = {
        "all_present": Fraction(-1, 480),
        "partly_absent": Fraction(-601, 3360),
        "all_absent": Fraction(-1327, 6720),
        "present_without_contribution": Fraction(-1327, 6720),
    }[case]
    assert first_trait["after_effect"] == pytest.approx(float(after_effect), rel=0, abs=1e-12)
    for row in effects:
        trait = row["trait_key"]
        assert row["removed_annotation_mass"] == expected_mass
        assert row["removed_trait_annotation_mass"] == sum(trait_masses.get(trait, (0, 0)))
        accounting = json.loads(row["lineage_accounting"])
        assert set(accounting) == set(lineages)
        for lineage in lineages:
            individual_mass = {LACHNO: 14, OSCILLO: 21}.get(lineage, 0)
            individual_trait_mass = {
                LACHNO: {TRAITS[0]: 7},
                OSCILLO: {TRAITS[0]: 8, TRAITS[1]: 5},
            }.get(lineage, {}).get(trait, 0)
            assert accounting[lineage] == {
                "graph_status": "absent_from_graph"
                if lineage in (ABSENT, ABSENT_OTHER)
                else "present_in_graph",
                "removed_annotation_mass": individual_mass,
                "removed_trait_annotation_mass": individual_trait_mass,
            }
        if expected_mass == 0:
            assert row["before_effect"] == row["after_effect"]
            assert row["lost_pairs"] == 0
            assert row["gained_pairs"] == 0
    for index, (direction, before) in enumerate((("inc", 21), ("dec", 32))):
        selected = [row for row in directional if row["direction"] == direction]
        assert {row["trait_key"] for row in selected} == set(base.columns)
        for row in selected:
            assert row["profile_count"] == 5
            assert row["before_annotation_mass"] == before
            assert row["removed_annotation_mass"] == removed_by_direction[index]
            assert row["after_annotation_mass"] == before - removed_by_direction[index]
            assert (
                row["removed_trait_annotation_mass"]
                == trait_masses.get(row["trait_key"], (0, 0))[index]
            )
            assert (
                row["before_trait_annotation_mass"] - row["after_trait_annotation_mass"]
                == row["removed_trait_annotation_mass"]
            )
        assert (
            sum(row["removed_trait_annotation_mass"] for row in selected)
            == (removed_by_direction[index])
        )


@pytest.mark.parametrize("include_absent", [False, True])
def test_absence_status_takes_precedence_even_when_one_direction_has_zero_removed_mass(
    joint_corpus: tuple[Base, Graph], include_absent: bool
) -> None:
    base, graph = joint_corpus
    # This is the only increased profile containing a Lachnospiraceae source.
    # Its other taxon remains, so all five paired experiments remain eligible.
    key = ("s4", "e1", "inc")
    del base.reported[key][SPECIES_L]
    for trait in base.columns:
        base.profiles[key][trait] -= int(trait in ANNOTATIONS[SPECIES_L])
    lineages = [LACHNO, ABSENT] if include_absent else [LACHNO]
    effects, directional = _rows(base, graph, [_scenario("zero_increased_mass", lineages)])
    expected_effect_status = (
        "partially_absent_lineage" if include_absent else "removed_annotation_mass"
    )
    assert {row["removal_status"] for row in effects} == {expected_effect_status}
    for row in effects:
        assert row["removed_annotation_mass"] == 12
        assert row["before_pairs"] == 5
        assert row["after_pairs"] == 5
        accounting = json.loads(row["lineage_accounting"])
        assert accounting[LACHNO] == {
            "graph_status": "present_in_graph",
            "removed_annotation_mass": 12,
            "removed_trait_annotation_mass": 6 if row["trait_key"] == TRAITS[0] else 0,
        }
        if include_absent:
            assert accounting[ABSENT] == {
                "graph_status": "absent_from_graph",
                "removed_annotation_mass": 0,
                "removed_trait_annotation_mass": 0,
            }
    for row in directional:
        increased = row["direction"] == "inc"
        expected_status = "no_removed_annotation_mass" if increased else "removed_annotation_mass"
        if include_absent:
            expected_status = "partially_absent_lineage"
        assert row["removal_status"] == expected_status
        assert row["profile_count"] == 5
        assert row["before_annotation_mass"] == (19 if increased else 32)
        assert row["removed_annotation_mass"] == (0 if increased else 12)
        assert row["after_annotation_mass"] == (19 if increased else 20)
        expected_trait_mass = (
            6 if not increased and row["trait_key"] in (TRAITS[0], BACKGROUND) else 0
        )
        assert row["removed_trait_annotation_mass"] == expected_trait_mass


def test_presence_is_resolved_per_scenario_when_all_statuses_are_computed_together(
    joint_corpus: tuple[Base, Graph],
) -> None:
    base, graph = joint_corpus
    _add_parent_only_empty_lineage(graph)
    scenarios = [_scenario(name, details[0]) for name, details in CASES.items()]
    effects, directional = _rows(base, graph, scenarios)
    assert len(effects) == 6 * len(CASES)
    assert len(directional) == 2 * len(base.columns) * len(CASES)
    for name, (_, expected_status, _) in CASES.items():
        for rows in (effects, directional):
            statuses = {
                row["removal_status"] for row in rows if row["scenario"] == f"status_{name}"
            }
            assert statuses == {expected_status}
