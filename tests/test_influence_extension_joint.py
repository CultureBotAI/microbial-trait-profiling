"""Exact joint-removal and correction-placement regressions for #960 (2026-10-23).

Adapt the lead's four-study/five-pair hand fixture, including a genus with two
family parents. Inputs are synthetic annotations; outputs are in-memory rows
and parsed Markdown cells. No production data or physiological claim is used.
"""

from __future__ import annotations

import json
from collections import defaultdict
from fractions import Fraction

import pytest
from scipy.stats import false_discovery_control

from mtp import influence_extension as extension
from mtp.contribution_inputs import Base, Graph
from mtp.taxon_contributions import Record, reconcile
from mtp.taxon_influence import comparison, removal

TRAITS = [f"METPO:t{i}" for i in range(6)]
BACKGROUND = "METPO:bg"
LACHNO, OSCILLO, OTHER = "NCBITaxon:186803", "NCBITaxon:216572", "NCBITaxon:900"
GENUS_L, GENUS_O, GENUS_X, GENUS_MULTI = (f"NCBITaxon:{i}" for i in (810, 820, 830, 840))
SPECIES_L, SPECIES_O, SPECIES_X, SPECIES_MULTI, SPECIES_B, SPECIES_C = (
    f"NCBITaxon:{i}" for i in range(1, 7)
)
ANNOTATIONS = {
    SPECIES_L: {TRAITS[0], BACKGROUND},
    SPECIES_O: {TRAITS[0], TRAITS[1], BACKGROUND},
    SPECIES_X: {TRAITS[1], BACKGROUND},
    SPECIES_MULTI: {TRAITS[0], BACKGROUND},
    SPECIES_B: {BACKGROUND},
    SPECIES_C: {TRAITS[0], BACKGROUND},
}
DESIGN = {
    ("s1", "e1"): ({SPECIES_X, SPECIES_B, SPECIES_C}, {SPECIES_L, SPECIES_O, SPECIES_B}),
    ("s2", "e1"): ({SPECIES_X}, {SPECIES_L, SPECIES_MULTI, SPECIES_B, SPECIES_C}),
    ("s3", "e1"): (
        {SPECIES_X, SPECIES_O, SPECIES_C},
        {SPECIES_L, SPECIES_O, SPECIES_MULTI, SPECIES_B},
    ),
    ("s3", "e2"): ({SPECIES_X, SPECIES_B}, {SPECIES_O, SPECIES_X}),
    ("s4", "e1"): ({SPECIES_L, SPECIES_X}, {SPECIES_MULTI, SPECIES_O, SPECIES_B}),
}
POOLED_JOINT = {
    "scope": "pooled",
    "rank": "family",
    "lineage": f"{LACHNO}+{OSCILLO}",
    "lineages": json.dumps([LACHNO, OSCILLO]),
    "name": "Lachnospiraceae + Oscillospiraceae",
    "scenario": "extension_exclude_joint_lachno_oscillo",
    "selection_kind": "declared_joint",
    "selection_reasons": "[]",
}


@pytest.fixture
def joint_corpus() -> tuple[Base, Graph]:
    parents = {
        SPECIES_L: {GENUS_L},
        SPECIES_O: {GENUS_O},
        SPECIES_X: {GENUS_X},
        SPECIES_MULTI: {GENUS_MULTI},
        SPECIES_B: {GENUS_X},
        SPECIES_C: {GENUS_X},
        GENUS_L: {LACHNO},
        GENUS_O: {OSCILLO},
        GENUS_X: {OTHER},
        GENUS_MULTI: {LACHNO, OSCILLO},
    }
    graph = Graph(parents, ANNOTATIONS, {}, {}, False)
    columns = [*TRAITS, BACKGROUND]
    matrix = {taxon: dict.fromkeys(values, 1) for taxon, values in ANNOTATIONS.items()}
    profiles, reported, metadata = {}, {}, {}
    for pair, (inc, dec) in DESIGN.items():
        metadata[pair] = {"condition": "Parkinson's disease", "body_site": "Feces"}
        for direction, taxa in (("inc", inc), ("dec", dec)):
            key = (*pair, direction)
            profiles[key] = {c: sum(matrix[t].get(c, 0) for t in taxa) for c in columns}
            reported[key] = {t: {"signature"} for t in taxa}
    ranks = dict.fromkeys(ANNOTATIONS, "species")
    ranks.update(dict.fromkeys((GENUS_L, GENUS_O, GENUS_X, GENUS_MULTI), "genus"))
    ranks.update(dict.fromkeys((LACHNO, OSCILLO, OTHER), "family"))
    base = Base(
        columns,
        dict(zip([f"k{i}" for i in range(6)], TRAITS, strict=True)),
        profiles,
        metadata,
        reported,
        matrix,
        {LACHNO: "Lachnospiraceae", OSCILLO: "Oscillospiraceae", OTHER: "Otheraceae"},
        ranks,
        [],
    )
    assert reconcile(base, graph).mismatches == []
    return base, graph


def _hand_effect(excluded: set[str]) -> Fraction:
    """Average exact pair fractions within studies, then equally across studies."""
    studies: defaultdict[str, list[Fraction]] = defaultdict(list)
    for (study, _), (inc, dec) in DESIGN.items():
        left, right = inc - excluded, dec - excluded
        left_total = sum(len(ANNOTATIONS[t]) for t in left)
        right_total = sum(len(ANNOTATIONS[t]) for t in right)
        if left_total and right_total:
            studies[study].append(
                Fraction(sum(TRAITS[0] in ANNOTATIONS[t] for t in left), left_total)
                - Fraction(sum(TRAITS[0] in ANNOTATIONS[t] for t in right), right_total)
            )
    means = [sum(values, Fraction(0)) / len(values) for values in studies.values()]
    return sum(means, Fraction(0)) / len(means)


def _joint_effects(base: Base, graph: Graph) -> list[Record]:
    effects, ambiguities = extension.lineage_effects(
        base, graph, reconcile(base, graph), [dict(POOLED_JOINT)], extension.load_protocol()
    )
    assert ambiguities == []
    return effects


def test_joint_effects_match_exact_hand_calculation(joint_corpus: tuple[Base, Graph]) -> None:
    base, graph = joint_corpus
    before = _hand_effect(set())
    after = _hand_effect({SPECIES_L, SPECIES_O, SPECIES_MULTI})
    only_l = _hand_effect({SPECIES_L, SPECIES_MULTI})
    only_o = _hand_effect({SPECIES_O, SPECIES_MULTI})
    assert before == Fraction(-1327, 6720)
    assert after == Fraction(-1, 480)
    assert only_l == Fraction(-601, 3360)
    assert only_o == Fraction(-13, 160)
    cell = next(row for row in _joint_effects(base, graph) if row["trait_key"] == TRAITS[0])
    for population, expected in (
        ("before", before),
        ("after", after),
        ("before_matched", before),
        ("after_matched", after),
        ("matched_change", Fraction(1313, 6720)),
    ):
        assert cell[f"{population}_effect"] == pytest.approx(float(expected), rel=0, abs=1e-12)
        assert cell[f"{population}_pairs"] == 5
        assert cell[f"{population}_studies"] == 4
    assert cell["lost_pairs"] == cell["gained_pairs"] == 0
    assert cell["after_effect"] != pytest.approx(float(before))
    assert cell["after_effect"] != pytest.approx(float(only_l + only_o - before))
    for lineage, expected in ((LACHNO, only_l), (OSCILLO, only_o)):
        scenario = {
            **POOLED_JOINT,
            "selection_kind": "ranked",
            "lineage": lineage,
            "lineages": json.dumps([lineage]),
        }
        rows, _ = extension.lineage_effects(
            base, graph, reconcile(base, graph), [scenario], extension.load_protocol()
        )
        assert rows[0]["after_effect"] == pytest.approx(float(expected), rel=0, abs=1e-12)


def test_joint_removal_deduplicates_multipath_sources(joint_corpus: tuple[Base, Graph]) -> None:
    base, graph = joint_corpus
    ledger = reconcile(base, graph)
    remaining, ambiguous = extension.joint_removal(base, graph, ledger, [LACHNO, OSCILLO])
    assert ambiguous == {}
    for pair, (inc, dec) in DESIGN.items():
        for direction, taxa in (("inc", inc), ("dec", dec)):
            kept = taxa - {SPECIES_L, SPECIES_O, SPECIES_MULTI}
            assert remaining[(*pair, direction)] == {
                trait: sum(trait in ANNOTATIONS[t] for t in kept) for trait in base.columns
            }
    for lineage in (LACHNO, OSCILLO, OTHER):
        joint, joint_ambiguous = extension.joint_removal(base, graph, ledger, [lineage])
        single, _, single_ambiguous = removal(base, graph, ledger, lineage)
        assert joint == single
        assert joint_ambiguous == single_ambiguous


def test_joint_effects_export_removed_mass_and_reduction(joint_corpus: tuple[Base, Graph]) -> None:
    base, graph = joint_corpus
    rows = _joint_effects(base, graph)
    for row in rows:
        assert row["removed_annotation_mass"] == 29
        assert row["removed_trait_annotation_mass"] == {TRAITS[0]: 12, TRAITS[1]: 5}.get(
            row["trait_key"], 0
        )
        assert row["removal_status"] == "removed_annotation_mass"
    cell = rows[0]
    for prefix in ("", "matched_"):
        assert cell[f"{prefix}effect_reduction_status"] == "defined"
        assert cell[f"{prefix}effect_reduction_fraction"] == pytest.approx(
            float(Fraction(1313, 1327))
        )
        assert cell[f"{prefix}effect_reduction_percent"] == pytest.approx(
            float(Fraction(131300, 1327))
        )
    for row in rows[2:]:
        for prefix in ("", "matched_"):
            assert row[f"{prefix}effect_reduction_status"] == "zero_before_effect"
            assert row[f"{prefix}effect_reduction_fraction"] == ""
            assert row[f"{prefix}effect_reduction_percent"] == ""


def test_absent_joint_lineage_reports_no_removed_mass(joint_corpus: tuple[Base, Graph]) -> None:
    base, graph = joint_corpus
    scenario = {**POOLED_JOINT, "lineage": "NCBITaxon:absent", "lineages": '["NCBITaxon:absent"]'}
    rows, _ = extension.lineage_effects(
        base, graph, reconcile(base, graph), [scenario], extension.load_protocol()
    )
    assert all(row["removed_annotation_mass"] == 0 for row in rows)
    assert all(row["removed_trait_annotation_mass"] == 0 for row in rows)
    assert all(row["removal_status"] == "no_removed_annotation_mass" for row in rows)
    assert rows[0]["before_effect"] == rows[0]["after_effect"]
    assert rows[0]["effect_reduction_fraction"] == 0


@pytest.mark.parametrize("add_zero_pair", [False, True])
def test_removed_mass_rows_reconcile_all_denominator_traits_and_directions(
    joint_corpus: tuple[Base, Graph], add_zero_pair: bool
) -> None:
    base, graph = joint_corpus
    if add_zero_pair:
        pair = ("ineligible", "e1")
        base.metadata[pair] = {"condition": "Parkinson's disease", "body_site": "Feces"}
        for direction, taxa in (("inc", set()), ("dec", {SPECIES_L, SPECIES_O})):
            key = (*pair, direction)
            base.reported[key] = {t: {"signature"} for t in taxa}
            base.profiles[key] = {
                trait: sum(trait in ANNOTATIONS[t] for t in taxa) for trait in base.columns
            }
    ledger = reconcile(base, graph)
    assert ledger.mismatches == []
    mass_rows: list[Record] = []
    effects, ambiguities = extension.lineage_effects(
        base,
        graph,
        ledger,
        [dict(POOLED_JOINT)],
        extension.load_protocol(),
        removed_mass_rows=mass_rows,
    )
    assert ambiguities == []
    assert len(effects) == 6
    assert len(mass_rows) == 2 * len(base.columns)
    lookup = {(row["direction"], row["trait_key"]): row for row in mass_rows}
    assert len(lookup) == len(mass_rows)
    for direction, before, after, removed in (("inc", 21, 16, 5), ("dec", 32, 8, 24)):
        rows = [lookup[(direction, trait)] for trait in base.columns]
        for row in rows:
            assert row["scope"] == "pooled"
            assert row["scenario"] == POOLED_JOINT["scenario"]
            assert row["population"] == "eligible_baseline_pairs"
            assert row["unit"] == "retained_trait_annotation_count"
            assert row["profile_count"] == 5
            assert row["before_annotation_mass"] == before
            assert row["after_annotation_mass"] == after
            assert row["removed_annotation_mass"] == removed
            assert row["removal_status"] == "removed_annotation_mass"
        for field, expected in (
            ("before_trait_annotation_mass", before),
            ("after_trait_annotation_mass", after),
            ("removed_trait_annotation_mass", removed),
        ):
            assert sum(row[field] for row in rows) == expected
    for trait, expected_mass in ((TRAITS[0], (2, 10)), (TRAITS[1], (1, 4)), (BACKGROUND, (2, 10))):
        assert (
            tuple(
                lookup[(direction, trait)]["removed_trait_annotation_mass"]
                for direction in ("inc", "dec")
            )
            == expected_mass
        )
    for row in effects:
        assert row["removed_annotation_mass"] == 29
        assert row["removed_trait_annotation_mass"] == sum(
            lookup[(direction, row["trait_key"])]["removed_trait_annotation_mass"]
            for direction in ("inc", "dec")
        )


def test_reduction_distinguishes_full_and_matched_support_without_clamping(
    joint_corpus: tuple[Base, Graph],
) -> None:
    base, graph = joint_corpus
    # The second study's increased profile now vanishes under the joint drop.
    # Remaining studies are unchanged, so all expectations stay exact by hand.
    key = ("s2", "e1", "inc")
    base.reported[key] = {SPECIES_L: {"signature"}}
    base.profiles[key] = {trait: int(trait in ANNOTATIONS[SPECIES_L]) for trait in base.columns}
    assert reconcile(base, graph).mismatches == []
    cell = _joint_effects(base, graph)[0]
    for population, expected_effect in (
        ("before", Fraction(-487, 6720)),
        ("after", Fraction(13, 120)),
        ("before_matched", Fraction(-607, 5040)),
        ("after_matched", Fraction(13, 120)),
        ("matched_change", Fraction(1153, 5040)),
    ):
        assert cell[f"{population}_effect"] == pytest.approx(
            float(expected_effect), rel=0, abs=1e-12
        )
    assert cell["lost_pairs"] == 1
    assert cell["before_pairs"] == 5
    assert cell["after_pairs"] == cell["before_matched_pairs"] == 4
    for prefix, expected in (("", Fraction(1215, 487)), ("matched_", Fraction(1153, 607))):
        assert cell[f"{prefix}effect_reduction_status"] == "defined"
        assert cell[f"{prefix}effect_reduction_fraction"] == pytest.approx(float(expected))
        assert cell[f"{prefix}effect_reduction_percent"] == pytest.approx(float(expected * 100))
        assert cell[f"{prefix}effect_reduction_fraction"] > 1


@pytest.mark.parametrize("empty_before", [False, True])
def test_reduction_marks_missing_support_explicitly(
    joint_corpus: tuple[Base, Graph], empty_before: bool
) -> None:
    base, graph = joint_corpus
    scenario = {
        **POOLED_JOINT,
        "scope": "saliva" if empty_before else "pooled",
        "lineages": json.dumps([LACHNO, OSCILLO, OTHER]),
    }
    rows, _ = extension.lineage_effects(
        base, graph, reconcile(base, graph), [scenario], extension.load_protocol()
    )
    cell = rows[0]
    assert cell["effect_reduction_status"] == (
        "no_before_effect" if empty_before else "no_after_effect"
    )
    assert cell["matched_effect_reduction_status"] == "no_before_effect"
    for prefix in ("", "matched_"):
        assert cell[f"{prefix}effect_reduction_fraction"] == ""
        assert cell[f"{prefix}effect_reduction_percent"] == ""


def test_lineage_bh_by_columns_match_independent_family_oracle(
    joint_corpus: tuple[Base, Graph], monkeypatch: pytest.MonkeyPatch
) -> None:
    base, graph = joint_corpus
    original = comparison
    calls = iter(range(6))
    pvalues = [0.001, 0.007, 0.019, 0.031, 0.12, None]
    populations = ("before", "after", "before_matched", "after_matched", "matched_change")

    def controlled_comparison(
        before: dict[tuple[str, str], float], after: dict[tuple[str, str], float]
    ) -> Record:
        row = original(before, after)
        p = pvalues[next(calls)]
        for multiplier, population in enumerate(populations, 1):
            row[f"{population}_p"] = "" if p is None else p * multiplier
        return row

    monkeypatch.setattr(extension, "comparison", controlled_comparison)
    rows = _joint_effects(base, graph)
    assert len(rows) == 6
    for multiplier, population in enumerate(populations, 1):
        complete_family = [1.0 if p is None else p * multiplier for p in pvalues]
        expected_bh = false_discovery_control(complete_family, method="bh")
        expected_by = false_discovery_control(complete_family, method="by")
        assert expected_bh[0] != pytest.approx(expected_by[0])
        for index, row in enumerate(rows):
            assert row["family_size"] == 6
            if pvalues[index] is None:
                assert row[f"{population}_q_bh_six"] == ""
                assert row[f"{population}_q_by_six"] == ""
            else:
                assert row[f"{population}_q_bh_six"] == pytest.approx(expected_bh[index])
                assert row[f"{population}_q_by_six"] == pytest.approx(expected_by[index])


def _markdown_rows(text: str, required_columns: set[str]) -> list[dict[str, str]]:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if not line.startswith("|"):
            continue
        headings = [cell.strip() for cell in line.strip("|").split("|")]
        if not required_columns <= set(headings):
            continue
        rows = []
        for following in lines[index + 2 :]:
            if not following.startswith("|"):
                break
            values = [cell.strip() for cell in following.strip("|").split("|")]
            rows.append(dict(zip(headings, values, strict=True)))
        return rows
    pytest.fail(f"No Markdown table has columns {sorted(required_columns)}")


def test_summary_bh_by_cells_are_in_exact_columns(joint_corpus: tuple[Base, Graph]) -> None:
    base, graph = joint_corpus
    effects = _joint_effects(base, graph)
    effects[0].update(
        before_q_bh_six=0.031,
        after_q_bh_six=0.047,
        before_q_by_six=0.076,
        after_q_by_six=0.115,
    )
    summaries = [
        {
            "scope": "pooled",
            "trait_key": TRAITS[0],
            "omissions": 4,
            "baseline_effect": 0.25,
            "baseline_q_bh_six": 0.031,
            "baseline_q_by_six": 0.076,
            "minimum_effect": -0.4,
            "maximum_effect": 0.8,
            "minimum_q_bh_six": 0.02,
            "maximum_q_bh_six": 0.07,
            "minimum_q_by_six": 0.05,
            "maximum_q_by_six": 0.2,
        }
    ]
    text = extension.combined_summary(None, [dict(POOLED_JOINT)], summaries, effects)
    effect_rows = _markdown_rows(text, {"BH-six before / after", "BY-six before / after"})
    focal = next(row for row in effect_rows if row["Trait"] == TRAITS[0])
    assert [float(value) for value in focal["BH-six before / after"].split(" / ")] == [0.031, 0.047]
    assert [float(value) for value in focal["BY-six before / after"].split(" / ")] == [0.076, 0.115]
    loso = _markdown_rows(text, {"BH-six q range", "BY-six q range"})
    assert len(loso) == 1
    assert [float(value) for value in loso[0]["BH-six q range"].split(" to ")] == [0.02, 0.07]
    assert [float(value) for value in loso[0]["BY-six q range"].split(" to ")] == [0.05, 0.2]
