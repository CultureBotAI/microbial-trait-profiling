"""Issue771: both-direction lineage selection and complete-family study LOSO."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from scipy.stats import binomtest, false_discovery_control

from mtp import influence_extension as extension
from mtp.contribution_inputs import Base, Graph
from mtp.influence_extension import fraction_family, lineage_effects, lineage_panel, study_loso
from mtp.paths import REPO_ROOT
from mtp.taxon_contributions import reconcile
from mtp.taxon_influence import effect_rows

TRAITS = [f"METPO:fixture{i}" for i in range(6)]
BACKGROUND = "METPO:background"


@pytest.fixture
def corpus() -> tuple[Base, Graph]:
    parents = {f"NCBITaxon:{i}": {f"NCBITaxon:{10 + i}"} for i in range(1, 5)}
    parents.update({f"NCBITaxon:{10 + i}": {f"NCBITaxon:{100 + i}"} for i in range(1, 5)})
    traits = {
        "NCBITaxon:1": {TRAITS[0], TRAITS[1], BACKGROUND},
        "NCBITaxon:2": {TRAITS[2], TRAITS[3], BACKGROUND},
        "NCBITaxon:3": {TRAITS[4], TRAITS[5], BACKGROUND},
        "NCBITaxon:4": {BACKGROUND},
    }
    graph = Graph(parents, traits, {}, {}, False)
    columns = [*TRAITS, BACKGROUND]
    matrix = {taxon: dict.fromkeys(values, 1) for taxon, values in traits.items()}
    profiles = {}
    reported = {}
    metadata = {}
    for study in range(1, 9):
        for experiment in range(1, 4 if study == 1 else 2):
            pair = (f"s{study}", f"e{experiment}")
            metadata[pair] = {"condition": "Parkinson's disease", "body_site": "Feces"}
            left = {"NCBITaxon:1", "NCBITaxon:3"}
            right = {"NCBITaxon:2", "NCBITaxon:4"}
            if study in {1, 2}:
                left, right = right, left
            for direction, taxa in (("inc", left), ("dec", right)):
                key = (*pair, direction)
                profiles[key] = {
                    trait: sum(matrix[t].get(trait, 0) for t in taxa) for trait in columns
                }
                reported[key] = {t: {"signature"} for t in taxa}
    ranks = {f"NCBITaxon:{i}": "species" for i in range(1, 5)}
    ranks.update({f"NCBITaxon:{10 + i}": "genus" for i in range(1, 5)})
    ranks.update({f"NCBITaxon:{100 + i}": "family" for i in range(1, 5)})
    base = Base(
        columns,
        dict(zip([f"trait{i}" for i in range(6)], TRAITS, strict=True)),
        profiles,
        metadata,
        reported,
        matrix,
        {},
        ranks,
        [],
    )
    return base, graph


def protocol() -> dict[str, object]:
    record: dict[str, object] = json.loads(
        (REPO_ROOT / "methods/table1_influence_extension_protocol.json").read_text()
    )
    return record


def test_panel_selects_both_directions_and_both_ranks(corpus: tuple[Base, Graph]) -> None:
    base, graph = corpus
    candidates, selected = lineage_panel(base, graph, reconcile(base, graph), protocol())
    assert {r["scope"] for r in selected if r["selection_kind"] == "ranked"} == {
        "pooled",
        "feces",
        "parkinsons_feces",
    }
    pooled = [r for r in candidates if r["scope"] == "pooled" and r["selected"]]
    assert {r["direction"] for r in pooled} == {"inc", "dec"}
    assert {r["rank"] for r in pooled} == {"family", "genus"}
    assert any(r["taxon"] == "NCBITaxon:102" and r["direction"] == "dec" for r in pooled)
    assert not any(r["scope"] == "saliva" and r["selection_kind"] == "ranked" for r in selected)


def test_panel_preserves_input_and_deterministic_tie_order(corpus: tuple[Base, Graph]) -> None:
    base, graph = corpus
    before = copy.deepcopy(base)
    ledger = reconcile(base, graph)
    first = lineage_panel(base, graph, ledger, protocol())
    base.profiles = dict(reversed(list(base.profiles.items())))
    assert lineage_panel(base, graph, ledger, protocol()) == first
    assert base.profiles == before.profiles
    assert base.traits == before.traits


def test_loso_omits_the_study_cluster_and_recorrects_every_trait(
    corpus: tuple[Base, Graph],
) -> None:
    base, _ = corpus
    keys = sorted(base.profiles)
    rows = study_loso(base, keys)
    assert len(rows) == 8 * 6
    omitted = [r for r in rows if r["omitted_study"] == "s1"]
    assert all(r["studies"] == 7 and r["pairs"] == 7 for r in omitted)
    expected_p = float(binomtest(6, 7).pvalue)
    assert all(r["p"] == pytest.approx(expected_p) for r in omitted)
    expected_by = false_discovery_control([expected_p] * 6, method="by")
    assert [r["q_by_six"] for r in omitted] == pytest.approx(expected_by)
    assert omitted[0]["q_bh_six"] != fraction_family(base, base.profiles, keys)[0]["q_bh_six"]


def test_zero_total_pairs_and_tied_traits_are_explicit(corpus: tuple[Base, Graph]) -> None:
    base, _ = corpus
    for key in base.profiles:
        if key[0] == "s8":
            base.profiles[key] = dict.fromkeys(base.columns, 0)
        else:
            base.profiles[key][TRAITS[5]] = 0
    result = fraction_family(base, base.profiles, sorted(base.profiles))
    tied = next(r for r in result if r["trait_key"] == TRAITS[5])
    assert tied["studies"] == 7
    assert tied["status"] == "all_tied"
    assert tied["p"] == ""
    assert tied["q_by_six"] == ""
    assert all(r["family_size"] == 6 for r in result)
    assert all(r["omitted_study"] != "s8" for r in study_loso(base, sorted(base.profiles)))


def test_incomplete_focal_family_is_rejected(corpus: tuple[Base, Graph]) -> None:
    base, _ = corpus
    base.traits.pop("trait0")
    with pytest.raises(ValueError, match="six"):
        study_loso(base, sorted(base.profiles))


def test_lineage_effects_recompute_full_denominator(corpus: tuple[Base, Graph]) -> None:
    base, graph = corpus
    ledger = reconcile(base, graph)
    assert ledger.mismatches == []
    _, panel = lineage_panel(base, graph, ledger, protocol())
    selected = [r for r in panel if r["scope"] == "pooled" and r["lineage"] == "NCBITaxon:11"]
    rows, _ = lineage_effects(base, graph, ledger, selected, protocol())
    assert len(rows) == 6
    third_source_trait = next(r for r in rows if r["trait_key"] == TRAITS[4])
    # Source1 carries another trait and background; its removal changes all shares.
    assert third_source_trait["before_effect"] == pytest.approx(1 / 12)
    assert third_source_trait["after_effect"] == pytest.approx(1 / 6)
    assert third_source_trait["lost_pairs"] == 0
    assert third_source_trait["before_matched_pairs"] == third_source_trait["after_matched_pairs"]
    assert all(r["family_size"] == 6 for r in rows)


def test_single_study_omission_keeps_six_undefined_members(corpus: tuple[Base, Graph]) -> None:
    base, _ = corpus
    keys = [key for key in base.profiles if key[0] == "s1"]
    rows = study_loso(base, keys)
    assert len(rows) == 6
    assert all(row["status"] == "no_eligible_pairs" for row in rows)
    assert all(row["p"] == "" for row in rows)
    assert all(row["q_by_six"] == "" for row in rows)


def test_combined_summary_retains_native_text_and_labels_dependent_diagnostics() -> None:
    from mtp.influence_extension import combined_summary

    text = combined_summary("# Declared baseline\nOriginal estimates.\n", [], [])
    assert "# Declared baseline\nOriginal estimates.\n" in text
    assert "Additional lineage exclusions" in text
    assert "Six-fraction study omissions" in text
    assert "dependent omissions" in text
    assert "not measurements of physiological activity" in text
    assert "Native baseline was not supplied" in combined_summary(None, [], [])


def test_native_summary_binding_rejects_an_arm_from_another_core(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from types import SimpleNamespace

    import mtp.influence_extension as extension

    core = tmp_path / "core.json"
    core.write_text("fixture core\n", encoding="utf-8")
    monkeypatch.setattr(
        extension,
        "verify_arm",
        lambda _path: SimpleNamespace(
            arm="noel",
            core_manifest_sha256="wrong",
            steps=[],
        ),
    )
    with pytest.raises(ValueError, match="this exact core"):
        extension.native_summary_source(core, tmp_path / "noel.json")


def ranking_corpus(*, half_empty: bool = False) -> tuple[Base, Graph]:
    """Six competing lineages; support and count order deliberately disagree."""
    columns = [*TRAITS, BACKGROUND]
    counts = {1: (4, 2), 2: (3, 2), 3: (3, 2), 4: (3, 2), 5: (2, 5)}
    if half_empty:
        counts[6] = (3, 3)
    traits = {f"NCBITaxon:{i}": {TRAITS[0], BACKGROUND} for i in counts}
    traits["NCBITaxon:99"] = {BACKGROUND}
    parents = {f"NCBITaxon:{i}": {f"NCBITaxon:{10 + i}"} for i in counts}
    parents.update({f"NCBITaxon:{10 + i}": {f"NCBITaxon:{100 + i}"} for i in counts})
    graph = Graph(parents, traits, {}, {}, False)
    matrix = {taxon: dict.fromkeys(keys, 1) for taxon, keys in traits.items()}
    profiles, reported, metadata = {}, {}, {}
    for source, (studies, experiments) in counts.items():
        for study in range(studies):
            for experiment in range(experiments):
                pair = (f"source{source}-study{study}", f"e{experiment}")
                metadata[pair] = {"condition": "Parkinson's disease", "body_site": "Feces"}
                for direction, taxon in (("inc", f"NCBITaxon:{source}"), ("dec", "NCBITaxon:99")):
                    key = (*pair, direction)
                    empty = source == 6 and study == 2 and direction == "dec"
                    profiles[key] = {
                        trait: 0 if empty else matrix[taxon].get(trait, 0) for trait in columns
                    }
                    reported[key] = {} if empty else {taxon: {"signature"}}
    ranks = {f"NCBITaxon:{i}": "species" for i in counts}
    ranks.update({f"NCBITaxon:{10 + i}": "genus" for i in counts})
    ranks.update({f"NCBITaxon:{100 + i}": "family" for i in counts})
    return Base(
        columns,
        dict(zip([f"trait{i}" for i in range(6)], TRAITS, strict=True)),
        profiles,
        metadata,
        reported,
        matrix,
        {},
        ranks,
        [],
    ), graph


def test_ranked_panel_has_top_three_supported_lineages_with_real_tie_at_cutoff() -> None:
    base, graph = ranking_corpus()
    ledger = reconcile(base, graph)
    assert ledger.mismatches == []
    candidates, _ = lineage_panel(base, graph, ledger, protocol())
    rows = {
        r["taxon"]: r
        for r in candidates
        if (r["scope"], r["trait_key"], r["direction"], r["rank"])
        == ("pooled", TRAITS[0], "inc", "genus")
    }
    assert {taxon for taxon, row in rows.items() if row["selected"]} == {
        "NCBITaxon:11",
        "NCBITaxon:12",
        "NCBITaxon:13",
    }
    assert [rows[f"NCBITaxon:{i}"]["count"] for i in range(11, 16)] == [8, 6, 6, 6, 10]
    assert [rows[f"NCBITaxon:{i}"]["eligible_rank"] for i in range(11, 15)] == [1, 2, 3, 4]
    assert rows["NCBITaxon:15"]["studies"] == 2
    assert rows["NCBITaxon:15"]["selection_status"] == "insufficient_studies"
    assert rows["NCBITaxon:14"]["selection_status"] == "outside_top_k"
    # Tie order cannot depend on graph, profile or reported-taxon insertion order.
    graph.parents = dict(reversed(list(graph.parents.items())))
    base.profiles = dict(reversed(list(base.profiles.items())))
    base.reported = dict(reversed(list(base.reported.items())))
    assert lineage_panel(base, graph, reconcile(base, graph), protocol())[0] == candidates


def test_ranking_uses_contributing_universe_including_half_empty_pairs() -> None:
    base, graph = ranking_corpus(half_empty=True)
    candidates, _ = lineage_panel(base, graph, reconcile(base, graph), protocol())
    row = next(
        r
        for r in candidates
        if (r["scope"], r["trait_key"], r["direction"], r["rank"], r["taxon"])
        == ("pooled", TRAITS[0], "inc", "genus", "NCBITaxon:16")
    )
    assert row["count"] == 9
    assert row["studies"] == 3
    assert row["selected"]
    assert row["eligible_rank"] == 1
    assert row["share"] == pytest.approx(9 / 45)


def test_declared_joint_panel_exists_even_in_scopes_without_support(
    corpus: tuple[Base, Graph],
) -> None:
    base, graph = corpus
    _, selected = lineage_panel(base, graph, reconcile(base, graph), protocol())
    joint = [row for row in selected if row["selection_kind"] == "declared_joint"]
    assert {row["scope"] for row in joint} == {"pooled", "feces", "saliva", "parkinsons_feces"}
    assert len(joint) == 4
    assert all(
        json.loads(row["lineages"]) == ["NCBITaxon:186803", "NCBITaxon:216572"] for row in joint
    )
    assert {row["scenario"] for row in joint} == {
        "extension_exclude_joint_cellulose_lachnospiraceae_oscillospiraceae"
    }


def test_joint_removal_unions_overlapping_descendants_without_double_subtraction(
    corpus: tuple[Base, Graph],
) -> None:
    base, graph = corpus
    graph.parents["NCBITaxon:11"].add("NCBITaxon:102")
    before = copy.deepcopy(base.profiles)
    ledger = reconcile(base, graph)
    remaining, ambiguous = extension.joint_removal(
        base, graph, ledger, ["NCBITaxon:101", "NCBITaxon:102", "NCBITaxon:101"]
    )
    for key, reported in base.reported.items():
        kept = set(reported) - {"NCBITaxon:1", "NCBITaxon:2"}
        expected = {
            trait: sum(base.taxa[taxon].get(trait, 0) for taxon in kept) for trait in base.columns
        }
        assert remaining[key] == expected
        assert all(value >= 0 for value in remaining[key].values())
    assert ambiguous == {}
    assert base.profiles == before


def test_loso_carries_declared_baseline_and_recomputes_complete_bh_family(
    corpus: tuple[Base, Graph],
) -> None:
    base, _ = corpus
    keys = sorted(base.profiles)
    baseline = {row["trait_key"]: row for row in fraction_family(base, base.profiles, keys)}
    rows = study_loso(base, keys)
    for study in {row["omitted_study"] for row in rows}:
        omission = [row for row in rows if row["omitted_study"] == study]
        expected = false_discovery_control([float(row["p"]) for row in omission], method="bh")
        assert [row["q_bh_six"] for row in omission] == pytest.approx(expected)
        for row in omission:
            original = baseline[row["trait_key"]]
            assert row["baseline_effect"] == original["effect"]
            assert row["baseline_q_bh_six"] == original["q_bh_six"]
            assert row["baseline_q_by_six"] == original["q_by_six"]


def test_fraction_baseline_matches_declared_shared_fraction_estimator(
    corpus: tuple[Base, Graph],
) -> None:
    base, _ = corpus
    keys = sorted(base.profiles)
    declared = {
        row["trait_key"]: row
        for row in effect_rows(base, base.profiles, keys)
        if row["representation"] == "shared_fraction"
    }
    assert len(declared) == 6
    for row in fraction_family(base, base.profiles, keys):
        native = declared[row["trait_key"]]
        assert native["correction"] == "BH6"
        for field in ("effect", "p", "studies", "pairs"):
            assert row[field] == native[f"before_{field}"]
        assert row["q_bh_six"] == native["before_q"]


def test_loso_summary_pins_actual_extrema_baseline_and_reversal_count() -> None:
    baseline = [
        {
            "scope": "pooled",
            "trait_key": TRAITS[0],
            "effect": 0.25,
            "q_bh_six": 0.031,
            "q_by_six": 0.076,
        }
    ]
    omissions = [
        {
            "scope": "pooled",
            "trait_key": TRAITS[0],
            "effect": effect,
            "q_bh_six": bh,
            "q_by_six": by,
        }
        for effect, bh, by in ((0.8, 0.07, 0.2), (-0.4, 0.04, 0.1), (0.1, 0.02, 0.05))
    ]
    row = extension.loso_summary(baseline, omissions)[0]
    expected = {
        "baseline_effect": 0.25,
        "baseline_q_bh_six": 0.031,
        "baseline_q_by_six": 0.076,
        "minimum_effect": -0.4,
        "maximum_effect": 0.8,
        "minimum_q_bh_six": 0.02,
        "maximum_q_bh_six": 0.07,
        "minimum_q_by_six": 0.05,
        "maximum_q_by_six": 0.2,
        "mean_shift_sign_reversal_count": 1,
    }
    for key, value in expected.items():
        assert row[key] == pytest.approx(value)
    assert "mean_sign_reversals" not in row


def test_combined_summary_includes_exclusion_effects_baseline_and_both_q_ranges(
    corpus: tuple[Base, Graph],
) -> None:
    base, graph = corpus
    _, panel = lineage_panel(base, graph, reconcile(base, graph), protocol())
    effects, _ = lineage_effects(base, graph, reconcile(base, graph), panel[:1], protocol())
    effects[0].update(
        before_effect=0.125, after_effect=0.375, before_q_bh_six=0.031, after_q_bh_six=0.047
    )
    summary = [
        {
            "scope": "pooled",
            "trait_key": TRAITS[0],
            "omissions": 3,
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
    text = extension.combined_summary(
        "Declared native baseline.", panel[:1], summary, effects=effects
    )
    assert "Declared native baseline." in text
    assert "BH" in text
    assert "BY" in text
    for number in (
        "0.125",
        "0.375",
        "0.031",
        "0.047",
        "0.076",
        "-0.4",
        "0.8",
        "0.02",
        "0.07",
        "0.05",
        "0.2",
    ):
        assert number in text
