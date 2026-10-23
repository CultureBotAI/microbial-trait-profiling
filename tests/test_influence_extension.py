"""Issue771: both-direction lineage selection and complete-family study LOSO."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from scipy.stats import binomtest, false_discovery_control

from mtp.contribution_inputs import Base, Graph
from mtp.influence_extension import fraction_family, lineage_effects, lineage_panel, study_loso
from mtp.paths import REPO_ROOT
from mtp.taxon_contributions import reconcile

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
    assert {r["scope"] for r in selected} == {"pooled", "feces", "parkinsons_feces"}
    pooled = [r for r in candidates if r["scope"] == "pooled" and r["selected"]]
    assert {r["direction"] for r in pooled} == {"inc", "dec"}
    assert {r["rank"] for r in pooled} == {"family", "genus"}
    assert any(r["taxon"] == "NCBITaxon:102" and r["direction"] == "dec" for r in pooled)
    assert not any(r["scope"] == "saliva" for r in selected)


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
