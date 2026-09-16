"""Gram-stain Sankey: paired unit contrasts, their tests, and the edge table (#222, #223)."""

from pathlib import Path

import pytest
from scipy.stats import binomtest, false_discovery_control

from bsdbng.gram_sankey import (
    EDGE_HEADER,
    build_edges,
    complementarity,
    drop_undefined_pairs,
    unit_contrasts,
    write_edge_table,
)
from bsdbng.trait_guilds import Profile


def _profile(
    study: str,
    direction: str,
    pos: int,
    neg: int,
    *,
    experiment: str = "1",
    condition: str = "Asthma",
    other: int = 0,
) -> Profile:
    return Profile(
        study_id=study,
        experiment_id=f"{study}-{experiment}",
        inc_dec=direction,
        condition=condition,
        body_site="Feces",
        guild_counts={"gram positive": pos, "gram negative": neg},
        total=pos + neg + other,
    )


def _paired(
    study: str,
    inc: tuple[int, int],
    dec: tuple[int, int],
    *,
    experiment: str = "1",
    condition: str = "Asthma",
) -> list[Profile]:
    return [
        _profile(study, "inc", *inc, experiment=experiment, condition=condition),
        _profile(study, "dec", *dec, experiment=experiment, condition=condition),
    ]


UP = ((8, 2), (2, 8))  # gram positive up, gram negative down
DOWN = ((2, 8), (8, 2))


def test_a_flow_is_drawn_once_in_its_net_direction() -> None:
    """#223: a 3-versus-1 split is one flow of width 2, not two flows of 3 and 1."""
    profiles: list[Profile] = []
    for study in ("bsdb:1", "bsdb:2", "bsdb:3"):
        profiles += _paired(study, *UP)
    profiles += _paired("bsdb:4", *DOWN)

    edges = {e.gram_class: e for e in build_edges(profiles, min_studies=1)}
    positive = edges["gram positive"]
    assert (positive.units_up, positive.units_down, positive.net) == (3, 1, 2)
    assert positive.direction == "inc"
    assert edges["gram negative"].direction == "dec"
    assert positive.sign_test_p == pytest.approx(float(binomtest(3, 4, p=0.5).pvalue))


def test_sign_test_q_is_corrected_within_each_gram_class() -> None:
    profiles: list[Profile] = []
    for index in range(6):
        profiles += _paired(f"a:{index}", *UP, condition="Asthma")
    for index, pattern in enumerate((UP, UP, DOWN)):
        profiles += _paired(f"b:{index}", *pattern, condition="Obesity")

    edges = build_edges(profiles, min_studies=1)
    positive = sorted(
        (e for e in edges if e.gram_class == "gram positive"), key=lambda e: e.condition
    )
    p_values = [e.sign_test_p for e in positive]
    expected = false_discovery_control(p_values, method="bh")
    assert [e.sign_test_q for e in positive] == pytest.approx([float(q) for q in expected])
    assert positive[0].sign_test_p == pytest.approx(float(binomtest(6, 6, p=0.5).pvalue))


def test_mcnemar_counts_discordant_gram_signs_per_condition_and_site() -> None:
    both_up = [_profile("bsdb:4", "inc", 6, 6), _profile("bsdb:4", "dec", 2, 2, other=8)]
    # bsdb:4 has both classes up, a concordant pair that McNemar leaves out.
    profiles = [*_paired("bsdb:1", *UP), *_paired("bsdb:2", *UP), *_paired("bsdb:3", *DOWN)]
    profiles += both_up

    edge = build_edges(profiles, min_studies=1)[0]
    assert (edge.mcnemar_pos_up_neg_down, edge.mcnemar_neg_up_pos_down) == (2, 1)
    assert edge.mcnemar_p == pytest.approx(float(binomtest(2, 3, p=0.5).pvalue))


def test_min_studies_counts_non_tied_units() -> None:
    profiles = _paired("bsdb:1", *UP) + _paired("bsdb:2", (5, 5), (5, 5))
    edges = build_edges(profiles, min_studies=1)
    assert edges
    assert all(e.units_tied == 1 for e in edges)
    assert build_edges(profiles, min_studies=2) == []


def test_an_even_split_is_recorded_without_a_direction(tmp_path: Path) -> None:
    profiles = _paired("bsdb:1", *UP) + _paired("bsdb:2", *DOWN)
    edges = build_edges(profiles, min_studies=1)
    assert all(e.direction is None for e in edges)
    path = write_edge_table(edges, tmp_path / "edges.tsv", "2026-04-24_00:41_UTC")
    header, *rows = path.read_text(encoding="utf-8").splitlines()
    assert tuple(header.split("\t")) == EDGE_HEADER
    assert len(rows) == 2
    assert all(row.split("\t")[3] == "" for row in rows)
    assert all(row.split("\t")[-1] == "2026-04-24_00:41_UTC" for row in rows)


def test_experiments_within_one_unit_are_averaged_not_counted_twice() -> None:
    profiles = _paired("bsdb:1", *UP, experiment="1") + _paired(
        "bsdb:1", (4, 6), (6, 4), experiment="2"
    )
    contrasts = unit_contrasts(profiles)
    assert list(contrasts) == [("bsdb:1", "Asthma", "Feces")]
    assert contrasts[("bsdb:1", "Asthma", "Feces")]["gram positive"] == pytest.approx(
        ((0.8 - 0.2) + (0.4 - 0.6)) / 2
    )


def test_an_unpaired_unit_raises_instead_of_dividing_by_a_wrong_count() -> None:
    """#222: the shared helper's paired guard applies here rather than a comment."""
    profiles = [*_paired("bsdb:1", *UP), _profile("bsdb:1", "inc", 5, 5, experiment="2")]
    with pytest.raises(ValueError, match="paired"):
        unit_contrasts(profiles)


def test_zero_total_and_unlabelled_experiments_are_removed_whole_and_counted() -> None:
    """#222: a zero-total side would otherwise dilute its study's contrast toward zero."""
    profiles = [
        *_paired("bsdb:1", *UP),
        _profile("bsdb:2", "inc", 0, 0),
        _profile("bsdb:2", "dec", 3, 1),
        *_paired("bsdb:3", *UP, condition=""),
    ]
    filtered = drop_undefined_pairs(profiles)
    assert filtered.zero_total_experiments == 1
    assert filtered.missing_metadata_experiments == 1
    assert {p.study_id for p in filtered.profiles} == {"bsdb:1"}


def test_complementarity_counts_opposite_and_same_signs() -> None:
    profiles = [
        *_paired("bsdb:1", *UP),
        _profile("bsdb:2", "inc", 6, 6),
        _profile("bsdb:2", "dec", 2, 2, other=8),
        *_paired("bsdb:3", (5, 5), (5, 5)),
    ]
    balance = complementarity(unit_contrasts(profiles))
    assert (balance.units, balance.opposite, balance.same, balance.either_tied) == (3, 1, 1, 1)
