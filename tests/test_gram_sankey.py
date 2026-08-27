"""Gram-stain Sankey edge construction."""

from pathlib import Path

import pytest

from bsdbng.gram_sankey import build_edges, write_edge_table
from bsdbng.trait_guilds import GRAM_GUILDS, Profile

GRAM_POS = GRAM_GUILDS["gram positive"][0]
GRAM_NEG = GRAM_GUILDS["gram negative"][0]


def _profile(
    study: str, direction: str, pos: int, neg: int, *, condition: str = "Asthma"
) -> Profile:
    return Profile(
        study_id=study,
        experiment_id=f"{study}-1",
        inc_dec=direction,
        condition=condition,
        body_site="Feces",
        guild_counts={"gram positive": pos, "gram negative": neg},
        total=pos + neg,
    )


def _paired(study: str, inc_pos: int, inc_neg: int, dec_pos: int, dec_neg: int) -> list[Profile]:
    return [
        _profile(study, "inc", inc_pos, inc_neg),
        _profile(study, "dec", dec_pos, dec_neg),
    ]


def test_edge_direction_follows_the_paired_contrast_sign() -> None:
    # Gram positive is a larger share on the increased side in all three studies.
    profiles: list[Profile] = []
    for study in ("bsdb:1", "bsdb:2", "bsdb:3"):
        profiles += _paired(study, inc_pos=8, inc_neg=2, dec_pos=2, dec_neg=8)

    edges = build_edges(profiles, min_studies=1)
    by_key = {(e.gram_class, e.direction): e for e in edges}

    assert by_key[("gram positive", "inc")].studies == 3
    assert by_key[("gram negative", "dec")].studies == 3
    assert ("gram positive", "dec") not in by_key
    assert ("gram negative", "inc") not in by_key

    # Contrast magnitude is 0.8 - 0.2 for each study.
    assert by_key[("gram positive", "inc")].mean_contrast == pytest.approx(0.6)


def test_min_studies_drops_thin_flows() -> None:
    profiles = _paired("bsdb:1", 9, 1, 1, 9)
    assert build_edges(profiles, min_studies=1)
    assert build_edges(profiles, min_studies=2) == []


def test_zero_contrast_contributes_no_edge() -> None:
    # Identical fractions on both sides: no directional claim to make.
    profiles = _paired("bsdb:1", 5, 5, 5, 5)
    assert build_edges(profiles, min_studies=1) == []


def test_profiles_without_metadata_are_skipped() -> None:
    profiles = [
        _profile("bsdb:1", "inc", 8, 2, condition=""),
        _profile("bsdb:1", "dec", 2, 8, condition=""),
    ]
    assert build_edges(profiles, min_studies=1) == []


def test_edge_table_records_direction_labels_and_vintage(tmp_path: Path) -> None:
    profiles = _paired("bsdb:1", 8, 2, 2, 8) + _paired("bsdb:2", 8, 2, 2, 8)
    edges = build_edges(profiles, min_studies=1)
    path = write_edge_table(edges, tmp_path / "edges.tsv", "2026-04-24_00:41_UTC")

    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == [
        "gram_class",
        "direction",
        "condition",
        "body_site",
        "studies",
        "mean_abs_contrast",
        "bugsigdb_vintage",
    ]
    body = [line.split("\t") for line in lines[1:]]
    assert {row[1] for row in body} == {"increased", "decreased"}
    assert all(row[6] == "2026-04-24_00:41_UTC" for row in body)
