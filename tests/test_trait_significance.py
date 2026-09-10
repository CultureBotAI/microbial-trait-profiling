"""Tests for the two arms' significance calls behind the heatmap borders (#172)."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest
from scipy.stats import hypergeom

from bsdbng.trait_guilds import Profile
from bsdbng.trait_significance import (
    BOTH,
    COMPOSITION,
    ENRICHMENT,
    EnrichmentCall,
    build_overlay,
    composition_calls,
    enrichment_calls,
    write_significance_tables,
)

TRAITS = {"trait_a": ("METPO:1",), "trait_b": ("METPO:2",)}


def write_profiles(path: Path, rows: list[tuple[str, str, str, int, int, int]]) -> Path:
    """Write a minimal profile matrix: two named traits plus one filler column."""
    lines = ["study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2\tMETPO:9"]
    lines.extend("\t".join(str(value) for value in row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def profile(study: str, experiment: str, direction: str, a: int, b: int, filler: int) -> Profile:
    return Profile(
        study_id=study,
        experiment_id=experiment,
        inc_dec=direction,
        condition="Colorectal cancer",
        body_site="Feces",
        guild_counts={"trait_a": a, "trait_b": b},
        total=a + b + filler,
    )


def test_enrichment_matches_the_hypergeometric_it_transcribes(tmp_path: Path) -> None:
    """The p-value is scipy's one-tailed survival function on the corpus pool."""
    rows = [
        ("s1", "e1", "inc", 90, 5, 5),
        ("s1", "e1", "dec", 5, 90, 5),
        ("s2", "e2", "inc", 10, 10, 80),
        ("s2", "e2", "dec", 10, 10, 80),
    ]
    path = write_profiles(tmp_path / "profiles.tsv", rows)
    keys = [(row[0], row[1], row[2]) for row in rows]
    calls, untested, profiles_tested, tests = enrichment_calls(keys, TRAITS, path)

    assert untested == []
    assert profiles_tested == 4
    # Two named traits plus the filler, in every one of the four profiles.
    assert tests == 12

    background_total = sum(sum(row[3:]) for row in rows)
    successes = sum(row[3] for row in rows)
    first = next(
        c for c in calls if c.study_id == "s1" and c.inc_dec == "inc" and c.trait == "trait_a"
    )
    expected_p = float(hypergeom.sf(90 - 1, background_total, successes, 100))
    assert first.p_value == pytest.approx(expected_p)
    assert first.fold == pytest.approx(90 / (100 * successes / background_total))


def test_a_profile_over_represented_in_a_trait_is_called_enriched() -> None:
    call = EnrichmentCall(
        study_id="s1",
        experiment_id="e1",
        inc_dec="inc",
        trait="trait_a",
        p_value=0.0,
        q_value=0.0,
        fold=4.0,
        count=90,
        expected=22.5,
    )
    assert call.verdict == "enriched"
    # Same evidence, but the fold sits between the two cut-offs.
    assert call.__class__(**{**call.__dict__, "fold": 1.0}).verdict is None
    assert call.__class__(**{**call.__dict__, "fold": 0.2}).verdict == "depleted"
    assert call.__class__(**{**call.__dict__, "q_value": 0.5}).verdict is None


def test_composition_needs_three_study_pairs_before_it_says_anything() -> None:
    """Two pairs is below the arm's own floor, so no contrast is emitted."""
    profiles = [
        profile("s1", "e1", "inc", 90, 5, 5),
        profile("s1", "e1", "dec", 5, 90, 5),
        profile("s2", "e2", "inc", 80, 10, 10),
        profile("s2", "e2", "dec", 10, 80, 10),
    ]
    bins = {"Colorectal cancer": "Colorectal cancer"}, {"Feces": "Feces"}
    calls = composition_calls(
        profiles, ["trait_a", "trait_b"], bin_condition=bins[0], bin_body_site=bins[1]
    )
    assert calls == []


def test_composition_reports_both_grainings_once_there_are_enough_studies() -> None:
    profiles = []
    # The contrasts have to vary across studies: four identical pairs give a
    # zero-variance column, which the arm's own t-test refuses.
    for index, weight in enumerate((90, 70, 80, 60)):
        study = f"s{index}"
        profiles.append(profile(study, f"e{index}", "inc", weight, 5, 100 - weight))
        profiles.append(profile(study, f"e{index}", "dec", 5, weight, 100 - weight))
    calls = composition_calls(
        profiles,
        ["trait_a", "trait_b"],
        bin_condition={"Colorectal cancer": "Colorectal cancer"},
        bin_body_site={"Feces": "Feces"},
    )
    kinds = {call.stratum_kind for call in calls}
    assert kinds == {"body_site", "condition_body_site"}
    site_call = next(c for c in calls if c.stratum_kind == "body_site" and c.trait == "trait_a")
    assert site_call.studies == 4
    assert site_call.condition is None


def test_a_profile_outside_the_drawn_labels_is_not_scored_against_a_cell(tmp_path: Path) -> None:
    """A body site the figure never draws contributes to no cell at all."""
    rows = [
        ("s1", "e1", "inc", 90, 5, 5),
        ("s1", "e1", "dec", 5, 90, 5),
    ]
    path = write_profiles(tmp_path / "profiles.tsv", rows)
    profiles = [
        profile("s1", "e1", "inc", 90, 5, 5),
        profile("s1", "e1", "dec", 5, 90, 5),
    ]
    overlay = build_overlay(
        profiles,
        TRAITS,
        bin_condition={"Colorectal cancer": "Colorectal cancer"},
        # Feces is absent, so nothing maps to a drawn body site.
        bin_body_site={},
        profile_path=path,
    )
    assert overlay.marks == {}
    assert overlay.support == {}


def test_one_study_never_speaks_for_a_cell(tmp_path: Path) -> None:
    """A cell backed by a single study is left unmarked however strong it is."""
    rows = [
        ("s1", "e1", "inc", 990, 5, 5),
        ("s1", "e1", "dec", 5, 990, 5),
        ("s2", "e2", "inc", 5, 5, 990),
        ("s2", "e2", "dec", 5, 5, 990),
    ]
    path = write_profiles(tmp_path / "profiles.tsv", rows)
    profiles = [
        profile("s1", "e1", "inc", 990, 5, 5),
        profile("s1", "e1", "dec", 5, 990, 5),
    ]
    overlay = build_overlay(
        profiles,
        TRAITS,
        bin_condition={"Colorectal cancer": "Colorectal cancer"},
        bin_body_site={"Feces": "Feces"},
        profile_path=path,
    )
    assert all(studies == 1 for studies, _called in overlay.support.values())
    assert overlay.marks == {}


def test_the_written_tables_carry_every_cell_and_its_support(tmp_path: Path) -> None:
    rows = [("s1", "e1", "inc", 90, 5, 5), ("s1", "e1", "dec", 5, 90, 5)]
    path = write_profiles(tmp_path / "profiles.tsv", rows)
    profiles = [
        profile("s1", "e1", "inc", 90, 5, 5),
        profile("s1", "e1", "dec", 5, 90, 5),
    ]
    overlay = build_overlay(
        profiles,
        TRAITS,
        bin_condition={"Colorectal cancer": "Colorectal cancer"},
        bin_body_site={"Feces": "Feces"},
        profile_path=path,
    )
    cells, contrasts = write_significance_tables(
        overlay, tmp_path / "cells.tsv", tmp_path / "contrasts.tsv"
    )
    header, *body = cells.read_text(encoding="utf-8").splitlines()
    assert header.split("\t")[:5] == [
        "condition",
        "body_site",
        "trait",
        "direction",
        "mark",
    ]
    # Two traits by two directions.
    assert len(body) == 4
    assert contrasts.read_text(encoding="utf-8").startswith("stratum_kind\t")


def test_the_three_border_styles_are_distinct() -> None:
    assert len({ENRICHMENT, COMPOSITION, BOTH}) == 3


@pytest.mark.integration
def test_reproduces_the_enrichment_arm_committed_counts() -> None:
    """Orthogonal check against the table this test's subject was transcribed from.

    The enrichment arm's own output on `origin/enrichment-analysis-heatmap-cluster`
    was produced by a different script, on a different machine, from the same
    corpus. If this module's hypergeometric is the arm's hypergeometric, the two
    have to agree on how many profiles each Table 1 trait is significant in.
    """
    from bsdbng.trait_guilds import (
        DEFAULT_EXPERIMENTS_PATH,
        DEFAULT_PROFILE_PATH,
        TABLE1_TRAITS,
        load_profiles,
    )

    committed = subprocess.run(
        [
            "git",
            "show",
            "origin/enrichment-analysis-heatmap-cluster:data/derived/trait_enrichment_results.tsv",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if committed.returncode != 0:
        pytest.skip("the enrichment arm's branch is not fetched")

    theirs: dict[str, int] = {}
    for line in committed.stdout.splitlines():
        if line.startswith(("#", "study_id")):
            continue
        trait = line.split("\t")[3]
        theirs[trait] = theirs.get(trait, 0) + 1

    profiles, _ = load_profiles(DEFAULT_PROFILE_PATH, DEFAULT_EXPERIMENTS_PATH, TABLE1_TRAITS)
    keys = [(p.study_id, p.experiment_id, p.inc_dec) for p in profiles]
    calls, _untested, _tested, _tests = enrichment_calls(keys, TABLE1_TRAITS, DEFAULT_PROFILE_PATH)

    ours: dict[str, int] = {}
    for call in calls:
        if call.verdict == "enriched":
            ours[call.trait] = ours.get(call.trait, 0) + 1

    for name, ids in TABLE1_TRAITS.items():
        term = ids[0]
        theirs_count = next(
            (count for key, count in theirs.items() if key.rsplit("|", 1)[-1] == term), None
        )
        assert theirs_count is not None, f"{term} absent from the committed table"
        assert ours.get(name, 0) == theirs_count, name
