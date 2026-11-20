"""Gram-stain directional Sankey: paired contrasts, sign tests, McNemar, edge table (#147, #223)."""

import csv
import json
from pathlib import Path

import pytest

from mtp.gram_sankey import (
    EDGES_FILENAME,
    GRAM_NEGATIVE,
    GRAM_POSITIVE,
    MCNEMAR_FAMILY_ID,
    STATUS_BELOW_MIN,
    STATUS_BELOW_MIN_DISCORDANT,
    STATUS_TESTED,
    SUMMARY_FILENAME,
    Complementarity,
    build_flows,
    build_parser,
    complementarity,
    drop_unplaced_experiments,
    figure_subtitle,
    main,
    sign_family_id,
    site_contrasts,
)
from mtp.guild_maps import exact_mcnemar_p
from mtp.trait_guilds import DEFAULT_PROFILE_PATH, GRAM_GUILDS, Profile

POSITIVE_TERM, NEGATIVE_TERM = "METPO:1000698", "METPO:1000699"
OTHER_TERM = "METPO:1000615"
GRAM_TERMS = (POSITIVE_TERM, NEGATIVE_TERM, OTHER_TERM)
#: Two-sided exact binomial p for 9 of 10 against one half: 2 * (1 + 10) / 2**10.
P_9_OF_10 = 22 / 1024
#: BH over [P_9_OF_10, 1.0]: the smaller p times 2 / 1.
BH_9_OF_10_IN_TWO = 2 * P_9_OF_10
#: BY multiplies BH by the harmonic number of the family size, 1 + 1/2.
BY_9_OF_10_IN_TWO = 1.5 * BH_9_OF_10_IN_TWO

#: One study per entry: (condition, body site, gram positive shifts up).
STUDIES = (
    *[("Colorectal cancer", "Feces", True)] * 9,
    ("Colorectal cancer", "Feces", False),
    *[("Asthma", "Feces", True)] * 2,
    ("Asthma", "Feces", False),
    *[("Autism", "Feces", True)] * 2,
)


def _sides(positive_up: bool) -> dict[str, list[int]]:
    """Gram positive, gram negative, other counts per side; totals are equal."""
    up, down = [6, 2, 2], [2, 6, 2]
    return {"inc": up, "dec": down} if positive_up else {"inc": down, "dec": up}


def _profiles() -> list[Profile]:
    profiles = []
    for index, (condition, body_site, positive_up) in enumerate(STUDIES):
        for inc_dec, (pos, neg, other) in _sides(positive_up).items():
            profiles.append(
                Profile(
                    study_id=f"bsdb:{index}",
                    experiment_id=f"bsdb:{index}-1",
                    inc_dec=inc_dec,
                    condition=condition,
                    body_site=body_site,
                    guild_counts={GRAM_POSITIVE: pos, GRAM_NEGATIVE: neg},
                    total=pos + neg + other,
                )
            )
    return profiles


def _write_corpus(tmp_path: Path) -> tuple[Path, Path]:
    profile_path = tmp_path / "profiles.tsv"
    experiments_path = tmp_path / "experiments.tsv"
    profile_lines = ["\t".join(["study_id", "experiment_id", "inc_dec", *GRAM_TERMS])]
    experiment_lines = [
        "\t".join(["study_id", "experiment_id", "condition", "body_site", "host_species"])
    ]
    for index, (condition, body_site, positive_up) in enumerate(STUDIES):
        study, experiment = f"bsdb:{index}", f"bsdb:{index}-1"
        for inc_dec, counts in _sides(positive_up).items():
            profile_lines.append("\t".join([study, experiment, inc_dec, *map(str, counts)]))
        experiment_lines.append(
            "\t".join([study, experiment, condition, body_site, "Homo sapiens"])
        )
    profile_path.write_text("\n".join(profile_lines) + "\n", encoding="utf-8")
    experiments_path.write_text("\n".join(experiment_lines) + "\n", encoding="utf-8")
    return profile_path, experiments_path


def test_gram_classes_read_the_two_gram_stain_terms() -> None:
    assert GRAM_GUILDS == {GRAM_POSITIVE: (POSITIVE_TERM,), GRAM_NEGATIVE: (NEGATIVE_TERM,)}


def test_default_input_is_the_shared_human_16s_shotgun_profile_table() -> None:
    args = build_parser().parse_args([])
    assert args.profile_input == DEFAULT_PROFILE_PATH


def test_exact_mcnemar_matches_the_binomial_on_discordant_pairs() -> None:
    assert exact_mcnemar_p(9, 1) == pytest.approx(P_9_OF_10)
    assert exact_mcnemar_p(6, 0) == pytest.approx(2 / 64)
    assert exact_mcnemar_p(0, 0) is None


def test_sign_tests_and_mcnemar_match_known_values() -> None:
    flows = {(f.gram_class, f.condition): f for f in build_flows(site_contrasts(_profiles()))}

    colorectal = flows[(GRAM_POSITIVE, "Colorectal cancer")]
    assert (colorectal.studies_up, colorectal.studies_down, colorectal.studies_tied) == (9, 1, 0)
    assert colorectal.sign_p_nominal == pytest.approx(P_9_OF_10)
    assert colorectal.sign_bh_q == pytest.approx(BH_9_OF_10_IN_TWO)
    assert colorectal.sign_by_q == pytest.approx(BY_9_OF_10_IN_TWO)
    assert colorectal.direction == "inc"
    assert colorectal.significant
    assert flows[(GRAM_NEGATIVE, "Colorectal cancer")].direction == "dec"

    asthma = flows[(GRAM_POSITIVE, "Asthma")]
    assert asthma.sign_p_nominal == pytest.approx(1.0)
    assert asthma.sign_bh_q == pytest.approx(1.0)
    assert not asthma.significant

    mcnemar = colorectal.mcnemar
    assert (mcnemar.positive_up_negative_down, mcnemar.negative_up_positive_down) == (9, 1)
    assert mcnemar.p_nominal == pytest.approx(P_9_OF_10)
    assert mcnemar.bh_q == pytest.approx(BH_9_OF_10_IN_TWO)
    assert (mcnemar.family_size, mcnemar.unfiltered_size) == (2, 3)


def test_flows_below_the_support_floor_are_reported_untested_outside_the_family() -> None:
    flows = build_flows(site_contrasts(_profiles()), min_studies=3)
    autism = next(f for f in flows if f.gram_class == GRAM_POSITIVE and f.condition == "Autism")
    assert autism.status == STATUS_BELOW_MIN
    assert autism.sign_p_nominal is None
    assert autism.sign_bh_q is None
    assert autism.direction is None
    assert (autism.sign_family_size, autism.sign_unfiltered_size) == (2, 3)
    assert autism.mcnemar.status == STATUS_BELOW_MIN_DISCORDANT
    assert autism.mcnemar.bh_q is None


def test_one_study_votes_once_per_flow_even_with_repeated_experiments() -> None:
    profiles = _profiles()
    repeat = [
        Profile(
            study_id=p.study_id,
            experiment_id=f"{p.experiment_id}-repeat",
            inc_dec=p.inc_dec,
            condition=p.condition,
            body_site=p.body_site,
            guild_counts=p.guild_counts,
            total=p.total,
        )
        for p in profiles
        if p.study_id == "bsdb:0"
    ]
    flows = build_flows(site_contrasts(profiles + repeat))
    colorectal = next(
        f for f in flows if f.gram_class == GRAM_POSITIVE and f.condition == "Colorectal cancer"
    )
    assert colorectal.studies_scored == 10


def test_an_unpaired_study_is_refused_rather_than_averaged() -> None:
    profiles = _profiles()
    extra = next(p for p in profiles if p.study_id == "bsdb:0" and p.inc_dec == "inc")
    unpaired = Profile(
        study_id=extra.study_id,
        experiment_id="bsdb:0-2",
        inc_dec="inc",
        condition=extra.condition,
        body_site=extra.body_site,
        guild_counts=extra.guild_counts,
        total=extra.total,
    )
    with pytest.raises(ValueError, match="paired"):
        site_contrasts([*profiles, unpaired])


def test_unplaced_experiments_are_dropped_whole() -> None:
    profiles = _profiles()
    first = profiles[0]
    blank = [
        Profile(
            study_id="bsdb:x",
            experiment_id="bsdb:x-1",
            inc_dec=side,
            condition="",
            body_site="Feces",
            guild_counts=first.guild_counts,
            total=first.total,
        )
        for side in ("inc", "dec")
    ]
    placed = drop_unplaced_experiments([*profiles, *blank])
    assert placed.missing_site_experiments == 1
    assert placed.profiles == profiles


def test_complementarity_counts_opposite_gram_shifts() -> None:
    balance = complementarity(site_contrasts(_profiles()))
    assert balance == Complementarity(
        units=len(STUDIES), opposite=len(STUDIES), same=0, either_tied=0
    )


def test_subtitle_says_when_no_flow_passes() -> None:
    flows = build_flows(site_contrasts(_profiles()))
    failing = [f for f in flows if f.gram_class == GRAM_POSITIVE and f.condition == "Asthma"]
    text = figure_subtitle(failing, GRAM_POSITIVE, "inc", Complementarity(3, 2, 1, 0), "v")
    assert "No flow passes its sign test" in text
    assert sign_family_id(GRAM_POSITIVE) in text
    assert "not abundance" in text


def test_main_writes_the_statistics_into_the_edge_table(tmp_path: Path) -> None:
    profile_path, experiments_path = _write_corpus(tmp_path)
    output_dir = tmp_path / "out"
    assert (
        main(
            [
                "--profile-input",
                str(profile_path),
                "--experiments-input",
                str(experiments_path),
                "--full-dump-input",
                str(tmp_path / "absent.csv"),
                "--output-dir",
                str(output_dir),
                "--no-figures",
            ]
        )
        == 0
    )
    with (output_dir / EDGES_FILENAME).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 6
    row = next(
        r
        for r in rows
        if r["gram_class"] == GRAM_POSITIVE and r["condition"] == "Colorectal cancer"
    )
    assert row["status"] == STATUS_TESTED
    assert row["drawn_direction"] == "increased"
    assert (row["studies_up"], row["studies_down"], row["net_studies"]) == ("9", "1", "8")
    assert float(row["sign_test_p_nominal"]) == pytest.approx(P_9_OF_10, rel=1e-5)
    assert float(row["sign_test_bh_q"]) == pytest.approx(BH_9_OF_10_IN_TWO, rel=1e-5)
    assert float(row["sign_test_by_q"]) == pytest.approx(BY_9_OF_10_IN_TWO, rel=1e-5)
    assert row["sign_test_family_id"] == sign_family_id(GRAM_POSITIVE)
    assert (row["sign_test_family_size"], row["sign_test_unfiltered_flows"]) == ("2", "3")
    assert (
        row["mcnemar_positive_up_negative_down"],
        row["mcnemar_negative_up_positive_down"],
    ) == ("9", "1")
    assert float(row["mcnemar_p_nominal"]) == pytest.approx(P_9_OF_10, rel=1e-5)
    assert float(row["mcnemar_bh_q"]) == pytest.approx(BH_9_OF_10_IN_TWO, rel=1e-5)
    assert row["mcnemar_family_id"] == MCNEMAR_FAMILY_ID
    assert row["primary_correction"] == "BH"

    untested = next(r for r in rows if r["condition"] == "Autism")
    assert untested["status"] == STATUS_BELOW_MIN
    assert untested["sign_test_p_nominal"] == untested["sign_test_bh_q"] == ""

    summary = json.loads((output_dir / SUMMARY_FILENAME).read_text(encoding="utf-8"))
    assert summary["sign_test"][GRAM_POSITIVE]["bh_q_passing"] == 1
    assert summary["complementarity"]["opposite"] == len(STUDIES)
