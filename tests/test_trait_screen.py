"""The widened directional family: its correction, its dependence, its missing states."""

import csv
import json
import random
from pathlib import Path

import pytest
from scipy.stats import false_discovery_control

from bsdbng.table1_trait_maps import compute_direction_shifts
from bsdbng.trait_guilds import TABLE1_TRAITS, Profile, load_profiles, retained_trait_columns
from bsdbng.trait_screen import (
    ALL_TIED,
    NO_ELIGIBLE_STUDIES,
    TESTED,
    Dependence,
    candidate_slots,
    cheverud_note,
    dependence,
    effective_family_size,
    main,
    run_screen,
    screen_candidates,
    significant,
    step_up,
    table1_slug,
)


def _profile(
    study: str,
    direction: str,
    counts: dict[str, int],
    *,
    experiment: str = "e1",
) -> Profile:
    return Profile(study, experiment, direction, "X", "Feces", counts, sum(counts.values()) or 1)


def _corpus(path: Path, columns: list[str], rows: list[tuple[str, str, str, list[int]]]) -> None:
    header = ["study_id", "experiment_id", "inc_dec", *columns, "mediadive.medium:X"]
    lines = ["\t".join(header)]
    for study, experiment, direction, values in rows:
        lines.append("\t".join([study, experiment, direction, *map(str, values), "7"]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# --- the correction --------------------------------------------------------


@pytest.mark.parametrize("method", ["bh", "by"])
def test_step_up_at_the_nominal_family_size_is_exactly_scipy(method: str) -> None:
    """The adjustable family size is the only thing this adds to scipy's procedure.

    #265 needs the same step-up run at an effective m, which
    ``false_discovery_control`` cannot express, so it is written here. It must
    reproduce scipy exactly whenever the family is the number of tests.
    """
    rng = random.Random(255)
    for size in (1, 2, 5, 40, 317):
        p_values = [rng.random() for _ in range(size)]
        ours = step_up(p_values, method=method, family_size=size)
        theirs = [float(q) for q in false_discovery_control(p_values, method=method)]
        assert ours == pytest.approx(theirs)


def test_a_smaller_family_is_less_conservative_but_never_below_the_raw_p() -> None:
    p_values = [0.001, 0.02, 0.2, 0.4, 0.9]
    nominal = step_up(p_values, method="by", family_size=len(p_values))
    effective = step_up(p_values, method="by", family_size=2)
    assert all(small <= large for small, large in zip(effective, nominal, strict=True))
    # An adjusted p-value below its own raw p-value would be meaningless.
    assert all(q >= p for q, p in zip(effective, p_values, strict=True))


@pytest.mark.parametrize(("method", "family_size"), [("bogus", 3), ("bh", 0), ("by", -1)])
def test_an_undeclarable_correction_is_rejected(method: str, family_size: int) -> None:
    with pytest.raises(ValueError, match=r"method|family_size"):
        step_up([0.1, 0.2, 0.3], method=method, family_size=family_size)


# --- the effective family size ---------------------------------------------


def test_effective_family_size_spans_independence_to_redundancy() -> None:
    """Every estimator must agree at the two limits, or none of them means anything."""
    independent = effective_family_size([1.0] * 10)
    assert independent == pytest.approx((10.0, 10.0, 10.0))
    # One perfectly correlated block: all the variance in a single component.
    redundant = effective_family_size([10.0, *[0.0] * 9])
    assert redundant == pytest.approx((1.0, 1.0, 1.0))


def test_cheverud_nyholt_is_computed_so_that_it_can_be_named_and_excluded() -> None:
    """#265 asks for the excluded estimator to be stated, not silently dropped."""
    eigenvalues = [5.0, 2.0, 1.5, 0.9, 0.4, 0.2]
    li_ji, galwey, cheverud = effective_family_size(eigenvalues)
    for estimate in (li_ji, galwey, cheverud):
        assert 1.0 <= estimate <= len(eigenvalues)
    # Galwey is the most aggressive of the three at discounting shared structure.
    assert galwey <= li_ji


# --- the candidate family --------------------------------------------------


def test_the_candidate_family_is_the_retained_universe_and_excludes_mediadive(
    tmp_path: Path,
) -> None:
    path = tmp_path / "profiles.tsv"
    _corpus(path, ["A", "B"], [("s1", "e1", "inc", [1, 2]), ("s1", "e1", "dec", [2, 1])])
    assert screen_candidates(path) == ("A", "B")
    header = path.open(encoding="utf-8").readline().rstrip("\n").split("\t")
    positions, columns, excluded = retained_trait_columns(header)
    assert columns == ["A", "B"]
    assert excluded == 1
    assert len(positions) == 2


def test_each_candidate_is_its_own_single_column_slot() -> None:
    assert candidate_slots(["A", "B"]) == {"A": ("A",), "B": ("B",)}


def test_the_six_table_1_traits_resolve_through_their_qualifiers() -> None:
    """Two of the six sit behind a METPO qualifier; a naive membership test loses them."""
    assert table1_slug("METPO:1000608") == "facultatively_aerobic"
    assert table1_slug("METPO:2000103|GO:0030245") == "cellulose_catabolic_process"
    assert table1_slug("METPO:1000602") is None


# --- the statistic and its missing states ----------------------------------


def _six_study_corpus() -> list[Profile]:
    """Six studies: A shifts up, B shifts down, C is tied at zero on both sides.

    Both sides carry the same total, because these are fractions: an unequal
    total moves every trait's fraction, so a trait with identical counts on both
    sides is only tied when the denominators match too.
    """
    profiles = []
    for index in range(6):
        profiles.append(_profile(f"s{index}", "inc", {"A": 8, "B": 5, "C": 0}))
        profiles.append(_profile(f"s{index}", "dec", {"A": 5, "B": 8, "C": 0}))
    return profiles


def test_the_widened_screen_votes_exactly_as_the_six_trait_arm_does() -> None:
    """#255 changes the family, not the statistic. The tallies must be identical."""
    profiles = _six_study_corpus()
    columns = ["A", "B", "C"]
    screen = run_screen(profiles, columns)
    legacy = compute_direction_shifts(profiles, columns)
    by_trait = {shift.trait: shift for shift in legacy.shifts}
    for row in screen.rows:
        shift = by_trait[row.trait]
        assert (row.studies_up, row.studies_down, row.studies_tied) == (
            shift.studies_up,
            shift.studies_down,
            shift.studies_tied,
        )
        if row.test_status == TESTED:
            assert row.p_value == pytest.approx(shift.p_value)


def test_the_six_study_corpus_gives_the_hand_computed_sign_test() -> None:
    """Pins the shared rule itself, which comparing the two arms cannot (#330).

    Both arms call the same tally, so agreement between them says nothing about
    whether that tally is right. By hand: A is 8/13 vs 5/13 in every study, so
    6 up, 0 down, and the two-sided exact binomial p is 2 * 0.5**6 = 1/32; B is
    the mirror; C is 0 vs 0 everywhere, so 6 tied and untested. With C held at
    p = 1 in a family of 3, BH gives A and B q = (1/32) * 3/2 = 3/64.
    """
    rows = {row.trait: row for row in run_screen(_six_study_corpus(), ["A", "B", "C"]).rows}
    assert (rows["A"].studies_up, rows["A"].studies_down, rows["A"].studies_tied) == (6, 0, 0)
    assert (rows["B"].studies_up, rows["B"].studies_down, rows["B"].studies_tied) == (0, 6, 0)
    assert (rows["C"].studies_up, rows["C"].studies_down, rows["C"].studies_tied) == (0, 0, 6)
    assert rows["A"].p_value == pytest.approx(1 / 32)
    assert rows["B"].p_value == pytest.approx(1 / 32)
    assert rows["A"].bh_q_nominal == pytest.approx(3 / 64)
    assert rows["B"].bh_q_nominal == pytest.approx(3 / 64)
    # Exact Clopper-Pearson for 6 of 6 is (0.025 ** (1 / 6), 1); testing the
    # down count instead of the up count would put it at (0, 0.459) (#415).
    assert rows["A"].ci_low == pytest.approx(0.025 ** (1 / 6))
    assert rows["A"].ci_high == pytest.approx(1.0)


def test_ties_are_left_out_of_the_binomial_n() -> None:
    """The manifest's test is of non-tied studies; a unanimous fixture cannot show it (#415).

    M shifts up in four studies, down in one, and is tied in one: p is the
    two-sided exact binomial of 4 of 5 = 0.375. Counting the tie as a trial
    would give 4 of 6 = 0.6875.
    """
    up = ({"M": 2, "Z": 1}, {"M": 1, "Z": 2})
    down = ({"M": 1, "Z": 2}, {"M": 2, "Z": 1})
    tied = ({"M": 1, "Z": 2}, {"M": 1, "Z": 2})
    profiles = []
    for index, (inc, dec) in enumerate([up, up, up, up, down, tied]):
        profiles.append(_profile(f"s{index}", "inc", inc))
        profiles.append(_profile(f"s{index}", "dec", dec))
    rows = {row.trait: row for row in run_screen(profiles, ["M", "Z"]).rows}
    assert (rows["M"].studies_up, rows["M"].studies_down, rows["M"].studies_tied) == (4, 1, 1)
    assert rows["M"].p_value == pytest.approx(0.375)
    assert rows["M"].share_up == pytest.approx(0.8)


def test_an_all_tied_trait_has_no_result_rather_than_a_share_of_zero() -> None:
    screen = run_screen(_six_study_corpus(), ["A", "B", "C"])
    rows = {row.trait: row for row in screen.rows}
    assert rows["C"].test_status == ALL_TIED
    assert rows["C"].studies_tied == 6
    assert rows["C"].share_up is None
    assert rows["C"].p_value is None
    assert rows["C"].by_q_nominal is rows["C"].by_q_effective is rows["C"].bh_q_nominal is None
    assert rows["A"].test_status == TESTED
    assert rows["A"].share_up == pytest.approx(1.0)
    assert rows["B"].share_up == pytest.approx(0.0)


def test_a_trait_zero_on_both_sides_is_all_tied_not_a_null_result() -> None:
    profiles = [
        _profile("s0", "inc", {"A": 1, "D": 0}),
        _profile("s0", "dec", {"A": 2, "D": 0}),
    ]
    # D is zero on both sides, so the study's contrast for it is exactly tied:
    # a paired study always contributes a contrast for every candidate (#327).
    screen = run_screen(profiles, ["A", "D"])
    rows = {row.trait: row for row in screen.rows}
    assert rows["D"].test_status == ALL_TIED
    assert rows["D"].studies_tied == 1
    assert rows["D"].p_value is None
    assert rows["D"].share_up is None


def test_no_eligible_studies_is_the_state_of_a_scope_with_no_paired_study() -> None:
    """Every paired study contrasts every candidate, so this state is per-scope (#327).

    It cannot single out one trait; it is what every candidate reads when the
    profiles passed in hold no paired study at all. The CLI has no subset
    option, so on the real pipeline this means load_profiles kept no paired
    human study.
    """
    screen = run_screen([], ["A", "D"])
    assert {row.test_status for row in screen.rows} == {NO_ELIGIBLE_STUDIES}
    assert all(row.p_value is None and row.share_up is None for row in screen.rows)


def test_untestable_candidates_still_count_toward_the_family_they_are_corrected_in() -> None:
    """Holding them at p = 1 is PR #181's convention; dropping them would shrink the family."""
    screen = run_screen(_six_study_corpus(), ["A", "B", "C"])
    assert screen.family_size == 3
    tested = [row for row in screen.rows if row.test_status == TESTED]
    assert len(tested) == 2
    assert all(row.by_q_nominal is not None for row in tested)
    # Corrected in a family of three, not of the two that happened to be testable.
    assert tested[0].by_q_nominal == pytest.approx(
        step_up(
            [tested[0].p_value or 1.0, tested[1].p_value or 1.0, 1.0], method="by", family_size=3
        )[0]
    )


# --- dependence ------------------------------------------------------------


def test_dependence_keeps_a_constant_candidate_in_the_family_but_out_of_the_correlation() -> None:
    matrix = [
        [1.0, 2.0, 3.0, 9.0, 5.0],
        [2.0, 4.0, 1.0, 3.0, 5.0],
        [3.0, 1.0, 4.0, 7.0, 5.0],
        [4.0, 3.0, 2.0, 1.0, 5.0],
    ]
    record = dependence(matrix, "test")
    assert record is not None
    assert record.constant_traits_excluded == 1
    assert record.traits == 4
    assert record.units == 4


def test_a_corpus_too_small_for_dependence_leaves_the_effective_size_undefined() -> None:
    """The estimate is a diagnostic; losing it must not invent a number or fail the screen."""
    assert dependence([[1.0, 2.0, 3.0]], "test") is None
    assert dependence([[1.0, 2.0], [2.0, 1.0], [3.0, 4.0]], "test") is None
    screen = run_screen(_six_study_corpus(), ["A", "B", "C"])
    assert screen.contrast_dependence is None
    assert screen.effective_family_size is None
    assert all(row.by_q_effective is None for row in screen.rows)
    assert screen.family_size == 3


def test_the_effective_set_contains_the_nominal_one() -> None:
    """A smaller family only ever adds traits, so the difference is gained, not swapped."""
    profiles = [
        profile
        for index in range(30)
        for profile in (
            _profile(f"s{index}", "inc", {"A": 9, "B": 5 + index % 3, "C": 4}),
            _profile(f"s{index}", "dec", {"A": 1, "B": 5, "C": 6}),
        )
    ]
    screen = run_screen(profiles, ["A", "B", "C"])
    nominal = significant(screen.rows, "by_q_nominal", 0.05)
    effective = significant(screen.rows, "by_q_effective", 0.05)
    assert nominal <= effective


# --- the pinned corpus and the CLI -----------------------------------------


def test_the_cli_writes_every_candidate_with_an_explicit_state(tmp_path: Path) -> None:
    profile_path = tmp_path / "profiles.tsv"
    experiments = tmp_path / "experiments.tsv"
    dump = tmp_path / "dump.csv"
    columns = [ids[0] for ids in TABLE1_TRAITS.values()] + ["METPO:9999999"]
    rows = []
    for index in range(4):
        rows.append((f"s{index}", "e1", "inc", [9, 1, 5, 5, 3, 7, 0]))
        rows.append((f"s{index}", "e1", "dec", [1, 9, 5, 5, 3, 2, 0]))
    _corpus(profile_path, columns, rows)
    experiments.write_text(
        "study_id\texperiment_id\tcondition\tbody_site\thost_species\n"
        + "".join(f"s{i}\te1\tX\tFeces\tHomo sapiens\n" for i in range(4)),
        encoding="utf-8",
    )
    dump.write_text("# BugSigDB pinned-test, License: test\n", encoding="utf-8")
    output = tmp_path / "out"
    assert (
        main(
            [
                "--profile-input",
                str(profile_path),
                "--experiments-input",
                str(experiments),
                "--full-dump-input",
                str(dump),
                "--output-dir",
                str(output),
                "--column-labels",
                str(tmp_path / "absent.tsv"),
            ]
        )
        == 0
    )
    with (output / "trait_direction_screen.tsv").open(encoding="utf-8") as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert len(written) == len(columns)
    assert {row["family_size_nominal"] for row in written} == {str(len(columns))}
    assert all(row["test_status"] in (TESTED, ALL_TIED, NO_ELIGIBLE_STUDIES) for row in written)
    # The six are flagged as a subset, never as the family.
    assert sum(row["in_table1"] == "true" for row in written) == 6
    # An untested candidate carries no statistic at all, not a zero.
    untested = [row for row in written if row["test_status"] != TESTED]
    assert untested
    assert all(
        row["share_up"] == row["binomial_p"] == row["by_q_nominal"] == "" for row in untested
    )
    # Without the label index every candidate falls back to its own column name.
    assert all(row["trait_label"] == row["trait"] for row in written)

    provenance = json.loads((output / "trait_screen_provenance.json").read_text(encoding="utf-8"))
    assert provenance["manifest"] == "trait-direction-screen-v1"
    assert provenance["correction_primary"] == "Benjamini-Yekutieli"
    assert provenance["effective_estimator_excluded"] == "cheverud_nyholt"
    assert provenance["bugsigdb_vintage"] == "pinned-test"
    assert provenance["candidates_with_a_resolved_label"] == 0

    family = list(
        csv.DictReader((output / "trait_screen_family.tsv").open(encoding="utf-8"), delimiter="\t")
    )
    quantities = {row["quantity"] for row in family}
    assert {"family_size_nominal", "effective_li_ji", "effective_galwey"} <= quantities
    assert any(
        row["quantity"] == "effective_cheverud_nyholt" and "EXCLUDED" in row["note"]
        for row in family
    )


@pytest.mark.integration
def test_the_pinned_corpus_screens_every_retained_trait_and_reproduces_the_quoted_dependence() -> (
    None
):
    """The profile-basis figures in #263/#265 must come back out of this code.

    Li & Ji 60.0, Galwey 48.9, Cheverud/Nyholt 198.9, PC1 60.4%, median rho
    +0.611 over 50,086 pairs, 63.5% above +0.5. They were computed independently
    of this module, so reproducing them checks the correlation and spectrum code
    rather than restating it. The values below are repinned to the 16S and shotgun
    corpus of #742 (5,629 profiles), where Li & Ji is 61.0 and Galwey 48.8.
    """
    columns = screen_candidates()
    slots = candidate_slots(columns)
    paired, stats = load_profiles(guilds=slots)
    unpaired, _ = load_profiles(guilds=slots, paired_only=False)
    assert len(columns) == stats.trait_columns_used == 317
    result = run_screen(paired, columns, dependence_profiles=unpaired)
    assert result.family_size == 317
    contrast = result.contrast_dependence
    assert contrast is not None
    assert contrast.li_ji == pytest.approx(112.0, abs=0.05)
    assert contrast.galwey == pytest.approx(109.1, abs=0.05)
    assert contrast.cheverud_nyholt == pytest.approx(294.3, abs=0.05)
    # Li & Ji is 111.99999999999996 here, so a floor instead of a round would
    # silently correct at 111; the 57 count below is unchanged from 110 to 113 (#416).
    assert result.effective_family_size_used == 112
    # The contrast-basis figures the methods note and report lead with (#326).
    assert contrast.median_rho == pytest.approx(0.088, abs=0.001)
    assert contrast.pc1_share == pytest.approx(0.206, abs=0.001)
    assert contrast.share_rho_above_half == pytest.approx(0.076, abs=0.001)
    assert (
        sum(row.by_q_nominal is not None and row.by_q_nominal <= 0.05 for row in result.rows) == 38
    )
    assert (
        sum(row.by_q_effective is not None and row.by_q_effective <= 0.05 for row in result.rows)
        == 57
    )
    # Tallies and BY q at m = 317 for the six Table 1 traits, published in #276 and
    # repinned to the 16S and shotgun corpus (#742).
    expected = {
        "facultatively_aerobic": (521, 377, 24, 0.00269326),
        "moderately_halophilic": (476, 353, 93, 0.00628741),
        "cellulose_catabolic_process": (330, 433, 159, 0.0407693),
        "arsenate_reduction": (341, 265, 316, 0.0957324),
        "sphere_shaped": (463, 370, 89, 0.0679272),
        "coccus_shaped": (488, 408, 26, 0.237474),
    }
    table1 = {row.table1_slug: row for row in result.rows if row.table1_slug}
    assert set(table1) == set(expected)
    for slug, (up, down, tied, by_q) in expected.items():
        row = table1[slug]
        assert (row.studies_up, row.studies_down, row.studies_tied) == (up, down, tied), slug
        assert row.by_q_nominal == pytest.approx(by_q, rel=1e-4), slug
    profiles = result.profile_dependence
    assert profiles is not None
    assert profiles.units == 5629
    assert profiles.li_ji == pytest.approx(61.0, abs=0.05)
    assert profiles.galwey == pytest.approx(48.8, abs=0.05)
    assert profiles.cheverud_nyholt == pytest.approx(198.85, abs=0.05)
    assert profiles.pc1_share == pytest.approx(0.604, abs=0.001)
    assert profiles.pairs == 50086
    assert profiles.median_rho == pytest.approx(0.6105, abs=0.001)
    assert profiles.share_rho_above_half == pytest.approx(0.634, abs=0.001)


def _dependence_record(*, li_ji: float, galwey: float, cheverud_nyholt: float) -> Dependence:
    return Dependence(
        basis="contrast",
        units=10,
        traits=5,
        constant_traits_excluded=0,
        pc1_share=0.5,
        eigenvalues=(),
        li_ji=li_ji,
        galwey=galwey,
        cheverud_nyholt=cheverud_nyholt,
        pairs=10,
        median_rho=0.1,
        mean_rho=0.1,
        share_rho_above_half=0.0,
        share_rho_above_nine_tenths=0.0,
    )


def test_cheverud_note_claims_strictness_only_when_it_is_the_largest() -> None:
    largest = cheverud_note(_dependence_record(li_ji=3.0, galwey=3.1, cheverud_nyholt=4.5))
    assert "Li & Ji (2005) found overly large" in largest
    assert "Here too it is the largest of the three" in largest

    smaller = cheverud_note(_dependence_record(li_ji=3.0, galwey=3.1, cheverud_nyholt=2.0))
    assert "Li & Ji (2005) found overly large" in smaller
    assert "largest" not in smaller
    assert "disagrees" not in smaller
