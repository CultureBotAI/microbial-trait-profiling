"""The integrated evidence table: states, families and what a join must not lose (#178)."""

from __future__ import annotations

import csv
import json
from dataclasses import fields, replace
from pathlib import Path
from typing import TypedDict

import pytest

from mtp.guild_maps import build_matrix, cell_status_rows, compute_correlation
from mtp.integrated_evidence import (
    ALL_TIED,
    ASSESSMENTS,
    CELL_STATES,
    DRAWN,
    GUILD_MIN_STUDIES,
    INSUFFICIENT,
    MASKED,
    NO_ELIGIBLE_PAIRS,
    NOT_APPLICABLE,
    NOT_DRAWN,
    NOT_SAMPLED,
    NOT_SIGNIFICANT,
    SIGNIFICANT,
    SOURCE_MISSING,
    TESTED,
    UNDECIDED,
    UNTESTED,
    ZERO_ANNOTATION,
    AssessmentStatus,
    Record,
    apply_statuses,
    assessment_rows,
    build_guild_matrix,
    cell_support_records,
    family_inventory,
    feature_column_rows,
    fraction_sign_records,
    guild_records,
    is_cell_row,
    log_ratio_records,
    main,
    placeholder_records,
    read_status,
    summarize,
    unmatched_terms,
    write_table,
    write_view,
)
from mtp.paths import file_sha256
from mtp.table1_trait_maps import compute_direction_shifts
from mtp.trait_composition import compute_log_ratio_shifts
from mtp.trait_guilds import (
    OXYGEN_GUILDS,
    TABLE1_TRAITS,
    Profile,
    drop_zero_total_pairs,
    resolve_trait_columns,
)


class _Thresholds(TypedDict):
    q_threshold: float
    strict_q_threshold: float


#: Both decision thresholds every record builder now takes (#294).
THRESHOLDS: _Thresholds = {"q_threshold": 0.05, "strict_q_threshold": 0.001}

TRAITS = tuple(TABLE1_TRAITS)
#: A per-assessment status as `read_status` now returns it. #182 has a verdict
#: for one procedure and excludes the rest; the others are pending.
STATUSES: dict[str, AssessmentStatus] = {
    "calibration_status": AssessmentStatus(
        issue="#182",
        verdicts={"fraction_sign_test": "pass"},
        default="not_calibrated_here #182",
    ),
    "sensitivity_status": AssessmentStatus(
        issue="#183",
        verdicts={"log_ratio_sign_test": "scenarios reported, no verdict #183"},
        default="pending #183",
    ),
    "coverage_status": AssessmentStatus(issue="#184", verdicts={}, default="pending #184"),
    "orientation_status": AssessmentStatus(
        issue="#188", verdicts={}, default="contrast population described, no verdict #188"
    ),
    "taxon_influence_status": AssessmentStatus(issue="#190", verdicts={}, default="pending #190"),
}

#: The matrix header, as the pinned corpus actually spells it: four of the six
#: Table 1 terms are bare CURIEs and two are METPO-qualified composites, so the
#: declared term is not the column name and cannot be joined on (#294).
MATRIX_COLUMNS: dict[str, str] = {
    "METPO:1000608": "METPO:1000608",
    "GO:0030245": "METPO:2000103|GO:0030245",
    "METPO:1000668": "METPO:1000668",
    "METPO:1000683": "METPO:1000683",
    "pathways:arsenate_reduction": "METPO:2000103|pathways:arsenate_reduction",
    "METPO:1000623": "METPO:1000623",
    "METPO:1000602": "METPO:1000602",
    "METPO:1000604": "METPO:1000604",
    "METPO:1000603": "METPO:1000603",
    "METPO:1000607": "METPO:1000607",
}


def write_matrix(tmp_path: Path, rows: list[list[str]]) -> Path:
    """Write a profile matrix whose header spells the columns as the corpus does."""
    header = ["study_id", "experiment_id", "inc_dec", *MATRIX_COLUMNS.values()]
    # A MediaDive growth-medium column, which the loader must exclude from both
    # the guild sums and the denominator.
    header.append("mediadive.medium:372")
    path = tmp_path / "study_trait_profiles_wide.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def resolved_columns(tmp_path: Path) -> dict[str, tuple[tuple[str, str], ...]]:
    """The real resolution of the real header, not a hand-written expectation."""
    return resolve_trait_columns(write_matrix(tmp_path, []), {**TABLE1_TRAITS, **OXYGEN_GUILDS})


#: Per-direction trait counts for one profile of the synthetic corpus. The
#: MediaDive column carries a count too, so the run has to exclude it from both
#: the guild sums and the denominator rather than merely tolerate it.
def _counts(direction: str, kind: str) -> dict[str, int]:
    high = direction == "inc"
    if kind == "zero_total":
        return dict.fromkeys(MATRIX_COLUMNS.values(), 0)
    shapes = {
        "METPO:2000103|GO:0030245": 9 if high else 4,
        "METPO:1000668": 9 if high else 4,
        "METPO:1000683": 9 if high else 4,
        "METPO:2000103|pathways:arsenate_reduction": 4 if high else 9,
        "METPO:1000623": 9 if high else 4,
    }
    if kind == "no_oxygen":
        # Annotated, and carrying no term of either oxygen guild: the cell is
        # sampled and its guild signal is genuinely absent.
        oxygen = dict.fromkeys(
            (
                "METPO:1000608",
                "METPO:1000602",
                "METPO:1000604",
                "METPO:1000603",
                "METPO:1000607",
            ),
            0,
        )
    else:
        oxygen = {
            "METPO:1000608": 30 if high else 5,
            "METPO:1000602": 20 if high else 5,
            "METPO:1000604": 20 if high else 5,
            "METPO:1000603": 2 if high else 30,
            "METPO:1000607": 2 if high else 30,
        }
    return {**shapes, **oxygen}


#: (studies, condition, body site, host, kind, directions) of the corpus below.
CORPUS_BLOCKS: tuple[tuple[int, str, str, str, str, tuple[str, ...]], ...] = (
    (6, "Colitis", "Feces", "Homo sapiens", "signal", ("inc", "dec")),
    (2, "Colitis", "Saliva", "Homo sapiens", "signal", ("inc", "dec")),
    (4, "IBS", "Saliva", "Homo sapiens", "signal", ("inc", "dec")),
    (3, "Crohn", "Feces", "Homo sapiens", "no_oxygen", ("inc", "dec")),
    (1, "Colitis", "Feces", "Homo sapiens", "zero_total", ("inc", "dec")),
    (1, "Colitis", "Feces", "Mus musculus", "signal", ("inc", "dec")),
    (1, "Colitis", "Feces", "Homo sapiens", "signal", ("inc",)),
)


def _corpus(tmp_path: Path) -> tuple[Path, Path]:
    """A matrix and an experiments table exercising every filter the run applies.

    Twelve human paired studies with an oxygen signal, three more that carry no
    oxygen term at all, one whose retained-annotation total is zero on both
    sides, one mouse study, and one experiment with only an increased side.
    """
    rows: list[list[str]] = []
    experiments: list[list[str]] = []
    index = 0
    for count, condition, site, host, kind, directions in CORPUS_BLOCKS:
        for _ in range(count):
            index += 1
            study = f"bsdb:{index}"
            experiment = f"{study}-1"
            experiments.append([study, experiment, condition, site, host])
            for direction in directions:
                counts = _counts(direction, kind)
                rows.append(
                    [
                        study,
                        experiment,
                        direction,
                        *(str(counts[column]) for column in MATRIX_COLUMNS.values()),
                        "77",
                    ]
                )

    profile_path = write_matrix(tmp_path, rows)
    experiments_path = tmp_path / "experiments.tsv"
    with experiments_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["study_id", "experiment_id", "condition", "body_site", "host_species"])
        writer.writerows(experiments)
    return profile_path, experiments_path


def _profiles(pattern: str = "signal") -> list[Profile]:
    """Twelve studies whose counts give a signal, a tie, or no annotation at all."""
    profiles: list[Profile] = []
    for study in range(12):
        for direction in ("inc", "dec"):
            counts = dict.fromkeys(TRAITS, 5)
            if pattern == "signal" and direction == "inc":
                counts[TRAITS[0]] = 30
            if pattern == "tied":
                counts = dict.fromkeys(TRAITS, 5)
            if pattern == "empty":
                counts = dict.fromkeys(TRAITS, 0)
            profiles.append(
                Profile(
                    study_id=f"s{study}",
                    experiment_id=f"s{study}-1",
                    inc_dec=direction,
                    condition="C",
                    body_site="Feces",
                    guild_counts=counts,
                    total=sum(counts.values()) + (0 if pattern == "empty" else 20),
                )
            )
    return profiles


def _ratio_profiles() -> list[Profile]:
    """Twelve studies of two experiments each, with uneven ratio eligibility.

    A log-ratio population is conditional: an experiment counts for a pair only
    when both traits are positive on both sides. Zeroing one trait in some
    experiments makes ``studies_total``, ``studies_eligible``,
    ``paired_experiments`` and ``experiments_eligible`` four different numbers.
    A fixture in which they coincide cannot tell those columns apart, so it would
    pass on code that exported one of them in another's place.
    """
    sparse = TRAITS[1]
    profiles: list[Profile] = []
    for study in range(12):
        for experiment in range(2):
            for direction in ("inc", "dec"):
                counts = dict.fromkeys(TRAITS, 5)
                if direction == "inc":
                    counts[TRAITS[0]] = 30
                # Studies 10 and 11 never carry the sparse trait, so every pair
                # involving it drops them; studies 8 and 9 carry it in their
                # second experiment only, so they are partially observed.
                if study >= 10 or (study >= 8 and experiment == 0):
                    counts[sparse] = 0
                profiles.append(
                    Profile(
                        study_id=f"s{study}",
                        experiment_id=f"s{study}-{experiment}",
                        inc_dec=direction,
                        condition="C",
                        body_site="Feces",
                        guild_counts=counts,
                        total=sum(counts.values()) + 20,
                    )
                )
    return profiles


def _guild_profiles(reverse: int = 0) -> list[Profile]:
    """Twelve studies with discordant oxygen guild signs, ``reverse`` of them flipped.

    The McNemar runs on the two pooled guilds, so the counts have to be keyed on
    ``aerobe`` and ``anaerobe`` rather than on the six Table 1 trait names. With
    ``reverse=0`` no study is discordant the other way and the discordant ratio
    has a zero denominator; with ``reverse>0`` the ratio is defined.
    """
    profiles: list[Profile] = []
    for study in range(12):
        flipped = study < reverse
        for direction in ("inc", "dec"):
            aerobe_high = (direction == "inc") != flipped
            counts = {
                "aerobe": 30 if aerobe_high else 5,
                "anaerobe": 5 if aerobe_high else 30,
            }
            profiles.append(
                Profile(
                    study_id=f"g{study}",
                    experiment_id=f"g{study}-1",
                    inc_dec=direction,
                    condition="C",
                    body_site="Feces",
                    guild_counts=counts,
                    total=sum(counts.values()) + 20,
                )
            )
    return profiles


def _built(
    pattern: str = "signal",
    *,
    q_threshold: float | None = 0.05,
    strict_q_threshold: float | None = 0.001,
) -> list[Record]:
    """Fraction-sign records with no resolved columns, for state-only assertions."""
    return fraction_sign_records(
        _profiles(pattern),
        source="x",
        columns={},
        q_threshold=q_threshold,
        strict_q_threshold=strict_q_threshold,
    )


def test_counts_reconcile_with_the_source_run() -> None:
    profiles = _profiles()
    built = fraction_sign_records(profiles, source="x", columns={}, **THRESHOLDS)
    records = {r.trait: r for r in built}
    for shift in compute_direction_shifts(profiles).shifts:
        record = records[shift.trait]
        assert record.testable_units == shift.studies_scored
        assert record.ties == shift.studies_tied
        assert record.effect_value == shift.share_up
        assert record.q_value == shift.q_value
        assert record.unique_studies == 12
        assert record.family_size == 6


def test_an_all_tied_trait_is_not_a_negative_result() -> None:
    records = fraction_sign_records(_profiles("tied"), source="x", columns={}, **THRESHOLDS)
    for record in records:
        assert record.test_status == ALL_TIED
        assert record.decision == UNDECIDED
        assert record.effect_value is None
        assert record.q_value is None
        assert record.interval_type == "unavailable"


def test_a_pair_with_no_eligible_experiment_is_not_absent_annotation() -> None:
    """An eligibility failure and absent annotation are different claims (#294).

    ``no_positive_pairs`` was mapped onto ``zero_annotation_signal``, which says
    the annotation is not there. The second fixture below shows why that is wrong:
    both traits are annotated on every profile, and the pair is still untestable
    because they are never both positive on both sides of one experiment.
    """
    empty = log_ratio_records(_profiles("empty"), source="x", columns={}, **THRESHOLDS)
    assert empty
    for record in empty:
        assert record.test_status == NO_ELIGIBLE_PAIRS
        assert record.test_status != ZERO_ANNOTATION
        assert record.decision == UNDECIDED
        assert record.q_value is None

    # Two traits that are each well annotated but never co-occur on one side.
    exclusive: list[Profile] = []
    for study in range(6):
        for direction in ("inc", "dec"):
            counts = dict.fromkeys(TRAITS, 4)
            counts[TRAITS[0]] = 9 if direction == "inc" else 0
            counts[TRAITS[1]] = 0 if direction == "inc" else 9
            exclusive.append(
                Profile(
                    study_id=f"e{study}",
                    experiment_id=f"e{study}-1",
                    inc_dec=direction,
                    condition="C",
                    body_site="Feces",
                    guild_counts=counts,
                    total=sum(counts.values()) + 20,
                )
            )
    pair = next(
        record
        for record in log_ratio_records(exclusive, source="x", columns={}, **THRESHOLDS)
        if (record.trait, record.reference) == (TRAITS[0], TRAITS[1])
    )
    assert pair.test_status == NO_ELIGIBLE_PAIRS
    assert pair.experiments_excluded_zero_count == 6
    # The annotation is present on every profile; only the pairing fails.
    assert all(
        profile.guild_counts[TRAITS[0]] + profile.guild_counts[TRAITS[1]] > 0
        for profile in exclusive
    )


def test_an_unsupplied_method_is_not_a_nonsignificant_result() -> None:
    records = placeholder_records("enrichment", "profile", "not supplied", TRAITS)
    assert len(records) == len(TRAITS)
    for record in records:
        assert record.test_status == SOURCE_MISSING
        assert record.decision == UNDECIDED
        assert record.q_value is None
        assert record.notes == "not supplied"


def test_each_method_keeps_its_own_family_and_units() -> None:
    profiles = _profiles()
    records = fraction_sign_records(
        profiles, source="x", columns={}, **THRESHOLDS
    ) + log_ratio_records(profiles, source="x", columns={}, **THRESHOLDS)
    families = {(r.method, r.family_id, r.correction_method) for r in records}
    assert (
        "directional_consistency_fraction_sign_test",
        "table1_six_fractions",
        "Benjamini-Hochberg",
    ) in families
    assert (
        "directional_consistency_log_ratio_sign_test",
        "table1_fifteen_ratios",
        "Benjamini-Yekutieli",
    ) in families
    assert len({r.effect_units for r in records}) == 2
    summary = {row["method"]: row for row in summarize(records)}
    assert summary["directional_consistency_fraction_sign_test"]["families"] == (
        "table1_six_fractions"
    )
    assert all(
        row["combined_guarantee"] == "none; each family is corrected on its own"
        for row in summary.values()
    )


def _support_profiles() -> list[Profile]:
    """A corpus whose figure cells land in each of the four support states.

    - Colitis in Feces: four studies, enough to be drawn.
    - Colitis in Saliva: two studies, below the min-studies rule, so masked.
    - Crohn in Feces: three studies in which no taxon carries either guild, so
      sampled with no signal.
    - IBS in Saliva: four studies, so the Saliva column block survives the rule
      that drops a body site no condition is drawn in - without it the masked
      Colitis/Saliva cell would be dropped with the whole block.
    - Crohn in Saliva and IBS in Feces: never sampled.
    """
    profiles: list[Profile] = []

    def add(study: str, condition: str, site: str, *, signal: bool) -> None:
        for direction in ("inc", "dec"):
            aerobe_high = direction == "inc"
            counts = (
                {"aerobe": 30 if aerobe_high else 5, "anaerobe": 5 if aerobe_high else 30}
                if signal
                else {"aerobe": 0, "anaerobe": 0}
            )
            profiles.append(
                Profile(
                    study_id=study,
                    experiment_id=f"{study}-{site}",
                    inc_dec=direction,
                    condition=condition,
                    body_site=site,
                    guild_counts=counts,
                    total=sum(counts.values()) + 20,
                )
            )

    for index in range(4):
        add(f"drawn{index}", "Colitis", "Feces", signal=True)
    for index in range(2):
        add(f"thin{index}", "Colitis", "Saliva", signal=True)
    for index in range(3):
        add(f"blank{index}", "Crohn", "Feces", signal=False)
    for index in range(4):
        add(f"keep{index}", "IBS", "Saliva", signal=True)
    return profiles


def _cells_by_stratum(records: list[Record]) -> dict[tuple[str, str, str, str], Record]:
    return {(r.condition, r.body_site, r.trait, r.direction): r for r in records if is_cell_row(r)}


def test_the_four_cell_support_states_come_from_the_real_matrix() -> None:
    """Drawn, masked, sampled-without-signal and not-sampled, each emitted (#294).

    Nothing used to set ``drawn`` or ``masked_below_min_studies``, and
    ``not_sampled`` was never emitted at all, so three of the enum's states named
    something the table could not produce. These rows come from the same
    ``cell_status`` classification the figure writer uses, over the matrix the
    McNemar row is computed on.
    """
    profiles = _support_profiles()
    matrix = build_guild_matrix(profiles)
    records = cell_support_records(matrix, source="x", columns={}, min_studies=GUILD_MIN_STUDIES)
    cells = _cells_by_stratum(records)

    drawn = cells[("Colitis", "Feces", "aerobe", "inc")]
    assert (drawn.display_status, drawn.test_status) == (DRAWN, UNTESTED)
    assert drawn.unique_studies == 4

    masked = cells[("Colitis", "Saliva", "aerobe", "inc")]
    assert (masked.display_status, masked.test_status) == (MASKED, INSUFFICIENT)
    assert masked.unique_studies == 2
    assert masked.support_threshold == GUILD_MIN_STUDIES == 3

    no_signal = cells[("Crohn", "Feces", "aerobe", "inc")]
    assert (no_signal.display_status, no_signal.test_status) == (DRAWN, ZERO_ANNOTATION)
    assert no_signal.unique_studies == 3, "sampled; the annotation is what is absent"

    not_sampled = cells[("Crohn", "Saliva", "aerobe", "inc")]
    assert (not_sampled.display_status, not_sampled.test_status) == (NOT_DRAWN, NOT_SAMPLED)
    assert not_sampled.unique_studies == 0

    # The four states are different claims, so no two of them collapse.
    assert len({(c.display_status, c.test_status) for c in cells.values()}) == 4
    # Every cell of the matrix is present, including the unsampled ones.
    assert len([r for r in records if is_cell_row(r)]) == len(matrix.conditions) * len(
        matrix.columns
    )


def test_cell_rows_line_up_with_the_figure_writers_cell_grid() -> None:
    """The rows must be the figure's own cells, in the figure's own order.

    This pins the alignment, not the mapping: it re-derives the expectation from
    ``CELL_STATES``, so a change to that table passes here and is caught by
    ``test_the_four_cell_support_states_come_from_the_real_matrix`` instead. What
    it does catch is a row emitted against the wrong cell - a transposed index,
    a column parsed into the wrong site or direction, or a skipped cell.
    """
    matrix = build_guild_matrix(_support_profiles())
    records = cell_support_records(matrix, source="x", columns={}, min_studies=GUILD_MIN_STUDIES)
    statuses = cell_status_rows(matrix)

    expected = [
        CELL_STATES[statuses[row][column]]
        for row in range(len(matrix.conditions))
        for column in range(len(matrix.columns))
    ]
    assert [(r.display_status, r.test_status) for r in records] == [
        (display, test) for test, display in expected
    ]
    # Each row names the cell it came from, so the grid can be rebuilt from it.
    assert [f"{r.body_site}|{r.trait}_{r.direction}" for r in records] == [
        column for _ in matrix.conditions for column in matrix.columns
    ]
    assert [r.condition for r in records] == [
        condition for condition in matrix.conditions for _ in matrix.columns
    ]


def test_a_cell_row_parses_a_trait_name_that_contains_an_underscore() -> None:
    """The column is ``site|feature_direction``, and five of the features have
    underscores in their own names.

    The oxygen guilds do not, so a fixture built on them cannot tell a split on
    the first underscore from a split on the last. The six Table 1 traits can,
    and this builder is given whichever matrix it is handed.
    """
    traits = ("facultatively_aerobic", "coccus_shaped")
    profiles = [
        Profile(
            study_id=f"u{study}",
            experiment_id=f"u{study}-1",
            inc_dec=direction,
            condition="Colitis",
            body_site="Feces",
            guild_counts={
                "facultatively_aerobic": 30 if direction == "inc" else 5,
                "coccus_shaped": 5,
            },
            total=60,
        )
        for study in range(4)
        for direction in ("inc", "dec")
    ]
    matrix = build_matrix(
        profiles,
        max_conditions=10,
        max_body_sites=4,
        min_studies=3,
        guild_order=traits,
    )
    records = cell_support_records(matrix, source="x", columns={}, min_studies=3)

    cells = _cells_by_stratum(records)
    assert ("Colitis", "Feces", "facultatively_aerobic", "inc") in cells
    assert ("Colitis", "Feces", "coccus_shaped", "dec") in cells
    assert {r.trait for r in records} == set(traits)
    assert {r.direction for r in records} == {"inc", "dec"}


def test_a_masked_cell_does_not_change_a_tested_result() -> None:
    """Display state never feeds back into a decision (#178).

    The previous version of this test built its own masked record with
    ``dataclasses.replace``, so it passed whatever the code did. This one takes
    a corpus in which the tested corpus-level result and a masked figure cell
    are both produced by the real builders, and checks that the masking is
    recorded without touching the result.
    """
    profiles = _support_profiles()
    matrix = build_guild_matrix(profiles)
    (tested,) = guild_records(profiles, matrix=matrix, source="x", columns={}, **THRESHOLDS)
    cells = _cells_by_stratum(
        cell_support_records(matrix, source="x", columns={}, min_studies=GUILD_MIN_STUDIES)
    )

    masked = cells[("Colitis", "Saliva", "aerobe", "inc")]
    assert masked.display_status == MASKED
    assert masked.decision == UNDECIDED
    assert masked.q_value is None

    # The corpus-level test ran over the same studies and is unaffected by what
    # the figure draws, so it is not a figure cell and has no display state.
    assert tested.test_status == TESTED
    assert tested.display_status == NOT_APPLICABLE
    assert not is_cell_row(tested)
    assert tested.p_value is not None


def test_statuses_are_attached_without_changing_any_result() -> None:
    records = _built()
    annotated = apply_statuses(records, STATUSES)
    for before, after in zip(records, annotated, strict=True):
        assert (after.decision, after.test_status, after.q_value) == (
            before.decision,
            before.test_status,
            before.q_value,
        )
        assert after.calibration_status == "pass"
        assert after.coverage_status == "pending #184"
    with pytest.raises(ValueError, match="no status supplied"):
        apply_statuses(records, {"calibration_status": STATUSES["calibration_status"]})


def test_each_row_carries_its_own_procedures_verdict_only() -> None:
    """A verdict belongs to the procedure it was reached for (#294).

    #182 reports per procedure and names the arms it does not cover. Copying a
    concatenation of every verdict onto every row attached a calibration result
    to the three arms that assessment had explicitly excluded.
    """
    profiles = _profiles()
    records = apply_statuses(
        fraction_sign_records(profiles, source="x", columns={}, **THRESHOLDS)
        + log_ratio_records(profiles, source="x", columns={}, **THRESHOLDS)
        + placeholder_records(
            "enrichment_arm_profile_hypergeometric", "profile", "not supplied", TRAITS
        ),
        STATUSES,
    )
    by_method: dict[str, Record] = {record.method: record for record in records}

    fraction = by_method["directional_consistency_fraction_sign_test"]
    ratio = by_method["directional_consistency_log_ratio_sign_test"]
    enrichment = by_method["enrichment_arm_profile_hypergeometric"]

    # Only the procedure #182 reached carries its verdict.
    assert fraction.calibration_status == "pass"
    assert ratio.calibration_status == "not_calibrated_here #182"
    assert enrichment.calibration_status == "not_calibrated_here #182"
    assert "pass" not in enrichment.calibration_status

    # #183 restricts the log-ratio population and assesses nothing else.
    assert ratio.sensitivity_status == "scenarios reported, no verdict #183"
    assert fraction.sensitivity_status == "pending #183"
    assert enrichment.sensitivity_status == "pending #183"

    # #188 describes the corpus, so it reaches every row the same way.
    assert {record.orientation_status for record in records} == {
        "contrast population described, no verdict #188"
    }


def test_a_sensitivity_scenario_count_is_not_a_verdict(tmp_path: Path) -> None:
    """#183 reports per scenario; a count of scenarios read as cleared (#294)."""
    path = tmp_path / "sensitivity_status.json"
    path.write_text(
        json.dumps(
            {
                "scenario_status": {
                    "zero_exclusion": {"tested": 9, "all_tied": 6},
                    "reference_choice": {"tested": 15},
                }
            }
        ),
        encoding="utf-8",
    )

    status = read_status(path, "sensitivity_status", "#183")

    applied = status.for_procedure("log_ratio_sign_test")
    assert "no verdict" in applied
    assert applied != "2 scenarios reported"
    assert not any(character.isdigit() for character in applied.replace("#183", ""))
    # A procedure the scenarios say nothing about stays pending.
    assert status.for_procedure("fraction_sign_test") == "pending #183"
    # The per-scenario numbers are kept, attached to the scenario they describe.
    assert ("scenario", "zero_exclusion", "tested", "9") in status.detail
    assert ("scenario", "reference_choice", "tested", "15") in status.detail


def test_a_calibration_file_keeps_its_uncovered_arms(tmp_path: Path) -> None:
    """`not_calibrated_here` is preserved rather than being read as a pass."""
    path = tmp_path / "calibration_status.json"
    path.write_text(
        json.dumps(
            {
                "calibration_status": {
                    "fraction_sign_test": {
                        "overall": "pass",
                        "criteria": {"per_test_level": "pass", "family_error": "pass"},
                    },
                    "guild_mcnemar": {
                        "overall": "fail",
                        "criteria": {"per_test_level": "pass", "family_error": "fail"},
                    },
                },
                "not_calibrated_here": ["the enrichment arm's hypergeometric urn (#131)"],
            }
        ),
        encoding="utf-8",
    )

    status = read_status(path, "calibration_status", "#182")

    assert status.for_procedure("fraction_sign_test") == (
        "pass #182 [per_test_level=pass, family_error=pass]"
    )
    # A failed family criterion is no longer dropped behind a per-test pass.
    assert status.for_procedure("guild_mcnemar") == (
        "fail #182 [per_test_level=pass, family_error=fail]"
    )
    assert ("procedure", "guild_mcnemar", "criterion:family_error", "fail") in status.detail
    assert status.for_procedure("log_ratio_sign_test") == "not_calibrated_here #182"
    assert (
        "note",
        "not_calibrated_here",
        "procedure",
        "the enrichment arm's hypergeometric urn (#131)",
    ) in status.detail


def test_absent_assessment_files_read_as_pending(tmp_path: Path) -> None:
    assert read_status(None, "calibration_status", "#182").default == "pending #182"
    missing = read_status(tmp_path / "missing.json", "coverage_status", "#184")
    assert missing.for_procedure("anything") == "pending #184"
    broken = tmp_path / "broken.json"
    broken.write_text("not json", encoding="utf-8")
    assert read_status(broken, "coverage_status", "#184").default == "unreadable #184"
    written = tmp_path / "calibration.json"
    written.write_text(
        '{"calibration_status": {"a": {"overall": "pass", "criteria": {"c": "pass"}}}}',
        encoding="utf-8",
    )
    assert read_status(written, "calibration_status", "#182").for_procedure("a") == (
        "pass #182 [c=pass]"
    )


@pytest.mark.parametrize(
    ("payload", "key", "why"),
    [
        ("[1, 2, 3]", "coverage_status", "JSON list"),
        ('"just a string"', "coverage_status", "JSON str"),
        ("17", "coverage_status", "JSON int"),
        ('{"coverage_status": {"nested": "object"}}', "coverage_status", "not a string"),
        ('{"coverage_status": [1]}', "coverage_status", "not a string"),
        ('{"calibration_status": ["a", "b"]}', "calibration_status", "not a mapping"),
        ('{"calibration_status": {"a": "pass"}}', "calibration_status", "single verdict"),
        ('{"calibration_status": {"a": {"b": "c"}}}', "calibration_status", "lacks an overall"),
        ('{"family_status": ["one"]}', "sensitivity_status", "not an object"),
        ('{"family_status": {"S4:x": {"tested": 6}}}', "sensitivity_status", "lacks one of"),
        ('{"scenario_status": ["one"]}', "sensitivity_status", "not an object"),
        ('{"class_counts": {"a": 1}}', "orientation_status", "not a list"),
        ('{"class_counts": ["a"]}', "orientation_status", "not an object"),
    ],
)
def test_a_status_file_of_the_wrong_shape_is_unreadable_and_says_why(
    tmp_path: Path, payload: str, key: str, why: str
) -> None:
    """A status file is validated, not trusted (#294).

    The previous version called `.get` on whatever JSON held, so an array raised
    `AttributeError` and killed the run, and a nested object was written into a
    status column verbatim.
    """
    path = tmp_path / "status.json"
    path.write_text(payload, encoding="utf-8")

    status = read_status(path, key, "#184")

    applied = status.for_procedure("fraction_sign_test")
    assert applied.startswith("unreadable #184")
    assert why in applied
    assert status.verdicts == {}


def test_every_field_survives_the_export(tmp_path: Path) -> None:
    records = apply_statuses(
        _built() + placeholder_records("enrichment", "profile", "not supplied", TRAITS[:1]),
        STATUSES,
    )
    path = write_table(records, tmp_path / "evidence.tsv")
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == len(records)
    assert {"test_status", "decision", "display_status", "family_id", "interval_type"} <= set(
        rows[0]
    )
    unsupplied = next(r for r in rows if r["method"] == "enrichment")
    assert unsupplied["q_value"] == ""
    assert unsupplied["test_status"] == SOURCE_MISSING
    assert unsupplied["decision"] == UNDECIDED
    signal = next(r for r in rows if r["trait"] == TRAITS[0] and r["method"] != "enrichment")
    assert signal["decision"] in (SIGNIFICANT, NOT_SIGNIFICANT)
    view = write_view(records, summarize(records), tmp_path / "evidence.md")
    text = view.read_text(encoding="utf-8")
    assert "no pooled p-value" in text
    assert "source_not_supplied" in text


def test_a_record_carries_the_selection_note() -> None:
    record = _built()[0]
    assert "selected on this corpus" in record.selection_note
    assert isinstance(record, Record)


def test_a_zero_total_pair_is_excluded_rather_than_counted_as_a_tie() -> None:
    """A pair with no retained annotation on a side has an undefined share (#294).

    ``Profile.fraction`` returns 0.0 for a zero total, so such a pair enters the
    contrast as a measured share of nothing, the contrast comes out exactly zero,
    and the sign test counts it as a tie. That merges all-tied studies with
    zero-annotation studies, which the shared comparison basis of #178 forbids.
    """
    scored = _profiles()
    blank = [
        Profile(
            study_id="bsdb:blank",
            experiment_id="bsdb:blank-1",
            inc_dec=direction,
            condition="c",
            body_site="s",
            guild_counts=dict.fromkeys(TRAITS, 0),
            total=0,
        )
        for direction in ("inc", "dec")
    ]

    exclusion = drop_zero_total_pairs(scored + blank)

    assert exclusion.zero_total_experiments == 1
    assert exclusion.studies_lost == 1
    assert exclusion.studies_after == len({p.study_id for p in scored})
    assert all(p.study_id != "bsdb:blank" for p in exclusion.profiles)

    # Kept as a tie under historical scoring, the blank study would inflate the tie count;
    # excluded under the shared comparison basis, the tallies match the corpus that has
    # a defined share (#294, #228).
    with_historical = fraction_sign_records(
        scored + blank, source="x", columns={}, exclude_zero_totals=False, **THRESHOLDS
    )
    without = fraction_sign_records(exclusion.profiles, source="x", columns={}, **THRESHOLDS)
    assert {r.trait: r.ties for r in with_historical} != {r.trait: r.ties for r in without}
    with_blank = fraction_sign_records(scored + blank, source="x", columns={}, **THRESHOLDS)
    assert {r.trait: r.ties for r in with_blank} == {r.trait: r.ties for r in without}


def test_each_row_is_decided_at_both_thresholds_with_their_sources() -> None:
    """One cutoff cannot speak for three families, so both are reported (#294)."""
    records = _built()

    for record in records:
        assert record.q_threshold == 0.05
        assert record.strict_q_threshold == 0.001
        # Each threshold says where it comes from rather than appearing as a bare number.
        assert record.q_threshold_source
        assert record.strict_q_threshold_source
        assert record.q_threshold_source != record.strict_q_threshold_source


def test_the_two_thresholds_decide_independently() -> None:
    """The strict column must read its own threshold, not copy the loose decision.

    The fixture's q values sit below both default thresholds, so a run at the
    defaults cannot tell the two columns apart - which is how a strict column
    that silently reuses ``q_threshold`` would pass. The thresholds here bracket
    the fixture's q instead, so the two decisions are required to disagree (#294).
    """
    records = _built(q_threshold=0.05, strict_q_threshold=1e-9)

    assert records
    for record in records:
        assert record.q_value is not None
        assert 1e-9 < record.q_value <= 0.05, "fixture must sit between the two thresholds"
        assert record.decision == SIGNIFICANT
        assert record.strict_decision == NOT_SIGNIFICANT


def test_the_stricter_threshold_can_only_shrink_the_significant_set() -> None:
    """A row significant at the stricter cutoff is significant at the looser one."""
    records = _built()

    assert records
    for record in records:
        if record.strict_decision == SIGNIFICANT:
            assert record.decision == SIGNIFICANT


def test_a_row_exports_the_raw_matrix_column_not_the_declared_term(tmp_path: Path) -> None:
    """The declared CURIE does not name a column, so it cannot be joined on (#294).

    Two of the six Table 1 traits are carried as METPO-qualified composites. A
    row exported under the bare term, and a ratio row that drops the reference's
    identifier altogether, cannot be joined back to the matrix at all.
    """
    columns = resolved_columns(tmp_path)
    profiles = _profiles()

    fractions = {
        r.trait: r
        for r in fraction_sign_records(profiles, source="x", columns=columns, **THRESHOLDS)
    }
    cellulose = fractions["cellulose_catabolic_process"]
    assert cellulose.trait_qualifier == "GO:0030245"
    assert cellulose.trait_columns == "METPO:2000103|GO:0030245"
    assert cellulose.trait_columns != cellulose.trait_qualifier
    arsenate = fractions["arsenate_reduction"]
    assert arsenate.trait_columns == "METPO:2000103|pathways:arsenate_reduction"

    ratios = log_ratio_records(profiles, source="x", columns=columns, **THRESHOLDS)
    pair = next(
        r
        for r in ratios
        if (r.trait, r.reference) == ("cellulose_catabolic_process", "arsenate_reduction")
    )
    # Both halves of the hypothesis are named, each with its own raw column.
    assert pair.reference_qualifier == "pathways:arsenate_reduction"
    assert pair.reference_columns == "METPO:2000103|pathways:arsenate_reduction"
    assert pair.trait_columns == "METPO:2000103|GO:0030245"


def test_the_feature_column_sidecar_is_atomic_and_covers_every_feature(tmp_path: Path) -> None:
    """A pooled guild lists several columns in one cell; the side-car splits them.

    Each of the six Table 1 traits resolves to exactly one column, so a side-car
    that simply copied the main table's packed cell would be indistinguishable on
    them. The aerobe guild pools three terms, which is what forces the split.
    """
    columns = resolved_columns(tmp_path)
    records = fraction_sign_records(
        _profiles(), source="x", columns=columns, **THRESHOLDS
    ) + placeholder_records(
        "oxygen_guild_placeholder", "corpus", "not supplied", ("aerobe",), columns=columns
    )
    rows = feature_column_rows(records, columns)

    assert {row["feature"] for row in rows} == {*TRAITS, "aerobe"}
    for row in rows:
        assert ";" not in row["matrix_column"]
        assert ";" not in row["term_id"]
    # One row per pooled term, not one packed row for the guild.
    aerobe = [row for row in rows if row["feature"] == "aerobe"]
    assert {row["matrix_column"] for row in aerobe} == set(OXYGEN_GUILDS["aerobe"])
    # Header order, not declaration order: the header puts METPO:1000608 first.
    assert [row["matrix_column"] for row in aerobe] == [
        "METPO:1000608",
        "METPO:1000602",
        "METPO:1000604",
    ]
    assert {
        "feature": "cellulose_catabolic_process",
        "term_id": "GO:0030245",
        "matrix_column": "METPO:2000103|GO:0030245",
    } in rows


def test_an_unresolved_declared_term_is_reported(tmp_path: Path) -> None:
    """A guild with no column at all is rejected upstream; one missing term is not.

    The aerobe guild pools three terms. If one of them matched no column the sum
    would silently cover two, so the run names the term rather than staying quiet.
    """
    columns = resolved_columns(tmp_path)
    assert unmatched_terms(columns, {**TABLE1_TRAITS, **OXYGEN_GUILDS}) == []
    assert unmatched_terms(columns, {"aerobe": ("METPO:1000602", "METPO:9999999")}) == [
        "aerobe/METPO:9999999"
    ]


def test_a_ratio_row_exports_its_own_population_and_effect_summary() -> None:
    """The ratio population is conditional, so each field is its own column (#294).

    The mean previously survived only as a string inside ``notes``, and the
    population fields the #178 spec lists were absent, so a reader could not tell
    how many studies a ratio hypothesis had dropped.
    """
    profiles = _ratio_profiles()
    records = {
        (r.trait, r.reference): r
        for r in log_ratio_records(profiles, source="x", columns={}, **THRESHOLDS)
    }
    summaries = compute_log_ratio_shifts(profiles).summaries
    assert summaries

    # The fixture must separate the four population counts, or an export that
    # put one of them in another's column would pass this test.
    assert any(r.studies_total != r.unique_studies for r in records.values())
    assert any(r.paired_experiments != r.experiments_eligible for r in records.values())
    assert any((r.studies_partially_observed or 0) > 0 for r in records.values())
    assert any((r.studies_excluded or 0) > 0 for r in records.values())

    for summary in summaries:
        record = records[(summary.numerator, summary.reference)]
        assert record.studies_total == summary.studies_total
        assert record.studies_excluded == summary.studies_excluded
        assert record.studies_partially_observed == summary.studies_partially_observed
        assert record.paired_experiments == summary.paired_experiments
        assert record.experiments_eligible == summary.eligible_experiments
        assert record.experiments_excluded_zero_count == summary.excluded_zero_experiments
        assert record.study_effect_mean == summary.mean_log2_ratio_change
        assert record.study_effect_median == summary.median_log2_ratio_change
        assert record.study_effect_units
        # The numbers are in numeric columns, so no text column carries one.
        assert record.notes == ""
        assert not any(character.isdigit() for character in record.zero_annotation_state)


def test_every_family_exports_its_size_and_the_size_is_checked() -> None:
    """A q is unreadable without its family and that family's size (#236, #289)."""
    profiles = _profiles()
    records = (
        fraction_sign_records(profiles, source="x", columns={}, **THRESHOLDS)
        + log_ratio_records(profiles, source="x", columns={}, **THRESHOLDS)
        + placeholder_records("enrichment", "profile", "not supplied", TRAITS)
    )

    inventory = {row["family_id"]: row for row in family_inventory(records)}
    assert inventory["table1_six_fractions"]["family_size"] == 6
    assert inventory["table1_six_fractions"]["rows_exported"] == 6
    assert inventory["table1_fifteen_ratios"]["family_size"] == 15
    assert inventory["table1_fifteen_ratios"]["rows_exported"] == 15
    # Unsupplied rows carry no q, so they belong to no correction family.
    assert "" not in inventory

    # A family whose declared size no longer matches the rows exported for it
    # means the q was corrected over a different set of tests than the table
    # shows, which makes every q in it unreadable.
    shrunk = [r for r in records if r.family_id != "table1_six_fractions"] + [
        r for r in records if r.family_id == "table1_six_fractions"
    ][:5]
    with pytest.raises(ValueError, match="declares 6 members but the table exports 5"):
        family_inventory(shrunk)


def test_the_mcnemar_row_keeps_its_p_and_leaves_q_blank() -> None:
    """A single uncorrected test has no q, so it gets no q-threshold decision.

    Writing the uncorrected p into ``q_value`` made the guild row a one-member
    FDR family, and both threshold columns then read as though a multiplicity
    correction had been applied to it (#294).
    """
    profiles = _guild_profiles()
    matrix = build_guild_matrix(profiles)
    (record,) = guild_records(profiles, matrix=matrix, source="x", columns={}, **THRESHOLDS)
    correlation = compute_correlation(profiles, matrix)
    up_down, down_up = correlation.discordant

    assert record.test_status == TESTED
    assert record.p_value == correlation.mcnemar_p
    assert record.p_value is not None
    assert record.q_value is None
    # The p clears both thresholds, so a decision on this row could only come
    # from reading the uncorrected p as though it were a q.
    assert record.p_value <= 0.001
    assert record.decision == UNDECIDED
    assert record.strict_decision == UNDECIDED
    assert record.testable_units == up_down + down_up == 12
    assert record.unique_studies == correlation.studies == 12
    assert "no q is computed" in record.correction_method


def test_an_undefined_odds_ratio_is_missing_rather_than_infinite() -> None:
    """Every study discordant the same way leaves the ratio undefined (#294).

    The discordant ratio has a zero denominator when no study is discordant in
    the reverse direction. ``float("inf")`` in a TSV cell is a sentinel value,
    which this repository does not use for missing data. The defined case is
    asserted alongside it, so code that simply never reported a ratio would
    fail here too.
    """
    undefined = _guild_profiles(reverse=0)
    matrix = build_guild_matrix(undefined)
    (blank,) = guild_records(undefined, matrix=matrix, source="x", columns={}, **THRESHOLDS)
    correlation = compute_correlation(undefined, matrix)

    assert correlation.discordant == (12, 0)
    assert correlation.odds_ratio is None
    assert blank.effect_value is None
    assert blank.effect_units == "odds ratio of discordant study signs"
    assert blank.test_status == TESTED, "the test ran; only the ratio is undefined"

    defined = _guild_profiles(reverse=3)
    defined_matrix = build_guild_matrix(defined)
    (ratio,) = guild_records(defined, matrix=defined_matrix, source="x", columns={}, **THRESHOLDS)
    defined_correlation = compute_correlation(defined, defined_matrix)

    assert defined_correlation.discordant == (9, 3)
    assert defined_correlation.odds_ratio == 3.0
    assert ratio.effect_value == 3.0


def test_a_q_value_without_a_family_is_rejected() -> None:
    """Every corrected value must name the family it was corrected in (#236)."""
    record = _built()[0]
    assert record.q_value is not None
    with pytest.raises(ValueError, match="no correction family"):
        family_inventory([replace(record, family_id="", family_size=None)])


def test_main_writes_every_artifact_and_its_provenance(tmp_path: Path) -> None:
    """The CLI end to end: the files it writes, and what they agree on (#294).

    Nothing covered `main()` or the provenance record, so the exclusion counts,
    the thresholds and the family sizes could drift from the table beside them
    without any test noticing.
    """
    profile_path, experiments_path = _corpus(tmp_path)
    output = tmp_path / "out"

    exit_code = main(
        [
            "--profile-input",
            str(profile_path),
            "--experiments-input",
            str(experiments_path),
            "--full-dump-input",
            str(tmp_path / "no_such_full_dump.csv"),
            "--calibration-status",
            str(tmp_path / "no_such_calibration_status.json"),
            "--sensitivity-status",
            str(tmp_path / "no_such_sensitivity_status.json"),
            "--orientation-provenance",
            str(tmp_path / "no_such_orientation_provenance.json"),
            "--output-dir",
            str(output),
            "--code-commit",
            "0123456789abcdef",
        ]
    )

    assert exit_code == 0
    written = {
        "integrated_evidence.tsv",
        "integrated_evidence_method_results.tsv",
        "integrated_evidence_features.tsv",
        "integrated_evidence_families.tsv",
        "integrated_evidence_assessments.tsv",
        "integrated_evidence.md",
        "integrated_evidence_provenance.json",
    }
    assert {path.name for path in output.iterdir()} == written

    with (output / "integrated_evidence.tsv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    provenance = json.loads(
        (output / "integrated_evidence_provenance.json").read_text(encoding="utf-8")
    )

    # Every declared field is exported, under its own name.
    assert list(rows[0]) == [field.name for field in fields(Record)]

    # The corpus cascade: one mouse study and one unpaired increased side are
    # dropped by the loader, and the zero-total pair by the shared basis.
    corpus = provenance["corpus"]
    assert corpus["dropped_non_human"] == 2
    assert corpus["dropped_unpaired"] == 1
    assert corpus["mediadive_columns_excluded"] == 1
    assert corpus["trait_columns_used"] == 10
    exclusion = provenance["pair_exclusion"]
    assert exclusion["zero_total_experiments"] == 1
    assert exclusion["studies_before"] - exclusion["studies_after"] == 1
    assert exclusion["studies_after"] == 15
    assert exclusion["experiments_kept"] == 15
    assert "undefined, not zero" in exclusion["rule"]

    # The exclusion is applied to the data, not merely reported beside it: the
    # tests run on the 15 remaining studies, and the only ties on the oxygen
    # trait are the three studies that carry no oxygen term.
    fractions = [
        row for row in rows if row["method"] == "directional_consistency_fraction_sign_test"
    ]
    assert {row["unique_studies"] for row in fractions} == {str(exclusion["studies_after"])}
    assert {row["studies_total"] for row in fractions} == {str(exclusion["studies_after"])}
    aerobic = next(row for row in fractions if row["trait"] == "facultatively_aerobic")
    assert aerobic["ties"] == "3"
    assert aerobic["testable_units"] == "12"

    # Both thresholds are recorded with their sources, not as bare numbers.
    assert provenance["thresholds"]["q_threshold"] == 0.05
    assert provenance["thresholds"]["strict_q_threshold"] == 0.001
    assert provenance["thresholds"]["q_threshold_source"]
    assert (
        provenance["thresholds"]["q_threshold_source"]
        != provenance["thresholds"]["strict_q_threshold_source"]
    )
    assert provenance["code_commit"] == "0123456789abcdef"
    assert provenance["bugsigdb_vintage"] is None
    assert provenance["declared_terms_without_a_matrix_column"] == []

    # The family export, the provenance record and the table all agree, and the
    # declared size is the number of rows exported for that family.
    with (output / "integrated_evidence_families.tsv").open(encoding="utf-8") as handle:
        families = {row["family_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert {name: row["family_size"] for name, row in families.items()} == {
        "table1_six_fractions": "6",
        "table1_fifteen_ratios": "15",
        "guild_mcnemar": "1",
    }
    for name, row in families.items():
        assert row["rows_exported"] == row["family_size"]
        assert sum(r["family_id"] == name for r in rows) == int(row["family_size"])
    assert {row["family_id"]: row["family_size"] for row in provenance["families"]} == {
        name: int(row["family_size"]) for name, row in families.items()
    }
    # No exported q belongs to no family.
    assert all(row["family_id"] for row in rows if row["q_value"])

    # The raw column is exported and is joinable against the matrix header.
    with (output / "integrated_evidence_features.tsv").open(encoding="utf-8") as handle:
        features = list(csv.DictReader(handle, delimiter="\t"))
    header = profile_path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert features
    for feature in features:
        assert feature["matrix_column"] in header
    assert {
        "feature": "cellulose_catabolic_process",
        "term_id": "GO:0030245",
        "matrix_column": "METPO:2000103|GO:0030245",
    } in features

    # The uncorrected guild test keeps its p and leaves q blank.
    mcnemar = next(row for row in rows if row["method"] == "oxygen_guild_mcnemar")
    assert mcnemar["p_value"]
    assert mcnemar["q_value"] == ""
    assert mcnemar["decision"] == UNDECIDED

    # The cell inventory covers all four support states, and nothing else does.
    cell_states = {
        (row["display_status"], row["test_status"]) for row in rows if row["support_threshold"]
    }
    assert cell_states == {
        (DRAWN, UNTESTED),
        (DRAWN, ZERO_ANNOTATION),
        (MASKED, INSUFFICIENT),
        (NOT_DRAWN, NOT_SAMPLED),
    }
    assert all(row["condition"] and row["body_site"] for row in rows if row["support_threshold"])

    # With no assessment output present, every status is pending, per procedure.
    with (output / "integrated_evidence_assessments.tsv").open(encoding="utf-8") as handle:
        assessments = list(csv.DictReader(handle, delimiter="\t"))
    assert assessments
    assert {row["state_value"] for row in assessments} == {
        f"pending {issue}" for issue in ASSESSMENTS.values() if issue != "#190"
    } | {"pending #190: influence evidence not supplied"}
    assert {row["calibration_status"] for row in rows} == {"pending #182"}

    view = (output / "integrated_evidence.md").read_text(encoding="utf-8")
    assert "no pooled p-value" in view
    assert "Figure cells, by what the figure does with them" in view
    # The cell inventory is counted in the view, never listed row by row.
    before, records_section = view.split("## Records", 1)
    assert "oxygen_guild_cell_support" not in records_section
    assert "oxygen_guild_cell_support" in before
    assert "oxygen_guild_mcnemar" in records_section


FIXTURES = Path(__file__).parent / "fixtures"


def test_the_sensitivity_producers_family_records_survive_the_export() -> None:
    """Cross-PR: #237's status at c9fb6f1 keeps S4's two families apart (#241 review).

    The producer and consumer tests each passed on their own while the consumer read
    only ``scenario_status``, which pools S4 into tested=21 with no family, size or
    correction. This feeds the producer's real output through the consumer.
    """
    status = read_status(FIXTURES / "sensitivity_status_c9fb6f1.json", "sensitivity_status", "#183")
    detail = set(status.detail)
    assert ("family", "S4:six_fractions", "correction", "Benjamini-Hochberg") in detail
    assert ("family", "S4:six_fractions", "family_size", "6") in detail
    assert ("family", "S4:six_fractions", "q_at_or_below_0.05", "2") in detail
    assert ("family", "S4:fifteen_ratios", "correction", "Benjamini-Yekutieli") in detail
    assert ("family", "S4:fifteen_ratios", "family_size", "15") in detail
    assert ("family", "S4:fifteen_ratios", "q_at_or_below_0.05", "4") in detail
    assert not any(value == "21" for _kind, _scope, _variable, value in detail)
    assert not any(kind == "scenario" for kind, _scope, _variable, _value in detail)
    rows = assessment_rows(
        apply_statuses(_built(), {**STATUSES, "sensitivity_status": status}),
        {**STATUSES, "sensitivity_status": status},
    )
    exported = {
        (str(r["scope"]), str(r["state_variable"]), str(r["state_value"]))
        for r in rows
        if r["assessment"] == "sensitivity_status"
    }
    assert ("S4:six_fractions", "correction", "Benjamini-Hochberg") in exported
    assert ("S4:fifteen_ratios", "family_size", "15") in exported


def test_a_quick_calibration_run_never_reaches_the_table_as_a_pass() -> None:
    """Cross-PR: #240's quick run (5 replicates, seed 1) at da8e1c1, end to end.

    Before the #240 fix this file said ``guild_mcnemar: "pass"`` and this table
    carried it as an unqualified pass on the McNemar row.
    """
    status = read_status(
        FIXTURES / "calibration_status_quick_run.json", "calibration_status", "#182"
    )
    for procedure in ("fraction_sign_test", "log_ratio_sign_test", "guild_mcnemar"):
        applied = status.for_procedure(procedure)
        assert applied.startswith("not_certified #182 ["), applied
        assert "insufficient_budget" in applied
    assert ("scenario", "complete_null_study_flip", "replicates_completed", "5") in status.detail
    records = apply_statuses(_built(), {**STATUSES, "calibration_status": status})
    assert not any(record.calibration_status.startswith("pass") for record in records)


def test_orientation_carries_class_orientation_and_validation() -> None:
    """Cross-PR: #238's provenance at 3d9d5cd, read for #188's fifth criterion."""
    status = read_status(
        FIXTURES / "contrast_orientation_provenance_3d9d5cd.json", "orientation_status", "#188"
    )
    detail = set(status.detail)
    assert ("contrast_class", "disease_vs_healthy:native", "experiments", "1003") in detail
    assert ("contrast_class", "disease_vs_healthy:reversed", "experiments", "1") in detail
    assert ("contrast_class", "disease_vs_healthy:native", "studies", "569") in detail
    assert ("validation", "hand_review", "disagreed", "17") in detail
    assert (
        "validation",
        "hand_review",
        "validation_status",
        "hand_review_disagreements_unresolved",
    ) in detail
    assert status.for_procedure("anything") == (
        "contrast population described, hand_review_disagreements_unresolved, no verdict #188"
    )


def _bound_status_payload(key: str, inputs: dict[str, str]) -> dict[str, object]:
    identity = {**inputs, "analysis_source_sha256": "a" * 64}
    if key == "calibration_status":
        return {
            "version": "calibration-v3",
            "code_commit": "historical-producer-revision",
            "inputs": {**identity, "protocol_sha256": "b" * 64},
            "protocol": "methods/calibration_protocol.json",
            "eligibility": {"id": "whole_pair_nonzero_retained_totals_v1"},
            "procedures": {
                "guild_mcnemar": {"eligibility_policy": "whole_pair_nonzero_retained_totals_v1"}
            },
            "calibration_status": {
                "fraction_sign_test": {"overall": "pass", "criteria": {"family_error": "pass"}},
                "guild_mcnemar": {"overall": "fail", "criteria": {"family_error": "fail"}},
            },
            "not_calibrated_here": ["native ORA and composition"],
        }
    return {
        **identity,
        "classifier_version": "contrast-orientation-v4",
        "code_commit": "historical-classifier-revision",
        "class_counts": [
            {"contrast_class": "unresolved", "orientation": "not_applicable", "experiments": 2}
        ],
        "review": {"reviewed_experiments": 60, "agreed": 58, "disagreed": 2},
    }


@pytest.mark.parametrize("key", ["calibration_status", "orientation_status"])
@pytest.mark.parametrize("bad_field", ["profile_sha256", "experiments_sha256", "full_dump_sha256"])
@pytest.mark.parametrize("corruption", ["mismatch", "missing", "malformed"])
def test_assessment_cannot_attach_verdicts_from_another_corpus(
    tmp_path: Path, key: str, bad_field: str, corruption: str
) -> None:
    expected = {
        "profile_sha256": "1" * 64,
        "experiments_sha256": "2" * 64,
        "full_dump_sha256": "3" * 64,
    }
    declared = dict(expected)
    if corruption == "missing":
        del declared[bad_field]
    else:
        declared[bad_field] = "4" * 64 if corruption == "mismatch" else "not-a-checksum"
    path = tmp_path / "assessment.json"
    path.write_text(json.dumps(_bound_status_payload(key, declared)), encoding="utf-8")
    status = read_status(path, key, "#729", expected_input_sha256=expected)
    assert status.verdicts == {}
    assert status.default == f"unreadable #729: assessment {bad_field} missing or mismatched"
    assert status.source["sha256"] == file_sha256(path)
    assert status.source["binding"] == "mismatch_or_missing_identity"


@pytest.mark.parametrize("key", ["calibration_status", "orientation_status"])
def test_bound_assessment_preserves_historical_identity_and_scoped_findings(
    tmp_path: Path, key: str
) -> None:
    expected = {
        "profile_sha256": "1" * 64,
        "experiments_sha256": "2" * 64,
        "full_dump_sha256": "3" * 64,
    }
    path = tmp_path / "assessment.json"
    payload = _bound_status_payload(key, expected)
    path.write_text(json.dumps(payload), encoding="utf-8")
    status = read_status(path, key, "#729", expected_input_sha256=expected)
    assert status.source["binding"] == "matched"
    assert status.source["input_sha256"] == expected
    assert status.source["code_commit"] == payload["code_commit"]
    assert status.source["sha256"] == file_sha256(path)
    assert "not an independent assessment" in str(status.source["validation_scope"])
    if key == "calibration_status":
        assert status.for_procedure("fraction_sign_test").startswith("pass #729")
        assert status.for_procedure("guild_mcnemar").startswith("fail #729")
        assert status.for_procedure("composition_arm_native_transformed_contrast") == (
            "not_calibrated_here #729"
        )
    else:
        assert "hand_review_disagreements_unresolved" in status.default
        assert "no verdict" in status.default


@pytest.mark.parametrize("key", ["calibration_status", "orientation_status"])
@pytest.mark.parametrize("missing", ["version", "code_commit", "analysis_source_sha256"])
def test_input_match_without_producer_identity_cannot_supply_a_verdict(
    tmp_path: Path, key: str, missing: str
) -> None:
    expected = {
        "profile_sha256": "1" * 64,
        "experiments_sha256": "2" * 64,
        "full_dump_sha256": "3" * 64,
    }
    payload = _bound_status_payload(key, expected)
    if missing == "analysis_source_sha256" and key == "calibration_status":
        inputs = payload["inputs"]
        assert isinstance(inputs, dict)
        del inputs[missing]
    else:
        field_name = missing
        if missing == "version" and key == "orientation_status":
            field_name = "classifier_version"
        del payload[field_name]
    path = tmp_path / "assessment.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    status = read_status(path, key, "#729", expected_input_sha256=expected)
    assert status.verdicts == {}
    assert status.default.startswith("unreadable #729: assessment lacks valid ")
    assert status.source["binding"] == "missing_producer_identity"


def test_calibration_without_protocol_identity_cannot_supply_a_verdict(tmp_path: Path) -> None:
    expected = dict.fromkeys(("profile_sha256", "experiments_sha256", "full_dump_sha256"), "1" * 64)
    payload = _bound_status_payload("calibration_status", expected)
    inputs = payload["inputs"]
    assert isinstance(inputs, dict)
    del inputs["protocol_sha256"]
    path = tmp_path / "assessment.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    status = read_status(path, "calibration_status", "#729", expected_input_sha256=expected)
    assert status.verdicts == {}
    assert status.default == "unreadable #729: assessment lacks valid protocol_sha256"


@pytest.mark.parametrize("missing", ["eligibility", "procedures"])
def test_matched_historical_calibration_does_not_certify_different_guild_eligibility(
    tmp_path: Path, missing: str
) -> None:
    expected = dict.fromkeys(("profile_sha256", "experiments_sha256", "full_dump_sha256"), "1" * 64)
    payload = _bound_status_payload("calibration_status", expected)
    payload["version"] = "calibration-v2"
    del payload[missing]
    path = tmp_path / "assessment.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    status = read_status(path, "calibration_status", "#182", expected_input_sha256=expected)
    assert status.source["binding"] == "matched"
    assert status.source["guild_eligibility_applicability"] == "not_established"
    assert status.for_procedure("guild_mcnemar").startswith("not_applicable #182")
    assert status.for_procedure("fraction_sign_test").startswith("pass #182")
    # The historical failure remains inspectable; it is not relabeled a pass.
    assert ("procedure", "guild_mcnemar", "overall", "fail") in status.detail


def test_cli_binds_imported_assessments_to_its_actual_inputs(tmp_path: Path) -> None:
    profiles, experiments = _corpus(tmp_path)
    dump = tmp_path / "dump.csv"
    dump.write_text("unversioned test dump\n", encoding="utf-8")
    expected = {
        "profile_sha256": file_sha256(profiles),
        "experiments_sha256": file_sha256(experiments),
        "full_dump_sha256": file_sha256(dump),
    }
    calibration = tmp_path / "calibration.json"
    calibration.write_text(
        json.dumps(_bound_status_payload("calibration_status", expected)), encoding="utf-8"
    )
    orientation = tmp_path / "orientation.json"
    orientation.write_text(
        json.dumps(
            _bound_status_payload("orientation_status", {**expected, "profile_sha256": "0" * 64})
        ),
        encoding="utf-8",
    )
    output = tmp_path / "integrated"
    assert (
        main(
            [
                "--profile-input",
                str(profiles),
                "--experiments-input",
                str(experiments),
                "--full-dump-input",
                str(dump),
                "--calibration-status",
                str(calibration),
                "--orientation-provenance",
                str(orientation),
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )
    provenance = json.loads((output / "integrated_evidence_provenance.json").read_text())
    statuses = provenance["assessment_status"]
    assert statuses["calibration_status"]["source"]["binding"] == "matched"
    assert statuses["orientation_status"]["source"]["binding"] == "mismatch_or_missing_identity"
    with (output / "integrated_evidence_method_results.tsv").open() as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    assert all(row["orientation_status"].startswith("unreadable #188") for row in records)
    fraction = next(
        row for row in records if row["method"] == "directional_consistency_fraction_sign_test"
    )
    assert fraction["calibration_status"].startswith("pass #182")
