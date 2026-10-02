"""The integrated evidence table: states, families and what a join must not lose (#178)."""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import pytest

from bsdbng.guild_maps import build_matrix, cell_status_rows, compute_correlation
from bsdbng.integrated_evidence import (
    ALL_TIED,
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
    Record,
    apply_statuses,
    build_guild_matrix,
    cell_support_records,
    family_inventory,
    feature_column_rows,
    fraction_sign_records,
    guild_records,
    is_cell_row,
    log_ratio_records,
    placeholder_records,
    read_status,
    summarize,
    unmatched_terms,
    write_table,
    write_view,
)
from bsdbng.table1_trait_maps import compute_direction_shifts
from bsdbng.trait_composition import compute_log_ratio_shifts
from bsdbng.trait_guilds import (
    OXYGEN_GUILDS,
    TABLE1_TRAITS,
    Profile,
    drop_zero_total_pairs,
    resolve_trait_columns,
)

#: Both decision thresholds every record builder now takes (#294).
THRESHOLDS = {"q_threshold": 0.05, "strict_q_threshold": 0.001}

TRAITS = tuple(TABLE1_TRAITS)
STATUSES = {
    "calibration_status": "fraction_sign_test: pass",
    "sensitivity_status": "8 scenarios reported",
    "coverage_status": "pending #184",
    "orientation_status": "1052 classified",
    "taxon_influence_status": "pending #190",
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
        assert after.calibration_status == "fraction_sign_test: pass"
        assert after.coverage_status == "pending #184"
    with pytest.raises(ValueError, match="no status supplied"):
        apply_statuses(records, {"calibration_status": "x"})


def test_absent_assessment_files_read_as_pending(tmp_path: Path) -> None:
    assert read_status(None, "calibration_status", "#182") == "pending #182"
    assert read_status(tmp_path / "missing.json", "coverage_status", "#184") == "pending #184"
    broken = tmp_path / "broken.json"
    broken.write_text("not json", encoding="utf-8")
    assert read_status(broken, "coverage_status", "#184") == "unreadable #184"
    written = tmp_path / "calibration.json"
    written.write_text('{"calibration_status": {"a": "pass"}}', encoding="utf-8")
    assert read_status(written, "calibration_status", "#182") == "a: pass"


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

    # Kept as a tie, the blank study would inflate the tie count; excluded, the
    # tallies match the corpus that has a defined share.
    with_blank = fraction_sign_records(scored + blank, source="x", columns={}, **THRESHOLDS)
    without = fraction_sign_records(exclusion.profiles, source="x", columns={}, **THRESHOLDS)
    assert {r.trait: r.ties for r in with_blank} != {r.trait: r.ties for r in without}


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
