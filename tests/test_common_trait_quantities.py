"""Common study-weighted quantities and denominator/reference invariants (#256)."""

from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

import pytest

from bsdbng.common_trait_quantities import (
    FRACTION_FULL,
    FRACTION_MATCHED,
    LOG_RATIO,
    Pair,
    QuantitySummary,
    audit_denominators,
    build_common_quantities,
    inspect_trait_columns,
    log_ratio_contrast,
)
from bsdbng.trait_guilds import Profile, load_profiles

TRAITS = ("A", "B")


def _profile(
    study: str,
    experiment: str,
    direction: str,
    a: int,
    b: int,
    *,
    total: int | None = None,
) -> Profile:
    counts = {"A": a, "B": b}
    return Profile(
        study_id=study,
        experiment_id=experiment,
        inc_dec=direction,
        condition="condition",
        body_site="Feces",
        guild_counts=counts,
        total=a + b if total is None else total,
    )


def _pair(
    study: str,
    experiment: str,
    inc: tuple[int, int],
    dec: tuple[int, int],
    *,
    inc_total: int | None = None,
    dec_total: int | None = None,
) -> list[Profile]:
    return [
        _profile(study, experiment, "inc", *inc, total=inc_total),
        _profile(study, experiment, "dec", *dec, total=dec_total),
    ]


def _summary(
    rows: tuple[QuantitySummary, ...], quantity: str, numerator: str, reference: str = ""
) -> QuantitySummary:
    return next(
        row
        for row in rows
        if row.quantity_type == quantity
        and row.numerator == numerator
        and row.reference == reference
    )


def test_denominator_only_change_moves_fixed_focal_fractions_but_not_ratio() -> None:
    profiles: list[Profile] = []
    for study in range(4):
        # A and B are fixed; an unrelated retained component changes the total.
        profiles += _pair(
            f"s{study}",
            "e1",
            (10, 20),
            (10, 20),
            inc_total=60 + 10 * study,
            dec_total=40,
        )
    result = build_common_quantities(profiles, TRAITS)
    a = _summary(result.summaries, FRACTION_FULL, "A")
    b = _summary(result.summaries, FRACTION_FULL, "B")
    ratio = _summary(result.summaries, LOG_RATIO, "A", "B")

    assert a.mean_study_effect is not None
    assert a.mean_study_effect < 0
    assert b.mean_study_effect is not None
    assert b.mean_study_effect < 0
    assert ratio.mean_study_effect == pytest.approx(0.0)
    assert ratio.tied_studies == 4
    assert ratio.proportion_positive_non_tied is None


def _write_loader_fixture(path: Path, media: tuple[int, int]) -> tuple[Path, Path]:
    path.mkdir(parents=True, exist_ok=True)
    profiles = path / "profiles.tsv"
    experiments = path / "experiments.tsv"
    profiles.write_text(
        "study_id\texperiment_id\tinc_dec\tA\tB\tother\tmediadive.medium:1\n"
        f"s1\te1\tinc\t4\t2\t5\t{media[0]}\n"
        f"s1\te1\tdec\t2\t2\t5\t{media[1]}\n",
        encoding="utf-8",
    )
    experiments.write_text(
        "study_id\texperiment_id\tcondition\tbody_site\thost_species\n"
        "s1\te1\tcondition\tFeces\tHomo sapiens\n",
        encoding="utf-8",
    )
    return profiles, experiments


def test_excluded_mediadive_changes_neither_denominator_fractions_nor_ratio(
    tmp_path: Path,
) -> None:
    left_paths = _write_loader_fixture(tmp_path / "left", (0, 0))
    right_paths = _write_loader_fixture(tmp_path / "right", (10_000, 1))
    guilds: dict[str, tuple[str, ...]] = {"A": ("A",), "B": ("B",)}
    left, left_stats = load_profiles(*left_paths, guilds, reject_missing_counts=True)
    right, right_stats = load_profiles(*right_paths, guilds, reject_missing_counts=True)
    left_result = build_common_quantities(left, TRAITS)
    right_result = build_common_quantities(right, TRAITS)

    assert [asdict(row) for row in left_result.summaries] == [
        asdict(row) for row in right_result.summaries
    ]
    assert [profile.total for profile in left] == [11, 9]
    assert [profile.total for profile in right] == [11, 9]
    assert left_stats.mediadive_columns_excluded == right_stats.mediadive_columns_excluded == 1


def test_reference_reversal_reverses_effect_without_changing_population() -> None:
    profiles = _pair("s1", "e1", (8, 2), (2, 2))
    forward = build_common_quantities(profiles, ("A", "B"))
    reverse = build_common_quantities(profiles, ("B", "A"))
    ab = _summary(forward.summaries, LOG_RATIO, "A", "B")
    reverse_ratio = _summary(reverse.summaries, LOG_RATIO, "B", "A")

    assert ab.mean_study_effect == pytest.approx(2.0)
    assert reverse_ratio.mean_study_effect == pytest.approx(-2.0)
    assert ab.eligible_experiment_pairs == reverse_ratio.eligible_experiment_pairs == 1
    pair = Pair("s1", "e1", "condition", "Feces", profiles[0], profiles[1])
    assert log_ratio_contrast(pair, "A", "B") == -log_ratio_contrast(pair, "B", "A")


def test_unequal_experiments_are_averaged_within_study_before_equal_weighting() -> None:
    profiles: list[Profile] = []
    for experiment in range(4):
        inc = (2, 1) if experiment < 3 else (1, 16)
        profiles += _pair("s1", f"e{experiment}", inc, (1, 1))
    profiles += _pair("s2", "e1", (2, 1), (1, 1))
    result = build_common_quantities(profiles, TRAITS)
    ratio = _summary(result.summaries, LOG_RATIO, "A", "B")
    ratio_studies = [row for row in result.study_effects if row.row_id == ratio.row_id]

    assert [row.mean_study_effect for row in ratio_studies] == pytest.approx([-0.25, 1.0])
    assert ratio.mean_study_effect == pytest.approx(0.375)
    assert (ratio.positive_studies, ratio.negative_studies) == (1, 1)


def test_zero_total_and_zero_named_count_policies_remain_distinct() -> None:
    profiles = [
        *_pair("s1", "e1", (0, 0), (1, 1), inc_total=0, dec_total=2),
        # An observed focal zero with positive retained totals is valid for fractions.
        *_pair("s2", "e1", (0, 1), (0, 1), inc_total=5, dec_total=5),
    ]
    result = build_common_quantities(profiles, TRAITS)
    fraction = _summary(result.summaries, FRACTION_FULL, "A")
    ratio = _summary(result.summaries, LOG_RATIO, "A", "B")

    assert fraction.eligible_experiment_pairs == 1
    assert fraction.excluded_zero_total_pairs == 1
    assert fraction.tied_studies == 1
    assert ratio.eligible_experiment_pairs == 0
    assert ratio.excluded_zero_named_count_pairs == 2
    assert ratio.status == "no_eligible_pairs"


def test_matched_fractions_and_ratio_point_to_the_same_exact_pair_manifest() -> None:
    profiles = [
        *_pair("s1", "e1", (4, 2), (2, 2)),
        *_pair("s2", "e1", (0, 2), (2, 2)),
    ]
    result = build_common_quantities(profiles, TRAITS)
    ratio = _summary(result.summaries, LOG_RATIO, "A", "B")
    matched = [
        row
        for row in result.summaries
        if row.quantity_type == FRACTION_MATCHED and row.population_id == ratio.population_id
    ]
    population = next(
        item for item in result.populations if item.population_id == ratio.population_id
    )

    assert len(matched) == 2
    assert {row.eligible_experiment_pairs for row in matched} == {1}
    assert {row.eligible_unique_studies for row in matched} == {1}
    assert [pair.pair_id for pair in population.pairs] == ["s1|e1"]


@pytest.mark.parametrize(
    ("raw", "message"),
    [("", "Missing trait count"), ("unknown", "Malformed trait count"), ("-1", "Negative trait")],
)
def test_missing_malformed_and_negative_counts_fail_before_eligibility(
    tmp_path: Path, raw: str, message: str
) -> None:
    paths = _write_loader_fixture(tmp_path, (0, 0))
    profile_path, experiments_path = paths
    text = profile_path.read_text(encoding="utf-8").replace("inc\t4\t", f"inc\t{raw}\t")
    profile_path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        load_profiles(
            profile_path,
            experiments_path,
            {"A": ("A",), "B": ("B",)},
            reject_missing_counts=True,
        )


def _write_reference_audit(path: Path, rows: list[tuple[str, int, int, bool, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "profile_id",
                "legacy_total",
                "harmonized_total",
                "mediadive_total",
                "ilyas_total",
                "noel_total",
                "paired",
            ]
        )
        for profile_id, harmonized, media, paired, noel in rows:
            writer.writerow(
                [profile_id, harmonized + media, harmonized, media, harmonized, noel, paired]
            )
    return path


def test_denominator_audit_separates_total_parity_from_population_mismatch(
    tmp_path: Path,
) -> None:
    profiles = _pair("s1", "e1", (2, 1), (1, 1))
    reference = _write_reference_audit(
        tmp_path / "audit.tsv",
        [
            ("s1|e1|inc", 3, 100, False, "3"),
            ("s1|e1|dec", 2, 20, True, "2"),
            ("s1|wrong|inc", 9, 1, True, "NA"),
        ],
    )
    audit = audit_denominators(profiles, reference, "revision")
    by_check = {str(row["check"]): row for row in audit.summary}

    assert audit.ok
    assert by_check["current_loader_vs_pr197_harmonized"]["mismatches"] == 0
    assert by_check["pr197_paired_membership_vs_current_shared_base"]["mismatches"] == 2
    assert by_check["pr197_paired_rows_missing_noel_total"]["mismatches"] == 1
    assert all(not bool(row["blocking"]) for row in audit.mismatches)


def test_denominator_value_mismatch_is_blocking(tmp_path: Path) -> None:
    profiles = _pair("s1", "e1", (2, 1), (1, 1))
    reference = _write_reference_audit(
        tmp_path / "audit.tsv",
        [
            ("s1|e1|inc", 99, 0, True, "99"),
            ("s1|e1|dec", 2, 0, True, "2"),
        ],
    )
    audit = audit_denominators(profiles, reference, "revision")
    assert not audit.ok
    assert any(
        row["check"] == "current_loader_vs_pr197_harmonized" and bool(row["blocking"])
        for row in audit.mismatches
    )


def test_column_manifest_excludes_media_and_resolves_qualified_terms(tmp_path: Path) -> None:
    profile = tmp_path / "profiles.tsv"
    profile.write_text(
        "study_id\texperiment_id\tinc_dec\tMETPO:1000608\tMETPO:2000103|GO:0030245\t"
        "other\tmediadive.medium:1\n",
        encoding="utf-8",
    )
    retained, identities = inspect_trait_columns(
        profile,
        {
            "facultatively_aerobic": ("METPO:1000608",),
            "cellulose_catabolic_process": ("GO:0030245",),
        },
    )
    assert retained == ["METPO:1000608", "METPO:2000103|GO:0030245", "other"]
    assert [identity.source_column for identity in identities] == retained[:2]


def test_two_trait_contract_emits_two_full_two_matched_and_one_ratio_rows() -> None:
    result = build_common_quantities(_pair("s1", "e1", (2, 1), (1, 1)), TRAITS)
    assert len(result.summaries) == 5
    assert [row.quantity_type for row in result.summaries].count(FRACTION_FULL) == 2
    assert [row.quantity_type for row in result.summaries].count(FRACTION_MATCHED) == 2
    assert [row.quantity_type for row in result.summaries].count(LOG_RATIO) == 1
