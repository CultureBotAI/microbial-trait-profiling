"""Tests for the evidence overlays behind the Table 1 heatmap marks (#172, #176, #177)."""

from __future__ import annotations

import csv
import math
import subprocess
from pathlib import Path

import pytest
from scipy.stats import false_discovery_control, hypergeom, ttest_1samp

from bsdbng.trait_guilds import Profile
from bsdbng.trait_significance import (
    CELL_COMPOSITION,
    INSUFFICIENT_STUDIES,
    MAJORITY_SUPPORT,
    NO_MAJORITY,
    NOT_SELECTED,
    SITE_BELOW_MIN_PAIRS,
    SITE_COMPOSITION,
    SUPPORT,
    CompositionResult,
    CountTable,
    EnrichmentCall,
    EnrichmentRun,
    EvidenceOverlay,
    aggregate_support,
    composition_contrasts,
    enrichment_calls,
    read_count_table,
    write_evidence_tables,
)

TRAITS = {"trait_a": ("METPO:1",), "trait_b": ("METPO:2",)}
HEADER = ["study_id", "experiment_id", "inc_dec", "METPO:1", "METPO:2", "METPO:9"]

# (study, experiment, condition, trait_a increased, trait_a decreased). trait_b is
# 10 and the filler 60 in every profile. Study s1 spans two conditions with unequal
# experiment counts; s2 opposes it. This is the shape that the old grouping by
# (condition, study) counted as two studies.
PAIRS = [
    ("s1", "e1", "X", 30, 5),
    ("s1", "e2", "X", 25, 10),
    ("s1", "e3", "X", 40, 5),
    ("s1", "e4", "Y", 35, 8),
    ("s2", "e5", "X", 5, 30),
    ("s3", "e6", "X", 20, 12),
    ("s4", "e7", "X", 0, 15),
]


def _profiles(pairs: list[tuple[str, str, str, int, int]], site: str = "Feces") -> list[Profile]:
    profiles = []
    for study, experiment, condition, up, down in pairs:
        for direction, a in (("inc", up), ("dec", down)):
            profiles.append(
                Profile(
                    study_id=study,
                    experiment_id=experiment,
                    inc_dec=direction,
                    condition=condition,
                    body_site=site,
                    guild_counts={"trait_a": a, "trait_b": 10},
                    total=a + 70,
                )
            )
    return profiles


def _table(profiles: list[Profile]) -> CountTable:
    return CountTable(
        columns=("METPO:1", "METPO:2", "METPO:9"),
        counts={
            (p.study_id, p.experiment_id, p.inc_dec): (
                p.guild_counts["trait_a"],
                p.guild_counts["trait_b"],
                p.total - p.guild_counts["trait_a"] - p.guild_counts["trait_b"],
            )
            for p in profiles
        },
    )


def _write(path: Path, rows: list[list[str]]) -> Path:
    path.write_text("\n".join("\t".join(row) for row in [HEADER, *rows]) + "\n", encoding="utf-8")
    return path


def test_count_reader_keeps_only_requested_profiles_and_drops_mediadive(tmp_path: Path) -> None:
    path = tmp_path / "profiles.tsv"
    path.write_text(
        "study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2|mediadive.medium:J1\n"
        "s1\te1\tinc\t4\t99\n"
        "s2\te2\tinc\t7\t99\n",
        encoding="utf-8",
    )
    table = read_count_table(path, [("s1", "e1", "inc")])
    assert table.columns == ("METPO:1",)
    assert table.counts == {("s1", "e1", "inc"): (4,)}
    assert table.total(("s1", "e1", "inc")) == 4


@pytest.mark.parametrize(
    ("cell", "message"),
    [("", "missing is not zero"), ("-1", "non-negative"), ("2.5", "non-negative"), ("x", "holds")],
)
def test_count_reader_rejects_values_that_would_become_silent_zeros(
    tmp_path: Path, cell: str, message: str
) -> None:
    path = _write(tmp_path / "profiles.tsv", [["s1", "e1", "inc", cell, "1", "1"]])
    with pytest.raises(ValueError, match=message):
        read_count_table(path, [("s1", "e1", "inc")])


def test_count_reader_rejects_ragged_rows_repeats_and_absent_profiles(tmp_path: Path) -> None:
    short = _write(tmp_path / "short.tsv", [["s1", "e1", "inc", "1", "1"]])
    with pytest.raises(ValueError, match="fields"):
        read_count_table(short, [("s1", "e1", "inc")])
    repeated = _write(
        tmp_path / "repeated.tsv",
        [["s1", "e1", "inc", "1", "1", "1"], ["s1", "e1", "inc", "2", "2", "2"]],
    )
    with pytest.raises(ValueError, match="repeats"):
        read_count_table(repeated, [("s1", "e1", "inc")])
    absent = _write(tmp_path / "absent.tsv", [["s1", "e1", "inc", "1", "1", "1"]])
    with pytest.raises(ValueError, match="absent"):
        read_count_table(absent, [("s1", "e1", "inc"), ("s9", "e9", "dec")])


def test_enrichment_matches_the_hypergeometric_it_transcribes() -> None:
    """The p-value is scipy's one-tailed survival function on the corpus pool."""
    counts: dict[tuple[str, str, str], tuple[int, ...]] = {
        ("s1", "e1", "inc"): (90, 5, 5),
        ("s1", "e1", "dec"): (5, 90, 5),
        ("s2", "e2", "inc"): (10, 10, 80),
        ("s2", "e2", "dec"): (10, 10, 80),
        # No retained annotation at all: in neither the pool nor the tests.
        ("s3", "e3", "inc"): (0, 0, 0),
    }
    run = enrichment_calls(CountTable(("METPO:1", "METPO:2", "METPO:9"), counts), TRAITS)
    assert run.untested_traits == []
    assert run.profiles_tested == 4
    assert run.tests_performed == 12

    background_total = 400
    successes = 115
    first = next(
        c for c in run.calls if c.study_id == "s1" and c.inc_dec == "inc" and c.trait == "trait_a"
    )
    assert first.p_value == pytest.approx(float(hypergeom.sf(89, background_total, successes, 100)))
    assert first.fold == pytest.approx(90 / (100 * successes / background_total))
    all_p = [
        float(hypergeom.sf(count - 1, background_total, pool, 100))
        for values in list(counts.values())[:4]
        for count, pool in zip(values, (115, 115, 170), strict=True)
    ]
    assert first.q_value == pytest.approx(float(false_discovery_control(all_p, method="bh")[0]))


def test_there_is_no_depletion_call() -> None:
    """The arm has no depletion test, so an under-represented profile is simply not called."""
    call = EnrichmentCall("s1", "e1", "inc", "trait_a", 1, 50.0, 0.02, 0.0, 0.0)
    assert not call.enriched
    assert EnrichmentCall("s1", "e1", "inc", "trait_a", 90, 20.0, 4.5, 0.0, 0.0).enriched
    assert not EnrichmentCall("s1", "e1", "inc", "trait_a", 90, 20.0, 4.5, 0.0, 0.01).enriched


def _site(results: list[CompositionResult], trait: str) -> CompositionResult:
    return next(r for r in results if r.evidence_type == SITE_COMPOSITION and r.trait == trait)


def test_a_study_with_several_conditions_is_one_site_observation() -> None:
    """#176: s1's two conditions and four experiments enter the site test as one value."""
    profiles = _profiles(PAIRS)
    results = composition_contrasts(
        profiles, _table(profiles), TRAITS, min_site_pairs=1, min_profile_support=1
    )
    site = _site(results, "trait_a")
    assert (site.studies, site.paired_experiments, site.conditions) == (4, 7, 2)

    # Independent check: the t-statistic does not depend on the z-score scale, so
    # it must equal a t-test on per-study means of the unscaled log1p contrast.
    def scaled(a: int) -> float:
        return math.log1p(1000.0 * a / (a + 70))

    by_study: dict[str, list[float]] = {}
    for study, _experiment, _condition, up, down in PAIRS:
        by_study.setdefault(study, []).append(scaled(up) - scaled(down))
    means = [sum(v) / len(v) for v in by_study.values()]
    expected = ttest_1samp(means, 0.0)
    assert site.t_statistic == pytest.approx(float(expected.statistic))
    assert site.p_value == pytest.approx(float(expected.pvalue))
    assert site.family == "site:Feces"


def test_conditions_keep_their_own_contrasts() -> None:
    profiles = _profiles(PAIRS)
    results = composition_contrasts(
        profiles, _table(profiles), TRAITS, min_site_pairs=1, min_profile_support=1
    )
    cells = {
        r.condition: r
        for r in results
        if r.evidence_type == CELL_COMPOSITION and r.trait == "trait_a"
    }
    assert cells["X"].studies == 4
    assert cells["X"].paired_experiments == 6
    # One study is below the t-test floor and is recorded, not dropped.
    assert cells["Y"].status == INSUFFICIENT_STUDIES
    assert cells["Y"].p_value is None
    assert {r.family for r in results if r.evidence_type == CELL_COMPOSITION and r.p_value} == {
        "all_cells"
    }


def test_ubiquitous_traits_and_small_sites_are_recorded_not_tested() -> None:
    """The arm's own selection: a trait in more than 95% of profiles is not selected."""
    profiles = _profiles(PAIRS)
    results = composition_contrasts(
        profiles, _table(profiles), TRAITS, min_site_pairs=1, min_profile_support=1
    )
    assert _site(results, "trait_b").status == NOT_SELECTED
    below = composition_contrasts(profiles, _table(profiles), TRAITS, min_site_pairs=50)
    assert {r.status for r in below} == {SITE_BELOW_MIN_PAIRS}
    assert {r.trait for r in below} == set(TRAITS)


def test_zero_total_pairs_are_excluded_and_counted() -> None:
    profiles = _profiles(PAIRS)
    empty = [
        Profile("s5", "e8", direction, "X", "Feces", {"trait_a": 0, "trait_b": 0}, 0)
        for direction in ("inc", "dec")
    ]
    table = _table(profiles + empty)
    results = composition_contrasts(
        profiles + empty, table, TRAITS, min_site_pairs=1, min_profile_support=1
    )
    site = _site(results, "trait_a")
    assert site.zero_total_pairs_excluded == 1
    assert site.studies == 4


def test_the_six_trait_family_is_a_separate_sensitivity_column() -> None:
    profiles = _profiles(PAIRS) + _profiles(
        [(s, f"{e}b", c, u, d) for s, e, c, u, d in PAIRS], site="Saliva"
    )
    results = composition_contrasts(
        profiles, _table(profiles), TRAITS, min_site_pairs=1, min_profile_support=1
    )
    feces, saliva = (
        next(r for r in results if r.evidence_type == SITE_COMPOSITION and r.body_site == site)
        for site in ("Feces", "Saliva")
    )
    # Each site is its own screen family; the named-trait family spans both sites.
    assert feces.family == "site:Feces"
    assert feces.family_size == 1
    assert feces.table1_family_size == 2
    assert saliva.family != feces.family


def _run(calls: list[tuple[str, str, str, bool]]) -> EnrichmentRun:
    return EnrichmentRun(
        calls=[
            EnrichmentCall(s, e, d, "trait_a", 1, 1.0, 9.0 if hit else 1.0, 0.0, 0.0)
            for s, e, d, hit in calls
        ],
        untested_traits=[],
        profiles_tested=len(calls),
        tests_performed=len(calls),
    )


def test_support_counts_unique_studies_and_the_alternative_rule_removes_extra_chances() -> None:
    """#177: s1 has three chances and uses one; the stricter rule does not count it."""
    pairs = [
        ("s1", "e1", "X", 5, 5),
        ("s1", "e2", "X", 5, 5),
        ("s1", "e3", "X", 5, 5),
        ("s2", "e4", "X", 5, 5),
        ("s3", "e5", "X", 5, 5),
    ]
    profiles = _profiles(pairs)
    run = _run(
        [
            ("s1", "e1", "inc", True),
            ("s1", "e2", "inc", False),
            ("s1", "e3", "inc", False),
            ("s2", "e4", "inc", True),
            ("s3", "e5", "inc", False),
        ]
    )
    records = aggregate_support(profiles, _table(profiles), run, ["trait_a"])
    cell = records[("X", "Feces", "trait_a", "inc")]
    assert (cell.eligible_studies, cell.supporting_studies) == (3, 2)
    assert cell.status == MAJORITY_SUPPORT
    assert cell.alternative_supporting_studies == 1
    assert cell.alternative_status == NO_MAJORITY
    assert (cell.eligible_profiles, cell.enriched_profiles) == (5, 2)
    assert (cell.studies_with_repeated_profiles, cell.max_profiles_per_study) == (1, 3)
    assert not hasattr(cell, "q_value")
    assert not hasattr(cell, "p_value")


def test_one_study_never_speaks_for_a_cell() -> None:
    profiles = _profiles([("s1", "e1", "X", 5, 5)])
    records = aggregate_support(
        profiles, _table(profiles), _run([("s1", "e1", "inc", True)]), ["trait_a"]
    )
    assert records[("X", "Feces", "trait_a", "inc")].status == INSUFFICIENT_STUDIES


def _overlay(results: list[CompositionResult], support_hit: bool = True) -> EvidenceOverlay:
    pairs = [("s1", "e1", "X", 5, 5), ("s2", "e2", "X", 5, 5)]
    profiles = _profiles(pairs)
    run = _run([("s1", "e1", "inc", support_hit), ("s2", "e2", "inc", support_hit)])
    return EvidenceOverlay(
        support=aggregate_support(profiles, _table(profiles), run, ["trait_a"]),
        composition=results,
        enrichment=run,
        traits=("trait_a",),
    )


def _result(evidence_type: str, condition: str | None, mean: float) -> CompositionResult:
    return CompositionResult(
        evidence_type=evidence_type,
        body_site="Feces",
        condition=condition,
        trait="trait_a",
        column="METPO:1",
        status="tested",
        studies=5,
        paired_experiments=5,
        conditions=1,
        zero_total_pairs_excluded=0,
        mean_contrast_z=mean,
        t_statistic=3.0,
        p_value=0.001,
        family="f",
        family_size=1,
        q_value=0.001,
        table1_family_size=1,
        table1_q_value=0.001,
    )


def test_site_results_are_never_promoted_to_cells() -> None:
    overlay = _overlay([_result(SITE_COMPOSITION, None, 0.4)], support_hit=False)
    assert overlay.cell_marks() == {}
    assert overlay.site_marks() == {("Feces", "trait_a"): "increased"}


def test_a_condition_contrast_marks_only_the_side_it_favours() -> None:
    overlay = _overlay([_result(CELL_COMPOSITION, "X", -0.4)])
    marks = overlay.cell_marks()
    assert marks[("X", "Feces", "trait_a", "dec")] == frozenset({CELL_COMPOSITION})
    assert marks[("X", "Feces", "trait_a", "inc")] == frozenset({SUPPORT})


def test_written_support_carries_rule_counts_and_no_p_or_q(tmp_path: Path) -> None:
    overlay = _overlay([_result(SITE_COMPOSITION, None, 0.4), _result(CELL_COMPOSITION, "X", 0.2)])
    paths = write_evidence_tables(
        overlay,
        tmp_path,
        "t",
        display_status=lambda cell: "drawn" if cell[3] == "inc" else "below_min_studies",
        site_drawn=lambda site: site == "Feces",
        provenance={"profile_sha256": "abc"},
    )
    support_path, contrasts_path, screen_path, provenance_path = paths
    with support_path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert not {"p_value", "q_value"} & set(rows[0])
    by_direction = {row["direction"]: row for row in rows}
    assert by_direction["inc"]["support_status"] == MAJORITY_SUPPORT
    assert by_direction["inc"]["display_status"] == "drawn"
    assert by_direction["dec"]["display_status"] == "below_min_studies"
    with contrasts_path.open(encoding="utf-8") as handle:
        kinds = {row["evidence_type"] for row in csv.DictReader(handle, delimiter="\t")}
    assert kinds == {SITE_COMPOSITION, CELL_COMPOSITION}
    assert screen_path.exists()
    text = provenance_path.read_text(encoding="utf-8")
    assert '"cell_p_or_q": null' in text
    assert "composition_departures" in text


# The enrichment arm's committed run this module reproduces.
ENRICHMENT_ARM_COMMIT = "c75744a"


@pytest.mark.integration
def test_reproduces_the_enrichment_arm_committed_counts() -> None:
    """Orthogonal check against the pinned table this test's subject was transcribed from.

    The arm's output at a fixed commit was produced by a different script from the
    same corpus. If this module's hypergeometric is the arm's, the two agree on how
    many profiles each Table 1 trait is enriched in. A missing commit fails rather
    than skips, so the check cannot pass by not running.
    """
    from bsdbng.trait_guilds import (
        DEFAULT_EXPERIMENTS_PATH,
        DEFAULT_PROFILE_PATH,
        TABLE1_TRAITS,
        load_profiles,
    )

    committed = subprocess.run(
        ["git", "show", f"{ENRICHMENT_ARM_COMMIT}:data/derived/trait_enrichment_results.tsv"],
        capture_output=True,
        text=True,
        check=False,
    )
    if committed.returncode != 0:
        pytest.fail(
            f"commit {ENRICHMENT_ARM_COMMIT} is not available; run "
            "`git fetch origin enrichment-analysis-heatmap-cluster`"
        )

    theirs: dict[str, int] = {}
    for line in committed.stdout.splitlines():
        if line.startswith(("#", "study_id")):
            continue
        trait = line.split("\t")[3]
        theirs[trait] = theirs.get(trait, 0) + 1

    profiles, _ = load_profiles(DEFAULT_PROFILE_PATH, DEFAULT_EXPERIMENTS_PATH, TABLE1_TRAITS)
    keys = [(p.study_id, p.experiment_id, p.inc_dec) for p in profiles]
    run = enrichment_calls(read_count_table(DEFAULT_PROFILE_PATH, keys), TABLE1_TRAITS)

    ours: dict[str, int] = {}
    for call in run.calls:
        if call.enriched:
            ours[call.trait] = ours.get(call.trait, 0) + 1

    for name, ids in TABLE1_TRAITS.items():
        term = ids[0]
        theirs_count = next(
            (count for key, count in theirs.items() if key.rsplit("|", 1)[-1] == term), None
        )
        assert theirs_count is not None, f"{term} absent from the committed table"
        assert ours.get(name, 0) == theirs_count, name
