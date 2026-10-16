"""The widened directional screen under the complete study-flip null (#263)."""

from __future__ import annotations

from mtp.calibration import SEED, SHARED_FLIP_REPLICATES, flip_assignments
from mtp.screen_calibration import (
    CORRECTIONS,
    STREAM,
    THRESHOLDS,
    ReplicateCount,
    count_calls,
    replicate_flipped_profiles,
    run_null_screens,
    shared_flip_artefact_provenance,
    summarize,
)
from mtp.trait_guilds import Profile
from mtp.trait_screen import run_screen

COLUMNS = ["A", "B", "C"]


def _profile(study: str, direction: str, counts: dict[str, int]) -> Profile:
    return Profile(study, f"{study}-1", direction, "X", "Feces", counts, sum(counts.values()))


def _shifted_corpus(studies: int = 30) -> list[Profile]:
    """Every study shifts the same way on all three traits, so the screen calls all three."""
    return [
        profile
        for index in range(studies)
        for profile in (
            _profile(f"s{index}", "inc", {"A": 9, "B": 5 + index % 3, "C": 4}),
            _profile(f"s{index}", "dec", {"A": 1, "B": 5, "C": 6}),
        )
    ]


def test_every_correction_and_threshold_is_counted_once() -> None:
    counts = count_calls(run_screen(_shifted_corpus(), COLUMNS), replicate=0)
    keys = [(c.q_column, c.threshold) for c in counts]
    assert len(keys) == len(set(keys)) == len(CORRECTIONS) * len(THRESHOLDS)
    by_key = {(c.q_column, c.threshold): c for c in counts}
    assert by_key[("by_q_nominal", 0.05)].significant == 3
    assert by_key[("by_q_nominal", 0.05)].family_size == len(COLUMNS)


def test_the_null_replicates_are_the_shared_artefacts_flips() -> None:
    """Replicate i here must flip exactly the studies replicate i of the artefact flips."""
    # Synthetic corpus with multiple experiments in some studies (#339)
    profiles = [
        _profile("s0", "inc", {"A": 9, "B": 5, "C": 4}),
        _profile("s0", "dec", {"A": 1, "B": 5, "C": 6}),
        Profile("s0", "s0-2", "inc", "X", "Feces", {"A": 8, "B": 4, "C": 3}, 15),
        Profile("s0", "s0-2", "dec", "X", "Feces", {"A": 2, "B": 4, "C": 5}, 11),
        _profile("s1", "inc", {"A": 7, "B": 6, "C": 5}),
        _profile("s1", "dec", {"A": 3, "B": 6, "C": 4}),
    ]
    artefact = flip_assignments(profiles, replicates=3, seed=SEED, stream=STREAM)
    for index in range(3):
        flipped = replicate_flipped_profiles(profiles, index)
        changed = {
            before.study_id
            for before, after in zip(profiles, flipped, strict=True)
            if before.inc_dec != after.inc_dec
        }
        listed = {
            str(row["study_id"])
            for row in artefact
            if row["replicate"] == index and row["flipped"] == 1
        }
        assert changed == listed


def test_flipping_whole_studies_destroys_a_perfectly_consistent_signal() -> None:
    """With every study shifting the same way, the null screens must call far less."""
    profiles = _shifted_corpus()
    observed = count_calls(run_screen(profiles, COLUMNS), replicate=-1)
    seen = {(c.q_column, c.threshold): c.significant for c in observed}
    counts = run_null_screens(profiles, COLUMNS, replicates=20, seed=1)
    assert {c.replicate for c in counts} == set(range(20))
    rows = summarize(counts, seen)
    nominal = next(r for r in rows if r["q_column"] == "by_q_nominal" and r["threshold"] == 0.05)
    assert nominal["observed_calls"] == 3
    assert nominal["replicates"] == 20
    assert float(str(nominal["mean_null_calls"])) < 3
    per_rep = len(CORRECTIONS) * len(THRESHOLDS)
    assert run_null_screens(profiles, COLUMNS, replicates=3, seed=1) == counts[: 3 * per_rep]


def test_the_summary_reports_family_error_with_an_exact_interval() -> None:
    counts = [
        ReplicateCount(index, "by_q_nominal", "Benjamini-Yekutieli", "nominal", 3, 0.05, value)
        for index, value in enumerate([0, 0, 1, 0, 3])
    ]
    (row,) = summarize(counts, {("by_q_nominal", 0.05): 2})
    assert row["replicates_with_any_call"] == 2
    assert row["family_error"] == 0.4
    assert float(str(row["family_error_ci_low"])) < 0.4 < float(str(row["family_error_ci_high"]))
    assert row["replicates_reaching_observed"] == 1
    assert row["mean_null_calls"] == 0.8
    assert row["max_null_calls"] == 3

    # When observed calls are 0, reaching observed is None (not applicable) (#342)
    (zero_row,) = summarize(counts, {("by_q_nominal", 0.05): 0})
    assert zero_row["replicates_reaching_observed"] is None


def test_shared_flip_artefact_provenance_records_seed_and_replicates() -> None:
    matched = shared_flip_artefact_provenance(SEED, SHARED_FLIP_REPLICATES)
    assert matched == {
        "shared_flip_artefact_seed": SEED,
        "shared_flip_artefact_replicates": SHARED_FLIP_REPLICATES,
        "matches_shared_flip_artefact": True,
    }
    diff_seed = shared_flip_artefact_provenance(SEED + 1, SHARED_FLIP_REPLICATES)
    assert diff_seed["matches_shared_flip_artefact"] is False
    too_many = shared_flip_artefact_provenance(SEED, SHARED_FLIP_REPLICATES + 1)
    assert too_many["matches_shared_flip_artefact"] is False
    zero_reps = shared_flip_artefact_provenance(SEED, 0)
    assert zero_reps["matches_shared_flip_artefact"] is False
