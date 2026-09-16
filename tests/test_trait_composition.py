"""Counterexamples to fraction inference and invariants of paired log-ratios."""

import csv
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

from bsdbng.table1_trait_maps import compute_direction_shifts
from bsdbng.trait_composition import compute_log_ratio_shifts, main
from bsdbng.trait_guilds import TABLE1_TRAITS, Profile


def _profile(
    study: str,
    direction: str,
    a: int,
    b: int,
    c: int = 1,
    *,
    experiment: str = "e1",
    condition: str = "X",
) -> Profile:
    return Profile(
        study, experiment, direction, condition, "Feces", {"A": a, "B": b, "C": c}, a + b + c
    )


def test_closure_moves_unchanged_traits_together_but_not_their_ratio() -> None:
    profiles = [
        profile
        for i in range(1, 9)
        for profile in (
            _profile(str(i), "inc", 10, 20, c=10 + 10 * i),
            _profile(str(i), "dec", 10, 20, c=10),
        )
    ]
    fractions = compute_direction_shifts(profiles, ("A", "B", "C"))
    a, b, c = fractions.shifts
    assert a.studies_down == b.studies_down == c.studies_up == 8
    assert a.p_value == pytest.approx(1 / 128)
    # Closure permits perfectly POSITIVE correlation of unchanged components.
    assert fractions.spearman["A~B"][0] == pytest.approx(1.0)

    (ratio,) = compute_log_ratio_shifts(profiles, ("A", "B")).summaries
    assert ratio.studies_tied == 8
    assert ratio.mean_log2_ratio_change == pytest.approx(0.0)
    assert ratio.status == "all_tied"
    assert ratio.share_ratio_up is ratio.sign_p is ratio.by_q is None


def test_true_ratio_shift_has_known_effect_interval_and_full_family_adjustment() -> None:
    profiles = [
        profile
        for i in range(8)
        for profile in (_profile(str(i), "inc", 4, 1), _profile(str(i), "dec", 1, 1))
    ]
    result = compute_log_ratio_shifts(profiles, ("A", "B", "C"))
    ab, ac, bc = result.summaries
    assert result.family_size == 3
    for shifted in (ab, ac):
        assert shifted.mean_log2_ratio_change == shifted.median_log2_ratio_change == 2.0
        assert shifted.sign_p == pytest.approx(2 / 2**8)
        assert shifted.share_ci_low == pytest.approx(0.025 ** (1 / 8))
        assert shifted.share_ci_high == 1.0
        # BY for [1/128, 1/128, 1]: p * (3/2) * (1 + 1/2 + 1/3).
        assert shifted.by_q == pytest.approx(11 / 512)
    assert bc.status == "all_tied"
    assert bc.sign_p is bc.by_q is None


def test_ratio_is_invariant_to_profile_scaling_and_unrelated_components() -> None:
    profiles = [
        profile
        for i in range(1, 5)
        for profile in (_profile(str(i), "inc", 2 + i, 2), _profile(str(i), "dec", 2, 2))
    ]
    changed = []
    for i, profile in enumerate(profiles):
        scale = 2 ** (i + 1)
        counts = {trait: count * scale for trait, count in profile.guild_counts.items()}
        counts["C"] += 1000 * (i + 1)
        changed.append(replace(profile, guild_counts=counts, total=sum(counts.values())))
    before = compute_log_ratio_shifts(profiles, ("A", "B"))
    after = compute_log_ratio_shifts(changed, ("A", "B"))
    assert [s.mean_log2_ratio_change for s in after.studies] == pytest.approx(
        [s.mean_log2_ratio_change for s in before.studies]
    )
    left, right = before.summaries[0], after.summaries[0]
    assert right.mean_log2_ratio_change == pytest.approx(left.mean_log2_ratio_change)
    assert (right.studies_up, right.studies_down, right.sign_p, right.by_q) == (
        left.studies_up,
        left.studies_down,
        left.sign_p,
        left.by_q,
    )


def test_reference_reversal_reverses_effect_and_interval_but_not_two_sided_p() -> None:
    profiles = [
        profile
        for i in range(8)
        for profile in (
            _profile(str(i), "inc", 4 if i < 6 else 1, 1),
            _profile(str(i), "dec", 1 if i < 6 else 2, 1),
        )
    ]
    (forward,) = compute_log_ratio_shifts(profiles, ("A", "B")).summaries
    (reverse,) = compute_log_ratio_shifts(profiles, ("B", "A")).summaries
    assert (forward.studies_up, forward.studies_down) == (6, 2)
    assert (reverse.studies_up, reverse.studies_down) == (2, 6)
    assert forward.mean_log2_ratio_change == 1.25
    assert reverse.mean_log2_ratio_change == -1.25
    assert forward.sign_p == reverse.sign_p == pytest.approx(74 / 256)
    assert forward.by_q == reverse.by_q
    assert forward.share_ci_low is not None
    assert forward.share_ci_high is not None
    assert reverse.share_ci_low == pytest.approx(1 - forward.share_ci_high)
    assert reverse.share_ci_high == pytest.approx(1 - forward.share_ci_low)


def test_each_study_votes_once_across_experiments_and_conditions() -> None:
    profiles = [
        profile
        for i in range(4)
        for profile in (
            _profile(
                "s1",
                "inc",
                2 if i < 3 else 1,
                1 if i < 3 else 16,
                experiment=str(i),
                condition=f"condition {i}",
            ),
            _profile("s1", "dec", 1, 1, experiment=str(i), condition=f"condition {i}"),
        )
    ]
    profiles.extend([_profile("s2", "inc", 2, 1), _profile("s2", "dec", 1, 1)])
    result = compute_log_ratio_shifts(profiles, ("A", "B"))
    assert result.studies[0].mean_log2_ratio_change == -0.25
    summary = result.summaries[0]
    # Four positive experiments and one negative; two studies split one each.
    assert (summary.studies_up, summary.studies_down, summary.studies_total) == (1, 1, 2)
    assert summary.mean_log2_ratio_change == 0.375
    assert summary.sign_p == 1.0


def test_zero_exclusions_preserve_whole_pairs_and_study_coverage() -> None:
    profiles = [
        _profile("s1", "inc", 4, 1),
        _profile("s1", "dec", 1, 1),
        _profile("s1", "inc", 0, 1, experiment="e2"),
        _profile("s1", "dec", 1024, 1, experiment="e2"),
        _profile("s2", "inc", 1, 0),
        _profile("s2", "dec", 1, 1),
    ]
    result = compute_log_ratio_shifts(profiles, ("A", "B"))
    summary = result.summaries[0]
    assert (
        summary.paired_experiments,
        summary.eligible_experiments,
        summary.excluded_zero_experiments,
    ) == (3, 1, 2)
    assert (
        summary.studies_total,
        summary.studies_eligible,
        summary.studies_excluded,
        summary.studies_partially_observed,
    ) == (2, 1, 1, 1)
    assert result.studies[0].mean_log2_ratio_change == 2.0
    assert result.studies[1].mean_log2_ratio_change is None

    (empty,) = compute_log_ratio_shifts(profiles[-2:], ("A", "B")).summaries
    assert empty.status == "no_positive_pairs"
    assert empty.mean_log2_ratio_change is None
    assert empty.share_ratio_up is empty.sign_p is empty.by_q is None


@pytest.mark.parametrize(
    ("defect", "message"),
    [
        ("duplicate", "Duplicate profile"),
        ("unpaired", "Unpaired experiment"),
        ("missing", "Missing trait"),
        ("negative", "Negative count"),
        ("metadata", "Mismatched experiment metadata"),
    ],
)
def test_invalid_input_does_not_silently_change_the_sample(defect: str, message: str) -> None:
    profiles = [_profile("s1", "inc", 2, 1), _profile("s1", "dec", 1, 1)]
    if defect == "duplicate":
        profiles.append(profiles[0])
    elif defect == "unpaired":
        profiles.pop()
    elif defect == "missing":
        profiles[0] = replace(profiles[0], guild_counts={"A": 2})
    elif defect == "negative":
        profiles[0] = replace(profiles[0], guild_counts={"A": -1, "B": 1})
    else:
        profiles[0] = replace(profiles[0], condition="different")
    with pytest.raises(ValueError, match=message):
        compute_log_ratio_shifts(profiles, ("A", "B"))


@pytest.mark.parametrize("traits", [("A",), ("A", "A"), ()])
def test_invalid_hypothesis_family_is_rejected(traits: tuple[str, ...]) -> None:
    with pytest.raises(ValueError, match="distinct traits"):
        compute_log_ratio_shifts([], traits)


def test_cli_writes_all_fifteen_hypotheses_missing_states_and_provenance(tmp_path: Path) -> None:
    profile_path = tmp_path / "profiles.tsv"
    experiments = tmp_path / "experiments.tsv"
    dump = tmp_path / "dump.csv"
    columns = [ids[0] for ids in TABLE1_TRAITS.values()]
    # No selected trait changes its ratio; enormous unrelated MediaDive change.
    profile_path.write_text(
        "\t".join(["study_id", "experiment_id", "inc_dec", *columns, "mediadive.medium:X"])
        + "\ns1\te1\tinc\t1\t1\t1\t1\t1\t0\t10000"
        + "\ns1\te1\tdec\t1\t1\t1\t1\t1\t0\t1\n",
        encoding="utf-8",
    )
    experiments.write_text(
        "study_id\texperiment_id\tcondition\tbody_site\thost_species\n"
        "s1\te1\tX\tFeces\tHomo sapiens\n",
        encoding="utf-8",
    )
    dump.write_text("# BugSigDB pinned-test, License: test\n", encoding="utf-8")
    output = tmp_path / "output"
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
                "--body-site",
                "Feces",
            ]
        )
        == 0
    )
    with (output / "trait_log_ratio_stats.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 15
    assert {r["family_size"] for r in rows} == {"15"}
    assert {r["inference_level"] for r in rows} == {"body_site"}
    assert {r["status"] for r in rows} == {"all_tied", "no_positive_pairs"}
    assert all(r["share_ratio_up"] == r["sign_p"] == r["by_q"] == "" for r in rows)
    provenance = json.loads((output / "trait_log_ratio_provenance.json").read_text())
    assert provenance["profile_sha256"] == hashlib.sha256(profile_path.read_bytes()).hexdigest()
    assert provenance["bugsigdb_vintage"] == "pinned-test"
    assert provenance["family_size"] == 15
    assert provenance["correction_method"] == "Benjamini-Yekutieli"
    assert "NA" in (output / "trait_log_ratio_report.md").read_text()

    # A missing count in the source must not turn into an observed zero/exclusion.
    profile_path.write_text(profile_path.read_text().replace("inc\t1\t", "inc\t\t"))
    with pytest.raises(ValueError, match="Missing trait count"):
        main(["--profile-input", str(profile_path), "--experiments-input", str(experiments)])
