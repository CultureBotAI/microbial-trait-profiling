"""Native kernels preserve the pinned tests, denominator, and missing results."""

import json
import math
from collections import Counter
from pathlib import Path
from statistics import mean, pstdev, stdev

import numpy as np
import pytest
from scipy.stats import hypergeom, ttest_1samp

from bsdbng.contribution_inputs import Base
from bsdbng.native_taxon_influence import (
    export_native,
    fold_pairs,
    hypergeom_probabilities,
    ora_cells,
    reconstruct,
    run,
    transformed,
)
from bsdbng.taxon_contributions import write_table
from bsdbng.taxon_influence import adjust, comparison

POLICY = {
    "minimum_profile_support": 1,
    "maximum_profile_support_fraction": 1.0,
    "minimum_pairs": 3,
    "transform_multiplier": 1000,
}


@pytest.mark.parametrize("tail", ["enrichment", "depletion"])
def test_parallel_native_probabilities_are_identical(tail: str) -> None:
    triples = np.asarray([(i % 9, 25 + i % 7, 10 + i % 3) for i in range(5000)], dtype=np.int64)
    serial = hypergeom_probabilities(triples, 60, tail, workers=1)
    parallel = hypergeom_probabilities(triples, 60, tail, workers=4)
    assert np.array_equal(serial, parallel)


def corpus() -> Base:
    profiles = {}
    for i, (inc, dec) in enumerate(((8, 2), (7, 3), (1, 4), (0, 0))):
        for direction, value in (("inc", inc), ("dec", dec)):
            profiles[("s1" if i < 2 else f"s{i}", f"e{i}", direction)] = {
                "A": value,
                "B": 10 - value if i < 3 else 0,
                "C": 0,
            }
    return Base(["A", "B", "C"], {"a": "A", "c": "C"}, profiles, {}, {}, {}, {}, {}, [])


def test_native_composition_matches_log_zscore_and_t_test() -> None:
    base = corpus()
    records, pairs, scales = transformed(base, base.profiles, list(base.profiles), POLICY)
    logs = [
        math.log1p(1000 * p["A"] / sum(p.values())) if sum(p.values()) else 0.0
        for p in base.profiles.values()
    ]
    sigma = pstdev(logs)
    expected = [(logs[i] - logs[i + 1]) / sigma for i in range(0, len(logs), 2)]
    assert scales["A"] == pytest.approx(sigma)
    assert list(pairs["A"].values()) == pytest.approx(expected)
    assert records["A"]["effect"] == pytest.approx(mean(expected))
    assert records["A"]["p"] == pytest.approx(ttest_1samp(expected, 0).pvalue)
    assert records["A"]["pairs"] == 4
    assert records["A"]["zero_total_profiles"] == 2
    assert records["A"]["family_tests"] == 2
    assert records["C"]["status"] == "not_selected"
    assert records["C"]["p"] == records["C"]["q"] == ""
    shared = comparison(pairs["A"], pairs["A"])
    assert shared["before_effect"] != pytest.approx(records["A"]["effect"])


def test_constant_paired_effect_is_not_a_significant_native_result() -> None:
    base = corpus()
    for key in base.profiles:
        base.profiles[key] = {"A": 8 if key[2] == "inc" else 2, "B": 1, "C": 0}
    records, _, _ = transformed(base, base.profiles, list(base.profiles), POLICY)
    assert records["A"]["selected"]
    assert records["A"]["status"] == "untestable"
    assert records["A"]["p"] == ""
    assert records["A"]["q"] == 1
    assert not records["A"]["significant_0.05"]


def test_ora_global_families_include_non_display_traits_and_depletion_zeros() -> None:
    base = corpus()
    background = Counter({t: sum(p[t] for p in base.profiles.values()) for t in base.columns})
    records, family = ora_cells(base, base.profiles, list(base.profiles), background)
    assert family["enrichment_tests"] == 12
    assert family["depletion_tests"] == 18
    observed = [p[t] for p in base.profiles.values() if sum(p.values()) for t in base.columns]
    successes = [
        background[t] for p in base.profiles.values() if sum(p.values()) for t in base.columns
    ]
    pvalues = hypergeom.cdf(observed, 60, successes, 10).tolist()
    expected_q = adjust(pvalues)
    first = records[0]
    assert first["depletion_p"] == pytest.approx(pvalues[0])
    assert first["depletion_q"] == pytest.approx(expected_q[0])
    assert first["enrichment_p"] == pytest.approx(hypergeom.sf(7, 60, 25, 10))
    assert first["null_count_low"] == hypergeom.ppf(0.025, 60, 25, 10)
    assert first["null_count_high"] == hypergeom.ppf(0.975, 60, 25, 10)
    assert all(r["trait_key"] in {"A", "C"} for r in records)
    zero = next(r for r in records if r["row_total"] == 0)
    assert zero["status"] == "zero_total"
    assert zero["depletion_p"] == zero["fold"] == ""
    missing = next(r for r in records if r["trait_key"] == "C" and r["row_total"])
    assert missing["status"] == "zero_background_trait"
    assert missing["depletion_p"] == missing["depletion_q"] == 1
    assert missing["fold"] == ""
    assert not missing["depletion_significant_0.05"]
    assert len(fold_pairs(records, "A")) == 3
    assert fold_pairs(records, "C") == {}


def test_fixed_and_recomputed_backgrounds_answer_different_questions() -> None:
    base = corpus()
    remaining = {k: {**v, "A": 0} for k, v in base.profiles.items()}
    fixed = Counter({t: sum(p[t] for p in base.profiles.values()) for t in base.columns})
    recomputed = Counter({t: sum(p[t] for p in remaining.values()) for t in base.columns})
    first, _ = ora_cells(base, remaining, list(remaining), fixed)
    second, _ = ora_cells(base, remaining, list(remaining), recomputed)
    assert first[0]["fold"] == 0
    assert second[0]["fold"] == ""
    assert first[0]["depletion_p"] < second[0]["depletion_p"] == 1
    assert first[0]["enrichment_p"] == second[0]["enrichment_p"] == ""


def test_empty_scope_and_all_removed_corpus_remain_explicit() -> None:
    base = corpus()
    records, pairs, _ = transformed(base, base.profiles, [], POLICY)
    assert all(r["p"] == "" and r["pairs"] == 0 for r in records.values())
    assert all(v == {} for v in pairs.values())
    cells, family = ora_cells(base, base.profiles, [], Counter())
    assert cells == []
    assert family["depletion_tests"] == 0
    remaining = {k: dict.fromkeys(base.columns, 0) for k in base.profiles}
    cells, family = ora_cells(base, remaining, list(remaining), Counter())
    assert family["enrichment_tests"] == family["depletion_tests"] == 0
    assert all(r["status"] == "zero_total" for r in cells)


def test_ora_rejects_a_noncontaining_background() -> None:
    base = corpus()
    with pytest.raises(ValueError, match="finite population"):
        ora_cells(base, base.profiles, list(base.profiles), Counter(A=1, B=1))


def test_sparse_reconstruction_verifies_full_denominator_and_ledger(tmp_path: Path) -> None:
    base = corpus()
    base.reported = {k: {"taxon": {"signature"}} for k in base.profiles}
    write_table(
        tmp_path / "removed_source_counts.tsv.gz",
        ["reported_taxon", "trait_key", "removed_count"],
        [],
    )
    entries = [
        {
            "study_id": k[0],
            "experiment_id": k[1],
            "direction": k[2],
            "trait_key": t,
            "original_count": p[t],
            "remaining_count": p[t],
            "original_total": sum(p.values()),
            "remaining_total": sum(p.values()),
        }
        for k, p in base.profiles.items()
        for t in base.traits.values()
    ]
    write_table(tmp_path / "profile_counts.tsv.gz", list(entries[0]), entries)
    assert reconstruct(base, tmp_path) == base.profiles
    entries[0]["remaining_total"] = 999
    write_table(tmp_path / "profile_counts.tsv.gz", list(entries[0]), entries)
    with pytest.raises(ValueError, match="reconstruction mismatch"):
        reconstruct(base, tmp_path)


def test_native_run_refuses_changed_audit_inputs(tmp_path: Path) -> None:
    (tmp_path / "manifest.json").write_text(
        json.dumps(
            {
                "completed": True,
                "mismatches": 0,
                "input_sha256": {"study_trait_profiles_wide.tsv": "old"},
            }
        )
    )
    (tmp_path / "study_trait_profiles_wide.tsv").write_text("changed\n")
    with pytest.raises(ValueError, match="Audit input changed"):
        run(tmp_path, tmp_path, tmp_path / "missing-protocol.json")


def test_native_mean_interval_targets_experiments_not_sign_share() -> None:
    base = corpus()
    records, pairs, _ = transformed(base, base.profiles, list(base.profiles), POLICY)
    from bsdbng.taxon_influence import critical_t

    values = list(pairs["A"].values())
    radius = critical_t(3) * stdev(values) / 2
    assert records["A"]["mean_ci_low"] == pytest.approx(mean(values) - radius)
    assert records["A"]["mean_ci_high"] == pytest.approx(mean(values) + radius)
    assert 0 <= records["A"]["share_ci_low"] <= records["A"]["share_ci_high"] <= 1


def test_native_exports_end_to_end_with_an_empty_stratum(tmp_path: Path) -> None:
    base = corpus()
    base.reported = {k: {"taxon": {"signature"}} for k in base.profiles}
    base.metadata = {k[:2]: {"condition": "c", "body_site": "site"} for k in base.profiles}
    directory = tmp_path / "baseline"
    directory.mkdir()
    write_table(
        directory / "removed_source_counts.tsv.gz",
        ["reported_taxon", "trait_key", "removed_count"],
        [],
    )
    entries = [
        {
            "study_id": k[0],
            "experiment_id": k[1],
            "direction": k[2],
            "trait_key": t,
            "original_count": p[t],
            "remaining_count": p[t],
            "original_total": sum(p.values()),
            "remaining_total": sum(p.values()),
        }
        for k, p in base.profiles.items()
        for t in base.traits.values()
    ]
    write_table(directory / "profile_counts.tsv.gz", list(entries[0]), entries)
    write_table(tmp_path / "scenarios.tsv", ["scenario", "status"], [])
    audit = {
        "protocol": {
            "scopes": [
                {"id": "pooled", "condition": "*", "body_site": "*"},
                {"id": "empty", "condition": "missing", "body_site": "*"},
            ]
        }
    }
    result = export_native(
        base,
        audit,
        {"composition": POLICY, "ora": {"backgrounds": ["fixed_baseline", "recomputed_remaining"]}},
        tmp_path,
    )
    assert result["ora_families"] == 4
    assert result["composition_rows"] == 6
    assert result["composition_comparisons"] == 8
    assert result["ora_comparisons"] == 8
    assert "Native lineage" in (tmp_path / "native_summary.md").read_text()
    from bsdbng.contribution_inputs import rows

    comparisons = list(rows(tmp_path / "native_composition_comparison.tsv"))
    assert {r["shared_correction"] for r in comparisons} == {"BH2_directional_sign_test"}
    assert all(r["before_p"] == "" for r in comparisons if r["scope"] == "empty")
    zero_reference = {
        r["reference"]: r for r in comparisons if r["scope"] == "pooled" and r["trait_key"] == "C"
    }
    refitted = zero_reference["refitted"]
    assert refitted["before_pairs"] == refitted["after_pairs"] == "4"
    assert refitted["lost_pairs"] == refitted["gained_pairs"] == "0"
    assert float(refitted["before_effect"]) == float(refitted["after_effect"]) == 0
    assert refitted["before_p"] == refitted["after_p"] == ""
    fixed = zero_reference["fixed_baseline_sd"]
    assert fixed["before_pairs"] == fixed["after_pairs"] == "0"
    assert fixed["before_effect"] == fixed["after_effect"] == ""
