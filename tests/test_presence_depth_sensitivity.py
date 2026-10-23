"""Study-weighted, population-bound annotation sensitivity regressions (#772)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from mtp.contribution_inputs import rows
from mtp.paths import REPO_ROOT
from mtp.presence_depth_sensitivity import (
    StudyEvidence,
    build_screen,
    estimate_hc3,
    load_core_evidence,
    make_study_observations,
    reconstruct_profiles,
    run_sensitivity,
)
from mtp.workflow.consumers import open_core_bundle, verify_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow


def test_hc3_matches_analytic_orthogonal_design() -> None:
    x = [-2.0, -1.0, 0.0, 1.0, 2.0]
    y = [0.2 + 0.1 * v + e for v, e in zip(x, [0.02, -0.04, 0.04, -0.04, 0.02], strict=True)]
    result = estimate_hc3(y, x)
    expected_se = math.sqrt((3 * 0.05**2 + 2 * (0.04 / 0.7) ** 2) / 25)
    assert result["effect"] == pytest.approx(0.2)
    assert result["slope"] == pytest.approx(0.1)
    assert result["standard_error"] == pytest.approx(expected_se)
    assert result["degrees_of_freedom"] == 3
    assert result["status"] == "tested"


def test_depth_only_artifact_and_real_shift_are_distinguished() -> None:
    x = [-1.0, 0.0, 1.0, 2.0, 3.0]
    noise = [0.02, -0.04, 0.04, -0.04, 0.02]
    artifact = [0.1 * value + error for value, error in zip(x, noise, strict=True)]
    assert estimate_hc3(artifact)["effect"] == pytest.approx(0.1)
    adjusted = estimate_hc3(artifact, x)
    assert adjusted["effect"] == pytest.approx(0, abs=1e-14)
    assert adjusted["p_value"] == pytest.approx(1.0)
    shifted = estimate_hc3([v + 0.2 for v in artifact], x)
    assert shifted["effect"] == pytest.approx(0.2)
    assert float(shifted["p_value"]) < 0.05


@pytest.mark.parametrize(
    ("x", "status"),
    [
        ([1.0, 2.0, 3.0, 4.0], "zero_depth_outside_support"),
        ([1.0, 1.0, 1.0, 1.0], "constant_nonzero_depth"),
        ([-1.0, 1.0], "insufficient_studies"),
    ],
)
def test_unsupported_adjustments_are_not_significant(x: list[float], status: str) -> None:
    result = estimate_hc3([0.1] * len(x), x)
    assert result["status"] == status
    assert result["p_value"] is None


def test_reconstruction_deduplicates_taxa_and_audits_unmatched_zero_rows() -> None:
    key = ("s", "e", "inc")
    records = reconstruct_profiles(
        {key: (10, 2)},
        {key: ["a", "a", "zero", "unmatched"]},
        {"a": (10, 2), "zero": (0, 0)},
        ("trait1", "trait2"),
    )
    profile = records[0]
    assert profile.raw == (10, 2)
    assert profile.presence == (1, 1)
    assert profile.reported_taxa == 3
    assert profile.matched_taxa == 2
    assert profile.unmatched_taxa == ("unmatched",)
    assert profile.unannotated_taxa == 1
    assert profile.depth_mean_total == 6
    assert profile.depth_primary == pytest.approx(math.log1p(6))
    assert profile.depth_secondary == pytest.approx(math.log10(11) / 2)


def test_inconsistent_reconstruction_is_an_error() -> None:
    key = ("s", "e", "inc")
    with pytest.raises(ValueError, match="reconstruct"):
        reconstruct_profiles({key: (11, 2)}, {key: ["a"]}, {"a": (10, 2)}, ("t1", "t2"))


def test_one_study_with_many_experiments_does_not_gain_weight() -> None:
    observed = {}
    reported = {}
    metadata = {}
    matrix = {"up": (3, 1), "down": (1, 3)}
    for study, experiments in [("many", 10), ("one", 1)]:
        for index in range(experiments):
            experiment = f"e{index}"
            metadata[(study, experiment)] = {
                "host_species": "Homo sapiens",
                "sequencing_type": "16S",
            }
            for direction in ("inc", "dec"):
                taxon = "up" if (study == "many") == (direction == "inc") else "down"
                key = (study, experiment, direction)
                observed[key] = matrix[taxon]
                reported[key] = [taxon]
    profiles = reconstruct_profiles(observed, reported, matrix, ("t1", "t2"))
    studies, ledger = make_study_observations(profiles, metadata)
    raw = [row for row in studies if row.representation == "raw"]
    assert len(raw) == 2
    assert sorted(row.experiments for row in raw) == [1, 10]
    assert np.mean([row.mean_differences[0] for row in raw]) == 0
    assert len(ledger) == 11
    assert all(row["status"] == "included" for row in ledger)
    assert all(
        row.mean_differences == (0.0, 0.0) for row in studies if row.representation == "presence"
    )


def test_zero_total_pair_is_excluded_not_converted_to_zero_effect() -> None:
    observed = {("s", "e", "inc"): (2, 1), ("s", "e", "dec"): (0, 0)}
    reported = {("s", "e", "inc"): ["a"], ("s", "e", "dec"): ["z"]}
    profiles = reconstruct_profiles(observed, reported, {"a": (2, 1), "z": (0, 0)}, ("t1", "t2"))
    studies, ledger = make_study_observations(
        profiles,
        {("s", "e"): {"host_species": "Homo sapiens", "sequencing_type": "16S"}},
    )
    assert studies == []
    assert ledger[0]["status"] == "zero_retained_total"


def test_depth_transform_is_applied_before_pairing_and_study_averaging() -> None:
    observed = {
        ("s", "e1", "inc"): (9, 0),
        ("s", "e1", "dec"): (0, 1),
        ("s", "e2", "inc"): (1, 0),
        ("s", "e2", "dec"): (0, 1),
    }
    matrix = {"a": (9, 0), "b": (0, 1), "c": (1, 0)}
    reported = {
        key: ["a" if value == (9, 0) else "c" if value == (1, 0) else "b"]
        for key, value in observed.items()
    }
    profiles = reconstruct_profiles(observed, reported, matrix, ("t1", "t2"))
    metadata = {
        ("s", e): {"host_species": "Homo sapiens", "sequencing_type": "WMS"} for e in ("e1", "e2")
    }
    studies, _ = make_study_observations(profiles, metadata)
    assert studies[0].depth_primary == pytest.approx(math.log(5) / 2)
    assert studies[0].depth_primary != pytest.approx(math.log(3))


def test_null_empty_and_untestable_candidates_remain_in_full_families() -> None:
    result = estimate_hc3([0.0] * 4)
    assert result["p_value"] is None
    screen, families = build_screen([], ["METPO:1", "METPO:2"])
    assert len(screen) == 16
    assert all(row["p_value"] is None and row["q_by"] is None for row in screen)
    assert all(row["family_size"] == 2 and row["tested"] == 0 for row in families)


def test_roundoff_ties_have_no_p_values_and_keep_family_slots() -> None:
    studies = [StudyEvidence(str(i), "presence", 1, (1e-14, 0.0), 0.0, 0.0) for i in range(5)]
    screen, families = build_screen(studies, ["METPO:1", "METPO:2"])
    sign = [row for row in screen if row["representation"] == "presence" and row["test"] == "sign"]
    assert all(row["n_positive"] == 0 and row["n_tied"] == 5 for row in sign)
    assert all(row["p_value"] is None and row["q_by"] is None for row in sign)
    family = next(row for row in families if row["id"].endswith(":presence:sign"))
    assert family["tested"] == 0
    assert family["family_size"] == 2


def _core(tmp_path: Path, subset: bool) -> Path:
    source = REPO_ROOT / "tests/fixtures/workflow"
    edges, nodes = tmp_path / "edges.tsv", tmp_path / "nodes.tsv"
    edges.write_text(
        (source / "merged-kg_edges.tsv").read_text()
        + "depth200\tNCBITaxon:200\tbiolink:has_phenotype\tMETPO:1000668\tfixture:depth\n"
        + "depth300\tNCBITaxon:300\tbiolink:has_phenotype\tMETPO:1000668\tfixture:depth\n"
    )
    nodes.write_text(
        (source / "merged-kg_nodes.tsv").read_text()
        + "METPO:1000668\tbiolink:PhenotypicQuality\tFixture coccus\n"
    )
    run = execute_workflow(
        REPO_ROOT / "workflows" / ("bugsigdb_core_subset.yaml" if subset else "bugsigdb_core.yaml"),
        {
            "raw_dump": source / "contrast_population.csv",
            "kg_edges": edges,
            "kg_nodes": nodes,
        },
        tmp_path / ("subset" if subset else "full"),
        OPERATIONS,
        parameter_overrides={
            "source_mode": "fixture",
            **(
                {"contrast_population": "disease_vs_healthy_tier1_2_subject_controls"}
                if subset
                else {}
            ),
        },
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    return Path(run.manifest_path)


@pytest.mark.parametrize("subset", [False, True])
def test_verified_core_run_records_actual_family_and_population(
    tmp_path: Path, subset: bool
) -> None:
    manifest = _core(tmp_path, subset)
    output = tmp_path / "sensitivity"
    binding = run_sensitivity(manifest, output)
    record = verify_consumer_record(binding)
    assert "reported_taxa" in record.product_ids
    provenance = json.loads((output / "provenance.json").read_text())
    bases = json.loads((output / "analysis_bases.json").read_text())
    assert provenance["retained_features"] == bases[0]["retained_feature_count"] == 2
    assert provenance["direction_basis"] == (
        "disease_vs_comparison" if subset else "group_1_vs_group_0"
    )
    audit = list(rows(output / "profile_audit.tsv"))
    if subset:
        assert {row["experiment_id"] for row in audit} == {"bsdb:1-1", "bsdb:1-2"}
    full_screen = list(rows(output / "trait_screen.tsv"))
    assert len(full_screen) == 16
    for row in rows(output / "focal_traits.tsv"):
        original = next(
            r
            for r in full_screen
            if r["family_id"] == row["family_id"] and r["feature_id"] == row["feature_id"]
        )
        assert row["q_by"] == original["q_by"]
    with pytest.raises(FileExistsError):
        run_sensitivity(manifest, output)
    (output / "trait_screen.tsv").write_text("tampered\n")
    with pytest.raises(ValueError, match="Consumer output changed"):
        verify_consumer_record(binding)


def test_reversed_subset_is_consumed_once_with_taxa_aligned(tmp_path: Path) -> None:
    full = open_core_bundle(_core(tmp_path, False))
    subset = open_core_bundle(_core(tmp_path, True))
    _, full_rows, _ = load_core_evidence(full)
    _, subset_rows, _ = load_core_evidence(subset)
    full_map = {row.key: row for row in full_rows}
    subset_map = {row.key: row for row in subset_rows}
    key = next(key for key in subset_map if key[1:] == ("bsdb:1-2", "inc"))
    opposite = (*key[:2], "dec")
    assert subset_map[key].raw == full_map[opposite].raw
    assert subset_map[key].presence == full_map[opposite].presence
    assert subset_map[key].raw != full_map[key].raw
    with pytest.raises(ValueError, match="overlaps"):
        run_sensitivity(subset.manifest, subset.manifest.parent / "unsafe")
