"""Study-weighted, population-bound annotation sensitivity regressions (#772)."""

from __future__ import annotations

import json
import math
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any
from unittest.mock import patch

import numpy as np
import pytest

import mtp.presence_depth_sensitivity as sensitivity
from mtp.contribution_inputs import rows
from mtp.paths import REPO_ROOT, file_sha256
from mtp.presence_depth_sensitivity import (
    StudyEvidence,
    build_screen,
    estimate_hc3,
    load_core_evidence,
    make_study_observations,
    reconstruct_profiles,
    run_sensitivity,
)
from mtp.workflow.consumers import CoreBundle, open_core_bundle, verify_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_execution_product
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
        ([-1.0, 1.0], "insufficient_study_support"),
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
    assert profile.depth_primary == pytest.approx(math.log10(11) / 2)
    assert profile.depth_secondary == pytest.approx(math.log1p(6))


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
    assert studies[0].depth_primary == pytest.approx(math.log10(5) / 2)
    assert studies[0].depth_primary != pytest.approx(math.log10(3))
    assert studies[0].depth_secondary == pytest.approx(math.log(5) / 2)


def test_presence_breadth_and_taxon_prevalence_are_distinct_estimands() -> None:
    observed = {("s", "e", "inc"): (1, 1, 1), ("s", "e", "dec"): (1, 1, 0)}
    reported = {("s", "e", "inc"): ["broad"], ("s", "e", "dec"): ["a", "b"]}
    profiles = reconstruct_profiles(
        observed,
        reported,
        {"broad": (1, 1, 1), "a": (1, 0, 0), "b": (0, 1, 0)},
        ("t1", "t2", "t3"),
    )
    studies, _ = make_study_observations(
        profiles, {("s", "e"): {"host_species": "Homo sapiens", "sequencing_type": "16S"}}
    )
    by_representation = {row.representation: row for row in studies}
    assert by_representation["presence"].mean_differences[0] == pytest.approx(-1 / 6)
    assert by_representation["prevalence"].mean_differences[0] == pytest.approx(0.5)


def test_focal_mapping_requires_the_actual_trait_identifiers() -> None:
    from mtp.trait_guilds import TABLE1_TRAITS

    columns = [f"METPO:relation|{terms[0]}" for terms in TABLE1_TRAITS.values()]
    mapped = sensitivity.focal_feature_map(columns)
    assert mapped == dict(zip(TABLE1_TRAITS, columns, strict=True))
    assert sensitivity.focal_feature_map(columns[:1]) == {next(iter(TABLE1_TRAITS)): columns[0]}


def test_focal_rows_keep_primary_corrections_and_exact_trait_identity() -> None:
    from mtp.trait_guilds import TABLE1_TRAITS, column_term_id

    columns = [f"METPO:relation|{terms[0]}" for terms in TABLE1_TRAITS.values()]
    screen = [
        {
            "family_id": "fixture:full",
            "representation": "raw",
            "test": "sign",
            "feature_id": column,
            "feature_ids": [column],
            "status": "tested",
            "p_value": 0.001,
            "q_by": 1.0,
            "q_bh": 0.317,
            "family_size": 317,
        }
        for column in columns
    ]
    for focal in sensitivity.focal_screen(screen, columns):
        assert column_term_id(focal["feature_id"]) in TABLE1_TRAITS[focal["focal_trait"]]
        assert focal["feature_ids"] == [focal["feature_id"]]
        assert focal["q_by"] == 1.0
        assert focal["q_bh"] == 0.317
        assert focal["family_size"] == 317

    # Missing focal classes retain their biological identity and explicit null calls.
    for focal in sensitivity.focal_screen(screen[:1], columns[:1]):
        assert column_term_id(focal["feature_id"]) in TABLE1_TRAITS[focal["focal_trait"]]
        if focal["status"] == "missing_feature":
            assert focal["q_by"] is None
            assert focal["q_bh"] is None


def test_native_sign_tally_mismatch_is_rejected() -> None:
    from mtp.trait_guilds import Profile
    from mtp.trait_screen import run_screen

    columns = ["METPO:1", "METPO:2"]
    studies = [
        StudyEvidence("s", name, 1, (0.5, -0.5), 0.0, 0.0) for name in sensitivity.REPRESENTATIONS
    ]
    native = {
        name: [
            Profile("s", "e", direction, "", "", dict(zip(columns, counts, strict=True)), 4)
            for direction, counts in (("inc", (3, 1)), ("dec", (1, 3)))
        ]
        for name in sensitivity.REPRESENTATIONS
    }
    # Every representation has a consistent nonempty cohort before corruption.
    baseline = build_screen(studies, columns, native)
    assert {row["representation"] for row in baseline} == set(sensitivity.REPRESENTATIONS)
    screened = [run_screen(native[name], columns) for name in sensitivity.REPRESENTATIONS]
    target = next(row for row in screened[0].rows if row.trait == "METPO:1")
    assert (target.studies_up, target.studies_down, target.studies_tied) == (1, 0, 0)
    # Keep sample size fixed but reverse a native sign: a count-only check is insufficient.
    bad = replace(target, studies_up=0, studies_down=1)
    assert bad.studies_up + bad.studies_down + bad.studies_tied == 1
    screened[0] = replace(
        screened[0], rows=[bad if row.trait == target.trait else row for row in screened[0].rows]
    )
    with (
        patch.object(sensitivity, "run_screen", side_effect=screened),
        pytest.raises(ValueError, match="study sign tallies disagree"),
    ):
        build_screen(studies, columns, native)


def test_native_screen_and_composition_are_called_with_reconstructed_inputs() -> None:
    from mtp.trait_composition import compute_log_ratio_shifts
    from mtp.trait_guilds import TABLE1_TRAITS
    from mtp.trait_screen import run_screen

    columns = [terms[0] for terms in TABLE1_TRAITS.values()]
    matrix = {
        "a": (4, 2, 0, 2, 2, 3),
        "b": (1, 0, 3, 2, 4, 3),
        "c": (0, 4, 2, 3, 1, 0),
        "d": (2, 1, 4, 0, 2, 1),
    }
    reported = {}
    metadata = {}
    for study, (inc, dec) in enumerate([("ac", "bd"), ("ab", "cd"), ("ad", "bc")]):
        metadata[(str(study), "e")] = {"host_species": "Homo sapiens", "sequencing_type": "16S"}
        reported[(str(study), "e", "inc")] = list(inc)
        reported[(str(study), "e", "dec")] = list(dec)
    observed = {
        key: tuple(sum(matrix[taxon][i] for taxon in taxa) for i in range(6))
        for key, taxa in reported.items()
    }
    profiles = reconstruct_profiles(observed, reported, matrix, columns)
    studies, ledger = make_study_observations(profiles, metadata)
    native = sensitivity.native_profile_inputs(profiles, metadata, columns, ledger)
    expected = run_screen(native["presence"], columns)
    with patch.object(sensitivity, "run_screen", wraps=run_screen) as called:
        screen = build_screen(studies, columns, native)
        assert called.call_count == 3
    for original in expected.rows:
        row = next(
            row
            for row in screen
            if row["representation"] == "presence"
            and row["test"] == "sign"
            and row["feature_id"] == original.trait
        )
        assert (row["p_value"], row["ci_low"], row["ci_high"]) == (
            original.p_value,
            original.ci_low,
            original.ci_high,
        )
        assert row["native_by_q_nominal"] == original.by_q_nominal
        assert row["native_by_q_effective"] == original.by_q_effective
    with patch.object(
        sensitivity, "compute_log_ratio_shifts", wraps=compute_log_ratio_shifts
    ) as called:
        ratios, per_study = sensitivity.composition_screen(native, columns)
        assert called.call_count == 2
    assert len(ratios) == 30
    assert len(per_study) == 90
    assert all(row["native_available_pair_family_size"] == 15 for row in ratios)
    for representation in ("raw", "presence"):
        mapped = sensitivity.focal_feature_map(columns)
        inputs = [
            replace(
                p, guild_counts={slug: p.guild_counts[column] for slug, column in mapped.items()}
            )
            for p in native[representation]
        ]
        original_ratios = compute_log_ratio_shifts(inputs)
        for summary in original_ratios.summaries:
            row = next(
                row
                for row in ratios
                if row["representation"] == representation
                and row["numerator"] == summary.numerator
                and row["reference"] == summary.reference
            )
            for key, value in asdict(summary).items():
                assert row["native_available_pair_by_q" if key == "by_q" else key] == value


def test_null_empty_and_untestable_candidates_remain_in_full_families() -> None:
    result = estimate_hc3([0.0] * 4)
    assert result["p_value"] is None
    screen = build_screen(
        [], ["METPO:1", "METPO:2"], {name: [] for name in sensitivity.REPRESENTATIONS}
    )
    assert len(screen) == 24
    assert all(row["p_value"] is None for row in screen)


def test_roundoff_ties_have_no_p_values_and_keep_family_slots() -> None:
    from mtp.trait_guilds import Profile

    columns = ["METPO:1", "METPO:2"]
    studies = [StudyEvidence(str(i), "presence", 1, (1e-14, -1e-14), 0.0, 0.0) for i in range(5)]
    native: dict[str, list[Profile]] = {name: [] for name in sensitivity.REPRESENTATIONS}
    native["presence"] = [
        Profile(str(i), "e", direction, "", "", dict(zip(columns, counts, strict=True)), 10**14)
        for i in range(5)
        for direction, counts in (("inc", (1, 10**14 - 1)), ("dec", (0, 10**14)))
    ]
    screen = build_screen(studies, columns, native)
    sign = [row for row in screen if row["representation"] == "presence" and row["test"] == "sign"]
    assert all(row["n_positive"] == 0 and row["n_tied"] == 5 for row in sign)
    assert all(row["p_value"] is None and row["native_by_q_nominal"] is None for row in sign)


def _core(tmp_path: Path, subset: bool, scope: str = "phenotype") -> Path:
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
            "strain_rollup": scope,
            **(
                {"contrast_population": "disease_vs_healthy_tier1_2_subject_controls"}
                if subset
                else {}
            ),
        },
        validator=validate_execution_product,
        repository_root=REPO_ROOT,
    )
    return Path(run.manifest_path)


@pytest.mark.parametrize("subset", [False, True])
@pytest.mark.parametrize("scope", ["off", "phenotype"])
def test_verified_core_run_records_actual_family_and_population(
    tmp_path: Path, subset: bool, scope: str
) -> None:
    from mtp.trait_guilds import TABLE1_TRAITS, column_term_id

    manifest = _core(tmp_path, subset, scope)
    output = tmp_path / "sensitivity"
    binding = run_sensitivity(manifest, output)
    record = verify_consumer_record(binding)
    assert "reported_taxa" in record.product_ids
    provenance = json.loads((output / "provenance.json").read_text())
    bases = json.loads((output / "analysis_bases.json").read_text())
    assert provenance["retained_features"] == bases[0]["retained_feature_count"] == 2
    assert provenance["strain_rollup"] == scope
    assert all(basis["rollup_policy"] == f"ncbi_taxonomy_closure_rollup_{scope}" for basis in bases)
    assert provenance["direction_basis"] == (
        "disease_vs_comparison" if subset else "group_1_vs_group_0"
    )
    audit = list(rows(output / "profile_audit.tsv"))
    if subset:
        assert {row["experiment_id"] for row in audit} == {"bsdb:1-1", "bsdb:1-2"}
    full_screen = list(rows(output / "trait_screen.tsv"))
    assert len(full_screen) == 24
    for row in rows(output / "focal_traits.tsv"):
        assert column_term_id(row["feature_id"]) in TABLE1_TRAITS[row["focal_trait"]]
        assert "q_by_six_trait" in row
        assert "q_bh_six_trait" in row
        if row["status"] == "missing_feature":
            assert row["p_value"] == ""
            assert row["q_by"] == row["q_bh"] == ""
            continue
        original = next(
            r
            for r in full_screen
            if r["family_id"] == row["full_family_id"] and r["feature_id"] == row["feature_id"]
        )
        assert row["q_by"] == original["q_by"]
        assert row["q_bh"] == original["q_bh"]
        assert row["family_size"] == original["family_size"]
        assert row["family_size_six_trait"] == "6"
    ora_record = verify_consumer_record(output / "ora_core_consumption.json")
    assert "ora_profiles" in ora_record.product_ids
    assert "count_profiles" not in ora_record.product_ids
    families = json.loads((output / "families.json").read_text())
    assert families
    assert {family["strain_rollup"] for family in families} == {scope}
    assert all(
        row["expected_family_size"] == 6
        for row in families
        if row["family_id"].endswith(":focal_six")
    )
    assert all(
        row["expected_family_size"] == 15
        for row in families
        if row["family_id"].endswith(":composition")
    )
    ratios = list(rows(output / "composition_presence.tsv"))
    assert len(ratios) == 30
    unavailable = [row for row in ratios if row["status"] == "missing_feature"]
    assert len(unavailable) == 28
    assert all(row["p_value"] == "" for row in unavailable)
    assert {
        (row["numerator"], row["reference"]) for row in ratios if row["status"] != "missing_feature"
    } == {("facultatively_aerobic", "coccus_shaped")}
    assert bases[-1]["denominator_policy"] == "positive_reference_trait_count_no_pseudocount"
    assert len(provenance["code_commit"]) == 40
    assert isinstance(provenance["code_dirty"], bool)
    assert {
        "src/mtp/trait_screen.py",
        "src/mtp/workflow/consumers.py",
        "src/mtp/guild_maps.py",
        "src/mtp/workflow/verification.py",
        "src/mtp/study_trait_profiles.py",
        "methods/multiple_testing_policy_v1.yaml",
    }.issubset(provenance["implementation_sources_sha256"])
    assert all(
        file_sha256(REPO_ROOT / path) == digest
        for path, digest in provenance["implementation_sources_sha256"].items()
    )
    eligibility = list(rows(output / "eligibility.tsv"))
    if not subset:
        assert (
            next(row for row in eligibility if row["experiment_id"] == "bsdb:2-1")["status"]
            == "non_human"
        )
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


def test_complete_family_keeps_untestable_slots_in_authoritative_correction(tmp_path: Path) -> None:
    from mtp.datamodel.multiple_testing_pydantic import (
        CorrectionRun,
        HypothesisFamily,
        RawHypothesis,
    )
    from mtp.multiple_testing import verify_correction
    from mtp.presence_depth_families import declare_and_correct

    manifest = _core(tmp_path, False)
    output = tmp_path / "sensitivity"
    run_sensitivity(manifest, output)
    bundle = open_core_bundle(manifest)
    columns, _, _ = load_core_evidence(bundle)
    bases = json.loads((output / "analysis_bases.json").read_text())
    basis_map = {basis["id"].removeprefix(f"{sensitivity.PROTOCOL}:"): basis for basis in bases}
    ledger = list(rows(output / "eligibility.tsv"))
    records: list[sensitivity.Result] = [
        {
            "family_id": "fixture:full",
            "representation": "raw",
            "test": "sign",
            "hypothesis_id": column,
            "feature_ids": [column],
            "p_value": 0.01 if index == 0 else None,
            "status": "tested" if index == 0 else "all_tied",
        }
        for index, column in enumerate(columns)
    ]
    families, corrections = declare_and_correct(records, basis_map, bundle, ledger, "off")
    assert families[0]["expected_family_size"] == 2
    assert records[0]["q_by"] == pytest.approx(0.03)
    assert records[0]["q_bh"] == pytest.approx(0.02)
    assert records[1]["q_by"] is None
    focal = sensitivity.focal_screen(
        [{**row, "feature_id": row["hypothesis_id"]} for row in records],
        columns,
    )
    focal_families, _ = declare_and_correct(focal, basis_map, bundle, ledger, "off")
    assert focal_families[0]["expected_family_size"] == 6
    tested_focal = [row for row in focal if row["p_value"] is not None]
    assert len(tested_focal) == 1
    assert tested_focal[0]["q_by"] == records[0]["q_by"]
    assert tested_focal[0]["q_bh"] == records[0]["q_bh"]
    assert tested_focal[0]["q_by_six_trait"] == pytest.approx(0.147)
    assert tested_focal[0]["q_bh_six_trait"] == pytest.approx(0.06)
    assert tested_focal[0]["family_size"] == 2
    assert tested_focal[0]["family_size_six_trait"] == 6
    corrected = CorrectionRun.model_validate(corrections[0])
    family = HypothesisFamily.model_validate(families[0])
    inputs = [
        RawHypothesis(
            hypothesis_id=row["hypothesis_id"],
            p_value=row["p_value"],
            status="tested" if row["p_value"] is not None else "untestable",
            reason=None if row["p_value"] is not None else row["status"],
        )
        for row in records
    ]
    verify_correction(corrected, family=family, hypotheses=inputs)
    corrected.results[0].q_by = 0.01
    with pytest.raises(ValueError, match="verification failed"):
        verify_correction(corrected, family=family, hypotheses=inputs)


def test_binary_ora_product_is_read_and_reconciled(tmp_path: Path) -> None:
    bundle = open_core_bundle(_core(tmp_path, True))
    _, profiles, _ = load_core_evidence(bundle)
    assert all(profile.ora_taxa_count == profile.matched_taxa for profile in profiles)
    path = bundle.path("ora_profiles")
    content = path.read_text().splitlines()
    header = content[0].split("\t")
    first = content[1].split("\t")
    first[header.index("taxa_count")] = str(int(first[header.index("taxa_count")]) + 1)
    content[1] = "\t".join(first)
    path.write_text("\n".join(content) + "\n")
    # Direct loader rechecks the consumed values, independently of opening the core again.
    with pytest.raises(ValueError, match="taxa_count fail source reconstruction"):
        load_core_evidence(bundle)


@pytest.mark.parametrize("change", ["missing", "extra"])
def test_binary_ora_population_requires_exact_equality(tmp_path: Path, change: str) -> None:
    bundle = open_core_bundle(_core(tmp_path, True))
    path = bundle.path("ora_profiles")
    lines = path.read_text().splitlines()
    if change == "missing":
        lines.pop(1)
    else:
        extra = lines[1].split("\t")
        extra[0] = "bsdb:outside-selected-population"
        lines.append("\t".join(extra))
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="population differs"):
        load_core_evidence(bundle)


def test_full_rollup_is_explicitly_unsupported_before_analysis(tmp_path: Path) -> None:
    manifest = _core(tmp_path, True, "full")
    output = tmp_path / "unsupported"
    with (
        patch.object(sensitivity, "load_core_evidence") as load,
        pytest.raises(ValueError, match=r"(?i)full.*not supported.*multiplicity"),
    ):
        run_sensitivity(manifest, output)
    load.assert_not_called()
    assert not output.exists()


def test_source_inventory_tracks_additionally_loaded_mtp_module() -> None:
    import mtp.coupling
    from mtp.multiple_testing import POLICY_PATH

    sources = sensitivity.implementation_sources()
    assert sources["src/mtp/coupling.py"] == file_sha256(Path(mtp.coupling.__file__))
    assert sources["src/mtp/workflow/verification.py"] == file_sha256(
        REPO_ROOT / "src/mtp/workflow/verification.py"
    )
    assert sources[POLICY_PATH.relative_to(REPO_ROOT).as_posix()] == file_sha256(POLICY_PATH)


@pytest.mark.parametrize("change", ["code", "policy", "lazy_code", "missing_loaded", "core_open"])
def test_changed_execution_identity_is_rejected_before_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, change: str
) -> None:
    manifest = _core(tmp_path, True)
    output = tmp_path / "changed-identity"
    hashes = sensitivity.implementation_sources()
    target = {
        "code": "src/mtp/presence_depth_sensitivity.py",
        "policy": "methods/multiple_testing_policy_v1.yaml",
        "lazy_code": "src/mtp/coupling.py",
        "missing_loaded": "src/mtp/presence_depth_sensitivity.py",
        "core_open": "src/mtp/presence_depth_sensitivity.py",
    }[change]
    if change == "lazy_code":
        hashes.pop(target, None)
    original = sensitivity.build_screen
    original_open = open_core_bundle

    def changed_open(*args: Any, **kwargs: Any) -> CoreBundle:
        bundle = original_open(*args, **kwargs)
        hashes[target] = "f" * 64
        return bundle

    def changed_analysis(*args: Any, **kwargs: Any) -> list[sensitivity.Result]:
        result = original(*args, **kwargs)
        if change == "missing_loaded":
            hashes.pop(target)
        else:
            # A lazy module's changed bytes precede its first inventory appearance.
            hashes[target] = "f" * 64
        return result

    monkeypatch.setattr(sensitivity, "implementation_sources", lambda: dict(hashes))
    if change == "core_open":
        monkeypatch.setattr(sensitivity, "open_core_bundle", changed_open)
    else:
        monkeypatch.setattr(sensitivity, "build_screen", changed_analysis)
    # Before the fix, even real final consumer verification accepts these outputs.
    with pytest.raises(ValueError, match="implementation or policy changed during analysis"):
        verify_consumer_record(run_sensitivity(manifest, output))
    assert not output.exists()


def test_unchanged_lazy_source_is_bound_to_its_preanalysis_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest = _core(tmp_path, True)
    output = tmp_path / "lazy-identity"
    hashes = sensitivity.implementation_sources()
    target = "src/mtp/coupling.py"
    hashes.pop(target, None)
    original = sensitivity.build_screen
    expected = file_sha256(REPO_ROOT / target)

    def lazy_analysis(*args: Any, **kwargs: Any) -> list[sensitivity.Result]:
        result = original(*args, **kwargs)
        hashes[target] = expected
        return result

    monkeypatch.setattr(sensitivity, "implementation_sources", lambda: dict(hashes))
    monkeypatch.setattr(sensitivity, "build_screen", lazy_analysis)
    verify_consumer_record(run_sensitivity(manifest, output))
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["implementation_sources_sha256"][target] == expected
