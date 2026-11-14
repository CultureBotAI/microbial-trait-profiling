"""Publication must preserve complete six-trait tables and reject changed runs."""

import json
from pathlib import Path

import pytest

from mtp.contribution_inputs import digest, rows
from mtp.native_influence_publication import interpretation_handoff, publish
from mtp.native_taxon_influence import compressed_rows
from mtp.taxon_contributions import write_table


@pytest.fixture
def run_directory(tmp_path: Path) -> Path:
    source = tmp_path / "run"
    source.mkdir()
    (source / "manifest.json").write_text(
        json.dumps({"completed": True, "mismatches": 0, "traits": {"a": "A"}})
    )
    output_hashes = {}
    for filename in (
        "native_composition.tsv.gz",
        "native_ora_fold_comparison.tsv",
        "native_ora_families.tsv",
    ):
        write_table(
            source / filename,
            ["trait_key", "value"],
            [{"trait_key": "A", "value": 1}, {"trait_key": "OTHER", "value": 2}],
        )
        output_hashes[filename] = digest(source / filename)
    write_table(
        source / "native_composition_comparison.tsv",
        [
            "scenario",
            "scope",
            "trait_key",
            "reference",
            "native_before_weighting",
            "native_before_unit",
            "native_before_family_tests",
            "native_before_q",
            "native_after_q",
        ],
        [
            {
                "scenario": "baseline",
                "scope": "pooled",
                "trait_key": "METPO:1000683",
                "reference": "refitted",
                "native_before_weighting": "equal_experiments",
                "native_before_unit": "within_scope_refitted_zscore_contrast",
                "native_before_family_tests": 315,
                "native_before_q": 0.00003,
                "native_after_q": 0.00003,
            },
            {
                "scenario": "baseline",
                "scope": "pooled",
                "trait_key": "METPO:2000103|pathways:arsenate_reduction",
                "reference": "refitted",
                "native_before_weighting": "equal_experiments",
                "native_before_unit": "within_scope_refitted_zscore_contrast",
                "native_before_family_tests": 315,
                "native_before_q": 0.00003,
                "native_after_q": 0.00003,
            },
            {
                "scenario": "exclude_genus_NCBITaxon:561",
                "scope": "pooled",
                "trait_key": "METPO:2000103|pathways:arsenate_reduction",
                "reference": "refitted",
                "native_before_weighting": "equal_experiments",
                "native_before_unit": "within_scope_refitted_zscore_contrast",
                "native_before_family_tests": 315,
                "native_before_q": 0.00003,
                "native_after_q": 0.67,
            },
        ],
    )
    output_hashes["native_composition_comparison.tsv"] = digest(
        source / "native_composition_comparison.tsv"
    )
    write_table(
        source / "effects.tsv.gz",
        [
            "scenario",
            "scope",
            "representation",
            "trait_key",
            "denominator_trait",
            "unit",
            "correction",
            "before_effect",
            "before_mean_ci_low",
            "before_mean_ci_high",
            "after_effect",
            "after_mean_ci_low",
            "after_mean_ci_high",
            "matched_change_effect",
            "matched_change_mean_ci_low",
            "matched_change_mean_ci_high",
            "before_q",
            "after_q",
        ],
        [
            {
                "scenario": "baseline",
                "scope": "pooled",
                "representation": "shared_fraction",
                "trait_key": "METPO:1000683",
                "denominator_trait": "",
                "unit": "fraction_difference",
                "correction": "BH6_directional_sign_test",
                "before_effect": 0.00021,
                "before_mean_ci_low": -0.00017,
                "before_mean_ci_high": 0.0006,
                "after_effect": 0.00021,
                "after_mean_ci_low": -0.00017,
                "after_mean_ci_high": 0.0006,
                "matched_change_effect": 0,
                "matched_change_mean_ci_low": 0,
                "matched_change_mean_ci_high": 0,
                "before_q": 0.0029,
                "after_q": 0.0029,
            },
            {
                "scenario": "exclude_genus_NCBITaxon:561",
                "scope": "pooled",
                "representation": "shared_fraction",
                "trait_key": "METPO:2000103|pathways:arsenate_reduction",
                "denominator_trait": "",
                "unit": "fraction_difference",
                "correction": "BH6_directional_sign_test",
                "before_effect": 0.000064,
                "before_mean_ci_low": 0.000033,
                "before_mean_ci_high": 0.000096,
                "after_effect": -0.000007,
                "after_mean_ci_low": -0.000019,
                "after_mean_ci_high": 0.000005,
                "matched_change_effect": -0.000065,
                "matched_change_mean_ci_low": -0.000092,
                "matched_change_mean_ci_high": -0.000037,
                "before_q": 0.0038,
                "after_q": 0.65,
            },
        ],
    )
    write_table(
        source / "contributing_universe.tsv",
        ["scope", "direction", "trait_key", "rank", "taxon", "name", "share"],
        [
            {
                "scope": "pooled",
                "direction": "inc",
                "trait_key": "METPO:2000103|pathways:arsenate_reduction",
                "rank": "genus",
                "taxon": "NCBITaxon:561",
                "name": "Escherichia",
                "share": 0.587,
            }
        ],
    )
    (source / "native_summary.md").write_text(
        "# Native summary\nSee native_manifest.json and native_composition_comparison.tsv.\n"
    )
    output_hashes["native_summary.md"] = digest(source / "native_summary.md")
    write_table(
        source / "scenarios.tsv",
        ["scenario", "name", "status"],
        [
            {"scenario": "exclude_A", "name": "Genus A", "status": "run"},
            {
                "scenario": "exclude_B",
                "name": "Genus B",
                "status": "outside_rule_declared_panel",
            },
        ],
    )
    write_table(
        source / "source_contributions.tsv.gz",
        ["reported_taxon", "trait_key", "count_contribution"],
        [{"reported_taxon": "taxon", "trait_key": "A", "count_contribution": 1}],
    )
    write_table(
        source / "profile_reported_taxa.tsv.gz",
        ["study_id", "experiment_id", "direction", "reported_taxon"],
        [
            {
                "study_id": "study",
                "experiment_id": "experiment",
                "direction": "inc",
                "reported_taxon": "taxon",
            }
        ],
    )
    write_table(
        source / "profile_counts.tsv.gz",
        ["study_id", "experiment_id", "direction", "trait_key", "profile_count"],
        [
            {
                "study_id": "study",
                "experiment_id": "experiment",
                "direction": "inc",
                "trait_key": "A",
                "profile_count": 1,
            }
        ],
    )
    (source / "native_manifest.json").write_text(
        json.dumps(
            {
                "completed": True,
                "code_dirty": False,
                "contribution_manifest_sha256": digest(source / "manifest.json"),
                "results": {
                    "input_sha256": {"scenarios.tsv": digest(source / "scenarios.tsv")},
                    "output_sha256": output_hashes,
                },
            }
        )
    )
    return source


def test_publication_preserves_comparisons_and_full_family_metadata(run_directory: Path) -> None:
    destination = run_directory.parent / "review"
    result = publish(run_directory, destination)
    tables = result["tables"]
    assert isinstance(tables, dict)
    assert len(tables) == 6
    assert (destination / "scenarios.tsv").read_bytes() == (
        run_directory / "scenarios.tsv"
    ).read_bytes()
    selected = list(compressed_rows(destination / "native_table1_stats.tsv.gz"))
    assert selected == [{"trait_key": "A", "value": "1"}]
    comparisons = list(compressed_rows(destination / "native_composition_comparison.tsv.gz"))
    assert len(comparisons) == 3
    handoff = list(rows(destination / "interpretation_handoff.tsv"))
    assert [row["finding_id"] for row in handoff] == [
        "arsenate_escherichia_sensitive",
        "sphere_raw_fraction_mean_interval_crosses_zero",
    ]
    assert handoff[0]["increased_arm_count_share"] == "0.587"
    assert handoff[0]["native_after_q"] == "0.67"
    assert handoff[1]["raw_fraction_baseline_ci_low"] == "-0.00017"
    reconciliation = result["source_reconciliation"]
    assert isinstance(reconciliation, dict)
    assert reconciliation["source_contribution_rows"] == 1
    assert reconciliation["mismatches"] == 0
    assert "FULL native families" in (destination / "README.md").read_text()
    summary = (destination / "README.md").read_text()
    assert "See run_manifest.json and native_composition_comparison.tsv.gz." in summary
    assert "Cross-scale interpretation flags (#213)" in summary
    assert "directional sign-test q is 0.0029" in summary


@pytest.mark.parametrize(
    ("version", "selected", "expected_size", "basis"),
    [
        ("190-native-v3", "4", "4", "all_selected_traits_padded_at_one"),
        ("190-native-v1", "4", "2", "historical_v1_finite_tests"),
        ("190-native-v1", None, "2", "historical_v1_finite_tests"),
        ("190-native-v2", "4", "2", "historical_finite_tests"),
        ("190-native-v2", None, "2", "historical_finite_tests"),
        (None, None, "2", "legacy_finite_tests_no_protocol_version"),
    ],
)
def test_handoff_distinguishes_correction_family_from_tested_count(
    run_directory: Path,
    version: str | None,
    selected: str | None,
    expected_size: str,
    basis: str,
) -> None:
    table = run_directory / "native_composition_comparison.tsv"
    records = list(rows(table))
    for record in records:
        record["native_before_family_tests"] = "2"
        if selected is not None:
            record["native_before_family_selected_traits"] = selected
    write_table(table, list(records[0]), records)
    manifest_path = run_directory / "native_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if version is not None:
        manifest["protocol"] = {"version": version}
    manifest["results"]["output_sha256"][table.name] = digest(table)
    manifest_path.write_text(json.dumps(manifest))
    original_manifest = manifest_path.read_bytes()

    destination = run_directory.parent / "review"
    refusal = ""
    try:
        publish(run_directory, destination)
    except ValueError as error:
        refusal = str(error)
    assert not refusal, f"Known historical/current protocol {version} refused: {refusal}"
    assert manifest_path.read_bytes() == original_manifest
    published_manifest = json.loads((destination / "run_manifest.json").read_text())
    assert published_manifest.get("protocol") == manifest.get("protocol")
    handoff = list(rows(destination / "interpretation_handoff.tsv"))
    assert len(handoff) == 2
    assert [record["native_family_size"] for record in handoff] == [expected_size] * 2
    assert [record.get("native_tested_count") for record in handoff] == ["2"] * 2
    assert [record.get("native_family_basis") for record in handoff] == [basis] * 2
    assert [record["native_baseline_q"] for record in handoff] == ["3e-05"] * 2


@pytest.mark.parametrize("version", ["190-native-v3", "unknown-native-protocol"])
def test_handoff_refuses_undeclared_correction_family(
    run_directory: Path,
    version: str,
) -> None:
    manifest_path = run_directory / "native_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["protocol"] = {"version": version}
    manifest_path.write_text(json.dumps(manifest))
    message = ""
    try:
        interpretation_handoff(run_directory)
    except ValueError as exc:
        message = str(exc)
    assert "family" in message or "protocol" in message


def test_publication_rejects_changed_outputs_before_writing(run_directory: Path) -> None:
    (run_directory / "native_summary.md").write_text("changed")
    destination = run_directory.parent / "review"
    with pytest.raises(ValueError, match="artifact changed"):
        publish(run_directory, destination)
    assert not destination.exists()


def test_publication_rejects_source_contribution_mismatch(run_directory: Path) -> None:
    write_table(
        run_directory / "source_contributions.tsv.gz",
        ["reported_taxon", "trait_key", "count_contribution"],
        [{"reported_taxon": "taxon", "trait_key": "A", "count_contribution": 2}],
    )
    with pytest.raises(ValueError, match="do not reconcile"):
        publish(run_directory, run_directory.parent / "review")


def test_publication_rejects_incomplete_runs(run_directory: Path) -> None:
    path = run_directory / "native_manifest.json"
    manifest = json.loads(path.read_text())
    manifest["completed"] = False
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Only completed"):
        publish(run_directory, run_directory.parent / "review")


def test_publication_rejects_dirty_or_legacy_provenance(run_directory: Path) -> None:
    path = run_directory / "native_manifest.json"
    manifest = json.loads(path.read_text())
    destination = run_directory.parent / "review"
    manifest["code_dirty"] = True
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="dirty or unrecorded"):
        publish(run_directory, destination)
    assert not destination.exists()

    del manifest["code_dirty"]
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="dirty or unrecorded"):
        publish(run_directory, destination)


def test_publication_rejects_paths_outside_the_run(run_directory: Path) -> None:
    path = run_directory / "native_manifest.json"
    manifest = json.loads(path.read_text())
    manifest["results"]["input_sha256"] = {"../outside": "unknown"}
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="escapes"):
        publish(run_directory, run_directory.parent / "review")
