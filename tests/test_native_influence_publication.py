"""Publication must preserve complete six-trait tables and reject changed runs."""

import json
from pathlib import Path

import pytest

from bsdbng.contribution_inputs import digest
from bsdbng.native_influence_publication import publish
from bsdbng.native_taxon_influence import compressed_rows
from bsdbng.taxon_contributions import write_table


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
        "native_composition_comparison.tsv",
        "native_ora_fold_comparison.tsv",
        "native_ora_families.tsv",
    ):
        write_table(
            source / filename,
            ["trait_key", "value"],
            [{"trait_key": "A", "value": 1}, {"trait_key": "OTHER", "value": 2}],
        )
        output_hashes[filename] = digest(source / filename)
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
    (source / "native_manifest.json").write_text(
        json.dumps(
            {
                "completed": True,
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
    assert len(tables) == 5
    assert (destination / "scenarios.tsv").read_bytes() == (
        run_directory / "scenarios.tsv"
    ).read_bytes()
    selected = list(compressed_rows(destination / "native_table1_stats.tsv.gz"))
    assert selected == [{"trait_key": "A", "value": "1"}]
    comparisons = list(compressed_rows(destination / "native_composition_comparison.tsv.gz"))
    assert len(comparisons) == 2
    assert "FULL native families" in (destination / "README.md").read_text()
    summary = (destination / "README.md").read_text()
    assert "See run_manifest.json and native_composition_comparison.tsv.gz." in summary


def test_publication_rejects_changed_outputs_before_writing(run_directory: Path) -> None:
    (run_directory / "native_summary.md").write_text("changed")
    destination = run_directory.parent / "review"
    with pytest.raises(ValueError, match="artifact changed"):
        publish(run_directory, destination)
    assert not destination.exists()


def test_publication_rejects_incomplete_runs(run_directory: Path) -> None:
    path = run_directory / "native_manifest.json"
    manifest = json.loads(path.read_text())
    manifest["completed"] = False
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Only completed"):
        publish(run_directory, run_directory.parent / "review")


def test_publication_rejects_paths_outside_the_run(run_directory: Path) -> None:
    path = run_directory / "native_manifest.json"
    manifest = json.loads(path.read_text())
    manifest["results"]["input_sha256"] = {"../outside": "unknown"}
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="escapes"):
        publish(run_directory, run_directory.parent / "review")
