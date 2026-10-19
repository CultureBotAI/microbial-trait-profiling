"""Bind integrated evidence to executed core inputs and preserve source files (#551)."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from mtp.integrated_evidence import _protect_source_files, main, ora_testable_profile_keys
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import OXYGEN_GUILDS, TABLE1_TRAITS
from mtp.workflow.consumers import open_core_bundle, verify_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow


@pytest.fixture
def core_manifest(tmp_path: Path) -> Path:
    sources = REPO_ROOT / "tests/fixtures/workflow"
    edges = tmp_path / "edges.tsv"
    nodes = tmp_path / "nodes.tsv"
    terms = sorted({term for group in (OXYGEN_GUILDS | TABLE1_TRAITS).values() for term in group})
    edge_rows = (sources / "merged-kg_edges.tsv").read_text()
    node_rows = (sources / "merged-kg_nodes.tsv").read_text()
    node_rows += "METPO:2000103\tbiolink:OntologyClass\tcapable of\n"
    for index, term in enumerate(terms):
        if term == "METPO:1000608":
            continue
        predicate = "biolink:has_phenotype" if term.startswith("METPO:") else "METPO:2000103"
        edge_rows += f"integrated{index}\tNCBITaxon:101\t{predicate}\t{term}\tfixture:annotation\n"
        node_rows += f"{term}\tbiolink:PhenotypicQuality\tFixture trait {index}\n"
    edges.write_text(edge_rows)
    nodes.write_text(node_rows)
    run = execute_workflow(
        REPO_ROOT / "workflows/bugsigdb_core.yaml",
        {
            "raw_dump": sources / "full_dump.csv",
            "kg_edges": edges,
            "kg_nodes": nodes,
        },
        tmp_path / "runs",
        OPERATIONS,
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    return Path(run.manifest_path)


def test_integrated_core_binds_outputs_and_does_not_read_checkout_assessments(
    core_manifest: Path, tmp_path: Path
) -> None:
    output = tmp_path / "integrated"
    assert main(["--core-run", str(core_manifest), "--output-dir", str(output)]) == 0
    record = verify_consumer_record(output / "core_consumption.json")
    provenance = json.loads((output / "integrated_evidence_provenance.json").read_text())
    assert provenance["core_run_id"] == record.run_id
    assert provenance["core_manifest_sha256"] == record.run_manifest_sha256
    assert record.consumer == "mtp.integrated_evidence"
    assert {Path(item.path).name for item in record.outputs} >= {
        "integrated_evidence.tsv",
        "integrated_evidence_method_results.tsv",
        "integrated_evidence_provenance.json",
    }
    assert all(
        not item["verdicts_by_procedure"] for item in provenance["assessment_status"].values()
    )
    assert all(source["sha256"] is None for source in provenance["native_sources"].values())


def test_integrated_core_rejects_conflicting_input_before_output(
    core_manifest: Path, tmp_path: Path
) -> None:
    output = tmp_path / "integrated"
    with pytest.raises(ValueError, match="conflicts"):
        main(
            [
                "--core-run",
                str(core_manifest),
                "--profile-input",
                str(tmp_path / "different.tsv"),
                "--output-dir",
                str(output),
            ]
        )
    assert not output.exists()


def test_explicit_missing_native_source_is_not_silently_marked_unsupplied(
    core_manifest: Path, tmp_path: Path
) -> None:
    output = tmp_path / "integrated"
    with pytest.raises(SystemExit):
        main(
            [
                "--core-run",
                str(core_manifest),
                "--arm2-enrichment-results",
                str(tmp_path / "missing.tsv"),
                "--output-dir",
                str(output),
            ]
        )
    assert not output.exists()


def test_integrated_core_cannot_write_inside_input_bundle(core_manifest: Path) -> None:
    bundle = open_core_bundle(core_manifest)
    with pytest.raises(ValueError, match="overlaps"):
        main(
            [
                "--core-run",
                str(core_manifest),
                "--output-dir",
                str(bundle.path("count_profiles").parent),
            ]
        )


@pytest.mark.parametrize("alias", ["inside", "hardlink", "symlink"])
def test_native_evidence_is_protected_from_output_aliases(tmp_path: Path, alias: str) -> None:
    output = tmp_path / "output"
    output.mkdir()
    source = output / "source.tsv" if alias == "inside" else tmp_path / "source.tsv"
    source.write_text("source evidence\n")
    if alias == "hardlink":
        (output / "integrated_evidence.tsv").hardlink_to(source)
    elif alias == "symlink":
        (output / "integrated_evidence.tsv").symlink_to(source)
    with pytest.raises(ValueError, match="input evidence"):
        _protect_source_files(output, [source])
    assert source.read_text() == "source evidence\n"


def test_ora_eligibility_pairs_before_excluding_zero_draws(tmp_path: Path) -> None:
    """N is matched reported taxa, independent of the retained annotation total."""
    source = tmp_path / "ora.tsv"
    source.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTRAIT\n"
        "s1\te1\tinc\t0\t0\n"
        "s1\te1\tdec\t3\t0\n"
        "s2\te1\tinc\t2\t1\n"
        "s3\te1\tinc\t2\t1\n"
        "s3\te1\tdec\t2\t1\n"
    )
    assert ora_testable_profile_keys(source) == {
        ("s1", "e1", "dec"),
        ("s3", "e1", "inc"),
        ("s3", "e1", "dec"),
    }


def test_integrated_core_reads_actual_ora_test_cohort(core_manifest: Path, tmp_path: Path) -> None:
    """Import complete production ORA output using the actual core ORA unit."""
    bundle = open_core_bundle(core_manifest)
    ora_path = bundle.path("ora_profiles")
    complete = tmp_path / "all_ora.tsv"
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts/trait_enrichment.py"),
            "--input",
            str(ora_path),
            "--taxa-matrix",
            str(bundle.path("taxa_matrix")),
            "--output",
            str(tmp_path / "ora-report.txt"),
            "--tsv-output",
            str(tmp_path / "ora-selected.tsv"),
            "--all-results-output",
            str(complete),
            "--p-threshold",
            "0.001",
            "--min-fold",
            "1.5",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    output = tmp_path / "integrated"
    assert (
        main(
            [
                "--core-run",
                str(core_manifest),
                "--arm2-enrichment-results",
                str(complete),
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )
    with (output / "integrated_evidence.tsv").open() as handle:
        records = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["method"] == "enrichment_arm_profile_hypergeometric"
        ]
    assert len(records) == 6
    expected = ora_testable_profile_keys(ora_path)
    assert expected
    assert all(int(row["n_tests"]) == len(expected) for row in records)
    assert all(row["test_status"] == "tested_profiles_summarized" for row in records)
    provenance = json.loads((output / "integrated_evidence_provenance.json").read_text())
    assert provenance["ora_eligibility_source"]["product_id"] == "ora_profiles"
    assert provenance["ora_eligibility_source"]["path"] == str(ora_path)
    assert provenance["ora_eligibility_source"]["sha256"] == file_sha256(ora_path)
    verify_consumer_record(output / "core_consumption.json")
