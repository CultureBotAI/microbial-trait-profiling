"""Downstream contracts reject plausible-looking but incompatible or stale inputs (#629)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT
from mtp.workflow.__main__ import describe, main
from mtp.workflow.consumers import (
    import_legacy_provenance,
    open_core_bundle,
    verify_consumer_record,
    write_consumer_record,
)
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow


@pytest.fixture
def core_manifest(tmp_path: Path) -> Path:
    source = REPO_ROOT / "tests" / "fixtures" / "workflow"
    run = execute_workflow(
        REPO_ROOT / "workflows" / "bugsigdb_core.yaml",
        {
            "raw_dump": source / "full_dump.csv",
            "kg_edges": source / "merged-kg_edges.tsv",
            "kg_nodes": source / "merged-kg_nodes.tsv",
        },
        tmp_path / "runs",
        OPERATIONS,
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    return Path(run.manifest_path)


def test_verified_consumption_and_changed_output(core_manifest: Path, tmp_path: Path) -> None:
    bundle = open_core_bundle(core_manifest)
    result = tmp_path / "result.tsv"
    result.write_text("statistic\tvalue\nexample\t1\n")
    provenance = write_consumer_record(
        bundle, tmp_path / "binding.json", consumer="test-analysis", outputs=[result]
    )
    record = verify_consumer_record(provenance)
    assert "count_profiles" in record.product_ids
    assert record.cohort == "all"
    result.write_text("statistic\tvalue\nexample\t2\n")
    with pytest.raises(ValueError, match="Consumer output changed"):
        verify_consumer_record(provenance)


def test_additional_consumed_product_is_bound_and_tamper_detected(
    core_manifest: Path, tmp_path: Path
) -> None:
    bundle = open_core_bundle(core_manifest)
    result = tmp_path / "result.tsv"
    result.write_text("value\n1\n")
    binding = write_consumer_record(
        bundle,
        tmp_path / "binding.json",
        consumer="reported-taxon-audit",
        outputs=[result],
        additional_products=("reported_taxa",),
    )
    record = verify_consumer_record(binding)
    assert record.product_ids.count("reported_taxa") == 1
    # Corrupt only the additional binding: the core and every artifact remain
    # valid, so core re-verification cannot accidentally make this test pass.
    original = json.loads(binding.read_text())
    index = original["product_ids"].index("reported_taxa")
    original["artifact_ids"][index] = original["artifact_ids"][0]
    binding.write_text(json.dumps(original))
    open_core_bundle(core_manifest)
    with pytest.raises(ValueError, match="unresolved product/artifact"):
        verify_consumer_record(binding)


@pytest.mark.parametrize(
    "extras",
    [
        ("reported_taxa", "reported_taxa"),
        ("count_profiles",),
        ("missing",),
        ("ora_profiles",),
    ],
)
def test_additional_product_errors_precede_sidecar_writes(
    core_manifest: Path, tmp_path: Path, extras: tuple[str, ...]
) -> None:
    bundle = open_core_bundle(core_manifest)
    result = tmp_path / "result.tsv"
    result.write_text("value\n1\n")
    binding = tmp_path / "binding.json"
    with pytest.raises(ValueError, match="Additional"):
        write_consumer_record(
            bundle, binding, consumer="audit", outputs=[result], additional_products=extras
        )
    assert not binding.exists()
    assert not binding.with_suffix(".coverage.tsv").exists()


def test_units_cohorts_and_changed_same_name_file(core_manifest: Path) -> None:
    count = open_core_bundle(core_manifest)
    ora = open_core_bundle(core_manifest, "ora_profiles")
    assert count.path("count_profiles") != ora.path("ora_profiles")
    with pytest.raises(ValueError, match="cohort"):
        open_core_bundle(core_manifest, "ora_profiles", expected_cohort="Mus musculus")
    path = count.path("count_profiles")
    path.write_text(ora.path("ora_profiles").read_text())
    with pytest.raises(ValueError, match="content or inventory changed"):
        open_core_bundle(core_manifest)


def test_old_provenance_never_becomes_an_execution(tmp_path: Path) -> None:
    source = tmp_path / "old.yaml"
    source.write_text("pipeline:\n  - command: mtp.cli nonexistent\n")
    destination = import_legacy_provenance(source, tmp_path / "legacy.json")
    record = json.loads(destination.read_text())
    assert record["verification_status"] == "imported_unverified"
    assert "nonexistent" in record["claims_json"]
    with pytest.raises(ValueError, match="validation errors for WorkflowRun"):
        open_core_bundle(destination)
    with pytest.raises(FileExistsError):
        import_legacy_provenance(source, destination)


def test_cli_validation_and_schema_derived_docs(
    core_manifest: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(["validate", str(core_manifest)]) == 0
    assert "scientific contracts" in capsys.readouterr().out
    assert main(["lineage", str(core_manifest), "ora_profiles"]) == 0
    lineage = json.loads(capsys.readouterr().out)
    assert {item["product_id"] for item in lineage["artifacts"]} >= {
        "raw_dump",
        "kg_edges",
        "kg_nodes",
        "studies",
        "taxa_matrix",
        "ora_profiles",
    }
    assert (
        describe(REPO_ROOT / "workflows" / "bugsigdb_core.yaml")
        == (REPO_ROOT / "docs" / "workflow-products.md").read_text()
    )


def test_ora_uses_catalog_media_filter_case_insensitively(tmp_path: Path) -> None:
    profiles, matrix = tmp_path / "profiles.tsv", tmp_path / "matrix.tsv"
    medium = "METPO:2000510|MediaDive.Medium:1"
    profiles.write_text(
        f"study_id\texperiment_id\tinc_dec\ttaxa_count\tMETPO:1\t{medium}\n"
        "bsdb:1\tbsdb:1-1\tinc\t1\t1\t1\n"
        "bsdb:1\tbsdb:1-1\tdec\t1\t0\t1\n"
    )
    matrix.write_text(f"taxon_id\tMETPO:1\t{medium}\nNCBITaxon:1\t1\t1\n")
    full = tmp_path / "full.tsv"
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "trait_enrichment.py"),
            "--input",
            str(profiles),
            "--taxa-matrix",
            str(matrix),
            "--output",
            str(tmp_path / "report.txt"),
            "--all-results-output",
            str(full),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert medium not in full.read_text()
    assert "METPO:1" in full.read_text()
