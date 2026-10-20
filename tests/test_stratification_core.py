"""Bind anatomical partitions and all downstream products to actual core inputs (#727)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from mtp import stratification
from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow.consumers import CoreBundle, open_core_bundle, verify_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow


@pytest.fixture
def core_bundle(tmp_path: Path, request: pytest.FixtureRequest) -> CoreBundle:
    source = REPO_ROOT / "tests" / "fixtures" / "workflow"
    raw = tmp_path / "full_dump.csv"
    with (source / "full_dump.csv").open() as handle:
        comment = handle.readline()
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        assert fields is not None
        rows = list(reader)
    if getattr(request, "param", "") == "zero_eligible":
        rows = [
            row
            for row in rows
            if row["Host species"] != "Homo sapiens" or row["Abundance in Group 1"] == "increased"
        ]
    for row in rows:
        if row["Study"] == "Study 3":
            row["Body site"] = "saliva"
    with raw.open("w", newline="") as handle:
        handle.write(comment)
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    run = execute_workflow(
        REPO_ROOT / "workflows" / "bugsigdb_core.yaml",
        {
            "raw_dump": raw,
            "kg_edges": source / "merged-kg_edges.tsv",
            "kg_nodes": source / "merged-kg_nodes.tsv",
        },
        tmp_path / "core",
        OPERATIONS,
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    return open_core_bundle(Path(run.manifest_path), "ora_profiles")


def test_core_strata_records_every_actual_product_and_source(
    core_bundle: CoreBundle, tmp_path: Path
) -> None:
    output = tmp_path / "analysis"
    assert (
        stratification.main(["--core-run", str(core_bundle.manifest), "--output-dir", str(output)])
        == 0
    )
    record = verify_consumer_record(output / "core_consumption.json")
    assert record.consumer == "mtp.stratification"
    assert record.cohort == "Homo sapiens"
    assert {"ora_profiles", "experiments", "taxa_matrix"}.issubset(record.product_ids)
    assert record.coverage_scope == "input_cohort_before_analysis_specific_exclusions"
    expected = {path.resolve() for path in output.iterdir() if path.name != "core_consumption.json"}
    assert {Path(item.path) for item in record.outputs} == expected
    assert len(expected) == 25
    for stratum in stratification.PRIMARY_STRATA:
        partition = output / f"study_trait_profiles_{stratum}.tsv"
        complete = output / f"trait_enrichment_all_results_{stratum}.tsv"
        assert f"# Profiles-SHA256: {file_sha256(partition)}" in complete.read_text()


@pytest.mark.parametrize("argument", ["profiles", "experiments", "taxa-matrix"])
def test_core_strata_rejects_explicit_conflicting_inputs(
    core_bundle: CoreBundle, tmp_path: Path, argument: str
) -> None:
    output = tmp_path / "analysis"
    with pytest.raises(ValueError, match=f"conflicts with explicit --{argument}"):
        stratification.main(
            [
                "--core-run",
                str(core_bundle.manifest),
                "--output-dir",
                str(output),
                f"--{argument}",
                str(tmp_path / "unrelated.tsv"),
            ]
        )
    assert not output.exists()


@pytest.mark.parametrize(
    "collision", ["core_directory", "input_symlink", "input_hardlink", "sidecar_alias"]
)
def test_core_strata_preflight_preserves_inputs_and_sidecars(
    core_bundle: CoreBundle, tmp_path: Path, collision: str
) -> None:
    output = tmp_path / "analysis"
    protected = core_bundle.path("experiments")
    before = protected.read_bytes()
    if collision == "core_directory":
        output = core_bundle.manifest.parent
    else:
        output.mkdir()
        partition = output / "study_trait_profiles_pooled.tsv"
        if collision == "input_symlink":
            partition.symlink_to(protected)
        elif collision == "input_hardlink":
            partition.hardlink_to(protected)
        else:
            partition.write_text("previous result\n")
            (output / "core_consumption.json").hardlink_to(partition)
    contents = {path.name: path.read_bytes() for path in output.iterdir() if path.is_file()}
    with pytest.raises(ValueError, match="Consumer output"):
        stratification.main(["--core-run", str(core_bundle.manifest), "--output-dir", str(output)])
    assert protected.read_bytes() == before
    assert {path.name: path.read_bytes() for path in output.iterdir() if path.is_file()} == contents


@pytest.mark.parametrize("changed", ["experiments", "partition", "result"])
def test_core_strata_binding_rejects_changed_inputs_or_outputs(
    core_bundle: CoreBundle, tmp_path: Path, changed: str
) -> None:
    output = tmp_path / "analysis"
    stratification.main(["--core-run", str(core_bundle.manifest), "--output-dir", str(output)])
    record = output / "core_consumption.json"
    verify_consumer_record(record)
    path = {
        "experiments": core_bundle.path("experiments"),
        "partition": output / "study_trait_profiles_feces.tsv",
        "result": output / "trait_enrichment_all_results_feces.tsv",
    }[changed]
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="changed"):
        verify_consumer_record(record)


@pytest.mark.parametrize(
    ("core_bundle", "status"),
    [("zero_eligible", "zero_eligible_tests"), ("zero_discoveries", "zero_discoveries")],
    indirect=["core_bundle"],
)
def test_core_strata_preserves_valid_empty_results(
    core_bundle: CoreBundle, tmp_path: Path, status: str
) -> None:
    output = tmp_path / "analysis"
    stratification.main(["--core-run", str(core_bundle.manifest), "--output-dir", str(output)])
    verify_consumer_record(output / "core_consumption.json")
    for stratum in stratification.PRIMARY_STRATA:
        complete = output / f"trait_enrichment_all_results_{stratum}.tsv"
        assert f"# Result-Status: {status}" in complete.read_text()
    record = json.loads((output / "core_consumption.json").read_text())
    assert record["outputs"]
