"""Fully verify a portable policy-less historical core without mocks (#1134)."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import tarfile
from collections.abc import Callable
from pathlib import Path

import pytest

from mtp.datamodel.workflow_pydantic import WorkflowRun
from mtp.workflow.coverage import write_feature_catalog
from mtp.workflow.runtime import artifact_path, inventory
from mtp.workflow.trait_policy import resolve_trait_policy
from mtp.workflow.verification import verify_core_run

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/trait_key_policy/historical_raw_core"
UNION_FIXTURE = ROOT / "tests/fixtures/trait_key_policy/historical_union_core.tar.gz"
NAMESPACE = Path("/fixture/trait-policy-511c891")
UNION_NAMESPACE = Path("/fixture/trait-policy-44897c6")
MEMBERS = {"METPO:1000605", "METPO:1000608", "METPO:1000612"}


def resolver(root: Path, namespace: Path = NAMESPACE) -> Callable[[Path], Path]:
    """Map only this committed test bundle; never consult original source paths."""
    root = root.resolve()

    def resolve(path: Path) -> Path:
        if path.is_relative_to(namespace):
            target = root / path.relative_to(namespace)
        elif path.is_relative_to(root):
            target = path
        else:
            raise ValueError(f"Path is outside the portable historical fixture: {path}")
        resolved = target.resolve(strict=True)
        if not resolved.is_relative_to(root):
            raise ValueError("Portable historical fixture path escapes its root")
        return resolved

    return resolve


def extract_union_fixture(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(UNION_FIXTURE, mode="r:gz") as tar:
        tar.extractall(destination, filter="data")


def payloads(root: Path) -> dict[str, bytes]:
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in root.rglob("*")
        if path.is_file()
    }


def test_real_verifier_accepts_relocated_policyless_raw_core_without_changing_bytes(
    tmp_path: Path,
) -> None:
    copied = tmp_path / "relocated historical core"
    shutil.copytree(FIXTURE, copied)
    before = payloads(copied)
    manifest = copied / "core/run.json"
    report = json.loads((copied / "generation.json").read_text(encoding="utf-8"))
    assert hashlib.sha256(manifest.read_bytes()).hexdigest() == report["normalized_manifest_sha256"]
    assert report["original_generation_verified"] is True
    assert report["normalized_historical_verification"] is True
    helper = ROOT / "scripts/adhoc/issue1134_generate_historical_core_fixture.py"
    assert hashlib.sha256(helper.read_bytes()).hexdigest() == report["helper_sha256"]
    record = json.loads(manifest.read_bytes())
    assert "trait_key_policy" not in record
    assert all(item["name"] != "trait_key_policy" for item in record["parameters"])
    assert all(item["name"] != "trait_key_policy" for item in record["definition"]["parameters"])
    run = verify_core_run(NAMESPACE / "core/run.json", resolver=resolver(copied))
    assert run.status == "completed"
    assert run.trait_key_policy is None
    assert run.git_commit == report["producer_commit"]
    assert not run.git_dirty
    policy = resolve_trait_policy(run)
    assert policy.policy == "raw_object_v1"
    assert policy.include_policy_metadata is False
    implementations = json.loads(
        (FIXTURE.parent / "historical_implementations.json").read_text(encoding="utf-8")
    )
    audited = next(row for row in implementations["records"] if row["policy"] == "raw_object_v1")
    assert report["producer_commit"] == audited["source_commit"]
    matrix_source = next(
        row for row in run.code_inventory if row.path == audited["code_inventory"]["path"]
    )
    assert matrix_source.model_dump(mode="json") == audited["code_inventory"]
    with artifact_path(run, "taxa_matrix", resolver(copied)).open(
        encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert MEMBERS.issubset(reader.fieldnames or [])
        rows = {row["taxon_id"]: row for row in reader}
    assert all(int(rows["NCBITaxon:100"][trait]) == 1 for trait in MEMBERS)
    assert int(rows["NCBITaxon:101"]["METPO:1000608"]) == 1
    assert int(rows["NCBITaxon:103"]["METPO:1000605"]) == 1
    with artifact_path(run, "feature_metadata", resolver(copied)).open(
        encoding="utf-8", newline=""
    ) as handle:
        metadata = csv.DictReader(handle, delimiter="\t")
        assert len(metadata.fieldnames or []) == 7
        assert MEMBERS.issubset({row["column_id"] for row in metadata})
    assert payloads(copied) == before


def test_historical_inventory_corruption_is_rejected(tmp_path: Path) -> None:
    copied = tmp_path / "corrupted historical core"
    shutil.copytree(FIXTURE, copied)
    manifest = copied / "core/run.json"
    run = WorkflowRun.model_validate_json(manifest.read_bytes())
    matrix = artifact_path(run, "taxa_matrix", resolver(copied))
    matrix.write_bytes(matrix.read_bytes() + b"corruption\n")
    with pytest.raises(ValueError, match="content or inventory changed"):
        verify_core_run(NAMESPACE / "core/run.json", resolver=resolver(copied))


def test_historical_scientific_matrix_reconciliation_is_not_bypassed(tmp_path: Path) -> None:
    copied = tmp_path / "scientifically wrong historical core"
    shutil.copytree(FIXTURE, copied)
    manifest = copied / "core/run.json"
    run = WorkflowRun.model_validate_json(manifest.read_bytes())
    matrix = artifact_path(run, "taxa_matrix", resolver(copied))
    with matrix.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        assert header is not None
        rows = list(reader)
    genus = next(row for row in rows if row["taxon_id"] == "NCBITaxon:100")
    genus["METPO:1000612"] = str(int(genus["METPO:1000612"]) + 1)
    with matrix.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    artifact = next(item for item in run.artifacts or [] if item.product_id == "taxa_matrix")
    artifact.sha256, artifact.size_bytes, artifact.inventory = inventory(matrix, "file")
    manifest.write_text(run.model_dump_json(indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Matrix count differs from consumed KG evidence"):
        verify_core_run(NAMESPACE / "core/run.json", resolver=resolver(copied))


def test_real_verifier_accepts_relocated_union_core_produced_by_main_without_changing_bytes(
    tmp_path: Path,
) -> None:
    """Verify a portable facultative_union_v1 core produced by main (44897c6) (#1231, #1235)."""
    copied = tmp_path / "relocated historical union core"
    extract_union_fixture(copied)
    before = payloads(copied)
    manifest = copied / "core/run.json"
    report = json.loads((copied / "generation.json").read_text(encoding="utf-8"))
    assert hashlib.sha256(manifest.read_bytes()).hexdigest() == report["normalized_manifest_sha256"]
    assert report["original_generation_verified"] is True
    assert report["normalized_historical_verification"] is True
    assert isinstance(report.get("helper_sha256"), str)
    assert len(report["helper_sha256"]) == 64
    helper_path = ROOT / "scripts/adhoc/issue1231_generate_historical_union_core_fixture.py"
    assert report["helper_sha256"] == hashlib.sha256(helper_path.read_bytes()).hexdigest()
    assert (FIXTURE.parent / "README.md").is_file()
    assert all("original_value_sha256" not in item for item in report["location_transformations"])
    run = verify_core_run(
        UNION_NAMESPACE / "core/run.json", resolver=resolver(copied, UNION_NAMESPACE)
    )
    assert run.status == "completed"
    assert run.trait_key_policy == "facultative_union_v1"
    assert run.git_commit == report["producer_commit"]
    assert not run.git_dirty
    policy = resolve_trait_policy(run)
    assert policy.policy == "facultative_union_v1"
    assert policy.include_policy_metadata is True
    with artifact_path(run, "feature_metadata", resolver(copied, UNION_NAMESPACE)).open(
        encoding="utf-8", newline=""
    ) as handle:
        metadata_rows = list(csv.DictReader(handle, delimiter="\t"))
    focal = next(row for row in metadata_rows if row["column_id"] == "METPO:1000608")
    assert focal["object_name"] == "facultatively aerobic"
    assert focal["column_label"] == "facultatively aerobic METPO:1000608"
    assert json.loads(focal["trait_member_ids_json"]) == sorted(MEMBERS)
    assert payloads(copied) == before


def test_union_core_verification_rejects_metadata_relabel_of_facultative_column(
    tmp_path: Path,
) -> None:
    """Mutating recorded object_name/column_label in a union core fails verification (#1231)."""
    copied = tmp_path / "relabelled union core"
    extract_union_fixture(copied)
    manifest = copied / "core/run.json"
    run = WorkflowRun.model_validate_json(manifest.read_bytes())
    resolve = resolver(copied, UNION_NAMESPACE)
    feature_meta = artifact_path(run, "feature_metadata", resolve)
    catalog_meta = artifact_path(run, "feature_catalog", resolve)
    text = feature_meta.read_text(encoding="utf-8").replace(
        "facultatively aerobic", "facultative (aerobe or anaerobe)"
    )
    feature_meta.write_text(text, encoding="utf-8")
    write_feature_catalog(
        feature_meta,
        catalog_meta,
        trait_key_policy="facultative_union_v1",
        include_policy_metadata=True,
    )
    for prod_id, path in (("feature_metadata", feature_meta), ("feature_catalog", catalog_meta)):
        artifact = next(item for item in run.artifacts or [] if item.product_id == prod_id)
        artifact.sha256, artifact.size_bytes, artifact.inventory = inventory(path, "file")
    manifest.write_text(run.model_dump_json(indent=2) + "\n", encoding="utf-8")
    with pytest.raises(
        ValueError, match="Feature metadata differs from consumed KG evidence and nodes"
    ):
        verify_core_run(UNION_NAMESPACE / "core/run.json", resolver=resolve)
