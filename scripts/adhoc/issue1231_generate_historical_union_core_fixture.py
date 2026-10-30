"""Generate and verify a portable historical union-core fixture (#1231, #1235, 2026-10-30).

Inputs: explicitly selected clean main checkout (44897c6) and synthetic workflow CSV.
Outputs: original run under a durable --output root, plus a normalized committed
test fixture at --fixture (directory or deterministic .tar.gz archive). Only
provenance locations and their inventory hashes change in the normalized fixture;
every scientific TSV/YAML byte stays unchanged.
The verify mode runs another checkout's real verifier without monkeypatching.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import Any

REPOSITORY = Path(__file__).resolve().parents[2]
NAMESPACE = Path("/fixture/trait-policy-44897c6")
MATRIX_SOURCE = "src/mtp/taxa_trait_matrix.py"
EXPECTED_MAIN_SHORT_COMMIT = "44897c6"
EXPECTED_MAIN_MATRIX_SIZE = 28555


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve_fixture(root: Path, path: Path) -> Path:
    root = root.resolve()
    if path.is_relative_to(NAMESPACE):
        candidate = root / path.relative_to(NAMESPACE)
    elif path.is_relative_to(root):
        candidate = path
    else:
        raise ValueError(f"Path outside portable fixture: {path}")
    target = candidate.resolve(strict=True)
    if not target.is_relative_to(root):
        raise ValueError("Fixture path escapes its root")
    return target


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_deterministic_tar_gz(source_dir: Path, archive_path: Path) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    files = sorted(
        path.relative_to(source_dir).as_posix() for path in source_dir.rglob("*") if path.is_file()
    )
    tar_buffer = io.BytesIO()
    with tarfile.open(fileobj=tar_buffer, mode="w") as tar:
        for rel_posix in files:
            data = (source_dir / rel_posix).read_bytes()
            info = tarfile.TarInfo(name=rel_posix)
            info.size = len(data)
            info.mtime = 0
            info.mode = 0o644
            info.uid = 0
            info.gid = 0
            info.uname = ""
            info.gname = ""
            tar.addfile(info, io.BytesIO(data))
    with (
        archive_path.open("wb") as raw_handle,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as gz_handle,
    ):
        gz_handle.write(tar_buffer.getvalue())


def extract_fixture_bundle(fixture: Path, destination: Path) -> None:
    if fixture.is_dir():
        shutil.copytree(fixture, destination)
        return
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fixture, mode="r:gz") as tar:
        tar.extractall(destination, filter="data")


def _inputs(output: Path) -> dict[str, Path]:
    inputs = output / "inputs"
    inputs.mkdir()
    shutil.copyfile(REPOSITORY / "tests/fixtures/workflow/full_dump.csv", inputs / "full_dump.csv")
    base = REPOSITORY / "tests/fixtures/workflow"
    edges = (base / "merged-kg_edges.tsv").read_text(encoding="utf-8").splitlines()
    # Use one sole carrier per raw term: direct608, direct612, strain605.
    edges = [row for row in edges if not row.startswith("annotation")]
    edges.extend(
        [
            "direct608\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000608\tfixture:annotation",
            "direct612\tNCBITaxon:102\tbiolink:has_phenotype\tMETPO:1000612\tBactoTraits",
            "strain_link\tkgmicrobe.strain:fixture\tbiolink:subclass_of\tNCBITaxon:103\tBacDive",
            "strain605\tkgmicrobe.strain:fixture\tbiolink:has_phenotype\tMETPO:1000605\tBacDive",
        ]
    )
    (inputs / "merged-kg_edges.tsv").write_text("\n".join(edges) + "\n", encoding="utf-8")
    nodes = (
        (base / "merged-kg_nodes.tsv")
        .read_text(encoding="utf-8")
        .replace(
            "METPO:1000608\tbiolink:PhenotypicQuality\tfacultatively anaerobic\n",
            "METPO:1000608\tbiolink:PhenotypicQuality\tfacultatively aerobic\n",
        )
    )
    nodes += (
        "METPO:1000605\tbiolink:PhenotypicQuality\tfacultatively anaerobic\n"
        "METPO:1000612\tbiolink:PhenotypicQuality\tfacultative\n"
        "kgmicrobe.strain:fixture\tbiolink:OrganismTaxon\tSynthetic BacDive strain\n"
    )
    (inputs / "merged-kg_nodes.tsv").write_text(nodes, encoding="utf-8")
    return {
        "raw_dump": inputs / "full_dump.csv",
        "kg_edges": inputs / "merged-kg_edges.tsv",
        "kg_nodes": inputs / "merged-kg_nodes.tsv",
    }


def generate(code: Path, output: Path, fixture: Path) -> None:
    from git import Repo

    from mtp.workflow.operations import OPERATIONS, validate_execution_product
    from mtp.workflow.runtime import execute_workflow, inventory
    from mtp.workflow.verification import verify_core_run

    if output.exists() or fixture.exists():
        raise FileExistsError("Generation output and committed fixture must both be fresh")
    with Repo(code) as repository:
        if repository.is_dirty(untracked_files=True):
            raise ValueError("Historical producer checkout must be clean")
        producer_commit = repository.head.commit.hexsha
    if not producer_commit.startswith(EXPECTED_MAIN_SHORT_COMMIT):
        raise ValueError(
            f"Generator requires the audited main producer commit {EXPECTED_MAIN_SHORT_COMMIT}"
        )
    matrix_file = code / MATRIX_SOURCE
    expected_matrix_sha256 = digest(REPOSITORY / MATRIX_SOURCE)
    if (
        digest(matrix_file) != expected_matrix_sha256
        or matrix_file.stat().st_size != EXPECTED_MAIN_MATRIX_SIZE
    ):
        raise ValueError("Producer matrix bytes differ from audited main source")
    output.mkdir(parents=True)
    inputs = _inputs(output)
    run = execute_workflow(
        code / "workflows/bugsigdb_core.yaml",
        inputs,
        output / "runs",
        OPERATIONS,
        parameter_overrides={
            "source_mode": "fixture",
            "strain_rollup": "phenotype",
            "trait_key_policy": "facultative_union_v1",
        },
        validator=validate_execution_product,
        repository_root=code,
    )
    original = Path(run.manifest_path)
    verify_core_run(original)
    original_manifest_hash = digest(original)
    before = {
        path.relative_to(original.parent).as_posix(): path.read_bytes()
        for path in original.parent.rglob("*")
        if path.is_file()
    }
    is_archive = fixture.name.endswith(".tar.gz")
    staging = output / "normalized_fixture" if is_archive else fixture
    staging.mkdir(parents=True)
    shutil.copytree(original.parent, staging / "core")
    shutil.copytree(output / "inputs", staging / "inputs")
    transformations: list[dict[str, str]] = []

    def relocate(field: str, replacement: str) -> str:
        transformations.append(
            {
                "field": field,
                "replacement": replacement,
            }
        )
        return replacement

    run.manifest_path = relocate("run.manifest_path", str(NAMESPACE / "core/run.json"))
    run.spec_path = relocate(
        "run.spec_path", str(NAMESPACE / "producer/workflows/bugsigdb_core.yaml")
    )
    runtime = json.loads(run.runtime_json)
    if "executable" in runtime:
        runtime["executable"] = relocate("runtime.executable", str(NAMESPACE / "runtime/python"))
    run.runtime_json = json.dumps(runtime, sort_keys=True)
    for artifact in run.artifacts or []:
        if artifact.product_id in inputs:
            replacement = str(NAMESPACE / "inputs" / inputs[artifact.product_id].name)
            artifact.path = relocate(f"artifact.{artifact.product_id}.path", replacement)
            if artifact.requested_path:
                artifact.requested_path = relocate(
                    f"artifact.{artifact.product_id}.requested_path",
                    replacement,
                )
    source_artifact = next(
        item for item in run.artifacts or [] if item.product_id == "source_identity"
    )
    source_path = staging / "core" / source_artifact.path
    source = json.loads(source_path.read_bytes())
    for name, item in source["consumed_files"].items():
        item["path"] = relocate(
            f"source_identity.{name}.path",
            str(NAMESPACE / "inputs" / inputs[name].name),
        )
    _write_json(source_path, source)
    old_source_hash = source_artifact.sha256
    source_artifact.sha256, source_artifact.size_bytes, source_artifact.inventory = inventory(
        source_path, "file"
    )
    normalized = staging / "core/run.json"
    normalized.write_text(run.model_dump_json(indent=2) + "\n", encoding="utf-8")
    verify_core_run(
        NAMESPACE / "core/run.json", resolver=lambda path: resolve_fixture(staging, path)
    )
    allowed = {"run.json", source_artifact.path}
    after = {
        path.relative_to(staging / "core").as_posix(): path.read_bytes()
        for path in (staging / "core").rglob("*")
        if path.is_file()
    }
    changed = {name for name in before if before[name] != after[name]}
    if changed != allowed or before.keys() != after.keys():
        raise ValueError(f"Normalization changed unexpected files: {sorted(changed)}")
    if digest(original) != original_manifest_hash:
        raise ValueError("Original generation manifest changed")
    original_after = {
        path.relative_to(original.parent).as_posix(): path.read_bytes()
        for path in original.parent.rglob("*")
        if path.is_file()
    }
    if original_after != before:
        raise ValueError("Original generation evidence changed during normalization")
    _write_json(
        staging / "generation.json",
        {
            "issue": 1231,
            "fixture_kind": "location-normalized copy, not byte-identical execution record",
            "producer_commit": producer_commit,
            "producer_matrix_inventory": {
                "path": MATRIX_SOURCE,
                "sha256": expected_matrix_sha256,
                "size_bytes": EXPECTED_MAIN_MATRIX_SIZE,
            },
            "original_manifest_sha256": original_manifest_hash,
            "normalized_manifest_sha256": digest(normalized),
            "original_generation_verified": True,
            "normalized_historical_verification": True,
            "location_transformations": transformations,
            "source_identity_inventory_rebinding": {
                "before_sha256": old_source_hash,
                "after_sha256": source_artifact.sha256,
            },
            "unchanged_scientific_files": sorted(set(before) - allowed),
        },
    )
    if is_archive:
        _write_deterministic_tar_gz(staging, fixture)
    print(f"Verified original union run and normalized fixture; {len(before)} retained files")


def verify(code: Path, fixture: Path, output: Path, expected_error: str | None) -> int:
    from mtp.workflow.verification import verify_core_run

    result: dict[str, Any] = {"verifier_checkout_name": code.name}
    with tempfile.TemporaryDirectory() as tmp:
        bundle_root = fixture if fixture.is_dir() else Path(tmp) / "bundle"
        if not fixture.is_dir():
            extract_fixture_bundle(fixture, bundle_root)
        try:
            run = verify_core_run(
                NAMESPACE / "core/run.json",
                resolver=lambda path: resolve_fixture(bundle_root, path),
            )
        except ValueError as error:
            result.update(status="rejected", error=str(error))
            success = expected_error is not None and expected_error in str(error)
        else:
            result.update(status="verified", producer_commit=run.git_commit)
            success = expected_error is None
    if output.exists():
        raise FileExistsError(output)
    _write_json(output, result)
    print(json.dumps(result, sort_keys=True))
    return int(not success)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("generate", "verify"))
    parser.add_argument("--code", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-error")
    args = parser.parse_args()
    os.environ.setdefault("GIT_CONFIG_GLOBAL", "/dev/null")
    os.environ.setdefault("GIT_CONFIG_NOSYSTEM", "1")
    code = args.code.resolve(strict=True)
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(code / "src"))
    if args.mode == "generate":
        generate(code, args.output.absolute(), args.fixture.absolute())
        return 0
    return verify(code, args.fixture.absolute(), args.output.absolute(), args.expected_error)


if __name__ == "__main__":
    raise SystemExit(main())
