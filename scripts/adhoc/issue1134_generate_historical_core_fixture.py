"""Generate and verify a portable historical core fixture (#1134, 2026-10-28).

Inputs: explicitly selected clean historical checkout and synthetic workflow CSV.
Outputs: original run under a durable --output root, plus a normalized committed
test fixture at --fixture. Only provenance locations and their inventory hashes
change in the copied fixture; every scientific TSV/YAML byte stays unchanged.
The verify mode runs another checkout's real verifier without monkeypatching.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path
from typing import Any, cast

from git import Repo

REPOSITORY = Path(__file__).resolve().parents[2]
NAMESPACE = Path("/fixture/trait-policy-511c891")
MATRIX_SOURCE = "src/mtp/taxa_trait_matrix.py"


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
    nodes = (base / "merged-kg_nodes.tsv").read_text(encoding="utf-8")
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
    from mtp.workflow.operations import OPERATIONS, validate_core_product
    from mtp.workflow.runtime import Validator, execute_workflow, inventory
    from mtp.workflow.verification import verify_core_run

    if output.exists() or fixture.exists():
        raise FileExistsError("Generation output and committed fixture must both be fresh")
    with Repo(code) as repository:
        if repository.is_dirty(untracked_files=True):
            raise ValueError("Historical producer checkout must be clean")
        producer_commit = repository.head.commit.hexsha
    audited = json.loads(
        (REPOSITORY / "tests/fixtures/trait_key_policy/historical_implementations.json").read_text()
    )
    identity = next(row for row in audited["records"] if row["policy"] == "raw_object_v1")
    if producer_commit != identity["source_commit"]:
        raise ValueError("Generator requires the audited pre-union producer commit")
    if digest(code / MATRIX_SOURCE) != identity["code_inventory"]["sha256"]:
        raise ValueError("Historical producer matrix bytes differ from audited source")
    output.mkdir(parents=True)
    inputs = _inputs(output)
    run = execute_workflow(
        code / "workflows/bugsigdb_core.yaml",
        inputs,
        output / "runs",
        OPERATIONS,
        parameter_overrides={"source_mode": "fixture", "strain_rollup": "phenotype"},
        # The authenticated historical checkout implements this older protocol;
        # the current checkout's validator additionally requires policy keywords.
        validator=cast(Validator, validate_core_product),
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
    fixture.mkdir(parents=True)
    shutil.copytree(original.parent, fixture / "core")
    shutil.copytree(output / "inputs", fixture / "inputs")
    transformations: list[dict[str, str]] = []

    def relocate(field: str, value: str, replacement: str) -> str:
        transformations.append(
            {
                "field": field,
                "original_value_sha256": hashlib.sha256(value.encode()).hexdigest(),
                "replacement": replacement,
            }
        )
        return replacement

    run.manifest_path = relocate(
        "run.manifest_path", run.manifest_path, str(NAMESPACE / "core/run.json")
    )
    run.spec_path = relocate(
        "run.spec_path", run.spec_path, str(NAMESPACE / "producer/workflows/bugsigdb_core.yaml")
    )
    runtime = json.loads(run.runtime_json)
    if "executable" in runtime:
        runtime["executable"] = relocate(
            "runtime.executable", runtime["executable"], str(NAMESPACE / "runtime/python")
        )
    run.runtime_json = json.dumps(runtime, sort_keys=True)
    for artifact in run.artifacts or []:
        if artifact.product_id in inputs:
            replacement = str(NAMESPACE / "inputs" / inputs[artifact.product_id].name)
            artifact.path = relocate(
                f"artifact.{artifact.product_id}.path", artifact.path, replacement
            )
            if artifact.requested_path:
                artifact.requested_path = relocate(
                    f"artifact.{artifact.product_id}.requested_path",
                    artifact.requested_path,
                    replacement,
                )
    source_artifact = next(
        item for item in run.artifacts or [] if item.product_id == "source_identity"
    )
    source_path = fixture / "core" / source_artifact.path
    source = json.loads(source_path.read_bytes())
    for name, item in source["consumed_files"].items():
        item["path"] = relocate(
            f"source_identity.{name}.path",
            item["path"],
            str(NAMESPACE / "inputs" / inputs[name].name),
        )
    _write_json(source_path, source)
    old_source_hash = source_artifact.sha256
    source_artifact.sha256, source_artifact.size_bytes, source_artifact.inventory = inventory(
        source_path, "file"
    )
    normalized = fixture / "core/run.json"
    normalized.write_text(run.model_dump_json(indent=2) + "\n", encoding="utf-8")
    verify_core_run(
        NAMESPACE / "core/run.json", resolver=lambda path: resolve_fixture(fixture, path)
    )
    allowed = {"run.json", source_artifact.path}
    after = {
        path.relative_to(fixture / "core").as_posix(): path.read_bytes()
        for path in (fixture / "core").rglob("*")
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
        fixture / "generation.json",
        {
            "issue": 1134,
            "fixture_kind": "location-normalized copy, not byte-identical execution record",
            "producer_commit": producer_commit,
            "producer_matrix_inventory": identity["code_inventory"],
            "helper_sha256": digest(Path(__file__)),
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
    print(f"Verified original historical run and normalized fixture; {len(before)} retained files")


def verify(code: Path, fixture: Path, output: Path, expected_error: str | None) -> int:
    from mtp.workflow.verification import verify_core_run

    result: dict[str, Any] = {"verifier_checkout_name": code.name}
    try:
        run = verify_core_run(
            NAMESPACE / "core/run.json", resolver=lambda path: resolve_fixture(fixture, path)
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
    code = args.code.resolve(strict=True)
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(code / "src"))
    if args.mode == "generate":
        generate(code, args.output.absolute(), args.fixture.absolute())
        return 0
    return verify(code, args.fixture.absolute(), args.output.absolute(), args.expected_error)


if __name__ == "__main__":
    raise SystemExit(main())
