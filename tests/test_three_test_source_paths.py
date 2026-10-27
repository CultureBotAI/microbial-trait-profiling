"""Read-only source-location protection, including symlink hops (#1078, 2026-10-26)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pytest

from mtp.datamodel.workflow_pydantic import (
    ArmRun,
    ArtifactInstance,
    ComparisonRun,
    ProductDefinition,
    StorageKind,
    WorkflowDefinition,
    WorkflowRun,
)
from mtp.multiplicity.sources import VerifiedSources
from mtp.three_test_table.source_paths import (
    declared_source_roots,
    manifest_roots,
    source_roots,
)


def file(path: Path, text: str = "{}\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def comparison(tmp_path: Path, core_record: dict[str, Any]) -> tuple[Path, Path]:
    core = file(tmp_path / "core" / "run.json", json.dumps(core_record))
    path = tmp_path / "comparison" / "comparison_run.json"
    record = {
        "manifest_path": str(path),
        "core_manifest_path": str(core),
        **{
            f"{arm}_manifest_path": str(file(tmp_path / arm / "arm_run.json"))
            for arm in ("noel", "jed", "ilyas")
        },
    }
    file(path, json.dumps(record))
    return path, core


def file_definition(product: str) -> dict[str, Any]:
    return WorkflowDefinition(
        id="synthetic-core",
        version="1",
        steps=[],
        products=[
            ProductDefinition(
                id=product,
                kind=product,
                storage_kind=StorageKind.file,
                format="tsv",
                filename=f"{product}.tsv",
            )
        ],
    ).model_dump(mode="json")


def verified(
    core_path: Path,
    artifacts: list[ArtifactInstance],
    *,
    storage: StorageKind = StorageKind.file,
) -> VerifiedSources:
    root = core_path.parent.parent
    manifests = {name: root / name / "arm_run.json" for name in ("noel", "jed", "ilyas")}
    record = ComparisonRun.model_construct(
        manifest_path=str(root / "comparison" / "comparison_run.json"),
        core_manifest_path=str(core_path),
        noel_manifest_path=str(manifests["noel"]),
        jed_manifest_path=str(manifests["jed"]),
        ilyas_manifest_path=str(manifests["ilyas"]),
    )
    return VerifiedSources(
        comparison=record,
        core=WorkflowRun.model_construct(
            manifest_path=str(core_path),
            artifacts=artifacts,
            definition=WorkflowDefinition(
                id="synthetic-core",
                version="1",
                steps=[],
                products=[
                    ProductDefinition(
                        id=artifact.product_id,
                        kind=artifact.product_id,
                        storage_kind=storage,
                        format="tsv",
                        filename=Path(artifact.path).name,
                    )
                    for artifact in artifacts
                ],
            ),
        ),
        arms={
            name: ArmRun.model_construct(manifest_path=str(path))
            for name, path in manifests.items()
        },
        bases={},
        features={},
        scope="off",
    )


def test_file_alias_chain_retains_every_middle_parent(tmp_path: Path) -> None:
    target = file(tmp_path / "target" / "run.json")
    middle = tmp_path / "middle" / "relay.json"
    alias = tmp_path / "alias" / "run.json"
    middle.parent.mkdir()
    alias.parent.mkdir()
    middle.symlink_to(Path("..") / "target" / target.name)
    alias.symlink_to(Path("..") / "middle" / middle.name)
    roots = manifest_roots([alias], multiplicity=None, label_support=None)
    assert {alias.parent, middle.parent, target.parent}.issubset(roots)


def test_component_link_does_not_hide_a_middle_file_alias(tmp_path: Path) -> None:
    target = file(tmp_path / "target" / "run.json")
    physical_middle = tmp_path / "physical-middle" / "relay.json"
    physical_middle.parent.mkdir()
    physical_middle.symlink_to(target)
    component = tmp_path / "bridge" / "route"
    component.parent.mkdir()
    component.symlink_to(physical_middle.parent, target_is_directory=True)
    alias = tmp_path / "alias" / "run.json"
    alias.parent.mkdir()
    alias.symlink_to(component / physical_middle.name)
    roots = manifest_roots([alias], multiplicity=None, label_support=None)
    assert {alias.parent, component, physical_middle.parent, target.parent}.issubset(roots)


@pytest.mark.parametrize("option", ["multiplicity", "label_support"])
def test_directory_alias_chain_protects_only_its_trees(tmp_path: Path, option: str) -> None:
    target = tmp_path / "target" / "evidence"
    target.mkdir(parents=True)
    middle = tmp_path / "middle" / "evidence"
    middle.parent.mkdir()
    middle.symlink_to(target, target_is_directory=True)
    alias = tmp_path / "alias" / "evidence"
    alias.parent.mkdir()
    alias.symlink_to(middle, target_is_directory=True)
    roots = manifest_roots(
        [],
        multiplicity=alias if option == "multiplicity" else None,
        label_support=alias if option == "label_support" else None,
    )
    assert set(roots) == {alias, middle, target}
    for sibling in (alias.parent / "output", middle.parent / "output", target.parent / "output"):
        assert not any(
            sibling.resolve().is_relative_to(root.resolve())
            or root.resolve().is_relative_to(sibling.resolve())
            for root in roots
        )


@pytest.mark.parametrize("kind", ["file", "component", "growing"])
def test_symlink_cycles_fail_closed(tmp_path: Path, kind: str) -> None:
    first, second = tmp_path / "first", tmp_path / "second"
    if kind == "growing":
        first.symlink_to("first/more")
    else:
        first.symlink_to(second, target_is_directory=kind == "component")
        second.symlink_to(first, target_is_directory=kind == "component")
    candidate = first if kind == "file" else first / "run.json"
    with pytest.raises(ValueError, match=r"[Ss]ymlink"):
        manifest_roots([candidate], multiplicity=None, label_support=None)


def test_symlink_then_parent_component_preserves_kernel_path_meaning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = file(tmp_path / "target" / "run.json")
    nested = target.parent / "nested"
    nested.mkdir()
    route = tmp_path / "route"
    route.mkdir()
    (route / "link").symlink_to(nested, target_is_directory=True)
    monkeypatch.chdir(tmp_path)
    supplied = Path("route/link/../run.json")
    roots = manifest_roots([supplied], multiplicity=None, label_support=None)
    assert target.parent in roots
    assert tmp_path / supplied.parent in roots
    assert route not in roots


def test_dangling_file_alias_still_protects_both_locations(tmp_path: Path) -> None:
    target = tmp_path / "future" / "run.json"
    alias = tmp_path / "alias" / "run.json"
    alias.parent.mkdir()
    alias.symlink_to(target)
    roots = manifest_roots([alias], multiplicity=None, label_support=None)
    assert {alias.parent, target.parent}.issubset(roots)
    assert not target.parent.exists()


def test_preflight_protects_self_declared_core_location_and_relative_artifacts(
    tmp_path: Path,
) -> None:
    self_declared = tmp_path / "original-core" / "run.json"
    requested = file(tmp_path / "requested" / "edges.tsv")
    path, core = comparison(
        tmp_path,
        {
            "manifest_path": str(self_declared),
            "definition": file_definition("kg_edges"),
            "artifacts": [
                {
                    "product_id": "kg_edges",
                    "path": "relative/edges.tsv",
                    "requested_path": str(requested),
                }
            ],
        },
    )
    roots = declared_source_roots(path, multiplicity=None, label_support=None)
    assert self_declared.parent in roots
    assert core.parent / "relative" in roots
    assert requested.parent in roots
    assert not self_declared.parent.exists()
    assert not (core.parent / "relative").exists()


def test_preflight_relative_artifact_uses_resolved_core_manifest_parent(tmp_path: Path) -> None:
    path, core = comparison(
        tmp_path,
        {
            "definition": file_definition("kg_edges"),
            "artifacts": [{"product_id": "kg_edges", "path": "products/edges.tsv"}],
        },
    )
    alias = tmp_path / "core-alias" / "run.json"
    alias.parent.mkdir()
    alias.symlink_to(core)
    record = json.loads(path.read_text())
    record["core_manifest_path"] = str(alias)
    path.write_text(json.dumps(record), encoding="utf-8")
    roots = declared_source_roots(path, multiplicity=None, label_support=None)
    assert core.parent / "products" in roots
    assert alias.parent in roots


def test_hard_linked_core_retains_its_self_declared_original_parent(tmp_path: Path) -> None:
    original = tmp_path / "original-core" / "run.json"
    path, supplied = comparison(tmp_path, {"manifest_path": str(original)})
    original.parent.mkdir()
    original.hardlink_to(supplied)
    assert original.samefile(supplied)
    roots = declared_source_roots(path, multiplicity=None, label_support=None)
    assert original.parent in roots
    assert supplied.parent in roots
    assert not (original.parent / "must-not-be-created").exists()


@pytest.mark.parametrize("product", ["kg_edges", "kg_nodes", "reported_taxa"])
def test_both_helpers_protect_every_external_core_file_and_requested_alias(
    tmp_path: Path, product: str
) -> None:
    target = file(tmp_path / "external" / product / "data.tsv")
    middle = tmp_path / "middle" / product / "data.tsv"
    middle.parent.mkdir(parents=True)
    middle.symlink_to(target)
    requested = tmp_path / "requested" / product / "data.tsv"
    requested.parent.mkdir(parents=True)
    requested.symlink_to(middle)
    declaration = {"product_id": product, "path": str(target), "requested_path": str(requested)}
    path, core = comparison(
        tmp_path,
        {
            "definition": file_definition(product),
            "artifacts": [declaration],
        },
    )
    sources = verified(
        core,
        [
            ArtifactInstance.model_construct(
                product_id=product, path=str(target), requested_path=str(requested)
            )
        ],
    )
    for roots in (
        declared_source_roots(path, multiplicity=None, label_support=None),
        source_roots(sources, multiplicity=None, label_support=None),
    ):
        assert {requested.parent, middle.parent, target.parent}.issubset(roots)


def test_verified_relative_artifact_resolves_against_its_core_bundle(tmp_path: Path) -> None:
    core = tmp_path / "core" / "run.json"
    sources = verified(
        core, [ArtifactInstance.model_construct(product_id="kg_edges", path="products/edges.tsv")]
    )
    roots = source_roots(sources, multiplicity=None, label_support=None)
    assert core.parent / "products" in roots


def test_synthetic_core_without_artifacts_is_only_a_location_preflight(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, core = comparison(tmp_path, {"synthetic": True})
    before = {entry: entry.read_bytes() for entry in tmp_path.rglob("*.json")}
    forbidden = Mock(side_effect=AssertionError("Read-only location helper attempted mkdir"))
    monkeypatch.setattr(Path, "mkdir", forbidden)
    roots = declared_source_roots(path, multiplicity=None, label_support=None)
    forbidden.assert_not_called()
    assert path.parent in roots
    assert core.parent in roots
    assert {entry: entry.read_bytes() for entry in tmp_path.rglob("*.json")} == before


def test_missing_core_retains_declared_roots_and_defers_execution_validation(
    tmp_path: Path,
) -> None:
    path, core = comparison(tmp_path, {})
    core.unlink()
    roots = declared_source_roots(path, multiplicity=None, label_support=None)
    assert core.parent in roots
    assert path.parent in roots
    assert not core.exists()


@pytest.mark.parametrize(
    "artifacts",
    [
        "invalid",
        ["invalid"],
        [{"product_id": "kg_edges", "path": " "}],
        [{"product_id": "kg_edges", "path": "data.tsv", "requested_path": 7}],
    ],
)
def test_malformed_location_declarations_fail_closed(tmp_path: Path, artifacts: Any) -> None:
    path, _ = comparison(
        tmp_path,
        {
            "definition": file_definition("kg_edges"),
            "artifacts": artifacts,
        },
    )
    with pytest.raises(ValueError, match=r"[Ll]ocation|[Aa]rtifact|[Pp]ath"):
        declared_source_roots(path, multiplicity=None, label_support=None)


@pytest.mark.parametrize("reader", ["declared", "verified"])
@pytest.mark.parametrize("actual", ["directory", "missing"])
def test_declared_file_kind_protects_every_alias_parent_when_live_type_changes(
    tmp_path: Path,
    reader: str,
    actual: str,
) -> None:
    target = tmp_path / "external-file" / "edges.tsv"
    if actual == "directory":
        target.mkdir(parents=True)
    middle = tmp_path / "middle-file" / "edges.tsv"
    requested = tmp_path / "requested-file" / "edges.tsv"
    middle.parent.mkdir()
    requested.parent.mkdir()
    middle.symlink_to(target)
    requested.symlink_to(middle)
    artifact = ArtifactInstance.model_construct(
        product_id="kg_edges",
        path=str(target),
        requested_path=str(requested),
    )
    path, core = comparison(tmp_path, {})
    sources = verified(core, [artifact])
    core.write_text(
        json.dumps(
            {
                "definition": sources.core.definition.model_dump(mode="json"),
                "artifacts": [artifact.model_dump(mode="json")],
            }
        ),
        encoding="utf-8",
    )
    roots = (
        declared_source_roots(path, multiplicity=None, label_support=None)
        if reader == "declared"
        else source_roots(sources, multiplicity=None, label_support=None)
    )
    assert {target.parent, middle.parent, requested.parent}.issubset(roots)
    assert target not in roots
    if actual == "missing":
        assert not target.parent.exists()


@pytest.mark.parametrize("reader", ["declared", "verified"])
@pytest.mark.parametrize("actual", ["directory", "missing", "file"])
def test_declared_collection_kind_protects_only_its_tree_despite_live_type(
    tmp_path: Path,
    reader: str,
    actual: str,
) -> None:
    target = tmp_path / "external-collection" / "studies"
    if actual == "directory":
        target.mkdir(parents=True)
    elif actual == "file":
        file(target)
    requested = tmp_path / "requested-collection" / "studies"
    requested.parent.mkdir()
    requested.symlink_to(target)
    artifact = ArtifactInstance.model_construct(
        product_id="studies",
        path=str(target),
        requested_path=str(requested),
    )
    path, core = comparison(tmp_path, {})
    sources = verified(core, [artifact], storage=StorageKind.collection)
    core.write_text(
        json.dumps(
            {
                "definition": sources.core.definition.model_dump(mode="json"),
                "artifacts": [artifact.model_dump(mode="json")],
            }
        ),
        encoding="utf-8",
    )
    roots = (
        declared_source_roots(path, multiplicity=None, label_support=None)
        if reader == "declared"
        else source_roots(sources, multiplicity=None, label_support=None)
    )
    assert {target, requested}.issubset(roots)
    assert target.parent not in roots
    assert requested.parent not in roots


@pytest.mark.parametrize(
    "definition",
    [
        None,
        "bad",
        {},
        {"products": "bad"},
        {"products": ["bad"]},
        {"products": [{"id": "kg_edges"}]},
        {"products": [{"id": "kg_edges", "storage_kind": "unknown"}]},
        {"products": [{"id": "kg_edges", "storage_kind": "file"}] * 2},
        {"products": [{"id": "other", "storage_kind": "file"}]},
    ],
)
def test_present_malformed_storage_declarations_fail_closed(
    tmp_path: Path,
    definition: Any,
) -> None:
    path, _ = comparison(
        tmp_path,
        {
            "definition": definition,
            "artifacts": [{"product_id": "kg_edges", "path": "external/edges.tsv"}],
        },
    )
    with pytest.raises(ValueError, match=r"(?i)definition|storage|product"):
        declared_source_roots(path, multiplicity=None, label_support=None)


def test_artifacts_without_storage_declarations_fail_closed(tmp_path: Path) -> None:
    target = tmp_path / "external" / "edges.tsv"
    target.mkdir(parents=True)
    path, _ = comparison(
        tmp_path,
        {
            "artifacts": [{"product_id": "kg_edges", "path": str(target)}],
        },
    )
    with pytest.raises(ValueError, match="storage"):
        declared_source_roots(path, multiplicity=None, label_support=None)
