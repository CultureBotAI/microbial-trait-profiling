"""Fail-before-write and once-only source contracts (#1078; 2026-10-26)."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import Mock

import pytest

from mtp.datamodel.workflow_pydantic import ArtifactInstance, ProductDefinition, StorageKind
from mtp.multiplicity import sources as source_module
from mtp.multiplicity import stage as multiplicity_stage
from mtp.multiplicity.sources import VerifiedSources
from mtp.paths import file_sha256
from mtp.three_test_table import stage
from mtp.workflow import output_roots
from mtp.workflow.runtime import inventory
from test_three_test_stage import sources_for_stage
from test_three_test_stage import tmp_path as tmp_path


def add_external_inputs(
    sources: VerifiedSources, root: Path, monkeypatch: pytest.MonkeyPatch
) -> dict[str, Path]:
    paths = {}
    assert sources.core.artifacts is not None
    original = {
        artifact.product_id: sources.product(artifact.product_id)
        for artifact in sources.core.artifacts
    }
    core_path = Path(sources.core.manifest_path)
    declared = json.loads(core_path.read_text())
    artifacts = declared.setdefault("artifacts", [])
    for name in ("kg_edges", "kg_nodes", "reported_taxa"):
        path = root / name / f"{name}.tsv"
        path.parent.mkdir(parents=True)
        path.write_text("synthetic bound source\n", encoding="utf-8")
        paths[name] = path
        digest, size, entries = inventory(path, "file")
        item = ArtifactInstance(
            id=f"fixture-{name}",
            run_id=sources.core.run_id,
            product_id=name,
            path=str(path),
            requested_path=str(path),
            sha256=digest,
            size_bytes=size,
            inventory=entries,
            validations=[],
        )
        sources.core.artifacts.append(item)
        artifacts.append(item.model_dump())
        sources.core.definition.products.append(
            ProductDefinition(
                id=name,
                kind=name,
                storage_kind=StorageKind.file,
                format="tsv",
                filename=path.name,
                external=True,
            )
        )
    declared["definition"] = sources.core.definition.model_dump(mode="json")
    core_path.write_text(json.dumps(declared), encoding="utf-8")
    sources.comparison.core_manifest_sha256 = file_sha256(core_path)
    for arm in sources.arms.values():
        arm.core_manifest_sha256 = file_sha256(core_path)
    Path(sources.comparison.manifest_path).write_text(sources.comparison.model_dump_json())
    monkeypatch.setattr(
        source_module,
        "artifact_path",
        lambda _run, name: paths[name] if name in paths else original[name],
    )
    return paths


def test_shared_temporary_root_guard_is_authoritative(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    shared = Mock(return_value=tmp_path / "new-shared-temporary-root")
    monkeypatch.setattr(output_roots, "temporary_root", shared)
    with pytest.raises(ValueError, match="durable"):
        stage.require_durable(tmp_path / "output")
    shared.assert_called_once()


def test_failed_verification_does_not_create_output_parents(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    output = tmp_path / "not-created" / "nested" / "output"
    monkeypatch.setattr(
        stage, "open_verified_sources", Mock(side_effect=ValueError("invalid source evidence"))
    )
    with pytest.raises(ValueError, match="invalid source evidence"):
        stage.run_table(
            Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
        )
    assert not (tmp_path / "not-created").exists()


def test_verify_checks_declared_roots_before_first_source_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    output = tmp_path / "completed"
    output.mkdir()
    (output / stage.MANIFEST).write_text("{}")
    forbidden = Mock(side_effect=AssertionError("source verification ran before preflight"))
    monkeypatch.setattr(stage, "open_verified_sources", forbidden)
    scratch = Path(sources.core.manifest_path).parent / "new" / "scratch"
    with pytest.raises(ValueError, match="disjoint"):
        stage.verify_table(output, Path(sources.comparison.manifest_path), scratch=scratch)
    forbidden.assert_not_called()
    assert not scratch.parent.exists()


def test_verify_opens_sources_once_for_the_complete_replay(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    output = tmp_path / "completed"
    record = stage.run_table(comparison, output, scratch=tmp_path / "initial-scratch")
    checked = Mock(return_value=sources)
    monkeypatch.setattr(stage, "open_verified_sources", checked)
    assert stage.verify_table(output, comparison, scratch=tmp_path / "replay-scratch") == record
    checked.assert_called_once()


@pytest.mark.parametrize("destination", ["output", "scratch"])
def test_middle_file_symlink_hop_is_protected_before_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, destination: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    original = Path(sources.comparison.manifest_path)
    middle = tmp_path / "middle" / "comparison.json"
    alias = tmp_path / "alias" / "comparison.json"
    middle.parent.mkdir()
    alias.parent.mkdir()
    middle.symlink_to(original)
    alias.symlink_to(middle)
    nested = middle.parent / "new" / destination
    output = nested if destination == "output" else tmp_path / "output"
    scratch = nested if destination == "scratch" else tmp_path / "scratch"
    forbidden = Mock(side_effect=AssertionError("source verification ran before alias preflight"))
    monkeypatch.setattr(stage, "open_verified_sources", forbidden)
    with pytest.raises(ValueError, match="disjoint"):
        stage.run_table(alias, output, scratch=scratch)
    forbidden.assert_not_called()
    assert not nested.parent.exists()


@pytest.mark.parametrize("role", ["kg_edges", "kg_nodes", "reported_taxa"])
@pytest.mark.parametrize("destination", ["output", "scratch"])
def test_external_core_file_parent_is_protected_before_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str, destination: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    paths = add_external_inputs(sources, tmp_path / "external", monkeypatch)
    nested = paths[role].parent / "new" / destination
    output = nested if destination == "output" else tmp_path / "output"
    scratch = nested if destination == "scratch" else tmp_path / "scratch"
    forbidden = Mock(
        side_effect=AssertionError("source verification ran before external preflight")
    )
    monkeypatch.setattr(stage, "open_verified_sources", forbidden)
    with pytest.raises(ValueError, match="disjoint"):
        stage.run_table(Path(sources.comparison.manifest_path), output, scratch=scratch)
    forbidden.assert_not_called()
    assert not nested.parent.exists()


def test_all_consumed_core_files_are_explicitly_hash_bound(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    paths = add_external_inputs(sources, tmp_path / "external", monkeypatch)
    record = stage.run_table(
        Path(sources.comparison.manifest_path), tmp_path / "output", scratch=tmp_path / "scratch"
    )
    for path in paths.values():
        digest = file_sha256(path)
        assert record["inputs"][f"{path.name}@sha256:{digest}"] == {
            "sha256": digest,
            "size_bytes": path.stat().st_size,
        }


def test_missing_storage_declarations_reject_before_probe_or_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    path = add_external_inputs(sources, tmp_path / "external", monkeypatch)["kg_edges"]
    path.unlink()
    path.mkdir()
    manifest = Path(sources.core.manifest_path)
    document = json.loads(manifest.read_text())
    del document["definition"]
    manifest.write_text(json.dumps(document))
    probe = Mock(side_effect=AssertionError("probe ran without source storage declarations"))
    verify = Mock(side_effect=AssertionError("verification ran before source preflight"))
    monkeypatch.setattr(stage, "_probe_publication", probe)
    monkeypatch.setattr(stage, "open_verified_sources", verify)
    output = path.parent / "new" / "output"
    with pytest.raises(ValueError, match="storage"):
        stage.run_table(
            Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
        )
    probe.assert_not_called()
    verify.assert_not_called()
    assert not output.parent.exists()


@pytest.mark.parametrize("role", ["kg_edges", "kg_nodes", "reported_taxa"])
@pytest.mark.parametrize("mutation", ["delete", "directory"])
def test_binding_keeps_declared_file_hash_after_live_type_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str, mutation: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    path = add_external_inputs(sources, tmp_path / "external", monkeypatch)[role]
    expected = file_sha256(path)
    path.unlink()
    if mutation == "directory":
        path.mkdir()
    bound = source_module.bind_recorded_hashes(
        sources, file_sha256(Path(sources.comparison.manifest_path))
    )
    assert bound.recorded_hashes.get(path.resolve()) == expected


def test_declared_collection_remains_outside_regular_file_hash_bindings(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    collection = Path(sources.core.manifest_path).parent / "studies"
    collection.mkdir()
    (collection / "study.yaml").write_text("study_id: synthetic-study\n", encoding="utf-8")
    digest, size, entries = inventory(collection, "collection")
    assert sources.core.artifacts is not None
    original = {
        artifact.product_id: sources.product(artifact.product_id)
        for artifact in sources.core.artifacts
    }
    sources.core.definition.products.append(
        ProductDefinition(
            id="studies",
            kind="study_collection",
            storage_kind=StorageKind.collection,
            format="yaml",
            filename="studies",
        )
    )
    artifact = ArtifactInstance(
        id="synthetic-studies",
        product_id="studies",
        run_id=sources.core.run_id,
        path=str(collection),
        sha256=digest,
        size_bytes=size,
        inventory=entries,
        validations=[],
    )
    sources.core.artifacts.append(artifact)
    monkeypatch.setattr(
        source_module,
        "artifact_path",
        lambda _run, name: collection if name == "studies" else original[name],
    )
    bound = source_module.bind_recorded_hashes(
        sources, file_sha256(Path(sources.comparison.manifest_path))
    )
    assert collection.resolve() not in bound.recorded_hashes
    for path in original.values():
        assert bound.expected_digest(path) == file_sha256(path)


@pytest.mark.parametrize("role", ["kg_edges", "kg_nodes", "reported_taxa"])
@pytest.mark.parametrize("mutation", ["delete", "directory"])
@pytest.mark.parametrize("timing", ["after_open", "after_bind"])
def test_declared_core_file_cannot_disappear_from_publication_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str, mutation: str, timing: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    path = add_external_inputs(sources, tmp_path / "external", monkeypatch)[role]
    original_bind = source_module.bind_recorded_hashes

    def change_type() -> None:
        path.unlink()
        if mutation == "directory":
            path.mkdir()

    def opened(_path: Path) -> VerifiedSources:
        if timing == "after_open":
            change_type()
        return sources

    def bound(verified: VerifiedSources, digest: str) -> VerifiedSources:
        result = original_bind(verified, digest)
        if timing == "after_bind":
            change_type()
        return result

    checked = Mock(side_effect=opened)
    monkeypatch.setattr(stage, "open_verified_sources", checked)
    monkeypatch.setattr(stage, "bind_recorded_hashes", bound)
    output = tmp_path / "uncreated" / "output"
    scratch = tmp_path / "scratch"
    with pytest.raises((FileNotFoundError, IsADirectoryError)):
        stage.run_table(Path(sources.comparison.manifest_path), output, scratch=scratch)
    checked.assert_called_once()
    assert not output.parent.exists()
    assert not scratch.exists()


def test_optional_multiplicity_reuses_verified_sources_and_declares_co_presence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    optional = tmp_path / "multiplicity"
    multiplicity_stage.run_post_arm(comparison, optional)
    forbidden = Mock(side_effect=AssertionError("multiplicity repeated source verification"))
    monkeypatch.setattr(multiplicity_stage, "open_verified_sources", forbidden)
    record = stage.run_table(
        comparison,
        tmp_path / "output",
        scratch=tmp_path / "scratch",
        multiplicity=optional,
    )
    forbidden.assert_not_called()
    assert record["multiplicity_binding"] == "verified_co_presence_only"


def test_failed_adaptation_does_not_create_output_parents(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    output = tmp_path / "not-created" / "nested" / "output"
    monkeypatch.setattr(stage, "_merge_tables", Mock(side_effect=ValueError("invalid adapter")))
    with pytest.raises(ValueError, match="invalid adapter"):
        stage.run_table(
            Path(sources.comparison.manifest_path), output, scratch=tmp_path / "scratch"
        )
    assert not (tmp_path / "not-created").exists()


@pytest.mark.parametrize("destination", ["output", "scratch"])
def test_symlink_then_parent_component_keeps_disjoint_siblings_available(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, destination: str
) -> None:
    sources = sources_for_stage(tmp_path, monkeypatch)
    comparison = Path(sources.comparison.manifest_path)
    nested = comparison.parent / "nested"
    nested.mkdir()
    route = tmp_path / "route"
    route.mkdir()
    (route / "link").symlink_to(nested, target_is_directory=True)
    supplied = route / "link" / ".." / comparison.name
    assert supplied.resolve() == comparison.resolve()
    sibling = route / "new" / destination
    output = sibling if destination == "output" else tmp_path / "output"
    scratch = sibling if destination == "scratch" else tmp_path / "scratch"
    checked = Mock(side_effect=ValueError("reached verification with disjoint sibling"))
    monkeypatch.setattr(stage, "open_verified_sources", checked)
    with pytest.raises(ValueError, match="reached verification with disjoint sibling"):
        stage.run_table(supplied, output, scratch=scratch)
    checked.assert_called_once()
    assert not sibling.parent.exists()
    roots = stage._manifest_roots([supplied], multiplicity=None, label_support=None)
    for protected in (comparison.parent / "blocked", route / "link" / ".." / "blocked"):
        with pytest.raises(ValueError, match="disjoint"):
            stage._disjoint(protected, roots)
