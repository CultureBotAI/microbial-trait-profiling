"""Declared-root source references retain byte identity and relocation (#1664)."""

from __future__ import annotations

import importlib
import json
import shutil
from pathlib import Path
from types import ModuleType

import pytest

from mtp.paths import file_sha256


@pytest.fixture
def api() -> ModuleType:
    assert importlib.util.find_spec("mtp.provenance_paths") is not None
    module = importlib.import_module("mtp.provenance_paths")
    assert module.__file__ is not None
    assert Path(module.__file__).resolve() == (
        Path(__file__).resolve().parents[1] / "src/mtp/provenance_paths.py"
    )
    return module


def test_relocated_reference_checks_the_named_bytes(tmp_path: Path, api: ModuleType) -> None:
    root = tmp_path / "original"
    (root / "source").mkdir(parents=True)
    (root / "records").mkdir()
    source = root / "source/input.tsv"
    source.write_text("first\n")
    record = root / "records/provenance.json"
    paths = api.PortablePaths.create(record, [source])
    assert paths.declaration == {"kind": "root-relative-v1", "root": ".."}
    ref = {"path": paths.encode(source), "sha256": file_sha256(source)}
    assert ref["path"] == "source/input.tsv"
    record.write_text(json.dumps({"path_basis": paths.declaration, "source": ref}))
    moved = tmp_path / "moved"
    shutil.copytree(root, moved)
    root.rename(tmp_path / "hidden")
    relocated = moved / "records/provenance.json"
    assert api.verify_reference(relocated, ref, path_basis=paths.declaration) == (
        moved / "source/input.tsv"
    )
    (moved / "source/input.tsv").write_text("other\n")
    with pytest.raises(ValueError, match=r"hash|SHA|digest"):
        api.verify_reference(relocated, ref, path_basis=paths.declaration)


@pytest.mark.parametrize(
    "name", ["../escape", "/outside", "a/../b", "a//b", "./b", "a\\b", "C:/b", ""]
)
def test_noncanonical_references_are_refused(tmp_path: Path, api: ModuleType, name: str) -> None:
    with pytest.raises(ValueError, match=r"relative|canonical|path"):
        api.resolve_reference(
            tmp_path / "record.json", name, path_basis={"kind": "root-relative-v1", "root": "."}
        )


@pytest.mark.parametrize(
    "basis",
    [
        {"kind": "unknown", "root": "."},
        {"kind": "root-relative-v1", "root": "/absolute"},
        {"kind": "root-relative-v1", "root": "../sibling"},
        {"kind": "root-relative-v1", "root": ".", "extra": True},
    ],
)
def test_ambiguous_root_declarations_are_refused(
    tmp_path: Path, api: ModuleType, basis: dict[str, object]
) -> None:
    with pytest.raises(ValueError, match=r"root|basis|declaration"):
        api.resolve_reference(tmp_path / "record.json", "data.tsv", path_basis=basis)


def test_declared_root_rejects_symlink_escape(tmp_path: Path, api: ModuleType) -> None:
    root = tmp_path / "bundle"
    root.mkdir()
    outside = tmp_path / "outside.tsv"
    outside.write_text("bytes\n")
    (root / "alias.tsv").symlink_to(outside)
    with pytest.raises(ValueError, match=r"root|escape"):
        api.resolve_reference(
            root / "record.json", "alias.tsv", path_basis={"kind": "root-relative-v1", "root": "."}
        )
    with pytest.raises(ValueError, match=r"root|outside"):
        api.PortablePaths.create(root / "record.json", [outside], root=root)


@pytest.mark.parametrize("absolute", [False, True])
def test_legacy_file_references_remain_readable(
    tmp_path: Path, api: ModuleType, absolute: bool
) -> None:
    source = tmp_path / "source.tsv"
    source.write_text("legacy bytes\n")
    ref = {"path": str(source) if absolute else source.name, "sha256": file_sha256(source)}
    assert api.verify_reference(tmp_path / "record.json", ref, path_basis=None) == source


@pytest.mark.parametrize("digest", ["", "g" * 64, "0" * 63, None])
def test_source_hash_must_be_well_formed(tmp_path: Path, api: ModuleType, digest: object) -> None:
    source = tmp_path / "source.tsv"
    source.write_text("bytes\n")
    with pytest.raises(ValueError, match=r"hash|SHA|digest"):
        api.verify_reference(
            tmp_path / "record.json",
            {"path": source.name, "sha256": digest},
            path_basis={"kind": "root-relative-v1", "root": "."},
        )
