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


def test_cross_top_level_inference_cannot_publish_machine_prefixes(api: ModuleType) -> None:
    anchor = Path(Path.cwd().anchor)
    record = anchor / "Users" / "example-user" / "records" / "record.json"
    source = anchor / "private" / "tmp" / "input.tsv"
    with pytest.raises(ValueError, match="--provenance-root"):
        api.PortablePaths.create(record, [source])


@pytest.mark.parametrize("root_parts", [("Users",), ("private",)])
def test_explicit_root_cannot_leave_local_prefix_components(
    api: ModuleType, root_parts: tuple[str, ...]
) -> None:
    anchor = Path(Path.cwd().anchor)
    root = anchor.joinpath(*root_parts)
    suffix = ("example-user", "source.tsv") if root_parts == ("Users",) else ("tmp", "source.tsv")
    with pytest.raises(ValueError, match=r"root|local prefix"):
        api.PortablePaths.create(root / "record.json", [root.joinpath(*suffix)], root=root)


def test_filesystem_root_is_never_an_explicit_publication_root(api: ModuleType) -> None:
    root = Path(Path.cwd().anchor)
    with pytest.raises(ValueError, match="--provenance-root"):
        api.PortablePaths.create(root / "record.json", [root / "data.tsv"], root=root)


@pytest.mark.parametrize(
    "components",
    [
        ("Users", "example-user", "data.tsv"),
        ("home", "example-user", "data.tsv"),
        ("private", "tmp", "run", "data.tsv"),
    ],
)
def test_text_backstop_detects_slashless_machine_prefixes(
    api: ModuleType, components: tuple[str, ...]
) -> None:
    text = json.dumps({"path": "/".join(components)})
    with pytest.raises(ValueError, match=r"local path|local prefix"):
        api.assert_portable_text(text)


def test_absolute_reference_inside_root_is_still_noncanonical(
    tmp_path: Path, api: ModuleType
) -> None:
    with pytest.raises(ValueError, match="canonical root-relative"):
        api.resolve_reference(
            tmp_path / "record.json",
            str(tmp_path / "source.tsv"),
            path_basis={"kind": "root-relative-v1", "root": "."},
        )


def test_explicit_source_root_supports_an_external_report_without_disclosure(
    tmp_path: Path, api: ModuleType
) -> None:
    root = tmp_path / "source-tree"
    root.mkdir()
    source = root / "source.tsv"
    source.write_text("retained bytes\n")
    record = tmp_path / "elsewhere" / "record.json"
    paths = api.PortablePaths.create(record, [source], root=root)
    assert paths.declaration == {"kind": "supplied-root-v1"}
    ref = {"path": paths.encode(source), "sha256": file_sha256(source)}
    assert ref["path"] == "source.tsv"
    assert str(tmp_path) not in json.dumps({"path_basis": paths.declaration, "source": ref})
    with pytest.raises(ValueError, match="--provenance-root"):
        api.verify_reference(record, ref, path_basis=paths.declaration)
    assert (
        api.verify_reference(record, ref, path_basis=paths.declaration, provenance_root=root)
        == source
    )
    relocated = tmp_path / "relocated-inputs"
    root.rename(relocated)
    assert (
        api.verify_reference(record, ref, path_basis=paths.declaration, provenance_root=relocated)
        == relocated / source.name
    )


def test_explicit_resolution_cannot_override_a_record_relative_root(
    tmp_path: Path, api: ModuleType
) -> None:
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    with pytest.raises(ValueError, match=r"match|conflict"):
        api.resolve_reference(
            tmp_path / "record.json",
            "source.tsv",
            path_basis={"kind": "root-relative-v1", "root": "."},
            provenance_root=elsewhere,
        )


def test_missing_reference_is_a_clear_verification_error(tmp_path: Path, api: ModuleType) -> None:
    with pytest.raises(ValueError, match=r"missing|read|unavailable"):
        api.verify_reference(
            tmp_path / "record.json",
            {"path": "missing.tsv", "sha256": file_sha256(__file__)},
            path_basis={"kind": "root-relative-v1", "root": "."},
        )
