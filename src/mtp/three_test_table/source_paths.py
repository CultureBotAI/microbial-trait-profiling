"""Read-only source locations for early publication guards (#1078, 2026-10-26).

These helpers discover locations, not scientific validity or byte identity. The
normal comparison/core/arm verifiers remain authoritative. File inputs protect
their parent directories; directory inputs protect only their own trees. Every
symlink expansion is retained so a middle file alias cannot disappear behind
Path.resolve(). No directories or outputs are created here.
"""

from __future__ import annotations

import json
import stat
from collections import deque
from pathlib import Path
from typing import Any

from mtp.datamodel.workflow_pydantic import StorageKind
from mtp.multiplicity.sources import VerifiedSources

_COMPARISON_LOCATIONS = (
    "manifest_path",
    "core_manifest_path",
    "noel_manifest_path",
    "jed_manifest_path",
    "ilyas_manifest_path",
)
_MAX_SYMLINK_EXPANSIONS = 40


def _absolute(path: Path) -> Path:
    # Do not collapse .. before expanding symlinks: that changes kernel meaning.
    return path if path.is_absolute() else Path.cwd() / path


def _spellings(path: Path) -> list[Path]:
    """Expand component links explicitly, retaining each complete path spelling."""
    original = _absolute(path)
    spellings = [original]
    visited = {original}
    pending = deque(original.parts[1:])
    resolved = Path(original.anchor)
    expansions = 0
    while pending:
        component = pending.popleft()
        if component == "..":
            resolved = resolved.parent
            continue
        candidate = resolved / component
        try:
            mode = candidate.lstat().st_mode
        except FileNotFoundError:
            resolved = candidate
            continue
        if not stat.S_ISLNK(mode):
            resolved = candidate
            continue
        expansions += 1
        if expansions > _MAX_SYMLINK_EXPANSIONS:
            raise ValueError("Symlink expansion limit reached while protecting source locations")
        # Include the physical spelling reached after earlier links and ..
        # components, as well as the untouched caller/target spellings.
        spellings.append(candidate.joinpath(*pending))
        target = candidate.readlink()
        expanded = target if target.is_absolute() else candidate.parent / target
        spelling = expanded.joinpath(*pending)
        if spelling in visited:
            raise ValueError("Symlink cycle while protecting source locations")
        visited.add(spelling)
        spellings.append(spelling)
        if target.is_absolute():
            resolved = Path(target.anchor)
            target_parts = target.parts[1:]
        else:
            resolved = candidate.parent
            target_parts = target.parts
        pending.extendleft(reversed(target_parts))
    spellings.append(resolved)
    return list(dict.fromkeys(spellings))


def _path_roots(path: Path, *, directory: bool | None) -> list[Path]:
    spellings = _spellings(path)
    is_directory = spellings[-1].is_dir() if directory is None else directory
    return spellings if is_directory else [spelling.parent for spelling in spellings]


def _location(value: object, description: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Invalid source path declaration: {description}")
    return Path(value)


def _record(path: Path) -> dict[str, Any]:
    record: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(record, dict):
        raise ValueError("Source location manifest must be a JSON object")
    return record


def _declared_storage(core: dict[str, Any]) -> dict[str, bool]:
    """Read only declared file/collection locations, not a workflow validity check."""
    if "definition" not in core:
        if core.get("artifacts"):
            raise ValueError("Core artifacts require declared product storage kinds")
        # Minimal location-only preflights may have no artifacts to classify.
        return {}
    definition = core["definition"]
    if not isinstance(definition, dict) or not isinstance(definition.get("products"), list):
        raise ValueError("Core definition must declare product storage kinds")
    kinds: dict[str, bool] = {}
    for product in definition["products"]:
        if not isinstance(product, dict):
            raise ValueError("Core product storage declaration must be an object")
        identity, kind = product.get("id"), product.get("storage_kind")
        if (
            not isinstance(identity, str)
            or not identity.strip()
            or identity in kinds
            or not isinstance(kind, str)
            or kind not in {"file", "collection"}
        ):
            raise ValueError("Core product storage declaration is missing, duplicate, or invalid")
        kinds[identity] = kind == StorageKind.collection.value
    return kinds


def manifest_roots(
    manifests: list[Path], *, multiplicity: Path | None, label_support: Path | None
) -> list[Path]:
    """Protect every file alias parent and each optional directory's own tree."""
    roots = [root for path in manifests for root in _path_roots(path, directory=False)]
    for optional in (multiplicity, label_support):
        if optional is not None:
            roots.extend(_path_roots(optional, directory=None))
    return list(dict.fromkeys(roots))


def declared_source_roots(
    comparison: Path, *, multiplicity: Path | None, label_support: Path | None
) -> list[Path]:
    """Read location declarations before any probe; never claim execution validity.

    A minimal synthetic core may omit artifacts. Real source verification still
    validates its complete model and inventory later. The actual core verifier
    resolves the supplied manifest before binding run-relative artifact paths.
    """
    roots = manifest_roots([comparison], multiplicity=multiplicity, label_support=label_support)
    record = _record(comparison)
    locations = {field: _location(record.get(field), field) for field in _COMPARISON_LOCATIONS}
    roots.extend(manifest_roots(list(locations.values()), multiplicity=None, label_support=None))
    core_path = locations["core_manifest_path"]
    try:
        core = _record(core_path)
    except FileNotFoundError:
        # Keep declared locations protected, but leave missing-source and run
        # status diagnostics to the authoritative execution verifier.
        return list(dict.fromkeys(roots))
    if "manifest_path" in core:
        roots.extend(
            _path_roots(_location(core["manifest_path"], "core manifest_path"), directory=False)
        )
    core_parent = _spellings(core_path)[-1].parent
    kinds = _declared_storage(core)
    artifacts = core.get("artifacts")
    if artifacts is None:
        artifacts = []
    if not isinstance(artifacts, list):
        raise ValueError("Core artifact location declarations must be a list")
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise ValueError("Core artifact location declaration must be an object")
        identity = artifact.get("product_id")
        if not isinstance(identity, str) or identity not in kinds:
            raise ValueError("Core artifact lacks a declared product storage kind")
        directory = kinds[identity]
        path = _location(artifact.get("path"), "core artifact path")
        roots.extend(
            _path_roots(path if path.is_absolute() else core_parent / path, directory=directory)
        )
        requested = artifact.get("requested_path")
        if requested is not None and requested != "":
            roots.extend(
                _path_roots(
                    _location(requested, "core artifact requested_path"), directory=directory
                )
            )
    return list(dict.fromkeys(roots))


def source_roots(
    sources: VerifiedSources, *, multiplicity: Path | None, label_support: Path | None
) -> list[Path]:
    """Protect verified manifests, all concrete core products and requested aliases."""
    core_path = Path(sources.core.manifest_path)
    manifests = [
        core_path,
        Path(sources.comparison.manifest_path),
        *(Path(arm.manifest_path) for arm in sources.arms.values()),
    ]
    for field in _COMPARISON_LOCATIONS:
        value = getattr(sources.comparison, field, None)
        if value is not None:
            manifests.append(_location(value, field))
    roots = manifest_roots(manifests, multiplicity=multiplicity, label_support=label_support)
    kinds = {product.id: product.storage_kind for product in sources.core.definition.products}
    for artifact in sources.core.artifacts or []:
        directory = kinds[artifact.product_id] == StorageKind.collection
        roots.extend(_path_roots(sources.product(artifact.product_id), directory=directory))
        path = _location(artifact.path, "core artifact path")
        roots.extend(
            _path_roots(
                path if path.is_absolute() else core_path.parent / path, directory=directory
            )
        )
        requested = artifact.requested_path
        if requested is not None and requested != "":
            roots.extend(
                _path_roots(
                    _location(requested, "core artifact requested_path"), directory=directory
                )
            )
    return list(dict.fromkeys(roots))


__all__ = ["declared_source_roots", "manifest_roots", "source_roots"]
