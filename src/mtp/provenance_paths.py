"""Portable source locations for standalone provenance records (#1664).

The declaration locates a common source root relative to the record's directory.
Every reference is a canonical POSIX path beneath that root. Relocation preserves
this tree, not an original machine prefix. Existing records without a declaration
retain their absolute or record-directory-relative path semantics. Hashes bind
file bytes independently of the location spelling.
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path, PurePosixPath, PureWindowsPath

from mtp.paths import file_sha256

_KIND = "root-relative-v1"
_LOCAL_PATH = re.compile(
    r"/(?:Users|home)/[A-Za-z0-9._-]+|[A-Za-z]:[\\/]+Users[\\/]+"
    r"|(?:/private)?/var/(?:folders)(?=/)|/private/(?:tmp)(?=/)"
)


def assert_portable_text(text: str) -> None:
    """Refuse leaked machine prefixes; never silently redact source metadata."""
    if _LOCAL_PATH.search(text):
        raise ValueError(
            "Provenance contains an absolute local path; declare portable source paths"
        )


def _canonical_reference(value: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if (
        not value
        or value == "."
        or path.is_absolute()
        or PureWindowsPath(value).drive
        or "\\" in value
        or ".." in path.parts
        or path.as_posix() != value
    ):
        raise ValueError("Source path must be a canonical root-relative file path")
    return path


@dataclass(frozen=True)
class PortablePaths:
    """Locations in one declared source tree, with no absolute root in the record."""

    record: Path
    root: Path

    @classmethod
    def create(
        cls,
        record: Path,
        sources: Iterable[Path],
        *,
        root: Path | None = None,
    ) -> PortablePaths:
        record = record.resolve()
        paths = tuple(source.resolve() for source in sources)
        base = (
            root.resolve()
            if root is not None
            else Path(os.path.commonpath([record.parent, *(path.parent for path in paths)]))
        )
        if not record.parent.is_relative_to(base) or any(
            not path.is_relative_to(base) for path in paths
        ):
            raise ValueError("Source or output is outside the declared provenance root")
        return cls(record, base)

    @property
    def declaration(self) -> dict[str, str]:
        return {
            "kind": _KIND,
            "root": os.path.relpath(self.root, self.record.parent).replace(os.sep, "/"),
        }

    def encode(self, path: Path) -> str:
        resolved = path.resolve()
        if not resolved.is_relative_to(self.root):
            raise ValueError("Source path is outside the declared provenance root")
        value = resolved.relative_to(self.root).as_posix()
        _canonical_reference(value)
        return value


def resolve_reference(
    record_path: Path,
    path: str,
    *,
    path_basis: Mapping[str, object] | None,
) -> Path:
    """Resolve legacy locations or strictly checked declared-root references."""
    if not isinstance(path, str) or not path:
        raise ValueError("Source path must be a nonempty string")
    parent = record_path.resolve().parent
    if path_basis is None:
        return (parent / path).resolve()
    root = path_basis.get("root")
    if (
        set(path_basis) != {"kind", "root"}
        or path_basis.get("kind") != _KIND
        or not isinstance(root, str)
        or (root != "." and re.fullmatch(r"(?:\.\./)*\.\.", root) is None)
    ):
        raise ValueError("Invalid portable path basis/root declaration")
    base = (parent / root).resolve()
    result = (base / _canonical_reference(path)).resolve()
    if not result.is_relative_to(base):
        raise ValueError("Source path escapes the declared provenance root")
    return result


def verify_reference(
    record_path: Path,
    reference: Mapping[str, object],
    *,
    path_basis: Mapping[str, object] | None,
) -> Path:
    """Resolve a recorded source and reject missing or different file bytes."""
    path, digest = reference.get("path"), reference.get("sha256")
    if not isinstance(path, str):
        raise ValueError("Source path must be a string")
    if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise ValueError("Source SHA-256 digest must be 64 lowercase hexadecimal characters")
    resolved = resolve_reference(record_path, path, path_basis=path_basis)
    if file_sha256(resolved) != digest:
        raise ValueError("Source hash does not match the recorded SHA-256 digest")
    return resolved
