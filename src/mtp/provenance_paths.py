"""Portable source locations for standalone provenance records (#1664).

The declaration locates a common source root relative to the record's directory.
Every reference is a canonical POSIX path beneath that root. Relocation preserves
this tree, not an original machine prefix. Existing records without a declaration
retain their absolute or record-directory-relative path semantics. Hashes bind
file bytes independently of the location spelling.

An explicitly supplied source root can also support a record outside that tree.
That record carries only ``supplied-root-v1``; readers must supply its root again.
No machine path is serialized to reconnect otherwise unrelated trees.
"""

from __future__ import annotations

import os
import pwd
import re
import tempfile
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path, PurePosixPath, PureWindowsPath

from mtp.paths import file_sha256

_KIND = "root-relative-v1"
_SUPPLIED_KIND = "supplied-root-v1"
_LOCAL_PATH = re.compile(
    r"/(?:Users|home)/[A-Za-z0-9._-]+|[A-Za-z]:[\\/]+Users[\\/]+"
    r"|(?:/private)?/var/(?:folders)(?=/)|/private/(?:tmp)(?=/)",
    re.IGNORECASE,
)
_SLASHLESS_LOCAL_PATH = re.compile(
    r"(?<![\w])(?:(?:Users|home)/[A-Za-z0-9._-]+"
    r"|(?:private/)?var/folders(?=/)|private/tmp(?=/))",
    re.IGNORECASE,
)
_ROOT_GUIDANCE = (
    "Choose --provenance-root below the home or temporary-directory prefix; "
    "all bound sources must be beneath it. Stage verified copies together if needed."
)


def assert_portable_text(text: str) -> None:
    """Refuse leaked machine prefixes; never silently redact source metadata."""
    if _LOCAL_PATH.search(text) or _SLASHLESS_LOCAL_PATH.search(text):
        raise ValueError("Provenance contains a local path or local prefix. " + _ROOT_GUIDANCE)
    folded = text.casefold()
    if any(name in folded for name in _account_names()):
        raise ValueError("Provenance contains a local account name. " + _ROOT_GUIDANCE)


def _account_names() -> set[str]:
    """Inspect every account authority, not just the first populated environment key."""
    names = {os.environ.get(key, "") for key in ("LOGNAME", "USER", "LNAME", "USERNAME")}
    names.add(Path.home().name)
    try:
        account = pwd.getpwuid(os.getuid())
    except KeyError:
        # A container may have a numeric UID without a password-database entry.
        pass
    else:
        names.update((account.pw_name, Path(account.pw_dir).name))
    return {name.casefold() for name in names if name}


def _private_roots() -> set[Path]:
    """Runtime home/temp boundaries and platform prefixes must stay above a bundle."""
    roots = {Path.home(), Path(tempfile.gettempdir())}
    roots.update(Path(value) for key in ("TMPDIR", "TEMP", "TMP") if (value := os.getenv(key)))
    try:
        account = pwd.getpwuid(os.getuid())
    except KeyError:
        pass
    else:
        if account.pw_dir:
            roots.add(Path(account.pw_dir))
    anchor = Path(Path.cwd().anchor)
    roots.update(
        anchor.joinpath(*parts)
        for parts in (
            ("tmp",),
            ("private", "tmp"),
            ("var", "tmp"),
            ("private", "var", "tmp"),
            ("var", "folders"),
            ("private", "var", "folders"),
        )
    )
    return {root.resolve() for root in roots}


def _require_root_directory(root: Path) -> None:
    if not root.exists():
        raise ValueError("--provenance-root does not exist; choose an existing source directory")
    if not root.is_dir():
        raise ValueError("--provenance-root is not a directory")


def _validate_root(root: Path, sources: Iterable[Path]) -> None:
    """Never expose a home/temp prefix by removing only its leading slash."""
    if root == root.parent:
        raise ValueError("The filesystem root is not a provenance root. " + _ROOT_GUIDANCE)
    # A case alias can name the same directory on macOS. Be conservative on
    # case-sensitive filesystems too rather than make publication depend on it.
    folded_root = PurePosixPath(root.as_posix().casefold())
    if any(
        PurePosixPath(boundary.as_posix().casefold()).is_relative_to(folded_root)
        for boundary in _private_roots()
    ):
        raise ValueError(
            "Provenance root must be below a home or temporary root. " + _ROOT_GUIDANCE
        )
    if any(match.end() == len(root.as_posix()) for match in _LOCAL_PATH.finditer(root.as_posix())):
        raise ValueError("Provenance root equals a private home/temp prefix. " + _ROOT_GUIDANCE)
    for source in sources:
        if not source.is_relative_to(root):
            raise ValueError("Source is outside the declared provenance root. " + _ROOT_GUIDANCE)
        if any(
            match.end() >= len(root.as_posix()) for match in _LOCAL_PATH.finditer(source.as_posix())
        ):
            raise ValueError(
                "Provenance root leaves a local prefix in a reference. " + _ROOT_GUIDANCE
            )
        assert_portable_text("/" + source.relative_to(root).as_posix())


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
        if root is not None:
            _require_root_directory(base)
        _validate_root(base, paths)
        return cls(record, base)

    @property
    def declaration(self) -> dict[str, str]:
        if not self.record.parent.is_relative_to(self.root):
            return {"kind": _SUPPLIED_KIND}
        return {
            "kind": _KIND,
            "root": os.path.relpath(self.root, self.record.parent).replace(os.sep, "/"),
        }

    def encode(self, path: Path) -> str:
        resolved = path.resolve()
        _validate_root(self.root, [resolved])
        value = resolved.relative_to(self.root).as_posix()
        _canonical_reference(value)
        return value


def resolve_provenance_root(
    record_path: Path,
    *,
    path_basis: Mapping[str, object] | None,
    provenance_root: Path | None = None,
) -> Path | None:
    """Validate a root declaration even when a record has no present sources."""
    parent = record_path.resolve().parent
    if path_basis is None:
        return None
    if not isinstance(path_basis, Mapping):
        raise ValueError("Invalid portable path basis declaration")
    if path_basis.get("kind") == _SUPPLIED_KIND:
        if set(path_basis) != {"kind"}:
            raise ValueError("Invalid supplied-root declaration")
        if provenance_root is None:
            raise ValueError("This record requires --provenance-root to resolve its sources")
        base = provenance_root.resolve()
        _require_root_directory(base)
        _validate_root(base, [])
        return base
    root = path_basis.get("root")
    if (
        set(path_basis) != {"kind", "root"}
        or path_basis.get("kind") != _KIND
        or not isinstance(root, str)
        or (root != "." and re.fullmatch(r"(?:\.\./)*\.\.", root) is None)
    ):
        raise ValueError("Invalid portable path basis/root declaration")
    base = (parent / root).resolve()
    _validate_root(base, [])
    if provenance_root is not None and provenance_root.resolve() != base:
        raise ValueError("--provenance-root conflicts with the record-relative root")
    return base


def resolve_reference(
    record_path: Path,
    path: str,
    *,
    path_basis: Mapping[str, object] | None,
    provenance_root: Path | None = None,
) -> Path:
    """Resolve legacy locations or strictly checked declared-root references."""
    if not isinstance(path, str) or not path:
        raise ValueError("Source path must be a nonempty string")
    base = resolve_provenance_root(
        record_path, path_basis=path_basis, provenance_root=provenance_root
    )
    if base is None:
        return (record_path.resolve().parent / path).resolve()
    result = (base / _canonical_reference(path)).resolve()
    if not result.is_relative_to(base):
        raise ValueError("Source path escapes the declared provenance root")
    _validate_root(base, [result])
    return result


def verify_reference(
    record_path: Path,
    reference: Mapping[str, object],
    *,
    path_basis: Mapping[str, object] | None,
    provenance_root: Path | None = None,
) -> Path:
    """Resolve a recorded source and reject missing or different file bytes."""
    path, digest = reference.get("path"), reference.get("sha256")
    if not isinstance(path, str):
        raise ValueError("Source path must be a string")
    if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise ValueError("Source SHA-256 digest must be 64 lowercase hexadecimal characters")
    resolved = resolve_reference(
        record_path, path, path_basis=path_basis, provenance_root=provenance_root
    )
    try:
        actual = file_sha256(resolved)
    except OSError as exc:
        raise ValueError("Recorded source is missing or unavailable to read") from exc
    if actual != digest:
        raise ValueError("Source hash does not match the recorded SHA-256 digest")
    return resolved
