"""Refuse run and archive output roots that the system can empty (#1049).

Run F's verified run directories were kept only under the system temporary directory
and were lost on a restart (#1040). Outputs a record cites must be written somewhere
durable; a throwaway run or a test opts in explicitly.

Symlink loops are invalid destinations even with that opt-in (#1175). Strict probes
before and after nonstrict normalization detect loops even when an initial missing
component masks a later loop (#1190). Fresh destinations retain nonstrict semantics.
Filesystem identity checks on existing ancestors also detect case and firmlink aliases
without casefold guesses (#1178).
These are point-in-time checks, not protection against concurrent filesystem changes.

Bandit's B108 check forbids a "/tmp" literal here, so /tmp is not a fixed root (#1051).
macOS's /tmp is still refused, because it resolves to /private/tmp. On Linux, /tmp is
refused only as the process temporary directory, which it is unless $TMPDIR points
elsewhere.
"""

from __future__ import annotations

import argparse
import errno
import os
import stat
import tempfile
from collections.abc import Iterable
from pathlib import Path

ALLOW_TEMPORARY_OUTPUT = "--allow-temporary-output"
# macOS links /var to /private/var; both spellings are listed so that the check does not
# depend on which links this system has.
FIXED_TEMPORARY_ROOTS = tuple(
    Path(item) for item in ("/private/tmp", "/var/folders", "/private/var/folders")
)


def system_temporary_directories() -> list[Path]:
    """The fixed roots, this process's temporary directory and ``$TMPDIR``."""
    directories = [*FIXED_TEMPORARY_ROOTS, Path(tempfile.gettempdir())]
    if environment := os.environ.get("TMPDIR"):
        directories.append(Path(environment))
    return directories


def _resolved_path(path: Path, *, role: str = "output path") -> Path:
    """Allow missing suffixes, but refuse other resolution failures on every platform."""
    try:
        return path.resolve(strict=True)
    except RuntimeError as error:
        # Python 3.11/3.12 use RuntimeError for a symlink loop.
        raise ValueError(f"Cannot resolve {role} {path}: symbolic link loop") from error
    except OSError as error:
        # Newer pathlib reports ELOOP only when strict=True. Only ENOENT can
        # represent a fresh output; permission and unsupported-path errors refuse.
        _require_missing_path(path, error, role=role)
    try:
        resolved = path.resolve(strict=False)
    except RuntimeError as error:
        raise ValueError(f"Cannot resolve {role} {path}: symbolic link loop") from error
    except OSError as error:
        raise ValueError(f"Cannot resolve {role} {path}: {error}") from error
    # A missing/../ prefix can stop the first strict probe before a real loop.
    # Nonstrict normalization removes that prefix but newer pathlib suppresses
    # the loop; probe the resulting spelling before accepting a fresh output.
    try:
        resolved.resolve(strict=True)
    except RuntimeError as error:
        raise ValueError(f"Cannot resolve {role} {path}: symbolic link loop") from error
    except OSError as error:
        _require_missing_path(path, error, role=role)
    return resolved


def _spellings(path: Path, *, role: str = "output path") -> set[Path]:
    """A path as written, made absolute, and with every symlink resolved."""
    return {Path(os.path.abspath(path)), _resolved_path(path, role=role)}


def _require_missing_path(path: Path, error: OSError, *, role: str = "output path") -> None:
    """Only ENOENT is evidence of a new path, never a failed inspection."""
    if error.errno == errno.ELOOP:
        raise ValueError(f"Cannot resolve {role} {path}: symbolic link loop") from error
    if error.errno != errno.ENOENT:
        raise ValueError(f"Cannot inspect {role} {path}: {error}") from error


def _filesystem_identity(
    path: Path, *, directory: bool = False, role: str = "output path"
) -> tuple[int, int] | None:
    """Follow aliases; refuse unusable identities rather than guessing containment."""
    try:
        information = path.stat()
    except OSError as error:
        _require_missing_path(path, error, role=role)
        return None
    device = getattr(information, "st_dev", None)
    inode = getattr(information, "st_ino", None)
    if (
        not isinstance(device, int)
        or isinstance(device, bool)
        or device < 0
        or not isinstance(inode, int)
        or isinstance(inode, bool)
        or inode <= 0
    ):
        raise ValueError(f"Cannot establish filesystem identity for {role} {path}")
    if directory:
        mode = getattr(information, "st_mode", None)
        if not isinstance(mode, int) or isinstance(mode, bool) or not stat.S_ISDIR(mode):
            raise ValueError(f"Cannot inspect {role} {path}: no valid directory mode")
    return device, inode


def containing_root(path: Path, directories: Iterable[Path]) -> Path | None:
    """Return a containing root by spelling or existing directory identity (#1178).

    Inspect existing candidate ancestors through the filesystem root in both
    written and resolved spellings. This detects case aliases and macOS
    firmlinks even when resolve() preserves their spelling. Identity includes the
    device: inode numbers alone are not unique across filesystems. Never promote
    a missing protected root to its existing parent; retain only its lexical bound.
    Failed inspections and unusable identities refuse; nothing is created. This
    checks the filesystem at inspection time and is not a race-proof write guard.
    """
    candidates = sorted(_spellings(path))
    roots = sorted(
        {
            spelling
            for directory in directories
            for spelling in _spellings(directory, role="protected root")
        }
    )
    for candidate in candidates:
        for root in roots:
            if candidate.is_relative_to(root):
                return root
    # An existing output may be a file (for example a pinned Madin table).
    # Missing suffixes are skipped, but all existing ancestors are inspected.
    ancestors = {
        ancestor for candidate in candidates for ancestor in (candidate, *candidate.parents)
    }
    identities = {
        identity
        for ancestor in sorted(ancestors)
        if (identity := _filesystem_identity(ancestor)) is not None
    }
    for root in roots:
        identity = _filesystem_identity(root, directory=True, role="protected root")
        if identity is not None and identity in identities:
            return root
    return None


def temporary_root(path: Path, directories: Iterable[Path] | None = None) -> Path | None:
    """The temporary root containing ``path`` by spelling or filesystem identity.

    Checking the written spelling also conservatively refuses a path beneath a
    temporary root whose symlink resolves out to durable storage.
    """
    return containing_root(
        path, system_temporary_directories() if directories is None else directories
    )


def add_temporary_output_option(parser: argparse.ArgumentParser) -> None:
    """Declare the explicit opt-in for a throwaway run in temporary storage."""
    parser.add_argument(
        ALLOW_TEMPORARY_OUTPUT,
        action="store_true",
        help="Permit an output root under the system temporary directory, for tests and "
        "throwaway runs that no record will cite",
    )


def require_durable_output_path(
    output: Path, *, allowed: bool, include_temporary_hint: bool = True
) -> Path:
    """Check the literal destination without writing; callers expand a home explicitly.

    The temporary-storage opt-in does not waive resolution or identity inspection.
    Callers without an opt-in suppress its hint without changing that policy.
    """
    root = containing_root(output, ()) if allowed else temporary_root(output)
    if root is not None:
        hint = (
            f", or pass {ALLOW_TEMPORARY_OUTPUT} for a throwaway run."
            if include_temporary_hint
            else "."
        )
        raise ValueError(
            f"output root {output} is under the temporary directory {root}, which the "
            "system can empty on restart (#1049). Write outputs that a record will cite "
            f"to a durable location{hint}"
        )
    return _resolved_path(output)


def require_durable_output(
    parser: argparse.ArgumentParser,
    output: Path,
    *,
    allowed: bool,
    include_temporary_hint: bool = True,
) -> None:
    """Exit with a usage error, before anything is written, if a restart could erase it."""
    try:
        require_durable_output_path(
            output, allowed=allowed, include_temporary_hint=include_temporary_hint
        )
    except ValueError as error:
        parser.error(str(error))


def main(argv: list[str] | None = None) -> None:
    """Check one output root for a recipe whose module the arms also run (#1087).

    Guarding such a module itself would refuse arm runs made with the opt-in, and adding
    the flag to its argv would change the arms' recorded plans, so the recipe checks first.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("output", type=Path)
    add_temporary_output_option(parser)
    args = parser.parse_args(argv)
    require_durable_output(parser, args.output, allowed=args.allow_temporary_output)


if __name__ == "__main__":
    main()
