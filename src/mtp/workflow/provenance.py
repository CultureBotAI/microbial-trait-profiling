"""Check supplied provenance-root layouts before workflow output creation (#1817)."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from mtp.provenance_paths import PortablePaths, _validate_root, assert_portable_path_text


def preflight_provenance_root(
    record: Path,
    sources: Iterable[Path],
    *,
    root: Path,
    planned_files: Iterable[Path] = (),
) -> None:
    """Apply producer policy without creating any output.

    Known paths use the producer's normal checks. Declared future files receive
    file context for suffix handling before they exist; their directory and stem
    tokens remain checked. The final producer rechecks actual files and metadata.
    A comparison record may sit outside the supplied source root.
    """
    portable = PortablePaths.create(record, sources, root=root)
    for path in planned_files:
        resolved = path.resolve()
        _validate_root(portable.root, (resolved,))
        assert_portable_path_text(resolved.relative_to(portable.root).as_posix(), file_locator=True)
