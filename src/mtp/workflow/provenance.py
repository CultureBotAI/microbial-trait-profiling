"""Check supplied provenance-root layouts before workflow output creation (#1817)."""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import NoReturn

from mtp.provenance_paths import PortablePaths, _validate_root, assert_portable_path_text


class _SourceArguments(argparse.ArgumentParser):
    def error(self, message: str) -> NoReturn:
        raise ValueError(message)


def preflight_integrated_sources(
    record: Path,
    argv: Sequence[str],
    *,
    count_profiles: Path,
    ora_profiles: Path,
    root: Path,
    planned_files: Iterable[Path] = (),
) -> None:
    """Select the producer's located sources, excluding hash-only core inputs.

    The caller has verified the core. Reuse the producer's parser, native defaults
    and source selector without repeating scientific validation. A missing optional
    assessment has no locator unless an earlier step will actually produce it.
    Metadata inside assessment files remains the final producer's responsibility.
    """
    from mtp.integrated_evidence import (
        _build_parser,
        _native_source_paths,
        _provenance_sources,
    )

    args = _build_parser(_SourceArguments).parse_args(argv)
    if args.profile_input is not None and args.profile_input.resolve() != count_profiles.resolve():
        raise ValueError("--core-run conflicts with explicit --profile-input")
    args.profile_input = count_profiles
    native = _native_source_paths(args)
    assessments = tuple(
        getattr(args, name)
        for name in (
            "calibration_status",
            "sensitivity_status",
            "orientation_provenance",
            "coverage_status",
            "influence_status",
        )
    )
    selected = {path.resolve() for path in (*native.values(), *assessments) if path is not None}
    declared_outputs = tuple(path.resolve() for path in planned_files)
    future = tuple(
        path for path in selected if any(path.is_relative_to(output) for output in declared_outputs)
    )
    sources = _provenance_sources(native, assessment_paths=assessments, ora_profiles=ora_profiles)
    preflight_provenance_root(
        record,
        sources,
        root=root,
        planned_files=future,
    )


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
