"""Refuse run and archive output roots that the system can empty (#1049).

Run F's verified run directories were kept only under the system temporary directory
and were lost on a restart (#1040). Outputs a record cites must be written somewhere
durable; a throwaway run or a test opts in explicitly.

Bandit's B108 check forbids a "/tmp" literal here, so /tmp is not a fixed root (#1051).
macOS's /tmp is still refused, because it resolves to /private/tmp. On Linux, /tmp is
refused only as the process temporary directory, which it is unless $TMPDIR points
elsewhere.
"""

from __future__ import annotations

import argparse
import os
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


def _spellings(path: Path) -> set[Path]:
    """A path as written, made absolute, and with every symlink resolved."""
    return {Path(os.path.abspath(path)), path.resolve()}


def temporary_root(path: Path, directories: Iterable[Path] | None = None) -> Path | None:
    """The temporary directory holding ``path`` in either spelling, or None.

    The resolved spelling refuses a durable-looking link into temporary storage. The
    executors resolve their output before recording it, so also checking the written
    spelling is a conservative over-refusal, not protection against a dangling record.
    """
    roots = sorted(
        {
            spelling
            for directory in (
                system_temporary_directories() if directories is None else directories
            )
            for spelling in _spellings(directory)
        }
    )
    for spelling in sorted(_spellings(path)):
        for root in roots:
            if spelling.is_relative_to(root):
                return root
    return None


def add_temporary_output_option(parser: argparse.ArgumentParser) -> None:
    """Declare the explicit opt-in for a throwaway run in temporary storage."""
    parser.add_argument(
        ALLOW_TEMPORARY_OUTPUT,
        action="store_true",
        help="Permit an output root under the system temporary directory, for tests and "
        "throwaway runs that no record will cite",
    )


def require_durable_output_path(output: Path, *, allowed: bool) -> Path:
    """Check the literal destination without writing; callers expand a home explicitly."""
    root = None if allowed else temporary_root(output)
    if root is not None:
        raise ValueError(
            f"output root {output} is under the temporary directory {root}, which the "
            "system can empty on restart (#1049). Write outputs that a record will cite "
            f"to a durable location, or pass {ALLOW_TEMPORARY_OUTPUT} for a throwaway run."
        )
    return output.resolve()


def require_durable_output(parser: argparse.ArgumentParser, output: Path, *, allowed: bool) -> None:
    """Exit with a usage error, before anything is written, if a restart could erase it."""
    try:
        require_durable_output_path(output, allowed=allowed)
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
