"""Explicit retained-file resolution without altering historical records (#726)."""

from collections.abc import Callable
from pathlib import Path

LocationResolver = Callable[[Path], Path]


def locate(path: Path, resolver: LocationResolver | None = None) -> Path:
    """Use the supplied resolver exclusively; archives must never fall back to live files."""
    return resolver(path) if resolver is not None else path
