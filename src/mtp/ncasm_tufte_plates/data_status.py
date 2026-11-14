"""Where each NCASM plate's numbers and data marks come from (#775, #776, #2264).

The manifest carries one record per plate, and the build acts on it: a plate whose
points come from a random generator is drawn only on request and carries a banner
saying so. A plate leaves that state when its builder reads run tables instead.
"""

from __future__ import annotations

from typing import Literal, TypedDict

__all__ = [
    "RANDOM_MARKS_BANNER",
    "DataStatus",
    "InputFile",
    "MarkSource",
    "NumberSource",
    "random_marks_status",
    "run_tables_status",
    "typed_literals_status",
]

NumberSource = Literal["typed_literals", "run_tables"]
MarkSource = Literal["typed_literals", "random_generator", "run_tables"]

RANDOM_MARKS_BANNER = (
    "ILLUSTRATION ONLY: the points on this plate are random draws, not data (#775)"
)


class InputFile(TypedDict):
    """A data file a plate was drawn from."""

    path: str
    sha256: str


class DataStatus(TypedDict):
    """The provenance of one plate's content.

    ``numbers`` covers every figure quoted in text; ``marks`` covers the points, bars
    and lines. ``typed_literals`` means the value is written in the plate's code with
    no recorded source. ``issues`` lists the open issues about the content of the plate
    set at this status; it is the same for every plate of a status and is not a
    plate-by-plate mapping (#2281)."""

    numbers: NumberSource
    marks: MarkSource
    inputs: list[InputFile]
    issues: list[int]


def typed_literals_status() -> DataStatus:
    """Plates I-XX as ported: numbers and marks are literals in the plate code (#776)."""
    return {
        "numbers": "typed_literals",
        "marks": "typed_literals",
        "inputs": [],
        "issues": [776, 777, 778, 779],
    }


def random_marks_status() -> DataStatus:
    """A data-story plot as ported: typed numbers and randomly drawn points (#775)."""
    return {
        "numbers": "typed_literals",
        "marks": "random_generator",
        "inputs": [],
        "issues": [775, 776, 777, 778, 779],
    }


def run_tables_status(inputs: tuple[InputFile, ...]) -> DataStatus:
    """A plate whose points and quoted numbers are all read from the listed run tables."""
    return {
        "numbers": "run_tables",
        "marks": "run_tables",
        "inputs": list(inputs),
        "issues": [],
    }
