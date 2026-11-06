"""Show that independent source pins reject changed displayed values (#1551, #1555)."""

from __future__ import annotations

import copy
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from functools import partial
from pathlib import Path
from typing import Any

import pytest
from matplotlib.figure import Figure

import test_cellulose_tufte_overlap_sources as overlap
import test_cellulose_tufte_statistics_sources as statistics
from test_cellulose_tufte_historical_labels import label_plates as label_plates
from test_cellulose_tufte_overlap_sources import reference as reference


@contextmanager
def _changed_text(figure: Figure, old: str, new: str) -> Iterator[None]:
    artists = [*figure.texts, *(text for axis in figure.axes for text in axis.texts)]
    matches = [artist for artist in artists if old in artist.get_text()]
    assert len(matches) == 1, (old, [artist.get_text() for artist in matches])
    artist = matches[0]
    original = artist.get_text()
    artist.set_text(original.replace(old, new))
    try:
        yield
    finally:
        artist.set_text(original)


@pytest.mark.parametrize(
    ("plate", "old", "new", "check"),
    [
        (
            "V",
            "matrix (2,335 experiments",
            "matrix (2,336 experiments",
            statistics.test_plate_five_shared_corpus_and_cascade_match_the_pinned_methods,
        ),
        (
            "VI",
            "single 1,912-row table",
            "single 1,913-row table",
            statistics.test_plate_six_total_rows_match_the_historical_methods_snapshot,
        ),
        (
            "V",
            "(6 traits) & BY (15 ratios, 317 traits)",
            "(6 traits) & BY (15 ratios, 344 traits)",
            statistics.test_plate_five_retained_trait_family_matches_the_matrix_source,
        ),
        (
            "V",
            "58.2% (n=917",
            "59.2% (n=917",
            partial(
                statistics.test_plate_five_fraction_labels_match_the_tracked_validation_rows,
                trait="FA",
            ),
        ),
        (
            "V",
            "59.5% up (n=467, BY q = 0.0011",
            "59.5% up (n=467, BY q = 0.0012",
            partial(
                statistics.test_plate_five_ratio_labels_match_the_tracked_validation_rows,
                ratio="FA/CEL",
            ),
        ),
        (
            "VI",
            "drawn (316 cells, 16.9%)",
            "drawn (317 cells, 16.9%)",
            partial(
                statistics.test_plate_six_cell_counts_match_the_historical_inventory,
                display="drawn",
                state="untested",
                prefix="drawn",
            ),
        ),
        (
            "V",
            "66,240 enrichments",
            "66,241 enrichments",
            statistics.test_plate_five_ora_numbers_match_the_historical_export_audit,
        ),
    ],
    ids=[
        "paired-corpus",
        "total-rows",
        "family-size",
        "fraction-share",
        "ratio-q",
        "cell-count",
        "ora-hits",
    ],
)
def test_display_number_mutations_are_rejected_by_source_pins(
    label_plates: dict[str, Figure],
    plate: str,
    old: str,
    new: str,
    check: Callable[[dict[str, Figure]], None],
) -> None:
    check(label_plates)
    with _changed_text(label_plates[plate], old, new), pytest.raises(AssertionError):
        check(label_plates)


def test_overlap_number_mutation_is_rejected_by_the_historical_source(
    label_plates: dict[str, Figure],
    reference: dict[str, Any],
) -> None:
    check = partial(
        overlap.test_plate_one_pair_labels_read_historical_intersections_and_unions,
        reference=reference,
        pair="PD vs. CRC",
    )
    check(label_plates)
    with (
        _changed_text(label_plates["I"], "PD vs. CRC (27/47 shared)", "PD vs. CRC (27/48 shared)"),
        pytest.raises(AssertionError),
    ):
        check(label_plates)


def test_historical_reference_change_cannot_be_ignored(
    label_plates: dict[str, Figure],
    reference: dict[str, Any],
) -> None:
    changed = copy.deepcopy(reference)
    changed["sets"]["intersections"]["PD vs. CRC"] += 1
    with pytest.raises(AssertionError):
        overlap.test_plate_one_pair_labels_read_historical_intersections_and_unions(
            label_plates, changed, "PD vs. CRC"
        )


@pytest.mark.parametrize(
    ("path", "old", "new", "check"),
    [
        (
            statistics.VALIDATION,
            "| FA | 534/383/42 | 917 |",
            "| FA | 534/383/42 | 918 |",
            partial(
                statistics.test_plate_five_fraction_labels_match_the_tracked_validation_rows,
                trait="FA",
            ),
        ),
        (
            statistics.ORA_AUDIT,
            "66,240",
            "66,241",
            statistics.test_plate_five_ora_numbers_match_the_historical_export_audit,
        ),
        (
            statistics.HISTORICAL / "integrated_evidence_a3c32cf.md",
            "| drawn | untested | 316 |",
            "| drawn | untested | 317 |",
            partial(
                statistics.test_plate_six_cell_counts_match_the_historical_inventory,
                display="drawn",
                state="untested",
                prefix="drawn",
            ),
        ),
    ],
    ids=["validation-row", "ora-audit", "historical-grid"],
)
def test_changed_source_rows_cannot_be_ignored(
    label_plates: dict[str, Figure],
    monkeypatch: pytest.MonkeyPatch,
    path: Path,
    old: str,
    new: str,
    check: Callable[[dict[str, Figure]], None],
) -> None:
    check(label_plates)
    original = path.read_text(encoding="utf-8")
    assert old in original
    changed = original.replace(old, new)
    reader = Path.read_text

    def substitute(source: Path, *args: Any, **kwargs: Any) -> str:
        return changed if source == path else reader(source, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", substitute)
    with pytest.raises(AssertionError):
        check(label_plates)


def test_threshold_column_swap_cannot_preserve_the_source_pin(
    label_plates: dict[str, Figure],
) -> None:
    panel = label_plates["VI"].axes[2]
    key = next(text for text in panel.texts if text.get_text() == "table1_Burst_direction")
    values = [
        text
        for text in panel.texts
        if text.get_position()[1] == key.get_position()[1] and text.get_text() in {"6 / 6", "3 / 6"}
    ]
    assert len(values) == 2
    originals = [text.get_text() for text in values]
    for text, value in zip(values, reversed(originals), strict=True):
        text.set_text(value)
    try:
        with pytest.raises(AssertionError):
            statistics.test_plate_six_method_counts_match_the_historical_methods_snapshot(
                label_plates,
                "directional consistency, fraction sign test",
                "table1_Burst_direction",
            )
    finally:
        for text, value in zip(values, originals, strict=True):
            text.set_text(value)
