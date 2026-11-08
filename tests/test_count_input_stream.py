"""Caller-owned audit streams preserve the strict count-table contract (#1357)."""

from __future__ import annotations

import csv
import io
from pathlib import Path
from types import GeneratorType
from unittest.mock import Mock

import pytest

from mtp.count_input import (
    LEGACY_SPARSE_ZERO_MODE,
    PROFILE_KEY_COLUMNS,
    CountTableReader,
    profile_count_reader,
    taxa_trait_count_reader,
)

HEADER = "study_id\texperiment_id\tinc_dec\ttrait_a\n"


class _NonseekableText(io.StringIO):
    def seek(self, offset: int, whence: int = 0) -> int:
        raise AssertionError(f"the provided stream must not be rewound: {offset=}, {whence=}")


def test_header_and_rows_share_reader_without_closing_or_reopening(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    csv_reader = Mock(wraps=csv.reader)
    monkeypatch.setattr(csv, "reader", csv_reader)
    stream = _NonseekableText(
        'study_id\texperiment_id\tinc_dec\t"trait\nname"\ns\te\tinc\t2\ns\te\tdec\t0\n'
    )
    reader = profile_count_reader(tmp_path / "absent.tsv.gz", text_stream=stream)
    assert reader.count_columns == ("trait\nname",)
    assert not stream.closed
    rows = list(reader)
    assert [(row.key, row.counts, row.line_number) for row in rows] == [
        (("s", "e", "inc"), (2,), 3),
        (("s", "e", "dec"), (0,), 4),
    ]
    assert reader.summary.rows == 2
    assert csv_reader.call_count == 1
    assert not stream.closed
    with pytest.raises(ValueError, match=r"(?:once|already|one.pass)"):
        list(reader)
    assert not stream.closed
    assert reader.summary.rows == 2


def test_abandoned_iteration_does_not_close_stream_or_allow_restart(tmp_path: Path) -> None:
    stream = io.StringIO(HEADER + "s\te\tinc\t2\ns\te\tdec\t1\n")
    reader = profile_count_reader(tmp_path / "absent.tsv", text_stream=stream)
    iterator = iter(reader)
    assert next(iterator).counts == (2,)
    assert isinstance(iterator, GeneratorType)
    iterator.close()
    assert not stream.closed
    with pytest.raises(ValueError, match=r"(?:once|already|one.pass)"):
        list(reader)


def test_taxon_stream_keeps_taxon_only_key_and_metadata(tmp_path: Path) -> None:
    stream = io.StringIO(
        "taxon_id\ttaxon_name\ttaxonomic_rank\ttrait_a\nNCBITaxon:1\tGenus species\tspecies\t2\n"
    )
    reader = taxa_trait_count_reader(tmp_path / "absent.tsv", text_stream=stream)
    row = next(iter(reader))
    assert row.key == ("NCBITaxon:1",)
    assert row.values[1:3] == ("Genus species", "species")
    assert row.counts == (2,)
    assert not stream.closed


def test_raw_buffer_and_text_stream_are_mutually_exclusive(tmp_path: Path) -> None:
    raw = (HEADER + "s\te\tinc\t2\n").encode()
    stream = io.StringIO(raw.decode())
    with pytest.raises(ValueError, match=r"raw_bytes.*text_stream|text_stream.*raw_bytes"):
        profile_count_reader(tmp_path / "absent.tsv", raw_bytes=raw, text_stream=stream)
    assert stream.tell() == 0
    assert not stream.closed


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("", "empty"),
        ("study_id\texperiment_id\tinc_dec\ttrait_a\ttrait_a\n", "duplicate column"),
        ("wrong\texperiment_id\tinc_dec\ttrait_a\n", "must start"),
        ("study_id\texperiment_id\tinc_dec\t\n", "blank header"),
        (HEADER + "s\te\tinc\t\n", "Missing count"),
        (HEADER + "s\te\tinc\t-1\n", "finite nonnegative integer"),
        (HEADER + "s\te\tinc\t1.5\n", "finite nonnegative integer"),
        (HEADER + "s\te\tinc\tNaN\n", "finite nonnegative integer"),
        (HEADER + "s\te\tinc\t\uff11\n", "finite nonnegative integer"),
        (HEADER + "s\te\tinc\n", "Wrong row width"),
        (HEADER + "s\te\tinc\t1\textra\n", "Wrong row width"),
        (HEADER + "s\t\tinc\t1\n", "Missing identifier"),
        (HEADER + "s\te\tinc\t1\ns\te\tinc\t2\n", "Duplicate key"),
    ],
)
def test_stream_mode_keeps_validation_and_caller_ownership(
    tmp_path: Path, text: str, message: str
) -> None:
    stream = io.StringIO(text)
    with pytest.raises(ValueError, match=message):
        list(profile_count_reader(tmp_path / "nominal.tsv", text_stream=stream))
    assert not stream.closed


def test_explicit_columns_and_sparse_zero_mode_work_on_the_same_stream(tmp_path: Path) -> None:
    stream = _NonseekableText(HEADER.rstrip("\n") + "\tnote\ns\te\tinc\t\tkept\n")
    reader = CountTableReader(
        tmp_path / "absent.tsv",
        leading_columns=PROFILE_KEY_COLUMNS,
        key_columns=PROFILE_KEY_COLUMNS,
        count_columns=("trait_a",),
        exact_columns=(*PROFILE_KEY_COLUMNS, "trait_a", "note"),
        mode=LEGACY_SPARSE_ZERO_MODE,
        text_stream=stream,
    )
    rows = list(reader)
    assert rows[0].counts == (0,)
    assert rows[0].values[-1] == "kept"
    assert reader.summary.blank_counts_interpreted_as_zero == 1
    assert not stream.closed


def test_stream_exact_columns_and_taxon_uniqueness_are_still_enforced(tmp_path: Path) -> None:
    stream = io.StringIO(HEADER + "s\te\tinc\t1\n")
    with pytest.raises(ValueError, match="exactly these columns"):
        CountTableReader(
            tmp_path / "absent.tsv",
            leading_columns=PROFILE_KEY_COLUMNS,
            key_columns=PROFILE_KEY_COLUMNS,
            exact_columns=(*PROFILE_KEY_COLUMNS, "other"),
            text_stream=stream,
        )
    taxa = io.StringIO(
        "taxon_id\ttaxon_name\ttaxonomic_rank\ttrait_a\n"
        "NCBITaxon:1\tOne\tspecies\t1\nNCBITaxon:1\tTwo\tgenus\t2\n"
    )
    with pytest.raises(ValueError, match="Duplicate key"):
        list(taxa_trait_count_reader(tmp_path / "absent.tsv", text_stream=taxa))
    assert not stream.closed
    assert not taxa.closed
