"""Shared count-table contract at every statistical input boundary."""

import gzip
from pathlib import Path

import pytest

from mtp.count_input import (
    LEGACY_SPARSE_ZERO_MODE,
    PROFILE_KEY_COLUMNS,
    CountTableReader,
    profile_count_reader,
    taxa_trait_count_reader,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_complete_profile_counts_preserve_explicit_zero(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "profiles.tsv",
        "study_id\texperiment_id\tinc_dec\ttrait_a\ttrait_b\ns1\te1\tinc\t0\t2\n",
    )
    reader = profile_count_reader(path)

    rows = list(reader)

    assert reader.count_columns == ("trait_a", "trait_b")
    assert rows[0].key == ("s1", "e1", "inc")
    assert rows[0].counts == (0, 2)
    assert reader.summary.blank_counts_interpreted_as_zero == 0


@pytest.mark.parametrize(
    ("body", "message"),
    [
        ("s1\te1\tinc\t1\n", "Wrong row width"),
        ("s1\te1\tinc\t1\t2\t3\n", "Wrong row width"),
        ("s1\te1\tinc\t\t2\n", "Missing count"),
        ("s1\te1\tinc\t-1\t2\n", "finite nonnegative integer"),
        ("s1\te1\tinc\t1.5\t2\n", "finite nonnegative integer"),
        ("s1\te1\tinc\tNaN\t2\n", "finite nonnegative integer"),
        ("s1\te1\tinc\tword\t2\n", "finite nonnegative integer"),
        ("s1\te1\tinc\t\uff11\uff12\t2\n", "finite nonnegative integer"),
        ("s1\te1\tinc\t\u0661\t2\n", "finite nonnegative integer"),
        ("s1\te1\tinc\t+1\t2\n", "finite nonnegative integer"),
        ("s1\te1\tinc\tinf\t2\n", "finite nonnegative integer"),
        ("s1\t\tinc\t1\t2\n", "Missing identifier"),
        ("\n", "Wrong row width"),
    ],
)
@pytest.mark.parametrize("buffered", [False, True])
def test_profile_count_defects_are_not_coerced(
    tmp_path: Path, body: str, message: str, buffered: bool
) -> None:
    path = _write(
        tmp_path / "profiles.tsv",
        "study_id\texperiment_id\tinc_dec\ttrait_a\ttrait_b\n" + body,
    )
    with pytest.raises(ValueError, match=message):
        list(
            profile_count_reader(path, raw_bytes=path.read_bytes())
            if buffered
            else profile_count_reader(path)
        )


@pytest.mark.parametrize("buffered", [False, True])
def test_duplicate_columns_and_profile_keys_are_rejected(tmp_path: Path, buffered: bool) -> None:
    duplicate_columns = _write(
        tmp_path / "duplicate_columns.tsv",
        "study_id\texperiment_id\tinc_dec\ttrait_a\ttrait_a\ns1\te1\tinc\t1\t2\n",
    )
    with pytest.raises(ValueError, match="duplicate column 'trait_a'"):
        list(
            profile_count_reader(duplicate_columns, raw_bytes=duplicate_columns.read_bytes())
            if buffered
            else profile_count_reader(duplicate_columns)
        )

    duplicate_keys = _write(
        tmp_path / "duplicate_keys.tsv",
        "study_id\texperiment_id\tinc_dec\ttrait_a\ns1\te1\tinc\t1\ns1\te1\tinc\t2\n",
    )
    with pytest.raises(ValueError, match=r"Duplicate key.*first seen at row 2"):
        list(
            profile_count_reader(duplicate_keys, raw_bytes=duplicate_keys.read_bytes())
            if buffered
            else profile_count_reader(duplicate_keys)
        )


def test_duplicate_taxon_id_is_rejected_even_when_metadata_differs(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "taxa.tsv",
        "taxon_id\ttaxon_name\ttaxonomic_rank\ttrait_a\n"
        "NCBITaxon:1\tFirst name\tspecies\t1\n"
        "NCBITaxon:1\tConflicting name\tgenus\t2\n",
    )
    with pytest.raises(ValueError, match=r"Duplicate key.*NCBITaxon:1"):
        list(taxa_trait_count_reader(path))


def test_legacy_sparse_mode_is_explicit_and_still_requires_full_width(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "legacy.tsv",
        "study_id\texperiment_id\tinc_dec\ttrait_a\ttrait_b\ns1\te1\tinc\t\t2\n",
    )
    reader = profile_count_reader(path, mode=LEGACY_SPARSE_ZERO_MODE)
    assert next(iter(reader)).counts == (0, 2)
    assert reader.summary.blank_counts_interpreted_as_zero == 1

    truncated = _write(
        tmp_path / "truncated.tsv",
        "study_id\texperiment_id\tinc_dec\ttrait_a\ttrait_b\ns1\te1\tinc\t1\n",
    )
    with pytest.raises(ValueError, match="Wrong row width"):
        list(profile_count_reader(truncated, mode=LEGACY_SPARSE_ZERO_MODE))


def test_gzip_profile_input_uses_the_same_contract(tmp_path: Path) -> None:
    path = tmp_path / "profiles.tsv.gz"
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("study_id\texperiment_id\tinc_dec\ttrait_a\ns1\te1\tinc\t3\n")
    assert next(iter(profile_count_reader(path))).counts == (3,)


@pytest.mark.parametrize("compressed", [False, True])
def test_profile_buffer_preserves_header_keys_and_counts_without_reopening_path(
    tmp_path: Path, compressed: bool
) -> None:
    raw = b"study_id\texperiment_id\tinc_dec\ttrait_a\ttrait_b\ns1\te1\tinc\t0\t2\n"
    decoy = b"study_id\texperiment_id\tinc_dec\tother\ns9\te9\tdec\t99\n"
    path = tmp_path / ("profiles.tsv.gz" if compressed else "profiles.tsv")
    encoded = gzip.compress(raw, mtime=0) if compressed else raw
    path.write_bytes(gzip.compress(decoy, mtime=0) if compressed else decoy)
    checked = profile_count_reader(path, raw_bytes=encoded)
    # A captured buffer remains usable when the original location has disappeared.
    path.unlink()
    result = list(checked)
    assert checked.count_columns == ("trait_a", "trait_b")
    assert [row.key for row in result] == [("s1", "e1", "inc")]
    assert [row.counts for row in result] == [(0, 2)]
    assert checked.summary.blank_counts_interpreted_as_zero == 0
    # Repeat iteration must reset its counters and reuse the same immutable input.
    assert [row.counts for row in checked] == [(0, 2)]
    assert checked.summary.rows == 1


def test_generic_count_reader_buffer_keeps_explicit_column_contract(tmp_path: Path) -> None:
    raw = b"study_id\texperiment_id\tinc_dec\ttrait_a\tnote\ns1\te1\tinc\t4\tkept\n"
    reader = CountTableReader(
        tmp_path / "no-file.tsv",
        leading_columns=PROFILE_KEY_COLUMNS,
        key_columns=PROFILE_KEY_COLUMNS,
        count_columns=("trait_a",),
        exact_columns=(*PROFILE_KEY_COLUMNS, "trait_a", "note"),
        raw_bytes=raw,
    )
    result = list(reader)
    assert result[0].counts == (4,)
    assert result[0].values[-1] == "kept"
    with pytest.raises(ValueError, match="exactly these columns"):
        CountTableReader(
            tmp_path / "no-file.tsv",
            leading_columns=PROFILE_KEY_COLUMNS,
            key_columns=PROFILE_KEY_COLUMNS,
            exact_columns=(*PROFILE_KEY_COLUMNS, "trait_a"),
            raw_bytes=raw,
        )


def test_buffer_keeps_explicit_legacy_sparse_mode(tmp_path: Path) -> None:
    raw = b"study_id\texperiment_id\tinc_dec\ttrait_a\ttrait_b\ns1\te1\tinc\t\t2\n"
    path = tmp_path / "no-file.tsv"
    with pytest.raises(ValueError, match="Missing count"):
        list(profile_count_reader(path, raw_bytes=raw))
    checked = profile_count_reader(path, raw_bytes=raw, mode=LEGACY_SPARSE_ZERO_MODE)
    assert [row.counts for row in checked] == [(0, 2)]
    assert checked.summary.blank_counts_interpreted_as_zero == 1
