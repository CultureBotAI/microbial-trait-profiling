"""Issue #1081: assembly preserves nullable counts without altering numerical estimates."""

import csv
import io
from pathlib import Path

import pandas as pd
import pytest

from mtp.three_test_table.carriers import carrier_table
from mtp.three_test_table.native import native_tables
from mtp.three_test_table.ora import ora_table
from mtp.three_test_table.stage import _merge_tables
from multiplicity_fixture import fixture_sources


def test_final_table_merge_preserves_nullable_count_types_and_tsv_blanks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path, monkeypatch)
    native, _ = native_tables(sources)
    ora, _ = ora_table(sources)
    carriers, _ = carrier_table(sources, None)
    result = _merge_tables(sources.features["count_profiles"], [native, ora, carriers])
    counts = (
        "ora_enriched_inc",
        "ora_enriched_dec",
        "ora_depleted_inc",
        "ora_depleted_dec",
        "native_selected_study_count",
        "native_selected_experiment_count",
        "native_selected_profile_support",
        "carrier_direct_reported_taxa",
        "carrier_inherited_reported_taxa",
        "carrier_unique_taxa",
        "carrier_source_record_count",
    )
    assert all(str(result[column].dtype) == "Int64" for column in counts)
    for component in (native, ora, carriers):
        pd.testing.assert_frame_equal(result[component.columns], component)
    exported = {
        row["trait"]: row
        for row in csv.DictReader(io.StringIO(result.to_csv(sep="\t", index=False)), delimiter="\t")
    }
    assert all(
        row[column] == "" or row[column].isdecimal()
        for row in exported.values()
        for column in counts
    )
    for column in counts[4:]:
        assert exported["qualified|extra"][column] == ""
