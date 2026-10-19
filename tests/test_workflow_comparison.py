"""Migration comparisons preserve row and qualified feature identity (#629)."""

from pathlib import Path

from mtp.workflow.comparison import compare_count_tables


def test_column_order_does_not_create_numeric_drift(tmp_path: Path) -> None:
    before, after = tmp_path / "old.tsv", tmp_path / "new.tsv"
    before.write_text("study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2\ns\te\tinc\t3\t4\n")
    after.write_text("study_id\texperiment_id\tinc_dec\tMETPO:2\tMETPO:1\ns\te\tinc\t4\t3\n")
    comparison = compare_count_tables(before, after)
    assert comparison["identical_counts"] is True
    assert comparison["compared_cells"] == 2


def test_comparison_exposes_new_rows_features_and_media_differences(tmp_path: Path) -> None:
    before, after = tmp_path / "old.tsv", tmp_path / "new.tsv"
    before.write_text(
        "study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2|mediadive.medium:1\n"
        "s\te\tinc\t3\t4\ns\te\tdec\t1\t0\n"
    )
    after.write_text(
        "study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2|mediadive.medium:1\tMETPO:3\n"
        "s\te\tinc\t3\t5\t1\ns2\te\tinc\t1\t0\t0\n"
    )
    comparison = compare_count_tables(before, after)
    assert comparison["identical_counts"] is False
    assert comparison["changed_cells"] == 1
    assert comparison["changed_nonmedia_cells"] == 0
    assert comparison["added_rows"] == [["s2", "e", "inc"]]
    assert comparison["removed_rows"] == [["s", "e", "dec"]]
    assert comparison["added_features"] == ["METPO:3"]
