"""Keep nominal ORA p-values separate from unavailable adjusted values (#1409)."""

from __future__ import annotations

import csv
import io
import math
from pathlib import Path

import pytest

from mtp.enrichment_results import RESULT_COLUMNS, MalformedResultsError, load_results_tsv


def _results(
    tmp_path: Path,
    *,
    fdr: str | None = "False",
    scope: str = "none",
    q_values: tuple[str, str] = ("", ""),
    family: str | None = None,
) -> Path:
    """Write two selected tests without relying on the production ORA writer."""
    path = tmp_path / "results.tsv"
    metadata = [
        f"# FDR-Scope: {scope}",
        "# Total-Tests-Performed: 2",
        "# Significant-Hits: 2",
        "# Result-Status: discoveries",
    ]
    if fdr is not None:
        metadata.insert(0, f"# FDR: {fdr}")
    columns = [*RESULT_COLUMNS, *(["correction_family"] if family is not None else [])]
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\n".join(metadata) + "\n")
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for direction, p_value, q_value in zip(
            ("inc", "dec"), ("0.01", "0.2"), q_values, strict=True
        ):
            row = {
                "study_id": "study",
                "experiment_id": "experiment",
                "inc_dec": direction,
                "trait": "TraitA",
                "p_value": p_value,
                "fdr_q": q_value,
                "fold": "2.0",
                "count": "2",
                "expected": "1.0",
                "testable_profiles": "2",
            }
            if family is not None:
                row["correction_family"] = family
            writer.writerow(row)
    return path


def _assert_nominal_values(path: Path) -> None:
    table = load_results_tsv(path)
    assert table.status == "discoveries"
    assert table.frame["p_value"].tolist() == [0.01, 0.2]
    assert all(isinstance(value, float) and math.isnan(value) for value in table.frame["fdr_q"])
    exported = list(
        csv.DictReader(io.StringIO(table.frame.to_csv(sep="\t", index=False)), delimiter="\t")
    )
    assert [row["fdr_q"] for row in exported] == ["", ""]
    assert [float(row["p_value"]) for row in exported] == [0.01, 0.2]


@pytest.mark.parametrize("family", [None, "none"])
def test_nominal_reader_accepts_explicit_blank_adjusted_values(
    tmp_path: Path, family: str | None
) -> None:
    _assert_nominal_values(_results(tmp_path, family=family))


@pytest.mark.parametrize("q_value", ["0.01", "0", "not-a-number", "nan", "NA", "inf", " "])
def test_nominal_reader_refuses_nonblank_adjusted_values(tmp_path: Path, q_value: str) -> None:
    path = _results(tmp_path, q_values=(q_value, ""), family="none")
    with pytest.raises(MalformedResultsError, match="fdr_q"):
        load_results_tsv(path)


@pytest.mark.parametrize("fdr", ["True", None])
def test_nominal_reader_does_not_allow_blank_q_without_no_fdr_declaration(
    tmp_path: Path, fdr: str | None
) -> None:
    path = _results(tmp_path, fdr=fdr, scope="global", family="global")
    with pytest.raises(MalformedResultsError):
        load_results_tsv(path)


@pytest.mark.parametrize("scope", ["global", "profile"])
@pytest.mark.parametrize("complete", [False, True])
def test_nominal_reader_legacy_copied_p_is_not_exposed_as_q(
    tmp_path: Path, scope: str, complete: bool
) -> None:
    _assert_nominal_values(
        _results(
            tmp_path,
            scope=scope,
            q_values=("0.01", "0.2"),
            family=scope if complete else None,
        )
    )


@pytest.mark.parametrize("scope", ["global", "profile"])
def test_nominal_reader_legacy_mode_requires_exact_copied_p(tmp_path: Path, scope: str) -> None:
    path = _results(tmp_path, scope=scope, q_values=("0.02", "0.2"), family=scope)
    with pytest.raises(MalformedResultsError, match="fdr_q"):
        load_results_tsv(path)


@pytest.mark.parametrize(
    ("scope", "q_values", "family"),
    [
        ("none", ("", ""), "global"),
        ("none", ("", ""), "profile"),
        ("none", ("", ""), "unknown"),
        ("global", ("0.01", "0.2"), "none"),
        ("profile", ("0.01", "0.2"), "global"),
    ],
)
def test_nominal_reader_refuses_inconsistent_family(
    tmp_path: Path, scope: str, q_values: tuple[str, str], family: str
) -> None:
    path = _results(tmp_path, scope=scope, q_values=q_values, family=family)
    with pytest.raises(MalformedResultsError, match="correction_family"):
        load_results_tsv(path)


@pytest.mark.parametrize("fdr", ["True", None])
def test_nominal_reader_preserves_existing_adjusted_results(
    tmp_path: Path, fdr: str | None
) -> None:
    path = _results(tmp_path, fdr=fdr, scope="global", q_values=("0.02", "0.4"), family="global")
    table = load_results_tsv(path)
    assert table.frame["p_value"].tolist() == [0.01, 0.2]
    assert table.frame["fdr_q"].tolist() == [0.02, 0.4]
