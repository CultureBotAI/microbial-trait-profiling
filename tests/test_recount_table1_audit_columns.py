"""Synthetic stdout and numerical-preservation checks for the #1430 recount labels."""

import csv
import importlib.util
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from scipy import stats

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts/adhoc/recount_table1_audit_columns_807.py"
)
SPEC = importlib.util.spec_from_file_location("recount_table1_audit_columns", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
recount_module = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = recount_module
SPEC.loader.exec_module(recount_module)

TRAITS = (
    ("METPO:1000608", "facultatively aerobic"),
    ("METPO:2000103|GO:0030245", "cellulose catabolic process"),
    ("METPO:1000668", "coccus shaped"),
    ("METPO:1000683", "sphere shaped"),
    ("METPO:2000103|pathways:arsenate_reduction", "arsenate reduction"),
    ("METPO:1000623", "moderately halophilic"),
)
COCCUS_VALUES = [0, 1, 2, 3, 4]
SPHERE_VALUES = [1, 0, 4, 2, 5]


def _write_tsv(path: Path, rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        csv.writer(stream, delimiter="\t").writerows(rows)


def _inputs(tmp_path: Path, *, constant: bool = False) -> tuple[Path, Path, Path]:
    """Four matrix taxa and five profiles; every required trait is explicit."""
    matrix = tmp_path / "matrix.tsv"
    profiles = tmp_path / "profiles.tsv"
    benchmark = tmp_path / "benchmark.tsv"
    _write_tsv(
        matrix,
        [
            ["taxon_id", *(trait for trait, _ in TRAITS)],
            ["taxon-a", 1, 0, 1, 1, 0, 0],
            ["taxon-b", 0, 1, 1, 0, 0, 0],
            ["taxon-c", 1, 1, 0, 1, 1, 0],
            ["taxon-d", 0, 0, 0, 0, 0, 1],
        ],
    )
    coccus = [1] * 5 if constant else COCCUS_VALUES
    sphere = [1] * 5 if constant else SPHERE_VALUES
    other_traits = [
        (0, 1, 0, 0),
        (1, 0, 0, 1),
        (0, 2, 0, 0),
        (1, 0, 0, 0),
        (1, 1, 1, 0),
    ]
    _write_tsv(
        profiles,
        [
            ["study_id", *(trait for trait, _ in TRAITS)],
            *[
                [
                    f"study-{index}",
                    facultative,
                    cellulose,
                    coccus[index],
                    sphere[index],
                    arsenate,
                    halophile,
                ]
                for index, (facultative, cellulose, arsenate, halophile) in enumerate(other_traits)
            ],
        ],
    )
    _write_tsv(
        benchmark,
        [
            ["trait_id", "matrix_taxon_count"],
            *[[trait, count] for (trait, _), count in zip(TRAITS, (2, 2, 2, 2, 1, 1), strict=True)],
        ],
    )
    return matrix, profiles, benchmark


def _count_lines(*, constant: bool = False) -> list[str]:
    """Expected counts are hand-counted from the rows above, not the producer."""
    lines = ["matrix taxa (rows): 4"]
    for (_, label), count in zip(TRAITS, (2, 2, 2, 2, 1, 1), strict=True):
        lines.append(
            f"  {label:28s} K={count:5d} ({25 * count:.1f}%)  benchmark matrix_taxon_count={count}"
        )
    lines.extend(
        [
            "coccus K=2 sphere K=2 intersection=1 union=3 sphere-in-coccus=50.0% "
            "coccus-in-sphere=50.0% Jaccard=0.3333",
            "profiles: 5",
        ]
    )
    counts = (3, 3, 5, 5, 1, 1) if constant else (3, 3, 4, 4, 1, 1)
    for (_, label), count in zip(TRAITS, counts, strict=True):
        lines.append(f"  {label:28s} nonzero={count:5d} ({20 * count:.1f}%)")
    return lines


def _summary(finite: int) -> str:
    return (
        f"Correlation comparisons: 2; finite nominal p-values: {finite}; "
        f"unavailable p-values: {2 - finite}; no multiplicity correction."
    )


def test_recount_preserves_counts_and_finite_statistics_with_nominal_labels(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths = _inputs(tmp_path)
    original_bytes = [path.read_bytes() for path in paths]
    pearson = stats.pearsonr(COCCUS_VALUES, SPHERE_VALUES)
    spearman = stats.spearmanr(COCCUS_VALUES, SPHERE_VALUES)
    pearsonr, spearmanr = stats.pearsonr, stats.spearmanr
    calls: list[tuple[str, list[int], list[int]]] = []

    def actual_pearson(left: pd.Series, right: pd.Series) -> Any:
        calls.append(("pearson", left.tolist(), right.tolist()))
        return pearsonr(left, right)

    def actual_spearman(left: pd.Series, right: pd.Series) -> Any:
        calls.append(("spearman", left.tolist(), right.tolist()))
        return spearmanr(left, right)

    monkeypatch.setattr(recount_module.stats, "pearsonr", actual_pearson)
    monkeypatch.setattr(recount_module.stats, "spearmanr", actual_spearman)
    recount_module.recount(*paths)

    assert calls == [
        ("pearson", COCCUS_VALUES, SPHERE_VALUES),
        ("spearman", COCCUS_VALUES, SPHERE_VALUES),
    ]
    assert [path.read_bytes() for path in paths] == original_bytes
    assert capsys.readouterr().out.splitlines() == [
        *_count_lines(),
        f"coccus vs sphere profile counts: Pearson r={pearson[0]:.3f} "
        f"(nominal p={pearson[1]:.3g}); Spearman rho={spearman[0]:.3f} "
        f"(nominal p={spearman[1]:.3g})",
        _summary(2),
    ]


def test_recount_constant_profiles_report_unavailable_p_without_test_status(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    paths = _inputs(tmp_path, constant=True)
    original_bytes = [path.read_bytes() for path in paths]
    with pytest.warns(stats.ConstantInputWarning) as caught:
        recount_module.recount(*paths)

    assert len(caught) == 2
    output = capsys.readouterr().out
    assert output.splitlines() == [
        *_count_lines(constant=True),
        "coccus vs sphere profile counts: Pearson r=nan (p unavailable); "
        "Spearman rho=nan (p unavailable)",
        _summary(0),
    ]
    assert "not tested" not in output.casefold()
    assert [path.read_bytes() for path in paths] == original_bytes


@pytest.mark.parametrize(
    "unavailable", [None, float("nan"), float("inf")], ids=["none", "nan", "infinity"]
)
@pytest.mark.parametrize("unavailable_method", ["pearson", "spearman"])
def test_recount_mixed_availability_counts_finite_zero(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
    unavailable: float | None,
    unavailable_method: str,
) -> None:
    """Injected result tuples test rendering, not a synthetic scientific finding."""
    paths = _inputs(tmp_path)
    original_bytes = [path.read_bytes() for path in paths]
    pearson_p = unavailable if unavailable_method == "pearson" else 0.0
    spearman_p = unavailable if unavailable_method == "spearman" else 0.0
    calls: list[tuple[str, list[int], list[int]]] = []

    def pearson_boundary(left: pd.Series, right: pd.Series) -> tuple[float, float | None]:
        calls.append(("pearson", left.tolist(), right.tolist()))
        return 0.25, pearson_p

    def spearman_boundary(left: pd.Series, right: pd.Series) -> tuple[float, float | None]:
        calls.append(("spearman", left.tolist(), right.tolist()))
        return -0.875, spearman_p

    monkeypatch.setattr(recount_module.stats, "pearsonr", pearson_boundary)
    monkeypatch.setattr(recount_module.stats, "spearmanr", spearman_boundary)
    try:
        recount_module.recount(*paths)
    except (TypeError, ValueError) as error:
        raise AssertionError("unavailable p-values must remain renderable") from error

    pearson_label = "p unavailable" if unavailable_method == "pearson" else "nominal p=0"
    spearman_label = "p unavailable" if unavailable_method == "spearman" else "nominal p=0"
    output = capsys.readouterr().out
    assert calls == [
        ("pearson", COCCUS_VALUES, SPHERE_VALUES),
        ("spearman", COCCUS_VALUES, SPHERE_VALUES),
    ]
    assert output.splitlines() == [
        *_count_lines(),
        f"coccus vs sphere profile counts: Pearson r=0.250 ({pearson_label}); "
        f"Spearman rho=-0.875 ({spearman_label})",
        _summary(1),
    ]
    assert "not tested" not in output.casefold()
    assert [path.read_bytes() for path in paths] == original_bytes
