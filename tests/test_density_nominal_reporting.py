"""Density reports label uncorrected p-values and count finite tests (#1430, R3).

The real module and consolidated writers run on 24 synthetic profiles and six
traits. Only the measured effects and twin summary are replaced to exercise
missing and nonfinite p-values without reading any derived scientific data.
"""

from __future__ import annotations

import importlib.util
import json
import math
from collections.abc import Callable
from dataclasses import asdict
from pathlib import Path
from types import ModuleType
from typing import Any
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from mtp import density_adjustment as density

ReportRenderer = Callable[[list[density.Table1TraitEffect], float], tuple[str, dict[str, Any]]]


@pytest.fixture(scope="module")
def audit_script() -> ModuleType:
    repo_root = Path(__file__).resolve().parents[1]
    assert density.__file__ is not None
    assert Path(density.__file__).resolve() == repo_root / "src/mtp/density_adjustment.py"
    path = repo_root / "scripts/audit_trait_density_correlation.py"
    spec = importlib.util.spec_from_file_location("density_nominal_under_test", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.__file__ is not None
    assert Path(module.__file__).resolve() == path
    return module


def _write_small_profiles(directory: Path) -> tuple[Path, Path]:
    """Write twelve paired experiments over the six required trait columns."""
    traits = list(density.TABLE_1_TRAITS)
    rng = np.random.default_rng(1430)
    rows = []
    for experiment in range(12):
        for direction in ("inc", "dec"):
            counts = rng.integers(1, 40, size=len(traits)).tolist()
            rows.append(
                {
                    "study_id": f"study-{experiment}",
                    "experiment_id": f"experiment-{experiment}",
                    "inc_dec": direction,
                    **dict(zip(traits, counts, strict=True)),
                }
            )
    profiles = directory / "profiles.tsv"
    columns = directory / "columns.tsv"
    pd.DataFrame(rows).to_csv(profiles, sep="\t", index=False)
    pd.DataFrame(
        {
            "column_id": traits,
            "column_label": list(density.TABLE_1_TRAITS.values()),
            "object_name": list(density.TABLE_1_TRAITS.values()),
            "object_category": ["biolink:PhenotypicQuality"] * len(traits),
        }
    ).to_csv(columns, sep="\t", index=False)
    return profiles, columns


def _effects(raw: list[float], adjusted: list[float]) -> list[density.Table1TraitEffect]:
    """Keep all six rows, including untestable ones, and preserve boundary q-values."""
    rows = []
    for index, ((trait_id, label), raw_p, adj_p) in enumerate(
        zip(density.TABLE_1_TRAITS.items(), raw, adjusted, strict=True)
    ):
        raw_q = 0.05004 if math.isfinite(raw_p) else float("nan")
        adj_q = 0.05 if math.isfinite(adj_p) else float("nan")
        raw_mean = 0.125 + index / 100
        adj_mean = raw_mean / 2
        rows.append(
            density.Table1TraitEffect(
                trait_id=trait_id,
                trait_label=label,
                raw_mean_diff=raw_mean,
                raw_wilcoxon_p=raw_p,
                adj_mean_diff=adj_mean,
                adj_wilcoxon_p=adj_p,
                relative_change_pct=-50.0,
                status=density.adjustment_status(raw_mean, raw_q, adj_mean, adj_q),
                magnitude_change_pct=-50.0,
                sign_flipped=False,
                raw_wilcoxon_bh_q=raw_q,
                adj_wilcoxon_bh_q=adj_q,
                bh_family_size=len(raw),
            )
        )
    return rows


@pytest.fixture(params=["module", "consolidated"])
def render_report(
    request: pytest.FixtureRequest,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    audit_script: ModuleType,
) -> ReportRenderer:
    profiles, columns = _write_small_profiles(tmp_path)

    def render(
        effects: list[density.Table1TraitEffect], slope_p: float
    ) -> tuple[str, dict[str, Any]]:
        twin = density.ConditionalTwinSummary(
            trait_a=density.TWIN_A,
            trait_b=density.TWIN_B,
            trait_a_label="coccus shaped",
            trait_b_label="sphere shaped",
            profile_correlation=0.4,
            profile_partial_r_depth=0.3,
            vif_profile=1.2,
            profile_slope_a_on_b=0.375,
            profile_intercept_a_on_b=0.1,
            n_paired_experiments=23,
            paired_correlation=-0.3,
            vif_paired=1.1,
            paired_slope_a_on_b=-0.3125,
            paired_slope_pvalue=slope_p,
        )
        expected = pd.DataFrame([asdict(effect) for effect in effects])
        monkeypatch.setattr(density, "evaluate_table1_paired_effects", Mock(return_value=effects))
        monkeypatch.setattr(density, "evaluate_conditional_twin_test", Mock(return_value=twin))
        output_dir = tmp_path / "report"
        result = density.run_density_audit(profiles, columns, output_dir)
        if request.param == "consolidated":
            ontology = {
                "n_redundant_pairs": 0,
                "n_clusters": 0,
                "n_collapsed_features": 0,
                "bounded_trait_space": len(effects),
                "native_metrics": result["raw_pearson"],
                "max_metrics": result["raw_pearson"],
                "adj_collapsed_metrics": result["adj_pc1_pearson"],
            }
            report_path = audit_script.write_consolidated_report(output_dir, result, ontology)
        else:
            report_path = result["report_txt"]
        # Rendering must preserve every effect, p, q, status and selected row/order.
        pd.testing.assert_frame_equal(
            pd.DataFrame([asdict(effect) for effect in result["table1_effects"]]), expected
        )
        pd.testing.assert_frame_equal(
            pd.read_csv(result["table1_tsv"], sep="\t"), expected, check_dtype=False
        )
        assert result["twin_summary"] is twin
        assert result["bh_family_size"] == len(effects)
        content = report_path.read_text(encoding="utf-8")
        assert content
        return content, result

    return render


@pytest.mark.parametrize(
    ("p", "q", "expected"),
    [
        (0.0, 0.0, "nominal p=0.00e+00, q=0.00e+00"),
        (0.01, 0.05004, "nominal p=1.00e-02, q=5.004e-02"),
        (1.0, 0.05, "nominal p=1.00e+00, q=5.00e-02"),
        (float("nan"), float("nan"), "nominal p=NA, q=NA"),
        (np.float64(0.01), np.float64(0.0496), "nominal p=1.00e-02, q=4.96e-02"),
    ],
)
def test_p_and_q_labels_nominal_without_changing_numeric_text(
    p: float, q: float, expected: str
) -> None:
    assert density.format_p_and_q(p, q) == expected


def _assert_count_line_placement(
    report: str,
    effects: list[density.Table1TraitEffect],
    raw_count: int,
    adjusted_count: int,
) -> None:
    """#1658: count line sits in section 3 right after the family statement and before row 1."""
    count_line = (
        f"Finite nominal Wilcoxon tests shown: raw={raw_count}; PC1-adjusted={adjusted_count}."
    )
    assert report.count("Finite nominal Wilcoxon tests shown:") == 1
    lines = report.splitlines()
    sec3_idx = next(
        idx
        for idx, line in enumerate(lines)
        if line.startswith(
            (
                "3. SURVIVAL OF CONVERGENT SIX TRAIT SIGNALS",
                "3. DECLARED TRAIT LIST SIGNALS",
            )
        )
    )
    family_end_idx = lines.index("the family at p = 1 and has no q.")
    count_idx = lines.index(count_line)
    first_trait_idx = next(
        idx for idx, line in enumerate(lines) if f"{effects[0].trait_label:<30s} |" in line
    )
    assert sec3_idx < family_end_idx
    assert count_idx == family_end_idx + 1
    assert first_trait_idx == count_idx + 1


@pytest.mark.parametrize(
    ("raw", "adjusted", "raw_count", "adjusted_count"),
    [
        (
            [0.0, 0.01, float("nan"), float("inf"), -float("inf"), float("nan")],
            [float("nan"), 0.0, 0.07, 0.4, float("inf"), float("nan")],
            2,
            3,
        ),
        (
            [0.0, 0.01, 1.0, float("nan"), float("inf"), -float("inf")],
            [0.0, float("nan"), float("inf"), -float("inf"), float("nan"), float("nan")],
            3,
            1,
        ),
        ([float("nan")] * 6, [float("nan")] * 6, 0, 0),
        ([float("nan")] * 5 + [0.03], [float("nan")] * 5 + [0.04], 1, 1),
    ],
)
def test_reports_count_finite_nominal_wilcoxon_values_in_each_displayed_column(
    render_report: ReportRenderer,
    raw: list[float],
    adjusted: list[float],
    raw_count: int,
    adjusted_count: int,
) -> None:
    effects = _effects(raw, adjusted)
    report, result = render_report(effects, 0.0125)
    _assert_count_line_placement(report, effects, raw_count, adjusted_count)
    assert result["q_values"]["p_value_type"] == "nominal"
    provenance = json.loads(result["provenance"].read_text(encoding="utf-8"))
    assert provenance["q_values"]["p_value_type"] == "nominal"
    assert "over all 6 reported traits" in report
    assert "PAIRED EXPERIMENTS N=23" in report
    for effect in effects:
        line = next(line for line in report.splitlines() if f"{effect.trait_label:<30s} |" in line)
        raw_text = "NA" if math.isnan(effect.raw_wilcoxon_p) else f"{effect.raw_wilcoxon_p:.2e}"
        adj_text = "NA" if math.isnan(effect.adj_wilcoxon_p) else f"{effect.adj_wilcoxon_p:.2e}"
        raw_q = "5.004e-02" if math.isfinite(effect.raw_wilcoxon_p) else "NA"
        adj_q = "5.00e-02" if math.isfinite(effect.adj_wilcoxon_p) else "NA"
        assert f"Raw mean: {effect.raw_mean_diff:+.3f} (nominal p={raw_text}, q={raw_q})" in line
        assert f"Adj mean: {effect.adj_mean_diff:+.3f} (nominal p={adj_text}, q={adj_q})" in line
        assert line.endswith(f" | {effect.status}")


@pytest.mark.parametrize(
    ("slope", "expected"),
    [
        (0.0, "0.0000"),
        (-0.3125, "-0.3125"),
        (1.5, "1.5000"),
        (float("nan"), "NA"),
        (float("inf"), "NA"),
        (-float("inf"), "NA"),
        (np.float64("nan"), "NA"),
    ],
)
def test_format_paired_regression_slope(slope: float, expected: str) -> None:
    """#1714: finite paired regression slopes format to 4 decimals; non-finite format as NA."""
    assert density.format_paired_regression_slope(slope) == expected


@pytest.mark.parametrize(
    ("slope_p", "count", "p_text", "test_word"),
    [
        (0.0, 1, "0.00e+00", "test"),
        (0.0125, 1, "1.25e-02", "test"),
        (1.0, 1, "1.00e+00", "test"),
        (float("nan"), 0, "NA", "tests"),
        (float("inf"), 0, "NA", "tests"),
        (-float("inf"), 0, "NA", "tests"),
    ],
)
def test_reports_label_paired_slope_and_count_only_its_finite_test(
    render_report: ReportRenderer, slope_p: float, count: int, p_text: str, test_word: str
) -> None:
    effects = _effects([0.01] * 6, [float("nan")] * 6)
    report, _ = render_report(effects, slope_p)
    _assert_count_line_placement(report, effects, 6, 0)
    line = next(line for line in report.splitlines() if line.startswith("Paired regression slope:"))
    assert line.startswith("Paired regression slope: -0.3125 ")
    assert f"nominal p={p_text}; {count} finite paired-regression {test_word} shown" in line
    assert "; diagnostic" in line
    assert ("untested" in line) is (count == 0)


@pytest.mark.parametrize("writer", ["module", "consolidated"])
@pytest.mark.parametrize("case", ["zero_pairs", "one_pair", "constant_delta_b", "zero_slope"])
def test_unmocked_undefined_vs_zero_paired_slope_report_and_readers(
    tmp_path: Path, audit_script: ModuleType, writer: str, case: str
) -> None:
    """#1714: <2 pairs or constant delta_b give NaN slope printed as NA;
    a true zero slope (constant delta_a with varying delta_b) still prints 0.0000;
    TSV and JSON outputs remain readable."""
    traits = list(density.TABLE_1_TRAITS)
    rows: list[dict[str, object]] = []
    if case == "zero_pairs":
        for idx in range(3):
            base = {trait_id: float(5 + idx) for trait_id in traits}
            rows.append(
                {
                    "study_id": f"study-{idx}",
                    "experiment_id": f"exp-{idx}",
                    "inc_dec": "inc",
                    **base,
                }
            )
    elif case == "one_pair":
        for direction, da, db in (("inc", 4.0, 2.0), ("dec", 0.0, 0.0)):
            base = dict.fromkeys(traits, 5.0)
            base[density.TWIN_A] += da
            base[density.TWIN_B] += db
            rows.append(
                {
                    "study_id": "study-0",
                    "experiment_id": "exp-0",
                    "inc_dec": direction,
                    **base,
                }
            )
    elif case == "constant_delta_b":
        for exp_idx, da in enumerate((2.0, 9.0, 15.0)):
            for direction in ("inc", "dec"):
                base = {trait_id: float(5 + exp_idx) for trait_id in traits}
                if direction == "inc":
                    base[density.TWIN_A] += da
                    base[density.TWIN_B] += 3.0
                rows.append(
                    {
                        "study_id": f"study-{exp_idx}",
                        "experiment_id": f"exp-{exp_idx}",
                        "inc_dec": direction,
                        **base,
                    }
                )
    else:
        for exp_idx, db in enumerate((1.0, 3.0, 6.0)):
            for direction in ("inc", "dec"):
                base = {trait_id: float(5 + exp_idx) for trait_id in traits}
                if direction == "inc":
                    base[density.TWIN_A] += 4.0
                    base[density.TWIN_B] += db
                rows.append(
                    {
                        "study_id": f"study-{exp_idx}",
                        "experiment_id": f"exp-{exp_idx}",
                        "inc_dec": direction,
                        **base,
                    }
                )

    profiles = tmp_path / f"profiles_{writer}_{case}.tsv"
    columns = tmp_path / f"columns_{writer}_{case}.tsv"
    pd.DataFrame(rows).to_csv(profiles, sep="\t", index=False)
    pd.DataFrame(
        {
            "column_id": traits,
            "column_label": list(density.TABLE_1_TRAITS.values()),
            "object_name": list(density.TABLE_1_TRAITS.values()),
            "object_category": ["biolink:PhenotypicQuality"] * len(traits),
        }
    ).to_csv(columns, sep="\t", index=False)

    out_dir = tmp_path / f"out_{writer}_{case}"
    result = density.run_density_audit(profiles, columns, out_dir)
    if writer == "consolidated":
        ontology = {
            "n_redundant_pairs": 0,
            "n_clusters": 0,
            "n_collapsed_features": 0,
            "bounded_trait_space": len(traits),
            "native_metrics": result["raw_pearson"],
            "max_metrics": result["raw_pearson"],
            "adj_collapsed_metrics": result["adj_pc1_pearson"],
        }
        report_path = audit_script.write_consolidated_report(out_dir, result, ontology)
    else:
        report_path = result["report_txt"]

    summary_df = pd.read_csv(result["summary_tsv"], sep="\t")
    table1_df = pd.read_csv(result["table1_tsv"], sep="\t")
    provenance = json.loads(result["provenance"].read_text(encoding="utf-8"))
    assert len(summary_df) == 4
    assert len(table1_df) == len(traits)
    assert provenance["q_values"]["p_value_type"] == "nominal"

    twin = result["twin_summary"]
    report = report_path.read_text(encoding="utf-8")
    slope_line = next(
        line for line in report.splitlines() if line.startswith("Paired regression slope:")
    )
    if case == "zero_slope":
        assert twin.paired_slope_a_on_b == pytest.approx(0.0)
        assert slope_line == (
            "Paired regression slope: 0.0000 "
            "(nominal p=NA; 0 finite paired-regression tests shown; diagnostic; untested)"
        )
    else:
        assert math.isnan(twin.paired_slope_a_on_b)
        assert math.isnan(twin.paired_slope_pvalue)
        assert slope_line == (
            "Paired regression slope: NA "
            "(nominal p=NA; 0 finite paired-regression tests shown; diagnostic; untested)"
        )


@pytest.mark.parametrize("writer", ["module", "consolidated"])
@pytest.mark.parametrize("n_pairs", [2, 3])
def test_unmocked_two_vs_three_pairs_paired_regression_p_and_report_line(
    tmp_path: Path, audit_script: ModuleType, writer: str, n_pairs: int
) -> None:
    """#1660: N=2 paired experiments yield NaN regression p (0 finite tests, NA/untested);
    N=3 paired experiments yield a finite regression p (1 finite test)."""
    traits = list(density.TABLE_1_TRAITS)
    rows = []
    # 3 distinct (delta_b, delta_a) pairs that are not collinear so N=3 has 0 < p < 1.
    pair_shifts = [(1.0, 2.0), (3.0, 9.0), (6.0, 11.0)]
    for exp_idx in range(n_pairs):
        db, da = pair_shifts[exp_idx]
        for direction in ("inc", "dec"):
            base = {trait_id: float(5 + exp_idx) for trait_id in traits}
            if direction == "inc":
                base[density.TWIN_A] += da
                base[density.TWIN_B] += db
            rows.append(
                {
                    "study_id": f"study-{exp_idx}",
                    "experiment_id": f"exp-{exp_idx}",
                    "inc_dec": direction,
                    **base,
                }
            )
    profiles = tmp_path / "profiles.tsv"
    columns = tmp_path / "columns.tsv"
    pd.DataFrame(rows).to_csv(profiles, sep="\t", index=False)
    pd.DataFrame(
        {
            "column_id": traits,
            "column_label": list(density.TABLE_1_TRAITS.values()),
            "object_name": list(density.TABLE_1_TRAITS.values()),
            "object_category": ["biolink:PhenotypicQuality"] * len(traits),
        }
    ).to_csv(columns, sep="\t", index=False)

    out_dir = tmp_path / f"out_{writer}_{n_pairs}"
    result = density.run_density_audit(profiles, columns, out_dir)
    if writer == "consolidated":
        ontology = {
            "n_redundant_pairs": 0,
            "n_clusters": 0,
            "n_collapsed_features": 0,
            "bounded_trait_space": len(traits),
            "native_metrics": result["raw_pearson"],
            "max_metrics": result["raw_pearson"],
            "adj_collapsed_metrics": result["adj_pc1_pearson"],
        }
        report_path = audit_script.write_consolidated_report(out_dir, result, ontology)
    else:
        report_path = result["report_txt"]

    twin = result["twin_summary"]
    assert twin.n_paired_experiments == n_pairs
    assert math.isfinite(twin.paired_slope_a_on_b)
    report = report_path.read_text(encoding="utf-8")
    _assert_count_line_placement(report, result["table1_effects"], 0, 0)
    slope_line = next(
        line for line in report.splitlines() if line.startswith("Paired regression slope:")
    )
    if n_pairs == 2:
        assert math.isnan(twin.paired_slope_pvalue)
        assert slope_line == (
            f"Paired regression slope: {twin.paired_slope_a_on_b:.4f} "
            "(nominal p=NA; 0 finite paired-regression tests shown; diagnostic; untested)"
        )
    else:
        assert math.isfinite(twin.paired_slope_pvalue)
        assert 0.0 < twin.paired_slope_pvalue < 1.0
        assert slope_line == (
            f"Paired regression slope: {twin.paired_slope_a_on_b:.4f} "
            f"(nominal p={twin.paired_slope_pvalue:.2e}; "
            "1 finite paired-regression test shown; diagnostic)"
        )
