"""Density reports label uncorrected p-values and count finite tests (#1430, R3).

The real module and consolidated writers run on 24 synthetic profiles and six
traits. Only the measured effects and twin summary are replaced to exercise
missing and nonfinite p-values without reading any derived scientific data.
"""

from __future__ import annotations

import importlib.util
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
    path = Path(__file__).resolve().parents[1] / "scripts/audit_trait_density_correlation.py"
    spec = importlib.util.spec_from_file_location("density_nominal_under_test", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
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
    report, _ = render_report(effects, 0.0125)
    assert (
        f"Finite nominal Wilcoxon tests shown: raw={raw_count}; PC1-adjusted={adjusted_count}."
    ) in report
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
    ("slope_p", "count"),
    [(0.0, 1), (0.0125, 1), (1.0, 1), (float("nan"), 0), (float("inf"), 0), (-float("inf"), 0)],
)
def test_reports_label_paired_slope_and_count_only_its_finite_test(
    render_report: ReportRenderer, slope_p: float, count: int
) -> None:
    effects = _effects([0.01] * 6, [float("nan")] * 6)
    report, _ = render_report(effects, slope_p)
    line = next(line for line in report.splitlines() if line.startswith("Paired regression slope:"))
    assert line.startswith("Paired regression slope: -0.3125 ")
    assert f"nominal p={slope_p:.2e}; {count} finite paired-regression tests shown" in line
    assert "; diagnostic" in line
    assert ("untested" in line) is (count == 0)
