"""Tests for trait coupling sensitivity, LOTO, LCO, and decision rule (Issue #260)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mtp.coupling import (
    COCCUS_ID,
    SPHERE_ID,
    compute_composition_arm,
    compute_directional_arm,
    compute_ora_arm,
    evaluate_decision_rule,
    evaluate_lco,
    evaluate_loto,
    run_coupling_audit,
)
from mtp.density_adjustment import (
    DEFAULT_COLUMNS_PATH,
    DEFAULT_PROFILES_PATH,
)


@pytest.fixture
def synthetic_coupling_data() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Create a minimal synthetic profile dataset for testing coupling mechanics."""
    traits = ["t1", "t2", "t3", "t4"]
    # 4 experiments across 2 studies
    rows = []
    for s_idx, study in enumerate(["study_A", "study_B"]):
        for e_idx, exp in enumerate([f"exp_{s_idx}_1", f"exp_{s_idx}_2"]):
            # inc condition
            rows.append(
                {
                    "study_id": study,
                    "experiment_id": exp,
                    "inc_dec": "inc",
                    "t1": 10.0 + e_idx,
                    "t2": 20.0 + e_idx,
                    "t3": 30.0,
                    "t4": 40.0,
                }
            )
            # dec condition
            rows.append(
                {
                    "study_id": study,
                    "experiment_id": exp,
                    "inc_dec": "dec",
                    "t1": 5.0 + e_idx,
                    "t2": 15.0 + e_idx,
                    "t3": 30.0,
                    "t4": 50.0,
                }
            )
    profiles_df = pd.DataFrame(rows)
    return profiles_df, profiles_df.copy(), traits


def test_compute_composition_arm(
    synthetic_coupling_data: tuple[pd.DataFrame, pd.DataFrame, list[str]],
) -> None:
    """Verify composition arm computation on synthetic data."""
    profiles_df, _, traits = synthetic_coupling_data
    focal = ["t1", "t2"]

    res_all = compute_composition_arm(profiles_df, traits, focal)
    assert len(res_all) == 2
    assert "t1" in res_all
    assert "t2" in res_all

    # t1 should have positive effect because inc > dec
    assert res_all["t1"].primary_stat > 0.0

    # Test LOTO: drop t3 from denominator
    traits_no_t3 = ["t1", "t2", "t4"]
    res_no_t3 = compute_composition_arm(profiles_df, traits_no_t3, focal)
    assert len(res_no_t3) == 2

    # Denominator shrunk, so relative fractions changed
    assert res_no_t3["t1"].primary_stat != pytest.approx(res_all["t1"].primary_stat, abs=1e-12)


def test_compute_directional_arm(
    synthetic_coupling_data: tuple[pd.DataFrame, pd.DataFrame, list[str]],
) -> None:
    """Verify directional sign test calculation."""
    profiles_df, _, traits = synthetic_coupling_data
    focal = ["t1", "t4"]

    res = compute_directional_arm(profiles_df, traits, focal)
    assert len(res) == 2

    # In synthetic data, t1 is always higher in inc -> share_up = 1.0
    assert res["t1"].primary_stat == 1.0
    # In synthetic data, t4 is lower in inc -> share_up = 0.0
    assert res["t4"].primary_stat == 0.0


def test_compute_ora_arm(
    synthetic_coupling_data: tuple[pd.DataFrame, pd.DataFrame, list[str]],
) -> None:
    """Verify ORA hypergeometric calculation."""
    _, paired_df, traits = synthetic_coupling_data
    focal = ["t1", "t2"]

    res = compute_ora_arm(paired_df, traits, focal)
    assert len(res) == 2
    assert res["t1"].primary_stat > 0.0
    assert res["t2"].primary_stat > 0.0


def test_evaluate_loto_structure(
    synthetic_coupling_data: tuple[pd.DataFrame, pd.DataFrame, list[str]],
) -> None:
    """Verify LOTO matrix dimensions and diagonal handling."""
    profiles_df, paired_df, traits = synthetic_coupling_data

    loto_rows, matrices = evaluate_loto(
        profiles_df=profiles_df,
        paired_profiles_df=paired_df,
        trait_universe=traits,
        focal_traits=traits,
    )

    # 4 traits * 3 dropped targets each = 12 rows per arm * 3 arms = 36 rows
    assert len(loto_rows) == 36

    # Verify composition delta matrix
    eff_mat = matrices["composition_delta_effect"]
    assert eff_mat.shape == (4, 4)

    # Diagonal must be NaN (cannot drop trait and evaluate itself)
    for col in eff_mat.columns:
        assert np.isnan(eff_mat.loc[col, col])


def test_evaluate_lco_structure(
    synthetic_coupling_data: tuple[pd.DataFrame, pd.DataFrame, list[str]],
) -> None:
    """Verify LCO cluster dropping mechanics."""
    profiles_df, paired_df, traits = synthetic_coupling_data
    cluster = {"t1": "Trait 1", "t2": "Trait 2"}

    lco_rows = evaluate_lco(
        profiles_df=profiles_df,
        paired_profiles_df=paired_df,
        trait_universe=traits,
        focal_traits=traits,
        cluster_to_drop=cluster,
    )

    # 2 traits remain (t3, t4) * 3 arms = 6 rows
    assert len(lco_rows) == 6
    remaining_ids = {r.trait_id for r in lco_rows}
    assert remaining_ids == {"t3", "t4"}


def test_decision_rule_criterion() -> None:
    """Verify decision rule pass/fail thresholding."""
    from mtp.coupling import ArmResult

    # Mock baseline composition results
    base_comp = {
        "t1": ArmResult("t1", "t1", -0.10, 0.01, 0.01),
        "t2": ArmResult("t2", "t2", 0.05, 0.01, 0.01),
        "t3": ArmResult("t3", "t3", 0.12, 0.01, 0.01),
    }
    # Sorted effects: [-0.10, 0.05, 0.12]
    # Distances: [0.15, 0.07]
    # min_adj_dist: 0.07

    # Case A: Perturbation is small (0.01 < 0.07) -> passes
    df_small = pd.DataFrame(
        [[np.nan, 0.01, 0.0], [0.01, np.nan, 0.0], [0.0, 0.0, np.nan]],
        index=["t1", "t2", "t3"],
        columns=["t1", "t2", "t3"],
    )
    res_pass = evaluate_decision_rule(
        base_comp, df_small, target_trait_id="t2", dropped_trait_id="t1"
    )
    assert res_pass.passes_stability_criterion is True
    assert res_pass.min_adjacent_distance == pytest.approx(0.07)
    assert res_pass.abs_delta_effect == pytest.approx(0.01)

    # Case B: Perturbation is large (0.08 >= 0.07) -> fails
    df_large = pd.DataFrame(
        [[np.nan, 0.08, 0.0], [0.08, np.nan, 0.0], [0.0, 0.0, np.nan]],
        index=["t1", "t2", "t3"],
        columns=["t1", "t2", "t3"],
    )
    res_fail = evaluate_decision_rule(
        base_comp, df_large, target_trait_id="t2", dropped_trait_id="t1"
    )
    assert res_fail.passes_stability_criterion is False
    assert res_fail.abs_delta_effect == pytest.approx(0.08)


def test_full_coupling_audit_on_production_data(tmp_path: Path) -> None:
    """Execute complete audit against production profile data and verify exact outputs."""
    if not DEFAULT_PROFILES_PATH.exists() or not DEFAULT_COLUMNS_PATH.exists():
        pytest.skip("Production data files not found.")

    res = run_coupling_audit(
        profiles_path=DEFAULT_PROFILES_PATH,
        columns_path=DEFAULT_COLUMNS_PATH,
        output_dir=tmp_path,
    )

    # Verify output files exist and are populated
    assert res["loto_tsv"].exists()
    assert res["lco_tsv"].exists()
    assert res["report_txt"].exists()

    loto_df = pd.read_csv(res["loto_tsv"], sep="\t")
    lco_df = pd.read_csv(res["lco_tsv"], sep="\t")

    # 6 traits * 5 targets = 30 rows per arm * 3 arms = 90 rows
    assert len(loto_df) == 90
    # 4 remaining traits * 3 arms = 12 rows
    assert len(lco_df) == 12

    # Verify decision rule result matches empirical finding
    rule = res["rule_result"]
    assert rule.target_trait_id == COCCUS_ID
    assert rule.dropped_trait_id == SPHERE_ID
    assert pytest.approx(rule.min_adjacent_distance, abs=1e-5) == 0.002040
    assert pytest.approx(rule.delta_effect, abs=1e-5) == 0.000237
    assert rule.passes_stability_criterion is True


def test_cli_writes_where_the_committed_outputs_live() -> None:
    """`just audit-trait-coupling` must refresh the committed files (#477)."""
    import importlib.util

    from mtp.coupling import DEFAULT_OUTPUT_DIR
    from mtp.paths import DERIVED_DATA_DIR

    script = Path(__file__).resolve().parents[1] / "scripts" / "audit_trait_coupling.py"
    spec = importlib.util.spec_from_file_location("audit_trait_coupling", script)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert module.parse_args([]).output_dir == DEFAULT_OUTPUT_DIR == DERIVED_DATA_DIR
    for name in (
        "trait_coupling_loto_matrix.tsv",
        "trait_coupling_lco_summary.tsv",
        "trait_coupling_report.txt",
    ):
        assert (DEFAULT_OUTPUT_DIR / name).is_file(), name
