"""Tests for trait coupling sensitivity, LOTO, LCO, and decision rule (Issue #260)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mtp.coupling import (
    COCCUS_ID,
    DEFAULT_OUTPUT_DIR,
    SPHERE_ID,
    ArmResult,
    DecisionRuleResult,
    LCOSummaryRow,
    build_coupling_report_text,
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
    TABLE_1_TRAITS,
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


def test_compute_ora_arm_unbiased_with_zeros() -> None:
    """Verify that ORA mean fold accounts for zero-count profiles without selection bias (#480)."""
    # 10 profiles of 100 total counts each (M = 1000)
    # t_focal has count 10 in profile 0, and 0 in profiles 1..9 (K = 10)
    # Expected in each profile: (100 * 10) / 1000 = 1.0
    # Profile 0 fold = 10 / 1.0 = 10.0
    # Profiles 1..9 fold = 0 / 1.0 = 0.0
    # Expected unbiased mean fold = (10.0 + 9 * 0.0) / 10 = 1.0
    rows = []
    for i in range(10):
        rows.append(
            {
                "study_id": f"s_{i}",
                "experiment_id": f"e_{i}",
                "inc_dec": "inc",
                "t_focal": 10 if i == 0 else 0,
                "t_bg": 90 if i == 0 else 100,
                "t_zero": 0,
            }
        )
    df = pd.DataFrame(rows)
    traits = ["t_focal", "t_bg", "t_zero"]

    res = compute_ora_arm(df, traits, ["t_focal", "t_zero"])
    # Mean fold must be 1.0 across all profiles, not 10.0 (which conditioned on k > 0)
    assert res["t_focal"].primary_stat == pytest.approx(1.0)
    # For a trait with 0 total counts (K = 0), mean fold is 0.0
    assert res["t_zero"].primary_stat == 0.0
    assert res["t_zero"].p_value == 1.0
    assert res["t_zero"].extra_stat == 0.0


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

    from mtp.paths import DERIVED_DATA_DIR

    script = Path(__file__).resolve().parents[1] / "scripts" / "audit_trait_coupling.py"
    spec = importlib.util.spec_from_file_location("audit_trait_coupling", script)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    args_empty = module.parse_args([])
    assert args_empty.output_dir == DEFAULT_OUTPUT_DIR == DERIVED_DATA_DIR
    assert args_empty.code_commit == ""

    args_commit = module.parse_args(["--code-commit", "deadbeef+dirty"])
    assert args_commit.code_commit == "deadbeef+dirty"

    for name in (
        "trait_coupling_loto_matrix.tsv",
        "trait_coupling_lco_summary.tsv",
        "trait_coupling_report.txt",
    ):
        assert (DEFAULT_OUTPUT_DIR / name).is_file(), name


def test_build_coupling_report_text_records_code_commit() -> None:
    """Report text must include Code commit: <hash> when provided (#477)."""
    dummy_df = pd.DataFrame()
    dummy_rule = DecisionRuleResult(
        target_trait_id="t1",
        target_trait_name="Trait 1",
        dropped_trait_id="t2",
        dropped_trait_name="Trait 2",
        delta_effect=0.001,
        abs_delta_effect=0.001,
        min_adjacent_distance=0.005,
        passes_stability_criterion=True,
        adjacent_distances=[0.005],
        sorted_effects=[0.01, 0.015],
        verdict="Passes",
    )
    text_without = build_coupling_report_text(dummy_df, dummy_df, [], dummy_rule, dummy_rule)
    assert "Code commit:" not in text_without

    text_with = build_coupling_report_text(
        dummy_df, dummy_df, [], dummy_rule, dummy_rule, code_commit="abc123commit"
    )
    assert "Code commit: abc123commit" in text_with


def test_build_coupling_report_text_content() -> None:
    """Verify build_coupling_report_text section contents and formatting (#480)."""
    t1_names = ["coccus shaped", "sphere shaped"]
    delta_eff = pd.DataFrame(
        [[np.nan, 0.0001], [0.0002, np.nan]],
        index=t1_names,
        columns=t1_names,
    )
    delta_q = pd.DataFrame(
        [[np.nan, -0.01], [0.02, np.nan]],
        index=t1_names,
        columns=t1_names,
    )
    lco_rows = [
        LCOSummaryRow(
            arm="composition",
            cluster_dropped="coccus + sphere",
            trait_id="METPO:1000608",
            trait_name="facultatively aerobic",
            base_stat=0.05,
            lco_stat=0.051,
            delta_stat=0.001,
            base_q=0.01,
            lco_q=0.01,
            delta_q=0.0,
            base_extra=2.5,
            lco_extra=2.55,
        )
    ]
    rule_res = DecisionRuleResult(
        target_trait_id=COCCUS_ID,
        target_trait_name="coccus shaped",
        dropped_trait_id=SPHERE_ID,
        dropped_trait_name="sphere shaped",
        delta_effect=0.0001,
        abs_delta_effect=0.0001,
        min_adjacent_distance=0.002,
        passes_stability_criterion=True,
        adjacent_distances=[0.002],
        sorted_effects=[0.01, 0.012],
        verdict=(
            "Passes stability criterion: |Delta effect| (0.000100) < "
            "min adjacent row distance (0.002000)."
        ),
    )
    rule_rev = DecisionRuleResult(
        target_trait_id=SPHERE_ID,
        target_trait_name="sphere shaped",
        dropped_trait_id=COCCUS_ID,
        dropped_trait_name="coccus shaped",
        delta_effect=0.0002,
        abs_delta_effect=0.0002,
        min_adjacent_distance=0.002,
        passes_stability_criterion=True,
        adjacent_distances=[0.002],
        sorted_effects=[0.01, 0.012],
        verdict=(
            "Passes stability criterion: |Delta effect| (0.000200) < "
            "min adjacent row distance (0.002000)."
        ),
    )

    text = build_coupling_report_text(
        delta_eff_mat=delta_eff,
        delta_q_mat=delta_q,
        lco_rows=lco_rows,
        rule_result=rule_res,
        rule_result_reverse=rule_rev,
        code_commit="commit_abc_test",
    )

    assert "TRAIT COUPLING SENSITIVITY AUDIT REPORT (Issue #260)" in text
    assert "Code commit: commit_abc_test" in text
    assert "1. LEAVE-ONE-TRAIT-OUT (LOTO) - COMPOSITION DELTA EFFECT MATRIX (6 x 6):" in text
    assert "2. LEAVE-ONE-TRAIT-OUT (LOTO) - COMPOSITION DELTA Q-VALUE MATRIX (6 x 6):" in text
    assert "3. LEAVE-CLUSTER-OUT (LCO) - DROPPING COCCUS + SPHERE SIMULTANEOUSLY:" in text
    assert "4. FORMAL DECISION RULE EVALUATION:" in text
    assert "5. METHODOLOGICAL & BIOLOGICAL SYNTHESIS:" in text
    assert "facultatively aerobic" in text
    assert "Direct Test (Drop Sphere -> Coccus Effect Shift):" in text
    assert "Reverse Test (Drop Coccus -> Sphere Effect Shift):" in text
    assert "Passes stability criterion" in text

    for line in text.splitlines():
        assert line == line.rstrip(), f"Line contains trailing whitespace: {line!r}"


@pytest.fixture
def synthetic_audit_files(tmp_path: Path) -> tuple[Path, Path]:
    """Create minimal synthetic wide profiles and columns TSVs for audit pipeline testing."""
    t1_ids = list(TABLE_1_TRAITS.keys())
    all_trait_ids = [*t1_ids, "METPO:bg_extra_trait"]

    columns_df = pd.DataFrame(
        [{"column_id": cid, "object_category": "METPO:trait"} for cid in all_trait_ids]
    )
    columns_path = tmp_path / "taxa_trait_matrix_columns.tsv"
    columns_df.to_csv(columns_path, sep="\t", index=False)

    rows: list[dict[str, object]] = []
    # 4 experiments across 2 studies with paired inc and dec profiles
    for s_idx in range(2):
        study_id = f"study_{s_idx}"
        for e_idx in range(2):
            exp_id = f"exp_{s_idx}_{e_idx}"
            row_inc: dict[str, object] = {
                "study_id": study_id,
                "experiment_id": exp_id,
                "inc_dec": "inc",
            }
            row_dec: dict[str, object] = {
                "study_id": study_id,
                "experiment_id": exp_id,
                "inc_dec": "dec",
            }
            for t_idx, tid in enumerate(all_trait_ids):
                base = 20.0 + t_idx * 5.0 + s_idx * 2.0 + e_idx * 3.0
                row_inc[tid] = base + (3.0 + e_idx if t_idx % 2 == 0 else -(2.0 + e_idx))
                row_dec[tid] = base + (1.0 - s_idx if t_idx % 2 == 0 else (2.0 - s_idx))
            rows.append(row_inc)
            rows.append(row_dec)

    profiles_df = pd.DataFrame(rows)
    profiles_path = tmp_path / "study_trait_profiles_wide.tsv"
    profiles_df.to_csv(profiles_path, sep="\t", index=False)

    return profiles_path, columns_path


def test_full_coupling_audit_on_synthetic_data(
    synthetic_audit_files: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    """Execute complete audit pipeline on synthetic data without skipping (#480)."""
    profiles_path, columns_path = synthetic_audit_files
    out_dir = tmp_path / "synthetic_audit_out"

    res = run_coupling_audit(
        profiles_path=profiles_path,
        columns_path=columns_path,
        output_dir=out_dir,
        code_commit="test_commit_sha_123",
    )

    # Output files must exist and be written to output_dir
    assert res["loto_tsv"] == out_dir / "trait_coupling_loto_matrix.tsv"
    assert res["lco_tsv"] == out_dir / "trait_coupling_lco_summary.tsv"
    assert res["report_txt"] == out_dir / "trait_coupling_report.txt"
    assert res["loto_tsv"].is_file()
    assert res["lco_tsv"].is_file()
    assert res["report_txt"].is_file()

    # 6 Table 1 traits * 5 dropped traits = 30 pairs * 3 arms = 90 rows
    loto_df = pd.read_csv(res["loto_tsv"], sep="\t")
    assert len(loto_df) == 90
    assert set(loto_df["arm"]) == {"composition", "directional", "ora"}

    # 4 remaining traits * 3 arms = 12 rows
    lco_df = pd.read_csv(res["lco_tsv"], sep="\t")
    assert len(lco_df) == 12

    # Decision rule results must be evaluated
    rule = res["rule_result"]
    rule_rev = res["rule_result_rev"]
    assert isinstance(rule, DecisionRuleResult)
    assert isinstance(rule_rev, DecisionRuleResult)
    assert rule.target_trait_id == COCCUS_ID
    assert rule.dropped_trait_id == SPHERE_ID
    assert rule_rev.target_trait_id == SPHERE_ID
    assert rule_rev.dropped_trait_id == COCCUS_ID

    # Report text must contain provenance and all section headers
    report = res["report_txt"].read_text(encoding="utf-8")
    assert "Code commit: test_commit_sha_123" in report
    assert "TRAIT COUPLING SENSITIVITY AUDIT REPORT (Issue #260)" in report
    assert "1. LEAVE-ONE-TRAIT-OUT (LOTO) - COMPOSITION DELTA EFFECT MATRIX" in report
    assert "2. LEAVE-ONE-TRAIT-OUT (LOTO) - COMPOSITION DELTA Q-VALUE MATRIX" in report
    assert "3. LEAVE-CLUSTER-OUT (LCO) - DROPPING COCCUS + SPHERE SIMULTANEOUSLY:" in report
    assert "4. FORMAL DECISION RULE EVALUATION:" in report
    assert "5. METHODOLOGICAL & BIOLOGICAL SYNTHESIS:" in report


def test_cli_main_invokes_audit_and_writes_output(
    synthetic_audit_files: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    """Verify CLI main() runs end-to-end, accepts arguments, and writes to output_dir (#480)."""
    import importlib.util

    script = Path(__file__).resolve().parents[1] / "scripts" / "audit_trait_coupling.py"
    spec = importlib.util.spec_from_file_location("audit_trait_coupling", script)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    profiles_path, columns_path = synthetic_audit_files
    cli_out = tmp_path / "cli_audit_out"

    exit_code = module.main(
        [
            f"--profiles={profiles_path}",
            f"--columns={columns_path}",
            f"--output-dir={cli_out}",
            "--code-commit=deadbeef_test",
        ]
    )

    assert exit_code == 0
    assert (cli_out / "trait_coupling_loto_matrix.tsv").is_file()
    assert (cli_out / "trait_coupling_lco_summary.tsv").is_file()
    assert (cli_out / "trait_coupling_report.txt").is_file()

    report_content = (cli_out / "trait_coupling_report.txt").read_text(encoding="utf-8")
    assert "Code commit: deadbeef_test" in report_content


def test_cli_main_handles_error(tmp_path: Path) -> None:
    """Verify CLI main() handles non-existent paths gracefully (#480)."""
    import importlib.util

    script = Path(__file__).resolve().parents[1] / "scripts" / "audit_trait_coupling.py"
    spec = importlib.util.spec_from_file_location("audit_trait_coupling", script)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    bad_profiles = tmp_path / "does_not_exist.tsv"
    exit_code = module.main([f"--profiles={bad_profiles}"])
    assert exit_code == 1
