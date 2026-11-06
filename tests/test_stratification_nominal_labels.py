"""Pin additive ORA statistic labels without changing calculations (#1503)."""

from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from mtp.stratification import (
    StratumOutputs,
    build_sensitivity_summary_table,
    generate_sensitivity_reports,
)


def _outputs(
    root: Path,
    stratum: str,
    stat_type: str | None,
    scope: str | None,
    *,
    hits: int = 1,
) -> StratumOutputs:
    directory = root / stratum
    directory.mkdir()
    tested = directory / "tested.tsv"
    ranked = directory / "ranked.tsv"
    frame = pd.DataFrame(
        {
            "trait": ["TraitA"],
            "trait_total": [20],
            "n_profiles_tested": [4],
            "prevalence": [0.25],
            "mean_fold": [2.5],
            "best_pvalue": [0.004],
            "best_qvalue": [0.02 if stat_type == "fdr_q" else float("nan")],
            "hit_rows": [hits],
            "test_status": ["tested_significant" if hits else "tested_nonsignificant"],
        }
    )
    frame.to_csv(tested, sep="\t", index=False)
    metadata = "# Total-Testable-Profiles: 4\n# Ranked-Traits: 1\n"
    if stat_type is not None:
        metadata += f"# Stat-Type: {stat_type}\n"
    if scope is not None:
        metadata += f"# FDR-Scope: {scope}\n"
    ranked.write_text(metadata + "trait\nTraitA\n", encoding="utf-8")
    return StratumOutputs(
        stratum=stratum,
        profiles_tsv=directory / "profiles.tsv",
        results_tsv=directory / "results.tsv",
        report_txt=directory / "report.txt",
        overall_ranked_tsv=ranked,
        pairs_ranked_tsv=directory / "pairs.tsv",
        tested_traits_tsv=tested,
    )


def test_mixed_strata_keep_values_and_declare_each_hit_and_statistic_basis(tmp_path: Path) -> None:
    outputs = {
        "pooled": _outputs(tmp_path, "pooled", "p_value", "none"),
        "feces": _outputs(tmp_path, "feces", "fdr_q", "profile"),
        "feces_excluded": _outputs(tmp_path, "feces_excluded", "fdr_q", "global", hits=0),
    }
    traits, strata = build_sensitivity_summary_table(outputs)
    row = traits.iloc[0]
    expected_legacy: dict[str, object] = {
        "trait": "TraitA",
        "trait_total_pooled": 20.0,
        "support_status": "feces_only_significant",
        "delta_mean_fold": 0.0,
        "delta_prevalence": 0.0,
    }
    for name, stat, significant in [
        ("pooled", 0.004, True),
        ("feces", 0.02, True),
        ("feces_excluded", 0.02, False),
    ]:
        expected_legacy.update(
            {
                f"{name}_n_tested": 4,
                f"{name}_prevalence": 0.25,
                f"{name}_best_stat": stat,
                f"{name}_mean_fold": 2.5,
                f"{name}_significant": significant,
                f"{name}_test_status": (
                    "tested_significant" if significant else "tested_nonsignificant"
                ),
            }
        )
    assert row[list(expected_legacy)].to_dict() == expected_legacy
    for name, basis, scope in [
        ("pooled", "nominal p", "none"),
        ("feces", "FDR q", "profile"),
        ("feces_excluded", "FDR q", "global"),
    ]:
        assert row[f"{name}_best_stat_basis"] == basis
        assert row[f"{name}_significance_basis"] == basis
        assert row[f"{name}_fdr_scope"] == scope
        meta = strata.set_index("stratum").loc[name]
        assert meta["threshold_statistic"] == basis
        assert meta["hit_statistic"] == basis
        assert meta["fdr_scope"] == scope
        assert meta["stat_type"] == ("p_value" if name == "pooled" else "fdr_q")


def test_missing_metadata_and_unavailable_analysis_are_not_called_corrected(tmp_path: Path) -> None:
    outputs = {"pooled": _outputs(tmp_path, "pooled", None, None, hits=0)}
    traits, strata = build_sensitivity_summary_table(outputs)
    row = traits.iloc[0]
    assert row["pooled_best_stat"] == 0.004
    assert row["pooled_best_stat_basis"] == "nominal p"
    assert row["pooled_significance_basis"] == "unknown"
    assert row["pooled_fdr_scope"] == "unknown"
    assert not row["pooled_significant"]
    assert row["support_status"] == "unavailable_comparison"
    assert row["feces_best_stat_basis"] == "unavailable"
    assert row["feces_significance_basis"] == "unavailable"
    assert row["feces_fdr_scope"] == "unavailable"
    assert row["feces_significant"] is None
    assert strata.iloc[0]["stat_type"] == "unknown"
    assert strata.iloc[0]["threshold_statistic"] == "unknown"
    assert strata.iloc[0]["hit_statistic"] == "unknown"


@pytest.mark.parametrize(
    ("stat_type", "scope", "basis", "hit_label"),
    [
        ("p_value", "none", "nominal p", "Nominal hits"),
        ("fdr_q", "profile", "FDR q", "FDR hits"),
        (None, None, "unknown", "Hits (statistic not declared)"),
    ],
)
def test_written_reports_and_density_counts_label_the_pooled_hit_basis(
    tmp_path: Path,
    stat_type: str | None,
    scope: str | None,
    basis: str,
    hit_label: str,
) -> None:
    outputs = {"pooled": _outputs(tmp_path, "pooled", stat_type, scope)}
    report_dir = tmp_path / "reports"
    tsv, txt = generate_sensitivity_reports(outputs, report_dir)
    summary = tsv.read_text(encoding="utf-8")
    report = txt.read_text(encoding="utf-8")
    assert f"# pooled-Threshold-Statistic: {basis}\n" in summary
    assert f"# pooled-Hit-Statistic: {basis}\n" in summary
    assert f"# pooled-FDR-Scope: {scope or 'unknown'}\n" in summary
    assert f"Hit statistic: {basis} | FDR scope: {scope or 'unknown'}" in report
    if stat_type is None:
        assert (
            "Hits (statistic not declared): traits with any test passing the analyzer's thresholds."
            in report
        )
    else:
        assert (
            f"{hit_label}: traits with any test passing the declared {basis} and fold thresholds."
            in report
        )
    assert "Sig Traits:" not in report
    assert "declared p/q" not in report
    assert (
        "Support-state names retain 'significant'; interpret each stratum using "
        "its hit statistic and FDR scope above." in report
    )
    assert f"{hit_label:<24} {'% Hits':<8}" in report
    density = pd.read_csv(report_dir / "annotation_density_deciles.tsv", sep="\t")
    assert density["threshold_statistic"].tolist() == [basis]
    assert density["hit_statistic"].tolist() == [basis]
    assert density["fdr_scope"].tolist() == [scope or "unknown"]
    assert density.iloc[0][
        [
            "decile",
            "min_trait_total",
            "max_trait_total",
            "total_traits",
            "significant_traits",
            "pct_significant",
        ]
    ].to_dict() == {
        "decile": 1,
        "min_trait_total": 20,
        "max_trait_total": 20,
        "total_traits": 1,
        "significant_traits": 1,
        "pct_significant": 100.0,
    }


@pytest.mark.parametrize(("stat_type", "scope"), [("p_value", "none"), ("fdr_q", "global")])
def test_untested_trait_does_not_inherit_declared_stratum_hit_basis(
    tmp_path: Path, stat_type: str, scope: str
) -> None:
    output = _outputs(tmp_path, "pooled", stat_type, scope)
    assert output.tested_traits_tsv is not None
    tested = pd.read_csv(output.tested_traits_tsv, sep="\t")
    untested = pd.DataFrame(
        {
            "trait": ["UntestedTrait"],
            "trait_total": [0],
            "n_profiles_tested": [0],
            "prevalence": [float("nan")],
            "mean_fold": [float("nan")],
            "best_pvalue": [float("nan")],
            "best_qvalue": [float("nan")],
            "hit_rows": [0],
            "test_status": ["untested_zero_trait_total"],
        }
    )
    pd.concat([tested, untested], ignore_index=True).to_csv(
        output.tested_traits_tsv, sep="\t", index=False
    )
    traits, strata = build_sensitivity_summary_table({"pooled": output})
    row = traits.set_index("trait").loc["UntestedTrait"]
    assert row["pooled_n_tested"] == 0
    assert row["pooled_significant"] is None
    assert row["pooled_test_status"] == "untested_zero_trait_total"
    assert row["pooled_best_stat_basis"] == "unavailable"
    assert row["pooled_significance_basis"] == "unavailable"
    assert row["pooled_fdr_scope"] == scope
    assert strata.iloc[0]["hit_statistic"] == ("nominal p" if stat_type == "p_value" else "FDR q")
    assert strata.iloc[0]["tested_traits"] == 1
    assert strata.iloc[0]["untested_traits"] == 1


@pytest.mark.parametrize("missing_path", [False, True], ids=["none", "absent-file"])
def test_missing_tested_table_resets_parsed_labels_and_explains_unavailable_hits(
    tmp_path: Path, missing_path: bool
) -> None:
    pooled = _outputs(tmp_path, "pooled", "fdr_q", "profile")
    pooled = replace(
        pooled,
        tested_traits_tsv=tmp_path / "absent-tested.tsv" if missing_path else None,
    )
    outputs = {
        "pooled": pooled,
        "feces": _outputs(tmp_path, "feces", "p_value", "none"),
    }
    traits, strata = build_sensitivity_summary_table(outputs)
    pooled_metadata = strata.set_index("stratum").loc["pooled"]
    for field in ("stat_type", "threshold_statistic", "hit_statistic", "fdr_scope"):
        assert pooled_metadata[field] == "unavailable"
    row = traits.set_index("trait").loc["TraitA"]
    assert row["pooled_test_status"] == "unavailable_analysis"
    assert row["pooled_significant"] is None
    assert pd.isna(row["pooled_best_stat"])
    assert row["pooled_best_stat_basis"] == "unavailable"
    assert row["pooled_significance_basis"] == "unavailable"
    assert row["pooled_fdr_scope"] == "unavailable"

    report_dir = tmp_path / "reports"
    tsv, txt = generate_sensitivity_reports(outputs, report_dir)
    summary = tsv.read_text(encoding="utf-8")
    report = txt.read_text(encoding="utf-8")
    assert "# pooled-Threshold-Statistic: unavailable\n" in summary
    assert "# pooled-Hit-Statistic: unavailable\n" in summary
    assert "# pooled-FDR-Scope: unavailable\n" in summary
    assert "Hit statistic: unavailable | FDR scope: unavailable" in report
    assert (
        "Hits (statistic not declared): traits with any test passing the analyzer's thresholds."
        in report
    )
    assert "declared unavailable" not in report
    assert pd.read_csv(report_dir / "annotation_density_deciles.tsv", sep="\t").empty
