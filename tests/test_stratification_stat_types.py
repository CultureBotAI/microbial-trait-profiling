"""Keep q/p identities, scopes, and call bases explicit without changing values (#1507-#1511)."""

import gzip
import importlib.util
from dataclasses import replace
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from mtp.stratification import (
    ANNOTATION_DENSITY_NOTE,
    BEST_STAT_NOTE,
    PRECAUTIONARY_NOTE,
    PRIMARY_STRATA,
    StratumOutputs,
    build_sensitivity_summary_table,
    generate_sensitivity_reports,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_issue727_compare() -> ModuleType:
    script_path = REPO_ROOT / "scripts" / "adhoc" / "issue727_compare_core_strata.py"
    spec = importlib.util.spec_from_file_location("issue727_compare_core_strata", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_outputs(
    root: Path,
    stratum: str,
    q: float | None,
    p: float | None,
    stat_type: str | None = None,
    trait: str = "trait-a",
    fdr_scope: str | None = None,
    stat_threshold: float | str | None = None,
) -> StratumOutputs:
    tested = root / f"tested-{stratum}.tsv"
    with tested.open("w", encoding="utf-8") as stream:
        if stat_type is not None:
            stream.write(f"# Stat-Type: {stat_type}\n")
        if fdr_scope is not None:
            stream.write(f"# FDR-Scope: {fdr_scope}\n")
        if stat_threshold is not None:
            stream.write(f"# Stat-Threshold: {stat_threshold}\n")
        pd.DataFrame(
            [
                {
                    "trait": trait,
                    "trait_total": 3,
                    "n_profiles_tested": 2,
                    "prevalence": 0.5,
                    "mean_fold": 1.5,
                    "hit_rows": 1,
                    "test_status": "tested_significant",
                    "best_qvalue": q,
                    "best_pvalue": p,
                }
            ]
        ).to_csv(stream, sep="\t", index=False)
    ranked = root / f"ranked-{stratum}.tsv"
    # Deliberately contradictory ranking metadata must not override the complete table.
    ranked.write_text(
        "# Stat-Type: fdr_q\n# FDR-Scope: profile\n"
        "# Total-Testable-Profiles: 2\n# Ranked-Traits: 1\n",
        encoding="utf-8",
    )
    return StratumOutputs(
        stratum=stratum,
        profiles_tsv=root / "profiles.tsv",
        results_tsv=root / "results.tsv",
        report_txt=root / "report.txt",
        overall_ranked_tsv=ranked,
        pairs_ranked_tsv=root / "pairs.tsv",
        tested_traits_tsv=tested,
    )


@pytest.mark.parametrize(
    ("metadata", "scope", "q", "p", "expected", "kind", "expected_scope"),
    [
        (None, None, 0.4, 0.02, 0.4, "fdr_q", "unknown"),
        (None, None, 0.4, None, 0.4, "fdr_q", "unknown"),
        (None, None, None, 0.02, 0.02, "p_value", "none"),
        (None, None, None, None, None, "unavailable", "unavailable"),
        ("p_value", None, 0.02, 0.02, 0.02, "p_value", "none"),
        ("p_value", "none", None, 0.02, 0.02, "p_value", "none"),
        ("fdr_q", None, 0.4, 0.02, 0.4, "fdr_q", "unknown"),
        ("fdr_q", "global", 0.4, 0.02, 0.4, "fdr_q", "global"),
        ("fdr_q", "profile", 0.4, 0.02, 0.4, "fdr_q", "profile"),
        ("fdr_q", "global", None, 0.02, 0.02, "p_value", "none"),
    ],
)
def test_statistic_type_tracks_the_preserved_value(
    tmp_path: Path,
    metadata: str | None,
    scope: str | None,
    q: float | None,
    p: float | None,
    expected: float | None,
    kind: str,
    expected_scope: str,
) -> None:
    outputs = write_outputs(tmp_path, "pooled", q, p, metadata, fdr_scope=scope)
    frame, strata = build_sensitivity_summary_table({"pooled": outputs})
    row = frame.iloc[0]
    assert row["pooled_best_stat_type"] == kind
    assert row["pooled_best_stat_scope"] == expected_scope
    assert (
        row["pooled_best_stat_basis"]
        == {
            "fdr_q": "FDR q",
            "p_value": "nominal p",
            "unavailable": "unavailable",
        }[kind]
    )
    if expected is None:
        assert pd.isna(row["pooled_best_stat"])
    else:
        assert row["pooled_best_stat"] == expected
    assert row["pooled_n_tested"] == 2
    assert row["pooled_mean_fold"] == 1.5
    assert row["support_status"] == "unavailable_comparison"
    assert strata.iloc[0]["tested_traits"] == 1
    for stratum in ("feces", "feces_excluded"):
        assert pd.isna(row[f"{stratum}_best_stat"])
        assert row[f"{stratum}_best_stat_type"] == "unavailable"
        assert row[f"{stratum}_best_stat_scope"] == "unavailable"
        assert row[f"{stratum}_best_stat_basis"] == "unavailable"


def test_types_are_per_value_and_persist_in_reports(tmp_path: Path) -> None:
    outputs = {
        "pooled": write_outputs(
            tmp_path, "pooled", 0.4, 0.02, "fdr_q", fdr_scope="global", stat_threshold=0.001
        ),
        "feces": write_outputs(
            tmp_path, "feces", 0.02, 0.02, "p_value", fdr_scope="none", stat_threshold=0.05
        ),
        "feces_excluded": write_outputs(tmp_path, "feces_excluded", None, None),
    }
    tsv, txt = generate_sensitivity_reports(outputs, tmp_path / "reports")
    tsv_text = tsv.read_text(encoding="utf-8")
    txt_text = txt.read_text(encoding="utf-8")
    row = pd.read_csv(tsv, sep="\t", comment="#").iloc[0]
    assert [row[f"{s}_best_stat_type"] for s in outputs] == ["fdr_q", "p_value", "unavailable"]
    assert [row[f"{s}_best_stat_scope"] for s in outputs] == ["global", "none", "unavailable"]
    assert row["pooled_best_stat"] == 0.4
    assert row["feces_best_stat"] == 0.02
    assert row["feces_best_stat_basis"] == "nominal p"
    # Tested table metadata is authoritative and must not be overridden by ranking metadata (#1508).
    assert row["feces_significance_basis"] == "nominal p"
    assert row["feces_fdr_scope"] == "none"
    assert pd.isna(row["feces_excluded_best_stat"])
    assert row["support_status"] == "both_strata_significant"
    assert row["delta_mean_fold"] == 0
    assert row["delta_prevalence"] == 0
    # Pin every pre-existing output field, not only the three statistic values.
    expected: dict[str, str | float | bool] = {"trait": "trait-a", "trait_total_pooled": 3.0}
    for stratum, value in zip(outputs, [0.4, 0.02, float("nan")], strict=True):
        expected.update(
            {
                f"{stratum}_n_tested": 2,
                f"{stratum}_prevalence": 0.5,
                f"{stratum}_best_stat": value,
                f"{stratum}_mean_fold": 1.5,
                f"{stratum}_significant": True,
                f"{stratum}_test_status": "tested_significant",
            }
        )
    expected.update(
        support_status="both_strata_significant", delta_mean_fold=0.0, delta_prevalence=0.0
    )
    exported = pd.read_csv(tsv, sep="\t", comment="#")
    additive_suffixes = (
        "best_stat_type",
        "best_stat_scope",
        "best_stat_basis",
        "significance_basis",
        "fdr_scope",
    )
    legacy = exported.drop(
        columns=[f"{stratum}_{suffix}" for stratum in outputs for suffix in additive_suffixes]
    )
    pd.testing.assert_frame_equal(legacy, pd.DataFrame([expected]))
    for stratum, best_basis, sig_basis, best_scope, fdr_scope in zip(
        outputs,
        ["FDR q", "nominal p", "unavailable"],
        ["FDR q", "nominal p", "FDR q"],
        ["global", "none", "unavailable"],
        ["global", "none", "profile"],
        strict=True,
    ):
        assert row[f"{stratum}_best_stat_basis"] == best_basis
        assert row[f"{stratum}_significance_basis"] == sig_basis
        assert row[f"{stratum}_best_stat_scope"] == best_scope
        assert row[f"{stratum}_fdr_scope"] == fdr_scope
    assert set(exported.columns) == set(expected) | {
        f"{stratum}_{suffix}" for stratum in outputs for suffix in additive_suffixes
    }
    assert f"# Best-Statistic: {BEST_STAT_NOTE}\n" in tsv_text
    assert "# pooled-Threshold-Statistic: FDR q\n" in tsv_text
    assert "# pooled-Hit-Statistic: FDR q\n" in tsv_text
    assert "# pooled-FDR-Scope: global\n" in tsv_text
    assert "# pooled-Stat-Threshold: 0.001\n" in tsv_text
    assert "# feces-Threshold-Statistic: nominal p\n" in tsv_text
    assert "# feces-Hit-Statistic: nominal p\n" in tsv_text
    assert "# feces-FDR-Scope: none\n" in tsv_text
    assert "# feces-Stat-Threshold: 0.05\n" in tsv_text
    assert BEST_STAT_NOTE in txt_text
    assert "Hit statistic: FDR q | FDR scope: global" in txt_text
    assert "Hit statistic: nominal p | FDR scope: none" in txt_text


def test_absent_trait_has_no_statistic_type(tmp_path: Path) -> None:
    outputs = {
        "pooled": write_outputs(tmp_path, "pooled", 0.4, 0.02, "fdr_q", fdr_scope="global"),
        "feces": write_outputs(
            tmp_path, "feces", 0.4, 0.02, "fdr_q", trait="other", fdr_scope="profile"
        ),
    }
    frame, _ = build_sensitivity_summary_table(outputs)
    row = frame.set_index("trait").loc["trait-a"]
    assert row["feces_best_stat_type"] == "unavailable"
    assert row["feces_best_stat_scope"] == "unavailable"
    assert row["feces_test_status"] == "unavailable_trait"


def test_one_stratum_can_contain_q_p_and_unavailable_values(tmp_path: Path) -> None:
    outputs = write_outputs(tmp_path, "pooled", 0.4, 0.02, "fdr_q", fdr_scope="profile")
    path = outputs.tested_traits_tsv
    assert path is not None
    with path.open("a", encoding="utf-8") as stream:
        stream.write("fallback\t3\t2\t0.5\t1.5\t1\ttested_significant\t\t0.1\n")
        stream.write("missing\t3\t0\t\t\t0\tzero_eligible_tests\t\t\n")
    frame, _ = build_sensitivity_summary_table({"pooled": outputs})
    rows = frame.set_index("trait")
    assert rows["pooled_best_stat_type"].to_dict() == {
        "trait-a": "fdr_q",
        "fallback": "p_value",
        "missing": "unavailable",
    }
    assert rows["pooled_best_stat_scope"].to_dict() == {
        "trait-a": "profile",
        "fallback": "none",
        "missing": "unavailable",
    }
    assert rows["pooled_significance_basis"].to_dict() == {
        "trait-a": "FDR q",
        "fallback": "FDR q",
        "missing": "unavailable",
    }
    assert rows["pooled_fdr_scope"].to_dict() == {
        "trait-a": "profile",
        "fallback": "profile",
        "missing": "profile",
    }
    assert rows.loc["trait-a", "pooled_best_stat"] == 0.4
    assert rows.loc["fallback", "pooled_best_stat"] == 0.1
    assert pd.isna(rows.loc["missing", "pooled_best_stat"])


@pytest.mark.parametrize("scope", ["global", "profile"])
def test_fdr_reports_from_complete_tested_tables_emit_call_basis_and_scope_labels(
    tmp_path: Path, scope: str
) -> None:
    outputs = {
        name: write_outputs(
            tmp_path, name, 0.02, 0.004, "fdr_q", fdr_scope=scope, stat_threshold=0.001
        )
        for name in PRIMARY_STRATA
    }
    report_dir = tmp_path / "reports"
    tsv, txt = generate_sensitivity_reports(outputs, report_dir)
    tsv_text = tsv.read_text(encoding="utf-8")
    txt_text = txt.read_text(encoding="utf-8")
    assert f"# Best-Statistic: {BEST_STAT_NOTE}\n" in tsv_text
    for name in PRIMARY_STRATA:
        assert f"# {name}-Threshold-Statistic: FDR q\n" in tsv_text
        assert f"# {name}-Hit-Statistic: FDR q\n" in tsv_text
        assert f"# {name}-FDR-Scope: {scope}\n" in tsv_text
        assert f"# {name}-Stat-Threshold: 0.001\n" in tsv_text
    exported = pd.read_csv(tsv, sep="\t", comment="#")
    for name in PRIMARY_STRATA:
        assert exported[f"{name}_best_stat_type"].tolist() == ["fdr_q"]
        assert exported[f"{name}_best_stat_scope"].tolist() == [scope]
        assert exported[f"{name}_best_stat_basis"].tolist() == ["FDR q"]
        assert exported[f"{name}_significance_basis"].tolist() == ["FDR q"]
        assert exported[f"{name}_fdr_scope"].tolist() == [scope]
    deciles = pd.read_csv(report_dir / "annotation_density_deciles.tsv", sep="\t")
    assert deciles["threshold_statistic"].tolist() == ["FDR q"]
    assert deciles["hit_statistic"].tolist() == ["FDR q"]
    assert deciles["fdr_scope"].tolist() == [scope]
    assert f"Hit statistic: FDR q | FDR scope: {scope}" in txt_text
    assert "FDR hits: traits with any test passing the declared FDR q and fold thresholds." in (
        txt_text
    )
    assert (
        "Support-state names retain 'significant'; interpret each stratum using "
        "its hit statistic and FDR scope above." in txt_text
    )


@pytest.mark.parametrize(("q", "p"), [(0.4, 0.02), (0.4, None)])
def test_nominal_metadata_cannot_relabel_a_different_q(
    tmp_path: Path,
    q: float,
    p: float | None,
) -> None:
    outputs = write_outputs(tmp_path, "pooled", q, p, "p_value")
    with pytest.raises(ValueError, match=r"nominal.*best_qvalue"):
        build_sensitivity_summary_table({"pooled": outputs})


@pytest.mark.parametrize("equal_value", [0.02, 1.0])
def test_headerless_table_with_equal_q_and_p_is_refused(tmp_path: Path, equal_value: float) -> None:
    outputs = write_outputs(tmp_path, "pooled", equal_value, equal_value, None)
    with pytest.raises(ValueError, match=r"Stat-Type.*best_qvalue|best_qvalue.*best_pvalue"):
        build_sensitivity_summary_table({"pooled": outputs})


@pytest.mark.parametrize("metadata", ["", "q", "BY", "p_value\n# Stat-Type: fdr_q"])
def test_invalid_statistic_metadata_is_refused(tmp_path: Path, metadata: str) -> None:
    outputs = write_outputs(tmp_path, "pooled", 0.02, 0.02, metadata)
    with pytest.raises(ValueError, match="Stat-Type"):
        build_sensitivity_summary_table({"pooled": outputs})


@pytest.mark.parametrize(
    ("stat_type", "fdr_scope", "stat_threshold", "match"),
    [
        ("fdr_q", "", None, "FDR-Scope"),
        ("fdr_q", "by", None, "FDR-Scope"),
        ("fdr_q", "global\n# FDR-Scope: profile", None, "FDR-Scope"),
        ("p_value", "global", None, "FDR-Scope"),
        ("p_value", "profile", None, "FDR-Scope"),
        ("fdr_q", "none", None, "FDR-Scope"),
        ("fdr_q", "global", "", "Stat-Threshold"),
        ("fdr_q", "global", "0", "Stat-Threshold"),
        ("fdr_q", "global", "-0.01", "Stat-Threshold"),
        ("fdr_q", "global", "1.5", "Stat-Threshold"),
        ("fdr_q", "global", "nan", "Stat-Threshold"),
        ("fdr_q", "global", "0.05\n# Stat-Threshold: 0.01", "Stat-Threshold"),
    ],
)
def test_invalid_or_contradictory_scope_and_threshold_metadata_are_refused(
    tmp_path: Path,
    stat_type: str,
    fdr_scope: str | None,
    stat_threshold: str | None,
    match: str,
) -> None:
    outputs = write_outputs(
        tmp_path,
        "pooled",
        0.02,
        0.02 if stat_type == "p_value" else 0.004,
        stat_type,
        fdr_scope=fdr_scope,
        stat_threshold=stat_threshold,
    )
    with pytest.raises(ValueError, match=match):
        build_sensitivity_summary_table({"pooled": outputs})


@pytest.mark.parametrize("compressed", [False, True])
def test_blank_header_lines_and_gzip_keep_nominal_identity(
    tmp_path: Path, compressed: bool
) -> None:
    outputs = write_outputs(tmp_path, "pooled", 0.02, 0.02, "p_value", fdr_scope="none")
    path = outputs.tested_traits_tsv
    assert path is not None
    content = "\n# Introductory comment\n\n" + path.read_text(encoding="utf-8")
    if compressed:
        path = path.with_suffix(".tsv.gz")
        with gzip.open(path, "wt", encoding="utf-8") as stream:
            stream.write(content)
    else:
        path.write_text(content, encoding="utf-8")
    frame, _ = build_sensitivity_summary_table({"pooled": replace(outputs, tested_traits_tsv=path)})
    assert frame.iloc[0]["pooled_best_stat"] == 0.02
    assert frame.iloc[0]["pooled_best_stat_type"] == "p_value"
    assert frame.iloc[0]["pooled_best_stat_scope"] == "none"


def _assemble_24_strata_products(
    directory: Path,
    *,
    include_stat_type_col: bool,
    include_full_labels: bool,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    header = "study_id\texperiment_id\tinc_dec\n"
    feces = header + "s1\te1\tinc\ns1\te1\tdec\n"
    other = header + "s2\te2\tinc\ns2\te2\tdec\n"
    (directory / "study_trait_profiles_feces.tsv").write_text(feces, encoding="utf-8")
    (directory / "study_trait_profiles_feces_excluded.tsv").write_text(other, encoding="utf-8")
    (directory / "study_trait_profiles_pooled.tsv").write_text(
        header + feces.split("\n", 1)[1] + other.split("\n", 1)[1],
        encoding="utf-8",
    )
    outputs: dict[str, StratumOutputs] = {}
    for stratum in PRIMARY_STRATA:
        tested_path = directory / f"trait_enrichment_traits_tested_{stratum}.tsv"
        ranked_path = directory / f"trait_enrichment_traits_overall_ranked_{stratum}.tsv"
        pairs_path = directory / f"trait_enrichment_trait_study_pairs_ranked_{stratum}.tsv"
        results_path = directory / f"trait_enrichment_results_{stratum}.tsv"
        all_results_path = directory / f"trait_enrichment_all_results_{stratum}.tsv"
        report_path = directory / f"trait_enrichment_report_{stratum}.txt"
        tested_header = (
            "# Stat-Type: fdr_q\n# FDR-Scope: global\n# Stat-Threshold: 0.001\n"
            if include_full_labels
            else ""
        )
        tested_df = pd.DataFrame(
            [
                {
                    "trait": "trait-a",
                    "trait_total": 3,
                    "n_profiles_tested": 2,
                    "prevalence": 0.5,
                    "mean_fold": 1.5,
                    "hit_rows": 1,
                    "test_status": "tested_significant",
                    "best_qvalue": 0.02,
                    "best_pvalue": 0.004,
                }
            ]
        )
        tested_path.write_text(
            tested_header + tested_df.to_csv(sep="\t", index=False),
            encoding="utf-8",
        )
        ranked_path.write_text(
            "# Total-Testable-Profiles: 2\n# Ranked-Traits: 1\n"
            "# Stat-Type: fdr_q\n# FDR-Scope: global\ntrait\ntrait-a\n",
            encoding="utf-8",
        )
        pairs_path.write_text("rank\ttrait\n1\ttrait-a\n", encoding="utf-8")
        results_path.write_text("study_id\ttrait\ns1\ttrait-a\n", encoding="utf-8")
        all_results_path.write_text("study_id\ttrait\ns1\ttrait-a\n", encoding="utf-8")
        report_path.write_text("FDR q <= 0.001\n", encoding="utf-8")
        outputs[stratum] = StratumOutputs(
            stratum=stratum,
            profiles_tsv=directory / f"study_trait_profiles_{stratum}.tsv",
            results_tsv=results_path,
            report_txt=report_path,
            overall_ranked_tsv=ranked_path,
            pairs_ranked_tsv=pairs_path,
            tested_traits_tsv=tested_path,
            all_results_tsv=all_results_path,
        )
    tsv_path, txt_path = generate_sensitivity_reports(outputs, directory)
    # Restore identical tested_traits_tsv bytes across before/after so only 3 summary files differ.
    for stratum in PRIMARY_STRATA:
        tested_path = directory / f"trait_enrichment_traits_tested_{stratum}.tsv"
        body = pd.read_csv(tested_path, sep="\t", comment="#").to_csv(sep="\t", index=False)
        tested_path.write_text(
            "# Stat-Type: fdr_q\n# FDR-Scope: global\n# Stat-Threshold: 0.001\n" + body,
            encoding="utf-8",
        )
    if not include_stat_type_col and not include_full_labels:
        # Simulate pre-#1489 summary TSV and text report (no *_best_stat_type or BEST_STAT_NOTE).
        df = pd.read_csv(tsv_path, sep="\t", comment="#")
        df = df.drop(columns=[c for c in df.columns if c.endswith("_best_stat_type")])
        with tsv_path.open("w", encoding="utf-8", newline="") as f:
            f.write(f"# Note: {PRECAUTIONARY_NOTE}\n")
            f.write(f"# Diagnostic: {ANNOTATION_DENSITY_NOTE}\n")
            f.write("# Evaluated-Strata: pooled,feces,feces_excluded\n")
            df.to_csv(f, sep="\t", index=False)
        txt_path.write_text(
            txt_path.read_text(encoding="utf-8").replace(f"\n{BEST_STAT_NOTE}\n", ""),
            encoding="utf-8",
        )


@pytest.mark.parametrize("before_has_1489", [False, True], ids=["pre-1489", "post-1489"])
def test_core_strata_comparison_accepts_additive_stat_and_scope_labels(
    tmp_path: Path, before_has_1489: bool
) -> None:
    issue727 = _load_issue727_compare()
    before = tmp_path / "before"
    after = tmp_path / "after"
    _assemble_24_strata_products(
        before, include_stat_type_col=before_has_1489, include_full_labels=False
    )
    _assemble_24_strata_products(after, include_stat_type_col=True, include_full_labels=True)
    result = issue727.compare(before, after)
    assert result["passed"] is True
    products = {entry["file"]: entry for entry in result["products"]}
    for report_name in (
        "stratification_sensitivity_summary.tsv",
        "annotation_density_deciles.tsv",
        "stratification_sensitivity_summary.txt",
    ):
        assert products[report_name]["byte_identical"] is False
        assert products[report_name]["scientifically_identical"] is True


@pytest.mark.parametrize(
    ("target_file", "old_substr", "new_substr"),
    [
        ("stratification_sensitivity_summary.tsv", "\t0.02\t", "\t0.03\t"),
        (
            "stratification_sensitivity_summary.tsv",
            "both_strata_significant",
            "neither_significant",
        ),
        ("stratification_sensitivity_summary.tsv", "\tglobal\t", "\tinvalid_scope\t"),
        (
            "stratification_sensitivity_summary.tsv",
            "\tglobal\tFDR q\tFDR q\tglobal\t",
            "\tglobal\tFDR q\tFDR q\tinvalid_scope\t",
        ),
        (
            "stratification_sensitivity_summary.tsv",
            "# pooled-FDR-Scope: global",
            "# pooled-FDR-Scope: bad",
        ),
        ("annotation_density_deciles.tsv", "100.0", "50.0"),
        ("annotation_density_deciles.tsv", "FDR q\tFDR q\tglobal", "FDR q\tFDR q\tinvalid_scope"),
        ("annotation_density_deciles.tsv", "FDR q\tFDR q\tglobal", "FDR q\tnominal p\tglobal"),
        ("stratification_sensitivity_summary.txt", "100.0   %", "50.0    %"),
        ("stratification_sensitivity_summary.txt", "FDR scope: global", "FDR scope: invalid"),
    ],
)
def test_core_strata_comparison_rejects_scientific_and_invalid_label_changes(
    tmp_path: Path, target_file: str, old_substr: str, new_substr: str
) -> None:
    issue727 = _load_issue727_compare()
    before = tmp_path / "before"
    after = tmp_path / "after"
    _assemble_24_strata_products(before, include_stat_type_col=False, include_full_labels=False)
    _assemble_24_strata_products(after, include_stat_type_col=True, include_full_labels=True)
    target = after / target_file
    content = target.read_text(encoding="utf-8")
    assert old_substr in content
    target.write_text(content.replace(old_substr, new_substr, 1), encoding="utf-8")
    with pytest.raises(ValueError, match=f"Scientific output differs: {target_file}"):
        issue727.compare(before, after)
