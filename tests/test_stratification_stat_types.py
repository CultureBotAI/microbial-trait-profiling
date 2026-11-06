"""Keep q/p identities explicit without changing summary values (#1418)."""

import gzip
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from mtp.stratification import (
    StratumOutputs,
    build_sensitivity_summary_table,
    generate_sensitivity_reports,
)


def write_outputs(
    root: Path,
    stratum: str,
    q: float | None,
    p: float | None,
    stat_type: str | None = None,
    trait: str = "trait-a",
) -> StratumOutputs:
    tested = root / f"tested-{stratum}.tsv"
    with tested.open("w") as stream:
        if stat_type is not None:
            stream.write(f"# Stat-Type: {stat_type}\n")
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
    ranked.write_text("# Stat-Type: fdr_q\n# Total-Testable-Profiles: 2\n# Ranked-Traits: 1\n")
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
    ("metadata", "q", "p", "expected", "kind"),
    [
        (None, 0.4, 0.02, 0.4, "fdr_q"),
        (None, 0.4, None, 0.4, "fdr_q"),
        (None, None, 0.02, 0.02, "p_value"),
        (None, None, None, None, "unavailable"),
        ("p_value", 0.02, 0.02, 0.02, "p_value"),
        ("p_value", None, 0.02, 0.02, "p_value"),
        ("fdr_q", 0.4, 0.02, 0.4, "fdr_q"),
        ("fdr_q", None, 0.02, 0.02, "p_value"),
    ],
)
def test_statistic_type_tracks_the_preserved_value(
    tmp_path: Path,
    metadata: str | None,
    q: float | None,
    p: float | None,
    expected: float | None,
    kind: str,
) -> None:
    outputs = write_outputs(tmp_path, "pooled", q, p, metadata)
    frame, strata = build_sensitivity_summary_table({"pooled": outputs})
    row = frame.iloc[0]
    assert row["pooled_best_stat_type"] == kind
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
        assert row[f"{stratum}_best_stat_basis"] == "unavailable"


def test_types_are_per_value_and_persist_in_reports(tmp_path: Path) -> None:
    outputs = {
        "pooled": write_outputs(tmp_path, "pooled", 0.4, 0.02, "fdr_q"),
        "feces": write_outputs(tmp_path, "feces", 0.02, 0.02, "p_value"),
        "feces_excluded": write_outputs(tmp_path, "feces_excluded", None, None),
    }
    tsv, txt = generate_sensitivity_reports(outputs, tmp_path / "reports")
    row = pd.read_csv(tsv, sep="\t", comment="#").iloc[0]
    assert [row[f"{s}_best_stat_type"] for s in outputs] == ["fdr_q", "p_value", "unavailable"]
    assert row["pooled_best_stat"] == 0.4
    assert row["feces_best_stat"] == 0.02
    assert row["feces_best_stat_basis"] == "nominal p"
    # The hit rule is read independently from the ranking metadata.
    assert row["feces_significance_basis"] == "FDR q"
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
    legacy = exported.drop(
        columns=[
            f"{stratum}_{suffix}"
            for stratum in outputs
            for suffix in ("best_stat_type", "best_stat_basis", "significance_basis", "fdr_scope")
        ]
    )
    pd.testing.assert_frame_equal(legacy, pd.DataFrame([expected]))
    for stratum, basis in zip(outputs, ["FDR q", "nominal p", "unavailable"], strict=True):
        assert row[f"{stratum}_best_stat_basis"] == basis
        assert row[f"{stratum}_significance_basis"] == "FDR q"
        assert row[f"{stratum}_fdr_scope"] == "unknown"
    assert set(exported.columns) == set(expected) | {
        f"{stratum}_{suffix}"
        for stratum in outputs
        for suffix in ("best_stat_type", "best_stat_basis", "significance_basis", "fdr_scope")
    }
    assert "best_stat_type" in tsv.read_text()
    assert "best_stat_type" in txt.read_text()


def test_absent_trait_has_no_statistic_type(tmp_path: Path) -> None:
    outputs = {
        "pooled": write_outputs(tmp_path, "pooled", 0.4, 0.02),
        "feces": write_outputs(tmp_path, "feces", 0.4, 0.02, trait="other"),
    }
    frame, _ = build_sensitivity_summary_table(outputs)
    row = frame.set_index("trait").loc["trait-a"]
    assert row["feces_best_stat_type"] == "unavailable"
    assert row["feces_test_status"] == "unavailable_trait"


def test_one_stratum_can_contain_q_p_and_unavailable_values(tmp_path: Path) -> None:
    outputs = write_outputs(tmp_path, "pooled", 0.4, 0.02, "fdr_q")
    path = outputs.tested_traits_tsv
    assert path is not None
    with path.open("a") as stream:
        stream.write("fallback\t3\t2\t0.5\t1.5\t1\ttested_significant\t\t0.1\n")
        stream.write("missing\t3\t0\t\t\t0\tzero_eligible_tests\t\t\n")
    frame, _ = build_sensitivity_summary_table({"pooled": outputs})
    rows = frame.set_index("trait")
    assert rows["pooled_best_stat_type"].to_dict() == {
        "trait-a": "fdr_q",
        "fallback": "p_value",
        "missing": "unavailable",
    }
    assert rows.loc["trait-a", "pooled_best_stat"] == 0.4
    assert rows.loc["fallback", "pooled_best_stat"] == 0.1
    assert pd.isna(rows.loc["missing", "pooled_best_stat"])


@pytest.mark.parametrize(("q", "p"), [(0.4, 0.02), (0.4, None)])
def test_nominal_metadata_cannot_relabel_a_different_q(
    tmp_path: Path,
    q: float,
    p: float | None,
) -> None:
    outputs = write_outputs(tmp_path, "pooled", q, p, "p_value")
    with pytest.raises(ValueError, match=r"nominal.*best_qvalue"):
        build_sensitivity_summary_table({"pooled": outputs})


@pytest.mark.parametrize("metadata", ["", "q", "BY", "p_value\n# Stat-Type: fdr_q"])
def test_invalid_statistic_metadata_is_refused(tmp_path: Path, metadata: str) -> None:
    outputs = write_outputs(tmp_path, "pooled", 0.02, 0.02, metadata)
    with pytest.raises(ValueError, match="Stat-Type"):
        build_sensitivity_summary_table({"pooled": outputs})


@pytest.mark.parametrize("compressed", [False, True])
def test_blank_header_lines_and_gzip_keep_nominal_identity(
    tmp_path: Path, compressed: bool
) -> None:
    outputs = write_outputs(tmp_path, "pooled", 0.02, 0.02, "p_value")
    path = outputs.tested_traits_tsv
    assert path is not None
    content = "\n# Introductory comment\n\n" + path.read_text()
    if compressed:
        path = path.with_suffix(".tsv.gz")
        with gzip.open(path, "wt") as stream:
            stream.write(content)
    else:
        path.write_text(content)
    frame, _ = build_sensitivity_summary_table({"pooled": replace(outputs, tested_traits_tsv=path)})
    assert frame.iloc[0]["pooled_best_stat"] == 0.02
    assert frame.iloc[0]["pooled_best_stat_type"] == "p_value"
