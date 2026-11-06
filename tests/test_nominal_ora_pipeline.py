"""Nominal ORA must never manufacture adjusted probabilities (#1409)."""

from __future__ import annotations

import csv
import importlib.util
import io
import sys
from pathlib import Path
from types import ModuleType

import pytest
from scipy.stats import false_discovery_control, hypergeom

from mtp.deck_visualizations import generate_figure2_calibrated_ora_infographic
from mtp.paths import REPO_ROOT
from mtp.quoted_results import OraCalls, refuse_nominal_ora
from mtp.stratification import run_stratum_enrichment
from mtp.workflow.downstream_contracts import get_standard_analysis_spec


def _script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / "scripts" / f"{name}.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def producers() -> tuple[ModuleType, ModuleType]:
    return _script("trait_enrichment"), _script("export_trait_enrichment_rankings")


@pytest.fixture
def inputs(tmp_path: Path) -> tuple[Path, Path]:
    profiles, matrix = tmp_path / "profiles.tsv", tmp_path / "matrix.tsv"
    profiles.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\tTraitB\n"
        "s1\te1\tinc\t2\t2\t1\n"
        "s1\te1\tdec\t2\t0\t2\n",
        encoding="utf-8",
    )
    matrix.write_text(
        "taxon_id\tTraitA\tTraitB\n"
        + "".join(f"NCBITaxon:{i}\t{int(i <= 3)}\t{int(i <= 5)}\n" for i in range(1, 11)),
        encoding="utf-8",
    )
    return profiles, matrix


def _table(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    metadata = dict(line[2:].split(": ", 1) for line in lines if line.startswith("# "))
    rows = list(
        csv.DictReader(
            io.StringIO("\n".join(line for line in lines if not line.startswith("#"))),
            delimiter="\t",
        )
    )
    return metadata, rows


@pytest.mark.parametrize("tail", ["upper", "lower"])
@pytest.mark.parametrize("mode", ["none", "global", "profile"])
def test_producer_declares_actual_correction_and_preserves_probabilities(
    tmp_path: Path,
    inputs: tuple[Path, Path],
    producers: tuple[ModuleType, ModuleType],
    tail: str,
    mode: str,
) -> None:
    analyzer, _ = producers
    profiles, matrix = inputs
    report, selected, complete = (
        tmp_path / name for name in ("report.txt", "selected.tsv", "complete.tsv")
    )
    analyzer.main(
        [
            "--input",
            str(profiles),
            "--taxa-matrix",
            str(matrix),
            "--output",
            str(report),
            "--tsv-output",
            str(selected),
            "--all-results-output",
            str(complete),
            "--tail",
            tail,
            "--p-threshold",
            "1",
            "--min-fold",
            "0",
            "--max-fold",
            "100",
            "--fdr-scope",
            mode if mode != "none" else "global",
            "--no-fdr" if mode == "none" else "--fdr",
        ]
    )
    metadata, rows = _table(complete)
    assert len(rows) == 4
    expected_p = {}
    for direction, counts in (("inc", (2, 1)), ("dec", (0, 2))):
        for trait, k, background in zip(("TraitA", "TraitB"), counts, (3, 5), strict=True):
            expected_p[(direction, trait)] = float(
                hypergeom.sf(k - 1, 10, background, 2)
                if tail == "upper"
                else hypergeom.cdf(k, 10, background, 2)
            )
    assert [float(row["p_value"]) for row in rows] == pytest.approx(
        [expected_p[(row["inc_dec"], row["trait"])] for row in rows]
    )
    assert metadata["FDR-Scope"] == mode
    assert {row["correction_family"] for row in rows} == {mode}
    selected_meta, selected_rows = _table(selected)
    assert selected_meta["FDR-Scope"] == mode
    assert len(selected_rows) == 4
    if mode == "none":
        assert all(row["fdr_q"] == "" for row in rows + selected_rows)
        assert [float(row["p_value"]) for row in rows] == pytest.approx(sorted(expected_p.values()))
        assert report.read_text().startswith("nominal p <= 1")
        assert "FDR" not in report.read_text()
        with pytest.raises(ValueError, match="not FDR-controlled"):
            OraCalls.from_header(complete.read_text().splitlines())
    else:
        for family in (
            [rows]
            if mode == "global"
            else [
                [row for row in rows if row["inc_dec"] == direction] for direction in ("inc", "dec")
            ]
        ):
            assert [float(row["fdr_q"]) for row in family] == pytest.approx(
                false_discovery_control([float(row["p_value"]) for row in family], method="bh")
            )
        assert report.read_text().startswith("FDR q <= 1")


def test_ora_spec_matches_independent_enrichment_and_depletion_families(
    tmp_path: Path, inputs: tuple[Path, Path], producers: tuple[ModuleType, ModuleType]
) -> None:
    """The arm's two real entry points each correct their own eligible tests (#1499)."""
    profiles, matrix = inputs
    enrichment, _ = producers
    depletion = _script("trait_depletion")
    probabilities: list[float] = []
    adjusted: list[float] = []
    for producer, tail in ((enrichment, "upper"), (depletion, "lower")):
        complete = tmp_path / f"{tail}_all_results.tsv"
        producer.main(
            [
                "--input",
                str(profiles),
                "--taxa-matrix",
                str(matrix),
                "--output",
                str(tmp_path / f"{tail}_report.txt"),
                "--tsv-output",
                str(tmp_path / f"{tail}_selected.tsv"),
                "--all-results-output",
                str(complete),
                "--fdr-scope",
                "global",
            ]
        )
        metadata, rows = _table(complete)
        assert len(rows) == 4
        assert metadata["FDR-Scope"] == "global"
        assert {row["correction_family"] for row in rows} == {"global"}
        p_values = [float(row["p_value"]) for row in rows]
        q_values = [float(row["fdr_q"]) for row in rows]
        # No --tail is supplied: exercise enrichment's default and the real
        # depletion forwarder, independently of the plan's requested labels.
        for row, observed in zip(rows, p_values, strict=True):
            count = int(row["count"])
            background = 3 if row["trait"] == "TraitA" else 5
            expected = (
                hypergeom.sf(count - 1, 10, background, 2)
                if tail == "upper"
                else hypergeom.cdf(count, 10, background, 2)
            )
            assert observed == pytest.approx(expected)
        assert q_values == pytest.approx(false_discovery_control(p_values, method="bh"))
        probabilities.extend(p_values)
        adjusted.extend(q_values)
    # This fixture distinguishes two tail families from one pooled family.
    assert adjusted != pytest.approx(false_discovery_control(probabilities, method="bh"))
    spec = get_standard_analysis_spec("ilyas_signature_hypergeometric_ora_v1")
    assert (
        spec.multiple_testing_policy == "global_BH_q_per_tail_across_eligible_profile_trait_tests"
    )


def test_nominal_stratum_rankings_and_report_fallback_keep_q_missing(
    tmp_path: Path, inputs: tuple[Path, Path], producers: tuple[ModuleType, ModuleType]
) -> None:
    profiles, matrix = inputs
    outputs = run_stratum_enrichment(
        "fixture",
        profiles,
        tmp_path / "stratum",
        p_threshold=1.0,
        min_fold=0.0,
        fdr=False,
        taxa_matrix=matrix,
    )
    for path in (outputs.overall_ranked_tsv, outputs.pairs_ranked_tsv, outputs.tested_traits_tsv):
        assert path is not None
        metadata, rows = _table(path)
        assert rows
        assert metadata["FDR-Scope"] == "none"
        assert metadata["Stat-Type"] == "p_value"
        for row in rows:
            assert all(
                value == ""
                for key, value in row.items()
                if key in {"fdr_q", "best_qvalue", "median_qvalue", "q_value"}
            )
    _, rankings = producers
    overall, pairs = tmp_path / "overall.tsv", tmp_path / "pairs.tsv"
    assert (
        rankings.main(
            [
                "--input",
                str(outputs.report_txt),
                "--results",
                str(tmp_path / "absent.tsv"),
                "--overall-output",
                str(overall),
                "--pair-output",
                str(pairs),
                "--fdr-scope",
                "global",
            ]
        )
        == 0
    )
    metadata, rows = _table(pairs)
    assert metadata["FDR-Scope"] == "none"
    assert rows
    assert all(row["fdr_q"] == "" for row in rows)
    assert all(0 <= float(row["p_value"]) <= 1 for row in rows)


@pytest.mark.parametrize("tail", ["enrichment", "depletion"])
@pytest.mark.parametrize("complete", [False, True])
@pytest.mark.parametrize("prefix", ["", "\n", "# Producer: fixture\n\n"])
def test_calibrated_figure_refuses_nominal_input_before_rendering(
    tmp_path: Path,
    tail: str,
    complete: bool,
    prefix: str,
) -> None:
    source = tmp_path / "nominal.tsv"
    source.write_text(
        prefix
        + "# FDR: False\n# FDR-Scope: none\ntrait\tfdr_q"
        + ("\tis_significant" if complete else "")
        + "\nTraitA\t"
        + ("\tTrue" if complete else "")
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="FDR"):
        generate_figure2_calibrated_ora_infographic(
            tmp_path / "matrix.tsv",
            tmp_path / "experiments.tsv",
            tmp_path / "figure.png",
            enrichment_results_path=source if tail == "enrichment" else None,
            depletion_results_path=source if tail == "depletion" else None,
        )
    assert not (tmp_path / "figure.png").exists()


@pytest.mark.parametrize("prefix", ["", "\n", "# Producer: fixture\n\n"])
def test_ora_figure_guard_refuses_explicit_none_scope(prefix: str) -> None:
    with pytest.raises(ValueError, match="FDR"):
        refuse_nominal_ora(io.StringIO(prefix + "# FDR-Scope: none\ntrait\tfdr_q\n"))


@pytest.mark.parametrize("metadata", ["# Stat-Type: p_value\n", "# FDR: False\n"])
def test_plate_rankings_refuse_explicit_nominal_hits(tmp_path: Path, metadata: str) -> None:
    plate = _script("generate_tufte_plates")
    path = tmp_path / "rankings.tsv"
    path.write_text(metadata + "trait\tinc_hits\tdec_hits\thit_rows\nTraitA\t3\t1\t4\n")
    with pytest.raises(ValueError, match="FDR"):
        plate.read_tsv(path, adjusted_ora=True)


@pytest.mark.parametrize("metadata", ["# FDR: True\n", "# Stat-Type: fdr_q\n"])
def test_plate_rankings_preserve_existing_adjusted_rows(tmp_path: Path, metadata: str) -> None:
    plate = _script("generate_tufte_plates")
    path = tmp_path / "rankings.tsv"
    path.write_text(metadata + "trait\tinc_hits\tdec_hits\thit_rows\nTraitA\t3\t1\t4\n")
    assert plate.read_tsv(path, adjusted_ora=True) == [
        {"trait": "TraitA", "inc_hits": "3", "dec_hits": "1", "hit_rows": "4"}
    ]


@pytest.mark.parametrize(
    "text",
    [
        "",
        "\n\n",
        "# Producer: fixture\n",
        "trait\tfdr_q\nTraitA\t\n",
        "\n# Producer: fixture\n\ntrait\tfdr_q\nTraitA\t\n",
        "# FDR-Scope: global\ntrait\tfdr_q\nTraitA\t0.01\n",
    ],
    ids=["empty", "blank", "comment-only", "headerless", "generic-header", "scope-only"],
)
def test_ora_figure_guard_requires_an_adjusted_declaration(text: str) -> None:
    """A table or a generic comment cannot establish FDR control (#1504)."""
    with pytest.raises(ValueError, match="FDR"):
        refuse_nominal_ora(io.StringIO(text))


@pytest.mark.parametrize("metadata", ["# FDR: True\n", "# Stat-Type: fdr_q\n"])
def test_ora_figure_guard_preserves_affirmative_adjusted_headers(metadata: str) -> None:
    refuse_nominal_ora(io.StringIO("\n# Producer: fixture\n\n" + metadata + "trait\tfdr_q\n"))


@pytest.mark.parametrize(
    "metadata", ["# FDR: False\n", "# Stat-Type: p_value\n", "# FDR-Scope: none\n"]
)
def test_ora_figure_guard_does_not_stop_after_an_adjusted_header(metadata: str) -> None:
    with pytest.raises(ValueError, match="Nominal ORA input is not FDR-controlled"):
        refuse_nominal_ora(io.StringIO("# FDR: True\n" + metadata + "trait\tfdr_q\n"))
