"""Readers of the older pipeline's tables check their strain rollup records (#847)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

from mtp import taxon_contributions
from mtp.paths import REPO_ROOT
from mtp.strain_rollup_record import (
    StrainRollupRecordError,
    recorded_strain_rollup,
    write_strain_rollup_record,
)
from mtp.stratification import run_stratification_pipeline
from mtp.taxa_trait_matrix import StrainRollup


def _recorded(path: Path, scope: StrainRollup) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("study_id\texperiment_id\tinc_dec\tMETPO:1\n", encoding="utf-8")
    write_strain_rollup_record(path, scope)
    return path


def _load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        f"strain_rollup_record_{name}", REPO_ROOT / "scripts" / f"{name}.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_ora_enrichment_refuses_profiles_and_matrix_of_different_scopes(
    tmp_path: Path,
) -> None:
    enrichment = _load_script("trait_enrichment")
    profiles = _recorded(tmp_path / "profiles.tsv", "off")
    matrix = _recorded(tmp_path / "matrix.tsv", "phenotype")
    output = tmp_path / "report.txt"

    with pytest.raises(StrainRollupRecordError, match="different strain rollup scopes"):
        enrichment.main(
            ["--input", str(profiles), "--taxa-matrix", str(matrix), "--output", str(output)]
        )
    assert not output.exists()


def test_ora_enrichment_refuses_an_unrecorded_table_read_with_a_recorded_one(
    tmp_path: Path,
) -> None:
    """The committed table decompressed by hand, without its record, next to a recorded
    matrix may not share its scope (#939)."""
    enrichment = _load_script("trait_enrichment")
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text("study_id\texperiment_id\tinc_dec\tMETPO:1\n", encoding="utf-8")
    matrix = _recorded(tmp_path / "matrix.tsv", "phenotype")
    argv = ["--input", str(profiles), "--taxa-matrix", str(matrix)]

    with pytest.raises(StrainRollupRecordError, match=r"profiles\.tsv \(no record\)"):
        enrichment.main([*argv, "--output", str(tmp_path / "report.txt")])


def test_ora_strata_record_each_subset_under_the_scope_of_its_tables(tmp_path: Path) -> None:
    """Each stratum's enrichment reads its subset with the matrix, so the subset carries
    the scope the full tables record (#939)."""
    matrix = tmp_path / "taxa_trait_matrix.tsv"
    matrix.write_text(
        "taxon_id\tTrait1\tTrait2\nNCBITaxon:1\t1\t1\nNCBITaxon:2\t1\t0\n"
        "NCBITaxon:3\t0\t1\nNCBITaxon:4\t1\t1\n",
        encoding="utf-8",
    )
    profiles = tmp_path / "study_trait_profiles_ora_wide.tsv"
    profiles.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTrait1\tTrait2\n"
        "s1\te1\tinc\t4\t3\t2\ns1\te1\tdec\t4\t1\t3\n"
        "s2\te2\tinc\t4\t3\t1\ns2\te2\tdec\t4\t2\t2\n",
        encoding="utf-8",
    )
    for table in (matrix, profiles):
        write_strain_rollup_record(table, "phenotype")
    experiments = tmp_path / "experiments.tsv"
    experiments.write_text(
        "study_id\texperiment_id\tbody_site\ns1\te1\tfeces\ns2\te2\tsaliva\n",
        encoding="utf-8",
    )

    outputs, _, _ = run_stratification_pipeline(
        profiles_path=profiles,
        experiments_path=experiments,
        output_dir=tmp_path / "strata",
        p_threshold=1.0,
        min_fold=0.0,
        fdr=False,
        strata=["feces"],
        taxa_matrix=matrix,
    )

    subset = tmp_path / "strata" / "study_trait_profiles_feces.tsv"
    assert outputs["feces"].overall_ranked_tsv.exists()
    assert recorded_strain_rollup(subset) == "phenotype"


def test_ora_strata_refuse_profiles_and_matrix_of_different_scopes(tmp_path: Path) -> None:
    """The strata write unrecorded subsets, so the full tables are checked first."""
    profiles = _recorded(tmp_path / "profiles.tsv", "off")
    matrix = _recorded(tmp_path / "matrix.tsv", "phenotype")
    output = tmp_path / "strata"

    with pytest.raises(StrainRollupRecordError, match="different strain rollup scopes"):
        run_stratification_pipeline(
            profiles_path=profiles,
            experiments_path=tmp_path / "experiments.tsv",
            output_dir=output,
            taxa_matrix=matrix,
        )
    assert not output.exists()


def test_the_ledger_refuses_inputs_built_under_another_scope(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """A mismatch is named before the KG is loaded, not after 32,220 count mismatches."""
    inputs = tmp_path / "derived"
    for name in ("study_trait_profiles_wide.tsv", "taxa_trait_matrix.tsv"):
        _recorded(inputs / name, "phenotype")
    output = tmp_path / "ledger"
    argv = ["--kg-dir", str(tmp_path / "kg"), "--inputs", str(inputs), "--output", str(output)]
    monkeypatch.setattr(sys, "argv", ["taxon_contributions", *argv, "--strain-rollup", "off"])

    with pytest.raises(SystemExit) as exited:
        taxon_contributions.main()

    assert exited.value.code == 2
    assert "phenotype strain rollup scope, not the declared off" in capsys.readouterr().err
    assert not output.exists()
