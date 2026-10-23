"""Readers of the older pipeline's tables check their strain rollup records (#847)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

from mtp import taxon_contributions
from mtp.paths import REPO_ROOT
from mtp.strain_rollup_record import StrainRollupRecordError, write_strain_rollup_record
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
