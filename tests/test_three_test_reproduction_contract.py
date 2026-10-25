"""Synthetic reproduction-interface regressions for issues #1060 and #1061."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

import test_three_test_reproduction as reproduction

_ENVIRONMENT = (
    "MTP_THREE_TEST_COMPARISON",
    "MTP_THREE_TEST_OUTPUT",
    "MTP_THREE_TEST_SCRATCH",
    "MTP_THREE_TEST_BASIS",
    "MTP_THREE_TEST_MULTIPLICITY",
    "MTP_THREE_TEST_LABEL_SUPPORT",
)


@pytest.fixture
def reproduction_environment(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> dict[str, Path]:
    for variable in _ENVIRONMENT:
        monkeypatch.delenv(variable, raising=False)
    paths = {
        "comparison": tmp_path / "comparison",
        "output": tmp_path / "output",
        "scratch": tmp_path / "scratch",
        "multiplicity": tmp_path / "multiplicity",
        "label_support": tmp_path / "label support",
    }
    for name in ("comparison", "output", "scratch"):
        monkeypatch.setenv(f"MTP_THREE_TEST_{name.upper()}", str(paths[name]))
    monkeypatch.setenv("MTP_THREE_TEST_BASIS", "511c891")
    return paths


def _reference_record() -> dict[str, Any]:
    """Synthetic values exercise assertion routing; no scientific run is represented."""
    return {
        "strain_rollup": "off",
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "family_size": 317,
        "source_commits": {"core": "511c891" + "0" * 33},
    }


def _reference_table(*args: Any, **kwargs: Any) -> pd.DataFrame:
    del args, kwargs
    return pd.DataFrame(
        [
            {
                "trait": "METPO:2000103|GO:0030245" if index == 0 else f"synthetic:{index}",
                "composition_studies": 548,
                "composition_pairs": 935,
                "study_t_status": "tested",
                "study_t_statistic": -5.24,
                "study_t_q_by": 0.00047,
            }
            for index in range(317)
        ]
    )


@pytest.mark.parametrize(
    "auxiliary", [(), ("multiplicity",), ("label_support",), ("multiplicity", "label_support")]
)
def test_reproduction_forwards_original_auxiliary_evidence(
    reproduction_environment: dict[str, Path],
    monkeypatch: pytest.MonkeyPatch,
    auxiliary: tuple[str, ...],
) -> None:
    paths = reproduction_environment
    for name in auxiliary:
        monkeypatch.setenv(f"MTP_THREE_TEST_{name.upper()}", str(paths[name]))
    calls: list[tuple[Path, Path, dict[str, Any]]] = []

    def verify(output: Path, comparison: Path, **kwargs: Any) -> dict[str, Any]:
        calls.append((output, comparison, kwargs))
        return _reference_record()

    monkeypatch.setattr(reproduction, "verify_table", verify)
    monkeypatch.setattr(pd, "read_csv", _reference_table)
    reproduction.test_verified_off_cellulose_target()
    assert calls == [
        (
            paths["output"],
            paths["comparison"],
            {
                "scratch": paths["scratch"],
                "multiplicity": paths["multiplicity"] if "multiplicity" in auxiliary else None,
                "label_support": paths["label_support"] if "label_support" in auxiliary else None,
            },
        )
    ]


@pytest.mark.parametrize("basis", [None, "", "abc", "511c89g", "0" * 41, " 511c891"])
def test_invalid_basis_is_rejected_before_replay(
    reproduction_environment: dict[str, Path],
    monkeypatch: pytest.MonkeyPatch,
    basis: str | None,
) -> None:
    del reproduction_environment
    if basis is None:
        monkeypatch.delenv("MTP_THREE_TEST_BASIS")
    else:
        monkeypatch.setenv("MTP_THREE_TEST_BASIS", basis)

    def fail_replay(*args: Any, **kwargs: Any) -> dict[str, Any]:
        del args, kwargs
        pytest.fail("Invalid basis reached expensive replay")

    monkeypatch.setattr(reproduction, "verify_table", fail_replay)
    with pytest.raises(ValueError, match="basis"):
        reproduction.test_verified_off_cellulose_target()


def test_missing_explicit_input_is_not_a_successful_skip(
    reproduction_environment: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    del reproduction_environment
    monkeypatch.delenv("MTP_THREE_TEST_SCRATCH")
    with pytest.raises(ValueError, match="MTP_THREE_TEST_SCRATCH"):
        reproduction.test_verified_off_cellulose_target()


def test_direct_integration_without_inputs_remains_an_explicit_skip(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    for variable in _ENVIRONMENT:
        monkeypatch.delenv(variable, raising=False)
    with pytest.raises(pytest.skip.Exception, match="verified archived off comparison"):
        reproduction.test_verified_off_cellulose_target()


@pytest.mark.parametrize("with_auxiliary", [False, True])
def test_reproduction_recipe_accepts_optional_evidence_paths(
    tmp_path: Path, with_auxiliary: bool
) -> None:
    arguments = [str(tmp_path / name) for name in ("comparison", "output", "scratch")]
    if with_auxiliary:
        arguments.extend(
            ["511c891", str(tmp_path / "multiplicity run"), str(tmp_path / "label support")]
        )
    result = subprocess.run(
        ["just", "--dry-run", "test-three-test-off-reproduction", *arguments],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
        check=False,
        timeout=10,
    )
    assert result.returncode == 0, result.stderr
    command = result.stdout + result.stderr
    assert 'MTP_THREE_TEST_MULTIPLICITY="$5"' in command
    assert 'MTP_THREE_TEST_LABEL_SUPPORT="$6"' in command
