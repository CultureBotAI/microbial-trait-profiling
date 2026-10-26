"""Declared inputs and command entry point for the #1052 derived analyses."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from unittest.mock import Mock

import pytest

from mtp.delta_z_derived import __main__ as cli


def arguments() -> list[str]:
    return [
        "--three-test-run",
        "parent/three_test_run.json",
        "--three-test-sha256",
        "a" * 64,
        "--comparison-run",
        "comparison/comparison_run.json",
        "--output",
        "derived",
    ]


def test_full_command_passes_explicit_trait_list(monkeypatch: pytest.MonkeyPatch) -> None:
    run = Mock(return_value={"run_id": "derived-fixture", "family_size": 10})
    monkeypatch.setattr(cli, "run_derived", run)
    cli.main([*arguments(), "--trait-list", "ten.txt", "--code-commit", "b" * 40])
    run.assert_called_once_with(
        Path("parent/three_test_run.json"),
        Path("comparison/comparison_run.json"),
        Path("derived"),
        expected_sha256="a" * 64,
        trait_list=Path("ten.txt"),
        summary_only=False,
        code_commit="b" * 40,
        top_diseases=7,
    )


def test_off_summary_does_not_require_or_infer_a_trait_list(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    run = Mock(return_value={"run_id": "derived-off", "family_size": 317})
    monkeypatch.setattr(cli, "run_derived", run)
    cli.main([*arguments(), "--summary-only"])
    assert run.call_args.kwargs["trait_list"] is None
    assert run.call_args.kwargs["summary_only"] is True


@pytest.mark.parametrize(
    "extra",
    [
        [],
        ["--summary-only", "--trait-list", "ten.txt"],
        ["--trait-list", "ten.txt", "--top-diseases", "0"],
        ["--trait-list", "ten.txt", "--top-diseases", "-1"],
    ],
)
def test_ambiguous_or_invalid_analysis_selection_fails_before_run(
    monkeypatch: pytest.MonkeyPatch,
    extra: list[str],
) -> None:
    run = Mock()
    monkeypatch.setattr(cli, "run_derived", run)
    with pytest.raises(SystemExit) as error:
        cli.main([*arguments(), *extra])
    assert error.value.code == 2
    run.assert_not_called()


def test_parent_manifest_hash_is_an_explicit_required_input(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    run = Mock()
    monkeypatch.setattr(cli, "run_derived", run)
    args = arguments()
    del args[2:4]
    with pytest.raises(SystemExit) as error:
        cli.main([*args, "--summary-only"])
    assert error.value.code == 2
    run.assert_not_called()


def test_real_module_help_names_inputs_and_summary_mode() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "mtp.delta_z_derived", "--help"],
        check=True,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert "--three-test-run" in result.stdout
    assert "--three-test-sha256" in result.stdout
    assert "--trait-list" in result.stdout
    assert "--summary-only" in result.stdout
