"""Declared inputs and command entry point for the #1052 derived analyses."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pytest

from mtp.delta_z_derived import __main__ as cli

REPOSITORY = Path(__file__).resolve().parents[1]
REPRODUCTION_ENV = (
    "MTP_DZ_DERIVED_OFF",
    "MTP_DZ_DERIVED_PHENOTYPE",
    "MTP_DZ_HISTORICAL_STATS",
    "MTP_DZ_DERIVED_COMMIT",
    "MTP_DZ_CORE_EXPERIMENTS",
    "MTP_DZ_HISTORICAL_SHA256",
)
CAPTURE_SOURCE = '''"""Capture recipe argv and selected environment for #1099 (2026-10-25).

Inputs: intercepted uv arguments and explicit test environment keys.
Output: the test's JSON capture file; no analysis or external commands run.
"""

import json
import os
import sys
from pathlib import Path

keys = json.loads(os.environ["MTP_DZ_TEST_ENV_KEYS"])
record = {
    "argv": sys.argv[1:],
    "environment": {name: os.environ.get(name) for name in keys},
}
Path(os.environ["MTP_DZ_TEST_CAPTURE"]).write_text(json.dumps(record), encoding="utf-8")
'''


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    """Keep generated command-interception fixtures under the repository (#1099)."""
    root = REPOSITORY / ".tmp" / "delta-z-recipe-tests"
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as directory:
        yield Path(directory)


@pytest.fixture
def just_capture(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    """Intercept only uv; execute the actual checked-out Justfile with its real shell."""
    tools = tmp_path / "intercepted tools"
    tools.mkdir()
    writer = tools / "capture recipe.py"
    writer.write_text(CAPTURE_SOURCE, encoding="utf-8")
    uv = tools / "uv"
    uv.write_text(
        "#!/bin/sh\n"
        "# Capture actual Justfile argument forwarding for #1099 (2026-10-25).\n"
        'exec "$MTP_DZ_TEST_PYTHON" "$MTP_DZ_TEST_WRITER" "$@"\n',
        encoding="utf-8",
    )
    uv.chmod(0o755)
    capture = tmp_path / "recipe capture.json"
    environment = {
        "PATH": os.pathsep.join([str(tools), os.environ.get("PATH", os.defpath)]),
        "MTP_DZ_TEST_PYTHON": sys.executable,
        "MTP_DZ_TEST_WRITER": str(writer),
        "MTP_DZ_TEST_CAPTURE": str(capture),
        "MTP_DZ_TEST_ENV_KEYS": json.dumps(REPRODUCTION_ENV),
    }
    return capture, environment


def run_just_recipe(
    capture: tuple[Path, dict[str, str]], recipe: str, arguments: list[str]
) -> dict[str, Any]:
    destination, environment = capture
    result = subprocess.run(
        ["just", "--justfile", str(REPOSITORY / "Justfile"), recipe, *arguments],
        cwd=REPOSITORY,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert destination.is_file(), "The actual recipe never reached the intercepted uv command"
    observed: dict[str, Any] = json.loads(destination.read_text(encoding="utf-8"))
    return observed


@pytest.mark.parametrize("summary_only", [False, True])
def test_just_recipe_forwards_delta_z_arguments_with_spaces(
    just_capture: tuple[Path, dict[str, str]], summary_only: bool
) -> None:
    arguments = [
        "--three-test-run",
        "parent with spaces/$(printf literal); three_test_run.json",
        "--three-test-sha256",
        "a" * 64,
        "--comparison-run",
        "comparison with spaces/comparison_run.json",
        "--output",
        "derived output with spaces",
        "--code-commit",
        "b" * 40,
    ]
    arguments.extend(
        ["--summary-only"]
        if summary_only
        else ["--trait-list", "declared traits/ten traits.txt", "--top-diseases", "3"]
    )
    observed = run_just_recipe(just_capture, "delta-z-derived", arguments)
    assert observed["argv"] == ["run", "python", "-m", "mtp.delta_z_derived", *arguments]


def test_just_recipe_forwards_all_six_reproduction_inputs_with_spaces(
    just_capture: tuple[Path, dict[str, str]],
) -> None:
    arguments = [
        "off summary with spaces",
        "phenotype bundle with spaces",
        "historical numbers/$(printf literal); reference file.json",
        "c" * 40,
        "verified core/experiments table.tsv",
        "d" * 64,
    ]
    observed = run_just_recipe(just_capture, "test-delta-z-reproduction", arguments)
    assert observed["argv"] == [
        "run",
        "pytest",
        "tests/test_delta_z_derived_reproduction.py",
        "-m",
        "integration",
        "-v",
    ]
    assert observed["environment"] == dict(zip(REPRODUCTION_ENV, arguments, strict=True))


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
