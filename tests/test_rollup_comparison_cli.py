"""Real recipe forwarding and required scope pins for #1041."""

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

from mtp.rollup_comparison import __main__ as cli

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    root = ROOT / ".tmp" / "rollup-cli-tests"
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as directory:
        yield Path(directory)


def arguments() -> list[str]:
    result: list[str] = []
    for scope in ("off", "phenotype"):
        result.extend(
            [
                f"--{scope}-three-test-run",
                f"{scope} parent/three_test_run.json",
                f"--{scope}-sha256",
                scope[0] * 64,
                f"--{scope}-comparison-run",
                f"{scope} source/comparison_run.json",
            ]
        )
    return [*result, "--output", "new output/$(printf literal); result", "--code-commit", "a" * 40]


def test_cli_passes_each_pin_to_the_correct_parent(monkeypatch: pytest.MonkeyPatch) -> None:
    run = Mock(return_value={"run_id": "fixture", "counts": {"shared_traits": 2}})
    monkeypatch.setattr(cli, "run_comparison", run)
    assert cli.main(arguments()) == 0
    run.assert_called_once_with(
        Path("off parent/three_test_run.json"),
        Path("off source/comparison_run.json"),
        Path("phenotype parent/three_test_run.json"),
        Path("phenotype source/comparison_run.json"),
        Path("new output/$(printf literal); result"),
        off_sha256="o" * 64,
        phenotype_sha256="p" * 64,
        code_commit="a" * 40,
    )


@pytest.mark.parametrize("missing", ["--off-sha256", "--phenotype-sha256"])
def test_both_pins_are_mandatory(missing: str) -> None:
    argv = arguments()
    at = argv.index(missing)
    del argv[at : at + 2]
    with pytest.raises(SystemExit) as caught:
        cli.main(argv)
    assert caught.value.code == 2


def test_actual_just_recipe_preserves_all_arguments(tmp_path: Path) -> None:
    tools = tmp_path / "intercepted tools"
    tools.mkdir()
    writer = tools / "capture.py"
    writer.write_text(
        '"""Capture recipe argv; inputs environment/output JSON; #1041, 2026-10-25."""\n'
        "import json, os, sys\nfrom pathlib import Path\n"
        'Path(os.environ["MTP_ROLLUP_CAPTURE"]).write_text(json.dumps(sys.argv[1:]))\n',
        encoding="utf-8",
    )
    uv = tools / "uv"
    uv.write_text(
        '#!/bin/sh\nexec "$MTP_ROLLUP_PYTHON" "$MTP_ROLLUP_WRITER" "$@"\n', encoding="utf-8"
    )
    uv.chmod(0o755)
    capture = tmp_path / "capture.json"
    environment = {
        "PATH": os.pathsep.join([str(tools), os.environ.get("PATH", os.defpath)]),
        "MTP_ROLLUP_PYTHON": sys.executable,
        "MTP_ROLLUP_WRITER": str(writer),
        "MTP_ROLLUP_CAPTURE": str(capture),
    }
    result = subprocess.run(
        ["just", "--justfile", str(ROOT / "Justfile"), "rollup-comparison", *arguments()],
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    observed: Any = json.loads(capture.read_text())
    assert observed == ["run", "python", "-m", "mtp.rollup_comparison", *arguments()]


def test_actual_qualification_recipe_forwards_all_seven_inputs(tmp_path: Path) -> None:
    """Use the real recipe/shell with explicit input captures rather than a dry run."""
    tools = tmp_path / "intercepted qualification tools"
    tools.mkdir()
    writer = tools / "capture.py"
    writer.write_text(
        '"""Capture qualification argv/env; #1041, 2026-10-25."""\n'
        "import json, os, sys\nfrom pathlib import Path\n"
        'keys = json.loads(os.environ["MTP_ROLLUP_KEYS"])\n'
        'record = {"argv": sys.argv[1:], "inputs": [os.environ.get(k) for k in keys]}\n'
        'Path(os.environ["MTP_ROLLUP_CAPTURE"]).write_text(json.dumps(record))\n',
        encoding="utf-8",
    )
    uv = tools / "uv"
    uv.write_text(
        '#!/bin/sh\nexec "$MTP_ROLLUP_PYTHON" "$MTP_ROLLUP_WRITER" "$@"\n', encoding="utf-8"
    )
    uv.chmod(0o755)
    capture = tmp_path / "capture.json"
    keys = [
        "MTP_ROLLUP_OFF_PARENT",
        "MTP_ROLLUP_PHENOTYPE_PARENT",
        "MTP_ROLLUP_RESULT",
        "MTP_ROLLUP_COMMIT",
        "MTP_ROLLUP_SOURCE_COMMIT",
        "MTP_ROLLUP_OFF_SHA256",
        "MTP_ROLLUP_PHENOTYPE_SHA256",
    ]
    values = [
        "off parent/file.json",
        "phenotype parent/file.json",
        "result $(printf literal); dir",
        "a" * 40,
        "b" * 40,
        "c" * 64,
        "d" * 64,
    ]
    result = subprocess.run(
        ["just", "--justfile", str(ROOT / "Justfile"), "test-rollup-reproduction", *values],
        cwd=ROOT,
        env={
            "PATH": os.pathsep.join([str(tools), os.environ.get("PATH", os.defpath)]),
            "MTP_ROLLUP_PYTHON": sys.executable,
            "MTP_ROLLUP_WRITER": str(writer),
            "MTP_ROLLUP_CAPTURE": str(capture),
            "MTP_ROLLUP_KEYS": json.dumps(keys),
        },
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    observed = json.loads(capture.read_text())
    assert observed["inputs"] == values
    assert observed["argv"] == [
        "run",
        "pytest",
        "tests/test_rollup_comparison_reproduction.py",
        "-m",
        "integration",
        "-s",
    ]
