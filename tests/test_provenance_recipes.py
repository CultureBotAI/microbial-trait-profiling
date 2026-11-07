"""Provenance recipes forward root and attestation options without a real run."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "recipe", ["audit-ora-refresh", "ora-refresh-inputs", "verify-integrated-provenance"]
)
def test_provenance_recipe_forwards_options(tmp_path: Path, recipe: str) -> None:
    capture = tmp_path / "arguments.txt"
    shim = tmp_path / "uv"
    shim.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "$PROVENANCE_CAPTURE"\n')
    shim.chmod(0o700)
    arguments = ["source record.json"]
    if recipe != "verify-integrated-provenance":
        arguments.append("output record.json")
    arguments += ["--provenance-root", "source root"]
    if recipe == "ora-refresh-inputs":
        arguments += ["--expected-audit", "previous audit.json"]
    result = subprocess.run(
        [
            "just",
            "--justfile",
            str(Path(__file__).resolve().parents[1] / "Justfile"),
            recipe,
            *arguments,
        ],
        env={
            **os.environ,
            "PATH": str(tmp_path) + os.pathsep + os.environ["PATH"],
            "PROVENANCE_CAPTURE": str(capture),
        },
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert capture.is_file()
    actual = capture.read_text().splitlines()
    assert actual[actual.index("--provenance-root") + 1] == "source root"
    if recipe == "ora-refresh-inputs":
        assert actual[actual.index("--expected-audit") + 1] == "previous audit.json"
        assert "--write-portable-manifest" in actual
    elif recipe == "verify-integrated-provenance":
        assert actual[:4] == ["run", "python", "-m", "mtp.integrated_provenance"]
        assert actual[4] == "source record.json"
