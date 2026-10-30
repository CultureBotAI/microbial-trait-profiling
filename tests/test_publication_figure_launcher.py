"""System-shell coverage of the pinned F figure launcher (#1226)."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from collections.abc import Iterator
from pathlib import Path

import pytest

REPOSITORY = Path(__file__).resolve().parents[1]
LAUNCHER = REPOSITORY / "scripts/adhoc/issue1046_slice2_f.sh"
FIGURES = ("rank-sensitivity", "candidate-triples", "trait-overlap")
CAPTURE_SOURCE = '''"""Capture launcher arguments for #1226, 2026-10-29.

Inputs: intercepted just arguments and an explicit capture destination.
Output: JSON arguments only; no scientific files or analysis commands are used.
"""
import json
import os
import sys
from pathlib import Path

Path(os.environ["MTP_LAUNCH_CAPTURE"]).write_text(json.dumps(sys.argv[1:]), encoding="utf-8")
'''


@pytest.fixture
def launcher_directory() -> Iterator[Path]:
    """Keep generated interception helpers beneath the repository."""
    root = REPOSITORY / ".tmp/publication-figure-launcher-tests"
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as directory:
        yield Path(directory)


@pytest.mark.parametrize("figure", ["all", *FIGURES])
def test_f_launcher_forwards_figures_under_system_bash(
    launcher_directory: Path, figure: str
) -> None:
    """Exercise actual nounset/array expansion, including macOS Bash 3.2."""
    bash = Path("/bin/bash")
    if not bash.is_file():
        pytest.skip("System /bin/bash is unavailable")
    tools = launcher_directory / "intercepted tools"
    tools.mkdir()
    writer = tools / "capture arguments.py"
    writer.write_text(CAPTURE_SOURCE, encoding="utf-8")
    intercepted = tools / "just"
    intercepted.write_text(
        "#!/bin/sh\n"
        "# Intercept just for the #1226 shell regression, 2026-10-29.\n"
        'exec "$MTP_LAUNCH_PYTHON" "$MTP_LAUNCH_WRITER" "$@"\n',
        encoding="utf-8",
    )
    intercepted.chmod(0o755)
    capture = launcher_directory / "captured arguments.json"
    evidence = launcher_directory / "evidence with spaces $(printf literal); source"
    output = launcher_directory / "fresh output with spaces"
    result = subprocess.run(
        [str(bash), str(LAUNCHER), str(evidence), str(output), figure],
        cwd=REPOSITORY,
        env={
            "PATH": os.pathsep.join((str(tools), os.defpath)),
            "MTP_LAUNCH_PYTHON": sys.executable,
            "MTP_LAUNCH_WRITER": str(writer),
            "MTP_LAUNCH_CAPTURE": str(capture),
        },
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert capture.is_file(), "Launcher did not reach the intercepted just command"
    observed = json.loads(capture.read_text(encoding="utf-8"))
    assert isinstance(observed, list)
    assert all(isinstance(item, str) for item in observed)
    assert observed[0] == "publication-figures-slice2"
    selected = [observed[index + 1] for index, value in enumerate(observed) if value == "--figure"]
    assert selected == (list(FIGURES) if figure == "all" else [figure])
    assert observed[-2:] == ["--output", str(output)]
    assert observed[observed.index("--scan-run") + 1] == str(
        evidence / "mtp-runs/F-511c891/additive/tier12-phenotype/scan835"
    )
    assert not evidence.exists()
    assert not output.exists()
