"""Real CLI publication applies inherited permissions without re-importing main (#1088)."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from mtp.three_test_table.stage import require_durable


@pytest.mark.parametrize("mask", [0o002, 0o022, 0o027, 0o077], ids=["002", "022", "027", "077"])
def test_cli_publication_inherits_umask_without_main_reimport(mask: int) -> None:
    root = Path(__file__).resolve().parents[1]
    scratch = require_durable(root / ".tmp" / "three-test-publication-tests")
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=scratch) as temporary:
        output_root = Path(temporary) / "run"
        process = subprocess.run(
            [
                sys.executable,
                str(root / "scripts/adhoc/issue1088_cli_publication.py"),
                "--output-root",
                str(output_root),
            ],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
            timeout=180,
            umask=mask,
        )
        assert process.returncode == 0, process.stderr
        assert "ISSUE1088_FORBIDDEN_MAIN_REIMPORT" not in process.stderr
        report = json.loads((output_root / "result.json").read_text(encoding="utf-8"))
        assert report["cli_exit_code"] == 0
        assert report["output_mode"] == 0o777 & ~mask
        assert report["parent_mode_before"] == report["parent_mode_after"]
        assert report["siblings"] == ["output"]
        assert report["trait_count"] == report["expected_trait_count"] == 7
        assert report["unique_traits"] == report["trait_count"]
        assert report["family_sizes"] == ["7"]
