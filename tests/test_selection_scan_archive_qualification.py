"""Guard the synthetic archive proof's destination and expected population (#1080)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest


@pytest.fixture
def helper() -> ModuleType:
    path = Path(__file__).resolve().parents[1] / (
        "scripts/adhoc/issue1080_scan_archive_qualification.py"
    )
    spec = importlib.util.spec_from_file_location("issue1080_proof_test", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("root", ["../durable", "~/durable"])
def test_proof_requires_absolute_root_before_preparation(
    helper: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    root: str,
) -> None:
    def forbidden(_root: Path) -> None:
        pytest.fail("Relative proof root reached preparation")

    monkeypatch.setattr(helper, "prepare", forbidden)
    monkeypatch.setattr(sys, "argv", ["proof", "prepare", root])
    with pytest.raises(SystemExit):
        helper.main()
    assert "absolute output root" in capsys.readouterr().err


def proof_outputs(root: Path, pool_studies: int) -> None:
    """Minimal already-produced evidence for comparator tests, not an archive verifier."""
    (root / "producer.json").write_text(json.dumps({"commit": "fixture", "files": {}}))
    report: dict[str, Any] = {
        "protocol": "cellulose-decr-selection-scan-v2",
        "follows_protocol": False,
        "provenance": {"generated_at_utc": "fixture"},
        "scopes": [],
    }
    for scope, role in (("phenotype", "primary"), ("off", "declared_sensitivity")):
        report["scopes"].append(
            {
                "scope": scope,
                "role": role,
                "counts": {
                    "eligible_conditions": 5,
                    "triples": 10,
                    "traits": 2,
                    "pool_studies": pool_studies,
                },
            }
        )
        archive = root / "restored" / scope
        data = archive / "data"
        data.mkdir(parents=True)
        (archive / "archive.json").write_text("{}")
        consumer = data / "core_consumption.json"
        consumer.write_text("{}")
        (root / f"resolved-{scope}.txt").write_text(str(consumer))
        retired = root / "retired-original"
        scan = retired / f"scan-{scope}"
        scan.mkdir(parents=True)
        (scan / "core_consumption.json").write_text("{}")
        core = retired / f"core-{scope}" / "run-fixture"
        core.mkdir(parents=True)
        (core / "run.json").write_text(
            json.dumps({"git_commit": "fixture", "git_dirty": False, "run_id": scope})
        )
    for name in ("unmoved", "relocated"):
        output = root / name
        output.mkdir()
        (output / "selection_scan_scopes.tsv").write_text("trait_id\nfixture\n")
        (output / "selection_scan_scopes.json").write_text(json.dumps(report))


@pytest.mark.parametrize("count", [17, 18, 19])
def test_proof_pins_pool_studies_even_when_relocation_parity_holds(
    helper: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, count: int
) -> None:
    monkeypatch.setattr(helper, "producer", lambda: {"commit": "fixture", "files": {}})
    proof_outputs(tmp_path, count)
    if count == 18:
        helper.compare(tmp_path)
        assert json.loads((tmp_path / "parity.json").read_text())["passed"] is True
    else:
        with pytest.raises(ValueError, match="Unexpected synthetic scan population"):
            helper.compare(tmp_path)
        assert not (tmp_path / "parity.json").exists()
