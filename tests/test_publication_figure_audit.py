"""Export-audit inventory and exact-byte regressions for #1198/#1199."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
import yaml
from PIL import Image

from mtp.workflow import output_roots

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/adhoc/issue1046_check_figure_bundle.py"
SPEC = importlib.util.spec_from_file_location("figure_audit_under_test", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load figure audit helper")
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


def identity(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


@dataclass
class Bundle:
    root: Path
    output: Path
    reports: Path
    record: dict[str, Any]

    @property
    def provenance(self) -> Path:
        return self.root / "forest/provenance.json"

    @property
    def report(self) -> Path:
        return self.reports / "forest-page-001.yaml"

    def seal(self) -> None:
        self.provenance.write_text(json.dumps(self.record), encoding="utf-8")
        root_record = {
            "version": "publication-figures-v1",
            "status": "completed",
            "figures": {"forest": {"path": "forest/provenance.json", **identity(self.provenance)}},
        }
        (self.root / "figures.json").write_text(json.dumps(root_record), encoding="utf-8")

    def audit(self, reports: bool = False) -> dict[str, Any]:
        result: dict[str, Any] = HELPER.audit(
            self.root, self.output, self.reports if reports else None
        )
        return result


@pytest.fixture
def bundle(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Bundle:
    monkeypatch.setattr(
        output_roots, "system_temporary_directories", lambda: [tmp_path / "volatile"]
    )
    root, reports = tmp_path / "bundle", tmp_path / "reports"
    figure = root / "forest"
    figure.mkdir(parents=True)
    reports.mkdir()
    for name in ("forest.pdf", "page-001.pdf", "page-001.svg"):
        (figure / name).write_bytes(f"synthetic export {name}".encode())
    with Image.new("RGB", (30, 30), "white") as image:
        image.save(figure / "page-001.png", dpi=(300, 300))
    record = {
        "version": "publication-figures-v1",
        "figure": "forest",
        "code_commit": "a" * 40,
        "code_dirty": False,
        "rendering": {
            "dpi": 300,
            "page_count": 1,
            "vector_formats": ["pdf", "svg"],
            "raster_formats": ["png"],
        },
        "pages": [{"page": 1, "row_ids": ["METPO:1"], "mark_count": 2, "details": {}}],
        "outputs": {path.name: identity(path) for path in figure.iterdir()},
    }
    result = Bundle(root, tmp_path / "audit", reports, record)
    result.seal()
    report = {
        "target": {"name": "journal_double_column", "final_width_in": 7.2},
        "verdict": "pass",
        "inputs": [
            {"format": suffix, "sha256": identity(figure / f"page-001.{suffix}")["sha256"]}
            for suffix in ("pdf", "png")
        ],
    }
    result.report.write_text(yaml.safe_dump(report), encoding="utf-8")
    return result


def test_valid_bundle_and_shell_nested_audit_contract(bundle: Bundle) -> None:
    result = bundle.audit()
    assert result["figures"]["forest"]["rows"] == 1
    assert result["figures"]["forest"]["marks"] == 2
    assert result["bundle_manifest"] == identity(bundle.root / "figures.json")
    (bundle.output / "design-check").mkdir()
    validated = HELPER.audit(bundle.root, bundle.output / "validated", bundle.reports)
    assert validated["figures"]["forest"]["design_check_verdicts"] == {"pass": 1}


@pytest.mark.parametrize("name", ["forest.pdf", "page-001.pdf", "page-001.svg", "page-001.png"])
def test_declared_and_disk_export_omission_refused(bundle: Bundle, name: str) -> None:
    del bundle.record["outputs"][name]
    (bundle.provenance.parent / name).unlink()
    bundle.seal()
    with pytest.raises(ValueError, match=r"inventory|export"):
        bundle.audit()


@pytest.mark.parametrize("kind", ["empty", "unknown", "wrong_path"])
def test_root_figure_inventory_refused(bundle: Bundle, kind: str) -> None:
    path = bundle.root / "figures.json"
    record = json.loads(path.read_bytes())
    if kind == "empty":
        record["figures"] = {}
    elif kind == "unknown":
        record["figures"]["unexpected"] = record["figures"].pop("forest")
    else:
        record["figures"]["forest"]["path"] = "forest/../forest/provenance.json"
    path.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match=r"inventory|figure|path"):
        bundle.audit()


@pytest.mark.parametrize("kind", ["empty", "duplicate", "noncontiguous", "page_count"])
def test_page_inventory_refused(bundle: Bundle, kind: str) -> None:
    if kind == "empty":
        bundle.record["pages"] = []
    elif kind == "duplicate":
        bundle.record["pages"].append(bundle.record["pages"][0].copy())
    elif kind == "noncontiguous":
        bundle.record["pages"][0]["page"] = 2
    else:
        bundle.record["rendering"]["page_count"] = 2
    bundle.seal()
    with pytest.raises(ValueError, match=r"page|inventory"):
        bundle.audit()


@pytest.mark.parametrize("location", ["root", "figure"])
def test_unrecorded_bundle_member_refused(bundle: Bundle, location: str) -> None:
    directory = bundle.root if location == "root" else bundle.provenance.parent
    (directory / "unrecorded.txt").write_text("extra", encoding="utf-8")
    with pytest.raises(ValueError, match="inventory"):
        bundle.audit()


def test_symlink_export_refused(bundle: Bundle) -> None:
    path = bundle.provenance.parent / "page-001.svg"
    alternate = bundle.root.parent / "external.svg"
    path.rename(alternate)
    path.symlink_to(alternate)
    with pytest.raises(ValueError, match=r"symlink|regular"):
        bundle.audit()


def test_provenance_hash_binds_parsed_bytes(
    bundle: Bundle, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = Path.read_bytes
    changed = {**bundle.record, "code_commit": "b" * 40}

    def read_bytes(path: Path) -> bytes:
        return json.dumps(changed).encode() if path == bundle.provenance else original(path)

    monkeypatch.setattr(Path, "read_bytes", read_bytes)
    with pytest.raises(ValueError, match=r"identity|hash|changed"):
        bundle.audit()


def test_root_change_after_parse_refused(bundle: Bundle, monkeypatch: pytest.MonkeyPatch) -> None:
    original = Path.read_bytes
    root_manifest = bundle.root / "figures.json"
    changed = False

    def read_bytes(path: Path) -> bytes:
        nonlocal changed
        raw = original(path)
        if path == root_manifest and not changed:
            changed = True
            path.write_bytes(raw + b"\n")
        return raw

    monkeypatch.setattr(Path, "read_bytes", read_bytes)
    with pytest.raises(ValueError, match=r"identity|hash|changed"):
        bundle.audit()


def test_export_change_after_validation_refused(
    bundle: Bundle, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = HELPER.Image.open
    target = bundle.provenance.parent / "page-001.svg"

    def open_image(*args: Any, **kwargs: Any) -> Any:
        image = original(*args, **kwargs)
        target.write_bytes(b"changed after export checks")
        return image

    monkeypatch.setattr(HELPER.Image, "open", open_image)
    with pytest.raises(ValueError, match=r"identity|hash|changed"):
        bundle.audit()


def test_report_change_after_parse_refused(bundle: Bundle, monkeypatch: pytest.MonkeyPatch) -> None:
    original = HELPER.yaml.safe_load

    def safe_load(*args: Any, **kwargs: Any) -> Any:
        result = original(*args, **kwargs)
        bundle.report.write_text("verdict: fail\n", encoding="utf-8")
        return result

    monkeypatch.setattr(HELPER.yaml, "safe_load", safe_load)
    with pytest.raises(ValueError, match=r"identity|hash|changed"):
        bundle.audit(reports=True)


def test_design_report_identity_is_retained(bundle: Bundle) -> None:
    result = bundle.audit(reports=True)
    assert result["design_reports"] == {bundle.report.name: identity(bundle.report)}


def test_validation_failure_does_not_create_output(bundle: Bundle) -> None:
    bundle.record["outputs"]["forest.pdf"]["sha256"] = "f" * 64
    bundle.seal()
    with pytest.raises(ValueError, match=r"identity|hash"):
        bundle.audit()
    assert not bundle.output.exists()


@pytest.mark.parametrize("bad_identity", [None, {}, {"sha256": "a" * 64, "size_bytes": None}])
def test_malformed_export_identity_refused(bundle: Bundle, bad_identity: Any) -> None:
    bundle.record["outputs"]["forest.pdf"] = bad_identity
    bundle.seal()
    with pytest.raises(ValueError, match=r"identity|inventory"):
        bundle.audit()
