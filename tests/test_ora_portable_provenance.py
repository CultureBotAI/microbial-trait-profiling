"""Portable ORA audit/input records retain verifiable byte identities (#1664)."""

from __future__ import annotations

import hashlib
import io
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from mtp.paths import file_sha256
from test_no_absolute_local_paths import ABSOLUTE_LOCAL_PATH_PATTERN
from test_ora_refresh_audit import _AUDIT, _nominal_refresh_pair, run_audit


@pytest.fixture(autouse=True)
def assigned_source() -> None:
    expected = Path(__file__).resolve().parents[1] / "scripts/adhoc/issue720_audit_ora_refresh.py"
    assert _AUDIT.__file__ is not None
    assert Path(_AUDIT.__file__).resolve() == expected


def legacy_manifest(root: Path, *, absolute: bool = False) -> Path:
    """Use the independent small statistical fixture, without running a producer."""
    inputs = root / "inputs"
    inputs.mkdir(parents=True)
    pair = _nominal_refresh_pair(inputs, "none")
    for key in ("historical", "current", "profiles", "taxa_matrix"):
        pair[key] = str(inputs / pair[key]) if absolute else f"inputs/{pair[key]}"
    manifest = root / "legacy.json"
    manifest.write_text(json.dumps({"pairs": [pair]}), encoding="utf-8")
    return manifest


def attest_sources(manifest: Path) -> Path:
    """Retain the source hashes established before a manifest conversion."""
    audit = manifest.with_name("expected-audit.json")
    if not audit.exists():
        result = _AUDIT.audit_manifest(manifest, audit)
        assert result["passed"]
        audit.write_text(json.dumps(result), encoding="utf-8")
    return audit


def assert_portable(path: Path, original_root: Path) -> dict[str, Any]:
    payload = path.read_text(encoding="utf-8")
    assert not ABSOLUTE_LOCAL_PATH_PATTERN.search(payload), payload
    assert str(original_root.resolve()) not in payload, payload
    raw = json.loads(payload)
    assert isinstance(raw, dict)
    assert raw["path_basis"]["kind"] == "root-relative-v1"
    return raw


@pytest.mark.parametrize("absolute", [False, True])
def test_legacy_manifest_still_audits_but_saved_report_is_portable(
    tmp_path: Path, absolute: bool
) -> None:
    manifest = legacy_manifest(tmp_path, absolute=absolute)
    output = tmp_path / "audit" / "report.json"
    result = run_audit(manifest, output)
    assert result.returncode == 0, result.stderr
    report = assert_portable(output, tmp_path)
    assert report["passed"]
    assert report["manifest_sha256"] == file_sha256(manifest)
    pair = report["pairs"][0]
    assert pair["table_body_equal"]
    assert pair["current"]["rows"] == 2
    assert pair["current"]["sha256"] == file_sha256(tmp_path / "inputs/current.tsv")
    assert pair["inputs"]["profiles"]["sha256"] == file_sha256(tmp_path / "inputs/profiles.tsv")


def test_failed_audit_diagnostic_does_not_serialize_absolute_source_path(tmp_path: Path) -> None:
    manifest = legacy_manifest(tmp_path)
    current = tmp_path / "inputs/current.tsv"
    current.write_bytes(current.read_bytes().replace(b"\tinc\t", b"\tunknown\t"))
    output = tmp_path / "audit" / "report.json"
    result = run_audit(manifest, output)
    assert result.returncode == 1, result.stderr
    report = assert_portable(output, tmp_path)
    assert not report["passed"]
    assert "invalid inc_dec" in report["pairs"][0]["errors"][0]
    assert "inputs/current.tsv" in report["pairs"][0]["errors"][0]


def test_unknown_absolute_metadata_is_refused_without_redacting_source(tmp_path: Path) -> None:
    manifest = legacy_manifest(tmp_path)
    external = Path("/") / "Users" / "outside" / "untracked.tsv"
    before: dict[Path, bytes] = {}
    for name in ("historical.tsv", "current.tsv"):
        source = tmp_path / "inputs" / name
        source.write_bytes(f"# External-Note: {external}\n".encode() + source.read_bytes())
        before[source] = source.read_bytes()
    output = tmp_path / "audit.json"
    result = run_audit(manifest, output)
    assert result.returncode != 0
    assert not output.exists()
    assert {source: source.read_bytes() for source in before} == before


def test_portable_manifest_pins_sources_and_survives_relocation(tmp_path: Path) -> None:
    root = tmp_path / "original"
    manifest = legacy_manifest(root, absolute=True)
    portable = root / "records" / "inputs.json"
    _AUDIT.write_portable_manifest(
        manifest, portable, root=root, expected_audit=attest_sources(manifest)
    )
    declared = assert_portable(portable, root)
    pair = declared["pairs"][0]
    for key, filename in (
        ("historical", "historical.tsv"),
        ("current", "current.tsv"),
        ("profiles", "profiles.tsv"),
        ("taxa_matrix", "matrix.tsv"),
    ):
        assert pair[key] == {
            "path": f"inputs/{filename}",
            "sha256": file_sha256(root / "inputs" / filename),
        }
    moved = tmp_path / "relocated"
    shutil.copytree(root, moved)
    root.rename(tmp_path / "hidden_original")
    output = moved / "results" / "audit.json"
    result = run_audit(moved / "records/inputs.json", output)
    assert result.returncode == 0, result.stderr
    report = assert_portable(output, tmp_path)
    assert report["passed"]
    assert report["pairs"][0]["table_body_equal"]
    assert report["manifest_sha256"] == file_sha256(moved / "records/inputs.json")


@pytest.mark.parametrize("source", ["historical.tsv", "current.tsv", "profiles.tsv", "matrix.tsv"])
def test_portable_manifest_refuses_same_size_changed_source(tmp_path: Path, source: str) -> None:
    manifest = legacy_manifest(tmp_path)
    portable = tmp_path / "portable.json"
    _AUDIT.write_portable_manifest(manifest, portable, expected_audit=attest_sources(manifest))
    target = tmp_path / "inputs" / source
    old = target.read_bytes()
    # Same-size byte substitution: the manifest's hash, not a size-only check,
    # must reject it before parsing or comparing the ORA population.
    target.write_bytes(old.replace(b"TraitA", b"TraitB", 1))
    assert target.stat().st_size == len(old)
    result = run_audit(portable, tmp_path / "audit.json")
    assert result.returncode != 0
    diagnostic = result.stderr + result.stdout
    if (tmp_path / "audit.json").exists():
        diagnostic += (tmp_path / "audit.json").read_text(encoding="utf-8")
    assert "sha256" in diagnostic.lower().replace("-", "")


def test_portable_manifest_existing_output_refusal(tmp_path: Path) -> None:
    manifest = legacy_manifest(tmp_path)
    portable = tmp_path / "records" / "portable.json"
    _AUDIT.write_portable_manifest(manifest, portable, expected_audit=attest_sources(manifest))
    before = portable.read_bytes()
    with pytest.raises(FileExistsError):
        _AUDIT.write_portable_manifest(manifest, portable, expected_audit=attest_sources(manifest))
    assert portable.read_bytes() == before
    assert_portable(portable, tmp_path)


def test_portable_manifest_conversion_cli(tmp_path: Path) -> None:
    manifest = legacy_manifest(tmp_path, absolute=True)
    portable = tmp_path / "records" / "portable.json"
    result = subprocess.run(
        [
            sys.executable,
            str(_AUDIT.__file__),
            "--manifest",
            str(manifest),
            "--output",
            str(portable),
            "--write-portable-manifest",
            "--expected-audit",
            str(attest_sources(manifest)),
            "--provenance-root",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert_portable(portable, tmp_path)
    audited = run_audit(portable, tmp_path / "audit.json")
    assert audited.returncode == 0, audited.stderr


@pytest.mark.parametrize("reference", ["../inputs/current.tsv", "/absolute/current.tsv"])
def test_portable_manifest_rejects_noncanonical_source_reference(
    tmp_path: Path, reference: str
) -> None:
    manifest = legacy_manifest(tmp_path)
    portable = tmp_path / "portable.json"
    _AUDIT.write_portable_manifest(manifest, portable, expected_audit=attest_sources(manifest))
    raw = json.loads(portable.read_text(encoding="utf-8"))
    raw["pairs"][0]["current"]["path"] = reference
    portable.write_text(json.dumps(raw), encoding="utf-8")
    result = run_audit(portable, tmp_path / "audit.json")
    assert result.returncode != 0
    assert not (tmp_path / "audit.json").exists()


def test_portable_reference_refuses_same_bytes_reached_through_outside_symlink(
    tmp_path: Path,
) -> None:
    root = tmp_path / "bundle"
    manifest = legacy_manifest(root)
    portable = root / "portable.json"
    _AUDIT.write_portable_manifest(
        manifest, portable, root=root, expected_audit=attest_sources(manifest)
    )
    source = root / "inputs/current.tsv"
    external = tmp_path / "outside.tsv"
    external.write_bytes(source.read_bytes())
    source.unlink()
    source.symlink_to(external)
    result = run_audit(portable, root / "audit.json")
    assert result.returncode != 0
    assert not (root / "audit.json").exists()


def test_manifest_writer_refuses_sources_outside_explicit_root(tmp_path: Path) -> None:
    manifest = legacy_manifest(tmp_path)
    narrow = tmp_path / "records"
    narrow.mkdir()
    output = narrow / "portable.json"
    with pytest.raises(ValueError, match="outside the declared provenance root"):
        _AUDIT.write_portable_manifest(
            manifest, output, root=narrow, expected_audit=attest_sources(manifest)
        )
    assert not output.exists()


@pytest.mark.parametrize("source", ["historical.tsv", "current.tsv", "profiles.tsv", "matrix.tsv"])
def test_audit_retains_manifest_hash_after_reference_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source: str
) -> None:
    manifest = legacy_manifest(tmp_path)
    portable = tmp_path / "portable.json"
    _AUDIT.write_portable_manifest(manifest, portable, expected_audit=attest_sources(manifest))
    target = tmp_path / "inputs" / source
    before = target.read_bytes()
    changed = (
        before.replace(b"# Tail: upper\n", b"#Tail: upper \n", 1)
        if source in {"historical.tsv", "current.tsv"}
        else before.replace(b"\n", b"\r", 1)
    )
    assert changed != before
    assert len(changed) == len(before)
    verified = _AUDIT.verify_reference
    replacements = 0

    def verify_then_change(*args: Any, **kwargs: Any) -> Path:
        nonlocal replacements
        resolved = verified(*args, **kwargs)
        assert isinstance(resolved, Path)
        if resolved == target:
            target.write_bytes(changed)
            replacements += 1
        return resolved

    monkeypatch.setattr(_AUDIT, "verify_reference", verify_then_change)
    report = _AUDIT.audit_manifest(portable, tmp_path / "audit.json")
    assert replacements == 1
    assert not report["passed"]
    diagnostic = " ".join(report["pairs"][0]["errors"]).lower().replace("-", "")
    assert "manifest" in diagnostic
    assert "sha256" in diagnostic


def test_table_sha256_is_from_the_parsed_stream_not_a_second_read(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    legacy_manifest(tmp_path)
    target = tmp_path / "inputs/current.tsv"
    original = target.read_bytes()
    changed = original.replace(b"# Tail: upper\n", b"# Tail: lower\n", 1)
    assert original != changed
    assert len(original) == len(changed)
    population = _AUDIT.input_population(
        tmp_path / "inputs/profiles.tsv", tmp_path / "inputs/matrix.tsv"
    )
    original_open = Path.open
    reads = 0

    def alternating_reads(path: Path, mode: str = "r", *args: Any, **kwargs: Any) -> Any:
        nonlocal reads
        if path == target and mode == "rb":
            reads += 1
            return io.BytesIO(original if reads == 1 else changed)
        return original_open(path, mode, *args, **kwargs)

    monkeypatch.setattr(Path, "open", alternating_reads)
    table = _AUDIT.read_table(target, population)
    assert table.metadata["Tail"] == "upper"
    assert table.sha256 == hashlib.sha256(original).hexdigest()
    assert reads == 1
