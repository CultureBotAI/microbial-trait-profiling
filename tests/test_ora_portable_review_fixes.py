"""Historical attestations and hash-check boundaries for PR #1749 (#1753/#1754/#1756)."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from mtp.paths import file_sha256
from test_no_absolute_local_paths import ABSOLUTE_LOCAL_PATH_PATTERN
from test_ora_portable_provenance import assigned_source as assigned_source
from test_ora_portable_provenance import attest_sources, legacy_manifest
from test_ora_refresh_audit import _AUDIT


def conversion_refusal(manifest: Path, output: Path, **kwargs: Any) -> str:
    try:
        _AUDIT.write_portable_manifest(manifest, output, **kwargs)
    except ValueError as exc:
        return str(exc)
    return ""


def test_legacy_conversion_requires_an_attesting_audit(tmp_path: Path) -> None:
    manifest = legacy_manifest(tmp_path)
    output = tmp_path / "portable.json"
    diagnostic = conversion_refusal(manifest, output)
    assert "expected-audit" in diagnostic
    assert not output.exists()


def test_conversion_binds_exact_manifest_and_audit_buffers(tmp_path: Path) -> None:
    manifest = legacy_manifest(tmp_path)
    audit = attest_sources(manifest)
    output = tmp_path / "portable.json"
    _AUDIT.write_portable_manifest(manifest, output, expected_audit=audit)
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["converted_from"] == {
        "manifest_sha256": file_sha256(manifest),
        "audit_sha256": file_sha256(audit),
    }
    verified = _AUDIT.audit_manifest(output, tmp_path / "new-audit.json")
    assert verified["passed"]


@pytest.mark.parametrize("source", ["historical.tsv", "current.tsv", "profiles.tsv", "matrix.tsv"])
def test_legacy_conversion_refuses_changed_attested_bytes(tmp_path: Path, source: str) -> None:
    manifest = legacy_manifest(tmp_path)
    audit = attest_sources(manifest)
    target = tmp_path / "inputs" / source
    before = target.read_bytes()
    target.write_bytes(before.replace(b"TraitA", b"TraitB", 1))
    assert target.stat().st_size == len(before)
    output = tmp_path / "portable.json"
    diagnostic = conversion_refusal(manifest, output, expected_audit=audit)
    assert "sha256" in diagnostic.lower().replace("-", "")
    assert not output.exists()


@pytest.mark.parametrize("mutation", ["manifest", "label", "duplicate", "missing_hash", "failed"])
def test_conversion_rejects_inapplicable_attestation(tmp_path: Path, mutation: str) -> None:
    manifest = legacy_manifest(tmp_path)
    audit = attest_sources(manifest)
    raw = json.loads(audit.read_text(encoding="utf-8"))
    if mutation == "manifest":
        raw["manifest_sha256"] = "0" * 64
    elif mutation == "label":
        raw["pairs"][0]["label"] = "other-pair"
    elif mutation == "duplicate":
        raw["pairs"].append(raw["pairs"][0])
    elif mutation == "missing_hash":
        del raw["pairs"][0]["inputs"]["profiles"]["sha256"]
    else:
        raw["passed"] = False
    audit.write_text(json.dumps(raw), encoding="utf-8")
    output = tmp_path / "portable.json"
    diagnostic = conversion_refusal(manifest, output, expected_audit=audit)
    assert diagnostic
    assert not output.exists()


def test_conversion_hashes_the_attestation_buffer_it_parsed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest = legacy_manifest(tmp_path)
    audit = attest_sources(manifest)
    original = audit.read_bytes()
    changed = original.replace(b'"passed": true', b'"passed":false', 1)
    assert original != changed
    read_bytes = Path.read_bytes
    reads = 0

    def alternate(path: Path) -> bytes:
        nonlocal reads
        if path == audit:
            reads += 1
            if reads == 1:
                # A later Path.open/hash read must also see the changed bytes,
                # not just a second call to this patched read_bytes method.
                audit.write_bytes(changed)
            return original if reads == 1 else changed
        return read_bytes(path)

    monkeypatch.setattr(Path, "read_bytes", alternate)
    output = tmp_path / "portable.json"
    _AUDIT.write_portable_manifest(manifest, output, expected_audit=audit)
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["converted_from"]["audit_sha256"] == hashlib.sha256(original).hexdigest()
    assert reads == 1


def test_portable_reconversion_rechecks_declared_hashes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest = legacy_manifest(tmp_path)
    portable = tmp_path / "portable.json"
    _AUDIT.write_portable_manifest(manifest, portable, expected_audit=attest_sources(manifest))
    target = tmp_path / "inputs/historical.tsv"
    original = target.read_bytes()
    changed = original.replace(b"# Tail: upper\n", b"#Tail: upper \n", 1)
    assert original != changed
    verified = _AUDIT.verify_reference

    def verify_then_change(*args: Any, **kwargs: Any) -> Path:
        path = verified(*args, **kwargs)
        assert isinstance(path, Path)
        if path == target:
            target.write_bytes(changed)
        return path

    monkeypatch.setattr(_AUDIT, "verify_reference", verify_then_change)
    output = tmp_path / "converted-again.json"
    diagnostic = conversion_refusal(portable, output)
    assert "sha256" in diagnostic.lower().replace("-", "")
    assert not output.exists()


def test_expected_hash_mismatch_refuses_before_parsing_sources(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest = legacy_manifest(tmp_path)
    pair = json.loads(manifest.read_text(encoding="utf-8"))["pairs"][0]
    expected = {key: file_sha256(tmp_path / value) for key, value in pair.items() if key != "label"}
    expected["historical"] = "0" * 64
    population = _AUDIT.input_population
    calls = 0

    def counted_population(*args: Any, **kwargs: Any) -> Any:
        nonlocal calls
        calls += 1
        return population(*args, **kwargs)

    monkeypatch.setattr(_AUDIT, "input_population", counted_population)
    diagnostic = ""
    try:
        _AUDIT.audit_pair(pair, tmp_path, expected_hashes=expected)
    except ValueError as exc:
        diagnostic = str(exc)
    assert "before auditing" in diagnostic
    assert calls == 0


def test_legacy_audit_refuses_table_changed_after_before_hash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest = legacy_manifest(tmp_path)
    target = tmp_path / "inputs/historical.tsv"
    original = target.read_bytes()
    changed = original.replace(b"# Tail: upper\n", b"#Tail: upper \n", 1)
    assert original != changed
    population = _AUDIT.input_population

    def change_after_before_hash(*args: Any, **kwargs: Any) -> Any:
        target.write_bytes(changed)
        return population(*args, **kwargs)

    monkeypatch.setattr(_AUDIT, "input_population", change_after_before_hash)
    result = _AUDIT.audit_manifest(manifest, tmp_path / "audit.json")
    assert not result["passed"]
    assert "changed while auditing" in " ".join(result["pairs"][0]["errors"])


def test_consumed_hash_must_match_still_current_caller_declaration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """M19 is redundant for immutable identities; pin its mutable public-API boundary."""
    manifest = legacy_manifest(tmp_path)
    pair = json.loads(manifest.read_text(encoding="utf-8"))["pairs"][0]
    expected = {key: file_sha256(tmp_path / value) for key, value in pair.items() if key != "label"}
    population = _AUDIT.input_population

    def change_declaration(*args: Any, **kwargs: Any) -> Any:
        expected["historical"] = "0" * 64
        return population(*args, **kwargs)

    monkeypatch.setattr(_AUDIT, "input_population", change_declaration)
    diagnostic = ""
    try:
        _AUDIT.audit_pair(pair, tmp_path, expected_hashes=expected)
    except ValueError as exc:
        diagnostic = str(exc)
    assert "Consumed source SHA-256 differs from manifest" in diagnostic


def test_audit_cli_accepts_explicit_provenance_root(tmp_path: Path) -> None:
    manifest = legacy_manifest(tmp_path)
    output = tmp_path / "record" / "audit.json"
    result = subprocess.run(
        [
            sys.executable,
            str(_AUDIT.__file__),
            "--manifest",
            str(manifest),
            "--output",
            str(output),
            "--provenance-root",
            str(tmp_path),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    record = json.loads(output.read_text(encoding="utf-8"))
    assert record["passed"]
    assert record["path_basis"] == {"kind": "root-relative-v1", "root": ".."}
    assert record["pairs"][0]["historical"]["path"] == "inputs/historical.tsv"


def test_audit_cli_external_output_uses_supplied_source_root_without_machine_prefixes(
    tmp_path: Path,
) -> None:
    source_root = tmp_path / "private" / "tmp" / "source-fixture"
    manifest = legacy_manifest(source_root)
    output = tmp_path / "Users" / "synthetic-owner" / "report.json"
    result = subprocess.run(
        [
            sys.executable,
            str(_AUDIT.__file__),
            "--manifest",
            str(manifest),
            "--output",
            str(output),
            "--provenance-root",
            str(source_root),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    text = output.read_text(encoding="utf-8")
    assert not ABSOLUTE_LOCAL_PATH_PATTERN.search(text)
    assert str(source_root) not in text
    assert str(source_root).lstrip("/") not in text
    record = json.loads(text)
    assert record["path_basis"] == {"kind": "supplied-root-v1"}
    assert record["manifest"] == "legacy.json"
    for key in ("historical", "current"):
        assert not ABSOLUTE_LOCAL_PATH_PATTERN.search("/" + record["pairs"][0][key]["path"])
    assert record["passed"]
