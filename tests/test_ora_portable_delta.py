"""ORA attestation boundaries and actionable source diagnostics (#1767/#1768)."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from mtp.paths import file_sha256
from mtp.provenance_paths import assert_portable_text
from test_ora_portable_provenance import assigned_source as assigned_source
from test_ora_portable_provenance import attest_sources, legacy_manifest
from test_ora_refresh_audit import _AUDIT

ROLES = ("historical", "current", "profiles", "taxa_matrix")


def refusal(manifest: Path, output: Path, *, convert: bool = False, **kwargs: Any) -> str:
    try:
        if convert:
            _AUDIT.write_portable_manifest(manifest, output, **kwargs)
        else:
            _AUDIT.audit_manifest(manifest, output, **kwargs)
    except (OSError, ValueError) as exc:
        return str(exc)
    return ""


def assert_no_location(message: str, root: Path) -> None:
    assert message
    assert str(root.resolve()) not in message
    assert_portable_text(message)


def test_attesting_audit_rejects_a_failed_pair_despite_overall_success(tmp_path: Path) -> None:
    manifest = legacy_manifest(tmp_path)
    audit = attest_sources(manifest)
    record = json.loads(audit.read_text())
    assert record["passed"] is True
    record["pairs"][0]["passed"] = False
    audit.write_text(json.dumps(record))
    output = tmp_path / "converted.json"
    diagnostic = refusal(manifest, output, convert=True, expected_audit=audit)
    assert "failed" in diagnostic.lower()
    assert not output.exists()


def test_expected_audit_cli_flag_requires_conversion(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    manifest = legacy_manifest(tmp_path)
    audit = attest_sources(manifest)
    output = tmp_path / "unused.json"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "ora-audit",
            "--manifest",
            str(manifest),
            "--output",
            str(output),
            "--expected-audit",
            str(audit),
        ],
    )
    with pytest.raises(SystemExit) as caught:
        _AUDIT.main()
    assert caught.value.code == 2
    assert "--expected-audit is only used with --write-portable-manifest" in capsys.readouterr().err
    assert not output.exists()


@pytest.mark.parametrize("convert", [False, True])
def test_manifest_outside_root_is_named_by_role(tmp_path: Path, convert: bool) -> None:
    manifest = legacy_manifest(tmp_path, absolute=True)
    audit = attest_sources(manifest)
    output = tmp_path / "report.json"
    kwargs: dict[str, Any] = {"root": tmp_path / "inputs"}
    if convert:
        kwargs["expected_audit"] = audit
    diagnostic = refusal(manifest, output, convert=convert, **kwargs)
    assert "manifest" in diagnostic.lower()
    assert "outside" in diagnostic.lower()
    assert_no_location(diagnostic, tmp_path)
    assert not output.exists()


@pytest.mark.parametrize("role", ROLES)
def test_source_outside_root_names_pair_and_role(tmp_path: Path, role: str) -> None:
    root = tmp_path / "source-root"
    manifest = legacy_manifest(root, absolute=True)
    record = json.loads(manifest.read_text())
    pair = record["pairs"][0]
    outside = tmp_path / "outside-source.tsv"
    shutil.copy2(pair[role], outside)
    pair[role] = str(outside)
    manifest.write_text(json.dumps(record))
    output = tmp_path / "report.json"
    diagnostic = refusal(manifest, output, root=root)
    assert pair["label"] in diagnostic
    assert role in diagnostic
    assert "outside" in diagnostic.lower()
    assert_no_location(diagnostic, tmp_path)
    assert not output.exists()


@pytest.mark.parametrize("convert", [False, True])
def test_absent_sources_are_counted_before_root_validation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, convert: bool
) -> None:
    pairs = []
    for number in range(5):
        source_manifest = legacy_manifest(tmp_path / f"unit-{number}", absolute=True)
        pair = json.loads(source_manifest.read_text())["pairs"][0]
        pair["label"] = f"unit-{number}"
        pairs.append(pair)
    manifest = tmp_path / "all-inputs.json"
    manifest.write_text(json.dumps({"pairs": pairs}))
    audit = attest_sources(manifest)
    for pair in pairs[:3]:
        for role in ROLES:
            Path(pair[role]).unlink()

    def root_validation_must_not_run(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("Source presence must be checked before provenance-root validation")

    monkeypatch.setattr(_AUDIT.PortablePaths, "create", root_validation_must_not_run)
    output = tmp_path / "report.json"
    kwargs: dict[str, Any] = {"root": tmp_path / "unrelated-missing-root"}
    if convert:
        kwargs["expected_audit"] = audit
    diagnostic = refusal(manifest, output, convert=convert, **kwargs)
    assert "12 of 20" in diagnostic
    assert "missing" in diagnostic.lower()
    for pair in pairs[:3]:
        for role in ROLES:
            assert f"{pair['label']}:{role}" in diagnostic
    assert_no_location(diagnostic, tmp_path)
    assert not output.exists()


@pytest.mark.parametrize("role", ROLES)
def test_missing_portable_source_names_pair_and_role(tmp_path: Path, role: str) -> None:
    manifest = legacy_manifest(tmp_path)
    portable = tmp_path / "portable.json"
    _AUDIT.write_portable_manifest(manifest, portable, expected_audit=attest_sources(manifest))
    record = json.loads(portable.read_text())
    pair = record["pairs"][0]
    (tmp_path / pair[role]["path"]).unlink()
    output = tmp_path / "report.json"
    diagnostic = refusal(portable, output)
    assert pair["label"] in diagnostic
    assert role in diagnostic
    assert "missing" in diagnostic.lower()
    assert_no_location(diagnostic, tmp_path)
    assert not output.exists()


@pytest.mark.parametrize("convert", [False, True])
def test_missing_manifest_is_named_without_disclosing_location(
    tmp_path: Path, convert: bool
) -> None:
    manifest = tmp_path / "missing-manifest.json"
    output = tmp_path / "report.json"
    diagnostic = refusal(manifest, output, convert=convert)
    assert "manifest" in diagnostic.lower()
    assert "missing" in diagnostic.lower()
    assert_no_location(diagnostic, tmp_path)
    assert not output.exists()


@pytest.mark.parametrize("role", ROLES)
def test_portable_digest_refusal_names_pair_and_role(tmp_path: Path, role: str) -> None:
    manifest = legacy_manifest(tmp_path)
    portable = tmp_path / "portable.json"
    _AUDIT.write_portable_manifest(manifest, portable, expected_audit=attest_sources(manifest))
    record = json.loads(portable.read_text())
    pair = record["pairs"][0]
    del pair[role]["sha256"]
    portable.write_text(json.dumps(record))
    output = tmp_path / "report.json"
    diagnostic = refusal(portable, output)
    assert pair["label"] in diagnostic
    assert role in diagnostic
    assert "sha" in diagnostic.lower() or "digest" in diagnostic.lower()
    assert_no_location(diagnostic, tmp_path)
    assert not output.exists()


def test_actual_cli_refuses_temp_boundary_but_accepts_deeper_source_root(tmp_path: Path) -> None:
    account = "syntheticaccount"
    boundary = tmp_path / "temporary-boundary"
    source_root = boundary / f"-Users-{account}-project"
    manifest = legacy_manifest(source_root)
    reports = tmp_path / "external-reports"
    env = {
        **os.environ,
        "USER": account,
        "TMPDIR": str(boundary),
        "PYTHONDONTWRITEBYTECODE": "1",
    }

    def invoke(root: Path, output: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(_AUDIT.__file__),
                "--manifest",
                str(manifest),
                "--output",
                str(output),
                "--provenance-root",
                str(root),
            ],
            env=env,
            check=False,
            capture_output=True,
            text=True,
        )

    refused = reports / "refused.json"
    broad = invoke(boundary, refused)
    assert broad.returncode == 2, broad.stderr
    assert "--provenance-root" in broad.stderr
    assert not refused.exists()

    accepted = reports / "accepted.json"
    deeper = invoke(source_root, accepted)
    assert deeper.returncode == 0, deeper.stderr
    assert accepted.is_file()
    assert accepted.stat().st_size > 0
    text = accepted.read_text()
    record = json.loads(text)
    assert record["passed"] is True
    assert record["path_basis"] == {"kind": "supplied-root-v1"}
    assert account.casefold() not in text.casefold()
    assert str(source_root) not in text
    assert_portable_text(text)


@pytest.mark.parametrize("account", ["root", "ci", "app"])
def test_actual_cli_accepts_root_and_short_accounts_in_nonpath_fields(
    tmp_path: Path, account: str
) -> None:
    """Both writers accept ordinary record text under container-style identities (#1773)."""
    root = tmp_path / "bundle"
    manifest = legacy_manifest(root)
    original = json.loads(manifest.read_text())
    original["pairs"][0]["label"] = account
    manifest.write_text(json.dumps(original))
    home = tmp_path / "runtime-home"
    temporary = tmp_path / "runtime-temp"
    home.mkdir()
    temporary.mkdir()
    env = {
        **os.environ,
        **dict.fromkeys(("USER", "LOGNAME", "LNAME", "USERNAME"), account),
        **{key: str(temporary) for key in ("TMPDIR", "TEMP", "TMP")},
        "HOME": str(home),
        "PYTHONPATH": str(Path(__file__).resolve().parents[1] / "src"),
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    audit = root / "audit" / "report.json"
    converted = root / "records" / "portable.json"
    for output, options in (
        (audit, []),
        (converted, ["--write-portable-manifest", "--expected-audit", str(audit)]),
    ):
        result = subprocess.run(
            [
                sys.executable,
                "-B",
                str(_AUDIT.__file__),
                "--manifest",
                str(manifest),
                "--output",
                str(output),
                *options,
            ],
            env=env,
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert output.is_file()
        assert output.stat().st_size > 0
        text = output.read_text()
        record = json.loads(text)
        assert record["path_basis"] == {"kind": "root-relative-v1", "root": ".."}
        assert record["pairs"][0]["label"] == account
        assert str(root) not in text
        assert_portable_text(text)
    report = json.loads(audit.read_text())
    assert report["passed"] is True
    assert report["pairs"][0]["current"]["sha256"] == file_sha256(root / "inputs/current.tsv")
    portable = json.loads(converted.read_text())
    assert portable["converted_from"] == {
        "manifest_sha256": file_sha256(manifest),
        "audit_sha256": file_sha256(audit),
    }
    for role, filename in (
        ("historical", "historical.tsv"),
        ("current", "current.tsv"),
        ("profiles", "profiles.tsv"),
        ("taxa_matrix", "matrix.tsv"),
    ):
        assert portable["pairs"][0][role] == {
            "path": "inputs/" + filename,
            "sha256": file_sha256(root / "inputs" / filename),
        }
