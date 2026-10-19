"""Reject incomplete or scientifically changed ORA refreshes (#720)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from mtp.paths import file_sha256


@pytest.fixture
def refresh_manifest(tmp_path: Path) -> Path:
    """Use actual producer output, including mixed metadata/table line endings."""
    profiles = tmp_path / "profiles.tsv"
    matrix = tmp_path / "matrix.tsv"
    profiles.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tTraitA\n"
        "s1\te1\tinc\t2\t2\n"
        "s1\te1\tdec\t2\t0\n",
        encoding="utf-8",
    )
    matrix.write_text("taxon_id\tTraitA\n1\t1\n2\t1\n3\t0\n4\t0\n", encoding="utf-8")
    current = tmp_path / "current.tsv"
    subprocess.run(
        [
            sys.executable,
            "scripts/trait_enrichment.py",
            "--input",
            str(profiles),
            "--taxa-matrix",
            str(matrix),
            "--output",
            str(tmp_path / "report.txt"),
            "--all-results-output",
            str(current),
            "--p-threshold",
            "1.0",
            "--min-fold",
            "0.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    (tmp_path / "historical.tsv").write_bytes(
        b"".join(
            line
            for line in current.read_bytes().splitlines(keepends=True)
            if not line.startswith((b"# Profiles-SHA256:", b"# Taxa-Matrix-SHA256:"))
        )
    )
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "pairs": [
                    {
                        "label": "fixture-upper",
                        "historical": "historical.tsv",
                        "current": "current.tsv",
                        "profiles": "profiles.tsv",
                        "taxa_matrix": "matrix.tsv",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    return manifest


def run_audit(manifest: Path, output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "scripts/adhoc/issue720_audit_ora_refresh.py",
            "--manifest",
            str(manifest),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        check=False,
    )


def test_refresh_audit_accepts_only_new_fingerprints(refresh_manifest: Path) -> None:
    output = refresh_manifest.parent / "audit.json"
    historical = refresh_manifest.parent / "historical.tsv"
    before = file_sha256(historical)
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 0, result.stderr
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["passed"]
    pair = report["pairs"][0]
    assert pair["table_body_equal"]
    assert pair["current"]["rows"] == 2
    assert pair["current"]["significant_rows"] == 2
    assert set(pair["metadata_changes"]) == {"Profiles-SHA256", "Taxa-Matrix-SHA256"}
    assert file_sha256(historical) == before
    assert pair["current"]["body_sha256"] == pair["historical"]["body_sha256"]


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("body", "Table body bytes differ"),
        ("truncated", "Total-Tests-Performed: declared 2, observed 1"),
        ("selected_count", "Significant-Hits: declared 3, observed 2"),
        ("metadata", "Unexpected metadata change: FDR-Scope"),
        ("profiles", "Profiles-SHA256: current export does not match supplied input"),
        ("matrix", "Taxa-Matrix-SHA256: current export does not match supplied input"),
        ("missing_fingerprint", "Profiles-SHA256: current export does not match supplied input"),
    ],
)
def test_refresh_audit_rejects_changed_or_incomplete_results(
    refresh_manifest: Path, mutation: str, message: str
) -> None:
    current = refresh_manifest.parent / "current.tsv"
    contents = current.read_bytes()
    if mutation == "body":
        current.write_bytes(contents.replace(b"s1\te1\tinc", b"s2\te1\tinc"))
    elif mutation == "truncated":
        current.write_bytes(b"".join(contents.splitlines(keepends=True)[:-1]))
    elif mutation == "selected_count":
        current.write_bytes(contents.replace(b"# Significant-Hits: 2", b"# Significant-Hits: 3"))
    elif mutation == "metadata":
        current.write_bytes(contents.replace(b"# FDR-Scope: global", b"# FDR-Scope: per-profile"))
    elif mutation in {"profiles", "matrix"}:
        path = refresh_manifest.parent / f"{mutation}.tsv"
        path.write_bytes(path.read_bytes() + b"\n")
    elif mutation == "missing_fingerprint":
        current.write_bytes(
            b"".join(
                line
                for line in contents.splitlines(keepends=True)
                if not line.startswith(b"# Profiles-SHA256:")
            )
        )
    output = refresh_manifest.parent / "audit.json"
    result = run_audit(refresh_manifest, output)
    assert result.returncode == 1, result.stderr
    report = json.loads(output.read_text(encoding="utf-8"))
    assert not report["passed"]
    assert any(message in error for error in report["pairs"][0]["errors"])


@pytest.mark.parametrize("hard_link", [False, True])
def test_refresh_audit_cannot_overwrite_historical_input(
    refresh_manifest: Path, hard_link: bool
) -> None:
    historical = refresh_manifest.parent / "historical.tsv"
    output = historical
    if hard_link:
        output = refresh_manifest.parent / "audit.json"
        output.hardlink_to(historical)
    before = file_sha256(historical)
    result = run_audit(refresh_manifest, output)
    assert result.returncode != 0
    assert "must not overwrite any input" in result.stderr
    assert file_sha256(historical) == before
