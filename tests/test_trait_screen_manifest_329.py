"""Regression tests for trait_screen_manifest.json self-attestation disclosure (#329)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[1]
REQUIRED_AMENDMENT_KEYS: set[str] = {
    "date",
    "change",
    "reason",
    "known_when_amended",
    "unchanged",
}


def validate_trait_screen_manifest_329(
    manifest_path: Path = ROOT / "methods" / "trait_screen_manifest.json",
) -> None:
    """Verify trait_screen_manifest.json preserves keys and discloses #329 self-attestation."""
    data: dict[str, Any] = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert isinstance(data, dict), "Expected top-level JSON object in trait_screen_manifest.json"

    assert data.get("protocol") == "trait-direction-screen-v1", (
        f"Expected protocol == 'trait-direction-screen-v1', got {data.get('protocol')!r}"
    )
    assert data.get("comparison_spec") == "three-framework-comparison-v1", (
        "Expected comparison_spec == 'three-framework-comparison-v1', "
        f"got {data.get('comparison_spec')!r}"
    )
    assert data.get("declared_before_first_run") is True, (
        "Expected declared_before_first_run to remain boolean True, "
        f"got {data.get('declared_before_first_run')!r}"
    )
    assert data.get("candidate_count_on_the_pinned_corpus") == 317, (
        "Expected candidate_count_on_the_pinned_corpus == 317, "
        f"got {data.get('candidate_count_on_the_pinned_corpus')!r}"
    )

    declared_by = data.get("declared_by")
    assert isinstance(declared_by, str), "Expected declared_by to be a string"
    assert "self-attestation" in declared_by, (
        "Expected declared_by to state that pre-run declaration rests on self-attestation (#329)"
    )
    assert "1f3d36c" in declared_by, (
        "Expected declared_by to cite commit 1f3d36c where manifest and methods doc landed (#329)"
    )

    verification = data.get("declared_before_first_run_verification")
    assert isinstance(verification, str), (
        "Expected declared_before_first_run_verification string in manifest (#329)"
    )
    assert "self-attestation" in verification, (
        "Expected declared_before_first_run_verification to name self-attestation (#329)"
    )
    assert "cannot be verified from Git history" in verification, (
        "Expected declared_before_first_run_verification to state that timing "
        "cannot be verified from Git history (#329)"
    )
    assert "1f3d36c" in verification, (
        "Expected declared_before_first_run_verification to cite commit 1f3d36c (#329)"
    )
    assert "Li & Ji" in verification, (
        "Expected declared_before_first_run_verification to cite Li & Ji (#329)"
    )
    assert "112.0" in verification, (
        "Expected declared_before_first_run_verification to cite effective family size 112.0 (#329)"
    )
    assert "committed on their own before" in verification, (
        "Expected declared_before_first_run_verification to require future pre-registration "
        "manifests to be committed on their own before results (#329)"
    )

    amendments = data.get("amendments")
    assert isinstance(amendments, list), "Expected amendments list in trait_screen_manifest.json"
    assert len(amendments) >= 3, (
        f"Expected at least 3 amendments (including #329), got {len(amendments)}"
    )
    for entry in amendments:
        assert isinstance(entry, dict), f"Expected dict amendment entry, got {entry!r}"
        assert set(entry) >= REQUIRED_AMENDMENT_KEYS, (
            f"Amendment missing required keys: {REQUIRED_AMENDMENT_KEYS - set(entry)}"
        )

    dates = [str(entry["date"]) for entry in amendments]
    assert dates == sorted(dates), f"Expected amendments in chronological date order, got {dates!r}"

    matching_329 = [
        entry
        for entry in amendments
        if "#329" in str(entry.get("reason", ""))
        and "1f3d36c" in str(entry.get("reason", ""))
        and "112.0" in str(entry.get("reason", ""))
        and "self-attestation" in str(entry.get("change", ""))
    ]
    assert len(matching_329) == 1, (
        "Expected exactly one #329 amendment entry naming 1f3d36c, 112.0, and self-attestation"
    )


def test_repo_trait_screen_manifest_discloses_self_attestation_329() -> None:
    """Repository methods/trait_screen_manifest.json must disclose #329 self-attestation."""
    validate_trait_screen_manifest_329()


def test_validate_trait_screen_manifest_329_rejects_missing_disclosures(
    tmp_path: Path,
) -> None:
    """Validator must fail closed when self-attestation, 1f3d36c, 112.0, or #329 is absent."""
    manifest_path = tmp_path / "trait_screen_manifest.json"
    base_manifest: dict[str, Any] = {
        "protocol": "trait-direction-screen-v1",
        "comparison_spec": "three-framework-comparison-v1",
        "declared_before_first_run": True,
        "candidate_count_on_the_pinned_corpus": 317,
        "declared_by": (
            "drafted by a coding agent before results (self-attestation only: "
            "cannot be verified from Git history because it landed in 1f3d36c; #329)"
        ),
        "declared_before_first_run_verification": (
            "self-attestation only (#329): cannot be verified from Git history alone "
            "because commit 1f3d36c already quoted Li & Ji 112.0; future pre-registration "
            "manifests must be committed on their own before results."
        ),
        "amendments": [
            {
                "date": "2026-10-08",
                "change": "c1",
                "reason": "r1",
                "known_when_amended": "k1",
                "unchanged": "u1",
            },
            {
                "date": "2026-10-08",
                "change": "c2",
                "reason": "r2",
                "known_when_amended": "k2",
                "unchanged": "u2",
            },
            {
                "date": "2026-11-09",
                "change": "disclose self-attestation for declared_before_first_run",
                "reason": "#329: commit 1f3d36c already quoted Li & Ji 112.0",
                "known_when_amended": "k3",
                "unchanged": "u3",
            },
        ],
    }
    manifest_path.write_text(json.dumps(base_manifest), encoding="utf-8")
    validate_trait_screen_manifest_329(manifest_path)

    bad = dict(base_manifest, declared_before_first_run=False)
    manifest_path.write_text(json.dumps(bad), encoding="utf-8")
    with pytest.raises(AssertionError, match="declared_before_first_run"):
        validate_trait_screen_manifest_329(manifest_path)

    bad2 = dict(base_manifest)
    bad2.pop("declared_before_first_run_verification")
    manifest_path.write_text(json.dumps(bad2), encoding="utf-8")
    with pytest.raises(AssertionError, match="declared_before_first_run_verification"):
        validate_trait_screen_manifest_329(manifest_path)
