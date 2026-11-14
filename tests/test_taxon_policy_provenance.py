"""Bind the contribution CLI's actual policy decision to its consumed bytes (#2307)."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

from mtp import taxon_contributions, taxon_influence
from mtp.trait_key_policy import TraitKeyPolicy
from test_contribution_trait_policy import _write_synthetic_contribution_inputs


def _run_cli(
    inputs: Path,
    protocol: Path,
    output: Path,
    monkeypatch: pytest.MonkeyPatch,
    explicit: TraitKeyPolicy | None = None,
) -> dict[str, object]:
    monkeypatch.setattr(taxon_influence, "run_influence", lambda *_args: {})
    arguments = [
        "ledger",
        "--kg-dir",
        str(inputs),
        "--inputs",
        str(inputs),
        "--protocol",
        str(protocol),
        "--output",
        str(output),
        "--strain-rollup",
        "phenotype",
    ]
    if explicit is not None:
        arguments.extend(["--trait-key-policy", explicit])
    monkeypatch.setattr(sys, "argv", arguments)
    try:
        taxon_contributions.main()
    except ValueError as exc:
        raise AssertionError(f"The captured policy must reconstruct the fixture: {exc}") from exc
    result = json.loads((output / "manifest.json").read_text())
    assert isinstance(result, dict)
    assert result["completed"] is True
    assert result["mismatches"] == 0
    return result


@pytest.mark.parametrize(
    ("policy", "explicit", "manifest", "columns", "basis"),
    [
        ("facultative_union_v1", False, False, False, "default"),
        ("raw_object_v1", True, False, False, "explicit"),
        ("raw_object_v1", False, False, True, "column_metadata"),
        ("raw_object_v1", False, True, False, "core_manifest"),
        ("raw_object_v1", False, True, True, "core_manifest"),
        ("raw_object_v1", True, True, True, "explicit"),
    ],
)
def test_taxon_manifest_records_policy_basis_and_consumed_sources(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    policy: TraitKeyPolicy,
    explicit: bool,
    manifest: bool,
    columns: bool,
    basis: str,
) -> None:
    inputs = tmp_path / "inputs"
    protocol = _write_synthetic_contribution_inputs(
        inputs,
        policy,
        with_run_manifest=manifest,
        with_column_metadata=columns,
    )
    sources = {}
    if manifest:
        sources["core_manifest"] = hashlib.sha256((inputs / "run.json").read_bytes()).hexdigest()
    if columns:
        sources["taxa_trait_matrix_columns.tsv"] = hashlib.sha256(
            (inputs / "taxa_trait_matrix_columns.tsv").read_bytes()
        ).hexdigest()

    result = _run_cli(
        inputs, protocol, tmp_path / "output", monkeypatch, policy if explicit else None
    )

    assert result["trait_key_policy"] == policy
    assert result.get("trait_key_policy_resolution") == {
        "policy": policy,
        "basis": basis,
        "source_sha256": sources,
    }


def test_taxon_manifest_hashes_the_policy_buffers_used_for_its_counts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    inputs = tmp_path / "inputs"
    protocol = _write_synthetic_contribution_inputs(
        inputs,
        "raw_object_v1",
        with_run_manifest=True,
        with_column_metadata=True,
    )
    source_paths = {
        "core_manifest": inputs / "run.json",
        "taxa_trait_matrix_columns.tsv": inputs / "taxa_trait_matrix_columns.tsv",
    }
    captured = {path: path.read_bytes() for path in source_paths.values()}
    calls = dict.fromkeys(captured, 0)
    original = Path.read_bytes

    def changed_after_capture(path: Path) -> bytes:
        value = original(path)
        if path in captured:
            calls[path] += 1
            changed = value.replace(b"raw_object_v1", b"facultative_union_v1")
            # A historical manifest can encode policy through producer identity,
            # without containing the policy name. Its bytes must still change.
            path.write_bytes(changed if changed != value else value + b"\n")
        return value

    monkeypatch.setattr(Path, "read_bytes", changed_after_capture)
    result = _run_cli(inputs, protocol, tmp_path / "output", monkeypatch)

    assert result["trait_key_policy"] == "raw_object_v1"
    assert result.get("trait_key_policy_resolution") == {
        "policy": "raw_object_v1",
        "basis": "core_manifest",
        "source_sha256": {
            role: hashlib.sha256(captured[path]).hexdigest() for role, path in source_paths.items()
        },
    }
    assert set(calls.values()) == {1}
    assert all(original(path) != value for path, value in captured.items())


def test_taxon_cli_describes_policy_metadata_and_the_callers_verification_duty() -> None:
    help_text = " ".join(taxon_contributions.argument_parser().format_help().split())
    assert "verify and inherit" not in help_text
    assert "policy metadata" in help_text
    assert "caller" in help_text.lower()
    assert "bind" in help_text
