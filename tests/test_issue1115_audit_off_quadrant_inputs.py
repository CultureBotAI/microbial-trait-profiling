"""Refuse audit report destinations inside completed input bundles (#1186)."""

from __future__ import annotations

import hashlib
import json
import runpy
import sys
from pathlib import Path
from typing import Any

import pytest

SCRIPT_PATH = Path(__file__).parents[1] / "scripts/adhoc/issue1115_audit_off_quadrant_inputs.py"


def _load_module() -> dict[str, Any]:
    return runpy.run_path(str(SCRIPT_PATH))


def _identity(raw: bytes) -> dict[str, str | int]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def _tree_state(root: Path) -> tuple[list[str], dict[str, bytes]]:
    entries = sorted(path.relative_to(root).as_posix() for path in root.rglob("*"))
    files = {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }
    return entries, files


@pytest.fixture
def audit_fixture(tmp_path: Path) -> dict[str, Any]:
    parent_bundle = tmp_path / "parent_bundle"
    parent_bundle.mkdir()
    pair_bytes = (
        b"study_id\texperiment_id\tMETPO:1000608\tMETPO:1000607\nstudy_1\texp_1\t0.25\t-0.25\n"
    )
    study_bytes = b"study_id\tMETPO:1000608\tMETPO:1000607\nstudy_1\t0.25\t-0.25\n"
    evidence_bytes = (
        b"trait\ttrait_label\n"
        b"METPO:1000608\tfacultatively aerobic\n"
        b"METPO:1000607\tobligately anaerobic\n"
    )
    (parent_bundle / "composition_pair_delta_z.tsv").write_bytes(pair_bytes)
    (parent_bundle / "composition_study_delta_z.tsv").write_bytes(study_bytes)
    (parent_bundle / "per_trait_evidence.tsv").write_bytes(evidence_bytes)
    parent_record = {
        "run_id": "three-test-off-fixture",
        "code_commit": "a" * 40,
        "strain_rollup": "off",
        "contrast_population": "all_retained",
        "direction_basis": "case_minus_control",
        "retained_traits": ["METPO:1000608", "METPO:1000607"],
        "outputs": {
            "composition_pair_delta_z.tsv": _identity(pair_bytes),
            "composition_study_delta_z.tsv": _identity(study_bytes),
            "per_trait_evidence.tsv": _identity(evidence_bytes),
        },
    }
    parent_manifest = parent_bundle / "three_test_run.json"
    parent_bytes = (json.dumps(parent_record, indent=2) + "\n").encode("utf-8")
    parent_manifest.write_bytes(parent_bytes)

    reference_bundle = tmp_path / "reference_bundle"
    reference_bundle.mkdir()
    reference_summary = b"trait\tsdz_active_at_least_50\nMETPO:1000608\tTrue\n"
    (reference_bundle / "trait_summary.tsv").write_bytes(reference_summary)
    reference_manifest = reference_bundle / "delta_z_derived_run.json"
    reference_manifest.write_bytes(
        (
            json.dumps(
                {
                    "run_id": "derived-off-reference",
                    "outputs": {"trait_summary.tsv": _identity(reference_summary)},
                },
                indent=2,
            )
            + "\n"
        ).encode("utf-8")
    )

    trait_dir = tmp_path / "declared_traits"
    trait_dir.mkdir()
    trait_list = trait_dir / "headline_ten_phenotype.txt"
    trait_list.write_text("METPO:1000608\nMETPO:1000607\n", encoding="utf-8")

    return {
        "tmp_path": tmp_path,
        "parent_bundle": parent_bundle,
        "parent_manifest": parent_manifest,
        "parent_sha256": hashlib.sha256(parent_bytes).hexdigest(),
        "reference_bundle": reference_bundle,
        "reference_manifest": reference_manifest,
        "reference_summary": reference_bundle / "trait_summary.tsv",
        "trait_list": trait_list,
    }


@pytest.mark.parametrize(
    "case",
    [
        "parent_lexical_file",
        "parent_lexical_nested",
        "reference_lexical_file",
        "reference_lexical_nested",
        "output_symlink_to_parent_bundle",
        "output_symlink_to_reference_bundle",
        "parent_input_symlink_output_in_alias",
        "parent_input_symlink_output_in_resolved",
        "reference_input_symlink_output_in_alias",
        "reference_input_symlink_output_in_resolved",
        "lexical_parent_traversal_outside",
        "lexical_reference_traversal_outside",
        "output_file_symlink_into_parent_bundle",
        "output_file_symlink_into_reference_bundle",
        "trait_list_collision",
    ],
)
def test_refuses_output_inside_either_completed_bundle_without_modifying_bundles(
    audit_fixture: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
    case: str,
) -> None:
    tmp_path: Path = audit_fixture["tmp_path"]
    parent_bundle: Path = audit_fixture["parent_bundle"]
    reference_bundle: Path = audit_fixture["reference_bundle"]
    parent_manifest: Path = audit_fixture["parent_manifest"]
    reference_manifest: Path = audit_fixture["reference_manifest"]
    reference_summary: Path = audit_fixture["reference_summary"]
    trait_list: Path = audit_fixture["trait_list"]

    parent_alias = tmp_path / "aliases" / "parent_link"
    reference_alias = tmp_path / "aliases" / "reference_link"
    parent_alias.parent.mkdir(parents=True, exist_ok=True)
    parent_alias.symlink_to(parent_bundle, target_is_directory=True)
    reference_alias.symlink_to(reference_bundle, target_is_directory=True)

    if case == "parent_lexical_file":
        output = parent_bundle / "audit_report.json"
    elif case == "parent_lexical_nested":
        output = parent_bundle / "nested" / "deep" / "audit_report.json"
    elif case == "reference_lexical_file":
        output = reference_bundle / "audit_report.json"
    elif case == "reference_lexical_nested":
        output = reference_bundle / "nested" / "deep" / "audit_report.json"
    elif case == "output_symlink_to_parent_bundle":
        output = parent_alias / "nested" / "audit_report.json"
    elif case == "output_symlink_to_reference_bundle":
        output = reference_alias / "nested" / "audit_report.json"
    elif case == "parent_input_symlink_output_in_alias":
        parent_manifest = parent_alias / "three_test_run.json"
        output = parent_alias / "nested" / "audit_report.json"
    elif case == "parent_input_symlink_output_in_resolved":
        parent_manifest = parent_alias / "three_test_run.json"
        output = parent_bundle / "nested" / "audit_report.json"
    elif case == "reference_input_symlink_output_in_alias":
        reference_manifest = reference_alias / "delta_z_derived_run.json"
        output = reference_alias / "nested" / "audit_report.json"
    elif case == "reference_input_symlink_output_in_resolved":
        reference_manifest = reference_alias / "delta_z_derived_run.json"
        output = reference_bundle / "nested" / "audit_report.json"
    elif case == "lexical_parent_traversal_outside":
        output = parent_bundle / ".." / "escaped_parent" / "audit_report.json"
    elif case == "lexical_reference_traversal_outside":
        output = reference_bundle / ".." / "escaped_reference" / "audit_report.json"
    elif case == "output_file_symlink_into_parent_bundle":
        output = tmp_path / "aliases" / "report_into_parent.json"
        output.symlink_to(parent_bundle / "nested" / "audit_report.json")
    elif case == "output_file_symlink_into_reference_bundle":
        output = tmp_path / "aliases" / "report_into_reference.json"
        output.symlink_to(reference_bundle / "nested" / "audit_report.json")
    else:
        output = trait_list

    parent_before = _tree_state(parent_bundle)
    reference_before = _tree_state(reference_bundle)

    module = _load_module()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT_PATH),
            "--parent-manifest",
            str(parent_manifest),
            "--parent-sha256",
            audit_fixture["parent_sha256"],
            "--trait-list",
            str(trait_list),
            "--reference-file",
            str(reference_manifest),
            "--reference-file",
            str(reference_summary),
            "--output",
            str(output),
        ],
    )
    with pytest.raises(ValueError, match=r"input bundle|disjoint|outside"):
        module["main"]()

    assert _tree_state(parent_bundle) == parent_before
    assert _tree_state(reference_bundle) == reference_before
    if case == "trait_list_collision":
        assert trait_list.read_text(encoding="utf-8") == "METPO:1000608\nMETPO:1000607\n"
    else:
        assert not output.exists()


def test_separate_output_directory_writes_report_and_leaves_bundles_unchanged(
    audit_fixture: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tmp_path: Path = audit_fixture["tmp_path"]
    parent_bundle: Path = audit_fixture["parent_bundle"]
    reference_bundle: Path = audit_fixture["reference_bundle"]
    output = tmp_path / "separate_output" / "nested" / "audit_report.json"

    parent_before = _tree_state(parent_bundle)
    reference_before = _tree_state(reference_bundle)

    module = _load_module()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT_PATH),
            "--parent-manifest",
            str(audit_fixture["parent_manifest"]),
            "--parent-sha256",
            audit_fixture["parent_sha256"],
            "--trait-list",
            str(audit_fixture["trait_list"]),
            "--reference-file",
            str(audit_fixture["reference_manifest"]),
            "--reference-file",
            str(audit_fixture["reference_summary"]),
            "--output",
            str(output),
        ],
    )
    module["main"]()

    assert _tree_state(parent_bundle) == parent_before
    assert _tree_state(reference_bundle) == reference_before
    assert output.is_file()
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["parent_run_id"] == "three-test-off-fixture"
    assert report["list_membership_blocks_off_quadrant"] is False
    assert report["input_bytes_unchanged"] is True
