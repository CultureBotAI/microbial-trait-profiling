"""Refuse audit report destinations inside completed input bundles (#1186, #1212)."""

from __future__ import annotations

import hashlib
import json
import runpy
import sys
import tempfile
from pathlib import Path
from typing import Any

import pytest

SCRIPT_PATH = Path(__file__).parents[1] / "scripts/adhoc/issue1115_audit_off_quadrant_inputs.py"
MACOS = pytest.mark.skipif(
    sys.platform != "darwin",
    reason="macOS APFS case-insensitivity and /System/Volumes/Data firmlinks",
)
BUNDLE_REFUSAL_MATCH = (
    r"^Report output must be outside completed input bundles and disjoint from audit inputs$"
)
SYMLINK_CYCLE_MATCH = r"^Symlink cycle while protecting input bundles$"
SYMLINK_LIMIT_MATCH = r"^Symlink expansion limit reached while protecting input bundles$"


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
        "output_is_parent_of_bundles",
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
    elif case == "output_is_parent_of_bundles":
        output = tmp_path
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
    with pytest.raises(ValueError, match=BUNDLE_REFUSAL_MATCH):
        module["main"]()

    assert _tree_state(parent_bundle) == parent_before
    assert _tree_state(reference_bundle) == reference_before
    if case == "trait_list_collision":
        assert trait_list.read_text(encoding="utf-8") == "METPO:1000608\nMETPO:1000607\n"
    elif case == "output_is_parent_of_bundles":
        assert output.is_dir()
    else:
        assert not output.exists()


@MACOS
@pytest.mark.parametrize(
    "case",
    [
        "parent_case_variant",
        "reference_case_variant_nested",
        "users_prefix_case_variant",
        "parent_firmlink",
        "reference_firmlink_nested",
        "case_variant_reference_file_input",
        "trait_list_case_variant",
    ],
)
def test_refuses_macos_case_variant_and_firmlink_bundle_spellings(
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

    if case == "parent_case_variant":
        output = tmp_path / "PARENT_BUNDLE" / "a.json"
    elif case == "reference_case_variant_nested":
        output = tmp_path / "Reference_Bundle" / "x" / "a.json"
    elif case == "users_prefix_case_variant":
        scratch_base = Path(__file__).resolve().parents[1] / ".tmp"
        scratch_base.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="audit-users-", dir=scratch_base) as raw_dir:
            users_root = Path(raw_dir).resolve()
            users_parent = users_root / "parent_bundle"
            users_parent.mkdir()
            for item in parent_bundle.iterdir():
                (users_parent / item.name).write_bytes(item.read_bytes())
            users_manifest = users_parent / "three_test_run.json"
            assert users_parent.parts[1] == "Users"
            upper_users_parent = Path("/USERS", *users_parent.parts[2:])
            output = upper_users_parent / "a.json"
            users_before = _tree_state(users_parent)
            module = _load_module()
            monkeypatch.setattr(
                sys,
                "argv",
                [
                    str(SCRIPT_PATH),
                    "--parent-manifest",
                    str(users_manifest),
                    "--parent-sha256",
                    audit_fixture["parent_sha256"],
                    "--trait-list",
                    str(trait_list),
                    "--reference-file",
                    str(reference_manifest),
                    "--output",
                    str(output),
                ],
            )
            with pytest.raises(ValueError, match=BUNDLE_REFUSAL_MATCH):
                module["main"]()
            assert _tree_state(users_parent) == users_before
            assert not (users_parent / "a.json").exists()
        if scratch_base.exists() and not any(scratch_base.iterdir()):
            scratch_base.rmdir()
        return
    elif case == "parent_firmlink":
        output = Path("/System/Volumes/Data" + str(parent_bundle.resolve())) / "a.json"
    elif case == "reference_firmlink_nested":
        output = Path("/System/Volumes/Data" + str(reference_bundle.resolve())) / "x" / "a.json"
    elif case == "case_variant_reference_file_input":
        reference_manifest = tmp_path / "REFERENCE_BUNDLE" / "delta_z_derived_run.json"
        reference_summary = tmp_path / "REFERENCE_BUNDLE" / "trait_summary.tsv"
        output = reference_bundle / "a.json"
    else:
        output = trait_list.parent / trait_list.name.upper()

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
    with pytest.raises(ValueError, match=BUNDLE_REFUSAL_MATCH):
        module["main"]()

    assert _tree_state(parent_bundle) == parent_before
    assert _tree_state(reference_bundle) == reference_before


def test_refuses_symlink_cycle_with_pinned_message(
    audit_fixture: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tmp_path: Path = audit_fixture["tmp_path"]
    loop_a = tmp_path / "loop_a"
    loop_b = tmp_path / "loop_b"
    loop_a.symlink_to(loop_b)
    loop_b.symlink_to(loop_a)

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
            "--output",
            str(loop_a / "audit_report.json"),
        ],
    )
    with pytest.raises(ValueError, match=SYMLINK_CYCLE_MATCH):
        module["main"]()


def test_refuses_symlink_expansion_limit_with_pinned_message(
    audit_fixture: dict[str, Any],
) -> None:
    tmp_path: Path = audit_fixture["tmp_path"]
    separate_dir = tmp_path / "separate_hop_target"
    separate_dir.mkdir()
    module = _load_module()
    chain_dir = tmp_path / "hop_chain"
    chain_dir.mkdir()
    target = separate_dir
    links: list[Path] = []
    for index in range(41):
        link = chain_dir / f"hop_{index}"
        link.symlink_to(target, target_is_directory=True)
        links.append(link)
        target = link

    validate_output_path = module["validate_output_path"]
    validate_output_path(
        links[39] / "audit_report.json",
        audit_fixture["parent_manifest"],
        audit_fixture["trait_list"],
        [audit_fixture["reference_manifest"]],
    )
    with pytest.raises(ValueError, match=SYMLINK_LIMIT_MATCH):
        validate_output_path(
            links[40] / "audit_report.json",
            audit_fixture["parent_manifest"],
            audit_fixture["trait_list"],
            [audit_fixture["reference_manifest"]],
        )


def test_accepts_symlinked_inputs_with_separate_output_directory(
    audit_fixture: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tmp_path: Path = audit_fixture["tmp_path"]
    parent_bundle: Path = audit_fixture["parent_bundle"]
    reference_bundle: Path = audit_fixture["reference_bundle"]
    aliases = tmp_path / "input_aliases"
    aliases.mkdir()
    parent_alias = aliases / "parent_bundle_link"
    reference_alias = aliases / "reference_bundle_link"
    trait_list_alias = aliases / "headline_ten_phenotype_link.txt"
    reference_summary_alias = aliases / "trait_summary_link.tsv"
    parent_alias.symlink_to(parent_bundle, target_is_directory=True)
    reference_alias.symlink_to(reference_bundle, target_is_directory=True)
    trait_list_alias.symlink_to(audit_fixture["trait_list"])
    reference_summary_alias.symlink_to(audit_fixture["reference_summary"])

    output = tmp_path / "separate_from_symlinks" / "audit_report.json"
    parent_before = _tree_state(parent_bundle)
    reference_before = _tree_state(reference_bundle)

    module = _load_module()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT_PATH),
            "--parent-manifest",
            str(parent_alias / "three_test_run.json"),
            "--parent-sha256",
            audit_fixture["parent_sha256"],
            "--trait-list",
            str(trait_list_alias),
            "--reference-file",
            str(reference_alias / "delta_z_derived_run.json"),
            "--reference-file",
            str(reference_summary_alias),
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
    assert report["input_bytes_unchanged"] is True


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
