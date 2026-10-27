"""Independent qualification-reference regressions (#1105/#1108/#1109).

Synthetic Git trees, caller pins and manifest bytes exercise the actual
qualification helpers. No completed scientific run is generated or replayed.
Date: 2026-10-25.
"""

from __future__ import annotations

import hashlib
import json
import tempfile
from collections.abc import Iterator
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest
from git import Actor, Repo

import test_rollup_comparison_reproduction as qualification

REPOSITORY = Path(__file__).resolve().parents[1]
FILES = {
    "src/mtp/__init__.py": '"""Synthetic package."""\n',
    "src/mtp/core.py": '"""Synthetic implementation."""\n',
    "src/mtp/nested/analysis.py": '"""Nested implementation must be bound."""\n',
    "Justfile": "fixture:\n    echo fixture\n",
    "pyproject.toml": "[project]\nname = 'fixture'\n",
    "uv.lock": "version = 1\n",
    "methods/three_test_table.md": "Parent methods fixture.\n",
    "methods/off_phenotype_comparison.md": "Comparison methods fixture.\n",
}
BOTH_OBSERVED = "phenotype_primary_selected_after_off_and_phenotype_results_were_seen"
OLD_DISCLOSURE = "phenotype_primary_selected_after_off_results_were_seen"


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    root = REPOSITORY / ".tmp" / "rollup-reference-tests"
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as directory:
        yield Path(directory)


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def inventory_digest(inventory: dict[str, str]) -> str:
    return digest(
        json.dumps(
            inventory, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
        ).encode("utf-8")
    )


@pytest.fixture
def producer(tmp_path: Path) -> tuple[Path, str, dict[str, Any]]:
    for name, content in {**FILES, "unrelated.txt": "Outside the producer contract.\n"}.items():
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    actor = Actor("Fixture", "fixture@example.invalid")
    with Repo.init(tmp_path) as repository:
        repository.index.add([*FILES, "unrelated.txt"])
        commit = repository.index.commit(
            "Synthetic qualification source", author=actor, committer=actor
        )
    inventory = {name: digest(content.encode("utf-8")) for name, content in FILES.items()}
    record = {
        "code_commit": commit.hexsha,
        "code_dirty": False,
        "code_inventory": inventory,
        "implementation_sha256": inventory_digest(inventory),
        "methods_sha256": inventory["methods/off_phenotype_comparison.md"],
        "lock_sha256": inventory["uv.lock"],
    }
    return tmp_path, commit.hexsha, record


def environment() -> dict[str, str]:
    return {
        "MTP_ROLLUP_OFF_PARENT": "off parent/three_test_run.json",
        "MTP_ROLLUP_PHENOTYPE_PARENT": "phenotype parent/three_test_run.json",
        "MTP_ROLLUP_RESULT": "completed result",
        "MTP_ROLLUP_COMMIT": "a" * 40,
        "MTP_ROLLUP_SOURCE_COMMIT": "b" * 40,
        "MTP_ROLLUP_OFF_SHA256": "c" * 64,
        "MTP_ROLLUP_PHENOTYPE_SHA256": "d" * 64,
    }


def test_qualification_requires_both_independent_pins_and_preserves_their_order() -> None:
    parsed = qualification.parse_inputs(environment())
    assert parsed is not None
    assert parsed.off_sha256 == "c" * 64
    assert parsed.phenotype_sha256 == "d" * 64
    assert parsed.commit == "a" * 40
    assert parsed.source_commit == "b" * 40
    assert qualification.parse_inputs({}) is None


@pytest.mark.parametrize("name", ["MTP_ROLLUP_OFF_SHA256", "MTP_ROLLUP_PHENOTYPE_SHA256"])
def test_missing_independent_parent_pin_is_rejected(name: str) -> None:
    values = environment()
    del values[name]
    with pytest.raises(ValueError, match="Missing"):
        qualification.parse_inputs(values)


@pytest.mark.parametrize("bad_pin", ["", " ", "not-a-sha256", "e" * 63, "g" * 64])
def test_malformed_independent_parent_pin_is_rejected(bad_pin: str) -> None:
    values = environment()
    values["MTP_ROLLUP_OFF_SHA256"] = bad_pin
    with pytest.raises(ValueError, match=r"Missing|SHA-256"):
        qualification.parse_inputs(values)


def test_exact_parent_pin_is_checked_against_the_bytes_that_are_parsed(tmp_path: Path) -> None:
    path = tmp_path / "three_test_run.json"
    raw = b'{"run_id":"original"}\n'
    path.write_bytes(raw)
    assert qualification.pinned_parent_manifest(path, digest(raw), digest(raw)) == {
        "run_id": "original"
    }
    with pytest.raises(AssertionError):
        qualification.pinned_parent_manifest(path, "f" * 64, digest(raw))


def test_mutated_parent_and_updated_record_do_not_replace_the_callers_pin(tmp_path: Path) -> None:
    path = tmp_path / "three_test_run.json"
    original = b'{"run_id":"original"}\n'
    changed = b'{"run_id":"different"}\n'
    path.write_bytes(changed)
    with pytest.raises(AssertionError):
        qualification.pinned_parent_manifest(path, digest(original), digest(changed))


def test_independent_pin_must_also_match_the_consumer_declared_parent(tmp_path: Path) -> None:
    path = tmp_path / "three_test_run.json"
    raw = b'{"run_id":"original"}\n'
    path.write_bytes(raw)
    with pytest.raises(AssertionError):
        qualification.pinned_parent_manifest(path, digest(raw), "e" * 64)


def test_complete_git_derived_inventory_and_digests_are_accepted(
    producer: tuple[Path, str, dict[str, Any]],
) -> None:
    root, commit, record = producer
    qualification.verify_producer_inventory(record, commit, root)


@pytest.mark.parametrize("change", ["empty", "drop_nested", "drop_recipe", "bonus"])
def test_inventory_membership_is_independently_derived_from_the_full_git_tree(
    producer: tuple[Path, str, dict[str, Any]], change: str
) -> None:
    root, commit, record = producer
    inventory = record["code_inventory"]
    if change == "empty":
        inventory.clear()
    elif change == "drop_nested":
        del inventory["src/mtp/nested/analysis.py"]
    elif change == "drop_recipe":
        del inventory["Justfile"]
    else:
        inventory["unrelated.txt"] = digest((root / "unrelated.txt").read_bytes())
    record["implementation_sha256"] = inventory_digest(inventory)
    with pytest.raises(AssertionError):
        qualification.verify_producer_inventory(record, commit, root)


@pytest.mark.parametrize("field", ["implementation_sha256", "methods_sha256", "lock_sha256"])
def test_aggregate_and_named_digests_are_bound_to_the_expected_git_blobs(
    producer: tuple[Path, str, dict[str, Any]], field: str
) -> None:
    root, commit, record = producer
    record[field] = "0" * 64
    with pytest.raises(AssertionError):
        qualification.verify_producer_inventory(record, commit, root)


def test_modified_member_digest_is_rejected(producer: tuple[Path, str, dict[str, Any]]) -> None:
    root, commit, record = producer
    record["code_inventory"]["src/mtp/core.py"] = "0" * 64
    record["implementation_sha256"] = inventory_digest(record["code_inventory"])
    with pytest.raises(AssertionError):
        qualification.verify_producer_inventory(record, commit, root)


def test_scope_disclosure_names_both_observed_results() -> None:
    qualification.assert_scope_selection(
        {"primary_scope": "phenotype", "sensitivity_scope": "off", "scope_selection": BOTH_OBSERVED}
    )


def test_stale_off_only_disclosure_is_rejected() -> None:
    with pytest.raises(AssertionError):
        qualification.assert_scope_selection(
            {
                "primary_scope": "phenotype",
                "sensitivity_scope": "off",
                "scope_selection": OLD_DISCLOSURE,
            }
        )


def test_native_noel_reference_is_four_per_scope() -> None:
    qualification.assert_native_noel_reference({"off": 4, "phenotype": 4})


def test_study_t_counts_are_not_accepted_as_the_noel_reference() -> None:
    with pytest.raises(AssertionError):
        qualification.assert_native_noel_reference({"off": 9, "phenotype": 10})


EXPORTED_PARENT_FIELDS = (
    "run_id",
    "code_commit",
    "code_dirty",
    "implementation_sha256",
    "environment",
    "strain_rollup",
    "scope_role",
    "contrast_population",
    "direction_basis",
    "family_size",
    "retained_traits",
    "source_commits",
    "core_run_id",
    "comparison_run_id",
    "source_run_ids",
    "analysis_bases",
    "input_inventory",
)
CONSUMED_TABLES = ("per_trait_evidence.tsv", "composition_study_delta_z.tsv")


@pytest.fixture
def parent_binding() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """A pinned parent and its independently copied consumer declaration (#1109)."""
    parent: dict[str, Any] = {
        "run_id": "parent-run",
        "code_commit": "a" * 40,
        "code_dirty": False,
        "implementation_sha256": "b" * 64,
        "environment": {"python": "3.11", "numpy": "2", "pandas": "3", "scipy": "1"},
        "strain_rollup": "off",
        "scope_role": "sensitivity",
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "direction_basis": "disease_vs_comparison",
        "family_size": 1,
        "retained_traits": ["TEST:trait"],
        "source_commits": dict.fromkeys(("core", "comparison", "noel", "jed", "ilyas"), "c" * 40),
        "core_run_id": "core-run",
        "comparison_run_id": "comparison-run",
        "source_run_ids": {"noel": "noel-run", "jed": "jed-run", "ilyas": "ilyas-run"},
        "analysis_bases": {"count_profiles": {"weighting_policy": "equal_study"}},
        "inputs": {"core.tsv@sha256:" + "d" * 64: {"sha256": "d" * 64, "size_bytes": 13}},
        "outputs": {
            name: {"sha256": letter * 64, "size_bytes": size}
            for name, letter, size in zip(CONSUMED_TABLES, ("e", "f"), (17, 23), strict=True)
        },
    }
    declared = {
        name: deepcopy(parent[name]) for name in EXPORTED_PARENT_FIELDS if name != "input_inventory"
    }
    declared.update({"manifest_sha256": "1" * 64, "input_inventory": deepcopy(parent["inputs"])})
    consumed = {
        f"{name}@sha256:{identity['sha256']}": deepcopy(identity)
        for name, identity in parent["outputs"].items()
    }
    return declared, parent, consumed


def test_complete_parent_declaration_and_keyed_inputs_are_accepted(
    parent_binding: tuple[dict[str, Any], dict[str, Any], dict[str, Any]],
) -> None:
    qualification.verify_parent_binding(*parent_binding)


@pytest.mark.parametrize("field", EXPORTED_PARENT_FIELDS)
def test_every_copied_parent_field_is_bound_to_the_pinned_parent(
    parent_binding: tuple[dict[str, Any], dict[str, Any], dict[str, Any]],
    field: str,
) -> None:
    declared, parent, consumed = parent_binding
    declared[field] = {"corrupted": True}
    with pytest.raises(AssertionError):
        qualification.verify_parent_binding(declared, parent, consumed)


@pytest.mark.parametrize("change", ["missing", "extra"])
def test_exported_parent_field_inventory_is_exact(
    parent_binding: tuple[dict[str, Any], dict[str, Any], dict[str, Any]],
    change: str,
) -> None:
    declared, parent, consumed = parent_binding
    if change == "missing":
        del declared["source_run_ids"]
    else:
        declared["unrecorded_source_claim"] = "unsupported"
    with pytest.raises(AssertionError):
        qualification.verify_parent_binding(declared, parent, consumed)


@pytest.mark.parametrize("name", CONSUMED_TABLES)
@pytest.mark.parametrize("change", ["wrong_name", "wrong_digest", "wrong_size", "missing"])
def test_consumed_parent_tables_require_exact_identity_keys_and_values(
    parent_binding: tuple[dict[str, Any], dict[str, Any], dict[str, Any]],
    name: str,
    change: str,
) -> None:
    declared, parent, consumed = parent_binding
    identity = parent["outputs"][name]
    key = f"{name}@sha256:{identity['sha256']}"
    if change == "wrong_name":
        consumed[f"different.tsv@sha256:{identity['sha256']}"] = consumed.pop(key)
    elif change == "wrong_digest":
        consumed[f"{name}@sha256:{'0' * 64}"] = consumed.pop(key)
    elif change == "wrong_size":
        consumed[key]["size_bytes"] += 1
    else:
        del consumed[key]
    with pytest.raises(AssertionError):
        qualification.verify_parent_binding(declared, parent, consumed)


def test_other_consumed_sources_remain_allowed_when_parent_table_keys_are_exact(
    parent_binding: tuple[dict[str, Any], dict[str, Any], dict[str, Any]],
) -> None:
    declared, parent, consumed = parent_binding
    consumed["other.tsv@sha256:" + "0" * 64] = {"sha256": "0" * 64, "size_bytes": 19}
    qualification.verify_parent_binding(declared, parent, consumed)
