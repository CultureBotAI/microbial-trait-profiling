"""Read-only qualification of completed F rollup comparisons (#1041, 2026-10-25).

Explicit environment inputs identify both parent manifests, the completed result,
its clean producer commit, archived source commit and two independent parent SHA-256
pins. No pipeline is replayed,
no Delta-z is transformed, and no significance test or correction is recomputed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
from git import Blob, Repo

from mtp.paths import file_sha256

REPOSITORY = Path(__file__).resolve().parents[1]
SCOPES = ("off", "phenotype")
ENVIRONMENT = (
    "MTP_ROLLUP_OFF_PARENT",
    "MTP_ROLLUP_PHENOTYPE_PARENT",
    "MTP_ROLLUP_RESULT",
    "MTP_ROLLUP_COMMIT",
    "MTP_ROLLUP_SOURCE_COMMIT",
    "MTP_ROLLUP_OFF_SHA256",
    "MTP_ROLLUP_PHENOTYPE_SHA256",
)
RANKS = tuple(
    f"{effect}_rank_{order}"
    for effect in ("experiment_mean_delta_z", "study_mean_delta_z")
    for order in ("ascending", "descending", "absolute_descending")
)
Membership = tuple[bool, bool, str, int | None, float | None, int]


@dataclass(frozen=True)
class Inputs:
    off: Path
    phenotype: Path
    result: Path
    commit: str
    source_commit: str
    off_sha256: str
    phenotype_sha256: str


@dataclass(frozen=True)
class Parent:
    manifest: Path
    record: dict[str, Any]
    evidence: dict[str, dict[str, str]]
    studies: list[dict[str, str]]


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        assert len(reader.fieldnames) == len(set(reader.fieldnames))
        rows = list(reader)
    assert all(all(isinstance(value, str) for value in row.values()) for row in rows)
    return rows


def keyed(rows: list[dict[str, str]], field: str) -> dict[str, dict[str, str]]:
    result = {row[field]: row for row in rows}
    assert len(result) == len(rows)
    assert all(result)
    return result


def same_value(observed: str, expected: str, column: str) -> None:
    """Permit numerical spelling changes only; require the identical binary float value."""
    if observed == expected:
        return
    assert observed
    assert expected
    assert column not in {"trait", "trait_label"}
    try:
        assert float(observed) == float(expected)
    except ValueError as error:
        raise AssertionError(f"Changed imported text: {column}") from error


def nullable_boolean(value: str) -> bool | None:
    assert value in {"", "True", "False"}
    return None if value == "" else value == "True"


def assert_scope_selection(result: dict[str, Any]) -> None:
    """Check the recorded scope choice independently of artifact qualification."""
    assert result["primary_scope"] == "phenotype"
    assert result["sensitivity_scope"] == "off"
    assert result["scope_selection"] == (
        "phenotype_primary_selected_after_off_and_phenotype_results_were_seen"
    )


def verify_producer_inventory(
    result: dict[str, Any],
    producer: str,
    repository_root: Path = REPOSITORY,
) -> None:
    """Derive the complete producer inventory from the caller-pinned Git tree."""
    assert result["code_commit"] == producer
    assert result["code_dirty"] is False
    required = {
        "Justfile",
        "pyproject.toml",
        "uv.lock",
        "methods/three_test_table.md",
        "methods/off_phenotype_comparison.md",
    }
    with Repo(repository_root) as repository:
        tree = repository.commit(producer).tree
        required.update(
            str(item.path)
            for item in tree.traverse()
            if isinstance(item, Blob)
            and str(item.path).startswith("src/mtp/")
            and str(item.path).endswith(".py")
        )
        assert any(name.startswith("src/mtp/") for name in required)
        expected = {}
        for name in sorted(required):
            blob = tree / name
            assert isinstance(blob, Blob)
            assert blob.mode in {0o100644, 0o100755}
            expected[name] = hashlib.sha256(blob.data_stream.read()).hexdigest()
    assert result["code_inventory"] == expected
    encoded = json.dumps(
        expected, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")
    assert result["implementation_sha256"] == hashlib.sha256(encoded).hexdigest()
    assert result["methods_sha256"] == expected["methods/off_phenotype_comparison.md"]
    assert result["lock_sha256"] == expected["uv.lock"]


def pinned_parent_manifest(
    path: Path, expected_sha256: str, declared_sha256: str
) -> dict[str, Any]:
    """Read parent bytes once, then check their declared identity before parsing."""
    assert re.fullmatch(r"[0-9a-f]{64}", expected_sha256)
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == expected_sha256
    assert declared_sha256 == expected_sha256
    result: dict[str, Any] = json.loads(raw)
    assert isinstance(result, dict)
    return result


def assert_native_noel_reference(counts: dict[str, int]) -> None:
    """Keep the native Noel nominal-BY reference distinct from study-t counts."""
    # Native F trait_screen_report.md reports four at nominal BY <= 0.05 per scope.
    assert counts == {"off": 4, "phenotype": 4}


def parse_inputs(environment: Mapping[str, str]) -> Inputs | None:
    """Parse explicit qualification inputs without consulting source metadata."""
    if not any(name in environment for name in ENVIRONMENT):
        return None
    missing = [name for name in ENVIRONMENT if not environment.get(name, "").strip()]
    if missing:
        raise ValueError(f"Missing rollup qualification inputs: {missing}")
    values = [environment[name] for name in ENVIRONMENT]
    assert re.fullmatch(r"[0-9a-f]{40}", values[3])
    assert re.fullmatch(r"[0-9a-f]{40}", values[4])
    for name, value in zip(ENVIRONMENT[5:], values[5:], strict=True):
        if not re.fullmatch(r"[0-9a-f]{64}", value):
            raise ValueError(f"Expected a lowercase SHA-256 pin for {name}")
    return Inputs(
        Path(values[0]),
        Path(values[1]),
        Path(values[2]),
        values[3],
        values[4],
        values[5],
        values[6],
    )


@pytest.fixture(scope="module")
def inputs() -> Inputs:
    parsed = parse_inputs(os.environ)
    if parsed is None:
        pytest.skip("Requires completed F rollup result, both parents and explicit commit pins")
    return parsed


def verify_parent_binding(
    declared: dict[str, Any], parent: dict[str, Any], consumed_inputs: dict[str, Any]
) -> None:
    """Bind every copied parent claim and exact consumed-table identity (#1109)."""
    fields = (
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
    )
    # The independently supplied manifest SHA was already checked against parsed bytes.
    assert set(declared) == {"manifest_sha256", "input_inventory", *fields}
    for name in fields:
        assert declared[name] == parent[name]
    assert declared["input_inventory"] == parent["inputs"]
    for name in ("per_trait_evidence.tsv", "composition_study_delta_z.tsv"):
        identity = parent["outputs"][name]
        assert consumed_inputs.get(f"{name}@sha256:{identity['sha256']}") == identity


@pytest.fixture(scope="module")
def completed(inputs: Inputs) -> tuple[dict[str, Any], dict[str, Parent]]:
    result: dict[str, Any] = json.loads((inputs.result / "rollup_comparison_run.json").read_text())
    assert result["version"] == "off-phenotype-comparison-v1"
    assert_scope_selection(result)
    assert set(result["outputs"]) == {
        "shared_traits.tsv",
        "trait_catalog.tsv",
        "rank_transitions.tsv",
        "off_per_trait_evidence.tsv",
        "phenotype_per_trait_evidence.tsv",
    }
    assert {path.name for path in inputs.result.iterdir()} == {
        "rollup_comparison_run.json",
        *result["outputs"],
    }
    for name, identity in result["outputs"].items():
        path = inputs.result / name
        assert not path.is_symlink()
        assert file_sha256(path) == identity["sha256"]
        assert path.stat().st_size == identity["size_bytes"]
    verify_producer_inventory(result, inputs.commit)
    parents = {}
    for scope, supplied, expected_sha256 in zip(
        SCOPES,
        (inputs.off, inputs.phenotype),
        (inputs.off_sha256, inputs.phenotype_sha256),
        strict=True,
    ):
        manifest = supplied / "three_test_run.json" if supplied.is_dir() else supplied
        declared = result["parents"][scope]
        parent = pinned_parent_manifest(manifest, expected_sha256, declared["manifest_sha256"])
        assert parent["code_dirty"] is False
        assert parent["strain_rollup"] == scope
        assert parent["contrast_population"] == "disease_vs_healthy_tier1_2_subject_controls"
        assert parent["direction_basis"] == "disease_vs_comparison"
        assert set(parent["source_commits"]) == {"comparison", "core", "noel", "jed", "ilyas"}
        assert set(parent["source_commits"].values()) == {inputs.source_commit}
        verify_parent_binding(declared, parent, result["inputs"][scope])
        for name in ("per_trait_evidence.tsv", "composition_study_delta_z.tsv"):
            identity = parent["outputs"][name]
            path = manifest.parent / name
            assert file_sha256(path) == identity["sha256"]
            assert path.stat().st_size == identity["size_bytes"]
        evidence = keyed(read_table(manifest.parent / "per_trait_evidence.tsv"), "trait")
        assert list(evidence) == parent["retained_traits"]
        studies = read_table(manifest.parent / "composition_study_delta_z.tsv")
        assert len(keyed(studies, "study_id")) == 548
        parents[scope] = Parent(manifest, parent, evidence, studies)
    assert parents["off"].record["code_commit"] == parents["phenotype"].record["code_commit"]
    return result, parents


def active_studies(parent: Parent) -> dict[str, int]:
    return {
        trait: sum(abs(float(row[trait])) > 1e-12 for row in parent.studies)
        for trait in parent.evidence
    }


@pytest.mark.integration
def test_f_identity_catalog_and_exact_parent_copies(
    inputs: Inputs, completed: tuple[dict[str, Any], dict[str, Parent]]
) -> None:
    result, parents = completed
    assert result["counts"] == {
        "off_traits": 317,
        "phenotype_traits": 344,
        "shared_traits": 317,
        "union_traits": 344,
        "comparisons": 935,
        "studies": 548,
    }
    off, phenotype = set(parents["off"].evidence), set(parents["phenotype"].evidence)
    assert len(off & phenotype) == 317
    assert len(phenotype - off) == 27
    assert not off - phenotype
    catalog = keyed(read_table(inputs.result / "trait_catalog.tsv"), "trait")
    assert set(catalog) == off | phenotype
    assert catalog["METPO:1000608"]["membership"] == "shared"
    assert catalog["METPO:1000605"]["membership"] == "phenotype_only"
    for trait, row in catalog.items():
        for scope in SCOPES:
            available = trait in parents[scope].evidence
            assert nullable_boolean(row[f"in_{scope}"]) == available
            expected = parents[scope].evidence[trait]["trait_label"] if available else ""
            assert row[f"{scope}_trait_label"] == expected
    for scope, parent in parents.items():
        original = parent.manifest.parent / "per_trait_evidence.tsv"
        assert (
            inputs.result / f"{scope}_per_trait_evidence.tsv"
        ).read_bytes() == original.read_bytes()


@pytest.mark.integration
def test_every_shared_parent_field_and_full_family_rank_is_preserved(
    inputs: Inputs, completed: tuple[dict[str, Any], dict[str, Parent]]
) -> None:
    _, parents = completed
    shared = keyed(read_table(inputs.result / "shared_traits.tsv"), "trait")
    assert set(shared) == set(parents["off"].evidence) & set(parents["phenotype"].evidence)
    for scope, parent in parents.items():
        support = active_studies(parent)
        for trait, observed in shared.items():
            original = parent.evidence[trait]
            for column, expected in original.items():
                if column != "trait":
                    same_value(observed[f"{scope}_{column}"], expected, column)
            for column in RANKS:
                assert int(observed[f"{scope}_{column}"]) == int(original[column])
            assert int(observed[f"{scope}_composition_family_size"]) == len(parent.evidence)
            assert int(observed[f"{scope}_noel_family_size"]) == len(parent.evidence)
            same_value(
                observed[f"{scope}_comparison_delta_z_mean"],
                original["experiment_mean_delta_z"],
                "experiment_mean_delta_z",
            )
            up = sum(float(row[trait]) > 1e-12 for row in parent.studies)
            down = sum(float(row[trait]) < -1e-12 for row in parent.studies)
            assert int(observed[f"{scope}_sdz_up"]) == up
            assert int(observed[f"{scope}_sdz_down"]) == down
            assert int(observed[f"{scope}_sdz_active_studies"]) == support[trait] == up + down
            assert int(observed[f"{scope}_sdz_tied"]) == len(parent.studies) - up - down


def expected_members(
    parent: Parent, quantity: str, direction: str, top: int, minimum: int
) -> dict[str, Membership]:
    """Independent scalar ranking and thresholds over the entire parent catalog."""
    support = active_studies(parent)
    result: dict[str, Membership] = {}
    if top:
        effects = {trait: float(row[quantity]) for trait, row in parent.evidence.items()}
        eligible_effects = {
            trait: value
            for trait, value in effects.items()
            if (value > 0 if direction == "positive" else value < 0) and support[trait] >= minimum
        }
        for trait, value in effects.items():
            rank = (
                1
                + sum(
                    other > value if direction == "positive" else other < value
                    for other in eligible_effects.values()
                )
                if trait in eligible_effects
                else None
            )
            selected = rank is not None and rank <= top
            status = (
                "selected" if selected else "not_selected" if rank is not None else "ineligible"
            )
            result[trait] = (
                trait in eligible_effects,
                selected,
                status,
                rank,
                value,
                len(eligible_effects),
            )
        return result
    status_column = "noel_status" if quantity == "noel_q_by_nominal" else "study_t_status"
    result = {}
    for trait, row in parent.evidence.items():
        tested = row[status_column] == "tested"
        probability = float(row[quantity]) if tested else None
        selected = probability is not None and probability <= 0.05
        status = "selected" if selected else "not_selected" if tested else "untestable"
        result[trait] = tested, selected, status, None, probability, len(parent.evidence)
    return result


@pytest.mark.integration
def test_all_transitions_against_independent_full_scope_selections(
    inputs: Inputs, completed: tuple[dict[str, Any], dict[str, Parent]]
) -> None:
    _, parents = completed
    views = [
        ("noel_by_0_05", "noel_q_by_nominal", "not_applicable", 0, 0),
        ("study_t_by_0_05", "study_t_q_by", "not_applicable", 0, 0),
        *(
            (
                f"comparison_mean_{direction}_top{top}_{label}",
                "experiment_mean_delta_z",
                direction,
                top,
                minimum,
            )
            for label, minimum in (("all", 0), ("active50", 50))
            for direction in ("positive", "negative")
            for top in (5, 10)
        ),
    ]
    rows = read_table(inputs.result / "rank_transitions.tsv")
    actual = {(row["view"], row["trait"]): row for row in rows}
    union = set(parents["off"].evidence) | set(parents["phenotype"].evidence)
    assert len(actual) == len(rows) == len(union) * 10
    assert set(actual) == {(view[0], trait) for view in views for trait in union}
    for name, quantity, direction, top, minimum in views:
        expected = {
            scope: expected_members(parent, quantity, direction, top, minimum)
            for scope, parent in parents.items()
        }
        counts = {
            scope: sum(value[1] for value in members.values())
            for scope, members in expected.items()
        }
        if name == "noel_by_0_05":
            assert_native_noel_reference(counts)
        if not top:
            print(json.dumps({"view": name, "selected": counts}, sort_keys=True))
        for trait in union:
            row = actual[name, trait]
            selected = {}
            for scope, parent in parents.items():
                available = trait in parent.evidence
                assert nullable_boolean(row[f"{scope}_available"]) == available
                if not available:
                    assert row[f"{scope}_status"] == "unavailable"
                    assert row[f"{scope}_selected"] == row[f"{scope}_eligible"] == ""
                    selected[scope] = False
                    continue
                eligible, chosen, status, rank, value, size = expected[scope][trait]
                selected[scope] = chosen
                assert nullable_boolean(row[f"{scope}_eligible"]) == eligible
                assert nullable_boolean(row[f"{scope}_selected"]) == chosen
                assert row[f"{scope}_status"] == status
                assert (
                    int(row[f"{scope}_view_rank"]) if row[f"{scope}_view_rank"] else None
                ) == rank
                assert (
                    float(row[f"{scope}_selection_value"])
                    if row[f"{scope}_selection_value"]
                    else None
                ) == value
                assert int(row[f"{scope}_retained_family_size"]) == len(parent.evidence)
                assert int(row[f"{scope}_selection_family_size"]) == size
            old, new = selected["off"], selected["phenotype"]
            state = (
                "retained" if old and new else "entered" if new else "left" if old else "neither"
            )
            assert row["membership_state"] == state
            assert row["selection_quantity"] == quantity
            assert row["direction"] == direction
            assert float(row["cutoff"]) == (top or 0.05)
            assert int(row["minimum_active_studies"]) == minimum


def test_qualification_numeric_comparison_preserves_adjacent_float_values() -> None:
    same_value("0.3", "0.3", "study_t_q_by")
    same_value("3.0", "3", "composition_family_size")
    with pytest.raises(AssertionError):
        same_value("0.3", "0.30000000000000004", "study_t_q_by")
