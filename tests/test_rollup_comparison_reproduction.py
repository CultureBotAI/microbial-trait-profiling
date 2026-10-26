"""Read-only qualification of completed F rollup comparisons (#1041, 2026-10-25).

Explicit environment inputs identify both parent manifests, the completed result,
its clean producer commit and the archived source commit. No pipeline is replayed,
no Delta-z is transformed, and no significance test or correction is recomputed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
from git import Blob, Repo

from mtp.paths import file_sha256

SCOPES = ("off", "phenotype")
ENVIRONMENT = (
    "MTP_ROLLUP_OFF_PARENT",
    "MTP_ROLLUP_PHENOTYPE_PARENT",
    "MTP_ROLLUP_RESULT",
    "MTP_ROLLUP_COMMIT",
    "MTP_ROLLUP_SOURCE_COMMIT",
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


@pytest.fixture(scope="module")
def inputs() -> Inputs:
    if not any(name in os.environ for name in ENVIRONMENT):
        pytest.skip("Requires completed F rollup result, both parents and explicit commit pins")
    missing = [name for name in ENVIRONMENT if not os.environ.get(name, "").strip()]
    if missing:
        raise ValueError(f"Missing rollup qualification inputs: {missing}")
    values = [os.environ[name] for name in ENVIRONMENT]
    assert re.fullmatch(r"[0-9a-f]{40}", values[3])
    assert re.fullmatch(r"[0-9a-f]{40}", values[4])
    return Inputs(Path(values[0]), Path(values[1]), Path(values[2]), values[3], values[4])


@pytest.fixture(scope="module")
def completed(inputs: Inputs) -> tuple[dict[str, Any], dict[str, Parent]]:
    result: dict[str, Any] = json.loads((inputs.result / "rollup_comparison_run.json").read_text())
    assert result["version"] == "off-phenotype-comparison-v1"
    assert result["code_commit"] == inputs.commit
    assert result["code_dirty"] is False
    assert result["primary_scope"] == "phenotype"
    assert result["sensitivity_scope"] == "off"
    assert result["scope_selection"] == "phenotype_primary_selected_after_off_results_were_seen"
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
    with Repo(Path(__file__).resolve().parents[1]) as repository:
        tree = repository.commit(inputs.commit).tree
        for name, digest in result["code_inventory"].items():
            blob = tree / name
            assert isinstance(blob, Blob)
            assert hashlib.sha256(blob.data_stream.read()).hexdigest() == digest
    parents = {}
    for scope, supplied in zip(SCOPES, (inputs.off, inputs.phenotype), strict=True):
        manifest = supplied / "three_test_run.json" if supplied.is_dir() else supplied
        raw = manifest.read_bytes()
        parent: dict[str, Any] = json.loads(raw)
        declared = result["parents"][scope]
        assert declared["manifest_sha256"] == hashlib.sha256(raw).hexdigest()
        assert parent["code_dirty"] is False
        assert parent["strain_rollup"] == scope
        assert parent["contrast_population"] == "disease_vs_healthy_tier1_2_subject_controls"
        assert parent["direction_basis"] == "disease_vs_comparison"
        assert set(parent["source_commits"]) == {"comparison", "core", "noel", "jed", "ilyas"}
        assert set(parent["source_commits"].values()) == {inputs.source_commit}
        for name in ("run_id", "code_commit", "family_size", "retained_traits", "source_commits"):
            assert declared[name] == parent[name]
        for name in ("per_trait_evidence.tsv", "composition_study_delta_z.tsv"):
            identity = parent["outputs"][name]
            path = manifest.parent / name
            assert file_sha256(path) == identity["sha256"]
            assert path.stat().st_size == identity["size_bytes"]
            assert identity in result["inputs"][scope].values()
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
            assert counts == {"off": 9, "phenotype": 10}
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
