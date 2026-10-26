"""Fail-closed paired-rollup contracts; synthetic persisted evidence (#1041, 2026-10-25)."""

from __future__ import annotations

import hashlib
import json
import tempfile
from collections.abc import Iterator
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from mtp.datamodel.workflow_pydantic import (
    AnalysisBasis,
    ArmRun,
    ArtifactInstance,
    ComparisonRun,
    FileInventoryEntry,
    ParameterValue,
    WorkflowRun,
)
from mtp.multiplicity.sources import VerifiedSources
from mtp.rollup_comparison.validation import validate_evidence, verify_comparable
from mtp.three_test_table import reader
from mtp.three_test_table.reader import CompletedTable

POPULATION = "disease_vs_healthy_tier1_2_subject_controls"
DIRECTION = "disease_vs_comparison"
SOURCE_COMMIT = "a" * 40
PRODUCER = "b" * 40
SCOPE_PRODUCTS = {
    "count_coverage",
    "count_profiles",
    "feature_catalog",
    "feature_metadata",
    "ora_coverage",
    "ora_profiles",
    "taxa_matrix",
    "source_identity",
}


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    root = Path(__file__).resolve().parents[1] / ".tmp" / "rollup-validation-tests"
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=root) as directory:
        yield Path(directory)


def completed(scope: str, root: Path) -> CompletedTable:
    """Small accepted-reader output; no inference, real source data or file IO."""
    traits = ["TEST:shift", "TEST:zero", "TEST:constant"]
    if scope == "phenotype":
        traits.append("TEST:added")
    pairs = pd.DataFrame(
        {
            trait: (
                [0.0, 0.0, 0.0]
                if trait.endswith("zero")
                else [1.0, 1.0, 1.0]
                if trait.endswith("constant")
                else [-1.0, 0.0, 2.0]
            )
            for trait in traits
        },
        index=pd.MultiIndex.from_tuples(
            [("s1", "e1"), ("s2", "e1"), ("s3", "e1")],
            names=["study_id", "experiment_id"],
        ),
    )
    studies = pairs.groupby(level="study_id").mean()
    records: list[dict[str, Any]] = []
    for trait in traits:
        reason = (
            "all_tied"
            if trait.endswith("zero")
            else ("insufficient_variance" if trait.endswith("constant") else "")
        )
        records.append(
            {
                "trait": trait,
                "trait_label": trait,
                "composition_family_size": len(traits),
                "noel_family_size": len(traits),
                "composition_pairs": 3,
                "composition_studies": 3,
                "composition_profiles": 6,
                "composition_input_profiles": 6,
                "composition_unpaired_profiles_excluded": 0,
                "composition_profile_support": 0 if trait.endswith("zero") else 5,
                "composition_profile_support_fraction": 0.0 if trait.endswith("zero") else 5 / 6,
                "composition_zero_total_profiles": 0,
                "experiment_mean_delta_z": pairs[trait].mean(),
                "study_mean_delta_z": studies[trait].mean(),
                "study_delta_z_sd": studies[trait].std(ddof=1),
                "study_t_status": "untestable" if reason else "tested",
                "study_t_reason": reason,
                "study_t_primary_correction": "BY",
                "study_t_statistic": "" if reason else 0.38,
                "study_t_p": "" if reason else 0.74,
                "study_t_q_bh": "" if reason else 0.8,
                "study_t_q_by": "" if reason else 1.0,
                "study_mean_delta_z_ci95_low": "" if reason else -3.0,
                "study_mean_delta_z_ci95_high": "" if reason else 3.7,
                "noel_status": "tested",
                "noel_studies_up": 2,
                "noel_studies_down": 1,
                "noel_studies_tied": 0,
                "noel_studies_scored": 3,
                "noel_p_value": 1.0,
                "noel_q_bh_nominal": 1.0,
                "noel_q_by_nominal": 1.0,
                "noel_share_up": 2 / 3,
                "noel_ci_low": 0.1,
                "noel_ci_high": 0.99,
            }
        )
    evidence = pd.DataFrame(records)
    for effect in ("experiment_mean_delta_z", "study_mean_delta_z"):
        for suffix, values, ascending in (
            ("ascending", evidence[effect], True),
            ("descending", evidence[effect], False),
            ("absolute_descending", evidence[effect].abs(), False),
        ):
            evidence[f"{effect}_rank_{suffix}"] = values.rank(method="min", ascending=ascending)
    inventory = [FileInventoryEntry(path="src/mtp/operation.py", sha256="1" * 64, size_bytes=3)]
    artifacts = [
        ArtifactInstance.model_construct(
            product_id=name,
            sha256=("2" if name in SCOPE_PRODUCTS and scope == "phenotype" else "1") * 64,
            size_bytes=5,
            path=f"{name}.tsv",
            unit=name,
            cohort="same",
            feature_policy="same",
            coverage_policy="same",
        )
        for name in sorted(
            {
                "raw_dump",
                "kg_edges",
                "kg_nodes",
                "experiments",
                "contrast_selection",
                *SCOPE_PRODUCTS,
            }
        )
    ]
    source_root = root / scope
    source_root.mkdir(parents=True, exist_ok=True)
    (source_root / "source_identity.tsv").write_text(
        json.dumps(
            {
                "recorded_at_utc": scope,
                "declared_source_mode": "pinned",
                "pinned_verified": True,
                "consumed_files": {
                    name: {
                        "path": f"{scope}/{name}",
                        "product_id": name,
                        "sha256": "1" * 64,
                        "size_bytes": 5,
                    }
                    for name in ("raw_dump", "kg_edges", "kg_nodes")
                },
                "verified_release": {"version": "fixture"},
            }
        ),
        encoding="utf-8",
    )
    core = WorkflowRun.model_construct(
        run_id=f"core-{scope}",
        artifacts=artifacts,
        git_commit=SOURCE_COMMIT,
        git_dirty=False,
        manifest_path=str(source_root / "run.json"),
        contrast_population=POPULATION,
        direction_basis=DIRECTION,
        code_sha256="1" * 64,
        code_inventory=inventory,
        spec_sha256="1" * 64,
        definition_sha256="1" * 64,
        lock_sha256="1" * 64,
        runtime_json='{"python":"3.11"}',
        parameters=[
            ParameterValue(name="strain_rollup", value=scope, source="override"),
            ParameterValue(name="count_host_species", value="all", source="default"),
        ],
    )
    arms = {
        name: ArmRun.model_construct(
            id=f"{name}-{scope}",
            arm=name,
            git_commit=SOURCE_COMMIT,
            git_dirty=False,
            code_inventory=deepcopy(inventory),
            runtime_json='{"python":"3.11"}',
            contrast_population=POPULATION,
            direction_basis=DIRECTION,
            profile_product="ora_profiles" if name == "ilyas" else "count_profiles",
            plan_sha256="1" * 64,
            parameters=[],
        )
        for name in ("noel", "jed", "ilyas")
    }
    bases = {
        name: AnalysisBasis.model_construct(
            id=f"basis:{name}",
            profile_unit=name,
            feature_catalog_sha256=("1" if scope == "off" else "2") * 64,
            feature_filter_id=f"sha256:{('1' if scope == 'off' else '2') * 64}",
            retained_feature_count=len(traits),
            rollup_policy=f"ncbi_taxonomy_closure_rollup_{scope}",
            aggregation_policy="sum",
            denominator_policy="sum",
            background_inventory_id="same",
            weighting_policy="equal",
        )
        for name in ("count_profiles", "ora_profiles")
    }
    sources = VerifiedSources(
        ComparisonRun.model_construct(
            id=f"comparison-{scope}",
            git_commit=SOURCE_COMMIT,
            git_dirty=False,
            contrast_population=POPULATION,
            direction_basis=DIRECTION,
        ),
        core,
        arms,
        bases,
        dict.fromkeys(bases, tuple(traits)),
        scope,
    )
    identity_path = source_root / "source_identity.tsv"
    raw = identity_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    identity_artifact = next(item for item in artifacts if item.product_id == "source_identity")
    identity_artifact.sha256 = digest
    identity_artifact.size_bytes = len(raw)
    sources.recorded_hashes[identity_path.resolve()] = digest
    record = {
        "code_commit": PRODUCER,
        "code_dirty": False,
        "code_inventory": {"src/mtp/three_test_table/stage.py": "1" * 64},
        "implementation_sha256": "1" * 64,
        "methods_sha256": "1" * 64,
        "lock_sha256": "1" * 64,
        "environment": {"python": "3.11", "numpy": "2.4", "pandas": "3.0", "scipy": "1.17"},
        "strain_rollup": scope,
        "scope_role": "sensitivity" if scope == "off" else "primary",
        "contrast_population": POPULATION,
        "direction_basis": DIRECTION,
        "source_commits": dict.fromkeys(["core", "comparison", *arms], SOURCE_COMMIT),
        "analysis_bases": {name: basis.model_dump() for name, basis in bases.items()},
        "retained_traits": traits,
        "family_size": len(traits),
    }
    return CompletedTable(
        Path(f"{scope}/three_test_run.json"), Path(scope), record, sources, evidence, pairs, studies
    )


@pytest.fixture
def parents(tmp_path: Path) -> tuple[CompletedTable, CompletedTable]:
    return completed("off", tmp_path), completed("phenotype", tmp_path)


def test_same_producer_allows_rollup_catalog_and_profile_differences(
    parents: tuple[CompletedTable, CompletedTable],
) -> None:
    off, phenotype = parents
    before = phenotype.evidence.copy(deep=True)
    verify_comparable(off, phenotype)
    pd.testing.assert_frame_equal(phenotype.evidence, before)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("code_commit", "c" * 40),
        ("code_dirty", True),
        ("code_inventory", {"src/mtp/three_test_table/stage.py": "9" * 64}),
        ("implementation_sha256", "9" * 64),
        ("methods_sha256", "9" * 64),
        ("lock_sha256", "9" * 64),
        ("environment", {"python": "different"}),
        ("source_commits", {"core": "9" * 40}),
        ("contrast_population", "all_paired"),
        ("direction_basis", "group_1_vs_group_0"),
        ("scope_role", "sensitivity"),
        ("strain_rollup", "full"),
    ],
)
def test_rejects_parent_contract_mismatch(
    parents: tuple[CompletedTable, CompletedTable], field: str, value: Any
) -> None:
    off, phenotype = parents
    phenotype.record[field] = value
    with pytest.raises(ValueError, match="Completed parent"):
        verify_comparable(off, phenotype)


def test_rejects_swapped_scope_arguments(parents: tuple[CompletedTable, CompletedTable]) -> None:
    with pytest.raises(ValueError, match="off sensitivity and phenotype primary"):
        verify_comparable(*reversed(parents))


@pytest.mark.parametrize("ledger", ["pair_delta_z", "study_delta_z"])
def test_equal_counts_do_not_establish_same_cohort(
    parents: tuple[CompletedTable, CompletedTable], ledger: str
) -> None:
    off, phenotype = parents
    frame = getattr(phenotype, ledger)
    frame.rename(index={"s1": "unrelated-study"}, inplace=True)
    with pytest.raises(ValueError, match="observational keys"):
        verify_comparable(off, phenotype)


@pytest.mark.parametrize(
    "role", ["raw_dump", "kg_edges", "kg_nodes", "experiments", "contrast_selection"]
)
def test_rejects_changed_source_content_by_role(
    parents: tuple[CompletedTable, CompletedTable], role: str
) -> None:
    off, phenotype = parents
    item = next(item for item in phenotype.sources.core.artifacts or [] if item.product_id == role)
    item.sha256 = "9" * 64
    with pytest.raises(ValueError, match="Core role"):
        verify_comparable(off, phenotype)


@pytest.mark.parametrize("change", ["missing", "duplicate", "unknown", "semantic"])
def test_rejects_changed_or_ambiguous_source_roles(
    parents: tuple[CompletedTable, CompletedTable], change: str
) -> None:
    off, phenotype = parents
    artifacts = phenotype.sources.core.artifacts
    assert artifacts is not None
    if change == "missing":
        artifacts.pop()
    elif change == "duplicate":
        artifacts.append(deepcopy(artifacts[0]))
    elif change == "unknown":
        item = deepcopy(artifacts[0])
        item.product_id = "another_input"
        artifacts.append(item)
    else:
        artifacts[0].cohort = "non-human"
    with pytest.raises(ValueError, match=r"Core (source|role)"):
        verify_comparable(off, phenotype)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("git_dirty", True),
        ("code_sha256", "9" * 64),
        ("runtime_json", '{"python":"different"}'),
        ("definition_sha256", "9" * 64),
    ],
)
def test_rejects_changed_core_execution(
    parents: tuple[CompletedTable, CompletedTable], field: str, value: Any
) -> None:
    off, phenotype = parents
    setattr(phenotype.sources.core, field, value)
    with pytest.raises(ValueError, match=r"(Core execution|Upstream producer)"):
        verify_comparable(off, phenotype)


@pytest.mark.parametrize("owner", ["comparison", "noel", "jed", "ilyas"])
def test_rejects_dirty_upstream_sources(
    parents: tuple[CompletedTable, CompletedTable], owner: str
) -> None:
    off, phenotype = parents
    source = (
        phenotype.sources.comparison if owner == "comparison" else phenotype.sources.arms[owner]
    )
    source.git_dirty = True
    with pytest.raises(ValueError, match="Upstream producer must be clean"):
        verify_comparable(off, phenotype)


def test_rejects_nonrollup_core_parameter(parents: tuple[CompletedTable, CompletedTable]) -> None:
    off, phenotype = parents
    assert phenotype.sources.core.parameters is not None
    phenotype.sources.core.parameters.append(
        ParameterValue(name="extra_filter", value="yes", source="override")
    )
    with pytest.raises(ValueError, match="Core parameters"):
        verify_comparable(off, phenotype)


def test_rejects_changed_analysis_basis_semantics(
    parents: tuple[CompletedTable, CompletedTable],
) -> None:
    off, phenotype = parents
    phenotype.record["analysis_bases"]["count_profiles"]["weighting_policy"] = "unequal"
    with pytest.raises(ValueError, match="Analysis basis semantics"):
        verify_comparable(off, phenotype)


@pytest.mark.parametrize(
    ("column", "value"),
    [
        ("study_t_statistic", np.inf),
        ("study_t_statistic", ""),
        ("study_t_statistic", -0.4),
        ("study_t_p", np.nan),
        ("study_t_p", -0.1),
        ("study_t_q_by", 1.01),
        ("study_t_q_by", 0.7),
        ("study_t_q_bh", 0.1),
        ("study_t_status", "unknown"),
        ("study_t_reason", "all_tied"),
        ("study_mean_delta_z_ci95_low", 1.0),
        ("study_mean_delta_z_ci95_high", -1.0),
        ("study_delta_z_sd", -1),
        ("composition_profile_support", 7),
        ("composition_profile_support", 0.5),
        ("composition_profile_support_fraction", 0.1),
        ("composition_pairs", 2),
        ("composition_studies", 4),
        ("composition_family_size", 2),
        ("experiment_mean_delta_z_rank_descending", 1.5),
        ("study_mean_delta_z_rank_absolute_descending", 1),
        ("noel_status", "unknown"),
        ("noel_q_by_nominal", 1.2),
        ("noel_studies_scored", 4),
    ],
)
def test_rejects_invalid_persisted_evidence(tmp_path: Path, column: str, value: Any) -> None:
    parent = completed("off", tmp_path)
    parent.evidence[column] = parent.evidence[column].astype(object)
    parent.evidence.at[0, column] = value
    with pytest.raises(ValueError, match="Inherited"):
        validate_evidence(parent)


@pytest.mark.parametrize(
    "column",
    [
        "study_t_statistic",
        "study_t_p",
        "study_t_q_by",
        "study_t_q_bh",
        "study_mean_delta_z_ci95_low",
        "study_mean_delta_z_ci95_high",
    ],
)
def test_untestable_rows_cannot_claim_inference(tmp_path: Path, column: str) -> None:
    parent = completed("off", tmp_path)
    parent.evidence.at[1, column] = 0.1
    with pytest.raises(ValueError, match="Untestable inherited evidence"):
        validate_evidence(parent)


def test_untestable_status_must_match_study_vector(tmp_path: Path) -> None:
    parent = completed("off", tmp_path)
    parent.evidence.at[1, "study_t_reason"] = "insufficient_variance"
    with pytest.raises(ValueError, match="t-test status/reason"):
        validate_evidence(parent)


def test_missing_required_evidence_is_rejected(tmp_path: Path) -> None:
    parent = completed("off", tmp_path)
    with pytest.raises(ValueError, match="required unique columns"):
        validate_evidence(replace(parent, evidence=parent.evidence.drop(columns="study_t_q_by")))


def test_rejects_changed_source_release_identity(
    parents: tuple[CompletedTable, CompletedTable],
) -> None:
    off, phenotype = parents
    path = phenotype.sources.product("source_identity")
    identity = json.loads(path.read_text(encoding="utf-8"))
    identity["verified_release"]["version"] = "different"
    path.write_text(json.dumps(identity), encoding="utf-8")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    artifact = next(
        item
        for item in phenotype.sources.core.artifacts or []
        if item.product_id == "source_identity"
    )
    artifact.sha256, artifact.size_bytes = digest, path.stat().st_size
    phenotype.sources.recorded_hashes[path.resolve()] = digest
    with pytest.raises(ValueError, match="Core source declarations"):
        verify_comparable(off, phenotype)


def test_source_identity_new_read_requires_verified_hash(
    parents: tuple[CompletedTable, CompletedTable],
) -> None:
    off, phenotype = parents
    path = phenotype.sources.product("source_identity")
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match=r"source.identity.*(hash|bytes)"):
        verify_comparable(off, phenotype)


def test_source_identity_new_read_requires_verified_size(
    parents: tuple[CompletedTable, CompletedTable],
) -> None:
    off, phenotype = parents
    artifact = next(
        item
        for item in phenotype.sources.core.artifacts or []
        if item.product_id == "source_identity"
    )
    artifact.size_bytes += 1
    with pytest.raises(ValueError, match=r"source.identity.*(size|bytes)"):
        verify_comparable(off, phenotype)


def test_source_identity_added_to_publication_recheck(
    parents: tuple[CompletedTable, CompletedTable],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    off, phenotype = parents
    verify_comparable(off, phenotype)
    path = phenotype.sources.product("source_identity")
    path.write_bytes(path.read_bytes() + b"\n")
    # Parent output inventories are synthetic; exercise the reader's actual source recheck.
    monkeypatch.setattr(reader, "_check_outputs", lambda _root, _record: None)
    with pytest.raises(ValueError, match="input hash changed"):
        reader.recheck_completed_table(phenotype)


@pytest.mark.parametrize("scope", ["off", "phenotype"])
def test_source_identity_alias_and_target_are_bound(
    parents: tuple[CompletedTable, CompletedTable],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    scope: str,
) -> None:
    off, phenotype = parents
    parent = off if scope == "off" else phenotype
    target = parent.sources.product("source_identity")
    alias_root = tmp_path / f"external-{scope}"
    alias_root.mkdir()
    alias = alias_root / "identity.json"
    alias.symlink_to(target)
    artifact = next(
        item for item in parent.sources.core.artifacts or [] if item.product_id == "source_identity"
    )
    artifact.requested_path = str(alias)
    verify_comparable(off, phenotype)
    assert alias in parent.snapshots
    assert target in parent.snapshots
    assert parent.resolved_targets[alias] == target
    identical = tmp_path / f"identical-{scope}.json"
    identical.write_bytes(target.read_bytes())
    alias.unlink()
    alias.symlink_to(identical)
    monkeypatch.setattr(reader, "_check_outputs", lambda _root, _record: None)
    with pytest.raises(ValueError, match="alias target changed"):
        reader.recheck_completed_table(parent)


def test_rank_validation_preserves_adjacent_float_decimal_strings(tmp_path: Path) -> None:
    parent = completed("off", tmp_path)
    for index, value in ((1, 0.3), (2, 0.1 + 0.2)):
        trait = parent.evidence.at[index, "trait"]
        parent.pair_delta_z[trait] = value
        parent.study_delta_z[trait] = value
        parent.evidence.at[index, "study_t_reason"] = "insufficient_variance"
        parent.evidence.at[index, "study_delta_z_sd"] = parent.study_delta_z[trait].std(ddof=1)
        parent.evidence.at[index, "composition_profile_support"] = 5
        parent.evidence.at[index, "composition_profile_support_fraction"] = 5 / 6
    for effect, ledger in (
        ("experiment_mean_delta_z", parent.pair_delta_z),
        ("study_mean_delta_z", parent.study_delta_z),
    ):
        values = ledger.mean().reindex(parent.evidence.trait).reset_index(drop=True)
        for suffix, ranked, ascending in (
            ("ascending", values, True),
            ("descending", values, False),
            ("absolute_descending", values.abs(), False),
        ):
            parent.evidence[f"{effect}_rank_{suffix}"] = ranked.rank(
                method="min", ascending=ascending
            )
        parent.evidence[effect] = values.map(repr)
    assert (
        parent.evidence.at[1, "study_mean_delta_z"] != parent.evidence.at[2, "study_mean_delta_z"]
    )
    validate_evidence(parent)
