"""Require comparable completed parents and well-formed inherited inference (#1041).

This layer consumes the pinned read-only reader's verified bundles. It does not
recompute delta-z, tests, confidence intervals or multiplicity corrections. The
checks below establish common inputs/producer and validate the inherited fields'
contracts, including descriptive ledger consistency and full-family effect ranks.
"""

from __future__ import annotations

import hashlib
import json
import math
from functools import partial
from pathlib import Path
from typing import Any

import numpy as np

from mtp.datamodel.workflow_pydantic import ArtifactInstance, ParameterValue
from mtp.three_test_table.reader import CompletedTable, _remember_path

_SCOPE_PRODUCTS = frozenset(
    {
        "count_coverage",
        "count_profiles",
        "feature_catalog",
        "feature_metadata",
        "ora_coverage",
        "ora_profiles",
        "taxa_matrix",
    }
)
_BASIS_SCOPE_FIELDS = frozenset(
    {
        "feature_catalog_sha256",
        "feature_filter_id",
        "retained_feature_count",
        "rollup_policy",
    }
)
_T_INFERENCE = (
    "study_t_statistic",
    "study_t_p",
    "study_t_q_bh",
    "study_t_q_by",
    "study_mean_delta_z_ci95_low",
    "study_mean_delta_z_ci95_high",
)
_NOEL_INFERENCE = (
    "noel_p_value",
    "noel_q_bh_nominal",
    "noel_q_by_nominal",
    "noel_share_up",
    "noel_ci_low",
    "noel_ci_high",
)
_EFFECTS = ("experiment_mean_delta_z", "study_mean_delta_z")
_RANKS = ("ascending", "descending", "absolute_descending")


def _number(value: Any, field: str) -> float:
    if isinstance(value, (bool, np.bool_)) or value == "":
        raise ValueError(f"Inherited evidence requires a finite number: {field}")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Inherited evidence requires a finite number: {field}") from error
    if not math.isfinite(number):
        raise ValueError(f"Inherited evidence requires a finite number: {field}")
    return number


def _count(value: Any, field: str) -> int:
    number = _number(value, field)
    if number < 0 or not number.is_integer():
        raise ValueError(f"Inherited evidence requires a nonnegative integer: {field}")
    return int(number)


def _close(actual: float, expected: float, field: str) -> None:
    if not math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-14):
        raise ValueError(f"Inherited evidence differs from its descriptive ledger: {field}")


def _probabilities(row: dict[str, Any], fields: tuple[str, str, str]) -> None:
    p, bh, by = (_number(row[name], name) for name in fields)
    if not all(0 <= value <= 1 for value in (p, bh, by)):
        raise ValueError("Inherited probabilities must be between zero and one")
    # Permit only serialization-scale roundoff, never a new correction or family.
    if bh + 1e-14 < p or by + 1e-14 < bh:
        raise ValueError("Inherited BY/BH probabilities violate their ordering")


def _blank(row: dict[str, Any], columns: tuple[str, ...]) -> None:
    if any(row[name] != "" for name in columns):
        raise ValueError("Untestable inherited evidence must leave inference fields blank")


def _t_fields(row: dict[str, Any], values: np.ndarray[Any, np.dtype[np.float64]]) -> None:
    count = len(values)
    sd = float(np.std(values, ddof=1)) if count >= 2 else None
    reason = (
        "insufficient_support"
        if count < 2
        else ("all_tied" if np.all(values == 0.0) else "insufficient_variance")
        if sd is not None and sd <= 1e-12
        else ""
    )
    if (
        row["study_t_status"] != ("untestable" if reason else "tested")
        or row["study_t_reason"] != reason
    ):
        raise ValueError("Inherited t-test status/reason differs from its study support")
    if row["study_t_primary_correction"] != "BY":
        raise ValueError("Inherited t-test correction differs")
    if sd is None:
        if row["study_delta_z_sd"] != "":
            raise ValueError("Insufficient-support study SD must be blank")
    else:
        _close(_number(row["study_delta_z_sd"], "study_delta_z_sd"), sd, "study_delta_z_sd")
    if reason:
        _blank(row, _T_INFERENCE)
        return
    _probabilities(row, ("study_t_p", "study_t_q_bh", "study_t_q_by"))
    statistic = _number(row["study_t_statistic"], "study_t_statistic")
    mean = _number(row["study_mean_delta_z"], "study_mean_delta_z")
    low = _number(row["study_mean_delta_z_ci95_low"], "study CI low")
    high = _number(row["study_mean_delta_z_ci95_high"], "study CI high")
    if not low <= mean <= high or np.sign(statistic) != np.sign(mean):
        raise ValueError("Inherited t statistic/CI does not agree with the study mean")


def _noel_fields(row: dict[str, Any], study_count: int) -> None:
    up, down, tied, scored = (
        _count(row[f"noel_studies_{name}"], f"noel_studies_{name}")
        for name in ("up", "down", "tied", "scored")
    )
    if up + down != scored or scored + tied > study_count:
        raise ValueError("Inherited Noel vote counts are inconsistent")
    status = "tested" if scored else "all_tied" if tied else "no_eligible_studies"
    if row["noel_status"] != status:
        raise ValueError("Inherited Noel status differs from its vote support")
    if status != "tested":
        _blank(row, _NOEL_INFERENCE)
        return
    _probabilities(row, ("noel_p_value", "noel_q_bh_nominal", "noel_q_by_nominal"))
    share, low, high = (_number(row[name], name) for name in _NOEL_INFERENCE[3:])
    if not 0 <= low <= share <= high <= 1:
        raise ValueError("Inherited Noel proportion/CI is invalid")


def validate_evidence(parent: CompletedTable) -> None:
    """Validate inherited summary fields without computing any inferential replacement."""
    evidence, pairs, studies = parent.evidence, parent.pair_delta_z, parent.study_delta_z
    traits = parent.record["retained_traits"]
    required = {
        "trait",
        "trait_label",
        "composition_family_size",
        "noel_family_size",
        "composition_pairs",
        "composition_studies",
        "composition_profiles",
        "composition_input_profiles",
        "composition_unpaired_profiles_excluded",
        "composition_profile_support",
        "composition_profile_support_fraction",
        "composition_zero_total_profiles",
        "study_delta_z_sd",
        "study_t_status",
        "study_t_reason",
        "study_t_primary_correction",
        "noel_status",
        *_T_INFERENCE,
        *_NOEL_INFERENCE,
        *_EFFECTS,
        *(f"noel_studies_{name}" for name in ("up", "down", "tied", "scored")),
        *(f"{effect}_rank_{suffix}" for effect in _EFFECTS for suffix in _RANKS),
    }
    if not required.issubset(evidence.columns) or not evidence.columns.is_unique:
        raise ValueError("Inherited evidence lacks required unique columns")
    if (
        evidence.trait.tolist() != traits
        or len(traits) != len(set(traits))
        or not traits
        or parent.record["family_size"] != len(traits)
        or list(pairs.columns) != traits
        or list(studies.columns) != traits
        or not pairs.index.is_unique
        or not studies.index.is_unique
    ):
        raise ValueError("Inherited evidence/ledgers differ from the full retained catalog")
    if (
        not np.isfinite(pairs.to_numpy(dtype=float)).all()
        or not np.isfinite(studies.to_numpy(dtype=float)).all()
    ):
        raise ValueError("Inherited descriptive ledgers must be finite")
    for record in evidence.to_dict(orient="records"):
        row = {str(name): value for name, value in record.items()}
        if not isinstance(row["trait_label"], str) or not row["trait_label"].strip():
            raise ValueError("Inherited trait label is empty")
        expected = {
            "composition_family_size": len(traits),
            "noel_family_size": len(traits),
            "composition_pairs": len(pairs),
            "composition_studies": len(studies),
            "composition_profiles": 2 * len(pairs),
        }
        for field, count in expected.items():
            if _count(row[field], field) != count:
                raise ValueError(f"Inherited evidence count differs: {field}")
        profiles = expected["composition_profiles"]
        support = _count(row["composition_profile_support"], "profile support")
        zero = _count(row["composition_zero_total_profiles"], "zero-total profiles")
        unpaired = _count(row["composition_unpaired_profiles_excluded"], "unpaired profiles")
        if (
            support > profiles
            or zero > profiles
            or _count(row["composition_input_profiles"], "input profiles") != profiles + unpaired
        ):
            raise ValueError("Inherited profile support/counts are inconsistent")
        if profiles:
            _close(
                _number(row["composition_profile_support_fraction"], "profile fraction"),
                support / profiles,
                "profile support fraction",
            )
        trait = row["trait"]
        _close(
            _number(row["experiment_mean_delta_z"], "pair mean"),
            float(pairs[trait].mean()),
            "pair mean",
        )
        _close(
            _number(row["study_mean_delta_z"], "study mean"),
            float(studies[trait].mean()),
            "study mean",
        )
        _t_fields(row, studies[trait].to_numpy(dtype=np.float64))
        _noel_fields(row, len(studies))
    for effect in _EFFECTS:
        values = evidence[effect].map(partial(_number, field=effect))
        for suffix in _RANKS:
            ranked = values.abs() if suffix == "absolute_descending" else values
            expected_ranks = ranked.rank(method="min", ascending=suffix == "ascending")
            rank_column = f"{effect}_rank_{suffix}"
            actual = evidence[rank_column].map(partial(_number, field=rank_column))
            if not np.array_equal(actual.to_numpy(), expected_ranks.to_numpy()):
                raise ValueError("Inherited full-family ranks differ from their effect estimates")


def _parameters(items: list[ParameterValue] | None, *, scope: str | None = None) -> dict[str, str]:
    values = {item.name: item.value for item in items or []}
    if len(values) != len(items or []):
        raise ValueError("Duplicate execution parameter names")
    if scope is not None and values.pop("strain_rollup", None) != scope:
        raise ValueError("Core strain-rollup parameter differs")
    return values


def _equal_fields(first: Any, second: Any, names: tuple[str, ...], description: str) -> None:
    for name in names:
        if getattr(first, name) != getattr(second, name):
            raise ValueError(f"{description} differs: {name}")


def _artifact_roles(parent: CompletedTable) -> dict[str, ArtifactInstance]:
    artifacts = parent.sources.core.artifacts or []
    roles = {artifact.product_id: artifact for artifact in artifacts}
    required = {
        "raw_dump",
        "kg_edges",
        "kg_nodes",
        "experiments",
        "contrast_selection",
        "source_identity",
    }
    if len(roles) != len(artifacts) or not required.issubset(roles):
        raise ValueError("Core source roles are missing or ambiguous")
    return roles


def _source_identity(parent: CompletedTable) -> dict[str, Any]:
    path = parent.sources.product("source_identity")
    target = path.resolve(strict=True)
    artifact = _artifact_roles(parent)["source_identity"]
    raw = path.read_bytes()
    snapshot = {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}
    if snapshot != {"sha256": artifact.sha256, "size_bytes": artifact.size_bytes} or snapshot[
        "sha256"
    ] != parent.sources.expected_digest(path):
        raise ValueError("Core source-identity bytes differ from recorded hash/size")
    aliases = {path}
    for value in (artifact.path, artifact.requested_path):
        if value:
            alias = Path(value)
            aliases.add(
                alias
                if alias.is_absolute()
                else Path(parent.sources.core.manifest_path).parent / alias
            )
    for alias in aliases:
        if alias.resolve(strict=True) != target:
            raise ValueError("Core source-identity alias target differs from the bound artifact")
        _remember_path(alias, snapshot, parent.snapshots, parent.resolved_targets)
    identity: Any = json.loads(raw)
    if not isinstance(identity, dict) or not isinstance(identity.get("consumed_files"), dict):
        raise ValueError("Core source-identity declaration is malformed")
    identity.pop("recorded_at_utc", None)
    roles = _artifact_roles(parent)
    files = identity["consumed_files"]
    if set(files) != {"raw_dump", "kg_edges", "kg_nodes"}:
        raise ValueError("Core source-identity input roles differ")
    for name, item in files.items():
        if not isinstance(item, dict):
            raise ValueError("Core source-identity input is malformed")
        item.pop("path", None)
        artifact = roles[name]
        if item != {
            "product_id": name,
            "sha256": artifact.sha256,
            "size_bytes": artifact.size_bytes,
        }:
            raise ValueError("Core source identity differs from its recorded artifact role")
    return identity


def verify_comparable(off: CompletedTable, phenotype: CompletedTable) -> None:
    """Require isolated rollup scope differences between two already verified parents."""
    for parent, scope, role in ((off, "off", "sensitivity"), (phenotype, "phenotype", "primary")):
        if (parent.record["strain_rollup"], parent.record["scope_role"], parent.sources.scope) != (
            scope,
            role,
            scope,
        ):
            raise ValueError("Completed parents must be off sensitivity and phenotype primary")
        if parent.record["code_dirty"] is not False:
            raise ValueError("Completed parent producer must be clean")
        sources = parent.sources
        if (
            sources.core.git_dirty is not False
            or sources.comparison.git_dirty is not False
            or any(arm.git_dirty is not False for arm in sources.arms.values())
        ):
            raise ValueError("Upstream producer must be clean for isolated scope comparison")
        validate_evidence(parent)
    for name in (
        "code_commit",
        "code_inventory",
        "implementation_sha256",
        "methods_sha256",
        "lock_sha256",
        "environment",
        "source_commits",
        "contrast_population",
        "direction_basis",
    ):
        if off.record[name] != phenotype.record[name]:
            raise ValueError(f"Completed parent contract differs: {name}")
    for ledger in ("pair_delta_z", "study_delta_z"):
        first, second = getattr(off, ledger), getattr(phenotype, ledger)
        if first.index.names != second.index.names or set(first.index) != set(second.index):
            raise ValueError(f"Completed parents have different observational keys: {ledger}")
    first_sources, second_sources = off.sources, phenotype.sources
    _equal_fields(
        first_sources.core,
        second_sources.core,
        (
            "git_commit",
            "code_sha256",
            "code_inventory",
            "spec_sha256",
            "definition_sha256",
            "lock_sha256",
            "runtime_json",
            "contrast_population",
            "direction_basis",
        ),
        "Core execution",
    )
    if _parameters(first_sources.core.parameters, scope="off") != _parameters(
        second_sources.core.parameters, scope="phenotype"
    ):
        raise ValueError("Core parameters differ beyond strain rollup")
    _equal_fields(
        first_sources.comparison,
        second_sources.comparison,
        (
            "git_commit",
            "contrast_population",
            "direction_basis",
        ),
        "Comparison execution",
    )
    if set(first_sources.arms) != set(second_sources.arms) or set(first_sources.arms) != {
        "noel",
        "jed",
        "ilyas",
    }:
        raise ValueError("Upstream arm identities differ")
    for name, arm in first_sources.arms.items():
        counterpart = second_sources.arms[name]
        _equal_fields(
            arm,
            counterpart,
            (
                "git_commit",
                "code_inventory",
                "runtime_json",
                "profile_product",
                "contrast_population",
                "direction_basis",
            ),
            f"{name} execution",
        )
        if _parameters(arm.parameters) != _parameters(counterpart.parameters):
            raise ValueError(f"{name} parameters differ")
    first_roles, second_roles = _artifact_roles(off), _artifact_roles(phenotype)
    if set(first_roles) != set(second_roles):
        raise ValueError("Core source role inventories differ")
    for name, artifact in first_roles.items():
        counterpart_artifact = second_roles[name]
        _equal_fields(
            artifact,
            counterpart_artifact,
            ("unit", "cohort", "feature_policy", "coverage_policy"),
            f"Core role {name}",
        )
        if name not in _SCOPE_PRODUCTS and name != "source_identity":
            _equal_fields(
                artifact, counterpart_artifact, ("sha256", "size_bytes"), f"Core role {name}"
            )
    if _source_identity(off) != _source_identity(phenotype):
        raise ValueError("Core source declarations differ beyond timestamps and locations")
    first_bases, second_bases = off.record["analysis_bases"], phenotype.record["analysis_bases"]
    if set(first_bases) != set(second_bases) or set(first_bases) != {
        "count_profiles",
        "ora_profiles",
    }:
        raise ValueError("Analysis basis identities differ")
    for name, basis in first_bases.items():
        counterpart = second_bases[name]
        if {key: value for key, value in basis.items() if key not in _BASIS_SCOPE_FIELDS} != {
            key: value for key, value in counterpart.items() if key not in _BASIS_SCOPE_FIELDS
        }:
            raise ValueError(f"Analysis basis semantics differ: {name}")
