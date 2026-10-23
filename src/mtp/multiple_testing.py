"""Issue #852: additive full-family correction, independent of native arm calls.

Adapters must derive the declared inventory and source identities from verified
execution before supplying p-values. The engine checks that declaration exactly;
it cannot establish study independence or repair post-selection inference.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import scipy
import yaml
from numpy.typing import NDArray
from scipy.stats import false_discovery_control

from mtp.datamodel import multiple_testing_pydantic as models
from mtp.datamodel.multiple_testing_pydantic import (
    AdjustedHypothesis,
    CorrectionRun,
    HypothesisDefinition,
    HypothesisFamily,
    HypothesisFamilyMetadata,
    MultipleTestingPolicy,
    RawHypothesis,
)
from mtp.workflow.runtime import content_hash

# Public compatibility name for adapters; use the shared canonical hashing implementation.
canonical_sha256 = content_hash

ENGINE_VERSION = "1.1.0"
POLICY_PATH = Path(__file__).resolve().parents[2] / "methods" / "multiple_testing_policy_v1.yaml"
# Pin the authoritative versioned artifact, without duplicating its policy values.
_POLICY_FILE_SHA256 = "d1d504783871b1847fb91bc2cb96cea707400708180e7b75ca98ceba1b92b1ba"


def load_policy(path: Path = POLICY_PATH) -> MultipleTestingPolicy:
    """Read and verify the authoritative policy; v1 values cannot be overridden."""
    policy = MultipleTestingPolicy.model_validate(yaml.safe_load(path.read_text(encoding="utf-8")))
    _check_policy(policy)
    return policy


def _check_policy(policy: MultipleTestingPolicy) -> None:
    policy_bytes = POLICY_PATH.read_bytes()
    if hashlib.sha256(policy_bytes).hexdigest() != _POLICY_FILE_SHA256:
        raise ValueError("unsupported or changed multiple-testing policy file")
    supported = MultipleTestingPolicy.model_validate(yaml.safe_load(policy_bytes))
    if policy.model_dump(mode="json") != supported.model_dump(mode="json"):
        raise ValueError("unsupported or changed multiple-testing policy")


def _code_sha256() -> str:
    return canonical_sha256(
        {
            "engine": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "models": hashlib.sha256(Path(models.__file__).read_bytes()).hexdigest(),
        }
    )


def _unique(values: Sequence[str], label: str) -> None:
    if any(not value.strip() for value in values):
        raise ValueError(f"blank {label}")
    if len(values) != len(set(values)):
        raise ValueError(f"duplicate {label}")


def scientific_hypothesis_sha256(member: HypothesisDefinition) -> str:
    """Identity of a scientific test; row aliases and eligibility hashes cannot split it."""
    return canonical_sha256(
        member.model_dump(mode="json", exclude={"hypothesis_id", "eligibility_sha256"})
    )


def validate_family_metadata(metadata: HypothesisFamilyMetadata) -> HypothesisFamilyMetadata:
    """Validate scientific scope without materializing a large member inventory."""
    declared = HypothesisFamilyMetadata.model_validate(metadata.model_dump(mode="json"))
    for key in (
        "family_id",
        "analysis_id",
        "population",
        "direction_basis",
        "profile_unit",
        "inference_unit",
        "hypothesis_unit",
        "estimand",
        "selection_rule",
    ):
        if not str(getattr(declared, key)).strip():
            raise ValueError(f"blank family metadata: {key}")
    if not declared.selection_history or any(not s.strip() for s in declared.selection_history):
        raise ValueError("selection history must be explicitly disclosed")
    if not declared.sources:
        raise ValueError("family must identify its upstream sources")
    _unique([s.identity for s in declared.sources], "source identity")
    if declared.strain_rollup not in {"off", "phenotype"}:
        raise ValueError("unsupported strain rollup scope")
    if re.fullmatch(r"sha256:[0-9a-f]{64}", declared.feature_filter_id) is None:
        raise ValueError("feature filter identity must be a complete SHA256 identifier")
    digests = [declared.analysis_basis_sha256, declared.source_membership_sha256]
    digests.extend(source.sha256 for source in declared.sources)
    if any(re.fullmatch(r"[0-9a-f]{64}", digest) is None for digest in digests):
        raise ValueError("family identity must be a complete lowercase SHA256 digest")
    declared.sources.sort(key=lambda s: s.identity)
    return declared


def validate_hypothesis_definition(member: HypothesisDefinition) -> HypothesisDefinition:
    """Validate one declared coordinate without changing its identity or order."""
    member = HypothesisDefinition.model_validate(member.model_dump(mode="json"))
    if not member.hypothesis_id.strip():
        raise ValueError("blank declared hypothesis identity")
    if not member.feature_ids or any(not f.strip() for f in member.feature_ids):
        raise ValueError("hypothesis feature identities must be nonempty")
    if not member.test.strip() or not member.tail.strip():
        raise ValueError("hypothesis test and tail must be nonempty")
    coordinates = member.coordinates
    for key in ("profile_id", "stratum_id", "split_id"):
        value = getattr(coordinates, key)
        if value is not None and not value.strip():
            raise ValueError(f"blank hypothesis coordinate: {key}")
    conditions = coordinates.condition_ids or []
    _unique(conditions, "condition coordinate")
    coordinates.condition_ids = sorted(conditions)
    if (
        member.eligibility_sha256 is not None
        and re.fullmatch(r"[0-9a-f]{64}", member.eligibility_sha256) is None
    ):
        raise ValueError("hypothesis eligibility must be a complete SHA256 digest")
    return member


def _family(family: HypothesisFamily) -> HypothesisFamily:
    declared = HypothesisFamily.model_validate(family.model_dump(mode="json"))
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    metadata = validate_family_metadata(metadata)
    _unique([m.hypothesis_id for m in declared.members], "declared hypothesis identity")
    if declared.expected_family_size != len(declared.members):
        raise ValueError("expected family size differs from declared inventory")
    members = [validate_hypothesis_definition(member) for member in declared.members]
    _unique([scientific_hypothesis_sha256(member) for member in members], "scientific hypothesis")
    return HypothesisFamily(
        **metadata.model_dump(), members=sorted(members, key=lambda m: m.hypothesis_id)
    )


def validate_raw_hypothesis(row: RawHypothesis) -> RawHypothesis:
    """Validate one raw test using the same rules in memory and disk-backed inventories."""
    row = RawHypothesis.model_validate(row.model_dump(mode="json"))
    if not row.hypothesis_id.strip():
        raise ValueError("blank raw hypothesis identity")
    if row.status not in {"tested", "untestable"}:
        raise ValueError("unsupported hypothesis test status")
    if row.status == "tested":
        if (
            row.p_value is None
            or not math.isfinite(row.p_value)
            or not 0 <= row.p_value <= 1
            or math.copysign(1, row.p_value) < 0
            or row.reason is not None
        ):
            raise ValueError("tested hypothesis needs a finite probability and no undefined reason")
    elif row.p_value is not None or row.reason is None or not row.reason.strip():
        raise ValueError("untestable hypothesis needs null p-value and an explicit reason")
    return row


def _raw(family: HypothesisFamily, hypotheses: Sequence[RawHypothesis]) -> list[RawHypothesis]:
    rows = [validate_raw_hypothesis(row) for row in hypotheses]
    _unique([row.hypothesis_id for row in rows], "raw hypothesis identity")
    if {row.hypothesis_id for row in rows} != {m.hypothesis_id for m in family.members}:
        raise ValueError("raw hypothesis membership differs from declared inventory")
    return sorted(rows, key=lambda row: row.hypothesis_id)


def adjusted_probabilities(
    values: Sequence[float] | NDArray[np.float64],
    *,
    policy: MultipleTestingPolicy,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Use one full-family numeric vector, without allocating per-member Python models."""
    probabilities = np.asarray(values, dtype=np.float64)
    _check_policy(policy)
    if (
        probabilities.ndim != 1
        or not np.all(np.isfinite(probabilities))
        or np.any((probabilities < 0) | (probabilities > 1))
        or np.any(np.signbit(probabilities))
    ):
        raise ValueError("probabilities must be a finite one-dimensional vector in [0, 1]")
    if not probabilities.size:
        return probabilities.copy(), probabilities.copy()
    return np.asarray(
        false_discovery_control(probabilities, method=policy.primary_method)
    ), np.asarray(false_discovery_control(probabilities, method=policy.comparison_method))


def adjusted_hypothesis(
    row: RawHypothesis, q_by: float, q_bh: float, *, policy: MultipleTestingPolicy
) -> AdjustedHypothesis:
    """Export one member using the policy validated by the correction entry point."""
    nominal, stringent = policy.thresholds
    result = AdjustedHypothesis(**row.model_dump())
    if row.status == "tested":
        result.q_by, result.q_bh = float(q_by), float(q_bh)
        result.reject_by_005 = q_by <= nominal
        result.reject_by_0001 = q_by <= stringent
        result.reject_bh_005 = q_bh <= nominal
        result.reject_bh_0001 = q_bh <= stringent
    return result


def _results(rows: list[RawHypothesis], policy: MultipleTestingPolicy) -> list[AdjustedHypothesis]:
    if not rows:
        return []
    p_values = [row.p_value if row.p_value is not None else 1.0 for row in rows]
    by, bh = adjusted_probabilities(p_values, policy=policy)
    return [
        adjusted_hypothesis(row, q_by, q_bh, policy=policy)
        for row, q_by, q_bh in zip(rows, by, bh, strict=True)
    ]


def correct_family(
    family: HypothesisFamily,
    hypotheses: Sequence[RawHypothesis],
    *,
    policy: MultipleTestingPolicy | None = None,
    recorded_at: datetime | None = None,
) -> CorrectionRun:
    """Correct the exact declared family, including reserved untestable members.

    The supplied expected family size is verified against both inventories. It
    is not assumed equal to retained_feature_count (e.g. ratios or joint ORA
    tails). Source byte/AnalysisBasis verification belongs to the adapter.
    """
    active_policy = load_policy() if policy is None else policy.model_copy(deep=True)
    _check_policy(active_policy)
    declared = _family(family)
    rows = _raw(declared, hypotheses)
    results = _results(rows, active_policy)
    timestamp = recorded_at if recorded_at is not None else datetime.now(UTC)
    offset = timestamp.utcoffset()
    if offset is None or offset.total_seconds() != 0:
        raise ValueError("recorded_at must be timezone-aware UTC")
    envelope = {
        "recorded_at": timestamp.astimezone(UTC).isoformat(),
        "engine_version": ENGINE_VERSION,
        "scipy_version": scipy.__version__,
        "numpy_version": np.__version__,
        "code_sha256": _code_sha256(),
        "policy_sha256": canonical_sha256(active_policy.model_dump(mode="json")),
        "family_sha256": canonical_sha256(declared.model_dump(mode="json")),
        "input_sha256": canonical_sha256([r.model_dump(mode="json") for r in rows]),
        "output_sha256": canonical_sha256([r.model_dump(mode="json") for r in results]),
    }
    return CorrectionRun(
        policy=active_policy,
        family=declared,
        hypotheses=rows,
        results=results,
        run_sha256=canonical_sha256(envelope),
        **envelope,
    )


def verify_correction(
    run: CorrectionRun,
    *,
    family: HypothesisFamily,
    hypotheses: Sequence[RawHypothesis],
    policy: MultipleTestingPolicy | None = None,
) -> None:
    """Replay against externally trusted inputs; self-reported hashes do not suffice.

    Verification requires the same recorded engine/model, NumPy and SciPy versions.
    Preserve that execution environment for archival verification after upgrades.
    No adjusted p-value, rejection, source identity or scientific metadata is
    trusted merely because its envelope is internally consistent.
    """
    actual = CorrectionRun.model_validate(run.model_dump(mode="json"))
    _check_policy(actual.policy)
    expected = correct_family(family, hypotheses, policy=policy, recorded_at=actual.recorded_at)
    actual.family = _family(actual.family)
    actual.hypotheses = _raw(actual.family, actual.hypotheses)
    _unique([r.hypothesis_id for r in actual.results], "adjusted hypothesis identity")
    actual.results.sort(key=lambda row: row.hypothesis_id)
    recorded_output = canonical_sha256([row.model_dump(mode="json") for row in actual.results])
    if recorded_output != actual.output_sha256:
        raise ValueError(
            "correction verification failed: output_sha256 differs from recorded results"
        )
    for field in CorrectionRun.model_fields:
        if getattr(actual, field) != getattr(expected, field):
            raise ValueError(f"correction verification failed: {field} differs from replay")
