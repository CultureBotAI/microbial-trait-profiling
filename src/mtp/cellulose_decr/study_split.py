"""Training-only selection and held-out U6 overlap tests (#854).

Repeated splits describe stability; their overlapping tests are not independent
replications and are never pooled into one inferential p-value.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass, replace
from datetime import datetime
from typing import Any

import numpy as np
import pandas as pd

from mtp.cellulose_decr import selection_scan as scanner
from mtp.cellulose_decr.analysis import FECES
from mtp.datamodel.multiple_testing_pydantic import (
    HypothesisCoordinates,
    HypothesisDefinition,
    HypothesisFamily,
    RawHypothesis,
    SourceArtifact,
)
from mtp.datamodel.workflow_pydantic import AnalysisBasis
from mtp.multiple_testing import canonical_sha256, correct_family, verify_correction
from mtp.trait_guilds import Profile
from mtp.trait_screen import TESTED
from mtp.workflow.coverage import feature_filter_id
from mtp.workflow.population import ORIENTED_DIRECTION_BASIS

STATISTICS = (scanner.MEAN_DISTANCE, scanner.CORE)


@dataclass(frozen=True)
class Protocol:
    """Fixture overrides are recorded separately from the committed defaults."""

    splits: int = 20
    seed: int = 20261023
    min_studies: int = 3
    top_traits: int = 10
    permutations: int = 100000
    tolerance: float = 1e-12

    def __post_init__(self) -> None:
        if min(self.splits, self.min_studies, self.top_traits, self.permutations) < 1:
            raise ValueError("Split, eligibility, trait and draw counts must be positive")
        if self.seed < 0 or not np.isfinite(self.tolerance) or self.tolerance < 0:
            raise ValueError("Invalid seed or distance tolerance")


@dataclass(frozen=True)
class StudySplit:
    index: int
    training: tuple[str, ...]
    testing: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.index < 0 or not self.training or not self.testing:
            raise ValueError("A split needs a nonnegative index and both study partitions")
        if len(set(self.training)) != len(self.training) or len(set(self.testing)) != len(
            self.testing
        ):
            raise ValueError("Duplicate study in split membership")
        if set(self.training) & set(self.testing):
            raise ValueError("Study leakage between training and testing")


@dataclass(frozen=True)
class Candidate:
    trait: str
    conditions: tuple[str, str, str]
    training_distance: float

    def __post_init__(self) -> None:
        if not self.trait or len(self.conditions) != 3 or len(set(self.conditions)) != 3:
            raise ValueError("A candidate needs a trait and three distinct conditions")
        if any(not name.strip() for name in self.conditions):
            raise ValueError("Blank candidate condition")
        if tuple(sorted(self.conditions)) != self.conditions:
            raise ValueError("Candidate conditions must be in canonical order")
        if not np.isfinite(self.training_distance) or not 0 <= self.training_distance <= 1:
            raise ValueError("Invalid training distance")

    def identity(self, statistic: str) -> str:
        return canonical_sha256([self.trait, list(self.conditions), statistic])


def make_splits(studies: Sequence[str], *, splits: int, seed: int) -> list[StudySplit]:
    """Split sorted distinct IDs without depending on input row ordering."""
    ordered = sorted(set(studies))
    if len(ordered) < 2 or any(not study.strip() for study in ordered):
        raise ValueError("At least two nonblank study IDs are required")
    if splits < 1 or seed < 0:
        raise ValueError("Invalid split count or seed")
    result = []
    for index in range(splits):
        permutation = np.random.default_rng([seed, index, 0]).permutation(len(ordered))
        boundary = len(ordered) // 2
        result.append(
            StudySplit(
                index,
                tuple(sorted(ordered[i] for i in permutation[:boundary])),
                tuple(sorted(ordered[i] for i in permutation[boundary:])),
            )
        )
    return result


def partition_frames(
    profiles: pd.DataFrame, taxa: pd.DataFrame, studies: Sequence[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = set(studies)
    return (
        profiles.loc[profiles.study_id.isin(selected)].sort_index().copy(),
        taxa.loc[taxa.study_id.isin(selected)].sort_values(["profile_id", "taxon_id"]).copy(),
    )


def _scanner_protocol(protocol: Protocol) -> scanner.Protocol:
    """Share scanner semantics while recording this sidecar's declared parameters."""
    return replace(
        scanner.load_protocol(),
        min_studies=protocol.min_studies,
        top_n=protocol.top_traits,
        permutations=protocol.permutations,
        tolerance=protocol.tolerance,
    )


def _training_selection(
    profiles: pd.DataFrame,
    taxa: pd.DataFrame,
    matrix: pd.DataFrame,
    counts: list[Profile],
    columns: Sequence[str],
    split: StudySplit,
    protocol: Protocol,
    *,
    strain_rollup: str = "off",
) -> tuple[list[Candidate], str]:
    """Select on training studies only, without forced historical candidates."""
    if set(split.training) & set(split.testing):
        raise ValueError("Study leakage between training and testing")
    ordered_counts = sorted(
        counts,
        key=lambda p: (p.study_id, p.experiment_id, p.inc_dec),
    )
    scan_protocol = _scanner_protocol(protocol)
    scope = strain_rollup
    role = scanner.scope_role(scope, scan_protocol)
    inputs = scanner.ScanInputs(
        profiles,
        taxa,
        matrix,
        tuple(ordered_counts),
        tuple(columns),
        {},
        contrast_population=scanner.TIER1_2,
        strain_rollup=scope,
        strain_rollup_role=role,
    )
    selection = scanner.select(
        inputs,
        split.training,
        protocol=scan_protocol,
        permutations=0,
        always_included=False,
    )
    # This sidecar declares exact ties; do not inherit the helper's 1e-12 tie window.
    chosen = scanner.closest_disjoint_triples(selection.triples, tolerance=0.0)
    candidates = [
        Candidate(
            str(row["trait_id"]),
            (str(row["condition_a"]), str(row["condition_b"]), str(row["condition_c"])),
            float(row[scanner.MEAN_DISTANCE]),
        )
        for row in chosen.to_dict(orient="records")
    ]
    status = selection.status
    if status == scanner.SELECTED and not candidates:
        status = "no_testable_study_disjoint_triple"
    return candidates, status


def select_candidates(
    profiles: pd.DataFrame,
    taxa: pd.DataFrame,
    matrix: pd.DataFrame,
    counts: list[Profile],
    columns: Sequence[str],
    split: StudySplit,
    protocol: Protocol,
    *,
    strain_rollup: str = "off",
) -> list[Candidate]:
    """Return the training-only candidates; execution also records the selection status."""
    candidates, _ = _training_selection(
        profiles, taxa, matrix, counts, columns, split, protocol, strain_rollup=strain_rollup
    )
    return candidates


def _condition_reason(profiles: pd.DataFrame, candidate: Candidate, minimum: int) -> str | None:
    stool = profiles.loc[profiles.compartment.eq(FECES) & profiles.dir.eq("dec")]
    members = [
        set(stool.loc[stool.condition.eq(name), "study_id"]) for name in candidate.conditions
    ]
    if any(len(studies) < minimum for studies in members):
        return "condition_below_minimum_studies"
    if members[0] & members[1] or members[0] & members[2] or members[1] & members[2]:
        return "shared_study_between_conditions"
    return None


def evaluate_candidates(
    candidates: Sequence[Candidate],
    heldout: pd.DataFrame,
    heldout_taxa: pd.DataFrame,
    matrix: pd.DataFrame,
    protocol: Protocol,
    *,
    seed: int,
    strain_rollup: str = "off",
) -> tuple[list[RawHypothesis], list[dict[str, Any]]]:
    """Test frozen selections with only held-out studies in observations and nulls."""
    rows = []
    details = []
    identities = [candidate.identity(stat) for candidate in candidates for stat in STATISTICS]
    if len(set(identities)) != len(identities):
        raise ValueError("Duplicate training-selected hypothesis")
    for index, candidate in enumerate(candidates):
        reason = _condition_reason(heldout, candidate, protocol.min_studies)
        result: dict[str, Any] | None = None
        summary: dict[str, Any] = {}
        if reason is None:
            scan_protocol = _scanner_protocol(protocol)
            scope = strain_rollup
            role = scanner.scope_role(scope, scan_protocol)
            inputs = scanner.ScanInputs(
                heldout,
                heldout_taxa,
                matrix,
                (),
                tuple(matrix.columns),
                {},
                contrast_population=scanner.TIER1_2,
                strain_rollup=scope,
                strain_rollup_role=role,
            )
            evaluated = scanner.evaluate(
                inputs,
                inputs.studies,
                [(*candidate.conditions, candidate.trait)],
                protocol=scan_protocol,
                seed=seed + index,
            )
            result = {str(key): value for key, value in evaluated.iloc[0].to_dict().items()}
            if result["test_status"] != TESTED:
                reason = str(result["test_status"])
            summary = result
        for statistic in STATISTICS:
            test_reason = reason
            if (
                test_reason is None
                and statistic == scanner.MEAN_DISTANCE
                and result is not None
                and result["u6_draws_without_distance"] > 0
            ):
                test_reason = "undefined_distance_in_null"
            p_value = (
                None if test_reason is not None else float(summary[scanner.P_COLUMNS[statistic]])
            )
            row = RawHypothesis(
                hypothesis_id=candidate.identity(statistic),
                status="tested" if test_reason is None else "untestable",
                p_value=p_value,
                reason=test_reason,
            )
            rows.append(row)
            observed = None
            if result is not None:
                value = result[statistic]
                observed = float(value) if pd.notna(value) and np.isfinite(value) else None
            details.append(
                {
                    **row.model_dump(mode="json"),
                    "trait_id": candidate.trait,
                    "conditions": list(candidate.conditions),
                    "statistic": statistic,
                    "observed": observed,
                    "training_distance": candidate.training_distance,
                    "pool_studies": None
                    if result is None
                    else int(
                        heldout.loc[
                            heldout.compartment.eq(FECES) & heldout.dir.eq("dec"), "study_id"
                        ].nunique()
                    ),
                    "draws": 0
                    if result is None or pd.isna(result["u6_draws"])
                    else int(result["u6_draws"]),
                }
            )
    return rows, details


def declare_family(
    candidates: Sequence[Candidate],
    split: StudySplit,
    basis: AnalysisBasis,
    sources: list[SourceArtifact],
    *,
    analysis_id: str,
    strain_rollup: str,
) -> HypothesisFamily:
    """Freeze both held-out statistics before evaluation, including unavailable tests."""
    if strain_rollup not in {"off", "phenotype"} or (
        basis.rollup_policy != f"ncbi_taxonomy_closure_rollup_{strain_rollup}"
        or basis.profile_unit != "count_profiles"
    ):
        raise ValueError("Split family scope/profile representation differs from AnalysisBasis")
    members = [
        HypothesisDefinition(
            hypothesis_id=candidate.identity(statistic),
            feature_ids=[candidate.trait],
            coordinates=HypothesisCoordinates(
                condition_ids=list(candidate.conditions),
                split_id=f"split-{split.index:03}",
                stratum_id=FECES,
            ),
            test=f"U6:{statistic}",
            tail="lower" if statistic == scanner.MEAN_DISTANCE else "upper",
            eligibility_sha256=canonical_sha256(
                {
                    "selected_conditions": list(candidate.conditions),
                    "split": asdict(split),
                }
            ),
        )
        for candidate in candidates
        for statistic in STATISTICS
    ]
    return HypothesisFamily(
        family_id=f"study-split-selection-v1:split-{split.index:03}",
        analysis_id=analysis_id,
        analysis_basis_sha256=canonical_sha256(basis.model_dump(mode="json")),
        source_membership_sha256=canonical_sha256(asdict(split)),
        retained_feature_count=basis.retained_feature_count,
        strain_rollup=strain_rollup,
        feature_filter_id=basis.feature_filter_id,
        expected_family_size=len(members),
        population=f"{scanner.TIER1_2}; exact Feces; held-out studies",
        direction_basis=ORIENTED_DIRECTION_BASIS,
        profile_unit=basis.profile_unit,
        inference_unit="whole study under U6 exchangeability",
        hypothesis_unit="training-selected trait and condition triple by overlap statistic",
        estimand="overlap of trait-annotated reported taxa across condition signatures",
        selection_rule="train-only decreasing trait screen then minimum observed distance "
        "study-disjoint triple; both held-out statistics reserved",
        selection_history=[
            "Historical same-corpus exploratory inspection preceded this protocol.",
            "This split's algorithm sees only training studies before fixing its family.",
            "Held-out tests are internal assessments, not external replication.",
            "No across-split p-value pooling or best-split inference.",
        ],
        sources=sources,
        members=members,
    )


def run_splits(
    profiles: pd.DataFrame,
    taxa: pd.DataFrame,
    matrix: pd.DataFrame,
    counts: list[Profile],
    columns: Sequence[str],
    basis: AnalysisBasis,
    sources: list[SourceArtifact],
    protocol: Protocol,
    *,
    analysis_id: str,
    strain_rollup: str,
    recorded_at: datetime,
) -> list[dict[str, Any]]:
    """Run every declared split, fixing its full family before held-out evaluation."""
    study_ids = set(profiles.study_id)
    if not set(taxa.study_id) <= study_ids or not {p.study_id for p in counts} <= study_ids:
        raise ValueError("Taxa or screen profiles contain studies outside the core population")
    if basis.retained_feature_count != len(columns) or len(set(columns)) != len(columns):
        raise ValueError("AnalysisBasis and screening candidate inventory disagree")
    if basis.feature_filter_id != feature_filter_id(columns) or tuple(matrix.columns) != tuple(
        columns
    ):
        raise ValueError("AnalysisBasis ordered features and matrix inventory disagree")
    results = []
    for split in make_splits(sorted(study_ids), splits=protocol.splits, seed=protocol.seed):
        candidates, selection_status = _training_selection(
            profiles, taxa, matrix, counts, columns, split, protocol, strain_rollup=strain_rollup
        )
        family = declare_family(
            candidates,
            split,
            basis,
            sources,
            analysis_id=analysis_id,
            strain_rollup=strain_rollup,
        )
        heldout, heldout_taxa = partition_frames(profiles, taxa, split.testing)
        null_seed = int(np.random.default_rng([protocol.seed, split.index, 1]).integers(0, 2**31))
        raw, details = evaluate_candidates(
            candidates,
            heldout,
            heldout_taxa,
            matrix,
            protocol,
            seed=null_seed,
            strain_rollup=strain_rollup,
        )
        correction = correct_family(family, raw, recorded_at=recorded_at)
        verify_correction(correction, family=family, hypotheses=raw)
        results.append(
            {
                "split": asdict(split),
                "selection_status": selection_status,
                "selected_candidates": [asdict(candidate) for candidate in candidates],
                "null_seed": null_seed,
                "correction": correction.model_dump(mode="json"),
                "test_details": details,
            }
        )
    return results


def summarize_stability(records: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Descriptive counts across overlapping splits, with no aggregate inference."""
    inventory: dict[str, dict[str, Any]] = {}
    for record in records:
        details = {row["hypothesis_id"]: row for row in record["test_details"]}
        for result in record["correction"]["results"]:
            identity = result["hypothesis_id"]
            detail = details[identity]
            if identity not in inventory:
                inventory[identity] = {
                    "hypothesis_id": identity,
                    "trait_id": detail["trait_id"],
                    "conditions": detail["conditions"],
                    "statistic": detail["statistic"],
                    "total_splits": len(records),
                    "selected_splits": 0,
                    "testable_splits": 0,
                    "by_rejections_005": 0,
                    "bh_rejections_005": 0,
                    "interpretation": "descriptive counts; overlapping splits are dependent",
                }
            row = inventory[identity]
            row["selected_splits"] += 1
            row["testable_splits"] += int(result["status"] == "tested")
            row["by_rejections_005"] += int(result["reject_by_005"] is True)
            row["bh_rejections_005"] += int(result["reject_bh_005"] is True)
    return [inventory[key] for key in sorted(inventory)]
