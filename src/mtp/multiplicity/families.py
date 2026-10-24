"""Declared full-screen and ordered-ratio families from verified native exports."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from collections.abc import Callable, Iterator
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

from mtp.datamodel.multiple_testing_pydantic import (
    HypothesisCoordinates,
    HypothesisDefinition,
    HypothesisFamilyMetadata,
    RawHypothesis,
    ReasonVocabulary,
)
from mtp.multiple_testing import canonical_sha256
from mtp.multiplicity.sources import VerifiedSources, recorded_arm_output
from mtp.paths import file_sha256
from mtp.table1_trait_maps import study_trait_contrasts, tally_direction_signs
from mtp.trait_composition import StudyLogRatio, log_ratio_support
from mtp.trait_guilds import TABLE1_TRAITS, drop_zero_total_pairs, load_profiles
from mtp.trait_screen import candidate_slots


@dataclass(frozen=True)
class NativeObservation:
    """Original p/q and row contents retained independently of sidecar decisions."""

    raw: RawHypothesis
    native_q: float | None
    native_method: str
    native_family: str
    native_threshold: float
    row: dict[str, str]


@dataclass(frozen=True)
class FamilyInput:
    metadata: HypothesisFamilyMetadata
    declarations: Callable[[], Iterator[HypothesisDefinition]]
    observations: Callable[[], Iterator[NativeObservation]]
    source_paths: tuple[Path, ...]


def rows(path: Path, expected_sha256: str | None = None) -> Iterator[dict[str, str]]:
    """Hash exactly the consumed bytes, including comments and original line endings."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:

        def lines() -> Iterator[str]:
            for raw_line in handle:
                digest.update(raw_line)
                line = raw_line.decode("utf-8")
                if not line.startswith("#"):
                    yield line

        reader = csv.DictReader(lines(), delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Missing or duplicate source columns")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed source table row")
            yield row
    if expected_sha256 is not None and digest.hexdigest() != expected_sha256:
        raise ValueError(f"Source consumed bytes differ from upstream recorded hash: {path.name}")


def require_fields(row: dict[str, Any], expected: dict[str, Any]) -> None:
    for key, value in expected.items():
        if row.get(key) != value:
            raise ValueError(f"Unsupported or mismatched source metadata: {key}")


def probability(value: str) -> float | None:
    if not value:
        return None
    number = float(value)
    if not math.isfinite(number) or not 0 <= number <= 1 or math.copysign(1, number) < 0:
        raise ValueError("Source probability must be finite and in [0,1]")
    return number


def support_count(row: dict[str, str], key: str) -> int:
    """Read the producer's nonnegative integer support, without rounding or coercion."""
    value = int(row[key])
    if value < 0:
        raise ValueError(f"Invalid producer support: {key}")
    return value


def validate_sign_status(row: dict[str, str], *, status: str, empty: str) -> None:
    """Check a sign-test reason against exported tallies, without recomputing p."""
    up, down, tied = (support_count(row, f"studies_{kind}") for kind in ("up", "down", "tied"))
    if support_count(row, "studies_scored") != up + down:
        raise ValueError("Producer sign-test support counts disagree")
    expected = "tested" if up + down else "all_tied" if tied else empty
    if row[status] != expected:
        raise ValueError("Producer sign-test status/reason disagrees with support")
    if "studies_eligible" in row:
        eligible = support_count(row, "studies_eligible")
        if eligible != up + down + tied:
            raise ValueError("Producer eligible-study support disagrees with sign tallies")
        if support_count(row, "studies_total") != eligible + support_count(row, "studies_excluded"):
            raise ValueError("Producer total-study support disagrees with exclusions")


def observation(
    identity: str,
    row: dict[str, str],
    *,
    p: str,
    q: str,
    status: str,
    undefined: set[str],
    method: str,
    native_family: str,
    threshold: float = 0.05,
    untestable_q: float | None = None,
) -> NativeObservation:
    test_status = row[status]
    if test_status != "tested" and test_status not in undefined:
        raise ValueError(f"Unsupported native test status: {test_status}")
    value = probability(row[p])
    if (value is not None) != (test_status == "tested"):
        raise ValueError("Native status and p-value disagree")
    native_q = probability(row[q])
    if test_status == "tested" and native_q is None:
        raise ValueError("Tested native row lacks its corrected probability")
    if test_status != "tested" and native_q != untestable_q:
        raise ValueError("Native undefined status disagrees with corrected probability convention")
    return NativeObservation(
        RawHypothesis(
            hypothesis_id=identity,
            status="tested" if value is not None else "untestable",
            p_value=value,
            reason=None if value is not None else test_status,
            reason_vocabulary=ReasonVocabulary.analysis_reason_v1,
        ),
        native_q,
        method,
        native_family,
        threshold,
        row,
    )


def input_paths(sources: VerifiedSources, product: str, paths: list[Path]) -> tuple[Path, ...]:
    common = [
        Path(sources.comparison.manifest_path),
        Path(sources.core.manifest_path),
        sources.product(product),
        sources.product("experiments"),
        sources.product("feature_catalog"),
    ]
    return tuple(sorted({*common, *paths}))


def metadata(
    sources: VerifiedSources,
    *,
    family_id: str,
    analysis_id: str,
    product: str,
    count: int,
    membership: Any,
    unit: str,
    inference: str,
    estimand: str,
    selection: str,
    history: list[str],
    paths: list[Path],
) -> HypothesisFamilyMetadata:
    basis = sources.bases[product]
    return HypothesisFamilyMetadata(
        family_id=family_id,
        analysis_id=analysis_id,
        analysis_basis_sha256=canonical_sha256(basis.model_dump(mode="json")),
        source_membership_sha256=canonical_sha256(membership),
        retained_feature_count=basis.retained_feature_count,
        expected_family_size=count,
        population=sources.comparison.contrast_population,
        direction_basis=sources.comparison.direction_basis,
        strain_rollup=sources.scope,
        feature_filter_id=basis.feature_filter_id,
        profile_unit=basis.profile_unit,
        inference_unit=inference,
        hypothesis_unit=unit,
        estimand=estimand,
        selection_rule=selection,
        selection_history=history,
        sources=sorted(
            {
                artifact.identity: artifact
                for path in input_paths(sources, product, paths)
                for artifact in [sources.artifact(path)]
            }.values(),
            key=lambda artifact: artifact.identity,
        ),
    )


def noel_family(sources: VerifiedSources) -> FamilyInput:
    arm = sources.arms["noel"]
    table = recorded_arm_output(arm, "directional_screen", "trait_direction_screen.tsv")
    provenance = recorded_arm_output(arm, "directional_screen", "trait_screen_provenance.json")
    features = sources.features["count_profiles"]
    details = json.loads(provenance.read_text(encoding="utf-8"))
    require_fields(
        details,
        {
            "analysis_version": "trait-direction-screen-v1",
            "profile_sha256": file_sha256(sources.product("count_profiles")),
            "experiments_sha256": file_sha256(sources.product("experiments")),
            "family_size_nominal": len(features),
            "correction_primary": "Benjamini-Yekutieli",
        },
    )
    paired, _ = load_profiles(
        sources.product("count_profiles"),
        sources.product("experiments"),
        candidate_slots(features),
        reject_missing_counts=True,
    )
    tallies = tally_direction_signs(study_trait_contrasts(paired, features), features)
    eligible_keys = sorted(
        (row.study_id, row.experiment_id, row.inc_dec)
        for row in drop_zero_total_pairs(paired).profiles
    )
    eligible_studies = len({key[0] for key in eligible_keys})
    family = metadata(
        sources,
        family_id="noel-full-screen",
        analysis_id="trait-direction-screen-v1",
        product="count_profiles",
        count=len(features),
        membership={"features": features, "positive_total_paired_profiles": eligible_keys},
        unit="qualified_trait_two_sided_sign_test",
        inference="study",
        estimand="study mean increased-minus-decreased fraction; two-sided sign test",
        selection="every retained non-MediaDive trait; zero-total pairs excluded; ties unscored",
        history=[
            "Full candidate screen; no focal selection or effect filter.",
            "Native p-values are imported at the TSV export precision (six significant digits).",
            "No focal q reuse: matching labels do not establish hypothesis/eligibility identity.",
        ],
        paths=[Path(arm.manifest_path), table, provenance],
    )

    def declarations() -> Iterator[HypothesisDefinition]:
        for feature in features:
            yield HypothesisDefinition(
                hypothesis_id=feature,
                feature_ids=[feature],
                test="binomial_sign",
                tail="two_sided",
                coordinates=HypothesisCoordinates(),
            )

    def observations() -> Iterator[NativeObservation]:
        for row in rows(table, sources.expected_digest(table)):
            require_fields(
                row,
                {
                    "family_size_nominal": str(len(features)),
                    "family_id": "trait-direction-screen-v1",
                },
            )
            validate_sign_status(row, status="test_status", empty="no_eligible_studies")
            if (
                sum(support_count(row, f"studies_{kind}") for kind in ("up", "down", "tied"))
                != eligible_studies
            ):
                raise ValueError("Screen support differs from independently loaded core cohort")
            actual = tuple(support_count(row, f"studies_{kind}") for kind in ("up", "down", "tied"))
            if row["trait"] not in tallies or actual != tallies[row["trait"]]:
                raise ValueError("Screen sign tallies disagree with core study contrasts")
            yield observation(
                row["trait"],
                row,
                p="binomial_p",
                q="by_q_nominal",
                status="test_status",
                undefined={"all_tied", "no_eligible_studies"},
                method="BY",
                native_family="noel-full-screen",
            )

    return FamilyInput(
        family,
        declarations,
        observations,
        input_paths(sources, "count_profiles", [Path(arm.manifest_path), table, provenance]),
    )


def jed_family(sources: VerifiedSources) -> FamilyInput:
    arm = sources.arms["jed"]
    table = recorded_arm_output(arm, "composition", "trait_log_ratio_stats.tsv")
    provenance = recorded_arm_output(arm, "composition", "trait_log_ratio_provenance.json")
    studies = recorded_arm_output(arm, "composition", "trait_log_ratio_studies.tsv")
    details = json.loads(provenance.read_text(encoding="utf-8"))
    require_fields(
        details,
        {
            "analysis_version": "paired-log-ratio-sign-v1",
            "profile_sha256": file_sha256(sources.product("count_profiles")),
            "experiments_sha256": file_sha256(sources.product("experiments")),
            "trait_ids": {key: list(value) for key, value in TABLE1_TRAITS.items()},
        },
    )
    pairs = list(combinations(TABLE1_TRAITS, 2))
    definitions = []
    ledger_support = {}
    study_rows = list(rows(studies, sources.expected_digest(studies)))
    paired, _ = load_profiles(
        sources.product("count_profiles"),
        sources.product("experiments"),
        {},
        reject_missing_counts=True,
    )
    paired_by_study = Counter(row.study_id for row in paired if row.inc_dec == "inc")
    for numerator, reference in pairs:
        groups = {
            slug: [
                feature
                for feature in sources.features["count_profiles"]
                if feature.rsplit("|", 1)[-1] in TABLE1_TRAITS[slug]
            ]
            for slug in (numerator, reference)
        }
        if any(len(group) != 1 for group in groups.values()):
            raise ValueError("Ratio trait must resolve to exactly one retained qualified feature")
        identifiers = [groups[numerator][0], groups[reference][0]]
        eligibility = [
            row
            for row in study_rows
            if row["numerator"] == numerator and row["reference"] == reference
        ]
        if (
            not eligibility
            or len({row["study_id"] for row in eligibility}) != len(eligibility)
            or {row["study_id"] for row in eligibility} != set(paired_by_study)
        ):
            raise ValueError("Missing or duplicate ratio eligibility records")
        for row in eligibility:
            expected = paired_by_study[row["study_id"]]
            if (
                support_count(row, "paired_experiments") != expected
                or support_count(row, "eligible_experiments")
                + support_count(row, "excluded_zero_experiments")
                != expected
            ):
                raise ValueError("Ratio eligibility support disagrees with core paired cohort")
        ledger_support[f"{numerator}/{reference}"] = log_ratio_support(
            [
                StudyLogRatio(
                    numerator=numerator,
                    reference=reference,
                    study_id=row["study_id"],
                    paired_experiments=support_count(row, "paired_experiments"),
                    eligible_experiments=support_count(row, "eligible_experiments"),
                    excluded_zero_experiments=support_count(row, "excluded_zero_experiments"),
                    mean_log2_ratio_change=float(row["mean_log2_ratio_change"])
                    if row["mean_log2_ratio_change"]
                    else None,
                )
                for row in eligibility
            ]
        )
        definitions.append(
            HypothesisDefinition(
                hypothesis_id=f"{numerator}/{reference}",
                feature_ids=identifiers,
                test="binomial_sign_log2_ratio",
                tail="two_sided",
                coordinates=HypothesisCoordinates(),
                eligibility_sha256=canonical_sha256(
                    sorted(eligibility, key=lambda row: row["study_id"])
                ),
            )
        )
    if len(study_rows) != sum(
        1 for row in study_rows if (row["numerator"], row["reference"]) in pairs
    ):
        raise ValueError("Unexpected ratio eligibility membership")
    family = metadata(
        sources,
        family_id="jed-15-log-ratios",
        analysis_id="paired-log-ratio-sign-v1",
        product="count_profiles",
        count=len(pairs),
        membership=[definition.model_dump(mode="json") for definition in definitions],
        unit="ordered_trait_ratio_two_sided_sign_test",
        inference="study",
        estimand="study mean log2(A_inc/B_inc)-log2(A_dec/B_dec); two-sided sign test",
        selection="15 focal pairs; four positive counts per eligible experiment; no pseudocount",
        history=[
            "Six focal annotations were selected historically; BY does not correct that selection.",
            "Ratio order and pair-specific study eligibility are part of the hypothesis.",
        ],
        paths=[Path(arm.manifest_path), table, provenance, studies],
    )

    def observations() -> Iterator[NativeObservation]:
        for row in rows(table, sources.expected_digest(table)):
            require_fields(
                row,
                {
                    "analysis_version": "paired-log-ratio-sign-v1",
                    "inference_level": "pooled",
                    "correction_method": "Benjamini-Yekutieli",
                    "family_size": "15",
                },
            )
            validate_sign_status(row, status="status", empty="no_positive_pairs")
            identity = f"{row['numerator']}/{row['reference']}"
            if identity not in ledger_support:
                raise ValueError("Ratio summary is absent from the recorded study ledger")
            expected_support = asdict(ledger_support[identity])
            if any(
                (row[key] if key == "status" else support_count(row, key)) != value
                for key, value in expected_support.items()
            ):
                raise ValueError("Ratio summary support/signs disagree with recorded study ledger")
            if support_count(row, "studies_total") != len(paired_by_study):
                raise ValueError("Ratio summary support disagrees with core paired cohort")
            yield observation(
                identity,
                row,
                p="sign_p",
                q="by_q",
                status="status",
                undefined={"all_tied", "no_positive_pairs"},
                method="BY",
                native_family="jed-15-log-ratios",
            )

    return FamilyInput(
        family,
        lambda: iter(definitions),
        observations,
        input_paths(
            sources, "count_profiles", [Path(arm.manifest_path), table, provenance, studies]
        ),
    )
