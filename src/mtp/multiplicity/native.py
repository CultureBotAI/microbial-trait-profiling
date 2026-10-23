"""Native pooled composition adapter: reconstruct selection, import primary p-values."""

from __future__ import annotations

import math
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd

from mtp.count_input import profile_count_reader
from mtp.datamodel.multiple_testing_pydantic import HypothesisCoordinates, HypothesisDefinition
from mtp.multiple_testing import canonical_sha256
from mtp.multiplicity.families import (
    FamilyInput,
    NativeObservation,
    input_paths,
    metadata,
    observation,
    require_fields,
    rows,
    support_count,
)
from mtp.multiplicity.sources import VerifiedSources, recorded_comparison_output
from mtp.paths import file_sha256
from mtp.trait_guilds import load_profiles
from mtp.workflow.consumers import verify_consumer_record


def validate_study_status(row: dict[str, str], *, expected_studies: int) -> None:
    """Check producer support/dispersion and undefined reason without retesting p."""
    count = support_count(row, "study_count")
    mean = float(row["study_mean_increased_minus_decreased"])
    std = float(row["study_std_increased_minus_decreased"])
    if count != expected_studies or not math.isfinite(mean) or not math.isfinite(std) or std < 0:
        raise ValueError("Native study support or dispersion disagrees with its producer contract")
    if count < 2:
        expected = "insufficient_support"
    elif np.isclose(std, 0.0):
        expected = "all_tied" if std == 0 and mean == 0 else "insufficient_variance"
    else:
        expected = "tested"
    if row["study_test_status"] != expected:
        raise ValueError("Native study status/reason disagrees with support and dispersion")


def native_candidates(
    profile: Path, experiments: Path, features: tuple[str, ...], labels: Path
) -> tuple[tuple[str, ...], list[tuple[str, str, str]]]:
    """Reconstruct the supported native default filter without reading selected results.

    Support 6..floor(.95*n), nonzero fraction variance, no maximum trait count.
    Fractions use all retained traits, and historical zero-total pairs remain.
    The scientific tests and transformed effects are never recomputed here.
    """
    profiles, _ = load_profiles(profile, experiments, {}, reject_missing_counts=True)
    keys = sorted((row.study_id, row.experiment_id, row.inc_dec) for row in profiles)
    if not keys:
        raise ValueError("Native composition has no paired human profiles")
    included = set(keys)
    reader = profile_count_reader(profile)
    indices = [reader.count_columns.index(feature) for feature in features]
    counts = pd.DataFrame(
        [[row.counts[index] for index in indices] for row in reader if row.key in included],
        columns=list(features),
        dtype=float,
    )
    if len(counts) != len(keys):
        raise ValueError("Native composition profile inventory mismatch")
    fractions = counts.div(counts.sum(axis=1).replace(0, float("nan")), axis=0).fillna(0)
    support, signal, variance = counts.gt(0).sum(), counts.sum(), fractions.var()
    label_rows = {row["column_id"]: row for row in rows(labels)}
    candidates = []
    for feature in features:
        label = label_rows.get(feature, {})
        name = label.get("object_name", "").strip().replace("_", " ")
        if not name or name == "unknown":
            name = label.get("column_label", "").strip().replace("_", " ")
        trait_label = f"{name} {feature.rsplit('|', 1)[-1]}"
        if (
            6 <= support[feature] <= max(2, int(len(counts) * 0.95))
            and variance[feature] > 0
            and "mediadive.medium" not in trait_label
        ):
            candidates.append(feature)
    candidates.sort(
        key=lambda feature: (
            -variance[feature] * math.log1p(signal[feature]),
            -variance[feature],
            -support[feature],
            feature,
        )
    )
    if len(candidates) < 2:
        raise ValueError("Fewer than two native traits satisfy the declared filter")
    return tuple(candidates), keys


def native_family(sources: VerifiedSources) -> FamilyInput:
    table = recorded_comparison_output(
        sources.comparison, "trait_profile_paired_contrast_trait_stats_pooled.tsv"
    )
    pairs = recorded_comparison_output(
        sources.comparison, "trait_profile_paired_contrast_pair_metadata.tsv"
    )
    consumer_path = table.parent / "core_consumption.json"
    recorded = {
        (Path(sources.comparison.manifest_path).parent / entry.path).resolve()
        for entry in sources.comparison.outputs
    }
    if consumer_path.resolve() not in recorded:
        raise ValueError("Native composition lacks a comparison-bound consumer record")
    consumer = verify_consumer_record(consumer_path)
    if (
        consumer.consumer != "native-pooled-composition-stats"
        or consumer.run_id != sources.core.run_id
        or consumer.run_manifest_sha256 != sources.comparison.core_manifest_sha256
        or not {str(table), str(pairs)}.issubset({entry.path for entry in consumer.outputs})
    ):
        raise ValueError("Native composition consumer provenance mismatch")
    features, keys = native_candidates(
        sources.product("count_profiles"),
        sources.product("experiments"),
        sources.features["count_profiles"],
        sources.product("feature_metadata"),
    )
    expected_hashes = {
        "profile_sha256": file_sha256(sources.product("count_profiles")),
        "experiments_sha256": file_sha256(sources.product("experiments")),
        "trait_labels_sha256": file_sha256(sources.product("feature_metadata")),
    }
    pair_keys = []
    for row in rows(pairs, sources.expected_digest(pairs)):
        require_fields(
            row,
            {**expected_hashes, "zero_total_policy": "included_in_historical_native_composition"},
        )
        pair_keys.extend(
            [(row["study_id"], row["experiment_id"], direction) for direction in ("inc", "dec")]
        )
    if sorted(pair_keys) != keys:
        raise ValueError("Native pair metadata differs from core paired population")
    family = metadata(
        sources,
        family_id="native-pooled-composition",
        analysis_id="native-pooled-study-t-v1",
        product="count_profiles",
        count=len(features),
        membership={"selected_traits": features, "profiles": keys},
        unit="filtered_qualified_trait_two_sided_t_test",
        inference="study",
        estimand="study mean paired shift of standardized log1p(1000 * retained trait fraction)",
        selection=(
            "paired human profiles; support>=6 and <=max(2,floor(.95*n)); "
            "fraction variance>0; no max_traits"
        ),
        history=[
            "Support/variance filtering defines this family; BY does not fix selection bias.",
            "Retained catalog count and filtered test family count are distinct.",
            "Historical zero-total profiles included; primary study p-values imported unchanged.",
        ],
        paths=[table, pairs, consumer_path, sources.product("feature_metadata")],
    )
    eligibility_hash = canonical_sha256(keys)
    expected_studies = len({key[0] for key in keys})

    def declarations() -> Iterator[HypothesisDefinition]:
        for feature in features:
            yield HypothesisDefinition(
                hypothesis_id=feature,
                feature_ids=[feature],
                test="two_sided_one_sample_t_test",
                tail="two_sided",
                coordinates=HypothesisCoordinates(),
                eligibility_sha256=eligibility_hash,
            )

    def observations() -> Iterator[NativeObservation]:
        for row in rows(table, sources.expected_digest(table)):
            require_fields(
                row,
                {
                    **expected_hashes,
                    "pair_metadata_sha256": file_sha256(pairs),
                    "stratum_kind": "pooled",
                    "stratum": "pooled",
                    "test_method": "two_sided_one_sample_t_test",
                    "correction_method": "BH",
                    "correction_family_size": str(len(features)),
                    "retained_trait_count": str(
                        sources.bases["count_profiles"].retained_feature_count
                    ),
                },
            )
            validate_study_status(row, expected_studies=expected_studies)
            yield observation(
                row["trait"],
                row,
                p="study_p_value",
                q="study_q_value",
                status="study_test_status",
                undefined={"all_tied", "insufficient_support", "insufficient_variance"},
                method="BH",
                native_family="native-pooled-composition",
                untestable_q=1.0,
            )

    return FamilyInput(
        family,
        declarations,
        observations,
        input_paths(
            sources,
            "count_profiles",
            [table, pairs, consumer_path, sources.product("feature_metadata")],
        ),
    )
