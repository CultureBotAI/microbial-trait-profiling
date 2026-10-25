"""Native ORA calls and declared study enrichment preference (#1044).

ora-study-enrichment-direction-v1 was declared before computation. Native hits
are upper-tail, separate global BH q <= .001 and fold >= 1.5. The directional
population contains only experiment pairs with both native profile N > 0 and
traits with background K > 0. Each study contributes the mean paired inc-hit
minus dec-hit difference. The integer numerator determines its sign; exact zero
is tied. A two-sided exact binomial sign test excludes ties. BY and BH span the
whole retained catalog, holding untestable members at p=1 internally and reporting
blank p/q. The exact 95% Clopper-Pearson interval describes the positive share.

The global BH call threshold depends on all traits and lists, so thresholded hits
are not independent across traits. This describes enrichment preference only;
there is no depletion-direction test and no ORA headline decision. Calibration
status is uncalibrated (methods/calibration.md excludes native ORA).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from mtp.count_input import parse_nonnegative_integer
from mtp.enrichment_stats import benjamini_hochberg
from mtp.multiplicity.families import rows
from mtp.multiplicity.ora import ora_family
from mtp.multiplicity.sources import VerifiedSources

DIRECTION_FAMILY = "ora-study-enrichment-direction-v1"
CALIBRATION_STATUS = "uncalibrated (methods/calibration.md excludes native ORA)"
ProfileKey = tuple[str, str, str]
PairKey = tuple[str, str]
_NULLABLE_COUNTS = ("ora_enriched_inc", "ora_enriched_dec", "ora_depleted_inc", "ora_depleted_dec")


@dataclass(frozen=True)
class _Population:
    """Independent paired inventory before and after the native positive-N rule."""

    testable: tuple[ProfileKey, ...]
    paired: tuple[PairKey, ...]
    common: tuple[PairKey, ...]
    inc_larger_pairs: int
    dec_larger_pairs: int
    equal_size_pairs: int


def _population(sources: VerifiedSources) -> _Population:
    path = sources.product("ora_profiles")
    draws: dict[ProfileKey, int] = {}
    directions: dict[PairKey, set[str]] = {}
    for line, row in enumerate(rows(path, sources.expected_digest(path)), start=2):
        key = (row["study_id"], row["experiment_id"], row["inc_dec"])
        if not all(key) or key[2] not in {"inc", "dec"} or key in draws:
            raise ValueError("ORA profiles have duplicate or invalid identifiers")
        draws[key] = parse_nonnegative_integer(
            row["taxa_count"], path=path, line_number=line, column_name="taxa_count"
        )
        directions.setdefault(key[:2], set()).add(key[2])
    paired = tuple(sorted(key for key, values in directions.items() if values == {"inc", "dec"}))
    paired_set = set(paired)
    testable = tuple(sorted(key for key, value in draws.items() if key[:2] in paired_set and value))
    common = tuple(key for key in paired if draws[(*key, "inc")] > 0 and draws[(*key, "dec")] > 0)
    return _Population(
        testable,
        paired,
        common,
        sum(draws[(*key, "inc")] > draws[(*key, "dec")] for key in common),
        sum(draws[(*key, "inc")] < draws[(*key, "dec")] for key in common),
        sum(draws[(*key, "inc")] == draws[(*key, "dec")] for key in common),
    )


def _background(sources: VerifiedSources, features: tuple[str, ...]) -> tuple[int, np.ndarray]:
    """Count the unchanged background, with a hash of exactly the consumed bytes."""
    path = sources.product("taxa_matrix")
    carriers = np.zeros(len(features), dtype=np.int64)
    taxa: set[str] = set()
    for line, row in enumerate(rows(path, sources.expected_digest(path)), start=2):
        taxon = row["taxon_id"]
        if not taxon or taxon in taxa:
            raise ValueError("ORA background has duplicate or blank taxon identifiers")
        taxa.add(taxon)
        for index, feature in enumerate(features):
            carriers[index] += (
                parse_nonnegative_integer(
                    row[feature], path=path, line_number=line, column_name=feature
                )
                > 0
            )
    if not taxa:
        raise ValueError("ORA background contains no taxa")
    return len(taxa), carriers


def _study_sums(population: _Population, upper: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    studies = sorted({key[0] for key in population.common})
    positions = {study: index for index, study in enumerate(studies)}
    profile_positions = {key: index for index, key in enumerate(population.testable)}
    numerators = np.zeros((len(studies), upper.shape[1]), dtype=np.int64)
    pair_counts = Counter(key[0] for key in population.common)
    for study, experiment in population.common:
        inc = upper[profile_positions[(study, experiment, "inc")]].astype(np.int64)
        dec = upper[profile_positions[(study, experiment, "dec")]].astype(np.int64)
        numerators[positions[study]] += inc - dec
    denominators = np.asarray([pair_counts[study] for study in studies], dtype=np.int64)
    return numerators, denominators


def _direction_row(
    numerator: np.ndarray,
    denominators: np.ndarray,
    *,
    background_k: int,
    population: _Population,
) -> dict[str, Any]:
    eligible = len(denominators) if background_k else 0
    eligible_pairs = len(population.common) if background_k else 0
    up = int(np.count_nonzero(numerator > 0)) if background_k else 0
    down = int(np.count_nonzero(numerator < 0)) if background_k else 0
    tied = eligible - up - down
    if not background_k:
        status = "no_background_taxa"
    elif not eligible_pairs:
        status = "no_eligible_pairs"
    else:
        status = "tested" if up + down else "all_tied"
    test = binomtest(up, up + down, p=0.5, alternative="two-sided") if up + down else None
    interval = test.proportion_ci(confidence_level=0.95, method="exact") if test else None
    preference = "none"
    if up > down:
        preference = "positive"
    elif down > up:
        preference = "negative"
    elif up:
        preference = "mixed"
    return {
        "ora_direction_status": status,
        "ora_direction_vote_preference": preference if eligible else "unavailable",
        "ora_direction_studies_up": up,
        "ora_direction_studies_down": down,
        "ora_direction_studies_tied": tied,
        "ora_direction_studies_scored": up + down,
        "ora_direction_studies_eligible": eligible,
        "ora_direction_studies_excluded": len({key[0] for key in population.paired}) - eligible,
        "ora_direction_studies_total": len({key[0] for key in population.paired}),
        "ora_direction_pairs_eligible": eligible_pairs,
        "ora_direction_pairs_n_inc_gt_dec": population.inc_larger_pairs if background_k else 0,
        "ora_direction_pairs_n_inc_lt_dec": population.dec_larger_pairs if background_k else 0,
        "ora_direction_pairs_n_inc_eq_dec": population.equal_size_pairs if background_k else 0,
        "ora_direction_pairs_excluded": len(population.paired) - eligible_pairs,
        "ora_direction_pairs_total": len(population.paired),
        "ora_direction_mean_study_delta": float(np.mean(numerator / denominators))
        if eligible
        else None,
        "ora_direction_share_up": up / (up + down) if test else None,
        "ora_direction_share_up_ci95_low": float(interval.low) if interval else None,
        "ora_direction_share_up_ci95_high": float(interval.high) if interval else None,
        "ora_direction_p": float(test.pvalue) if test else None,
        "ora_direction_family_id": DIRECTION_FAMILY,
        "ora_direction_inference_unit": "study",
        "ora_direction_primary_correction": "BY",
    }


def ora_table(sources: VerifiedSources) -> tuple[pd.DataFrame, tuple[Path, ...]]:
    """Reuse complete validated native calls, then summarize the declared study test.

    Native validation uses disk-backed SQLite and exact per-tail BH comparisons.
    Only boolean profile-by-trait calls and small integer count arrays remain in
    memory here. The caller must route tempfile scratch through a durable TMPDIR.
    An exception anywhere in either complete tail aborts before a table returns.
    """
    features = sources.features["ora_profiles"]
    if not features or len(features) != len(set(features)):
        raise ValueError("ORA retained trait catalog must be nonempty and unique")
    family = ora_family(sources)
    population = _population(sources)
    universe, carriers = _background(sources, features)
    feature_positions = {feature: index for index, feature in enumerate(features)}
    profile_positions = {key: index for index, key in enumerate(population.testable)}
    upper = np.zeros((len(population.testable), len(features)), dtype=bool)
    counts = np.zeros((2, 2, len(features)), dtype=np.int64)
    tails = {"ilyas-upper-global": 0, "ilyas-lower-global": 1}
    directions = {"inc": 0, "dec": 1}
    for observation in family.observations():
        row = observation.row
        feature_index = feature_positions[row["trait"]]
        key = (row["study_id"], row["experiment_id"], row["inc_dec"])
        profile = profile_positions[key]
        tail = tails[observation.native_family]
        if row["is_significant"].lower() == "true":
            counts[tail, directions[key[2]], feature_index] += 1
            if tail == 0:
                upper[profile, feature_index] = True
    numerators, denominators = _study_sums(population, upper)
    testable_directions = Counter(key[2] for key in population.testable)
    records = []
    for index, feature in enumerate(features):
        background_k = int(carriers[index])
        records.append(
            {
                "trait": feature,
                "ora_enriched_inc": int(counts[0, 0, index]) if background_k else None,
                "ora_enriched_dec": int(counts[0, 1, index]) if background_k else None,
                "ora_depleted_inc": int(counts[1, 0, index]) if background_k else None,
                "ora_depleted_dec": int(counts[1, 1, index]) if background_k else None,
                "ora_native_testable_inc_profiles": testable_directions["inc"]
                if background_k
                else 0,
                "ora_native_testable_dec_profiles": testable_directions["dec"]
                if background_k
                else 0,
                "ora_background_m": universe,
                "ora_background_k": background_k,
                "ora_native_tests_per_tail": family.metadata.expected_family_size // 2,
                "ora_native_correction": "BH_separate_global_tails",
                "ora_native_q_threshold": 0.001,
                "ora_enrichment_min_fold": 1.5,
                "ora_depletion_max_fold": 0.5,
                "ora_calibration_status": CALIBRATION_STATUS,
                "ora_direction_family_size": len(features),
                **_direction_row(
                    numerators[:, index],
                    denominators,
                    background_k=background_k,
                    population=population,
                ),
            }
        )
    held = [
        record["ora_direction_p"] if record["ora_direction_status"] == "tested" else 1.0
        for record in records
    ]
    for method in ("by", "bh"):
        adjusted = benjamini_hochberg(held, family_size=len(features), method=method)
        for record, q in zip(records, adjusted, strict=True):
            record[f"ora_direction_{method}_q"] = (
                q if record["ora_direction_status"] == "tested" else None
            )
    table = pd.DataFrame(records).astype(dict.fromkeys(_NULLABLE_COUNTS, "Int64"))
    return table, family.source_paths
