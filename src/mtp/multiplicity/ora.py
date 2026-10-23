"""One joint upper/lower ORA family with independently declared coordinates (#853)."""

from __future__ import annotations

import json
import math
from collections.abc import Iterator
from pathlib import Path

from mtp.count_input import taxa_trait_count_reader
from mtp.datamodel.multiple_testing_pydantic import (
    HypothesisCoordinates,
    HypothesisDefinition,
    RawHypothesis,
)
from mtp.integrated_evidence import ora_testable_profile_keys
from mtp.multiple_testing import canonical_sha256
from mtp.multiplicity.families import (
    FamilyInput,
    NativeObservation,
    input_paths,
    metadata,
    probability,
    require_fields,
    rows,
)
from mtp.multiplicity.sources import VerifiedSources, recorded_arm_output
from mtp.paths import file_sha256


def ora_identity(key: tuple[str, str, str], trait: str, tail: str) -> str:
    """Unambiguous tuple encoding, even if names contain separator characters."""
    return json.dumps([*key, trait, tail], ensure_ascii=False, separators=(",", ":"))


def ora_comments(path: Path) -> dict[str, str]:
    comments: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            key, separator, value = line[1:].partition(":")
            if separator:
                if key.strip() in comments:
                    raise ValueError("Duplicate ORA metadata")
                comments[key.strip()] = value.strip()
    return comments


def ora_family(sources: VerifiedSources) -> FamilyInput:
    arm = sources.arms["ilyas"]
    exports = {
        "upper": recorded_arm_output(arm, "enrichment", "trait_enrichment_all_results.tsv"),
        "lower": recorded_arm_output(arm, "depletion", "trait_depletion_all_results.tsv"),
    }
    profiles = sources.product("ora_profiles")
    matrix = sources.product("taxa_matrix")
    keys = sorted(ora_testable_profile_keys(profiles))
    reader = taxa_trait_count_reader(matrix)
    positions = {column: index for index, column in enumerate(reader.count_columns)}
    retained = sources.features["ora_profiles"]
    if not set(retained).issubset(positions):
        raise ValueError("ORA background lacks retained features")
    prevalence = dict.fromkeys(retained, 0)
    universe = 0
    for row in reader:
        universe += 1
        for feature in retained:
            prevalence[feature] += row.counts[positions[feature]] > 0
    candidates = tuple(feature for feature in retained if prevalence[feature] > 0)
    candidate_set = set(candidates)
    tests_per_tail = len(keys) * len(candidates)
    expected = {
        "Profiles-SHA256": file_sha256(profiles),
        "Taxa-Matrix-SHA256": file_sha256(matrix),
        "Testable-Set": "K>0 & N>0",
        "FDR": "True",
        "FDR-Scope": "global",
        "P-Threshold": "0.001",
        "Testable-Profiles": str(len(keys)),
        "Background-Taxa-Universe-M": str(universe),
        "Background-Mode": f"TaxaMatrix ({matrix.name})",
        "Total-Tests-Performed": str(tests_per_tail),
    }
    for tail, path in exports.items():
        comments = ora_comments(path)
        require_fields(comments, {**expected, "Tail": tail})
        fold_key = "Min-Fold" if tail == "upper" else "Max-Fold"
        if float(comments[fold_key]) != (1.5 if tail == "upper" else 0.5):
            raise ValueError("Unsupported native ORA fold threshold")
    eligibility_hash = canonical_sha256(keys)
    family = metadata(
        sources,
        family_id="ilyas-joint-ora",
        analysis_id="hypergeometric-joint-tails-v1",
        product="ora_profiles",
        count=2 * tests_per_tail,
        membership={
            "profiles": keys,
            "traits": candidates,
            "tails": ["upper", "lower"],
            "background_positive_counts": prevalence,
            "background_universe": universe,
        },
        unit="profile_qualified_trait_tail",
        inference="signature_profile",
        estimand="trait annotation over/under-representation against unchanged taxon matrix",
        selection="paired profiles then N>0; retained traits with background K>0; joint tails",
        history=[
            "Native legacy BH families were separate global upper/lower tails.",
            "Additive BY/BH use one complete profile x eligible trait x two-tail family.",
            "Fold cutoffs remain display/effect filters, separate from multiplicity rejection.",
        ],
        paths=[Path(arm.manifest_path), matrix, *exports.values()],
    )

    def declarations() -> Iterator[HypothesisDefinition]:
        for key in keys:
            for feature in candidates:
                for tail in ("upper", "lower"):
                    yield HypothesisDefinition(
                        hypothesis_id=ora_identity(key, feature, tail),
                        feature_ids=[feature],
                        test="hypergeometric",
                        tail=tail,
                        coordinates=HypothesisCoordinates(
                            profile_id=json.dumps(
                                list(key), ensure_ascii=False, separators=(",", ":")
                            )
                        ),
                        eligibility_sha256=eligibility_hash,
                    )

    def observations() -> Iterator[NativeObservation]:
        for tail, path in exports.items():
            count = 0
            for row in rows(path, sources.expected_digest(path)):
                count += 1
                feature = row["trait"]
                require_fields(
                    row,
                    {
                        "eligibility_status": "testable",
                        "correction_family": "global",
                        "background_m": str(universe),
                        "testable_profiles": str(len(keys)),
                    },
                )
                if feature not in candidate_set or int(row["background_k"]) != prevalence[feature]:
                    raise ValueError("ORA background feature prevalence mismatch")
                key = (row["study_id"], row["experiment_id"], row["inc_dec"])
                p, q = probability(row["p_value"]), probability(row["fdr_q"])
                if p is None or q is None:
                    raise ValueError("ORA testable member lacks raw or adjusted probability")
                fold = float(row["fold"])
                if not math.isfinite(fold) or fold < 0:
                    raise ValueError("Invalid ORA effect")
                called = q <= 0.001 and (fold >= 1.5 if tail == "upper" else fold <= 0.5)
                if row["is_significant"].lower() != str(called).lower():
                    raise ValueError("ORA native hit flag disagrees with q/effect thresholds")
                yield NativeObservation(
                    RawHypothesis(
                        hypothesis_id=ora_identity(key, feature, tail), status="tested", p_value=p
                    ),
                    q,
                    "BH",
                    f"ilyas-{tail}-global",
                    0.001,
                    row,
                )
            if count != tests_per_tail:
                raise ValueError("ORA all-results export is truncated or has unexpected members")

    return FamilyInput(
        family,
        declarations,
        observations,
        input_paths(sources, "ora_profiles", [Path(arm.manifest_path), matrix, *exports.values()]),
    )
