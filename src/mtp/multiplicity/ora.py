"""One joint upper/lower ORA family with independently declared coordinates (#853)."""

from __future__ import annotations

import json
import math
import sqlite3
import tempfile
from collections.abc import Iterator
from contextlib import closing
from pathlib import Path

import numpy as np

from mtp.count_input import taxa_trait_count_reader
from mtp.datamodel.multiple_testing_pydantic import (
    HypothesisCoordinates,
    HypothesisDefinition,
    RawHypothesis,
    ReasonVocabulary,
)
from mtp.enrichment_stats import benjamini_hochberg
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


def native_bh_rows(
    path: Path,
    *,
    expected_sha256: str | None,
    profile_keys: set[tuple[str, str, str]],
    features: set[str],
) -> Iterator[dict[str, str]]:
    """Check one complete native tail on disk before yielding its original rows.

    Unique declared profile/feature coordinates plus the exact Cartesian size
    establish membership. Only numeric correction arrays reside in memory;
    source rows stay in SQLite and are never reread after their stream hash is
    checked. Native probabilities are compared, never substituted.
    """
    expected_count = len(profile_keys) * len(features)
    with (
        tempfile.TemporaryDirectory(prefix="mtp-native-ora-") as scratch,
        closing(sqlite3.connect(Path(scratch) / "native.sqlite")) as db,
        db,
    ):
        db.execute(
            "CREATE TABLE native (ordinal INTEGER PRIMARY KEY, identity TEXT UNIQUE NOT NULL, "
            "record TEXT NOT NULL, p REAL NOT NULL, q REAL NOT NULL)"
        )
        count = 0
        for row in rows(path, expected_sha256):
            key = (row["study_id"], row["experiment_id"], row["inc_dec"])
            if key not in profile_keys or row["trait"] not in features:
                raise ValueError("ORA native BH family has an undeclared coordinate")
            p, q = probability(row["p_value"]), probability(row["fdr_q"])
            if p is None or q is None:
                raise ValueError("ORA testable member lacks raw or adjusted probability")
            try:
                db.execute(
                    "INSERT INTO native VALUES (?, ?, ?, ?, ?)",
                    (
                        count,
                        json.dumps([*key, row["trait"]], separators=(",", ":")),
                        json.dumps(row, ensure_ascii=False, separators=(",", ":")),
                        p,
                        q,
                    ),
                )
            except sqlite3.IntegrityError as error:
                raise ValueError("ORA native BH family has duplicate coordinates") from error
            count += 1
        if count != expected_count:
            raise ValueError("ORA all-results export is truncated or has unexpected members")
        if count:
            p_values = np.fromiter(
                (entry[0] for entry in db.execute("SELECT p FROM native ORDER BY ordinal")),
                dtype=float,
                count=count,
            )
            adjusted = benjamini_hochberg(p_values, family_size=count)
            for (exported,), expected in zip(
                db.execute("SELECT q FROM native ORDER BY ordinal"), adjusted, strict=True
            ):
                if exported != expected:
                    raise ValueError("ORA native per-tail BH disagrees with complete raw-p family")
            del p_values, adjusted
        for (record,) in db.execute("SELECT record FROM native ORDER BY ordinal"):
            yield json.loads(record)


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
            for row in native_bh_rows(
                path,
                expected_sha256=sources.expected_digest(path),
                profile_keys=set(keys),
                features=candidate_set,
            ):
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
                        hypothesis_id=ora_identity(key, feature, tail),
                        status="tested",
                        p_value=p,
                        reason_vocabulary=ReasonVocabulary.analysis_reason_v1,
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
