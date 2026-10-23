"""Bind issue #772 sensitivity families to the authoritative multiplicity contract (#892)."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from mtp.datamodel.multiple_testing_pydantic import (
    HypothesisCoordinates,
    HypothesisDefinition,
    HypothesisFamily,
    RawHypothesis,
    SourceArtifact,
)
from mtp.multiple_testing import canonical_sha256, correct_family, verify_correction
from mtp.paths import file_sha256
from mtp.workflow.consumers import CoreBundle


def declare_and_correct(
    records: Sequence[dict[str, Any]],
    bases: Mapping[str, dict[str, Any]],
    bundle: CoreBundle,
    ledger: Sequence[dict[str, Any]],
    strain_rollup: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Correct complete declared inventories; native corrections remain separate fields."""
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in records:
        grouped[row["family_id"]].append(row)
    sources = [
        SourceArtifact(identity=f"core:{name}", sha256=file_sha256(bundle.path(name)))
        for name in (
            "count_profiles",
            "ora_profiles",
            "reported_taxa",
            "experiments",
            "taxa_matrix",
            "feature_catalog",
        )
    ]
    families, corrections = [], []
    for identity, rows in sorted(grouped.items()):
        representation, test = rows[0]["representation"], rows[0]["test"]
        focal = identity.endswith(":focal_six")
        ratio = test == "paired_log2_ratio_sign"
        basis = bases[f"{representation}:composition" if ratio else representation]
        family = HypothesisFamily(
            family_id=identity,
            analysis_id="presence-depth-study-v2",
            analysis_basis_sha256=canonical_sha256(basis),
            source_membership_sha256=canonical_sha256(ledger),
            retained_feature_count=basis["retained_feature_count"],
            expected_family_size=len(rows),
            population=bundle.run.contrast_population or "all_paired",
            direction_basis=bundle.run.direction_basis or "group_1_vs_group_0",
            strain_rollup=strain_rollup,
            feature_filter_id=basis["feature_filter_id"],
            profile_unit=basis["profile_unit"],
            inference_unit="study",
            hypothesis_unit="ordered_trait_pair" if ratio else "trait",
            estimand=(
                "study sign of paired log2 annotation ratio change"
                if ratio
                else f"{test} of paired {representation} fraction differences at equal study weight"
            ),
            selection_rule=(
                "all fifteen declared focal trait pairs; missing features remain unavailable"
                if ratio
                else "six declared Table 1 traits; missing features remain unavailable"
                if focal
                else "every retained non-MediaDive core feature, including untestable slots"
            ),
            selection_history=[
                "additive sensitivity declared before rerun; focal traits selected in earlier "
                "exploratory work",
                "comparisons across representations and test families are exploratory, "
                "not replication",
                "six-trait BH is a native-method comparison; full-family BY remains "
                "the primary screen",
            ],
            sources=sources,
            members=[
                HypothesisDefinition(
                    hypothesis_id=row["hypothesis_id"],
                    feature_ids=row["feature_ids"],
                    coordinates=HypothesisCoordinates(),
                    test=test,
                    tail="two_sided",
                )
                for row in rows
            ],
        )
        raw = [
            RawHypothesis(
                hypothesis_id=row["hypothesis_id"],
                status="tested" if row["p_value"] is not None else "untestable",
                p_value=row["p_value"],
                reason=None if row["p_value"] is not None else row["status"],
            )
            for row in rows
        ]
        correction = correct_family(family, raw)
        verify_correction(correction, family=family, hypotheses=raw)
        by_id = {row.hypothesis_id: row for row in correction.results}
        for row in rows:
            adjusted = by_id[row["hypothesis_id"]]
            row.update(q_by=adjusted.q_by, q_bh=adjusted.q_bh, family_size=len(rows))
            if focal:
                row["focal_six_bh_comparison"] = adjusted.q_bh
        families.append(family.model_dump(mode="json"))
        corrections.append(correction.model_dump(mode="json"))
    return families, corrections
