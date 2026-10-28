"""Explicit annotation coverage and feature identity for core products (#629).

Coverage refers to matrix lookup and retained annotations, never in vivo activity.
The aggregate matrix does not preserve the source edges needed to distinguish
direct from descendant assertions; its evidence origin is explicitly unresolved.
Under a strain rollup scope, the features it admits may also count BacDive strain
records at their linked NCBITaxon nodes, a third unresolved origin (#788).
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict
from pathlib import Path

from mtp.count_input import COVERAGE_COLUMNS, ProfileCoverage, annotation_coverage_status
from mtp.study_trait_profiles import ProfileKey, TaxaTraitMatrix
from mtp.taxa_trait_matrix import StrainRollup, strain_trait_rolls_up
from mtp.trait_key_policy import (
    DEFAULT_TRAIT_KEY_POLICY,
    TRAIT_POLICY_METADATA_COLUMNS,
    TraitKeyPolicy,
    parse_trait_key_policy,
)

FEATURE_METADATA_COLUMNS = (
    "column_id",
    "column_label",
    "predicate_id",
    "predicate_name",
    "object_id",
    "object_name",
    "object_category",
)
FEATURE_CATALOG_COLUMNS = (*FEATURE_METADATA_COLUMNS, "feature_class", "annotation_evidence")
UNRESOLVED_EVIDENCE = "direct_or_descendant_unresolved"
STRAIN_ATTRIBUTED_EVIDENCE = "direct_descendant_or_strain_unresolved"


def feature_metadata_columns(*, include_policy_metadata: bool) -> tuple[str, ...]:
    """Select an explicit serialization contract independently of scientific policy."""
    return (
        (*FEATURE_METADATA_COLUMNS, *TRAIT_POLICY_METADATA_COLUMNS)
        if include_policy_metadata
        else FEATURE_METADATA_COLUMNS
    )


def feature_catalog_columns(*, include_policy_metadata: bool) -> tuple[str, ...]:
    """The catalog retains every field of its source metadata format."""
    return (
        *feature_metadata_columns(include_policy_metadata=include_policy_metadata),
        "feature_class",
        "annotation_evidence",
    )


def feature_filter_id(feature_columns: Sequence[str]) -> str:
    """Fingerprint the ordered retained feature vocabulary, including MediaDive policy."""
    if len(set(feature_columns)) != len(feature_columns) or any(not c for c in feature_columns):
        raise ValueError("Feature vocabulary must contain unique, nonempty IDs")
    payload = json.dumps(list(feature_columns), ensure_ascii=True, separators=(",", ":"))
    return "sha256:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def annotation_evidence(feature_id: str, strain_rollup: StrainRollup) -> str:
    """The origins a feature's counts may combine; the aggregate never resolves which.

    A feature the scope rolls up from BacDive strains may count strain records at their
    linked NCBITaxon nodes as well as direct and descendant assertions (#788).
    """
    if strain_trait_rolls_up(feature_id, strain_rollup):
        return STRAIN_ATTRIBUTED_EVIDENCE
    return UNRESOLVED_EVIDENCE


def product_evidence(feature_ids: Iterable[str], strain_rollup: StrainRollup) -> str:
    """The origins any count in a product over these features may combine (#788)."""
    if any(
        annotation_evidence(feature, strain_rollup) == STRAIN_ATTRIBUTED_EVIDENCE
        for feature in feature_ids
    ):
        return STRAIN_ATTRIBUTED_EVIDENCE
    return UNRESOLVED_EVIDENCE


def feature_class(feature_id: str) -> str:
    """Use the native downstream MediaDive identifier rule, without dropping columns."""
    return "culture_medium" if "mediadive.medium" in feature_id.lower() else "trait_annotation"


def build_profile_coverage(
    profile_taxa: Mapping[ProfileKey, Sequence[str]],
    matrix: TaxaTraitMatrix,
    *,
    feature_columns: Sequence[str] | None = None,
) -> dict[ProfileKey, ProfileCoverage]:
    """Reconcile distinct reported IDs against matrix rows and retained positive counts.

    A KG-absent taxon emitted as an all-zero matrix row is lookup-matched and
    unannotated. Only an absent matrix row is lookup-unmatched. Mixed ranks are
    retained as reported IDs; a genus and its reported species remain two IDs.
    """
    columns = tuple(matrix.trait_columns if feature_columns is None else feature_columns)
    unknown = set(columns).difference(matrix.trait_columns)
    if unknown:
        raise ValueError(f"Coverage filter contains unknown matrix features: {sorted(unknown)}")
    filter_id = feature_filter_id(columns)
    result: dict[ProfileKey, ProfileCoverage] = {}
    for key, observations in profile_taxa.items():
        if len(key) != 3 or any(not value.strip() for value in key) or key[2] not in {"inc", "dec"}:
            raise ValueError(f"Invalid profile key: {key!r}")
        reported = set(observations)
        if any(not taxon.strip() for taxon in reported):
            raise ValueError(f"Missing reported taxon ID in profile {key!r}")
        matched = reported.intersection(matrix.counts_by_taxon)
        unmatched = reported.difference(matched)
        unannotated = {
            taxon
            for taxon in matched
            if not any(matrix.counts_by_taxon[taxon].get(column, 0) > 0 for column in columns)
        }
        annotated = matched.difference(unannotated)
        status = annotation_coverage_status(
            reported=len(reported),
            matched=len(matched),
            unmatched=len(unmatched),
            with_annotations=len(annotated),
            without_annotations=len(unannotated),
        )
        result[key] = ProfileCoverage(
            *key,
            len(reported),
            len(matched),
            len(unmatched),
            len(annotated),
            len(unannotated),
            status,
            tuple(sorted(unmatched)),
            tuple(sorted(unannotated)),
            filter_id,
        )
    return result


def write_profile_coverage(path: Path, coverage: Mapping[ProfileKey, ProfileCoverage]) -> None:
    """Write the count-input-v1 sidecar; numeric zero and unmatched IDs stay distinct."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=COVERAGE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(coverage):
            record = coverage[key]
            if record.key != key:
                raise ValueError(f"Coverage mapping key disagrees with record: {key!r}")
            values = asdict(record)
            for column in ("unmatched_taxon_ids", "matched_without_annotation_taxon_ids"):
                values[column] = json.dumps(values[column], separators=(",", ":"))
            writer.writerow(values)


def write_feature_catalog(
    metadata_path: Path,
    output_path: Path,
    *,
    strain_rollup: StrainRollup = "off",
    trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY,
    include_policy_metadata: bool = True,
) -> None:
    """Extend native feature metadata with explicit culture-medium and origin semantics."""
    # Local import avoids a module cycle while sharing strict TSV validation.
    from mtp.workflow.contracts import read_table, validate_product

    trait_key_policy = parse_trait_key_policy(trait_key_policy)
    columns = feature_metadata_columns(include_policy_metadata=include_policy_metadata)
    validate_product(
        metadata_path,
        "feature_metadata",
        trait_key_policy=trait_key_policy,
        include_policy_metadata=include_policy_metadata,
    )
    _, records = read_table(metadata_path, required=columns, keys=("column_id",))
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=feature_catalog_columns(include_policy_metadata=include_policy_metadata),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    **{column: record[column] for column in columns},
                    "feature_class": feature_class(record["column_id"]),
                    "annotation_evidence": annotation_evidence(record["column_id"], strain_rollup),
                }
            )
