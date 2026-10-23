"""Scientific contracts for existing wide core products (#629).

These checks describe annotation alignment, not measured physiology. Native
matrix cells count distinct annotated source taxa at or below a reported taxon;
count profiles sum those cells, whereas ORA profiles count positive reported
taxon IDs and retain all lookup-matched IDs in their denominator.
"""

from __future__ import annotations

import csv
import hashlib
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path

from pydantic import BaseModel

from mtp.count_input import (
    CountTableReader,
    open_input_text,
    parse_nonnegative_integer,
    profile_count_reader,
    read_profile_coverage,
    taxa_trait_count_reader,
)
from mtp.datamodel.core_products_pydantic import (
    CountProfileRow,
    ExperimentOntologyTermRow,
    ExperimentTableRow,
    FeatureDefinition,
    ORAProfileRow,
    ProfileCoverageRecord,
    SignatureTableRow,
    StudyTableRow,
    TaxonOccurrenceRow,
)
from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES
from mtp.study_trait_profiles import ProfileKey, direction_to_inc_dec, read_taxa_trait_matrix
from mtp.taxa_trait_matrix import StrainRollup
from mtp.workflow.coverage import (
    FEATURE_CATALOG_COLUMNS,
    FEATURE_METADATA_COLUMNS,
    STRAIN_ATTRIBUTED_EVIDENCE,
    annotation_evidence,
    build_profile_coverage,
    feature_class,
    feature_filter_id,
    product_evidence,
)
from mtp.workflow.matrix_contracts import validate_feature_node_metadata, validate_matrix_sources

CONTRACT_VERSION = "core-products-v1"
_CURIE = re.compile(r"[A-Za-z][A-Za-z0-9_.-]*:[^\s|]+\Z")
_ALIASES = {
    "studies": "study_table",
    "experiments": "experiment_table",
    "signatures": "signature_table",
    "taxa": "reported_taxa",
    "experiment_ontology_terms": "ontology_term_table",
    "taxon_trait_matrix": "taxa_matrix",
}
_TABLE_KEYS = {
    "study_table": ("study_id",),
    "experiment_table": ("study_id", "experiment_id"),
    "signature_table": ("study_id", "experiment_id", "signature_id"),
    "reported_taxa": (),  # Repeated observations are retained, not silently deduplicated.
    "ontology_term_table": ("study_id", "experiment_id", "slot", "term_id"),
}
_TABLE_MODELS: dict[str, type[BaseModel]] = {
    "study_table": StudyTableRow,
    "experiment_table": ExperimentTableRow,
    "signature_table": SignatureTableRow,
    "reported_taxa": TaxonOccurrenceRow,
    "ontology_term_table": ExperimentOntologyTermRow,
}


def read_table(
    path: Path, *, required: Sequence[str] = (), keys: Sequence[str] = ()
) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    """Read a rectangular TSV without DictReader's duplicate/extra-field coercions."""
    with open_input_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = tuple(next(reader, ()))
        if not header or any(not name or name != name.strip() for name in header):
            raise ValueError(f"{path}: missing or whitespace-padded column name")
        if len(set(header)) != len(header):
            raise ValueError(f"{path}: duplicate column names")
        missing = set(required).union(keys).difference(header)
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        records: list[dict[str, str]] = []
        seen: set[tuple[str, ...]] = set()
        for values in reader:
            if len(values) != len(header):
                raise ValueError(f"{path}:{reader.line_num}: wrong row width")
            row = dict(zip(header, values, strict=True))
            if keys:
                key = tuple(row[name].strip() for name in keys)
                if any(not value for value in key):
                    raise ValueError(f"{path}:{reader.line_num}: missing key identifier")
                if key in seen:
                    raise ValueError(f"{path}:{reader.line_num}: duplicate key {key!r}")
                seen.add(key)
            records.append(row)
    return header, records


def _path(context: Mapping[str, object], *names: str) -> Path | None:
    for name in names:
        value = context.get(name)
        if isinstance(value, Path):
            return value
    return None


def _records(context: Mapping[str, object], *names: str) -> list[dict[str, str]] | None:
    path = _path(context, *names)
    return read_table(path)[1] if path is not None else None


def _qualified(value: str, *, label: str) -> None:
    if not _CURIE.fullmatch(value):
        raise ValueError(f"{label}: expected a qualified identifier, got {value!r}")


def _feature_id(value: str) -> None:
    parts = value.split("|")
    if len(parts) not in {1, 2} or not parts[0].startswith("METPO:"):
        raise ValueError(f"Invalid native trait feature ID: {value!r}")
    for part in parts:
        _qualified(part, label="Feature")


def _check_foreign_key(
    rows: Sequence[Mapping[str, str]],
    parent: Sequence[Mapping[str, str]] | None,
    columns: tuple[str, ...],
    *,
    label: str,
) -> None:
    if parent is None:
        return
    parent_keys = {tuple(row[column] for column in columns) for row in parent}
    for row in rows:
        key = tuple(row[column] for column in columns)
        if key not in parent_keys:
            raise ValueError(f"{label}: missing parent for {key!r}")


def _table_contract(path: Path, kind: str, context: Mapping[str, object]) -> dict[str, object]:
    row_model = _TABLE_MODELS[kind]
    expected_header = tuple(row_model.model_fields)
    header, rows = read_table(path, required=expected_header, keys=_TABLE_KEYS[kind])
    if header != expected_header:
        raise ValueError(f"{path}: columns differ from the generated {row_model.__name__} contract")
    identifiers = [
        name for name in header if name in {"study_id", "experiment_id", "signature_id", "taxon_id"}
    ]
    for number, row in enumerate(rows, start=2):
        for column in identifiers:
            _qualified(row[column], label=f"{path}:{number} {column}")
        for column in header:
            if column.endswith("_count"):
                parse_nonnegative_integer(
                    row[column], path=path, line_number=number, column_name=column
                )
        if "direction" in row and row["direction"] not in {"increased", "decreased"}:
            raise ValueError(f"{path}:{number}: invalid signature direction {row['direction']!r}")
        if "taxon_id" in row and not row["taxon_id"].startswith("NCBITaxon:"):
            raise ValueError(f"{path}:{number}: expected NCBITaxon ID")
        if "term_id" in row:
            _qualified(row["term_id"], label="Ontology term")
        row_model.model_validate(row)
    studies = _records(context, "study_table")
    experiments = _records(context, "experiments")
    signatures = _records(context, "signatures")
    taxa = _records(context, "reported_taxa")
    if kind != "study_table":
        _check_foreign_key(rows, studies, ("study_id",), label=kind)
    if kind in {"signature_table", "reported_taxa", "ontology_term_table"}:
        _check_foreign_key(rows, experiments, ("study_id", "experiment_id"), label=kind)
    if kind == "reported_taxa":
        _check_foreign_key(
            rows, signatures, ("study_id", "experiment_id", "signature_id", "direction"), label=kind
        )
    child_sets = {
        "experiment_count": experiments,
        "signature_count": signatures,
        "taxon_count": taxa,
    }
    grouping = _TABLE_KEYS[kind]
    if kind in {"study_table", "experiment_table", "signature_table"}:
        for column, children in child_sets.items():
            if column not in header or children is None:
                continue
            counts = Counter(tuple(row[key] for key in grouping) for row in children)
            for row in rows:
                key = tuple(row[name] for name in grouping)
                if int(row[column]) != counts[key]:
                    raise ValueError(f"{path}: {column} does not reconcile for {key!r}")
    return {
        "rows": len(rows),
        "columns": len(header),
        "duplicate_observations_retained": kind == "reported_taxa",
    }


def _feature_contract(
    path: Path, kind: str, context: Mapping[str, object], strain_rollup: StrainRollup
) -> dict[str, object]:
    required = FEATURE_CATALOG_COLUMNS if kind == "feature_catalog" else FEATURE_METADATA_COLUMNS
    _, rows = read_table(path, required=required, keys=("column_id",))
    for row in rows:
        feature = row["column_id"]
        _feature_id(feature)
        expected = (
            f"{row['predicate_id']}|{row['object_id']}" if row["predicate_id"] else row["object_id"]
        )
        if feature != expected:
            raise ValueError(f"Feature ID does not match predicate/object: {feature!r}")
        if kind == "feature_catalog":
            FeatureDefinition.model_validate(row)
            if row["feature_class"] != feature_class(feature):
                raise ValueError(f"Incorrect MediaDive classification for {feature!r}")
            if row["annotation_evidence"] != annotation_evidence(feature, strain_rollup):
                raise ValueError(
                    "Aggregate feature catalog evidence must be the unresolved origins the "
                    f"strain rollup scope {strain_rollup!r} admits for {feature!r}"
                )
    features = tuple(row["column_id"] for row in rows)
    matrix = _path(context, "taxa_matrix", "taxon_trait_matrix")
    if matrix is not None and taxa_trait_count_reader(matrix).count_columns != features:
        raise ValueError("Feature catalog vocabulary/order does not match matrix columns")
    nodes = _path(context, "kg_nodes")
    if nodes is not None:
        validate_feature_node_metadata(path, nodes)
    metadata = _path(context, "feature_metadata")
    if kind == "feature_catalog" and metadata is not None:
        _, source_rows = read_table(metadata, required=FEATURE_METADATA_COLUMNS)
        observed = [tuple(row[column] for column in FEATURE_METADATA_COLUMNS) for row in rows]
        expected_metadata = [
            tuple(row[column] for column in FEATURE_METADATA_COLUMNS) for row in source_rows
        ]
        if observed != expected_metadata:
            raise ValueError("Feature catalog metadata differs from its source feature metadata")
    return {
        "rows": len(rows),
        "feature_filter_id": feature_filter_id(features),
        "culture_medium_features": sum(
            feature_class(feature) == "culture_medium" for feature in features
        ),
        "annotation_evidence": product_evidence(features, strain_rollup),
        "strain_attributed_features": sum(
            annotation_evidence(feature, strain_rollup) == STRAIN_ATTRIBUTED_EVIDENCE
            for feature in features
        ),
        "strain_rollup": strain_rollup,
        "kg_node_metadata_checked": nodes is not None,
        "source_metadata_checked": kind == "feature_catalog" and metadata is not None,
    }


def _matrix_contract(
    path: Path, context: Mapping[str, object], strain_rollup: StrainRollup
) -> dict[str, object]:
    reader = taxa_trait_count_reader(path)
    for feature in reader.count_columns:
        _feature_id(feature)
    ids: set[str] = set()
    zero_rows = 0
    for row in reader:
        _qualified(row.key[0], label="Matrix taxon")
        if not row.key[0].startswith("NCBITaxon:"):
            raise ValueError(f"Matrix requires NCBITaxon IDs: {row.key[0]!r}")
        ids.add(row.key[0])
        zero_rows += not any(row.counts)
    taxa = _records(context, "reported_taxa")
    if taxa is not None:
        expected = {
            row["taxon_id"] for row in taxa if row["direction"] in {"increased", "decreased"}
        }
        if ids != expected:
            raise ValueError(
                "Matrix taxon vocabulary differs from distinct reported directional IDs"
            )
    catalog = _path(context, "feature_catalog", "feature_metadata")
    if catalog is not None:
        _, records = read_table(catalog, required=("column_id",), keys=("column_id",))
        if tuple(row["column_id"] for row in records) != reader.count_columns:
            raise ValueError("Matrix feature vocabulary/order differs from catalog")
    source_facts: dict[str, object] = {"kg_source_aggregation_checked": False}
    edges = _path(context, "kg_edges")
    if edges is not None:
        reported = _path(context, "reported_taxa")
        nodes = _path(context, "kg_nodes")
        if reported is None or nodes is None:
            raise ValueError("KG matrix reconciliation requires reported taxa and KG nodes")
        source_facts = validate_matrix_sources(
            path,
            reported,
            edges,
            nodes,
            _path(context, "feature_metadata"),
            strain_rollup=strain_rollup,
        )
    return {
        "rows": reader.summary.rows,
        "features": len(reader.count_columns),
        "all_zero_rows": zero_rows,
        "feature_filter_id": feature_filter_id(reader.count_columns),
        "unit": "distinct_annotated_source_taxa_at_or_below_reported_taxon",
        "annotation_evidence": product_evidence(reader.count_columns, strain_rollup),
        "strain_rollup": strain_rollup,
        **source_facts,
    }


def _profile_taxa(context: Mapping[str, object], cohort: str) -> dict[ProfileKey, list[str]] | None:
    taxa = _records(context, "reported_taxa")
    if taxa is None:
        return None
    experiments = _records(context, "experiments")
    if experiments is None:
        raise ValueError("Cohort validation requires the experiments product")
    if experiments and "sequencing_type" not in experiments[0]:
        raise ValueError(
            "Cohort validation requires sequencing_type in the experiments product; "
            "a core made before #742 lacks it and must be rebuilt (#762)"
        )
    # Profiles keep 16S and shotgun experiments only, whatever the host cohort (#742)
    allowed = {
        (r["study_id"], r["experiment_id"])
        for r in experiments
        if r["sequencing_type"].strip() in CORPUS_SEQUENCING_TYPES
        and (cohort == "all" or r["host_species"].strip() == cohort)
    }
    result: dict[ProfileKey, list[str]] = {}
    # Signature records preserve valid profiles with zero reported taxa.
    for row in _records(context, "signatures") or taxa:
        direction = direction_to_inc_dec(row["direction"])
        if direction and (row["study_id"], row["experiment_id"]) in allowed:
            result.setdefault((row["study_id"], row["experiment_id"], direction), [])
    for row in taxa:
        direction = direction_to_inc_dec(row["direction"])
        key = (row["study_id"], row["experiment_id"], direction or "")
        if direction and key in result:
            result[key].append(row["taxon_id"])
    return result


def _profile_contract(
    path: Path, kind: str, context: Mapping[str, object], cohort: str
) -> dict[str, object]:
    reader = profile_count_reader(path)
    ora = kind == "ora_profiles"
    if ora and (not reader.count_columns or reader.count_columns[0] != "taxa_count"):
        raise ValueError("ORA profiles require taxa_count immediately after profile identifiers")
    if not ora and "taxa_count" in reader.count_columns:
        raise ValueError("Annotation count profiles must not contain an ORA denominator")
    features = reader.count_columns[1:] if ora else reader.count_columns
    for feature in features:
        _feature_id(feature)
    matrix_path = _path(context, "taxa_matrix", "taxon_trait_matrix")
    matrix = read_taxa_trait_matrix(matrix_path) if matrix_path is not None else None
    if matrix is not None and tuple(matrix.trait_columns) != features:
        raise ValueError("Profile feature vocabulary/order differs from taxon matrix")
    source = _profile_taxa(context, cohort)
    expected_coverage = (
        build_profile_coverage(source, matrix)
        if source is not None and matrix is not None
        else None
    )
    coverage_path = _path(
        context,
        "ora_coverage" if ora else "count_coverage",
        "profile_coverage_ora" if ora else "profile_coverage_count",
    )
    coverage = read_profile_coverage(coverage_path) if coverage_path is not None else None
    seen: set[ProfileKey] = set()
    filter_id = feature_filter_id(features)
    for row in reader:
        key = (row.key[0], row.key[1], row.key[2])
        if key[2] not in {"inc", "dec"}:
            raise ValueError(f"Invalid profile direction: {key!r}")
        _qualified(key[0], label="Profile study")
        _qualified(key[1], label="Profile experiment")
        seen.add(key)
        values = row.counts[1:] if ora else row.counts
        logical: dict[str, object] = {
            "study_id": key[0],
            "experiment_id": key[1],
            "inc_dec": key[2],
            "trait_filter_id": filter_id,
            # The complete-count reader has already checked every zero cell.
            "feature_counts": [
                {"feature_id": feature, "count": count}
                for feature, count in zip(features, values, strict=True)
                if count
            ],
        }
        if ora:
            logical["sample_size"] = row.counts[0]
            ORAProfileRow.model_validate(logical)
        else:
            CountProfileRow.model_validate(logical)
        if ora and any(value > row.counts[0] for value in values):
            raise ValueError(f"ORA trait count exceeds matched-taxon denominator for {key!r}")
        if source is not None and matrix is not None:
            if key not in source:
                raise ValueError(f"Unexpected profile outside declared cohort: {key!r}")
            matched = set(source[key]).intersection(matrix.counts_by_taxon)
            expected = tuple(
                sum(
                    (matrix.counts_by_taxon[taxon].get(feature, 0) > 0)
                    if ora
                    else matrix.counts_by_taxon[taxon].get(feature, 0)
                    for taxon in matched
                )
                for feature in features
            )
            if values != expected:
                raise ValueError(f"{kind} values do not match source aggregation for {key!r}")
            if ora and row.counts[0] != len(matched):
                raise ValueError(f"ORA denominator is not distinct lookup-matched taxa for {key!r}")
        if coverage is not None:
            if key not in coverage:
                raise ValueError(f"Missing coverage for profile {key!r}")
            if coverage[key].trait_filter_id != feature_filter_id(features):
                raise ValueError(f"Coverage feature-filter identity mismatch for {key!r}")
            if ora and coverage[key].lookup_matched_taxa != row.counts[0]:
                raise ValueError(f"Coverage and ORA denominator disagree for {key!r}")
    if source is not None and seen != set(source):
        raise ValueError("Profile membership differs from declared source cohort")
    if coverage is not None and set(coverage) != seen:
        raise ValueError("Coverage and profile membership differ")
    if coverage is not None and expected_coverage is not None and coverage != expected_coverage:
        raise ValueError("Coverage counts or identifier lists differ from source annotations")
    return {
        "rows": len(seen),
        "features": len(features),
        "cohort": cohort,
        "feature_filter_id": feature_filter_id(features),
        "unit": "distinct_lookup_matched_reported_taxa_with_feature"
        if ora
        else "summed_annotation_counts",
        "denominator": "distinct_lookup_matched_reported_taxa_including_unannotated"
        if ora
        else "none",
        "source_aggregation_checked": source is not None and matrix is not None,
        "mixed_rank_overlap_policy": "retain_distinct_reported_ids",
    }


def _coverage_contract(
    path: Path, context: Mapping[str, object], cohort: str | None
) -> dict[str, object]:
    coverage = read_profile_coverage(path)
    for row in read_table(path)[1]:
        ProfileCoverageRecord.model_validate(row)
    filters = {record.trait_filter_id for record in coverage.values()}
    if len(filters) > 1:
        raise ValueError("Coverage table mixes feature-filter identities")
    ora_path = _path(context, "ora_coverage", "profile_coverage_ora")
    is_ora = ora_path is not None and ora_path.resolve() == path.resolve()
    effective_cohort = cohort or ("Homo sapiens" if is_ora else "all")
    matrix_path = _path(context, "taxa_matrix", "taxon_trait_matrix")
    source = _profile_taxa(context, effective_cohort)
    if matrix_path is not None:
        matrix = read_taxa_trait_matrix(matrix_path)
        if filters and filters != {feature_filter_id(matrix.trait_columns)}:
            raise ValueError("Coverage feature-filter identity differs from matrix")
        if source is not None and coverage != build_profile_coverage(source, matrix):
            raise ValueError("Coverage does not reconcile to reported taxa and matrix annotations")
    return {
        "rows": len(coverage),
        "cohort": effective_cohort,
        "feature_filter_ids": sorted(filters),
    }


def _taxon_profiles(path: Path, kind: str, context: Mapping[str, object]) -> dict[str, object]:
    if kind == "taxon_profile_wide":
        reader = profile_count_reader(path)
        taxa = _records(context, "reported_taxa")
        expected: dict[tuple[str, ...], Counter[str]] = {}
        if taxa is not None:
            for record in taxa:
                key = (
                    record["study_id"],
                    record["experiment_id"],
                    direction_to_inc_dec(record["direction"]) or "",
                )
                feature = record["taxon_name"].strip() or record["taxon_id"].strip()
                expected.setdefault(key, Counter())[feature] += 1
            vocabulary = sorted({feature for counts in expected.values() for feature in counts})
            if tuple(vocabulary) != reader.count_columns:
                raise ValueError("Taxon name vocabulary differs from reported observations")
        seen_profiles: set[tuple[str, ...]] = set()
        for row in reader:
            if row.key[2] not in {"inc", "dec"}:
                raise ValueError("Invalid taxon profile direction")
            seen_profiles.add(row.key)
            if taxa is not None and (
                row.key not in expected
                or row.counts
                != tuple(expected[row.key][feature] for feature in reader.count_columns)
            ):
                raise ValueError("Taxon profile counts differ from reported observations")
        if taxa is not None and seen_profiles != set(expected):
            raise ValueError("Taxon profile membership differs from reported observations")
        return {
            "rows": reader.summary.rows,
            "features": len(reader.count_columns),
            "unit": "reported_taxon_observations",
            "source_aggregation_checked": taxa is not None,
        }
    reader = CountTableReader(path, leading_columns=("taxon",), key_columns=("taxon",))
    wide_path = _path(context, "taxon_profiles_wide")
    wide = profile_count_reader(wide_path) if wide_path is not None else None
    wide_rows = list(wide) if wide is not None else None
    if wide_rows is not None and reader.count_columns != tuple(
        "_".join(row.key) for row in wide_rows
    ):
        raise ValueError("Transposed taxon profile column identity differs from wide table")
    seen: list[str] = []
    for index, row in enumerate(reader):
        seen.append(row.key[0])
        if wide is not None and wide_rows is not None:
            if index >= len(wide.count_columns) or row.key[0] != wide.count_columns[index]:
                raise ValueError(
                    "Transposed taxon profile feature identity differs from wide table"
                )
            if row.counts != tuple(record.counts[index] for record in wide_rows):
                raise ValueError("Transposed taxon profile values differ from wide table")
    if wide is not None and tuple(seen) != wide.count_columns:
        raise ValueError("Transposed taxon profile feature membership differs from wide table")
    return {
        "rows": len(seen),
        "profiles": len(reader.count_columns),
        "layout": "transpose_of_wide_not_tidy_long",
    }


def validate_product(
    path: Path,
    kind: str,
    *,
    context: Mapping[str, object] | None = None,
    cohort: str | None = None,
    strain_rollup: StrainRollup = "off",
) -> dict[str, object]:
    """Validate a scientific product, raising ValueError before publication on failure.

    ``context`` maps workflow product IDs to resolved paths. When related inputs
    are supplied, keys, vocabulary, denominators, and aggregation are checked
    against them. The returned facts explicitly state whether source aggregation
    was available. ``cohort='all'`` means no host filter; other values match the
    experiments' host_species field exactly. Either way the cohort keeps only
    experiments whose sequencing_type is one of ``CORPUS_SEQUENCING_TYPES`` (#742).
    ``strain_rollup`` is the scope the producing run declared: the matrix and its feature
    catalog are checked against the BacDive strain records that scope admits (#788).
    """
    kind = _ALIASES.get(kind, kind)
    related = context or {}
    if cohort is not None:
        cohort = cohort.strip()
        if not cohort:
            raise ValueError("A declared cohort must be nonempty; use 'all' for no host filter")
    if kind in _TABLE_KEYS:
        facts = _table_contract(path, kind, related)
    elif kind in {"feature_metadata", "feature_catalog"}:
        facts = _feature_contract(path, kind, related, strain_rollup)
    elif kind == "taxa_matrix":
        facts = _matrix_contract(path, related, strain_rollup)
    elif kind in {"count_profiles", "ora_profiles"}:
        effective_cohort = cohort or ("Homo sapiens" if kind == "ora_profiles" else "all")
        facts = _profile_contract(path, kind, related, effective_cohort)
    elif kind == "profile_coverage":
        facts = _coverage_contract(path, related, cohort)
    elif kind in {"taxon_profile_wide", "taxon_profile_long"}:
        facts = _taxon_profiles(path, kind, related)
    else:
        raise ValueError(f"Unknown scientific product contract: {kind!r}")
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"contract_version": CONTRACT_VERSION, "kind": kind, "sha256": digest, **facts}
