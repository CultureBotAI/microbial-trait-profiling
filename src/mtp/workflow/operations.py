"""Explicit-path raw-to-core operations and source tracing for issue #629.

The workflow graph lives in workflows/bugsigdb_core.yaml. These adapters call
the native builders, retaining their numerical definitions and host defaults.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from mtp.datamodel.core_products_pydantic import SourceIdentity
from mtp.download import (
    BUGSIGDB_EXPORTS_DOI,
    BUGSIGDB_EXPORTS_TAG_SHA,
    BUGSIGDB_EXPORTS_VERSION,
    PINNED_FULL_DUMP_BYTES,
    PINNED_FULL_DUMP_URL,
    IntegrityError,
    verify_pinned_dump,
)
from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES
from mtp.ingest import ingest
from mtp.paths import file_sha256
from mtp.study_trait_profiles import (
    ProfileKey,
    collect_profile_taxa,
    create_study_trait_profiles,
    read_taxa_trait_matrix,
)
from mtp.study_tsv import GeneratedTables, generate_tables
from mtp.taxa_trait_matrix import (
    STRAIN_PREFIX,
    SUBCLASS_OF,
    StrainRollup,
    assertion_trait_key,
    create_taxa_trait_matrix,
    parse_strain_rollup,
    read_assertion_term_metadata,
    read_node_metadata,
    scoped_taxon_traits,
    strain_trait_rolls_up,
)
from mtp.workflow.coverage import (
    build_profile_coverage,
    write_feature_catalog,
    write_profile_coverage,
)
from mtp.workflow.kg_inputs import audit_kg_edges, iter_kg_edges
from mtp.workflow.models import (
    ProductDefinition,
    ValidationResult,
    ValidationStatus,
    WorkflowRun,
)
from mtp.workflow.population import select_contrasts, validate_population
from mtp.workflow.source_contracts import validate_export_sources, validate_study_sources
from mtp.workflow.source_identity import validate_source_identity

Operation = Callable[[dict[str, Path], dict[str, Path], dict[str, object]], None]
RAW_COLUMNS = ("BSDB ID", "NCBI Taxonomy IDs", "MetaPhlAn taxon names", "Abundance in Group 1")


@dataclass(frozen=True)
class _ScientificSemantics:
    unit: str
    cohort: str | None = None
    feature_policy: str | None = None
    coverage_policy: str | None = None
    schema_ref: str | None = None


# These are implementation contracts, not a second workflow graph. A declarative
# label cannot change the numerical meaning implemented by the native builders.
_SCHEMA = "schema/core_products.yaml#"
SOURCE_EXPORTS = (
    "study_table",
    "experiments",
    "experiment_ontology_terms",
    "signatures",
    "reported_taxa",
    "taxon_profiles_wide",
    "taxon_profiles_long",
)
_CORE_SEMANTICS = {
    "raw_csv": _ScientificSemantics("literature_signature_record", "as_supplied"),
    "kg_edges": _ScientificSemantics("knowledge_graph_edge", "as_supplied"),
    "kg_nodes": _ScientificSemantics("knowledge_graph_node", "as_supplied"),
    "source_identity": _ScientificSemantics(
        "consumed_source_file", "as_supplied", schema_ref=_SCHEMA + "SourceIdentity"
    ),
    "study_collection": _ScientificSemantics(
        "study",
        "all_hosts",
        coverage_policy="exact_current_ingest_collection",
        schema_ref="schema/mtp.yaml#Study",
    ),
    "contrast_study_collection": _ScientificSemantics(
        "study",
        coverage_policy="contrast_population_v2",
        schema_ref="schema/mtp.yaml#Study",
    ),
    "contrast_selection": _ScientificSemantics(
        "experiment_selection_decision", coverage_policy="contrast_population_v2"
    ),
    "ingest_log": _ScientificSemantics("ingest_event", "all_hosts"),
    "study_table": _ScientificSemantics("study", "all_hosts", schema_ref=_SCHEMA + "StudyTableRow"),
    "experiment_table": _ScientificSemantics(
        "experiment", "all_hosts", schema_ref=_SCHEMA + "ExperimentTableRow"
    ),
    "ontology_term_table": _ScientificSemantics(
        "experiment_ontology_annotation",
        "all_hosts",
        schema_ref=_SCHEMA + "ExperimentOntologyTermRow",
    ),
    "signature_table": _ScientificSemantics(
        "literature_signature", "all_hosts", schema_ref=_SCHEMA + "SignatureTableRow"
    ),
    "reported_taxa": _ScientificSemantics(
        "reported_taxon_occurrence",
        "all_hosts",
        coverage_policy="preserve_repeated_occurrences",
        schema_ref=_SCHEMA + "TaxonOccurrenceRow",
    ),
    "taxon_profile_wide": _ScientificSemantics("taxon_name_occurrence", "all_hosts"),
    "taxon_profile_long": _ScientificSemantics("taxon_name_occurrence", "all_hosts"),
    "taxa_matrix": _ScientificSemantics(
        "trait_bearing_source_taxa_at_or_below_reported_taxon",
        "all_reported_taxa",
        "include_mediadive",
        "all_reported_taxa_including_zero_annotation_rows",
    ),
    "feature_metadata": _ScientificSemantics(
        "annotation_feature", "all_reported_taxa", "include_mediadive"
    ),
    "feature_catalog": _ScientificSemantics(
        "annotation_feature",
        "all_reported_taxa",
        "include_mediadive",
        schema_ref=_SCHEMA + "FeatureDefinition",
    ),
    "count_profiles": _ScientificSemantics(
        "rolled_up_source_taxon_count",
        feature_policy="include_mediadive",
        coverage_policy="count_input_16s_wms_v1",
        schema_ref=_SCHEMA + "CountProfileRow",
    ),
    "ora_profiles": _ScientificSemantics(
        "distinct_matched_reported_taxon_count",
        feature_policy="include_mediadive",
        coverage_policy="count_input_16s_wms_v1",
        schema_ref=_SCHEMA + "ORAProfileRow",
    ),
    "profile_coverage": _ScientificSemantics(
        "distinct_reported_taxon_lookup",
        feature_policy="include_mediadive",
        coverage_policy="count_input_16s_wms_v1",
        schema_ref=_SCHEMA + "ProfileCoverageRecord",
    ),
}


def _validate_semantic_descriptor(product: ProductDefinition) -> None:
    expected = _CORE_SEMANTICS.get(product.kind)
    if expected is None:
        raise ValueError(f"{product.id}: no implemented scientific contract for {product.kind}")
    for field in ("unit", "feature_policy", "coverage_policy", "schema_ref", "cohort"):
        value = getattr(expected, field)
        if field == "cohort" and value is None:
            # Runtime resolves host overrides before validation. Row recomputation
            # below independently checks that the data match the effective cohort.
            if not product.cohort or not product.cohort.strip():
                raise ValueError(
                    f"{product.id}: declared cohort must name the effective host policy"
                )
        elif getattr(product, field) != value:
            raise ValueError(
                f"{product.id}: declared {field} does not match implemented "
                f"{product.kind} semantics (expected {value!r})"
            )


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _host(parameters: Mapping[str, object], name: str) -> str | None:
    value = str(parameters[name]).strip()
    if not value:
        raise ValueError(f"{name} cannot be empty; use 'all' for an unfiltered cohort")
    return None if value == "all" else value


def identify_sources(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    """Describe consumed bytes; a release label alone never establishes the pin."""
    if set(parameters) != {"source_mode"}:
        raise ValueError("identify_sources requires exactly the source_mode parameter")
    declared = str(parameters["source_mode"])
    if declared not in {"unknown", "fixture", "live", "pinned"}:
        raise ValueError("source_mode must be unknown, fixture, live, or pinned")
    raw = inputs["raw_dump"]
    pinned = False
    if raw.stat().st_size == PINNED_FULL_DUMP_BYTES:
        try:
            verify_pinned_dump(raw.read_bytes())
        except IntegrityError:
            pinned = False
        else:
            pinned = True
    if declared == "pinned" and not pinned:
        raise ValueError("Declared pinned BugSigDB source does not match verified release bytes")
    with raw.open(encoding="utf-8") as handle:
        first_line = handle.readline().rstrip("\n")
    facts: dict[str, object] = {
        "recorded_at_utc": datetime.now(UTC).isoformat(),
        "declared_source_mode": declared,
        "pinned_verified": pinned,
        "verification_status": "pinned_bytes_verified" if pinned else "origin_unverified",
        "snapshot_header": first_line if first_line.startswith("#") else None,
        "consumed_files": {
            name: {
                "product_id": name,
                "path": str(path.resolve()),
                "sha256": file_sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for name, path in sorted(inputs.items())
        },
        "interpretation": "Taxon reports and organism annotations; no physiological measurement.",
    }
    if pinned:
        facts["verified_release"] = {
            "version": BUGSIGDB_EXPORTS_VERSION,
            "doi": BUGSIGDB_EXPORTS_DOI,
            "commit": BUGSIGDB_EXPORTS_TAG_SHA,
            "url": PINNED_FULL_DUMP_URL,
        }
    SourceIdentity.model_validate(facts)
    _write_json(outputs["source_identity"], facts)


def ingest_operation(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    """Ingest into an empty attempt-owned collection, never an existing study set."""
    if parameters:
        raise ValueError("ingest has no configurable scientific parameters")
    studies = outputs["studies"]
    if studies.exists() and any(studies.iterdir()):
        raise ValueError("Ingest collection must be empty; refusing to reuse stale study YAML")
    ingest(
        raw_dir=inputs["raw_dump"].parent,
        study_dir=studies,
        derived_dir=outputs["ingest_log"].parent,
        full_dump_path=inputs["raw_dump"],
        log_path=outputs["ingest_log"],
    )


def export_tables_operation(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    """Flatten only this execution's study collection to declared filenames."""
    if parameters:
        raise ValueError("export_tables has no configurable scientific parameters")
    paths = GeneratedTables(
        studies=outputs["study_table"],
        experiments=outputs["experiments"],
        experiment_ontology_terms=outputs["experiment_ontology_terms"],
        signatures=outputs["signatures"],
        taxa=outputs["reported_taxa"],
        taxon_profiles_wide=outputs["taxon_profiles_wide"],
        taxon_profiles_long=outputs["taxon_profiles_long"],
    )
    generate_tables(inputs["studies"], outputs["study_table"].parent, output_paths=paths)


def ingest_source_operation(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    """Retain the complete native ingest before the explicit selection step."""
    ingest_operation(
        inputs,
        {"studies": outputs["source_studies"], "ingest_log": outputs["ingest_log"]},
        parameters,
    )


def export_source_tables_operation(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    """Keep full-source exports for provenance and the fixed ORA background."""
    export_tables_operation(
        {"studies": inputs["source_studies"]},
        {name: outputs[f"source_{name}"] for name in SOURCE_EXPORTS},
        parameters,
    )


def source_taxa_matrix_operation(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    """Build the same full-source taxon matrix regardless of contrast selection."""
    taxa_matrix_operation(
        {
            "reported_taxa": inputs["source_reported_taxa"],
            "kg_edges": inputs["kg_edges"],
            "kg_nodes": inputs["kg_nodes"],
        },
        outputs,
        parameters,
    )


def taxa_matrix_operation(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    """Use native direct-or-descendant source-taxon counting without feature filtering.

    BacDive strain records count at their linked NCBITaxon nodes as far as the declared
    strain_rollup scope admits (#788).
    """
    if set(parameters) != {"strain_rollup"}:
        raise ValueError("taxa_matrix requires exactly the strain_rollup parameter")
    create_taxa_trait_matrix(
        taxa_path=inputs["reported_taxa"],
        edges_path=inputs["kg_edges"],
        nodes_path=inputs["kg_nodes"],
        output_path=outputs["taxa_matrix"],
        column_metadata_path=outputs["feature_metadata"],
        strain_rollup=parse_strain_rollup(parameters["strain_rollup"]),
    )


def feature_catalog_operation(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    """Classify medium annotations explicitly while preserving every matrix feature.

    The evidence column marks the features whose counts the declared strain_rollup scope
    lets strain records enter (#788).
    """
    if set(parameters) != {"strain_rollup"}:
        raise ValueError("feature_catalog requires exactly the strain_rollup parameter")
    outputs["feature_catalog"].parent.mkdir(parents=True, exist_ok=True)
    write_feature_catalog(
        inputs["feature_metadata"],
        outputs["feature_catalog"],
        strain_rollup=parse_strain_rollup(parameters["strain_rollup"]),
    )


def _profiles_operation(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object], *, ora: bool
) -> None:
    prefix = "ora" if ora else "count"
    if set(parameters) != {f"{prefix}_host_species"}:
        raise ValueError(f"{prefix}_profiles requires exactly the {prefix}_host_species parameter")
    host = _host(parameters, f"{prefix}_host_species")
    # Both profile products keep 16S and shotgun experiments only (#742)
    create_study_trait_profiles(
        study_dir=inputs["studies"],
        taxa_trait_matrix_path=inputs["taxa_matrix"],
        output_path=outputs[f"{prefix}_profiles"],
        host_species=host,
        binary=ora,
        include_taxa_count=ora,
        deduplicate_taxa=True,
        sequencing_types=CORPUS_SEQUENCING_TYPES,
    )
    matrix = read_taxa_trait_matrix(inputs["taxa_matrix"])
    profile_taxa, _ = collect_profile_taxa(
        inputs["studies"], host_species=host, sequencing_types=CORPUS_SEQUENCING_TYPES
    )
    coverage = build_profile_coverage(profile_taxa, matrix)
    output = outputs[f"{prefix}_coverage"]
    output.parent.mkdir(parents=True, exist_ok=True)
    write_profile_coverage(output, coverage)


def count_profiles_operation(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    """Sum native source-taxon counts over distinct reported IDs, all hosts by default."""
    _profiles_operation(inputs, outputs, parameters, ora=False)


def ora_profiles_operation(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    """Count positive reported IDs with matched-ID N, human hosts by default."""
    _profiles_operation(inputs, outputs, parameters, ora=True)


OPERATIONS: dict[str, Operation] = {
    "identify_sources": identify_sources,
    "ingest": ingest_operation,
    "export_tables": export_tables_operation,
    "taxa_matrix": taxa_matrix_operation,
    "feature_catalog": feature_catalog_operation,
    "count_profiles": count_profiles_operation,
    "ora_profiles": ora_profiles_operation,
    "ingest_source": ingest_source_operation,
    "export_source_tables": export_source_tables_operation,
    "source_taxa_matrix": source_taxa_matrix_operation,
    "select_contrasts": select_contrasts,
}


def _validate_delimited(path: Path, required: tuple[str, ...], *, raw_csv: bool = False) -> int:
    with path.open(encoding="utf-8", newline="") as handle:
        lines = (line for line in handle if not raw_csv or not line.startswith("#"))
        reader = csv.reader(lines, delimiter="," if raw_csv else "\t")
        header = next(reader, [])
        if len(set(header)) != len(header) or not set(required).issubset(header):
            raise ValueError(f"{path}: missing or duplicate required columns {required}")
        if not raw_csv and tuple(header[: len(required)]) != required:
            raise ValueError(f"{path}: native KG parser requires leading columns {required}")
        rows = 0
        for line_number, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(f"{path}:{line_number}: row width differs from header")
            rows += 1
    if not rows:
        raise ValueError(f"{path}: empty input table")
    return rows


def run_strain_rollup(run: WorkflowRun) -> StrainRollup:
    """The strain rollup scope a core run recorded, ``off`` for a run that predates it (#788)."""
    recorded = {item.name: item.value for item in run.parameters or []}
    return parse_strain_rollup(recorded.get("strain_rollup", "off"))


def validate_core_product(
    product: ProductDefinition,
    path: Path,
    context: Mapping[str, Path],
    *,
    recorded_sources: Mapping[str, Path] | None = None,
    parameters: Mapping[str, str] | None = None,
) -> list[ValidationResult]:
    """Validate one scientific product using the execution's declared related products.

    ``parameters`` are the run's resolved parameters. The matrix and its feature catalog
    are checked under the declared strain_rollup scope, ``off`` for a run that predates
    the parameter (#788).
    """
    from mtp.workflow.contracts import validate_product

    strain_rollup = parse_strain_rollup((parameters or {}).get("strain_rollup", "off"))

    _validate_semantic_descriptor(product)
    kind = product.kind
    if product.id in {f"source_{name}" for name in SOURCE_EXPORTS}:
        context = {
            **context,
            "studies": context["source_studies"],
            **{name: context[f"source_{name}"] for name in SOURCE_EXPORTS},
        }
    elif kind == "taxa_matrix" and "source_reported_taxa" in context:
        context = {**context, "reported_taxa": context["source_reported_taxa"]}
    details: dict[str, object]
    if kind == "raw_csv":
        details = {"rows": _validate_delimited(path, RAW_COLUMNS, raw_csv=True)}
    elif kind == "kg_edges":
        details = audit_kg_edges(path)
    elif kind == "kg_nodes":
        details = {"rows": _validate_delimited(path, ("id", "category", "name"))}
    elif kind == "study_collection":
        details = validate_study_sources(context["raw_dump"], path, context["ingest_log"])
    elif kind == "source_identity":
        details = validate_source_identity(path, context, recorded_sources=recorded_sources)
    elif kind == "ingest_log":
        value = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(value, list) or any(not isinstance(item, dict) for item in value):
            raise ValueError("Ingest log must be a list of event mappings")
        details = {"events": len(value)}
    elif kind in {"contrast_study_collection", "contrast_selection"}:
        details = validate_population(
            context["studies"],
            context["contrast_selection"],
            context["source_studies"],
            context["raw_dump"],
            product.cohort or "",
        )
    else:
        details = validate_product(
            path, kind, context=context, cohort=product.cohort, strain_rollup=strain_rollup
        )
        if kind == "study_table":
            details.update(validate_export_sources(context["studies"], context))
    return [
        ValidationResult(
            validator=f"core:{kind}:v1",
            status=ValidationStatus.passed,
            message="Scientific product contract passed",
            details_json=json.dumps(details, sort_keys=True),
        )
    ]


def trace_profile_feature(
    study_dir: Path,
    matrix_path: Path,
    edges_path: Path,
    nodes_path: Path,
    profile_key: ProfileKey,
    feature_id: str,
    *,
    host_species: str | None = None,
    strain_rollup: StrainRollup = "off",
) -> dict[str, object]:
    """Trace a chosen profile feature to distinct KG source taxa and original edges.

    Only the selected feature's evidence edges are retained. Taxonomy links are
    traversed upward with cycle detection. Duplicate evidence edges do not add
    source taxa, and genus/species overlap remains visible rather than corrected.
    Under the matrix's strain_rollup scope, the BacDive strain edges it admits are
    evidence for the strains' linked NCBITaxon nodes, by the builder's rule (#788).
    """
    matrix = read_taxa_trait_matrix(matrix_path)
    if feature_id not in matrix.trait_columns:
        raise ValueError(f"Feature is absent from the matrix: {feature_id}")
    profiles, _ = collect_profile_taxa(
        study_dir, host_species=host_species, sequencing_types=CORPUS_SEQUENCING_TYPES
    )
    if profile_key not in profiles:
        raise ValueError(f"Profile is absent from selected study cohort: {profile_key}")
    reported = set(profiles[profile_key])
    matched = reported.intersection(matrix.counts_by_taxon)
    parents: defaultdict[str, set[str]] = defaultdict(set)
    direct: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    strain_taxa: defaultdict[str, set[str]] = defaultdict(set)
    strain_evidence: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    kg_taxa: set[str] = set()
    terms = read_assertion_term_metadata(nodes_path, edges_path)
    for row in iter_kg_edges(edges_path):
        subject, predicate, obj = (row[name].strip() for name in ("subject", "predicate", "object"))
        kg_taxa.update(item for item in (subject, obj) if item.startswith("NCBITaxon:"))
        is_link = predicate == SUBCLASS_OF and obj.startswith("NCBITaxon:")
        is_evidence = (
            subject.startswith(("NCBITaxon:", STRAIN_PREFIX))
            and (obj.startswith("METPO:") or predicate.startswith("METPO:"))
            and assertion_trait_key(row, terms) == feature_id
        )
        if subject.startswith("NCBITaxon:"):
            if is_link:
                parents[subject].add(obj)
            elif is_evidence:
                direct[subject].append(row)
        elif strain_rollup != "off" and subject.startswith(STRAIN_PREFIX):
            if is_link:
                strain_taxa[subject].add(obj)
            elif is_evidence:
                strain_evidence[subject].append(row)
    sources = scoped_taxon_traits(
        {source: {feature_id} for source in direct},
        dict(strain_taxa),
        {strain: {feature_id} for strain in strain_evidence},
        strain_rollup,
    )
    linked: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    if strain_trait_rolls_up(feature_id, strain_rollup):
        for strain, rows in sorted(strain_evidence.items()):
            for source in sorted(strain_taxa.get(strain, ())):
                linked[source].extend(rows)
    evidence = {source: [*direct.get(source, []), *linked.get(source, [])] for source in sources}
    needed = reported | set(evidence)
    metadata = read_node_metadata(nodes_path, needed)
    kg_taxa.update(metadata)
    ledger: list[dict[str, object]] = []
    source_multiplicity: Counter[str] = Counter()
    reconstructed: Counter[str] = Counter()
    for source, edges in sorted(evidence.items()):
        visited: set[str] = set()
        stack = [source]
        contributors: set[str] = set()
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            if current in matched:
                contributors.add(current)
            stack.extend(parents.get(current, ()))
        for reported_taxon in sorted(contributors):
            reconstructed[reported_taxon] += 1
            source_multiplicity[source] += 1
            ledger.append(
                {
                    "reported_taxon_id": reported_taxon,
                    "source_taxon_id": source,
                    "source_taxon_name": metadata[source].name if source in metadata else None,
                    "annotation_origin": "direct" if reported_taxon == source else "descendant",
                    "source_count_contribution": 1,
                    "evidence_edges": edges,
                }
            )
    for taxon in matched:
        if matrix.counts_by_taxon[taxon].get(feature_id, 0) != reconstructed[taxon]:
            raise ValueError(f"KG source ledger does not reconcile matrix count for {taxon}")
    count = sum(matrix.counts_by_taxon[t].get(feature_id, 0) for t in matched)
    ora = sum(matrix.counts_by_taxon[t].get(feature_id, 0) > 0 for t in matched)
    return {
        "profile_key": list(profile_key),
        "feature_id": feature_id,
        "count_profile_value": count,
        "ora_profile_value": ora,
        "ora_taxa_count": len(matched),
        "reported_taxon_ids": sorted(reported),
        "matrix_matched_taxon_ids": sorted(matched),
        "matrix_unmatched_taxon_ids": sorted(reported - matched),
        "reported_taxa_absent_from_kg": sorted(reported - kg_taxa),
        "unique_annotation_source_taxa": len(source_multiplicity),
        "source_taxon_multiplicity": dict(sorted(source_multiplicity.items())),
        "overlapping_source_taxa": sorted(t for t, n in source_multiplicity.items() if n > 1),
        "ledger": ledger,
        "interpretation": "Annotation alignment; count contributions overlap across reported ranks "
        "and do not measure cells, in vivo activity, or independent organisms.",
    }
