"""Inventory and compare committed LinkML records for #1316 (2026-11-08).

The inventory pass reads tracked path metadata and at most 4 KiB of each candidate
text stream; it does not import models or execute validation. The comparison pass
reports baseline/new schema and generated-model outcomes without data edits.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib
import importlib.util
import io
import json
import stat
import sys
import tarfile
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, BinaryIO, Protocol, TextIO, cast

import yaml
from git import Repo
from pydantic import BaseModel, ValidationError

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[2]
TEXT_ENDINGS = (".yaml", ".yml", ".json", ".jsonl", ".ndjson")
SPECIAL_JSON = {".secrets.baseline"}
HEADER_BYTES = 4096
FIXTURE_ROOTS = ("tests/fixtures/", "examples/", "workflows/")
MAX_MEMBER_BYTES = 8 * 1024 * 1024
MAX_ARCHIVE_BYTES = 64 * 1024 * 1024
MAX_ARCHIVE_MEMBERS = 1000
MAX_RECORD_CHARS = 1024 * 1024
MAX_COLUMNS = 16384
SCHEMAS = ("mtp", "bsdbng", "core_products", "workflow", "multiple_testing")
MODEL_PATHS = {name: f"src/mtp/datamodel/{name}_pydantic.py" for name in SCHEMAS}
MODEL_PATHS["bsdbng"] = "src/bsdbng/datamodel/bsdbng_pydantic.py"
CORE_TABLES = {
    "studies.tsv": "StudyTableRow",
    "experiments.tsv": "ExperimentTableRow",
    "signatures.tsv": "SignatureTableRow",
    "bugsigdb_taxa.tsv": "TaxonOccurrenceRow",
    "experiment_ontology_terms.tsv": "ExperimentOntologyTermRow",
    "feature_catalog.tsv": "FeatureDefinition",
    "study_trait_profiles_coverage.tsv": "ProfileCoverageRecord",
    "study_trait_profiles_ora_coverage.tsv": "ProfileCoverageRecord",
}
VISUAL_JSON = {
    "trait-visualization/trait-profile-investigation/PD-cellulose/pd_exploratory_taxa_six_traits_influence_heatmap_manifest.json",
    "trait-visualization/trait-profile-investigation/PD-cellulose/pd_six_trait_taxon_evidence_provenance.json",
    "trait-visualization/trait-profile-investigation/PD-cellulose/pd_species_six_traits_influence_heatmap_manifest.json",
    "trait-visualization/trait-profile-investigation/PD-cellulose/pd_trait_contributors_manifest.json",
    "trait-visualization/trait-profile-investigation/six-shared-traits/oral/oral_six_trait_provenance.json",
    "trait-visualization/trait-profile-investigation/six-shared-traits/six_shared_traits_profile_graph_data.json",
    "trait-visualization/trait-profile-investigation/stratified_stats_provenance.json",
    "trait-visualization/trait-profile-investigation/taxa-coverage/provenance.json",
}
FIXTURE_DOCS = {
    "tests/fixtures/cellulose_tufte_historical_sources/README.md",
    "tests/fixtures/cellulose_tufte_historical_sources/integrated_evidence_a3c32cf.md",
    "tests/fixtures/cellulose_tufte_historical_sources/methods_integrated_evidence_a3c32cf.md",
    "tests/fixtures/comparison_specs/README.md",
    "tests/fixtures/ora_figure_1505/README.md",
    "tests/fixtures/trait_key_policy/README.md",
    "tests/fixtures/trait_key_policy/historical_raw_core/README.md",
    "tests/fixtures/workflow/README.md",
}
PHYSICAL_TABLES = {
    "core/steps/count_profiles/derived/study_trait_profiles_wide.tsv",
    "core/steps/ora_profiles/derived/study_trait_profiles_ora_wide.tsv",
    "core/steps/export_tables/derived/taxon_profiles_long.tsv",
    "core/steps/export_tables/derived/taxon_profiles_wide.tsv",
    "core/steps/taxa_matrix/derived/taxa_trait_matrix.tsv",
    "core/steps/taxa_matrix/derived/taxa_trait_matrix_columns.tsv",
}
CORE_ROOTS = (
    "tests/fixtures/trait_key_policy/historical_raw_core/",
    "tests/fixtures/trait_key_policy/historical_union_core.tar.gz!",
)
NATIVE_TABLES = {
    "methods/table1-native-results/native_composition_comparison.tsv.gz",
    "methods/table1-native-results/native_ora_fold_comparison.tsv.gz",
    "methods/table1-native-results/native_table1_stats.tsv.gz",
}
CONFIG_FILES = {
    ".github/dependabot.yml",
    ".github/workflows/ci.yaml",
    ".github/workflows/codeql.yaml",
    ".github/workflows/release.yaml",
    ".pre-commit-config.yaml",
    ".yamllint.yaml",
    ".yamllint-historical.yaml",
    ".secrets.baseline",
    "schema/.linkmllint.yaml",
}
CUSTOM_JSON = {
    "data/analysis/cellulose_decr/cellulose_decr_provenance.json": (
        "src/mtp/cellulose_decr/analysis.py"
    ),
    "data/analysis/cellulose_decr/cellulose_decr_summary.json": (
        "src/mtp/cellulose_decr/analysis.py"
    ),
    "data/analysis/cellulose_decr/four_highlighted_pairs_summary.json": (
        "src/mtp/cellulose_decr/tufte_plates.py"
    ),
    "data/derived/calibration/composition/composition_calibration_summary.json": (
        "src/mtp/composition_calibration/composition.py"
    ),
    "data/derived/integrated_evidence/integrated_evidence_provenance.json": (
        "src/mtp/integrated_evidence.py"
    ),
    "data/derived/ora_calibration/calibration_summary.json": "src/mtp/ora_calibration/harness.py",
    "data/kg-microbe/provenance.json": "Justfile",
    "data/raw/download_provenance.json": "Justfile",
    "docs/reviews/2026-10-16-composition-archive-manifest.json": (
        "docs/reviews/2026-10-16-composition-archive-manifest.json"
    ),
    "docs/reviews/2026-10-19-ora-provenance-refresh-audit.json": (
        "scripts/adhoc/issue720_audit_ora_refresh.py"
    ),
    "docs/reviews/2026-10-19-ora-provenance-refresh-inputs.json": (
        "scripts/adhoc/issue720_audit_ora_refresh.py"
    ),
    "methods/calibration_protocol.json": "src/mtp/calibration.py",
    "methods/cellulose_decr_selection_scan_protocol.json": (
        "src/mtp/cellulose_decr/selection_scan.py"
    ),
    "methods/cellulose_decr_selection_scan_protocol_v2.json": (
        "src/mtp/cellulose_decr/selection_scan.py"
    ),
    "methods/study_split_selection_protocol.json": "src/mtp/cellulose_decr/study_split.py",
    "methods/table1-native-results/contribution_manifest.json": (
        "src/mtp/native_influence_publication.py"
    ),
    "methods/table1-native-results/publication_manifest.json": (
        "src/mtp/native_influence_publication.py"
    ),
    "methods/table1-native-results/run_manifest.json": "src/mtp/native_influence_publication.py",
    "methods/table1_influence_extension_protocol.json": "src/mtp/influence_extension.py",
    "methods/table1_log_ratio_sensitivity_protocol.json": "src/mtp/log_ratio_sensitivity.py",
    "methods/table1_native_influence_protocol.json": "src/mtp/native_taxon_influence.py",
    "methods/table1_taxon_influence_protocol.json": "src/mtp/taxon_influence.py",
    "methods/trait_screen_manifest.json": "src/mtp/trait_screen.py",
    "tests/fixtures/calibration/issue784_before_matrix_removal.json": (
        "tests/test_calibration_matrix_cleanup.py:31"
    ),
    "tests/fixtures/calibration_status_quick_run.json": "tests/test_integrated_evidence.py:1813",
    "tests/fixtures/cellulose_historical_numeric_snapshot.json": (
        "tests/test_cellulose_historical_numeric_snapshot.py:39"
    ),
    "tests/fixtures/cellulose_tufte_historical_sources/study_modal_7970cba.json": (
        "tests/test_cellulose_tufte_overlap_sources.py:27"
    ),
    "tests/fixtures/contrast_orientation_guard_digests.json": (
        "tests/test_contrast_orientation.py:2460"
    ),
    "tests/fixtures/contrast_orientation_provenance_3d9d5cd.json": (
        "tests/test_integrated_evidence.py:1827"
    ),
    "tests/fixtures/sensitivity_status_c9fb6f1.json": "tests/test_integrated_evidence.py:1783",
}


def classification(name: str) -> dict[str, Any]:
    """Explicit format mapping; unknown additions are never silently non-instances."""
    core_path = next(
        (name.removeprefix(prefix) for prefix in CORE_ROOTS if name.startswith(prefix)), None
    )
    leaf = Path(core_path if core_path is not None else name).name
    core = core_path is not None
    result: dict[str, Any] = {"kind": "unknown", "reason": "No reviewed format mapping"}

    def typed(schema: str, target: str, reader: str, encoding: str = "document") -> dict[str, Any]:
        return {
            "kind": "instance",
            "schema": schema,
            "class": target,
            "encoding": encoding,
            "reader": reader,
        }

    if name in {f"schema/{schema}.yaml" for schema in SCHEMAS}:
        return {"kind": "schema", "reader": "Justfile:993", "schema": Path(name).stem}
    if name.startswith("workflows/") and name.endswith(".yaml"):
        return typed("workflow", "WorkflowDefinition", "Justfile:55")
    if name == "methods/multiple_testing_policy_v1.yaml":
        return typed("multiple_testing", "MultipleTestingPolicy", "Justfile:59")
    if name == "tests/fixtures/multiple_testing_family.yaml":
        return typed("multiple_testing", "HypothesisFamily", "Justfile:60")
    if name == "tests/fixtures/comparison_specs/analysis_specs_f0cb5b7.json":
        return typed(
            "workflow",
            "AnalysisSpec",
            "tests/test_comparison_spec_fixture_extractor.py:16",
            "array",
        )
    if name.endswith("historical_union_core.tar.gz"):
        return {"kind": "archive", "reader": "tests/fixtures/trait_key_policy/README.md"}
    if core and name.endswith("core/run.json"):
        return typed("workflow", "WorkflowRun", "src/mtp/workflow/runtime.py:654")
    if core and leaf == "_step.json":
        return typed("workflow", "StepExecution", "src/mtp/workflow/runtime.py:719")
    if core and leaf == "source_identity.json":
        return typed("core_products", "SourceIdentity", "src/mtp/workflow/operations.py:246")
    if core and "/studies/" in name and leaf.endswith(".yaml"):
        result = typed("mtp", "Study", "Justfile:1024")
        result["also_schema"] = "bsdbng"
        return result
    if core and leaf in CORE_TABLES:
        return typed("core_products", CORE_TABLES[leaf], "src/mtp/workflow/contracts.py:78", "tsv")
    if leaf == "merged-kg_nodes.tsv" and (core or name.startswith("tests/fixtures/workflow/")):
        result = typed("mtp", "Node", "schema/mtp.yaml:175", "tsv")
        result["also_schema"] = "bsdbng"
        return result
    if name == "tests/fixtures/trait_key_policy/historical_implementations.json":
        return typed(
            "workflow",
            "FileInventoryEntry",
            "tests/test_historical_trait_policy_core.py:95",
            "code_inventories",
        )
    if core and leaf == "generation.json":
        return typed(
            "workflow",
            "FileInventoryEntry",
            "tests/fixtures/trait_key_policy/README.md",
            "producer_inventory",
        )
    if core_path == "core/steps/ingest/derived/ingest_log.json":
        return {
            "kind": "noninstance",
            "format": "JSON event list",
            "reader": "src/mtp/workflow/operations.py:494",
            "reason": "Ingest diagnostics are event dictionaries, not a Study or run contract",
        }
    if name in CUSTOM_JSON:
        return {
            "kind": "noninstance",
            "format": "producer-specific JSON",
            "reader": CUSTOM_JSON[name],
            "reason": (
                "Custom summary, protocol or provenance object; no corresponding generated class"
            ),
        }
    if name in {
        "tests/fixtures/ora_figure_1505/combined.json",
        "tests/fixtures/ora_figure_1505/separate.json",
    }:
        return {
            "kind": "noninstance",
            "format": "Plotly figure JSON",
            "reader": "tests/test_ora_figure_single_read.py:17",
            "reason": "Plotly data/layout, not a workflow or core product record",
        }
    if name in VISUAL_JSON:
        return {
            "kind": "noninstance",
            "format": "committed visualization JSON",
            "reader": "trait-visualization/figure_methods.md",
            "reason": (
                "Custom figure manifest, graph or provenance; "
                "five owned schemas do not declare this format"
            ),
        }
    if name == "methods/multiple_testing_families.yaml":
        return {
            "kind": "noninstance",
            "format": "registry YAML",
            "reader": "tests/test_multiple_testing_registry.py",
            "reason": "Executable-site registry, not the per-analysis HypothesisFamily record",
        }
    if name in FIXTURE_DOCS:
        return {
            "kind": "noninstance",
            "format": "Markdown documentation",
            "reader": name,
            "reason": (
                "Historical narrative or fixture instructions, not a serialized domain record"
            ),
        }
    if leaf.endswith((".tsv", ".tsv.gz", ".csv")):
        if core_path == "inputs/full_dump.csv" or name in {
            "tests/fixtures/workflow/full_dump.csv",
            "tests/fixtures/workflow/contrast_population.csv",
        }:
            return {
                "kind": "noninstance",
                "format": "raw BugSigDB CSV",
                "reader": "tests/fixtures/workflow/README.md",
                "reason": (
                    "Raw export includes comment/header rows; "
                    "ingestion constructs nested Study records"
                ),
            }
        if core_path == "inputs/merged-kg_edges.tsv" or name in {
            "tests/fixtures/workflow/merged-kg_edges.tsv",
            "tests/fixtures/coverage/taxonomy_rank_class.tsv",
        }:
            return {
                "kind": "noninstance",
                "format": "raw KGX edge TSV",
                "reader": "src/mtp/taxa_trait_matrix.py:292",
                "reason": (
                    "Raw KGX id/primary_knowledge_source fields differ from "
                    "Edge edge_id/provided_by; no invented remapping"
                ),
            }
        if core_path in PHYSICAL_TABLES or name == "data/derived/study_trait_profiles_wide.tsv.gz":
            return {
                "kind": "noninstance",
                "format": "physical wide/metadata TSV",
                "reader": "src/mtp/workflow/contracts.py",
                "reason": (
                    "Dynamic wide counts or physical feature metadata; not direct logical "
                    "CountProfileRow, ORAProfileRow or TaxonTraitMatrixRow. "
                    "No missing units, feature hashes or nested counts fabricated"
                ),
            }
        if name in NATIVE_TABLES:
            return {
                "kind": "noninstance",
                "format": "native result TSV",
                "reader": "src/mtp/native_influence_publication.py:227",
                "reason": (
                    "Native result columns are a custom statistical output, not a core-product row"
                ),
            }
    if name in CONFIG_FILES:
        return {
            "kind": "noninstance",
            "format": "tool configuration",
            "reader": ".pre-commit-config.yaml",
            "reason": (
                "Vendor configuration; syntax checked without claiming domain-schema validation"
            ),
        }
    return result


def _git(*args: str) -> bytes:
    return str(Repo(ROOT).git.execute(["git", *args])).encode()


def tracked_candidates() -> tuple[list[str], list[str]]:
    """Enumerate HEAD and index paths so staged deletions cannot omit committed inputs."""
    names = {
        value.decode("utf-8")
        for command in (("ls-files", "-z"), ("ls-tree", "-r", "--name-only", "-z", "HEAD"))
        for value in _git(*command).split(b"\0")
        if value
    }
    candidates = [
        name
        for name in names
        if name.lower().endswith((*TEXT_ENDINGS, ".gz", ".gzip", ".tar", ".tgz", ".zip", ".zst"))
        or name in SPECIAL_JSON
        or name.startswith(FIXTURE_ROOTS)
    ]
    other_archives = [
        name for name in names if name.lower().endswith((".zip", ".tar", ".tgz", ".zst"))
    ]
    return sorted(set(candidates)), sorted(set(other_archives))


def archive_members(path: Path) -> list[dict[str, Any]]:
    """Read bounded regular tar members without extracting or following links."""
    rows: list[dict[str, Any]] = []
    names: set[str] = set()
    total = 0
    with tarfile.open(path, "r|*") as archive:
        for member in archive:
            if len(rows) >= MAX_ARCHIVE_MEMBERS:
                raise ValueError("Archive member-count limit exceeded")
            name = member.name
            if name in names or Path(name).is_absolute() or ".." in Path(name).parts:
                raise ValueError(f"Unsafe or duplicate archive member: {name!r}")
            names.add(name)
            if not member.isfile():
                raise ValueError(f"Archive member is not a regular file: {name!r}")
            total += member.size
            if member.size > MAX_MEMBER_BYTES or total > MAX_ARCHIVE_BYTES:
                raise ValueError("Archive decompressed-size limit exceeded")
            handle = archive.extractfile(member)
            if handle is None:
                raise ValueError(f"Archive member unavailable: {name!r}")
            with handle:
                header = handle.read(HEADER_BYTES)
            rows.append(
                {
                    "member": name,
                    "size_bytes": member.size,
                    "header_preview": header.decode("utf-8"),
                    "header_bytes_read": len(header),
                    "mapping": classification(f"{path.relative_to(ROOT)}!{name}"),
                }
            )
    return rows


def inventory() -> dict[str, Any]:
    """Collect bounded header previews; preserve unknowns rather than infer validity."""
    candidates, other_archives = tracked_candidates()
    rows = []
    for name in candidates:
        path = ROOT / name
        stat = path.lstat()
        compressed = name.lower().endswith((".gz", ".gzip"))
        row: dict[str, Any] = {
            "path": name,
            "size_bytes": stat.st_size,
            "compressed": compressed,
            "suffixes": path.suffixes,
            "mapping": classification(name),
        }
        if path.is_symlink() or not path.is_file():
            row["header_error"] = "tracked candidate is not a regular non-symlink file"
        elif name.endswith((".tar.gz", ".tar", ".tgz")):
            row["archive_members"] = archive_members(path)
        else:
            try:
                opener = gzip.open if compressed else Path.open
                with opener(path, "rb") as handle:
                    raw = handle.read(HEADER_BYTES)
                row["header_preview"] = raw.decode("utf-8")
                row["header_bytes_read"] = len(raw)
            except (OSError, UnicodeError, EOFError) as error:
                row["header_error"] = f"{type(error).__name__}: {error}"
        rows.append(row)
    return {
        "created_at": datetime.now(UTC).isoformat(),
        "head": _git("rev-parse", "HEAD").decode().strip(),
        "mode": "bounded_tracked_inventory_only",
        "header_limit_bytes": HEADER_BYTES,
        "candidate_files": len(rows),
        "fixture_example_workflow_files": sum(
            row["path"].startswith(FIXTURE_ROOTS) for row in rows
        ),
        "archive_member_files": sum(len(row.get("archive_members", [])) for row in rows),
        "mapping_counts": dict(sorted(Counter(row["mapping"]["kind"] for row in rows).items())),
        "physical_bytes": sum(row["size_bytes"] for row in rows),
        "top_directory_counts": dict(
            sorted(Counter(Path(row["path"]).parts[0] for row in rows).items())
        ),
        "other_tracked_archives_not_expanded": other_archives,
        "files": rows,
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_bytes(value: Any) -> bytes:
    """Serialize without omitting model defaults, nulls or nested fields."""
    return json.dumps(
        value, sort_keys=True, ensure_ascii=False, separators=(",", ":"), default=str
    ).encode()


class JsonError(Protocol):
    absolute_path: Any
    absolute_schema_path: Any
    validator: Any
    validator_value: Any
    message: str
    instance: Any
    context: list[Any]


class JsonValidator(Protocol):
    def iter_errors(self, instance: Any) -> Iterator[JsonError]: ...


class Context(Protocol):
    def json_schema_validator(
        self, *, closed: bool, include_range_class_descendants: bool
    ) -> JsonValidator: ...


class LinkResult(Protocol):
    message: str


class LinkValidator(Protocol):
    def iter_results(self, instance: Any, target_class: str) -> Iterator[LinkResult]: ...


class ValidatorFactory(Protocol):
    def __call__(self, schema: str, *, validation_plugins: list[object]) -> LinkValidator: ...


class PluginFactory(Protocol):
    def __call__(self, *, closed: bool) -> object: ...


class SchemaLinter(Protocol):
    @staticmethod
    def validate_schema(schema_path: str) -> Iterator[LinkResult]: ...


class AllConstraints:
    """Observe the actual LinkML context before best_match hides sibling errors.

    The stock closed JsonschemaValidationPlugin still runs. This observer retains
    every generated-JSON-Schema error, including nested anyOf branches, so a new
    pattern failure cannot disappear behind an inherited unrelated failure.
    """

    def __init__(self) -> None:
        self.errors: list[dict[str, Any]] = []

    def pre_process(self, _context: Context) -> None:
        self.errors = []

    def post_process(self, _context: Context) -> None:
        return None

    def process(self, instance: Any, context: Context) -> Iterator[LinkResult]:
        validator = context.json_schema_validator(closed=True, include_range_class_descendants=True)
        for error in validator.iter_errors(instance):
            self._collect(error)
        return iter(())

    def _collect(self, error: JsonError) -> None:
        if len(self.errors) >= 10000:
            raise ValueError("Per-record error bound exceeded; no errors silently truncated")
        self.errors.append(
            {
                "path": list(error.absolute_path),
                "schema_path": list(error.absolute_schema_path),
                "validator": error.validator,
                "constraint": error.validator_value,
                "message": error.message,
                "raw_value": error.instance,
            }
        )
        for child in error.context:
            self._collect(child)


class ValidationSide:
    """Load each side's schemas and self-contained generated modules separately."""

    def __init__(self, root: Path, label: str) -> None:
        self.root = root
        self.label = label
        self.models: dict[str, Any] = {}
        self.validators: dict[str, tuple[LinkValidator, AllConstraints]] = {}
        self.schema_documents: dict[str, Any] = {}
        self.validator_factory = cast(
            ValidatorFactory, importlib.import_module("linkml.validator").Validator
        )
        self.plugin_factory = cast(
            PluginFactory,
            importlib.import_module("linkml.validator.plugins").JsonschemaValidationPlugin,
        )
        for schema in SCHEMAS:
            path = root / MODEL_PATHS[schema]
            module_name = f"_issue1316_{label}_{schema}"
            spec = importlib.util.spec_from_file_location(module_name, path)
            if spec is None or spec.loader is None:
                raise ValueError(f"Cannot load generated module: {path}")
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            spec.loader.exec_module(module)
            self.models[schema] = module
            schema_path = root / "schema" / f"{schema}.yaml"
            self.schema_documents[schema] = yaml.safe_load(schema_path.read_text())

    def typed_table_row(self, schema: str, target: str, row: dict[str, str]) -> dict[str, Any]:
        """Parse only schema-declared integer TSV cells; retain empty strings.

        Core contracts validate original strings with Pydantic. LinkML receives
        the corresponding primitive integers, independently of either generated
        model. All string cells, especially patterned identifiers, stay exact.
        This is not the LinkML CLI TSV loader, which drops empty cells.
        """
        attributes = self.schema_documents[schema]["classes"][target].get("attributes", {})
        result: dict[str, Any] = dict(row)
        for key, definition in attributes.items():
            if definition.get("range") == "integer" and key in row:
                try:
                    result[key] = int(row[key])
                except ValueError:
                    result[key] = row[key]
        return result

    def validate(
        self, schema: str, target: str, value: Any, *, table: bool = False
    ) -> tuple[dict[str, Any], bytes | None]:
        if schema not in self.validators:
            observer = AllConstraints()
            validator = self.validator_factory(
                str(self.root / "schema" / f"{schema}.yaml"),
                validation_plugins=[self.plugin_factory(closed=True), observer],
            )
            self.validators[schema] = validator, observer
        validator, observer = self.validators[schema]
        linked_value = self.typed_table_row(schema, target, value) if table else value
        cli_messages = [result.message for result in validator.iter_results(linked_value, target)]
        result: dict[str, Any] = {
            "linkml": {
                "valid": not observer.errors,
                "errors": observer.errors,
                "stock_plugin_messages": cli_messages,
            },
        }
        model = cast(type[BaseModel], getattr(self.models[schema], target))
        serialized: bytes | None = None
        try:
            instance = model.model_validate(value)
            serialized = stable_bytes(
                instance.model_dump(
                    mode="json", exclude_none=False, exclude_defaults=False, exclude_unset=False
                )
            )
            result["model"] = {
                "valid": True,
                "errors": [],
                "serialized_sha256": hashlib.sha256(serialized).hexdigest(),
                "serialized_size_bytes": len(serialized),
            }
        except ValidationError as error:
            errors = []
            for detail in error.errors(include_url=False):
                errors.append(
                    {
                        "path": list(detail["loc"]),
                        "validator": detail["type"],
                        "message": detail["msg"],
                        "raw_value": detail.get("input"),
                        "context": detail.get("ctx"),
                    }
                )
            result["model"] = {"valid": False, "errors": errors}
        return result, serialized


def error_identity(error: dict[str, Any]) -> bytes:
    return stable_bytes([error["path"], error["validator"], error.get("schema_path")])


def constraint_delta(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    old = {error_identity(error): error for error in before["errors"]}
    new = {error_identity(error): error for error in after["errors"]}
    return {
        "newly_failing_constraints": [new[key] for key in sorted(new.keys() - old.keys())],
        "no_longer_failing_constraints": [old[key] for key in sorted(old.keys() - new.keys())],
        "inherited_failure_constraints": len(old.keys() & new.keys()),
        "valid_to_invalid": before["valid"] and not after["valid"],
    }


def document_records(mapping: dict[str, Any], data: Any) -> Iterator[tuple[str, Any]]:
    encoding = mapping.get("encoding")
    if encoding == "array":
        if not isinstance(data, list):
            raise ValueError("Expected top-level record array")
        yield from ((f"/{index}", item) for index, item in enumerate(data))
    elif encoding == "code_inventories":
        for index, item in enumerate(data["records"]):
            yield f"/records/{index}/code_inventory", item["code_inventory"]
    elif encoding == "producer_inventory":
        yield "/producer_matrix_inventory", data["producer_matrix_inventory"]
    else:
        yield "", data


@contextmanager
def input_stream(path: Path) -> Iterator[BinaryIO]:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Input is not a regular non-symlink file: {path}")
    if path.name.endswith((".gz", ".gzip")):
        with gzip.open(path, "rb") as handle:
            yield cast(BinaryIO, handle)
    else:
        with path.open("rb") as handle:
            yield handle


def bounded_document(handle: BinaryIO, name: str) -> Any:
    raw = handle.read(MAX_MEMBER_BYTES + 1)
    if len(raw) > MAX_MEMBER_BYTES:
        raise ValueError(f"Serialized document exceeds {MAX_MEMBER_BYTES} bytes: {name}")
    text = raw.decode("utf-8")
    if name.endswith((".yaml", ".yml")):
        return yaml.safe_load(text)
    return json.loads(text)


class BoundedCSVLines:
    """Bound the whole logical record, including quoted multiline cells."""

    def __init__(self, source: TextIO) -> None:
        self.source = source
        self.characters = 0

    def __iter__(self) -> BoundedCSVLines:
        return self

    def __next__(self) -> str:
        line = self.source.readline(MAX_RECORD_CHARS + 1)
        if not line:
            raise StopIteration
        self.characters += len(line)
        if self.characters > MAX_RECORD_CHARS:
            raise ValueError(f"CSV record exceeds {MAX_RECORD_CHARS} characters")
        return line


def bounded_csv_rows(text: TextIO, *, delimiter: str) -> Iterator[list[str]]:
    lines = BoundedCSVLines(text)
    for row in csv.reader(lines, delimiter=delimiter):
        if len(row) > MAX_COLUMNS:
            raise ValueError(f"CSV record exceeds {MAX_COLUMNS} columns")
        lines.characters = 0
        yield row


def table_records(handle: BinaryIO, name: str) -> Iterator[tuple[int, dict[str, str]]]:
    """Stream rectangular direct row records; no strip, field drops or filling."""
    with io.TextIOWrapper(handle, encoding="utf-8", newline="") as text:
        reader = bounded_csv_rows(text, delimiter="\t")
        header = next(reader)
        if len(set(header)) != len(header):
            raise ValueError(f"Duplicate table header: {name}")
        for number, values in enumerate(reader, 2):
            if len(values) != len(header):
                raise ValueError(f"Ragged direct-instance table {name}, row {number}")
            yield number, dict(zip(header, values, strict=True))


class Comparison:
    def __init__(self, baseline: Path, sink: TextIO) -> None:
        self.before = ValidationSide(baseline, "baseline")
        self.after = ValidationSide(ROOT, "current")
        self.sink = sink
        self.counts: Counter[str] = Counter()
        self.files: list[dict[str, Any]] = []

    def record(
        self, name: str, locator: str, mapping: dict[str, Any], value: Any, schema: str
    ) -> None:
        target = mapping["class"]
        table = mapping.get("encoding") == "tsv"
        before, old_bytes = self.before.validate(schema, target, value, table=table)
        after, new_bytes = self.after.validate(schema, target, value, table=table)
        row: dict[str, Any] = {
            "path": name,
            "record": locator,
            "schema": schema,
            "class": target,
            "raw_record_sha256": hashlib.sha256(stable_bytes(value)).hexdigest(),
            "baseline": before,
            "current": after,
            "delta": {
                layer: constraint_delta(before[layer], after[layer])
                for layer in ("linkml", "model")
            },
        }
        if old_bytes is not None and new_bytes is not None:
            equal = old_bytes == new_bytes
            row["model_serialized_bytes_equal"] = equal
            self.counts["both_models_valid"] += 1
            self.counts[
                "model_serialization_equal" if equal else "model_serialization_changed"
            ] += 1
        else:
            row["model_serialized_bytes_equal"] = None
        for layer in ("linkml", "model"):
            for side in ("baseline", "current"):
                self.counts[
                    f"{side}_{layer}_{'valid' if row[side][layer]['valid'] else 'invalid'}"
                ] += 1
            self.counts[f"{layer}_newly_failing_constraints"] += len(
                row["delta"][layer]["newly_failing_constraints"]
            )
            self.counts[f"{layer}_valid_to_invalid"] += int(row["delta"][layer]["valid_to_invalid"])
        self.counts["validated_record_schema_pairs"] += 1
        self.sink.write(stable_bytes(row).decode() + "\n")

    def consume(self, name: str, handle: BinaryIO, size: int) -> None:
        mapping = classification(name)
        result: dict[str, Any] = {
            "path": name,
            "size_bytes": size,
            "mapping": mapping,
            "records": 0,
        }
        kind = mapping["kind"]
        if kind == "unknown":
            self.counts["unclassified_files"] += 1
            self.files.append(result)
            return
        if kind == "instance":
            schema_names = [mapping["schema"]]
            if "also_schema" in mapping:
                schema_names.append(mapping["also_schema"])
            if mapping.get("encoding") == "tsv":
                for number, value in table_records(handle, name):
                    for schema in schema_names:
                        self.record(name, f"/row/{number}", mapping, value, schema)
                    result["records"] += 1
            else:
                data = bounded_document(handle, name)
                for locator, value in document_records(mapping, data):
                    for schema in schema_names:
                        self.record(name, locator, mapping, value, schema)
                    result["records"] += 1
            result["schema_validations_per_record"] = len(schema_names)
            if mapping.get("encoding") in {"producer_inventory", "code_inventories"}:
                result["projection_limitation"] = (
                    "Only the named inventory subrecord has a generated contract; "
                    "enclosing helper receipt is syntax checked, "
                    "not claimed to be a WorkflowRun"
                )
        elif kind == "schema":
            linter = cast(SchemaLinter, importlib.import_module("linkml.linter.linter").Linter)
            result["metamodel"] = {
                label: [problem.message for problem in linter.validate_schema(str(root / name))]
                for label, root in (("baseline", self.before.root), ("current", ROOT))
            }
            result["model_validation"] = (
                "Not applicable: schema definitions are validated against the LinkML "
                "metamodel, not the five domain models"
            )
            self.counts["schema_documents"] += 1
            self.counts["current_schema_errors"] += len(result["metamodel"]["current"])
        elif kind == "noninstance":
            if name.endswith((*TEXT_ENDINGS, ".baseline")):
                data = bounded_document(handle, name)
                result["syntax"] = "parsed"
                result["top_level_type"] = type(data).__name__
                result["records"] = len(data) if isinstance(data, list) else 1
                if isinstance(data, dict):
                    result["top_level_keys"] = sorted(str(key) for key in data)
            elif name.endswith((".tsv", ".tsv.gz", ".csv")):
                # Raw exports may have comment rows or deliberate malformed KGX
                # headers; count physical records without claiming domain validity.
                with io.TextIOWrapper(handle, encoding="utf-8", newline="") as text:
                    reader = bounded_csv_rows(
                        text, delimiter="," if name.endswith(".csv") else "\t"
                    )
                    widths: Counter[int] = Counter()
                    for number, row in enumerate(reader):
                        widths[len(row)] += 1
                        if number == 0:
                            result["physical_first_row"] = row
                    result["physical_csv_records_including_header_comments"] = sum(widths.values())
                    result["row_width_counts"] = dict(sorted(widths.items()))
                result["records"] = None
            else:
                raw = handle.read(MAX_MEMBER_BYTES + 1)
                if len(raw) > MAX_MEMBER_BYTES:
                    raise ValueError(f"Document-size limit exceeded: {name}")
                raw.decode("utf-8")
                result["syntax"] = "UTF-8 text, no serialized record schema"
            self.counts["noninstance_files"] += 1
        else:
            raise ValueError(f"Unexpected standalone mapping {kind}: {name}")
        self.files.append(result)


def source_inventory(root: Path) -> dict[str, str]:
    paths = [f"schema/{name}.yaml" for name in SCHEMAS] + list(MODEL_PATHS.values())
    return {name: sha256_file(root / name) for name in paths}


def sha256_stream(handle: BinaryIO) -> str:
    digest = hashlib.sha256()
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(chunk)
    return digest.hexdigest()


def committed_input_bindings(names: list[str]) -> dict[str, dict[str, Any]]:
    """Require HEAD, index and working records to agree in bytes and Git mode.

    Only the five deliberately edited schema source files permit byte changes.
    They still require ordinary regular-file modes, and all three hashes are
    retained. Generated models are bound separately as changed source products.
    """
    repo = Repo(ROOT)
    bindings: dict[str, dict[str, Any]] = {}
    allowed_schemas = {f"schema/{name}.yaml" for name in SCHEMAS}
    for name in names:
        path = ROOT / name
        mode = path.lstat().st_mode
        if not stat.S_ISREG(mode):
            raise ValueError(f"Committed input must be a regular file: {name}")
        try:
            head = repo.head.commit.tree / name
        except KeyError as error:
            raise ValueError(f"Input is absent from committed HEAD: {name}") from error
        entry = repo.index.entries.get((name, 0))
        if entry is None:
            raise ValueError(f"Input lacks an unconflicted index entry: {name}")
        working_mode = 0o100755 if mode & stat.S_IXUSR else 0o100644
        if (
            head.mode not in {0o100644, 0o100755}
            or entry.mode != head.mode
            or working_mode != head.mode
        ):
            raise ValueError(f"Committed/index/working Git file mode differs: {name}")
        head_hash = sha256_stream(cast(BinaryIO, head.data_stream))
        index_hash = sha256_stream(cast(BinaryIO, repo.odb.stream(entry.binsha)))
        working_hash = sha256_file(path)
        if name not in allowed_schemas and len({head_hash, index_hash, working_hash}) != 1:
            raise ValueError(f"Committed HEAD, index and working bytes differ: {name}")
        bindings[name] = {
            "head_sha256": head_hash,
            "index_sha256": index_hash,
            "working_sha256": working_hash,
            "git_mode": oct(head.mode),
            "intentional_schema_source_edit": name in allowed_schemas
            and len({head_hash, index_hash, working_hash}) != 1,
        }
    return bindings


def compare(baseline: Path, output: Path) -> dict[str, Any]:
    """Compare immutable records; never run builders, normalize or rewrite data."""
    candidates, other_archives = tracked_candidates()
    before_sources = source_inventory(baseline)
    current_sources = source_inventory(ROOT)
    helper_before = sha256_file(Path(__file__))
    baseline_commit = (baseline / "COMMIT").read_text().strip()
    baseline_tree = Repo(ROOT).commit(baseline_commit).tree
    for name, expected_hash in before_sources.items():
        if sha256_stream(cast(BinaryIO, (baseline_tree / name).data_stream)) != expected_hash:
            raise ValueError(f"Baseline source differs from its declared committed blob: {name}")
    committed_bindings = committed_input_bindings(candidates)
    before_inputs = {name: row["working_sha256"] for name, row in committed_bindings.items()}
    records_path = output.with_suffix(".records.jsonl")
    with records_path.open("x", encoding="utf-8") as sink:
        comparison = Comparison(baseline, sink)
        for name in candidates:
            path = ROOT / name
            if classification(name)["kind"] == "archive":
                members = archive_members(path)
                expected = {row["member"]: row for row in members}
                seen: set[str] = set()
                with tarfile.open(path, "r|*") as archive:
                    for member in archive:
                        if (
                            member.name not in expected
                            or member.name in seen
                            or not member.isfile()
                        ):
                            raise ValueError("Archive changed between inventory and validation")
                        seen.add(member.name)
                        extracted = archive.extractfile(member)
                        if extracted is None:
                            raise ValueError("Missing archive member stream")
                        with extracted:
                            raw = extracted.read(MAX_MEMBER_BYTES + 1)
                        if (
                            len(raw) != expected[member.name]["size_bytes"]
                            or len(raw) > MAX_MEMBER_BYTES
                        ):
                            raise ValueError("Archive member size changed or exceeded limit")
                        comparison.consume(f"{name}!{member.name}", io.BytesIO(raw), len(raw))
                        comparison.files[-1]["member_sha256"] = hashlib.sha256(raw).hexdigest()
                if seen != expected.keys():
                    raise ValueError("Archive member set changed")
                comparison.files.append(
                    {"path": name, "mapping": classification(name), "members": len(seen)}
                )
            else:
                with input_stream(path) as handle:
                    comparison.consume(name, handle, path.stat().st_size)
    after_bindings = committed_input_bindings(candidates)
    if committed_bindings != after_bindings or helper_before != sha256_file(Path(__file__)):
        raise ValueError("Committed/index/working input binding or helper changed during audit")
    if baseline_commit != (baseline / "COMMIT").read_text().strip():
        raise ValueError("Baseline commit declaration changed during audit")
    after_inputs = {name: row["working_sha256"] for name, row in after_bindings.items()}
    if (
        before_inputs != after_inputs
        or before_sources != source_inventory(baseline)
        or current_sources != source_inventory(ROOT)
    ):
        raise ValueError("Inputs or schema/model sources changed during audit")
    count = comparison.counts
    refused = any(
        count[key]
        for key in (
            "unclassified_files",
            "model_serialization_changed",
            "linkml_newly_failing_constraints",
            "model_newly_failing_constraints",
            "current_schema_errors",
        )
    )
    return {
        "created_at": datetime.now(UTC).isoformat(),
        "head": _git("rev-parse", "HEAD").decode().strip(),
        "baseline_commit": baseline_commit,
        "helper_sha256": helper_before,
        "status": "REVIEW_REQUIRED" if refused else "NO_NEW_FAILURES_OR_SERIALIZATION_CHANGES",
        "tracked_files": len(candidates),
        "counts": dict(sorted(count.items())),
        "other_archive_candidates": other_archives,
        "files": comparison.files,
        "baseline_sources": before_sources,
        "current_sources": current_sources,
        "input_sha256": before_inputs,
        "committed_input_bindings": committed_bindings,
        "input_and_source_hashes_unchanged": True,
        "records": {"path": str(records_path), "sha256": sha256_file(records_path)},
        "scope": (
            "All HEAD/index YAML/JSON and compressed formats plus every HEAD/index "
            "fixture/example/workflow file, including archive members. Non-schema "
            "inputs bind HEAD, index and working bytes; the five intentionally edited "
            "schema sources are separately identified. No execution records outside "
            "committed HEAD and the git index."
        ),
        "limitations": [
            "Inherited validation failures remain failures and are reported without suppression.",
            (
                "Core TSV string cells and empty values are preserved; only declared "
                "integer cells are typed independently for LinkML. "
                "Pydantic receives original strings."
            ),
            (
                "Wide physical tables and raw exports are inventoried, "
                "not fabricated as nested domain records."
            ),
            (
                "Schema documents use the actual LinkML metamodel validator; "
                "recipes and lint are separate evidence."
            ),
            (
                "Exact model serialization comparisons include all defaults, "
                "null values and nested objects."
            ),
            "No data regenerated; no scientific comparison or historical run verification claimed.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("inventory", "compare"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline-root", type=Path)
    args = parser.parse_args()
    if not args.output.parent.is_dir():
        parser.error("Output parent must already exist")
    if args.output.exists() or args.output.is_symlink():
        parser.error("Refusing to overwrite output")
    if args.mode == "compare" and args.baseline_root is None:
        parser.error("compare requires --baseline-root")
    record = inventory() if args.mode == "inventory" else compare(args.baseline_root, args.output)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    print(
        json.dumps(
            {
                key: record[key]
                for key in (
                    "mode",
                    "status",
                    "head",
                    "candidate_files",
                    "fixture_example_workflow_files",
                    "archive_member_files",
                    "tracked_files",
                    "counts",
                )
                if key in record
            },
            indent=2,
        )
    )
    if record.get("status") == "REVIEW_REQUIRED":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
