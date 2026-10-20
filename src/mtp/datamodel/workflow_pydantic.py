from __future__ import annotations

import re
import sys
from datetime import (
    date,
    datetime,
    time
)
from decimal import Decimal
from enum import Enum
from typing import (
    Any,
    ClassVar,
    Literal,
    Optional,
    Union
)

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    RootModel,
    SerializationInfo,
    SerializerFunctionWrapHandler,
    field_validator,
    model_serializer
)


metamodel_version = "1.7.0"
version = "1.0.0"


class ConfiguredBaseModel(BaseModel):
    model_config = ConfigDict(
        serialize_by_alias = True,
        validate_by_name = True,
        validate_assignment = True,
        validate_default = True,
        extra = "forbid",
        arbitrary_types_allowed = True,
        use_enum_values = True,
        strict = False,
    )





class LinkMLMeta(RootModel):
    root: dict[str, Any] = {}
    model_config = ConfigDict(frozen=True)

    def __getattr__(self, key:str):
        return getattr(self.root, key)

    def __getitem__(self, key:str):
        return self.root[key]

    def __setitem__(self, key:str, value):
        self.root[key] = value

    def __contains__(self, key:str) -> bool:
        return key in self.root


linkml_meta = LinkMLMeta({'default_prefix': 'mtpw',
     'default_range': 'string',
     'description': 'Versioned product and workflow contracts, execution attempts, '
                    'and verified artifact instances. Structural validation is '
                    'supplemented by graph, filesystem, checksum, and scientific '
                    'product checks.',
     'id': 'https://w3id.org/mtp/workflow',
     'imports': ['linkml:types'],
     'license': 'MIT',
     'name': 'mtp_workflow',
     'prefixes': {'linkml': {'prefix_prefix': 'linkml',
                             'prefix_reference': 'https://w3id.org/linkml/'},
                  'mtpw': {'prefix_prefix': 'mtpw',
                           'prefix_reference': 'https://w3id.org/mtp/workflow/'},
                  'prov': {'prefix_prefix': 'prov',
                           'prefix_reference': 'http://www.w3.org/ns/prov#'}},
     'source_file': 'schema/workflow.yaml',
     'title': 'Execution-checked microbial trait workflow'} )

class StorageKind(str, Enum):
    """
    Permitted storage kind values.
    """
    file = "file"
    collection = "collection"


class ExecutionStatus(str, Enum):
    """
    Permitted execution status values.
    """
    running = "running"
    completed = "completed"
    failed = "failed"
    imported_unverified = "imported_unverified"


class ValidationStatus(str, Enum):
    """
    Permitted validation status values.
    """
    passed = "passed"
    failed = "failed"



class ProductDefinition(ConfiguredBaseModel):
    """
    A stable scientific product contract, independent of its location.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Product Definition'})

    id: str = Field(default=..., title="Id", description="""Identifier unique within its containing document.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive']} })
    description: Optional[str] = Field(default=None, title="Description", description="""Human-readable product meaning.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'StepDefinition', 'WorkflowDefinition']} })
    kind: str = Field(default=..., title="Kind", description="""Scientific validator contract name.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition']} })
    storage_kind: StorageKind = Field(default=..., title="Storage kind", description="""Whether the artifact is a regular file or an inventoried directory.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArchivedLocation']} })
    format: str = Field(default=..., title="Format", description="""Serialization format.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition']} })
    filename: str = Field(default=..., title="Filename", description="""Relative path inside a producer bundle.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition']} })
    external: Optional[bool] = Field(default=False, title="External", description="""Whether this product is supplied explicitly rather than produced by a step.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition'], 'ifabsent': 'false'} })
    schema_ref: Optional[str] = Field(default=None, title="Schema ref", description="""Versioned scientific schema reference.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition']} })
    unit: Optional[str] = Field(default=None, title="Unit", description="""Meaning of one row or observation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArtifactInstance']} })
    cohort: Optional[str] = Field(default=None, title="Cohort", description="""Expected population selection policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArtifactInstance', 'ConsumerRecord']} })
    cohort_parameter: Optional[str] = Field(default=None, title="Cohort parameter", description="""Parameter determining effective cohort.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition']} })
    feature_policy: Optional[str] = Field(default=None, title="Feature policy", description="""Feature inclusion and identity policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArtifactInstance', 'ConsumerRecord']} })
    coverage_policy: Optional[str] = Field(default=None, title="Coverage policy", description="""Required coverage and exclusion accounting.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArtifactInstance']} })

    @field_validator('id')
    def pattern_id(cls, v):
        pattern=re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid id format: {v}"
            raise ValueError(err_msg)
        return v


class ParameterDefinition(ConfiguredBaseModel):
    """
    An explicitly resolved string parameter; undeclared overrides are rejected.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'Parameter Definition'})

    name: str = Field(default=..., title="Name", description="""Declared parameter identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ParameterDefinition', 'ParameterValue']} })
    default: Optional[str] = Field(default=None, title="Default", description="""Default value before explicit environment and API overrides.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ParameterDefinition']} })
    environment: Optional[str] = Field(default=None, title="Environment", description="""Optional explicitly permitted environment variable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ParameterDefinition']} })
    required: Optional[bool] = Field(default=False, title="Required", description="""Whether execution requires an effective parameter value.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ParameterDefinition'], 'ifabsent': 'false'} })

    @field_validator('name')
    def pattern_name(cls, v):
        pattern=re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid name format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid name format: {v}"
            raise ValueError(err_msg)
        return v


class ParameterValue(ConfiguredBaseModel):
    """
    The effective value and source of one declared parameter.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Parameter Value'})

    name: str = Field(default=..., title="Name", description="""Declared parameter identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ParameterDefinition', 'ParameterValue']} })
    value: str = Field(default=..., title="Value", description="""Effective resolved parameter value passed to the operation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ParameterValue']} })
    source: str = Field(default=..., title="Source", description="""Origin of the effective value: default, declared environment variable, or override.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ParameterValue']} })


class StepDefinition(ConfiguredBaseModel):
    """
    An operation registry key and its typed input and output products.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Step Definition'})

    id: str = Field(default=..., title="Id", description="""Identifier unique within its containing document.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive']} })
    operation: str = Field(default=..., title="Operation", description="""Key resolving to an explicitly registered Python callable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition']} })
    description: Optional[str] = Field(default=None, title="Description", description="""Purpose of this step.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'StepDefinition', 'WorkflowDefinition']} })
    inputs: list[str] = Field(default=..., title="Inputs", description="""Product references whose verified artifact paths are passed to the operation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition']} })
    outputs: list[str] = Field(default=..., title="Outputs", description="""Declared outputs; references are resolved and checked against their producers.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition', 'WorkflowRun', 'ConsumerRecord']} })
    depends_on: Optional[list[str]] = Field(default=None, title="Depends on", description="""Additional prerequisite step references beyond input-producing dependencies.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition']} })
    parameters: Optional[list[str]] = Field(default=None, title="Parameters", description="""Declared or effective parameters scoped to this definition or execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowDefinition',
                       'StepExecution',
                       'WorkflowRun']} })

    @field_validator('id')
    def pattern_id(cls, v):
        pattern=re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid id format: {v}"
            raise ValueError(err_msg)
        return v


class WorkflowDefinition(ConfiguredBaseModel):
    """
    Versioned graph of typed products and registered Python operations.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'Workflow Definition',
         'tree_root': True})

    id: str = Field(default=..., title="Id", description="""Identifier unique within its containing document.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive']} })
    version: str = Field(default=..., title="Version", description="""Version of the declared workflow contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowDefinition']} })
    description: Optional[str] = Field(default=None, title="Description", description="""Workflow purpose and scope.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'StepDefinition', 'WorkflowDefinition']} })
    products: list[ProductDefinition] = Field(default=..., title="Products", description="""Scientific product contracts available to this workflow.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowDefinition']} })
    steps: list[StepDefinition] = Field(default=..., title="Steps", description="""Ordered or dependency-resolved workflow stages or their attempted executions.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowDefinition', 'WorkflowRun']} })
    parameters: Optional[list[ParameterDefinition]] = Field(default=None, title="Parameters", description="""Declared or effective parameters scoped to this definition or execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowDefinition',
                       'StepExecution',
                       'WorkflowRun']} })

    @field_validator('id')
    def pattern_id(cls, v):
        pattern=re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid id format: {v}"
            raise ValueError(err_msg)
        return v


class FileInventoryEntry(ConfiguredBaseModel):
    """
    One regular file in an artifact; symlinks and special files are rejected.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'File Inventory Entry'})

    path: str = Field(default=..., title="Path", description="""Relative member path, or the file basename.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance']} })
    sha256: str = Field(default=..., title="Sha256", description="""Lowercase SHA-256 of file bytes or the canonical collection inventory.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance', 'ArchivedLocation']} })
    size_bytes: int = Field(default=..., title="Size bytes", description="""Total number of regular-file bytes represented by this record.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance', 'ArchivedLocation']} })

    @field_validator('sha256')
    def pattern_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class ValidationResult(ConfiguredBaseModel):
    """
    A recorded structural, semantic, or scientific validation outcome.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Validation Result'})

    validator: str = Field(default=..., title="Validator", description="""Versioned validation procedure identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult']} })
    status: ValidationStatus = Field(default=..., title="Status", description="""Explicit lifecycle or validation outcome; only completed verified runs are consumable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult', 'StepExecution', 'WorkflowRun']} })
    message: str = Field(default=..., title="Message", description="""Human-readable validation result.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult']} })
    details_json: Optional[str] = Field(default=None, title="Details json", description="""Validator-specific facts encoded as canonical JSON.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult']} })


class ArtifactInstance(ConfiguredBaseModel):
    """
    Immutable product content bound to a particular execution attempt.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'class_uri': 'prov:Entity',
         'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'Artifact Instance'})

    id: str = Field(default=..., title="Id", description="""Identifier unique within its containing document.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive']} })
    product_id: str = Field(default=..., title="Product id", description="""Reference to the embedded product contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance']} })
    run_id: str = Field(default=..., title="Run id", description="""Unique execution-attempt identifier, distinct from the deterministic fingerprint.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance', 'WorkflowRun', 'ConsumerRecord']} })
    producer_step: Optional[str] = Field(default=None, title="Producer step", description="""Reference to the execution that produced this artifact; absent for inputs.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance'], 'slot_uri': 'prov:wasGeneratedBy'} })
    path: str = Field(default=..., title="Path", description="""Absolute external path or run-relative output path.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance']} })
    requested_path: Optional[str] = Field(default=None, title="Requested path", description="""Original external binding before symlink resolution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance']} })
    sha256: str = Field(default=..., title="Sha256", description="""Lowercase SHA-256 of file bytes or the canonical collection inventory.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance', 'ArchivedLocation']} })
    size_bytes: int = Field(default=..., title="Size bytes", description="""Total number of regular-file bytes represented by this record.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance', 'ArchivedLocation']} })
    inventory: list[FileInventoryEntry] = Field(default=..., title="Inventory", description="""Complete sorted list of regular-file members and their byte hashes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance', 'ArchivedLocation']} })
    validations: list[ValidationResult] = Field(default=..., title="Validations", description="""Recorded inventory and scientific product validation outcomes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance']} })
    unit: Optional[str] = Field(default=None, title="Unit", description="""Effective observational unit.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArtifactInstance']} })
    cohort: Optional[str] = Field(default=None, title="Cohort", description="""Effective population selection policy or parameter value.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArtifactInstance', 'ConsumerRecord']} })
    feature_policy: Optional[str] = Field(default=None, title="Feature policy", description="""Effective feature inclusion policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArtifactInstance', 'ConsumerRecord']} })
    coverage_policy: Optional[str] = Field(default=None, title="Coverage policy", description="""Effective exclusion accounting policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArtifactInstance']} })

    @field_validator('sha256')
    def pattern_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class StepExecution(ConfiguredBaseModel):
    """
    One attempted operation with concrete inputs and staged outputs.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'class_uri': 'prov:Activity',
         'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'Step Execution'})

    id: str = Field(default=..., title="Id", description="""Identifier unique within its containing document.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive']} })
    step_id: str = Field(default=..., title="Step id", description="""Reference to the declared step contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution']} })
    status: ExecutionStatus = Field(default=..., title="Status", description="""Explicit lifecycle or validation outcome; only completed verified runs are consumable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult', 'StepExecution', 'WorkflowRun']} })
    started_at: datetime  = Field(default=..., title="Started at", description="""UTC timestamp when this execution attempt began.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution', 'WorkflowRun']} })
    finished_at: Optional[datetime ] = Field(default=None, title="Finished at", description="""UTC timestamp when this execution attempt finished or failed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution', 'WorkflowRun']} })
    input_artifacts: Optional[list[str]] = Field(default=None, title="Input artifacts", description="""Concrete artifact references actually bound as operation inputs.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution'], 'slot_uri': 'prov:used'} })
    output_artifacts: Optional[list[str]] = Field(default=None, title="Output artifacts", description="""Concrete artifact references published by this completed operation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution'], 'slot_uri': 'prov:generated'} })
    parameters: Optional[list[ParameterValue]] = Field(default=None, title="Parameters", description="""Declared or effective parameters scoped to this definition or execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowDefinition',
                       'StepExecution',
                       'WorkflowRun']} })
    error: Optional[str] = Field(default=None, title="Error", description="""Failure diagnostic; failed artifacts are never published.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution', 'WorkflowRun']} })


class WorkflowRun(ConfiguredBaseModel):
    """
    Unique execution attempt with an embedded contract and content fingerprint.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'class_uri': 'prov:Activity',
         'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'Workflow Run'})

    run_id: str = Field(default=..., title="Run id", description="""Unique execution-attempt identifier, distinct from the deterministic fingerprint.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance', 'WorkflowRun', 'ConsumerRecord']} })
    definition: WorkflowDefinition = Field(default=..., title="Definition", description="""Embedded validated workflow definition used by this attempt.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    status: ExecutionStatus = Field(default=..., title="Status", description="""Explicit lifecycle or validation outcome; only completed verified runs are consumable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult', 'StepExecution', 'WorkflowRun']} })
    started_at: datetime  = Field(default=..., title="Started at", description="""UTC timestamp when this execution attempt began.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution', 'WorkflowRun']} })
    finished_at: Optional[datetime ] = Field(default=None, title="Finished at", description="""UTC timestamp when this execution attempt finished or failed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution', 'WorkflowRun']} })
    manifest_path: str = Field(default=..., title="Manifest path", description="""Location of this run manifest; relocation does not change execution identity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    spec_path: str = Field(default=..., title="Spec path", description="""Source location of the workflow specification at execution time.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    spec_sha256: str = Field(default=..., title="Spec sha256", description="""SHA-256 of the exact source specification bytes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    definition_sha256: str = Field(default=..., title="Definition sha256", description="""SHA-256 of the canonical embedded validated workflow definition.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    execution_fingerprint: str = Field(default=..., title="Execution fingerprint", description="""Deterministic hash of contracts, inputs, code, environment, and resolved values.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    code_sha256: str = Field(default=..., title="Code sha256", description="""SHA-256 of the complete recorded code inventory.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    code_inventory: list[FileInventoryEntry] = Field(default=..., title="Code inventory", description="""Source and schema files, including untracked code, used to identify implementation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    lock_sha256: str = Field(default=..., title="Lock sha256", description="""SHA-256 of the dependency lockfile used for this attempt.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    git_commit: str = Field(default=..., title="Git commit", description="""Repository commit at execution start; code hashes also capture local changes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    git_dirty: bool = Field(default=..., title="Git dirty", description="""Whether tracked or untracked repository changes existed at execution start.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    runtime_json: str = Field(default=..., title="Runtime json", description="""Canonical JSON describing Python, installed package versions, platform, and operations.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    parameters: Optional[list[ParameterValue]] = Field(default=None, title="Parameters", description="""Declared or effective parameters scoped to this definition or execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowDefinition',
                       'StepExecution',
                       'WorkflowRun']} })
    steps: Optional[list[StepExecution]] = Field(default=None, title="Steps", description="""Ordered or dependency-resolved workflow stages or their attempted executions.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowDefinition', 'WorkflowRun']} })
    artifacts: Optional[list[ArtifactInstance]] = Field(default=None, title="Artifacts", description="""Concrete external input and published output artifact instances.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    outputs: Optional[list[str]] = Field(default=None, title="Outputs", description="""Declared outputs; references are resolved and checked against their producers.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition', 'WorkflowRun', 'ConsumerRecord']} })
    error: Optional[str] = Field(default=None, title="Error", description="""Failure diagnostic for the attempt.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution', 'WorkflowRun']} })

    @field_validator('spec_sha256')
    def pattern_spec_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid spec_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid spec_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('definition_sha256')
    def pattern_definition_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid definition_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid definition_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('execution_fingerprint')
    def pattern_execution_fingerprint(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid execution_fingerprint format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid execution_fingerprint format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('code_sha256')
    def pattern_code_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid code_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid code_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('lock_sha256')
    def pattern_lock_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid lock_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid lock_sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class ConsumerRecord(ConfiguredBaseModel):
    """
    Verified binding of downstream analysis products to a completed workflow run.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Consumer Record'})

    id: str = Field(default=..., title="Id", description="""Identifier unique within its containing document.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive']} })
    run_id: str = Field(default=..., title="Run id", description="""Unique execution-attempt identifier, distinct from the deterministic fingerprint.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance', 'WorkflowRun', 'ConsumerRecord']} })
    run_manifest_path: str = Field(default=..., title="Run manifest path", description="""Explicit path to the verified producing run manifest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord']} })
    run_manifest_sha256: str = Field(default=..., title="Run manifest sha256", description="""SHA-256 of the producing run manifest consumed by the analysis.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord']} })
    consumer: str = Field(default=..., title="Consumer", description="""Identifier of the downstream analysis or presentation operation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord']} })
    product_ids: list[str] = Field(default=..., title="Product ids", description="""Product contracts explicitly consumed by this analysis.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord']} })
    artifact_ids: list[str] = Field(default=..., title="Artifact ids", description="""Concrete verified artifacts consumed by this analysis.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord']} })
    feature_policy: str = Field(default=..., title="Feature policy", description="""Feature inclusion and identity policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArtifactInstance', 'ConsumerRecord']} })
    cohort: str = Field(default=..., title="Cohort", description="""Expected population selection policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArtifactInstance', 'ConsumerRecord']} })
    created_at: datetime  = Field(default=..., title="Created at", description="""UTC timestamp when this consumer record was created.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord', 'EvidenceArchive']} })
    coverage_path: str = Field(default=..., title="Coverage path", description="""Path to reconciled annotation coverage under the actual consumer feature filter.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord']} })
    feature_filter_id: str = Field(default=..., title="Feature filter id", description="""SHA-256 identity of the exact ordered retained analysis features.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord']} })
    coverage_scope: str = Field(default=..., title="Coverage scope", description="""Coverage describes the input cohort before analysis-specific eligibility exclusions.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord']} })
    outputs: list[FileInventoryEntry] = Field(default=..., title="Outputs", description="""Declared outputs; references are resolved and checked against their producers.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition', 'WorkflowRun', 'ConsumerRecord']} })

    @field_validator('run_manifest_sha256')
    def pattern_run_manifest_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid run_manifest_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid run_manifest_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('feature_filter_id')
    def pattern_feature_filter_id(cls, v):
        pattern=re.compile(r"^sha256:[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid feature_filter_id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid feature_filter_id format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('coverage_scope')
    def pattern_coverage_scope(cls, v):
        pattern=re.compile(r"^input_cohort_before_analysis_specific_exclusions$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid coverage_scope format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid coverage_scope format: {v}"
            raise ValueError(err_msg)
        return v


class LegacyImportRecord(ConfiguredBaseModel):
    """
    Historical claims preserved without asserting a verified producing execution.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'Legacy Import Record'})

    id: str = Field(default=..., title="Id", description="""Identifier unique within its containing document.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive']} })
    source_path: str = Field(default=..., title="Source path", description="""Location of the historical source whose claims are being preserved.""", json_schema_extra = { "linkml_meta": {'domain_of': ['LegacyImportRecord']} })
    source_sha256: str = Field(default=..., title="Source sha256", description="""SHA-256 of the historical source bytes at import time.""", json_schema_extra = { "linkml_meta": {'domain_of': ['LegacyImportRecord']} })
    imported_at: datetime  = Field(default=..., title="Imported at", description="""UTC timestamp of legacy import, not the unknown historical execution time.""", json_schema_extra = { "linkml_meta": {'domain_of': ['LegacyImportRecord']} })
    verification_status: str = Field(default=..., title="Verification status", description="""Explicitly unverified historical execution; byte integrity alone is insufficient.""", json_schema_extra = { "linkml_meta": {'domain_of': ['LegacyImportRecord']} })
    claims_json: str = Field(default=..., title="Claims json", description="""Original historical claims preserved as JSON without inventing executed steps.""", json_schema_extra = { "linkml_meta": {'domain_of': ['LegacyImportRecord']} })

    @field_validator('source_sha256')
    def pattern_source_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid source_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid source_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('verification_status')
    def pattern_verification_status(cls, v):
        pattern=re.compile(r"^imported_unverified$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid verification_status format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid verification_status format: {v}"
            raise ValueError(err_msg)
        return v


class ArchivedLocation(ConfiguredBaseModel):
    """
    Original file or collection locations mapped to one immutable retained copy.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Archived Location'})

    original_paths: list[str] = Field(default=..., title="Original paths", description="""Absolute recorded locations, including any declared input alias; never rewritten in source records.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArchivedLocation']} })
    archive_path: str = Field(default=..., title="Archive path", description="""Safe relative location inside the archive, independent of its current root directory.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArchivedLocation']} })
    storage_kind: StorageKind = Field(default=..., title="Storage kind", description="""File or collection storage contract for the retained copy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArchivedLocation']} })
    sha256: str = Field(default=..., title="Sha256", description="""File digest or canonical collection inventory digest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance', 'ArchivedLocation']} })
    size_bytes: int = Field(default=..., title="Size bytes", description="""Total bytes of retained regular files.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance', 'ArchivedLocation']} })
    inventory: list[FileInventoryEntry] = Field(default=..., title="Inventory", description="""Complete regular-file inventory, including hidden members.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance', 'ArchivedLocation']} })

    @field_validator('sha256')
    def pattern_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class EvidenceArchive(ConfiguredBaseModel):
    """
    Retained bytes and explicit location mappings; archival verification does not invent historical execution or authorship.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Evidence Archive'})

    id: str = Field(default=..., title="Id", description="""Unique archive identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive']} })
    archive_version: str = Field(default=..., title="Archive version", description="""Version of the retention and verification contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['EvidenceArchive']} })
    created_at: datetime  = Field(default=..., title="Created at", description="""UTC archival time, distinct from original execution times.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord', 'EvidenceArchive']} })
    core_manifests: list[str] = Field(default=..., title="Core manifests", description="""Original absolute locations of retained core runs whose full scientific contracts are revalidated.""", json_schema_extra = { "linkml_meta": {'domain_of': ['EvidenceArchive']} })
    consumer_records: Optional[list[str]] = Field(default=None, title="Consumer records", description="""Original absolute locations of consumer records whose output and core bindings are revalidated.""", json_schema_extra = { "linkml_meta": {'domain_of': ['EvidenceArchive']} })
    retained_paths: Optional[list[str]] = Field(default=None, title="Retained paths", description="""Explicit supplemental evidence preserved by byte identity; preservation alone does not validate its scientific claims.""", json_schema_extra = { "linkml_meta": {'domain_of': ['EvidenceArchive']} })
    locations: list[ArchivedLocation] = Field(default=..., title="Locations", description="""Complete hash-checked source-to-retained-location mapping.""", json_schema_extra = { "linkml_meta": {'domain_of': ['EvidenceArchive']} })

    @field_validator('archive_version')
    def pattern_archive_version(cls, v):
        pattern=re.compile(r"^1$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid archive_version format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid archive_version format: {v}"
            raise ValueError(err_msg)
        return v


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
ProductDefinition.model_rebuild()
ParameterDefinition.model_rebuild()
ParameterValue.model_rebuild()
StepDefinition.model_rebuild()
WorkflowDefinition.model_rebuild()
FileInventoryEntry.model_rebuild()
ValidationResult.model_rebuild()
ArtifactInstance.model_rebuild()
StepExecution.model_rebuild()
WorkflowRun.model_rebuild()
ConsumerRecord.model_rebuild()
LegacyImportRecord.model_rebuild()
ArchivedLocation.model_rebuild()
EvidenceArchive.model_rebuild()
