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

class ContrastPopulation(str, Enum):
    """
    Executed experiment selection policy, separate from downstream trait eligibility.
    """
    all_paired = "all_paired"
    """
    Preserve the unfiltered core; downstream analyses retain their human pairing rules.
    """
    disease_vs_healthy_tier1_2_subject_controls = "disease_vs_healthy_tier1_2_subject_controls"
    """
    Paired human disease contrasts with Tier 1 or Tier 2 subject controls.
    """
    case_vs_control_subject = "case_vs_control_subject"
    """
    Paired human disease contrasts with subject controls, including Tier 3.
    """
    disease_vs_healthy_oriented = "disease_vs_healthy_oriented"
    """
    All paired human contrasts classified as disease versus healthy, including non-subject comparisons.
    """


class DirectionBasis(str, Enum):
    """
    Meaning of increased and decreased in the consumed core products.
    """
    group_1_vs_group_0 = "group_1_vs_group_0"
    """
    Original BugSigDB Group 1 relative to Group 0.
    """
    disease_vs_comparison = "disease_vs_comparison"
    """
    Disease group relative to the selected comparison, after explicit reversal where required.
    """


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


class CalibrationOutcome(str, Enum):
    """
    Standardized empirical calibration classification across analytical arms.
    """
    calibrated = "calibrated"
    """
    Null type-I error and FDR meet declared tolerances and power meets minimum on the exact production AnalysisSpec.
    """
    partially_calibrated = "partially_calibrated"
    """
    Some calibration criteria pass while others fail or remain unassessed on the evaluated AnalysisSpec.
    """
    anti_conservative = "anti_conservative"
    """
    Empirical type-I error or null FDR exceeds the declared tolerance upper bound.
    """
    null_controlled_low_power = "null_controlled_low_power"
    """
    Empirical null control passes, but power under the declared injected effect falls below minimum.
    """
    insufficient_budget = "insufficient_budget"
    """
    Replicate budget is below the protocol minimum required to certify calibration.
    """
    different_estimator_uncalibrated = "different_estimator_uncalibrated"
    """
    Calibration harness evaluated a different AnalysisSpec than the production estimator.
    """



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
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    description: Optional[str] = Field(default=None, title="Description", description="""Human-readable product meaning.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'StepDefinition', 'WorkflowDefinition']} })
    kind: str = Field(default=..., title="Kind", description="""Scientific validator contract name.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition']} })
    storage_kind: StorageKind = Field(default=..., title="Storage kind", description="""Whether the artifact is a regular file or an inventoried directory.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArchivedLocation', 'RecordedPath']} })
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
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    operation: str = Field(default=..., title="Operation", description="""Key resolving to an explicitly registered Python callable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition']} })
    description: Optional[str] = Field(default=None, title="Description", description="""Purpose of this step.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'StepDefinition', 'WorkflowDefinition']} })
    inputs: list[str] = Field(default=..., title="Inputs", description="""Product references whose verified artifact paths are passed to the operation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition', 'ArmCommandDefinition', 'ArmCommandExecution']} })
    outputs: list[str] = Field(default=..., title="Outputs", description="""Declared outputs; references are resolved and checked against their producers.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowRun',
                       'ConsumerRecord',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    depends_on: Optional[list[str]] = Field(default=None, title="Depends on", description="""Additional prerequisite step references beyond input-producing dependencies.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition']} })
    parameters: Optional[list[str]] = Field(default=None, title="Parameters", description="""Declared or effective parameters scoped to this definition or execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowDefinition',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmRun']} })

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
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    version: str = Field(default=..., title="Version", description="""Version of the declared workflow contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowDefinition']} })
    description: Optional[str] = Field(default=None, title="Description", description="""Workflow purpose and scope.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'StepDefinition', 'WorkflowDefinition']} })
    products: list[ProductDefinition] = Field(default=..., title="Products", description="""Scientific product contracts available to this workflow.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowDefinition']} })
    steps: list[StepDefinition] = Field(default=..., title="Steps", description="""Ordered or dependency-resolved workflow stages or their attempted executions.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowDefinition', 'WorkflowRun', 'ArmRun']} })
    parameters: Optional[list[ParameterDefinition]] = Field(default=None, title="Parameters", description="""Declared or effective parameters scoped to this definition or execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowDefinition',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmRun']} })

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

    path: str = Field(default=..., title="Path", description="""Relative member path, or the file basename.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance', 'RecordedPath']} })
    sha256: str = Field(default=..., title="Sha256", description="""Lowercase SHA-256 of file bytes or the canonical collection inventory.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry',
                       'ArtifactInstance',
                       'ArchivedLocation',
                       'RecordedPath']} })
    size_bytes: int = Field(default=..., title="Size bytes", description="""Total number of regular-file bytes represented by this record.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry',
                       'ArtifactInstance',
                       'ArchivedLocation',
                       'RecordedPath']} })

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
    status: ValidationStatus = Field(default=..., title="Status", description="""Explicit lifecycle or validation outcome; only completed verified runs are consumable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
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
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    product_id: str = Field(default=..., title="Product id", description="""Reference to the embedded product contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance']} })
    run_id: str = Field(default=..., title="Run id", description="""Unique execution-attempt identifier, distinct from the deterministic fingerprint.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance', 'WorkflowRun', 'ConsumerRecord']} })
    producer_step: Optional[str] = Field(default=None, title="Producer step", description="""Reference to the execution that produced this artifact; absent for inputs.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance'], 'slot_uri': 'prov:wasGeneratedBy'} })
    path: str = Field(default=..., title="Path", description="""Absolute external path or run-relative output path.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance', 'RecordedPath']} })
    requested_path: Optional[str] = Field(default=None, title="Requested path", description="""Original external binding before symlink resolution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance']} })
    sha256: str = Field(default=..., title="Sha256", description="""Lowercase SHA-256 of file bytes or the canonical collection inventory.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry',
                       'ArtifactInstance',
                       'ArchivedLocation',
                       'RecordedPath']} })
    size_bytes: int = Field(default=..., title="Size bytes", description="""Total number of regular-file bytes represented by this record.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry',
                       'ArtifactInstance',
                       'ArchivedLocation',
                       'RecordedPath']} })
    inventory: list[FileInventoryEntry] = Field(default=..., title="Inventory", description="""Complete sorted list of regular-file members and their byte hashes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance', 'ArchivedLocation', 'RecordedPath']} })
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
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    step_id: str = Field(default=..., title="Step id", description="""Reference to the declared step contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution']} })
    status: ExecutionStatus = Field(default=..., title="Status", description="""Explicit lifecycle or validation outcome; only completed verified runs are consumable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    started_at: datetime  = Field(default=..., title="Started at", description="""UTC timestamp when this execution attempt began.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    finished_at: Optional[datetime ] = Field(default=None, title="Finished at", description="""UTC timestamp when this execution attempt finished or failed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution', 'WorkflowRun']} })
    input_artifacts: Optional[list[str]] = Field(default=None, title="Input artifacts", description="""Concrete artifact references actually bound as operation inputs.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution'], 'slot_uri': 'prov:used'} })
    output_artifacts: Optional[list[str]] = Field(default=None, title="Output artifacts", description="""Concrete artifact references published by this completed operation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution'], 'slot_uri': 'prov:generated'} })
    parameters: Optional[list[ParameterValue]] = Field(default=None, title="Parameters", description="""Declared or effective parameters scoped to this definition or execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowDefinition',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmRun']} })
    error: Optional[str] = Field(default=None, title="Error", description="""Failure diagnostic; failed artifacts are never published.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })


class WorkflowRun(ConfiguredBaseModel):
    """
    Unique execution attempt with an embedded contract and content fingerprint.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'class_uri': 'prov:Activity',
         'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'Workflow Run'})

    run_id: str = Field(default=..., title="Run id", description="""Unique execution-attempt identifier, distinct from the deterministic fingerprint.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance', 'WorkflowRun', 'ConsumerRecord']} })
    contrast_population: Optional[ContrastPopulation] = Field(default=ContrastPopulation.all_paired, title="Contrast population", description="""Executed selection policy; historical records without this field retain the unfiltered policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun',
                       'ArmRun',
                       'EligibilityLedger',
                       'ScientificResultRecord',
                       'ComparisonRun'],
         'ifabsent': 'string(all_paired)'} })
    direction_basis: Optional[DirectionBasis] = Field(default=DirectionBasis.group_1_vs_group_0, title="Direction basis", description="""Direction semantics of selected Study signatures and all products derived from them.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun',
                       'ArmRun',
                       'ScientificResultRecord',
                       'ComparisonRun'],
         'ifabsent': 'string(group_1_vs_group_0)'} })
    definition: WorkflowDefinition = Field(default=..., title="Definition", description="""Embedded validated workflow definition used by this attempt.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    status: ExecutionStatus = Field(default=..., title="Status", description="""Explicit lifecycle or validation outcome; only completed verified runs are consumable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    started_at: datetime  = Field(default=..., title="Started at", description="""UTC timestamp when this execution attempt began.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    finished_at: Optional[datetime ] = Field(default=None, title="Finished at", description="""UTC timestamp when this execution attempt finished or failed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution', 'WorkflowRun']} })
    manifest_path: str = Field(default=..., title="Manifest path", description="""Location of this run manifest; relocation does not change execution identity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun', 'ComparisonRun']} })
    spec_path: str = Field(default=..., title="Spec path", description="""Source location of the workflow specification at execution time.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    spec_sha256: str = Field(default=..., title="Spec sha256", description="""SHA-256 of the exact source specification bytes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    definition_sha256: str = Field(default=..., title="Definition sha256", description="""SHA-256 of the canonical embedded validated workflow definition.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    execution_fingerprint: str = Field(default=..., title="Execution fingerprint", description="""Deterministic hash of contracts, inputs, code, environment, and resolved values.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    code_sha256: str = Field(default=..., title="Code sha256", description="""SHA-256 of the complete recorded code inventory.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    code_inventory: list[FileInventoryEntry] = Field(default=..., title="Code inventory", description="""Source and schema files, including untracked code, used to identify implementation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun']} })
    lock_sha256: str = Field(default=..., title="Lock sha256", description="""SHA-256 of the dependency lockfile used for this attempt.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    git_commit: str = Field(default=..., title="Git commit", description="""Repository commit at execution start; code hashes also capture local changes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun', 'ComparisonRun']} })
    git_dirty: bool = Field(default=..., title="Git dirty", description="""Whether tracked or untracked repository changes existed at execution start.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun', 'ComparisonRun']} })
    runtime_json: str = Field(default=..., title="Runtime json", description="""Canonical JSON describing Python, installed package versions, platform, and operations.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun']} })
    parameters: Optional[list[ParameterValue]] = Field(default=None, title="Parameters", description="""Declared or effective parameters scoped to this definition or execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowDefinition',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmRun']} })
    steps: Optional[list[StepExecution]] = Field(default=None, title="Steps", description="""Ordered or dependency-resolved workflow stages or their attempted executions.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowDefinition', 'WorkflowRun', 'ArmRun']} })
    artifacts: Optional[list[ArtifactInstance]] = Field(default=None, title="Artifacts", description="""Concrete external input and published output artifact instances.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun']} })
    outputs: Optional[list[str]] = Field(default=None, title="Outputs", description="""Declared outputs; references are resolved and checked against their producers.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowRun',
                       'ConsumerRecord',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    error: Optional[str] = Field(default=None, title="Error", description="""Failure diagnostic for the attempt.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })

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
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
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
    feature_filter_id: str = Field(default=..., title="Feature filter id", description="""SHA-256 identity of the exact ordered retained analysis features.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord', 'AnalysisBasis']} })
    coverage_scope: str = Field(default=..., title="Coverage scope", description="""Coverage describes the input cohort before analysis-specific eligibility exclusions.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord']} })
    outputs: list[FileInventoryEntry] = Field(default=..., title="Outputs", description="""Declared outputs; references are resolved and checked against their producers.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowRun',
                       'ConsumerRecord',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })

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
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
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
    storage_kind: StorageKind = Field(default=..., title="Storage kind", description="""File or collection storage contract for the retained copy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArchivedLocation', 'RecordedPath']} })
    sha256: str = Field(default=..., title="Sha256", description="""File digest or canonical collection inventory digest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry',
                       'ArtifactInstance',
                       'ArchivedLocation',
                       'RecordedPath']} })
    size_bytes: int = Field(default=..., title="Size bytes", description="""Total bytes of retained regular files.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry',
                       'ArtifactInstance',
                       'ArchivedLocation',
                       'RecordedPath']} })
    inventory: list[FileInventoryEntry] = Field(default=..., title="Inventory", description="""Complete regular-file inventory, including hidden members.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance', 'ArchivedLocation', 'RecordedPath']} })

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
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
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


class RecordedPath(ConfiguredBaseModel):
    """
    Observed file or directory snapshot; member bytes remain independently checkable.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Recorded Path'})

    path: str = Field(default=..., title="Path", description="""Absolute location observed during execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry', 'ArtifactInstance', 'RecordedPath']} })
    storage_kind: StorageKind = Field(default=..., title="Storage kind", description="""File or directory snapshot.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition', 'ArchivedLocation', 'RecordedPath']} })
    sha256: str = Field(default=..., title="Sha256", description="""File checksum or canonical inventory checksum.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry',
                       'ArtifactInstance',
                       'ArchivedLocation',
                       'RecordedPath']} })
    size_bytes: int = Field(default=..., title="Size bytes", description="""Observed total file size.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['FileInventoryEntry',
                       'ArtifactInstance',
                       'ArchivedLocation',
                       'RecordedPath']} })
    inventory: list[FileInventoryEntry] = Field(default=..., title="Inventory", description="""Exact members observed at this stage; later stages may add directory members.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArtifactInstance', 'ArchivedLocation', 'RecordedPath']} })

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


class ArmCommandDefinition(ConfiguredBaseModel):
    """
    Explicit command and required file contracts for one downstream stage.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'Arm Command Definition'})

    id: str = Field(default=..., title="Id", description="""Unique stage identifier within the arm.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    argv: list[str] = Field(default=..., title="Argv", description="""Exact argument vector executed without a shell.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmCommandDefinition']} })
    inputs: Optional[list[str]] = Field(default=None, title="Inputs", description="""Required input paths, explicitly bound to core or earlier stage products.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition', 'ArmCommandDefinition', 'ArmCommandExecution']} })
    outputs: list[str] = Field(default=..., title="Outputs", description="""Required output files or inventoried directories.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowRun',
                       'ConsumerRecord',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })

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


class ArmCommandExecution(ConfiguredBaseModel):
    """
    An attempted command with observed bytes and exit status.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'Arm Command Execution'})

    id: str = Field(default=..., title="Id", description="""Corresponding command definition identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    status: ExecutionStatus = Field(default=..., title="Status", description="""Actual execution state; failed or running stages are never complete.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    started_at: datetime  = Field(default=..., title="Started at", description="""UTC start time.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    ended_at: Optional[datetime ] = Field(default=None, title="Ended at", description="""UTC completion or failure time.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmCommandExecution', 'ArmRun', 'ComparisonRun']} })
    exit_code: Optional[int] = Field(default=None, title="Exit code", description="""Actual process exit code when a process was started.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmCommandExecution']} })
    inputs: Optional[list[RecordedPath]] = Field(default=None, title="Inputs", description="""Observed input snapshots before execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition', 'ArmCommandDefinition', 'ArmCommandExecution']} })
    outputs: Optional[list[RecordedPath]] = Field(default=None, title="Outputs", description="""Observed required outputs after successful execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowRun',
                       'ConsumerRecord',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    log: Optional[RecordedPath] = Field(default=None, title="Log", description="""Captured process output and its checksum.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmCommandExecution']} })
    error: Optional[str] = Field(default=None, title="Error", description="""Failure diagnostic, if any.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })


class ArmRun(ConfiguredBaseModel):
    """
    Execution-generated downstream rerun bound to a verified core bundle.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Arm Run'})

    id: str = Field(default=..., title="Id", description="""Unique attempt identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    arm: str = Field(default=..., title="Arm", description="""Owner of the declared standard analysis suite.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun']} })
    contrast_population: Optional[ContrastPopulation] = Field(default=ContrastPopulation.all_paired, title="Contrast population", description="""Selection inherited from the verified upstream core, never independently overridden.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun',
                       'ArmRun',
                       'EligibilityLedger',
                       'ScientificResultRecord',
                       'ComparisonRun'],
         'ifabsent': 'string(all_paired)'} })
    direction_basis: Optional[DirectionBasis] = Field(default=DirectionBasis.group_1_vs_group_0, title="Direction basis", description="""Upstream direction semantics; disease-oriented inputs must not be reversed a second time.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun',
                       'ArmRun',
                       'ScientificResultRecord',
                       'ComparisonRun'],
         'ifabsent': 'string(group_1_vs_group_0)'} })
    scope: str = Field(default=..., title="Scope", description="""Explicit scope; historical PD and presentation investigations are separate.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun']} })
    status: ExecutionStatus = Field(default=..., title="Status", description="""Actual attempt status.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    started_at: datetime  = Field(default=..., title="Started at", description="""UTC start time.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    ended_at: Optional[datetime ] = Field(default=None, title="Ended at", description="""UTC completion or failure time.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmCommandExecution', 'ArmRun', 'ComparisonRun']} })
    manifest_path: str = Field(default=..., title="Manifest path", description="""Absolute location of this execution record.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun', 'ComparisonRun']} })
    core_manifest_path: str = Field(default=..., title="Core manifest path", description="""Verified upstream run manifest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun', 'ComparisonRun']} })
    core_manifest_sha256: str = Field(default=..., title="Core manifest sha256", description="""Upstream manifest checksum before execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun', 'ComparisonRun']} })
    core_run_id: str = Field(default=..., title="Core run id", description="""Upstream execution identity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun', 'ComparisonRun']} })
    profile_product: str = Field(default=..., title="Profile product", description="""Scientific count unit used by the primary analysis.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun']} })
    plan_sha256: str = Field(default=..., title="Plan sha256", description="""Canonical checksum of the complete ordered command definitions.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun']} })
    git_commit: str = Field(default=..., title="Git commit", description="""Actual checkout revision.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun', 'ComparisonRun']} })
    git_dirty: bool = Field(default=..., title="Git dirty", description="""Whether the checkout contained local changes at start.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun', 'ComparisonRun']} })
    code_inventory: list[FileInventoryEntry] = Field(default=..., title="Code inventory", description="""Actual source, schema, methods and workflow configuration bytes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun']} })
    repository_path: str = Field(default=..., title="Repository path", description="""Source checkout used to execute the commands.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun']} })
    runtime_json: str = Field(default=..., title="Runtime json", description="""Python, platform, package versions and explicitly controlled process settings.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun']} })
    parameters: list[ParameterValue] = Field(default=..., title="Parameters", description="""Resolved arm options; absent overrides preserve stage-specific defaults.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowDefinition',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmRun']} })
    prepared_inputs: Optional[list[RecordedPath]] = Field(default=None, title="Prepared inputs", description="""Fresh byte-identical copies of verified core products under canonical filenames.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun']} })
    plan: list[ArmCommandDefinition] = Field(default=..., title="Plan", description="""Complete ordered command plan recorded before execution.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun']} })
    steps: Optional[list[ArmCommandExecution]] = Field(default=None, title="Steps", description="""Commands actually attempted in order.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowDefinition', 'WorkflowRun', 'ArmRun']} })
    outputs: Optional[list[FileInventoryEntry]] = Field(default=None, title="Outputs", description="""Exact final output-tree inventory excluding this manifest and its temporary write.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowRun',
                       'ConsumerRecord',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    error: Optional[str] = Field(default=None, title="Error", description="""Failure diagnostic; no failed attempt can verify as completed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })

    @field_validator('arm')
    def pattern_arm(cls, v):
        pattern=re.compile(r"^(ilyas|noel|jed)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid arm format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid arm format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('scope')
    def pattern_scope(cls, v):
        pattern=re.compile(r"^standard_analysis_calibration_and_optional_figures$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid scope format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid scope format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('core_manifest_sha256')
    def pattern_core_manifest_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid core_manifest_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid core_manifest_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('profile_product')
    def pattern_profile_product(cls, v):
        pattern=re.compile(r"^(count_profiles|ora_profiles)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid profile_product format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid profile_product format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('plan_sha256')
    def pattern_plan_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid plan_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid plan_sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class AnalysisSpec(ConfiguredBaseModel):
    """
    Shared executable specification of estimand, kernel, transformation, sampling unit, weighting, and eligibility.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Analysis Spec'})

    id: str = Field(default=..., title="Id", description="""Stable analysis specification identifier shared by production, calibration, and sensitivity evaluations.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    estimand: str = Field(default=..., title="Estimand", description="""Scientific quantity estimated by the analysis.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisSpec']} })
    kernel: str = Field(default=..., title="Kernel", description="""Statistical test or summary kernel applied to observations.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisSpec']} })
    transformation: str = Field(default=..., title="Transformation", description="""Scale transformation applied before aggregation or testing.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisSpec']} })
    sampling_unit: str = Field(default=..., title="Sampling unit", description="""Observational unit at which independent or paired records are evaluated.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisSpec']} })
    weighting: str = Field(default=..., title="Weighting", description="""Observation or study weighting rule.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisSpec']} })
    eligibility_rule: str = Field(default=..., title="Eligibility rule", description="""Rule determining which contrasts, studies, or signatures enter the denominator.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisSpec']} })
    multiple_testing_policy: str = Field(default=..., title="Multiple testing policy", description="""Multiplicity correction procedure and family scope.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisSpec']} })

    @field_validator('id')
    def pattern_id(cls, v):
        pattern=re.compile(r"^[A-Za-z][A-Za-z0-9_:.-]*$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid id format: {v}"
            raise ValueError(err_msg)
        return v


class EligibilityLedger(ConfiguredBaseModel):
    """
    Explicit accounting of included sampling units, exclusion reasons, and common-support counts.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Eligibility Ledger'})

    id: str = Field(default=..., title="Id", description="""Unique eligibility ledger identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    analysis_spec_id: str = Field(default=..., title="Analysis spec id", description="""Analysis specification whose eligibility rules produced this ledger.""", json_schema_extra = { "linkml_meta": {'domain_of': ['EligibilityLedger', 'ScientificResultRecord']} })
    contrast_population: ContrastPopulation = Field(default=..., title="Contrast population", description="""Upstream contrast population before analysis-specific eligibility filtering.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun',
                       'ArmRun',
                       'EligibilityLedger',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    included_count: int = Field(default=..., title="Included count", description="""Number of sampling units retained after eligibility filtering.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['EligibilityLedger']} })
    excluded_count: int = Field(default=..., title="Excluded count", description="""Number of candidate sampling units excluded by the eligibility rule.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['EligibilityLedger']} })
    included_ids: Optional[list[str]] = Field(default=None, title="Included ids", description="""Deterministic sorted identifiers of included sampling units.""", json_schema_extra = { "linkml_meta": {'domain_of': ['EligibilityLedger']} })
    exclusion_counts_json: str = Field(default=..., title="Exclusion counts json", description="""Canonical JSON mapping each exclusion reason code to its count.""", json_schema_extra = { "linkml_meta": {'domain_of': ['EligibilityLedger']} })
    common_support_count: Optional[int] = Field(default=None, title="Common support count", description="""Number of sampling units shared with a paired comparator specification when evaluated on common support.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['EligibilityLedger']} })


class AnalysisBasis(ConfiguredBaseModel):
    """
    Feature catalog, vocabulary filter, rollup policy, denominator, background universe, and weighting basis.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Analysis Basis'})

    id: str = Field(default=..., title="Id", description="""Unique analysis basis identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    feature_catalog_sha256: str = Field(default=..., title="Feature catalog sha256", description="""SHA-256 of the authoritative trait feature catalog.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisBasis']} })
    feature_filter_id: str = Field(default=..., title="Feature filter id", description="""SHA-256 identity of the exact ordered retained analysis features.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ConsumerRecord', 'AnalysisBasis']} })
    retained_feature_count: int = Field(default=..., title="Retained feature count", description="""Number of retained qualified trait features in the analysis family.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisBasis']} })
    profile_unit: str = Field(default=..., title="Profile unit", description="""Underlying profile representation unit.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisBasis']} })
    aggregation_policy: str = Field(default=..., title="Aggregation policy", description="""How taxon-level annotations are aggregated into signature, contrast, or study summaries.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisBasis']} })
    rollup_policy: str = Field(default=..., title="Rollup policy", description="""Taxonomic trait inheritance and strain rollup policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisBasis']} })
    denominator_policy: str = Field(default=..., title="Denominator policy", description="""Denominator definition for fractions, ratios, or hypergeometric backgrounds.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisBasis']} })
    background_inventory_id: str = Field(default=..., title="Background inventory id", description="""Identifier or digest of the background taxon universe when applicable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisBasis']} })
    weighting_policy: str = Field(default=..., title="Weighting policy", description="""Study, contrast, or taxon weighting policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AnalysisBasis']} })

    @field_validator('feature_catalog_sha256')
    def pattern_feature_catalog_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid feature_catalog_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid feature_catalog_sha256 format: {v}"
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

    @field_validator('profile_unit')
    def pattern_profile_unit(cls, v):
        pattern=re.compile(r"^(count_profiles|ora_profiles)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid profile_unit format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid profile_unit format: {v}"
            raise ValueError(err_msg)
        return v


class CalibrationVerdict(ConfiguredBaseModel):
    """
    Standardized calibration assessment binding the evaluated harness specification to the production estimator.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Calibration Verdict'})

    id: str = Field(default=..., title="Id", description="""Unique calibration verdict identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    production_analysis_spec_id: str = Field(default=..., title="Production analysis spec id", description="""Analysis specification of the production estimator being assessed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    evaluated_analysis_spec_id: str = Field(default=..., title="Evaluated analysis spec id", description="""Analysis specification actually executed by the calibration harness.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    estimator_match: bool = Field(default=..., title="Estimator match", description="""Whether the calibration harness evaluated the exact production AnalysisSpec.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    outcome: CalibrationOutcome = Field(default=..., title="Outcome", description="""Standardized calibration classification across arms.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    null_control_passed: Optional[bool] = Field(default=None, title="Null control passed", description="""Whether empirical type-I error and null FDR stayed within declared tolerances.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    power_passed: Optional[bool] = Field(default=None, title="Power passed", description="""Whether empirical power under the declared injected effect met the minimum threshold.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    replicate_budget_sufficient: bool = Field(default=..., title="Replicate budget sufficient", description="""Whether the executed replicate count met the protocol minimum.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    summary_message: str = Field(default=..., title="Summary message", description="""Human-readable explanation of the calibration verdict and any estimator mismatch.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    negative_control_procedure: Optional[str] = Field(default=None, title="Negative control procedure", description="""Harness procedure declared in advance to fail the acceptance criteria, recorded apart from the production procedure's own error control.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    negative_control_outcome: Optional[str] = Field(default=None, title="Negative control outcome", description="""Overall status the calibration harness reported for the negative control procedure.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    negative_control_as_expected: Optional[bool] = Field(default=None, title="Negative control as expected", description="""Whether the negative control failed, as declared, every criterion the production verdict relies on. False means it passed one, which is then too weak to certify the production procedure; unset means its run did not settle it.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })
    selection_uncertainty: Optional[str] = Field(default=None, title="Selection uncertainty", description="""Uncertainty from selecting the analysed features on the same corpus, which calibrating the fixed-feature procedure does not resolve in either direction.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CalibrationVerdict']} })


class ScientificResultRecord(ConfiguredBaseModel):
    """
    Standardized scientific result row bound to explicit analysis specification, eligibility ledger, and basis contracts.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow',
         'title': 'Scientific Result Record'})

    id: str = Field(default=..., title="Id", description="""Unique scientific result record identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    contrast_population: ContrastPopulation = Field(default=..., title="Contrast population", description="""Executed experiment selection policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun',
                       'ArmRun',
                       'EligibilityLedger',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    direction_basis: DirectionBasis = Field(default=..., title="Direction basis", description="""Meaning of increased and decreased in the consumed core products.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun',
                       'ArmRun',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    analysis_spec_id: str = Field(default=..., title="Analysis spec id", description="""Reference to the exact AnalysisSpec used to compute this result.""", json_schema_extra = { "linkml_meta": {'domain_of': ['EligibilityLedger', 'ScientificResultRecord']} })
    eligibility_ledger_id: str = Field(default=..., title="Eligibility ledger id", description="""Reference to the EligibilityLedger recording included and excluded units.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    analysis_basis_id: str = Field(default=..., title="Analysis basis id", description="""Reference to the AnalysisBasis recording feature vocabulary, rollup, and denominator rules.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    calibration_verdict_id: Optional[str] = Field(default=None, title="Calibration verdict id", description="""Reference to the CalibrationVerdict for this estimator, if assessed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    feature_id: str = Field(default=..., title="Feature id", description="""Qualified trait identifier or pair identifier, never a bare display label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    feature_label: str = Field(default=..., title="Feature label", description="""Human-readable display label associated with the qualified feature identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    effect_estimate: Optional[float] = Field(default=None, title="Effect estimate", description="""Primary point estimate on the scale declared by the AnalysisSpec.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    effect_scale: str = Field(default=..., title="Effect scale", description="""Scale and orientation of the point estimate.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    p_value: Optional[float] = Field(default=None, title="P value", description="""Raw two-sided or one-sided p-value as declared by the AnalysisSpec.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    q_value: Optional[float] = Field(default=None, title="Q value", description="""Multiplicity-adjusted q-value within the declared family scope.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    effective_q_value: Optional[float] = Field(default=None, title="Effective q value", description="""Eigenvalue-effective family-size adjusted q-value when applicable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    n_eligible: int = Field(default=..., title="N eligible", description="""Number of eligible sampling units contributing to this test.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    positive_count: Optional[int] = Field(default=None, title="Positive count", description="""Count of positive shifts or enriched signatures.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    negative_count: Optional[int] = Field(default=None, title="Negative count", description="""Count of negative shifts or depleted signatures.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })
    tie_or_zero_count: Optional[int] = Field(default=None, title="Tie or zero count", description="""Count of exact ties or zero-total excluded units.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ScientificResultRecord']} })


class ComparisonRun(ConfiguredBaseModel):
    """
    Complete cross-arm integration record consuming verified core, Noel, Jed, and Ilyas arm manifests without missing-source placeholders.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/workflow', 'title': 'Comparison Run'})

    id: str = Field(default=..., title="Id", description="""Unique cross-arm comparison attempt identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProductDefinition',
                       'StepDefinition',
                       'WorkflowDefinition',
                       'ArtifactInstance',
                       'StepExecution',
                       'ConsumerRecord',
                       'LegacyImportRecord',
                       'EvidenceArchive',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'AnalysisSpec',
                       'EligibilityLedger',
                       'AnalysisBasis',
                       'CalibrationVerdict',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    contrast_population: ContrastPopulation = Field(default=..., title="Contrast population", description="""Shared contrast population verified across core and all three arms.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun',
                       'ArmRun',
                       'EligibilityLedger',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    direction_basis: DirectionBasis = Field(default=..., title="Direction basis", description="""Shared direction basis verified across core and all three arms.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun',
                       'ArmRun',
                       'ScientificResultRecord',
                       'ComparisonRun']} })
    status: ExecutionStatus = Field(default=..., title="Status", description="""Explicit lifecycle outcome of the cross-arm comparison run.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ValidationResult',
                       'StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    started_at: datetime  = Field(default=..., title="Started at", description="""UTC timestamp when the cross-arm integration began.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    ended_at: Optional[datetime ] = Field(default=None, title="Ended at", description="""UTC timestamp when the cross-arm integration finished or failed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmCommandExecution', 'ArmRun', 'ComparisonRun']} })
    manifest_path: str = Field(default=..., title="Manifest path", description="""Absolute location of this ComparisonRun manifest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun', 'ComparisonRun']} })
    core_manifest_path: str = Field(default=..., title="Core manifest path", description="""Verified upstream core run manifest path.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun', 'ComparisonRun']} })
    core_manifest_sha256: str = Field(default=..., title="Core manifest sha256", description="""SHA-256 of the verified upstream core run manifest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun', 'ComparisonRun']} })
    core_run_id: str = Field(default=..., title="Core run id", description="""Verified upstream core run identifier shared by all three arms.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ArmRun', 'ComparisonRun']} })
    noel_manifest_path: str = Field(default=..., title="Noel manifest path", description="""Verified Noel ArmRun manifest path.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ComparisonRun']} })
    noel_manifest_sha256: str = Field(default=..., title="Noel manifest sha256", description="""SHA-256 of the verified Noel ArmRun manifest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ComparisonRun']} })
    jed_manifest_path: str = Field(default=..., title="Jed manifest path", description="""Verified Jed ArmRun manifest path.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ComparisonRun']} })
    jed_manifest_sha256: str = Field(default=..., title="Jed manifest sha256", description="""SHA-256 of the verified Jed ArmRun manifest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ComparisonRun']} })
    ilyas_manifest_path: str = Field(default=..., title="Ilyas manifest path", description="""Verified Ilyas ArmRun manifest path.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ComparisonRun']} })
    ilyas_manifest_sha256: str = Field(default=..., title="Ilyas manifest sha256", description="""SHA-256 of the verified Ilyas ArmRun manifest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ComparisonRun']} })
    git_commit: str = Field(default=..., title="Git commit", description="""Repository commit at execution start.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun', 'ComparisonRun']} })
    git_dirty: bool = Field(default=..., title="Git dirty", description="""Whether tracked or untracked changes existed at execution start.""", json_schema_extra = { "linkml_meta": {'domain_of': ['WorkflowRun', 'ArmRun', 'ComparisonRun']} })
    analysis_specs: list[AnalysisSpec] = Field(default=..., title="Analysis specs", description="""Declared AnalysisSpec contracts evaluated or compared in this run.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ComparisonRun']} })
    eligibility_ledgers: list[EligibilityLedger] = Field(default=..., title="Eligibility ledgers", description="""Eligibility ledgers for each evaluated analysis specification.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ComparisonRun']} })
    analysis_bases: list[AnalysisBasis] = Field(default=..., title="Analysis bases", description="""Feature, rollup, denominator, and background bases for each arm.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ComparisonRun']} })
    calibration_verdicts: list[CalibrationVerdict] = Field(default=..., title="Calibration verdicts", description="""Standardized calibration verdicts with explicit estimator-match checks.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ComparisonRun']} })
    outputs: list[FileInventoryEntry] = Field(default=..., title="Outputs", description="""Complete output file inventory produced by the cross-arm comparison stage.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepDefinition',
                       'WorkflowRun',
                       'ConsumerRecord',
                       'ArmCommandDefinition',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })
    error: Optional[str] = Field(default=None, title="Error", description="""Failure diagnostic if the comparison run failed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StepExecution',
                       'WorkflowRun',
                       'ArmCommandExecution',
                       'ArmRun',
                       'ComparisonRun']} })

    @field_validator('core_manifest_sha256')
    def pattern_core_manifest_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid core_manifest_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid core_manifest_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('noel_manifest_sha256')
    def pattern_noel_manifest_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid noel_manifest_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid noel_manifest_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('jed_manifest_sha256')
    def pattern_jed_manifest_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid jed_manifest_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid jed_manifest_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('ilyas_manifest_sha256')
    def pattern_ilyas_manifest_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid ilyas_manifest_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid ilyas_manifest_sha256 format: {v}"
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
RecordedPath.model_rebuild()
ArmCommandDefinition.model_rebuild()
ArmCommandExecution.model_rebuild()
ArmRun.model_rebuild()
AnalysisSpec.model_rebuild()
EligibilityLedger.model_rebuild()
AnalysisBasis.model_rebuild()
CalibrationVerdict.model_rebuild()
ScientificResultRecord.model_rebuild()
ComparisonRun.model_rebuild()
