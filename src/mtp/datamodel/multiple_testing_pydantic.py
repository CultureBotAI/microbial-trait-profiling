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
version = "1.1.0"


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


linkml_meta = LinkMLMeta({'default_prefix': 'mtpmt',
     'default_range': 'string',
     'description': 'Scientific family inventories and replayable corrections. '
                    'These records do not change native arm calls or certify the '
                    'validity of their input p-values.',
     'id': 'https://w3id.org/mtp/multiple-testing',
     'imports': ['linkml:types'],
     'license': 'MIT',
     'name': 'mtp_multiple_testing',
     'prefixes': {'linkml': {'prefix_prefix': 'linkml',
                             'prefix_reference': 'https://w3id.org/linkml/'},
                  'mtpmt': {'prefix_prefix': 'mtpmt',
                            'prefix_reference': 'https://w3id.org/mtp/multiple-testing/'}},
     'source_file': 'schema/multiple_testing.yaml',
     'title': 'Additive full-family multiple-testing contracts'} )


class MultipleTestingPolicy(ConfiguredBaseModel):
    """
    Versioned immutable correction policy checked against the supported policy file.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'Multiple testing policy'})

    policy_id: str = Field(default=..., title="Policy id", description="""Stable versioned policy identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['MultipleTestingPolicy']} })
    version: str = Field(default=..., title="Version", description="""Semantic policy version.""", json_schema_extra = { "linkml_meta": {'domain_of': ['MultipleTestingPolicy', 'PostArmMultiplicityRun']} })
    role: str = Field(default=..., title="Role", description="""Additive status; native arm calls remain unchanged.""", json_schema_extra = { "linkml_meta": {'domain_of': ['MultipleTestingPolicy']} })
    primary_method: str = Field(default=..., title="Primary method", description="""Correction used as primary within this additive sidecar.""", json_schema_extra = { "linkml_meta": {'domain_of': ['MultipleTestingPolicy']} })
    comparison_method: str = Field(default=..., title="Comparison method", description="""Companion correction for comparison only.""", json_schema_extra = { "linkml_meta": {'domain_of': ['MultipleTestingPolicy']} })
    thresholds: list[float] = Field(default=..., title="Thresholds", description="""Inclusive adjusted-p thresholds reported by the policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['MultipleTestingPolicy']} })
    missing_test_policy: str = Field(default=..., title="Missing test policy", description="""Explicit untestable members reserve internal p one and export undefined p and q.""", json_schema_extra = { "linkml_meta": {'domain_of': ['MultipleTestingPolicy']} })
    family_policy: str = Field(default=..., title="Family policy", description="""Full declared membership is mandatory, without significance or effect filtering.""", json_schema_extra = { "linkml_meta": {'domain_of': ['MultipleTestingPolicy']} })


class SourceArtifact(ConfiguredBaseModel):
    """
    Logical identity and byte hash of an upstream source, verified by the calling adapter.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'Source artifact'})

    identity: str = Field(default=..., title="Identity", description="""Stable upstream product or execution identifier, not a machine-specific path.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceArtifact']} })
    sha256: str = Field(default=..., title="Sha256", description="""Lowercase SHA256 of the referenced source bytes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceArtifact']} })

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


class HypothesisCoordinates(ConfiguredBaseModel):
    """
    Explicit non-feature coordinates. Empty coordinates identify a pooled test. Coordinates identify scientific tests independently of arbitrary row IDs and eligibility ledger hashes.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'Hypothesis coordinates'})

    profile_id: Optional[str] = Field(default=None, title="Profile id", description="""Stable source profile identity, for example study, experiment and direction.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisCoordinates']} })
    stratum_id: Optional[str] = Field(default=None, title="Stratum id", description="""Source stratum identity or label, for example an anatomical site.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisCoordinates']} })
    condition_ids: Optional[list[str]] = Field(default=None, title="Condition ids", description="""Distinct source condition identifiers or exact labels, canonicalized in sorted order. These need not be ontology-mapped identifiers. Adapters must verify the expected cardinality, for example three for a triple.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisCoordinates']} })
    split_id: Optional[str] = Field(default=None, title="Split id", description="""Declared training or held-out split identity, where applicable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisCoordinates']} })


class HypothesisDefinition(ConfiguredBaseModel):
    """
    A declared member of a scientific family, before its observed p-value is read.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'Hypothesis definition'})

    hypothesis_id: str = Field(default=..., title="Hypothesis id", description="""Unique stable hypothesis identifier within a family.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisDefinition', 'RawHypothesis']} })
    feature_ids: list[str] = Field(default=..., title="Feature ids", description="""Ordered qualified feature identities; ratios retain numerator then denominator.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisDefinition']} })
    coordinates: HypothesisCoordinates = Field(default=..., title="Coordinates", description="""Explicit profile, stratum, condition-set and split coordinates; empty for pooled tests.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisDefinition']} })
    test: str = Field(default=..., title="Test", description="""Statistical test that produced this hypothesis p-value.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisDefinition']} })
    tail: str = Field(default=..., title="Tail", description="""Test alternative; independent of a curated abundance direction.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisDefinition']} })
    eligibility_sha256: Optional[str] = Field(default=None, title="Eligibility sha256", description="""Optional identity of this hypothesis-specific eligibility ledger.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisDefinition']} })

    @field_validator('eligibility_sha256')
    def pattern_eligibility_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid eligibility_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid eligibility_sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class HypothesisFamilyMetadata(ConfiguredBaseModel):
    """
    Complete declared inventory and scientific scope. Expected size is supplied from verified execution by an adapter; the engine does not infer it from observed results or use a fixed retained feature count.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'Hypothesis family metadata'})

    family_id: str = Field(default=..., title="Family id", description="""Stable scientific family identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    analysis_id: str = Field(default=..., title="Analysis id", description="""Analysis or source arm execution identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    analysis_basis_sha256: str = Field(default=..., title="Analysis basis sha256", description="""Canonical identity of the verified upstream AnalysisBasis.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    source_membership_sha256: str = Field(default=..., title="Source membership sha256", description="""Exact cohort or eligibility membership identity supplied by the verified adapter.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    retained_feature_count: int = Field(default=..., title="Retained feature count", description="""Feature count from AnalysisBasis; distinct from the number of hypotheses.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    strain_rollup: str = Field(default=..., title="Strain rollup", description="""Verified annotation scope; full rollup is not supported by this contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata', 'PostArmMultiplicityRun']} })
    feature_filter_id: str = Field(default=..., title="Feature filter id", description="""Ordered retained-feature identity copied from the verified AnalysisBasis.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    expected_family_size: int = Field(default=..., title="Expected family size", description="""Complete hypothesis count including explicitly untestable members.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    population: str = Field(default=..., title="Population", description="""Verified contrast population and any anatomical or condition stratum.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    direction_basis: str = Field(default=..., title="Direction basis", description="""Curated group direction or verified disease-control orientation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata', 'PostArmMultiplicityRun']} })
    profile_unit: str = Field(default=..., title="Profile unit", description="""Source profile representation copied verbatim from AnalysisBasis, such as count_profiles or ora_profiles.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    inference_unit: str = Field(default=..., title="Inference unit", description="""Independent unit assumed by the source test, for example study or taxon.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    hypothesis_unit: str = Field(default=..., title="Hypothesis unit", description="""Family coordinate, for example trait, ordered ratio, or profile by trait by tail.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    estimand: str = Field(default=..., title="Estimand", description="""Scientific quantity tested, preserving annotation versus physiological measurement.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    selection_rule: str = Field(default=..., title="Selection rule", description="""Declared membership rule, independent of display and effect-size filters.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    selection_history: list[str] = Field(default=..., title="Selection history", description="""Chronological disclosures of data-driven selection; correction alone does not remove selection bias.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    sources: list[SourceArtifact] = Field(default=..., title="Sources", description="""Verified upstream artifacts; identities must be unique.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata', 'PostArmMultiplicityRun']} })

    @field_validator('analysis_basis_sha256')
    def pattern_analysis_basis_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid analysis_basis_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid analysis_basis_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('source_membership_sha256')
    def pattern_source_membership_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid source_membership_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid source_membership_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('strain_rollup')
    def pattern_strain_rollup(cls, v):
        pattern=re.compile(r"^(off|phenotype)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid strain_rollup format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid strain_rollup format: {v}"
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


class HypothesisFamily(HypothesisFamilyMetadata):
    """
    Complete in-memory declaration used for small scientific families.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'Hypothesis family'})

    members: list[HypothesisDefinition] = Field(default=..., title="Members", description="""Complete declared hypothesis inventory; empty is allowed only with expected size zero.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamily']} })
    family_id: str = Field(default=..., title="Family id", description="""Stable scientific family identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    analysis_id: str = Field(default=..., title="Analysis id", description="""Analysis or source arm execution identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    analysis_basis_sha256: str = Field(default=..., title="Analysis basis sha256", description="""Canonical identity of the verified upstream AnalysisBasis.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    source_membership_sha256: str = Field(default=..., title="Source membership sha256", description="""Exact cohort or eligibility membership identity supplied by the verified adapter.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    retained_feature_count: int = Field(default=..., title="Retained feature count", description="""Feature count from AnalysisBasis; distinct from the number of hypotheses.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    strain_rollup: str = Field(default=..., title="Strain rollup", description="""Verified annotation scope; full rollup is not supported by this contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata', 'PostArmMultiplicityRun']} })
    feature_filter_id: str = Field(default=..., title="Feature filter id", description="""Ordered retained-feature identity copied from the verified AnalysisBasis.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    expected_family_size: int = Field(default=..., title="Expected family size", description="""Complete hypothesis count including explicitly untestable members.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    population: str = Field(default=..., title="Population", description="""Verified contrast population and any anatomical or condition stratum.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    direction_basis: str = Field(default=..., title="Direction basis", description="""Curated group direction or verified disease-control orientation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata', 'PostArmMultiplicityRun']} })
    profile_unit: str = Field(default=..., title="Profile unit", description="""Source profile representation copied verbatim from AnalysisBasis, such as count_profiles or ora_profiles.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    inference_unit: str = Field(default=..., title="Inference unit", description="""Independent unit assumed by the source test, for example study or taxon.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    hypothesis_unit: str = Field(default=..., title="Hypothesis unit", description="""Family coordinate, for example trait, ordered ratio, or profile by trait by tail.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    estimand: str = Field(default=..., title="Estimand", description="""Scientific quantity tested, preserving annotation versus physiological measurement.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    selection_rule: str = Field(default=..., title="Selection rule", description="""Declared membership rule, independent of display and effect-size filters.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    selection_history: list[str] = Field(default=..., title="Selection history", description="""Chronological disclosures of data-driven selection; correction alone does not remove selection bias.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata']} })
    sources: list[SourceArtifact] = Field(default=..., title="Sources", description="""Verified upstream artifacts; identities must be unique.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata', 'PostArmMultiplicityRun']} })

    @field_validator('analysis_basis_sha256')
    def pattern_analysis_basis_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid analysis_basis_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid analysis_basis_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('source_membership_sha256')
    def pattern_source_membership_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid source_membership_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid source_membership_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('strain_rollup')
    def pattern_strain_rollup(cls, v):
        pattern=re.compile(r"^(off|phenotype)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid strain_rollup format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid strain_rollup format: {v}"
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


class RawHypothesis(ConfiguredBaseModel):
    """
    One executed test or an explicitly untestable reserved member.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'Raw hypothesis'})

    hypothesis_id: str = Field(default=..., title="Hypothesis id", description="""Identifier matching exactly one declared family member.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisDefinition', 'RawHypothesis']} })
    status: str = Field(default=..., title="Status", description="""Tested or explicitly untestable, never an inferred missing row.""", json_schema_extra = { "linkml_meta": {'domain_of': ['RawHypothesis']} })
    p_value: Optional[float] = Field(default=None, title="P value", description="""Finite probability for a tested member; undefined for an untestable member.""", ge=0, le=1, json_schema_extra = { "linkml_meta": {'domain_of': ['RawHypothesis']} })
    reason: Optional[str] = Field(default=None, title="Reason", description="""Nonempty explanation required only for an untestable member.""", json_schema_extra = { "linkml_meta": {'domain_of': ['RawHypothesis']} })

    @field_validator('status')
    def pattern_status(cls, v):
        pattern=re.compile(r"^(tested|untestable)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid status format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid status format: {v}"
            raise ValueError(err_msg)
        return v


class AdjustedHypothesis(RawHypothesis):
    """
    Recomputed full-family corrections and threshold decisions; no effect filter is applied.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'Adjusted hypothesis'})

    q_by: Optional[float] = Field(default=None, title="Q by", description="""BY adjusted p-value; undefined for untestable members.""", ge=0, le=1, json_schema_extra = { "linkml_meta": {'domain_of': ['AdjustedHypothesis']} })
    q_bh: Optional[float] = Field(default=None, title="Q bh", description="""BH adjusted p-value; undefined for untestable members.""", ge=0, le=1, json_schema_extra = { "linkml_meta": {'domain_of': ['AdjustedHypothesis']} })
    reject_by_005: Optional[bool] = Field(default=None, title="Reject by 005", description="""BY adjusted p at most 0.05; undefined for untestable members.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AdjustedHypothesis']} })
    reject_by_0001: Optional[bool] = Field(default=None, title="Reject by 0001", description="""BY adjusted p at most 0.001; undefined for untestable members.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AdjustedHypothesis']} })
    reject_bh_005: Optional[bool] = Field(default=None, title="Reject bh 005", description="""BH adjusted p at most 0.05; undefined for untestable members.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AdjustedHypothesis']} })
    reject_bh_0001: Optional[bool] = Field(default=None, title="Reject bh 0001", description="""BH adjusted p at most 0.001; undefined for untestable members.""", json_schema_extra = { "linkml_meta": {'domain_of': ['AdjustedHypothesis']} })
    hypothesis_id: str = Field(default=..., title="Hypothesis id", description="""Identifier matching exactly one declared family member.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisDefinition', 'RawHypothesis']} })
    status: str = Field(default=..., title="Status", description="""Tested or explicitly untestable, never an inferred missing row.""", json_schema_extra = { "linkml_meta": {'domain_of': ['RawHypothesis']} })
    p_value: Optional[float] = Field(default=None, title="P value", description="""Finite probability for a tested member; undefined for an untestable member.""", ge=0, le=1, json_schema_extra = { "linkml_meta": {'domain_of': ['RawHypothesis']} })
    reason: Optional[str] = Field(default=None, title="Reason", description="""Nonempty explanation required only for an untestable member.""", json_schema_extra = { "linkml_meta": {'domain_of': ['RawHypothesis']} })

    @field_validator('status')
    def pattern_status(cls, v):
        pattern=re.compile(r"^(tested|untestable)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid status format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid status format: {v}"
            raise ValueError(err_msg)
        return v


class CorrectionRun(ConfiguredBaseModel):
    """
    In-memory replayable execution with all input and result identities. Large callers should store inventories and result tables as separate hashed artifacts rather than copying them into a compact provenance summary. Verification requires external trusted family and raw inputs.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'Correction run',
         'tree_root': True})

    policy: MultipleTestingPolicy = Field(default=..., title="Policy", description="""Immutable supported policy executed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'PostArmMultiplicityRun']} })
    family: HypothesisFamily = Field(default=..., title="Family", description="""Complete canonical declaration executed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun']} })
    hypotheses: list[RawHypothesis] = Field(default=..., title="Hypotheses", description="""Complete raw result inventory, canonicalized by hypothesis identity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun']} })
    results: list[AdjustedHypothesis] = Field(default=..., title="Results", description="""Complete adjusted result inventory, canonicalized by hypothesis identity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun']} })
    recorded_at: datetime  = Field(default=..., title="Recorded at", description="""UTC execution timestamp, timezone required by the runtime verifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun',
                       'ExternalFamilyCorrection',
                       'PostArmMultiplicityRun']} })
    engine_version: str = Field(default=..., title="Engine version", description="""Version of the correction engine contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun']} })
    scipy_version: str = Field(default=..., title="Scipy version", description="""Version of the established numerical implementation used.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'PostArmMultiplicityRun']} })
    numpy_version: str = Field(default=..., title="Numpy version", description="""Version of the array implementation used by SciPy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'PostArmMultiplicityRun']} })
    code_sha256: str = Field(default=..., title="Code sha256", description="""Digest of engine and generated model source bytes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'PostArmMultiplicityRun']} })
    policy_sha256: str = Field(default=..., title="Policy sha256", description="""Canonical executed policy digest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'ExternalFamilyCorrection']} })
    family_sha256: str = Field(default=..., title="Family sha256", description="""Canonical declared inventory and scientific metadata digest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'ExternalFamilyCorrection']} })
    input_sha256: str = Field(default=..., title="Input sha256", description="""Canonical raw hypothesis digest; includes test status and undefined reasons.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'ExternalFamilyCorrection']} })
    output_sha256: str = Field(default=..., title="Output sha256", description="""Canonical adjusted hypothesis digest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'ExternalFamilyCorrection']} })
    run_sha256: str = Field(default=..., title="Run sha256", description="""Execution envelope digest including scientific hashes, software and timestamp.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun']} })

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

    @field_validator('policy_sha256')
    def pattern_policy_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid policy_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid policy_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('family_sha256')
    def pattern_family_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid family_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid family_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('input_sha256')
    def pattern_input_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid input_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid input_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('output_sha256')
    def pattern_output_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid output_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid output_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('run_sha256')
    def pattern_run_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid run_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid run_sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class ExternalFamilyCorrection(ConfiguredBaseModel):
    """
    A complete correction with inventories stored once as canonical JSON Lines artifacts.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'External family correction'})

    metadata: HypothesisFamilyMetadata = Field(default=..., title="Metadata", description="""Scientific family metadata independent of observed p-values.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ExternalFamilyCorrection']} })
    declared: SourceArtifact = Field(default=..., title="Declared", description="""Sorted complete HypothesisDefinition inventory with one record per line.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ExternalFamilyCorrection']} })
    raw: SourceArtifact = Field(default=..., title="Raw", description="""Sorted complete RawHypothesis inventory with one record per line.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ExternalFamilyCorrection']} })
    corrected: SourceArtifact = Field(default=..., title="Corrected", description="""Sorted complete AdjustedHypothesis inventory with one record per line.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ExternalFamilyCorrection']} })
    recorded_at: datetime  = Field(default=..., title="Recorded at", description="""UTC correction timestamp.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun',
                       'ExternalFamilyCorrection',
                       'PostArmMultiplicityRun']} })
    policy_sha256: str = Field(default=..., title="Policy sha256", description="""Canonical policy identity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'ExternalFamilyCorrection']} })
    family_sha256: str = Field(default=..., title="Family sha256", description="""Canonical full family identity, equal to the in-memory contract digest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'ExternalFamilyCorrection']} })
    input_sha256: str = Field(default=..., title="Input sha256", description="""Canonical raw hypothesis identity, equal to the in-memory contract digest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'ExternalFamilyCorrection']} })
    output_sha256: str = Field(default=..., title="Output sha256", description="""Canonical adjusted hypothesis identity, equal to the in-memory contract digest.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'ExternalFamilyCorrection']} })
    tested: int = Field(default=..., title="Tested", description="""Number of members with defined tests.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ExternalFamilyCorrection']} })
    untestable: int = Field(default=..., title="Untestable", description="""Number of declared members reserving internal p one.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ExternalFamilyCorrection']} })

    @field_validator('policy_sha256')
    def pattern_policy_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid policy_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid policy_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('family_sha256')
    def pattern_family_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid family_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid family_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('input_sha256')
    def pattern_input_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid input_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid input_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('output_sha256')
    def pattern_output_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid output_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid output_sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class PostArmMultiplicityRun(ConfiguredBaseModel):
    """
    Verified additive stage consuming one shared core and completed three-arm comparison.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/multiple-testing',
         'title': 'Post arm multiplicity run'})

    version: str = Field(default=..., title="Version", description="""Supported sidecar execution format.""", json_schema_extra = { "linkml_meta": {'domain_of': ['MultipleTestingPolicy', 'PostArmMultiplicityRun']} })
    recorded_at: datetime  = Field(default=..., title="Recorded at", description="""Time the additive stage executed, in UTC.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun',
                       'ExternalFamilyCorrection',
                       'PostArmMultiplicityRun']} })
    comparison_manifest_path: str = Field(default=..., title="Comparison manifest path", description="""Actual verified upstream comparison record.""", json_schema_extra = { "linkml_meta": {'domain_of': ['PostArmMultiplicityRun']} })
    comparison_manifest_sha256: str = Field(default=..., title="Comparison manifest sha256", description="""Exact consumed comparison record bytes.""", json_schema_extra = { "linkml_meta": {'domain_of': ['PostArmMultiplicityRun']} })
    core_run_id: str = Field(default=..., title="Core run id", description="""Shared verified core execution identity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['PostArmMultiplicityRun']} })
    contrast_population: str = Field(default=..., title="Contrast population", description="""Shared verified population.""", json_schema_extra = { "linkml_meta": {'domain_of': ['PostArmMultiplicityRun']} })
    direction_basis: str = Field(default=..., title="Direction basis", description="""Shared verified orientation; never reapplied by this stage.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata', 'PostArmMultiplicityRun']} })
    strain_rollup: str = Field(default=..., title="Strain rollup", description="""Verified source core strain scope.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata', 'PostArmMultiplicityRun']} })
    policy: MultipleTestingPolicy = Field(default=..., title="Policy", description="""Exact supported additive policy.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'PostArmMultiplicityRun']} })
    code_sha256: str = Field(default=..., title="Code sha256", description="""Correction, adapter and generated model source identity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'PostArmMultiplicityRun']} })
    scipy_version: str = Field(default=..., title="Scipy version", description="""Numerical implementation version used by the stage.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'PostArmMultiplicityRun']} })
    numpy_version: str = Field(default=..., title="Numpy version", description="""Array implementation version used by the stage.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CorrectionRun', 'PostArmMultiplicityRun']} })
    sources: list[SourceArtifact] = Field(default=..., title="Sources", description="""Consumed source artifact identities and hashes, independently checked during replay.""", json_schema_extra = { "linkml_meta": {'domain_of': ['HypothesisFamilyMetadata', 'PostArmMultiplicityRun']} })
    families: list[ExternalFamilyCorrection] = Field(default=..., title="Families", description="""Complete external corrections; each family remains a distinct scientific unit.""", json_schema_extra = { "linkml_meta": {'domain_of': ['PostArmMultiplicityRun']} })
    reports: list[SourceArtifact] = Field(default=..., title="Reports", description="""Hashed migration and interpretation reports.""", json_schema_extra = { "linkml_meta": {'domain_of': ['PostArmMultiplicityRun']} })

    @field_validator('version')
    def pattern_version(cls, v):
        pattern=re.compile(r"^1\.1\.0$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid version format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid version format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('comparison_manifest_sha256')
    def pattern_comparison_manifest_sha256(cls, v):
        pattern=re.compile(r"^[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid comparison_manifest_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid comparison_manifest_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('strain_rollup')
    def pattern_strain_rollup(cls, v):
        pattern=re.compile(r"^(off|phenotype)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid strain_rollup format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid strain_rollup format: {v}"
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


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
MultipleTestingPolicy.model_rebuild()
SourceArtifact.model_rebuild()
HypothesisCoordinates.model_rebuild()
HypothesisDefinition.model_rebuild()
HypothesisFamilyMetadata.model_rebuild()
HypothesisFamily.model_rebuild()
RawHypothesis.model_rebuild()
AdjustedHypothesis.model_rebuild()
CorrectionRun.model_rebuild()
ExternalFamilyCorrection.model_rebuild()
PostArmMultiplicityRun.model_rebuild()
