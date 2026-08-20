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


metamodel_version = "1.11.0"
version = "None"


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


linkml_meta = LinkMLMeta({'default_prefix': 'bsdbng',
     'default_range': 'string',
     'description': 'LinkML schema for normalized BugSigDB studies, experiments, '
                    'signatures, and taxa.',
     'id': 'https://w3id.org/bsdbng/schema',
     'imports': ['linkml:types'],
     'license': 'MIT',
     'name': 'bsdbng',
     'prefixes': {'EFO': {'prefix_prefix': 'EFO',
                          'prefix_reference': 'http://identifiers.org/efo/'},
                  'NCBITaxon': {'prefix_prefix': 'NCBITaxon',
                                'prefix_reference': 'http://purl.obolibrary.org/obo/NCBITaxon_'},
                  'UBERON': {'prefix_prefix': 'UBERON',
                             'prefix_reference': 'http://purl.obolibrary.org/obo/UBERON_'},
                  'bsdb': {'prefix_prefix': 'bsdb',
                           'prefix_reference': 'https://bugsigdb.org/Study_'},
                  'bsdbng': {'prefix_prefix': 'bsdbng',
                             'prefix_reference': 'https://w3id.org/bsdbng/'},
                  'linkml': {'prefix_prefix': 'linkml',
                             'prefix_reference': 'https://w3id.org/linkml/'}},
     'source_file': 'schema/bsdbng.yaml',
     'title': 'BugSigDB Normalized Graph Schema'} )

class DirectionEnum(str, Enum):
    """
    Allowed values for reported abundance direction.
    """
    increased = "increased"
    decreased = "decreased"


class DiversityDirectionEnum(str, Enum):
    """
    Reported direction of change in an alpha diversity metric between the two sample groups.
    """
    increased = "increased"
    decreased = "decreased"
    unchanged = "unchanged"


class TaxonomicRankEnum(str, Enum):
    """
    Allowed taxonomic rank values for normalized BugSigDB taxa.
    """
    strain = "strain"
    species = "species"
    subspecies = "subspecies"
    isolate = "isolate"
    genus = "genus"
    family = "family"
    order = "order"
    class_ = "class"
    phylum = "phylum"
    superkingdom = "superkingdom"



class NamedThing(ConfiguredBaseModel):
    """
    An entity with a stable identifier. All identifiable classes in this schema inherit from NamedThing.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True,
         'from_schema': 'https://w3id.org/bsdbng/schema',
         'title': 'Named Thing'})

    id: str = Field(default=..., title="ID", description="""Stable identifier for this entity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing']} })


class Study(NamedThing):
    """
    One published microbiome study from which BugSigDB collects and standardizes microbial signatures.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/bsdbng/schema',
         'id_prefixes': ['bsdb'],
         'slot_usage': {'id': {'description': 'BugSigDB study identifier. Usually '
                                              'numeric (e.g. bsdb:10202341) but some '
                                              'entries use DOIs or PMC IDs (e.g. '
                                              'bsdb:PMC11017998).',
                               'name': 'id',
                               'pattern': '^bsdb:.+$'}},
         'title': 'Study',
         'tree_root': True})

    source_record_id: str = Field(default=..., title="Source Record ID", description="""Study identifier exactly as exported by BugSigDB.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    pmid: Optional[int] = Field(default=None, title="PMID", description="""PubMed identifier for the study.""", ge=1, json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    title: Optional[str] = Field(default=None, title="Title", description="""Study title.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    publication_year: Optional[int] = Field(default=None, title="Publication Year", description="""Publication year.""", ge=1900, le=2100, json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    doi: Optional[str] = Field(default=None, title="DOI", description="""Digital object identifier for the study.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    url: Optional[str] = Field(default=None, title="URL", description="""Stable URL for the study when no DOI or PMID is available.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    study_design: Optional[str] = Field(default=None, title="Study Design", description="""Study design type (e.g. case-control, cross-sectional, meta-analysis).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    authors_list: Optional[str] = Field(default=None, title="Authors List", description="""Author list as provided by BugSigDB.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    journal: Optional[str] = Field(default=None, title="Journal", description="""Journal of publication.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    keywords: Optional[str] = Field(default=None, title="Keywords", description="""Keywords associated with the study.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    experiments: list[Experiment] = Field(default=..., title="Experiments", description="""Semantic experiment units recorded for the study.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    id: str = Field(default=..., title="ID", description="""BugSigDB study identifier. Usually numeric (e.g. bsdb:10202341) but some entries use DOIs or PMC IDs (e.g. bsdb:PMC11017998).""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing']} })

    @field_validator('doi')
    def pattern_doi(cls, v):
        pattern=re.compile(r"^10\..+$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid doi format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid doi format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('id')
    def pattern_id(cls, v):
        pattern=re.compile(r"^bsdb:.+$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid id format: {v}"
            raise ValueError(err_msg)
        return v


class Experiment(NamedThing):
    """
    One semantic unit within a study that records contrasted sample groups and key metadata relevant to differential abundance findings.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/bsdbng/schema',
         'id_prefixes': ['bsdb'],
         'slot_usage': {'id': {'description': 'BugSigDB experiment identifier, e.g. '
                                              'bsdb:10202341-1.',
                               'name': 'id',
                               'pattern': '^bsdb:.+-[0-9]+$'}},
         'title': 'Experiment'})

    experiment_name: Optional[str] = Field(default=None, title="Experiment Name", description="""Label for the experiment.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    location_of_subjects: Optional[str] = Field(default=None, title="Location of Subjects", description="""Geographic location of study subjects.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    host_species: Optional[str] = Field(default=None, title="Host Species", description="""Host organism species (e.g. Homo sapiens, Mus musculus).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    body_site: Optional[str] = Field(default=None, title="Body Site", description="""Anatomical site of sample collection (e.g. stool, skin, oral cavity).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    body_site_ontology_id: Optional[list[str]] = Field(default=None, title="Body Site Ontology ID", description="""UBERON ontology term(s) for the body site (e.g. UBERON:0001988).""", json_schema_extra = { "linkml_meta": {'aliases': ['UBERON ID'], 'domain_of': ['Experiment']} })
    condition: Optional[str] = Field(default=None, title="Condition", description="""Disease or condition under study (e.g. obesity, colorectal cancer).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    condition_ontology_id: Optional[list[str]] = Field(default=None, title="Condition Ontology ID", description="""Ontology term(s) for the condition. BugSigDB labels this column \"EFO ID\" but it contains CURIEs from multiple ontologies including EFO, MONDO, HP, CHEBI, GO, and others.""", json_schema_extra = { "linkml_meta": {'aliases': ['EFO ID'], 'domain_of': ['Experiment']} })
    group_0_name: Optional[str] = Field(default=None, title="Group 0 Name", description="""Label for the control or reference sample group.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    group_1_name: Optional[str] = Field(default=None, title="Group 1 Name", description="""Label for the case or experimental sample group.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    group_1_definition: Optional[str] = Field(default=None, title="Group 1 Definition", description="""Detailed definition of the case group.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    group_0_sample_size: Optional[int] = Field(default=None, title="Group 0 Sample Size", description="""Number of subjects in the control group.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    group_1_sample_size: Optional[int] = Field(default=None, title="Group 1 Sample Size", description="""Number of subjects in the case group.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    antibiotics_exclusion: Optional[str] = Field(default=None, title="Antibiotics Exclusion", description="""Whether subjects on antibiotics were excluded.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    sequencing_type: Optional[str] = Field(default=None, title="Sequencing Type", description="""Sequencing approach (e.g. 16S, shotgun, metatranscriptomics).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    variable_region_16s: Optional[str] = Field(default=None, title="16S Variable Region", description="""Targeted 16S rRNA variable region (e.g. V4, V3V4).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    sequencing_platform: Optional[str] = Field(default=None, title="Sequencing Platform", description="""Sequencing platform used (e.g. Illumina, Ion Torrent).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    data_transformation: Optional[str] = Field(default=None, title="Data Transformation", description="""Data normalization or transformation applied (e.g. relative abundances).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    statistical_test: Optional[str] = Field(default=None, title="Statistical Test", description="""Statistical test used for differential abundance (e.g. Mann-Whitney, DESeq2).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    significance_threshold: Optional[float] = Field(default=None, title="Significance Threshold", description="""P-value or FDR threshold used for significance.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    mht_correction: Optional[str] = Field(default=None, title="MHT Correction", description="""Multiple hypothesis testing correction method (e.g. BH, Bonferroni).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    lda_score_above: Optional[float] = Field(default=None, title="LDA Score Above", description="""LEfSe LDA score threshold if applicable.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    matched_on: Optional[str] = Field(default=None, title="Matched On", description="""Variables on which groups were matched (e.g. age, BMI).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    confounders_controlled_for: Optional[str] = Field(default=None, title="Confounders Controlled For", description="""Confounders accounted for in the analysis.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    pielou: Optional[DiversityDirectionEnum] = Field(default=None, title="Pielou Evenness", description="""Direction of change in Pielou's evenness index between groups.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    shannon: Optional[DiversityDirectionEnum] = Field(default=None, title="Shannon Diversity", description="""Direction of change in Shannon diversity index between groups.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    chao1: Optional[DiversityDirectionEnum] = Field(default=None, title="Chao1 Richness", description="""Direction of change in Chao1 richness estimator between groups.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    simpson: Optional[DiversityDirectionEnum] = Field(default=None, title="Simpson Index", description="""Direction of change in Simpson diversity index between groups.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    inverse_simpson: Optional[DiversityDirectionEnum] = Field(default=None, title="Inverse Simpson", description="""Direction of change in inverse Simpson index between groups.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    richness: Optional[DiversityDirectionEnum] = Field(default=None, title="Richness", description="""Direction of change in observed species richness between groups.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    signatures: list[Signature] = Field(default=..., title="Signatures", description="""Signatures associated with the experiment.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    id: str = Field(default=..., title="ID", description="""BugSigDB experiment identifier, e.g. bsdb:10202341-1.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing']} })

    @field_validator('body_site_ontology_id')
    def pattern_body_site_ontology_id(cls, v):
        pattern=re.compile(r"^UBERON:[0-9]+$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid body_site_ontology_id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid body_site_ontology_id format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('id')
    def pattern_id(cls, v):
        pattern=re.compile(r"^bsdb:.+-[0-9]+$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid id format: {v}"
            raise ValueError(err_msg)
        return v


class Signature(NamedThing):
    """
    A microbial signature represented as an unordered set of taxa sharing a common property or response to a study condition.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/bsdbng/schema',
         'id_prefixes': ['bsdb'],
         'slot_usage': {'direction': {'name': 'direction', 'required': True},
                        'id': {'description': 'BugSigDB signature identifier, e.g. '
                                              'bsdb:10202341-1-1.',
                               'name': 'id',
                               'pattern': '^bsdb:.+-[0-9]+-[0-9]+$'}},
         'title': 'Signature'})

    direction: DirectionEnum = Field(default=..., title="Direction", description="""Reported direction of abundance change for taxa in the signature.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Signature']} })
    signature_source: Optional[str] = Field(default=None, title="Source", description="""Source of the signature within the publication (e.g. Table S6, Figure 3).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Signature']} })
    signature_description: Optional[str] = Field(default=None, title="Description", description="""Free-text description of the signature from BugSigDB curators.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Signature']} })
    taxa: list[Taxon] = Field(default=..., title="Taxa", description="""Taxa listed in the signature.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Signature']} })
    id: str = Field(default=..., title="ID", description="""BugSigDB signature identifier, e.g. bsdb:10202341-1-1.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing']} })

    @field_validator('id')
    def pattern_id(cls, v):
        pattern=re.compile(r"^bsdb:.+-[0-9]+-[0-9]+$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid id format: {v}"
            raise ValueError(err_msg)
        return v


class Taxon(NamedThing):
    """
    A taxonomic unit of any rank designating a microbial organism or group of microbial organisms.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/bsdbng/schema',
         'id_prefixes': ['NCBITaxon'],
         'slot_usage': {'id': {'description': 'NCBITaxon CURIE, e.g. NCBITaxon:9606.',
                               'name': 'id',
                               'pattern': '^NCBITaxon:[0-9]+$',
                               'recommended': True,
                               'required': False},
                        'taxon_name': {'name': 'taxon_name', 'required': True},
                        'taxonomic_rank': {'name': 'taxonomic_rank',
                                           'recommended': True,
                                           'required': False}},
         'title': 'Taxon'})

    taxon_name: str = Field(default=..., title="Taxon Name", description="""Reported taxon label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon']} })
    taxonomic_rank: Optional[TaxonomicRankEnum] = Field(default=None, title="Taxonomic Rank", description="""Taxonomic rank of the reported taxon.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon'], 'recommended': True} })
    id: str = Field(default=..., title="ID", description="""NCBITaxon CURIE, e.g. NCBITaxon:9606.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing'], 'recommended': True} })

    @field_validator('id')
    def pattern_id(cls, v):
        pattern=re.compile(r"^NCBITaxon:[0-9]+$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid id format: {v}"
            raise ValueError(err_msg)
        return v


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
NamedThing.model_rebuild()
Study.model_rebuild()
Experiment.model_rebuild()
Signature.model_rebuild()
Taxon.model_rebuild()
