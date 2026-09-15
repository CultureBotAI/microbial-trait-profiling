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
     'prefixes': {'CHEBI': {'prefix_prefix': 'CHEBI',
                            'prefix_reference': 'http://purl.obolibrary.org/obo/CHEBI_'},
                  'EFO': {'prefix_prefix': 'EFO',
                          'prefix_reference': 'http://identifiers.org/efo/'},
                  'ENVO': {'prefix_prefix': 'ENVO',
                           'prefix_reference': 'http://purl.obolibrary.org/obo/ENVO_'},
                  'METPO': {'prefix_prefix': 'METPO',
                            'prefix_reference': 'http://purl.obolibrary.org/obo/METPO_'},
                  'NCBITaxon': {'prefix_prefix': 'NCBITaxon',
                                'prefix_reference': 'http://purl.obolibrary.org/obo/NCBITaxon_'},
                  'UBERON': {'prefix_prefix': 'UBERON',
                             'prefix_reference': 'http://purl.obolibrary.org/obo/UBERON_'},
                  'biolink': {'prefix_prefix': 'biolink',
                              'prefix_reference': 'https://w3id.org/biolink/vocab/'},
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


class BiolinkCategoryEnum(str, Enum):
    """
    Common Biolink categories used by KG-Microbe nodes.
    """
    biolinkCOLONOrganismTaxon = "biolink:OrganismTaxon"
    """
    A taxonomic grouping of organisms.
    """
    biolinkCOLONPhenotype = "biolink:Phenotype"
    """
    An observable characteristic or trait.
    """
    biolinkCOLONPhenotypicQuality = "biolink:PhenotypicQuality"
    """
    A phenotype quality term.
    """
    biolinkCOLONOntologyClass = "biolink:OntologyClass"
    """
    A general ontology class.
    """
    biolinkCOLONChemicalEntity = "biolink:ChemicalEntity"
    """
    A chemical substance.
    """
    biolinkCOLONSmallMolecule = "biolink:SmallMolecule"
    """
    A small molecule.
    """
    biolinkCOLONEnvironment = "biolink:Environment"
    """
    An environment or habitat.
    """
    biolinkCOLONBiologicalProcess = "biolink:BiologicalProcess"
    """
    A biological process.
    """
    biolinkCOLONMolecularActivity = "biolink:MolecularActivity"
    """
    A molecular function or activity.
    """
    biolinkCOLONMacromolecularComplex = "biolink:MacromolecularComplex"
    """
    A macromolecular complex.
    """
    biolinkCOLONAttribute = "biolink:Attribute"
    """
    An attribute or quality.
    """
    biolinkCOLONNamedThing = "biolink:NamedThing"
    """
    A generic named entity.
    """


class BiolinkPredicateEnum(str, Enum):
    """
    Common relationship predicates used by KG-Microbe edges.
    """
    biolinkCOLONsubclass_of = "biolink:subclass_of"
    """
    Taxonomic hierarchy relationship.
    """
    biolinkCOLONhas_phenotype = "biolink:has_phenotype"
    """
    Taxon has a phenotype or trait.
    """
    biolinkCOLONcapable_of = "biolink:capable_of"
    """
    Taxon is capable of an activity.
    """
    biolinkCOLONinteracts_with = "biolink:interacts_with"
    """
    Taxon interacts with an entity.
    """
    biolinkCOLONproduces = "biolink:produces"
    """
    Taxon produces a substance.
    """
    biolinkCOLONlocation_of = "biolink:location_of"
    """
    Environment is location of taxon.
    """
    biolinkCOLONlocated_in = "biolink:located_in"
    """
    Taxon is located in an environment.
    """
    biolinkCOLONaffects = "biolink:affects"
    """
    Taxon affects another entity.
    """
    biolinkCOLONparticipates_in = "biolink:participates_in"
    """
    Taxon participates in a process.
    """



class NamedThing(ConfiguredBaseModel):
    """
    An entity with a stable identifier. All identifiable classes in this schema inherit from NamedThing.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True,
         'from_schema': 'https://w3id.org/bsdbng/schema',
         'title': 'Named Thing'})

    id: str = Field(default=..., title="ID", description="""Stable identifier for this entity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node']} })


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
    id: str = Field(default=..., title="ID", description="""BugSigDB study identifier. Usually numeric (e.g. bsdb:10202341) but some entries use DOIs or PMC IDs (e.g. bsdb:PMC11017998).""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node']} })

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
    id: str = Field(default=..., title="ID", description="""BugSigDB experiment identifier, e.g. bsdb:10202341-1.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node']} })

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
    id: str = Field(default=..., title="ID", description="""BugSigDB signature identifier, e.g. bsdb:10202341-1-1.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node']} })

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

    taxon_name: str = Field(default=..., title="Taxon Name", description="""Reported taxon label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon', 'TaxonTraitMatrixRow']} })
    taxonomic_rank: Optional[TaxonomicRankEnum] = Field(default=None, title="Taxonomic Rank", description="""Taxonomic rank of the reported taxon.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon', 'TaxonTraitMatrixRow'], 'recommended': True} })
    id: str = Field(default=..., title="ID", description="""NCBITaxon CURIE, e.g. NCBITaxon:9606.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node'], 'recommended': True} })

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


class Node(ConfiguredBaseModel):
    """
    A node in the KG-Microbe knowledge graph. Nodes can represent taxa, phenotypes, traits, chemicals, environments, or other entities.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/bsdbng/schema',
         'slot_usage': {'category': {'name': 'category', 'required': True},
                        'id': {'identifier': True, 'name': 'id', 'required': True}},
         'title': 'KG-Microbe Node'})

    id: str = Field(default=..., title="ID", description="""Stable identifier for this entity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node']} })
    category: str = Field(default=..., title="Category", description="""Biolink category of a KG-Microbe node.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node']} })
    name: Optional[str] = Field(default=None, title="Name", description="""Human-readable node label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node']} })
    xref: Optional[list[str]] = Field(default=None, title="Cross Reference", description="""Cross-reference to another database.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node']} })
    provided_by: Optional[list[str]] = Field(default=None, title="Provided By", description="""Data source or provenance for a node or edge.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node', 'Edge']} })
    synonym: Optional[list[str]] = Field(default=None, title="Synonym", description="""Alternative node label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node']} })


class Edge(ConfiguredBaseModel):
    """
    An edge or relationship in the KG-Microbe knowledge graph connecting two nodes.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/bsdbng/schema',
         'slot_usage': {'object': {'name': 'object', 'required': True},
                        'predicate': {'name': 'predicate', 'required': True},
                        'subject': {'name': 'subject', 'required': True}},
         'title': 'KG-Microbe Edge'})

    edge_id: Optional[str] = Field(default=None, title="Edge ID", description="""Stable identifier for a KG-Microbe edge.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Edge']} })
    subject: str = Field(default=..., title="Subject", description="""Subject node of a KG-Microbe relationship.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Edge']} })
    predicate: str = Field(default=..., title="Predicate", description="""Relationship predicate.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Edge']} })
    object: str = Field(default=..., title="Object", description="""Object node of a KG-Microbe relationship.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Edge']} })
    relation: Optional[str] = Field(default=None, title="Relation", description="""Ontology relation term supporting the predicate.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Edge']} })
    provided_by: Optional[list[str]] = Field(default=None, title="Provided By", description="""Data source or provenance for a node or edge.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node', 'Edge']} })


class Phenotype(ConfiguredBaseModel):
    """
    A microbial phenotype or trait from the microbial phenotypes ontology associated with a taxon in KG-Microbe.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/bsdbng/schema',
         'slot_usage': {'phenotype_id': {'name': 'phenotype_id', 'required': True}},
         'title': 'KG-Microbe Phenotype'})

    phenotype_id: str = Field(default=..., title="Phenotype ID", description="""METPO phenotype or trait identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype']} })
    phenotype_name: Optional[str] = Field(default=None, title="Phenotype Name", description="""Human-readable phenotype or trait label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype']} })
    count: Optional[int] = Field(default=None, title="Count", description="""Number of taxa at or below a row taxon that have the trait assertion.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype', 'KgTraitAssertion']} })
    from_descendants: Optional[bool] = Field(default=None, title="From Descendants", description="""Whether the trait assertion count includes descendant taxa.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype', 'KgTraitAssertion']} })

    @field_validator('phenotype_id')
    def pattern_phenotype_id(cls, v):
        pattern=re.compile(r"^METPO:[0-9]+$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid phenotype_id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid phenotype_id format: {v}"
            raise ValueError(err_msg)
        return v


class KgTraitAssertion(ConfiguredBaseModel):
    """
    Normalized representation of one dynamic trait column in the taxa-by-traits matrix. Direct phenotype columns use a METPO object such as METPO:1000615. Relation-value columns combine a METPO predicate with a specific KG object, such as METPO:2000006 plus a CHEBI carbon source or METPO:2000517 plus a MediaDive growth medium.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/bsdbng/schema',
         'slot_usage': {'count': {'name': 'count', 'required': True},
                        'object_id': {'name': 'object_id', 'required': True},
                        'trait_assertion_id': {'name': 'trait_assertion_id',
                                               'required': True}},
         'title': 'KG-Microbe Trait Assertion'})

    trait_assertion_id: str = Field(default=..., title="Trait Assertion ID", description="""Stable key for a generated trait matrix column. Direct METPO phenotype columns use the METPO object CURIE; relation-value columns combine the predicate and object identifiers.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion']} })
    trait_label: Optional[str] = Field(default=None, title="Trait Label", description="""Human-readable label for a generated trait matrix column.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion']} })
    predicate_id: Optional[str] = Field(default=None, title="Predicate ID", description="""Predicate used by the KG-Microbe assertion. Direct phenotype assertions may use biolink:has_phenotype, while relation-value assertions use METPO predicates such as METPO:2000006, METPO:2000103, or METPO:2000517.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion']} })
    predicate_name: Optional[str] = Field(default=None, title="Predicate Name", description="""Human-readable predicate label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion']} })
    object_id: str = Field(default=..., title="Object ID", description="""KG object identifier for a trait assertion. This may be a METPO phenotype, CHEBI chemical, GO process, pathway identifier, MediaDive medium, or other KG-Microbe object.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion']} })
    object_name: Optional[str] = Field(default=None, title="Object Name", description="""Human-readable object label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion']} })
    object_category: Optional[str] = Field(default=None, title="Object Category", description="""Biolink category of the KG object.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion']} })
    count: int = Field(default=..., title="Count", description="""Number of taxa at or below a row taxon that have the trait assertion.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype', 'KgTraitAssertion']} })
    from_descendants: Optional[bool] = Field(default=None, title="From Descendants", description="""Whether the trait assertion count includes descendant taxa.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype', 'KgTraitAssertion']} })


class TaxonTraitMatrixRow(ConfiguredBaseModel):
    """
    One BugSigDB-reported taxon row in a taxa-by-traits matrix. Trait columns in the generated TSV are dynamic KG-Microbe phenotype or trait assertion labels; the trait_counts slot provides the normalized representation of the same counts.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/bsdbng/schema',
         'slot_usage': {'taxon_id': {'name': 'taxon_id', 'required': True},
                        'taxon_name': {'name': 'taxon_name', 'required': True}},
         'title': 'Taxon Trait Matrix Row'})

    taxon_id: str = Field(default=..., title="Taxon ID", description="""NCBITaxon identifier for a matrix row.""", json_schema_extra = { "linkml_meta": {'domain_of': ['TaxonTraitMatrixRow']} })
    taxon_name: str = Field(default=..., title="Taxon Name", description="""Reported taxon label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon', 'TaxonTraitMatrixRow']} })
    taxonomic_rank: Optional[TaxonomicRankEnum] = Field(default=None, title="Taxonomic Rank", description="""Taxonomic rank of the reported taxon.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon', 'TaxonTraitMatrixRow']} })
    trait_counts: Optional[list[KgTraitAssertion]] = Field(default=None, title="Trait Counts", description="""Normalized KG-Microbe trait assertion counts for a matrix row.""", json_schema_extra = { "linkml_meta": {'domain_of': ['TaxonTraitMatrixRow']} })

    @field_validator('taxon_id')
    def pattern_taxon_id(cls, v):
        pattern=re.compile(r"^NCBITaxon:[0-9]+$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid taxon_id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid taxon_id format: {v}"
            raise ValueError(err_msg)
        return v


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
NamedThing.model_rebuild()
Study.model_rebuild()
Experiment.model_rebuild()
Signature.model_rebuild()
Taxon.model_rebuild()
Node.model_rebuild()
Edge.model_rebuild()
Phenotype.model_rebuild()
KgTraitAssertion.model_rebuild()
TaxonTraitMatrixRow.model_rebuild()
