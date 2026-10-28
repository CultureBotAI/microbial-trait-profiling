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


linkml_meta = LinkMLMeta({'default_prefix': 'mtpc',
     'default_range': 'string',
     'description': 'Scientific payload identities and coverage accounting for '
                    'core pipeline products. Canonical study and KG assertion '
                    'types are reused from mtp.',
     'id': 'https://w3id.org/mtp/core-products',
     'imports': ['linkml:types', 'mtp'],
     'license': 'MIT',
     'name': 'mtp_core_products',
     'prefixes': {'linkml': {'prefix_prefix': 'linkml',
                             'prefix_reference': 'https://w3id.org/linkml/'},
                  'mtpc': {'prefix_prefix': 'mtpc',
                           'prefix_reference': 'https://w3id.org/mtp/core-products/'}},
     'source_file': 'schema/core_products.yaml',
     'title': 'Core microbial trait product contracts'} )

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


class FeatureClass(str, Enum):
    """
    Permitted feature class values.
    """
    trait_annotation = "trait_annotation"
    culture_medium = "culture_medium"


class CoverageStatus(str, Enum):
    """
    Permitted coverage status values.
    """
    no_reported_taxa = "no_reported_taxa"
    no_lookup_matches = "no_lookup_matches"
    partial_lookup_coverage = "partial_lookup_coverage"
    matched_without_annotations = "matched_without_annotations"
    partial_annotation_coverage = "partial_annotation_coverage"
    complete_annotation_coverage = "complete_annotation_coverage"



class NamedThing(ConfiguredBaseModel):
    """
    An entity with a stable identifier. All identifiable classes in this schema inherit from NamedThing.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'abstract': True,
         'from_schema': 'https://w3id.org/mtp/schema',
         'title': 'Named Thing'})

    id: str = Field(default=..., title="ID", description="""Stable identifier for this entity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node', 'CoreProductCatalog']} })


class Study(NamedThing):
    """
    One published microbiome study from which BugSigDB collects and standardizes microbial signatures.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/schema',
         'id_prefixes': ['bsdb'],
         'slot_usage': {'id': {'description': 'BugSigDB study identifier. Usually '
                                              'numeric (e.g. bsdb:10202341) but some '
                                              'entries use DOIs or PMC IDs (e.g. '
                                              'bsdb:PMC11017998).',
                               'name': 'id',
                               'pattern': '^bsdb:.+$'}},
         'title': 'Study',
         'tree_root': True})

    source_record_id: str = Field(default=..., title="Source Record ID", description="""Study identifier exactly as exported by BugSigDB.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'StudyTableRow']} })
    pmid: Optional[int] = Field(default=None, title="PMID", description="""PubMed identifier for the study.""", ge=1, json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'StudyTableRow']} })
    title: Optional[str] = Field(default=None, title="Title", description="""Study title.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'StudyTableRow']} })
    publication_year: Optional[int] = Field(default=None, title="Publication Year", description="""Publication year.""", ge=1900, le=2100, json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'StudyTableRow']} })
    doi: Optional[str] = Field(default=None, title="DOI", description="""Digital object identifier for the study.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'VerifiedRelease', 'StudyTableRow']} })
    url: Optional[str] = Field(default=None, title="URL", description="""Stable URL for the study when no DOI or PMID is available.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'VerifiedRelease', 'StudyTableRow']} })
    study_design: Optional[str] = Field(default=None, title="Study Design", description="""Study design type (e.g. case-control, cross-sectional, meta-analysis).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    authors_list: Optional[str] = Field(default=None, title="Authors List", description="""Author list as provided by BugSigDB.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    journal: Optional[str] = Field(default=None, title="Journal", description="""Journal of publication.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    keywords: Optional[str] = Field(default=None, title="Keywords", description="""Keywords associated with the study.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    experiments: list[Experiment] = Field(default=..., title="Experiments", description="""Semantic experiment units recorded for the study.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study']} })
    id: str = Field(default=..., title="ID", description="""BugSigDB study identifier. Usually numeric (e.g. bsdb:10202341) but some entries use DOIs or PMC IDs (e.g. bsdb:PMC11017998).""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node', 'CoreProductCatalog']} })

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
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/schema',
         'id_prefixes': ['bsdb'],
         'slot_usage': {'id': {'description': 'BugSigDB experiment identifier, e.g. '
                                              'bsdb:10202341-1.',
                               'name': 'id',
                               'pattern': '^bsdb:.+-[0-9]+$'}},
         'title': 'Experiment'})

    experiment_name: Optional[str] = Field(default=None, title="Experiment Name", description="""Label for the experiment.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    location_of_subjects: Optional[str] = Field(default=None, title="Location of Subjects", description="""Geographic location of study subjects.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    host_species: Optional[str] = Field(default=None, title="Host Species", description="""Host organism species (e.g. Homo sapiens, Mus musculus).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    body_site: Optional[str] = Field(default=None, title="Body Site", description="""Anatomical site of sample collection (e.g. stool, skin, oral cavity).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    body_site_ontology_id: Optional[list[str]] = Field(default=None, title="Body Site Ontology ID", description="""UBERON ontology term(s) for the body site (e.g. UBERON:0001988).""", json_schema_extra = { "linkml_meta": {'aliases': ['UBERON ID'], 'domain_of': ['Experiment']} })
    condition: Optional[str] = Field(default=None, title="Condition", description="""Disease or condition under study (e.g. obesity, colorectal cancer).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    condition_ontology_id: Optional[list[str]] = Field(default=None, title="Condition Ontology ID", description="""Ontology term(s) for the condition. BugSigDB labels this column \"EFO ID\" but it contains CURIEs from multiple ontologies including EFO, MONDO, HP, CHEBI, GO, and others.""", json_schema_extra = { "linkml_meta": {'aliases': ['EFO ID'], 'domain_of': ['Experiment']} })
    group_0_name: Optional[str] = Field(default=None, title="Group 0 Name", description="""Label for the control or reference sample group.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    group_1_name: Optional[str] = Field(default=None, title="Group 1 Name", description="""Label for the case or experimental sample group.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    group_1_definition: Optional[str] = Field(default=None, title="Group 1 Definition", description="""Detailed definition of the case group.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    group_0_sample_size: Optional[int] = Field(default=None, title="Group 0 Sample Size", description="""Number of subjects in the control group.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    group_1_sample_size: Optional[int] = Field(default=None, title="Group 1 Sample Size", description="""Number of subjects in the case group.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    antibiotics_exclusion: Optional[str] = Field(default=None, title="Antibiotics Exclusion", description="""Whether subjects on antibiotics were excluded.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment']} })
    sequencing_type: Optional[str] = Field(default=None, title="Sequencing Type", description="""Sequencing approach (e.g. 16S, shotgun, metatranscriptomics).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
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
    id: str = Field(default=..., title="ID", description="""BugSigDB experiment identifier, e.g. bsdb:10202341-1.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node', 'CoreProductCatalog']} })

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
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/schema',
         'id_prefixes': ['bsdb'],
         'slot_usage': {'direction': {'name': 'direction', 'required': True},
                        'id': {'description': 'BugSigDB signature identifier, e.g. '
                                              'bsdb:10202341-1-1.',
                               'name': 'id',
                               'pattern': '^bsdb:.+-[0-9]+-[0-9]+$'}},
         'title': 'Signature'})

    direction: DirectionEnum = Field(default=..., title="Direction", description="""Reported direction of abundance change for taxa in the signature.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Signature', 'SignatureTableRow', 'TaxonOccurrenceRow']} })
    signature_source: Optional[str] = Field(default=None, title="Source", description="""Source of the signature within the publication (e.g. Table S6, Figure 3).""", json_schema_extra = { "linkml_meta": {'domain_of': ['Signature']} })
    signature_description: Optional[str] = Field(default=None, title="Description", description="""Free-text description of the signature from BugSigDB curators.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Signature']} })
    taxa: list[Taxon] = Field(default=..., title="Taxa", description="""Taxa listed in the signature.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Signature']} })
    id: str = Field(default=..., title="ID", description="""BugSigDB signature identifier, e.g. bsdb:10202341-1-1.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node', 'CoreProductCatalog']} })

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
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/schema',
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

    taxon_name: str = Field(default=..., title="Taxon Name", description="""Reported taxon label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon', 'TaxonTraitMatrixRow', 'TaxonOccurrenceRow']} })
    taxonomic_rank: Optional[TaxonomicRankEnum] = Field(default=None, title="Taxonomic Rank", description="""Taxonomic rank of the reported taxon.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon', 'TaxonTraitMatrixRow', 'TaxonOccurrenceRow'],
         'recommended': True} })
    id: str = Field(default=..., title="ID", description="""NCBITaxon CURIE, e.g. NCBITaxon:9606.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node', 'CoreProductCatalog'], 'recommended': True} })

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
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/schema',
         'slot_usage': {'category': {'name': 'category', 'required': True},
                        'id': {'identifier': True, 'name': 'id', 'required': True}},
         'title': 'KG-Microbe Node'})

    id: str = Field(default=..., title="ID", description="""Stable identifier for this entity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node', 'CoreProductCatalog']} })
    category: str = Field(default=..., title="Category", description="""Biolink category of a KG-Microbe node.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node']} })
    name: Optional[str] = Field(default=None, title="Name", description="""Human-readable node label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node']} })
    xref: Optional[list[str]] = Field(default=None, title="Cross Reference", description="""Cross-reference to another database.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node']} })
    provided_by: Optional[list[str]] = Field(default=None, title="Provided By", description="""Data source or provenance for a node or edge.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node', 'Edge']} })
    synonym: Optional[list[str]] = Field(default=None, title="Synonym", description="""Alternative node label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Node']} })


class Edge(ConfiguredBaseModel):
    """
    An edge or relationship in the KG-Microbe knowledge graph connecting two nodes.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/schema',
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
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/schema',
         'slot_usage': {'phenotype_id': {'name': 'phenotype_id', 'required': True}},
         'title': 'KG-Microbe Phenotype'})

    phenotype_id: str = Field(default=..., title="Phenotype ID", description="""METPO phenotype or trait identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype']} })
    phenotype_name: Optional[str] = Field(default=None, title="Phenotype Name", description="""Human-readable phenotype or trait label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype']} })
    count: Optional[int] = Field(default=None, title="Count", description="""Number of taxa at or below a row taxon that have the trait assertion.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype', 'KgTraitAssertion', 'FeatureCount']} })
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
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/schema',
         'slot_usage': {'count': {'name': 'count', 'required': True},
                        'object_id': {'name': 'object_id', 'required': True},
                        'trait_assertion_id': {'name': 'trait_assertion_id',
                                               'required': True}},
         'title': 'KG-Microbe Trait Assertion'})

    trait_assertion_id: str = Field(default=..., title="Trait Assertion ID", description="""Stable key for a generated trait matrix column. Direct METPO phenotype columns use the METPO object CURIE; relation-value columns combine the predicate and object identifiers.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion']} })
    trait_label: Optional[str] = Field(default=None, title="Trait Label", description="""Human-readable label for a generated trait matrix column.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion']} })
    predicate_id: Optional[str] = Field(default=None, title="Predicate ID", description="""Predicate used by the KG-Microbe assertion. Direct phenotype assertions may use biolink:has_phenotype, while relation-value assertions use METPO predicates such as METPO:2000006, METPO:2000103, or METPO:2000517.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion', 'FeatureDefinition']} })
    predicate_name: Optional[str] = Field(default=None, title="Predicate Name", description="""Human-readable predicate label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion', 'FeatureDefinition']} })
    object_id: str = Field(default=..., title="Object ID", description="""KG object identifier for a trait assertion. This may be a METPO phenotype, CHEBI chemical, GO process, pathway identifier, MediaDive medium, or other KG-Microbe object.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion', 'FeatureDefinition']} })
    object_name: Optional[str] = Field(default=None, title="Object Name", description="""Human-readable object label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion', 'FeatureDefinition']} })
    object_category: Optional[str] = Field(default=None, title="Object Category", description="""Biolink category of the KG object.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion', 'FeatureDefinition']} })
    count: int = Field(default=..., title="Count", description="""Number of taxa at or below a row taxon that have the trait assertion.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype', 'KgTraitAssertion', 'FeatureCount']} })
    from_descendants: Optional[bool] = Field(default=None, title="From Descendants", description="""Whether the trait assertion count includes descendant taxa.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype', 'KgTraitAssertion']} })


class TaxonTraitMatrixRow(ConfiguredBaseModel):
    """
    One BugSigDB-reported taxon row in a taxa-by-traits matrix. Trait columns in the generated TSV are dynamic KG-Microbe phenotype or trait assertion labels; the trait_counts slot provides the normalized representation of the same counts.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/schema',
         'slot_usage': {'taxon_id': {'name': 'taxon_id', 'required': True},
                        'taxon_name': {'name': 'taxon_name', 'required': True}},
         'title': 'Taxon Trait Matrix Row'})

    taxon_id: str = Field(default=..., title="Taxon ID", description="""NCBITaxon identifier for a matrix row.""", json_schema_extra = { "linkml_meta": {'domain_of': ['TaxonTraitMatrixRow', 'TaxonOccurrenceRow']} })
    taxon_name: str = Field(default=..., title="Taxon Name", description="""Reported taxon label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon', 'TaxonTraitMatrixRow', 'TaxonOccurrenceRow']} })
    taxonomic_rank: Optional[TaxonomicRankEnum] = Field(default=None, title="Taxonomic Rank", description="""Taxonomic rank of the reported taxon.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon', 'TaxonTraitMatrixRow', 'TaxonOccurrenceRow']} })
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


class FeatureCount(ConfiguredBaseModel):
    """
    One nonnegative count attached to a catalog feature identity.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products', 'title': 'Feature Count'})

    feature_id: str = Field(default=..., title="Feature id", description="""Reference to the feature catalog column identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FeatureCount']} })
    count: int = Field(default=..., title="Count", description="""Nonnegative count in the profile's explicitly declared unit.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['Phenotype', 'KgTraitAssertion', 'FeatureCount']} })


class CountProfileRow(ConfiguredBaseModel):
    """
    Logical row of the raw count profile before presentation filtering or pairing.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Count Profile Row'})

    study_id: str = Field(default=..., title="Study id", description="""Study component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'StudyTableRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    experiment_id: str = Field(default=..., title="Experiment id", description="""Experiment component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    inc_dec: str = Field(default=..., title="Inc dec", description="""Reported increased or decreased direction component of the profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow', 'ORAProfileRow', 'ProfileCoverageRecord']} })
    feature_counts: list[FeatureCount] = Field(default=..., title="Feature counts", description="""Nonzero counts normalized after validating every wide TSV cell; omitted features are verified zeros.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow', 'ORAProfileRow']} })
    trait_filter_id: str = Field(default=..., title="Trait filter id", description="""Hash of the complete ordered feature vocabulary including verified zero coordinates.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow', 'ORAProfileRow', 'ProfileCoverageRecord']} })
    unit: str = Field(default="rolled_up_source_taxon_count", title="Unit", description="""Counts sum descendant-rollup source counts over distinct reported taxon IDs. Repeated reported IDs are deduplicated; the same source annotation can contribute through distinct reported ancestor and descendant IDs.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow', 'ORAProfileRow'],
         'ifabsent': 'string(rolled_up_source_taxon_count)'} })

    @field_validator('inc_dec')
    def pattern_inc_dec(cls, v):
        pattern=re.compile(r"^(inc|dec)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid inc_dec format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid inc_dec format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('trait_filter_id')
    def pattern_trait_filter_id(cls, v):
        pattern=re.compile(r"^sha256:[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid trait_filter_id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid trait_filter_id format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('unit')
    def pattern_unit(cls, v):
        pattern=re.compile(r"^rolled_up_source_taxon_count$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid unit format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid unit format: {v}"
            raise ValueError(err_msg)
        return v


class ORAProfileRow(ConfiguredBaseModel):
    """
    Logical ORA row with explicit distinct matched taxon denominator.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products', 'title': 'ORAProfile Row'})

    study_id: str = Field(default=..., title="Study id", description="""Study component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'StudyTableRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    experiment_id: str = Field(default=..., title="Experiment id", description="""Experiment component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    inc_dec: str = Field(default=..., title="Inc dec", description="""Reported increased or decreased direction component of the profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow', 'ORAProfileRow', 'ProfileCoverageRecord']} })
    sample_size: int = Field(default=..., title="Sample size", description="""Number of distinct reported taxa matched to the lookup.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ORAProfileRow']} })
    feature_counts: list[FeatureCount] = Field(default=..., title="Feature counts", description="""Nonzero presence counts across distinct matched reported taxa; omitted features are verified zeros.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow', 'ORAProfileRow']} })
    trait_filter_id: str = Field(default=..., title="Trait filter id", description="""Hash of the complete ordered feature vocabulary including verified zero coordinates.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow', 'ORAProfileRow', 'ProfileCoverageRecord']} })
    unit: str = Field(default="distinct_matched_reported_taxon_count", title="Unit", description="""Counts distinct matched reported taxa without descendant multiplicity.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow', 'ORAProfileRow'],
         'ifabsent': 'string(distinct_matched_reported_taxon_count)'} })

    @field_validator('inc_dec')
    def pattern_inc_dec(cls, v):
        pattern=re.compile(r"^(inc|dec)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid inc_dec format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid inc_dec format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('trait_filter_id')
    def pattern_trait_filter_id(cls, v):
        pattern=re.compile(r"^sha256:[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid trait_filter_id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid trait_filter_id format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('unit')
    def pattern_unit(cls, v):
        pattern=re.compile(r"^distinct_matched_reported_taxon_count$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid unit format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid unit format: {v}"
            raise ValueError(err_msg)
        return v


class SourceFileIdentity(ConfiguredBaseModel):
    """
    Actual consumed source bytes with an explicit product binding.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Source File Identity'})

    product_id: str = Field(default=..., title="Product id", description="""Reference to the embedded product contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceFileIdentity']} })
    path: str = Field(default=..., title="Path", description="""Resolved source path actually consumed.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceFileIdentity']} })
    sha256: str = Field(default=..., title="Sha256", description="""Lowercase SHA-256 of file bytes or the canonical collection inventory.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceFileIdentity']} })
    size_bytes: int = Field(default=..., title="Size bytes", description="""Total number of regular-file bytes represented by this record.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['SourceFileIdentity']} })

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


class VerifiedRelease(ConfiguredBaseModel):
    """
    Release identity assigned only after checking the pinned source bytes.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Verified Release'})

    version: str = Field(default=..., title="Version", description="""Version of the declared workflow contract.""", json_schema_extra = { "linkml_meta": {'domain_of': ['VerifiedRelease']} })
    doi: str = Field(default=..., title="Doi", description="""Verified release DOI.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'VerifiedRelease', 'StudyTableRow']} })
    commit: str = Field(default=..., title="Commit", description="""Verified upstream release commit.""", json_schema_extra = { "linkml_meta": {'domain_of': ['VerifiedRelease']} })
    url: str = Field(default=..., title="Url", description="""Verified source release URL.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'VerifiedRelease', 'StudyTableRow']} })


class SourceIdentity(ConfiguredBaseModel):
    """
    Source classification with verified byte identities and explicit uncertainty.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Source Identity'})

    recorded_at_utc: datetime  = Field(default=..., title="Recorded at utc", description="""Time source identity was recorded.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceIdentity']} })
    declared_source_mode: str = Field(default=..., title="Declared source mode", description="""User-declared source mode; declaration alone does not verify origin.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceIdentity']} })
    pinned_verified: bool = Field(default=..., title="Pinned verified", description="""Whether exact pinned BugSigDB bytes were verified.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceIdentity']} })
    verification_status: str = Field(default=..., title="Verification status", description="""Whether exact pinned source bytes were verified or source origin remains unknown.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceIdentity']} })
    snapshot_header: Optional[str] = Field(default=None, title="Snapshot header", description="""Source snapshot comment if present; otherwise unknown.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceIdentity']} })
    consumed_files: dict[str, SourceFileIdentity] = Field(default=..., title="Consumed files", description="""File identities keyed by their explicit product identifiers.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceIdentity']} })
    interpretation: str = Field(default=..., title="Interpretation", description="""Scientific scope and annotation-versus-measurement interpretation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceIdentity']} })
    verified_release: Optional[VerifiedRelease] = Field(default=None, title="Verified release", description="""Verified release identity absent when source origin is unverified.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SourceIdentity']} })

    @field_validator('declared_source_mode')
    def pattern_declared_source_mode(cls, v):
        pattern=re.compile(r"^(unknown|fixture|live|pinned)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid declared_source_mode format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid declared_source_mode format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('verification_status')
    def pattern_verification_status(cls, v):
        pattern=re.compile(r"^(pinned_bytes_verified|origin_unverified)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid verification_status format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid verification_status format: {v}"
            raise ValueError(err_msg)
        return v


class StudyTableRow(ConfiguredBaseModel):
    """
    Exact logical row contract for the flattened studies TSV.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Study Table Row'})

    study_id: str = Field(default=..., title="Study id", description="""Study component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'StudyTableRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    source_record_id: str = Field(default=..., title="Source record id", description="""Original BugSigDB source record identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'StudyTableRow']} })
    pmid: str = Field(default=..., title="Pmid", description="""PubMed identifier or an empty TSV value when absent.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'StudyTableRow']} })
    title: str = Field(default=..., title="Title", description="""Reported publication title or an empty TSV value.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'StudyTableRow']} })
    publication_year: str = Field(default=..., title="Publication year", description="""Publication year serialized as text with empty missing values.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'StudyTableRow']} })
    doi: str = Field(default=..., title="Doi", description="""Publication DOI or an empty TSV value.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'VerifiedRelease', 'StudyTableRow']} })
    url: str = Field(default=..., title="Url", description="""Publication URL or an empty TSV value.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Study', 'VerifiedRelease', 'StudyTableRow']} })
    experiment_count: int = Field(default=..., title="Experiment count", description="""Number of exported experiments in this study.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['StudyTableRow']} })
    signature_count: int = Field(default=..., title="Signature count", description="""Number of exported signatures in this study.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['StudyTableRow', 'ExperimentTableRow']} })
    taxon_count: int = Field(default=..., title="Taxon count", description="""Number of reported taxon occurrences preserving repeated occurrences.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['StudyTableRow', 'ExperimentTableRow', 'SignatureTableRow']} })
    input_file: str = Field(default=..., title="Input file", description="""Source study YAML member name.""", json_schema_extra = { "linkml_meta": {'domain_of': ['StudyTableRow']} })


class ExperimentTableRow(ConfiguredBaseModel):
    """
    Exact logical row contract for the flattened experiments TSV.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Experiment Table Row'})

    study_id: str = Field(default=..., title="Study id", description="""Study component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'StudyTableRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    experiment_id: str = Field(default=..., title="Experiment id", description="""Experiment component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    experiment_name: str = Field(default=..., title="Experiment name", description="""Source experiment label or an empty TSV value.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    condition: str = Field(default=..., title="Condition", description="""Reported experimental condition.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    body_site: str = Field(default=..., title="Body site", description="""Reported body site.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    host_species: str = Field(default=..., title="Host species", description="""Reported host species used by explicit cohort selection.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    sequencing_type: str = Field(default=..., title="Sequencing type", description="""Reported sequencing type used by explicit cohort selection (#742), or an empty TSV value.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    location_of_subjects: str = Field(default=..., title="Location of subjects", description="""Reported subject location.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    group_0_name: str = Field(default=..., title="Group 0 name", description="""Reported reference group label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    group_1_name: str = Field(default=..., title="Group 1 name", description="""Reported comparison group label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Experiment', 'ExperimentTableRow']} })
    signature_count: int = Field(default=..., title="Signature count", description="""Number of exported signatures in this experiment.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['StudyTableRow', 'ExperimentTableRow']} })
    taxon_count: int = Field(default=..., title="Taxon count", description="""Number of reported taxon occurrences in this experiment.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['StudyTableRow', 'ExperimentTableRow', 'SignatureTableRow']} })


class SignatureTableRow(ConfiguredBaseModel):
    """
    Exact logical row contract for the flattened signatures TSV.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Signature Table Row'})

    study_id: str = Field(default=..., title="Study id", description="""Study component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'StudyTableRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    experiment_id: str = Field(default=..., title="Experiment id", description="""Experiment component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    signature_id: str = Field(default=..., title="Signature id", description="""Source signature identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SignatureTableRow', 'TaxonOccurrenceRow']} })
    direction: str = Field(default=..., title="Direction", description="""Curated reported taxon abundance direction in Group 1 relative to Group 0. Clinical case-control orientation is not established by this core field.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Signature', 'SignatureTableRow', 'TaxonOccurrenceRow']} })
    taxon_count: int = Field(default=..., title="Taxon count", description="""Number of taxon occurrences reported by this signature.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['StudyTableRow', 'ExperimentTableRow', 'SignatureTableRow']} })

    @field_validator('direction')
    def pattern_direction(cls, v):
        pattern=re.compile(r"^(increased|decreased)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid direction format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid direction format: {v}"
            raise ValueError(err_msg)
        return v


class TaxonOccurrenceRow(ConfiguredBaseModel):
    """
    Exact logical row contract preserving every reported taxon occurrence.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Taxon Occurrence Row'})

    study_id: str = Field(default=..., title="Study id", description="""Study component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'StudyTableRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    experiment_id: str = Field(default=..., title="Experiment id", description="""Experiment component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    signature_id: str = Field(default=..., title="Signature id", description="""Source signature identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['SignatureTableRow', 'TaxonOccurrenceRow']} })
    direction: str = Field(default=..., title="Direction", description="""Reported abundance direction.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Signature', 'SignatureTableRow', 'TaxonOccurrenceRow']} })
    taxon_id: str = Field(default=..., title="Taxon id", description="""Reported taxon identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['TaxonTraitMatrixRow', 'TaxonOccurrenceRow']} })
    taxon_name: str = Field(default=..., title="Taxon name", description="""Reported real taxon name.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon', 'TaxonTraitMatrixRow', 'TaxonOccurrenceRow']} })
    taxonomic_rank: str = Field(default=..., title="Taxonomic rank", description="""Reported taxonomic rank when available.""", json_schema_extra = { "linkml_meta": {'domain_of': ['Taxon', 'TaxonTraitMatrixRow', 'TaxonOccurrenceRow']} })

    @field_validator('direction')
    def pattern_direction(cls, v):
        pattern=re.compile(r"^(increased|decreased)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid direction format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid direction format: {v}"
            raise ValueError(err_msg)
        return v


class ExperimentOntologyTermRow(ConfiguredBaseModel):
    """
    Exact logical row contract for multivalued experimental ontology annotations.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Experiment Ontology Term Row'})

    study_id: str = Field(default=..., title="Study id", description="""Study component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'StudyTableRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    experiment_id: str = Field(default=..., title="Experiment id", description="""Experiment component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    slot: str = Field(default=..., title="Slot", description="""Name of the source scientific ontology field.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ExperimentOntologyTermRow']} })
    term_id: str = Field(default=..., title="Term id", description="""One ontology term identifier retained from the source field.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ExperimentOntologyTermRow']} })


class CoreProductCatalog(ConfiguredBaseModel):
    """
    Scientific identity and coverage metadata accompanying file products.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Core Product Catalog',
         'tree_root': True})

    id: str = Field(default=..., title="Id", description="""Identifier unique within its containing document.""", json_schema_extra = { "linkml_meta": {'domain_of': ['NamedThing', 'Node', 'CoreProductCatalog']} })
    features: Optional[list[FeatureDefinition]] = Field(default=None, title="Features", description="""Predicate-object feature definitions for the matrix vocabulary.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CoreProductCatalog']} })
    coverage: Optional[list[ProfileCoverageRecord]] = Field(default=None, title="Coverage", description="""Per-profile accounting of reported, matched, and annotated taxa.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CoreProductCatalog']} })
    studies: Optional[list[Study]] = Field(default=None, title="Studies", description="""Canonical normalized BugSigDB study records.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CoreProductCatalog']} })
    assertions: Optional[list[KgTraitAssertion]] = Field(default=None, title="Assertions", description="""Canonical KG-Microbe trait assertion records.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CoreProductCatalog']} })
    taxon_trait_rows: Optional[list[TaxonTraitMatrixRow]] = Field(default=None, title="Taxon trait rows", description="""Canonical taxon-to-feature matrix rows.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CoreProductCatalog']} })


class FeatureDefinition(ConfiguredBaseModel):
    """
    Predicate-object feature identity preserving trait and culture-medium semantics.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Feature Definition'})

    column_id: str = Field(default=..., title="Column id", description="""Stable predicate-object column identifier matching the matrix header.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FeatureDefinition']} })
    column_label: str = Field(default=..., title="Column label", description="""Human-readable matrix feature label.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FeatureDefinition']} })
    predicate_id: str = Field(default=..., title="Predicate id", description="""Relation identifier in the KG-Microbe assertion.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion', 'FeatureDefinition']} })
    predicate_name: str = Field(default=..., title="Predicate name", description="""Human-readable relation name.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion', 'FeatureDefinition']} })
    object_id: str = Field(default=..., title="Object id", description="""Trait or culture-medium object identifier.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion', 'FeatureDefinition']} })
    object_name: str = Field(default=..., title="Object name", description="""Human-readable trait or medium name.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion', 'FeatureDefinition']} })
    object_category: str = Field(default=..., title="Object category", description="""KG category of the assertion object.""", json_schema_extra = { "linkml_meta": {'domain_of': ['KgTraitAssertion', 'FeatureDefinition']} })
    trait_key_policy: Optional[str] = Field(default=None, title="Trait key policy", description="""Versioned trait identity rule; required by the current scientific contract and absent only in verified historical formats.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FeatureDefinition']} })
    trait_member_ids_json: Optional[str] = Field(default=None, title="Trait member identifiers JSON", description="""Canonical JSON array of declared raw object CURIEs represented by this feature, not an inventory of observed contributors.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FeatureDefinition']} })
    trait_source_vocabularies_json: Optional[str] = Field(default=None, title="Trait source vocabularies JSON", description="""Canonical JSON array of declared object_id, source_id and source_label mappings; does not imply those sources contributed to every matrix cell.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FeatureDefinition']} })
    feature_class: FeatureClass = Field(default=..., title="Feature class", description="""Whether the feature is a trait annotation or a culture-medium relation.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FeatureDefinition']} })
    annotation_evidence: str = Field(default=..., title="Annotation evidence", description="""The assertion origins the feature's counts may combine, unresolved by this matrix; direct or descendant, and BacDive strain records when the run's strain rollup scope admits the feature.""", json_schema_extra = { "linkml_meta": {'domain_of': ['FeatureDefinition']} })

    @field_validator('trait_key_policy')
    def pattern_trait_key_policy(cls, v):
        pattern=re.compile(r"^(raw_object_v1|facultative_union_v1)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid trait_key_policy format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid trait_key_policy format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('annotation_evidence')
    def pattern_annotation_evidence(cls, v):
        pattern=re.compile(r"^(direct_or_descendant_unresolved|direct_descendant_or_strain_unresolved)$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid annotation_evidence format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid annotation_evidence format: {v}"
            raise ValueError(err_msg)
        return v


class ProfileCoverageRecord(ConfiguredBaseModel):
    """
    Reconciled membership and exclusion accounting for a signature profile.
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://w3id.org/mtp/core-products',
         'title': 'Profile Coverage Record'})

    study_id: str = Field(default=..., title="Study id", description="""Study component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'StudyTableRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    experiment_id: str = Field(default=..., title="Experiment id", description="""Experiment component of the unique profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow',
                       'ORAProfileRow',
                       'ExperimentTableRow',
                       'SignatureTableRow',
                       'TaxonOccurrenceRow',
                       'ExperimentOntologyTermRow',
                       'ProfileCoverageRecord']} })
    inc_dec: str = Field(default=..., title="Inc dec", description="""Reported increased or decreased direction component of the profile key.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow', 'ORAProfileRow', 'ProfileCoverageRecord']} })
    reported_distinct_taxa: int = Field(default=..., title="Reported distinct taxa", description="""Number of distinct taxon identifiers reported in the source profile.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ProfileCoverageRecord']} })
    lookup_matched_taxa: int = Field(default=..., title="Lookup matched taxa", description="""Reported taxon identifiers present in the explicit taxon-trait lookup.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ProfileCoverageRecord']} })
    lookup_unmatched_taxa: int = Field(default=..., title="Lookup unmatched taxa", description="""Reported taxon identifiers absent from the explicit taxon-trait lookup.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ProfileCoverageRecord']} })
    matched_taxa_with_annotations: int = Field(default=..., title="Matched taxa with annotations", description="""Matched taxa with at least one retained feature annotation.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ProfileCoverageRecord']} })
    matched_taxa_without_annotations: int = Field(default=..., title="Matched taxa without annotations", description="""Matched taxa with no retained feature annotations.""", ge=0, json_schema_extra = { "linkml_meta": {'domain_of': ['ProfileCoverageRecord']} })
    coverage_status: CoverageStatus = Field(default=..., title="Coverage status", description="""Explicit coverage category determined by the reconciled counts.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProfileCoverageRecord']} })
    unmatched_taxon_ids: str = Field(default=..., title="Unmatched taxon ids", description="""Sorted unique taxon identifiers encoded as JSON.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProfileCoverageRecord']} })
    matched_without_annotation_taxon_ids: str = Field(default=..., title="Matched without annotation taxon ids", description="""Sorted unique taxon identifiers encoded as JSON.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProfileCoverageRecord']} })
    trait_filter_id: str = Field(default=..., title="Trait filter id", description="""SHA-256 identifying the exact ordered retained feature vocabulary.""", json_schema_extra = { "linkml_meta": {'domain_of': ['CountProfileRow', 'ORAProfileRow', 'ProfileCoverageRecord']} })
    lineage_resolution: str = Field(default=..., title="Lineage resolution", description="""Unresolved direct-versus-descendant lineage; presence does not establish direct evidence.""", json_schema_extra = { "linkml_meta": {'domain_of': ['ProfileCoverageRecord']} })

    @field_validator('trait_filter_id')
    def pattern_trait_filter_id(cls, v):
        pattern=re.compile(r"^sha256:[0-9a-f]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid trait_filter_id format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid trait_filter_id format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('lineage_resolution')
    def pattern_lineage_resolution(cls, v):
        pattern=re.compile(r"^unresolved$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid lineage_resolution format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid lineage_resolution format: {v}"
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
FeatureCount.model_rebuild()
CountProfileRow.model_rebuild()
ORAProfileRow.model_rebuild()
SourceFileIdentity.model_rebuild()
VerifiedRelease.model_rebuild()
SourceIdentity.model_rebuild()
StudyTableRow.model_rebuild()
ExperimentTableRow.model_rebuild()
SignatureTableRow.model_rebuild()
TaxonOccurrenceRow.model_rebuild()
ExperimentOntologyTermRow.model_rebuild()
CoreProductCatalog.model_rebuild()
FeatureDefinition.model_rebuild()
ProfileCoverageRecord.model_rebuild()
