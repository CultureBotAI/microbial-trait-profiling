"""Semantic audit, ontology mapping, and source evidence for Table 1 traits.

This module provides data models, reference constants, and audit functions
evaluating the semantic basis, ontology overlap, and evidence provenance for
the six canonical microbial traits analyzed in BugSigDB profiles.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from bsdbng.paths import DERIVED_DATA_DIR, REPO_ROOT

KG_MICROBE_DIR = REPO_ROOT / "data" / "kg-microbe"
PROVENANCE_FILE = KG_MICROBE_DIR / "provenance.json"
NODES_FILE = KG_MICROBE_DIR / "merged-kg_nodes.tsv"
EDGES_FILE = KG_MICROBE_DIR / "merged-kg_edges.tsv"
TAXA_TRAIT_MATRIX_FILE = DERIVED_DATA_DIR / "taxa_trait_matrix.tsv"
STUDY_TRAIT_PROFILES_FILE = DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv"
DEFAULT_OUTPUT_TSV = DERIVED_DATA_DIR / "table1_trait_semantic_audit.tsv"

# Pinned KG-Microbe checksums from data/kg-microbe/provenance.json
PINNED_KG_PROVENANCE = {
    "merged-kg_edges.tsv": {
        "md5": "555ecf5ba5d89288cdc44e368375374f",  # pragma: allowlist secret
        "size": 1034256091,
        "file_id": "1iIVd-FEleeRIv7CPS0l1F732p3J9JO5d",  # pragma: allowlist secret
        "fetched_at": "2026-08-30T18:06:23+00:00",
    },
    "merged-kg_nodes.tsv": {
        "md5": "dd76b4656f9341db25a7c4afa96cd58e",  # pragma: allowlist secret
        "size": 209545218,
        "file_id": "1lDVjE67E82LMUdNtqWGvyczrHvNu3BdY",  # pragma: allowlist secret
        "fetched_at": "2026-08-30T18:05:46+00:00",
    },
}

# Six canonical Table 1 traits
TABLE_1_TRAIT_KEYS: list[str] = [
    "METPO:1000608",
    "METPO:2000103|GO:0030245",
    "METPO:1000668",
    "METPO:1000683",
    "METPO:2000103|pathways:arsenate_reduction",
    "METPO:1000623",
]


@dataclass(frozen=True)
class TraitSemanticAuditRecord:
    """Detailed semantic audit record for a single Table 1 trait feature."""

    trait_key: str
    feature_name: str
    ontology_id: str
    ontology_label: str
    ontology_definition: str
    kg_edge_predicate: str
    kg_edge_relation: str
    source_knowledge_source: str
    source_file_hashes: str
    original_source_term: str
    evidence_type: str
    taxonomic_resolution: str
    direct_kg_taxa_count: int
    bugsigdb_matrix_taxa_count: int
    profile_support_count: int
    profile_support_pct: float
    evidence_status: str
    decision: str
    approved_wording: str
    unsupported_interpretations: str


@dataclass(frozen=True)
class CoccusSphereComparison:
    """Quantitative comparison of coccus-shaped vs sphere-shaped morphology features."""

    coccus_trait_key: str = "METPO:1000668"
    sphere_trait_key: str = "METPO:1000683"
    ontology_relationship: str = "Sibling classes under METPO:1000666 (cell shape)"
    coccus_direct_taxa: int = 5614
    sphere_direct_taxa: int = 381
    direct_taxa_intersection: int = 266
    direct_taxa_jaccard: float = 0.0464
    coccus_bugsigdb_taxa: int = 675
    sphere_bugsigdb_taxa: int = 283
    bugsigdb_taxa_intersection: int = 262
    bugsigdb_taxa_jaccard: float = 0.3764
    coccus_profile_nonzero_count: int = 4409
    sphere_profile_nonzero_count: int = 3541
    pearson_correlation: float = 0.9110
    spearman_correlation: float = 0.8468
    decision: str = "keep_separate"
    rationale: str = (
        "Keep separate in Table 1 and primary outputs to preserve provenance "
        "traceability to distinct source databases (Madin et al. vs Bactotraits). "
        "Any future collapse requires unique-assertion deduplication before "
        "taxonomic rollup, denominator recalculation, and a separate testing family."
    )


# Curated semantic specifications and approved interpretations for Table 1 traits
TRAIT_SEMANTIC_SPECS: dict[str, dict[str, Any]] = {
    "METPO:1000608": {
        "feature_name": "Facultatively aerobic",
        "ontology_id": "METPO:1000608",
        "ontology_label": "facultatively aerobic",
        "ontology_definition": (
            "An oxygen preference in which growth can occur without oxygen but is "
            "capable of aerobic growth."
        ),
        "kg_edge_predicate": "biolink:has_phenotype",
        "kg_edge_relation": "RO:0002215 (capable of) / RO:0002200 (has phenotype)",
        "source_knowledge_source": (
            "infores:madin_etal (5,496 taxa); infores:bacdive (53 strains)"
        ),
        "original_source_term": "oxygen tolerance = 'facultative' / 'facultatively aerobic'",
        "evidence_type": "observation (curated isolate growth preference)",
        "taxonomic_resolution": (
            "NCBI species/genus isolates; propagated to BugSigDB taxa via NCBI hierarchy"
        ),
        "direct_kg_taxa_count": 5496,
        "bugsigdb_matrix_taxa_count": 861,
        "profile_support_count": 4491,
        "profile_support_pct": 77.05,
        "evidence_status": "measured (isolate phenotype) + propagated (taxonomic rollup)",
        "decision": "keep_separate",
        "approved_wording": (
            "Taxa annotated as facultatively aerobic (capable of growth in both the presence "
            "and absence of oxygen)."
        ),
        "unsupported_interpretations": (
            "Do not infer anaerobic respiration using nitrate or tetrathionate. Do not infer "
            "specific respiratory electron transport chains (e.g., nar, nap, ttr). Do not "
            "infer mucosal inflammation or host-derived oxidant availability."
        ),
    },
    "METPO:2000103|GO:0030245": {
        "feature_name": "Cellulose catabolic process",
        "ontology_id": "GO:0030245",
        "ontology_label": "cellulose catabolic process",
        "ontology_definition": (
            "The chemical reactions and pathways resulting in the breakdown of cellulose, "
            "a linear beta1-4 glucan of molecular mass 50-400 kDa with the pyranose units "
            "in the -4C1 conformation."
        ),
        "kg_edge_predicate": "METPO:2000103 (capable of)",
        "kg_edge_relation": "RO:0002215 (capable of)",
        "source_knowledge_source": "infores:madin_etal (186 taxa)",
        "original_source_term": "pathway/function = 'cellulose catabolic process' [GO:0030245]",
        "evidence_type": "observation (curated pathway capability assertion)",
        "taxonomic_resolution": (
            "NCBI species/genus isolates; propagated to BugSigDB taxa via NCBI hierarchy"
        ),
        "direct_kg_taxa_count": 186,
        "bugsigdb_matrix_taxa_count": 164,
        "profile_support_count": 2963,
        "profile_support_pct": 50.83,
        "evidence_status": "measured (pathway annotation) + propagated (taxonomic rollup)",
        "decision": "keep_separate",
        "approved_wording": (
            "Taxa annotated with pathway capability for cellulose catabolic processes."
        ),
        "unsupported_interpretations": (
            "Do not infer demonstrated in vivo cellulose degradation in patients. Do not "
            "equate cellulose catabolism with butyrate or short-chain fatty acid (SCFA) "
            "production (end-products frequently comprise acetate, succinate, lactate, or ethanol)."
        ),
    },
    "METPO:1000668": {
        "feature_name": "Coccus shaped",
        "ontology_id": "METPO:1000668",
        "ontology_label": "coccus shaped",
        "ontology_definition": (
            "A cell shape in which an organism has a spherical or nearly spherical morphology, "
            "with roughly equal dimensions in all directions."
        ),
        "kg_edge_predicate": "biolink:has_phenotype",
        "kg_edge_relation": "RO:0002200 (has phenotype)",
        "source_knowledge_source": (
            "infores:madin_etal (5,614 taxa); infores:bacdive (1,476 strains)"
        ),
        "original_source_term": "cell_shape = 'coccus'",
        "evidence_type": "observation (microscopic isolate morphology)",
        "taxonomic_resolution": (
            "NCBI species/genus isolates; propagated to BugSigDB taxa via NCBI hierarchy"
        ),
        "direct_kg_taxa_count": 5614,
        "bugsigdb_matrix_taxa_count": 675,
        "profile_support_count": 4409,
        "profile_support_pct": 75.64,
        "evidence_status": "measured (isolate morphology) + propagated (taxonomic rollup)",
        "decision": "keep_separate",
        "approved_wording": (
            "Taxa annotated with coccus-shaped cell morphology (primarily derived from "
            "Madin et al. isolate compilation)."
        ),
        "unsupported_interpretations": (
            "Do not merge with sphere-shaped by naive string rewriting or profile summation. "
            "Do not claim biological non-equivalence without noting that differences reflect "
            "source database curation boundaries (Madin et al. vs Bactotraits)."
        ),
    },
    "METPO:1000683": {
        "feature_name": "Sphere shaped",
        "ontology_id": "METPO:1000683",
        "ontology_label": "sphere shaped",
        "ontology_definition": (
            "A cell shape defined in METPO as a subclass of cell shape (METPO:1000666); "
            "unpopulated narrative definition in KG nodes (synonyms: S_sphere, sphere-shaped)."
        ),
        "kg_edge_predicate": "biolink:has_phenotype",
        "kg_edge_relation": "RO:0002200 (has phenotype)",
        "source_knowledge_source": "infores:bactotraits (381 taxa); infores:bacdive (11 strains)",
        "original_source_term": "shape code = 'S_sphere' (Bactotraits / BacDive)",
        "evidence_type": "observation (microscopic isolate morphology)",
        "taxonomic_resolution": (
            "NCBI species/genus isolates; propagated to BugSigDB taxa via NCBI hierarchy"
        ),
        "direct_kg_taxa_count": 381,
        "bugsigdb_matrix_taxa_count": 283,
        "profile_support_count": 3541,
        "profile_support_pct": 60.75,
        "evidence_status": "measured (isolate morphology) + propagated (taxonomic rollup)",
        "decision": "keep_separate",
        "approved_wording": (
            "Taxa annotated with sphere-shaped cell morphology (derived from Bactotraits "
            "and BacDive morphology records)."
        ),
        "unsupported_interpretations": (
            "Do not treat as a biologically distinct morphology from coccus shaped. "
            "Do not sum profile counts with coccus shaped without deduplicating shared taxa."
        ),
    },
    "METPO:2000103|pathways:arsenate_reduction": {
        "feature_name": "Arsenate reduction",
        "ontology_id": "pathways:arsenate_reduction",
        "ontology_label": "arsenate_reduction",
        "ontology_definition": (
            "Pathway annotation in KG-Microbe for microbial reduction of arsenate; "
            "lacks narrative definition in KG nodes; mapped under METPO:2000103 (capable of)."
        ),
        "kg_edge_predicate": "METPO:2000103 (capable of)",
        "kg_edge_relation": "RO:0002215 (capable of)",
        "source_knowledge_source": "infores:madin_etal (51 taxa)",
        "original_source_term": "pathway = 'arsenate_reduction' (Madin et al. / FAPROTAX)",
        "evidence_type": "observation (pathway compilation assertion)",
        "taxonomic_resolution": (
            "NCBI species isolates; propagated to BugSigDB taxa via NCBI hierarchy"
        ),
        "direct_kg_taxa_count": 51,
        "bugsigdb_matrix_taxa_count": 58,
        "profile_support_count": 1770,
        "profile_support_pct": 30.37,
        "evidence_status": (
            "unresolved mechanism (0/51 KG edges resolve gene or pathway mechanism; "
            "unstratified mixture of arrAB dissimilatory respiration and arsC detoxification) "
            "+ propagated (taxonomic rollup)"
        ),
        "decision": "unresolved",
        "approved_wording": (
            "Taxa annotated with pathway capability for arsenate reduction (derived from "
            "Madin et al. isolate compilation; mechanism unresolved between dissimilatory "
            "respiration and cytoplasmic detoxification)."
        ),
        "unsupported_interpretations": (
            "Do not classify as arsC detoxification. Do not classify as arrAB dissimilatory "
            "respiration. Do not infer gene presence, plasmid carriage, patient arsenic exposure, "
            "or antimicrobial co-selection from this taxon-level pathway label."
        ),
    },
    "METPO:1000623": {
        "feature_name": "Moderately halophilic",
        "ontology_id": "METPO:1000623",
        "ontology_label": "moderately halophilic",
        "ontology_definition": (
            "A halophily preference where growth and proliferation requires high levels "
            "of sodium chloride, usually above or about 0.2 M."
        ),
        "kg_edge_predicate": "biolink:has_phenotype",
        "kg_edge_relation": "RO:0002200 (has phenotype)",
        "source_knowledge_source": "infores:madin_etal (426 taxa); infores:bacdive (4 strains)",
        "original_source_term": "range_salinity / halophily = 'moderate-halophilic'",
        "evidence_type": "observation (curated phenotypic salinity preference)",
        "taxonomic_resolution": (
            "NCBI species/genus isolates; propagated to BugSigDB taxa via NCBI hierarchy"
        ),
        "direct_kg_taxa_count": 426,
        "bugsigdb_matrix_taxa_count": 337,
        "profile_support_count": 3424,
        "profile_support_pct": 58.74,
        "evidence_status": (
            "measured (qualitative halophily category in Madin et al.) + unresolved "
            "(missing quantitative salt concentration, optimum, range, medium, "
            "and culture conditions) + propagated (taxonomic rollup)"
        ),
        "decision": "keep_separate",
        "approved_wording": (
            "Taxa annotated as moderately halophilic (defined in METPO as requiring high levels "
            "of sodium chloride, typically above or about 0.2 M NaCl)."
        ),
        "unsupported_interpretations": (
            "Do not read as 'halotolerant' (halophily requires salt for growth, whereas "
            "halotolerance permits growth without requiring salt). Do not claim specific "
            "quantitative salt ranges, growth media, or culture conditions from KG edges."
        ),
    },
}


def build_semantic_audit_records() -> list[TraitSemanticAuditRecord]:
    """Build structured audit records for the six Table 1 traits."""
    edges_provenance = PINNED_KG_PROVENANCE["merged-kg_edges.tsv"]
    nodes_provenance = PINNED_KG_PROVENANCE["merged-kg_nodes.tsv"]
    hashes_str = f"edges_md5:{edges_provenance['md5']}|nodes_md5:{nodes_provenance['md5']}"

    records: list[TraitSemanticAuditRecord] = []
    for trait_key in TABLE_1_TRAIT_KEYS:
        spec = TRAIT_SEMANTIC_SPECS[trait_key]
        records.append(
            TraitSemanticAuditRecord(
                trait_key=trait_key,
                feature_name=spec["feature_name"],
                ontology_id=spec["ontology_id"],
                ontology_label=spec["ontology_label"],
                ontology_definition=spec["ontology_definition"],
                kg_edge_predicate=spec["kg_edge_predicate"],
                kg_edge_relation=spec["kg_edge_relation"],
                source_knowledge_source=spec["source_knowledge_source"],
                source_file_hashes=hashes_str,
                original_source_term=spec["original_source_term"],
                evidence_type=spec["evidence_type"],
                taxonomic_resolution=spec["taxonomic_resolution"],
                direct_kg_taxa_count=spec["direct_kg_taxa_count"],
                bugsigdb_matrix_taxa_count=spec["bugsigdb_matrix_taxa_count"],
                profile_support_count=spec["profile_support_count"],
                profile_support_pct=spec["profile_support_pct"],
                evidence_status=spec["evidence_status"],
                decision=spec["decision"],
                approved_wording=spec["approved_wording"],
                unsupported_interpretations=spec["unsupported_interpretations"],
            )
        )
    return records


def write_semantic_audit_tsv(
    records: list[TraitSemanticAuditRecord],
    output_path: Path = DEFAULT_OUTPUT_TSV,
) -> Path:
    """Export semantic audit records to a tab-separated values file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        raise ValueError("Cannot write empty semantic audit records")

    fieldnames = list(asdict(records[0]).keys())
    with output_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for record in records:
            writer.writerow(asdict(record))

    return output_path


def read_semantic_audit_tsv(tsv_path: Path = DEFAULT_OUTPUT_TSV) -> list[TraitSemanticAuditRecord]:
    """Read and validate semantic audit records from a TSV file."""
    if not tsv_path.exists():
        raise FileNotFoundError(f"Semantic audit TSV not found: {tsv_path}")

    records: list[TraitSemanticAuditRecord] = []
    with tsv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            records.append(
                TraitSemanticAuditRecord(
                    trait_key=row["trait_key"],
                    feature_name=row["feature_name"],
                    ontology_id=row["ontology_id"],
                    ontology_label=row["ontology_label"],
                    ontology_definition=row["ontology_definition"],
                    kg_edge_predicate=row["kg_edge_predicate"],
                    kg_edge_relation=row["kg_edge_relation"],
                    source_knowledge_source=row["source_knowledge_source"],
                    source_file_hashes=row["source_file_hashes"],
                    original_source_term=row["original_source_term"],
                    evidence_type=row["evidence_type"],
                    taxonomic_resolution=row["taxonomic_resolution"],
                    direct_kg_taxa_count=int(row["direct_kg_taxa_count"]),
                    bugsigdb_matrix_taxa_count=int(row["bugsigdb_matrix_taxa_count"]),
                    profile_support_count=int(row["profile_support_count"]),
                    profile_support_pct=float(row["profile_support_pct"]),
                    evidence_status=row["evidence_status"],
                    decision=row["decision"],
                    approved_wording=row["approved_wording"],
                    unsupported_interpretations=row["unsupported_interpretations"],
                )
            )
    return records


def verify_provenance_checksums() -> bool:
    """Verify that data/kg-microbe/provenance.json matches pinned checksums."""
    if not PROVENANCE_FILE.exists():
        return False

    with PROVENANCE_FILE.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    files = data.get("files", {})
    for filename, pinned in PINNED_KG_PROVENANCE.items():
        entry = files.get(filename)
        if entry is None:
            return False
        if entry.get("md5") != pinned["md5"] or entry.get("size") != pinned["size"]:
            return False

    return True
