"""Build a BugSigDB taxon by KG-Microbe trait count matrix."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from mtp.paths import KG_MICROBE_DIR, REPO_ROOT
from mtp.trait_key_policy import (
    DEFAULT_TRAIT_KEY_POLICY,
    TRAIT_KEY_POLICIES,
    TRAIT_POLICY_METADATA_COLUMNS,
    TraitKeyPolicy,
    canonical_object,
    canonical_object_name,
    parse_trait_key_policy,
    policy_metadata,
)
from mtp.trait_key_policy import (
    FACULTATIVE_TRAIT_ID as FACULTATIVE_TRAIT_ID,
)

TAXON_PREFIX = "NCBITaxon:"
STRAIN_PREFIX = "kgmicrobe.strain:"
TRAIT_PREFIX = "METPO:"
SUBCLASS_OF = "biolink:subclass_of"
TRAIT_KEY_SEPARATOR = "|"

DEFAULT_TAXA_PATH = REPO_ROOT / "data" / "derived" / "bugsigdb_taxa.tsv"
DEFAULT_KG_DIR = KG_MICROBE_DIR
DEFAULT_OUTPUT_PATH = REPO_ROOT / "data" / "derived" / "taxa_trait_matrix.tsv"
EDGES_FILENAME = "merged-kg_edges.tsv"
NODES_FILENAME = "merged-kg_nodes.tsv"

# Which BacDive strain trait edges count at the strains' NCBITaxon nodes (#788): none
# (off), bare METPO phenotype classes (phenotype), or every edge that does not assert an
# absence (full).
StrainRollup = Literal["off", "phenotype", "full"]
STRAIN_ROLLUPS: tuple[StrainRollup, ...] = ("off", "phenotype", "full")
# The owner's primary scope (#788). The core workflow and the legacy matrix builder and
# ledger default to it; off is the sensitivity scope and full the exploratory one.
PRIMARY_STRAIN_ROLLUP: StrainRollup = "phenotype"

# Handling and clinical classes that BacDive strain records assert, not ecophysiological
# traits, so the phenotype scope never rolls them up (#788). Labels as in the pinned
# node table; `just verify-trait-terms` checks them.
PHENOTYPE_EXCLUDED_TRAITS = {
    "METPO:1001102": "biosafety level 1",
    "METPO:1001103": "biosafety level 2",
    "METPO:1001104": "biosafety level 3",
    "METPO:1004002": "animal pathogen",
    "METPO:1004003": "plant pathogen",
    "METPO:1004004": "human pathogen",
}

# Every METPO relation in the pinned node table that asserts an activity or capability
# is absent (#788). `just verify-trait-terms` checks the labels and that none is missing.
NEGATED_PREDICATES = {
    "METPO:2000021": "does not use for aerobic catabolization",
    "METPO:2000022": "does not use for aerobic growth",
    "METPO:2000023": "does not use for anaerobic catabolization",
    "METPO:2000024": "does not use for anaerobic growth",
    "METPO:2000025": "does not use for anaerobic growth in the dark",
    "METPO:2000026": "does not use for anaerobic growth with light",
    "METPO:2000027": "does not assimilate",
    "METPO:2000028": "does not build acid from",
    "METPO:2000029": "does not build base from",
    "METPO:2000030": "does not build gas from",
    "METPO:2000031": "does not use as carbon source",
    "METPO:2000033": "does not degrade",
    "METPO:2000034": "does not use as electron acceptor",
    "METPO:2000035": "does not use as electron donor",
    "METPO:2000036": "does not use as energy source",
    "METPO:2000037": "does not ferment",
    "METPO:2000038": "does not use for growth",
    "METPO:2000039": "does not hydrolyze",
    "METPO:2000040": "does not use as nitrogen source",
    "METPO:2000041": "does not use in other way",
    "METPO:2000042": "does not oxidize",
    "METPO:2000044": "does not reduce",
    "METPO:2000045": "is not required for growth",
    "METPO:2000046": "does not use for respiration",
    "METPO:2000047": "does not use as sulfur source",
    "METPO:2000220": "does not disproportionate",
    "METPO:2000222": "does not produce",
    "METPO:2000227": "does not transport",
    "METPO:2000228": "does not import",
    "METPO:2000229": "does not export",
    "METPO:2000230": "does not accumulate",
    "METPO:2000231": "does not sequester",
    "METPO:2000232": "does not compartmentalize",
    "METPO:2000303": "does not show activity of",
    "METPO:2000518": "does not grow in",
}


@dataclass(frozen=True)
class ReportedTaxon:
    """A unique BugSigDB-reported taxon used as an output matrix row."""

    taxon_id: str
    taxon_name: str
    taxonomic_rank: str


@dataclass(frozen=True)
class MatrixStats:
    """Summary of a generated taxon-trait matrix."""

    reported_taxa: int
    traits: int
    output_path: Path
    column_metadata_path: Path
    edges_path: Path
    nodes_path: Path


@dataclass(frozen=True)
class NodeMetadata:
    """Label metadata for a KG node used in matrix column descriptions."""

    name: str
    category: str


@dataclass(frozen=True)
class KgTraitEdges:
    """Taxonomy links and trait sets read from KG edges.

    ``strain_taxa`` maps each BacDive strain to the NCBITaxon nodes it is a
    ``biolink:subclass_of``, and ``strain_traits`` holds the trait keys of every strain
    with trait edges, only those a scope admits when one was given (#774, #788). Both are
    empty unless strains were read.
    """

    child_to_parents: dict[str, set[str]]
    taxon_traits: dict[str, set[str]]
    strain_taxa: dict[str, set[str]]
    strain_traits: dict[str, set[str]]


def _is_taxon(value: str) -> bool:
    return value.startswith(TAXON_PREFIX)


def _is_strain(value: str) -> bool:
    return value.startswith(STRAIN_PREFIX)


def _is_trait(value: str) -> bool:
    return value.startswith(TRAIT_PREFIX)


def _clean_label(value: str) -> str:
    return " ".join(value.split())


def resolve_kg_file(filename: str, kg_dir: Path, explicit_path: Path | None = None) -> Path:
    """Resolve a KG file from an explicit path or from kg_dir, failing clearly if missing."""
    if explicit_path is not None:
        if explicit_path.exists():
            return explicit_path
        raise FileNotFoundError(f"KG file not found at explicit path: {explicit_path}")

    kg_candidate = kg_dir / filename
    if kg_candidate.exists():
        return kg_candidate

    raise FileNotFoundError(
        f"KG file not found: expected {filename} in kg_dir {kg_dir}. "
        f"Pass an explicit path or ensure the file exists in that directory."
    )


def read_reported_taxa(taxa_path: Path) -> list[ReportedTaxon]:
    """Read unique increased/decreased BugSigDB taxa from data/derived/bugsigdb_taxa.tsv."""
    taxa_by_id: dict[str, ReportedTaxon] = {}
    required_columns = {"direction", "taxon_id", "taxon_name", "taxonomic_rank"}

    with taxa_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = required_columns.difference(reader.fieldnames or [])
        if missing:
            missing_text = ", ".join(sorted(missing))
            raise ValueError(f"{taxa_path} is missing required columns: {missing_text}")

        for row in reader:
            direction = (row.get("direction") or "").strip().lower()
            if direction not in {"increased", "decreased"}:
                continue

            taxon_id = (row.get("taxon_id") or "").strip()
            if not _is_taxon(taxon_id):
                continue

            taxon_name = _clean_label((row.get("taxon_name") or "").strip()) or taxon_id
            taxonomic_rank = _clean_label((row.get("taxonomic_rank") or "").strip())
            taxa_by_id.setdefault(taxon_id, ReportedTaxon(taxon_id, taxon_name, taxonomic_rank))

    return sorted(
        taxa_by_id.values(),
        key=lambda taxon: (taxon.taxon_name.casefold(), taxon.taxon_id),
    )


def read_node_metadata(nodes_path: Path, node_ids: set[str]) -> dict[str, NodeMetadata]:
    """Read labels and categories for KG nodes needed by output trait columns."""
    node_metadata: dict[str, NodeMetadata] = {}
    with nodes_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 3:
                continue
            node_id = row[0].strip()
            if node_id in node_ids:
                category = row[1].strip() if len(row) > 1 else ""
                node_metadata[node_id] = NodeMetadata(
                    name=_clean_label(row[2]),
                    category=category,
                )
    return node_metadata


def parse_strain_rollup(value: object) -> StrainRollup:
    """Return a declared strain rollup scope, refusing any other value."""
    for scope in STRAIN_ROLLUPS:
        if value == scope:
            return scope
    raise ValueError(f"strain_rollup must be one of {', '.join(STRAIN_ROLLUPS)}: {value!r}")


def make_trait_key(
    predicate: str, obj: str, *, trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY
) -> str:
    """Create the output trait key for a KG taxon assertion.

    A METPO object is its own key, except under a negated relation, which keeps its
    predicate so that the absence of a phenotype is never keyed as the phenotype (#788).
    The explicit policy either preserves source object IDs or unifies the three
    declared facultative oxygen-preference terms as METPO:1000608 (#1129).
    """
    obj = canonical_object(obj, trait_key_policy=trait_key_policy)
    if _is_trait(obj) and predicate not in NEGATED_PREDICATES:
        return obj
    return f"{predicate}{TRAIT_KEY_SEPARATOR}{obj}"


def split_trait_key(trait_key: str) -> tuple[str | None, str]:
    """Return optional predicate and object IDs encoded in a trait key."""
    if TRAIT_KEY_SEPARATOR not in trait_key:
        return None, trait_key
    predicate, obj = trait_key.split(TRAIT_KEY_SEPARATOR, 1)
    return predicate, obj


def is_negated(trait_key: str) -> bool:
    """Whether a trait key asserts that an activity or capability is absent (#788)."""
    predicate, _ = split_trait_key(trait_key)
    return predicate in NEGATED_PREDICATES


def strain_trait_rolls_up(trait_key: str, strain_rollup: StrainRollup) -> bool:
    """Whether a BacDive strain's trait key counts at its NCBITaxon node under a scope.

    ``phenotype`` takes bare METPO phenotype classes, without the handling and clinical
    classes of PHENOTYPE_EXCLUDED_TRAITS; ``full`` takes every key that does not assert
    an absence; ``off`` takes none (#788).
    """
    if strain_rollup == "phenotype":
        return TRAIT_KEY_SEPARATOR not in trait_key and trait_key not in PHENOTYPE_EXCLUDED_TRAITS
    if strain_rollup == "full":
        return not is_negated(trait_key)
    return False


def node_ids_for_trait_keys(trait_keys: set[str]) -> set[str]:
    """Return KG node IDs needed to label trait keys."""
    node_ids: set[str] = set()
    for trait_key in trait_keys:
        predicate, obj = split_trait_key(trait_key)
        node_ids.add(obj)
        if predicate is not None:
            node_ids.add(predicate)
    return node_ids


def read_kg_trait_edges(
    edges_path: Path,
    *,
    strains: bool = False,
    strain_rollup: StrainRollup | None = None,
    trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY,
) -> KgTraitEdges:
    """Read taxonomy parent links and direct trait assertions from KG edges.

    Trait sets are keyed by NCBITaxon subject. With ``strains``, the trait edges of
    BacDive strain subjects (``kgmicrobe.strain:``) and the strains' NCBITaxon
    ``subclass_of`` links are also kept, apart from the taxon sets (#774); a
    ``strain_rollup`` scope keeps only the strain keys that scope rolls up (#788).

    Memory: materializes the full taxonomy graph plus per-taxon trait sets in
    Python dicts. On the current ~1 GB merged-kg_edges.tsv this peaks around
    2-3 GB of resident memory, and around 6 GB with strains; switch to a
    streaming / two-pass implementation if that becomes a constraint.
    """
    trait_key_policy = parse_trait_key_policy(trait_key_policy)
    child_to_parents: defaultdict[str, set[str]] = defaultdict(set)
    taxon_traits: defaultdict[str, set[str]] = defaultdict(set)
    strain_taxa: defaultdict[str, set[str]] = defaultdict(set)
    strain_traits: defaultdict[str, set[str]] = defaultdict(set)

    with edges_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 4:
                continue

            subject = row[1].strip()
            predicate = row[2].strip()
            obj = row[3].strip()

            if predicate == SUBCLASS_OF and _is_taxon(subject) and _is_taxon(obj):
                child_to_parents[subject].add(obj)
            elif _is_taxon(subject) and (_is_trait(obj) or _is_trait(predicate)):
                taxon_traits[subject].add(
                    make_trait_key(predicate, obj, trait_key_policy=trait_key_policy)
                )
            elif strains and _is_strain(subject):
                if predicate == SUBCLASS_OF and _is_taxon(obj):
                    strain_taxa[subject].add(obj)
                elif _is_trait(obj) or _is_trait(predicate):
                    key = make_trait_key(predicate, obj, trait_key_policy=trait_key_policy)
                    keys = strain_traits[subject]
                    if strain_rollup is None or strain_trait_rolls_up(key, strain_rollup):
                        keys.add(key)

    return KgTraitEdges(
        dict(child_to_parents), dict(taxon_traits), dict(strain_taxa), dict(strain_traits)
    )


def roll_up_strain_traits(
    taxon_traits: dict[str, set[str]],
    strain_taxa: dict[str, set[str]],
    strain_traits: dict[str, set[str]],
) -> dict[str, set[str]]:
    """Count each strain's trait keys as assertions on the NCBITaxon nodes it links to.

    A node carries a trait once, whether from its own edges or from any number of its
    strains, so a strain adds traits to its species without becoming another annotated
    taxon: the matrix still counts trait-bearing NCBITaxon nodes. A strain linked to
    several nodes adds its traits to each, and one linked to none adds nothing (#774).
    The arguments are left unchanged.
    """
    rolled = {taxon_id: set(traits) for taxon_id, traits in taxon_traits.items()}
    for strain_id, traits in strain_traits.items():
        for taxon_id in strain_taxa.get(strain_id, ()):
            rolled.setdefault(taxon_id, set()).update(traits)
    return rolled


def scoped_taxon_traits(
    taxon_traits: dict[str, set[str]],
    strain_taxa: dict[str, set[str]],
    strain_traits: dict[str, set[str]],
    strain_rollup: StrainRollup,
) -> dict[str, set[str]]:
    """The trait set each NCBITaxon node carries under a strain rollup scope (#788).

    The strain keys the scope admits count at the strains' linked nodes, once per node
    (roll_up_strain_traits). ``full`` also drops every negated key, the nodes' own
    included, so no node carries both an activity and its absence. ``off`` returns the
    direct sets unchanged. Every reader of strain records shares this rule.
    """
    if strain_rollup == "off":
        return taxon_traits
    admitted = {
        strain: kept
        for strain, keys in strain_traits.items()
        if (kept := {key for key in keys if strain_trait_rolls_up(key, strain_rollup)})
    }
    rolled = roll_up_strain_traits(taxon_traits, strain_taxa, admitted)
    if strain_rollup == "phenotype":
        return rolled
    affirmed = {node: {key for key in keys if not is_negated(key)} for node, keys in rolled.items()}
    return {node: keys for node, keys in affirmed.items() if keys}


def read_kg_edges(
    edges_path: Path,
    *,
    strain_rollup: StrainRollup = "off",
    trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY,
) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Read taxonomy parent links and per-taxon trait sets from KG edges.

    By default only NCBITaxon subjects carry traits, so this reader returns the direct
    annotations. A ``strain_rollup`` scope adds the strain keys it admits at the strains'
    NCBITaxon nodes (scoped_taxon_traits); the matrix builders pass the primary scope
    unless told otherwise (#774, #788).
    """
    strains = strain_rollup != "off"
    edges = read_kg_trait_edges(
        edges_path, strains=strains, strain_rollup=strain_rollup, trait_key_policy=trait_key_policy
    )
    return edges.child_to_parents, scoped_taxon_traits(
        edges.taxon_traits, edges.strain_taxa, edges.strain_traits, strain_rollup
    )


def _target_ancestors(
    taxon_id: str,
    target_ids: set[str],
    child_to_parents: dict[str, set[str]],
    cache: dict[str, frozenset[str]],
) -> frozenset[str]:
    cached = cache.get(taxon_id)
    if cached is not None:
        return cached

    found: set[str] = set()
    seen: set[str] = set()
    stack = [taxon_id]
    while stack:
        current = stack.pop()
        if current in seen:
            continue
        seen.add(current)

        cached_current = cache.get(current)
        if cached_current is not None:
            found.update(cached_current)
            continue

        if current in target_ids:
            found.add(current)
        stack.extend(child_to_parents.get(current, ()))

    result = frozenset(found)
    cache[taxon_id] = result
    return result


def compute_trait_counts(
    reported_taxa: list[ReportedTaxon],
    child_to_parents: dict[str, set[str]],
    taxon_traits: dict[str, set[str]],
) -> dict[str, Counter[str]]:
    """Count trait-bearing taxa at or below each reported taxon."""
    target_ids = {taxon.taxon_id for taxon in reported_taxa}
    counts = {taxon_id: Counter[str]() for taxon_id in target_ids}
    ancestor_cache: dict[str, frozenset[str]] = {}

    for trait_taxon_id, traits in taxon_traits.items():
        ancestors = _target_ancestors(
            trait_taxon_id,
            target_ids,
            child_to_parents,
            ancestor_cache,
        )
        for target_id in ancestors:
            counts[target_id].update(traits)

    return counts


def trait_column_id(trait_key: str) -> str:
    """Return the compact stable TSV column identifier for a trait key."""
    predicate, obj = split_trait_key(trait_key)
    if predicate is None:
        return obj
    return f"{predicate}|{obj}"


def _node_label(
    node_id: str,
    node_metadata: dict[str, NodeMetadata],
    *,
    trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY,
) -> str:
    node_name = canonical_object_name(
        node_id,
        node_metadata.get(node_id, NodeMetadata("", "")).name,
        trait_key_policy=trait_key_policy,
    )
    if node_name:
        return f"{node_name} {node_id}"
    return node_id


def trait_column_label(
    trait_key: str,
    node_metadata: dict[str, NodeMetadata],
    *,
    trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY,
) -> str:
    """Return the human-readable label for a KG trait."""
    predicate, obj = split_trait_key(trait_key)
    if predicate is None:
        return _node_label(obj, node_metadata, trait_key_policy=trait_key_policy)

    predicate_name = node_metadata.get(predicate, NodeMetadata(predicate, "")).name or predicate
    return (
        f"{predicate_name}: "
        f"{_node_label(obj, node_metadata, trait_key_policy=trait_key_policy)} [{predicate}]"
    )


def matrix_rows(
    reported_taxa: list[ReportedTaxon],
    counts: dict[str, Counter[str]],
    node_metadata: dict[str, NodeMetadata],
    *,
    trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY,
) -> tuple[list[str], list[list[str]]]:
    """Build TSV headers and rows for the wide taxa-by-traits matrix."""

    def trait_sort_key(trait_id: str) -> tuple[str, str]:
        trait_label = trait_column_label(trait_id, node_metadata, trait_key_policy=trait_key_policy)
        return trait_label.casefold(), trait_id

    trait_ids: list[str] = sorted(
        {
            trait_id
            for taxon_counts in counts.values()
            for trait_id, count in taxon_counts.items()
            if count > 0
        },
        key=trait_sort_key,
    )

    headers = ["taxon_id", "taxon_name", "taxonomic_rank"]
    headers.extend(trait_column_id(trait_id) for trait_id in trait_ids)

    rows: list[list[str]] = []
    for taxon in reported_taxa:
        taxon_counts = counts.get(taxon.taxon_id, Counter[str]())
        row = [taxon.taxon_id, taxon.taxon_name, taxon.taxonomic_rank]
        row.extend(str(taxon_counts.get(trait_id, 0)) for trait_id in trait_ids)
        rows.append(row)

    return headers, rows


def column_metadata_rows(
    counts: dict[str, Counter[str]],
    node_metadata: dict[str, NodeMetadata],
    *,
    trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY,
    include_policy_metadata: bool = True,
) -> tuple[list[str], list[list[str]]]:
    """Describe columns, preserving literal node labels and declaring feature semantics.

    Historical verification selects the old seven-column layout explicitly; its
    scientific policy is independent of its serialization format (#1129).
    """
    trait_key_policy = parse_trait_key_policy(trait_key_policy)

    def trait_sort_key(trait_id: str) -> tuple[str, str]:
        trait_label = trait_column_label(trait_id, node_metadata, trait_key_policy=trait_key_policy)
        return trait_label.casefold(), trait_id

    trait_ids: list[str] = sorted(
        {
            trait_id
            for taxon_counts in counts.values()
            for trait_id, count in taxon_counts.items()
            if count > 0
        },
        key=trait_sort_key,
    )

    headers = [
        "column_id",
        "column_label",
        "predicate_id",
        "predicate_name",
        "object_id",
        "object_name",
        "object_category",
    ]
    if include_policy_metadata:
        headers.extend(TRAIT_POLICY_METADATA_COLUMNS)
    rows: list[list[str]] = []
    for trait_id in trait_ids:
        predicate, obj = split_trait_key(trait_id)
        predicate_metadata = (
            node_metadata.get(predicate, NodeMetadata("", "")) if predicate else None
        )
        object_metadata = node_metadata.get(obj, NodeMetadata("", ""))
        rows.append(
            [
                trait_column_id(trait_id),
                trait_column_label(trait_id, node_metadata, trait_key_policy=trait_key_policy),
                predicate or "",
                predicate_metadata.name if predicate_metadata else "",
                obj,
                canonical_object_name(obj, object_metadata.name, trait_key_policy=trait_key_policy),
                object_metadata.category,
            ]
        )
        if include_policy_metadata:
            rows[-1].extend(policy_metadata(obj, trait_key_policy=trait_key_policy))
    return headers, rows


def write_matrix(output_path: Path, headers: list[str], rows: list[list[str]]) -> None:
    """Write the matrix TSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(rows)


def create_taxa_trait_matrix(
    taxa_path: Path = DEFAULT_TAXA_PATH,
    kg_dir: Path = DEFAULT_KG_DIR,
    output_path: Path = DEFAULT_OUTPUT_PATH,
    column_metadata_path: Path | None = None,
    edges_path: Path | None = None,
    nodes_path: Path | None = None,
    strain_rollup: StrainRollup = PRIMARY_STRAIN_ROLLUP,
    trait_key_policy: TraitKeyPolicy = DEFAULT_TRAIT_KEY_POLICY,
    include_policy_metadata: bool = True,
) -> MatrixStats:
    """Create data/derived/taxa_trait_matrix.tsv from BugSigDB taxa and KG-Microbe.

    ``strain_rollup`` also counts the BacDive strain trait edges its scope admits
    (read_kg_edges); the default is the primary scope, as in the core workflow (#774, #788).
    """
    trait_key_policy = parse_trait_key_policy(trait_key_policy)
    resolved_edges_path = resolve_kg_file(EDGES_FILENAME, kg_dir, edges_path)
    resolved_nodes_path = resolve_kg_file(NODES_FILENAME, kg_dir, nodes_path)
    resolved_column_metadata_path = column_metadata_path or output_path.with_name(
        f"{output_path.stem}_columns{output_path.suffix}"
    )

    if not taxa_path.exists():
        raise FileNotFoundError(f"BugSigDB taxa TSV not found: {taxa_path}")

    reported_taxa = read_reported_taxa(taxa_path)
    child_to_parents, taxon_traits = read_kg_edges(
        resolved_edges_path, strain_rollup=strain_rollup, trait_key_policy=trait_key_policy
    )
    counts = compute_trait_counts(reported_taxa, child_to_parents, taxon_traits)
    trait_keys = {
        trait_key
        for taxon_counts in counts.values()
        for trait_key, count in taxon_counts.items()
        if count > 0
    }
    node_metadata = read_node_metadata(resolved_nodes_path, node_ids_for_trait_keys(trait_keys))
    headers, rows = matrix_rows(
        reported_taxa, counts, node_metadata, trait_key_policy=trait_key_policy
    )
    metadata_headers, metadata_rows = column_metadata_rows(
        counts,
        node_metadata,
        trait_key_policy=trait_key_policy,
        include_policy_metadata=include_policy_metadata,
    )
    write_matrix(output_path, headers, rows)
    write_matrix(resolved_column_metadata_path, metadata_headers, metadata_rows)

    return MatrixStats(
        reported_taxa=len(reported_taxa),
        traits=len(headers) - 3,
        output_path=output_path,
        column_metadata_path=resolved_column_metadata_path,
        edges_path=resolved_edges_path,
        nodes_path=resolved_nodes_path,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a BugSigDB taxa x KG-Microbe trait count matrix. Rows are "
            "unique NCBITaxon taxa reported as increased or decreased in BugSigDB "
            "signatures. Trait columns are METPO phenotypes sorted alphabetically."
        )
    )
    parser.add_argument(
        "--taxa",
        type=Path,
        default=DEFAULT_TAXA_PATH,
        help="Input BugSigDB taxa TSV.",
    )
    parser.add_argument(
        "--kg-dir",
        type=Path,
        default=DEFAULT_KG_DIR,
        help="Directory containing merged-kg_edges.tsv and merged-kg_nodes.tsv.",
    )
    parser.add_argument(
        "--edges",
        type=Path,
        help="Explicit KG-Microbe merged edges TSV path.",
    )
    parser.add_argument(
        "--nodes",
        type=Path,
        help="Explicit KG-Microbe merged nodes TSV path.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help="Output taxa x traits TSV path.",
    )
    parser.add_argument(
        "--column-metadata-output",
        type=Path,
        default=None,
        help=(
            "Output TSV path describing compact trait matrix columns. "
            "Defaults to a path derived from --output."
        ),
    )
    parser.add_argument(
        "--trait-key-policy",
        choices=TRAIT_KEY_POLICIES,
        default=DEFAULT_TRAIT_KEY_POLICY,
        help="Versioned feature mapping: source object IDs or the declared facultative union.",
    )
    parser.add_argument(
        "--strain-rollup",
        choices=STRAIN_ROLLUPS,
        default=PRIMARY_STRAIN_ROLLUP,
        help=(
            "Count BacDive strain (kgmicrobe.strain:) trait edges at the NCBITaxon node "
            "each strain is a subclass of: bare METPO phenotype classes without biosafety "
            "and pathogen labels (phenotype, the primary and default), none (off, the "
            "sensitivity scope) or every edge that does not assert an absence (full, "
            "exploratory). The analyses read the default output as the primary, so off and "
            "full need another --output (#774, #788)."
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    # The analyses read the default path as the primary matrix, so another scope must not
    # replace it (#797, #788)
    secondary = args.strain_rollup != PRIMARY_STRAIN_ROLLUP
    if secondary and args.output.resolve() == DEFAULT_OUTPUT_PATH.resolve():
        parser.error(
            f"--strain-rollup {args.strain_rollup} needs an --output other than the default "
            f"matrix path, which holds the primary {PRIMARY_STRAIN_ROLLUP} matrix"
        )
    stats = create_taxa_trait_matrix(
        taxa_path=args.taxa,
        kg_dir=args.kg_dir,
        output_path=args.output,
        column_metadata_path=args.column_metadata_output,
        edges_path=args.edges,
        nodes_path=args.nodes,
        strain_rollup=args.strain_rollup,
        trait_key_policy=parse_trait_key_policy(args.trait_key_policy),
    )
    print(f"Wrote matrix: {stats.output_path}")
    print(f"Wrote column metadata: {stats.column_metadata_path}")
    print(f"Reported taxa rows: {stats.reported_taxa}")
    print(f"Trait columns: {stats.traits}")
    if args.strain_rollup != "off":
        print(f"BacDive strain trait edges counted at their NCBITaxon nodes: {args.strain_rollup}")
    print(f"KG edges: {stats.edges_path}")
    print(f"KG nodes: {stats.nodes_path}")
    return 0
