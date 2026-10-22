"""Build a BugSigDB taxon by KG-Microbe trait count matrix."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from mtp.paths import KG_MICROBE_DIR, REPO_ROOT

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
    ``biolink:subclass_of``, and ``strain_traits`` holds the strain's own trait keys
    (#774). Both are empty unless strains were read.
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


def make_trait_key(predicate: str, obj: str) -> str:
    """Create the output trait key for a KG taxon assertion."""
    if _is_trait(obj):
        return obj
    return f"{predicate}{TRAIT_KEY_SEPARATOR}{obj}"


def split_trait_key(trait_key: str) -> tuple[str | None, str]:
    """Return optional predicate and object IDs encoded in a trait key."""
    if TRAIT_KEY_SEPARATOR not in trait_key:
        return None, trait_key
    predicate, obj = trait_key.split(TRAIT_KEY_SEPARATOR, 1)
    return predicate, obj


def node_ids_for_trait_keys(trait_keys: set[str]) -> set[str]:
    """Return KG node IDs needed to label trait keys."""
    node_ids: set[str] = set()
    for trait_key in trait_keys:
        predicate, obj = split_trait_key(trait_key)
        node_ids.add(obj)
        if predicate is not None:
            node_ids.add(predicate)
    return node_ids


def read_kg_trait_edges(edges_path: Path, *, strains: bool = False) -> KgTraitEdges:
    """Read taxonomy parent links and direct trait assertions from KG edges.

    Trait sets are keyed by NCBITaxon subject. With ``strains``, the trait edges of
    BacDive strain subjects (``kgmicrobe.strain:``) and the strains' NCBITaxon
    ``subclass_of`` links are also kept, apart from the taxon sets (#774).

    Memory: materializes the full taxonomy graph plus per-taxon trait sets in
    Python dicts. On the current ~1 GB merged-kg_edges.tsv this peaks around
    2-3 GB of resident memory, and around 6 GB with strains; switch to a
    streaming / two-pass implementation if that becomes a constraint.
    """
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
                taxon_traits[subject].add(make_trait_key(predicate, obj))
            elif strains and _is_strain(subject):
                if predicate == SUBCLASS_OF and _is_taxon(obj):
                    strain_taxa[subject].add(obj)
                elif _is_trait(obj) or _is_trait(predicate):
                    strain_traits[subject].add(make_trait_key(predicate, obj))

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


def read_kg_edges(
    edges_path: Path, *, strain_rollup: bool = False
) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Read taxonomy parent links and per-taxon trait sets from KG edges.

    Only NCBITaxon subjects carry traits by default, so BacDive strain records never
    reach the matrix. ``strain_rollup`` adds them at the strains' NCBITaxon nodes through
    roll_up_strain_traits; it stays off by default because it adds traits to the
    vocabulary of every downstream analysis (#774).
    """
    edges = read_kg_trait_edges(edges_path, strains=strain_rollup)
    if not strain_rollup:
        return edges.child_to_parents, edges.taxon_traits
    return edges.child_to_parents, roll_up_strain_traits(
        edges.taxon_traits, edges.strain_taxa, edges.strain_traits
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


def _node_label(node_id: str, node_metadata: dict[str, NodeMetadata]) -> str:
    node_name = node_metadata.get(node_id, NodeMetadata("", "")).name
    if node_name:
        return f"{node_name} {node_id}"
    return node_id


def trait_column_label(trait_key: str, node_metadata: dict[str, NodeMetadata]) -> str:
    """Return the human-readable label for a KG trait."""
    predicate, obj = split_trait_key(trait_key)
    if predicate is None:
        return _node_label(obj, node_metadata)

    predicate_name = node_metadata.get(predicate, NodeMetadata(predicate, "")).name or predicate
    return f"{predicate_name}: {_node_label(obj, node_metadata)} [{predicate}]"


def matrix_rows(
    reported_taxa: list[ReportedTaxon],
    counts: dict[str, Counter[str]],
    node_metadata: dict[str, NodeMetadata],
) -> tuple[list[str], list[list[str]]]:
    """Build TSV headers and rows for the wide taxa-by-traits matrix."""

    def trait_sort_key(trait_id: str) -> tuple[str, str]:
        trait_label = trait_column_label(trait_id, node_metadata)
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
) -> tuple[list[str], list[list[str]]]:
    """Build metadata rows describing compact matrix columns."""

    def trait_sort_key(trait_id: str) -> tuple[str, str]:
        trait_label = trait_column_label(trait_id, node_metadata)
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
                trait_column_label(trait_id, node_metadata),
                predicate or "",
                predicate_metadata.name if predicate_metadata else "",
                obj,
                object_metadata.name,
                object_metadata.category,
            ]
        )
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
    strain_rollup: bool = False,
) -> MatrixStats:
    """Create data/derived/taxa_trait_matrix.tsv from BugSigDB taxa and KG-Microbe.

    ``strain_rollup`` also counts BacDive strain trait edges (read_kg_edges, #774).
    """
    resolved_edges_path = resolve_kg_file(EDGES_FILENAME, kg_dir, edges_path)
    resolved_nodes_path = resolve_kg_file(NODES_FILENAME, kg_dir, nodes_path)
    resolved_column_metadata_path = column_metadata_path or output_path.with_name(
        f"{output_path.stem}_columns{output_path.suffix}"
    )

    if not taxa_path.exists():
        raise FileNotFoundError(f"BugSigDB taxa TSV not found: {taxa_path}")

    reported_taxa = read_reported_taxa(taxa_path)
    child_to_parents, taxon_traits = read_kg_edges(resolved_edges_path, strain_rollup=strain_rollup)
    counts = compute_trait_counts(reported_taxa, child_to_parents, taxon_traits)
    trait_keys = {
        trait_key
        for taxon_counts in counts.values()
        for trait_key, count in taxon_counts.items()
        if count > 0
    }
    node_metadata = read_node_metadata(resolved_nodes_path, node_ids_for_trait_keys(trait_keys))
    headers, rows = matrix_rows(reported_taxa, counts, node_metadata)
    metadata_headers, metadata_rows = column_metadata_rows(counts, node_metadata)
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
        "--strain-rollup",
        action="store_true",
        help=(
            "Also count the trait edges of BacDive strains (kgmicrobe.strain:) at the "
            "NCBITaxon node each strain is a subclass of. Off by default: it adds traits "
            "to the matrix vocabulary (#774). Needs an --output other than the default."
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    # The analyses read the default path, so a rolled-up matrix must not replace it (#797)
    if args.strain_rollup and args.output.resolve() == DEFAULT_OUTPUT_PATH.resolve():
        parser.error("--strain-rollup needs an --output other than the default matrix path")
    stats = create_taxa_trait_matrix(
        taxa_path=args.taxa,
        kg_dir=args.kg_dir,
        output_path=args.output,
        column_metadata_path=args.column_metadata_output,
        edges_path=args.edges,
        nodes_path=args.nodes,
        strain_rollup=args.strain_rollup,
    )
    print(f"Wrote matrix: {stats.output_path}")
    print(f"Wrote column metadata: {stats.column_metadata_path}")
    print(f"Reported taxa rows: {stats.reported_taxa}")
    print(f"Trait columns: {stats.traits}")
    if args.strain_rollup:
        print("BacDive strain trait edges counted at their NCBITaxon nodes (--strain-rollup)")
    print(f"KG edges: {stats.edges_path}")
    print(f"KG nodes: {stats.nodes_path}")
    return 0
