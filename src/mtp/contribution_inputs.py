"""Strict, provenance-preserving inputs for the Table 1 contribution audit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

from mtp.paths import REPO_ROOT
from mtp.taxa_trait_matrix import (
    PRIMARY_STRAIN_ROLLUP,
    STRAIN_PREFIX,
    StrainRollup,
    make_trait_key,
    scoped_taxon_traits,
    strain_trait_rolls_up,
)
from mtp.trait_guilds import MEDIADIVE_MARKER, TABLE1_TRAITS, column_term_id
from mtp.trait_key_policy import FACULTATIVE_MEMBER_IDS, FACULTATIVE_TRAIT_ID, TraitKeyPolicy

Key = tuple[str, str, str]
Pair = tuple[str, str]
RANKS = ("species", "genus", "family")


def digest(path: Path, algorithm: str = "sha256") -> str:
    result = hashlib.new(algorithm, usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def rows(path: Path) -> Iterator[dict[str, str]]:
    """Reject duplicate columns and malformed records instead of padding counts."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if not header or len(set(header)) != len(header):
            raise ValueError(f"Empty or duplicate header: {path}")
        for number, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(f"Wrong field count: {path}:{number}")
            yield dict(zip(header, row, strict=True))


def count(value: str) -> int:
    if not value or not value.isascii() or not value.isdecimal():
        raise ValueError(f"Expected nonnegative integer count, got {value!r}")
    return int(value)


@dataclass
class Base:
    columns: list[str]
    traits: dict[str, str]
    profiles: dict[Key, dict[str, int]]
    metadata: dict[Pair, dict[str, str]]
    reported: dict[Key, dict[str, set[str]]]
    taxa: dict[str, dict[str, int]]
    names: dict[str, str]
    ranks: dict[str, str]
    eligibility: list[dict[str, str]]


def load_base(
    profile_path: Path, experiment_path: Path, taxa_path: Path, matrix_path: Path
) -> Base:
    metadata: dict[Pair, dict[str, str]] = {}
    for row in rows(experiment_path):
        for required in ("condition", "body_site", "host_species", "study_id", "experiment_id"):
            if required not in row:
                raise ValueError(f"Missing experiment field {required}")
        key = (row["study_id"], row["experiment_id"])
        if key in metadata:
            raise ValueError(f"Duplicate experiment metadata: {key}")
        metadata[key] = row

    profiles: dict[Key, dict[str, int]] = {}
    columns: list[str] = []
    for row in rows(profile_path):
        key3 = (row.pop("study_id"), row.pop("experiment_id"), row.pop("inc_dec"))
        if not columns:
            columns = [c for c in row if MEDIADIVE_MARKER not in c.lower()]
        values = {c: count(v) for c, v in row.items()}
        if key3 in profiles or key3[2] not in {"inc", "dec"}:
            raise ValueError(f"Duplicate/invalid profile key: {key3}")
        profiles[key3] = {c: values[c] for c in columns}
    if not columns:
        raise ValueError("No profile trait columns")
    traits: dict[str, str] = {}
    for slug, terms in TABLE1_TRAITS.items():
        matches = [c for c in columns if c in terms or column_term_id(c) in terms]
        if len(matches) != 1:
            raise ValueError(f"Trait {slug} resolves to {matches}")
        traits[slug] = matches[0]

    eligibility = []
    retained: dict[Key, dict[str, int]] = {}
    for key3, values in profiles.items():
        pair = key3[:2]
        meta = metadata.get(pair)
        reason = (
            "missing_metadata"
            if meta is None
            else "non_human"
            if meta["host_species"].strip() != "Homo sapiens"
            else "unpaired"
            if (*pair, "inc") not in profiles or (*pair, "dec") not in profiles
            else "retained"
        )
        eligibility.append(
            {
                "study_id": pair[0],
                "experiment_id": pair[1],
                "direction": key3[2],
                "status": reason,
                "orientation_status": "unreviewed_188",
            }
        )
        if reason == "retained":
            retained[key3] = values

    taxa: dict[str, dict[str, int]] = {}
    names: dict[str, str] = {}
    ranks: dict[str, str] = {}
    for row in rows(matrix_path):
        taxon = row.pop("taxon_id")
        names[taxon] = row.pop("taxon_name")
        ranks[taxon] = row.pop("taxonomic_rank").strip().lower()
        if taxon in taxa or not set(columns).issubset(row):
            raise ValueError(f"Duplicate taxon or missing matrix trait columns: {taxon}")
        values = {c: count(v) for c, v in row.items()}
        taxa[taxon] = {c: values[c] for c in columns if values[c]}

    reported: dict[Key, dict[str, set[str]]] = {k: {} for k in retained}
    for row in rows(taxa_path):
        direction = {"increased": "inc", "decreased": "dec"}.get(row["direction"])
        if direction is None:
            raise ValueError(f"Invalid reported direction: {row['direction']}")
        key3 = (row["study_id"], row["experiment_id"], direction)
        taxon = row["taxon_id"]
        rank = row["taxonomic_rank"].strip().lower()
        if taxon in ranks and rank and ranks[taxon] and ranks[taxon] != rank:
            raise ValueError(f"Conflicting ranks for {taxon}")
        ranks[taxon] = rank or ranks.get(taxon, "")
        names.setdefault(taxon, row["taxon_name"])
        if key3 in reported:
            reported[key3].setdefault(taxon, set()).add(row["signature_id"])
    return Base(columns, traits, retained, metadata, reported, taxa, names, ranks, eligibility)


@dataclass
class Assertion:
    edge_id: str
    predicate: str
    object_id: str
    primary_source: str
    knowledge_source: str
    knowledge_level: str
    source_row: int


@dataclass
class Graph:
    parents: dict[str, set[str]]
    traits: dict[str, set[str]]
    assertions: dict[tuple[str, str], list[Assertion]]
    names: dict[str, str]
    header_repaired: bool
    continued_records: int = 0
    ancestor_cache: dict[str, frozenset[str]] = field(default_factory=dict)

    def ancestors(self, taxon: str) -> frozenset[str]:
        """All ancestors including self; reject cycles, retain multiple parents."""
        cached = self.ancestor_cache.get(taxon)
        if cached is not None:
            return cached
        active: set[str] = set()

        def visit(node: str) -> frozenset[str]:
            if node in active:
                raise ValueError(f"Taxonomy cycle at {node}")
            if node in self.ancestor_cache:
                return self.ancestor_cache[node]
            active.add(node)
            found = {node}
            for parent in sorted(self.parents.get(node, set())):
                found.update(visit(parent))
            active.remove(node)
            value = frozenset(found)
            self.ancestor_cache[node] = value
            return value

        return visit(taxon)

    def path(self, source: str, target: str) -> tuple[list[str], bool]:
        """A deterministic path plus whether another valid parent branch exists."""
        if source == target:
            return [source], False
        candidates = sorted(
            p for p in self.parents.get(source, set()) if target in self.ancestors(p)
        )
        if not candidates:
            raise ValueError(f"No rollup path from {source} to {target}")
        tail, multiple = self.path(candidates[0], target)
        return [source, *tail], multiple or len(candidates) > 1

    def ranked(self, taxon: str, rank: str, ranks: dict[str, str]) -> tuple[str, str]:
        found = sorted(a for a in self.ancestors(taxon) if ranks.get(a) == rank)
        if len(found) == 1:
            return found[0], "resolved_from_reported_rank_and_pinned_taxonomy"
        return "", "ambiguous_rank" if found else "unresolved_rank"


def load_graph(
    edges: Path,
    nodes: Path,
    columns: set[str],
    six: set[str],
    strain_rollup: StrainRollup = PRIMARY_STRAIN_ROLLUP,
) -> Graph:
    """Read the taxonomy, retained traits and six-trait edge provenance of the KG.

    Under the matrix's strain rollup scope, the BacDive strain records it admits count
    at the strains' linked NCBITaxon nodes by the builder's rule, and their edges are
    assertions of those nodes. The default is the primary scope, which the legacy matrix
    builder also defaults to (#788).
    """
    parents: defaultdict[str, set[str]] = defaultdict(set)
    traits: defaultdict[str, set[str]] = defaultdict(set)
    assertions: defaultdict[tuple[str, str], list[Assertion]] = defaultdict(list)
    strain_taxa: defaultdict[str, set[str]] = defaultdict(set)
    strain_traits: defaultdict[str, set[str]] = defaultdict(set)
    strain_assertions: defaultdict[tuple[str, str], list[Assertion]] = defaultdict(list)
    repaired = False
    continued = 0
    policy: TraitKeyPolicy = (
        "raw_object_v1"
        if any(m in columns for m in FACULTATIVE_MEMBER_IDS if m != FACULTATIVE_TRAIT_ID)
        else "facultative_union_v1"
    )
    with edges.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:4] != ["id", "subject", "predicate", "object"]:
            raise ValueError("Unexpected KG edge header")
        for row in reader:
            if row[:2] == ["", "knowledge_level"] and header[-1] == "agent_type":
                # The pinned 20260120 export has a physical newline in its header.
                if repaired:
                    raise ValueError("Repeated KG header continuation")
                header += row[1:]
                repaired = True
                continue
            source_row = reader.line_num
            if len(row) == 7 and row[6] == "manual_validation_of_automated_agent":
                continuation = next(reader, [])
                if len(continuation) != 5 or continuation[0] != "":
                    raise ValueError(f"Invalid agent-field continuation at line {source_row}")
                row += continuation[1:]
                continued += 1
            if len(row) != len(header):
                raise ValueError(f"Malformed KG record at line {reader.line_num}")
            edge, subject, predicate, obj = row[:4]
            strain = strain_rollup != "off" and subject.startswith(STRAIN_PREFIX)
            if not (strain or subject.startswith("NCBITaxon:")):
                continue
            if predicate == "biolink:subclass_of" and obj.startswith("NCBITaxon:"):
                (strain_taxa if strain else parents)[subject].add(obj)
            elif predicate.startswith("METPO:") or obj.startswith("METPO:"):
                key = make_trait_key(predicate, obj, trait_key_policy=policy)
                if key not in columns or (strain and not strain_trait_rolls_up(key, strain_rollup)):
                    continue
                (strain_traits if strain else traits)[subject].add(key)
                if key in six:
                    meta = dict(zip(header, row, strict=True))
                    (strain_assertions if strain else assertions)[(subject, key)].append(
                        Assertion(
                            edge,
                            predicate,
                            obj,
                            meta.get("primary_knowledge_source", ""),
                            meta.get("knowledge_source", ""),
                            meta.get("knowledge_level", ""),
                            source_row,
                        )
                    )
    for (source, key), evidence in strain_assertions.items():
        for taxon in sorted(strain_taxa.get(source, ())):
            assertions[(taxon, key)].extend(evidence)
    scoped = scoped_taxon_traits(
        dict(traits), dict(strain_taxa), dict(strain_traits), strain_rollup
    )
    names = {r["id"]: r["name"] for r in rows(nodes) if r["id"].startswith("NCBITaxon:")}
    return Graph(dict(parents), scoped, dict(assertions), names, repaired, continued)


def verify_kg(kg_dir: Path, expected: dict[str, str]) -> dict[str, str]:
    hashes = {}
    for filename, wanted in expected.items():
        path = kg_dir / filename
        actual = digest(path, "md5")
        if actual != wanted:
            raise ValueError(f"Pinned KG mismatch for {path}: {actual} != {wanted}")
        hashes[filename] = digest(path)
    return hashes


def encoded(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def main(argv: list[str] | None = None) -> None:
    """Check data/kg-microbe against the md5s pinned in the Table 1 influence protocol."""
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("--kg-dir", type=Path, default=REPO_ROOT / "data/kg-microbe")
    parser.add_argument(
        "--protocol", type=Path, default=REPO_ROOT / "methods/table1_taxon_influence_protocol.json"
    )
    args = parser.parse_args(argv)
    expected = json.loads(args.protocol.read_text(encoding="utf-8"))["kg_md5"]
    missing = sorted(name for name in expected if not (args.kg_dir / name).is_file())
    if missing:
        msg = (
            f"{args.kg_dir} lacks the pinned KG files {', '.join(missing)}; run `just download-kg`"
        )
        raise SystemExit(msg)
    try:
        verify_kg(args.kg_dir, expected)
    except ValueError as error:
        # `just download-kg` keeps a local file whose size matches, so it would not refetch it.
        msg = f"{error}; refetch the KG with `just download-kg-force`"
        raise SystemExit(msg) from error
    print(f"{', '.join(sorted(expected))} in {args.kg_dir} match the md5s in {args.protocol.name}")


__all__ = [
    "RANKS",
    "Assertion",
    "Base",
    "Graph",
    "Key",
    "Pair",
    "count",
    "digest",
    "encoded",
    "load_base",
    "load_graph",
    "main",
    "rows",
    "verify_kg",
]


if __name__ == "__main__":
    main()
