"""Check the trait identifiers the figures group against the pinned KG-Microbe nodes.

The oxygen guilds, the gram-stain classes and the six Table 1 traits are fixed
identifiers chosen a priori. Their methods text says the labels were verified
against ``merged-kg_nodes.tsv``; that file is downloaded rather than committed, so
this module makes the claim re-checkable on any machine that has run
``just download-kg`` (#165). It checks that every identifier exists in the node
table with the expected label. It does not check that a grouping is biologically
right: placing microaerophilic and facultatively aerobic in the aerobe guild is an
a priori decision recorded in the guild methods text, not something a label can
confirm.

The strain rollup rule names identifiers too (#788): the handling and clinical classes
the phenotype scope excludes, and every negated METPO relation. Their labels are
checked the same way, and a negated relation the rule does not list fails the check.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from mtp.paths import KG_MICROBE_DIR
from mtp.taxa_trait_matrix import (
    NEGATED_PREDICATES,
    PHENOTYPE_EXCLUDED_TRAITS,
    negative_relation_metadata,
    read_assertion_term_metadata,
    read_metpo_term_metadata,
)
from mtp.taxa_trait_matrix import NEGATION_PREFIXES as NEGATION_PREFIXES
from mtp.trait_guilds import GRAM_GUILDS, OXYGEN_GUILDS, TABLE1_TRAITS

DEFAULT_NODES_PATH = KG_MICROBE_DIR / "merged-kg_nodes.tsv"
DEFAULT_PROVENANCE_PATH = KG_MICROBE_DIR / "provenance.json"
#: Read in blocks so a 200 MB node table is not held in memory to be hashed.
_HASH_BLOCK_BYTES = 1 << 20

#: The label each grouped identifier, and each one the strain rollup rule names (#788),
#: carries in the pinned node table.
EXPECTED_LABELS: dict[str, str] = {
    "METPO:1000602": "aerobic",
    "METPO:1000604": "microaerophilic",
    "METPO:1000608": "facultatively aerobic",
    "METPO:1000603": "anaerobic",
    "METPO:1000607": "obligately anaerobic",
    "METPO:1000698": "gram positive",
    "METPO:1000699": "gram negative",
    "METPO:1000668": "coccus shaped",
    "METPO:1000683": "sphere shaped",
    "METPO:1000623": "moderately halophilic",
    "GO:0030245": "cellulose catabolic process",
    "pathways:arsenate_reduction": "arsenate_reduction",
    **PHENOTYPE_EXCLUDED_TRAITS,
    **NEGATED_PREDICATES,
}
#: How the node table labels a relation that asserts an absence.


def grouped_identifiers() -> set[str]:
    """Every identifier any figure groups, so an unlisted one cannot slip through."""
    groups: Iterable[Mapping[str, Sequence[str]]] = (OXYGEN_GUILDS, GRAM_GUILDS, TABLE1_TRAITS)
    return {term for group in groups for terms in group.values() for term in terms}


@dataclass(frozen=True)
class TermCheck:
    identifier: str
    expected: str
    found: str | None

    @property
    def ok(self) -> bool:
        return self.found == self.expected


@dataclass(frozen=True)
class NodeTableIdentity:
    """Which node table was read, and whether it is the pinned one."""

    path: Path
    declared_md5: str | None
    actual_md5: str
    #: Why the identity could not be established, or ``None`` when it could.
    undetermined: str | None = None

    @property
    def pinned(self) -> bool:
        return self.declared_md5 is not None and self.declared_md5 == self.actual_md5


def file_md5(path: Path) -> str:
    """The MD5 of ``path``, read in blocks."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while block := handle.read(_HASH_BLOCK_BYTES):
            digest.update(block)
    return digest.hexdigest()


def declared_node_md5(provenance_path: Path, filename: str) -> str | None:
    """The MD5 ``provenance.json`` records for ``filename``, or ``None``.

    Returns ``None`` rather than raising when the sidecar is absent or does not
    describe the file: that is "unknown", which the caller reports differently
    from "known and wrong".
    """
    if not provenance_path.is_file():
        return None
    try:
        recorded = json.loads(provenance_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    entry = recorded.get("files", {}).get(filename) if isinstance(recorded, dict) else None
    if not isinstance(entry, dict):
        return None
    declared = entry.get("md5")
    return declared if isinstance(declared, str) else None


def identify_node_table(nodes_path: Path, provenance_path: Path) -> NodeTableIdentity:
    """Establish whether ``nodes_path`` is the KG-Microbe build the sidecar pins.

    A label mismatch means one of two very different things: the pinned build
    renamed a term, or this is simply not the pinned build. Hashing first keeps
    those apart, so a check run against a newer KG cannot be read as evidence
    that a term is missing (#289).
    """
    declared = declared_node_md5(provenance_path, nodes_path.name)
    actual = file_md5(nodes_path)
    undetermined = None if declared is not None else f"{provenance_path} records no md5"
    return NodeTableIdentity(nodes_path, declared, actual, undetermined)


def check_term_labels(nodes_path: Path, expected: Mapping[str, str]) -> list[TermCheck]:
    """Look every expected identifier up in a KG-Microbe nodes TSV."""
    wanted = dict(expected)
    found: dict[str, str] = {}
    with nodes_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_columns = {"id", "name"}.difference(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(f"{nodes_path} lacks columns {sorted(missing_columns)}")
        for row in reader:
            if row["id"] in wanted:
                found[row["id"]] = row["name"]
                if len(found) == len(wanted):
                    break
    return [TermCheck(term, label, found.get(term)) for term, label in sorted(wanted.items())]


def unlisted_negations(nodes_path: Path, edges_path: Path | None = None) -> list[str]:
    """Unreviewed negative relations; an edge snapshot makes the role census authoritative.

    Nodes alone cannot establish assertion roles. That limited census inspects METPO
    terms except explicit phenotype classes, whose absence labels may be affirmed
    phenotype annotations (#914). Observed predicates are checked regardless of category.
    """
    terms = (
        read_assertion_term_metadata(nodes_path, edges_path)
        if edges_path is not None
        else read_metpo_term_metadata(nodes_path)
    )
    return sorted(
        term
        for term, row in terms.items()
        if term not in NEGATED_PREDICATES
        and (edges_path is not None or "PhenotypicQuality" not in row.get("category", ""))
        and negative_relation_metadata(row)
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=Path, default=DEFAULT_NODES_PATH)
    parser.add_argument(
        "--edges",
        type=Path,
        help="Read actual predicate roles for an authoritative negative-relation census.",
    )
    parser.add_argument("--provenance", type=Path, default=DEFAULT_PROVENANCE_PATH)
    parser.add_argument(
        "--allow-unpinned-nodes",
        action="store_true",
        help="Report the labels even when the node table is not the pinned build.",
    )
    args = parser.parse_args(argv)
    unlisted = sorted(grouped_identifiers().difference(EXPECTED_LABELS))
    if unlisted:
        print(f"grouped identifiers without an expected label: {', '.join(unlisted)}")
        return 1

    identity = identify_node_table(args.nodes, args.provenance)
    print(f"node table {identity.path}")
    print(f"  md5 {identity.actual_md5}")
    print(f"  pinned md5 {identity.declared_md5 or '(not recorded)'}")
    if not identity.pinned:
        reason = identity.undetermined or "it does not match the pinned md5"
        print(
            f"this is not the pinned KG-Microbe node table: {reason}. "
            "A label reported as missing here is evidence about this build, not about "
            "the pinned one.",
            file=sys.stderr,
        )
        if not args.allow_unpinned_nodes:
            return 1

    checks = check_term_labels(args.nodes, EXPECTED_LABELS)
    for check in checks:
        status = "ok" if check.ok else "MISMATCH"
        print(
            f"{status:<9} {check.identifier:<30} expected {check.expected!r}, found {check.found!r}"
        )
    failed = [check for check in checks if not check.ok]
    if failed:
        print(
            f"{len(failed)} of {len(checks)} identifiers do not match {args.nodes}", file=sys.stderr
        )
        return 1
    unlisted = unlisted_negations(args.nodes, args.edges)
    if unlisted:
        print(f"negated relations the strain rollup rule does not list: {', '.join(unlisted)}")
        return 1
    print(f"all {len(checks)} identifiers match {args.nodes}")
    if args.edges is None:
        print("Nodes-only negation census: pass --edges to check actual predicate roles.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
