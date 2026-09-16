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
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from bsdbng.paths import KG_MICROBE_DIR
from bsdbng.trait_guilds import GRAM_GUILDS, OXYGEN_GUILDS, TABLE1_TRAITS

DEFAULT_NODES_PATH = KG_MICROBE_DIR / "merged-kg_nodes.tsv"

#: The label each grouped identifier carries in the pinned node table.
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
}


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


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=Path, default=DEFAULT_NODES_PATH)
    args = parser.parse_args(argv)
    unlisted = sorted(grouped_identifiers().difference(EXPECTED_LABELS))
    if unlisted:
        print(f"grouped identifiers without an expected label: {', '.join(unlisted)}")
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
    print(f"all {len(checks)} identifiers match {args.nodes}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
