"""Declared trait lists for analyses whose default is the six Table 1 traits (#1054).

A trait-list file holds one trait identifier per line, in report order. Blank lines
and lines starting with ``#`` are ignored. An identifier is a CURIE or a rollup
composite such as ``METPO:2000103|GO:0030245``. Outputs record the file name and its
SHA-256, so each result names the exact list it was computed for.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

_CURIE = r"[A-Za-z][\w.-]*:[^\s|]+"
TRAIT_ID_PATTERN = re.compile(rf"{_CURIE}(?:\|{_CURIE})*")


@dataclass(frozen=True)
class DeclaredTraitList:
    """Trait identifiers in declared order and the identity of the file declaring them."""

    trait_ids: tuple[str, ...]
    source_name: str
    sha256: str

    def provenance(self) -> dict[str, object]:
        """Return the fields every trait-list output records."""
        return {
            "trait_list_file": self.source_name,
            "trait_list_sha256": self.sha256,
            "trait_count": len(self.trait_ids),
            "trait_ids": list(self.trait_ids),
        }


def parse_trait_list(text: str) -> tuple[str, ...]:
    """Parse declared identifiers; reject malformed, repeated or absent entries."""
    trait_ids: list[str] = []
    for number, line in enumerate(text.splitlines(), start=1):
        value = line.strip()
        if not value or value.startswith("#"):
            continue
        if not TRAIT_ID_PATTERN.fullmatch(value):
            raise ValueError(f"Trait list line {number} is not a trait identifier: {value!r}")
        trait_ids.append(value)
    if not trait_ids:
        raise ValueError("Trait list declares no traits")
    repeated = sorted({value for value in trait_ids if trait_ids.count(value) > 1})
    if repeated:
        raise ValueError(f"Trait list repeats {repeated}")
    return tuple(trait_ids)


def load_trait_list(path: Path) -> DeclaredTraitList:
    """Read a trait-list file and hash the exact bytes that were parsed."""
    data = path.read_bytes()
    return DeclaredTraitList(
        trait_ids=parse_trait_list(data.decode("utf-8")),
        source_name=path.name,
        sha256=hashlib.sha256(data).hexdigest(),
    )


def trait_labels(columns: pd.DataFrame, trait_ids: Sequence[str]) -> dict[str, str]:
    """Label declared traits from a taxa_trait_matrix_columns table's ``object_name``.

    Every declared trait must be a column of the table and carry a label; a declared
    list never falls back to an identifier or a placeholder.
    """
    if not {"column_id", "object_name"} <= set(columns.columns):
        raise ValueError("Trait columns table needs column_id and object_name")
    names = columns.drop_duplicates("column_id").set_index("column_id")["object_name"]
    missing = [trait_id for trait_id in trait_ids if trait_id not in names.index]
    if missing:
        raise ValueError(f"Declared traits absent from the trait columns table: {missing}")
    labels = {trait_id: names[trait_id] for trait_id in trait_ids}
    unlabeled = [
        trait_id
        for trait_id, label in labels.items()
        if not isinstance(label, str) or not label.strip()
    ]
    if unlabeled:
        raise ValueError(f"Declared traits without an object_name label: {unlabeled}")
    return {trait_id: str(label).strip() for trait_id, label in labels.items()}


def file_sha256(path: Path) -> str:
    """Hash one input file for trait-list provenance."""
    return hashlib.sha256(path.read_bytes()).hexdigest()
