"""Import bounded, reconciled carrier evidence without inventing source records (#1044, #1055)."""

from __future__ import annotations

import csv
import hashlib
import json
import sqlite3
import tempfile
from collections import Counter
from collections.abc import Iterator
from contextlib import closing
from pathlib import Path
from typing import Any

import pandas as pd

from mtp.count_input import parse_nonnegative_integer
from mtp.label_support_audit import AUDIT_TRAITS
from mtp.multiplicity.families import rows
from mtp.multiplicity.sources import VerifiedSources
from mtp.paths import file_sha256

CARRIER_POPULATION = "reported_taxa_observed_in_supplied_count_profiles"
SCOPE_DESCRIPTION = "exact supplied profile keys; no independent host or orientation claim"
Cell = tuple[str, str]


def _empty(traits: tuple[str, ...], status: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "trait": traits,
            "carrier_status": status,
            "carrier_population": CARRIER_POPULATION,
            "carrier_population_note": SCOPE_DESCRIPTION,
            "carrier_catalog_coverage": "focal_and_selected_oxygen_traits_only",
            "carrier_direct_reported_taxa": None,
            "carrier_inherited_reported_taxa": None,
            "carrier_unique_taxa": None,
            "carrier_source_record_count": None,
            "carrier_source_record_status": "unresolved_original_record_identifiers_unavailable",
            "carrier_madin_status": "not_supplied",
        }
    )


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict) or any(not isinstance(key, str) for key in value):
        raise ValueError(f"Malformed label-support {label}")
    return value


def _checked_rows(path: Path, digest: str | None, required: set[str]) -> Iterator[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    if not required <= set(header):
        raise ValueError(f"Label-support {path.name} lacks required columns")
    yield from rows(path, digest)


def _count(row: dict[str, str], column: str, path: Path) -> int:
    if column not in row:
        raise ValueError(f"Label-support {path.name} lacks required columns: {column}")
    return parse_nonnegative_integer(row[column], path=path, line_number=0, column_name=column)


def _valid_taxon(value: str) -> bool:
    prefix, separator, suffix = value.partition(":")
    return prefix == "NCBITaxon" and bool(separator) and suffix.isascii() and suffix.isdecimal()


def _core_inventory(
    sources: VerifiedSources,
) -> tuple[dict[Cell, int], Counter[str], dict[str, Any]]:
    """Derive positive audit cells and the exact supplied-profile occurrence population."""
    matrix = sources.product("taxa_matrix")
    expected: dict[Cell, int] = {}
    matrix_taxa: set[str] = set()
    for row in _checked_rows(
        matrix, sources.expected_digest(matrix), {"taxon_id", *sources.features["count_profiles"]}
    ):
        taxon = row["taxon_id"]
        if not _valid_taxon(taxon) or taxon in matrix_taxa:
            raise ValueError("Invalid or duplicate label-support matrix taxon")
        matrix_taxa.add(taxon)
        for trait in AUDIT_TRAITS.keys() & row.keys():
            value = _count(row, trait, matrix)
            if value:
                expected[(taxon, trait)] = value
    if not matrix_taxa:
        raise ValueError("Label-support matrix contains no taxon rows")
    profile = sources.product("count_profiles")
    keys: set[tuple[str, str, str]] = set()
    for row in _checked_rows(
        profile, sources.expected_digest(profile), {"study_id", "experiment_id", "inc_dec"}
    ):
        key = row["study_id"], row["experiment_id"], row["inc_dec"]
        if not all(key) or key[2] not in {"inc", "dec"} or key in keys:
            raise ValueError("Label-support scope has duplicate or invalid profile identifiers")
        keys.add(key)
    occurrences = sources.product("reported_taxa")
    bound: set[tuple[str, str, str, str]] = set()
    input_rows = 0
    for row in _checked_rows(
        occurrences,
        sources.expected_digest(occurrences),
        {"study_id", "experiment_id", "direction", "taxon_id"},
    ):
        input_rows += 1
        direction = {"increased": "inc", "decreased": "dec"}.get(row["direction"].strip().lower())
        if direction is None:
            continue
        key = row["study_id"], row["experiment_id"], direction
        if key in keys:
            if not _valid_taxon(row["taxon_id"]):
                raise ValueError("Label-support scope requires explicit NCBITaxon identifiers")
            bound.add((*key, row["taxon_id"]))
    if {row[:3] for row in bound} != keys:
        raise ValueError("Label-support scope is missing a supplied profile side")
    scope = {
        "input_occurrence_rows": input_rows,
        "bound_distinct_occurrences": len(bound),
        "bound_profiles": len(keys),
        "profile_counts_reconciled": True,
        "unmapped_occurrences": sum(row[3] not in matrix_taxa for row in bound),
        "population": SCOPE_DESCRIPTION,
    }
    return expected, Counter(row[3] for row in bound), scope


def _validate_support(
    path: Path, digest: str, expected: dict[Cell, int], occurrences: Counter[str]
) -> dict[Cell, tuple[bool, int]]:
    support: dict[Cell, tuple[bool, int]] = {}
    for row in _checked_rows(path, digest, {"taxon_id", "trait_key", "carrier_count"}):
        key = row["taxon_id"], row["trait_key"]
        if key in support or key not in expected:
            raise ValueError("Duplicate or unsupported label-support audit cell")
        count = _count(row, "carrier_count", path)
        inherited = _count(row, "inherited_carrier_count", path)
        observed = _count(row, "signature_occurrences", path)
        direct_text = row.get("carrier_at_reported_taxon", "").lower()
        if direct_text not in {"true", "false"}:
            raise ValueError("Invalid direct-carrier boolean")
        direct = direct_text == "true"
        if count != expected[key] or count != inherited + int(direct):
            raise ValueError("Label-support carrier counts differ from the core matrix")
        if observed != occurrences[key[0]]:
            raise ValueError("Label-support signature occurrence count differs from bound scope")
        support[key] = direct, inherited
    if set(support) != set(expected):
        raise ValueError("Label-support audit is missing positive core matrix cells")
    return support


def _validate_evidence(
    path: Path,
    digest: str,
    expected: dict[Cell, int],
    support: dict[Cell, tuple[bool, int]],
    occurrences: Counter[str],
) -> dict[str, int]:
    """Check unique evidence per cell with bounded memory; scratch follows durable TMPDIR."""
    with (
        tempfile.TemporaryDirectory(prefix="mtp-carriers-") as scratch,
        closing(sqlite3.connect(Path(scratch) / "carriers.sqlite")) as db,
        db,
    ):
        db.execute(
            "CREATE TABLE evidence (taxon TEXT, trait TEXT, carrier TEXT, selected INTEGER, "
            "PRIMARY KEY (taxon, trait, carrier)) WITHOUT ROWID"
        )
        for row in _checked_rows(
            path, digest, {"reported_taxon_id", "trait_key", "carrier_taxon_id"}
        ):
            key = row["reported_taxon_id"], row["trait_key"]
            carrier = row["carrier_taxon_id"]
            if key not in expected or not _valid_taxon(carrier):
                raise ValueError("Orphan or malformed label-support carrier evidence")
            try:
                db.execute(
                    "INSERT INTO evidence VALUES (?, ?, ?, ?)",
                    (*key, carrier, int(occurrences[key[0]] > 0)),
                )
            except sqlite3.IntegrityError as error:
                raise ValueError("Duplicate label-support carrier evidence") from error
        seen: set[Cell] = set()
        for taxon, trait, count, direct in db.execute(
            "SELECT taxon, trait, COUNT(*), MAX(taxon = carrier) "
            "FROM evidence GROUP BY taxon, trait"
        ):
            key = taxon, trait
            if count != expected[key] or bool(direct) != support[key][0]:
                raise ValueError("Label-support carrier evidence disagrees with support counts")
            seen.add(key)
        if seen != set(expected):
            raise ValueError("Label-support carrier evidence is missing positive cells")
        return dict(
            db.execute(
                "SELECT trait, COUNT(DISTINCT carrier) FROM evidence "
                "WHERE selected = 1 GROUP BY trait"
            )
        )


def carrier_table(
    sources: VerifiedSources, audit: Path | None
) -> tuple[pd.DataFrame, tuple[Path, ...]]:
    """Bind audit bytes and reconcile full focal coverage before reporting any observed zeros."""
    features = sources.features["count_profiles"]
    if audit is None:
        return _empty(features, "not_supplied"), ()
    manifest = audit / "provenance.json" if audit.is_dir() else audit
    root = manifest.parent.resolve()
    encoded = manifest.read_bytes()
    manifest_digest = hashlib.sha256(encoded).hexdigest()
    record = _mapping(json.loads(encoded), "provenance")
    if (
        record.get("schema_version") != "1.0"
        or record.get("analysis") != "label_support_sidecar"
        or record.get("strain_rollup") != sources.scope
        or record.get("verification_status")
        != "hashes_and_counts_reconciled_not_core_execution_verified"
        or record.get("support_table_population") != "all_supplied_matrix_taxa"
    ):
        raise ValueError("Unsupported label-support audit provenance or scope")
    inputs = _mapping(record.get("inputs"), "inputs")
    for role, product in (
        ("profiles", "count_profiles"),
        ("matrix", "taxa_matrix"),
        ("taxa", "reported_taxa"),
        ("kg_edges", "kg_edges"),
        ("kg_nodes", "kg_nodes"),
    ):
        expected_hash = sources.artifact(sources.product(product)).sha256
        if _mapping(inputs.get(role), f"input {role}").get("sha256") != expected_hash:
            raise ValueError(f"Label-support audit input differs from verified core: {role}")
    outputs = _mapping(record.get("outputs"), "outputs")
    required = {"taxon_trait_label_support.tsv", "carrier_source_evidence.tsv"}
    if not required <= set(outputs):
        raise ValueError("Label-support audit lacks required output inventory")
    paths = [manifest]
    for name, value in outputs.items():
        metadata = _mapping(value, "output inventory")
        path = root / name
        if (
            Path(name).name != name
            or path.is_symlink()
            or not path.resolve().is_relative_to(root)
            or not path.is_file()
            or file_sha256(path) != metadata.get("sha256")
            or path.stat().st_size != metadata.get("size_bytes")
        ):
            raise ValueError("Label-support audit output inventory mismatch")
        paths.append(path)
    expected, occurrences, scope = _core_inventory(sources)
    if record.get("scope") != scope:
        raise ValueError("Label-support audit scope facts differ from core inputs")
    support_path = root / "taxon_trait_label_support.tsv"
    support = _validate_support(
        support_path, outputs[support_path.name]["sha256"], expected, occurrences
    )
    if len(support) != record.get("support_cells"):
        raise ValueError("Label-support audit cell count differs from inventory")
    evidence_path = root / "carrier_source_evidence.tsv"
    carriers = _validate_evidence(
        evidence_path, outputs[evidence_path.name]["sha256"], expected, support, occurrences
    )
    table = _empty(features, "not_in_focal_audit")
    for index, feature in enumerate(features):
        if feature not in AUDIT_TRAITS:
            continue
        cells = [
            value
            for (taxon, trait), value in support.items()
            if trait == feature and occurrences[taxon] > 0
        ]
        table.loc[index, "carrier_status"] = "audit_hash_bound_focal_coverage"
        table.loc[index, "carrier_direct_reported_taxa"] = sum(direct for direct, _ in cells)
        table.loc[index, "carrier_inherited_reported_taxa"] = sum(
            inherited > 0 for _, inherited in cells
        )
        table.loc[index, "carrier_unique_taxa"] = carriers.get(feature, 0)
    table["carrier_madin_status"] = record.get("madin_status", "unresolved")
    if file_sha256(manifest) != manifest_digest:
        raise ValueError("Label-support provenance changed while reading the audit")
    return table, tuple(sorted(paths))
