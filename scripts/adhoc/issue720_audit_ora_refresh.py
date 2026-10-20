"""Audit regenerated complete ORA exports for #720 (2026-10-19).

Inputs: a JSON manifest of historical/current TSV pairs and their explicit profile
and taxa-matrix inputs. Output: a JSON audit; nonzero exit on unexpected changes.
Only new input fingerprint metadata is allowed. Historical exports are read-only.
Input eligibility, taxon/background counts and profile/trait completeness are
checked independently (#725); p/q values and effect estimates are not recalculated.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from itertools import chain
from pathlib import Path

from mtp.count_input import PROFILE_KEY_COLUMNS, open_input_text
from mtp.enrichment_results import (
    COUNT_COLUMNS,
    FLOAT_COLUMNS,
    ID_COLUMNS,
    META_RESULT_STATUS,
    META_SIGNIFICANT_HITS,
    META_TESTS_PERFORMED,
    PROBABILITY_COLUMNS,
    RESULT_COLUMNS,
    classify_result_status,
)
from mtp.paths import file_sha256
from mtp.trait_guilds import MEDIADIVE_MARKER

FINGERPRINTS = {"Profiles-SHA256": "profiles", "Taxa-Matrix-SHA256": "taxa_matrix"}
PAIR_FIELDS = {"label", "historical", "current", "profiles", "taxa_matrix"}
BACKGROUND_COLUMNS = ("background_m", "background_k")
COMPLETE_COLUMNS = {
    *RESULT_COLUMNS,
    *BACKGROUND_COLUMNS,
    "is_significant",
    "correction_family",
    "eligibility_status",
}
ProfileKey = tuple[str, str, str]


def table_rows(path: Path) -> Iterator[tuple[list[str], list[str]]]:
    """Stream input rows without accepting duplicate columns or ragged records."""
    with open_input_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        header = next(reader, [])
        if not header or len(header) != len(set(header)) or any(not col for col in header):
            raise ValueError(f"{path}: missing or duplicate input columns")
        for row in reader:
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(f"{path}: malformed input row {reader.line_num}")
            yield header, row


@dataclass(frozen=True)
class InputPopulation:
    """Input-defined profile/trait family; no result rows or statistics are retained."""

    profile_indices: dict[ProfileKey, int]
    paired_profiles: int
    background_m: int
    background_k: dict[str, int]

    def check_metadata(self, metadata: dict[str, str]) -> None:
        for key, expected in (
            ("Total-Profiles-Input", self.paired_profiles),
            ("Testable-Profiles", len(self.profile_indices)),
            ("Background-Taxa-Universe-M", self.background_m),
        ):
            if number(metadata[key], key, integer=True) != expected:
                raise ValueError(f"{key} differs from input population")
        if not metadata.get("Background-Mode", "").startswith("TaxaMatrix"):
            raise ValueError("Background-Mode must declare the supplied taxa matrix")

    def observe(self, row: list[str], indices: dict[str, int], masks: dict[str, int]) -> None:
        key = (row[indices["study_id"]], row[indices["experiment_id"]], row[indices["inc_dec"]])
        if key not in self.profile_indices:
            raise ValueError(f"profile absent from eligible inputs: {key}")
        trait = row[indices["trait"]]
        if trait not in self.background_k:
            raise ValueError(f"trait absent from eligible inputs: {trait}")
        background = number(row[indices["background_k"]], "background_k", integer=True)
        if background != self.background_k[trait]:
            raise ValueError(f"background_k differs from input matrix: {trait}")
        bit = 1 << self.profile_indices[key]
        if masks[trait] & bit:
            raise ValueError(f"duplicate profile/trait test: {key}, {trait}")
        masks[trait] |= bit


def input_population(profiles: Path, matrix: Path) -> InputPopulation:
    """Recover producer eligibility: paired profiles with N>0, non-MediaDive traits with K>0."""
    draws: dict[ProfileKey, int] = {}
    directions: dict[tuple[str, str], set[str]] = {}
    traits: list[str] = []
    positions: dict[str, int] = {}
    for header, row in table_rows(profiles):
        if not positions:
            if not {*PROFILE_KEY_COLUMNS, "taxa_count"}.issubset(header):
                raise ValueError(f"{profiles}: missing profile identifiers or taxa_count")
            traits = [
                col
                for col in header
                if col not in {*PROFILE_KEY_COLUMNS, "taxa_count"}
                and MEDIADIVE_MARKER not in col.lower()
            ]
            positions = {col: index for index, col in enumerate(header)}
        key = tuple(row[positions[col]].strip() for col in PROFILE_KEY_COLUMNS)
        profile_key = (key[0], key[1], key[2])
        if not all(profile_key) or profile_key[2] not in {"inc", "dec"}:
            raise ValueError(f"{profiles}: invalid profile identifier")
        if profile_key in draws:
            raise ValueError(f"{profiles}: duplicate input profile: {profile_key}")
        draws[profile_key] = int(number(row[positions["taxa_count"]], "taxa_count", integer=True))
        directions.setdefault(profile_key[:2], set()).add(profile_key[2])
    if not draws or not traits:
        raise ValueError(f"{profiles}: no input profiles or retained trait columns")
    paired = {key for key, observed in directions.items() if observed == {"inc", "dec"}}
    eligible = [key for key, count in draws.items() if key[:2] in paired and count > 0]
    background = dict.fromkeys(traits, 0)
    taxa: set[str] = set()
    positions = {}
    trait_positions: dict[str, int] = {}
    for header, row in table_rows(matrix):
        if not positions:
            if "taxon_id" not in header:
                raise ValueError(f"{matrix}: missing taxon_id")
            positions = {col: index for index, col in enumerate(header)}
            trait_positions = {trait: positions[trait] for trait in traits if trait in positions}
        taxon = row[positions["taxon_id"]].strip()
        if not taxon or taxon in taxa:
            raise ValueError(f"{matrix}: missing or duplicate taxon_id: {taxon!r}")
        taxa.add(taxon)
        for trait, position in trait_positions.items():
            # The producer interprets empty trait cells as zero and accepts integer floats.
            count = number(row[position].strip() or "0", trait, integer=True)
            background[trait] += count > 0
    if not taxa:
        raise ValueError(f"{matrix}: no taxon rows")
    if any(draws[key] > len(taxa) for key in draws if key[:2] in paired):
        raise ValueError("profile taxa_count exceeds input matrix population")
    return InputPopulation(
        profile_indices={key: index for index, key in enumerate(eligible)},
        paired_profiles=sum(key[:2] in paired for key in draws),
        background_m=len(taxa),
        background_k={trait: count for trait, count in background.items() if count > 0},
    )


def number(
    value: str,
    field: str,
    *,
    integer: bool = False,
    probability: bool = False,
    positive_infinity: bool = False,
) -> float:
    """Apply the shared numeric domains without materializing an ORA table (#724)."""
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"invalid {field}: {value!r}") from exc
    if (
        (not math.isfinite(parsed) and not (positive_infinity and parsed == math.inf))
        or parsed < 0
        or (integer and not parsed.is_integer())
        or (probability and parsed > 1)
    ):
        raise ValueError(f"invalid {field}: {value!r}")
    return parsed


@dataclass(frozen=True)
class RowContract:
    """Validate each complete row against the producer's declared analysis settings."""

    indices: dict[str, int]
    upper_tail: bool
    fdr: bool
    family: str
    threshold: float
    fold_cutoff: float
    background_m: float
    testable_profiles: float

    def selected(self, row: list[str]) -> bool:
        indices = self.indices
        for column in ID_COLUMNS:
            if not row[indices[column]].strip():
                raise ValueError(f"empty {column}")
        if row[indices["inc_dec"]] not in {"inc", "dec"}:
            raise ValueError("invalid inc_dec")
        values = {}
        for column in PROBABILITY_COLUMNS + FLOAT_COLUMNS + COUNT_COLUMNS + BACKGROUND_COLUMNS:
            values[column] = number(
                row[indices[column]],
                column,
                integer=column in COUNT_COLUMNS + BACKGROUND_COLUMNS,
                probability=column in PROBABILITY_COLUMNS,
            )
        if (
            values["background_m"] != self.background_m
            or values["testable_profiles"] != self.testable_profiles
        ):
            raise ValueError("row background_m or testable_profiles differs from metadata")
        if not 0 < values["background_k"] <= values["background_m"]:
            raise ValueError("background_k must be positive and not exceed background_m")
        if values["count"] > values["background_k"]:
            raise ValueError("count exceeds background_k")
        if values["expected"] <= 0 or values["testable_profiles"] <= 0:
            raise ValueError("testable rows require positive expected and testable_profiles")
        if row[indices["eligibility_status"]] != "testable":
            raise ValueError("invalid eligibility_status")
        if row[indices["correction_family"]] != self.family:
            raise ValueError("correction_family differs from FDR-Scope")
        if not self.fdr and values["fdr_q"] != values["p_value"]:
            raise ValueError("fdr_q differs from p_value with FDR disabled")
        selected = row[indices["is_significant"]]
        if selected not in {"True", "False"}:
            raise ValueError("invalid is_significant")
        statistic = values["fdr_q" if self.fdr else "p_value"]
        fold_passes = (
            values["fold"] >= self.fold_cutoff
            if self.upper_tail
            else values["fold"] <= self.fold_cutoff
        )
        if (selected == "True") != (statistic <= self.threshold and fold_passes):
            raise ValueError("is_significant disagrees with declared thresholds")
        return selected == "True"


def row_contract(header: list[str], metadata: dict[str, str]) -> RowContract:
    """Read settings once, accepting both tails, correction scopes and FDR modes."""
    for key, allowed in (
        ("Tail", {"upper", "lower"}),
        ("FDR", {"True", "False"}),
        ("FDR-Scope", {"global", "profile"}),
        ("Testable-Set", {"K>0 & N>0"}),
    ):
        if metadata.get(key) not in allowed:
            raise ValueError(f"missing or invalid {key}")
    fold_key = "Min-Fold" if metadata["Tail"] == "upper" else "Max-Fold"
    try:
        total_profiles = number(
            metadata["Total-Profiles-Input"], "Total-Profiles-Input", integer=True
        )
        testable_profiles = number(metadata["Testable-Profiles"], "Testable-Profiles", integer=True)
        if testable_profiles > total_profiles:
            raise ValueError("Testable-Profiles exceeds Total-Profiles-Input")
        return RowContract(
            indices={column: index for index, column in enumerate(header)},
            upper_tail=metadata["Tail"] == "upper",
            fdr=metadata["FDR"] == "True",
            family=metadata["FDR-Scope"],
            threshold=number(metadata["P-Threshold"], "P-Threshold", probability=True),
            fold_cutoff=number(metadata[fold_key], fold_key, positive_infinity=True),
            background_m=number(
                metadata["Background-Taxa-Universe-M"], "Background-Taxa-Universe-M", integer=True
            ),
            testable_profiles=testable_profiles,
        )
    except KeyError as exc:
        raise ValueError(f"missing metadata: {exc.args[0]}") from exc


def same_path(left: Path, right: Path) -> bool:
    """Include existing hard links when protecting input files from overwrite."""
    return left.resolve() == right.resolve() or (
        left.exists() and right.exists() and left.samefile(right)
    )


@dataclass(frozen=True)
class TableAudit:
    path: str
    sha256: str
    body_sha256: str
    body_bytes: int
    rows: int
    significant_rows: int
    metadata: dict[str, str]
    errors: list[str]


def read_table(path: Path, population: InputPopulation) -> TableAudit:
    """Stream a complete export and reconcile its declared and observed counts."""
    metadata: dict[str, str] = {}
    body_digest = hashlib.sha256()
    body_bytes = 0
    rows = 0
    significant = 0
    errors: list[str] = []
    masks = dict.fromkeys(population.background_k, 0)
    with path.open("rb") as handle:
        first_line = handle.readline()
        while first_line.startswith(b"#"):
            key, separator, value = first_line.decode("utf-8").strip()[1:].strip().partition(": ")
            if not separator or not key or key in metadata:
                raise ValueError(f"{path}: malformed or duplicate metadata: {first_line!r}")
            metadata[key] = value
            first_line = handle.readline()

        def lines() -> Iterator[str]:
            nonlocal body_bytes
            for line in chain((first_line,), handle):
                body_digest.update(line)
                body_bytes += len(line)
                yield line.decode("utf-8")

        reader = csv.reader(lines(), delimiter="\t", strict=True)
        header = next(reader, [])
        if len(set(header)) != len(header) or not COMPLETE_COLUMNS.issubset(header):
            raise ValueError(f"{path}: missing or duplicate complete-result columns")
        try:
            contract = row_contract(header, metadata)
            population.check_metadata(metadata)
        except ValueError as exc:
            raise ValueError(f"{path}: {exc}") from exc
        for row in reader:
            rows += 1
            if len(row) != len(header):
                raise ValueError(f"{path}: malformed table row {rows}")
            try:
                significant += contract.selected(row)
                population.observe(row, contract.indices, masks)
            except ValueError as exc:
                raise ValueError(f"{path}: row {rows}: {exc}") from exc
    for key, observed in ((META_TESTS_PERFORMED, rows), (META_SIGNIFICANT_HITS, significant)):
        try:
            declared = int(metadata[key])
        except (KeyError, ValueError) as exc:
            raise ValueError(f"{path}: missing or invalid {key}") from exc
        if declared != observed:
            errors.append(f"{key}: declared {declared}, observed {observed}")
    status = classify_result_status(rows, significant)
    if metadata.get(META_RESULT_STATUS) != status:
        errors.append(
            f"{META_RESULT_STATUS}: declared {metadata.get(META_RESULT_STATUS)!r}, "
            f"observed {status!r}"
        )
    complete = (1 << len(population.profile_indices)) - 1
    if any(mask != complete for mask in masks.values()):
        errors.append("incomplete eligible profile/trait family")
    return TableAudit(
        path=str(path),
        sha256=file_sha256(path),
        body_sha256=body_digest.hexdigest(),
        body_bytes=body_bytes,
        rows=rows,
        significant_rows=significant,
        metadata=metadata,
        errors=errors,
    )


def audit_pair(pair: dict[str, str], base: Path) -> dict[str, object]:
    paths = {key: (base / value).resolve() for key, value in pair.items() if key != "label"}
    if same_path(paths["historical"], paths["current"]):
        raise ValueError(f"{pair['label']}: historical and current paths must differ")
    population = input_population(paths["profiles"], paths["taxa_matrix"])
    historical = read_table(paths["historical"], population)
    current = read_table(paths["current"], population)
    old_meta = historical.metadata
    new_meta = current.metadata
    inputs = {
        key: {"path": str(paths[key]), "sha256": file_sha256(paths[key])}
        for key in FINGERPRINTS.values()
    }
    changes = {
        key: {"historical": old_meta.get(key), "current": new_meta.get(key)}
        for key in sorted(old_meta.keys() | new_meta.keys())
        if old_meta.get(key) != new_meta.get(key)
    }
    errors = [
        f"{name}: {error}"
        for name, table in (("historical", historical), ("current", current))
        for error in table.errors
    ]
    body_equal = (
        historical.body_sha256 == current.body_sha256
        and historical.body_bytes == current.body_bytes
    )
    if not body_equal:
        errors.append("Table body bytes differ")
    for key in changes:
        if key not in FINGERPRINTS or key in old_meta:
            errors.append(f"Unexpected metadata change: {key}")
    for key, input_key in FINGERPRINTS.items():
        if new_meta.get(key) != inputs[input_key]["sha256"]:
            errors.append(f"{key}: current export does not match supplied input")
    return {
        "label": pair["label"],
        "historical": asdict(historical),
        "current": asdict(current),
        "inputs": inputs,
        "metadata_changes": changes,
        "table_body_equal": body_equal,
        "passed": not errors,
        "errors": errors,
    }


def audit_manifest(manifest: Path, output: Path) -> dict[str, object]:
    raw = json.loads(manifest.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or set(raw) != {"pairs"}:
        raise ValueError("Manifest must contain exactly the 'pairs' key")
    pairs = raw["pairs"]
    if not isinstance(pairs, list) or not pairs:
        raise ValueError("Manifest 'pairs' must be a nonempty list")
    labels: set[str] = set()
    for pair in pairs:
        if (
            not isinstance(pair, dict)
            or set(pair) != PAIR_FIELDS
            or any(not isinstance(value, str) or not value.strip() for value in pair.values())
        ):
            raise ValueError(f"Each pair must supply nonempty strings for {sorted(PAIR_FIELDS)}")
        if pair["label"] in labels:
            raise ValueError(f"Duplicate pair label: {pair['label']}")
        labels.add(pair["label"])
        for key in PAIR_FIELDS - {"label"}:
            if same_path(manifest.parent / pair[key], output):
                raise ValueError("Audit output must not overwrite any input")
    if same_path(manifest, output):
        raise ValueError("Audit output must not overwrite manifest")
    results = []
    for pair in pairs:
        try:
            results.append(audit_pair(pair, manifest.parent))
        except (OSError, ValueError, csv.Error) as exc:
            results.append({"label": pair["label"], "passed": False, "errors": [str(exc)]})
    return {
        "issue": 720,
        "manifest": str(manifest.resolve()),
        "manifest_sha256": file_sha256(manifest),
        "passed": all(result["passed"] for result in results),
        "pairs": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit_manifest(args.manifest, args.output)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"ORA refresh audit: {'PASS' if result['passed'] else 'FAIL'} ({args.output})")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
