"""Audit signed-rank inference on pinned persisted study contrasts (#1073).

This is a parent-output sensitivity audit, not an upstream execution verifier.
Manifest metadata is preserved as declarations. No count transform, delta-z
calculation, native arm execution, or historical Git verification occurs here.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.stats import false_discovery_control, wilcoxon

from mtp.three_test_table.signed_rank import (
    LEGACY_WILCOXON_CONTRACT,
    WILCOXON_CONTRACT,
    study_wilcoxon,
)

_MANIFEST = "three_test_run.json"
_OUTPUTS = (
    "composition_pair_delta_z.tsv",
    "composition_study_delta_z.tsv",
    "per_trait_evidence.tsv",
)
_BASELINE_TOLERANCE = {"relative": 1e-12, "absolute": 1e-15}
_CHANGE_TOLERANCE = {"relative": 1e-12, "absolute": 1e-15}
_METRICS = ("p", "q_bh", "q_by")


@dataclass(frozen=True)
class AuditResult:
    """Read-only audit tables and hashes for the publishing caller to recheck."""

    table: pd.DataFrame
    summary: dict[str, Any]
    input_hashes: dict[Path, str]


def _object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for name, value in pairs:
        if name in result:
            raise ValueError(f"Duplicate manifest key: {name}")
        result[name] = value
    return result


def _digest(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _consume(path: Path, hashes: dict[Path, str], *, regular: bool = False) -> bytes:
    if not path.is_file() or (regular and path.is_symlink()):
        raise ValueError(f"Parent output must be a regular file: {path.name}")
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    for spelling in (path, path.resolve(strict=True)):
        if spelling in hashes and hashes[spelling] != digest:
            raise ValueError("Parent input identity changed while reading")
        hashes[spelling] = digest
    return raw


def _records(raw: bytes, keys: list[str]) -> tuple[list[str], list[dict[str, str]]]:
    reader = csv.reader(io.StringIO(raw.decode("utf-8"), newline=""), delimiter="\t")
    header = next(reader, [])
    if (
        not header
        or header[: len(keys)] != keys
        or len(set(header)) != len(header)
        or any(not name.strip() for name in header)
    ):
        raise ValueError("Parent TSV has missing, duplicate, or reordered identity columns")
    records: list[dict[str, str]] = []
    identities: set[tuple[str, ...]] = set()
    for fields in reader:
        if len(fields) != len(header):
            raise ValueError("Parent TSV row width differs from its header")
        row = dict(zip(header, fields, strict=True))
        identity = tuple(row[key] for key in keys)
        if any(not value.strip() for value in identity) or identity in identities:
            raise ValueError("Parent TSV has blank or duplicate identities")
        identities.add(identity)
        records.append(row)
    if not records:
        raise ValueError("Parent TSV has no rows")
    return header, records


def _number(value: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError("Parent numerical field is invalid") from error
    if not math.isfinite(result):
        raise ValueError("Parent numerical field must be finite")
    return result


def _integer(value: str, expected: int) -> None:
    if not value.isdecimal() or int(value) != expected:
        raise ValueError("Parent support or family count differs from its ledger")


def _metadata(record: dict[str, Any], scope: str) -> list[str]:
    if (
        scope not in {"off", "phenotype"}
        or record.get("version") != "three-test-per-trait-v1"
        or record.get("strain_rollup") != scope
        or record.get("scope_role") != ("primary" if scope == "phenotype" else "sensitivity")
        or record.get("code_dirty") is not False
        or not isinstance(record.get("code_commit"), str)
        or re.fullmatch(r"[0-9a-f]{40}", record["code_commit"]) is None
    ):
        raise ValueError("Parent declared version, scope, or producer identity is invalid")
    for field in (
        "run_id",
        "core_run_id",
        "comparison_run_id",
        "contrast_population",
        "direction_basis",
    ):
        if not isinstance(record.get(field), str) or not record[field].strip():
            raise ValueError(f"Parent declaration is missing: {field}")
    traits = record.get("retained_traits")
    if (
        not isinstance(traits, list)
        or not traits
        or any(not isinstance(trait, str) or not trait.strip() for trait in traits)
        or len(set(traits)) != len(traits)
        or type(record.get("family_size")) is not int
        or record["family_size"] != len(traits)
    ):
        raise ValueError("Parent retained trait catalog is invalid")
    return list(traits)


def _load(
    parent: Path, expected_sha256: str, scope: str
) -> tuple[dict[str, Any], list[dict[str, str]], NDArray[np.float64], int, dict[Path, str]]:
    supplied = parent / _MANIFEST if parent.is_dir() else parent
    supplied = supplied if supplied.is_absolute() else Path.cwd() / supplied
    hashes: dict[Path, str] = {}
    raw = _consume(supplied, hashes)
    if not _digest(expected_sha256) or hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise ValueError("Parent manifest SHA256 pin differs")
    record = json.loads(raw, object_pairs_hook=_object)
    if not isinstance(record, dict):
        raise ValueError("Parent manifest must be an object")
    traits = _metadata(record, scope)
    manifest = supplied.resolve(strict=True)
    root = manifest.parent
    inventory = record.get("outputs")
    if (
        not isinstance(inventory, dict)
        or set(inventory) != set(_OUTPUTS)
        or {path.name for path in root.iterdir()} != {_MANIFEST, *_OUTPUTS}
    ):
        raise ValueError("Parent output inventory must contain exactly the three completed tables")
    tables: dict[str, bytes] = {}
    for name in _OUTPUTS:
        item = inventory[name]
        path = supplied.parent / name if supplied.parent.resolve() == root else root / name
        data = _consume(path, hashes, regular=True)
        if (
            not isinstance(item, dict)
            or set(item) != {"sha256", "size_bytes"}
            or not _digest(item["sha256"])
            or type(item["size_bytes"]) is not int
            or item["sha256"] != hashlib.sha256(data).hexdigest()
            or item["size_bytes"] != len(data)
        ):
            raise ValueError(f"Parent output hash or size differs: {name}")
        tables[name] = data
    pair_header, pairs = _records(tables[_OUTPUTS[0]], ["study_id", "experiment_id"])
    study_header, studies = _records(tables[_OUTPUTS[1]], ["study_id"])
    _, evidence = _records(tables[_OUTPUTS[2]], ["trait"])
    if (
        pair_header != ["study_id", "experiment_id", *traits]
        or study_header != ["study_id", *traits]
        or [row["trait"] for row in evidence] != traits
        or {row["study_id"] for row in pairs} != {row["study_id"] for row in studies}
    ):
        raise ValueError("Parent ledger population or ordered trait catalog differs")
    # Python float conversion preserves the serialized float64 values and their ties.
    values = np.asarray([[_number(row[trait]) for trait in traits] for row in studies], dtype=float)
    for row in pairs:
        for trait in traits:
            _number(row[trait])
    return record, evidence, values, len(pairs), hashes


def _legacy(values: NDArray[np.float64]) -> dict[str, Any]:
    nonzero = int(np.count_nonzero(values))
    if nonzero:
        result = wilcoxon(
            values.copy(),
            alternative="two-sided",
            zero_method="wilcox",
            method="approx",
            correction=False,
        )
        statistic, probability = float(result.statistic), float(result.pvalue)
    else:
        statistic, probability = math.nan, 1.0
    return {
        "statistic": statistic,
        "p": probability,
        "nonzero_studies": nonzero,
        "zero_studies": len(values) - nonzero,
        "status": "tested" if nonzero else "untestable",
        "reason": "" if nonzero else "all_tied",
        "method": "approx" if nonzero else "not_tested",
    }


def _correct(results: list[dict[str, Any]]) -> None:
    for correction in ("bh", "by"):
        corrected = false_discovery_control([row["p"] for row in results], method=correction)
        for row, probability in zip(results, corrected, strict=True):
            row[f"q_{correction}"] = float(probability) if row["status"] == "tested" else math.nan
    for row in results:
        if row["status"] == "untestable":
            row["p"] = math.nan


def _archived(
    row: dict[str, str],
    reproduced: dict[str, Any],
    record: dict[str, Any],
    n_studies: int,
    n_pairs: int,
) -> dict[str, Any]:
    expected = {
        "composition_contract": "full-catalog-study-delta-z-v1",
        "composition_test_unit": "study",
        "composition_effect_units": "standardized_log1p_1000_retained_fraction_difference",
        "composition_zero_total_policy": "retained_as_zero_fraction",
        "study_wilcoxon_contract": LEGACY_WILCOXON_CONTRACT,
        **{
            name: record[name]
            for name in ("strain_rollup", "contrast_population", "direction_basis", "core_run_id")
        },
    }
    for field, value in expected.items():
        if row.get(field) != value:
            raise ValueError(f"Parent evidence contract differs: {field}")
    for field, count in (
        ("composition_family_size", record["family_size"]),
        ("composition_studies", n_studies),
        ("composition_pairs", n_pairs),
    ):
        _integer(row.get(field, ""), count)
    result = dict(reproduced)
    for field in ("status", "reason"):
        if row.get(f"study_wilcoxon_{field}") != reproduced[field]:
            raise ValueError(f"Legacy Wilcoxon {field} does not reproduce")
    for field in ("nonzero_studies", "zero_studies"):
        _integer(row.get(f"study_wilcoxon_{field}", ""), reproduced[field])
    for field in ("statistic", *_METRICS):
        value = row.get(f"study_wilcoxon_{field}")
        if value is None:
            raise ValueError(f"Missing legacy Wilcoxon field: {field}")
        if reproduced["status"] == "untestable":
            if value != "":
                raise ValueError("Untestable legacy inference must be blank")
            continue
        number = _number(value)
        if field != "statistic" and not 0 <= number <= 1:
            raise ValueError("Legacy probability is outside [0, 1]")
        same = (
            number == reproduced[field]
            if field == "statistic"
            else math.isclose(
                number,
                reproduced[field],
                rel_tol=_BASELINE_TOLERANCE["relative"],
                abs_tol=_BASELINE_TOLERANCE["absolute"],
            )
        )
        if not same:
            raise ValueError(f"Legacy Wilcoxon {field} does not reproduce")
        result[field] = number
    return result


def audit_parent(parent: Path, expected_sha256: str, *, scope: str) -> AuditResult:
    """Reproduce v1 before comparing v2 on exactly the pinned stored study vectors."""
    record, evidence, values, n_pairs, hashes = _load(parent, expected_sha256, scope)
    reproduced = [_legacy(values[:, index]) for index in range(values.shape[1])]
    _correct(reproduced)
    archived = [
        _archived(row, old, record, len(values), n_pairs)
        for row, old in zip(evidence, reproduced, strict=True)
    ]
    updated: list[dict[str, Any]] = []
    for index in range(values.shape[1]):
        result = study_wilcoxon(values[:, index].copy())
        updated.append(
            {
                "statistic": result.statistic,
                "p": result.pvalue,
                "nonzero_studies": result.nonzero_studies,
                "zero_studies": result.zero_studies,
                "status": "untestable" if result.reason else "tested",
                "reason": result.reason,
                "method": result.method,
            }
        )
    _correct(updated)
    rows: list[dict[str, Any]] = []
    for source, saved, old, new in zip(evidence, archived, reproduced, updated, strict=True):
        row: dict[str, Any] = {
            "trait": source["trait"],
            "trait_label": source.get("trait_label", ""),
            "old_contract": LEGACY_WILCOXON_CONTRACT,
            "new_contract": WILCOXON_CONTRACT,
            "status_changed": old["status"] != new["status"],
        }
        for prefix, data in (("archived", saved), ("reproduced", old), ("new", new)):
            row.update({f"{prefix}_{name}": value for name, value in data.items()})
        for metric in _METRICS:
            comparable = old["status"] == new["status"] == "tested"
            row[f"{metric}_changed_exact"] = comparable and old[metric] != new[metric]
            row[f"{metric}_changed_material"] = comparable and not math.isclose(
                old[metric],
                new[metric],
                rel_tol=_CHANGE_TOLERANCE["relative"],
                abs_tol=_CHANGE_TOLERANCE["absolute"],
            )
            row[f"{metric}_absolute_delta"] = (
                abs(new[metric] - old[metric]) if comparable else math.nan
            )
            if metric != "p":
                row[f"{metric}_call_changed_005"] = (
                    old["status"] == "tested" and old[metric] <= 0.05
                ) != (new["status"] == "tested" and new[metric] <= 0.05)
        rows.append(row)
    table = pd.DataFrame.from_records(rows)
    changes = {}
    for metric in _METRICS:
        deltas = table[f"{metric}_absolute_delta"].dropna()
        changes[metric] = {
            "exact": int(table[f"{metric}_changed_exact"].sum()),
            "material": int(table[f"{metric}_changed_material"].sum()),
            "comparable_tested": len(deltas),
            "max_absolute_delta": float(deltas.max()) if len(deltas) else None,
        }
        if metric != "p":
            changes[metric]["call_changes_005"] = int(table[f"{metric}_call_changed_005"].sum())
    summary = {
        "verification_scope": "pinned_parent_outputs_and_legacy_reproduction_only",
        "parent_declarations": {
            name: record[name]
            for name in (
                "run_id",
                "code_commit",
                "code_dirty",
                "core_run_id",
                "comparison_run_id",
                "strain_rollup",
                "contrast_population",
                "direction_basis",
            )
        },
        "family_size": len(evidence),
        "n_studies": len(values),
        "n_pairs": n_pairs,
        "old_contract": LEGACY_WILCOXON_CONTRACT,
        "new_contract": WILCOXON_CONTRACT,
        "baseline_tolerance": dict(_BASELINE_TOLERANCE),
        "change_tolerance": dict(_CHANGE_TOLERANCE),
        "changes": changes,
        "status_changes": int(table.status_changed.sum()),
        "method_counts": {
            str(method): int(count) for method, count in table.new_method.value_counts().items()
        },
    }
    for path, digest in hashes.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError("Parent input changed during audit")
    return AuditResult(table, summary, hashes)


__all__ = ["AuditResult", "audit_parent"]
