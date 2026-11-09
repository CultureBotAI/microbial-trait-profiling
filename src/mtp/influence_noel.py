"""Hash-bound Noel six-trait baseline comparisons for the influence sidecar (#1029)."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, BinaryIO

from scipy.stats import false_discovery_control

from mtp.taxon_contributions import Record

RELATIVE_TOLERANCE = 1e-10
ABSOLUTE_TOLERANCE = 0.0
REQUIRED = {
    "scenario",
    "scope",
    "representation",
    "trait_key",
    "denominator_trait",
    "unit",
    "correction",
    "before_p",
    "before_q",
}


class _HashedReader(io.RawIOBase):
    """Hash compressed bytes as the gzip parser consumes them, with bounded buffers."""

    def __init__(self, source: BinaryIO) -> None:
        super().__init__()
        self.source = source
        self.digest = hashlib.sha256()

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: Any) -> int:
        block = self.source.read(len(buffer))
        buffer[: len(block)] = block
        self.digest.update(block)
        return len(block)


def _read_baseline(
    path: Path, expected_sha256: str, expected: set[tuple[str, str]]
) -> dict[tuple[str, str], dict[str, str]]:
    selected = {}
    with path.open("rb") as source, _HashedReader(source) as hashed:
        with (
            gzip.GzipFile(fileobj=hashed, mode="rb") as compressed,
            io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text,
        ):
            reader = csv.reader(text, delimiter="\t")
            header = next(reader, [])
            if len(header) != len(set(header)) or not REQUIRED.issubset(header):
                raise ValueError("Noel effects table has missing or duplicate columns")
            for number, values in enumerate(reader, 2):
                if len(values) != len(header):
                    raise ValueError(f"Noel effects row {number} has the wrong field count")
                row = dict(zip(header, values, strict=True))
                if row["scenario"] != "baseline" or row["representation"] != "shared_fraction":
                    continue
                key = row["scope"], row["trait_key"]
                if key not in expected or key in selected:
                    raise ValueError(f"Unexpected or duplicate Noel baseline identity: {key}")
                selected[key] = row
        if hashed.digest.hexdigest() != expected_sha256:
            raise ValueError("Noel effects consumed-byte hash differs from the verified arm digest")
    if set(selected) != expected:
        raise ValueError(f"Missing Noel baseline family rows: {sorted(expected - set(selected))}")
    return selected


def _probability(value: object, name: str) -> float | None:
    if value == "":
        return None
    if isinstance(value, bool) or not isinstance(value, str | int | float):
        raise ValueError(f"Invalid {name}: expected a probability or an empty untestable value")
    try:
        result = float(value)
    except ValueError as error:
        raise ValueError(f"Invalid {name}: {value!r}") from error
    if not math.isfinite(result) or not 0 <= result <= 1:
        raise ValueError(f"Invalid {name}: {value!r}")
    return result


def _match(left: float | None, right: float | None, identity: str) -> None:
    equal = (
        left is right
        if left is None or right is None
        else math.isclose(left, right, rel_tol=RELATIVE_TOLERANCE, abs_tol=ABSOLUTE_TOLERANCE)
    )
    if not equal:
        raise ValueError(f"Noel baseline mismatch for {identity}: sidecar={left!r}, Noel={right!r}")


def check_noel_baseline(path: Path, expected_sha256: str, baseline: list[Record]) -> Record:
    """Compare stored p/BH and independently derive BY from each complete Noel family.

    There is no native BY-six column. Undefined members retain their family slot
    at internal p=1 and remain undefined in comparisons. No absolute tolerance
    can turn an extremely small nonzero p/q into zero or a much larger value.
    """
    grouped: defaultdict[str, list[Record]] = defaultdict(list)
    expected: set[tuple[str, str]] = set()
    for row in baseline:
        key = row["scope"], row["trait_key"]
        if key in expected or row["family_size"] != 6:
            raise ValueError("Duplicate sidecar baseline or non-six-member family")
        expected.add(key)
        grouped[row["scope"]].append(row)
    if not grouped or any(len(family) != 6 for family in grouped.values()):
        raise ValueError("Every baseline scope must contain all six focal traits")
    trait_sets = {frozenset(row["trait_key"] for row in family) for family in grouped.values()}
    if len(trait_sets) != 1:
        raise ValueError("Baseline scopes must contain the same six focal traits")
    native = _read_baseline(path, expected_sha256, expected)
    comparisons = []
    for scope, family in sorted(grouped.items()):
        original = [native[(scope, row["trait_key"])] for row in family]
        for row in original:
            if (
                row["unit"] != "fraction_difference"
                or row["correction"] != "BH6"
                or row["denominator_trait"]
            ):
                raise ValueError("Noel baseline must be the shared_fraction BH6 family")
        pvalues = [_probability(row["before_p"], "Noel before_p") for row in original]
        bh = [_probability(row["before_q"], "Noel before_q") for row in original]
        by = false_discovery_control([1.0 if p is None else p for p in pvalues], method="by")
        for row, p, q, derived in zip(family, pvalues, bh, by, strict=True):
            if (p is None) != (q is None):
                raise ValueError("Noel baseline p/BH missingness differs")
            derived_q = None if p is None else float(derived)
            sidecar_values = {
                field: _probability(row[field], field) for field in ("p", "q_bh_six", "q_by_six")
            }
            for field, value in (("p", p), ("q_bh_six", q), ("q_by_six", derived_q)):
                _match(sidecar_values[field], value, f"{scope}/{row['trait_key']}/{field}")
            comparisons.append(
                {
                    "scope": scope,
                    "trait_key": row["trait_key"],
                    "family_size": 6,
                    "sidecar_p": sidecar_values["p"],
                    "noel_p": p,
                    "sidecar_q_bh_six": sidecar_values["q_bh_six"],
                    "noel_q_bh_six": q,
                    "sidecar_q_by_six": sidecar_values["q_by_six"],
                    "noel_q_by_six_derived": derived_q,
                    "p_status": "matched",
                    "bh_status": "matched",
                    "by_status": "recomputed_from_noel_p",
                }
            )
    return {
        "status": "matched",
        "rows": len(comparisons),
        "families": len(grouped),
        "representation": "shared_fraction",
        "scenario": "baseline",
        "family_size": 6,
        "source_sha256": expected_sha256,
        "by_source": "derived_from_noel_p_values",
        "relative_tolerance": RELATIVE_TOLERANCE,
        "absolute_tolerance": ABSOLUTE_TOLERANCE,
        "comparisons": comparisons,
    }
