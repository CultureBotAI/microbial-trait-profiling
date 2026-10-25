"""Import unchanged native statistics and retain each producer's declared family (#1044)."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from mtp.datamodel.multiple_testing_pydantic import AdjustedHypothesis, HypothesisFamily
from mtp.multiple_testing import correct_family
from mtp.multiplicity.families import (
    FamilyInput,
    NativeObservation,
    noel_family,
    probability,
    support_count,
)
from mtp.multiplicity.native import native_family
from mtp.multiplicity.sources import VerifiedSources

__all__ = ["native_tables"]

_NOEL_COUNTS = ("studies_up", "studies_down", "studies_tied", "studies_scored")
_NOEL_PROBABILITIES = {
    "share_up": "share_up",
    "ci_low": "ci_low",
    "ci_high": "ci_high",
    "p_value": "binomial_p",
    "q_by_nominal": "by_q_nominal",
    "q_bh_nominal": "bh_q_nominal",
}
_NATIVE_VALUES = {
    "study_mean_dz": "study_mean_increased_minus_decreased",
    "study_ci95_lower": "study_ci95_lower",
    "study_ci95_upper": "study_ci95_upper",
    "study_t_statistic": "study_t_statistic",
    "experiment_mean_dz": "mean_increased_minus_decreased",
}
_NATIVE_COLUMNS = (
    "study_count",
    *_NATIVE_VALUES,
    "study_p_value",
    "study_q_bh",
    "study_q_by",
    "experiment_count",
)


def _observations(
    family: FamilyInput,
) -> tuple[dict[str, NativeObservation], dict[str, AdjustedHypothesis]]:
    """Check exact membership before indexing rows, including duplicate scientific identities."""
    declared = HypothesisFamily(**family.metadata.model_dump(), members=list(family.declarations()))
    observations = list(family.observations())
    correction = correct_family(declared, [item.raw for item in observations])
    return (
        {item.raw.hypothesis_id: item for item in observations},
        {item.hypothesis_id: item for item in correction.results},
    )


def _number(row: dict[str, str], key: str) -> float | None:
    value = row[key]
    if not value:
        return None
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Native numerical field must be finite or blank: {key}")
    return number


def native_tables(sources: VerifiedSources) -> tuple[pd.DataFrame, tuple[Path, ...]]:
    """Join full Noel votes and selected native study statistics by qualified trait.

    Sources must already be execution-verified and bound to recorded hashes.
    The existing adapters check producer contracts and consumed bytes. The shared
    correction engine additionally requires the exact declared membership before
    rows are indexed, so a duplicate cannot overwrite another row in the join.

    Noel p/BY/BH are imported at their original export precision. Native study
    BH is also unchanged; the separately named BY uses only the independently
    declared selected family, including its untestable members at internal p=1.
    Unselected catalog traits have explicit status and no invented statistics.
    """
    noel_input = noel_family(sources)
    noel, _ = _observations(noel_input)
    native_input = native_family(sources)
    native, adjusted = _observations(native_input)
    records: list[dict[str, str | int | float | None]] = []
    for trait in sources.features["count_profiles"]:
        screen = noel[trait]
        record: dict[str, str | int | float | None] = {
            "trait": trait,
            "noel_family_id": noel_input.metadata.family_id,
            "noel_family_size": noel_input.metadata.expected_family_size,
            "noel_status": screen.row["test_status"],
            **{f"noel_{name}": support_count(screen.row, name) for name in _NOEL_COUNTS},
            **{
                f"noel_{target}": probability(screen.row[source])
                for target, source in _NOEL_PROBABILITIES.items()
            },
            "native_selected_family_id": native_input.metadata.family_id,
            "native_selected_family_size": native_input.metadata.expected_family_size,
            "native_selected_status": "not_selected",
            **dict.fromkeys(f"native_selected_{name}" for name in _NATIVE_COLUMNS),
        }
        if screen.raw.status == "tested" and any(
            record[f"noel_{name}"] is None for name in _NOEL_PROBABILITIES
        ):
            raise ValueError("Tested Noel trait lacks an exported probability or interval")
        if trait in native:
            observed = native[trait]
            record.update(
                {
                    "native_selected_status": observed.row["study_test_status"],
                    "native_selected_study_count": support_count(observed.row, "study_count"),
                    "native_selected_experiment_count": support_count(observed.row, "pair_count"),
                    "native_selected_study_p_value": observed.raw.p_value,
                    "native_selected_study_q_bh": observed.native_q,
                    "native_selected_study_q_by": adjusted[trait].q_by,
                    **{
                        f"native_selected_{target}": _number(observed.row, source)
                        for target, source in _NATIVE_VALUES.items()
                    },
                }
            )
            if observed.raw.status == "tested" and any(
                record[f"native_selected_{name}"] is None for name in _NATIVE_VALUES
            ):
                raise ValueError("Tested native trait lacks an exported statistic or interval")
        records.append(record)
    consumed = tuple(sorted(set(noel_input.source_paths) | set(native_input.source_paths)))
    return pd.DataFrame.from_records(records), consumed
