"""Import unchanged native statistics and retain each producer's declared family (#1044)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
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

__all__ = ["native_tables", "validate_native_overlap"]

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
    "study_sd": "study_std_increased_minus_decreased",
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
    "profile_support",
)
_OVERLAP_COUNTS = {
    "composition_studies": "study_count",
    "composition_pairs": "experiment_count",
    "composition_profile_support": "profile_support",
}
_OVERLAP_EFFECTS = {
    "study_mean_delta_z": "study_mean_dz",
    "experiment_mean_delta_z": "experiment_mean_dz",
}
_OVERLAP_INFERENCE = {
    "study_t_statistic": "study_t_statistic",
    "study_t_p": "study_p_value",
    "study_mean_delta_z_ci95_low": "study_ci95_lower",
    "study_mean_delta_z_ci95_high": "study_ci95_upper",
}


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
                    "native_selected_profile_support": support_count(
                        observed.row, "profile_support"
                    ),
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


def _overlap_number(row: pd.Series, column: str, trait: str) -> float:
    value = row[column]
    if not isinstance(value, (float, int, np.floating, np.integer)) or not math.isfinite(value):
        raise ValueError(f"Native overlap requires a finite number: {trait}, {column}")
    return float(value)


def _overlap_equal(
    full: pd.Series, observed: pd.Series, full_column: str, native_column: str, trait: str
) -> None:
    first = _overlap_number(full, full_column, trait)
    second = _overlap_number(observed, f"native_selected_{native_column}", trait)
    # Native TSVs use default, unrounded float export. Accommodate transform and
    # reduction roundoff, never display precision; tiny p values need relative
    # agreement rather than being swallowed by an absolute tolerance.
    absolute = 0.0 if full_column == "study_t_p" else 1e-14
    if not math.isclose(first, second, rel_tol=1e-10, abs_tol=absolute):
        raise ValueError(f"Native overlap differs: {trait}, {full_column}, {native_column}")


def _overlap_status(full: pd.Series, observed: pd.Series, trait: str) -> tuple[str, str]:
    """Validate overlap without changing either producer's degeneracy threshold."""
    for full_column, native_column in _OVERLAP_COUNTS.items():
        first = _overlap_number(full, full_column, trait)
        second = _overlap_number(observed, f"native_selected_{native_column}", trait)
        if first < 0 or first != int(first) or first != second:
            raise ValueError(f"Native overlap count differs: {trait}, {full_column}")
    for full_column, native_column in _OVERLAP_EFFECTS.items():
        _overlap_equal(full, observed, full_column, native_column, trait)
    count = _overlap_number(full, "composition_studies", trait)
    std = _overlap_number(observed, "native_selected_study_sd", trait)
    mean = _overlap_number(observed, "native_selected_study_mean_dz", trait)
    if std < 0:
        raise ValueError(f"Native overlap has negative dispersion: {trait}")
    if count < 2:
        expected = "insufficient_support"
    else:
        _overlap_equal(full, observed, "study_delta_z_sd", "study_sd", trait)
        expected = (
            ("all_tied" if std == 0 and mean == 0 else "insufficient_variance")
            if np.isclose(std, 0.0)
            else "tested"
        )
    if observed["native_selected_status"] != expected:
        raise ValueError(f"Native overlap has unexplained native testability: {trait}")
    if expected == "tested":
        if full["study_t_status"] != "tested":
            raise ValueError(f"Native overlap has unexplained full-family testability: {trait}")
        for full_column, native_column in _OVERLAP_INFERENCE.items():
            _overlap_equal(full, observed, full_column, native_column, trait)
        return "tested_match", ""
    if any(
        not pd.isna(observed[f"native_selected_{column}"]) for column in _OVERLAP_INFERENCE.values()
    ):
        raise ValueError(f"Native overlap has defined inference for an untestable trait: {trait}")
    if full["study_t_status"] == "untestable" and full["study_t_reason"] == expected:
        return "untestable_match", expected
    if (
        full["study_t_status"] == "tested"
        and expected == "insufficient_variance"
        and _overlap_number(full, "study_delta_z_sd", trait) > 1e-12
    ):
        return (
            "native_sd_threshold_difference",
            "native np.isclose(SD, 0) excludes SD <= 1e-8; full-family SD > 1e-12 remains tested",
        )
    raise ValueError(f"Native overlap has unexplained testability difference: {trait}")


def validate_native_overlap(composition: pd.DataFrame, native: pd.DataFrame) -> pd.DataFrame:
    """Return native rows annotated with an enforced full-family overlap invariant.

    Both inputs must contain the same unique full catalog. For native-selected
    traits, require identical study/pair/support counts, both means and study SD.
    For jointly tested traits, also require matching raw t, p and study 95% CI.
    Q values are deliberately excluded: native-selected and full-catalog families
    differ. Native TSVs have no decimal rounding; numerical tolerance is relative
    1e-10 plus absolute 1e-14, except p values use zero absolute tolerance.

    Native ``study_std_increased_minus_decreased`` is imported as
    ``native_selected_study_sd``. Its producer uses ``np.isclose(SD, 0)``
    (atol=1e-8), whereas the full-family adapter uses SD <= 1e-12. Preserve and
    label that legitimate native exclusion, comparing descriptive quantities but
    leaving its undefined t/p/CI blank. Other status or population mismatches fail.
    """
    required_full = {
        "trait",
        "study_t_status",
        "study_t_reason",
        "study_delta_z_sd",
        *_OVERLAP_COUNTS,
        *_OVERLAP_EFFECTS,
        *_OVERLAP_INFERENCE,
    }
    required_native = {
        "trait",
        "native_selected_status",
        "native_selected_study_sd",
        *(
            f"native_selected_{name}"
            for name in (
                *_OVERLAP_COUNTS.values(),
                *_OVERLAP_EFFECTS.values(),
                *_OVERLAP_INFERENCE.values(),
            )
        ),
    }
    for frame, required in ((composition, required_full), (native, required_native)):
        if (
            not frame.columns.is_unique
            or not required.issubset(frame.columns)
            or not frame.trait.is_unique
            or any(not isinstance(trait, str) or not trait.strip() for trait in frame.trait)
        ):
            raise ValueError("Native overlap requires unique traits and complete columns")
    if set(composition.trait) != set(native.trait):
        raise ValueError("Native overlap catalog membership differs")
    full_rows = {str(row["trait"]): row for _, row in composition.iterrows()}
    result = native.copy()
    statuses, reasons = [], []
    for _, observed in native.iterrows():
        trait = str(observed["trait"])
        status, reason = (
            ("not_selected", "outside native support/variance selection")
            if observed["native_selected_status"] == "not_selected"
            else _overlap_status(full_rows[trait], observed, trait)
        )
        statuses.append(status)
        reasons.append(reason)
    result["native_selected_overlap_status"] = statuses
    result["native_selected_overlap_reason"] = reasons
    return result
