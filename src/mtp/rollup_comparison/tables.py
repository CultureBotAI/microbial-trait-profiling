"""Compare completed scope-specific evidence without recomputing inference (#1041).

Inputs are full parent per-trait tables and persisted study-mean Delta-z ledgers.
Every imported column remains separately labelled by scope. Selection uses native
within-scope BY values or comparison-weighted means in declared signed rank views;
neither correction nor ranking is recalculated on the shared-trait intersection.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

ACTIVITY_EPSILON = 1e-12
MINIMUM_ACTIVE_STUDIES = 50
BY_THRESHOLD = 0.05
RANK_COLUMNS = tuple(
    f"{effect}_rank_{order}"
    for effect in ("experiment_mean_delta_z", "study_mean_delta_z")
    for order in ("ascending", "descending", "absolute_descending")
)
DERIVED_COLUMNS = (
    "comparison_delta_z_mean",
    "sdz_up",
    "sdz_down",
    "sdz_active_studies",
    "sdz_tied",
    "sdz_studies_total",
)


@dataclass(frozen=True)
class RollupComparisonTables:
    """Shared evidence, explicit catalog membership and all declared selection changes."""

    shared_traits: pd.DataFrame
    trait_catalog: pd.DataFrame
    rank_transitions: pd.DataFrame


@dataclass(frozen=True)
class _View:
    name: str
    quantity: str
    status_column: str = ""
    direction: str = "not_applicable"
    top: int = 0
    minimum_active_studies: int = 0


VIEWS = (
    _View("noel_by_0_05", "noel_q_by_nominal", status_column="noel_status"),
    _View("study_t_by_0_05", "study_t_q_by", status_column="study_t_status"),
    *(
        _View(
            f"comparison_mean_{direction}_top{top}_{support}",
            "experiment_mean_delta_z",
            direction=direction,
            top=top,
            minimum_active_studies=minimum,
        )
        for support, minimum in (("all", 0), ("active50", MINIMUM_ACTIVE_STUDIES))
        for direction in ("positive", "negative")
        for top in (5, 10)
    ),
)
VIEW_NAMES = tuple(view.name for view in VIEWS)


def _identities(values: pd.Index | pd.Series, description: str) -> None:
    if any(not isinstance(value, str) or not value.strip() for value in values):
        raise ValueError(f"Blank or invalid {description}")
    if not values.is_unique:
        raise ValueError(f"Duplicate {description}")


def _scope(evidence: pd.DataFrame, studies: pd.DataFrame) -> pd.DataFrame:
    required = {
        "trait",
        "trait_label",
        "composition_family_size",
        "noel_family_size",
        "noel_status",
        "noel_q_by_nominal",
        "study_t_status",
        "study_t_q_by",
        "experiment_mean_delta_z",
        "study_mean_delta_z",
        *RANK_COLUMNS,
    }
    if evidence.empty or not evidence.columns.is_unique or not required.issubset(evidence):
        raise ValueError("Missing or duplicate parent evidence columns")
    if set(DERIVED_COLUMNS).intersection(evidence.columns):
        raise ValueError("Parent evidence collides with derived comparison columns")
    _identities(evidence.trait, "parent trait identity")
    if any(not isinstance(value, str) or not value.strip() for value in evidence.trait_label):
        raise ValueError("Blank or invalid scope-specific trait label")
    result = evidence.set_index("trait").copy(deep=True).sort_index()
    for name in ("composition_family_size", "noel_family_size"):
        if not result[name].map(float).eq(len(result)).all():
            raise ValueError("Parent family size differs from the full scope catalog")
    if (
        studies.empty
        or isinstance(studies.index, pd.MultiIndex)
        or studies.index.name != "study_id"
    ):
        raise ValueError("Study ledger requires a nonempty study_id index")
    _identities(studies.index, "study identity")
    if not studies.columns.is_unique or set(studies.columns) != set(result.index):
        raise ValueError("Study ledger trait catalog differs from the complete evidence family")
    if any(
        not pd.api.types.is_numeric_dtype(dtype)
        or pd.api.types.is_bool_dtype(dtype)
        or pd.api.types.is_complex_dtype(dtype)
        for dtype in studies.dtypes
    ):
        raise ValueError("Study ledger requires finite numeric contrasts")
    values = studies.loc[:, result.index].sort_index()
    if not np.isfinite(values.to_numpy(dtype=float)).all():
        raise ValueError("Study ledger requires finite numeric contrasts")
    for name in ("study_mean_delta_z", "experiment_mean_delta_z"):
        if not np.isfinite(result[name].map(float).to_numpy(dtype=float)).all():
            raise ValueError("Parent means must be finite")
    if not np.allclose(
        result.study_mean_delta_z.to_numpy(dtype=float),
        values.mean().to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-14,
    ):
        raise ValueError("Parent study mean differs from the supplied study ledger")
    result["comparison_delta_z_mean"] = result.experiment_mean_delta_z
    result["sdz_up"] = values.gt(ACTIVITY_EPSILON).sum()
    result["sdz_down"] = values.lt(-ACTIVITY_EPSILON).sum()
    result["sdz_active_studies"] = result.sdz_up + result.sdz_down
    result["sdz_tied"] = len(values) - result.sdz_active_studies
    result["sdz_studies_total"] = len(values)
    return result


def _members(data: pd.DataFrame, view: _View) -> pd.DataFrame:
    result = pd.DataFrame(index=data.index)
    result["available"] = True
    result["selected"] = False
    result["status"] = "not_selected"
    result["view_rank"] = pd.Series(pd.NA, index=data.index, dtype="Int64")
    result["selection_value"] = np.nan
    result["retained_family_size"] = len(data)
    result["active_studies"] = data.sdz_active_studies
    if view.status_column:
        eligible = data[view.status_column].eq("tested")
        # Python float preserves round-trip decimal values at the BY boundary (#1106).
        probabilities = data.loc[eligible, view.quantity].map(float)
        if not (
            np.isfinite(probabilities.to_numpy(dtype=float)).all()
            and probabilities.between(0, 1).all()
        ):
            raise ValueError("Tested parent BY values must be finite probabilities")
        result.loc[~eligible, "status"] = "untestable"
        result.loc[eligible, "selection_value"] = probabilities
        result.loc[eligible, "selected"] = probabilities.le(BY_THRESHOLD)
        result["selection_family_size"] = len(data)
    else:
        means = data[view.quantity].map(float)
        eligible = means.gt(0) if view.direction == "positive" else means.lt(0)
        eligible &= data.sdz_active_studies.ge(view.minimum_active_studies)
        ranks = (
            means.loc[eligible]
            .rank(method="min", ascending=view.direction == "negative")
            .astype("Int64")
        )
        result.loc[~eligible, "status"] = "ineligible"
        result.loc[eligible, "view_rank"] = ranks
        result.loc[eligible, "selected"] = ranks.le(view.top)
        result["selection_value"] = means
        result["selection_family_size"] = int(eligible.sum())
    result["eligible"] = eligible
    result.loc[result.selected, "status"] = "selected"
    return result


def _transitions(off: pd.DataFrame, phenotype: pd.DataFrame, catalog: pd.DataFrame) -> pd.DataFrame:
    frames = []
    union = pd.Index(catalog.trait, name="trait")
    for view in VIEWS:
        selections = {
            scope: _members(data, view).reindex(union)
            for scope, data in (("off", off), ("phenotype", phenotype))
        }
        for selected in selections.values():
            selected["available"] = selected.available.eq(True)
            selected["selected"] = selected.selected.astype("boolean")
            selected["eligible"] = selected.eligible.astype("boolean")
            selected.loc[~selected.available, "status"] = "unavailable"
            for column in (
                "view_rank",
                "retained_family_size",
                "selection_family_size",
                "active_studies",
            ):
                selected[column] = selected[column].astype("Int64")
        frame = (
            selections["off"]
            .add_prefix("off_")
            .join(selections["phenotype"].add_prefix("phenotype_"))
        )
        old = frame.off_selected.fillna(False).to_numpy(dtype=bool)
        new = frame.phenotype_selected.fillna(False).to_numpy(dtype=bool)
        frame["membership_state"] = np.select(
            [old & new, ~old & new, old & ~new], ["retained", "entered", "left"], default="neither"
        )
        frame["view"] = view.name
        frame["selection_quantity"] = view.quantity
        frame["direction"] = view.direction
        frame["cutoff"] = view.top if view.top else BY_THRESHOLD
        frame["minimum_active_studies"] = view.minimum_active_studies
        frame["tie_policy"] = "minimum_rank_includes_ties" if view.top else "not_applicable"
        frames.append(frame.reset_index())
    return pd.concat(frames, ignore_index=True)


def build_tables(
    off_evidence: pd.DataFrame,
    phenotype_evidence: pd.DataFrame,
    off_studies: pd.DataFrame,
    phenotype_studies: pd.DataFrame,
) -> RollupComparisonTables:
    """Retain each full family, then compare identities and declared selections.

    Parent execution identity and deeper native-statistic contracts are checked by
    the shared verifier. This pure boundary additionally rejects lossy joins and
    inconsistent study support. Missing scope membership remains unavailable.
    """
    off, phenotype = (
        _scope(off_evidence, off_studies),
        _scope(phenotype_evidence, phenotype_studies),
    )
    if set(off_studies.index) != set(phenotype_studies.index):
        raise ValueError("Scope study identities differ")
    union = sorted(set(off.index) | set(phenotype.index))
    catalog = pd.DataFrame({"trait": union})
    catalog["in_off"] = catalog.trait.isin(off.index)
    catalog["in_phenotype"] = catalog.trait.isin(phenotype.index)
    catalog["off_trait_label"] = catalog.trait.map(off.trait_label)
    catalog["phenotype_trait_label"] = catalog.trait.map(phenotype.trait_label)
    catalog["membership"] = np.select(
        [catalog.in_off & catalog.in_phenotype, catalog.in_off],
        ["shared", "off_only"],
        default="phenotype_only",
    )
    shared = sorted(set(off.index) & set(phenotype.index))
    joined = off.loc[shared].add_prefix("off_").join(phenotype.loc[shared].add_prefix("phenotype_"))
    return RollupComparisonTables(
        joined.reset_index(), catalog, _transitions(off, phenotype, catalog)
    )


__all__ = ["VIEW_NAMES", "RollupComparisonTables", "build_tables"]
