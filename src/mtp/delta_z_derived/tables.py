"""Derive rows 1/2/3/5 from supplied Delta-z ledgers (#1052, 2026-10-25).

Inputs are published pair/study Delta-z matrices, their per-trait evidence and
experiment metadata. No fractions, logarithms or z scores are recomputed here.
Study votes use a strict 1e-12 threshold; oxygen quadrants use exact zero. Stool
means the exact experiments.body_site value 'Feces'. Inferential quantities here
are nominal comparison-level summaries, not study-independent headline tests.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy import stats

FACULTATIVE = "METPO:1000608"
OBLIGATE = "METPO:1000607"
ACTIVE_THRESHOLD = 1e-12
MEAN_RTOL = 1e-12
MEAN_ATOL = 1e-14
PAIR_KEYS = ["study_id", "experiment_id"]
METADATA_FIELDS = [*PAIR_KEYS, "body_site", "condition"]
CELL_COLUMNS = [
    "body_site",
    "condition",
    "trait",
    "trait_label",
    "mean_delta_z",
    "n_comparisons",
    "n_studies",
]
STOOL_COLUMNS = [
    *CELL_COLUMNS,
    "disease_rank",
    "wilcoxon_statistic",
    "wilcoxon_p",
    "wilcoxon_status",
    "wilcoxon_nonzero_comparisons",
    "wilcoxon_zero_comparisons",
    "wilcoxon_zero_method",
    "wilcoxon_method",
    "test_unit",
    "p_adjustment",
]


@dataclass(frozen=True)
class DerivedTables:
    """Descriptive rows and explicitly nominal tests from immutable parent ledgers."""

    trait_summary: pd.DataFrame
    quadrant_pairs: pd.DataFrame
    quadrant_summary: pd.DataFrame
    stool_disease_shifts: pd.DataFrame
    site_disease_shifts: pd.DataFrame


def _require_columns(frame: pd.DataFrame, columns: Sequence[str], description: str) -> None:
    if not frame.columns.is_unique:
        raise ValueError(f"Duplicate {description} columns")
    if not set(columns).issubset(frame.columns):
        raise ValueError(f"Missing required {description} columns")


def _identifiers(values: Sequence[object], description: str) -> None:
    if any(not isinstance(value, str) or not value.strip() for value in values):
        raise ValueError(f"Blank or invalid {description} identity")


def _values(frame: pd.DataFrame, description: str) -> NDArray[np.float64]:
    if any(
        not pd.api.types.is_numeric_dtype(dtype)
        or pd.api.types.is_bool_dtype(dtype)
        or pd.api.types.is_complex_dtype(dtype)
        for dtype in frame.dtypes
    ):
        raise ValueError(f"{description} must contain finite numeric values")
    try:
        values = frame.to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{description} must contain finite numeric values") from error
    if not np.isfinite(values).all():
        raise ValueError(f"{description} must contain finite numeric values")
    return values


def _matrix(frame: pd.DataFrame, *, paired: bool) -> pd.DataFrame:
    expected = PAIR_KEYS if paired else ["study_id"]
    if list(frame.index.names) != expected or isinstance(frame.index, pd.MultiIndex) != paired:
        raise ValueError(f"Delta-z identity requires index names {expected}")
    if not frame.index.is_unique or not frame.columns.is_unique:
        raise ValueError("Duplicate Delta-z identity or trait columns")
    if frame.empty:
        raise ValueError("Delta-z ledgers must contain studies and traits")
    for name in expected:
        _identifiers(list(frame.index.get_level_values(name)), f"Delta-z {name}")
    _identifiers(list(frame.columns), "trait")
    _values(frame, "Delta-z ledger")
    return frame.astype(float).sort_index()


def _evidence(frame: pd.DataFrame, traits: Sequence[str]) -> pd.DataFrame:
    required = [
        "trait",
        "trait_label",
        "experiment_mean_delta_z",
        "study_mean_delta_z",
        "noel_studies_up",
        "noel_studies_down",
        "noel_studies_tied",
    ]
    _require_columns(frame, required, "per-trait evidence")
    _identifiers(list(frame.trait), "evidence trait")
    if frame.trait.duplicated().any():
        raise ValueError("Duplicate per-trait evidence identity")
    if set(frame.trait) != set(traits):
        raise ValueError("Evidence trait catalog differs from Delta-z ledger")
    _identifiers(list(frame.trait_label), "trait label")
    _values(frame[["experiment_mean_delta_z", "study_mean_delta_z"]], "evidence means")
    return frame.copy(deep=True).reset_index(drop=True)


def _check_mean(observed: pd.Series, expected: pd.Series, description: str) -> None:
    if not np.allclose(
        observed.to_numpy(dtype=float),
        expected.to_numpy(dtype=float),
        rtol=MEAN_RTOL,
        atol=MEAN_ATOL,
    ):
        raise ValueError(f"{description} mean differs from supplied Delta-z ledger")


def trait_summary(study_delta_z: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    """Append study-mean Delta-z votes beside unchanged native raw-fraction votes.

    This row-1-only interface needs neither oxygen traits nor a headline list,
    allowing a retained-catalog summary under off as well as phenotype rollup.
    """
    studies = _matrix(study_delta_z, paired=False)
    result = _evidence(evidence, list(studies.columns))
    traits = list(result.trait)
    _check_mean(result.study_mean_delta_z, studies[traits].mean(), "Evidence study")
    values = studies[traits]
    up, down = values.gt(ACTIVE_THRESHOLD).sum(), values.lt(-ACTIVE_THRESHOLD).sum()
    result["sdz_up"] = up.to_numpy(dtype=int)
    result["sdz_down"] = down.to_numpy(dtype=int)
    result["sdz_active_studies"] = (up + down).to_numpy(dtype=int)
    result["sdz_tied"] = len(studies) - result.sdz_active_studies
    result["sdz_studies_total"] = len(studies)
    result["sdz_active_threshold"] = ACTIVE_THRESHOLD
    result["sdz_active_at_least_50"] = result.sdz_active_studies.ge(50)
    result["sdz_quantity"] = "sign_of_study_mean_delta_z"
    result["noel_quantity"] = "native_raw_fraction_direction"
    return result


def _metadata(experiments: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    _require_columns(experiments, METADATA_FIELDS, "experiment metadata")
    for key in PAIR_KEYS:
        _identifiers(list(experiments[key]), f"metadata {key}")
    if experiments.duplicated(PAIR_KEYS).any():
        raise ValueError("Duplicate experiment metadata identity")
    metadata = experiments[METADATA_FIELDS].copy()
    for column in ("body_site", "condition"):
        metadata[column] = metadata[column].fillna("")
        if any(not isinstance(value, str) for value in metadata[column]):
            raise ValueError(f"Experiment metadata {column} must be text or blank")
    metadata = metadata.set_index(PAIR_KEYS)
    pair_index = pairs.index
    if not pair_index.isin(metadata.index).all():
        raise ValueError("Experiment metadata is missing a paired identity")
    return metadata.reindex(pair_index)


def _validate_pairs(pairs: pd.DataFrame, studies: pd.DataFrame, evidence: pd.DataFrame) -> None:
    if set(pairs.columns) != set(studies.columns):
        raise ValueError("Pair and study trait catalogs differ")
    means = pairs.groupby(level="study_id").mean().reindex(columns=studies.columns)
    if set(means.index) != set(studies.index):
        raise ValueError("Pair and study identity sets differ")
    if not np.allclose(
        means.to_numpy(dtype=float),
        studies.reindex(means.index).to_numpy(dtype=float),
        rtol=MEAN_RTOL,
        atol=MEAN_ATOL,
    ):
        raise ValueError("Supplied study means differ from pair Delta-z ledger")
    _check_mean(
        evidence.experiment_mean_delta_z, pairs[list(evidence.trait)].mean(), "Evidence comparison"
    )


def _quadrant_rows(pairs: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    for trait in (FACULTATIVE, OBLIGATE):
        if trait not in pairs.columns:
            raise ValueError(f"Required oxygen trait absent from catalog: {trait}")
    result = metadata.reset_index()
    result["fa_trait"] = FACULTATIVE
    result["oa_trait"] = OBLIGATE
    result["fa_delta_z"] = pairs[FACULTATIVE].to_numpy()
    result["oa_delta_z"] = pairs[OBLIGATE].to_numpy()
    fa, oa = result.fa_delta_z, result.oa_delta_z
    result["quadrant"] = np.select(
        [(fa > 0) & (oa < 0), (fa < 0) & (oa > 0), (fa > 0) & (oa > 0), (fa < 0) & (oa < 0)],
        ["fa_up_oa_down", "fa_down_oa_up", "both_up", "both_down"],
        default="axis",
    )
    result["is_stool"] = result.body_site.eq("Feces")
    return result


def _pearson(frame: pd.DataFrame) -> dict[str, object]:
    if len(frame) < 2:
        status = "insufficient_support"
    elif frame.fa_delta_z.nunique() == 1 or frame.oa_delta_z.nunique() == 1:
        status = "constant_input"
    else:
        result = stats.pearsonr(frame.fa_delta_z.to_numpy(), frame.oa_delta_z.to_numpy())
        return {
            "pearson_r": float(result.statistic),
            "pearson_p": float(result.pvalue),
            "pearson_status": "tested",
        }
    return {"pearson_r": np.nan, "pearson_p": np.nan, "pearson_status": status}


def _quadrant_summary(per_pair: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for population, frame in (("all", per_pair), ("stool", per_pair.loc[per_pair.is_stool])):
        counts = frame.quadrant.value_counts()
        q4 = int(counts.get("fa_up_oa_down", 0))
        q2 = int(counts.get("fa_down_oa_up", 0))
        rows.append(
            {
                "population": population,
                "fa_trait": FACULTATIVE,
                "oa_trait": OBLIGATE,
                "fa_label": "facultatively aerobic",
                "oa_label": "obligately anaerobic",
                "fa_up_oa_down": q4,
                "fa_down_oa_up": q2,
                "both_up": int(counts.get("both_up", 0)),
                "both_down": int(counts.get("both_down", 0)),
                "axis": int(counts.get("axis", 0)),
                "n_comparisons": len(frame),
                "n_studies": int(frame.study_id.nunique()),
                "sign_test_n": q4 + q2,
                "sign_test_p": float(stats.binomtest(q4, q4 + q2, alternative="two-sided").pvalue)
                if q4 + q2
                else np.nan,
                "sign_test_status": "tested" if q4 + q2 else "no_discordant_pairs",
                "sign_threshold": 0.0,
                "test_unit": "comparison",
                "p_adjustment": "none_nominal",
                "stool_field": "experiments.body_site",
                "stool_value": "Feces",
                **_pearson(frame),
            }
        )
    return pd.DataFrame.from_records(rows)


def _cell(
    pairs: pd.DataFrame,
    group: pd.DataFrame,
    trait: str,
    label: str,
    body_site: str,
    condition: str,
) -> dict[str, object]:
    return {
        "body_site": body_site,
        "condition": condition,
        "trait": trait,
        "trait_label": label,
        "mean_delta_z": float(pairs.loc[group.index, trait].mean()),
        "n_comparisons": len(group),
        "n_studies": int(group.index.get_level_values("study_id").nunique()),
    }


def _wilcoxon(values: NDArray[np.float64]) -> dict[str, object]:
    nonzero = int(np.count_nonzero(values))
    statistic, probability = np.nan, np.nan
    if nonzero:
        result = stats.wilcoxon(
            values, zero_method="wilcox", alternative="two-sided", correction=False, method="auto"
        )
        statistic, probability = float(result.statistic), float(result.pvalue)
    return {
        "wilcoxon_statistic": statistic,
        "wilcoxon_p": probability,
        "wilcoxon_status": "tested" if nonzero else "all_tied",
        "wilcoxon_nonzero_comparisons": nonzero,
        "wilcoxon_zero_comparisons": len(values) - nonzero,
        "wilcoxon_zero_method": "wilcox",
        "wilcoxon_method": "auto",
        "test_unit": "comparison",
        "p_adjustment": "none_nominal",
    }


def _stool_shifts(
    pairs: pd.DataFrame,
    metadata: pd.DataFrame,
    labels: dict[str, str],
    traits: Sequence[str],
    top_diseases: int,
) -> pd.DataFrame:
    stool = metadata.loc[metadata.body_site.eq("Feces")]
    groups = {
        str(condition): group
        for condition, group in stool.groupby("condition", sort=True, dropna=False)
    }
    selected = sorted(groups, key=lambda condition: (-len(groups[condition]), condition))[
        :top_diseases
    ]
    rows: list[dict[str, object]] = []
    for rank, condition in enumerate(selected, start=1):
        group = groups[condition]
        for trait in traits:
            rows.append(
                {
                    **_cell(pairs, group, trait, labels[trait], "Feces", condition),
                    "disease_rank": rank,
                    **_wilcoxon(pairs.loc[group.index, trait].to_numpy(dtype=float)),
                }
            )
    return pd.DataFrame.from_records(rows, columns=STOOL_COLUMNS)


def _site_shifts(
    pairs: pd.DataFrame, metadata: pd.DataFrame, labels: dict[str, str], traits: Sequence[str]
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (body_site, condition), group in metadata.groupby(
        ["body_site", "condition"], sort=True, dropna=False
    ):
        rows.extend(
            _cell(pairs, group, trait, labels[trait], str(body_site), str(condition))
            for trait in traits
        )
    return pd.DataFrame.from_records(rows, columns=CELL_COLUMNS)


def build_tables(
    pair_delta_z: pd.DataFrame,
    study_delta_z: pd.DataFrame,
    evidence: pd.DataFrame,
    experiments: pd.DataFrame,
    headline_traits: Sequence[str],
    *,
    top_diseases: int = 7,
) -> DerivedTables:
    """Build nominal derived analyses without recomputing or editing supplied Delta-z.

    Metadata may contain additional unpaired experiments, but may not omit a pair
    or contain duplicate identities. Blank body-site/condition cells are retained.
    Disease selection is comparison count descending, then exact condition name
    ascending. Every declared headline trait must be retained; none is substituted.
    """
    pairs, studies = _matrix(pair_delta_z, paired=True), _matrix(study_delta_z, paired=False)
    summary = trait_summary(studies, evidence)
    _validate_pairs(pairs, studies, summary)
    traits = tuple(headline_traits)
    if not traits or len(set(traits)) != len(traits):
        raise ValueError("Headline trait list must be nonempty and unique")
    _identifiers(list(traits), "headline trait")
    if not set(traits).issubset(pairs.columns):
        raise ValueError("Declared headline trait is absent from retained catalog")
    if isinstance(top_diseases, bool) or not isinstance(top_diseases, int) or top_diseases < 1:
        raise ValueError("The number of disease groups must be a positive integer")
    metadata = _metadata(experiments, pairs)
    labels = dict(zip(summary.trait, summary.trait_label, strict=True))
    per_pair = _quadrant_rows(pairs, metadata)
    return DerivedTables(
        trait_summary=summary,
        quadrant_pairs=per_pair,
        quadrant_summary=_quadrant_summary(per_pair),
        stool_disease_shifts=_stool_shifts(pairs, metadata, labels, traits, top_diseases),
        site_disease_shifts=_site_shifts(pairs, metadata, labels, list(summary.trait)),
    )
