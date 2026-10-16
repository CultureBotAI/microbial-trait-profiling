#!/usr/bin/env python3
"""Export ranked enrichment trait summaries from a trait enrichment report.

Both outputs always carry their declared column headers plus ``# Key: Value``
metadata recording the result status and counts, so an empty table is
interpretable (#296):

* ``Result-Status``: ``discoveries``, ``zero_discoveries`` or
  ``zero_eligible_tests`` (see :mod:`bsdbng.enrichment_results`).
* ``Total-Tests-Performed`` / ``Significant-Hits``: counts from the analyzer.
* ``Exported-Rows`` (pair table) / ``Ranked-Traits`` (overall table).

Malformed or missing input exits with status 1 and an actionable message
instead of writing an empty successful report.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any

import pandas as pd

from mtp.enrichment_results import (
    META_RESULT_STATUS,
    META_SIGNIFICANT_HITS,
    META_TESTS_PERFORMED,
    RESULT_STATUSES,
    STATUS_DISCOVERIES,
    MalformedResultsError,
    load_results_tsv,
)
from mtp.paths import DERIVED_DATA_DIR

DEFAULT_INPUT_CANDIDATES = [
    DERIVED_DATA_DIR / "trait_enrichment_report_p0.001.txt",
    DERIVED_DATA_DIR / "trait_enrichment_report.txt",
]
DEFAULT_RESULTS_CANDIDATES = [
    DERIVED_DATA_DIR / "trait_enrichment_results_p0.001.tsv",
    DERIVED_DATA_DIR / "trait_enrichment_results.tsv",
]
DEFAULT_OVERALL_OUTPUT_STEM = "trait_enrichment_traits_overall_ranked"
DEFAULT_PAIR_OUTPUT_STEM = "trait_enrichment_trait_study_pairs_ranked"
DEFAULT_MIN_FOLD = 1.5
DEFAULT_MAX_FOLD = 0.67
# #115: traits must be significant in at least this many profiles to be ranked.
MIN_HIT_ROWS_FOR_RANKING = 2


@dataclass(frozen=True)
class ParsedRow:
    study_id: str
    experiment_id: str
    inc_dec: str
    profile_id: str
    trait: str
    p_value: float
    fdr_q: float
    fold: float
    obs: float
    exp: float
    testable_profiles: int
    hit_count: int | None = None
    trait_total: int = 0
    is_significant: bool = True


PARSED_ROW_COLUMNS: list[str] = [f.name for f in fields(ParsedRow)]


def overall_columns(
    stat_type: str, censored: bool, body_site_group: str | None = None
) -> list[str]:
    """Declared column order of the overall trait ranking table."""
    is_fdr = stat_type == "fdr_q"
    best_stat_col = "best_qvalue" if is_fdr else "best_pvalue"
    median_stat_col = "median_qvalue" if is_fdr else "median_pvalue"
    other_best_col = "best_pvalue" if is_fdr else "best_qvalue"
    columns = [
        "rank",
        "trait",
        "net_direction",
        "inc_hits",
        "dec_hits",
        "prevalence_summary",
        "stat_type",
        best_stat_col,
        median_stat_col,
        "hit_rows",
        "unique_profiles",
        "mean_fold",
        "max_fold",
        "testable_profiles",
        "n_profiles_tested",
        "trait_total",
        other_best_col,
        "min_p_direction",
        "min_p_study_id",
        "min_p_experiment_id",
        "min_p_profile_id",
    ]
    if censored:
        columns.append("n_censored_profiles")
    columns.append("prevalence")
    if body_site_group is not None:
        columns.append("body_site_group")
    return columns


def pair_columns(
    stat_type: str,
    body_site_group: str | None = None,
    *,
    unfiltered: bool = False,
) -> list[str]:
    """Declared column order of the trait-study pair ranking table."""
    renamed = [
        "q_value" if (stat_type == "fdr_q" and c == "fdr_q") else c
        for c in PARSED_ROW_COLUMNS
        if c != "trait_total" and (unfiltered or c != "is_significant")
    ]
    cols = ["rank", *renamed, "is_censored_hit"]
    if body_site_group is not None:
        cols.append("body_site_group")
    return cols


def resolve_default_input() -> Path:
    for candidate in DEFAULT_INPUT_CANDIDATES:
        if candidate.exists():
            return candidate
    return DEFAULT_INPUT_CANDIDATES[-1]


def resolve_default_results(report_path: Path | None = None) -> Path:
    if report_path:
        if report_path.suffix == ".tsv":
            return report_path
        # Try to find a .tsv with same name/prefix as report
        stem = report_path.stem
        if stem.startswith("trait_enrichment_report"):
            tsv_stem = stem.replace("report", "results")
            candidate = report_path.parent / f"{tsv_stem}.tsv"
            if candidate.exists():
                return candidate

    for candidate in DEFAULT_RESULTS_CANDIDATES:
        if candidate.exists():
            return candidate
    return DEFAULT_RESULTS_CANDIDATES[-1]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export ranked trait enrichment TSV outputs from a text report."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=resolve_default_input(),
        help="Input trait enrichment report path (used for header detection).",
    )
    parser.add_argument(
        "--results",
        type=Path,
        help=(
            "Input trait enrichment results TSV (used for data aggregation). "
            "Defaults to companion of --input."
        ),
    )
    parser.add_argument(
        "--overall-output",
        type=Path,
        help="Output TSV for overall trait ranking. Defaults to dynamic path based on stats.",
    )
    parser.add_argument(
        "--pair-output",
        type=Path,
        help="Output TSV for trait-study pair ranking. Defaults to dynamic path based on stats.",
    )
    parser.add_argument(
        "--body-site-group",
        type=str,
        default=None,
        help="Optional body site group / stratum label to record in headers and columns.",
    )
    parser.add_argument(
        "--profiles",
        type=Path,
        default=None,
        help="Optional path to profiles TSV to compute empirical background K (trait totals).",
    )
    parser.add_argument(
        "--tail",
        choices=["upper", "lower"],
        default=None,
        help=(
            "Direction of tail: upper for enrichment, lower for depletion "
            "(detected from metadata if omitted)."
        ),
    )
    parser.add_argument(
        "--unfiltered",
        "--all-results",
        action="store_true",
        dest="unfiltered",
        help="Process an unfiltered results TSV containing all tested pairs.",
    )
    parser.add_argument(
        "--fdr-scope",
        choices=["global", "profile"],
        default=None,
        help="Explicit FDR multiplicity scope ('global' or 'profile'; overrides metadata if set).",
    )
    return parser.parse_args(argv)


def compute_trait_totals_from_profiles(profiles_path: Path) -> dict[str, int]:
    """Compute empirical background sum K for each trait from a profiles TSV."""
    if not profiles_path.exists():
        return {}
    try:
        df = pd.read_csv(profiles_path, sep="\t", low_memory=False)
        if "study_id" in df.columns and "experiment_id" in df.columns and "inc_dec" in df.columns:
            exp_dirs = df.groupby(["study_id", "experiment_id"])["inc_dec"].apply(
                lambda s: set(s.str.lower())
            )
            paired = {k for k, dirs in exp_dirs.items() if "inc" in dirs and "dec" in dirs}
            mask = df.set_index(["study_id", "experiment_id"]).index.isin(paired)
            df = df[mask]

        non_trait = {"study_id", "experiment_id", "inc_dec", "taxa_count", "body_site"}
        trait_cols = [c for c in df.columns if c not in non_trait and "mediadive.medium:" not in c]
        if not trait_cols:
            return {}
        numeric_df = df[trait_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        return {col: int(total) for col, total in numeric_df.sum().items()}
    except Exception as exc:
        print(
            f"Warning: failed to compute trait totals from {profiles_path}: {exc}",
            file=sys.stderr,
        )
        return {}


def _parse_summary_values(summary: str) -> tuple[float | None, float | None, float | None]:
    obs_match = re.search(r"obs=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", summary)
    exp_match = re.search(r"exp=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", summary)
    fold_match = re.search(r"fold=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", summary)
    obs = float(obs_match.group(1)) if obs_match else None
    exp = float(exp_match.group(1)) if exp_match else None
    fold = float(fold_match.group(1)) if fold_match else None
    return obs, exp, fold


def detect_metadata(path: Path) -> dict[str, Any]:
    """Detect metadata from the report or results TSV."""
    meta: dict[str, Any] = {
        "stat_type": "p_value",
        "threshold": 0.001,
        "min_fold": None,
        "max_fold": None,
        "tail": "upper",
        "testable_profiles": None,
        "max_traits": None,
        "fdr_scope": "global",
        "tests_performed": None,
        "significant_hits": None,
        "result_status": None,
    }
    if not path.exists():
        return meta

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            # Text report format
            if line.startswith("p <= "):
                meta["stat_type"] = "p_value"
                val = line.split("<=")[1]
                if "(Depletion)" in val:
                    meta["tail"] = "lower"
                    val = val.replace("(Depletion)", "")
                meta["threshold"] = float(val.strip())
            elif line.startswith("p < "):
                meta["stat_type"] = "p_value"
                val = line.split("<")[1]
                if "(Depletion)" in val:
                    meta["tail"] = "lower"
                    val = val.replace("(Depletion)", "")
                meta["threshold"] = float(val.strip())
            elif line.startswith("FDR q <= "):
                meta["stat_type"] = "fdr_q"
                val = line.split("<=")[1]
                if "(Depletion)" in val:
                    meta["tail"] = "lower"
                    val = val.replace("(Depletion)", "")
                meta["threshold"] = float(val.strip())
            elif line.startswith("FDR q < "):
                meta["stat_type"] = "fdr_q"
                val = line.split("<")[1]
                if "(Depletion)" in val:
                    meta["tail"] = "lower"
                    val = val.replace("(Depletion)", "")
                meta["threshold"] = float(val.strip())
            elif line.startswith("Min Fold >= "):
                meta["min_fold"] = float(line.split(">=")[1].strip())
            elif line.startswith("Max Fold <= "):
                meta["max_fold"] = float(line.split("<=")[1].strip())
            elif line.startswith("Tail: "):
                meta["tail"] = line.split(": ")[1].strip().lower()
            elif line.startswith("Tests performed: "):
                meta["tests_performed"] = int(line.split(": ")[1].strip())
            elif line.startswith("Significant hits: "):
                meta["significant_hits"] = int(line.split(": ")[1].strip())
            elif line.startswith("Result status: "):
                meta["result_status"] = line.split(": ")[1].strip()
            # TSV metadata header format
            elif line.startswith("# Tail: "):
                meta["tail"] = line.split(": ")[1].strip().lower()
            elif line.startswith("# FDR: "):
                is_fdr = line.split(": ")[1].strip().lower() == "true"
                meta["stat_type"] = "fdr_q" if is_fdr else "p_value"
            elif line.startswith(("# P-Threshold: ", "# Stat-Threshold: ")):
                meta["threshold"] = float(line.split(": ")[1].strip())
            elif line.startswith("# Min-Fold: "):
                meta["min_fold"] = float(line.split(": ")[1].strip())
            elif line.startswith("# Max-Fold: "):
                meta["max_fold"] = float(line.split(": ")[1].strip())
            elif line.startswith("# FDR-Scope: "):
                meta["fdr_scope"] = line.split(": ")[1].strip().lower()
            elif line.startswith("# Testable-Profiles: "):
                meta["testable_profiles"] = int(line.split(": ")[1].strip())
            elif line.startswith("# Max-Traits-Per-Row-Report: "):
                val = line.split(": ")[1].strip()
                meta["max_traits"] = None if val in ("0", "-1") else int(val)
            elif line.startswith("Max traits per profile: "):
                val = line.split(": ")[1].strip()
                meta["max_traits"] = None if val == "unlimited" else int(val)

    return meta


def parse_results_tsv(
    path: Path,
    stat_type: str,
    threshold: float,
    min_fold: float = DEFAULT_MIN_FOLD,
    max_fold: float = DEFAULT_MAX_FOLD,
    tail: str = "upper",
    trait_totals: dict[str, int] | None = None,
    unfiltered: bool = False,
) -> list[ParsedRow]:
    """Validate the analyzer results TSV and return the selected rows.

    Raises ``FileNotFoundError`` or :class:`MalformedResultsError` for invalid
    input. A valid file with zero selected rows returns an empty list.
    """
    table = load_results_tsv(path, allow_unfiltered=unfiltered)
    df = table.frame

    # col used for filtering
    col = "p_value" if stat_type == "p_value" else "fdr_q"

    if "is_significant" in df.columns:
        sig_mask = df["is_significant"].astype(str).str.lower().isin(("true", "1"))
    else:
        if tail == "lower":
            sig_mask = (df[col] <= threshold) & (df["fold"] <= max_fold)
        else:
            sig_mask = (df[col] <= threshold) & (df["fold"] >= min_fold)

    fold_desc = f"max_fold={max_fold}" if tail == "lower" else f"min_fold={min_fold}"
    if unfiltered:
        print(
            f"Loaded {len(df)} all-result rows ({sig_mask.sum()} significant hits) "
            f"(source={path.name}, threshold={threshold}, {fold_desc}, tail={tail})"
        )
    else:
        df = df[sig_mask]
        print(
            f"Filtered to {len(df)} significant hits "
            f"(source={path.name}, threshold={threshold}, {fold_desc}, tail={tail})"
        )

    rows = []
    for idx, r in df.iterrows():
        trait = str(r["trait"])
        # Exclude MediaDive growth-medium traits for consistency
        if "mediadive.medium:" in trait:
            continue

        study_id = str(r["study_id"])
        experiment_id = str(r["experiment_id"])
        inc_dec = str(r["inc_dec"])
        trait_total_raw = r["background_k"] if "background_k" in r else r.get("trait_total", 0)
        trait_total_val = int(trait_total_raw) if pd.notna(trait_total_raw) else 0
        if trait_total_val == 0 and trait_totals and trait in trait_totals:
            trait_total_val = trait_totals[trait]
        is_sig = bool(sig_mask.loc[idx]) if unfiltered else True
        rows.append(
            ParsedRow(
                study_id=study_id,
                experiment_id=experiment_id,
                inc_dec=inc_dec,
                profile_id=f"{study_id}|{experiment_id}|{inc_dec}",
                trait=trait,
                p_value=float(r["p_value"]),
                fdr_q=float(r["fdr_q"]),
                fold=float(r["fold"]),
                obs=float(r["count"]),
                exp=float(r["expected"]),
                testable_profiles=int(r["testable_profiles"]),
                trait_total=trait_total_val,
                is_significant=is_sig,
            )
        )
    return rows


def parse_report(
    path: Path,
    stat_type: str = "p_value",
    allow_empty: bool = False,
    trait_totals: dict[str, int] | None = None,
) -> list[ParsedRow]:
    rows: list[ParsedRow] = []
    current_header: tuple[str, str, str] | None = None

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if (
                not line
                or line == "---"
                or line.startswith(
                    (
                        "p < ",
                        "p <= ",
                        "FDR q < ",
                        "FDR q <= ",
                        "Min Fold >= ",
                        "Max Fold <= ",
                        "Tail: ",
                        "Tests performed: ",
                        "Significant hits: ",
                        "Result status: ",
                    )
                )
            ):
                continue

            if line.startswith("- "):
                if current_header is None:
                    raise ValueError(f"Trait row found before profile header: {line!r}")
                payload = line[2:]
                parts = payload.split(", ")
                if len(parts) < 3:
                    raise ValueError(f"Unparsable trait row: {line!r}")

                trait = parts[0].strip()
                val = float(parts[1])
                hit_count: int | None = None
                if re.fullmatch(r"\d+", parts[2]):
                    hit_count = int(parts[2])

                summary = parts[3] if len(parts) >= 4 else parts[2]
                obs, exp, fold = _parse_summary_values(summary)

                study_id, experiment_id, inc_dec = current_header
                profile_id = f"{study_id}|{experiment_id}|{inc_dec}"
                # Exclude MediaDive growth-medium traits for consistency
                if "mediadive.medium:" in trait:
                    continue

                # Set whichever statistic we found, and default the other to 1.0
                p_value = val if stat_type == "p_value" else 1.0
                fdr_q = val if stat_type == "fdr_q" else 1.0

                trait_tot = trait_totals.get(trait, 0) if trait_totals else 0
                rows.append(
                    ParsedRow(
                        study_id=study_id,
                        experiment_id=experiment_id,
                        inc_dec=inc_dec,
                        profile_id=profile_id,
                        trait=trait,
                        p_value=p_value,
                        fdr_q=fdr_q,
                        fold=fold or 0.0,
                        obs=obs or 0.0,
                        exp=exp or 0.0,
                        testable_profiles=0,  # Not available in text report
                        hit_count=hit_count,
                        trait_total=trait_tot,
                    )
                )
                continue

            header_parts = line.split("\t")
            if len(header_parts) != 3:
                continue
            current_header = (header_parts[0], header_parts[1], header_parts[2])

    if not rows and not allow_empty:
        raise MalformedResultsError(
            f"No trait rows parsed from {path}, and the report does not declare a "
            "'Result status:' of zero_discoveries or zero_eligible_tests. Regenerate it "
            "with `scripts/trait_enrichment.py`."
        )

    return rows


def rows_to_frame(rows: list[ParsedRow], max_traits: int | None) -> pd.DataFrame:
    """Build the parsed-row frame with a fixed schema (even when empty)."""
    frame = pd.DataFrame([row.__dict__ for row in rows], columns=PARSED_ROW_COLUMNS)
    # If we are parsing from a truncated source, identify potentially censored profiles.
    # A profile is censored if it hit the cap.
    if max_traits and not frame.empty:
        counts_per_profile = frame.groupby("profile_id").size()
        censored_profiles = set(counts_per_profile[counts_per_profile >= max_traits].index)
        frame["is_censored_hit"] = frame["profile_id"].isin(censored_profiles)
    else:
        frame["is_censored_hit"] = pd.Series(False, index=frame.index, dtype=bool)
    return frame


def build_overall_table(
    frame: pd.DataFrame,
    stat_type: str,
    tail: str = "upper",
    total_profiles: int | None = None,
    max_traits: int | None = None,
    fdr_scope: str = "global",
    body_site_group: str | None = None,
    trait_totals: dict[str, int] | None = None,
) -> pd.DataFrame:
    _ = total_profiles  # Reserved for future weighting
    if fdr_scope not in {"global", "profile"}:
        raise ValueError(f"Unsupported fdr_scope: {fdr_scope!r}; expected 'global' or 'profile'.")
    columns = overall_columns(stat_type, censored=bool(max_traits), body_site_group=body_site_group)
    if frame.empty:
        return pd.DataFrame(columns=columns)

    if "is_censored_hit" not in frame.columns:
        frame = frame.assign(is_censored_hit=False)

    sig_frame = frame[frame["is_significant"]] if "is_significant" in frame.columns else frame
    if sig_frame.empty:
        return pd.DataFrame(columns=columns)

    grouped = sig_frame.groupby("trait", dropna=False)

    # Use the relevant stat for choosing the 'best' row metadata
    stat_col = "fdr_q" if stat_type == "fdr_q" else "p_value"
    best_stat_col = "best_qvalue" if stat_type == "fdr_q" else "best_pvalue"
    median_stat_col = "median_qvalue" if stat_type == "fdr_q" else "median_pvalue"

    # Multiplicity control for the per-trait "best" statistic:
    # Since we are taking the minimum over all profiles for a trait, we should
    # ideally apply a correction here. For now, we keep the naming explicit.

    best_idx = grouped[stat_col].idxmin()
    best_rows = sig_frame.loc[
        best_idx, ["trait", "inc_dec", "study_id", "experiment_id", "profile_id"]
    ]
    best_rows = best_rows.rename(
        columns={
            "inc_dec": "min_p_direction",
            "study_id": "min_p_study_id",
            "experiment_id": "min_p_experiment_id",
            "profile_id": "min_p_profile_id",
        }
    )

    agg_map = {
        best_stat_col: (stat_col, "min"),
        median_stat_col: (stat_col, "median"),
        "hit_rows": ("trait", "size"),
        "unique_profiles": ("profile_id", "nunique"),
        "mean_fold": ("fold", "mean"),
        "max_fold": ("fold", "max"),
        "testable_profiles": ("testable_profiles", "first"),
    }
    if "trait_total" in sig_frame.columns:
        agg_map["trait_total"] = ("trait_total", "max")

    # Also include the other statistic if it's available (non-1.0)
    if stat_type == "fdr_q":
        agg_map["best_pvalue"] = ("p_value", "min")
    else:
        agg_map["best_qvalue"] = ("fdr_q", "min")

    summary = grouped.agg(**agg_map)

    inc_counts = (
        sig_frame[sig_frame["inc_dec"] == "inc"].groupby("trait")["trait"].size().rename("inc_hits")
    )
    dec_counts = (
        sig_frame[sig_frame["inc_dec"] == "dec"].groupby("trait")["trait"].size().rename("dec_hits")
    )

    overall = summary.join(inc_counts, how="left").join(dec_counts, how="left")
    overall[["inc_hits", "dec_hits"]] = overall[["inc_hits", "dec_hits"]].fillna(0).astype(int)
    overall = overall.reset_index().merge(best_rows, on="trait", how="left")
    overall.insert(1, "stat_type", stat_type)
    if "trait_total" not in overall.columns:
        overall["trait_total"] = 0
    else:
        overall["trait_total"] = overall["trait_total"].fillna(0).astype(int)
    if trait_totals:
        for t, tot in trait_totals.items():
            mask = overall["trait"] == t
            if mask.any():
                overall.loc[mask, "trait_total"] = tot

    if "is_significant" in frame.columns and frame["profile_id"].nunique() > 0:
        tested_per_trait = frame.groupby("trait")["profile_id"].nunique()
        overall["n_profiles_tested"] = (
            overall["trait"].map(tested_per_trait).fillna(overall["testable_profiles"]).astype(int)
        )
    else:
        overall["n_profiles_tested"] = overall["testable_profiles"]

    # Calculate net direction based on vote count
    def get_net_direction(row: pd.Series) -> str:
        if row["inc_hits"] > row["dec_hits"]:
            return "inc"
        if row["dec_hits"] > row["inc_hits"]:
            return "dec"
        return "mixed"

    overall["net_direction"] = overall.apply(get_net_direction, axis=1)

    if max_traits:
        censored_counts = (
            sig_frame[sig_frame["is_censored_hit"]]
            .groupby("trait")
            .size()
            .rename("n_censored_profiles")
        )
        overall = overall.merge(censored_counts, on="trait", how="left")
        overall["n_censored_profiles"] = overall["n_censored_profiles"].fillna(0).astype(int)

    # #119: Calculate prevalence fraction using the per-trait tested profile denominator
    # (n_profiles_tested) so it matches the per-trait ORA testing scope.
    overall["prevalence"] = overall["hit_rows"] / overall["n_profiles_tested"]

    # Add human-readable prevalence summary
    overall["prevalence_summary"] = (
        overall["hit_rows"].astype(str) + " / " + overall["n_profiles_tested"].astype(str)
    )

    # #115: Minimum hit thresholding for overall ranking to reduce noise.
    # We require traits to be significant in at least 2 profiles to appear in rankings.
    overall = overall[overall["hit_rows"] >= MIN_HIT_ROWS_FOR_RANKING]

    # #119: Rank by prevalence (significant / testable) descending, then by significance.
    # This provides a directly interpretable ranking that avoids extreme-value inflation.
    sort_stat = "best_qvalue" if stat_type == "fdr_q" else "best_pvalue"
    if tail == "lower":
        overall = overall.sort_values(
            ["prevalence", sort_stat, "mean_fold"],
            ascending=[False, True, True],
        )
    else:
        overall = overall.sort_values(
            ["prevalence", sort_stat, "max_fold"],
            ascending=[False, True, False],
        )
    overall.insert(0, "rank", range(1, len(overall) + 1))
    if body_site_group is not None:
        overall["body_site_group"] = body_site_group

    # Enforce the declared schema (raises if the computed table drifts from it).
    return overall[columns]


def build_pair_table(
    frame: pd.DataFrame,
    stat_type: str,
    tail: str = "upper",
    body_site_group: str | None = None,
    *,
    unfiltered: bool = False,
) -> pd.DataFrame:
    cols = pair_columns(stat_type, body_site_group=body_site_group, unfiltered=unfiltered)
    if "is_censored_hit" not in frame.columns:
        frame = frame.assign(is_censored_hit=False)
    if unfiltered:
        pairs = frame.copy()
    elif "is_significant" in frame.columns:
        pairs = frame[frame["is_significant"]].copy()
    else:
        pairs = frame.copy()
    if pairs.empty:
        return pd.DataFrame(columns=cols)

    sort_col = "fdr_q" if stat_type == "fdr_q" else "p_value"
    if tail == "lower":
        pairs = pairs.sort_values([sort_col, "fold"], ascending=[True, True])
    else:
        pairs = pairs.sort_values([sort_col, "fold"], ascending=[True, False])

    # Rename stat column for self-description
    if stat_type == "fdr_q":
        pairs = pairs.rename(columns={"fdr_q": "q_value"})

    pairs.insert(0, "rank", range(1, len(pairs) + 1))
    if body_site_group is not None:
        pairs["body_site_group"] = body_site_group
    return pairs[cols]


def _write_table(
    path: Path,
    table: pd.DataFrame,
    header: list[tuple[str, Any]],
    warnings: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        for key, value in header:
            f.write(f"# {key}: {value}\n")
        for warning in warnings:
            f.write(f"# WARNING: {warning}\n")
        table.to_csv(f, sep="\t", index=False)


def _fail(message: str) -> int:
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    meta = detect_metadata(args.input)
    tail = args.tail or meta.get("tail", "upper")
    stat_type = meta["stat_type"]
    stat_threshold = meta["threshold"]
    min_fold = meta["min_fold"] if meta["min_fold"] is not None else DEFAULT_MIN_FOLD
    max_fold = meta["max_fold"] if meta["max_fold"] is not None else DEFAULT_MAX_FOLD
    testable_profiles = meta["testable_profiles"]
    fdr_scope = args.fdr_scope or meta["fdr_scope"]
    max_traits = meta.get("max_traits")
    tests_performed: int | None = meta["tests_performed"]
    analyzer_hits: int | None = meta["significant_hits"]
    result_status: str | None = meta["result_status"]
    is_unfiltered: bool = bool(args.unfiltered)

    results_path = args.results or resolve_default_results(args.input)

    trait_totals: dict[str, int] = {}
    if args.profiles and args.profiles.is_file():
        trait_totals = compute_trait_totals_from_profiles(args.profiles)

    try:
        if results_path.exists():
            print(f"Loading data from results TSV: {results_path}")
            results_meta = detect_metadata(results_path)

            # Prefer TSV metadata if available
            if args.tail is None and results_meta.get("tail"):
                tail = results_meta["tail"]
            stat_type = results_meta.get("stat_type", stat_type)
            stat_threshold = results_meta.get("threshold", stat_threshold)
            if results_meta.get("min_fold") is not None:
                min_fold = results_meta["min_fold"]
            if results_meta.get("max_fold") is not None:
                max_fold = results_meta["max_fold"]
            fdr_scope = args.fdr_scope or results_meta.get("fdr_scope", fdr_scope)
            testable_profiles = results_meta.get("testable_profiles") or testable_profiles

            if args.unfiltered or "all_results" in str(results_path):
                is_unfiltered = True
            table = load_results_tsv(results_path, allow_unfiltered=is_unfiltered)
            tests_performed = table.tests_performed
            analyzer_hits = table.significant_hits
            result_status = table.status

            rows = parse_results_tsv(
                results_path,
                stat_type,
                stat_threshold,
                min_fold=min_fold,
                max_fold=max_fold,
                tail=tail,
                trait_totals=trait_totals,
                unfiltered=is_unfiltered,
            )
            # Results TSV contains all hits regardless of the report truncation cap
            is_source_truncated = False
            max_traits = None
        else:
            if not args.input.is_file():
                return _fail(
                    f"Neither the results TSV ({results_path}) nor the report ({args.input}) "
                    "exists. Run `scripts/trait_enrichment.py` (or `just analyze`) first."
                )
            print(
                f"Warning: Results TSV not found at {results_path}. "
                f"Falling back to report parser ({args.input.name})."
            )
            if result_status is not None and result_status not in RESULT_STATUSES:
                raise MalformedResultsError(
                    f"{args.input}: unknown 'Result status: {result_status}' "
                    f"(expected one of {list(RESULT_STATUSES)})."
                )
            allow_empty = result_status is not None and result_status != STATUS_DISCOVERIES
            rows = parse_report(
                args.input,
                stat_type=stat_type,
                allow_empty=allow_empty,
                trait_totals=trait_totals,
            )
            if result_status is None:
                result_status = STATUS_DISCOVERIES

            # Check for truncation in the report metadata
            if max_traits is not None:
                is_source_truncated = True
                print(
                    f"\nCRITICAL WARNING: The source report '{args.input.name}' "
                    f"is truncated to {max_traits} traits per profile."
                )
                print(
                    "Downstream rankings (hit_rows, unique_profiles, mean_fold) "
                    "will be systematically biased."
                )
                print(
                    "Please regenerate the report with --max-traits-per-row 0 "
                    "or ensure the results TSV is available.\n"
                )
            else:
                is_source_truncated = False
    except (FileNotFoundError, MalformedResultsError) as exc:
        return _fail(str(exc))

    frame = rows_to_frame(rows, max_traits)

    overall = build_overall_table(
        frame,
        stat_type,
        tail=tail,
        total_profiles=testable_profiles,
        max_traits=max_traits,
        fdr_scope=fdr_scope,
        body_site_group=args.body_site_group,
        trait_totals=trait_totals,
    )
    pairs = build_pair_table(
        frame,
        stat_type,
        tail=tail,
        body_site_group=args.body_site_group,
        unfiltered=is_unfiltered,
    )

    # Resolve dynamic filenames
    # Include 'profile' in suffix if using per-profile FDR
    method_suffix = "q" if stat_type == "fdr_q" else "p"
    if stat_type == "fdr_q" and fdr_scope == "profile":
        method_suffix = "q_profile"

    suffix = f"_{method_suffix}{stat_threshold}"

    overall_stem = (
        DEFAULT_OVERALL_OUTPUT_STEM.replace("enrichment", "depletion")
        if tail == "lower"
        else DEFAULT_OVERALL_OUTPUT_STEM
    )
    pair_stem = (
        DEFAULT_PAIR_OUTPUT_STEM.replace("enrichment", "depletion")
        if tail == "lower"
        else DEFAULT_PAIR_OUTPUT_STEM
    )

    write_overall = True
    write_pair = True
    if args.overall_output is not None and args.pair_output is None:
        write_pair = False
    elif args.pair_output is not None and args.overall_output is None:
        write_overall = False

    overall_out = args.overall_output
    if overall_out is None and write_overall:
        overall_out = DERIVED_DATA_DIR / f"{overall_stem}{suffix}.tsv"

    pair_out = args.pair_output
    if pair_out is None and write_pair:
        pair_out = DERIVED_DATA_DIR / f"{pair_stem}{suffix}.tsv"

    def na(value: int | None) -> str:
        return "NA" if value is None else str(value)

    common: list[tuple[str, Any]] = [
        ("Tail", tail),
        ("Stat-Type", stat_type),
        ("FDR-Scope", fdr_scope),
        ("Stat-Threshold", stat_threshold),
    ]
    if tail == "lower":
        common.append(("Max-Fold", max_fold))
    else:
        common.append(("Min-Fold", min_fold))
    common.extend(
        [
            (META_RESULT_STATUS, result_status),
            (META_TESTS_PERFORMED, na(tests_performed)),
            (META_SIGNIFICANT_HITS, na(analyzer_hits)),
        ]
    )
    if args.body_site_group is not None:
        common.append(("Body-Site-Group", args.body_site_group))

    overall_header = list(common)
    if testable_profiles:
        overall_header.append(("Total-Testable-Profiles", testable_profiles))
    overall_header.append(("Min-Profiles-Per-Ranked-Trait", MIN_HIT_ROWS_FOR_RANKING))
    overall_header.append(("Ranked-Traits", len(overall)))
    pair_header = [*common, ("Exported-Rows", len(pairs))]
    overall_warnings: list[str] = []
    pair_warnings: list[str] = []
    if is_source_truncated:
        overall_header.append(("Max-Traits-Per-Profile-Cap", max_traits))
        pair_header.append(("Max-Traits-Per-Profile-Cap", max_traits))
        overall_warnings.append("Totals are censored due to source truncation.")
        pair_warnings.append("Data is censored due to source truncation.")

    if write_overall and overall_out is not None:
        _write_table(overall_out, overall, overall_header, overall_warnings)
    if write_pair and pair_out is not None:
        _write_table(pair_out, pairs, pair_header, pair_warnings)

    print(f"Result status: {result_status}")
    print(f"Tests performed: {na(tests_performed)}")
    print(f"Analyzer significant hits: {na(analyzer_hits)}")
    print(f"Input rows parsed: {len(frame)}")
    print(f"Unique traits: {frame['trait'].nunique()}")
    print(f"Unique profiles: {frame['profile_id'].nunique()}")
    print(f"Ranked traits: {len(overall)}")
    if write_overall and overall_out is not None:
        print(f"Wrote overall ranking: {overall_out}")
    if write_pair and pair_out is not None:
        print(f"Wrote pair ranking: {pair_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
