#!/usr/bin/env python3
"""Export ranked enrichment trait summaries from a trait enrichment report."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

DEFAULT_INPUT_CANDIDATES = [
    Path("data/derived/trait_enrichment_report_p0.001.txt"),
    Path("data/derived/trait_enrichment_report.txt"),
]
DEFAULT_RESULTS_CANDIDATES = [
    Path("data/derived/trait_enrichment_results_p0.001.tsv"),
    Path("data/derived/trait_enrichment_results.tsv"),
]
DEFAULT_OVERALL_OUTPUT_STEM = "trait_enrichment_traits_overall_ranked"
DEFAULT_PAIR_OUTPUT_STEM = "trait_enrichment_trait_study_pairs_ranked"


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


def resolve_default_input() -> Path:
    for candidate in DEFAULT_INPUT_CANDIDATES:
        if candidate.exists():
            return candidate
    return DEFAULT_INPUT_CANDIDATES[-1]


def resolve_default_results(report_path: Path | None = None) -> Path:
    if report_path:
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


def parse_args() -> argparse.Namespace:
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
        help="Input trait enrichment results TSV (used for data aggregation). Defaults to companion of --input.",
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
    return parser.parse_args()


def _parse_summary_values(summary: str) -> tuple[float | None, float | None, float | None]:
    obs_match = re.search(r"obs=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", summary)
    exp_match = re.search(r"exp=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", summary)
    fold_match = re.search(r"fold=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", summary)
    obs = float(obs_match.group(1)) if obs_match else None
    exp = float(exp_match.group(1)) if exp_match else None
    fold = float(fold_match.group(1)) if fold_match else None
    return obs, exp, fold


def detect_metadata(path: Path) -> dict[str, any]:
    """Detect metadata from the report or results TSV."""
    meta = {
        "stat_type": "p_value",
        "threshold": 0.001,
        "testable_profiles": None,
        "max_traits": None,
        "fdr_scope": "global",
    }
    if not path.exists():
        return meta

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            # Text report format
            if line.startswith("p < "):
                meta["stat_type"] = "p_value"
                meta["threshold"] = float(line.split("<")[1].strip())
            elif line.startswith("FDR q < "):
                meta["stat_type"] = "fdr_q"
                meta["threshold"] = float(line.split("<")[1].strip())
            # TSV metadata header format
            elif line.startswith("# FDR: "):
                is_fdr = line.split(": ")[1].strip().lower() == "true"
                meta["stat_type"] = "fdr_q" if is_fdr else "p_value"
            elif line.startswith("# P-Threshold: "):
                meta["threshold"] = float(line.split(": ")[1].strip())
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
    path: Path, stat_type: str, threshold: float, min_fold: float = 1.5
) -> list[ParsedRow]:
    if not path.exists():
        raise FileNotFoundError(f"Results TSV not found: {path}")

    # Use comment='#' to skip the metadata header emitted by trait_enrichment.py
    df = pd.read_csv(path, sep="\t", comment="#")

    # col used for filtering
    col = "p_value" if stat_type == "p_value" else "fdr_q"

    # Apply filtering to match report logic
    df = df[(df[col] <= threshold) & (df["fold"] >= min_fold)]
    print(
        f"Filtered to {len(df)} significant hits (source={path.name}, threshold={threshold}, fold={min_fold})"
    )

    rows = []
    for _, r in df.iterrows():
        trait = str(r["trait"])
        # Exclude MediaDive growth-medium traits for consistency
        if "mediadive.medium:" in trait:
            continue

        study_id = str(r["study_id"])
        experiment_id = str(r["experiment_id"])
        inc_dec = str(r["inc_dec"])
        rows.append(
            ParsedRow(
                study_id=study_id,
                experiment_id=experiment_id,
                inc_dec=inc_dec,
                profile_id=f"{study_id}|{experiment_id}|{inc_dec}",
                trait=str(r["trait"]),
                p_value=float(r["p_value"]),
                fdr_q=float(r["fdr_q"]),
                fold=float(r["fold"]),
                obs=float(r["count"]),
                exp=float(r["expected"]),
                testable_profiles=int(r["testable_profiles"]),
            )
        )
    return rows


def parse_report(path: Path, stat_type: str = "p_value") -> list[ParsedRow]:
    rows: list[ParsedRow] = []
    current_header: tuple[str, str, str] | None = None

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line == "---" or line.startswith("p < ") or line.startswith("FDR q < "):
                continue

            if line.startswith("- "):
                if current_header is None:
                    raise ValueError(f"Trait row found before profile header: {line!r}")
                payload = line[2:]
                parts = payload.split(", ")
                if len(parts) < 3:
                    raise ValueError(f"Unparseable trait row: {line!r}")

                trait = parts[0].strip()
                val = float(parts[1])
                hit_count: int | None = None
                if re.fullmatch(r"\d+", parts[2]):
                    hit_count = int(parts[2])

                summary = parts[3] if len(parts) >= 4 else parts[2]
                obs, exp, fold = _parse_summary_values(summary)

                # Pylint: current_header is checked at top of block
                assert current_header is not None
                study_id, experiment_id, inc_dec = current_header
                profile_id = f"{study_id}|{experiment_id}|{inc_dec}"
                # Exclude MediaDive growth-medium traits for consistency
                if "mediadive.medium:" in trait:
                    continue

                # Set whichever statistic we found, and default the other to 1.0
                p_value = val if stat_type == "p_value" else 1.0
                fdr_q = val if stat_type == "fdr_q" else 1.0

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
                    )
                )
                continue

            header_parts = line.split("\t")
            if len(header_parts) != 3:
                continue
            current_header = (header_parts[0], header_parts[1], header_parts[2])

    if not rows:
        raise ValueError(f"No trait rows parsed from {path}")

    return rows


def build_overall_table(
    frame: pd.DataFrame,
    stat_type: str,
    total_profiles: int | None = None,
    max_traits: int | None = None,
    fdr_scope: str = "global",
) -> pd.DataFrame:
    grouped = frame.groupby("trait", dropna=False)

    # If we are parsing from a truncated source, identify potentially censored profiles.
    # A profile is censored if it hit the cap.
    if max_traits:
        counts_per_profile = frame.groupby("profile_id").size()
        censored_profiles = set(counts_per_profile[counts_per_profile >= max_traits].index)
        frame["is_censored_hit"] = frame["profile_id"].isin(censored_profiles)
    else:
        frame["is_censored_hit"] = False

    # Use the relevant stat for choosing the 'best' row metadata
    stat_col = "fdr_q" if stat_type == "fdr_q" else "p_value"
    best_stat_col = "best_qvalue" if stat_type == "fdr_q" else "best_pvalue"
    median_stat_col = "median_qvalue" if stat_type == "fdr_q" else "median_pvalue"

    # Multiplicity control for the per-trait "best" statistic:
    # Since we are taking the minimum over all profiles for a trait, we should
    # ideally apply a correction here. For now, we keep the naming explicit.
    
    best_idx = grouped[stat_col].idxmin()
    best_rows = frame.loc[best_idx, ["trait", "inc_dec", "study_id", "experiment_id", "profile_id"]]
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
    
    # Also include the other statistic if it's available (non-1.0)
    if stat_type == "fdr_q":
        agg_map["best_pvalue"] = ("p_value", "min")
    else:
        agg_map["best_qvalue"] = ("fdr_q", "min")

    summary = grouped.agg(**agg_map)

    inc_counts = (
        frame[frame["inc_dec"] == "inc"].groupby("trait")["trait"].size().rename("inc_hits")
    )
    dec_counts = (
        frame[frame["inc_dec"] == "dec"].groupby("trait")["trait"].size().rename("dec_hits")
    )

    overall = summary.join(inc_counts, how="left").join(dec_counts, how="left")
    overall[["inc_hits", "dec_hits"]] = overall[["inc_hits", "dec_hits"]].fillna(0).astype(int)
    overall = overall.reset_index().merge(best_rows, on="trait", how="left")
    overall.insert(1, "stat_type", stat_type)

    # Calculate net direction based on vote count
    def get_net_direction(row):
        if row["inc_hits"] > row["dec_hits"]:
            return "inc"
        if row["dec_hits"] > row["inc_hits"]:
            return "dec"
        return "mixed"

    overall["net_direction"] = overall.apply(get_net_direction, axis=1)

    if max_traits:
        censored_counts = (
            frame[frame["is_censored_hit"]].groupby("trait").size().rename("n_censored_profiles")
        )
        overall = overall.join(censored_counts, how="left")
        overall["n_censored_profiles"] = overall["n_censored_profiles"].fillna(0).astype(int)

    # #119: Calculate prevalence fraction.
    # Denominator is the number of profiles where the trait was testable
    # (i.e., appeared at least once).
    overall["prevalence"] = overall["hit_rows"] / overall["testable_profiles"]
    
    # Add human-readable prevalence summary
    overall["prevalence_summary"] = (
        overall["hit_rows"].astype(str) + " / " + overall["testable_profiles"].astype(str)
    )

    # #115: Minimum hit thresholding for overall ranking to reduce noise.
    # We require traits to be significant in at least 2 profiles to appear in rankings.
    overall = overall[overall["hit_rows"] >= 2]

    # #119: Rank by prevalence (significant / testable) descending, then by significance.
    # This provides a directly interpretable ranking that avoids extreme-value inflation.
    sort_stat = "best_qvalue" if stat_type == "fdr_q" else "best_pvalue"
    overall = overall.sort_values(
        ["prevalence", sort_stat, "max_fold"],
        ascending=[False, True, False],
    )
    overall.insert(0, "rank", range(1, len(overall) + 1))

    # Reorder columns to put consensus metrics upfront
    cols = list(overall.columns)
    # Move net_direction, inc_hits, dec_hits after trait
    front_cols = ["rank", "trait", "net_direction", "inc_hits", "dec_hits", "prevalence_summary"]
    remaining = [c for c in cols if c not in front_cols]
    overall = overall[front_cols + remaining]

    return overall


def build_pair_table(frame: pd.DataFrame, stat_type: str) -> pd.DataFrame:
    pairs = frame.copy()
    sort_col = "fdr_q" if stat_type == "fdr_q" else "p_value"
    pairs = pairs.sort_values([sort_col, "fold"], ascending=[True, False])
    
    # Rename stat column for self-description
    if stat_type == "fdr_q":
        pairs = pairs.rename(columns={"fdr_q": "q_value"})
    else:
        pairs = pairs.rename(columns={"p_value": "p_value"})

    pairs.insert(0, "rank", range(1, len(pairs) + 1))
    return pairs


def main() -> int:
    args = parse_args()

    meta = detect_metadata(args.input)
    stat_type = meta["stat_type"]
    stat_threshold = meta["threshold"]
    testable_profiles = meta["testable_profiles"]
    fdr_scope = meta["fdr_scope"]
    max_traits = meta.get("max_traits")
    is_source_truncated = False

    results_path = args.results or resolve_default_results(args.input)

    if results_path.exists():
        print(f"Loading data from results TSV: {results_path}")
        results_meta = detect_metadata(results_path)
        
        # Prefer TSV metadata if available
        stat_type = results_meta.get("stat_type", stat_type)
        stat_threshold = results_meta.get("threshold", stat_threshold)
        fdr_scope = results_meta.get("fdr_scope", fdr_scope)
        testable_profiles = results_meta.get("testable_profiles") or testable_profiles
        
        rows = parse_results_tsv(results_path, stat_type, stat_threshold)
        # Results TSV contains all hits regardless of the report truncation cap
        is_source_truncated = False
        max_traits = None
    else:
        print(
            f"Warning: Results TSV not found at {results_path}. Falling back to report parser ({args.input.name})."
        )
        rows = parse_report(args.input, stat_type=stat_type)
        
        # Check for truncation in the report metadata
        if max_traits is not None:
            is_source_truncated = True
            print(
                f"\nCRITICAL WARNING: The source report '{args.input.name}' is truncated to {max_traits} traits per profile."
            )
            print(
                "Downstream rankings (hit_rows, unique_profiles, mean_fold) will be systematically biased."
            )
            print(
                "Please regenerate the report with --max-traits-per-row 0 or ensure the results TSV is available.\n"
            )

    frame = pd.DataFrame([row.__dict__ for row in rows])

    overall = build_overall_table(
        frame,
        stat_type,
        total_profiles=testable_profiles,
        max_traits=max_traits,
        fdr_scope=fdr_scope,
    )
    pairs = build_pair_table(frame, stat_type)

    # Resolve dynamic filenames
    # Include 'profile' in suffix if using per-profile FDR
    method_suffix = "q" if stat_type == "fdr_q" else "p"
    if stat_type == "fdr_q" and fdr_scope == "profile":
        method_suffix = "q_profile"
        
    suffix = f"_{method_suffix}{stat_threshold}"
    
    overall_out = args.overall_output
    if overall_out is None:
        overall_out = Path("data/derived") / f"{DEFAULT_OVERALL_OUTPUT_STEM}{suffix}.tsv"
        
    pair_out = args.pair_output
    if pair_out is None:
        pair_out = Path("data/derived") / f"{DEFAULT_PAIR_OUTPUT_STEM}{suffix}.tsv"

    overall_out.parent.mkdir(parents=True, exist_ok=True)
    pair_out.parent.mkdir(parents=True, exist_ok=True)

    # Write overall table with metadata header
    with overall_out.open("w", encoding="utf-8", newline="") as f:
        f.write(f"# Stat-Type: {stat_type}\n")
        f.write(f"# FDR-Scope: {fdr_scope}\n")
        f.write(f"# Stat-Threshold: {stat_threshold}\n")
        if testable_profiles:
            f.write(f"# Total-Testable-Profiles: {testable_profiles}\n")
        if is_source_truncated:
            f.write(f"# Max-Traits-Per-Profile-Cap: {max_traits}\n")
            f.write("# WARNING: Totals are censored due to source truncation.\n")
        overall.to_csv(f, sep="\t", index=False)

    # Write pair table with metadata header
    with pair_out.open("w", encoding="utf-8", newline="") as f:
        f.write(f"# Stat-Type: {stat_type}\n")
        f.write(f"# FDR-Scope: {fdr_scope}\n")
        f.write(f"# Stat-Threshold: {stat_threshold}\n")
        if is_source_truncated:
            f.write(f"# Max-Traits-Per-Profile-Cap: {max_traits}\n")
            f.write("# WARNING: Data is censored due to source truncation.\n")
        pairs.to_csv(f, sep="\t", index=False)

    print(f"Input rows parsed: {len(frame)}")
    print(f"Unique traits: {frame['trait'].nunique()}")
    print(f"Unique profiles: {frame['profile_id'].nunique()}")
    print(f"Wrote overall ranking: {overall_out}")
    print(f"Wrote pair ranking: {pair_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
