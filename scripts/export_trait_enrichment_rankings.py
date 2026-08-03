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
    meta = {"stat_type": "p_value", "threshold": 0.001, "testable_profiles": None}
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
            elif line.startswith("# Testable-Profiles: "):
                meta["testable_profiles"] = int(line.split(": ")[1].strip())
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


def build_overall_table(frame: pd.DataFrame, stat_type: str, total_profiles: int | None = None) -> pd.DataFrame:
    grouped = frame.groupby("trait", dropna=False)

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
            "inc_dec": "best_direction",
            "study_id": "best_study_id",
            "experiment_id": "best_experiment_id",
            "profile_id": "best_profile_id",
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

    # #119: Calculate prevalence fraction.
    # Denominator is the number of profiles where the trait was testable
    # (i.e., appeared at least once).
    overall["prevalence"] = overall["hit_rows"] / overall["testable_profiles"]

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

    results_path = args.results or resolve_default_results(args.input)

    if results_path.exists():
        print(f"Loading data from results TSV: {results_path}")
        # If we just detected metadata from the report, check the results TSV too
        if testable_profiles is None:
            results_meta = detect_metadata(results_path)
            stat_type = results_meta["stat_type"]
            stat_threshold = results_meta["threshold"]
            testable_profiles = results_meta["testable_profiles"]
        
        rows = parse_results_tsv(results_path, stat_type, stat_threshold)
    else:
        print(
            f"Warning: Results TSV not found at {results_path}. Falling back to report parser ({args.input.name})."
        )
        rows = parse_report(args.input, stat_type=stat_type)
        
        # Check for truncation in the report metadata
        max_traits = meta.get("max_traits")
        if max_traits is not None:
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

    overall = build_overall_table(frame, stat_type, total_profiles=testable_profiles)
    pairs = build_pair_table(frame, stat_type)

    # Resolve dynamic filenames
    suffix = f"_{'q' if stat_type == 'fdr_q' else 'p'}{stat_threshold}"
    
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
        f.write(f"# Stat-Threshold: {stat_threshold}\n")
        if testable_profiles:
            f.write(f"# Total-Testable-Profiles: {testable_profiles}\n")
        overall.to_csv(f, sep="\t", index=False)

    # Write pair table with metadata header
    with pair_out.open("w", encoding="utf-8", newline="") as f:
        f.write(f"# Stat-Type: {stat_type}\n")
        f.write(f"# Stat-Threshold: {stat_threshold}\n")
        pairs.to_csv(f, sep="\t", index=False)

    print(f"Input rows parsed: {len(frame)}")
    print(f"Unique traits: {frame['trait'].nunique()}")
    print(f"Unique profiles: {frame['profile_id'].nunique()}")
    print(f"Wrote overall ranking: {overall_out}")
    print(f"Wrote pair ranking: {pair_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
