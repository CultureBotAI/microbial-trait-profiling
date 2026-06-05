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
DEFAULT_OVERALL_OUTPUT = Path("data/derived/trait_enrichment_traits_overall_ranked.tsv")
DEFAULT_PAIR_OUTPUT = Path("data/derived/trait_enrichment_trait_study_pairs_ranked.tsv")


@dataclass(frozen=True)
class ParsedRow:
    study_id: str
    experiment_id: str
    inc_dec: str
    profile_id: str
    trait: str
    pvalue: float
    fold: float | None
    obs: float | None
    exp: float | None
    hit_count: int | None


def resolve_default_input() -> Path:
    for candidate in DEFAULT_INPUT_CANDIDATES:
        if candidate.exists():
            return candidate
    return DEFAULT_INPUT_CANDIDATES[-1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export ranked trait enrichment TSV outputs from a text report."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=resolve_default_input(),
        help="Input trait enrichment report path.",
    )
    parser.add_argument(
        "--overall-output",
        type=Path,
        default=DEFAULT_OVERALL_OUTPUT,
        help="Output TSV for overall trait ranking.",
    )
    parser.add_argument(
        "--pair-output",
        type=Path,
        default=DEFAULT_PAIR_OUTPUT,
        help="Output TSV for trait-study pair ranking.",
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


def parse_report(path: Path) -> list[ParsedRow]:
    rows: list[ParsedRow] = []
    current_header: tuple[str, str, str] | None = None

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line == "---" or line.startswith("p < "):
                continue

            if line.startswith("- "):
                if current_header is None:
                    raise ValueError(f"Trait row found before profile header: {line!r}")
                payload = line[2:]
                parts = payload.split(", ")
                if len(parts) < 3:
                    raise ValueError(f"Unparseable trait row: {line!r}")

                trait = parts[0].strip()
                pvalue = float(parts[1])
                hit_count: int | None = None
                if re.fullmatch(r"\d+", parts[2]):
                    hit_count = int(parts[2])

                summary = parts[3] if len(parts) >= 4 else parts[2]
                obs, exp, fold = _parse_summary_values(summary)

                study_id, experiment_id, inc_dec = current_header
                profile_id = f"{study_id}|{experiment_id}|{inc_dec}"
                rows.append(
                    ParsedRow(
                        study_id=study_id,
                        experiment_id=experiment_id,
                        inc_dec=inc_dec,
                        profile_id=profile_id,
                        trait=trait,
                        pvalue=pvalue,
                        fold=fold,
                        obs=obs,
                        exp=exp,
                        hit_count=hit_count,
                    )
                )
                continue

            header_parts = line.split("\t")
            if len(header_parts) != 3:
                raise ValueError(f"Unparseable profile header: {line!r}")
            current_header = (header_parts[0], header_parts[1], header_parts[2])

    if not rows:
        raise ValueError(f"No trait rows parsed from {path}")

    return rows


def build_overall_table(frame: pd.DataFrame) -> pd.DataFrame:
    grouped = frame.groupby("trait", dropna=False)

    best_idx = grouped["pvalue"].idxmin()
    best_rows = frame.loc[best_idx, ["trait", "inc_dec", "study_id", "experiment_id", "profile_id"]]
    best_rows = best_rows.rename(
        columns={
            "inc_dec": "best_direction",
            "study_id": "best_study_id",
            "experiment_id": "best_experiment_id",
            "profile_id": "best_profile_id",
        }
    )

    summary = grouped.agg(
        best_pvalue=("pvalue", "min"),
        median_pvalue=("pvalue", "median"),
        hit_rows=("trait", "size"),
        unique_profiles=("profile_id", "nunique"),
        mean_fold=("fold", "mean"),
        max_fold=("fold", "max"),
    )

    inc_counts = (
        frame[frame["inc_dec"] == "inc"].groupby("trait")["trait"].size().rename("inc_hits")
    )
    dec_counts = (
        frame[frame["inc_dec"] == "dec"].groupby("trait")["trait"].size().rename("dec_hits")
    )

    overall = summary.join(inc_counts, how="left").join(dec_counts, how="left")
    overall[["inc_hits", "dec_hits"]] = overall[["inc_hits", "dec_hits"]].fillna(0).astype(int)
    overall = overall.reset_index().merge(best_rows, on="trait", how="left")
    overall = overall.sort_values(
        ["best_pvalue", "hit_rows", "max_fold"],
        ascending=[True, False, False],
    )
    overall.insert(0, "rank", range(1, len(overall) + 1))
    return overall


def build_pair_table(frame: pd.DataFrame) -> pd.DataFrame:
    pairs = frame.copy()
    pairs = pairs.sort_values(["pvalue", "fold"], ascending=[True, False])
    pairs.insert(0, "rank", range(1, len(pairs) + 1))
    return pairs


def main() -> int:
    args = parse_args()
    rows = parse_report(args.input)
    frame = pd.DataFrame([row.__dict__ for row in rows])

    overall = build_overall_table(frame)
    pairs = build_pair_table(frame)

    args.overall_output.parent.mkdir(parents=True, exist_ok=True)
    args.pair_output.parent.mkdir(parents=True, exist_ok=True)
    overall.to_csv(args.overall_output, sep="\t", index=False)
    pairs.to_csv(args.pair_output, sep="\t", index=False)

    print(f"Input rows parsed: {len(frame)}")
    print(f"Unique traits: {frame['trait'].nunique()}")
    print(f"Unique profiles: {frame['profile_id'].nunique()}")
    print(f"Wrote overall ranking: {args.overall_output}")
    print(f"Wrote pair ranking: {args.pair_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
