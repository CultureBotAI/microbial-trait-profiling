#!/usr/bin/env python3
"""Trait enrichment using a hypergeometric test.

Method summary:
1. Parse trait count rows grouped by study_id, experiment_id, and inc_dec
2. Compute per-row hypergeometric enrichment p-values for observed trait counts
3. Apply global Benjamini-Hochberg FDR correction
4. Keep only traits passing a significance threshold (default FDR q <= 0.05)
5. Emit a compact manuscript-style bullet list and optional structured TSV
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

from scipy.stats import hypergeom as _scipy_hypergeom

from bsdbng.enrichment_results import (
    META_RESULT_STATUS,
    META_SIGNIFICANT_HITS,
    META_TESTS_PERFORMED,
    classify_result_status,
)
from bsdbng.enrichment_stats import (
    benjamini_hochberg,
)
from bsdbng.enrichment_stats import (
    hypergeom_sf as _core_hypergeom_sf,
)
from bsdbng.enrichment_stats import (
    hypergeom_sf_many as _core_hypergeom_sf_many,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Trait enrichment report with hypergeometric p-values"
    )
    parser.add_argument("--input", type=Path, required=True, help="Path to trait TSV file")
    parser.add_argument("--output", type=Path, required=True, help="Path to write text report")
    parser.add_argument(
        "--tsv-output",
        type=Path,
        help="Path to write structured TSV results",
    )
    parser.add_argument(
        "--all-results-output",
        type=Path,
        help="Path to write complete unfiltered TSV results",
    )
    parser.add_argument(
        "--taxa-matrix",
        type=Path,
        default=None,
        help="Optional path to taxa_trait_matrix.tsv to define M and K.",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.001,
        help="Significance threshold (applied to FDR if used, else p-value)",
    )
    parser.add_argument(
        "--fdr",
        action="store_true",
        default=True,
        help="Apply Benjamini-Hochberg FDR correction (default: True)",
    )
    parser.add_argument(
        "--no-fdr",
        action="store_false",
        dest="fdr",
        help="Disable FDR correction",
    )
    parser.add_argument(
        "--fdr-scope",
        choices=["global", "profile"],
        default="global",
        help="Scope for FDR correction (default: global)",
    )
    parser.add_argument(
        "--min-fold",
        type=float,
        default=1.5,
        help="Minimum fold enrichment (count / expected)",
    )
    parser.add_argument(
        "--max-traits-per-row",
        type=int,
        default=0,
        help="Maximum number of bullet traits to show per row in text report (0 or -1 for all)",
    )
    return parser.parse_args(argv)


def _to_float(value: str) -> float:
    cleaned = value.strip()
    if cleaned == "":
        return 0.0
    return float(cleaned)


def _to_int(value: str) -> int:
    return int(_to_float(value))


def hypergeom_sf(x: int, population_size: int, population_successes: int, draws: int) -> float:
    """Return P[X >= x] using the required SciPy implementation."""
    return _core_hypergeom_sf(
        x, population_size, population_successes, draws, scipy_hypergeom=_scipy_hypergeom
    )


def hypergeom_sf_many(
    counts: list[int],
    population_size: int,
    population_successes: list[int],
    draws: int,
) -> list[float]:
    """Compute each trait's tail probability without silently replacing failed calculations."""
    return _core_hypergeom_sf_many(
        counts, population_size, population_successes, draws, scipy_hypergeom=_scipy_hypergeom
    )


def format_pvalue(value: float) -> str:
    if value < 0.0001:
        return f"{value:.6g}"
    return f"{value:.4f}"


def format_trait_name(name: str) -> str:
    words = name.split()
    if not words:
        return name
    return " ".join(word if word.isupper() else word.capitalize() for word in words)


def hypergeom_summary(count: int, row_total: int, trait_total: int, background_total: int) -> str:
    if background_total <= 0 or row_total <= 0:
        return "obs=0 exp=0 fold=NA"
    expected = row_total * (trait_total / background_total)
    fold = (count / expected) if expected > 0 else float("inf")
    fold_str = "inf" if math.isinf(fold) else f"{fold:.2f}"
    return f"obs={count} exp={expected:.2f} fold={fold_str}"


def open_input_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    with open_input_text(args.input) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        required = {"study_id", "experiment_id", "inc_dec"}
        missing = required - set(fieldnames)
        if missing:
            raise ValueError(f"Input TSV missing required columns: {sorted(missing)}")

        exclude_cols = {"study_id", "experiment_id", "inc_dec", "taxa_count"}
        trait_columns = [
            col for col in fieldnames if col not in exclude_cols and "mediadive.medium:" not in col
        ]
        if not trait_columns:
            raise ValueError("No trait columns found in input TSV")

        # 1. Read all rows to identify paired experiments
        all_raw_rows = list(reader)
        if not all_raw_rows:
            raise ValueError(
                f"Input TSV '{args.input}' has a header but no profile rows; "
                "regenerate it with `just study-trait-profiles`."
            )
        experiment_directions = defaultdict(set)
        for row in all_raw_rows:
            sid = row.get("study_id", "").strip()
            eid = row.get("experiment_id", "").strip()
            direction = row.get("inc_dec", "").strip().lower()
            if sid and eid and direction in {"inc", "dec"}:
                experiment_directions[(sid, eid)].add(direction)

        # Identify experiments with both 'inc' and 'dec'
        paired_experiments = {
            key for key, dirs in experiment_directions.items() if "inc" in dirs and "dec" in dirs
        }

        print(f"Total experiments found: {len(experiment_directions)}")
        print(f"Paired experiments identified: {len(paired_experiments)}")
        print(
            f"Dropping {len(experiment_directions) - len(paired_experiments)} unpaired experiments."
        )

        if not paired_experiments:
            # Valid but untestable input: still write outputs (status
            # zero_eligible_tests) so downstream steps never read stale files (#296).
            print("No paired experiments found; no tests are eligible.")

        rows = []
        # Background counts for hypergeometric test:
        # K(T) = number of unique taxa carrying trait T in reference background
        # M = total number of unique taxa in reference background
        background_counts = dict.fromkeys(trait_columns, 0)
        experiment_traits: dict[tuple[str, str], set[str]] = defaultdict(set)
        profile_taxa_counts: list[int] = []

        # 2. Process only paired rows
        for row_idx, row in enumerate(all_raw_rows, start=2):  # line 1 is header
            sid = row.get("study_id", "").strip()
            eid = row.get("experiment_id", "").strip()
            if (sid, eid) not in paired_experiments:
                continue

            row_traits_data = {}
            raw_taxa_count = row.get("taxa_count")
            direction = row.get("inc_dec", "").strip()

            for t in trait_columns:
                val = row.get(t)
                if val and val != "0":
                    try:
                        c = int(val)
                    except ValueError:
                        try:
                            c = int(float(val))
                        except ValueError:
                            raise ValueError(
                                f"File '{args.input}', line {row_idx}: invalid non-integer "
                                f"trait count '{val}' for trait '{t}' "
                                f"(study_id='{sid}', experiment_id='{eid}', inc_dec='{direction}')"
                            ) from None
                    if c < 0:
                        raise ValueError(
                            f"File '{args.input}', line {row_idx}: negative trait count '{c}' "
                            f"for trait '{t}' (study_id='{sid}', experiment_id='{eid}', "
                            f"inc_dec='{direction}')"
                        )
                    if c > 0:
                        row_traits_data[t] = c
                        experiment_traits[(sid, eid)].add(t)

            # Validate required draw count taxa_count (N) strictly without fallback
            if raw_taxa_count is None or not raw_taxa_count.strip():
                raise ValueError(
                    f"File '{args.input}', line {row_idx}: missing required 'taxa_count' column "
                    f"or value for profile (study_id='{sid}', experiment_id='{eid}', "
                    f"inc_dec='{direction}')"
                )

            try:
                n_taxa = int(raw_taxa_count.strip())
            except ValueError:
                raise ValueError(
                    f"File '{args.input}', line {row_idx}: malformed 'taxa_count' value "
                    f"'{raw_taxa_count}' (study_id='{sid}', experiment_id='{eid}', "
                    f"inc_dec='{direction}')"
                ) from None

            if n_taxa < 0:
                raise ValueError(
                    f"File '{args.input}', line {row_idx}: negative 'taxa_count' {n_taxa} "
                    f"(study_id='{sid}', experiment_id='{eid}', inc_dec='{direction}')"
                )

            max_k_row = max(row_traits_data.values(), default=0)
            if n_taxa < max_k_row:
                raise ValueError(
                    f"File '{args.input}', line {row_idx}: profile draw 'taxa_count'={n_taxa} "
                    f"is less than observed trait count {max_k_row} "
                    f"(study_id='{sid}', experiment_id='{eid}', inc_dec='{direction}')"
                )

            rows.append(
                {
                    "study_id": sid,
                    "experiment_id": eid,
                    "inc_dec": row.get("inc_dec", ""),
                    "traits": row_traits_data,
                    "n_taxa": n_taxa,
                }
            )
            profile_taxa_counts.append(n_taxa)

    # Determine reference population M and background successes K
    if args.taxa_matrix is not None:
        if not args.taxa_matrix.exists():
            raise FileNotFoundError(
                f"Specified taxa matrix file does not exist: {args.taxa_matrix}"
            )
        matrix_taxa_total = 0
        with open_input_text(args.taxa_matrix) as m_fh:
            m_reader = csv.DictReader(m_fh, delimiter="\t")
            for m_row in m_reader:
                matrix_taxa_total += 1
                for t in trait_columns:
                    val = m_row.get(t)
                    if val and val != "0":
                        try:
                            if int(val) > 0:
                                background_counts[t] += 1
                        except ValueError:
                            pass
        if matrix_taxa_total == 0:
            raise ValueError(f"Taxa matrix file {args.taxa_matrix} contains no taxon rows")
        background_taxa_universe = matrix_taxa_total
        bg_mode = f"TaxaMatrix ({args.taxa_matrix})"
    else:
        # Background empirical approximation:
        # Sum of profile taxa N across all testable profiles, and sum of trait presences k
        for r in rows:
            for t, c in r["traits"].items():
                background_counts[t] += c
        background_taxa_universe = sum(profile_taxa_counts)
        bg_mode = "EmpiricalProfileSum"

    # Validate profile draws and background counts
    max_k_trait = max(background_counts.values(), default=0)
    max_n_profile = max(profile_taxa_counts, default=0)
    if background_taxa_universe < max(max_k_trait, max_n_profile):
        raise ValueError(
            f"Background universe M ({background_taxa_universe}) is smaller than "
            f"maximum trait count K ({max_k_trait}) or maximum profile draw N ({max_n_profile})"
        )

    # Aligned testability denominators:
    # A profile is testable for trait t if it has n_taxa > 0, t is observed in the
    # experiment's paired arm (experiment_traits), and background_counts[t] > 0.
    testable_counts = dict.fromkeys(trait_columns, 0)
    for r in rows:
        if r["n_taxa"] <= 0:
            continue
        exp_key = (r["study_id"], r["experiment_id"])
        for t in experiment_traits[exp_key]:
            if background_counts[t] > 0:
                testable_counts[t] += 1

    # Second pass: Collect all P-values
    all_hits: list[dict] = []
    for row in rows:
        n_draws = row["n_taxa"]
        if n_draws == 0:
            continue

        traits_dict = row["traits"]
        if not traits_dict:
            continue

        row_traits = list(traits_dict.keys())
        row_counts = list(traits_dict.values())
        row_successes = [int(background_counts[t]) for t in row_traits]

        p_values = hypergeom_sf_many(row_counts, background_taxa_universe, row_successes, n_draws)

        # Per-profile FDR correction if requested
        if args.fdr and args.fdr_scope == "profile":
            q_values = benjamini_hochberg(p_values)
        else:
            q_values = p_values  # Will be corrected globally or left as p-values

        for trait, p, q, count, trait_total in zip(
            row_traits, p_values, q_values, row_counts, row_successes, strict=True
        ):
            expected = (
                n_draws * (trait_total / background_taxa_universe)
                if background_taxa_universe > 0
                else 0.0
            )
            fold = (count / expected) if expected > 0 else float("inf")
            all_hits.append(
                {
                    "study_id": row["study_id"],
                    "experiment_id": row["experiment_id"],
                    "inc_dec": row["inc_dec"],
                    "trait": trait,
                    "p_value": p,
                    "fdr_q": q,  # Preliminary q-value
                    "fold": fold,
                    "count": count,
                    "expected": expected,
                    "row_total": n_draws,
                    "trait_total": trait_total,
                    "testable_profiles": testable_counts[trait],
                }
            )

    # Apply global FDR correction if requested
    if args.fdr:
        if args.fdr_scope == "global":
            p_vals = [h["p_value"] for h in all_hits]
            q_vals = benjamini_hochberg(p_vals)
            for h, q in zip(all_hits, q_vals, strict=True):
                h["fdr_q"] = q
                h["is_significant"] = (q <= args.p_threshold) and (h["fold"] >= args.min_fold)
        else:
            # Already corrected per-profile
            for h in all_hits:
                h["is_significant"] = (h["fdr_q"] <= args.p_threshold) and (
                    h["fold"] >= args.min_fold
                )
    else:
        for h in all_hits:
            h["fdr_q"] = h["p_value"]
            h["is_significant"] = (h["p_value"] <= args.p_threshold) and (
                h["fold"] >= args.min_fold
            )

    # Filter significant hits
    significant_hits = [h for h in all_hits if h["is_significant"]]
    result_status = classify_result_status(len(all_hits), len(significant_hits))

    # Write TSV output if requested
    if args.tsv_output:
        args.tsv_output.parent.mkdir(parents=True, exist_ok=True)
        # Sort globally by significance (FDR q), then by study/experiment/direction/trait
        significant_hits.sort(
            key=lambda h: (h["fdr_q"], h["study_id"], h["experiment_id"], h["trait"])
        )
        with args.tsv_output.open("w", encoding="utf-8", newline="") as out_fh:
            # Metadata header
            out_fh.write(f"# FDR: {args.fdr}\n")
            out_fh.write(f"# FDR-Scope: {args.fdr_scope}\n")
            out_fh.write(f"# P-Threshold: {args.p_threshold}\n")
            out_fh.write(f"# Min-Fold: {args.min_fold}\n")
            out_fh.write(f"# Max-Traits-Per-Row-Report: {args.max_traits_per_row}\n")
            out_fh.write(f"# Background-Mode: {bg_mode}\n")
            testable_profile_count = len([n for n in profile_taxa_counts if n > 0])
            mean_taxa = (
                sum(profile_taxa_counts) / testable_profile_count
                if testable_profile_count > 0
                else 0
            )
            out_fh.write(f"# Total-Profiles-Input: {len(rows)}\n")
            out_fh.write(f"# Testable-Profiles: {testable_profile_count}\n")
            out_fh.write(f"# Background-Taxa-Universe-M: {background_taxa_universe}\n")
            out_fh.write(f"# Mean-Taxa-Per-Testable-Profile-N: {mean_taxa:.2f}\n")
            out_fh.write(f"# Mean-Assertions-Per-Testable-Profile: {mean_taxa:.2f}\n")
            out_fh.write(f"# {META_TESTS_PERFORMED}: {len(all_hits)}\n")
            out_fh.write(f"# {META_SIGNIFICANT_HITS}: {len(significant_hits)}\n")
            out_fh.write(f"# {META_RESULT_STATUS}: {result_status}\n")

            writer = csv.DictWriter(
                out_fh,
                fieldnames=[
                    "study_id",
                    "experiment_id",
                    "inc_dec",
                    "trait",
                    "p_value",
                    "fdr_q",
                    "fold",
                    "count",
                    "expected",
                    "testable_profiles",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            for h in significant_hits:
                writer.writerow(
                    {
                        "study_id": h["study_id"],
                        "experiment_id": h["experiment_id"],
                        "inc_dec": h["inc_dec"],
                        "trait": h["trait"],
                        "p_value": h["p_value"],
                        "fdr_q": h["fdr_q"],
                        "fold": h["fold"],
                        "count": h["count"],
                        "expected": h["expected"],
                        "testable_profiles": h["testable_profiles"],
                    }
                )

    # Write all-results output if requested
    if args.all_results_output:
        args.all_results_output.parent.mkdir(parents=True, exist_ok=True)
        # Sort globally by significance (FDR q), then by study/experiment/direction/trait
        all_hits.sort(key=lambda h: (h["fdr_q"], h["study_id"], h["experiment_id"], h["trait"]))
        with args.all_results_output.open("w", encoding="utf-8", newline="") as out_fh:
            # Metadata header
            out_fh.write(f"# FDR: {args.fdr}\n")
            out_fh.write(f"# FDR-Scope: {args.fdr_scope}\n")
            out_fh.write(f"# P-Threshold: {args.p_threshold}\n")
            out_fh.write(f"# Min-Fold: {args.min_fold}\n")
            out_fh.write(f"# Background-Mode: {bg_mode}\n")
            out_fh.write(f"# {META_TESTS_PERFORMED}: {len(all_hits)}\n")
            out_fh.write(f"# {META_SIGNIFICANT_HITS}: {len(significant_hits)}\n")
            out_fh.write(f"# {META_RESULT_STATUS}: {result_status}\n")

            writer = csv.DictWriter(
                out_fh,
                fieldnames=[
                    "study_id",
                    "experiment_id",
                    "inc_dec",
                    "trait",
                    "p_value",
                    "fdr_q",
                    "fold",
                    "count",
                    "expected",
                    "testable_profiles",
                    "background_m",
                    "background_k",
                    "is_significant",
                    "correction_family",
                    "eligibility_status",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            for h in all_hits:
                writer.writerow(
                    {
                        "study_id": h["study_id"],
                        "experiment_id": h["experiment_id"],
                        "inc_dec": h["inc_dec"],
                        "trait": h["trait"],
                        "p_value": h["p_value"],
                        "fdr_q": h["fdr_q"],
                        "fold": h["fold"],
                        "count": h["count"],
                        "expected": h["expected"],
                        "testable_profiles": h["testable_profiles"],
                        "background_m": background_taxa_universe,
                        "background_k": h["trait_total"],
                        "is_significant": h["is_significant"],
                        "correction_family": args.fdr_scope,
                        "eligibility_status": (
                            "testable" if h["testable_profiles"] > 0 else "ineligible"
                        ),
                    }
                )

    # Generate text report
    report_lines: list[str] = []
    if args.fdr:
        report_lines.append(f"FDR q < {args.p_threshold:g}")
    else:
        report_lines.append(f"p < {args.p_threshold:g}")
    if args.min_fold > 0:
        report_lines.append(f"Min Fold >= {args.min_fold:g}")
    if args.max_traits_per_row > 0:
        report_lines.append(f"Max traits per profile: {args.max_traits_per_row}")
    else:
        report_lines.append("Max traits per profile: unlimited")
    report_lines.append(f"Tests performed: {len(all_hits)}")
    report_lines.append(f"Significant hits: {len(significant_hits)}")
    report_lines.append(f"Result status: {result_status}")
    report_lines.append("")

    # Group significant hits for report
    report_groups: dict[tuple, list[dict]] = {}
    for h in significant_hits:
        key = (h["study_id"], h["experiment_id"], h["inc_dec"])
        if key not in report_groups:
            report_groups[key] = []
        report_groups[key].append(h)

    rows_with_hits = 0
    for key, group in report_groups.items():
        rows_with_hits += 1
        # Sort by significance then fold
        group.sort(key=lambda h: (h["fdr_q"], -h["fold"], h["trait"]))

        report_lines.append(f"{key[0]}\t{key[1]}\t{key[2]}")
        report_lines.append("")

        limit = args.max_traits_per_row
        subset = group if limit <= 0 else group[:limit]

        for h in subset:
            report_lines.append(
                "- "
                f"{format_trait_name(h['trait'])}, "
                f"{format_pvalue(h['fdr_q'])}, "
                f"fold={h['fold']:.2f}"
            )
        report_lines.append("")
        report_lines.append("---")
        report_lines.append("")

    if rows_with_hits == 0:
        report_lines.append("No traits passed the threshold.")
        report_lines.append("")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(report_lines).rstrip() + "\n", encoding="utf-8")

    print(f"Input rows: {len(rows)}")
    print(f"Trait columns: {len(trait_columns)}")
    testable_profile_count = len([n for n in profile_taxa_counts if n > 0])
    print(f"Total testable profiles: {testable_profile_count}")
    print(f"Background taxa universe M: {background_taxa_universe}")
    print(f"Total tests performed: {len(all_hits)}")
    print(f"Significant hits: {len(significant_hits)}")
    print(f"Result status: {result_status}")
    print(f"Rows with hits: {rows_with_hits}")
    print(f"Wrote report: {args.output}")
    if args.tsv_output:
        print(f"Wrote results TSV: {args.tsv_output}")
    if args.all_results_output:
        print(f"Wrote all-results TSV: {args.all_results_output}")


if __name__ == "__main__":
    main()
