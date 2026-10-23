#!/usr/bin/env python3
"""Trait enrichment or depletion using a hypergeometric test.

Method summary:
1. Parse trait count rows grouped by study_id, experiment_id, and inc_dec
2. Pre-specify testable family: all candidate traits with background K > 0
   tested against all profiles with taxa_count N > 0 (Testable-Set: K>0 & N>0)
3. Compute hypergeometric enrichment (upper tail) or depletion (lower tail) p-values
4. Apply global Benjamini-Hochberg FDR correction
5. Keep only traits passing significance thresholds (default FDR q <= 0.001)
6. Emit a compact manuscript-style bullet list and optional structured TSVs
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, TextIO

from scipy.stats import hypergeom as _scipy_hypergeom

from mtp.enrichment_results import (
    META_RESULT_STATUS,
    META_SIGNIFICANT_HITS,
    META_TESTS_PERFORMED,
    classify_result_status,
)
from mtp.enrichment_stats import (
    benjamini_hochberg,
)
from mtp.enrichment_stats import (
    hypergeom_cdf as _core_hypergeom_cdf,
)
from mtp.enrichment_stats import (
    hypergeom_cdf_many as _core_hypergeom_cdf_many,
)
from mtp.enrichment_stats import (
    hypergeom_sf as _core_hypergeom_sf,
)
from mtp.enrichment_stats import (
    hypergeom_sf_many as _core_hypergeom_sf_many,
)
from mtp.paths import file_sha256
from mtp.strain_rollup_record import shared_strain_rollup
from mtp.trait_guilds import MEDIADIVE_MARKER


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Trait enrichment or depletion report with hypergeometric p-values"
    )
    parser.add_argument("--core-run", type=Path, help="Verified core workflow run manifest")
    parser.add_argument("--input", type=Path, help="Path to trait TSV file")
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
        help=(
            "Path to taxa_trait_matrix.tsv defining reference background M "
            "(BugSigDB-reported taxa rows; mixed ranks, inherited positivity) and K."
        ),
    )
    parser.add_argument(
        "--tail",
        choices=["upper", "lower"],
        default="upper",
        help=(
            "Direction of hypergeometric tail: upper for enrichment, "
            "lower for depletion (default: upper)"
        ),
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.001,
        help="Significance threshold (applied to FDR if used, else p-value; default: 0.001)",
    )
    parser.add_argument(
        "--fdr",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Apply Benjamini-Hochberg FDR correction (default: True)",
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
        help="Minimum fold enrichment (count / expected) for upper tail (default: 1.5)",
    )
    parser.add_argument(
        "--max-fold",
        type=float,
        default=0.67,
        help="Maximum fold depletion (count / expected) for lower tail (default: 0.67)",
    )
    parser.add_argument(
        "--max-traits-per-row",
        type=int,
        default=0,
        help="Maximum number of bullet traits to show per row in text report (0 or -1 for all)",
    )
    args = parser.parse_args(argv)
    if args.core_run is None and (args.input is None or args.taxa_matrix is None):
        parser.error("Supply --core-run or both --input and --taxa-matrix")
    return args


def _parse_count(value: Any, trait: str = "", row_idx: int = 0, source: str = "") -> int:
    """Parse count strictly.

    Accepts integer strings and integer-valued floats (e.g. '1.0').
    Rejects malformed strings, non-integer floats (e.g. '1.5'), and negative numbers.
    """
    if value is None or value == "":
        return 0
    if isinstance(value, int):
        if value < 0:
            raise ValueError(f"Negative count {value} for trait '{trait}' in {source}:{row_idx}")
        return value
    val_str = str(value).strip()
    if not val_str or val_str == "0":
        return 0
    try:
        val_float = float(val_str)
    except ValueError as err:
        raise ValueError(
            f"Invalid count value {value!r} for trait '{trait}' in {source}:{row_idx}"
        ) from err
    if not val_float.is_integer():
        raise ValueError(
            f"Count value must be an integer, got {value!r} for trait '{trait}' "
            f"in {source}:{row_idx}"
        )
    count = int(val_float)
    if count < 0:
        raise ValueError(f"Negative count {count} for trait '{trait}' in {source}:{row_idx}")
    return count


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


def hypergeom_cdf(x: int, population_size: int, population_successes: int, draws: int) -> float:
    """Return P[X <= x] using the required SciPy implementation."""
    return _core_hypergeom_cdf(
        x, population_size, population_successes, draws, scipy_hypergeom=_scipy_hypergeom
    )


def hypergeom_cdf_many(
    counts: list[int],
    population_size: int,
    population_successes: list[int],
    draws: int,
) -> list[float]:
    """Compute each trait's lower tail probability without silently replacing failures."""
    return _core_hypergeom_cdf_many(
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
    fold = count / expected if expected > 0 else (float("inf") if count > 0 else 0.0)
    fold_str = "inf" if math.isinf(fold) else f"{fold:.2f}"
    return f"obs={count} exp={expected:.2f} fold={fold_str}"


def open_input_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def main(argv: list[str] | None = None) -> None:
    from mtp.workflow.consumers import (
        open_core_bundle,
        preflight_consumer_outputs,
        write_consumer_record,
    )

    args = parse_args(argv)
    bundle = None
    if args.core_run is not None:
        bundle = open_core_bundle(args.core_run, "ora_profiles")
        for argument, product in (("input", "ora_profiles"), ("taxa_matrix", "taxa_matrix")):
            path = bundle.path(product)
            supplied = getattr(args, argument)
            if supplied is not None and supplied.resolve() != path.resolve():
                raise ValueError(f"--core-run conflicts with --{argument.replace('_', '-')}")
            setattr(args, argument, path)
        preflight_consumer_outputs(
            bundle,
            args.output.with_suffix(".core_consumption.json"),
            [
                path
                for path in (args.output, args.tsv_output, args.all_results_output)
                if path is not None
            ],
        )

    # A profile table and matrix recorded under different scopes are refused (#847)
    shared_strain_rollup(args.input, args.taxa_matrix)
    # Bind exported tests to actual input bytes, including non-core invocations (#705).
    profiles_sha256 = file_sha256(args.input)
    with open_input_text(args.input) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        required = {"study_id", "experiment_id", "inc_dec"}
        missing = required - set(fieldnames)
        if missing:
            raise ValueError(f"Input TSV missing required columns: {sorted(missing)}")

        exclude_cols = {"study_id", "experiment_id", "inc_dec", "taxa_count"}
        trait_columns = [
            col
            for col in fieldnames
            if col not in exclude_cols and MEDIADIVE_MARKER not in col.lower()
        ]
        if not trait_columns:
            raise ValueError("No trait columns found in input TSV")

        # 1. Read all rows to identify paired experiments
        all_raw_rows = list(reader)
        if not all_raw_rows:
            raise ValueError(
                f"Input TSV '{args.input}' has a header but no profile rows; "
                "regenerate it with `just study-trait-profiles-ora`."
            )
        experiment_directions: dict[tuple[str, str], set[str]] = defaultdict(set)
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

        rows: list[dict[str, Any]] = []
        profile_taxa_counts: list[int] = []

        # 2. Process only paired rows
        for row_idx, row in enumerate(all_raw_rows, start=2):  # line 1 is header
            sid = row.get("study_id", "").strip()
            eid = row.get("experiment_id", "").strip()
            if (sid, eid) not in paired_experiments:
                continue

            row_traits_data: dict[str, int] = {}
            raw_taxa_count = row.get("taxa_count")
            direction = row.get("inc_dec", "").strip()

            for t in trait_columns:
                val = row.get(t)
                c = _parse_count(val, trait=t, row_idx=row_idx, source=str(args.input))
                if c > 0:
                    row_traits_data[t] = c

            # Validate required draw count taxa_count (N) strictly without fallback
            if raw_taxa_count is None or not raw_taxa_count.strip():
                raise ValueError(
                    f"File '{args.input}', line {row_idx}: missing required 'taxa_count' column "
                    f"or value for profile (study_id='{sid}', experiment_id='{eid}', "
                    f"inc_dec='{direction}')"
                )

            try:
                val_float = float(raw_taxa_count.strip())
            except ValueError:
                raise ValueError(
                    f"File '{args.input}', line {row_idx}: malformed 'taxa_count' value "
                    f"'{raw_taxa_count}' (study_id='{sid}', experiment_id='{eid}', "
                    f"inc_dec='{direction}')"
                ) from None

            if not val_float.is_integer():
                raise ValueError(
                    f"File '{args.input}', line {row_idx}: non-integer 'taxa_count' value "
                    f"'{raw_taxa_count}' (study_id='{sid}', experiment_id='{eid}', "
                    f"inc_dec='{direction}')"
                )
            n_taxa = int(val_float)

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
                    "inc_dec": direction,
                    "traits": row_traits_data,
                    "n_taxa": n_taxa,
                }
            )
            profile_taxa_counts.append(n_taxa)

    # Determine reference population M and background successes K from taxa_matrix
    if not args.taxa_matrix.exists():
        raise FileNotFoundError(f"Specified taxa matrix file does not exist: {args.taxa_matrix}")
    taxa_matrix_sha256 = file_sha256(args.taxa_matrix)
    background_counts: dict[str, int] = dict.fromkeys(trait_columns, 0)
    matrix_taxa_total = 0
    with open_input_text(args.taxa_matrix) as m_fh:
        m_reader = csv.DictReader(m_fh, delimiter="\t")
        m_fieldnames = set(m_reader.fieldnames or [])
        for row_idx, m_row in enumerate(m_reader, start=2):
            matrix_taxa_total += 1
            for t in trait_columns:
                if t not in m_fieldnames:
                    continue
                val = m_row.get(t)
                c = _parse_count(val, trait=t, row_idx=row_idx, source=str(args.taxa_matrix))
                if c > 0:
                    background_counts[t] += 1

    if matrix_taxa_total == 0:
        raise ValueError(f"Taxa matrix file {args.taxa_matrix} contains no taxon rows")

    background_taxa_universe = matrix_taxa_total
    bg_mode = f"TaxaMatrix ({args.taxa_matrix.name})"

    # Validate profile draws and background counts against reference universe M
    max_k_trait = max(background_counts.values(), default=0)
    max_n_profile = max(profile_taxa_counts, default=0)
    if background_taxa_universe < max(max_k_trait, max_n_profile):
        raise ValueError(
            f"Background universe M ({background_taxa_universe}) is smaller than "
            f"maximum trait count K ({max_k_trait}) or maximum profile draw N ({max_n_profile})"
        )

    # Issue #397: Validate profile counts k against background universe K
    # Reject k > K and k > 0 when K == 0 as hard errors.
    for r in rows:
        for t, k in r["traits"].items():
            K = background_counts.get(t, 0)
            if k > 0 and K == 0:
                raise ValueError(
                    f"Profile ({r['study_id']}, {r['experiment_id']}, {r['inc_dec']}) "
                    f"has count k={k} > 0 for trait '{t}', but reference background K is 0."
                )
            if k > K:
                raise ValueError(
                    f"Profile ({r['study_id']}, {r['experiment_id']}, {r['inc_dec']}) "
                    f"has count k={k} > reference background K={K} for trait '{t}'."
                )

    # Issue #398: Define one pre-specified family before testing:
    # All candidate traits with K > 0 and all paired profiles with N > 0.
    candidate_traits = [t for t in trait_columns if background_counts[t] > 0]
    testable_profiles_list = [r for r in rows if r["n_taxa"] > 0]
    testable_profile_count = len(testable_profiles_list)

    # Hypergeometric testing
    all_hits: list[dict[str, Any]] = []
    if candidate_traits and testable_profiles_list:
        trait_totals = [background_counts[t] for t in candidate_traits]
        for row in testable_profiles_list:
            n_draws = row["n_taxa"]
            row_counts = [row["traits"].get(t, 0) for t in candidate_traits]

            if args.tail == "upper":
                p_values = hypergeom_sf_many(
                    row_counts, background_taxa_universe, trait_totals, n_draws
                )
            else:
                p_values = hypergeom_cdf_many(
                    row_counts, background_taxa_universe, trait_totals, n_draws
                )

            # Per-profile FDR correction if requested
            if args.fdr and args.fdr_scope == "profile":
                q_values = benjamini_hochberg(p_values)
            else:
                q_values = p_values

            for trait, p, q, count, trait_total in zip(
                candidate_traits, p_values, q_values, row_counts, trait_totals, strict=True
            ):
                expected = (
                    n_draws * (trait_total / background_taxa_universe)
                    if background_taxa_universe > 0
                    else 0.0
                )
                fold = count / expected if expected > 0 else (float("inf") if count > 0 else 0.0)

                all_hits.append(
                    {
                        "study_id": row["study_id"],
                        "experiment_id": row["experiment_id"],
                        "inc_dec": row["inc_dec"],
                        "trait": trait,
                        "p_value": p,
                        "fdr_q": q,
                        "fold": fold,
                        "count": count,
                        "expected": expected,
                        "row_total": n_draws,
                        "trait_total": trait_total,
                        "testable_profiles": testable_profile_count,
                    }
                )

    def is_hit(stat: float, fold_val: float) -> bool:
        if args.tail == "upper":
            return (stat <= args.p_threshold) and (fold_val >= args.min_fold)
        return (stat <= args.p_threshold) and (fold_val <= args.max_fold)

    # Apply global FDR correction if requested
    if args.fdr:
        if args.fdr_scope == "global":
            p_vals = [h["p_value"] for h in all_hits]
            q_vals = benjamini_hochberg(p_vals)
            for h, q in zip(all_hits, q_vals, strict=True):
                h["fdr_q"] = q
                h["is_significant"] = is_hit(q, h["fold"])
        else:
            for h in all_hits:
                h["is_significant"] = is_hit(h["fdr_q"], h["fold"])
    else:
        for h in all_hits:
            h["fdr_q"] = h["p_value"]
            h["is_significant"] = is_hit(h["p_value"], h["fold"])

    # Filter significant hits
    significant_hits = [h for h in all_hits if h["is_significant"]]
    result_status = classify_result_status(len(all_hits), len(significant_hits))

    mean_taxa = (
        sum(profile_taxa_counts) / testable_profile_count if testable_profile_count > 0 else 0.0
    )

    # Write TSV output if requested
    if args.tsv_output:
        args.tsv_output.parent.mkdir(parents=True, exist_ok=True)
        # Sort globally by significance (FDR q), then by study/experiment/direction/trait
        significant_hits.sort(
            key=lambda h: (h["fdr_q"], h["study_id"], h["experiment_id"], h["trait"])
        )
        with args.tsv_output.open("w", encoding="utf-8", newline="") as out_fh:
            # Metadata header
            out_fh.write(f"# Tail: {args.tail}\n")
            out_fh.write("# Testable-Set: K>0 & N>0\n")
            out_fh.write(f"# FDR: {args.fdr}\n")
            out_fh.write(f"# FDR-Scope: {args.fdr_scope}\n")
            out_fh.write(f"# P-Threshold: {args.p_threshold}\n")
            if args.tail == "upper":
                out_fh.write(f"# Min-Fold: {args.min_fold}\n")
            else:
                out_fh.write(f"# Max-Fold: {args.max_fold}\n")
            out_fh.write(f"# Max-Traits-Per-Row-Report: {args.max_traits_per_row}\n")
            out_fh.write(f"# Background-Mode: {bg_mode}\n")
            out_fh.write(f"# Total-Profiles-Input: {len(rows)}\n")
            out_fh.write(f"# Testable-Profiles: {testable_profile_count}\n")
            out_fh.write(f"# Background-Taxa-Universe-M: {background_taxa_universe}\n")
            out_fh.write(
                "# Background-Universe-Definition: BugSigDB-reported taxa rows of "
                f"{args.taxa_matrix.name} (mixed ranks, inherited positivity)\n"
            )
            out_fh.write(f"# Mean-Taxa-Per-Testable-Profile-N: {mean_taxa:.2f}\n")
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
            out_fh.write(f"# Profiles-SHA256: {profiles_sha256}\n")
            out_fh.write(f"# Taxa-Matrix-SHA256: {taxa_matrix_sha256}\n")
            out_fh.write(f"# Tail: {args.tail}\n")
            out_fh.write("# Testable-Set: K>0 & N>0\n")
            out_fh.write(f"# FDR: {args.fdr}\n")
            out_fh.write(f"# FDR-Scope: {args.fdr_scope}\n")
            out_fh.write(f"# P-Threshold: {args.p_threshold}\n")
            if args.tail == "upper":
                out_fh.write(f"# Min-Fold: {args.min_fold}\n")
            else:
                out_fh.write(f"# Max-Fold: {args.max_fold}\n")
            out_fh.write(f"# Background-Mode: {bg_mode}\n")
            out_fh.write(f"# Total-Profiles-Input: {len(rows)}\n")
            out_fh.write(f"# Testable-Profiles: {testable_profile_count}\n")
            out_fh.write(f"# Background-Taxa-Universe-M: {background_taxa_universe}\n")
            out_fh.write(
                "# Background-Universe-Definition: BugSigDB-reported taxa rows of "
                f"{args.taxa_matrix.name} (mixed ranks, inherited positivity)\n"
            )
            out_fh.write(f"# Mean-Taxa-Per-Testable-Profile-N: {mean_taxa:.2f}\n")
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
                        "eligibility_status": "testable",
                    }
                )

    # Generate text report
    report_lines: list[str] = []
    tail_suffix = " (Depletion)" if args.tail == "lower" else ""
    if args.fdr:
        report_lines.append(f"FDR q <= {args.p_threshold:g}{tail_suffix}")
    else:
        report_lines.append(f"p <= {args.p_threshold:g}{tail_suffix}")
    if args.tail == "upper":
        if args.min_fold > 0:
            report_lines.append(f"Min Fold >= {args.min_fold:g}")
    else:
        if args.max_fold < float("inf"):
            report_lines.append(f"Max Fold <= {args.max_fold:g}")
    if args.max_traits_per_row > 0:
        report_lines.append(f"Max traits per profile: {args.max_traits_per_row}")
    else:
        report_lines.append("Max traits per profile: unlimited")
    report_lines.append(f"Tests performed: {len(all_hits)}")
    report_lines.append(f"Significant hits: {len(significant_hits)}")
    report_lines.append(f"Result status: {result_status}")
    report_lines.append("")

    # Group significant hits for report
    report_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for h in significant_hits:
        key = (h["study_id"], h["experiment_id"], h["inc_dec"])
        if key not in report_groups:
            report_groups[key] = []
        report_groups[key].append(h)

    rows_with_hits = 0
    for key, group in report_groups.items():
        rows_with_hits += 1
        # Sort by significance then fold
        if args.tail == "upper":
            group.sort(key=lambda h: (h["fdr_q"], -h["fold"], h["trait"]))
        else:
            group.sort(key=lambda h: (h["fdr_q"], h["fold"], h["trait"]))

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
        if args.tail == "lower":
            report_lines.append("No traits passed the depletion threshold.")
        else:
            report_lines.append("No traits passed the threshold.")
        report_lines.append("")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(report_lines).rstrip() + "\n", encoding="utf-8")

    print(f"Input rows: {len(rows)}")
    print(f"Trait columns: {len(trait_columns)}")
    print(f"Candidate traits (K > 0): {len(candidate_traits)}")
    print(f"Total testable profiles (N > 0): {testable_profile_count}")
    if bundle is not None:
        write_consumer_record(
            bundle,
            args.output.with_suffix(".core_consumption.json"),
            consumer="scripts/trait_enrichment.py",
            outputs=[
                path
                for path in (args.output, args.tsv_output, args.all_results_output)
                if path is not None
            ],
        )
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
