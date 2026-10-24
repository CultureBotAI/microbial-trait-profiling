"""Measure #976 ORA validation on a synthetic complete tail (2026-10-23).

Inputs: deterministic --profiles x --traits synthetic probabilities, never real
analysis inputs. Outputs: tail.tsv plus generate/parse/validate JSON evidence.
Run each phase in a fresh process; large phases require the shared heavy slot.
The parse baseline measures streamed probability/coordinate checks, not an old
whole-adapter replay. Peak RSS is process-wide, not incremental allocation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

from mtp.enrichment_stats import benjamini_hochberg
from mtp.multiplicity.families import probability, rows
from mtp.multiplicity.ora import native_bh_rows
from mtp.paths import file_sha256


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("generate", "parse", "validate"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--profiles", type=int, default=1870)
    parser.add_argument("--traits", type=int, default=317)
    parser.add_argument("--seed", type=int, default=976)
    args = parser.parse_args()
    if args.profiles < 1 or args.traits < 1:
        raise ValueError("Synthetic family dimensions must be positive")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    path = output / "tail.tsv"
    keys = [
        (f"study-{index // 2}", "experiment", "inc" if index % 2 else "dec")
        for index in range(args.profiles)
    ]
    features = [
        f"biolink:has_phenotype|METPO:synthetic-{index:07d}" for index in range(args.traits)
    ]
    count = args.profiles * args.traits
    if args.phase == "generate":
        if path.exists():
            raise FileExistsError(path)
        p_values = np.random.default_rng(args.seed).permutation(
            (np.arange(count, dtype=float) + 1) / (count + 1)
        )
        q_values = benjamini_hochberg(p_values)
        with path.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
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
                ]
            )
            for index, (p, q) in enumerate(zip(p_values, q_values, strict=True)):
                writer.writerow(
                    [
                        *keys[index // args.traits],
                        features[index % args.traits],
                        p,
                        q,
                        1.0,
                        12,
                        12.0,
                        args.profiles,
                        3221,
                        101,
                        False,
                        "global",
                        "testable",
                    ]
                )
        record = {
            "synthetic": True,
            "profiles": args.profiles,
            "traits": args.traits,
            "native_tests": count,
            "columns": 15,
            "scope": "synthetic native-export-width rows; p values are not hypergeometric tests",
            "p_order": "seeded permutation of equally spaced positive probabilities",
            "seed": args.seed,
            "sha256": file_sha256(path),
            "bytes": path.stat().st_size,
        }
        (output / "generate.json").write_text(json.dumps(record, indent=2) + "\n")
        return 0
    generated = json.loads((output / "generate.json").read_text())
    expected_hash = generated["sha256"]
    if (generated["profiles"], generated["traits"], generated["seed"]) != (
        args.profiles,
        args.traits,
        args.seed,
    ):
        raise ValueError("Synthetic family dimensions or seed differ from generation")
    if file_sha256(path) != expected_hash:
        raise ValueError("Synthetic input changed since generation")
    profile_keys, feature_set = set(keys), set(features)
    consumed = hashlib.sha256()
    before = resource.getrusage(resource.RUSAGE_SELF)
    start = time.perf_counter()
    iterator = (
        native_bh_rows(
            path, expected_sha256=expected_hash, profile_keys=profile_keys, features=feature_set
        )
        if args.phase == "validate"
        else rows(path, expected_hash)
    )
    observed = 0
    for row in iterator:
        if args.phase == "parse":
            if (row["study_id"], row["experiment_id"], row["inc_dec"]) not in profile_keys or row[
                "trait"
            ] not in feature_set:
                raise ValueError("Synthetic coordinate mismatch")
            if probability(row["p_value"]) is None or probability(row["fdr_q"]) is None:
                raise ValueError("Synthetic probability is missing")
        consumed.update(json.dumps(row, sort_keys=True, separators=(",", ":")).encode())
        observed += 1
    elapsed = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_SELF)
    if observed != count or file_sha256(path) != expected_hash:
        raise ValueError("Measured run changed input or lost records")
    record = {
        "phase": args.phase,
        "synthetic": True,
        "profiles": args.profiles,
        "traits": args.traits,
        "native_tests": observed,
        "seed": args.seed,
        "p_order": generated["p_order"],
        "source_sha256": expected_hash,
        "source_bytes": path.stat().st_size,
        "parsed_row_sha256": consumed.hexdigest(),
        "elapsed_seconds": elapsed,
        "cpu_seconds": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
        "max_rss_bytes": int(after.ru_maxrss) * (1 if sys.platform == "darwin" else 1024),
        "rss_scope": "whole fresh process, not incremental adapter allocation",
        "baseline_scope": "streamed probability and coordinate checks, not whole old adapter",
        "helper_sha256": file_sha256(Path(__file__)),
    }
    (output / f"{args.phase}.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
