"""Report pipeline statistics for comparison against BugSigDB live counts."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

import yaml

from bsdbng.paths import DERIVED_DATA_DIR, RAW_DATA_DIR, STUDY_DATA_DIR


def print_stats(
    raw_dir: Path | None = None,
    study_dir: Path | None = None,
    derived_dir: Path | None = None,
) -> None:
    """Print data source and pipeline output statistics."""
    raw_dir = raw_dir or RAW_DATA_DIR
    study_dir = study_dir or STUDY_DATA_DIR
    derived_dir = derived_dir or DERIVED_DATA_DIR

    print("=== Data sources ===")

    # Raw CSV stats
    studies_csv = raw_dir / "studies.csv"
    if studies_csv.exists():
        with studies_csv.open(encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh))
        states = Counter(r.get("State", "").strip() for r in rows)
        print(f"  studies.csv: {len(rows)} studies")
        for state, count in states.most_common():
            print(f"    {count:>5}  State={state!r}")

    full_dump = raw_dir / "full_dump.csv"
    if full_dump.exists():
        with full_dump.open(encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
        study_ids = {r["BSDB ID"].split("/")[0] for r in rows if r.get("BSDB ID", "").strip()}
        print(f"  full_dump.csv: {len(rows)} rows, {len(study_ids)} unique studies")

    # Pipeline output stats
    print()
    print("=== Pipeline output ===")

    yaml_files = sorted(study_dir.glob("*.yaml"))
    print(f"  Study YAML files: {len(yaml_files)}")

    if yaml_files:
        total_experiments = 0
        total_signatures = 0
        total_taxa = 0
        conditions: Counter[str] = Counter()
        body_sites: Counter[str] = Counter()
        host_species: Counter[str] = Counter()

        for path in yaml_files:
            study = yaml.safe_load(path.read_text())
            for exp in study.get("experiments", []):
                total_experiments += 1
                cond = exp.get("condition")
                if cond:
                    conditions[cond] += 1
                site = exp.get("body_site")
                if site:
                    body_sites[site] += 1
                host = exp.get("host_species")
                if host:
                    host_species[host] += 1
                for sig in exp.get("signatures", []):
                    total_signatures += 1
                    total_taxa += len(sig.get("taxa", []))

        print(f"  Experiments: {total_experiments}")
        print(f"  Signatures: {total_signatures}")
        print(f"  Taxa (total, with duplicates): {total_taxa}")
        print(f"  Unique conditions: {len(conditions)}")
        print(f"  Unique body sites: {len(body_sites)}")
        print(f"  Unique host species: {len(host_species)}")
        print()
        print("  Top 5 conditions:")
        for cond, count in conditions.most_common(5):
            print(f"    {count:>5}  {cond}")
        print()
        print("  Top 5 body sites:")
        for site, count in body_sites.most_common(5):
            print(f"    {count:>5}  {site}")

    # Ingest log stats
    log_path = derived_dir / "ingest_log.json"
    if log_path.exists():
        log = json.loads(log_path.read_text())
        skips = [e for e in log if e.get("level") == "skip"]
        infos = [e for e in log if e.get("level") == "info"]
        print()
        print(f"  Ingest log: {len(skips)} skips, {len(infos)} info entries")

    print()
    print("=== BugSigDB live (for comparison) ===")
    print("  Check https://bugsigdb.org/Main_Page for authoritative counts")
