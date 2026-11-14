"""Produce synthetic verified Noel artifacts for #1029/#1030 (2026-11-08).

Inputs: prepared synthetic core tables and the repository's contribution protocol.
Outputs: baseline effects.tsv.gz plus a summary in a separate recorded arm step.
Optional altered p/q or omitted effects are deliberate, recorded test fixtures.
No KG download, live corpus read, or production arm rerun is performed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from mtp.contribution_inputs import load_base
from mtp.paths import file_sha256
from mtp.taxon_contributions import Record, write_table
from mtp.taxon_influence import effect_rows, scope_keys

INPUT_NAMES = (
    "study_trait_profiles_wide.tsv",
    "experiments.tsv",
    "bugsigdb_taxa.tsv",
    "taxa_trait_matrix.tsv",
)
PROTOCOL = Path(__file__).resolve().parents[2] / "methods/table1_taxon_influence_protocol.json"
SUMMARY = "# Synthetic native baseline\nVerified fixture binding marker.\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    alteration = parser.add_mutually_exclusive_group()
    alteration.add_argument("--alter-p", action="store_true")
    alteration.add_argument("--alter-q", action="store_true")
    parser.add_argument("--omit-effects", action="store_true")
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()
    if args.summary_only:
        if args.alter_p or args.alter_q or args.omit_effects:
            parser.error("Summary-only mode cannot alter the contribution fixture")
        with (args.output / "native_summary.md").open("x", encoding="utf-8") as handle:
            handle.write(SUMMARY)
        return

    base = load_base(*(args.inputs / name for name in INPUT_NAMES))
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    records: list[Record] = [
        {"scenario": "baseline", "scope": scope["id"], **row}
        for scope in protocol["scopes"]
        for row in effect_rows(base, base.profiles, scope_keys(base, scope))
    ]
    altered_field = "before_p" if args.alter_p else "before_q" if args.alter_q else None
    if altered_field is not None:
        eligible = [
            row
            for row in records
            if row["representation"] == "shared_fraction" and row[altered_field] != ""
        ]
        if not eligible:
            raise ValueError("Synthetic fixture has no eligible value to alter")
        row = eligible[0]
        row[altered_field] = 0.9 if float(row[altered_field]) < 0.5 else 0.1

    args.output.mkdir(parents=True, exist_ok=False)
    if not args.omit_effects:
        write_table(args.output / "effects.tsv.gz", list(records[0]), records)
    manifest = {
        "completed": True,
        "protocol": protocol,
        "traits": base.traits,
        "denominator_columns": base.columns,
        "input_sha256": {name: file_sha256(args.inputs / name) for name in INPUT_NAMES},
        "fixture_altered_field": altered_field,
        "fixture_omitted_effects": args.omit_effects,
    }
    with (args.output / "manifest.json").open("x", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
