"""Verify a completed native run and publish compact, lossless review tables."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from mtp.contribution_inputs import count, digest, rows
from mtp.native_taxon_influence import compressed_rows
from mtp.paths import REPO_ROOT as ROOT
from mtp.taxon_contributions import write_table


def reconcile_exported_source_counts(source: Path) -> dict[str, object]:
    """Re-sum the source ledger and require equality with every profile count."""
    source_path = source / "source_contributions.tsv.gz"
    links_path = source / "profile_reported_taxa.tsv.gz"
    profile_path = source / "profile_counts.tsv.gz"
    source_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    source_rows = 0
    for row in compressed_rows(source_path):
        source_rows += 1
        source_counts[row["reported_taxon"]][row["trait_key"]] += count(row["count_contribution"])

    reconstructed: Counter[tuple[str, str, str, str]] = Counter()
    link_rows = 0
    for row in compressed_rows(links_path):
        link_rows += 1
        prefix = (row["study_id"], row["experiment_id"], row["direction"])
        for trait, contribution in source_counts[row["reported_taxon"]].items():
            reconstructed[(*prefix, trait)] += contribution

    expected: dict[tuple[str, str, str, str], int] = {}
    profile_rows = 0
    for row in compressed_rows(profile_path):
        profile_rows += 1
        key = (row["study_id"], row["experiment_id"], row["direction"], row["trait_key"])
        if key in expected:
            raise ValueError(f"Duplicate exported profile count: {key}")
        expected[key] = count(row["profile_count"])
    extra = set(reconstructed).difference(expected)
    mismatches = [
        key for key, expected_count in expected.items() if reconstructed[key] != expected_count
    ]
    if extra or mismatches:
        raise ValueError(
            "Exported source contributions do not reconcile with profile counts: "
            f"{len(mismatches)} mismatches, {len(extra)} unexpected keys"
        )
    return {
        "source_contribution_rows": source_rows,
        "profile_reported_taxon_rows": link_rows,
        "profile_count_rows": profile_rows,
        "mismatches": 0,
        "sha256": {path.name: digest(path) for path in (source_path, links_path, profile_path)},
    }


def publish(source: Path, destination: Path) -> dict[str, object]:
    native = json.loads((source / "native_manifest.json").read_text(encoding="utf-8"))
    audit = json.loads((source / "manifest.json").read_text(encoding="utf-8"))
    if not native["completed"] or not audit["completed"] or audit["mismatches"]:
        raise ValueError("Only completed, reconciled runs may be published")
    if digest(source / "manifest.json") != native["contribution_manifest_sha256"]:
        raise ValueError("The parent contribution manifest has changed")
    source_reconciliation = reconcile_exported_source_counts(source)
    for group in ("input_sha256", "output_sha256"):
        for filename, expected in native["results"][group].items():
            path = (source / filename).resolve()
            if not path.is_relative_to(source.resolve()):
                raise ValueError("Native artifact path escapes the run directory")
            if digest(path) != expected:
                raise ValueError(f"Native run artifact changed: {filename}")
    destination.mkdir(parents=True, exist_ok=True)
    six = set(audit["traits"].values())
    published: dict[str, object] = {}
    for filename, records in (
        (
            "native_table1_stats.tsv.gz",
            [
                r
                for r in compressed_rows(source / "native_composition.tsv.gz")
                if r["trait_key"] in six
            ],
        ),
        (
            "native_composition_comparison.tsv.gz",
            list(rows(source / "native_composition_comparison.tsv")),
        ),
        (
            "native_ora_fold_comparison.tsv.gz",
            list(rows(source / "native_ora_fold_comparison.tsv")),
        ),
        ("native_ora_families.tsv", list(rows(source / "native_ora_families.tsv"))),
        ("scenarios.tsv", list(rows(source / "scenarios.tsv"))),
    ):
        if not records:
            raise ValueError(f"Missing review table: {filename}")
        path = destination / filename
        number = write_table(path, list(records[0]), records)
        published[filename] = {"rows": number, "sha256": digest(path)}
    summary = (source / "native_summary.md").read_text(encoding="utf-8")
    for original, review in (
        ("native_manifest.json", "run_manifest.json"),
        ("native_composition_comparison.tsv", "native_composition_comparison.tsv.gz"),
    ):
        summary = summary.replace(original, review)
    summary += (
        "\n## Review Bundle\n\n"
        "The `.tsv.gz` comparison files contain all six traits, scenarios and scopes, "
        "including null and untestable results. `native_table1_stats.tsv.gz` is the "
        "six-trait subset of the native composition table; its q-values were computed "
        "over the FULL native families, whose sizes remain attached. The separate ORA "
        "family table includes both backgrounds and both tail-family sizes. "
        "`scenarios.tsv` supplies lineage names, support and selection status for "
        "the full genus/family inventory, including lineages outside the panel.\n\n"
        "The full 317-trait composition table, profile-level ORA results and source "
        "ledger are generated under `data/derived/table1_taxon_influence/`, following "
        "the repository's untracked generated-data convention. Regenerate them with "
        "`just table1-taxon-influence <pinned-KG-directory>`. This compact review bundle "
        "does not discard any six-trait influence comparison. `run_manifest.json` "
        "records hashes for the full data; `publication_manifest.json` identifies "
        "these review copies. All full-run file hashes were checked before publication.\n\n"
        "Publication also re-sums the complete source-contribution ledger through the "
        "profile/reported-taxon links and requires exact equality with every exported "
        "profile/trait count. The publication manifest records those file hashes and row "
        "counts.\n"
    )
    (destination / "README.md").write_text(summary, encoding="utf-8")
    (destination / "run_manifest.json").write_text(
        json.dumps(native, indent=2) + "\n", encoding="utf-8"
    )
    (destination / "contribution_manifest.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    result: dict[str, object] = {
        "native_manifest_sha256": digest(source / "native_manifest.json"),
        "publisher_sha256": digest(Path(__file__)),
        "source_reconciliation": source_reconciliation,
        "tables": published,
    }
    (destination / "publication_manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "data/derived/table1_taxon_influence")
    parser.add_argument("--destination", type=Path, default=ROOT / "methods/table1-native-results")
    args = parser.parse_args()
    publish(args.source, args.destination)


if __name__ == "__main__":
    main()


__all__ = [
    "main",
    "publish",
    "reconcile_exported_source_counts",
]
