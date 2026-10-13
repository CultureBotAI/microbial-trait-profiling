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

ARSENATE = "METPO:2000103|pathways:arsenate_reduction"
SPHERE = "METPO:1000683"
ESCHERICHIA_SCENARIO = "exclude_genus_NCBITaxon:561"


def _one(records: list[dict[str, str]], **criteria: str) -> dict[str, str]:
    matches = [
        record
        for record in records
        if all(record.get(field) == value for field, value in criteria.items())
    ]
    if len(matches) != 1:
        description = ", ".join(f"{field}={value!r}" for field, value in criteria.items())
        raise ValueError(f"Expected one row for {description}, found {len(matches)}")
    return matches[0]


def interpretation_handoff(source: Path) -> list[dict[str, str]]:
    """Extract the two #213 interpretation flags without merging their estimands."""
    effects = list(compressed_rows(source / "effects.tsv.gz"))
    native = list(rows(source / "native_composition_comparison.tsv"))
    contributors = list(rows(source / "contributing_universe.tsv"))

    shared_criteria = {
        "scope": "pooled",
        "representation": "shared_fraction",
        "denominator_trait": "",
    }
    sphere = _one(effects, scenario="baseline", trait_key=SPHERE, **shared_criteria)
    arsenate = _one(
        effects,
        scenario=ESCHERICHIA_SCENARIO,
        trait_key=ARSENATE,
        **shared_criteria,
    )
    sphere_native = _one(
        native,
        scenario="baseline",
        scope="pooled",
        trait_key=SPHERE,
        reference="refitted",
    )
    arsenate_native = _one(
        native,
        scenario=ESCHERICHIA_SCENARIO,
        scope="pooled",
        trait_key=ARSENATE,
        reference="refitted",
    )
    escherichia = _one(
        contributors,
        scope="pooled",
        direction="inc",
        trait_key=ARSENATE,
        rank="genus",
        taxon="NCBITaxon:561",
        name="Escherichia",
    )

    common = {
        "scope": "pooled",
        "raw_fraction_weighting": "equal experiments within study; equal studies",
        "raw_fraction_unit": "fraction_difference",
        "raw_fraction_interval": "conditional Student t 95% interval for the study mean",
        "native_weighting": sphere_native["native_before_weighting"],
        "native_unit": sphere_native["native_before_unit"],
        "q_threshold": "0.05",
        "integration_rule": (
            "interpretation flag only; preserve raw-fraction and native results as separate "
            "method rows"
        ),
    }
    return [
        {
            "finding_id": "arsenate_escherichia_sensitive",
            "trait_key": ARSENATE,
            "trait_label": "arsenate reduction",
            "scenario": ESCHERICHIA_SCENARIO,
            "excluded_taxon": "Escherichia (NCBITaxon:561)",
            "increased_arm_count_share": escherichia["share"],
            "raw_fraction_baseline_mean": arsenate["before_effect"],
            "raw_fraction_baseline_ci_low": arsenate["before_mean_ci_low"],
            "raw_fraction_baseline_ci_high": arsenate["before_mean_ci_high"],
            "raw_fraction_after_mean": arsenate["after_effect"],
            "raw_fraction_after_ci_low": arsenate["after_mean_ci_low"],
            "raw_fraction_after_ci_high": arsenate["after_mean_ci_high"],
            "raw_fraction_matched_change": arsenate["matched_change_effect"],
            "raw_fraction_matched_change_ci_low": arsenate["matched_change_mean_ci_low"],
            "raw_fraction_matched_change_ci_high": arsenate["matched_change_mean_ci_high"],
            "directional_correction": arsenate["correction"],
            "directional_baseline_q": arsenate["before_q"],
            "directional_after_q": arsenate["after_q"],
            "native_family_size": arsenate_native["native_before_family_tests"],
            "native_baseline_q": arsenate_native["native_before_q"],
            "native_after_q": arsenate_native["native_after_q"],
            "flag_status": "lineage_sensitive",
            "manuscript_boundary": (
                "Describe an Escherichia-sensitive annotation signal on this corpus, not a "
                "lineage-independent arsenate mechanism; taxon removal is not causal attribution."
            ),
            **common,
        },
        {
            "finding_id": "sphere_raw_fraction_mean_interval_crosses_zero",
            "trait_key": SPHERE,
            "trait_label": "sphere shaped",
            "scenario": "baseline",
            "excluded_taxon": "",
            "increased_arm_count_share": "",
            "raw_fraction_baseline_mean": sphere["before_effect"],
            "raw_fraction_baseline_ci_low": sphere["before_mean_ci_low"],
            "raw_fraction_baseline_ci_high": sphere["before_mean_ci_high"],
            "raw_fraction_after_mean": "",
            "raw_fraction_after_ci_low": "",
            "raw_fraction_after_ci_high": "",
            "raw_fraction_matched_change": "",
            "raw_fraction_matched_change_ci_low": "",
            "raw_fraction_matched_change_ci_high": "",
            "directional_correction": sphere["correction"],
            "directional_baseline_q": sphere["before_q"],
            "directional_after_q": "",
            "native_family_size": sphere_native["native_before_family_tests"],
            "native_baseline_q": sphere_native["native_before_q"],
            "native_after_q": "",
            "flag_status": "cross_scale_non_reproduction",
            "manuscript_boundary": (
                "The study-weighted raw-fraction mean interval crosses zero, while the separate "
                "directional sign test and native transformed test remain significant; report the "
                "estimand, weighting and test rather than treating either result as replication or "
                "refutation of the other."
            ),
            **common,
        },
    ]


def _publish_table(path: Path, fields: list[str], records: list[dict[str, str]]) -> int:
    """Keep an identical existing table byte-for-byte, including its gzip header."""
    if path.exists():
        existing = list(compressed_rows(path)) if path.suffix == ".gz" else list(rows(path))
        if existing == records:
            return len(records)
    return write_table(path, fields, records)


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
    if native.get("code_dirty", True):
        raise ValueError("Review bundles cannot be published from a dirty or unrecorded tree")
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
    review_tables: list[tuple[str, list[dict[str, str]]]] = [
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
        ("interpretation_handoff.tsv", interpretation_handoff(source)),
    ]
    for filename, records in review_tables:
        if not records:
            raise ValueError(f"Missing review table: {filename}")
        path = destination / filename
        number = _publish_table(path, list(records[0]), records)
        published[filename] = {"rows": number, "sha256": digest(path)}
    summary = (source / "native_summary.md").read_text(encoding="utf-8")
    for original, review in (
        ("native_manifest.json", "run_manifest.json"),
        ("native_composition_comparison.tsv", "native_composition_comparison.tsv.gz"),
    ):
        summary = summary.replace(original, review)
    flags = {row["finding_id"]: row for row in review_tables[-1][1]}
    arsenate = flags["arsenate_escherichia_sensitive"]
    sphere = flags["sphere_raw_fraction_mean_interval_crosses_zero"]
    summary += (
        "\n## Cross-scale interpretation flags (#213)\n\n"
        "The study-weighted raw-fraction mean and the native transformed analysis are "
        "different estimands. The raw-fraction directional q-value is a third statistic: "
        "it tests the share of non-tied study directions, not whether the mean-effect "
        "interval excludes zero. The machine-readable `interpretation_handoff.tsv` keeps "
        "those fields separate for #178.\n\n"
        "| Trait | Study-weighted raw-fraction finding | Native transformed finding | "
        "Interpretation |\n"
        "| --- | --- | --- | --- |\n"
        f"| Arsenate reduction | Excluding Escherichia changes the matched mean by "
        f"{float(arsenate['raw_fraction_matched_change']) * 100:.4f} pp "
        f"[{float(arsenate['raw_fraction_matched_change_ci_low']) * 100:.4f}, "
        f"{float(arsenate['raw_fraction_matched_change_ci_high']) * 100:.4f}]; the remaining "
        f"mean is {float(arsenate['raw_fraction_after_mean']) * 100:.4f} pp "
        f"[{float(arsenate['raw_fraction_after_ci_low']) * 100:.4f}, "
        f"{float(arsenate['raw_fraction_after_ci_high']) * 100:.4f}]. Escherichia supplies "
        f"{float(arsenate['increased_arm_count_share']) * 100:.1f}% of increased-arm counts. | "
        f"q changes from {float(arsenate['native_baseline_q']):.6g} to "
        f"{float(arsenate['native_after_q']):.6g} after the genus exclusion. | "
        "Lineage-sensitive on this corpus; do not describe a lineage-independent mechanism "
        "or treat removal as causal attribution. |\n"
        f"| Sphere shaped | Baseline mean {float(sphere['raw_fraction_baseline_mean']) * 100:.4f} "
        f"pp [{float(sphere['raw_fraction_baseline_ci_low']) * 100:.4f}, "
        f"{float(sphere['raw_fraction_baseline_ci_high']) * 100:.4f}] crosses zero; the separate "
        f"directional sign-test q is {float(sphere['directional_baseline_q']):.6g}. | "
        f"Baseline native q is {float(sphere['native_baseline_q']):.6g}. | "
        "Cross-scale non-reproduction of the mean interval, not a refutation or independent "
        "replication of the native result. |\n"
    )
    summary += (
        "\n## Review Bundle\n\n"
        "The `.tsv.gz` comparison files contain all six traits, scenarios and scopes, "
        "including null and untestable results. `native_table1_stats.tsv.gz` is the "
        "six-trait subset of the native composition table; its q-values were computed "
        "over the FULL native families, whose sizes remain attached. The separate ORA "
        "family table includes both backgrounds and both tail-family sizes. "
        "`interpretation_handoff.tsv` carries the two measured #213 flags without "
        "collapsing their scales or tests. "
        "`scenarios.tsv` supplies lineage names, support and selection status for "
        "the full genus/family inventory, including lineages outside the panel.\n\n"
        "The full 317-trait composition table, profile-level ORA results and source "
        "ledger are generated under `data/derived/table1_taxon_influence/`, following "
        "the repository's untracked generated-data convention. Regenerate the default "
        "input chain with `just table1-taxon-influence`, or use `just table1-taxon-influence-from "
        "<pinned-KG-directory>` for already materialized inputs. The recipe writes a candidate "
        "review bundle under the untracked run directory; this committed directory is "
        "an explicitly reviewed snapshot. This compact review bundle "
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
        "interpretation_source_sha256": {
            filename: digest(source / filename)
            for filename in (
                "effects.tsv.gz",
                "contributing_universe.tsv",
                "native_composition_comparison.tsv",
            )
        },
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
    parser.add_argument(
        "--destination",
        type=Path,
        help="Output directory. Defaults to <source>/review under untracked derived data.",
    )
    args = parser.parse_args()
    publish(args.source, args.destination or args.source / "review")


if __name__ == "__main__":
    main()


__all__ = [
    "interpretation_handoff",
    "main",
    "publish",
    "reconcile_exported_source_counts",
]
