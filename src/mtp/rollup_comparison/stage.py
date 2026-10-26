"""Publish the off-versus-phenotype comparison from pinned completed tables (#1041)."""

from __future__ import annotations

import json
import shutil
import uuid
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from git import Repo

from mtp.delta_z_derived.stage import _declared_roots, _durable, _input_inventory
from mtp.multiple_testing import canonical_sha256
from mtp.paths import file_sha256
from mtp.rollup_comparison.tables import build_tables
from mtp.rollup_comparison.validation import verify_comparable
from mtp.three_test_table.reader import (
    CompletedTable,
    load_completed_table,
    recheck_completed_table,
)
from mtp.three_test_table.stage import (
    _disjoint,
    _inventory,
    _manifest_roots,
    _probe_publication,
    _publish_directory,
    _source_roots,
)
from mtp.three_test_table.stage import _producer as _parent_producer

VERSION = "off-phenotype-comparison-v1"
MANIFEST = "rollup_comparison_run.json"
METHODS = "methods/off_phenotype_comparison.md"
EVIDENCE = "per_trait_evidence.tsv"


def _producer(requested_commit: str | None) -> dict[str, Any]:
    record = _parent_producer(requested_commit)
    package = Path(__file__).resolve().parent
    root = package.parents[2]
    inventory = dict(record["code_inventory"])
    for path in [*package.glob("*.py"), root / METHODS]:
        inventory[path.relative_to(root).as_posix()] = file_sha256(path)
    with Repo(root) as repository:
        record["code_dirty"] = repository.is_dirty(untracked_files=True)
    record["code_inventory"] = dict(sorted(inventory.items()))
    record["implementation_sha256"] = canonical_sha256(record["code_inventory"])
    record["methods_sha256"] = file_sha256(root / METHODS)
    return record


def _recheck(
    bundles: tuple[CompletedTable, CompletedTable],
    producer: dict[str, Any],
    code_commit: str | None,
) -> None:
    for bundle in bundles:
        recheck_completed_table(bundle)
    if _producer(code_commit) != producer:
        raise ValueError("Rollup comparison producer changed")


def _parent_record(bundle: CompletedTable, expected_sha256: str) -> dict[str, Any]:
    return {
        "manifest_sha256": expected_sha256,
        **{
            name: bundle.record[name]
            for name in (
                "run_id",
                "code_commit",
                "code_dirty",
                "implementation_sha256",
                "environment",
                "strain_rollup",
                "scope_role",
                "contrast_population",
                "direction_basis",
                "family_size",
                "retained_traits",
                "source_commits",
                "core_run_id",
                "comparison_run_id",
                "source_run_ids",
                "analysis_bases",
            )
        },
        "input_inventory": bundle.record["inputs"],
    }


def run_comparison(
    off: Path,
    off_comparison: Path,
    phenotype: Path,
    phenotype_comparison: Path,
    output: Path,
    *,
    off_sha256: str,
    phenotype_sha256: str,
    code_commit: str | None = None,
) -> dict[str, Any]:
    """Compare verified same-basis parents without recomputing inference or Delta-z."""
    if output.exists() or output.is_symlink():
        raise FileExistsError("Choose a new rollup-comparison output directory")
    output = _durable(output)
    declared = (
        (
            off / "three_test_run.json" if off.is_dir() else off,
            off_comparison / "comparison_run.json" if off_comparison.is_dir() else off_comparison,
        ),
        (
            phenotype / "three_test_run.json" if phenotype.is_dir() else phenotype,
            phenotype_comparison / "comparison_run.json"
            if phenotype_comparison.is_dir()
            else phenotype_comparison,
        ),
    )
    roots = [
        path for parent, comparison in declared for path in _declared_roots(parent, comparison)
    ]
    _disjoint(output, roots)
    producer = _producer(code_commit)
    bundles = (
        load_completed_table(off, off_comparison, expected_sha256=off_sha256),
        load_completed_table(phenotype, phenotype_comparison, expected_sha256=phenotype_sha256),
    )
    for bundle in bundles:
        roots.extend(
            [bundle.root, *_source_roots(bundle.sources, multiplicity=None, label_support=None)]
        )
        roots.extend(_manifest_roots(list(bundle.snapshots), multiplicity=None, label_support=None))
    _disjoint(output, roots)
    verify_comparable(*bundles)
    # Cross-parent validation binds additional source-identity artifacts.
    for bundle in bundles:
        roots.extend(_manifest_roots(list(bundle.snapshots), multiplicity=None, label_support=None))
    _disjoint(output, roots)
    _probe_publication(output, roots)
    off_bundle, phenotype_bundle = bundles
    tables = build_tables(
        off_bundle.evidence,
        phenotype_bundle.evidence,
        off_bundle.study_delta_z,
        phenotype_bundle.study_delta_z,
    )
    _recheck(bundles, producer, code_commit)
    record: dict[str, Any] = {
        "version": VERSION,
        "run_id": f"rollup-comparison-{uuid.uuid4().hex}",
        "recorded_at": datetime.now(UTC).isoformat(),
        **producer,
        "primary_scope": "phenotype",
        "sensitivity_scope": "off",
        "scope_selection": "phenotype_primary_selected_after_off_and_phenotype_results_were_seen",
        "parents": {
            "off": _parent_record(off_bundle, off_sha256),
            "phenotype": _parent_record(phenotype_bundle, phenotype_sha256),
        },
        "policies": {
            "inference": "unchanged_parent_full_family_within_each_scope",
            "shared_family_recorrection": False,
            "native_ranks": "preserved_parent_full_family",
            "view_ranks": "within_scope_direction_and_support_filter_min_ties",
            "significance_views": ["noel_by_0_05", "study_t_by_0_05"],
            "effect_views": {
                "sizes": [5, 10],
                "directions": ["positive", "negative"],
                "minimum_active_studies": [0, 50],
            },
            "activity": "abs(study_mean_delta_z) > 1e-12",
            "interval": "study_mean_delta_z_Student_t_95_percent",
            "comparison_mean": "equal_weight_comparisons_without_interval",
            "ora": "directional_only_not_a_headline_criterion",
        },
        "counts": {
            "off_traits": len(off_bundle.evidence),
            "phenotype_traits": len(phenotype_bundle.evidence),
            "shared_traits": len(tables.shared_traits),
            "union_traits": len(tables.trait_catalog),
            "comparisons": len(off_bundle.pair_delta_z),
            "studies": len(off_bundle.study_delta_z),
        },
        "inputs": {
            "off": _input_inventory(off_bundle.snapshots),
            "phenotype": _input_inventory(phenotype_bundle.snapshots),
        },
    }
    staging = output.parent / f".{output.name}-staging-{uuid.uuid4().hex}"
    _durable(staging)
    _disjoint(staging, roots)
    staging.mkdir()
    published = False
    try:
        for scope, bundle in zip(("off", "phenotype"), bundles, strict=True):
            destination = staging / f"{scope}_{EVIDENCE}"
            shutil.copyfile(bundle.root / EVIDENCE, destination)
            expected = bundle.snapshots.get((bundle.root / EVIDENCE).resolve())
            if expected is None or _inventory(destination) != expected:
                raise ValueError(f"Copied parent evidence identity differs: {scope}")
        for name, table in (
            ("shared_traits.tsv", tables.shared_traits),
            ("trait_catalog.tsv", tables.trait_catalog),
            ("rank_transitions.tsv", tables.rank_transitions),
        ):
            table.to_csv(staging / name, sep="\t", index=False, lineterminator="\n")
        record["outputs"] = {path.name: _inventory(path) for path in sorted(staging.iterdir())}
        _recheck(bundles, producer, code_commit)
        (staging / MANIFEST).write_text(
            json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
        )
        _publish_directory(staging, output)
        published = True
    finally:
        if not published:
            shutil.rmtree(staging)
    return record
