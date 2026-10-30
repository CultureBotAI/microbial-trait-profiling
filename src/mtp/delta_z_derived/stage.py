"""Publish derived analyses from a pinned completed Delta-z bundle (#1052)."""

from __future__ import annotations

import errno
import hashlib
import io
import json
import shutil
import stat
import uuid
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
from git import Repo

from mtp.delta_z_derived.tables import build_tables, trait_summary
from mtp.multiple_testing import canonical_sha256
from mtp.paths import file_sha256
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
from mtp.trait_lists import DeclaredTraitList, load_trait_list
from mtp.workflow.output_roots import temporary_root

VERSION = "delta-z-derived-v1"
MANIFEST = "delta_z_derived_run.json"
METHODS = "methods/delta_z_derived.md"
PAIR_FILE = "composition_pair_delta_z.tsv"
STUDY_FILE = "composition_study_delta_z.tsv"


def _producer(requested_commit: str | None) -> dict[str, Any]:
    """Extend the existing producer inventory to include this consumer before staging."""
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


def _durable(path: Path) -> Path:
    # Resolving first would silently erase a missing component followed by '..',
    # or accept a dangling link as a directory that this invocation may create.
    absolute = path if path.is_absolute() else Path.cwd() / path
    prefix = Path(absolute.anchor)
    for component in absolute.parts[1:]:
        if component == ".." and not prefix.is_dir():
            raise ValueError(f"Output path '..' requires an existing directory: {prefix}")
        prefix /= component
        if prefix.is_symlink():
            try:
                prefix.resolve(strict=True)
            except (OSError, RuntimeError) as error:
                raise ValueError(f"Output path contains an invalid symlink: {prefix}") from error
    if temporary_root(path) is not None:
        raise ValueError("Use a durable output directory outside system temporary roots")
    return path.resolve()


def _directory_identity(path: Path) -> tuple[int, int] | None:
    """Identify a real directory without treating a symlink as owned storage."""
    try:
        metadata = path.lstat()
    except FileNotFoundError:
        return None
    if not stat.S_ISDIR(metadata.st_mode):
        return None
    return metadata.st_dev, metadata.st_ino


def _create_parents(parent: Path, created: list[tuple[Path, tuple[int, int]]]) -> None:
    """Record each exclusive mkdir immediately, including a partially created chain."""
    missing = []
    cursor = parent
    while _directory_identity(cursor) is None:
        if cursor.exists() or cursor.is_symlink():
            raise FileExistsError(f"Output ancestor is not a directory: {cursor}")
        missing.append(cursor)
        cursor = cursor.parent
    for directory in reversed(missing):
        try:
            directory.mkdir()
        except FileExistsError:
            if _directory_identity(directory) is None:
                raise
        else:
            identity = _directory_identity(directory)
            if identity is None:
                raise FileNotFoundError(f"Created output ancestor disappeared: {directory}")
            created.append((directory, identity))


def _still_owned(path: Path, created: list[tuple[Path, tuple[int, int]]]) -> bool:
    return all(
        _directory_identity(directory) == identity
        for directory, identity in created
        if path.is_relative_to(directory)
    )


def _cleanup_failed_output(
    staging: Path,
    staging_identity: tuple[int, int] | None,
    created: list[tuple[Path, tuple[int, int]]],
) -> None:
    """Remove owned staging and empty owned ancestors; retain replacements and siblings."""
    if (
        staging_identity is not None
        and _still_owned(staging, created)
        and _directory_identity(staging) == staging_identity
    ):
        shutil.rmtree(staging)
    for directory, _identity in reversed(created):
        if not _still_owned(directory, created):
            continue
        try:
            directory.rmdir()
        except OSError as error:
            if error.errno not in (errno.ENOENT, errno.ENOTEMPTY, errno.EEXIST):
                raise


def _declared_roots(three_test: Path, comparison: Path) -> list[Path]:
    """Protect supplied aliases and declared source bundles before any filesystem write."""
    declared = json.loads(comparison.read_text(encoding="utf-8"))
    if not isinstance(declared, dict):
        raise ValueError("Comparison source record must be an object")
    manifests = [three_test, comparison]
    for field in (
        "manifest_path",
        "core_manifest_path",
        "noel_manifest_path",
        "jed_manifest_path",
        "ilyas_manifest_path",
    ):
        value = declared.get(field)
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"Comparison source record lacks {field}")
        manifests.append(Path(value))
    return _manifest_roots(manifests, multiplicity=None, label_support=None)


def _policies(top_diseases: int) -> dict[str, Any]:
    return {
        "delta_z_source": "parent_completed_ledgers_without_recomputation",
        "row1_activity_epsilon": 1e-12,
        "row1_activity_rule": "abs(study_mean_delta_z) > epsilon",
        "row1_minimum_active_studies": 50,
        "row1_noel_votes": "native_study_mean_raw_fraction_difference_signs",
        "row2_x_trait": "METPO:1000608",
        "row2_y_trait": "METPO:1000607",
        "row2_zero_policy": "exact_zero",
        "row2_sign_test": "two_sided_exact_binomial_Q4_vs_Q2",
        "row2_pearson_population": "all_unclipped_pairs_including_axes",
        "stool_metadata_field": "body_site",
        "stool_metadata_value": "Feces",
        "disease_metadata_field": "condition",
        "row3_selection": "comparison_count_descending_then_condition_ascending",
        "row3_top_diseases": top_diseases,
        "row3_test": "two_sided_wilcoxon_signed_rank_zero_method_wilcox_method_auto",
        "row3_effect_weighting": "equal_comparison",
        "row3_zero_method": "wilcox",
        "row3_correction": False,
        "row3_method_requested": "auto",
        "row3_resolved_method_field": "wilcoxon_method",
        "row3_study_sensitivity_test": "two_sided_wilcoxon_signed_rank",
        "row3_study_sensitivity_effect_weighting": "equal_weight_within_cell_study_means",
        "row3_study_sensitivity_inference_unit": "study",
        "row3_study_sensitivity_correction": "none_nominal",
        "row5_grouping": ["body_site", "condition"],
        "row5_effect_weighting": "equal_comparison",
        "row5_trait_family": "all_parent_retained_traits",
        "row2_row3_inference_unit": "comparison",
        "multiplicity_correction": "none_nominal",
        "dependence_limit": "comparisons_can_share_studies_and_controls",
    }


def _recheck(
    bundle: CompletedTable,
    additional: dict[Path, dict[str, Any]],
    producer: dict[str, Any],
    code_commit: str | None,
) -> None:
    recheck_completed_table(bundle)
    if any(_inventory(path) != identity for path, identity in additional.items()):
        raise ValueError("Derived-analysis input identity changed")
    if _producer(code_commit) != producer:
        raise ValueError("Derived-analysis producer changed")


def _input_inventory(snapshots: dict[Path, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        f"{path.name}@sha256:{identity['sha256']}": identity
        for path, identity in sorted(snapshots.items())
    }


def run_derived(
    three_test: Path,
    comparison: Path,
    output: Path,
    *,
    expected_sha256: str,
    trait_list: Path | None = None,
    summary_only: bool = False,
    code_commit: str | None = None,
    top_diseases: int = 7,
) -> dict[str, Any]:
    """Consume verified parent bytes and publish one new bundle without replaying Delta-z."""
    if output.exists() or output.is_symlink():
        raise FileExistsError("Choose a new derived-analysis output directory")
    output = _durable(output)
    if isinstance(top_diseases, bool) or not isinstance(top_diseases, int) or top_diseases < 1:
        raise ValueError("top_diseases must be a positive integer")
    if summary_only and trait_list is not None:
        raise ValueError("A summary-only run does not consume a trait list")
    if not summary_only and trait_list is None:
        raise ValueError("A full run requires an explicitly declared trait list")
    parent_path = three_test / "three_test_run.json" if three_test.is_dir() else three_test
    comparison_path = comparison / "comparison_run.json" if comparison.is_dir() else comparison
    roots = _declared_roots(parent_path, comparison_path)
    _disjoint(output, roots)
    if trait_list is not None:
        _disjoint(output, [trait_list, trait_list.resolve()])
    producer = _producer(code_commit)
    bundle = load_completed_table(three_test, comparison, expected_sha256=expected_sha256)
    roots.extend(
        [bundle.root, *_source_roots(bundle.sources, multiplicity=None, label_support=None)]
    )
    # Snapshots include auxiliary source files whose roots the comparison may not declare.
    roots.extend(_manifest_roots(list(bundle.snapshots), multiplicity=None, label_support=None))
    _disjoint(output, roots)
    experiments_path = bundle.sources.product("experiments")
    experiments_raw = experiments_path.read_bytes()
    additional = {
        experiments_path: {
            "sha256": hashlib.sha256(experiments_raw).hexdigest(),
            "size_bytes": len(experiments_raw),
        }
    }
    if bundle.snapshots.get(experiments_path.resolve()) != additional[experiments_path]:
        raise ValueError("Experiments identity differs from the verified parent source")
    declared: DeclaredTraitList | None = None
    if trait_list is not None:
        declared = load_trait_list(trait_list)
        additional[trait_list] = _inventory(trait_list)
        if additional[trait_list]["sha256"] != declared.sha256:
            raise ValueError("Declared trait list identity changed while reading")
    _probe_publication(output, roots)
    if summary_only:
        tables = {"trait_summary.tsv": trait_summary(bundle.study_delta_z, bundle.evidence)}
    else:
        if declared is None:
            raise ValueError("A full run requires an explicitly declared trait list")
        experiments = pd.read_csv(
            io.BytesIO(experiments_raw), sep="\t", dtype=str, keep_default_na=False
        )
        derived = build_tables(
            bundle.pair_delta_z,
            bundle.study_delta_z,
            bundle.evidence,
            experiments,
            declared.trait_ids,
            top_diseases=top_diseases,
        )
        tables = {
            "trait_summary.tsv": derived.trait_summary,
            "quadrant_pairs.tsv": derived.quadrant_pairs,
            "quadrant_summary.tsv": derived.quadrant_summary,
            "stool_disease_shifts.tsv": derived.stool_disease_shifts,
            "stool_study_sensitivity.tsv": derived.stool_study_sensitivity,
            "site_disease_shifts.tsv": derived.site_disease_shifts,
        }
    _recheck(bundle, additional, producer, code_commit)
    record: dict[str, Any] = {
        "version": VERSION,
        "run_id": f"delta-z-derived-{uuid.uuid4().hex}",
        "recorded_at": datetime.now(UTC).isoformat(),
        **producer,
        "parent_run_id": bundle.record["run_id"],
        "parent_manifest_sha256": expected_sha256,
        "parent_producer": {name: bundle.record[name] for name in ("code_commit", "code_dirty")},
        **{
            name: bundle.record[name]
            for name in (
                "core_run_id",
                "comparison_run_id",
                "source_run_ids",
                "source_commits",
                "strain_rollup",
                "scope_role",
                "contrast_population",
                "direction_basis",
                "analysis_bases",
                "family_size",
                "retained_traits",
            )
        },
        "summary_only": summary_only,
        "counts": {
            "traits": len(bundle.pair_delta_z.columns),
            "comparisons": len(bundle.pair_delta_z),
            "studies": len(bundle.study_delta_z),
            "traits_with_at_least_50_active_studies": int(
                tables["trait_summary.tsv"]["sdz_active_at_least_50"].sum()
            ),
        },
        "trait_list": declared.provenance() if declared is not None else None,
        "policies": _policies(top_diseases),
        "parent_input_inventory": bundle.record["inputs"],
        "inputs": _input_inventory({**bundle.snapshots, **additional}),
    }
    staging = output.parent / f".{output.name}-staging-{uuid.uuid4().hex}"
    _durable(staging)
    _disjoint(staging, roots)
    created: list[tuple[Path, tuple[int, int]]] = []
    staging_identity: tuple[int, int] | None = None
    published = False
    try:
        _create_parents(output.parent, created)
        # An exclusive mkdir never grants ownership of a colliding path.
        staging.mkdir()
        staging_identity = _directory_identity(staging)
        if staging_identity is None:
            raise FileNotFoundError(f"Created staging directory disappeared: {staging}")
        for name in (PAIR_FILE, STUDY_FILE):
            shutil.copyfile(bundle.root / name, staging / name)
        for name, table in tables.items():
            table.to_csv(staging / name, sep="\t", index=False, lineterminator="\n")
        record["outputs"] = {path.name: _inventory(path) for path in sorted(staging.iterdir())}
        for name in (PAIR_FILE, STUDY_FILE):
            expected = bundle.snapshots.get((bundle.root / name).resolve())
            if expected is None or record["outputs"][name] != expected:
                raise ValueError(f"Copied parent ledger identity differs: {name}")
        _recheck(bundle, additional, producer, code_commit)
        (staging / MANIFEST).write_text(
            json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        _publish_directory(staging, output)
        published = True
    finally:
        if not published:
            _cleanup_failed_output(staging, staging_identity, created)
    return record
