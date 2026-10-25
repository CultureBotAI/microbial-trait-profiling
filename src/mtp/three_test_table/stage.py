"""Verified additive per-trait table, portable provenance and strict replay (#1044)."""

from __future__ import annotations

import ctypes
import json
import os
import platform
import shutil
import sys
import tempfile
import uuid
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy
from git import Repo

from mtp.count_input import parse_nonnegative_integer
from mtp.datamodel.multiple_testing_pydantic import PostArmMultiplicityRun
from mtp.datamodel.workflow_pydantic import ComparisonRun, ContrastPopulation, DirectionBasis
from mtp.multiple_testing import canonical_sha256
from mtp.multiplicity.families import rows
from mtp.multiplicity.sources import VerifiedSources, bind_recorded_hashes, open_verified_sources
from mtp.multiplicity.stage import verify_post_arm
from mtp.paths import file_sha256
from mtp.three_test_table.carriers import carrier_table
from mtp.three_test_table.composition import composition_table
from mtp.three_test_table.native import native_tables, validate_native_overlap
from mtp.three_test_table.ora import ora_table

VERSION = "three-test-per-trait-v1"
MANIFEST = "three_test_run.json"
METHODS = "methods/three_test_table.md"
CALIBRATION_STATUS = "uncalibrated (methods/calibration.md excludes native ORA)"


def require_durable(path: Path) -> Path:
    """Reject declared temporary roots and aliases under the #1049 contract.

    This matches the shared guard proposed in #1063 until it can be imported.
    Bandit B108 prevents some fixed POSIX temporary-root literals (#1051).
    An unlisted root is not refused unless gettempdir() or TMPDIR identifies it.
    Explicit scratch must also be outside those inherited temporary roots.
    """
    resolved = path.resolve()
    directories = [
        Path("/private/tmp"),
        Path("/var/folders"),
        Path("/private/var/folders"),
        Path(tempfile.gettempdir()),
    ]
    if environment := os.environ.get("TMPDIR"):
        directories.append(Path(environment))
    roots = {
        spelling
        for directory in directories
        for spelling in (Path(os.path.abspath(directory)), directory.resolve())
    }
    if any(
        spelling.is_relative_to(root)
        for spelling in (Path(os.path.abspath(path)), resolved)
        for root in roots
    ):
        raise ValueError("Use a durable directory for outputs and scratch")
    return resolved


def _disjoint(path: Path, roots: list[Path]) -> None:
    resolved = path.resolve()
    for root in roots:
        root = root.resolve()
        if resolved.is_relative_to(root) or root.is_relative_to(resolved):
            raise ValueError("Sidecar outputs and scratch must be disjoint from source bundles")


def _manifest_roots(
    manifests: list[Path], *, multiplicity: Path | None, label_support: Path | None
) -> list[Path]:
    """Protect file aliases and their targets; directory inputs protect only themselves."""
    files = list(manifests)
    roots: list[Path] = []
    for optional in (multiplicity, label_support):
        if optional is not None:
            if optional.is_dir():
                roots.append(optional)
            else:
                files.append(optional)
    for manifest in files:
        # Resolve the file before taking its parent: a file symlink can point
        # outside its alias directory, and both locations need write protection.
        roots.extend((manifest.parent, manifest.resolve().parent))
    return roots


def _source_roots(
    sources: VerifiedSources, *, multiplicity: Path | None, label_support: Path | None
) -> list[Path]:
    return _manifest_roots(
        [
            Path(sources.core.manifest_path),
            Path(sources.comparison.manifest_path),
            *(Path(arm.manifest_path) for arm in sources.arms.values()),
        ],
        multiplicity=multiplicity,
        label_support=label_support,
    )


@contextmanager
def durable_scratch(path: Path) -> Iterator[None]:
    """Route existing adapters' tempfile use into the declared durable scratch root."""
    path = require_durable(path)
    path.mkdir(parents=True, exist_ok=True)
    previous_env, previous_cache = os.environ.get("TMPDIR"), tempfile.tempdir
    os.environ["TMPDIR"] = str(path)
    tempfile.tempdir = str(path)
    try:
        yield
    finally:
        tempfile.tempdir = previous_cache
        if previous_env is None:
            os.environ.pop("TMPDIR", None)
        else:
            os.environ["TMPDIR"] = previous_env


def _producer(requested_commit: str | None) -> dict[str, Any]:
    package = Path(__file__).resolve().parents[1]
    with Repo(package, search_parent_directories=True) as repository:
        if repository.working_tree_dir is None:
            raise ValueError("The table producer requires a Git worktree")
        root = Path(repository.working_tree_dir).resolve()
        commit = repository.head.commit.hexsha
        if requested_commit is not None and requested_commit != commit:
            raise ValueError("Requested code commit differs from the executing checkout")
        paths: set[Path] = set()
        for name, stage in repository.index.entries:
            if str(name).startswith("src/mtp/") and str(name).endswith(".py"):
                if stage:
                    raise ValueError("Resolve tracked source conflicts before execution")
                paths.add(root / name)
        # New package files must be bound even before their first commit.
        paths.update(Path(__file__).parent.glob("*.py"))
        paths.update(root / name for name in (METHODS, "uv.lock", "pyproject.toml", "Justfile"))
        inventory = {path.relative_to(root).as_posix(): file_sha256(path) for path in sorted(paths)}
        return {
            "code_commit": commit,
            "code_dirty": repository.is_dirty(untracked_files=False),
            "implementation_sha256": canonical_sha256(inventory),
            "code_inventory": inventory,
            "lock_sha256": file_sha256(root / "uv.lock"),
            "methods_sha256": file_sha256(root / METHODS),
            "environment": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "pandas": pd.__version__,
            },
        }


def _inventory(path: Path) -> dict[str, Any]:
    return {"sha256": file_sha256(path), "size_bytes": path.stat().st_size}


def _publish_directory(staging: Path, destination: Path) -> None:
    """Atomically publish without replacing even an empty competing directory.

    Darwin renamex_np(RENAME_EXCL) and Linux renameat2(RENAME_NOREPLACE)
    provide this guarantee; ordinary POSIX rename does not. An unavailable
    platform, libc symbol or filesystem operation fails closed, with no
    nonexclusive fallback. Staging and destination must share a filesystem.
    """
    library = ctypes.CDLL(None, use_errno=True)
    try:
        if sys.platform == "darwin":
            rename = library.renamex_np
            rename.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
            rename.restype = ctypes.c_int
            # Darwin sys/stdio.h: RENAME_EXCL = 0x00000004.
            result = rename(os.fsencode(staging), os.fsencode(destination), 4)
        elif sys.platform == "linux":
            rename = library.renameat2
            rename.argtypes = [
                ctypes.c_int,
                ctypes.c_char_p,
                ctypes.c_int,
                ctypes.c_char_p,
                ctypes.c_uint,
            ]
            rename.restype = ctypes.c_int
            # Linux uapi/linux/{fcntl,fs}.h: AT_FDCWD = -100; RENAME_NOREPLACE = 1.
            result = rename(-100, os.fsencode(staging), -100, os.fsencode(destination), 1)
        else:
            raise NotImplementedError("Atomic exclusive directory rename requires Darwin or Linux")
    except AttributeError as error:
        raise NotImplementedError(
            "Atomic exclusive directory rename is unavailable in libc"
        ) from error
    if result != 0:
        error_number = ctypes.get_errno()
        raise OSError(error_number, os.strerror(error_number), str(destination))


def _declared_source_roots(
    comparison: Path, *, multiplicity: Path | None, label_support: Path | None
) -> list[Path]:
    """Read only declared locations for early write protection, never source verification."""
    record = ComparisonRun.model_validate_json(comparison.read_text(encoding="utf-8"))
    return _manifest_roots(
        [
            comparison,
            Path(record.manifest_path),
            Path(record.core_manifest_path),
            Path(record.noel_manifest_path),
            Path(record.jed_manifest_path),
            Path(record.ilyas_manifest_path),
        ],
        multiplicity=multiplicity,
        label_support=label_support,
    )


def _probe_publication(output: Path, protected: list[Path]) -> None:
    """Probe both no-replace refusal and a real directory move on the output filesystem."""
    output.parent.mkdir(parents=True, exist_ok=True)
    probe = require_durable(output.parent / f".{output.name}-probe-{uuid.uuid4().hex}")
    _disjoint(probe, protected)
    # An exclusive mkdir reserves only our private probe; the requested output
    # and any existing user directory are never used as a disposable target.
    probe.mkdir(mode=0o700)
    try:
        source, target = probe / "source", probe / "target"
        source.mkdir()
        target.mkdir()
        original_source, original_target = source.stat(), target.stat()
        try:
            _publish_directory(source, target)
        except FileExistsError:
            pass
        else:
            raise RuntimeError("Filesystem probe did not enforce exclusive rename")
        if (
            not source.is_dir()
            or not target.is_dir()
            or (source.stat().st_dev, source.stat().st_ino)
            != (original_source.st_dev, original_source.st_ino)
            or (target.stat().st_dev, target.stat().st_ino)
            != (original_target.st_dev, original_target.st_ino)
        ):
            raise RuntimeError("Filesystem probe changed an existing rename target")
        target.rmdir()
        _publish_directory(source, target)
        if (
            source.exists()
            or not target.is_dir()
            or (target.stat().st_dev, target.stat().st_ino)
            != (original_source.st_dev, original_source.st_ino)
        ):
            raise RuntimeError("Filesystem probe did not complete the exclusive rename")
    finally:
        # All paths are private, empty directories created above. No recursive
        # deletion can reach a pre-existing directory or a source bundle.
        for name in ("source", "target"):
            member = probe / name
            if member.exists():
                member.rmdir()
        probe.rmdir()


def _inputs(paths: set[Path], sources: VerifiedSources) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for path in sorted(paths):
        item = _inventory(path)
        expected = sources.recorded_hashes.get(path.resolve())
        if expected is not None and item["sha256"] != expected:
            raise ValueError(f"Source bytes differ from recorded hash: {path.name}")
        identity = f"{path.name}@sha256:{item['sha256']}"
        result[identity] = item
    return result


def _check_multiplicity_inputs(root: Path, bindings: dict[Path, str]) -> None:
    """Require the verified optional inventory, including after other adapters run (#1059)."""
    members = list(root.rglob("*"))
    if any(path.is_symlink() for path in members) or {
        path for path in members if path.is_file()
    } != set(bindings):
        raise ValueError("Multiplicity input inventory changed after verification")
    for path, expected in bindings.items():
        if file_sha256(path) != expected:
            raise ValueError("Multiplicity input hash changed after verification")


def _bind_multiplicity_inputs(path: Path, comparison_sha256: str) -> tuple[Path, dict[Path, str]]:
    """Bind the verified record's hashes, never a later unverified disk snapshot (#1059)."""
    manifest = (path / "multiplicity_run.json" if path.is_dir() else path).resolve()
    digest = file_sha256(manifest)
    record = verify_post_arm(path)
    if file_sha256(manifest) != digest:
        raise ValueError("Multiplicity manifest changed during verification")
    if PostArmMultiplicityRun.model_validate_json(manifest.read_text(encoding="utf-8")) != record:
        raise ValueError("Multiplicity verified record differs from the bound manifest")
    if record.comparison_manifest_sha256 != comparison_sha256:
        raise ValueError("Multiplicity sidecar belongs to a different comparison")
    root = manifest.parent
    bindings = {manifest: digest}
    artifacts = [
        *record.reports,
        *(
            artifact
            for family in record.families
            for artifact in (family.declared, family.raw, family.corrected)
        ),
    ]
    for artifact in artifacts:
        relative = Path(artifact.identity)
        target = root / relative
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or relative.as_posix() != artifact.identity
            or not target.resolve().is_relative_to(root)
            or target in bindings
        ):
            raise ValueError("Multiplicity artifact identity is duplicate or escapes its bundle")
        bindings[target] = artifact.sha256
    _check_multiplicity_inputs(root, bindings)
    return root, bindings


def composition_profiles(sources: VerifiedSources) -> pd.DataFrame:
    """Read exact bound count bytes and reproduce the human-count arm input selection."""
    profile = sources.product("count_profiles")
    experiments = sources.product("experiments")
    features = sources.features["count_profiles"]
    hosts = {}
    for row in rows(experiments, sources.expected_digest(experiments)):
        key = row["study_id"], row["experiment_id"]
        if key in hosts or not all(key):
            raise ValueError("Missing or duplicate experiment identity")
        hosts[key] = row["host_species"]
    records = []
    identifiers = ("study_id", "experiment_id", "inc_dec")
    for index, row in enumerate(rows(profile, sources.expected_digest(profile)), start=2):
        key = row["study_id"], row["experiment_id"]
        if key not in hosts:
            raise ValueError("Count profile lacks experiment host metadata")
        if row["inc_dec"] not in {"inc", "dec"}:
            raise ValueError("Count profile has a noncanonical native population direction")
        if hosts[key].strip() != "Homo sapiens":
            continue
        record: dict[str, Any] = {name: row[name] for name in identifiers}
        for feature in features:
            record[feature] = parse_nonnegative_integer(
                row[feature], path=profile, line_number=index, column_name=feature
            )
        records.append(record)
    return pd.DataFrame.from_records(records, columns=[*identifiers, *features])


def _labels(sources: VerifiedSources) -> pd.DataFrame:
    path = sources.product("feature_metadata")
    labels = {}
    for row in rows(path, sources.expected_digest(path)):
        identity = row["column_id"]
        if identity in labels:
            raise ValueError("Duplicate feature metadata identity")
        labels[identity] = row["object_name"] or row["column_label"] or identity
    records = []
    for feature in sources.features["count_profiles"]:
        if feature not in labels:
            raise ValueError("Retained trait lacks feature metadata")
        label = labels[feature]
        object_curie = feature.rsplit("|", 1)[-1]
        if object_curie == "METPO:1000603" and label != "anaerobic":
            raise ValueError("METPO:1000603 must retain its anaerobic identity")
        if object_curie == "METPO:1000607" and label != "obligately anaerobic":
            raise ValueError("METPO:1000607 must retain its obligately anaerobic identity")
        records.append({"trait": feature, "trait_label": label})
    return pd.DataFrame(records)


def _merge_tables(features: tuple[str, ...], tables: list[pd.DataFrame]) -> pd.DataFrame:
    result = pd.DataFrame({"trait": features})
    for table in tables:
        if table.trait.duplicated().any() or set(table.trait) != set(features):
            raise ValueError("Per-trait component does not cover the complete catalog exactly once")
        if (set(table.columns) & set(result.columns)) - {"trait"}:
            raise ValueError("Per-trait components have ambiguous overlapping columns")
        result = result.merge(table, on="trait", how="left", validate="one_to_one", sort=False)
    return result


def run_table(
    comparison: Path,
    output: Path,
    *,
    scratch: Path,
    code_commit: str | None = None,
    multiplicity: Path | None = None,
    label_support: Path | None = None,
    recorded_at: str | None = None,
    run_id: str | None = None,
) -> dict[str, Any]:
    """Build a new full-catalog table without modifying any input or native call."""
    output, scratch = require_durable(output), require_durable(scratch)
    if output.exists() or output.is_symlink():
        raise FileExistsError("Choose a new sidecar output directory")
    _disjoint(output, [scratch])
    producer = _producer(code_commit)
    comparison_path = comparison / "comparison_run.json" if comparison.is_dir() else comparison
    declared_roots = _declared_source_roots(
        comparison_path, multiplicity=multiplicity, label_support=label_support
    )
    _disjoint(output, declared_roots)
    _disjoint(scratch, declared_roots)
    _probe_publication(output, [scratch, *declared_roots])
    comparison_hash = file_sha256(comparison_path)
    sources = bind_recorded_hashes(open_verified_sources(comparison_path), comparison_hash)
    roots = _source_roots(sources, multiplicity=multiplicity, label_support=label_support)
    _disjoint(output, roots)
    _disjoint(scratch, roots)
    features = sources.features["count_profiles"]
    if set(features) != set(sources.features["ora_profiles"]):
        raise ValueError("Count and ORA retained catalogs differ")
    paths = {
        Path(sources.core.manifest_path),
        comparison_path,
        *(Path(arm.manifest_path) for arm in sources.arms.values()),
        *(
            sources.product(name)
            for name in (
                "count_profiles",
                "ora_profiles",
                "experiments",
                "taxa_matrix",
                "feature_metadata",
                "feature_catalog",
            )
        ),
    }
    baseline = _inputs(paths, sources)
    multiplicity_inputs: tuple[Path, dict[Path, str]] | None = None
    with durable_scratch(scratch):
        if multiplicity is not None:
            multiplicity_inputs = _bind_multiplicity_inputs(multiplicity, comparison_hash)
            _, bindings = multiplicity_inputs
            for path, expected in bindings.items():
                previous = sources.recorded_hashes.get(path)
                if previous is not None and previous != expected:
                    raise ValueError("Multiplicity input has conflicting recorded hashes")
            sources = replace(sources, recorded_hashes={**sources.recorded_hashes, **bindings})
            paths.update(bindings)
        native, native_paths = native_tables(sources)
        ora, ora_paths = ora_table(sources)
        carriers, carrier_paths = carrier_table(sources, label_support)
        paths.update((*native_paths, *ora_paths, *carrier_paths))
        consumed = _inputs(paths, sources)
        composition = composition_table(composition_profiles(sources), features)
        native = validate_native_overlap(composition.table, native)
        table = _merge_tables(
            features, [_labels(sources), composition.table, native, ora, carriers]
        )
    if any(consumed.get(key) != value for key, value in baseline.items()):
        raise ValueError("Source input identity changed during adaptation")
    if _inputs(paths, sources) != consumed or _producer(code_commit) != producer:
        raise ValueError("Source or producer changed during table generation")
    if multiplicity_inputs is not None:
        _check_multiplicity_inputs(*multiplicity_inputs)
    population = ContrastPopulation(sources.comparison.contrast_population).value
    direction = DirectionBasis(sources.comparison.direction_basis).value
    table["strain_rollup"] = sources.scope
    table["contrast_population"] = population
    table["direction_basis"] = direction
    table["core_run_id"] = sources.core.run_id
    table["primary_ranking"] = "undecided"
    table["ora_calibration_status"] = CALIBRATION_STATUS
    record: dict[str, Any] = {
        "version": VERSION,
        "run_id": run_id or f"three-test-{uuid.uuid4().hex}",
        "recorded_at": recorded_at or datetime.now(UTC).isoformat(),
        **producer,
        "core_run_id": sources.core.run_id,
        "comparison_run_id": sources.comparison.id,
        "source_run_ids": {name: arm.id for name, arm in sources.arms.items()},
        "source_commits": {
            "core": sources.core.git_commit,
            "comparison": sources.comparison.git_commit,
            **{name: arm.git_commit for name, arm in sources.arms.items()},
        },
        "comparison_manifest_sha256": comparison_hash,
        "strain_rollup": sources.scope,
        "scope_role": "primary" if sources.scope == "phenotype" else "sensitivity",
        "contrast_population": population,
        "direction_basis": direction,
        "analysis_bases": {
            name: basis.model_dump(mode="json") for name, basis in sources.bases.items()
        },
        "retained_traits": list(features),
        "family_size": len(features),
        "primary_ranking": "undecided",
        "ora_declaration": "ora-study-enrichment-direction-v1",
        "ora_calibration_status": CALIBRATION_STATUS,
        "test_families": [
            "full-retained-study-t",
            "secondary-full-retained-study-wilcoxon",
            "ora-study-enrichment-direction-v1",
            "noel-native-full-screen",
            "native-selected-pooled-composition",
        ],
        "inputs": consumed,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = output.parent / f".{output.name}-staging-{uuid.uuid4().hex}"
    require_durable(staging)
    _disjoint(staging, [scratch, *roots])
    # Keep the lexical name: resolving a dangling symlink before mkdir would
    # create its target instead of refusing the existing staging-path collision.
    # Plain mkdir applies the inherited umask without querying or changing it.
    staging.mkdir()
    published = False
    try:
        table.to_csv(staging / "per_trait_evidence.tsv", sep="\t", index=False, lineterminator="\n")
        composition.pair_contrasts.to_csv(staging / "composition_pair_delta_z.tsv", sep="\t")
        composition.study_contrasts.to_csv(staging / "composition_study_delta_z.tsv", sep="\t")
        record["outputs"] = {path.name: _inventory(path) for path in sorted(staging.iterdir())}
        if _inputs(paths, sources) != consumed or _producer(code_commit) != producer:
            raise ValueError("Source or producer changed before provenance publication")
        if multiplicity_inputs is not None:
            _check_multiplicity_inputs(*multiplicity_inputs)
        (staging / MANIFEST).write_text(
            json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
        )
        # The completed bundle becomes visible in one no-replace rename. Any
        # preceding failure removes staging and leaves this destination retryable.
        _publish_directory(staging, output)
        published = True
    finally:
        # Own staging only after exclusive creation and until publication. A
        # replacement at its former name after a successful rename is not ours.
        if not published:
            shutil.rmtree(staging)
    return record


def verify_table(
    output: Path,
    comparison: Path,
    *,
    scratch: Path,
    multiplicity: Path | None = None,
    label_support: Path | None = None,
) -> dict[str, Any]:
    """Verify inventory and regenerate the table from verified original sources."""
    manifest = output / MANIFEST if output.is_dir() else output
    root = manifest.parent
    scratch = require_durable(scratch)
    _disjoint(scratch, _manifest_roots([manifest], multiplicity=None, label_support=None))
    sources = open_verified_sources(comparison)
    # Verified models may use canonical paths; retain the caller's file alias
    # and the comparison's declared aliases before creating any replay scratch.
    comparison_path = comparison / "comparison_run.json" if comparison.is_dir() else comparison
    _disjoint(
        scratch,
        [
            *_declared_source_roots(
                comparison_path, multiplicity=multiplicity, label_support=label_support
            ),
            *_source_roots(sources, multiplicity=multiplicity, label_support=label_support),
        ],
    )
    record: dict[str, Any] = json.loads(manifest.read_text(encoding="utf-8"))
    if record.get("version") != VERSION:
        raise ValueError("Unsupported per-trait evidence version")
    expected = {MANIFEST, *record["outputs"]}
    if {path.name for path in root.iterdir()} != expected:
        raise ValueError("Per-trait output inventory differs")
    for name, item in record["outputs"].items():
        path = root / name
        if Path(name).name != name or path.is_symlink() or _inventory(path) != item:
            raise ValueError("Per-trait output hash differs")
    # Keep the caller's temporary roots unchanged until the nested run has
    # independently checked its output and scratch destinations.
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="three-test-replay-", dir=scratch) as temp:
        replay = run_table(
            comparison,
            Path(temp) / "result",
            scratch=Path(temp) / "scratch",
            code_commit=record["code_commit"],
            multiplicity=multiplicity,
            label_support=label_support,
            recorded_at=record["recorded_at"],
            run_id=record["run_id"],
        )
    if replay != record:
        raise ValueError("Per-trait replay differs from recorded evidence")
    return record
