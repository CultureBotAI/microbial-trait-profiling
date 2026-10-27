"""Read pinned completed three-test bundles without replay or delta-z recomputation.

The caller pins the parent manifest. Its producer is checked against that historical
Git tree, independently of the current consumer commit. Upstream execution checks,
exact artifact identities and persisted-table consistency remain mandatory. This
reader verifies the completed evidence contract; it does not claim a fresh replay.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from git import Blob, Repo
from git.exc import BadName, BadObject, InvalidGitRepositoryError

from mtp.count_input import profile_count_reader
from mtp.multiple_testing import canonical_sha256
from mtp.multiplicity.families import rows
from mtp.multiplicity.sources import (
    VerifiedSources,
    bind_recorded_hashes,
    open_verified_sources,
    recorded_arm_output,
)
from mtp.paths import file_sha256
from mtp.three_test_table.composition import COMPOSITION_CONTRACT
from mtp.three_test_table.signed_rank import (
    LEGACY_WILCOXON_CONTRACT,
    WILCOXON_CONTRACT,
)
from mtp.three_test_table.stage import CALIBRATION_STATUS, MANIFEST, METHODS, VERSION

SUPPORTED_WILCOXON_CONTRACTS = frozenset({LEGACY_WILCOXON_CONTRACT, WILCOXON_CONTRACT})

_OUTPUTS = (
    "composition_pair_delta_z.tsv",
    "composition_study_delta_z.tsv",
    "per_trait_evidence.tsv",
)
_SHA256 = re.compile(r"[0-9a-f]{64}")
_COMMIT = re.compile(r"[0-9a-f]{40}")
_NOEL_FIELDS = {
    "studies_up": "studies_up",
    "studies_down": "studies_down",
    "studies_tied": "studies_tied",
    "studies_scored": "studies_scored",
    "share_up": "share_up",
    "ci_low": "ci_low",
    "ci_high": "ci_high",
    "p_value": "binomial_p",
    "q_by_nominal": "by_q_nominal",
    "q_bh_nominal": "bh_q_nominal",
}


@dataclass(frozen=True)
class CompletedTable:
    """Pinned parent evidence and verified sources for read-only downstream consumers."""

    manifest: Path
    root: Path
    record: dict[str, Any]
    sources: VerifiedSources
    evidence: pd.DataFrame
    pair_delta_z: pd.DataFrame
    study_delta_z: pd.DataFrame
    snapshots: dict[Path, dict[str, Any]] = field(default_factory=dict)
    resolved_targets: dict[Path, Path] = field(default_factory=dict)


def _snapshot(path: Path) -> dict[str, Any]:
    return {"sha256": file_sha256(path), "size_bytes": path.stat().st_size}


def _absolute(path: Path) -> Path:
    """Keep the supplied spelling instead of collapsing a file or directory alias."""
    return path if path.is_absolute() else Path.cwd() / path


def _remember_path(
    path: Path,
    identity: dict[str, Any],
    snapshots: dict[Path, dict[str, Any]],
    targets: dict[Path, Path],
) -> None:
    alias = _absolute(path)
    target = alias.resolve(strict=True)
    if alias in targets and targets[alias] != target:
        raise ValueError(f"Completed-table source alias target changed: {alias.name}")
    if _snapshot(alias) != identity:
        raise ValueError(f"Completed-table input hash changed: {alias.name}")
    for spelling in (alias, target):
        if spelling in snapshots and snapshots[spelling] != identity:
            raise ValueError("Completed-table source aliases have conflicting identities")
        snapshots[spelling] = dict(identity)
        targets[spelling] = target


def _declared_aliases(comparison: dict[str, Any]) -> set[Path]:
    """Read path declarations before verification can canonicalize their spellings."""
    aliases = {
        _absolute(Path(comparison[name]))
        for name in (
            "manifest_path",
            "core_manifest_path",
            "noel_manifest_path",
            "jed_manifest_path",
            "ilyas_manifest_path",
        )
        if isinstance(comparison.get(name), str) and comparison[name].strip()
    }
    core_path = Path(comparison["core_manifest_path"])
    core = json.loads(core_path.read_text(encoding="utf-8"))
    if not isinstance(core, dict):
        raise ValueError("Completed-table declared core is not an object")
    core_root = Path(core.get("manifest_path", str(core_path))).parent
    for artifact in core.get("artifacts") or []:
        if artifact.get("product_id") not in {"count_profiles", "experiments", "feature_metadata"}:
            continue
        for name in ("path", "requested_path"):
            value = artifact.get(name)
            if isinstance(value, str) and value.strip():
                path = Path(value)
                aliases.add(_absolute(path if path.is_absolute() else core_root / path))
    noel_path = Path(comparison["noel_manifest_path"])
    noel = json.loads(noel_path.read_text(encoding="utf-8"))
    if not isinstance(noel, dict):
        raise ValueError("Completed-table declared Noel run is not an object")
    for step in noel.get("steps") or []:
        if step.get("id") != "directional_screen":
            continue
        for output in step.get("outputs") or []:
            path = Path(output["path"])
            if output.get("storage_kind") == "file":
                if path.name == "trait_direction_screen.tsv":
                    aliases.add(_absolute(path))
            else:
                for entry in output.get("inventory") or []:
                    if Path(entry["path"]).name == "trait_direction_screen.tsv":
                        aliases.add(_absolute(path / entry["path"]))
    return aliases


def _verified_product_aliases(sources: VerifiedSources) -> set[Path]:
    aliases: set[Path] = set()
    for artifact in sources.core.artifacts or []:
        if artifact.product_id not in {"count_profiles", "experiments", "feature_metadata"}:
            continue
        for name in ("path", "requested_path"):
            value = getattr(artifact, name, None)
            if isinstance(value, str) and value.strip():
                path = Path(value)
                aliases.add(
                    _absolute(
                        path
                        if path.is_absolute()
                        else Path(sources.core.manifest_path).parent / path
                    )
                )
    return aliases


def _digest(value: object) -> bool:
    return isinstance(value, str) and _SHA256.fullmatch(value) is not None


def _inventory_item(item: object) -> bool:
    return (
        isinstance(item, dict)
        and set(item) == {"sha256", "size_bytes"}
        and _digest(item["sha256"])
        and type(item["size_bytes"]) is int
        and item["size_bytes"] >= 0
    )


def _check_outputs(root: Path, record: dict[str, Any]) -> None:
    outputs = record.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != set(_OUTPUTS):
        raise ValueError("Completed-table output inventory differs")
    if {path.name for path in root.iterdir()} != {MANIFEST, *_OUTPUTS}:
        raise ValueError("Completed-table bundle inventory differs")
    for name in (MANIFEST, *_OUTPUTS):
        path = root / name
        if path.is_symlink() or not path.is_file():
            raise ValueError("Completed-table output must be a regular, non-symlink file")
        if name != MANIFEST and (
            not _inventory_item(outputs[name]) or _snapshot(path) != outputs[name]
        ):
            raise ValueError("Completed-table output hash or size differs")


def _check_producer(record: dict[str, Any]) -> None:
    commit_id = record.get("code_commit")
    if (
        record.get("code_dirty") is not False
        or not isinstance(commit_id, str)
        or _COMMIT.fullmatch(commit_id) is None
    ):
        raise ValueError("Completed-table producer must be a clean pinned Git commit")
    inventory = record.get("code_inventory")
    if not isinstance(inventory, dict) or not all(
        isinstance(name, str) and _digest(digest) for name, digest in inventory.items()
    ):
        raise ValueError("Completed-table producer inventory is malformed")
    try:
        with Repo(Path(__file__).resolve().parent, search_parent_directories=True) as repository:
            commit = repository.commit(commit_id)
            expected = {
                str(item.path): hashlib.sha256(item.data_stream.read()).hexdigest()
                for item in commit.tree.traverse()
                if isinstance(item, Blob)
                and (
                    (str(item.path).startswith("src/mtp/") and str(item.path).endswith(".py"))
                    or str(item.path) in {METHODS, "uv.lock", "pyproject.toml", "Justfile"}
                )
            }
    except (BadName, BadObject, InvalidGitRepositoryError, ValueError) as error:
        raise ValueError("Completed-table producer Git commit is unavailable") from error
    if (
        inventory != expected
        or record.get("implementation_sha256") != canonical_sha256(expected)
        or record.get("methods_sha256") != expected.get(METHODS)
        or record.get("lock_sha256") != expected.get("uv.lock")
    ):
        raise ValueError("Completed-table producer inventory differs from its Git tree")
    environment = record.get("environment")
    if (
        not isinstance(environment, dict)
        or set(environment) != {"python", "pandas", "numpy", "scipy"}
        or not all(isinstance(value, str) and value.strip() for value in environment.values())
    ):
        raise ValueError("Completed-table numerical environment is incomplete")


def _check_contract(record: dict[str, Any], sources: VerifiedSources, comparison_hash: str) -> None:
    traits = list(sources.features["count_profiles"])
    expected = {
        "version": VERSION,
        "core_run_id": sources.core.run_id,
        "comparison_run_id": sources.comparison.id,
        "comparison_manifest_sha256": comparison_hash,
        "source_run_ids": {name: arm.id for name, arm in sources.arms.items()},
        "source_commits": {
            "core": sources.core.git_commit,
            "comparison": sources.comparison.git_commit,
            **{name: arm.git_commit for name, arm in sources.arms.items()},
        },
        "analysis_bases": {
            name: basis.model_dump(mode="json") for name, basis in sources.bases.items()
        },
        "strain_rollup": sources.scope,
        "scope_role": "primary" if sources.scope == "phenotype" else "sensitivity",
        "contrast_population": sources.comparison.contrast_population,
        "direction_basis": sources.comparison.direction_basis,
        "retained_traits": traits,
        "family_size": len(traits),
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
    }
    if (
        not traits
        or len(set(traits)) != len(traits)
        or tuple(traits) != sources.features["ora_profiles"]
        or type(record.get("family_size")) is not int
        or not isinstance(record.get("run_id"), str)
        or not record["run_id"].strip()
        or not isinstance(record.get("recorded_at"), str)
        or not record["recorded_at"].strip()
    ):
        raise ValueError("Completed-table contract has an invalid catalog or identity")
    for name, value in expected.items():
        if record.get(name) != value:
            raise ValueError(f"Completed-table source contract identity differs: {name}")


def _bound_snapshot(path: Path, sources: VerifiedSources, record: dict[str, Any]) -> dict[str, Any]:
    snapshot = _snapshot(path)
    if snapshot["sha256"] != sources.expected_digest(path):
        raise ValueError(f"Completed-table source differs from verified input: {path.name}")
    identity = f"{path.name}@sha256:{snapshot['sha256']}"
    if record["inputs"].get(identity) != snapshot:
        raise ValueError(f"Completed-table source input is not bound: {path.name}")
    return snapshot


def _read_table(path: Path, keys: list[str]) -> pd.DataFrame:
    with path.open(encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    if not header or len(set(header)) != len(header) or any(not name for name in header):
        raise ValueError("Completed-table TSV header has missing or duplicate columns")
    if header[: len(keys)] != keys:
        raise ValueError("Completed-table ledger identity columns differ")
    frame: pd.DataFrame = pd.read_csv(
        path,
        sep="\t",
        dtype=dict.fromkeys(keys, str),
        keep_default_na=False,
        float_precision="round_trip",
    )
    if (
        frame.empty
        or frame.duplicated(keys).any()
        or any(
            not isinstance(value, str) or not value.strip() for key in keys for value in frame[key]
        )
    ):
        raise ValueError("Completed-table ledger identity is empty or duplicate")
    return frame


def _ledger(path: Path, keys: list[str], traits: list[str]) -> pd.DataFrame:
    frame: pd.DataFrame = _read_table(path, keys)
    if list(frame.columns) != [*keys, *traits]:
        raise ValueError("Completed-table ledger columns differ from the retained catalog")
    frame = frame.set_index(keys)
    try:
        frame = frame.apply(pd.to_numeric, errors="raise")
        values = frame.to_numpy(dtype=float)
    except (ValueError, TypeError) as error:
        raise ValueError("Completed-table ledger requires finite numerical contrasts") from error
    if not np.isfinite(values).all():
        raise ValueError("Completed-table ledger requires finite numerical contrasts")
    return frame


def _pair_population(sources: VerifiedSources) -> tuple[set[tuple[str, str]], int]:
    metadata: dict[tuple[str, str], str] = {}
    for row in rows(
        sources.product("experiments"), sources.expected_digest(sources.product("experiments"))
    ):
        key = row["study_id"], row["experiment_id"]
        if key in metadata or not all(value.strip() for value in key):
            raise ValueError("Completed-table source experiment identity is missing or duplicate")
        metadata[key] = row["host_species"]
    directions: dict[tuple[str, str], set[str]] = {}
    profiles = 0
    for profile_row in profile_count_reader(sources.product("count_profiles")):
        study, experiment, direction = profile_row.key
        key = study, experiment
        if key not in metadata or direction not in {"inc", "dec"}:
            raise ValueError("Completed-table source profile population is invalid")
        if metadata[key].strip() == "Homo sapiens":
            directions.setdefault(key, set()).add(direction)
            profiles += 1
    return {key for key, directions in directions.items() if directions == {"inc", "dec"}}, profiles


def _same_number(actual: Any, expected: Any, description: str) -> None:
    if actual == "" or expected == "":
        if actual != expected:
            raise ValueError(f"Completed-table {description} differs")
        return
    try:
        first, second = float(actual), float(expected)
    except (ValueError, TypeError) as error:
        raise ValueError(f"Completed-table {description} is not numerical") from error
    if (
        not math.isfinite(first)
        or not math.isfinite(second)
        or not math.isclose(first, second, rel_tol=1e-12, abs_tol=1e-14)
    ):
        raise ValueError(f"Completed-table {description} differs")


def _same_noel(actual: Any, expected: str, field_name: str) -> None:
    if field_name.startswith("studies_"):
        try:
            count = float(actual)
            reference = int(expected)
        except (TypeError, ValueError) as error:
            raise ValueError("Completed-table Noel count is not an integer") from error
        if not math.isfinite(count) or not count.is_integer() or count < 0 or count != reference:
            raise ValueError("Completed-table Noel count differs from the native count")
    else:
        _same_number(actual, expected, f"Noel {field_name}")
        if actual != "" and not 0 <= float(actual) <= 1:
            raise ValueError("Completed-table Noel probability is outside zero to one")


def _check_summary(
    evidence: pd.DataFrame,
    pairs: pd.DataFrame,
    studies: pd.DataFrame,
    sources: VerifiedSources,
    input_profiles: int,
) -> None:
    traits = list(sources.features["count_profiles"])
    if evidence.trait.tolist() != traits:
        raise ValueError("Completed-table summary catalog differs")
    labels: dict[str, str] = {}
    for row in rows(
        sources.product("feature_metadata"),
        sources.expected_digest(sources.product("feature_metadata")),
    ):
        trait = row["column_id"]
        if trait in labels:
            raise ValueError("Completed-table source label identity is duplicate")
        labels[trait] = row["object_name"] or row["column_label"] or trait
    screen_path = recorded_arm_output(
        sources.arms["noel"], "directional_screen", "trait_direction_screen.tsv"
    )
    screen_rows = list(rows(screen_path, sources.expected_digest(screen_path)))
    screen = {row["trait"]: row for row in screen_rows}
    if len(screen) != len(screen_rows) or set(screen) != set(traits):
        raise ValueError("Completed-table native Noel catalog differs")
    contracts = {
        "composition_contract": COMPOSITION_CONTRACT,
        "composition_family_size": len(traits),
        "composition_input_profiles": input_profiles,
        "composition_profiles": 2 * len(pairs),
        "composition_pairs": len(pairs),
        "composition_studies": len(studies),
        "composition_unpaired_profiles_excluded": input_profiles - 2 * len(pairs),
        "composition_effect_units": "standardized_log1p_1000_retained_fraction_difference",
        "composition_test_unit": "study",
        "composition_zero_total_policy": "retained_as_zero_fraction",
        "study_t_primary_correction": "BY",
        "noel_family_id": "noel-full-screen",
        "noel_family_size": len(traits),
        "strain_rollup": sources.scope,
        "contrast_population": sources.comparison.contrast_population,
        "direction_basis": sources.comparison.direction_basis,
        "core_run_id": sources.core.run_id,
        "primary_ranking": "undecided",
        "ora_calibration_status": CALIBRATION_STATUS,
    }
    required = {
        "trait_label",
        "experiment_mean_delta_z",
        "study_mean_delta_z",
        "study_wilcoxon_contract",
        "noel_status",
        *contracts,
        *(f"noel_{name}" for name in _NOEL_FIELDS),
    }
    if not required.issubset(evidence.columns):
        raise ValueError("Completed-table summary contract columns are missing")
    wilcoxon_contracts = set(evidence["study_wilcoxon_contract"])
    if len(wilcoxon_contracts) != 1 or not wilcoxon_contracts <= SUPPORTED_WILCOXON_CONTRACTS:
        raise ValueError("Completed-table summary contract differs: study_wilcoxon_contract")
    for summary in evidence.to_dict(orient="records"):
        trait = summary["trait"]
        if summary["trait_label"] != labels.get(trait):
            raise ValueError("Completed-table trait label differs from verified metadata")
        for name, expected in contracts.items():
            if summary[name] != expected:
                raise ValueError(f"Completed-table summary contract differs: {name}")
        _same_number(summary["experiment_mean_delta_z"], pairs[trait].mean(), "pair mean")
        _same_number(summary["study_mean_delta_z"], studies[trait].mean(), "study mean")
        if summary["noel_status"] != screen[trait]["test_status"]:
            raise ValueError("Completed-table native Noel status differs")
        for output, native in _NOEL_FIELDS.items():
            _same_noel(summary[f"noel_{output}"], screen[trait][native], output)
    for effect in ("experiment_mean_delta_z", "study_mean_delta_z"):
        values = pd.to_numeric(evidence[effect], errors="raise")
        for name, ranked, ascending in (
            ("ascending", values, True),
            ("descending", values, False),
            ("absolute_descending", values.abs(), False),
        ):
            column = f"{effect}_rank_{name}"
            if column not in evidence or not np.array_equal(
                pd.to_numeric(evidence[column], errors="raise").to_numpy(),
                ranked.rank(method="min", ascending=ascending).to_numpy(),
            ):
                raise ValueError("Completed-table summary ranks differ from exported means")


def recheck_completed_table(bundle: CompletedTable) -> None:
    """Check the same parent bytes and bound sources again before publishing a consumer."""
    for alias, target in bundle.resolved_targets.items():
        if alias.resolve() != target:
            raise ValueError(f"Completed-table source alias target changed: {alias.name}")
    _check_outputs(bundle.root, bundle.record)
    for path, expected in bundle.snapshots.items():
        if not path.is_file() or _snapshot(path) != expected:
            raise ValueError(f"Completed-table input hash changed: {path.name}")


def load_completed_table(
    manifest: Path,
    comparison: Path,
    *,
    expected_sha256: str,
) -> CompletedTable:
    """Read one pinned completed bundle while verifying its original source bindings."""
    manifest = _absolute(manifest / MANIFEST if manifest.is_dir() else manifest)
    supplied_manifest = manifest
    resolved_targets = {manifest: manifest.resolve(strict=True)}
    raw = manifest.read_bytes()
    if not _digest(expected_sha256) or hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise ValueError("Completed-table manifest SHA256 pin differs")
    manifest = manifest.resolve()
    record = json.loads(raw)
    if not isinstance(record, dict) or record.get("version") != VERSION:
        raise ValueError("Unsupported completed-table contract version")
    _check_outputs(manifest.parent, record)
    _check_producer(record)
    inputs = record.get("inputs")
    if (
        not isinstance(inputs, dict)
        or not inputs
        or any(
            not isinstance(identity, str)
            or not _inventory_item(item)
            or identity != f"{identity.split('@sha256:', 1)[0]}@sha256:{item['sha256']}"
            for identity, item in inputs.items()
        )
    ):
        raise ValueError("Completed-table input inventory is malformed")
    comparison = _absolute(
        comparison / "comparison_run.json" if comparison.is_dir() else comparison
    )
    resolved_targets[comparison] = comparison.resolve(strict=True)
    comparison_hash = file_sha256(comparison)
    declared_comparison = json.loads(comparison.read_text(encoding="utf-8"))
    if not isinstance(declared_comparison, dict):
        raise ValueError("Completed-table declared comparison is not an object")
    declared_aliases = _declared_aliases(declared_comparison)
    for path in declared_aliases:
        target = path.resolve(strict=True)
        if path in resolved_targets and resolved_targets[path] != target:
            raise ValueError(f"Completed-table source alias target changed: {path.name}")
        resolved_targets[path] = target
    sources = bind_recorded_hashes(open_verified_sources(comparison), comparison_hash)
    _check_contract(record, sources, comparison_hash)
    snapshots: dict[Path, dict[str, Any]] = {}
    _remember_path(
        supplied_manifest,
        {"sha256": expected_sha256, "size_bytes": len(raw)},
        snapshots,
        resolved_targets,
    )
    for name, item in record["outputs"].items():
        _remember_path(manifest.parent / name, item, snapshots, resolved_targets)
    source_paths = {
        comparison,
        Path(sources.core.manifest_path),
        *(Path(arm.manifest_path) for arm in sources.arms.values()),
        *(sources.product(name) for name in ("count_profiles", "experiments", "feature_metadata")),
        recorded_arm_output(
            sources.arms["noel"], "directional_screen", "trait_direction_screen.tsv"
        ),
    }
    for path in source_paths:
        _remember_path(path, _bound_snapshot(path, sources, record), snapshots, resolved_targets)
    for alias in declared_aliases | _verified_product_aliases(sources):
        target = alias.resolve(strict=True)
        if target not in snapshots:
            raise ValueError("Declared source alias is not bound to a consumed input")
        _remember_path(alias, snapshots[target], snapshots, resolved_targets)
    traits = list(sources.features["count_profiles"])
    pairs = _ledger(manifest.parent / _OUTPUTS[0], ["study_id", "experiment_id"], traits)
    studies = _ledger(manifest.parent / _OUTPUTS[1], ["study_id"], traits)
    expected_pairs, input_profiles = _pair_population(sources)
    if set(pairs.index) != expected_pairs or set(studies.index) != {
        key[0] for key in expected_pairs
    }:
        raise ValueError("Completed-table ledger population differs from the verified core")
    means = pairs.groupby(level="study_id", sort=False).mean().reindex(studies.index)
    if not np.allclose(
        studies.to_numpy(dtype=float), means.to_numpy(dtype=float), rtol=1e-12, atol=1e-14
    ):
        raise ValueError("Completed-table study means differ from its pair ledger")
    evidence = _read_table(manifest.parent / _OUTPUTS[2], ["trait"])
    _check_summary(evidence, pairs, studies, sources, input_profiles)
    bundle = CompletedTable(
        manifest,
        manifest.parent,
        record,
        sources,
        evidence,
        pairs,
        studies,
        snapshots,
        resolved_targets,
    )
    recheck_completed_table(bundle)
    return bundle


__all__ = ["CompletedTable", "load_completed_table", "recheck_completed_table"]
