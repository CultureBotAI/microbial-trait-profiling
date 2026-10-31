"""Read pinned completed figure inputs without rerunning analysis (#1046).

The parent reader verifies historical upstream execution. This reader binds the
completed derived bundle to that parent, checks descriptive identities used by
figures, and preserves original inference fields and metadata labels.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from git import Blob, Repo
from git.exc import BadName, BadObject, InvalidGitRepositoryError

from mtp.multiple_testing import canonical_sha256
from mtp.rollup_comparison.validation import validate_evidence
from mtp.three_test_table.reader import (
    CompletedTable,
    load_completed_table,
    recheck_completed_table,
)
from mtp.three_test_table.source_paths import _spellings

MANIFEST = "delta_z_derived_run.json"
METHODS = "methods/delta_z_derived.md"
PAIR = "composition_pair_delta_z.tsv"
STUDY = "composition_study_delta_z.tsv"
OUTPUTS = (
    PAIR,
    STUDY,
    "trait_summary.tsv",
    "quadrant_pairs.tsv",
    "quadrant_summary.tsv",
    "stool_disease_shifts.tsv",
    "stool_study_sensitivity.tsv",
    "site_disease_shifts.tsv",
)
FA, OA = "METPO:1000608", "METPO:1000607"
KEYS = ["study_id", "experiment_id"]
_SHA = re.compile(r"[0-9a-f]{64}")


@dataclass(frozen=True)
class FigureInputs:
    """Verified parent, original derived tables and all publication recheck identities."""

    parent: CompletedTable
    record: dict[str, Any]
    tables: dict[str, pd.DataFrame]
    labels: dict[str, str]
    snapshots: dict[Path, dict[str, Any]]
    targets: dict[Path, Path]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _identity(raw: bytes) -> dict[str, Any]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def _inventory_item(value: Any) -> bool:
    return (
        isinstance(value, dict)
        and set(value) == {"sha256", "size_bytes"}
        and isinstance(value["sha256"], str)
        and _SHA.fullmatch(value["sha256"]) is not None
        and type(value["size_bytes"]) is int
        and value["size_bytes"] >= 0
    )


def _read(
    path: Path,
    snapshots: dict[Path, dict[str, Any]],
    targets: dict[Path, Path],
    *,
    expected: dict[str, Any] | None = None,
    digest: str | None = None,
) -> bytes:
    """Bind aliases before reading, then hash and parse this one byte buffer."""
    spellings = _spellings(path)
    resolved = {alias: alias.resolve(strict=True) for alias in spellings}
    _require(len(set(resolved.values())) == 1, "Input aliases resolve inconsistently")
    _require(all(alias.is_file() for alias in spellings), "Input must be a regular file")
    raw = path.read_bytes()
    actual = _identity(raw)
    if digest is not None:
        _require(
            _SHA.fullmatch(digest) is not None and actual["sha256"] == digest,
            "Figure input manifest SHA256 pin differs",
        )
    if expected is not None:
        _require(
            _inventory_item(expected) and actual == expected,
            f"Figure input hash or size differs: {path.name}",
        )
    for alias, target in resolved.items():
        _require(alias.resolve(strict=True) == target, "Input alias target changed while reading")
        _require(alias not in targets or targets[alias] == target, "Input alias target changed")
        _require(alias not in snapshots or snapshots[alias] == actual, "Input identity changed")
        snapshots[alias], targets[alias] = actual.copy(), target
    return raw


def _json(raw: bytes) -> dict[str, Any]:
    def unique(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require(key not in result, f"Duplicate manifest field: {key}")
            result[key] = value
        return result

    value = json.loads(raw, object_pairs_hook=unique)
    _require(isinstance(value, dict), "Figure input manifest must be an object")
    result: dict[str, Any] = value
    return result


def _table(raw: bytes) -> pd.DataFrame:
    reader = csv.reader(io.StringIO(raw.decode("utf-8")), delimiter="\t")
    header = next(reader, [])
    _require(
        bool(header) and all(header) and len(header) == len(set(header)),
        "Missing or duplicate figure input columns",
    )
    _require(all(len(row) == len(header) for row in reader), "Malformed figure input row width")
    identifiers = {
        "study_id",
        "experiment_id",
        "trait",
        "trait_label",
        "body_site",
        "condition",
        "fa_trait",
        "oa_trait",
        "fa_label",
        "oa_label",
        "population",
        "column_id",
        "object_name",
        "column_label",
    }
    frame: pd.DataFrame = pd.read_csv(
        io.BytesIO(raw),
        sep="\t",
        keep_default_na=False,
        low_memory=False,
        dtype={name: str for name in header if name in identifiers},
        float_precision="round_trip",
    )
    return frame


def _producer(record: dict[str, Any]) -> None:
    commit_id = record.get("code_commit")
    _require(
        record.get("code_dirty") is False
        and isinstance(commit_id, str)
        and re.fullmatch(r"[0-9a-f]{40}", commit_id) is not None,
        "Derived producer must be a clean pinned Git commit",
    )
    try:
        with Repo(Path(__file__).resolve().parent, search_parent_directories=True) as repository:
            commit = repository.commit(commit_id)
            expected = {
                str(item.path): hashlib.sha256(item.data_stream.read()).hexdigest()
                for item in commit.tree.traverse()
                if isinstance(item, Blob)
                and (
                    (str(item.path).startswith("src/mtp/") and str(item.path).endswith(".py"))
                    or str(item.path)
                    in {
                        METHODS,
                        "methods/three_test_table.md",
                        "uv.lock",
                        "pyproject.toml",
                        "Justfile",
                    }
                )
            }
    except (BadName, BadObject, InvalidGitRepositoryError, ValueError) as error:
        raise ValueError("Derived producer Git commit is unavailable") from error
    _require(
        record.get("code_inventory") == expected
        and record.get("implementation_sha256") == canonical_sha256(expected)
        and record.get("methods_sha256") == expected.get(METHODS)
        and record.get("lock_sha256") == expected.get("uv.lock"),
        "Derived producer inventory differs from its historical Git tree",
    )
    environment = record.get("environment")
    _require(
        isinstance(environment, dict)
        and set(environment) == {"python", "numpy", "pandas", "scipy"}
        and all(isinstance(value, str) and value.strip() for value in environment.values()),
        "Derived producer numerical environment is incomplete",
    )


def _outputs(root: Path, record: dict[str, Any]) -> None:
    outputs = record.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != set(OUTPUTS):
        raise ValueError("Derived output inventory must contain the full eight tables")
    _require(
        {path.name for path in root.iterdir()} == {MANIFEST, *OUTPUTS},
        "Derived bundle inventory differs",
    )
    for name in (MANIFEST, *OUTPUTS):
        path = root / name
        _require(
            path.is_file() and not path.is_symlink(),
            "Derived bundle requires regular nonsymlink files",
        )
        if name != MANIFEST:
            _require(_inventory_item(outputs[name]), "Malformed derived output identity")


def _parent_binding(record: dict[str, Any], parent: CompletedTable, pin: str) -> None:
    expected = {
        name: parent.record[name]
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
    }
    expected.update(
        parent_run_id=parent.record["run_id"],
        parent_manifest_sha256=pin,
        parent_producer={name: parent.record[name] for name in ("code_commit", "code_dirty")},
        parent_input_inventory=parent.record["inputs"],
    )
    for name, value in expected.items():
        _require(record.get(name) == value, f"Derived parent identity differs: {name}")
    consumed = record.get("inputs")
    if not isinstance(consumed, dict) or not consumed:
        raise ValueError("Derived consumed input inventory is missing")
    for key, item in consumed.items():
        _require(
            isinstance(key, str)
            and _inventory_item(item)
            and key == f"{key.split('@sha256:', 1)[0]}@sha256:{item['sha256']}",
            "Malformed derived consumed input identity",
        )
    for item in parent.snapshots.values():
        _require(
            item in consumed.values(), "Derived input inventory does not bind the verified parent"
        )
    for name in (PAIR, STUDY):
        _require(
            record["outputs"][name] == parent.record["outputs"][name],
            f"Copied parent ledger identity differs: {name}",
        )


def _indexed(frame: pd.DataFrame, keys: list[str], description: str) -> pd.DataFrame:
    _require(set(keys).issubset(frame.columns), f"{description} identity columns are missing")
    _require(
        not frame.empty
        and not frame.duplicated(keys).any()
        and all(
            isinstance(value, str) and bool(value.strip()) for key in keys for value in frame[key]
        ),
        f"{description} identity keys are empty or duplicate",
    )
    return frame.set_index(keys)


def _number(value: Any) -> float:
    _require(not isinstance(value, bool), "Boolean is not a figure statistic")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError("Figure statistic is not numeric") from error
    _require(math.isfinite(result), "Figure statistic must be finite")
    return result


def _ora_evidence(parent: CompletedTable) -> None:
    """Check inherited ORA support/status arithmetic without recomputing inference.

    Directional counts are always present in the native contract, including
    zero-support rows. Nullable enrichment/depletion counts and p/q fields are
    retained unchanged; they are not inputs to these descriptive checks.
    """
    names = (
        "studies_up",
        "studies_down",
        "studies_tied",
        "studies_scored",
        "studies_eligible",
        "studies_excluded",
        "studies_total",
        "pairs_eligible",
        "pairs_excluded",
        "pairs_total",
    )
    columns = {name: f"ora_direction_{name}" for name in names}
    columns["background_k"] = "ora_background_k"
    _require(
        {*columns.values(), "ora_direction_status"}.issubset(parent.evidence.columns),
        "ORA directional support or status fields are missing",
    )
    for row in parent.evidence.to_dict(orient="records"):
        counts: dict[str, int] = {}
        for name, column in columns.items():
            value = row[column]
            _require(not isinstance(value, (bool, np.bool_)), "ORA count cannot be boolean")
            try:
                numeric = _number(value)
            except ValueError as error:
                raise ValueError(f"ORA count must be a nonnegative integer: {column}") from error
            _require(
                numeric >= 0 and numeric.is_integer(),
                f"ORA count must be a nonnegative integer: {column}",
            )
            counts[name] = int(numeric)
        _require(
            counts["studies_up"] + counts["studies_down"] == counts["studies_scored"]
            and counts["studies_scored"] + counts["studies_tied"] == counts["studies_eligible"]
            and counts["studies_eligible"] + counts["studies_excluded"]
            == counts["studies_total"]
            == len(parent.study_delta_z)
            and counts["pairs_eligible"] + counts["pairs_excluded"]
            == counts["pairs_total"]
            == len(parent.pair_delta_z)
            and counts["studies_eligible"] <= counts["pairs_eligible"],
            "ORA directional support counts are inconsistent",
        )
        if not counts["background_k"]:
            expected_status = "no_background_taxa"
            _require(counts["pairs_eligible"] == 0, "ORA absent background has eligible pairs")
        elif not counts["pairs_eligible"]:
            expected_status = "no_eligible_pairs"
        else:
            _require(counts["studies_eligible"] > 0, "ORA eligible pairs have no eligible studies")
            expected_status = "tested" if counts["studies_scored"] else "all_tied"
        _require(
            row["ora_direction_status"] == expected_status,
            "ORA directional status differs from its support counts",
        )


def _labels(frame: pd.DataFrame, parent: CompletedTable) -> dict[str, str]:
    _require(
        {"column_id", "object_name", "column_label"}.issubset(frame.columns),
        "Verified feature metadata lacks label columns",
    )
    rows = _indexed(frame, ["column_id"], "Feature metadata")
    labels: dict[str, str] = {}
    for trait in parent.record["retained_traits"]:
        _require(trait in rows.index, "Retained trait lacks verified metadata")
        row = rows.loc[trait]
        label = row["object_name"] or row["column_label"]
        _require(isinstance(label, str) and bool(label.strip()), "Verified trait label is empty")
        labels[trait] = label
    for row in parent.evidence.to_dict(orient="records"):
        _require(
            row["trait_label"] == labels[row["trait"]],
            "Parent label differs from verified metadata",
        )
    return labels


def _summary(frame: pd.DataFrame, parent: CompletedTable) -> int:
    indexed = _indexed(frame, ["trait"], "Trait summary")
    _require(
        frame.trait.tolist() == parent.evidence.trait.tolist(), "Trait summary catalog differs"
    )
    _require(
        set(parent.evidence.columns).issubset(frame.columns), "Trait summary lost parent evidence"
    )
    original = parent.evidence.reset_index(drop=True)
    try:
        pd.testing.assert_frame_equal(
            frame[list(original.columns)].reset_index(drop=True),
            original,
            check_dtype=False,
            check_exact=True,
        )
    except AssertionError as error:
        raise ValueError("Derived trait summary differs from parent evidence") from error
    active_count = 0
    for trait in parent.record["retained_traits"]:
        values = parent.study_delta_z[trait]
        up, down = int(values.gt(1e-12).sum()), int(values.lt(-1e-12).sum())
        expected = {
            "sdz_up": up,
            "sdz_down": down,
            "sdz_active_studies": up + down,
            "sdz_tied": len(values) - up - down,
            "sdz_studies_total": len(values),
            "sdz_active_threshold": 1e-12,
            "sdz_active_at_least_50": up + down >= 50,
            "sdz_quantity": "sign_of_study_mean_delta_z",
            "noel_quantity": "native_raw_fraction_direction",
        }
        row = indexed.loc[trait]
        for name, value in expected.items():
            _require(name in row and row[name] == value, f"Trait summary activity differs: {name}")
        active_count += int(up + down >= 50)
    return active_count


def _quadrants(
    tables: dict[str, pd.DataFrame],
    parent: CompletedTable,
    experiments: pd.DataFrame,
    labels: dict[str, str],
) -> None:
    pairs = _indexed(tables["quadrant_pairs.tsv"], KEYS, "Quadrant pair")
    reference = parent.pair_delta_z
    _require(
        set(pairs.index) == set(reference.index), "Quadrant pair population differs from parent"
    )
    pairs = pairs.reindex(reference.index)
    metadata = _indexed(experiments, KEYS, "Experiments")
    _require(set(reference.index).issubset(metadata.index), "Quadrant pair metadata is missing")
    metadata = metadata.reindex(reference.index)
    required = {
        "fa_trait",
        "oa_trait",
        "fa_delta_z",
        "oa_delta_z",
        "quadrant",
        "is_stool",
        "body_site",
        "condition",
    }
    _require(required.issubset(pairs.columns), "Quadrant pair columns are missing")
    _require(
        {"body_site", "condition"}.issubset(metadata.columns),
        "Experiment metadata columns are missing",
    )
    for column, trait in (("fa", FA), ("oa", OA)):
        _require(
            trait in reference and bool(pairs[f"{column}_trait"].eq(trait).all()),
            "Quadrant trait identity differs",
        )
        values = pairs[f"{column}_delta_z"].map(_number).to_numpy()
        _require(
            np.array_equal(values, reference[trait].to_numpy(dtype=float)),
            "Quadrant values differ from parent pairs",
        )
    for name in ("body_site", "condition"):
        _require(
            pairs[name].equals(metadata[name]), "Quadrant pair metadata differs from verified core"
        )
    stool = metadata.body_site.eq("Feces")
    _require(
        all(isinstance(value, (bool, np.bool_)) for value in pairs.is_stool)
        and np.array_equal(pairs.is_stool.to_numpy(), stool.to_numpy()),
        "Quadrant stool membership differs",
    )
    x, y = pairs.fa_delta_z, pairs.oa_delta_z
    expected = np.select(
        [(x > 0) & (y < 0), (x < 0) & (y > 0), (x > 0) & (y > 0), (x < 0) & (y < 0)],
        ["fa_up_oa_down", "fa_down_oa_up", "both_up", "both_down"],
        default="axis",
    )
    _require(np.array_equal(pairs.quadrant.to_numpy(), expected), "Quadrant signs or axes differ")
    summaries = _indexed(tables["quadrant_summary.tsv"], ["population"], "Quadrant summary")
    _require(set(summaries.index) == {"all", "stool"}, "Quadrant summary population differs")
    for population, subset in (("all", pairs), ("stool", pairs.loc[stool])):
        row = summaries.loc[population]
        counts = subset.quadrant.value_counts()
        expected_counts: dict[str, Any] = {
            name: int(counts.get(name, 0))
            for name in ("fa_up_oa_down", "fa_down_oa_up", "both_up", "both_down", "axis")
        }
        expected_counts.update(
            n_comparisons=len(subset),
            n_studies=subset.index.get_level_values("study_id").nunique(),
            sign_test_n=expected_counts["fa_up_oa_down"] + expected_counts["fa_down_oa_up"],
            fa_trait=FA,
            oa_trait=OA,
            fa_label=labels[FA],
            oa_label=labels[OA],
            sign_threshold=0.0,
            test_unit="comparison",
            p_adjustment="none_nominal",
            stool_field="experiments.body_site",
            stool_value="Feces",
        )
        for name, value in expected_counts.items():
            _require(name in row and row[name] == value, f"Quadrant summary differs: {name}")
        for status, fields in (
            ("sign_test_status", ("sign_test_p",)),
            ("pearson_status", ("pearson_r", "pearson_p")),
        ):
            _require(status in row, "Quadrant inference status is missing")
            if row[status] == "tested":
                for name in fields:
                    value = _number(row[name])
                    _require(
                        (-1 if name == "pearson_r" else 0) <= value <= 1,
                        "Quadrant inference value is outside its range",
                    )


def load_inputs(
    three_test: Path,
    comparison: Path,
    derived: Path,
    *,
    three_test_sha256: str,
    derived_sha256: str,
) -> FigureInputs:
    """Load completed inputs read-only; never call an analysis producer or replay."""
    manifest = derived / MANIFEST if derived.is_dir() else derived
    snapshots: dict[Path, dict[str, Any]] = {}
    targets: dict[Path, Path] = {}
    record = _json(_read(manifest, snapshots, targets, digest=derived_sha256))
    _require(
        record.get("version") == "delta-z-derived-v1" and record.get("summary_only") is False,
        "Figures require the supported full derived bundle",
    )
    _require(
        isinstance(record.get("run_id"), str) and bool(record["run_id"].strip()),
        "Derived run ID is missing",
    )
    root = manifest.resolve(strict=True).parent
    _outputs(root, record)
    _producer(record)
    parent = load_completed_table(three_test, comparison, expected_sha256=three_test_sha256)
    validate_evidence(parent)
    _ora_evidence(parent)
    _parent_binding(record, parent, three_test_sha256)
    snapshots.update(parent.snapshots)
    targets.update(parent.resolved_targets)
    tables: dict[str, pd.DataFrame] = {}
    for name in OUTPUTS:
        path = manifest.parent / name if manifest.parent.resolve() == root else root / name
        tables[name] = _table(_read(path, snapshots, targets, expected=record["outputs"][name]))
    source_tables = {}
    for name in ("feature_metadata", "experiments"):
        path = parent.sources.product(name)
        expected = parent.snapshots.get(path.resolve())
        _require(expected is not None, f"Parent does not bind {name}")
        source_tables[name] = _table(_read(path, snapshots, targets, expected=expected))
    labels = _labels(source_tables["feature_metadata"], parent)
    active = _summary(tables["trait_summary.tsv"], parent)
    _quadrants(tables, parent, source_tables["experiments"], labels)
    _require(
        record.get("counts")
        == {
            "traits": len(parent.evidence),
            "comparisons": len(parent.pair_delta_z),
            "studies": len(parent.study_delta_z),
            "traits_with_at_least_50_active_studies": active,
        },
        "Derived count inventory differs from persisted parent",
    )
    policies = record.get("policies", {})
    for name, value in {
        "row2_x_trait": FA,
        "row2_y_trait": OA,
        "row2_zero_policy": "exact_zero",
        "stool_metadata_field": "body_site",
        "stool_metadata_value": "Feces",
        "row1_activity_epsilon": 1e-12,
        "row1_minimum_active_studies": 50,
    }.items():
        _require(policies.get(name) == value, f"Derived figure policy differs: {name}")
    result = FigureInputs(parent, record, tables, labels, snapshots, targets)
    recheck_inputs(result)
    return result


def recheck_inputs(bundle: FigureInputs) -> None:
    """Refuse changed input bytes, bundle membership or alias targets before publication."""
    recheck_completed_table(bundle.parent)
    for alias, target in bundle.targets.items():
        _require(alias.resolve(strict=True) == target, "Figure input alias target changed")
    for path, expected in bundle.snapshots.items():
        _require(
            path.is_file() and _identity(path.read_bytes()) == expected,
            f"Figure input hash changed: {path.name}",
        )
    roots = {target.parent for path, target in bundle.targets.items() if path.name == MANIFEST}
    _require(len(roots) == 1, "Derived manifest identity is ambiguous")
    _outputs(roots.pop(), bundle.record)
