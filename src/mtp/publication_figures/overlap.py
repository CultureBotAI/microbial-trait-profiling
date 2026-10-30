"""Read and display a complete declared carrier-overlap table (#1046).

The native consumer verifier authenticates upstream execution. Figure validation
checks saved descriptive quantities and identities; it never recomputes the
matrix or selects pairs by their overlap.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections.abc import Iterator
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

import matplotlib
import pandas as pd
from git import Blob, Repo
from git.exc import BadName, BadObject, InvalidGitRepositoryError
from matplotlib import colormaps
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

from mtp.datamodel.workflow_pydantic import ConsumerRecord, WorkflowRun
from mtp.publication_figures.inputs import _json, _read, _require, _table
from mtp.publication_figures.plots import STYLE, Page, _figure, _wrap, validate_authored_wording
from mtp.trait_lists import parse_trait_list, trait_labels
from mtp.workflow.consumers import verify_consumer_record
from mtp.workflow.trait_policy import resolve_trait_policy

MANIFEST = "core_consumption.json"
TABLE = "matrix_trait_overlap.tsv"
METHOD = "methods/matrix_coverage_and_pd_overlap.md"
SCRIPTS = (
    "src/mtp/matrix_overlap_coverage.py",
    "src/mtp/trait_lists.py",
    "src/mtp/matrix_overlap_inputs.py",
    "src/mtp/human_host_corpus.py",
    "src/mtp/workflow/coverage.py",
    "src/mtp/workflow/output_roots.py",
)
OUTPUTS = (
    TABLE,
    "tier12_coverage_by_rank.tsv",
    "filter_cascade.tsv",
    "parkinsons_stool_study_overlap_summary.tsv",
    "parkinsons_stool_study_pairs.tsv",
    "summary.json",
    "provenance.json",
    "core_consumption.coverage.tsv",
)
COLUMNS = (
    "trait_a_id",
    "trait_a_label",
    "trait_b_id",
    "trait_b_label",
    "matrix_taxa_total",
    "trait_a_taxa",
    "trait_b_taxa",
    "intersection_taxa",
    "union_taxa",
    "jaccard_fraction",
    "jaccard_pct",
)
POPULATION = "disease_vs_healthy_tier1_2_subject_controls"
# Native tables use %.12g. These tolerances cover decimal serialization only.
RATIO_REL_TOLERANCE = 1e-11
RATIO_ABS_TOLERANCE = 1e-12


@dataclass(frozen=True)
class OverlapInputs:
    """Pinned original records, verified labels and unchanged saved pair rows."""

    record: dict[str, Any]
    summary: dict[str, Any]
    consumer: ConsumerRecord
    context: dict[str, str]
    trait_ids: tuple[str, ...]
    labels: dict[str, str]
    snapshots: dict[Path, dict[str, Any]]
    targets: dict[Path, Path]
    tables: dict[str, pd.DataFrame]
    source_roots: list[Path]
    root: Path


def _producer(record: dict[str, Any]) -> None:
    commit = record.get("code_commit")
    _require(
        isinstance(commit, str) and re.fullmatch(r"[0-9a-f]{40}", commit) is not None,
        "Matrix producer must name a clean Git commit",
    )
    try:
        with Repo(Path(__file__).resolve().parent, search_parent_directories=True) as repository:
            expected = {}
            for name in (*SCRIPTS, METHOD):
                item = repository.commit(commit).tree / name
                if not isinstance(item, Blob):
                    raise ValueError("Matrix producer dependency is not a Git blob")
                expected[name] = {"sha256": hashlib.sha256(item.data_stream.read()).hexdigest()}
    except (BadName, BadObject, InvalidGitRepositoryError, ValueError, KeyError) as error:
        raise ValueError("Matrix producer Git inventory is unavailable") from error
    _require(
        record.get("scripts") == {name: expected[name] for name in SCRIPTS}
        and record.get("methods") == {METHOD: expected[METHOD]}
        and record.get("script_sha256") == expected[SCRIPTS[0]]["sha256"]
        and record.get("methods_sha256") == expected[METHOD]["sha256"],
        "Matrix producer inventory differs from its Git tree",
    )


def _membership(root: Path) -> None:
    _require(
        {path.name for path in root.iterdir()} == {MANIFEST, *OUTPUTS},
        "Matrix bundle output inventory differs",
    )
    _require(
        all(
            (root / name).is_file() and not (root / name).is_symlink()
            for name in (MANIFEST, *OUTPUTS)
        ),
        "Matrix bundle has nonregular output members",
    )


def _finite(value: Any) -> float:
    _require(not isinstance(value, bool) and value != "", "Missing or boolean overlap quantity")
    try:
        number = float(value)
    except (ValueError, TypeError) as error:
        raise ValueError("Non-numeric overlap quantity") from error
    _require(math.isfinite(number), "Nonfinite overlap quantity")
    return number


def _count(value: Any) -> int:
    number = _finite(value)
    _require(number >= 0 and number.is_integer(), "Overlap counts must be nonnegative integers")
    return int(number)


def _validate_pairs(
    frame: pd.DataFrame, traits: tuple[str, ...], labels: dict[str, str], total: Any
) -> None:
    _require(tuple(frame.columns) == COLUMNS, "Matrix overlap column inventory differs")
    _require(
        list(frame[["trait_a_id", "trait_b_id"]].itertuples(index=False, name=None))
        == list(combinations(traits, 2)),
        "Matrix overlap pair inventory or order differs",
    )
    carriers: dict[str, int] = {}
    expected_total = _count(total)
    for row in frame.to_dict(orient="records"):
        a, b = str(row["trait_a_id"]), str(row["trait_b_id"])
        _require(
            row["trait_a_label"] == labels[a] and row["trait_b_label"] == labels[b],
            "Overlap label differs from verified metadata",
        )
        n, na, nb, shared, union = (
            _count(row[name])
            for name in (
                "matrix_taxa_total",
                "trait_a_taxa",
                "trait_b_taxa",
                "intersection_taxa",
                "union_taxa",
            )
        )
        _require(
            n == expected_total
            and shared <= min(na, nb)
            and max(na, nb) <= union <= n
            and union == na + nb - shared,
            "Inconsistent overlap count identities",
        )
        for trait, count in ((a, na), (b, nb)):
            _require(
                trait not in carriers or carriers[trait] == count,
                "Trait carrier count changes between pairs",
            )
            carriers[trait] = count
        if union == 0:
            _require(
                row["jaccard_fraction"] == row["jaccard_pct"] == "",
                "Empty union must preserve undefined overlap",
            )
        else:
            for name, expected in (
                ("jaccard_fraction", shared / union),
                ("jaccard_pct", 100.0 * shared / union),
            ):
                value = _finite(row[name])
                upper = 1.0 if name == "jaccard_fraction" else 100.0
                _require(0.0 <= value <= upper, "Saved Jaccard is outside its domain")
                _require(
                    math.isclose(
                        value, expected, rel_tol=RATIO_REL_TOLERANCE, abs_tol=RATIO_ABS_TOLERANCE
                    ),
                    "Saved Jaccard differs from its descriptive counts",
                )


def load_overlap_inputs(path: Path, *, manifest_sha256: str, trait_list: Path) -> OverlapInputs:
    """Verify a complete matrix sidecar, preserving the declared nine-key family."""
    manifest = path / MANIFEST if path.is_dir() else path
    snapshots: dict[Path, dict[str, Any]] = {}
    targets: dict[Path, Path] = {}
    raw = _read(manifest, snapshots, targets, digest=manifest_sha256)
    consumer = ConsumerRecord.model_validate_json(raw)
    _require(
        consumer.consumer == "matrix-overlap-coverage"
        and consumer.feature_policy == "exclude_mediadive",
        "Wrong overlap consumer contract",
    )
    root = manifest.resolve(strict=True).parent
    _membership(root)
    _require(
        len(consumer.outputs) == len(OUTPUTS)
        and {Path(item.path).name for item in consumer.outputs} == set(OUTPUTS)
        and all(
            Path(item.path).resolve(strict=True) == root / Path(item.path).name
            for item in consumer.outputs
        ),
        "Consumer output inventory differs",
    )
    buffers = {}
    for item in consumer.outputs:
        buffers[Path(item.path).name] = _read(
            Path(item.path),
            snapshots,
            targets,
            expected={"sha256": item.sha256, "size_bytes": item.size_bytes},
        )
    core_path = Path(consumer.run_manifest_path)
    core_raw = _read(core_path, snapshots, targets, digest=consumer.run_manifest_sha256)
    core = _json(core_raw)
    core_run = WorkflowRun.model_validate_json(core_raw)
    verified = verify_consumer_record(manifest)
    _require(verified == consumer, "Consumer native parse differs from pinned buffer")
    record, summary = _json(buffers["provenance.json"]), _json(buffers["summary.json"])
    _producer(record)
    _require(
        core.get("status") == "completed"
        and core.get("git_dirty") is False
        and core.get("run_id") == consumer.run_id,
        "Core is not the completed clean consumer source",
    )
    params = core.get("parameters", [])
    parameters = {item["name"]: item["value"] for item in params}
    _require(len(parameters) == len(params), "Duplicate core parameter identity")
    for name in ("contrast_population", "direction_basis"):
        _require(
            name not in parameters or parameters[name] == core.get(name),
            "Core parameter contradicts its authoritative field",
        )
    _require(
        core.get("contrast_population") == POPULATION
        and core.get("direction_basis") == "disease_vs_comparison",
        "Overlap requires disease-oriented Tier 1+2 source",
    )
    scope = parameters.get("strain_rollup")
    _require(scope in {"off", "phenotype", "full"}, "Unknown matrix strain rollup")
    policy = resolve_trait_policy(core_run).policy
    context = {
        "core_run_id": consumer.run_id,
        "strain_rollup": str(scope),
        "contrast_population": POPULATION,
        "direction_basis": "disease_vs_comparison",
        "source_commit": str(core["git_commit"]),
        "consumer_commit": str(record["code_commit"]),
        "trait_key_policy": str(policy),
    }
    for value in (record, summary):
        _require(
            value.get("schema_version") == "1.0.0"
            and value.get("consumer") == consumer.consumer
            and value.get("core_run_id") == consumer.run_id
            and value.get("workflow_id") == core["definition"]["id"]
            and value.get("strain_rollup") == scope
            and value.get("contrast_population") == POPULATION,
            "Matrix record scope differs from verified core",
        )
    _require(record.get("core_git_commit_sha") == core["git_commit"], "Matrix core commit differs")
    expected_outputs = {
        name: {"relative_path": name, "sha256": hashlib.sha256(data).hexdigest()}
        for name, data in buffers.items()
        if name not in {"provenance.json", "core_consumption.coverage.tsv"}
    }
    _require(
        record.get("outputs") == expected_outputs, "Matrix provenance output inventory differs"
    )
    list_raw = _read(trait_list, snapshots, targets)
    traits = parse_trait_list(list_raw.decode("utf-8"))
    _require(len(traits) == 9, "Overlap figure requires the explicit nine-trait list")
    digest = hashlib.sha256(list_raw).hexdigest()
    for value in (record, summary):
        declared = value.get("declared_traits", {})
        _require(
            declared.get("trait_ids") == list(traits)
            and declared.get("trait_count") == len(traits)
            and declared.get("trait_list_sha256") == digest
            and declared.get("sha256") == digest
            and declared.get("trait_list_file") == trait_list.name
            and declared.get("source_name") == trait_list.name,
            "Declared trait list differs from pinned bytes or order",
        )
    _require(
        summary["declared_traits"].get("pair_count") == len(list(combinations(traits, 2))),
        "Declared pair count differs",
    )
    artifacts = core.get("artifacts", [])
    metadata_items = [item for item in artifacts if item["product_id"] == "feature_metadata"]
    _require(len(metadata_items) == 1, "Core feature metadata identity is missing or repeated")
    item = metadata_items[0]
    metadata_path = Path(item["path"])
    if not metadata_path.is_absolute():
        metadata_path = core_path.parent / metadata_path
    metadata_raw = _read(
        metadata_path,
        snapshots,
        targets,
        expected={"sha256": item["sha256"], "size_bytes": item["size_bytes"]},
    )
    metadata_binding = record.get("consumed_products", {}).get("feature_metadata", {})
    _require(
        metadata_binding
        == {
            "relative_path": metadata_path.resolve()
            .relative_to(core_path.resolve().parent)
            .as_posix(),
            "sha256": hashlib.sha256(metadata_raw).hexdigest(),
        },
        "Metadata differs from consumed product binding",
    )
    metadata = _table(metadata_raw)
    _require(
        "column_id" in metadata and not metadata.column_id.duplicated().any(),
        "Missing or duplicate feature metadata IDs",
    )
    labels = trait_labels(metadata, traits)
    frame = _table(buffers[TABLE])
    _validate_pairs(
        frame, traits, labels, summary.get("coverage_and_cascade", {}).get("matrix_taxa_total")
    )
    roots = {root, core_path.parent, *(item.parent for item in snapshots)}
    for artifact in artifacts:
        for key in ("path", "requested_path"):
            if artifact.get(key):
                source = Path(artifact[key])
                if not source.is_absolute():
                    source = core_path.parent / source
                roots.update((source.parent, source.resolve().parent))
    result = OverlapInputs(
        record,
        summary,
        consumer,
        context,
        traits,
        labels,
        snapshots,
        targets,
        {TABLE: frame, "feature_metadata": metadata},
        sorted(roots),
        root,
    )
    recheck_overlap_inputs(result)
    return result


def recheck_overlap_inputs(bundle: OverlapInputs) -> None:
    """Check every parsed buffer and alias before publication, without replay."""
    for path, target in bundle.targets.items():
        _require(path.resolve(strict=True) == target, "Overlap input alias changed")
    for path, expected in bundle.snapshots.items():
        _read(path, {}, {}, expected=expected)
    _membership(bundle.root)


def overlap_audit(bundle: OverlapInputs) -> pd.DataFrame:
    """Preserve every original pair and quantity, adding only display identity/status."""
    result = bundle.tables[TABLE].copy()
    for column in ("jaccard_fraction", "jaccard_pct"):
        result[column] = [None if value == "" else float(value) for value in result[column]]
    result["audit_row_id"] = [f"{a}::{b}" for a, b in combinations(bundle.trait_ids, 2)]
    result["overlap_status"] = [
        "undefined_empty_union" if _count(value) == 0 else "available"
        for value in result.union_taxa
    ]
    return result


def overlap_pages(bundle: OverlapInputs) -> Iterator[Page]:
    """Display all unordered pairs, omitting the unobserved diagonal and lower triangle."""
    audit = overlap_audit(bundle)
    count = len(bundle.trait_ids)
    labels = [
        _wrap(bundle.labels[trait], 43) + "\n" + _wrap(trait, 45) for trait in bundle.trait_ids
    ]
    # Two columns of verified labels; expand the page if long source labels wrap.
    lines = [
        max(
            len(labels[index].splitlines()),
            len(labels[index + 5].splitlines()) if index + 5 < count else 0,
        )
        for index in range(5)
    ]
    legend_height = sum(lines) * 0.14 + len(lines) * 0.16
    height = 6.55 + legend_height
    with matplotlib.rc_context(STYLE):
        figure = _figure(height)
        axis = figure.add_axes((0.15, (height - 5.65) / height, 0.66, 4.75 / height))
        palette, norm = colormaps["Blues"], Normalize(0, 100)
        positions = {trait: index for index, trait in enumerate(bundle.trait_ids)}
        for row in audit.to_dict(orient="records"):
            a, b = positions[str(row["trait_a_id"])], positions[str(row["trait_b_id"])]
            missing = row["overlap_status"] == "undefined_empty_union"
            value = 0.0 if missing else _finite(row["jaccard_pct"])
            patch = Rectangle(
                (b, a),
                1,
                1,
                facecolor="#e5e5e5" if missing else palette(norm(value)),
                edgecolor="white",
                linewidth=0.6,
            )
            axis.add_patch(patch)
            axis.text(
                b + 0.5,
                a + 0.5,
                validate_authored_wording("—" if missing else f"{value:.1f}"),
                ha="center",
                va="center",
                fontsize=8,
                color="white" if value > 65 else "#263746",
            )
        axis.set(xlim=(0, count), ylim=(count, 0), aspect="equal")
        axis.set_xticks(
            [index + 0.5 for index in range(count)],
            [validate_authored_wording(str(index + 1)) for index in range(count)],
        )
        axis.set_yticks(
            [index + 0.5 for index in range(count)],
            [validate_authored_wording(str(index + 1)) for index in range(count)],
        )
        axis.xaxis.tick_top()
        axis.tick_params(length=0)
        for spine in axis.spines.values():
            spine.set_visible(False)
        color_axis = figure.add_axes((0.85, (height - 5.65) / height, 0.025, 4.75 / height))
        bar = figure.colorbar(ScalarMappable(norm=norm, cmap=palette), cax=color_axis)
        if bar.solids is not None:
            bar.solids.set_rasterized(False)
        bar.set_label(validate_authored_wording("Taxon-carrier Jaccard (%)"))
        figure.text(
            0.05,
            1 - 0.23 / height,
            validate_authored_wording(f"Taxon-carrier overlap among {count} declared traits"),
            fontsize=10,
            weight="bold",
            va="top",
        )
        figure.text(
            0.05,
            1 - 0.48 / height,
            validate_authored_wording(
                f"Rollup: {bundle.context['strain_rollup']}; "
                f"trait keys: {bundle.context['trait_key_policy']}"
            ),
            va="top",
        )
        top = height - 6.0
        for index in range(5):
            for position, x in ((index, 0.05), (index + 5, 0.53)):
                if position < count:
                    figure.text(
                        x,
                        top / height,
                        validate_authored_wording(f"{position + 1}. ") + labels[position],
                        va="top",
                        fontsize=8,
                    )
            top -= lines[index] * 0.14 + 0.16
        figure.text(
            0.05,
            0.13 / height,
            validate_authored_wording(
                "Positive annotation counts define carriers across all matrix taxa. "
                "—: empty union; diagonal not shown."
            ),
            fontsize=8,
        )
        yield Page(
            figure,
            list(audit.audit_row_id),
            len(audit),
            {
                "pair_count": len(audit),
                "trait_ids": list(bundle.trait_ids),
                "undefined_pairs": int(audit.overlap_status.eq("undefined_empty_union").sum()),
                "matrix_taxa_total": _count(audit.matrix_taxa_total.iloc[0]),
                "units": "taxon-carrier Jaccard percent",
                "context": bundle.context.copy(),
                "selection": "all unordered pairs in declared trait-list order; no diagonal",
            },
        )
