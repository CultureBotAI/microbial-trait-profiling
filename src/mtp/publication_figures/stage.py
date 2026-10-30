"""Publish figure files and their input/rendering provenance without rerunning analysis."""

from __future__ import annotations

import json
import shutil
import uuid
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from git import Repo

from mtp.multiple_testing import canonical_sha256
from mtp.paths import file_sha256
from mtp.publication_figures.inputs import FigureInputs, load_inputs, recheck_inputs
from mtp.publication_figures.plots import forest_pages, quadrant_pages, table_pages
from mtp.three_test_table import source_paths
from mtp.three_test_table.stage import _disjoint, _inventory, _publish_directory
from mtp.three_test_table.stage import _producer as _parent_producer
from mtp.workflow.output_roots import require_durable_output_path

VERSION = "publication-figures-v1"
METHODS = "methods/publication_figures.md"
FIGURES = ("forest", "quadrant", "three-test")


def _producer(requested_commit: str | None) -> dict[str, Any]:
    """Bind tracked dependencies and new figure code, including uncommitted development."""
    import matplotlib

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
    record["environment"]["matplotlib"] = matplotlib.__version__
    return record


def _source_roots(sources: Any) -> list[Path]:
    return source_paths.source_roots(sources, multiplicity=None, label_support=None)


def _manifest(path: Path, name: str) -> Path:
    return path / name if path.is_dir() else path


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )


def _figure(
    name: str,
    bundle: FigureInputs,
    directory: Path,
    producer: dict[str, Any],
    *,
    ranking: str,
    top_n: int,
    dpi: int,
) -> dict[str, Any]:
    """Stream pages so a complete trait table does not retain every figure in memory."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    directory.mkdir()
    pages = (
        forest_pages(bundle, ranking=ranking, top_n=top_n)
        if name == "forest"
        else quadrant_pages(bundle)
        if name == "quadrant"
        else table_pages(bundle)
    )
    page_records: list[dict[str, Any]] = []
    metadata = {"Creator": VERSION, "CreationDate": None, "ModDate": None}
    with (
        mpl.rc_context({"pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": VERSION}),
        PdfPages(directory / f"{name}.pdf", metadata=metadata) as combined,
    ):
        for number, page in enumerate(pages, 1):
            try:
                stem = directory / f"page-{number:03d}"
                page.figure.savefig(stem.with_suffix(".pdf"), metadata=metadata)
                page.figure.savefig(
                    stem.with_suffix(".svg"), metadata={"Date": None, "Creator": VERSION}
                )
                page.figure.savefig(
                    stem.with_suffix(".png"), dpi=dpi, metadata={"Software": VERSION}
                )
                combined.savefig(page.figure)
                page_records.append(
                    {
                        "page": number,
                        "row_ids": page.row_ids,
                        "mark_count": page.mark_count,
                        "details": page.details,
                    }
                )
            finally:
                plt.close(page.figure)
    if not page_records:
        raise ValueError(f"Figure {name} produced no pages")
    parent = bundle.parent.record
    record = {
        "version": VERSION,
        "figure": name,
        **producer,
        **{
            key: parent[key]
            for key in (
                "core_run_id",
                "comparison_run_id",
                "source_run_ids",
                "source_commits",
                "strain_rollup",
                "scope_role",
                "contrast_population",
                "direction_basis",
                "analysis_bases",
            )
        },
        "parent_run_id": parent["run_id"],
        "derived_run_id": bundle.record["run_id"],
        "parent_code_commit": parent["code_commit"],
        "derived_code_commit": bundle.record["code_commit"],
        "input_files": [
            {"path": str(path), **identity} for path, identity in sorted(bundle.snapshots.items())
        ],
        "selection": {
            "forest": {"ranking": ranking, "top_n": top_n},
            "quadrant": {
                "populations": ["all", "stool"],
                "sampling": "none",
                "stool_rule": 'body_site == "Feces"',
            },
            "three-test": {
                "family": "complete_retained_family",
                "ordering": "qualified_trait_id_ascending",
            },
        }[name],
        "rendering": {
            "dpi": dpi,
            "page_count": len(page_records),
            "vector_formats": ["pdf", "svg"],
            "raster_formats": ["png"],
        },
        "pages": page_records,
        "outputs": {path.name: _inventory(path) for path in sorted(directory.iterdir())},
    }
    _write_json(directory / "provenance.json", record)
    return record


def build_figures(
    three_test: Path,
    comparison: Path,
    derived: Path,
    output: Path,
    *,
    three_test_sha256: str,
    derived_sha256: str,
    figures: Sequence[str] = FIGURES,
    ranking: str = "comparison_mean",
    top_n: int = 5,
    dpi: int = 300,
    code_commit: str | None = None,
    allow_temporary_output: bool = False,
) -> dict[str, Any]:
    """Read verified table buffers and atomically publish a fresh, disjoint figure bundle."""
    if output.exists() or output.is_symlink():
        raise FileExistsError("Choose a new figure output directory")
    if isinstance(dpi, bool) or not isinstance(dpi, int) or dpi < 300:
        raise ValueError("PNG dpi must be an integer of at least 300")
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
        raise ValueError("top_n must be a positive integer")
    if ranking not in {"comparison_mean", "study_mean"}:
        raise ValueError("ranking must be comparison_mean or study_mean")
    if (
        not figures
        or len(set(figures)) != len(figures)
        or any(name not in FIGURES for name in figures)
    ):
        raise ValueError("Choose distinct figure names: forest, quadrant, three-test")
    original = output
    output = require_durable_output_path(output, allowed=allow_temporary_output)
    roots = [Path(__file__).resolve().parents[3]]
    roots.extend(
        source_paths.manifest_roots(
            [
                _manifest(three_test, "three_test_run.json"),
                _manifest(derived, "delta_z_derived_run.json"),
            ],
            multiplicity=None,
            label_support=None,
        )
    )
    roots.extend(
        source_paths.declared_source_roots(
            _manifest(comparison, "comparison_run.json"),
            multiplicity=None,
            label_support=None,
        )
    )
    _disjoint(original, roots)
    producer = _producer(code_commit)
    bundle = load_inputs(
        three_test,
        comparison,
        derived,
        three_test_sha256=three_test_sha256,
        derived_sha256=derived_sha256,
    )
    roots.extend(_source_roots(bundle.parent.sources))
    roots.extend(
        source_paths.manifest_roots(list(bundle.snapshots), multiplicity=None, label_support=None)
    )
    _disjoint(original, roots)
    recheck_inputs(bundle)
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = output.parent / f".{output.name}.staging-{uuid.uuid4().hex}"
    staging.mkdir()
    try:
        entries = {}
        for name in figures:
            _figure(name, bundle, staging / name, producer, ranking=ranking, top_n=top_n, dpi=dpi)
            path = staging / name / "provenance.json"
            entries[name] = {"path": f"{name}/provenance.json", **_inventory(path)}
        record = {"version": VERSION, "status": "completed", "figures": entries}
        _write_json(staging / "figures.json", record)
        recheck_inputs(bundle)
        if _producer(code_commit) != producer:
            raise ValueError("Figure rendering producer changed during execution")
        if original.resolve() != output:
            raise ValueError("Figure destination alias changed during execution")
        _disjoint(original, roots)
        _publish_directory(staging, output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return record
