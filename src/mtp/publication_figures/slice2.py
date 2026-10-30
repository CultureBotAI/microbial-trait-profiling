"""Render declared rank, candidate-triple and carrier-overlap figures (#1046).

All estimates, ranks and inference come from completed, bound tables. This entry
point adds figure outputs only; it neither changes the forest selection nor runs
the native statistical producers.
"""

from __future__ import annotations

import argparse
import shutil
import uuid
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from mtp.cellulose_decr.selection_scan_scopes import pair_scopes
from mtp.publication_figures.inputs import load_inputs, recheck_inputs
from mtp.publication_figures.overlap import (
    load_overlap_inputs,
    overlap_audit,
    overlap_pages,
    recheck_overlap_inputs,
)
from mtp.publication_figures.plots import Page
from mtp.publication_figures.rank_sensitivity import rank_sensitivity_pages, rank_sensitivity_table
from mtp.publication_figures.scan import (
    load_scan_inputs,
    recheck_scan_inputs,
    scan_audit,
    scan_pages,
)
from mtp.publication_figures.stage import _manifest, _producer, _source_roots, _write_json
from mtp.three_test_table import source_paths
from mtp.three_test_table.stage import _disjoint, _inventory, _publish_directory
from mtp.workflow.output_roots import add_temporary_output_option, require_durable_output_path

VERSION = "publication-figures-slice2-v1"
FIGURES = ("rank-sensitivity", "candidate-triples", "trait-overlap")
CONTEXT = ("core_run_id", "strain_rollup", "contrast_population", "direction_basis")


def _inputs(bundles: Sequence[Any]) -> list[dict[str, Any]]:
    snapshots: dict[Path, dict[str, Any]] = {}
    targets: dict[Path, Path] = {}
    for bundle in bundles:
        for path, identity in bundle.snapshots.items():
            if path in snapshots and (
                snapshots[path] != identity or targets[path] != bundle.targets[path]
            ):
                raise ValueError("Figure inputs disagree about shared source bytes or aliases")
            snapshots[path] = identity
            targets[path] = bundle.targets[path]
    return [{"path": str(path), **identity} for path, identity in sorted(snapshots.items())]


def _export(
    name: str,
    directory: Path,
    pages: Iterable[Page],
    audit: pd.DataFrame,
    producer: dict[str, Any],
    context: dict[str, Any],
    sources: dict[str, Any],
    selection: dict[str, Any],
    *,
    dpi: int,
) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    directory.mkdir()
    audit.to_csv(directory / "audit.tsv", sep="\t", index=False, float_format="%.17g")
    page_records = []
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
    _write_json(
        directory / "provenance.json",
        {
            "version": VERSION,
            "figure": name,
            **producer,
            "context": context,
            "source_records": {role: bundle.record for role, bundle in sources.items()},
            "source_contexts": {
                role: bundle.context
                for role, bundle in sources.items()
                if hasattr(bundle, "context")
            },
            "input_files": _inputs(list(sources.values())),
            "selection": selection,
            "audit_rows": len(audit),
            "pages": page_records,
            "rendering": {
                "dpi": dpi,
                "page_count": len(page_records),
                "vector_formats": ["pdf", "svg"],
                "raster_formats": ["png"],
            },
            "outputs": {path.name: _inventory(path) for path in sorted(directory.iterdir())},
        },
    )


def _scan_pages(bundles: Sequence[Any]) -> Iterator[Page]:
    for bundle in bundles:
        for page in scan_pages(bundle):
            page.details["source_context"] = bundle.context
            yield page


def build_slice2_figures(
    three_test: Path,
    comparison: Path,
    derived: Path,
    output: Path,
    *,
    three_test_sha256: str,
    derived_sha256: str,
    scan_run: Path | None = None,
    scan_sha256: str | None = None,
    scan_companion_run: Path | None = None,
    scan_companion_sha256: str | None = None,
    matrix_run: Path | None = None,
    matrix_sha256: str | None = None,
    trait_list: Path | None = None,
    figures: Sequence[str] = FIGURES,
    top_n: int = 5,
    dpi: int = 300,
    code_commit: str | None = None,
    allow_temporary_output: bool = False,
) -> dict[str, Any]:
    """Publish one fresh, atomic bundle, requiring one verified source unit."""
    if output.exists() or output.is_symlink():
        raise FileExistsError("Choose a new figure output directory")
    if type(dpi) is not int or dpi < 300:
        raise ValueError("PNG dpi must be an integer of at least 300")
    if type(top_n) is not int or top_n < 1:
        raise ValueError("top_n must be a positive integer")
    if (
        not figures
        or len(set(figures)) != len(figures)
        or any(name not in FIGURES for name in figures)
    ):
        raise ValueError("Choose distinct slice 2 figures")
    if "candidate-triples" in figures and (scan_run is None or scan_sha256 is None):
        raise ValueError("Candidate-triple figure requires a scan run and manifest pin")
    if (scan_companion_run is None) != (scan_companion_sha256 is None):
        raise ValueError("Companion scan requires both path and manifest pin")
    if "trait-overlap" in figures and (
        matrix_run is None or matrix_sha256 is None or trait_list is None
    ):
        raise ValueError(
            "Trait overlap requires a matrix run, manifest pin and explicit trait list"
        )
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
            _manifest(comparison, "comparison_run.json"), multiplicity=None, label_support=None
        )
    )
    roots.extend(
        path if path.is_dir() else path.parent
        for path in (scan_run, scan_companion_run, matrix_run, trait_list)
        if path is not None
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
    context = {key: bundle.parent.record[key] for key in CONTEXT}
    core_manifest_sha256 = bundle.parent.sources.comparison.core_manifest_sha256
    context["core_manifest_sha256"] = core_manifest_sha256
    # Include the original parent record: this binds comparison, arm and basis identities.
    context["three_test_record"] = bundle.parent.record
    roots.extend(_source_roots(bundle.parent.sources))
    sources: dict[str, Any] = {"derived": bundle}
    scans = []
    if "candidate-triples" in figures:
        if scan_run is None or scan_sha256 is None:
            raise ValueError("Candidate-triple figure requires a scan run and manifest pin")
        scan = load_scan_inputs(scan_run, manifest_sha256=scan_sha256)
        scans.append(scan)
        sources["scan"] = scan
    if "trait-overlap" in figures:
        if matrix_run is None or matrix_sha256 is None or trait_list is None:
            raise ValueError("Trait overlap requires a matrix run, pin and explicit trait list")
        sources["matrix"] = load_overlap_inputs(
            matrix_run, manifest_sha256=matrix_sha256, trait_list=trait_list
        )
    for role, sidecar in sources.items():
        if role == "derived":
            continue
        if any(sidecar.context.get(key) != context[key] for key in CONTEXT):
            raise ValueError(f"Figure source context differs for {role}")
        sidecar_core_pin = (
            sidecar.core_manifest_sha256 if role == "scan" else sidecar.consumer.run_manifest_sha256
        )
        if sidecar_core_pin != core_manifest_sha256:
            raise ValueError(f"Figure source core manifest differs for {role}")
        source_commit = bundle.parent.record.get("source_commits", {}).get("core")
        if (
            source_commit is not None
            and sidecar.context.get("source_commit", source_commit) != source_commit
        ):
            raise ValueError(f"Figure source context commit differs for {role}")
    if scans:
        companion_scopes = scans[0].context.get("report_beside", [])
        if companion_scopes and scan_companion_run is None:
            raise ValueError("Recorded scan protocol requires the companion scope beside this scan")
        if scan_companion_run is not None:
            if scan_companion_sha256 is None:
                raise ValueError("Companion scan requires both path and manifest pin")
            companion = load_scan_inputs(scan_companion_run, manifest_sha256=scan_companion_sha256)
            if (
                companion.context["strain_rollup"] not in companion_scopes
                or context["strain_rollup"] not in companion.context.get("report_beside", [])
                or any(
                    companion.context.get(key) != scans[0].context.get(key)
                    for key in (
                        "contrast_population",
                        "direction_basis",
                        "source_commit",
                        "protocol_version",
                    )
                )
            ):
                raise ValueError("Companion scan source context differs from protocol")
            main = scans[0]
            if main.native_run is None or companion.native_run is None:
                raise ValueError("Companion scan lacks native comparability records")
            try:
                ordered = pair_scopes(
                    [main.native_run, companion.native_run], main.native_run.protocol
                )
            except ValueError as error:
                raise ValueError(f"Companion scan native comparability failed: {error}") from error
            by_core = {item.context["core_run_id"]: item for item in (main, companion)}
            scans = [by_core[item.core_run_id] for item in ordered]
            sources["scan_companion"] = companion
    for role, sidecar in sources.items():
        if role != "derived":
            roots.extend(sidecar.source_roots)
    _inputs(list(sources.values()))
    _disjoint(original, roots)

    def recheck() -> None:
        recheck_inputs(bundle)
        for scan in scans:
            recheck_scan_inputs(scan)
        if "matrix" in sources:
            recheck_overlap_inputs(sources["matrix"])

    recheck()
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = output.parent / f".{output.name}.staging-{uuid.uuid4().hex}"
    staging.mkdir()
    try:
        entries = {}
        for name in figures:
            if name == "rank-sensitivity":
                audit = rank_sensitivity_table(bundle, top_n=top_n)
                pages = rank_sensitivity_pages(bundle, top_n=top_n)
                selected_sources = {"derived": bundle}
                selection = {
                    "headline_ranking": "experiment_mean_delta_z",
                    "sensitivity_ranking": "study_mean_delta_z",
                    "top_n": top_n,
                    "ties": "retained",
                    "audit_family": "complete_retained_family",
                }
            elif name == "candidate-triples":
                audit = pd.concat(
                    [
                        scan_audit(item).assign(strain_rollup=item.context["strain_rollup"])
                        for item in scans
                    ],
                    ignore_index=True,
                )
                pages = _scan_pages(scans)
                selected_sources = {
                    key: value
                    for key, value in sources.items()
                    if key in {"derived", "scan", "scan_companion"}
                }
                selection = {
                    "reference_set": "all_tested_triples",
                    "sampling": "none",
                    "quantity": "candidate_triple_distribution",
                    "scopes": [item.context["strain_rollup"] for item in scans],
                }
            else:
                matrix = sources["matrix"]
                audit, pages = overlap_audit(matrix), overlap_pages(matrix)
                selected_sources = {"derived": bundle, "matrix": matrix}
                selection = {
                    "traits": "explicit_nine_key_list",
                    "pairs": "all_unordered_pairs",
                    "sampling": "none",
                }
            _export(
                name,
                staging / name,
                pages,
                audit,
                producer,
                context,
                selected_sources,
                selection,
                dpi=dpi,
            )
            path = staging / name / "provenance.json"
            entries[name] = {"path": f"{name}/provenance.json", **_inventory(path)}
        record = {"version": VERSION, "status": "completed", "figures": entries}
        _write_json(staging / "figures.json", record)
        recheck()
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("three-test-run", "comparison-run", "derived-run", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    for name in ("three-test-sha256", "derived-sha256"):
        parser.add_argument(f"--{name}", required=True)
    for name in ("scan-run", "scan-companion-run", "matrix-run", "trait-list"):
        parser.add_argument(f"--{name}", type=Path)
    for name in ("scan-sha256", "scan-companion-sha256", "matrix-sha256", "code-commit"):
        parser.add_argument(f"--{name}")
    parser.add_argument("--figure", action="append", choices=FIGURES, dest="figures")
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--dpi", type=int, default=300)
    add_temporary_output_option(parser)
    args = parser.parse_args(argv)
    build_slice2_figures(
        args.three_test_run,
        args.comparison_run,
        args.derived_run,
        args.output,
        three_test_sha256=args.three_test_sha256,
        derived_sha256=args.derived_sha256,
        scan_run=args.scan_run,
        scan_sha256=args.scan_sha256,
        scan_companion_run=args.scan_companion_run,
        scan_companion_sha256=args.scan_companion_sha256,
        matrix_run=args.matrix_run,
        matrix_sha256=args.matrix_sha256,
        trait_list=args.trait_list,
        figures=tuple(args.figures) if args.figures else FIGURES,
        top_n=args.top_n,
        dpi=args.dpi,
        code_commit=args.code_commit,
        allow_temporary_output=args.allow_temporary_output,
    )
    print(f"Slice 2 figures published: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
