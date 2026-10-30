"""Publication and source-boundary tests for #1046 slice 2; synthetic inputs only."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import Mock

import pandas as pd
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from PIL import Image

from mtp.cellulose_decr import selection_scan as native_scan
from mtp.cellulose_decr.selection_scan_scopes import ScanRun
from mtp.paths import file_sha256
from mtp.publication_figures import slice2
from mtp.publication_figures.plots import Page
from mtp.three_test_table import source_paths


@pytest.fixture
def sources(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    root = tmp_path / "sources"
    root.mkdir()
    source = root / "input.tsv"
    source.write_text("trait\tvalue\nMETPO:1\t0.25\n")
    context = {
        "core_run_id": "core-1",
        "strain_rollup": "phenotype",
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "direction_basis": "disease_vs_comparison",
    }
    snapshots = {source: {"sha256": file_sha256(source), "size_bytes": source.stat().st_size}}
    targets = {source: source.resolve()}
    parent = SimpleNamespace(
        record={
            **context,
            "run_id": "parent-1",
            "code_commit": "a" * 40,
            "source_commits": {"core": "b" * 40},
        },
        sources=SimpleNamespace(comparison=SimpleNamespace(core_manifest_sha256="e" * 64)),
    )
    bundle = SimpleNamespace(
        parent=parent, record={"run_id": "derived-1"}, snapshots=snapshots, targets=targets
    )
    scan = SimpleNamespace(
        record={"run_id": "scan-1"},
        context=dict(context),
        snapshots=snapshots,
        targets=targets,
        source_roots=[root],
        core_manifest_sha256="e" * 64,
        pairing_identity=None,
        native_run=None,
    )
    overlap = SimpleNamespace(
        record={"run_id": "matrix-1"},
        context=dict(context),
        snapshots=snapshots,
        targets=targets,
        source_roots=[root],
        consumer=SimpleNamespace(run_manifest_sha256="e" * 64),
    )
    loaders = {}
    for name, value in (
        ("load_inputs", bundle),
        ("load_scan_inputs", scan),
        ("load_overlap_inputs", overlap),
    ):
        loaders[name] = Mock(return_value=value)
        monkeypatch.setattr(slice2, name, loaders[name])

    def recheck(value: Any) -> None:
        for path, item in value.snapshots.items():
            if file_sha256(path) != item["sha256"] or path.resolve() != value.targets[path]:
                raise ValueError("Input changed")

    for name in ("recheck_inputs", "recheck_scan_inputs", "recheck_overlap_inputs"):
        monkeypatch.setattr(slice2, name, recheck)
    monkeypatch.setattr(slice2, "_source_roots", lambda _sources: [root])
    monkeypatch.setattr(source_paths, "manifest_roots", lambda *_args, **_kwargs: [root])
    monkeypatch.setattr(source_paths, "declared_source_roots", lambda *_args, **_kwargs: [root])
    producer = {"code_commit": "c" * 40, "code_dirty": False, "implementation_sha256": "d" * 64}
    monkeypatch.setattr(slice2, "_producer", lambda _commit: producer.copy())
    audit = pd.DataFrame({"trait": ["METPO:1", "METPO:2"], "value": [0.25, 0.5]})
    for name in ("rank_sensitivity_table", "scan_audit", "overlap_audit"):
        monkeypatch.setattr(slice2, name, Mock(return_value=audit))

    def pages(_bundle: Any, **_kwargs: Any) -> list[Page]:
        figure = Figure(figsize=(2, 1))
        FigureCanvasAgg(figure)
        figure.text(0.1, 0.5, "Recorded annotation labels")
        return [Page(figure, ["METPO:1", "METPO:2"], 2, {"sampling": "none"})]

    for name in ("rank_sensitivity_pages", "scan_pages", "overlap_pages"):
        monkeypatch.setattr(slice2, name, pages)
    args = {
        "three_test": root / "three_test_run.json",
        "comparison": root / "comparison_run.json",
        "derived": root / "delta_z_derived_run.json",
        "three_test_sha256": "a" * 64,
        "derived_sha256": "b" * 64,
        "scan_run": root / "scan",
        "scan_sha256": "c" * 64,
        "matrix_run": root / "matrix",
        "matrix_sha256": "d" * 64,
        "trait_list": root / "traits.txt",
        "allow_temporary_output": True,
    }
    return SimpleNamespace(
        args=args,
        source=source,
        bundle=bundle,
        scan=scan,
        overlap=overlap,
        scan_loader=loaders["load_scan_inputs"],
        overlap_loader=loaders["load_overlap_inputs"],
        pages=pages,
        producer=producer,
    )


def native_run(scan: Any) -> ScanRun:
    """Construct real native records; publication tests exercise the native checker."""
    identity = scan.pairing_identity
    protocol = native_scan.load_protocol()
    return ScanRun(
        scan.source_roots[0] / scan.context["strain_rollup"],
        {
            "rule": identity["rule"],
            "strain_rollup": {
                "scope": scan.context["strain_rollup"],
                "role": "primary"
                if scan.context["strain_rollup"] == "phenotype"
                else "sensitivity",
            },
        },
        {
            **{
                key: identity[key]
                for key in ("protocol", "code_commit", "scripts", "environment", "methods", "chunk")
            },
            "core_run": {
                "contrast_population": identity["core population"][0],
                "direction_basis": identity["core population"][1],
                "git_commit": identity["core git_commit"][0],
                "git_dirty": identity["core git_commit"][1],
                "parameters": {
                    **identity["core parameters"],
                    "strain_rollup": scan.context["strain_rollup"],
                },
                "sources": identity["core sources"],
            },
        },
        pd.DataFrame(),
        protocol,
        {},
        scan.context["core_run_id"],
    )


def companion_scan(sources: Any, monkeypatch: pytest.MonkeyPatch) -> Any:
    """Model the complete native comparable contract, with distinct scoped cores."""
    sources.scan.context.update(protocol_version="v2", report_beside=["off"])
    sources.scan.pairing_identity = {
        "protocol": {"version": "v2", "sha256": "1" * 64},
        "code_commit": "c" * 40,
        "scripts": {"selection_scan.py": "2" * 64},
        "environment": {"numpy": "2.4.4"},
        "methods": {"selection_scan.md": "3" * 64},
        "chunk": 32,
        "rule": {"seed": 0, "permutations": 999},
        "core population": (sources.scan.context["contrast_population"], "disease_vs_comparison"),
        "core git_commit": ("b" * 40, False),
        "core parameters": {"contrast_population": sources.scan.context["contrast_population"]},
        "core sources": {"raw_dump": "4" * 64, "kg_edges": "5" * 64},
    }
    companion = SimpleNamespace(**vars(sources.scan))
    companion.pairing_identity = deepcopy(sources.scan.pairing_identity)
    companion.core_manifest_sha256 = "f" * 64
    companion.context = {
        **sources.scan.context,
        "strain_rollup": "off",
        "core_run_id": "core-off",
        "report_beside": ["phenotype"],
    }
    sources.scan.native_run = native_run(sources.scan)
    companion.native_run = native_run(companion)
    monkeypatch.setattr(slice2, "load_scan_inputs", Mock(side_effect=[sources.scan, companion]))
    return companion


def test_slice2_exports_complete_audits_pages_and_bound_provenance(
    sources: Any, tmp_path: Path
) -> None:
    output = tmp_path / "figures"
    result = slice2.build_slice2_figures(output=output, **sources.args)
    assert set(result["figures"]) == set(slice2.FIGURES)
    for name, entry in result["figures"].items():
        path = output / entry["path"]
        assert file_sha256(path) == entry["sha256"]
        record = json.loads(path.read_bytes())
        assert record["context"]["core_run_id"] == "core-1"
        assert record["context"]["direction_basis"] == "disease_vs_comparison"
        assert record["code_commit"] == "c" * 40
        assert record["pages"][0]["row_ids"] == ["METPO:1", "METPO:2"]
        assert record["pages"][0]["mark_count"] == 2
        assert len(pd.read_csv(path.parent / "audit.tsv", sep="\t")) == 2
        assert record["input_files"][0]["sha256"] == file_sha256(sources.source)
        for filename in (
            f"{name}.pdf",
            "page-001.pdf",
            "page-001.svg",
            "page-001.png",
            "audit.tsv",
        ):
            artifact = path.parent / filename
            assert artifact.stat().st_size > 0
            assert record["outputs"][filename]["sha256"] == file_sha256(artifact)
        with Image.open(path.parent / "page-001.png") as image:
            assert image.info["dpi"][0] >= 299.9


@pytest.mark.parametrize("role", ["scan", "overlap"])
@pytest.mark.parametrize(
    "field", ["core_run_id", "strain_rollup", "contrast_population", "direction_basis"]
)
def test_slice2_refuses_mixed_source_units(
    sources: Any, tmp_path: Path, role: str, field: str
) -> None:
    getattr(sources, role).context[field] = "different"
    output = tmp_path / "figures"
    with pytest.raises(ValueError, match="source context"):
        slice2.build_slice2_figures(output=output, **sources.args)
    assert not output.exists()


def test_slice2_rechecks_inputs_after_rendering(
    sources: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = sources.pages

    def changing(bundle: Any, **kwargs: Any) -> Any:
        pages = original(bundle, **kwargs)
        sources.source.write_text("trait\tvalue\nMETPO:1\t9.25\n")
        return pages

    monkeypatch.setattr(slice2, "rank_sensitivity_pages", changing)
    output = tmp_path / "figures"
    with pytest.raises(ValueError, match="Input changed"):
        slice2.build_slice2_figures(output=output, **sources.args)
    assert not output.exists()
    assert not list(tmp_path.glob(".figures.staging-*"))


def test_slice2_refuses_existing_output(sources: Any, tmp_path: Path) -> None:
    output = tmp_path / "figures"
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("existing evidence")
    with pytest.raises(FileExistsError):
        slice2.build_slice2_figures(output=output, **sources.args)
    assert marker.read_text() == "existing evidence"


def test_slice2_refuses_output_inside_source(sources: Any) -> None:
    with pytest.raises(ValueError, match=r"overlap|source"):
        slice2.build_slice2_figures(output=sources.source.parent / "figures", **sources.args)


def test_rank_only_does_not_require_other_sidecars(sources: Any, tmp_path: Path) -> None:
    args = {
        key: value
        for key, value in sources.args.items()
        if key not in {"scan_run", "scan_sha256", "matrix_run", "matrix_sha256", "trait_list"}
    }
    result = slice2.build_slice2_figures(
        output=tmp_path / "rank", figures=("rank-sensitivity",), **args
    )
    assert set(result["figures"]) == {"rank-sensitivity"}
    sources.scan_loader.assert_not_called()
    sources.overlap_loader.assert_not_called()


def test_slice2_refuses_low_resolution(sources: Any, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="300"):
        slice2.build_slice2_figures(output=tmp_path / "figures", dpi=299, **sources.args)


def test_v2_scan_requires_recorded_companion(sources: Any, tmp_path: Path) -> None:
    sources.scan.context.update(protocol_version="v2", report_beside=["off"])
    with pytest.raises(ValueError, match="requires the companion"):
        slice2.build_slice2_figures(output=tmp_path / "figures", **sources.args)


def test_v2_scan_exports_both_scopes(
    sources: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    companion_scan(sources, monkeypatch)
    output = tmp_path / "figures"
    slice2.build_slice2_figures(
        output=output,
        scan_companion_run=sources.source.parent / "companion",
        scan_companion_sha256="f" * 64,
        **sources.args,
    )
    record = json.loads((output / "candidate-triples/provenance.json").read_bytes())
    assert record["selection"]["scopes"] == ["phenotype", "off"]
    assert [page["details"]["source_context"]["strain_rollup"] for page in record["pages"]] == [
        "phenotype",
        "off",
    ]
    audit = pd.read_csv(output / "candidate-triples/audit.tsv", sep="\t")
    assert audit.groupby("strain_rollup").size().to_dict() == {"off": 2, "phenotype": 2}


@pytest.mark.parametrize("role", ["scan", "overlap"])
def test_sidecar_same_run_id_with_different_core_manifest_pin_refuses(
    sources: Any, tmp_path: Path, role: str
) -> None:
    if role == "scan":
        sources.scan.core_manifest_sha256 = "f" * 64
    else:
        sources.overlap.consumer.run_manifest_sha256 = "f" * 64
    output = tmp_path / "fresh" / "figures"
    with pytest.raises(ValueError, match="core manifest"):
        slice2.build_slice2_figures(output=output, **sources.args)
    assert not output.parent.exists()


@pytest.mark.parametrize(
    "field",
    [
        "protocol",
        "code_commit",
        "scripts",
        "environment",
        "methods",
        "chunk",
        "rule",
        "core population",
        "core git_commit",
        "core parameters",
        "core sources",
    ],
)
def test_v2_companion_must_share_every_native_comparability_field(
    sources: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, field: str
) -> None:
    companion = companion_scan(sources, monkeypatch)
    previous = companion.pairing_identity[field]
    if isinstance(previous, dict):
        companion.pairing_identity[field] = {**previous, "changed": True}
    elif isinstance(previous, tuple):
        companion.pairing_identity[field] = ("changed", *previous[1:])
    elif isinstance(previous, int):
        companion.pairing_identity[field] = previous + 1
    else:
        companion.pairing_identity[field] = "d" * 40
    companion.native_run = native_run(companion)
    output = tmp_path / "fresh" / "figures"
    with pytest.raises(ValueError, match="native comparability"):
        slice2.build_slice2_figures(
            output=output,
            scan_companion_run=sources.source.parent / "companion",
            scan_companion_sha256="f" * 64,
            **sources.args,
        )
    assert not output.parent.exists()


@pytest.mark.parametrize("role", ["scan", "companion"])
def test_v2_companion_requires_a_native_pairing_identity(
    sources: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, role: str
) -> None:
    companion = companion_scan(sources, monkeypatch)
    (sources.scan if role == "scan" else companion).native_run = None
    output = tmp_path / "fresh" / "figures"
    with pytest.raises(ValueError, match="native comparability"):
        slice2.build_slice2_figures(
            output=output,
            scan_companion_run=sources.source.parent / "companion",
            scan_companion_sha256="f" * 64,
            **sources.args,
        )
    assert not output.parent.exists()


def test_v2_companion_requires_distinct_scoped_core_runs(
    sources: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    companion = companion_scan(sources, monkeypatch)
    companion.context["core_run_id"] = sources.scan.context["core_run_id"]
    companion.native_run = native_run(companion)
    output = tmp_path / "fresh" / "figures"
    with pytest.raises(ValueError, match="same core run"):
        slice2.build_slice2_figures(
            output=output,
            scan_companion_run=sources.source.parent / "companion",
            scan_companion_sha256="f" * 64,
            **sources.args,
        )
    assert not output.parent.exists()


def test_historical_v1_scan_remains_standalone_without_pairing_identity(
    sources: Any, tmp_path: Path
) -> None:
    sources.scan.context.update(protocol_version="v1", report_beside=[])
    sources.scan.pairing_identity = None
    result = slice2.build_slice2_figures(
        output=tmp_path / "historical", figures=("candidate-triples",), **sources.args
    )
    assert set(result["figures"]) == {"candidate-triples"}


def test_slice2_checks_producer_again_before_publishing(
    sources: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    producer = sources.producer.copy()
    monkeypatch.setattr(
        slice2, "_producer", Mock(side_effect=[producer, {**producer, "code_dirty": True}])
    )
    output = tmp_path / "figures"
    with pytest.raises(ValueError, match="producer changed"):
        slice2.build_slice2_figures(output=output, **sources.args)
    assert not output.exists()


def test_slice2_cli_forwards_explicit_pins_and_selection(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    builder = Mock(return_value={})
    monkeypatch.setattr(slice2, "build_slice2_figures", builder)
    args = ["--figure", "rank-sensitivity"]
    for flag in ("three-test-run", "comparison-run", "derived-run", "output"):
        args.extend([f"--{flag}", str(tmp_path / flag)])
    args.extend(["--three-test-sha256", "a" * 64, "--derived-sha256", "b" * 64])
    assert slice2.main(args) == 0
    assert builder.call_args.kwargs["figures"] == ("rank-sensitivity",)
    assert builder.call_args.kwargs["three_test_sha256"] == "a" * 64
    assert builder.call_args.kwargs["derived_sha256"] == "b" * 64
