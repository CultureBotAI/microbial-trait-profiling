"""Synthetic figure publication/provenance regressions for issue #1046.

Inputs: tiny verified-reader stand-ins and real Matplotlib pages, 2026-10-28.
Outputs: disposable PDF/SVG/PNG bundles under pytest scratch; no F/P analysis.
"""

from __future__ import annotations

import json
import re
from collections.abc import Iterator
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import Mock

import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.figure import Figure
from PIL import Image

from mtp.paths import file_sha256
from mtp.publication_figures import stage
from mtp.publication_figures.plots import Page
from mtp.three_test_table.stage import _publish_directory
from mtp.workflow import output_roots


@pytest.fixture(autouse=True)
def close_test_figures() -> Iterator[None]:
    """Own only tiny figures made by this test; publication must close them first."""
    existing = set(plt.get_fignums())
    yield
    for number in set(plt.get_fignums()) - existing:
        plt.close(number)


def identity(path: Path) -> dict[str, Any]:
    return {"sha256": file_sha256(path), "size_bytes": path.stat().st_size}


def synthetic_inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    """Stub the separately tested reader, retaining real files and identity checks."""
    source_root = tmp_path / "sources"
    records: dict[str, Path] = {}
    for role in ("core", "comparison", "noel", "jed", "ilyas"):
        directory = source_root / role
        directory.mkdir(parents=True)
        records[role] = directory / "run.json"
        records[role].write_text("{}\n", encoding="utf-8")
    records["comparison"].write_text(
        json.dumps(
            {
                "manifest_path": str(records["comparison"]),
                "core_manifest_path": str(records["core"]),
                **{f"{arm}_manifest_path": str(records[arm]) for arm in ("noel", "jed", "ilyas")},
            }
        ),
        encoding="utf-8",
    )
    parent_root = tmp_path / "three-test"
    parent_root.mkdir()
    parent_record = {
        "run_id": "three-test-fixture",
        "code_commit": "a" * 40,
        "code_dirty": False,
        "core_run_id": "core-fixture",
        "comparison_run_id": "comparison-fixture",
        "source_run_ids": {arm: f"{arm}-fixture" for arm in ("noel", "jed", "ilyas")},
        "source_commits": dict.fromkeys(records, "b" * 40),
        "strain_rollup": "phenotype",
        "scope_role": "primary",
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "direction_basis": "disease_vs_comparison",
        "analysis_bases": {"count_profiles": {"id": "basis-fixture"}},
    }
    parent_manifest = parent_root / "three_test_run.json"
    parent_manifest.write_text(json.dumps(parent_record), encoding="utf-8")
    parent_table = parent_root / "per_trait_evidence.tsv"
    parent_table.write_text(
        "trait\ttrait_label\nMETPO:1\tRecorded source label\n", encoding="utf-8"
    )
    derived_root = tmp_path / "derived"
    derived_root.mkdir()
    derived_record = {"run_id": "derived-fixture", "code_commit": "c" * 40}
    derived_manifest = derived_root / "delta_z_derived_run.json"
    derived_manifest.write_text(json.dumps(derived_record), encoding="utf-8")
    derived_table = derived_root / "quadrant_pairs.tsv"
    derived_table.write_text("study_id\texperiment_id\n001\t01\n", encoding="utf-8")
    paths = [parent_manifest, parent_table, derived_manifest, derived_table, *records.values()]
    snapshots = {path: identity(path) for path in paths}
    targets = {path: path.resolve() for path in paths}
    parent = SimpleNamespace(
        manifest=parent_manifest,
        root=parent_root,
        record=parent_record,
        sources=SimpleNamespace(),
        snapshots={path: item for path, item in snapshots.items() if path.parent != derived_root},
        resolved_targets=targets,
        evidence=pd.DataFrame({"trait": ["METPO:1"], "trait_label": ["Recorded source label"]}),
    )
    bundle = SimpleNamespace(
        parent=parent,
        record=derived_record,
        tables={},
        labels={"METPO:1": "Recorded source label"},
        snapshots=snapshots,
        targets=targets,
    )
    loader = Mock(return_value=bundle)
    monkeypatch.setattr(stage, "load_inputs", loader)

    def recheck(_bundle: Any) -> None:
        for path, expected in snapshots.items():
            if identity(path) != expected or path.resolve() != targets[path]:
                raise ValueError("Figure source identity changed")

    checker = Mock(side_effect=recheck)
    monkeypatch.setattr(stage, "recheck_inputs", checker)
    monkeypatch.setattr(
        stage, "_source_roots", lambda _sources: [path.parent for path in records.values()]
    )
    producer = {
        "code_commit": "d" * 40,
        "code_dirty": False,
        "code_inventory": {"src/mtp/publication_figures/stage.py": "e" * 64},
        "implementation_sha256": "f" * 64,
    }
    monkeypatch.setattr(stage, "_producer", lambda _commit: producer.copy())

    def pages(_bundle: Any, **_kwargs: Any) -> Iterator[Page]:
        for row_ids in (["first"], ["second", "third"]):
            figure, axis = plt.subplots(figsize=(1.6, 1.2))
            axis.plot(range(len(row_ids)), range(len(row_ids)), "o")
            axis.set_xlabel("Recorded annotation shift")
            yield Page(figure, list(row_ids), len(row_ids), {"units": "fixture annotation shift"})

    renderers = {}
    for name in ("forest_pages", "quadrant_pages", "table_pages"):
        renderer = Mock(side_effect=pages)
        monkeypatch.setattr(stage, name, renderer)
        renderers[name] = renderer
    return SimpleNamespace(
        bundle=bundle,
        three_test=parent_manifest,
        comparison=records["comparison"],
        derived=derived_manifest,
        source=derived_table,
        loader=loader,
        checker=checker,
        producer=producer,
        renderers=renderers,
        three_test_sha256=file_sha256(parent_manifest),
        derived_sha256=file_sha256(derived_manifest),
    )


def build(fixture: Any, output: Path, **kwargs: Any) -> dict[str, Any]:
    kwargs.setdefault("allow_temporary_output", True)
    return stage.build_figures(
        fixture.three_test,
        fixture.comparison,
        fixture.derived,
        output,
        three_test_sha256=fixture.three_test_sha256,
        derived_sha256=fixture.derived_sha256,
        **kwargs,
    )


def test_formats_pages_marks_and_provenance_are_all_published(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    caller_figure = plt.figure()
    initial_figures = set(plt.get_fignums())
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    output = tmp_path / "fresh" / "nested" / "figures"
    record = build(fixture, output, ranking="comparison_mean", top_n=5, dpi=300)
    assert record == json.loads((output / "figures.json").read_text())
    assert record["version"] == "publication-figures-v1"
    assert record["status"] == "completed"
    assert set(record["figures"]) == {"forest", "quadrant", "three-test"}
    for name, entry in record["figures"].items():
        provenance_path = output / name / "provenance.json"
        assert entry == {"path": f"{name}/provenance.json", **identity(provenance_path)}
        provenance = json.loads(provenance_path.read_text())
        assert provenance["figure"] == name
        assert provenance["parent_run_id"] == "three-test-fixture"
        assert provenance["derived_run_id"] == "derived-fixture"
        assert provenance["parent_code_commit"] == "a" * 40
        assert provenance["derived_code_commit"] == "c" * 40
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
        ):
            assert provenance[key] == fixture.bundle.parent.record[key]
        for key in ("code_commit", "code_dirty", "code_inventory", "implementation_sha256"):
            assert provenance[key] == fixture.producer[key]
        assert {
            item["path"]: {key: item[key] for key in ("sha256", "size_bytes")}
            for item in provenance["input_files"]
        } == {str(path): expected for path, expected in fixture.bundle.snapshots.items()}
        assert provenance["rendering"] == {
            "dpi": 300,
            "page_count": 2,
            "vector_formats": ["pdf", "svg"],
            "raster_formats": ["png"],
        }
        assert (
            provenance["selection"]
            == {
                "forest": {"ranking": "comparison_mean", "top_n": 5},
                "quadrant": {
                    "populations": ["all", "stool"],
                    "sampling": "none",
                    "stool_rule": 'body_site == "Feces"',
                },
                "three-test": {
                    "family": "complete_retained_family",
                    "ordering": "qualified_trait_id_ascending",
                },
            }[name]
        )
        assert provenance["pages"] == [
            {
                "page": 1,
                "row_ids": ["first"],
                "mark_count": 1,
                "details": {"units": "fixture annotation shift"},
            },
            {
                "page": 2,
                "row_ids": ["second", "third"],
                "mark_count": 2,
                "details": {"units": "fixture annotation shift"},
            },
        ]
        expected_files = {f"{name}.pdf"} | {
            f"page-{page:03d}.{suffix}" for page in (1, 2) for suffix in ("pdf", "svg", "png")
        }
        assert set(provenance["outputs"]) == expected_files
        assert {path.name for path in provenance_path.parent.iterdir()} == expected_files | {
            "provenance.json"
        }
        for filename, expected in provenance["outputs"].items():
            assert identity(provenance_path.parent / filename) == expected
        combined = (provenance_path.parent / f"{name}.pdf").read_bytes()
        assert combined.startswith(b"%PDF-")
        assert len(re.findall(rb"/Type\s*/Page\b", combined)) == 2
        assert "<svg" in (provenance_path.parent / "page-001.svg").read_text()
        with Image.open(provenance_path.parent / "page-001.png") as png:
            assert png.info["dpi"] == pytest.approx((300, 300), abs=0.01)
    assert set(plt.get_fignums()) == initial_figures
    assert plt.figure(caller_figure.number) is caller_figure
    fixture.loader.assert_called_once_with(
        fixture.three_test,
        fixture.comparison,
        fixture.derived,
        three_test_sha256=fixture.three_test_sha256,
        derived_sha256=fixture.derived_sha256,
    )
    assert fixture.checker.call_count >= 1


def test_only_requested_figure_is_rendered(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    record = build(fixture, tmp_path / "output", figures=("quadrant",))
    assert set(record["figures"]) == {"quadrant"}
    fixture.renderers["forest_pages"].assert_not_called()
    fixture.renderers["table_pages"].assert_not_called()
    fixture.renderers["quadrant_pages"].assert_called_once_with(fixture.bundle)


@pytest.mark.parametrize("kind", ["directory", "file", "dangling_symlink"])
def test_existing_output_is_preserved_before_loading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    output = tmp_path / "output"
    if kind == "directory":
        output.mkdir()
    elif kind == "file":
        output.write_text("unowned", encoding="utf-8")
    else:
        output.symlink_to(tmp_path / "missing")
    before = output.lstat()
    with pytest.raises(FileExistsError):
        build(fixture, output)
    assert output.lstat().st_ino == before.st_ino
    fixture.loader.assert_not_called()


@pytest.mark.parametrize(
    "location", ["three-test", "derived", "sources/core", "sources/comparison", "sources/noel"]
)
def test_source_bundle_parents_are_protected_before_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, location: str
) -> None:
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    output = tmp_path / location / "new" / "output"
    with pytest.raises(ValueError, match=r"source|disjoint|overlap"):
        build(fixture, output)
    assert not output.parent.exists()
    fixture.loader.assert_not_called()


def test_symlink_into_source_bundle_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    alias = tmp_path / "alias"
    alias.symlink_to(fixture.derived.parent, target_is_directory=True)
    with pytest.raises(ValueError, match=r"source|disjoint|overlap"):
        build(fixture, alias / "new" / "output")
    assert not (fixture.derived.parent / "new").exists()
    fixture.loader.assert_not_called()


def test_current_checkout_is_not_an_output_destination(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    output = Path(__file__).resolve().parents[1] / ".tmp" / f"forbidden-figures-{tmp_path.name}"
    assert not output.exists()
    with pytest.raises(ValueError, match=r"checkout|working tree|repository|source|disjoint"):
        build(fixture, output)
    assert not output.exists()
    fixture.loader.assert_not_called()


def test_shared_temporary_guard_requires_explicit_opt_in(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    monkeypatch.setattr(output_roots, "system_temporary_directories", lambda: [tmp_path])
    output = tmp_path / "new" / "output"
    with pytest.raises(ValueError, match=r"temporary|durable"):
        build(fixture, output, allow_temporary_output=False)
    fixture.loader.assert_not_called()
    assert not output.parent.exists()
    assert (
        build(fixture, output, figures=("forest",), allow_temporary_output=True)["status"]
        == "completed"
    )


@pytest.mark.parametrize(
    "options",
    [
        {"dpi": 299},
        {"dpi": 0},
        {"figures": ()},
        {"figures": ("unknown",)},
        {"top_n": 0},
        {"ranking": "unknown"},
    ],
)
def test_invalid_render_parameters_fail_before_source_loading(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, options: dict[str, Any]
) -> None:
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    output = tmp_path / "new" / "output"
    with pytest.raises(ValueError, match=r"dpi|DPI|figure|top_n|top.n|ranking"):
        build(fixture, output, **options)
    fixture.loader.assert_not_called()
    assert not output.parent.exists()


def test_failed_reader_does_not_create_missing_output_ancestors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    fixture.loader.side_effect = ValueError("Source verification failed")
    with pytest.raises(ValueError, match="verification"):
        build(fixture, tmp_path / "new" / "nested" / "output")
    assert not (tmp_path / "new").exists()


def test_source_mutation_during_render_refuses_complete_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    caller_figure = plt.figure()
    initial_figures = set(plt.get_fignums())
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    original = Figure.savefig

    def save_and_mutate(self: Figure, *args: Any, **kwargs: Any) -> Any:
        result = original(self, *args, **kwargs)
        fixture.source.write_text("changed source\n", encoding="utf-8")
        return result

    monkeypatch.setattr(Figure, "savefig", save_and_mutate)
    before = set(tmp_path.iterdir())
    with pytest.raises(ValueError, match=r"identity|changed"):
        build(fixture, tmp_path / "output", figures=("forest",))
    assert set(tmp_path.iterdir()) == before
    assert set(plt.get_fignums()) == initial_figures
    assert plt.figure(caller_figure.number) is caller_figure


def test_failed_page_render_cleans_owned_staging_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    caller_figure = plt.figure()
    initial_figures = set(plt.get_fignums())
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    monkeypatch.setattr(Figure, "savefig", Mock(side_effect=RuntimeError("render failed")))
    sentinel = tmp_path / "unrelated"
    sentinel.write_bytes(b"keep these bytes")
    before = set(tmp_path.iterdir())
    with pytest.raises(RuntimeError, match="render failed"):
        build(fixture, tmp_path / "output", figures=("forest",))
    assert set(tmp_path.iterdir()) == before
    assert sentinel.read_bytes() == b"keep these bytes"
    assert set(plt.get_fignums()) == initial_figures
    assert plt.figure(caller_figure.number) is caller_figure


def test_changed_rendering_producer_cannot_publish(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    calls = 0

    def producer(_commit: str | None) -> dict[str, Any]:
        nonlocal calls
        calls += 1
        return {**fixture.producer, "code_dirty": calls > 1}

    monkeypatch.setattr(stage, "_producer", producer)
    before = set(tmp_path.iterdir())
    with pytest.raises(ValueError, match=r"producer|changed"):
        build(fixture, tmp_path / "output", figures=("forest",))
    assert set(tmp_path.iterdir()) == before


@pytest.mark.parametrize("kind", ["directory", "dangling_symlink"])
def test_racing_destination_is_not_replaced(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    fixture = synthetic_inputs(tmp_path, monkeypatch)
    output = tmp_path / "output"

    def race(staging: Path, destination: Path) -> None:
        if kind == "directory":
            destination.mkdir()
            (destination / "sentinel").write_bytes(b"concurrent owner")
        else:
            destination.symlink_to(tmp_path / "unowned-missing")
        _publish_directory(staging, destination)

    monkeypatch.setattr(stage, "_publish_directory", race)
    before = set(tmp_path.iterdir())
    with pytest.raises((FileExistsError, OSError)):
        build(fixture, output, figures=("forest",))
    assert set(tmp_path.iterdir()) == before | {output}
    if kind == "directory":
        assert (output / "sentinel").read_bytes() == b"concurrent owner"
        assert not (output / "figures.json").exists()
    else:
        assert output.is_symlink()
        assert output.readlink() == tmp_path / "unowned-missing"


def test_cli_forwards_pins_selection_and_temporary_opt_in(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from mtp.publication_figures import __main__ as cli

    builder = Mock(return_value={"figures": {"forest": {}, "quadrant": {}}})
    monkeypatch.setattr(cli, "build_figures", builder)
    three_test, comparison, derived, output = (
        tmp_path / name for name in ("three-test.json", "comparison.json", "derived.json", "out")
    )
    assert (
        cli.main(
            [
                "--three-test-run",
                str(three_test),
                "--three-test-sha256",
                "1" * 64,
                "--comparison-run",
                str(comparison),
                "--derived-run",
                str(derived),
                "--derived-sha256",
                "2" * 64,
                "--output",
                str(output),
                "--figure",
                "forest",
                "--figure",
                "quadrant",
                "--ranking",
                "study_mean",
                "--top-n",
                "3",
                "--dpi",
                "400",
                "--code-commit",
                "3" * 40,
                "--allow-temporary-output",
            ]
        )
        == 0
    )
    builder.assert_called_once_with(
        three_test,
        comparison,
        derived,
        output,
        three_test_sha256="1" * 64,
        derived_sha256="2" * 64,
        figures=("forest", "quadrant"),
        ranking="study_mean",
        top_n=3,
        dpi=400,
        code_commit="3" * 40,
        allow_temporary_output=True,
    )
