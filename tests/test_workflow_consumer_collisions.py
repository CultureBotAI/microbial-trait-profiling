"""Consumer output preflight preserves results and verified core inputs (#654)."""

from __future__ import annotations

import subprocess
import sys
from functools import partial
from pathlib import Path

import pytest

from mtp import table1_trait_maps, trait_composition, trait_screen
from mtp.paths import REPO_ROOT
from mtp.workflow.consumers import CoreBundle, open_core_bundle, write_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow


@pytest.fixture
def bundle(tmp_path: Path) -> CoreBundle:
    source = REPO_ROOT / "tests" / "fixtures" / "workflow"
    inputs = {}
    for product, filename in (
        ("raw_dump", "full_dump.csv"),
        ("kg_edges", "merged-kg_edges.tsv"),
        ("kg_nodes", "merged-kg_nodes.tsv"),
    ):
        target = tmp_path / filename
        target.write_bytes((source / filename).read_bytes())
        inputs[product] = target
    run = execute_workflow(
        REPO_ROOT / "workflows" / "bugsigdb_core.yaml",
        inputs,
        tmp_path / "runs",
        OPERATIONS,
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    return open_core_bundle(Path(run.manifest_path))


@pytest.mark.parametrize("collision", ["coverage", "record", "symlink", "hardlink"])
def test_writer_rejects_reserved_path_aliases_before_overwrite(
    bundle: CoreBundle, tmp_path: Path, collision: str
) -> None:
    record = tmp_path / "binding.json"
    reserved = record if collision == "record" else record.with_suffix(".coverage.tsv")
    reserved.write_bytes(b"previous analysis result\n")
    output = reserved
    if collision in {"symlink", "hardlink"}:
        output = tmp_path / "alias.tsv"
        if collision == "symlink":
            output.symlink_to(reserved)
        else:
            output.hardlink_to(reserved)
    with pytest.raises(ValueError, match="collid"):
        write_consumer_record(bundle, record, consumer="test", outputs=[output])
    assert reserved.read_bytes() == b"previous analysis result\n"


@pytest.mark.parametrize(
    "collision",
    [
        "coverage",
        "record",
        "input",
        "product",
        "collection",
        "input_symlink",
        "input_hardlink",
        "publication_hardlink",
    ],
)
def test_ora_cli_preflight_preserves_existing_bytes(
    bundle: CoreBundle, tmp_path: Path, collision: str
) -> None:
    report = tmp_path / "analysis" / "report.txt"
    report.parent.mkdir()
    record = report.with_suffix(".core_consumption.json")
    candidates = {
        "coverage": record.with_suffix(".coverage.tsv"),
        "record": record,
        "input": bundle.path("raw_dump"),
        "product": bundle.path("ora_profiles"),
        "collection": next(bundle.path("studies").glob("*.yaml")),
    }
    if collision in {"input_symlink", "input_hardlink", "publication_hardlink"}:
        output = report.parent / "aliased_input.tsv"
        protected = bundle.path("raw_dump")
        if collision == "publication_hardlink":
            protected = bundle.manifest.parent / "steps" / "ingest" / "_step.json"
        if collision == "input_symlink":
            output.symlink_to(protected)
        else:
            output.hardlink_to(protected)
    else:
        output = candidates[collision]
    if not output.exists():
        output.write_bytes(b"previous analysis result\n")
    before = output.read_bytes()
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "trait_enrichment.py"),
            "--core-run",
            str(bundle.manifest),
            "--output",
            str(report),
            "--tsv-output",
            str(output),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "Consumer output" in result.stderr
    assert output.read_bytes() == before
    assert not report.exists()


@pytest.mark.parametrize("consumer", ["screen", "composition", "table1"])
def test_directory_consumers_reject_immutable_run_destination(
    bundle: CoreBundle, monkeypatch: pytest.MonkeyPatch, consumer: str
) -> None:
    arguments = ["--core-run", str(bundle.manifest), "--output-dir", str(bundle.manifest.parent)]
    before = set(bundle.manifest.parent.rglob("*"))
    if consumer == "table1":
        monkeypatch.setattr(sys, "argv", ["table1_trait_maps", *arguments])
        run_consumer = table1_trait_maps.main
    elif consumer == "composition":
        run_consumer = partial(trait_composition.main, arguments)
    else:
        run_consumer = partial(trait_screen.main, arguments)
    with pytest.raises(ValueError, match="Consumer output"):
        run_consumer()
    assert set(bundle.manifest.parent.rglob("*")) == before


@pytest.mark.parametrize("alias", ["sidecar_symlink", "sidecar_hardlink", "input_hardlink"])
def test_directory_preflight_rejects_existing_output_aliases(
    bundle: CoreBundle, tmp_path: Path, alias: str
) -> None:
    directory = tmp_path / "analysis"
    directory.mkdir()
    report = directory / "trait_screen_report.md"
    if alias == "input_hardlink":
        report.hardlink_to(bundle.path("raw_dump"))
    else:
        report.write_bytes(b"previous analysis result\n")
        sidecar = directory / "core_consumption.coverage.tsv"
        if alias == "sidecar_symlink":
            sidecar.symlink_to(report)
        else:
            sidecar.hardlink_to(report)
    before = report.read_bytes()
    with pytest.raises(ValueError, match="Consumer output"):
        trait_screen.main(["--core-run", str(bundle.manifest), "--output-dir", str(directory)])
    assert report.read_bytes() == before
    assert not (directory / "trait_direction_screen.tsv").exists()


def test_directory_preflight_rejects_dangling_output_alias(
    bundle: CoreBundle, tmp_path: Path
) -> None:
    directory = tmp_path / "analysis"
    directory.mkdir()
    report = directory / "trait_screen_report.md"
    statistics = directory / "trait_direction_screen.tsv"
    report.symlink_to(statistics)
    with pytest.raises(ValueError, match="Consumer output"):
        trait_screen.main(["--core-run", str(bundle.manifest), "--output-dir", str(directory)])
    assert report.is_symlink()
    assert report.readlink() == statistics
    assert not statistics.exists()
    assert list(directory.iterdir()) == [report]
