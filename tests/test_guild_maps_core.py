"""Guild figures bind actual KG provenance and protect verified inputs (#658)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from mtp import guild_maps
from mtp.corpus_identity import CorpusIdentity
from mtp.paths import REPO_ROOT
from mtp.workflow.consumers import CoreBundle, open_core_bundle, verify_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow


@pytest.fixture
def bundle(tmp_path: Path) -> CoreBundle:
    fixture = REPO_ROOT / "tests/fixtures/workflow"
    nodes = tmp_path / "nodes.tsv"
    edges = tmp_path / "edges.tsv"
    nodes.write_text(
        (fixture / "merged-kg_nodes.tsv").read_text()
        + "METPO:1000603\tbiolink:PhenotypicQuality\tanaerobic\n"
    )
    edges.write_text(
        (fixture / "merged-kg_edges.tsv").read_text()
        + "anaerobe\tNCBITaxon:300\tbiolink:has_phenotype\tMETPO:1000603\tfixture:annotation\n"
    )
    run = execute_workflow(
        REPO_ROOT / "workflows/bugsigdb_core.yaml",
        {
            "raw_dump": fixture / "full_dump.csv",
            "kg_edges": edges,
            "kg_nodes": nodes,
        },
        tmp_path / "core",
        OPERATIONS,
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    return open_core_bundle(Path(run.manifest_path))


def test_core_guild_identity_never_uses_checkout_sidecar(
    bundle: CoreBundle, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "guilds"

    def reject_default_identity(*_args: object, **_kwargs: object) -> CorpusIdentity:
        raise AssertionError("The core adapter consulted the default KG sidecar")

    def write_render_stub(
        _matrix: guild_maps.GuildMatrix, path: Path, _vintage: str | None
    ) -> Path:
        path.write_text("renderer excluded from numerical provenance test\n")
        return path

    monkeypatch.setattr(guild_maps, "corpus_identity", reject_default_identity)
    monkeypatch.setattr(guild_maps, "write_figure", write_render_stub)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "guild_maps",
            "--core-run",
            str(bundle.manifest),
            "--output-dir",
            str(output),
            "--min-studies",
            "1",
        ],
    )
    assert guild_maps.main() == 0
    record = verify_consumer_record(output / "core_consumption.json")
    assert record.consumer == "mtp.guild_maps"
    identity = json.loads((output / "corpus_identity.json").read_text())
    assert identity["matrix"]["path"] == str(bundle.path("count_profiles"))
    kg = identity["kg_microbe"]
    assert kg["release_tag"] is None
    assert kg["provenance_path"] == str(bundle.path("source_identity"))
    assert {entry["path"] for entry in kg["inputs"]} == {
        str(bundle.path("kg_edges")),
        str(bundle.path("kg_nodes")),
    }


@pytest.mark.parametrize("failure", ["input_conflict", "output_collision"])
def test_core_guild_rejects_conflicts_before_analysis(
    bundle: CoreBundle, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    output = bundle.manifest.parent if failure == "output_collision" else tmp_path / "guilds"
    arguments = ["guild_maps", "--core-run", str(bundle.manifest), "--output-dir", str(output)]
    if failure == "input_conflict":
        arguments.extend(("--full-dump-input", str(tmp_path / "unrelated.csv")))
    monkeypatch.setattr(sys, "argv", arguments)
    with pytest.raises(ValueError, match=r"conflicts|core|input"):
        guild_maps.main()
