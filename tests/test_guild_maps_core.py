"""Guild figures bind actual KG provenance and protect verified inputs (#658)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from mtp import guild_maps
from mtp.corpus_identity import CorpusIdentity
from mtp.paths import REPO_ROOT
from mtp.trait_guilds import CorpusStats
from mtp.workflow.consumers import CoreBundle, open_core_bundle, verify_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_execution_product
from mtp.workflow.runtime import execute_workflow


@pytest.mark.parametrize(
    ("include_non_human", "has_mouse_profiles", "expected_profiles"),
    [
        pytest.param(False, False, 2, id="default-human"),
        pytest.param(False, True, 2, id="default-mixed"),
        pytest.param(True, True, 4, id="include-mixed"),
        pytest.param(True, False, 0, id="refuse-human-with-mouse-metadata"),
    ],
)
def test_guild_cli_requires_nonhuman_rows_when_explicitly_requested(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    include_non_human: bool,
    has_mouse_profiles: bool,
    expected_profiles: int,
) -> None:
    """Use the real parser and loader; stop before identity/plotting after loading."""
    profiles = tmp_path / "profiles.tsv"
    experiments = tmp_path / "experiments.tsv"
    output = tmp_path / "output"
    studies = ("human", "mouse") if has_mouse_profiles else ("human",)
    profiles.write_text(
        "study_id\texperiment_id\tinc_dec\tMETPO:1000602\tMETPO:1000603\n"
        + "".join(
            f"{study}\t{study}-1\t{direction}\t3\t2\n"
            for study in studies
            for direction in ("inc", "dec")
        ),
        encoding="utf-8",
    )
    experiments.write_text(
        "study_id\texperiment_id\tcondition\tbody_site\thost_species\n"
        "human\thuman-1\tColitis\tFeces\tHomo sapiens\n"
        "mouse\tmouse-1\tColitis\tFeces\tMus musculus\n",
        encoding="utf-8",
    )
    before = (profiles.read_bytes(), experiments.read_bytes())
    observed: list[CorpusStats] = []

    class LoadedBoundary(Exception):
        pass

    def stop_after_loading(stats: CorpusStats, _vintage: str | None) -> None:
        assert expected_profiles > 0, "Explicit non-human request crossed the loader guard"
        assert stats.kept_profiles == expected_profiles
        observed.append(stats)
        raise LoadedBoundary

    monkeypatch.setattr(guild_maps, "write_stats_report", stop_after_loading)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "guild-maps",
            "--profile-input",
            str(profiles),
            "--experiments-input",
            str(experiments),
            "--full-dump-input",
            str(tmp_path / "absent-dump.csv"),
            "--output-dir",
            str(output),
            *(["--include-non-human"] if include_non_human else []),
        ],
    )
    if expected_profiles == 0:
        with pytest.raises(ValueError, match=r"no .*non-human.*all-host"):
            guild_maps.main()
        assert observed == []
    else:
        with pytest.raises(LoadedBoundary):
            guild_maps.main()
        assert len(observed) == 1
    assert not output.exists()
    assert (profiles.read_bytes(), experiments.read_bytes()) == before


def test_guild_cli_help_states_the_all_host_input_requirement() -> None:
    help_text = " ".join(guild_maps.build_parser().format_help().split()).lower()
    assert "--include-non-human" in help_text
    assert "all-host profile table" in help_text


@pytest.fixture
def bundle(tmp_path: Path) -> CoreBundle:
    fixture = REPO_ROOT / "tests/fixtures/workflow"
    nodes = tmp_path / "nodes.tsv"
    edges = tmp_path / "edges.tsv"
    # Each oxygen guild needs a matrix column. The fixture KG's own oxygen term is
    # METPO:1000608, which is in neither guild since #1237, so one aerobic and one
    # anaerobic annotation are added.
    nodes.write_text(
        (fixture / "merged-kg_nodes.tsv").read_text()
        + "METPO:1000603\tbiolink:PhenotypicQuality\tanaerobic\n"
        + "METPO:1000602\tbiolink:PhenotypicQuality\taerobic\n"
    )
    edges.write_text(
        (fixture / "merged-kg_edges.tsv").read_text()
        + "anaerobe\tNCBITaxon:300\tbiolink:has_phenotype\tMETPO:1000603\tfixture:annotation\n"
        + "aerobe\tNCBITaxon:300\tbiolink:has_phenotype\tMETPO:1000602\tfixture:annotation\n"
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
        validator=validate_execution_product,
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
        _matrix: guild_maps.GuildMatrix, path: Path, _vintage: str | None, **_options: object
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
