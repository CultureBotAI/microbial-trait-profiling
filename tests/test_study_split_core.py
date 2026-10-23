"""Real synthetic core/basis binding and sidecar replay for #854.

The compact replay cases stub comparison verification while exercising real core
execution. The integration case also executes the actual consumed arm stages and
comparison, with the comparison fixture's explicit synthetic influence assessment.
"""

from __future__ import annotations

import csv
import json
import random
import subprocess
from pathlib import Path

import pytest

from mtp.cellulose_decr import selection_scan as scanner
from mtp.cellulose_decr import study_split_run as runner
from mtp.cellulose_decr.study_split import Protocol
from mtp.datamodel.workflow_pydantic import ComparisonRun
from mtp.multiplicity.sources import VerifiedSources, bind_analysis_basis
from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow.comparison import core_analysis_bases
from mtp.workflow.consumers import open_core_bundle
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow


def sources(root: Path) -> dict[str, Path]:
    fixture = REPO_ROOT / "tests/fixtures/workflow"
    root.mkdir()
    paths = {name: root / name for name in ("raw_dump", "kg_edges", "kg_nodes")}
    paths["kg_edges"].write_text(
        (fixture / "merged-kg_edges.tsv").read_text()
        + "shape\tNCBITaxon:102\tbiolink:has_phenotype\tMETPO:1000678\tfixture\n"
        + "strain\tkgmicrobe.strain:1\tbiolink:subclass_of\tNCBITaxon:300\tfixture\n"
        + "strainshape\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1000682\tfixture\n"
    )
    paths["kg_nodes"].write_text(
        (fixture / "merged-kg_nodes.tsv").read_text()
        + "METPO:1000678\tbiolink:PhenotypicQuality\toval shaped\n"
        + "METPO:1000682\tbiolink:PhenotypicQuality\tfixture phenotype\n"
    )
    with (fixture / "contrast_population.csv").open() as handle:
        template = next(csv.DictReader(line for line in handle if not line.startswith("#")))
    with paths["raw_dump"].open("w", newline="") as handle:
        handle.write("# BugSigDB fixture-854, License: CC0 synthetic test data\n")
        writer = csv.DictWriter(handle, fieldnames=list(template))
        writer.writeheader()
        for index in range(18):
            condition = ("Parkinson's disease", "HIV infection", "colorectal cancer")[index % 3]
            for number, direction in enumerate(("increased", "decreased"), 1):
                reverse = index == 0
                disease_higher = (direction == "increased") != reverse
                row = dict(template)
                row.update(
                    {
                        "BSDB ID": f"bsdb:{index + 1}/1/{number}",
                        "Study": f"Study {index + 1}",
                        "Condition": condition,
                        "Body site": "Feces",
                        "Group 0 name": condition if reverse else "healthy controls",
                        "Group 1 name": "healthy controls" if reverse else condition,
                        "Group 1 definition": "",
                        "Abundance in Group 1": direction,
                        "NCBI Taxonomy IDs": "102;300" if disease_higher else "101;100",
                        "MetaPhlAn taxon names": "s__Examplegenus_beta,s__Unannotated_bacterium"
                        if disease_higher
                        else "s__Examplegenus_alpha,g__Examplegenus",
                    }
                )
                writer.writerow(row)
    return paths


@pytest.mark.parametrize("scope", ["off", "phenotype"])
@pytest.mark.parametrize(
    "tampering",
    ["contents", "line_endings", "extra_file", "host_newlines", "helper_source", "git_context"],
)
def test_real_scoped_core_binds_features_and_replays_sidecar(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scope: str, tampering: str
) -> None:
    core = execute_workflow(
        REPO_ROOT / "workflows/bugsigdb_core_subset.yaml",
        sources(tmp_path / "inputs"),
        tmp_path / "core",
        OPERATIONS,
        parameter_overrides={
            "source_mode": "fixture",
            "strain_rollup": scope,
            "contrast_population": scanner.TIER1_2,
        },
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    core_path = Path(core.manifest_path)
    bundle = open_core_bundle(core_path)
    bases = {basis.profile_unit: basis for basis in core_analysis_bases(core)}
    _, features = bind_analysis_basis(bundle, bases["count_profiles"])
    assert ("METPO:1000682" in features) == (scope == "phenotype")
    comparison_path = tmp_path / "comparison" / "comparison_run.json"
    comparison_path.parent.mkdir()
    comparison_path.write_text('{"fixture_boundary": "comparison verifier stub"}\n')
    comparison = ComparisonRun.model_construct(id="fixture-comparison")
    verified = VerifiedSources(comparison, core, {}, bases, {"count_profiles": features}, scope)
    paths = {f"core:{name}": bundle.path(name) for name in runner.PRODUCTS}
    paths.update({"comparison:manifest": comparison_path, "core:manifest": core_path})
    monkeypatch.setattr(runner, "_inputs", lambda _path: (verified, paths))
    originals = {name: file_sha256(path) for name, path in paths.items()}
    protocol = Protocol(splits=2, min_studies=1, permutations=19)
    output = tmp_path / "split-sidecar"
    if tampering == "host_newlines":
        original_write_text = Path.write_text

        def translated_write_text(
            path: Path,
            data: str,
            encoding: str | None = None,
            errors: str | None = None,
            newline: str | None = None,
        ) -> int:
            if path.is_relative_to(output):
                data = data.replace("\n", "\r\n")
            return original_write_text(
                path, data, encoding=encoding, errors=errors, newline=newline
            )

        monkeypatch.setattr(Path, "write_text", translated_write_text)
    manifest = runner.execute(comparison_path, output, protocol)
    records = json.loads((output / "split_results.json").read_text())
    assert len(records) == 2
    assert all(record["selected_candidates"] for record in records)
    for record in records:
        assert len(record["split"]["training"]) == len(record["split"]["testing"]) == 9
        assert not set(record["split"]["training"]) & set(record["split"]["testing"])
        assert record["correction"]["family"]["retained_feature_count"] == len(features)
        assert record["correction"]["family"]["strain_rollup"] == scope
        assert record["correction"]["family"]["feature_filter_id"] == (
            bases["count_profiles"].feature_filter_id
        )
    assert runner.execute(comparison_path, output, protocol, verify=True) == manifest
    assert originals == {name: file_sha256(path) for name, path in paths.items()}
    if tampering == "git_context":
        saved = manifest.read_bytes()
        monkeypatch.setattr(runner, "repository_state", lambda _root: ("e" * 40, False))
        assert runner.execute(comparison_path, output, protocol, verify=True) == manifest
        monkeypatch.setattr(runner, "repository_state", lambda _root: ("d" * 40, True))
        assert runner.execute(comparison_path, output, protocol, verify=True) == manifest
        assert manifest.read_bytes() == saved
        return
    if tampering == "host_newlines":
        assert all(b"\r\n" not in path.read_bytes() for path in output.iterdir())
        return
    if tampering == "helper_source":
        # A comment-only helper change leaves numerical results untouched, but
        # its source identity must still invalidate replay (#902).
        original_hash = file_sha256

        def changed_helper(path: Path) -> str:
            return "f" * 64 if path.name == "table1_trait_maps.py" else original_hash(path)

        monkeypatch.setattr(runner, "file_sha256", changed_helper)
        with pytest.raises(ValueError, match="sources, protocol, code or environment changed"):
            runner.execute(comparison_path, output, protocol, verify=True)
        return
    if tampering == "contents":
        (output / "split_stability.tsv").write_text("corrupt\n")
    elif tampering == "line_endings":
        table = output / "split_stability.tsv"
        table.write_bytes(table.read_bytes().replace(b"\n", b"\r\n"))
    else:
        (output / "stale_results.tsv").write_text("unrecorded\n")
    with pytest.raises(ValueError, match=r"replay differs|inventory"):
        runner.execute(comparison_path, output, protocol, verify=True)


@pytest.mark.parametrize("scope", ["off", "phenotype"])
def test_real_arm_comparison_and_split_source_verification(tmp_path: Path, scope: str) -> None:
    """Exercise the public verified source boundary without substituting its verifier."""
    from mtp.workflow.arms.runtime import execute_arm
    from mtp.workflow.comparison import execute_comparison_run, verify_comparison_run
    from test_comparison_run import _consumed_plan, _write_inputs

    inputs = _write_inputs(tmp_path / "inputs")
    with inputs["kg_edges"].open(newline="") as handle:
        edges = list(csv.DictReader(handle, delimiter="\t"))
    with inputs["kg_nodes"].open(newline="") as handle:
        names = {row["id"]: row["name"] for row in csv.DictReader(handle, delimiter="\t")}
    annotated = sorted(
        {row["subject"] for row in edges if row["subject"].startswith("NCBITaxon:10")}
        - {"NCBITaxon:1999"}
    )
    cellulose = sorted({row["subject"] for row in edges if row["object"] == "GO:0030245"})
    other = sorted(set(annotated) - set(cellulose))
    assert len(cellulose) >= 3
    assert len(other) >= 3
    with inputs["raw_dump"].open(newline="") as handle:
        comment = handle.readline()
        template = next(csv.DictReader(handle))
    rng = random.Random(854)
    with inputs["raw_dump"].open("w", newline="") as handle:
        handle.write(comment)
        writer = csv.DictWriter(handle, fieldnames=list(template))
        writer.writeheader()
        for study in range(1, 19):
            condition = ("Parkinson's disease", "HIV infection", "colorectal cancer")[study % 3]
            for signature, direction in ((1, "increased"), (2, "decreased")):
                # Keep the decreasing signal below the native 95% support ceiling.
                candidates = other if direction == "increased" or study == 1 else cellulose
                chosen = sorted(rng.sample(candidates, rng.randint(3, min(6, len(candidates)))))
                row = dict(template)
                row.update(
                    {
                        "BSDB ID": f"bsdb:{study}/1/{signature}",
                        "Study": f"Study {study}",
                        "Condition": condition,
                        "Body site": "Feces",
                        "Abundance in Group 1": direction,
                        "NCBI Taxonomy IDs": ";".join(
                            taxon.removeprefix("NCBITaxon:") for taxon in chosen
                        ),
                        "MetaPhlAn taxon names": ",".join(
                            "s__" + names[taxon].replace(" ", "_") for taxon in chosen
                        ),
                    }
                )
                writer.writerow(row)
    core = execute_workflow(
        REPO_ROOT / "workflows/bugsigdb_core_subset.yaml",
        inputs,
        tmp_path / "core",
        OPERATIONS,
        parameter_overrides={
            "source_mode": "fixture",
            "strain_rollup": scope,
            "contrast_population": scanner.TIER1_2,
        },
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    core_path = Path(core.manifest_path)
    arms = {
        arm: Path(
            execute_arm(
                arm,
                core_path,
                tmp_path / arm,
                _consumed_plan(arm),
                replicates=2,
                seed=0,
                figures=False,
            ).manifest_path
        )
        for arm in ("noel", "jed", "ilyas")
    }
    comparison = execute_comparison_run(
        core_path, arms["noel"], arms["jed"], arms["ilyas"], tmp_path / "comparison"
    )
    comparison_path = Path(comparison.manifest_path)
    assert verify_comparison_run(comparison_path).id == comparison.id
    protocol = Protocol(splits=2, min_studies=1, permutations=19)
    output = tmp_path / "split-sidecar"
    manifest = runner.execute(comparison_path, output, protocol)
    records = json.loads((output / "split_results.json").read_text())
    assert all(record["selected_candidates"] for record in records)
    assert all(record["correction"]["family"]["strain_rollup"] == scope for record in records)
    assert runner.execute(comparison_path, output, protocol, verify=True) == manifest
    cli_output = tmp_path / "cli split sidecar"
    command = [
        "just",
        "selection-study-splits",
        str(comparison_path.parent),
        str(cli_output),
        "--splits",
        "2",
        "--permutations",
        "19",
    ]
    subprocess.run(command, cwd=REPO_ROOT, check=True, capture_output=True, text=True)
    recorded = json.loads((cli_output / "study_split_run.json").read_text())
    assert recorded["provenance"]["parameters"]["splits"] == 2
    assert recorded["provenance"]["parameters"]["permutations"] == 19
    subprocess.run(
        [*command, "--verify"], cwd=REPO_ROOT, check=True, capture_output=True, text=True
    )
