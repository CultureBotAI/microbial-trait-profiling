"""Output isolation and replay boundary regressions for the split sidecar (#854)."""

import csv
import hashlib
import io
import json
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from mtp.cellulose_decr import selection_scan as scanner
from mtp.cellulose_decr import study_split as split_module
from mtp.cellulose_decr import study_split_run as runner
from mtp.cellulose_decr.study_split import Protocol
from mtp.multiplicity import sources
from mtp.multiplicity import stage as multiplicity_stage
from mtp.paths import file_sha256


def test_protocol_v2_preserves_and_binds_v1() -> None:
    historical = runner.PROTOCOL_FILE.with_name("study_split_selection_protocol_v1.json")
    assert historical.is_file(), "The original study-split v1 declaration must remain available"
    old_bytes = historical.read_bytes()
    old = json.loads(old_bytes)
    assert old["protocol"] == "study-split-selection-v1"
    assert old["primary_strain_rollup"] == "off"
    assert old["sensitivity_strain_rollup"] == "phenotype"
    new = json.loads(runner.PROTOCOL_FILE.read_bytes())
    assert new["protocol"] == "study-split-selection-v2"
    assert new["amends"] == {
        "path": "methods/study_split_selection_protocol_v1.json",
        "sha256": hashlib.sha256(old_bytes).hexdigest(),
        "protocol": "study-split-selection-v1",
    }
    assert new["primary_strain_rollup"] == "phenotype"
    assert new["sensitivity_strain_rollups"] == ["off"]
    assert new["absent_means"] == "off"
    assert new["scope_reporting"] == "side_by_side"
    assert new["known_results"]["declared_after_results"] == ["off", "phenotype"]
    assert new["known_results"]["source"] == (
        "methods/cellulose_decr_selection_scan_protocol_v2.json#/amendments/0/known_when_amended"
    )
    assert "historical corpus" in new["known_results"]["limitations"].lower()
    assert "methods/study_split_selection_protocol_v1.json" in runner.IMPLEMENTATIONS


@pytest.mark.parametrize("mutation", ["bytes", "path", "version", "primary", "sensitivity"])
def test_protocol_amendment_refuses_unbound_history(tmp_path: Path, mutation: str) -> None:
    load = getattr(split_module, "_declaration", None)
    assert callable(load), "The execution must validate its versioned protocol amendment"
    methods = tmp_path / "methods"
    methods.mkdir()
    declaration = json.loads(runner.PROTOCOL_FILE.read_bytes())
    historical = runner.PROTOCOL_FILE.with_name("study_split_selection_protocol_v1.json")
    target = methods / historical.name
    target.write_bytes(historical.read_bytes())
    if mutation == "bytes":
        target.write_bytes(target.read_bytes() + b"\n")
    elif mutation == "path":
        declaration["amends"]["path"] = "methods/unrelated.json"
    elif mutation == "version":
        declaration["amends"]["protocol"] = "study-split-selection-v0"
    elif mutation == "primary":
        declaration["primary_strain_rollup"] = "off"
    else:
        declaration["sensitivity_strain_rollups"] = ["phenotype"]
    current = methods / runner.PROTOCOL_FILE.name
    current.write_text(json.dumps(declaration), encoding="utf-8")
    refused = False
    try:
        load(current)
    except ValueError:
        refused = True
    assert refused, f"Protocol accepted an invalid amendment: {mutation}"


def _verified_scope(scope: str) -> Any:
    core = SimpleNamespace(
        run_id=f"core-{scope}",
        contrast_population=scanner.TIER1_2,
        direction_basis="disease_vs_comparison",
        trait_key_policy="fixture-policy",
        definition_sha256=hashlib.sha256(b"core-definition").hexdigest(),
        code_sha256=hashlib.sha256(b"core-code").hexdigest(),
        lock_sha256=hashlib.sha256(b"lock").hexdigest(),
        runtime_json='{"python": "fixture"}',
        git_commit="fixture-commit",
        git_dirty=False,
        parameters=[SimpleNamespace(name="strain_rollup", value=scope)],
        artifacts=[
            SimpleNamespace(product_id=name, sha256=hashlib.sha256(name.encode()).hexdigest())
            for name in scanner.CORE_SOURCES
        ],
    )
    return SimpleNamespace(core=core, scope=scope)


@pytest.mark.parametrize("scope", ["phenotype", "off"])
def test_scope_binding_uses_verified_sources_and_declared_roles(scope: str) -> None:
    bind = getattr(runner, "_scope_binding", None)
    assert callable(bind), "Study-split runs need a bound scope-comparison contract"
    value = bind(_verified_scope(scope))
    assert value["strain_rollup"] == scope
    assert value["role"] == (
        "primary_phenotype" if scope == "phenotype" else "declared_sensitivity"
    )
    assert value["comparison_contract"]["parameters"] == {}
    assert set(value["comparison_contract"]["sources"]) == set(scanner.CORE_SOURCES)


def test_absent_core_scope_means_historical_off_only() -> None:
    bind = getattr(runner, "_scope_binding", None)
    assert callable(bind)
    historical = _verified_scope("off")
    historical.core.parameters = []
    assert bind(historical)["role"] == "declared_sensitivity"
    historical.scope = "phenotype"
    refused = False
    try:
        bind(historical)
    except ValueError as error:
        refused = "contradicts" in str(error)
    assert refused, "A legacy core without a scope must not be silently reinterpreted as phenotype"


def _saved_scope(tmp_path: Path, scope: str) -> Path:
    """Synthetic empty-selection run, with independent hash-bound saved bytes."""
    directory = tmp_path / scope
    directory.mkdir()
    record = {
        "split": {"index": 0, "training": ["s1"], "testing": ["s2"]},
        "selection_status": "no_traits_selected",
        "selected_candidates": [],
        "correction": {
            "family": {
                "family_id": "study-split-selection-v2:split-000",
                "strain_rollup": scope,
                "profile_unit": "count_profiles",
                "direction_basis": "disease_vs_comparison",
            },
            "results": [],
        },
        "test_details": [],
    }
    documents = {
        "split_results.json": json.dumps([record]) + "\n",
        "split_stability.tsv": "\t".join(runner.STABILITY_FIELDS) + "\n",
        "README.md": f"Synthetic {scope} empty selection fixture.\n",
    }
    for name, contents in documents.items():
        (directory / name).write_text(contents, encoding="utf-8")
    manifest = {
        "protocol_version": "study-split-selection-v2",
        "provenance": {
            "protocol": json.loads(runner.PROTOCOL_FILE.read_bytes()),
            "parameters": asdict(Protocol(splits=1, permutations=3)),
            "inputs": {"comparison:manifest": hashlib.sha256(scope.encode()).hexdigest()},
            "implementation": {"fixture": hashlib.sha256(b"same-code").hexdigest()},
            "code_sha256": hashlib.sha256(b"same-transitive-code").hexdigest(),
            "environment": {"fixture": "same-environment"},
        },
        "scope_binding": {
            "core_run_id": f"core-{scope}",
            "strain_rollup": scope,
            "role": "primary_phenotype" if scope == "phenotype" else "declared_sensitivity",
            "comparison_contract": {
                "population": scanner.TIER1_2,
                "direction_basis": "disease_vs_comparison",
                "trait_key_policy": "fixture-policy",
                "producer": {
                    "definition_sha256": hashlib.sha256(b"core-definition").hexdigest(),
                    "code_sha256": hashlib.sha256(b"core-code").hexdigest(),
                    "lock_sha256": hashlib.sha256(b"lock").hexdigest(),
                    "runtime_sha256": hashlib.sha256(b"runtime").hexdigest(),
                    "git_commit": "fixture-commit",
                    "git_dirty": False,
                },
                "parameters": {"contrast_population": scanner.TIER1_2},
                "sources": {
                    name: hashlib.sha256(name.encode()).hexdigest() for name in scanner.CORE_SOURCES
                },
            },
        },
        "output_sha256": {
            name: hashlib.sha256(content.encode()).hexdigest()
            for name, content in documents.items()
        },
    }
    (directory / "study_split_run.json").write_text(json.dumps(manifest), encoding="utf-8")
    return directory


def test_scope_report_keeps_both_roles_and_replays_without_nulls(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    report = getattr(runner, "report_scopes", None)
    assert callable(report), "Both saved scope runs must be reported side by side"
    phenotype = _saved_scope(tmp_path, "phenotype")
    off = _saved_scope(tmp_path, "off")

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("A scope report must not resample or execute a split run")

    monkeypatch.setattr(runner, "run_splits", forbidden)
    monkeypatch.setattr(runner, "execute", forbidden)
    output = tmp_path / "paired"
    manifest_path = report(phenotype, off, output)
    assert manifest_path == output / "study_split_scopes.json"
    before = {path.name: path.read_bytes() for path in output.iterdir()}
    saved = json.loads(before["study_split_scopes.json"])
    assert saved["scope_order"] == ["phenotype", "off"]
    for scope, directory in (("phenotype", phenotype), ("off", off)):
        assert saved["scopes"][scope]["manifest_sha256"] == file_sha256(
            directory / "study_split_run.json"
        )
        assert saved["scopes"][scope]["splits"] == 1
    assert saved["scopes"]["phenotype"]["role"] == "primary_phenotype"
    assert saved["scopes"]["off"]["role"] == "declared_sensitivity"
    assert "dependent" in before["README.md"].decode()
    assert report(phenotype, off, output, verify=True) == manifest_path
    assert before == {path.name: path.read_bytes() for path in output.iterdir()}


@pytest.mark.parametrize(
    "mutation",
    [
        "hash",
        "role",
        "population",
        "source",
        "producer",
        "parameters",
        "partition",
        "family",
        "stability",
        "same_core",
        "missing_source",
        "core_parameters",
        "code",
    ],
)
def test_scope_report_refuses_noncomparable_or_inconsistent_runs(
    tmp_path: Path, mutation: str
) -> None:
    report = getattr(runner, "report_scopes", None)
    assert callable(report), "Scope pairing must check comparable hash-bound records"
    phenotype = _saved_scope(tmp_path, "phenotype")
    off = _saved_scope(tmp_path, "off")
    path = off / "study_split_run.json"
    manifest = json.loads(path.read_bytes())
    if mutation == "hash":
        (off / "README.md").write_text("changed\n", encoding="utf-8")
    elif mutation == "role":
        manifest["scope_binding"]["role"] = "primary_phenotype"
    elif mutation == "population":
        manifest["scope_binding"]["comparison_contract"]["population"] = "all_paired"
    elif mutation == "source":
        manifest["scope_binding"]["comparison_contract"]["sources"]["raw_dump"] = hashlib.sha256(
            b"different-corpus"
        ).hexdigest()
    elif mutation == "producer":
        manifest["scope_binding"]["comparison_contract"]["producer"]["code_sha256"] = (
            hashlib.sha256(b"different-producer").hexdigest()
        )
    elif mutation == "same_core":
        manifest["scope_binding"]["core_run_id"] = "core-phenotype"
    elif mutation == "missing_source":
        del manifest["scope_binding"]["comparison_contract"]["sources"]["raw_dump"]
    elif mutation == "core_parameters":
        manifest["scope_binding"]["comparison_contract"]["parameters"]["other_filter"] = "yes"
    elif mutation == "parameters":
        manifest["provenance"]["parameters"]["seed"] += 1
    elif mutation == "code":
        manifest["provenance"]["code_sha256"] = hashlib.sha256(b"different-code").hexdigest()
    else:
        name = "split_stability.tsv" if mutation == "stability" else "split_results.json"
        if mutation == "stability":
            contents = (off / name).read_text() + "invented\trow\n"
        else:
            records = json.loads((off / name).read_bytes())
            if mutation == "partition":
                records[0]["split"]["testing"] = ["s3"]
            else:
                records[0]["correction"]["family"]["strain_rollup"] = "phenotype"
            contents = json.dumps(records)
        (off / name).write_text(contents, encoding="utf-8")
        manifest["output_sha256"][name] = hashlib.sha256(contents.encode()).hexdigest()
    path.write_text(json.dumps(manifest), encoding="utf-8")
    refused = False
    try:
        report(phenotype, off, tmp_path / "paired")
    except ValueError:
        refused = True
    assert refused, f"Pairing accepted invalid saved evidence: {mutation}"
    assert not (tmp_path / "paired").exists()


def test_scope_report_preserves_nonempty_stability_in_primary_first_order(tmp_path: Path) -> None:
    report = getattr(runner, "report_scopes", None)
    assert callable(report)
    runs = {scope: _saved_scope(tmp_path, scope) for scope in ("phenotype", "off")}
    for scope, directory in runs.items():
        records = json.loads((directory / "split_results.json").read_bytes())
        record = records[0]
        record["selected_candidates"] = [{"trait": f"trait-{scope}"}]
        record["test_details"] = [
            {
                "hypothesis_id": "fixture-hypothesis",
                "trait_id": f"trait-{scope}",
                "conditions": ["A", "B", "C"],
                "statistic": "three_way_core",
            }
        ]
        record["correction"]["results"] = [
            {
                "hypothesis_id": "fixture-hypothesis",
                "status": "tested",
                "reject_by_005": scope == "phenotype",
                "reject_bh_005": True,
            }
        ]
        result_bytes = json.dumps(records).encode()
        row = {
            "hypothesis_id": "fixture-hypothesis",
            "trait_id": f"trait-{scope}",
            "conditions": json.dumps(["A", "B", "C"]),
            "statistic": "three_way_core",
            "total_splits": 1,
            "selected_splits": 1,
            "testable_splits": 1,
            "by_rejections_005": int(scope == "phenotype"),
            "bh_rejections_005": 1,
            "interpretation": "descriptive counts; overlapping splits are dependent",
        }
        output = io.StringIO(newline="")
        writer = csv.DictWriter(
            output, fieldnames=runner.STABILITY_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)
        manifest_path = directory / "study_split_run.json"
        manifest = json.loads(manifest_path.read_bytes())
        for name, raw in (
            ("split_results.json", result_bytes),
            ("split_stability.tsv", output.getvalue().encode()),
        ):
            (directory / name).write_bytes(raw)
            manifest["output_sha256"][name] = hashlib.sha256(raw).hexdigest()
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    output_dir = tmp_path / "paired"
    report(runs["phenotype"], runs["off"], output_dir)
    with (output_dir / "study_split_scopes.tsv").open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert [row["strain_rollup"] for row in rows] == ["phenotype", "off"]
    assert [row["trait_id"] for row in rows] == ["trait-phenotype", "trait-off"]
    assert [row["by_rejections_005"] for row in rows] == ["1", "0"]
    refused = False
    try:
        report(runs["phenotype"], runs["off"], output_dir)
    except ValueError as error:
        refused = "new or empty" in str(error)
    assert refused, "The scope report must not overwrite a previous report"


@pytest.mark.parametrize("field", ["trait_key_policy", "lock_sha256"])
def test_scope_report_requires_complete_metadata_in_both_runs(tmp_path: Path, field: str) -> None:
    report = getattr(runner, "report_scopes", None)
    assert callable(report)
    runs = {scope: _saved_scope(tmp_path, scope) for scope in ("phenotype", "off")}
    for directory in runs.values():
        path = directory / "study_split_run.json"
        manifest = json.loads(path.read_bytes())
        contract = manifest["scope_binding"]["comparison_contract"]
        target = contract if field == "trait_key_policy" else contract["producer"]
        del target[field]
        path.write_text(json.dumps(manifest), encoding="utf-8")
    refused = False
    try:
        report(runs["phenotype"], runs["off"], tmp_path / "paired")
    except ValueError:
        refused = True
    assert refused, f"Deleting {field} from both runs must not make them comparable"
    assert not (tmp_path / "paired").exists()


def test_scope_report_cli_routes_both_scopes_without_running_splits(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    report = getattr(runner, "report_scopes", None)
    assert callable(report)
    calls: list[tuple[Path, Path, Path, bool]] = []

    def capture(phenotype: Path, off: Path, output: Path, *, verify: bool = False) -> Path:
        calls.append((phenotype, off, output, verify))
        return output / "study_split_scopes.json"

    monkeypatch.setattr(runner, "report_scopes", capture)
    paths = [tmp_path / name for name in ("phenotype", "off", "paired")]
    assert (
        runner.main(
            [
                "--phenotype-run",
                str(paths[0]),
                "--off-run",
                str(paths[1]),
                "--output-dir",
                str(paths[2]),
                "--verify",
            ]
        )
        == 0
    )
    assert calls == [(*paths, True)]


def test_consumed_scanner_protocol_is_bound_in_provenance() -> None:
    # The scanner reads v2, which amends v1 by its hash (#1050); both are bound.
    assert "methods/cellulose_decr_selection_scan_protocol.json" in runner.IMPLEMENTATIONS
    assert "methods/cellulose_decr_selection_scan_protocol_v2.json" in runner.IMPLEMENTATIONS


def test_active_scientific_helper_changes_provenance(monkeypatch: pytest.MonkeyPatch) -> None:
    before = runner._provenance({}, Protocol())
    original_hash = file_sha256

    def changed_helper(path: Path) -> str:
        if path.name == "table1_trait_maps.py":
            return "f" * 64
        return original_hash(path)

    monkeypatch.setattr(multiplicity_stage, "file_sha256", changed_helper)
    assert runner._provenance({}, Protocol()) != before


@pytest.mark.parametrize("location", ["child", "same", "ancestor", "symlink"])
def test_outputs_cannot_modify_an_upstream_bundle(tmp_path: Path, location: str) -> None:
    source = tmp_path / "verified" / "arms" / "noel"
    source.mkdir(parents=True)
    manifest = source / "arm_run.json"
    manifest.write_text("{}\n")
    outputs = {
        "child": source / "sidecar",
        "same": source,
        "ancestor": source.parent,
        "symlink": tmp_path / "alias",
    }
    outputs["symlink"].symlink_to(source, target_is_directory=True)
    with pytest.raises(ValueError, match="immutable upstream bundle"):
        runner.preflight_bundle_roots(outputs[location], {"arm:noel:manifest": manifest})
    assert manifest.read_text() == "{}\n"
    assert not (source / "sidecar").exists()


def test_sibling_sidecar_is_allowed_without_writing(tmp_path: Path) -> None:
    output = tmp_path / "sidecars" / "splits"
    runner.preflight_bundle_roots(
        output,
        {
            "comparison:manifest": tmp_path / "comparison" / "comparison_run.json",
            "core:manifest": tmp_path / "core" / "run.json",
            "arm:noel:manifest": tmp_path / "arms" / "noel" / "arm_run.json",
        },
    )
    assert not output.exists()


def test_directory_input_uses_the_verified_manifest_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    comparison_dir = tmp_path / "comparison"
    comparison_dir.mkdir()
    manifest = comparison_dir / "comparison_run.json"
    manifest.write_text('{"routing_fixture": true}\n')
    original_digest = file_sha256(manifest)
    upstream = SimpleNamespace(
        core=SimpleNamespace(
            contrast_population="disease_vs_healthy_tier1_2_subject_controls",
            direction_basis="disease_vs_comparison",
        ),
        scope="off",
        comparison=SimpleNamespace(
            manifest_path=str(manifest),
            core_manifest_path=str(tmp_path / "core" / "run.json"),
            noel_manifest_path=str(tmp_path / "noel" / "arm_run.json"),
            jed_manifest_path=str(tmp_path / "jed" / "arm_run.json"),
            ilyas_manifest_path=str(tmp_path / "ilyas" / "arm_run.json"),
        ),
        product=lambda name: tmp_path / "core" / name,
    )

    observed: list[str | None] = []

    def open_then_change_comparison(
        _path: Path, *, expected_sha256: str | None = None
    ) -> SimpleNamespace:
        observed.append(expected_sha256)
        assert expected_sha256 == original_digest, "Existing comparison digest was not forwarded"
        manifest.write_text(manifest.read_text() + "\n", encoding="utf-8")
        assert file_sha256(manifest) != original_digest
        return upstream

    monkeypatch.setattr(sources, "open_verified_sources", open_then_change_comparison)
    bindings: list[str] = []

    def bind(verified: SimpleNamespace, comparison_sha256: str) -> SimpleNamespace:
        assert verified is upstream
        bindings.append(comparison_sha256)
        return verified

    monkeypatch.setattr(sources, "bind_recorded_hashes", bind)
    _, paths = runner._inputs(comparison_dir)
    assert observed == [original_digest]
    assert paths["comparison:manifest"] == manifest
    assert bindings == [original_digest]
    runner.preflight_bundle_roots(tmp_path / "sidecar", paths)


@pytest.mark.parametrize(
    ("population", "direction", "scope"),
    [
        ("all_paired", "group_1_vs_group_0", "off"),
        ("disease_vs_healthy_tier1_2_subject_controls", "group_1_vs_group_0", "off"),
        ("disease_vs_healthy_tier1_2_subject_controls", "disease_vs_comparison", "full"),
    ],
)
def test_source_boundary_refuses_other_populations_directions_and_scopes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    population: str,
    direction: str,
    scope: str,
) -> None:
    upstream = SimpleNamespace(
        core=SimpleNamespace(contrast_population=population, direction_basis=direction),
        scope=scope,
    )
    comparison = tmp_path / "comparison_run.json"
    comparison.write_text('{"invalid_population_fixture": true}\n')

    def verified(path: Path, *, expected_sha256: str | None = None) -> SimpleNamespace:
        assert path == comparison
        assert expected_sha256 == file_sha256(comparison)
        return upstream

    monkeypatch.setattr(sources, "open_verified_sources", verified)
    with pytest.raises(ValueError, match=r"Tier 1\+2|phenotype sensitivity"):
        runner.execute(comparison, tmp_path / "out", Protocol())
    assert not (tmp_path / "out").exists()


def test_runner_ignores_untracked_file_and_dangling_editor_lock(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import shutil

    from git import Actor, Repo

    from mtp.paths import REPO_ROOT

    checkout = tmp_path / "checkout"
    for rel_path in {*runner.IMPLEMENTATIONS, "src/mtp/multiplicity/stage.py"}:
        src_file = REPO_ROOT / rel_path
        dst_file = checkout / rel_path
        dst_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src_file, dst_file)

    actor = Actor("Test", "test@example.com")
    with Repo.init(checkout) as repository:
        repository.index.add([*runner.IMPLEMENTATIONS, "src/mtp/multiplicity/stage.py"])
        repository.index.commit("initial tracked tree", author=actor, committer=actor)

    stage_path = checkout / "src/mtp/multiplicity/stage.py"
    monkeypatch.setattr(multiplicity_stage, "__file__", str(stage_path))
    monkeypatch.setattr(runner, "REPO_ROOT", checkout)
    monkeypatch.setattr(
        runner, "PROTOCOL_FILE", checkout / "methods/study_split_selection_protocol.json"
    )
    monkeypatch.setattr(
        split_module, "PROTOCOL_PATH", checkout / "methods/study_split_selection_protocol.json"
    )

    upstream = tmp_path / "upstream"
    upstream.mkdir()
    comparison = upstream / "comparison" / "comparison_run.json"
    comparison.parent.mkdir()
    comparison.write_text("{}\n", encoding="utf-8")
    core_manifest = upstream / "core" / "run.json"
    core_manifest.parent.mkdir()
    core_manifest.write_text("{}\n", encoding="utf-8")
    verified = _verified_scope("off")
    verified.comparison = SimpleNamespace(id="comparison-fixture")
    paths = {"core:manifest": core_manifest, "comparison:manifest": comparison}
    digests = {path: file_sha256(path) for path in paths.values()}
    verified.expected_digest = digests.get
    monkeypatch.setattr(runner, "_inputs", lambda _comparison: (verified, paths))
    monkeypatch.setattr(runner, "open_core_bundle", lambda _manifest: SimpleNamespace())
    monkeypatch.setattr(runner, "preflight_consumer_outputs", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        runner,
        "_documents",
        lambda _verified, _paths, _protocol, _ts: {
            "split_results.json": "[]\n",
            "split_stability.tsv": "hypothesis_id\n",
            "README.md": "# Study split fixture\n",
        },
    )

    clean_out = tmp_path / "sidecars" / "clean"
    manifest_path = runner.execute(comparison, clean_out, Protocol(), verify=False)
    clean_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert clean_manifest["repository_context"]["code_dirty"] is False

    package_dir = checkout / "src" / "mtp"
    (package_dir / "untracked_scratch.py").write_text("# unrelated scratch\n", encoding="utf-8")
    (package_dir / ".#study_split_run.py").symlink_to(package_dir / "missing.py")

    assert runner.execute(comparison, clean_out, Protocol(), verify=True) == manifest_path
    dirty_out = tmp_path / "sidecars" / "with_lock"
    dirty_manifest_path = runner.execute(comparison, dirty_out, Protocol(), verify=False)
    dirty_manifest = json.loads(dirty_manifest_path.read_text(encoding="utf-8"))
    assert (
        dirty_manifest["provenance"]["code_sha256"] == clean_manifest["provenance"]["code_sha256"]
    )
    assert (
        dirty_manifest["provenance"]["implementation"]
        == clean_manifest["provenance"]["implementation"]
    )
    assert dirty_manifest["repository_context"]["code_dirty"] is True


@pytest.mark.parametrize(
    ("protocol", "expected_follows", "expected_seed"),
    [
        (Protocol(), "follows_protocol=True", "seed=20261023"),
        (Protocol(splits=2, seed=99, permutations=19), "follows_protocol=False", "seed=99"),
    ],
)
def test_readme_states_follows_protocol_and_seed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    protocol: Protocol,
    expected_follows: str,
    expected_seed: str,
) -> None:
    from datetime import UTC, datetime

    verified: Any = SimpleNamespace(
        scope="phenotype",
        features={"count_profiles": ("trait",)},
        bases={"count_profiles": SimpleNamespace()},
        comparison=SimpleNamespace(id="comparison-fixture"),
    )
    paths = {
        "core:taxa_matrix": tmp_path / "taxa_matrix.tsv",
        "core:experiments": tmp_path / "experiments.tsv",
        "core:reported_taxa": tmp_path / "reported_taxa.tsv",
        "core:count_profiles": tmp_path / "count_profiles.tsv",
    }
    monkeypatch.setattr(scanner, "read_matrix", lambda *_args: SimpleNamespace(index=[]))
    monkeypatch.setattr(scanner, "load_population", lambda *_args: (None, None))
    monkeypatch.setattr(runner, "load_profiles", lambda *_args, **_kwargs: ([], None))
    monkeypatch.setattr(runner, "_recorded_inputs", lambda *_args: {})
    monkeypatch.setattr(runner, "run_splits", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(runner, "summarize_stability", lambda _records: [])

    documents = runner._documents(verified, paths, protocol, datetime(2026, 11, 13, tzinfo=UTC))
    readme = documents["README.md"]
    assert expected_follows in readme
    assert expected_seed in readme


def test_cli_defaults_match_protocol_defaults(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from dataclasses import dataclass

    @dataclass(frozen=True)
    class CustomDefaultsProtocol:
        splits: int = 3
        seed: int = 777
        permutations: int = 17

    recorded: list[Any] = []

    def fake_execute(
        _comparison: Path, output_dir: Path, protocol: Any, *, verify: bool = False
    ) -> Path:
        recorded.append((protocol, verify))
        return output_dir / "study_split_run.json"

    monkeypatch.setattr(runner, "Protocol", CustomDefaultsProtocol)
    monkeypatch.setattr(runner, "execute", fake_execute)

    comparison = tmp_path / "comparison_run.json"
    output = tmp_path / "out"
    assert (
        runner.main(["--comparison-run", str(comparison), "--output-dir", str(output), "--verify"])
        == 0
    )
    assert recorded == [(CustomDefaultsProtocol(splits=3, seed=777, permutations=17), True)]

    with pytest.raises(SystemExit):
        runner.main(
            [
                "--phenotype-run",
                str(tmp_path / "phenotype"),
                "--off-run",
                str(tmp_path / "off"),
                "--output-dir",
                str(tmp_path / "paired"),
                "--splits",
                "20",
            ]
        )
