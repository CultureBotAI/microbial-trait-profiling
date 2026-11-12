"""Pre-publication provenance recheck and untestable-family reason tests (#997)."""

from __future__ import annotations

import hashlib
from collections.abc import Iterator
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import pandas as pd
import pytest

from mtp.cellulose_decr import selection_scan as scanner
from mtp.cellulose_decr import study_split as split_module
from mtp.cellulose_decr import study_split_run as runner
from mtp.cellulose_decr.study_split import (
    Candidate,
    Protocol,
    StudySplit,
    _training_selection,
    run_splits,
    summarize_stability,
)
from mtp.datamodel.multiple_testing_pydantic import SourceArtifact
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import Profile
from mtp.workflow.coverage import feature_filter_id
from mtp.workflow.downstream_contracts import build_analysis_basis


def _synthetic_corpus(
    *, shared_conditions: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[Profile]]:
    profile_rows = []
    taxa_rows = []
    counts = []
    for index in range(18):
        study = f"s{index:02}"
        condition = ("A", "B", "C")[index % 3]
        conditions = ("A", "B", "C") if shared_conditions else (condition, condition)
        for experiment_index, cond in enumerate(conditions, start=1):
            experiment = f"e{experiment_index}"
            for direction in ("inc", "dec"):
                key = f"{study}|{experiment}|{direction}"
                row = {
                    "study_id": study,
                    "experiment_id": experiment,
                    "condition": cond,
                    "body_site": "Feces",
                    "compartment": "Feces",
                    "dir": direction,
                }
                profile_rows.append(dict(profile_id=key, **row))
                taxa_rows.append(dict(profile_id=key, taxon_id="t1", **row))
                counts.append(
                    Profile(
                        study,
                        experiment,
                        direction,
                        cond,
                        "Feces",
                        {"trait": 1 if direction == "inc" else 3, "other": 3},
                        4 if direction == "inc" else 6,
                    )
                )
    return (
        pd.DataFrame(profile_rows).set_index("profile_id"),
        pd.DataFrame(taxa_rows),
        pd.DataFrame({"trait": [1, 0], "other": [0, 1]}, index=["t1", "t2"]),
        counts,
    )


def _setup_execute_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, dict[str, Path], dict[Path, str]]:
    upstream_root = tmp_path / "upstream"
    upstream_root.mkdir(parents=True)
    paths: dict[str, Path] = {}
    for product in runner.PRODUCTS:
        product_path = upstream_root / "core" / f"{product}.tsv"
        product_path.parent.mkdir(parents=True, exist_ok=True)
        product_path.write_text(f"fixture:{product}\n", encoding="utf-8")
        paths[f"core:{product}"] = product_path

    comparison_manifest = upstream_root / "comparison" / "comparison_run.json"
    comparison_manifest.parent.mkdir(parents=True, exist_ok=True)
    comparison_manifest.write_text('{"id": "comparison-fixture"}\n', encoding="utf-8")
    paths["comparison:manifest"] = comparison_manifest

    core_manifest = upstream_root / "core" / "run.json"
    core_manifest.write_text('{"id": "core-fixture"}\n', encoding="utf-8")
    paths["core:manifest"] = core_manifest

    for arm in ("jed", "noel", "ilyas"):
        arm_manifest = upstream_root / "arms" / arm / "arm_run.json"
        arm_manifest.parent.mkdir(parents=True, exist_ok=True)
        arm_manifest.write_text(f'{{"arm": "{arm}"}}\n', encoding="utf-8")
        paths[f"arm:{arm}:manifest"] = arm_manifest

    expected_digests: dict[Path, str] = {path: file_sha256(path) for path in paths.values()}
    verified = SimpleNamespace(
        scope="off",
        # Real VerifiedSources always carries the verified core used by scope binding.
        core=SimpleNamespace(
            run_id="core-fixture",
            contrast_population=scanner.TIER1_2,
            direction_basis="disease_vs_comparison",
            trait_key_policy="raw_object_v1",
            definition_sha256=hashlib.sha256(b"core-definition").hexdigest(),
            code_sha256=hashlib.sha256(b"core-code").hexdigest(),
            lock_sha256=hashlib.sha256(b"core-lock").hexdigest(),
            runtime_json='{"python": "fixture"}',
            git_commit="0adcceed",
            git_dirty=False,
            parameters=[SimpleNamespace(name="strain_rollup", value="off")],
            artifacts=[
                SimpleNamespace(product_id=name, sha256=hashlib.sha256(name.encode()).hexdigest())
                for name in scanner.CORE_SOURCES
            ],
        ),
        comparison=SimpleNamespace(id="comparison-fixture"),
        expected_digest=lambda path: expected_digests.get(path),
    )
    monkeypatch.setattr(runner, "_inputs", lambda _comparison: (verified, paths))
    monkeypatch.setattr(runner, "open_core_bundle", lambda _manifest: SimpleNamespace())
    monkeypatch.setattr(runner, "preflight_consumer_outputs", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(runner, "repository_state", lambda _root: ("0adcceed", False))
    monkeypatch.setattr(runner, "code_sha256", lambda: "a" * 64)
    return comparison_manifest, paths, expected_digests


def test_execute_publishes_and_verifies_when_provenance_stays_constant(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    comparison_manifest, _, _ = _setup_execute_fixture(tmp_path, monkeypatch)
    output_dir = tmp_path / "sidecars" / "splits"
    monkeypatch.setattr(
        runner,
        "_documents",
        lambda _verified, _paths, _protocol, _ts: {
            "split_results.json": "[]\n",
            "split_stability.tsv": "hypothesis_id\n",
            "README.md": "# Study split fixture\n",
        },
    )
    manifest_path = runner.execute(comparison_manifest, output_dir, Protocol(), verify=False)
    assert manifest_path == output_dir / "study_split_run.json"
    assert manifest_path.is_file()
    assert (output_dir / "split_results.json").is_file()
    assert (output_dir / "split_stability.tsv").is_file()
    assert (output_dir / "README.md").is_file()

    replayed = runner.execute(comparison_manifest, output_dir, Protocol(), verify=True)
    assert replayed == manifest_path


def test_execute_refuses_publication_when_code_sha256_drifts_during_execution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    comparison_manifest, _, _ = _setup_execute_fixture(tmp_path, monkeypatch)
    output_dir = tmp_path / "sidecars" / "splits_code_drift"

    def drift_code_during_documents(
        _verified: object, _paths: dict[str, Path], _protocol: Protocol, _ts: datetime
    ) -> dict[str, str]:
        monkeypatch.setattr(runner, "code_sha256", lambda: "b" * 64)
        return {
            "split_results.json": "[]\n",
            "split_stability.tsv": "hypothesis_id\n",
            "README.md": "# Study split fixture\n",
        }

    monkeypatch.setattr(runner, "_documents", drift_code_during_documents)
    caught: ValueError | None = None
    try:
        runner.execute(comparison_manifest, output_dir, Protocol(), verify=False)
    except ValueError as exc:
        caught = exc
    assert caught is not None
    assert str(caught) == "Study-split inputs or implementation changed during execution"
    assert not output_dir.exists()


def test_execute_refuses_publication_when_implementation_file_drifts_during_execution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    comparison_manifest, _, _ = _setup_execute_fixture(tmp_path, monkeypatch)
    output_dir = tmp_path / "sidecars" / "splits_impl_drift"
    real_file_sha256 = file_sha256

    def drift_impl_during_documents(
        _verified: object, _paths: dict[str, Path], _protocol: Protocol, _ts: datetime
    ) -> dict[str, str]:
        target_impl = (REPO_ROOT / "src/mtp/cellulose_decr/study_split.py").resolve()

        def mutated_file_sha256(path: Path | str) -> str:
            if Path(path).resolve() == target_impl:
                return "f" * 64
            return real_file_sha256(path)

        monkeypatch.setattr(runner, "file_sha256", mutated_file_sha256)
        return {
            "split_results.json": "[]\n",
            "split_stability.tsv": "hypothesis_id\n",
            "README.md": "# Study split fixture\n",
        }

    monkeypatch.setattr(runner, "_documents", drift_impl_during_documents)
    caught: ValueError | None = None
    try:
        runner.execute(comparison_manifest, output_dir, Protocol(), verify=False)
    except ValueError as exc:
        caught = exc
    assert caught is not None
    assert str(caught) == "Study-split inputs or implementation changed during execution"
    assert not output_dir.exists()


def test_execute_refuses_publication_when_consumed_input_drifts_during_execution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    comparison_manifest, paths, expected_digests = _setup_execute_fixture(tmp_path, monkeypatch)
    output_dir = tmp_path / "sidecars" / "splits_input_drift"

    def drift_input_during_documents(
        _verified: object, _paths: dict[str, Path], _protocol: Protocol, _ts: datetime
    ) -> dict[str, str]:
        target = paths["core:count_profiles"]
        target.write_text("mutated:count_profiles\n", encoding="utf-8")
        expected_digests[target] = file_sha256(target)
        return {
            "split_results.json": "[]\n",
            "split_stability.tsv": "hypothesis_id\n",
            "README.md": "# Study split fixture\n",
        }

    monkeypatch.setattr(runner, "_documents", drift_input_during_documents)
    caught: ValueError | None = None
    try:
        runner.execute(comparison_manifest, output_dir, Protocol(), verify=False)
    except ValueError as exc:
        caught = exc
    assert caught is not None
    assert str(caught) == "Study-split inputs or implementation changed during execution"
    assert not output_dir.exists()


def test_execute_verify_mode_refuses_mid_replay_provenance_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    comparison_manifest, _, _ = _setup_execute_fixture(tmp_path, monkeypatch)
    output_dir = tmp_path / "sidecars" / "splits_verify_drift"
    monkeypatch.setattr(
        runner,
        "_documents",
        lambda _verified, _paths, _protocol, _ts: {
            "split_results.json": "[]\n",
            "split_stability.tsv": "hypothesis_id\n",
            "README.md": "# Study split fixture\n",
        },
    )
    runner.execute(comparison_manifest, output_dir, Protocol(), verify=False)

    def drift_on_replay(
        _verified: object, _paths: dict[str, Path], _protocol: Protocol, _ts: datetime
    ) -> dict[str, str]:
        monkeypatch.setattr(runner, "code_sha256", lambda: "c" * 64)
        return {
            "split_results.json": "[]\n",
            "split_stability.tsv": "hypothesis_id\n",
            "README.md": "# Study split fixture\n",
        }

    monkeypatch.setattr(runner, "_documents", drift_on_replay)
    caught: ValueError | None = None
    try:
        runner.execute(comparison_manifest, output_dir, Protocol(), verify=True)
    except ValueError as exc:
        caught = exc
    assert caught is not None
    assert str(caught) == "Study-split inputs or implementation changed during execution"


@pytest.mark.parametrize(
    "case",
    [
        "no_traits",
        "no_traits_zero_counts",
        "no_testable_study_disjoint_triple",
        "fewer_than_three_eligible_conditions",
        "selected",
    ],
)
def test_training_selection_and_run_splits_untestable_reasons(case: str) -> None:
    profiles, taxa, matrix, counts = _synthetic_corpus(
        shared_conditions=case == "no_testable_study_disjoint_triple"
    )
    columns: tuple[str, ...] = ("trait", "other")
    if case == "no_traits":
        counts = [
            replace(profile, guild_counts={"trait": 1, "other": 3}, total=4) for profile in counts
        ]
    elif case == "no_traits_zero_counts":
        counts = [
            replace(profile, guild_counts={"trait": 0, "other": 0}, total=4) for profile in counts
        ]

    min_studies = 100 if case == "fewer_than_three_eligible_conditions" else 1
    protocol = Protocol(splits=2, min_studies=min_studies, permutations=19)
    split = split_module.make_splits(sorted(set(profiles.study_id)), splits=1, seed=protocol.seed)[
        0
    ]

    candidates, status = _training_selection(
        profiles,
        taxa,
        matrix,
        counts,
        columns,
        split,
        protocol,
        strain_rollup="off",
    )

    expected_status = "no_traits" if case == "no_traits_zero_counts" else case
    assert status == expected_status
    if case == "selected":
        assert candidates == [Candidate("trait", ("A", "B", "C"), 0.0)]
    else:
        assert candidates == []

    basis = build_analysis_basis(
        basis_id="fixture",
        feature_catalog_sha256="a" * 64,
        feature_filter_id=feature_filter_id(columns),
        retained_feature_count=len(columns),
        profile_unit="count_profiles",
        aggregation_policy="sum",
        rollup_policy="ncbi_taxonomy_closure_rollup_off",
        denominator_policy="retained",
        background_inventory_id="none",
        weighting_policy="equal_study",
    )
    records = run_splits(
        profiles,
        taxa,
        matrix,
        counts,
        columns,
        basis,
        [SourceArtifact(identity="fixture", sha256="c" * 64)],
        protocol,
        analysis_id="fixture",
        strain_rollup="off",
        recorded_at=datetime(2026, 10, 22, tzinfo=UTC),
    )
    assert len(records) == 2
    summary = summarize_stability(records)
    if case == "selected":
        assert all(record["selection_status"] == "selected" for record in records)
        assert len(summary) == 2
    else:
        assert summary == []
        for record in records:
            assert record["selection_status"] == expected_status
            assert record["selected_candidates"] == []
            assert record["test_details"] == []
            assert record["correction"]["family"]["expected_family_size"] == 0
            assert record["correction"]["results"] == []


def test_training_selection_converts_selected_without_disjoint_triples() -> None:
    """Direct unit check on study_split.py:182-185 when scanner.select returns SELECTED."""
    profiles, taxa, matrix, counts = _synthetic_corpus(shared_conditions=True)
    protocol = Protocol(splits=1, min_studies=1, permutations=19)
    split = StudySplit(
        0,
        tuple(f"s{index:02}" for index in range(9)),
        tuple(f"s{index:02}" for index in range(9, 18)),
    )
    candidates, status = _training_selection(
        profiles,
        taxa,
        matrix,
        counts,
        ("trait", "other"),
        split,
        protocol,
        strain_rollup="off",
    )
    assert candidates == []
    assert status == "no_testable_study_disjoint_triple"


@pytest.fixture
def verifier_997() -> Iterator[ModuleType]:
    """Load the helper by file path in both the live and copied test suites."""
    import importlib.util
    import sys

    path = REPO_ROOT / "scripts/adhoc/verify_997_mutants_agy_20261111.py"
    spec = importlib.util.spec_from_file_location("verifier_997_regression", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Nested verifier pytest runs evict these modules; preserve the outer suite's
    # identities so later tests patch the same objects their consumers import.
    outer_modules = {
        name: value
        for name, value in sys.modules.items()
        if name in {"mtp", "tests"} or name.startswith(("mtp.", "tests.", "test_"))
    }
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
        yield module
    finally:
        module._evict_project_modules()
        sys.modules.update(outer_modules)
        sys.modules.pop(spec.name, None)


@pytest.mark.parametrize("case", ["root", "nested", "symlink"])
def test_verifier_rejects_repository_scratch_before_mkdir(
    case: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, verifier_997: ModuleType
) -> None:
    scratch = REPO_ROOT
    if case == "nested":
        scratch = REPO_ROOT / "methods" / "nonexistent_scratch"
    elif case == "symlink":
        alias = tmp_path / "methods_alias"
        alias.symlink_to(REPO_ROOT / "methods", target_is_directory=True)
        scratch = alias / "sub"
    before = set((REPO_ROOT / "methods").iterdir())

    def forbidden_mkdir(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("containment must be checked before any mkdir")

    monkeypatch.setattr(Path, "mkdir", forbidden_mkdir)
    with pytest.raises(ValueError, match=r"^scratch_dir must be outside the repository root$"):
        verifier_997.verify_mutants(REPO_ROOT, scratch, canary_only=True)
    assert set((REPO_ROOT / "methods").iterdir()) == before
    if case != "root":
        assert not scratch.exists()


def test_verifier_rejects_lexically_internal_outward_symlink(
    tmp_path: Path, verifier_997: ModuleType
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    external = tmp_path / "external"
    external.mkdir()
    (repo / "alias").symlink_to(external, target_is_directory=True)
    with pytest.raises(ValueError, match="scratch_dir must be outside the repository root"):
        verifier_997.verify_mutants(repo, repo / "alias" / "new", canary_only=True)
    assert list(external.iterdir()) == []


def test_verifier_rejects_post_call_internal_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, verifier_997: ModuleType
) -> None:
    """Real pytest returns 3 after an assertion when logfinish raises INTERNALERROR."""
    import json

    source = tmp_path / "source"
    source.mkdir()

    def populate(_source: Path, destination: Path) -> None:
        destination.mkdir()
        (destination / "probe.py").write_text("VALUE = 1\n")
        (destination / "test_probe.py").write_text(
            "from pathlib import Path\n"
            "def test_probe():\n"
            "    assert Path('probe.py').read_text() == 'VALUE = 1\\n'\n"
        )
        (destination / "conftest.py").write_text(
            "failed_call = False\n"
            "def pytest_runtest_logreport(report):\n"
            "    global failed_call\n"
            "    if report.when == 'call':\n"
            "        failed_call = report.failed\n"
            "def pytest_runtest_logfinish(nodeid, location):\n"
            "    if failed_call:\n"
            "        raise RuntimeError('post-call internal error')\n"
        )

    original_run = verifier_997._run_pytest_node

    def run(source_root: Path, scratch: Path, node: str, junit: Path) -> tuple[int, Any, str]:
        target = "test_probe.py" if "::" not in node else node
        code, inspector, output = original_run(source_root, scratch, target, junit)
        return code, inspector, output

    monkeypatch.setattr(verifier_997, "_populate_scratch_repo", populate)
    monkeypatch.setattr(verifier_997, "_run_pytest_node", run)
    monkeypatch.setattr(
        verifier_997,
        "MUTANTS",
        (
            verifier_997.MutantSpec(
                "probe", "internal error", "probe.py", "1", "2", "test_probe.py::test_probe"
            ),
        ),
    )
    parent = tmp_path / "scratch"
    with pytest.raises(RuntimeError, match="not killed by call-phase AssertionError"):
        verifier_997.verify_mutants(source, parent, canary_only=True)
    receipts = list(parent.glob("*/verdicts.json"))
    assert len(receipts) == 1
    verdicts = json.loads(receipts[0].read_text())
    assert len(verdicts) == 1
    verdict = verdicts[0]
    assert verdict["status"] == "SURVIVED_OR_INVALID"
    assert verdict["exit_code"] == 3
    assert verdict["call_exc_types"] == ["AssertionError"]
    assert verdict["non_call_failures"] == []
    assert verdict["killed"] is False
    assert "INTERNALERROR" in verdict["detail"]
