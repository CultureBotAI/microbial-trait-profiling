"""Unit tests for Issue #1992: refuse missing trait keys in study_direction_means."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from mtp.contrast_orientation import comparison_rows
from mtp.guild_maps import study_direction_means
from mtp.table1_trait_maps import compute_direction_shifts
from mtp.trait_guilds import TABLE1_TRAITS, Profile


@pytest.fixture
def verifier_1992() -> Iterator[ModuleType]:
    """Load the verifier module by file path and restore outer project module identities."""
    repo_root = Path(__file__).resolve().parents[1]
    script_path = repo_root / "scripts/adhoc/verify_1992_mutants_agy_20261111.py"
    spec = importlib.util.spec_from_file_location("verify_1992_mutants_agy_20261111", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
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


def test_study_direction_means_rejects_missing_trait_in_later_profile_and_later_guild() -> None:
    """Direct study_direction_means call refuses a Profile missing a guild in guild_order."""
    profiles = [
        Profile(
            study_id="study_1",
            experiment_id="exp_1",
            inc_dec="inc",
            condition="Condition A",
            body_site="Feces",
            guild_counts={"aerobe": 4, "anaerobe": 6},
            total=10,
        ),
        Profile(
            study_id="study_1",
            experiment_id="exp_1",
            inc_dec="dec",
            condition="Condition A",
            body_site="Feces",
            guild_counts={"aerobe": 2, "anaerobe": 8},
            total=10,
        ),
        Profile(
            study_id="study_2",
            experiment_id="exp_2",
            inc_dec="inc",
            condition="Condition A",
            body_site="Feces",
            guild_counts={"aerobe": 5},
            total=10,
        ),
        Profile(
            study_id="study_2",
            experiment_id="exp_2",
            inc_dec="dec",
            condition="Condition A",
            body_site="Feces",
            guild_counts={"aerobe": 3, "anaerobe": 7},
            total=10,
        ),
    ]
    caught: ValueError | None = None
    try:
        study_direction_means(profiles, guild_order=("aerobe", "anaerobe"))
    except ValueError as exc:
        caught = exc

    assert caught is not None
    assert re.search(r"Missing trait .* for \('.*', '.*'\); missing is not zero", str(caught))
    assert str(caught) == ("Missing trait anaerobe for ('study_2', 'exp_2'); missing is not zero")


def test_compute_direction_shifts_rejects_missing_trait() -> None:
    """compute_direction_shifts refuses a Profile missing a trait in trait_order."""
    profiles = [
        Profile(
            study_id="study_1",
            experiment_id="exp_1",
            inc_dec="inc",
            condition="Condition A",
            body_site="Feces",
            guild_counts={"facultatively_aerobic": 6},
            total=10,
        ),
        Profile(
            study_id="study_1",
            experiment_id="exp_1",
            inc_dec="dec",
            condition="Condition A",
            body_site="Feces",
            guild_counts={
                "facultatively_aerobic": 2,
                "cellulose_catabolic_process": 0,
            },
            total=10,
        ),
    ]
    caught: ValueError | None = None
    try:
        compute_direction_shifts(
            profiles,
            trait_order=("facultatively_aerobic", "cellulose_catabolic_process"),
        )
    except ValueError as exc:
        caught = exc

    assert caught is not None
    assert str(caught) == (
        "Missing trait cellulose_catabolic_process for ('study_1', 'exp_1'); missing is not zero"
    )


def test_comparison_rows_rejects_absent_table1_trait_and_accepts_explicit_zero() -> None:
    """Exact #1992 scenario: comparison_rows rejects missing Table 1 trait, accepts explicit 0."""
    counts_without_cellulose = {
        trait: (4 if trait == "facultatively_aerobic" else 2)
        for trait in TABLE1_TRAITS
        if trait != "cellulose_catabolic_process"
    }
    profiles_missing = [
        Profile(
            study_id="study_1",
            experiment_id="exp_1",
            inc_dec="inc",
            condition="Condition A",
            body_site="Feces",
            guild_counts=dict(counts_without_cellulose),
            total=10,
        ),
        Profile(
            study_id="study_1",
            experiment_id="exp_1",
            inc_dec="dec",
            condition="Condition A",
            body_site="Feces",
            guild_counts=dict(counts_without_cellulose),
            total=10,
        ),
    ]
    caught: ValueError | None = None
    try:
        comparison_rows(
            [("pop", profiles_missing)],
            ratio_traits=("facultatively_aerobic", "coccus_shaped"),
        )
    except ValueError as exc:
        caught = exc

    assert caught is not None
    assert str(caught) == (
        "Missing trait cellulose_catabolic_process for ('study_1', 'exp_1'); missing is not zero"
    )

    counts_with_explicit_zero = {
        **counts_without_cellulose,
        "cellulose_catabolic_process": 0,
    }
    profiles_explicit_zero = [
        Profile(
            study_id="study_1",
            experiment_id="exp_1",
            inc_dec="inc",
            condition="Condition A",
            body_site="Feces",
            guild_counts=dict(counts_with_explicit_zero),
            total=10,
        ),
        Profile(
            study_id="study_1",
            experiment_id="exp_1",
            inc_dec="dec",
            condition="Condition A",
            body_site="Feces",
            guild_counts=dict(counts_with_explicit_zero),
            total=10,
        ),
    ]
    unexpected: ValueError | None = None
    rows: list[dict[str, object]] = []
    try:
        rows = comparison_rows(
            [("pop", profiles_explicit_zero)],
            ratio_traits=("facultatively_aerobic", "coccus_shaped"),
        )
    except ValueError as exc:
        unexpected = exc

    assert unexpected is None
    cellulose_rows = [
        row
        for row in rows
        if row["comparison"] == "fraction" and row["numerator"] == "cellulose_catabolic_process"
    ]
    assert len(cellulose_rows) == 1
    cellulose_row = cellulose_rows[0]
    assert cellulose_row["studies_scored"] == 0
    assert cellulose_row["studies_tied"] == 1


def test_study_direction_means_and_comparison_rows_reject_missing_trait_on_decreased_profile() -> (
    None
):
    """F3 regression: complete increased counts + decreased profile missing a trait must raise."""
    counts_complete = {
        trait: (4 if trait == "facultatively_aerobic" else 2) for trait in TABLE1_TRAITS
    }
    counts_dec_missing_cellulose = {
        trait: val
        for trait, val in counts_complete.items()
        if trait != "cellulose_catabolic_process"
    }
    profiles_dec_missing = [
        Profile(
            study_id="study_1",
            experiment_id="exp_1",
            inc_dec="inc",
            condition="Condition A",
            body_site="Feces",
            guild_counts=dict(counts_complete),
            total=10,
        ),
        Profile(
            study_id="study_1",
            experiment_id="exp_1",
            inc_dec="dec",
            condition="Condition A",
            body_site="Feces",
            guild_counts=dict(counts_dec_missing_cellulose),
            total=10,
        ),
    ]

    direct_caught: ValueError | None = None
    try:
        study_direction_means(profiles_dec_missing, guild_order=tuple(TABLE1_TRAITS))
    except ValueError as exc:
        direct_caught = exc

    assert direct_caught is not None
    assert str(direct_caught) == (
        "Missing trait cellulose_catabolic_process for ('study_1', 'exp_1'); missing is not zero"
    )

    comparison_caught: ValueError | None = None
    try:
        comparison_rows(
            [("pop", profiles_dec_missing)],
            ratio_traits=("facultatively_aerobic", "coccus_shaped"),
        )
    except ValueError as exc:
        comparison_caught = exc

    assert comparison_caught is not None
    assert str(comparison_caught) == (
        "Missing trait cellulose_catabolic_process for ('study_1', 'exp_1'); missing is not zero"
    )


def test_mutant_verifier_rejects_setup_error_mentioning_assertion_error(
    tmp_path: Path, verifier_1992: ModuleType
) -> None:
    """Verifier rejects a setup-phase ValueError even when its text says AssertionError."""
    scratch_repo = tmp_path / "repo"
    (scratch_repo / "src").mkdir(parents=True)
    tests_dir = scratch_repo / "tests"
    tests_dir.mkdir(parents=True)
    test_file = tests_dir / "test_setup_boom.py"
    test_file.write_text(
        "import pytest\n\n"
        "@pytest.fixture\n"
        "def bad_fixture() -> None:\n"
        '    raise ValueError("AssertionError in setup")\n\n'
        "def test_boom(bad_fixture: None) -> None:\n"
        "    assert True\n",
        encoding="utf-8",
    )
    junit_xml = tmp_path / "junit_setup_boom.xml"
    killed, detail = verifier_1992._run_targeted_pytest(
        scratch_repo, "tests/test_setup_boom.py::test_boom", junit_xml
    )
    assert killed is False
    assert "setup" in detail
    assert "ValueError" in detail


@pytest.mark.parametrize("case", ["root", "nested", "symlink"])
def test_verifier_rejects_repository_scratch_before_mkdir(
    case: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, verifier_1992: ModuleType
) -> None:
    """F1 regression: reject repository root, nested path, or inbound symlink before any mkdir."""
    repo_root = Path(__file__).resolve().parents[1]
    scratch = repo_root
    if case == "nested":
        scratch = repo_root / "src" / "nonexistent_scratch"
    elif case == "symlink":
        alias = tmp_path / "src_alias"
        alias.symlink_to(repo_root / "src", target_is_directory=True)
        scratch = alias / "sub"
    before = set((repo_root / "src").iterdir())

    def forbidden_mkdir(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("containment must be checked before any mkdir")

    monkeypatch.setattr(Path, "mkdir", forbidden_mkdir)
    with pytest.raises(ValueError, match=r"^scratch_dir must be outside the repository root$"):
        verifier_1992.verify_mutants(repo_root, scratch, canary_only=True)
    assert set((repo_root / "src").iterdir()) == before
    if case != "root":
        assert not scratch.exists()


def test_verifier_rejects_lexically_internal_outward_symlink_and_accepts_external(
    tmp_path: Path, verifier_1992: ModuleType
) -> None:
    """F1 regression: reject lexically internal outward symlink and accept external scratch dir."""
    repo = tmp_path / "repo"
    repo.mkdir()
    external = tmp_path / "external"
    external.mkdir()
    (repo / "alias").symlink_to(external, target_is_directory=True)
    with pytest.raises(ValueError, match=r"^scratch_dir must be outside the repository root$"):
        verifier_1992.verify_mutants(repo, repo / "alias" / "new", canary_only=True)
    assert list(external.iterdir()) == []

    verifier_1992._validate_scratch_parent(repo, external / "allowed")
    parsed_all = verifier_1992._parse_args(["--scratch-dir", str(external / "allowed"), "--all"])
    assert parsed_all.all is True
    assert parsed_all.canary is False


def test_verifier_rejects_post_call_internal_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, verifier_1992: ModuleType
) -> None:
    """F2 regression: real pytest exit 3 after call-phase AssertionError is rejected."""
    source = tmp_path / "source"
    source.mkdir()

    def populate(_source: Path, destination: Path) -> None:
        destination.mkdir()
        (destination / "probe.py").write_text("VALUE = 1\n", encoding="utf-8")
        (destination / "test_probe.py").write_text(
            "from pathlib import Path\n"
            "def test_probe() -> None:\n"
            "    assert Path('probe.py').read_text(encoding='utf-8') == 'VALUE = 1\\n'\n",
            encoding="utf-8",
        )
        (destination / "conftest.py").write_text(
            "failed_call = False\n"
            "def pytest_runtest_logreport(report):\n"
            "    global failed_call\n"
            "    if report.when == 'call':\n"
            "        failed_call = report.failed\n"
            "def pytest_runtest_logfinish(nodeid, location):\n"
            "    if failed_call:\n"
            "        raise RuntimeError('post-call internal error')\n",
            encoding="utf-8",
        )

    original_run = verifier_1992._run_pytest_node

    def run(source_root: Path, scratch: Path, node: str, junit: Path) -> tuple[int, Any, str]:
        target = "test_probe.py" if "::" not in node else node
        code, inspector, output = original_run(source_root, scratch, target, junit)
        return code, inspector, output

    monkeypatch.setattr(verifier_1992, "_populate_scratch_repo", populate)
    monkeypatch.setattr(verifier_1992, "_run_pytest_node", run)
    monkeypatch.setattr(
        verifier_1992,
        "MUTANTS",
        (
            verifier_1992.MutantSpec(
                "probe", "internal error", "probe.py", "1", "2", "test_probe.py::test_probe"
            ),
        ),
    )
    parent = tmp_path / "scratch"
    with pytest.raises(RuntimeError, match="not killed by call-phase AssertionError"):
        verifier_1992.verify_mutants(source, parent, canary_only=True)
    receipts = list(parent.glob("*/verdicts.json"))
    assert len(receipts) == 1
    verdicts = json.loads(receipts[0].read_text(encoding="utf-8"))
    assert len(verdicts) == 1
    verdict = verdicts[0]
    assert verdict["status"] == "SURVIVED_OR_INVALID"
    assert verdict["exit_code"] == 3
    assert verdict["call_exc_types"] == ["AssertionError"]
    assert verdict["non_call_failures"] == []
    assert verdict["killed"] is False
    assert "INTERNALERROR" in verdict["detail"]
