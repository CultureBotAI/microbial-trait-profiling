"""Unit tests for Issue #1992: refuse missing trait keys in study_direction_means."""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path
from typing import Protocol, cast

from mtp.contrast_orientation import comparison_rows
from mtp.guild_maps import study_direction_means
from mtp.table1_trait_maps import compute_direction_shifts
from mtp.trait_guilds import TABLE1_TRAITS, Profile


class _VerifierModule(Protocol):
    def _run_targeted_pytest(
        self,
        scratch_repo: Path,
        nodeid: str,
        junit_xml: Path,
        *,
        repo_root: Path | None = None,
    ) -> tuple[bool, str]: ...


def _load_verifier_module() -> _VerifierModule:
    repo_root = Path(__file__).resolve().parents[1]
    script_path = repo_root / "scripts/adhoc/verify_1992_mutants_agy_20261111.py"
    spec = importlib.util.spec_from_file_location("verify_1992_mutants_agy_20261111", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return cast(_VerifierModule, module)


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


def test_mutant_verifier_rejects_setup_error_mentioning_assertion_error(tmp_path: Path) -> None:
    """Verifier rejects a setup-phase ValueError even when its text says AssertionError."""
    verifier = _load_verifier_module()
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
    killed, detail = verifier._run_targeted_pytest(
        scratch_repo, "tests/test_setup_boom.py::test_boom", junit_xml
    )
    assert killed is False
    assert "setup" in detail
    assert "ValueError" in detail
