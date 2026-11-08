"""Cross-check real recorded Noel outputs against a synthetic verified core (#1029)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from mtp import influence_extension as extension
from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow.arms.runtime import execute_arm, verify_arm
from mtp.workflow.arms.types import ArmContext, Command
from mtp.workflow.consumers import verify_consumer_record
from test_influence_extension_core import extension_core as extension_core

PREPARED_NAMES = (
    "study_trait_profiles_wide.tsv",
    "experiments.tsv",
    "bugsigdb_taxa.tsv",
    "taxa_trait_matrix.tsv",
)


def noel_fixture_plan(
    context: ArmContext,
    *,
    alter_p: bool = False,
    alter_q: bool = False,
    omit_effects: bool = False,
) -> list[Command]:
    """Produce normal or deliberately mismatched artifacts before arm hashes are sealed."""
    producer = REPO_ROOT / "scripts/adhoc/issue1030_noel_fixture.py"
    protocol = REPO_ROOT / "methods/table1_taxon_influence_protocol.json"
    output = context.output / "table1_taxon_influence"
    arguments = (
        sys.executable,
        str(producer),
        "--inputs",
        str(context.inputs),
        "--output",
        str(output),
    )
    changes = tuple(
        flag
        for enabled, flag in (
            (alter_p, "--alter-p"),
            (alter_q, "--alter-q"),
            (omit_effects, "--omit-effects"),
        )
        if enabled
    )
    return [
        Command(
            "taxon_contributions",
            (*arguments, *changes),
            (*(context.inputs / name for name in PREPARED_NAMES), protocol),
            (output,),
        ),
        Command(
            "native_taxon_influence",
            (*arguments, "--summary-only"),
            (output,),
            (output / "native_summary.md",),
        ),
    ]


def _produce_noel(
    core: Path,
    output: Path,
    *,
    alter_p: bool = False,
    alter_q: bool = False,
    omit_effects: bool = False,
) -> Path:
    def plan(context: ArmContext) -> list[Command]:
        return noel_fixture_plan(
            context, alter_p=alter_p, alter_q=alter_q, omit_effects=omit_effects
        )

    run = execute_arm("noel", core, output, plan, figures=False)
    manifest = Path(run.manifest_path)
    verified = verify_arm(manifest)
    assert verified.arm == "noel"
    assert verified.core_manifest_sha256 == file_sha256(core)
    assert [step.id for step in verified.steps or []] == [
        "taxon_contributions",
        "native_taxon_influence",
    ]
    return manifest


def test_verified_noel_baseline_crosscheck_is_recorded_without_changing_arm(tmp_path: Path) -> None:
    core = extension_core(tmp_path)
    noel = _produce_noel(core, tmp_path / "noel")
    effects = noel.parent / "table1_taxon_influence/effects.tsv.gz"
    summary = noel.parent / "table1_taxon_influence/native_summary.md"
    before = {path: file_sha256(path) for path in (noel, effects, summary)}
    output = tmp_path / "extension"

    result = extension.run_extension(core, output, noel)

    provenance = json.loads(result.read_text(encoding="utf-8"))
    assert "noel_baseline_check" in provenance
    check = provenance["noel_baseline_check"]
    assert check["status"] == "matched"
    assert check["rows"] == 24
    assert check["families"] == 4
    assert check["by_source"] == "derived_from_noel_p_values"
    assert provenance["input_sha256"]["noel_manifest"] == before[noel]
    assert provenance["input_sha256"]["native_summary"] == before[summary]
    assert summary.read_text(encoding="utf-8") in (
        output / "influence_extension_summary.md"
    ).read_text(encoding="utf-8")
    record = verify_consumer_record(output / "core_consumption.json")
    assert record.run_id == provenance["core_run_id"]
    assert {path: file_sha256(path) for path in before} == before
    verify_arm(noel)


@pytest.mark.parametrize("field", ["p", "q"])
def test_verified_but_mismatched_noel_probability_refuses_before_publication(
    tmp_path: Path, field: str
) -> None:
    core = extension_core(tmp_path)
    noel = _produce_noel(core, tmp_path / "noel", alter_p=field == "p", alter_q=field == "q")
    before = file_sha256(noel)
    # Verification has already succeeded: the altered probability is a produced,
    # hash-bound artifact, so only the statistical comparison can reject it.
    output = tmp_path / "must-not-publish"
    with pytest.raises(ValueError, match=r"(?i)Noel.*(baseline|mismatch|differ)"):
        extension.run_extension(core, output, noel)
    assert not output.exists()
    assert file_sha256(noel) == before
    verify_arm(noel)


def test_verified_noel_without_declared_effects_refuses_before_publication(tmp_path: Path) -> None:
    core = extension_core(tmp_path)
    noel = _produce_noel(core, tmp_path / "noel", omit_effects=True)
    verified = verify_arm(noel)
    assert verified.steps
    assert verified.steps[0].outputs
    assert all(
        member.path != "effects.tsv.gz"
        for artifact in verified.steps[0].outputs
        for member in artifact.inventory
    )
    output = tmp_path / "must-not-publish"
    with pytest.raises(ValueError, match=r"effects\.tsv\.gz"):
        extension.run_extension(core, output, noel)
    assert not output.exists()
    verify_arm(noel)
