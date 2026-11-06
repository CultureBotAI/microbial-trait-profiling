"""The shared comparison ledger must implement its spec's eligibility rule (#1611)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from mtp.datamodel.workflow_pydantic import AnalysisSpec, ArtifactInstance, WorkflowRun
from mtp.trait_guilds import HUMAN_HOST_SPECIES, OXYGEN_GUILDS, TABLE1_TRAITS
from mtp.workflow import comparison
from mtp.workflow.downstream_contracts import STANDARD_ANALYSIS_SPECS
from mtp.workflow.models import ContrastPopulation

TABLE1 = STANDARD_ANALYSIS_SPECS["noel_table1_six_trait_sign_test_v1"]
POOLED = STANDARD_ANALYSIS_SPECS["jed_pooled_contrast_fraction_diff_v1"]


@pytest.fixture(autouse=True)
def _local_import() -> None:
    assert Path(comparison.__file__).resolve() == (
        Path(__file__).resolve().parents[1] / "src/mtp/workflow/comparison.py"
    )


def _core(root: Path) -> WorkflowRun:
    """Minimal path bindings for this unit test, not a verified execution receipt."""
    return WorkflowRun.model_construct(
        manifest_path=str(root / "run.json"),
        artifacts=[
            ArtifactInstance.model_construct(product_id=name, path=f"{name}.tsv")
            for name in ("count_profiles", "experiments")
        ],
    )


@pytest.mark.parametrize(
    "spec",
    [
        pytest.param(POOLED, id="pooled-study-retains-zero"),
        pytest.param(
            TABLE1.model_copy(update={"eligibility_rule": POOLED.eligibility_rule}),
            id="noel-id-cannot-bypass-eligibility",
        ),
        pytest.param(
            STANDARD_ANALYSIS_SPECS["jed_table1_paired_log_ratio_sign_test_v1"],
            id="positive-ratio-is-a-different-rule",
        ),
        pytest.param(
            TABLE1.model_copy(update={"eligibility_rule": "unknown_rule"}),
            id="unknown-rule-refused",
        ),
    ],
)
def test_incompatible_eligibility_is_refused_before_reading_inputs(
    tmp_path: Path, spec: AnalysisSpec
) -> None:
    assert spec.sampling_unit == "study"
    # Neither input exists. Refusal must be for eligibility, before any file read.
    with pytest.raises(ValueError, match=r"shared-basis ledger.*eligibility") as error:
        comparison.shared_basis_ledger(spec, _core(tmp_path), ContrastPopulation.all_paired)
    assert spec.id in str(error.value)
    assert spec.eligibility_rule in str(error.value)


def test_supported_eligibility_does_not_bypass_the_sampling_unit(tmp_path: Path) -> None:
    spec = TABLE1.model_copy(update={"sampling_unit": "directional_signature"})
    with pytest.raises(ValueError, match="directional_signature"):
        comparison.shared_basis_ledger(spec, _core(tmp_path), ContrastPopulation.all_paired)


@pytest.mark.parametrize(
    "spec",
    [
        TABLE1,
        STANDARD_ANALYSIS_SPECS["noel_trait_direction_screen_sign_test_v1"],
        TABLE1.model_copy(update={"id": "same_eligibility_new_id"}),
    ],
    ids=["production-six-trait", "full-screen", "rule-not-id"],
)
def test_supported_rule_keeps_the_shared_basis_exclusion(
    tmp_path: Path, spec: AnalysisSpec
) -> None:
    terms = sorted({term for ids in {**TABLE1_TRAITS, **OXYGEN_GUILDS}.values() for term in ids})
    pairs = [
        ("bsdb:1", "1", 1, 1),  # A nonzero-total tie stays on the shared basis.
        ("bsdb:2", "1", 1, 0),  # The only experiment: this study is lost.
        ("bsdb:3", "1", 0, 1),  # Whole experiment excluded, but study survives.
        ("bsdb:3", "2", 2, 1),
    ]
    with (tmp_path / "count_profiles.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["study_id", "experiment_id", "inc_dec", *terms])
        for study, experiment, inc, dec in pairs:
            writer.writerow([study, experiment, "inc", *([inc] * len(terms))])
            writer.writerow([study, experiment, "dec", *([dec] * len(terms))])
    with (tmp_path / "experiments.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["study_id", "experiment_id", "condition", "body_site", "host_species"])
        for study, experiment, _, _ in pairs:
            writer.writerow([study, experiment, "example", "Feces", HUMAN_HOST_SPECIES])
    ledger = comparison.shared_basis_ledger(spec, _core(tmp_path), ContrastPopulation.all_paired)
    assert ledger.analysis_spec_id == spec.id
    assert ledger.contrast_population == ContrastPopulation.all_paired
    assert ledger.included_ids == ["bsdb:1", "bsdb:3"]
    assert ledger.included_count == 2
    assert ledger.excluded_count == 1
    assert json.loads(ledger.exclusion_counts_json) == {"zero_retained_total": 1}
