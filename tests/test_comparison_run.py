"""The compare stage reads its inputs and verdicts from verified arm runs (#814, #815, #844).

The end-to-end fixture runs the real arm steps whose outputs the stage consumes on a
small synthetic corpus that carries every Table 1 and oxygen-guild term. Study bsdb:10
has no annotated taxon on its decreased side, so the shared-basis exclusion drops it.
"""

from __future__ import annotations

import csv
import json
import random
import shutil
import sys
from collections.abc import Callable
from pathlib import Path

import pytest

from mtp.datamodel.workflow_pydantic import (
    ArmRun,
    CalibrationOutcome,
    CalibrationVerdict,
    ComparisonRun,
    EligibilityLedger,
)
from mtp.integrated_evidence import main as integrated_evidence_main
from mtp.paths import REPO_ROOT
from mtp.workflow import comparison
from mtp.workflow.arms import ilyas, jed, noel
from mtp.workflow.arms.__main__ import main as arms_main
from mtp.workflow.arms.runtime import INPUT_NAMES, Planner, execute_arm, verify_arm
from mtp.workflow.arms.types import ArmContext, Command
from mtp.workflow.comparison import (
    ARM_EVIDENCE_INPUTS,
    CALIBRATION_SOURCES,
    execute_comparison_run,
    jed_calibration_verdict,
    noel_calibration_verdict,
    reconcile_ledger,
    recorded_output,
    shared_basis_ledger,
    verify_comparison_run,
)
from mtp.workflow.downstream_contracts import STANDARD_ANALYSIS_SPECS, build_eligibility_ledger
from mtp.workflow.models import ContrastPopulation
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow
from mtp.workflow.verification import verify_core_run

INFLUENCE_FIXTURE = REPO_ROOT / "scripts" / "adhoc" / "issue844_influence_assessment_fixture.py"
TABLE1_SPEC = STANDARD_ANALYSIS_SPECS["noel_table1_six_trait_sign_test_v1"]

#: Every Table 1 and oxygen-guild term, with the predicate that attaches it to a taxon.
TRAITS = {
    "METPO:1000602": "biolink:has_phenotype",
    "METPO:1000603": "biolink:has_phenotype",
    "METPO:1000604": "biolink:has_phenotype",
    "METPO:1000607": "biolink:has_phenotype",
    "METPO:1000608": "biolink:has_phenotype",
    "METPO:1000623": "biolink:has_phenotype",
    "METPO:1000668": "biolink:has_phenotype",
    "METPO:1000683": "biolink:has_phenotype",
    "GO:0030245": "METPO:2000103",
    "pathways:arsenate_reduction": "METPO:2000102",
}
DUMP_HEADER = [
    "BSDB ID",
    "Study",
    "Title",
    "Year",
    "Host species",
    "Body site",
    "Condition",
    "Group 0 name",
    "Group 1 name",
    "Abundance in Group 1",
    "NCBI Taxonomy IDs",
    "MetaPhlAn taxon names",
    "State",
    "Reviewer",
    "Sequencing type",
    "PMID",
    "DOI",
    "Study design",
    "EFO ID",
    "UBERON ID",
    "Group 1 definition",
    "Group 0 sample size",
    "Group 1 sample size",
    "Description",
]

#: The real arm steps whose outputs the compare stage reads, and their prerequisites.
CONSUMED_STEPS = {
    "noel": {"directional_screen", "empirical_calibration", "contrast_orientation"},
    "jed": {"human_counts", "composition_sensitivity", "composition_calibration"},
    "ilyas": {"enrichment", "depletion"},
}
PLANNERS: dict[str, Planner] = {
    "noel": noel.build_commands,
    "jed": jed.build_commands,
    "ilyas": ilyas.build_commands,
}


def _write_inputs(root: Path) -> dict[str, Path]:
    """A raw export and KG whose taxa carry every guild term (deterministic)."""
    rng = random.Random(7)
    root.mkdir(parents=True)
    species = {
        f"NCBITaxon:{1000 + index}": (
            f"NCBITaxon:{2000 + index % 4}",
            f"Synthgenus{index % 4} species{index}",
            sorted(rng.sample(sorted(TRAITS), rng.randint(2, 5))),
        )
        for index in range(1, 25)
    }
    nodes = [["id", "category", "name"], ["NCBITaxon:3000", "biolink:OrganismTaxon", "Family"]]
    nodes += [["NCBITaxon:1999", "biolink:OrganismTaxon", "Unannotated bacterium"]]
    edges = [["id", "subject", "predicate", "object", "primary_knowledge_source"]]
    edges += [["t1999", "NCBITaxon:1999", "biolink:subclass_of", "NCBITaxon:3000", "fixture"]]
    for genus in sorted({genus for genus, _, _ in species.values()}):
        nodes.append([genus, "biolink:OrganismTaxon", f"Synthgenus{genus[-1]}"])
        edges.append([f"t{genus}", genus, "biolink:subclass_of", "NCBITaxon:3000", "fixture"])
    for taxon, (genus, name, traits) in species.items():
        nodes.append([taxon, "biolink:OrganismTaxon", name])
        edges.append([f"t{taxon}", taxon, "biolink:subclass_of", genus, "fixture"])
        edges += [[f"a{taxon}{term}", taxon, TRAITS[term], term, "fixture"] for term in traits]
    nodes += [[term, "biolink:PhenotypicQuality", term] for term in TRAITS]
    nodes += [["METPO:2000102", "biolink:OntologyClass", "has pathway"]]
    nodes += [["METPO:2000103", "biolink:OntologyClass", "capable of"]]
    paths = {"kg_edges": root / "merged-kg_edges.tsv", "kg_nodes": root / "merged-kg_nodes.tsv"}
    for key, rows in (("kg_edges", edges), ("kg_nodes", nodes)):
        with paths[key].open("w", encoding="utf-8", newline="") as handle:
            csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)

    conditions = ("Crohn's disease", "colorectal cancer", "Parkinson's disease", "obesity")
    sites = ("feces", "feces", "feces", "saliva")
    taxa = sorted(species)
    paths["raw_dump"] = root / "full_dump.csv"
    with paths["raw_dump"].open("w", encoding="utf-8", newline="") as handle:
        handle.write("# BugSigDB fixture-814, License: CC0 synthetic test data\n")
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(DUMP_HEADER)
        for study, experiment, (signature, direction) in (
            (study, experiment, side)
            for study in range(1, 11)
            for experiment in (1, 2)
            for side in ((1, "increased"), (2, "decreased"))
        ):
            chosen = (
                ["NCBITaxon:1999"]
                if study == 10 and direction == "decreased"
                else sorted(rng.sample(taxa, rng.randint(3, 6)))
            )
            names = [
                "s__" + species[taxon][1].replace(" ", "_")
                if taxon in species
                else "s__Unannotated_bacterium"
                for taxon in chosen
            ]
            writer.writerow(
                [
                    f"bsdb:{study}/{experiment}/{signature}",
                    f"Study {study}",
                    f"Synthetic study {study}",
                    "2026",
                    "Homo sapiens",
                    sites[study % 4],
                    conditions[study % 4],
                    "healthy controls",
                    "patients",
                    direction,
                    ";".join(taxon.removeprefix("NCBITaxon:") for taxon in chosen),
                    ",".join(names),
                    "Complete",
                    "Fixture",
                    "WMS" if study % 3 == 0 else "16S",
                    str(1000 + study),
                    f"10.0/fixture.{study}",
                    "case-control",
                    "",
                    "",
                    "patients with the condition",
                    "20",
                    "20",
                    "",
                ]
            )
    return paths


def _influence_stand_in(context: ArmContext) -> Command:
    """The influence protocol only runs on the production KG bytes, so an explicit
    assessment bound to this core's inputs takes the native manifest's place."""
    profiles, experiments = context.products["count_profiles"], context.products["experiments"]
    target = context.output / "table1_taxon_influence" / "native_manifest.json"
    argv = (
        sys.executable,
        str(INFLUENCE_FIXTURE),
        "--profiles",
        str(profiles),
        "--experiments",
        str(experiments),
        "--output",
        str(target),
    )
    return Command("native_taxon_influence", argv, (profiles, experiments), (target,))


def _consumed_plan(arm: str) -> Planner:
    def plan(context: ArmContext) -> list[Command]:
        commands = [item for item in PLANNERS[arm](context) if item.id in CONSUMED_STEPS[arm]]
        return [*commands, _influence_stand_in(context)] if arm == "noel" else commands

    return plan


@pytest.fixture(scope="module")
def stage(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """One verified core, three verified arm runs and the comparison run over them."""
    root = tmp_path_factory.mktemp("comparison-stage")
    core = execute_workflow(
        REPO_ROOT / "workflows" / "bugsigdb_core.yaml",
        _write_inputs(root / "raw"),
        root / "core",
        OPERATIONS,
        parameter_overrides={"source_mode": "fixture"},
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    paths = {"core": Path(core.manifest_path), "comparison": root / "comparison"}
    for arm in ("noel", "jed", "ilyas"):
        record = execute_arm(
            arm, paths["core"], root / arm, _consumed_plan(arm), replicates=2, seed=0, figures=False
        )
        paths[arm] = Path(record.manifest_path)
    arguments = [str(paths[name]) for name in ("core", "noel", "jed", "ilyas", "comparison")]
    # The stage lives in temporary storage, which compare refuses without the opt-in
    assert arms_main(["compare", *arguments, "--allow-temporary-output"]) == 0
    return paths


def _comparison(stage: dict[str, Path]) -> ComparisonRun:
    path = stage["comparison"] / "comparison_run.json"
    return ComparisonRun.model_validate_json(path.read_text(encoding="utf-8"))


def _provenance(stage: dict[str, Path]) -> dict[str, object]:
    path = stage["comparison"] / "integrated_evidence" / "integrated_evidence_provenance.json"
    content = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(content, dict)
    return content


def test_every_file_the_stage_reads_is_declared_by_its_arm_plan(tmp_path: Path) -> None:
    # The fixture stands in for Noel's influence step, so the real plans are checked here.
    products = {product: tmp_path / "core" / name for product, name in INPUT_NAMES.items()}
    context = ArmContext(
        REPO_ROOT,
        tmp_path / "core" / "run.json",
        tmp_path / "arm",
        tmp_path / "arm" / "inputs",
        products,
        "a" * 40,
        None,
        None,
        True,
    )
    declared = {
        (arm, command.id, path.name)
        for arm, planner in PLANNERS.items()
        for command in planner(context)
        for path in command.outputs
    }
    read = {(arm, step, name) for _, arm, step, name in ARM_EVIDENCE_INPUTS}
    assert read | set(CALIBRATION_SOURCES) <= declared


def test_compare_stage_runs_on_arm_runs_and_verifies(stage: dict[str, Path]) -> None:
    assert arms_main(["verify-comparison", str(stage["comparison"])]) == 0


def test_every_evidence_input_is_the_file_its_arm_step_recorded(stage: dict[str, Path]) -> None:
    runs = {arm: verify_arm(stage[arm]) for arm in ("noel", "jed", "ilyas")}
    expected = {
        option: str(recorded_output(runs[arm], step, name))
        for option, arm, step, name in ARM_EVIDENCE_INPUTS
    }
    assert expected["--arm2-enrichment-results"] == str(
        stage["ilyas"].parent / "trait_enrichment_all_results.tsv"
    )
    assert expected["--full-screen"] == str(
        stage["noel"].parent / "trait_direction_screen" / "trait_direction_screen.tsv"
    )
    assert expected["--sensitivity-status"] == str(
        stage["jed"].parent / "composition_sensitivity" / "sensitivity_status.json"
    )
    provenance = _provenance(stage)
    native = provenance["native_sources"]
    assert isinstance(native, dict)
    for key, option in (
        ("arm2_enrichment", "--arm2-enrichment-results"),
        ("arm2_depletion", "--arm2-depletion-results"),
        ("full_screen", "--full-screen"),
        ("full_screen_provenance", "--full-screen-provenance"),
    ):
        assert native[key]["path"] == expected[option]
    statuses = provenance["assessment_status"]
    assert isinstance(statuses, dict)
    for key, option in (
        ("calibration_status", "--calibration-status"),
        ("orientation_status", "--orientation-provenance"),
        ("sensitivity_status", "--sensitivity-status"),
    ):
        assert statuses[key]["source"]["path"] == expected[option]
    # Nothing produces the #184 coverage status, so it is recorded as not supplied.
    assert statuses["coverage_status"]["source"] == {"path": None, "binding": "not_supplied"}


def test_the_ledger_lists_the_studies_integrated_evidence_kept(stage: dict[str, Path]) -> None:
    (ledger,) = _comparison(stage).eligibility_ledgers
    assert ledger.analysis_spec_id == TABLE1_SPEC.id
    assert TABLE1_SPEC.sampling_unit == "study"
    assert ledger.included_ids == [f"bsdb:{study}" for study in range(1, 10)]
    assert ledger.included_count == 9
    assert ledger.excluded_count == 1
    assert json.loads(ledger.exclusion_counts_json) == {"zero_retained_total": 1}


def test_the_verdicts_are_the_arms_own_calibration_results(stage: dict[str, Path]) -> None:
    verdicts = {verdict.id: verdict for verdict in _comparison(stage).calibration_verdicts}
    # Two replicates leave every criterion of Noel's harness below its declared budget.
    table1 = verdicts["verdict:noel_table1_six_trait_sign_test_v1"]
    assert table1.estimator_match is True
    assert table1.outcome == CalibrationOutcome.insufficient_budget
    assert table1.replicate_budget_sufficient is False
    assert table1.null_control_passed is None
    assert table1.power_passed is None
    assert "fraction_sign_test not_certified [family_error=insufficient_budget" in (
        table1.summary_message
    )
    assert "screened_fraction_sign_test not_certified" in table1.summary_message
    assert table1.negative_control_procedure == "screened_fraction_sign_test"
    assert table1.negative_control_outcome == "not_certified"
    # Two replicates cannot settle whether the control fails as declared.
    assert table1.negative_control_as_expected is None
    assert table1.selection_uncertainty is not None
    assert "p-values are exploratory" in table1.selection_uncertainty
    assert table1.selection_uncertainty in table1.summary_message
    composition = verdicts["verdict:jed_log_ratio_vs_transformed_t_calibration"]
    assert composition.selection_uncertainty is None
    assert composition.outcome == CalibrationOutcome.different_estimator_uncalibrated
    assert composition.replicate_budget_sufficient is False
    assert composition.null_control_passed is None
    assert composition.power_passed is None
    assert "2 of 2 replicates completed" in composition.summary_message


def test_the_influence_status_is_read_and_bound(stage: dict[str, Path]) -> None:
    statuses = _provenance(stage)["assessment_status"]
    assert isinstance(statuses, dict)
    influence = statuses["taxon_influence_status"]
    assert influence["verdicts_by_procedure"] == {"fraction_sign_test": "fixture_assessment #190"}
    assert not influence["status_for_a_procedure_it_is_silent_about"].startswith("unreadable")


def _evidence_argv(
    stage: dict[str, Path], output: Path, statuses: dict[str, str | None]
) -> list[str]:
    """The stage's own integrated_evidence arguments with some status inputs replaced."""
    runs = {arm: verify_arm(stage[arm]) for arm in ("noel", "jed", "ilyas")}
    inputs: dict[str, str | None] = {
        option: str(recorded_output(runs[arm], step, name))
        for option, arm, step, name in ARM_EVIDENCE_INPUTS
    }
    inputs.update(statuses)
    pooled = stage["comparison"] / "arm1_pooled_stats"
    return [
        "--core-run",
        str(stage["core"]),
        "--output-dir",
        str(output),
        "--arm1-pooled-stats",
        str(pooled / "trait_profile_paired_contrast_trait_stats_pooled.tsv"),
        *(part for option, path in inputs.items() if path is not None for part in (option, path)),
        "--require-complete-sources",
    ]


@pytest.mark.parametrize(
    ("option", "supplied"),
    [
        ("--calibration-status", False),
        ("--sensitivity-status", True),
        ("--orientation-provenance", True),
        ("--influence-status", False),
        ("--coverage-status", True),
    ],
)
def test_complete_sources_reject_a_missing_status_input(
    stage: dict[str, Path], tmp_path: Path, option: str, supplied: bool
) -> None:
    missing = str(tmp_path / "no_such_status.json") if supplied else None
    argv = _evidence_argv(stage, tmp_path / "integrated_evidence", {option: missing})
    status = option.removeprefix("--").replace("-", "_")
    with pytest.raises(ValueError, match=f"missing source files: {status}"):
        integrated_evidence_main(argv)


def test_complete_sources_reject_an_unreadable_status(
    stage: dict[str, Path], tmp_path: Path
) -> None:
    # The influence file the fixture used before #844: valid JSON, but not a status.
    wrong = tmp_path / "native_manifest.json"
    shutil.copyfile(REPO_ROOT / "methods" / "table1_native_influence_protocol.json", wrong)
    argv = _evidence_argv(stage, tmp_path / "out", {"--influence-status": str(wrong)})
    with pytest.raises(ValueError, match=r"unbound statuses: taxon_influence_status \(unreadable"):
        integrated_evidence_main(argv)


def test_complete_sources_reject_a_status_from_other_input_bytes(
    stage: dict[str, Path], tmp_path: Path
) -> None:
    source = recorded_output(
        verify_arm(stage["noel"]), "empirical_calibration", "calibration_status.json"
    )
    content = json.loads(source.read_text(encoding="utf-8"))
    content["inputs"]["profile_sha256"] = "0" * 64
    other = tmp_path / "calibration_status.json"
    other.write_text(json.dumps(content), encoding="utf-8")
    output = tmp_path / "out"
    argv = _evidence_argv(stage, output, {"--calibration-status": str(other)})
    with pytest.raises(ValueError, match="unbound statuses: calibration_status"):
        integrated_evidence_main(argv)
    # Without the flag the same status is recorded rather than rejected.
    integrated_evidence_main([item for item in argv if item != "--require-complete-sources"])
    recorded = json.loads((output / "integrated_evidence_provenance.json").read_text())
    binding = recorded["assessment_status"]["calibration_status"]["source"]["binding"]
    assert binding == "mismatch_or_missing_identity"


@pytest.mark.parametrize(
    ("fault", "message"),
    [
        ("calibration", "absent_control"),
        ("ledger", "unreadable count profiles"),
        ("bases", "unreadable count profiles"),
    ],
)
def test_a_bad_verdict_ledger_or_basis_input_writes_nothing(
    stage: dict[str, Path],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fault: str,
    message: str,
) -> None:
    def unreadable(*_args: object, **_kwargs: object) -> object:
        raise ValueError("unreadable count profiles")

    if fault == "calibration":
        monkeypatch.setattr(comparison, "NOEL_NEGATIVE_CONTROL", "absent_control")
    elif fault == "ledger":
        monkeypatch.setattr(comparison, "load_profiles", unreadable)
    else:
        monkeypatch.setattr(comparison, "retained_profile_features", unreadable)
    output = tmp_path / "comparison"
    with pytest.raises(ValueError, match=message):
        execute_comparison_run(stage["core"], stage["noel"], stage["jed"], stage["ilyas"], output)
    assert not output.exists()


def _noel_status(path: Path, table1: dict[str, str], control: dict[str, str]) -> Path:
    """A calibration status whose overall statuses follow Noel's rule for these criteria."""
    from mtp.calibration import calibration_status

    rows = [
        {"procedure": procedure, "criterion": criterion, "verdict": verdict}
        for procedure, criteria in (
            ("fraction_sign_test", table1),
            ("screened_fraction_sign_test", control),
        )
        for criterion, verdict in criteria.items()
    ]
    path.write_text(json.dumps({"calibration_status": calibration_status(rows)}))
    return path


#: The production procedure's criteria as Noel's harness reports a clean run of it.
CLEAN = {"per_test_level": "pass (sanity check)", "family_error": "pass", "power_reported": "pass"}
#: The negative control failing, as the protocol declares and as on the Tier 1+2 run.
CONTROL_FAILS = {
    "per_test_level": "fail",
    "family_error": "fail",
    "selection_inflation_measured": "reported",
}
CONTROL_PASSES = {
    "per_test_level": "pass",
    "family_error": "pass",
    "selection_inflation_measured": "reported",
}


def test_clean_production_with_the_expected_control_failure(tmp_path: Path) -> None:
    status = _noel_status(tmp_path / "calibration_status.json", CLEAN, CONTROL_FAILS)
    verdict = noel_calibration_verdict(status, TABLE1_SPEC)
    # The control failing as declared says nothing against the production procedure (#855).
    assert verdict.outcome != CalibrationOutcome.anti_conservative
    assert verdict.null_control_passed is True
    assert verdict.replicate_budget_sufficient is True
    # Power has no threshold in Noel's protocol, so this is the best outcome there (#850).
    assert verdict.power_passed is None
    assert verdict.outcome == CalibrationOutcome.partially_calibrated
    assert verdict.negative_control_procedure == "screened_fraction_sign_test"
    assert verdict.negative_control_outcome == "fail"
    assert verdict.negative_control_as_expected is True
    assert verdict.selection_uncertainty is not None
    assert "selected on this corpus" in verdict.selection_uncertainty
    assert "p-values are exploratory" in verdict.selection_uncertainty
    assert "does not validate the post-hoc six-trait analysis" in verdict.selection_uncertainty
    # The qualifier travels with the summary too, so no reader sees the label without it.
    assert verdict.selection_uncertainty in verdict.summary_message
    assert "fraction_sign_test pass [" in verdict.summary_message
    assert "screened_fraction_sign_test fail [family_error=fail" in verdict.summary_message
    assert "selecting and then testing exceeds the tolerance" in verdict.summary_message


@pytest.mark.parametrize(
    "control",
    [CONTROL_PASSES, {**CONTROL_FAILS, "family_error": "pass"}],
    ids=["control_passes", "control_passes_family_error"],
)
def test_an_unexpected_control_pass_withholds_certification(
    tmp_path: Path, control: dict[str, str]
) -> None:
    status = _noel_status(tmp_path / "calibration_status.json", CLEAN, control)
    verdict = noel_calibration_verdict(status, TABLE1_SPEC)
    assert verdict.null_control_passed is None
    assert verdict.outcome != CalibrationOutcome.calibrated
    assert verdict.negative_control_as_expected is False
    assert "passed family_error, which the protocol declares it must fail" in (
        verdict.summary_message
    )
    assert "certification is withheld" in verdict.summary_message


def test_an_actual_production_null_failure_is_anti_conservative(tmp_path: Path) -> None:
    table1 = {**CLEAN, "family_error": "fail"}
    status = _noel_status(tmp_path / "calibration_status.json", table1, CONTROL_FAILS)
    verdict = noel_calibration_verdict(status, TABLE1_SPEC)
    assert verdict.outcome == CalibrationOutcome.anti_conservative
    assert verdict.null_control_passed is False
    assert verdict.negative_control_as_expected is True


def test_labelling_the_selection_is_one_switch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The owner may label the selected traits anti-conservative; that is the only change.
    monkeypatch.setattr(comparison, "SELECTION_LABELS_ANTI_CONSERVATIVE", True)
    status = _noel_status(tmp_path / "calibration_status.json", CLEAN, CONTROL_FAILS)
    verdict = noel_calibration_verdict(status, TABLE1_SPEC)
    assert verdict.outcome == CalibrationOutcome.anti_conservative
    assert "The label counts the selection of the features" in verdict.summary_message


def test_a_failed_noel_criterion_is_recorded_whatever_the_budget(tmp_path: Path) -> None:
    table1 = {**CLEAN, "family_error": "fail", "power_reported": "insufficient_budget"}
    control = {**CONTROL_FAILS, "family_error": "insufficient_budget"}
    status = _noel_status(tmp_path / "calibration_status.json", table1, control)
    verdict = noel_calibration_verdict(status, TABLE1_SPEC)
    assert verdict.outcome == CalibrationOutcome.anti_conservative
    assert verdict.null_control_passed is False
    assert verdict.replicate_budget_sufficient is False


def test_a_failed_power_criterion_is_failed_power(tmp_path: Path) -> None:
    status = _noel_status(
        tmp_path / "calibration_status.json", {**CLEAN, "power_reported": "fail"}, CONTROL_FAILS
    )
    verdict = noel_calibration_verdict(status, TABLE1_SPEC)
    assert verdict.outcome == CalibrationOutcome.null_controlled_low_power
    assert verdict.null_control_passed is True
    assert verdict.power_passed is False
    assert "fraction_sign_test fail [" in verdict.summary_message


@pytest.mark.parametrize("unevaluated", ["insufficient_budget", "not_evaluated"])
@pytest.mark.parametrize("procedure", ["fraction_sign_test", "screened_fraction_sign_test"])
def test_a_not_certified_procedure_is_an_insufficient_budget(
    tmp_path: Path, procedure: str, unevaluated: str
) -> None:
    table1, control = CLEAN, CONTROL_FAILS
    if procedure == "fraction_sign_test":
        table1 = {**CLEAN, "family_error": unevaluated}
    else:
        control = {**CONTROL_FAILS, "per_test_level": unevaluated, "family_error": unevaluated}
    status = _noel_status(tmp_path / "calibration_status.json", table1, control)
    verdict = noel_calibration_verdict(status, TABLE1_SPEC)
    assert verdict.outcome == CalibrationOutcome.insufficient_budget
    assert verdict.replicate_budget_sufficient is False
    assert verdict.null_control_passed is None
    assert f"{procedure} not_certified" in verdict.summary_message
    if procedure == "screened_fraction_sign_test":
        assert verdict.negative_control_as_expected is None


@pytest.mark.parametrize(
    ("procedure", "overall"),
    [("fraction_sign_test", "pass"), ("screened_fraction_sign_test", "not_certified")],
)
def test_an_overall_status_its_criteria_contradict_is_an_error(
    tmp_path: Path, procedure: str, overall: str
) -> None:
    path = _noel_status(
        tmp_path / "calibration_status.json", {**CLEAN, "family_error": "fail"}, CONTROL_FAILS
    )
    content = json.loads(path.read_text())
    content["calibration_status"][procedure]["overall"] = overall
    path.write_text(json.dumps(content))
    with pytest.raises(ValueError, match=f"records '{procedure}' as '{overall}'"):
        noel_calibration_verdict(path, TABLE1_SPEC)


@pytest.mark.parametrize("procedure", ["fraction_sign_test", "screened_fraction_sign_test"])
def test_a_status_without_either_procedure_is_an_error(tmp_path: Path, procedure: str) -> None:
    path = _noel_status(tmp_path / "calibration_status.json", CLEAN, CONTROL_FAILS)
    content = json.loads(path.read_text())
    del content["calibration_status"][procedure]
    path.write_text(json.dumps(content))
    with pytest.raises(ValueError, match=procedure):
        noel_calibration_verdict(path, TABLE1_SPEC)


COMPOSITION_COLUMNS = (
    "seed",
    "stream",
    "unit",
    "replicate",
    "correction",
    "threshold",
    "n_tests",
    "significant",
)


def _composition_rows(replicates: int) -> list[dict[str, object]]:
    """Rows in the producer's order: each threshold, then BH and BY, for each replicate."""
    return [
        {
            "seed": 182,
            "stream": "complete_null_study_flip",
            "unit": "study",
            "replicate": replicate,
            "correction": correction,
            "threshold": threshold,
            "n_tests": 317,
            "significant": int(replicate % 7 == 0 and correction == "BH"),
        }
        for replicate in range(replicates)
        for threshold in (0.05, 0.01, 0.001)
        for correction in ("BH", "BY")
    ]


def _composition_summary(replicates: int) -> dict[str, object]:
    called = len(range(0, replicates, 7))
    return {
        "seed": 182,
        "stream": "complete_null_study_flip",
        "unit": "study",
        "n_traits": 317,
        "replicates_requested": replicates,
        "replicates_completed": replicates,
        "summary": [
            {
                "correction": correction,
                "threshold": threshold,
                "replicates": replicates,
                "replicates_with_any_call": called if correction == "BH" else 0,
            }
            for correction in ("BH", "BY")
            for threshold in (0.05, 0.01, 0.001)
        ],
    }


def _composition_verdict(
    root: Path, rows: list[dict[str, object]], summary: dict[str, object]
) -> CalibrationVerdict:
    summary_path = root / "composition_calibration_summary.json"
    summary_path.write_text(json.dumps(summary))
    table = root / "composition_calibration_replicates.tsv"
    lines = ["\t".join(COMPOSITION_COLUMNS)]
    lines += ["\t".join(str(row[column]) for column in COMPOSITION_COLUMNS) for row in rows]
    table.write_text("\n".join(lines) + "\n")
    return jed_calibration_verdict(
        summary_path,
        table,
        STANDARD_ANALYSIS_SPECS["jed_table1_paired_log_ratio_sign_test_v1"],
        STANDARD_ANALYSIS_SPECS["jed_composition_arcsine_sqrt_ttest_calibration_v1"],
    )


def test_a_complete_composition_table_meets_the_budget(tmp_path: Path) -> None:
    verdict = _composition_verdict(tmp_path, _composition_rows(100), _composition_summary(100))
    assert verdict.replicate_budget_sufficient is True
    assert verdict.null_control_passed is None
    assert verdict.negative_control_procedure is None
    assert "100 of 100 replicates completed" in verdict.summary_message


@pytest.mark.parametrize(
    ("case", "message"),
    [
        ("final_replicate_truncated", "truncated"),
        ("duplicate_cell", "duplicate row"),
        ("replicate_out_of_range", "outside 0..99"),
        ("other_stream", "inconsistent"),
        ("summary_cell_count", "per-cell replicate counts disagree"),
        ("summary_completed", "records 100 complete"),
    ],
)
def test_an_incomplete_or_inconsistent_composition_table_is_an_error(
    tmp_path: Path, case: str, message: str
) -> None:
    rows, summary = _composition_rows(100), _composition_summary(100)
    if case == "final_replicate_truncated":
        # The producer's last five rows: replicate 99 keeps only its first cell (#856).
        rows = rows[:-5]
    elif case == "duplicate_cell":
        rows = [*rows, rows[0]]
    elif case == "replicate_out_of_range":
        rows = [{**row, "replicate": 100} if row["replicate"] == 99 else row for row in rows]
    elif case == "other_stream":
        rows = [*rows[:3], {**rows[3], "stream": "complete_null_experiment_flip"}, *rows[4:]]
    elif case == "summary_cell_count":
        groups = [{"correction": "BH", "threshold": 0.05, "replicates": 99}]
        summary = {**summary, "summary": groups}
    else:
        summary = {**summary, "replicates_completed": 99}
    with pytest.raises(ValueError, match=message):
        _composition_verdict(tmp_path, rows, summary)


def _ledger(included: list[str], excluded: int) -> EligibilityLedger:
    return build_eligibility_ledger(
        ledger_id=f"ledger:{TABLE1_SPEC.id}",
        analysis_spec_id=TABLE1_SPEC.id,
        contrast_population=ContrastPopulation.all_paired,
        included_ids=included,
        exclusion_counts={"zero_retained_total": excluded},
    )


@pytest.mark.parametrize(
    "exclusion",
    [
        None,
        {"studies_before": 3, "studies_after": 2, "studies_lost": 1},
        {"study_ids_after": ["bsdb:1", "bsdb:2"], "studies_after": 2, "studies_lost": 1},
        {
            "study_ids_after": ["bsdb:1", "bsdb:2"],
            "studies_before": 4,
            "studies_after": 2,
            "studies_lost": 1,
        },
        {
            "study_ids_after": ["bsdb:1", "bsdb:3"],
            "studies_before": 3,
            "studies_after": 2,
            "studies_lost": 1,
        },
    ],
)
def test_a_ledger_must_reconcile_with_the_recorded_exclusion(
    exclusion: dict[str, object] | None,
) -> None:
    ledger = _ledger(["bsdb:1", "bsdb:2"], 1)
    provenance = {} if exclusion is None else {"pair_exclusion": exclusion}
    with pytest.raises(ValueError, match="pair_exclusion"):
        reconcile_ledger(ledger, provenance)
    complete = {
        "study_ids_after": ["bsdb:1", "bsdb:2"],
        "studies_before": 3,
        "studies_after": 2,
        "studies_lost": 1,
    }
    reconcile_ledger(ledger, {"pair_exclusion": complete})


def test_a_ledger_needs_a_study_level_spec(stage: dict[str, Path]) -> None:
    pooled = STANDARD_ANALYSIS_SPECS["jed_pooled_contrast_fraction_diff_v1"]
    with pytest.raises(ValueError, match="paired_experiment_contrast"):
        shared_basis_ledger(pooled, verify_core_run(stage["core"]), ContrastPopulation.all_paired)


def test_an_input_no_step_recorded_is_an_error(stage: dict[str, Path]) -> None:
    run: ArmRun = verify_arm(stage["noel"])
    with pytest.raises(ValueError, match="records no 'taxon_influence' step"):
        recorded_output(run, "taxon_influence", "taxon_influence_status.json")
    with pytest.raises(ValueError, match=r"0 files named 'coverage_status\.json'"):
        recorded_output(run, "contrast_orientation", "coverage_status.json")


Bases = list[dict[str, object]]
COUNT_BASIS, ORA_BASIS = "basis:count_profiles_non_mediadive", "basis:ora_profiles_binary"


def _comparison_copy(
    stage: dict[str, Path], tmp_path: Path, edit: Callable[[Bases], Bases] | None = None
) -> Path:
    """A copy of the verified comparison run, with its recorded analysis bases edited."""
    copy = tmp_path / "comparison"
    shutil.copytree(stage["comparison"], copy)
    manifest = copy / "comparison_run.json"
    record = json.loads(manifest.read_text(encoding="utf-8"))
    record["manifest_path"] = str(manifest.resolve())
    bases = record["analysis_bases"]
    assert [basis["id"] for basis in bases] == [COUNT_BASIS, ORA_BASIS]
    if edit is not None:
        record["analysis_bases"] = edit(bases)
    manifest.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return copy


def _replace(index: int, field: str, value: str | int) -> Callable[[Bases], Bases]:
    """Replace one field of one recorded basis, which must change its value."""

    def edit(bases: Bases) -> Bases:
        assert bases[index][field] != value
        return [{**basis, field: value} if at == index else basis for at, basis in enumerate(bases)]

    return edit


def _as_merged(basis: dict[str, object]) -> dict[str, object]:
    """A basis as a bad merge of #828 and #829 could bring it back."""
    return {
        **basis,
        "retained_feature_count": 317,
        "rollup_policy": "ncbi_taxonomy_closure_rollup_off",
    }


def test_an_unchanged_copy_of_the_comparison_verifies(
    stage: dict[str, Path], tmp_path: Path
) -> None:
    """The copy the tampering tests start from is itself a valid comparison run."""
    copy = _comparison_copy(stage, tmp_path)
    assert verify_comparison_run(copy).analysis_bases == _comparison(stage).analysis_bases


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("retained_feature_count", 317),
        ("rollup_policy", "ncbi_taxonomy_closure_rollup_off"),
        ("feature_filter_id", "sha256:" + "0" * 64),
        ("feature_catalog_sha256", "0" * 64),
    ],
)
def test_verification_rederives_the_analysis_bases_from_the_core(
    stage: dict[str, Path], tmp_path: Path, field: str, value: str | int
) -> None:
    """A record whose bases its verified core does not derive fails: another trait count,
    scope, retained-feature identity or catalog (#817, #860)."""
    copy = _comparison_copy(stage, tmp_path, _replace(0, field, value))
    with pytest.raises(ValueError, match=f"{COUNT_BASIS} {field}"):
        verify_comparison_run(copy)


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        pytest.param(
            lambda bases: [_as_merged(bases[0]), _as_merged(bases[1]), *bases],
            "records analysis bases",
            id="both_sides_of_a_merge",
        ),
        pytest.param(lambda bases: [*bases, *bases], "records analysis bases", id="repeated"),
        pytest.param(lambda bases: bases[::-1], "records analysis bases", id="swapped"),
        pytest.param(lambda bases: bases[:1], "records analysis bases", id="missing_ora"),
        pytest.param(
            lambda bases: [*bases, {**bases[1], "id": "basis:extra"}],
            "records analysis bases",
            id="extra",
        ),
        pytest.param(
            _replace(1, "retained_feature_count", 317),
            f"{ORA_BASIS} retained_feature_count",
            id="ora_field",
        ),
    ],
)
def test_verification_refuses_any_other_list_of_bases(
    stage: dict[str, Path], tmp_path: Path, edit: Callable[[Bases], Bases], message: str
) -> None:
    """The recorded bases must be exactly the derived ones, in order, each once (#941, #942)."""
    copy = _comparison_copy(stage, tmp_path, edit)
    with pytest.raises(ValueError, match=message):
        verify_comparison_run(copy)
