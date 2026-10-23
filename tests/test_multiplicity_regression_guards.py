"""Independent weakening-guard regressions for additive multiplicity (#974)."""

from pathlib import Path

import pytest

from mtp.datamodel.multiple_testing_pydantic import HypothesisFamilyMetadata
from mtp.multiplicity import native, ora, stage
from mtp.multiplicity.families import FamilyInput, NativeObservation, rows, validate_sign_status
from mtp.multiplicity.inventory import correct_inventory
from mtp.multiplicity.sources import recorded_arm_output
from multiplicity_fixture import fixture_sources, write_table
from test_multiple_testing import family, raw


@pytest.mark.parametrize("guard", ["support_floor", "fraction_variance"])
def test_native_selection_guards_are_independent(tmp_path: Path, guard: str) -> None:
    # These deliberately extreme, finite integer counts isolate the zero-variance
    # guard: tiny fractions underflow when squared. They are numerical fixtures,
    # not plausible organism counts or data-derived biological expectations.
    large = 10**170 if guard == "fraction_variance" else 10
    features = ("qualified|a", "qualified|b", "qualified|excluded")
    profile, experiments, labels = (tmp_path / name for name in ("p.tsv", "e.tsv", "l.tsv"))
    records = []
    for index in range(16):
        records.append(
            {
                "study_id": f"s{index // 2}",
                "experiment_id": "e1",
                "inc_dec": "inc" if index % 2 else "dec",
                features[0]: large if index % 2 else 0,
                features[1]: 0 if index % 2 else large,
                features[2]: (large if index < 5 else 0)
                if guard == "support_floor"
                else int(index < 8),
            }
        )
    write_table(profile, records)
    write_table(
        experiments,
        [
            {
                "study_id": f"s{i}",
                "experiment_id": "e1",
                "condition": "fixture",
                "body_site": "Feces",
                "host_species": "Homo sapiens",
            }
            for i in range(8)
        ],
    )
    write_table(labels, [{"column_id": feature, "object_name": feature} for feature in features])
    selected, keys = native.native_candidates(profile, experiments, features, labels)
    assert len(keys) == 16
    assert set(selected) == set(features[:2]), f"Native {guard} admitted an excluded candidate"


@pytest.mark.parametrize("tail", ["upper", "lower"])
def test_ora_background_prevalence_is_checked(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, tail: str
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    prefix = "enrichment" if tail == "upper" else "depletion"
    path = recorded_arm_output(sources.arms["ilyas"], prefix, f"trait_{prefix}_all_results.tsv")
    imported = ora.ora_family(sources)
    records = list(rows(path))
    # Keep metadata and every other value fixed: only K disagrees with the
    # independently read background matrix (one carrier in two taxa).
    records[0]["background_k"] = "2"
    write_table(path, records, ora.ora_comments(path))
    with pytest.raises(ValueError, match="ORA background feature prevalence mismatch"):
        list(imported.observations())


@pytest.mark.parametrize(("tail", "fold"), [("upper", "1.0"), ("lower", "1.0")])
def test_ora_native_hit_needs_fold_as_well_as_q(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, tail: str, fold: str
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    prefix = "enrichment" if tail == "upper" else "depletion"
    path = recorded_arm_output(sources.arms["ilyas"], prefix, f"trait_{prefix}_all_results.tsv")
    imported = ora.ora_family(sources)
    records = list(rows(path))
    records[0].update(p_value="0.0001", fdr_q="0.0001", fold=fold, is_significant="True")
    write_table(path, records, ora.ora_comments(path))
    with pytest.raises(ValueError, match="ORA native hit flag disagrees with q/effect thresholds"):
        list(imported.observations())


def test_sign_eligible_count_reconciles_with_tallies() -> None:
    row = {
        "studies_up": "3",
        "studies_down": "2",
        "studies_tied": "1",
        "studies_scored": "5",
        "studies_eligible": "7",
        "studies_total": "8",
        "studies_excluded": "1",
        "status": "tested",
    }
    with pytest.raises(ValueError, match="eligible-study support disagrees with sign tallies"):
        validate_sign_status(row, status="status", empty="no_eligible_studies")


def test_migration_native_q005_is_not_native_reported_threshold(tmp_path: Path) -> None:
    declared = family(4)
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    observations = [
        NativeObservation(
            value, q, "BH", "native-upper-tail", 0.001, {"is_significant": str(q <= 0.001)}
        )
        for value, q in zip(raw([0.0001] * 4), [0.001, 0.04, 0.05, 0.050001], strict=True)
    ]
    imported = FamilyInput(metadata, lambda: iter(declared.members), lambda: iter(observations), ())
    correction = correct_inventory(
        metadata, declared.members, [item.raw for item in observations], tmp_path / "corrected"
    )
    summary = stage.migration_report(
        imported, tmp_path / "corrected" / correction.corrected.identity, tmp_path / "migration.tsv"
    )
    assert summary["native_reject_005"] == 3, "Native q<=0.05 comparison changed its threshold"
    assert summary["native_reported_hits"] == 1
    assert {row["native_threshold"] for row in rows(tmp_path / "migration.tsv")} == {"0.001"}


def test_final_source_recheck_blocks_published_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    original = stage.migration_report
    count = 0

    def report_then_change_source(
        imported: FamilyInput, corrected: Path, destination: Path
    ) -> dict[str, int]:
        nonlocal count
        result = original(imported, corrected, destination)
        count += 1
        if count == 4:
            experiments = sources.product("experiments")
            experiments.write_text(experiments.read_text() + "\n", encoding="utf-8")
        return result

    monkeypatch.setattr(stage, "migration_report", report_then_change_source)
    output = tmp_path / "sidecar"
    with pytest.raises(ValueError, match=r"Source changed .*sidecars"):
        stage.run_post_arm(Path(sources.comparison.manifest_path), output)
    assert count == 4, "The source change must occur after all scientific import/correction work"
    assert not (output / "multiplicity_run.json").exists()
