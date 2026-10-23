"""Issue #853: independent inventories and disk-backed full-family corrections."""

import json
import runpy
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path

import pytest

from mtp.datamodel.multiple_testing_pydantic import HypothesisCoordinates, HypothesisFamilyMetadata
from mtp.datamodel.workflow_pydantic import AnalysisBasis
from mtp.multiple_testing import canonical_sha256, correct_family, load_policy
from mtp.multiplicity import sources as source_module
from mtp.multiplicity import stage as stage_module
from mtp.multiplicity.families import (
    FamilyInput,
    jed_family,
    metadata,
    noel_family,
    observation,
    rows,
)
from mtp.multiplicity.inventory import correct_inventory, read_corrected_inventory
from mtp.multiplicity.native import native_family
from mtp.multiplicity.ora import ora_family
from mtp.multiplicity.sources import (
    VerifiedSources,
    recorded_arm_output,
    recorded_comparison_output,
    retained_features,
    validate_basis,
)
from mtp.multiplicity.stage import run_post_arm, verify_post_arm
from mtp.paths import file_sha256
from mtp.workflow import runtime as workflow_runtime
from mtp.workflow.coverage import feature_filter_id
from multiplicity_fixture import fixture_sources, write_table
from test_multiple_testing import family, raw


def test_adapter_code_identity_binds_the_shared_serializer(monkeypatch: pytest.MonkeyPatch) -> None:
    before = stage_module.code_sha256()
    serializer = Path(workflow_runtime.__file__)

    def changed_digest(path: Path) -> str:
        return "f" * 64 if path == serializer else file_sha256(path)

    monkeypatch.setattr(stage_module, "file_sha256", changed_digest)
    assert stage_module.code_sha256() != before


@pytest.mark.parametrize("kind", ["code", "policy"])
def test_post_arm_refuses_identity_changes_during_execution(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    """The final manifest cannot name new disk code after running cached old code."""
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    fingerprint = "a" * 64
    policy = load_policy()
    monkeypatch.setattr(stage_module, "code_sha256", lambda: fingerprint)
    monkeypatch.setattr(stage_module, "load_policy", lambda: policy)
    original_report = stage_module.migration_report

    def report_then_edit_code(
        family: FamilyInput, corrected: Path, destination: Path
    ) -> dict[str, int]:
        nonlocal fingerprint
        counts = original_report(family, corrected, destination)
        if kind == "code":
            fingerprint = "b" * 64
        else:
            policy.thresholds.append(0.02)
        return counts

    monkeypatch.setattr(stage_module, "migration_report", report_then_edit_code)
    destination = tmp_path / "sidecars"
    with pytest.raises(ValueError, match=rf"{kind}.*changed"):
        run_post_arm(Path(sources.comparison.manifest_path), destination)
    assert not (destination / "multiplicity_run.json").exists()


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_adapter_membership_hash_rejects_nonfinite_numbers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, value: float
) -> None:
    sources = fixture_sources(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="Out of range float"):
        metadata(
            sources,
            family_id="fixture",
            analysis_id="fixture",
            product="count_profiles",
            count=1,
            membership={"eligibility": [{"invalid_number": value}]},
            unit="trait",
            inference="study",
            estimand="fixture",
            selection="all fixture traits",
            history=["fixture"],
            paths=[],
        )


def test_disk_correction_matches_foundation_without_changing_family(tmp_path: Path) -> None:
    declared = family(4)
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    values = raw([0.01, 0.04, 0.04, None])
    instant = datetime(2026, 10, 23, tzinfo=UTC)
    expected = correct_family(declared, values, recorded_at=instant)
    actual = correct_inventory(
        metadata, iter(declared.members), iter(values), tmp_path, recorded_at=instant
    )
    assert actual.family_sha256 == expected.family_sha256
    assert actual.input_sha256 == expected.input_sha256
    assert actual.output_sha256 == expected.output_sha256
    assert list(read_corrected_inventory(tmp_path / actual.corrected.identity)) == expected.results


def test_disk_inventory_rejects_missing_and_duplicate_rows(tmp_path: Path) -> None:
    declared = family(2)
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    values = raw([0.01, 0.2])
    with pytest.raises(ValueError, match="membership"):
        correct_inventory(metadata, iter(declared.members), iter(values[:1]), tmp_path / "missing")
    with pytest.raises(ValueError, match="duplicate"):
        correct_inventory(
            metadata, iter(declared.members), iter([values[0], values[0]]), tmp_path / "duplicate"
        )


def test_disk_inventory_rejects_renamed_scientific_duplicates(tmp_path: Path) -> None:
    declared = family(2)
    declared.members[1] = declared.members[0].model_copy(
        update={"hypothesis_id": "h1", "eligibility_sha256": "2" * 64}
    )
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    with pytest.raises(ValueError, match="duplicate scientific hypothesis"):
        correct_inventory(metadata, iter(declared.members), iter(raw([0.1, 0.1])), tmp_path)


def test_disk_inventory_rejects_negative_zero(tmp_path: Path) -> None:
    declared = family(1)
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    with pytest.raises(ValueError, match="finite probability"):
        correct_inventory(metadata, iter(declared.members), iter(raw([-0.0])), tmp_path)


def test_disk_inventory_distinguishes_profiles_and_canonicalizes_conditions(tmp_path: Path) -> None:
    declared = family(2)
    for index, member in enumerate(declared.members):
        member.feature_ids = ["same-feature"]
        member.coordinates = HypothesisCoordinates(
            profile_id=f"profile-{index}", condition_ids=["condition-b", "condition-a"]
        )
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    values = raw([0.1, 0.2])
    actual = correct_inventory(metadata, iter(declared.members), iter(values), tmp_path)
    expected = correct_family(declared, values, recorded_at=actual.recorded_at)
    assert actual.family_sha256 == expected.family_sha256
    assert actual.output_sha256 == expected.output_sha256


def test_two_tails_are_one_family_not_two_rescaled_exports(tmp_path: Path) -> None:
    declared = family(4)
    for index, member in enumerate(declared.members):
        member.tail = "upper" if index < 2 else "lower"
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    corrected = correct_inventory(
        metadata, iter(declared.members), iter(raw([0.01, 0.9, 0.02, 0.8])), tmp_path
    )
    rows = list(read_corrected_inventory(tmp_path / corrected.corrected.identity))
    assert [row.q_bh for row in rows] == pytest.approx([0.04, 0.9, 0.04, 0.9])
    assert [row.q_by for row in rows] == pytest.approx([1 / 12, 1.0, 1 / 12, 1.0])


@pytest.mark.parametrize("values", [[], [None], [None, None], [0.0, 1.0]])
def test_disk_boundary_families(tmp_path: Path, values: list[float | None]) -> None:
    declared = family(len(values))
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    result = correct_inventory(metadata, iter(declared.members), iter(raw(values)), tmp_path)
    assert result.tested + result.untestable == len(values)
    expected = correct_family(declared, raw(values), recorded_at=result.recorded_at)
    assert result.output_sha256 == expected.output_sha256


def test_disk_order_independence_and_duplicate_declarations(tmp_path: Path) -> None:
    declared = family(2)
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    instant = datetime(2026, 10, 23, tzinfo=UTC)
    before = correct_inventory(
        metadata,
        iter(declared.members),
        iter(raw([0.01, 0.2])),
        tmp_path / "a",
        recorded_at=instant,
    )
    after = correct_inventory(
        metadata,
        iter(reversed(declared.members)),
        iter(reversed(raw([0.01, 0.2]))),
        tmp_path / "b",
        recorded_at=instant,
    )
    assert before == after
    with pytest.raises(ValueError, match="duplicate"):
        correct_inventory(
            metadata,
            iter([declared.members[0]] * 2),
            iter(raw([0.01, 0.2])),
            tmp_path / "duplicate",
        )


def make_basis(
    catalog: Path, features: tuple[str, ...], product: str = "count_profiles", scope: str = "off"
) -> AnalysisBasis:
    return AnalysisBasis(
        id=f"basis:{product}",
        feature_catalog_sha256=file_sha256(catalog),
        feature_filter_id=feature_filter_id(features),
        retained_feature_count=len(features),
        profile_unit=product,
        aggregation_policy="fixture",
        denominator_policy="fixture",
        background_inventory_id="fixture",
        weighting_policy="fixture",
        rollup_policy=f"ncbi_taxonomy_closure_rollup_{scope}",
    )


def test_basis_uses_profile_order_and_excludes_ora_metadata(tmp_path: Path) -> None:
    profile, catalog = tmp_path / "ora.tsv", tmp_path / "catalog.tsv"
    profile.write_text(
        "study_id\texperiment_id\tinc_dec\ttaxa_count\tB\tA\tmediadive.medium:X\ns\te\tinc\t2\t1\t0\t0\n"
    )
    catalog.write_text(
        "column_id\tfeature_class\nA\ttrait_annotation\nB\ttrait_annotation\nmediadive.medium:X\tculture_medium\n"
    )
    features = retained_features(profile, catalog)
    assert features == ("B", "A")
    basis = make_basis(catalog, features, "ora_profiles", "phenotype")
    validate_basis(
        basis, product="ora_profiles", features=features, catalog=catalog, scope="phenotype"
    )
    for name, value, message in (
        ("retained_feature_count", 317, "count"),
        ("feature_filter_id", "sha256:" + file_sha256(catalog), "identity"),
        ("feature_filter_id", feature_filter_id(tuple(reversed(features))), "identity"),
        ("feature_catalog_sha256", "a" * 64, "catalog"),
        ("rollup_policy", "ncbi_taxonomy_closure_rollup_off", "scope"),
        ("profile_unit", "count_profiles", "unit"),
    ):
        altered = basis.model_copy(update={name: value})
        with pytest.raises(ValueError, match=message):
            validate_basis(
                altered,
                product="ora_profiles",
                features=features,
                catalog=catalog,
                scope="phenotype",
            )


@pytest.mark.parametrize(
    ("status", "p"),
    [("tested", ""), ("all_tied", ".1"), ("tested", "nan"), ("tested", "-0.0"), ("made_up", "")],
)
def test_adapter_never_coerces_undefined_or_invalid_tests(status: str, p: str) -> None:
    with pytest.raises(ValueError, match=r"status|probability"):
        observation(
            "h",
            {"p": p, "q": "0.2", "status": status},
            p="p",
            q="q",
            status="status",
            undefined={"all_tied"},
            method="BH",
            native_family="fixture",
        )


@pytest.mark.parametrize("arm", ["noel", "jed", "native"])
def test_untestable_reasons_must_match_producer_support(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, arm: str
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    if arm == "native":
        path = recorded_comparison_output(
            sources.comparison, "trait_profile_paired_contrast_trait_stats_pooled.tsv"
        )
        status, p, q = "study_test_status", "study_p_value", "study_q_value"
        family_input = native_family(sources)
    else:
        step, filename, status, p, q = {
            "noel": (
                "directional_screen",
                "trait_direction_screen.tsv",
                "test_status",
                "binomial_p",
                "by_q_nominal",
            ),
            "jed": ("composition", "trait_log_ratio_stats.tsv", "status", "sign_p", "by_q"),
        }[arm]
        path = recorded_arm_output(sources.arms[arm], step, filename)
        family_input = noel_family(sources) if arm == "noel" else jed_family(sources)
    records = list(rows(path))
    if arm == "jed":
        selected = records[0]
        reason = "no_positive_pairs" if int(selected["studies_eligible"]) else "all_tied"
    else:
        selected = next(row for row in records if row[status] == "tested")
        reason = "all_tied"
    selected[status], selected[p], selected[q] = reason, "", "1" if arm == "native" else ""
    write_table(path, records)
    with pytest.raises(ValueError, match=r"status|support|reason"):
        list(family_input.observations())


def test_undefined_native_q_must_follow_producer_convention() -> None:
    with pytest.raises(ValueError, match="status"):
        observation(
            "h",
            {"p": "", "q": ".2", "status": "all_tied"},
            p="p",
            q="q",
            status="status",
            undefined={"all_tied"},
            method="BY",
            native_family="noel-full-screen",
        )


def test_ratio_eligibility_cannot_drop_a_core_study(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    path = recorded_arm_output(sources.arms["jed"], "composition", "trait_log_ratio_studies.tsv")
    records = list(rows(path))
    write_table(path, records[1:])
    with pytest.raises(ValueError, match="eligibility"):
        jed_family(sources)


def test_screen_support_must_match_independently_loaded_core_cohort(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    family_input = noel_family(sources)
    path = recorded_arm_output(
        sources.arms["noel"], "directional_screen", "trait_direction_screen.tsv"
    )
    records = list(rows(path))
    selected = next(row for row in records if row["test_status"] == "tested")
    for key in ("studies_up", "studies_scored"):
        selected[key] = str(int(selected[key]) + 1)
    write_table(path, records)
    with pytest.raises(ValueError, match="cohort"):
        list(family_input.observations())


@pytest.mark.parametrize("subset", [False, True])
def test_all_adapters_replay_full_and_oriented_subset_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, subset: bool
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch, subset=subset)
    for product, basis in sources.bases.items():
        validate_basis(
            basis,
            product=product,
            features=sources.features[product],
            catalog=sources.product("feature_catalog"),
            scope=sources.scope,
        )
    result = run_post_arm(Path(sources.comparison.manifest_path), tmp_path / "sidecars")
    assert verify_post_arm(tmp_path / "sidecars") == result
    assert result.strain_rollup == ("phenotype" if subset else "off")
    for family_record in result.families:
        family_metadata = family_record.metadata
        basis = sources.bases[family_metadata.profile_unit]
        assert family_metadata.strain_rollup == sources.scope
        assert family_metadata.feature_filter_id == basis.feature_filter_id
    assert result.families[0].metadata.retained_feature_count == (6 if subset else 7)
    assert result.families[0].metadata.source_membership_sha256 == canonical_sha256(
        {
            "features": sources.features["count_profiles"],
            "positive_total_paired_profiles": [
                (f"s{i}", "e1", direction) for i in range(8) for direction in ("dec", "inc")
            ],
        }
    )
    assert result.families[1].metadata.expected_family_size == 15
    ratio_definitions = list(jed_family(sources).declarations())
    assert {row.raw.status for row in jed_family(sources).observations()} == {
        "tested",
        "untestable",
    }
    assert all(len(row.feature_ids) == 2 for row in ratio_definitions)
    assert all(
        feature in sources.features["count_profiles"]
        for row in ratio_definitions
        for feature in row.feature_ids
    )
    assert result.families[2].metadata.expected_family_size == 6
    assert result.families[3].metadata.expected_family_size == 2 * 16 * (6 if subset else 7)
    assert {
        tuple(json.loads(member.coordinates.profile_id or "null"))
        for member in ora_family(sources).declarations()
    } == {(f"s{i}", "e1", direction) for i in range(8) for direction in ("inc", "dec")}
    assert all(
        not Path(artifact.identity).is_absolute()
        for family in result.families
        for artifact in family.metadata.sources
    )
    rows_path = tmp_path / "sidecars" / result.families[3].corrected.identity
    rows = list(read_corrected_inventory(rows_path))
    assert {row.p_value for row in rows} == {0.1, 0.9}
    assert {row.q_bh for row in rows} == {0.2, 0.9}


@pytest.mark.parametrize(
    "kind", ["q", "metadata", "policy", "source", "extra", "numpy", "scope", "filter"]
)
def test_strict_post_arm_verifier_rejects_tampering(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    destination = tmp_path / "sidecars"
    run_post_arm(Path(sources.comparison.manifest_path), destination)
    manifest = destination / "multiplicity_run.json"
    document = json.loads(manifest.read_text())
    if kind == "q":
        path = destination / document["families"][0]["corrected"]["identity"]
        records = [json.loads(line) for line in path.read_text().splitlines()]
        records[0]["q_by"] = 0.0001
        path.write_text("".join(json.dumps(record) + "\n" for record in records))
        document["families"][0]["corrected"]["sha256"] = file_sha256(path)
    elif kind == "metadata":
        document["families"][0]["metadata"]["estimand"] = "different biological claim"
    elif kind == "policy":
        document["policy"]["primary_method"] = "BH"
    elif kind == "source":
        path = sources.product("count_profiles")
        path.write_text(path.read_text() + "\n")
    elif kind == "numpy":
        document["numpy_version"] = "unsupported-array-implementation"
    elif kind == "scope":
        document["families"][0]["metadata"]["strain_rollup"] = "phenotype"
    elif kind == "filter":
        document["families"][0]["metadata"]["feature_filter_id"] = "sha256:" + "1" * 64
    else:
        (destination / "unrecorded.txt").write_text("unexpected")
    manifest.write_text(json.dumps(document))
    with pytest.raises(ValueError, match=r"replay|policy|metadata|inventory"):
        verify_post_arm(destination)


@pytest.mark.parametrize("adapter", [noel_family, jed_family, native_family, ora_family])
def test_each_adapter_rejects_truncated_family(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    adapter: Callable[[VerifiedSources], FamilyInput],
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    # Full declared inventory is independent of imported result membership.
    family_input = adapter(sources)
    values = list(family_input.observations())
    with pytest.raises(ValueError, match="membership"):
        correct_inventory(
            family_input.metadata,
            family_input.declarations(),
            iter(observation.raw for observation in values[1:]),
            tmp_path / "truncated",
        )


def test_ora_adapter_consumes_actual_both_tail_producer_exports(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    profiles = sources.product("ora_profiles")
    records = list(rows(profiles))
    for record in records:
        # Both fixture background taxa are drawn; exactly one has each trait.
        for feature in sources.features["ora_profiles"]:
            record[feature] = "1"
    write_table(profiles, records)
    producer = runpy.run_path(str(Path(__file__).parents[1] / "scripts/trait_enrichment.py"))[
        "main"
    ]
    for tail, stage_name in (("upper", "enrichment"), ("lower", "depletion")):
        export = recorded_arm_output(
            sources.arms["ilyas"], stage_name, f"trait_{stage_name}_all_results.tsv"
        )
        producer(
            [
                "--input",
                str(profiles),
                "--taxa-matrix",
                str(sources.product("taxa_matrix")),
                "--output",
                str(export.with_suffix(".txt")),
                "--all-results-output",
                str(export),
                "--tail",
                tail,
                "--max-fold",
                "0.5",
            ]
        )
    declared = ora_family(sources)
    corrected = correct_inventory(
        declared.metadata,
        declared.declarations(),
        (observation.raw for observation in declared.observations()),
        tmp_path / "corrected",
    )
    assert corrected.tested == 2 * 16 * len(sources.features["ora_profiles"])
    assert {
        row.p_value for row in read_corrected_inventory(tmp_path / "corrected/adjusted.jsonl")
    } == {1.0}


def test_source_loader_checks_actual_basis_after_upstream_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = fixture_sources(tmp_path / "inputs", monkeypatch)
    monkeypatch.setattr(source_module, "verify_comparison_run", lambda _path: fixture.comparison)
    monkeypatch.setattr(source_module, "verify_core_run", lambda _path: fixture.core)
    monkeypatch.setattr(
        source_module,
        "verify_arm",
        lambda path: next(arm for arm in fixture.arms.values() if Path(arm.manifest_path) == path),
    )
    opened = source_module.open_verified_sources(Path(fixture.comparison.manifest_path))
    assert opened.features == fixture.features
    fixture.bases["count_profiles"].feature_filter_id = "sha256:" + file_sha256(
        fixture.product("feature_catalog")
    )
    with pytest.raises(ValueError, match="ordered retained feature identity"):
        source_module.open_verified_sources(Path(fixture.comparison.manifest_path))


def test_source_loader_rejects_mixed_arm_core(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = fixture_sources(tmp_path / "inputs", monkeypatch)
    monkeypatch.setattr(source_module, "verify_comparison_run", lambda _path: fixture.comparison)
    monkeypatch.setattr(source_module, "verify_core_run", lambda _path: fixture.core)
    monkeypatch.setattr(
        source_module,
        "verify_arm",
        lambda path: next(arm for arm in fixture.arms.values() if Path(arm.manifest_path) == path),
    )
    fixture.arms["jed"].core_manifest_sha256 = "b" * 64
    with pytest.raises(ValueError, match="core binding"):
        source_module.open_verified_sources(Path(fixture.comparison.manifest_path))


@pytest.mark.parametrize("product", ["count_profiles", "ora_profiles"])
@pytest.mark.parametrize(
    "field",
    ["aggregation_policy", "denominator_policy", "background_inventory_id", "weighting_policy"],
)
def test_source_loader_rejects_changed_producer_semantics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, product: str, field: str
) -> None:
    fixture = fixture_sources(tmp_path / "inputs", monkeypatch)
    monkeypatch.setattr(source_module, "verify_comparison_run", lambda _path: fixture.comparison)
    monkeypatch.setattr(source_module, "verify_core_run", lambda _path: fixture.core)
    monkeypatch.setattr(
        source_module,
        "verify_arm",
        lambda path: next(arm for arm in fixture.arms.values() if Path(arm.manifest_path) == path),
    )
    setattr(fixture.bases[product], field, "different-scientific-basis")
    with pytest.raises(ValueError, match="producer contract"):
        source_module.open_verified_sources(Path(fixture.comparison.manifest_path))


def test_sidecars_cannot_overwrite_or_live_within_source_bundles(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = fixture_sources(tmp_path / "inputs", monkeypatch)
    with pytest.raises(ValueError, match="disjoint"):
        run_post_arm(
            Path(fixture.comparison.manifest_path),
            Path(fixture.core.manifest_path).parent / "sidecars",
        )
    with pytest.raises(FileExistsError, match="new directory"):
        run_post_arm(Path(fixture.comparison.manifest_path), tmp_path)
