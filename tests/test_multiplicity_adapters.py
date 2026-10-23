"""Issue #853: independent inventories and disk-backed full-family corrections."""

import json
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path

import pytest

from mtp.datamodel.multiple_testing_pydantic import HypothesisFamilyMetadata
from mtp.datamodel.workflow_pydantic import AnalysisBasis
from mtp.multiple_testing import canonical_sha256, correct_family
from mtp.multiplicity import sources as source_module
from mtp.multiplicity.families import FamilyInput, jed_family, noel_family, observation
from mtp.multiplicity.inventory import correct_inventory, read_corrected_inventory
from mtp.multiplicity.native import native_family
from mtp.multiplicity.ora import ora_family
from mtp.multiplicity.sources import VerifiedSources, retained_features, validate_basis
from mtp.multiplicity.stage import run_post_arm, verify_post_arm
from mtp.paths import file_sha256
from mtp.workflow.coverage import feature_filter_id
from multiplicity_fixture import fixture_sources
from test_multiple_testing import family, raw


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
    ("status", "p"), [("tested", ""), ("all_tied", ".1"), ("tested", "nan"), ("made_up", "")]
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
    assert all(len(row.feature_ids) == 2 for row in ratio_definitions)
    assert all(
        feature in sources.features["count_profiles"]
        for row in ratio_definitions
        for feature in row.feature_ids
    )
    assert result.families[2].metadata.expected_family_size == 6
    assert result.families[3].metadata.expected_family_size == 2 * 16 * (6 if subset else 7)
    assert all(
        not Path(artifact.identity).is_absolute()
        for family in result.families
        for artifact in family.metadata.sources
    )
    rows_path = tmp_path / "sidecars" / result.families[3].corrected.identity
    rows = list(read_corrected_inventory(rows_path))
    assert {row.p_value for row in rows} == {0.1, 0.9}
    assert {row.q_bh for row in rows} == {0.2, 0.9}


@pytest.mark.parametrize("kind", ["q", "metadata", "policy", "source", "extra"])
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
