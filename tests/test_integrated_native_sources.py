"""Guard native evidence joins against changed tests, populations and families (#551)."""

from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path

import pytest

from mtp.integrated_evidence import (
    PROFILE_SUMMARY,
    SIGNIFICANT,
    TESTED,
    UNDECIDED,
    UNTESTED,
    Record,
    attach_full_screen,
    family_inventory,
    native_composition_records,
    ora_profile_records,
    placeholder_records,
)
from mtp.trait_guilds import TABLE1_TRAITS, Profile

PROFILE_HASH = "profile-fixture-sha256"
EXPERIMENTS_HASH = "experiments-fixture-sha256"
RAW_COLUMNS = {
    "facultatively_aerobic": "METPO:1000608",
    "cellulose_catabolic_process": "METPO:2000103|GO:0030245",
    "coccus_shaped": "METPO:1000668",
    "sphere_shaped": "METPO:1000683",
    "arsenate_reduction": "METPO:2000103|pathways:arsenate_reduction",
    "moderately_halophilic": "METPO:1000623",
}
COLUMNS = {slug: ((TABLE1_TRAITS[slug][0], column),) for slug, column in RAW_COLUMNS.items()}


def write_source(
    path: Path, rows: list[dict[str, str]], metadata: dict[str, str] | None = None
) -> Path:
    """Write a complete miniature native export, including its declared metadata."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        for key, value in (metadata or {}).items():
            handle.write(f"# {key}: {value}\n")
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def composition_rows() -> list[dict[str, str]]:
    """Include one non-Table-1 trait so the native family exceeds the exported subset."""
    return [
        {
            "trait": trait,
            "stratum_kind": "pooled",
            "stratum": "pooled",
            "profile_sha256": PROFILE_HASH,
            "experiments_sha256": EXPERIMENTS_HASH,
            "correction_family_size": "7",
            "correction_method": "BH",
            "test_method": "two_sided_one_sample_t_test",
            "pair_count": "4",
            "study_count": "3",
            "profile_count": "8",
            "n_nonzero_experiments": "3",
            "n_zero_total_experiments": "1",
            "n_zero_total_excluded": "0",
            "zero_total_policy": "historical zero-total pairs retained",
            "study_mean_increased_minus_decreased": "0.2",
            "study_median_increased_minus_decreased": "0.15",
            "study_ci95_lower": "0.1",
            "study_ci95_upper": "0.3",
            "study_p_value": "0.02",
            "study_q_value": "0.02",
            "p_value": "0.00001",
            "q_value": "0.00001",
            "mean_increased_minus_decreased": "-0.7",
        }
        for trait in [*RAW_COLUMNS.values(), "METPO:1000603"]
    ]


def import_composition(path: Path) -> list[Record]:
    return native_composition_records(
        path,
        columns=COLUMNS,
        q_threshold=0.05,
        strict_q_threshold=0.001,
        profile_sha256=PROFILE_HASH,
        experiments_sha256=EXPERIMENTS_HASH,
    )


def test_native_study_test_keeps_pair_test_secondary_and_actual_family(tmp_path: Path) -> None:
    records = import_composition(write_source(tmp_path / "native.tsv", composition_rows()))
    assert len(records) == 6
    assert {record.trait: record.trait_columns for record in records} == RAW_COLUMNS
    assert all(record.p_value == 0.02 and record.q_value == 0.02 for record in records)
    assert all(record.native_pair_q_value == 0.00001 for record in records)
    assert all(record.native_pair_p_value == 0.00001 for record in records)
    assert all(
        record.effect_value == 0.2 and record.study_effect_median == 0.15 for record in records
    )
    assert all(record.interval_low == 0.1 and record.interval_high == 0.3 for record in records)
    assert all(record.testable_units == 3 and record.paired_experiments == 4 for record in records)
    assert all(
        record.test_status == TESTED and record.decision == SIGNIFICANT for record in records
    )
    assert all(record.n_nonzero_experiments == 3 for record in records)
    assert all(record.n_zero_total_experiments == 1 for record in records)
    assert all(record.n_zero_total_excluded == 0 for record in records)
    family = family_inventory(records)[0]
    assert family["family_size"] == 7
    assert family["rows_exported"] == 6
    assert family["family_scope"] == "source_subset"
    assert str(family["family_source"]).startswith("native.tsv@")


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("profile_sha256", "another-profile", "checksum"),
        ("experiments_sha256", "another-experiment-cohort", "experiments checksum"),
        ("correction_family_size", "317", "family size"),
        ("n_nonzero_experiments", "2", "reconcile"),
        ("stratum", "feces", "pooled"),
        ("study_q_value", "nan", "non-finite"),
        ("trait", "METPO:1000662", "lacks"),
        ("correction_method", "Benjamini-Yekutieli", "correction"),
        ("test_method", "Wilcoxon", "test"),
    ],
)
def test_native_rejects_incompatible_sources(
    tmp_path: Path, field: str, value: str, message: str
) -> None:
    rows = composition_rows()
    rows[2][field] = value
    with pytest.raises(ValueError, match=message):
        import_composition(write_source(tmp_path / "native.tsv", rows))


def test_core_native_composition_requires_both_source_hashes(tmp_path: Path) -> None:
    rows = composition_rows()
    for row in rows:
        del row["experiments_sha256"]
    source = write_source(tmp_path / "legacy-native.tsv", rows)
    assert len(import_composition(source)) == 6
    with pytest.raises(ValueError, match="requires profile and experiments hashes"):
        native_composition_records(
            source,
            columns=COLUMNS,
            q_threshold=0.05,
            strict_q_threshold=0.001,
            profile_sha256=PROFILE_HASH,
            experiments_sha256=EXPERIMENTS_HASH,
            require_source_hashes=True,
        )


def test_native_rejects_more_excluded_pairs_than_zero_pairs(tmp_path: Path) -> None:
    rows = composition_rows()
    for row in rows:
        row.update(pair_count="2", study_count="2", profile_count="4", n_zero_total_excluded="2")
    with pytest.raises(ValueError, match="zero-total"):
        import_composition(write_source(tmp_path / "native.tsv", rows))


def test_native_duplicate_trait_is_not_silently_overwritten(tmp_path: Path) -> None:
    rows = composition_rows()
    rows[-1]["trait"] = rows[0]["trait"]
    with pytest.raises(ValueError, match="unique trait"):
        import_composition(write_source(tmp_path / "native.tsv", rows))


def ora_metadata(tail: str = "upper") -> dict[str, str]:
    metadata = {
        "Tail": tail,
        "Testable-Set": "K>0 & N>0",
        "FDR": "True",
        "FDR-Scope": "global",
        "P-Threshold": "0.001",
        "Total-Profiles-Input": "2",
        "Testable-Profiles": "2",
        "Background-Taxa-Universe-M": "100",
        "Background-Mode": "TaxaMatrix (fixture.tsv)",
        "Background-Universe-Definition": "fixture taxa",
        "Total-Tests-Performed": "12",
    }
    metadata["Min-Fold" if tail == "upper" else "Max-Fold"] = "1.5" if tail == "upper" else "0.5"
    return metadata


def ora_rows(*, enriched: bool = True) -> list[dict[str, str]]:
    return [
        {
            "study_id": "study1",
            "experiment_id": "experiment1",
            "inc_dec": direction,
            "trait": trait,
            "fdr_q": "0.0001" if enriched and direction == "inc" else "0.5",
            "fold": "2" if enriched and direction == "inc" else "1",
            "is_significant": str(enriched and direction == "inc"),
            "correction_family": "global",
            "eligibility_status": "testable",
            "background_m": "100",
            "background_k": "10",
            "testable_profiles": "2",
        }
        for trait in RAW_COLUMNS.values()
        for direction in ("inc", "dec")
    ]


def test_ora_summary_preserves_profile_calls_without_fabricated_aggregate_q(tmp_path: Path) -> None:
    enrichment = write_source(tmp_path / "enrichment.tsv", ora_rows(), ora_metadata())
    records = ora_profile_records(enrichment, tmp_path / "missing.tsv", columns=COLUMNS)
    assert len(records) == 12
    for record in records:
        assert record.p_value is None
        assert record.q_value is None
        assert record.decision == UNDECIDED
        assert record.strict_decision == UNDECIDED
        assert record.n_tests == 2
        assert record.n_significant_inc == 1
        assert record.n_significant_dec == 0
        assert record.inc_dec_ratio is None
        assert record.best_profile_q == record.median_profile_q == 0.0001
        assert record.median_hit_fold == record.max_hit_fold == 2
        assert record.n_depleted_inc is None
        assert record.n_depleted_dec is None
        assert record.depletion_source == ""
    assert {r.test_status for r in records} == {PROFILE_SUMMARY, UNTESTED}
    family = family_inventory(records)[0]
    assert family["family_size"] == 12
    assert family["rows_exported"] == 6
    assert family["tested"] == 0
    assert family["with_q_value"] == 0


def test_ora_full_depletion_source_proves_zero_hits(tmp_path: Path) -> None:
    enrichment = write_source(tmp_path / "enrichment.tsv", ora_rows(), ora_metadata())
    depletion = write_source(
        tmp_path / "depletion.tsv", ora_rows(enriched=False), ora_metadata("lower")
    )
    records = ora_profile_records(enrichment, depletion, columns=COLUMNS)
    assert all(r.n_depleted_inc == 0 and r.n_depleted_dec == 0 for r in records)
    assert all(r.depletion_q_threshold == 0.001 and r.depletion_max_fold == 0.5 for r in records)
    assert all(r.depletion_source.startswith("depletion.tsv@") for r in records)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("P-Threshold", "0.05", "production q/fold"),
        ("Min-Fold", "1.1", "production q/fold"),
        ("Tail", "lower", "tail and global FDR"),
        ("FDR-Scope", "profile", "tail and global FDR"),
        ("Total-Tests-Performed", "13", "correction family"),
        ("Testable-Profiles", "4", "profile"),
        ("Background-Mode", "AnnotationRows", "background"),
    ],
)
def test_ora_rejects_inconsistent_metadata(
    tmp_path: Path, field: str, value: str, message: str
) -> None:
    metadata = ora_metadata()
    metadata[field] = value
    enrichment = write_source(tmp_path / "enrichment.tsv", ora_rows(), metadata)
    with pytest.raises(ValueError, match=message):
        ora_profile_records(enrichment, tmp_path / "missing.tsv", columns=COLUMNS)


def test_ora_rejects_call_flag_that_disagrees_with_threshold(tmp_path: Path) -> None:
    rows = ora_rows()
    rows[0]["is_significant"] = "False"
    enrichment = write_source(tmp_path / "enrichment.tsv", rows, ora_metadata())
    with pytest.raises(ValueError, match="hit flag"):
        ora_profile_records(enrichment, tmp_path / "missing.tsv", columns=COLUMNS)


def test_ora_rejects_duplicate_profile_trait_test(tmp_path: Path) -> None:
    rows = ora_rows()
    rows[1] = rows[0].copy()
    enrichment = write_source(tmp_path / "enrichment.tsv", rows, ora_metadata())
    with pytest.raises(ValueError, match="duplicate ORA"):
        ora_profile_records(enrichment, tmp_path / "missing.tsv", columns=COLUMNS)


def test_ora_rejects_truncated_all_results(tmp_path: Path) -> None:
    enrichment = write_source(tmp_path / "enrichment.tsv", ora_rows()[:-1], ora_metadata())
    with pytest.raises(ValueError, match="correction family"):
        ora_profile_records(enrichment, tmp_path / "missing.tsv", columns=COLUMNS)


@pytest.mark.parametrize(
    ("field", "value"),
    [("background_m", "101"), ("background_k", "101"), ("fold", "-1")],
)
def test_ora_rejects_invalid_background_or_fold(tmp_path: Path, field: str, value: str) -> None:
    rows = ora_rows()
    rows[0][field] = value
    enrichment = write_source(tmp_path / "enrichment.tsv", rows, ora_metadata())
    with pytest.raises(ValueError, match=r"background|fold"):
        ora_profile_records(enrichment, tmp_path / "missing.tsv", columns=COLUMNS)


def test_ora_same_counts_but_different_profiles_do_not_join(tmp_path: Path) -> None:
    rows = ora_rows()
    rows[0]["study_id"] = "another-study"
    enrichment = write_source(tmp_path / "enrichment.tsv", rows, ora_metadata())
    with pytest.raises(ValueError, match="profile"):
        ora_profile_records(enrichment, tmp_path / "missing.tsv", columns=COLUMNS)


def test_ora_depletion_profile_identity_must_match_enrichment(tmp_path: Path) -> None:
    enrichment = write_source(tmp_path / "enrichment.tsv", ora_rows(), ora_metadata())
    depleted = ora_rows(enriched=False)
    for row in depleted:
        row["study_id"] = "another-study"
    depletion = write_source(tmp_path / "depletion.tsv", depleted, ora_metadata("lower"))
    with pytest.raises(ValueError, match="profile"):
        ora_profile_records(enrichment, depletion, columns=COLUMNS)


def test_ora_rejects_profiles_outside_current_integrated_population(tmp_path: Path) -> None:
    enrichment = write_source(tmp_path / "enrichment.tsv", ora_rows(), ora_metadata())
    expected = [
        Profile("another-study", "experiment1", direction, "", "", {}, 1)
        for direction in ("inc", "dec")
    ]
    with pytest.raises(ValueError, match="profile"):
        ora_profile_records(
            enrichment, tmp_path / "missing.tsv", columns=COLUMNS, expected_profiles=expected
        )


def test_ora_family_duplicate_is_rejected_even_outside_table1(tmp_path: Path) -> None:
    rows = ora_rows()
    extra = dict(rows[0], trait="METPO:1000603")
    rows.extend([extra, extra.copy()])
    metadata = ora_metadata()
    metadata["Total-Tests-Performed"] = "14"
    enrichment = write_source(tmp_path / "enrichment.tsv", rows, metadata)
    with pytest.raises(ValueError, match="duplicate"):
        ora_profile_records(enrichment, tmp_path / "missing.tsv", columns=COLUMNS)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("correction_family", "profile"),
        ("background_m", "101"),
        ("fdr_q", "nan"),
        ("fold", "-1"),
        ("is_significant", "unknown"),
        ("eligibility_status", "unknown"),
    ],
)
def test_ora_validates_family_rows_outside_selected_traits(
    tmp_path: Path, field: str, value: str
) -> None:
    rows = ora_rows()
    extras = [dict(row, trait="METPO:1000603") for row in rows[:2]]
    extras[1][field] = value
    metadata = ora_metadata()
    metadata["Total-Tests-Performed"] = "14"
    source = write_source(tmp_path / "all.tsv", rows + extras, metadata)
    with pytest.raises(ValueError, match=r"ORA|non-finite"):
        ora_profile_records(source, tmp_path / "missing.tsv", columns=COLUMNS)


def test_ora_rejects_partial_unselected_trait_even_if_total_is_changed(tmp_path: Path) -> None:
    rows = [*ora_rows(), dict(ora_rows()[0], trait="METPO:1000603")]
    metadata = ora_metadata()
    metadata["Total-Tests-Performed"] = "13"
    source = write_source(tmp_path / "all.tsv", rows, metadata)
    with pytest.raises(ValueError, match="per-trait profile populations"):
        ora_profile_records(source, tmp_path / "missing.tsv", columns=COLUMNS)


def test_ora_rejects_background_prevalence_that_changes_by_profile(tmp_path: Path) -> None:
    rows = ora_rows()
    rows[1]["background_k"] = "11"
    source = write_source(tmp_path / "all.tsv", rows, ora_metadata())
    with pytest.raises(ValueError, match="prevalence changes"):
        ora_profile_records(source, tmp_path / "missing.tsv", columns=COLUMNS)


def fraction_records() -> list[Record]:
    return [
        replace(
            record,
            test_status=TESTED,
            unique_studies=11,
            testable_units=10,
            ties=1,
            effect_value=0.8,
            p_value=0.109375,
            q_value=0.109375,
        )
        for record in placeholder_records(
            "directional_consistency_fraction_sign_test",
            "study",
            "fixture",
            TABLE1_TRAITS,
            columns=COLUMNS,
        )
    ]


def screen_sources(tmp_path: Path) -> tuple[Path, Path]:
    rows = [
        {
            "trait": trait,
            "studies_scored": "10",
            "studies_tied": "1",
            "family_size_nominal": "7",
            "family_id": "trait-direction-screen-v1",
            "test_status": TESTED,
            "binomial_p": "0.109375",
            "share_up": "0.8",
            "bh_q_nominal": "0.2",
            "by_q_nominal": "0.5",
        }
        for trait in [*RAW_COLUMNS.values(), "METPO:1000603"]
    ]
    path = write_source(tmp_path / "screen.tsv", rows)
    provenance = tmp_path / "screen.json"
    provenance.write_text(
        json.dumps(
            {
                "profile_sha256": PROFILE_HASH,
                "experiments_sha256": EXPERIMENTS_HASH,
                "family_size_nominal": 7,
            }
        ),
        encoding="utf-8",
    )
    return path, provenance


def test_full_screen_attaches_both_corrections_without_replacing_six_trait_q(
    tmp_path: Path,
) -> None:
    records = attach_full_screen(
        fraction_records(),
        *screen_sources(tmp_path),
        profile_sha256=PROFILE_HASH,
        experiments_sha256=EXPERIMENTS_HASH,
    )
    assert {r.trait: r.trait_columns for r in records} == RAW_COLUMNS
    assert all(r.q_value == 0.109375 for r in records)
    assert all(r.full_screen_bh_q == 0.2 and r.full_screen_by_q == 0.5 for r in records)
    assert all(r.full_screen_family_size == 7 for r in records)
    assert all(r.full_screen_source.startswith("screen.tsv@") for r in records)


@pytest.mark.parametrize("wrong_hash", ["profile", "experiments"])
def test_full_screen_rejects_incompatible_input_hash(tmp_path: Path, wrong_hash: str) -> None:
    with pytest.raises(ValueError, match="checksum mismatch"):
        attach_full_screen(
            fraction_records(),
            *screen_sources(tmp_path),
            profile_sha256="different" if wrong_hash == "profile" else PROFILE_HASH,
            experiments_sha256="different" if wrong_hash == "experiments" else EXPERIMENTS_HASH,
        )


@pytest.mark.parametrize("mismatch", ["population", "contrast"])
def test_full_screen_rejects_different_population_or_effect(tmp_path: Path, mismatch: str) -> None:
    records = fraction_records()
    records[0] = (
        replace(records[0], testable_units=9)
        if mismatch == "population"
        else replace(records[0], effect_value=0.2)
    )
    with pytest.raises(ValueError, match=r"population|contrast"):
        attach_full_screen(
            records,
            *screen_sources(tmp_path),
            profile_sha256=PROFILE_HASH,
            experiments_sha256=EXPERIMENTS_HASH,
        )


def test_external_family_subset_requires_source_provenance(tmp_path: Path) -> None:
    records = import_composition(write_source(tmp_path / "native.tsv", composition_rows()))
    with pytest.raises(ValueError, match="complete source"):
        family_inventory([replace(record, family_source="") for record in records])
