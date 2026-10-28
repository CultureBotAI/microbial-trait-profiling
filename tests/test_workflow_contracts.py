"""Scientific product contracts retain rank overlap, missingness, and count units."""

import csv
from collections import Counter
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path

import pytest

from mtp.count_input import annotation_coverage_status, read_profile_coverage
from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES
from mtp.study_trait_profiles import TaxaTraitMatrix, read_taxa_trait_matrix
from mtp.workflow.contracts import read_table, validate_product
from mtp.workflow.coverage import (
    FEATURE_METADATA_COLUMNS,
    build_profile_coverage,
    feature_filter_id,
    write_feature_catalog,
    write_profile_coverage,
)
from mtp.workflow.operations import _CORE_SEMANTICS

TRAIT = "METPO:1000001"
MEDIUM = "METPO:2000001|mediadive.medium:1"
KEY = ("bsdb:1", "bsdb:1-1", "inc")


def _write(path: Path, header: Sequence[str], rows: Sequence[Sequence[str | int]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def _matrix(tmp_path: Path) -> Path:
    return _write(
        tmp_path / "matrix.tsv",
        ["taxon_id", "taxon_name", "taxonomic_rank", TRAIT, MEDIUM],
        [
            ["NCBITaxon:10", "Genus", "genus", 3, 0],
            ["NCBITaxon:1", "Species", "species", 1, 0],
            ["NCBITaxon:20", "Medium only", "family", 0, 1],
            ["NCBITaxon:30", "Unannotated", "species", 0, 0],
        ],
    )


def _sources(tmp_path: Path) -> dict[str, Path]:
    matrix = _matrix(tmp_path)
    observations = [
        "NCBITaxon:10",
        "NCBITaxon:10",
        "NCBITaxon:1",
        "NCBITaxon:20",
        "NCBITaxon:30",
        "NCBITaxon:99",
    ]
    taxa = _write(
        tmp_path / "taxa.tsv",
        [
            "study_id",
            "experiment_id",
            "signature_id",
            "direction",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
        ],
        [
            ["bsdb:1", "bsdb:1-1", "bsdb:1-1-1", "increased", taxon, "Named taxon", "species"]
            for taxon in observations
        ]
        + [
            [
                "bsdb:1",
                "bsdb:1-2",
                "bsdb:1-2-1",
                "increased",
                "NCBITaxon:1",
                "Mouse species",
                "species",
            ]
        ],
    )
    experiments = _write(
        tmp_path / "experiments.tsv",
        [
            "study_id",
            "experiment_id",
            "host_species",
            "sequencing_type",
            "signature_count",
            "taxon_count",
        ],
        [
            ["bsdb:1", "bsdb:1-1", "Homo sapiens", "16S", 2, 6],
            ["bsdb:1", "bsdb:1-2", "Mus musculus", "WMS", 1, 1],
        ],
    )
    signatures = _write(
        tmp_path / "signatures.tsv",
        ["study_id", "experiment_id", "signature_id", "direction", "taxon_count"],
        [
            ["bsdb:1", "bsdb:1-1", "bsdb:1-1-1", "increased", 6],
            ["bsdb:1", "bsdb:1-1", "bsdb:1-1-2", "decreased", 0],
            ["bsdb:1", "bsdb:1-2", "bsdb:1-2-1", "increased", 1],
        ],
    )
    return {
        "taxa_matrix": matrix,
        "reported_taxa": taxa,
        "experiments": experiments,
        "signatures": signatures,
    }


def test_count_and_ora_units_with_mixed_ranks_duplicates_and_unmatched(tmp_path: Path) -> None:
    context = _sources(tmp_path)
    counts = _write(
        tmp_path / "counts.tsv",
        ["study_id", "experiment_id", "inc_dec", TRAIT, MEDIUM],
        [
            [*KEY, 4, 1],
            ["bsdb:1", "bsdb:1-1", "dec", 0, 0],
            ["bsdb:1", "bsdb:1-2", "inc", 1, 0],
        ],
    )
    ora = _write(
        tmp_path / "ora.tsv",
        ["study_id", "experiment_id", "inc_dec", "taxa_count", TRAIT, MEDIUM],
        [
            [*KEY, 4, 2, 1],
            ["bsdb:1", "bsdb:1-1", "dec", 0, 0, 0],
        ],
    )
    count_report = validate_product(counts, "count_profiles", context=context)
    ora_report = validate_product(ora, "ora_profiles", context=context)
    assert count_report["rows"] == 3
    assert count_report["source_aggregation_checked"] is True
    assert ora_report["rows"] == 2
    assert ora_report["cohort"] == "Homo sapiens"
    # A genus plus a reported descendant remains two reported IDs, although
    # their annotation source lineages overlap. Duplicate genus observations
    # do not multiply either profile unit.
    assert count_report["unit"] == "summed_annotation_counts"
    assert (
        ora_report["denominator"] == "distinct_lookup_matched_reported_taxa_including_unannotated"
    )


@pytest.mark.parametrize(
    ("n", "trait", "message"),
    [
        (2, 2, "denominator is not distinct"),
        (4, 4, "values do not match source aggregation"),
        (1, 2, "exceeds matched-taxon denominator"),
    ],
)
def test_ora_rejects_annotation_counts_or_annotated_only_denominator(
    tmp_path: Path, n: int, trait: int, message: str
) -> None:
    context = _sources(tmp_path)
    path = _write(
        tmp_path / "ora.tsv",
        ["study_id", "experiment_id", "inc_dec", "taxa_count", TRAIT, MEDIUM],
        [[*KEY, n, trait, 1]],
    )
    with pytest.raises(ValueError, match=message):
        validate_product(path, "ora_profiles", context=context)


def test_profile_cohort_membership_is_checked_not_inferred_from_rows(tmp_path: Path) -> None:
    context = _sources(tmp_path)
    path = _write(
        tmp_path / "ora.tsv",
        ["study_id", "experiment_id", "inc_dec", "taxa_count", TRAIT, MEDIUM],
        [[*KEY, 4, 2, 1]],
    )
    with pytest.raises(ValueError, match="membership differs"):
        validate_product(path, "ora_profiles", context=context)


def test_profile_cohort_keeps_16s_and_shotgun_experiments_only(tmp_path: Path) -> None:
    """A count profile of a PCR experiment is outside the cohort, whatever the host (#742)."""
    context = _sources(tmp_path)
    experiments = context["experiments"]
    text = experiments.read_text(encoding="utf-8")
    experiments.write_text(text.replace("\tWMS\t", "\tPCR\t"), encoding="utf-8")
    rows: list[list[str | int]] = [[*KEY, 4, 1], ["bsdb:1", "bsdb:1-1", "dec", 0, 0]]
    header = ["study_id", "experiment_id", "inc_dec", TRAIT, MEDIUM]
    with_pcr = _write(
        tmp_path / "with_pcr.tsv", header, [*rows, ["bsdb:1", "bsdb:1-2", "inc", 1, 0]]
    )
    with pytest.raises(ValueError, match="outside declared cohort"):
        validate_product(with_pcr, "count_profiles", context=context)
    without_pcr = _write(tmp_path / "without_pcr.tsv", header, rows)
    assert validate_product(without_pcr, "count_profiles", context=context)["rows"] == 2


def test_experiments_from_before_the_sequencing_rule_are_named(tmp_path: Path) -> None:
    """An experiments table without sequencing_type says the core must be rebuilt (#762)."""
    context = _sources(tmp_path)
    experiments = context["experiments"]
    lines = experiments.read_text(encoding="utf-8").splitlines()
    rows = [line.split("\t") for line in lines]
    column = rows[0].index("sequencing_type")
    experiments.write_text(
        "".join("\t".join(r[:column] + r[column + 1 :]) + "\n" for r in rows), encoding="utf-8"
    )
    counts = _write(
        tmp_path / "counts.tsv",
        ["study_id", "experiment_id", "inc_dec", TRAIT, MEDIUM],
        [[*KEY, 4, 1]],
    )
    with pytest.raises(ValueError, match="before #742"):
        validate_product(counts, "count_profiles", context=context)


def test_policy_label_names_the_corpus_sequencing_types() -> None:
    """Editing the constant must change the declared policy label too (#764)."""
    label = "count_input_" + "_".join(t.lower() for t in CORPUS_SEQUENCING_TYPES) + "_v1"
    for kind in ("count_profiles", "ora_profiles", "profile_coverage"):
        assert _CORE_SEMANTICS[kind].coverage_policy == label
    for spec in ("bugsigdb_core.yaml", "bugsigdb_core_subset.yaml"):
        text = (Path(__file__).parents[1] / "workflows" / spec).read_text(encoding="utf-8")
        assert text.count(f"coverage_policy: {label}") == 4


def test_coverage_keeps_lookup_zero_and_media_only_distinct(tmp_path: Path) -> None:
    matrix = read_taxa_trait_matrix(_matrix(tmp_path))
    taxa = {KEY: ["NCBITaxon:10", "NCBITaxon:10", "NCBITaxon:20", "NCBITaxon:30", "NCBITaxon:99"]}
    all_features = build_profile_coverage(taxa, matrix)
    record = all_features[KEY]
    assert (
        record.reported_distinct_taxa,
        record.lookup_matched_taxa,
        record.lookup_unmatched_taxa,
    ) == (4, 3, 1)
    assert (record.matched_taxa_with_annotations, record.matched_taxa_without_annotations) == (2, 1)
    assert record.unmatched_taxon_ids == ("NCBITaxon:99",)
    assert record.matched_without_annotation_taxon_ids == ("NCBITaxon:30",)
    retained = build_profile_coverage(taxa, matrix, feature_columns=[TRAIT])[KEY]
    assert retained.matched_without_annotation_taxon_ids == ("NCBITaxon:20", "NCBITaxon:30")
    assert retained.trait_filter_id != record.trait_filter_id
    path = tmp_path / "coverage.tsv"
    write_profile_coverage(path, all_features)
    assert read_profile_coverage(path) == all_features
    assert validate_product(path, "profile_coverage")["rows"] == 1


@pytest.mark.parametrize(
    ("reported", "matched", "unmatched", "annotated", "unannotated", "expected"),
    [
        (0, 0, 0, 0, 0, "no_reported_taxa"),
        (1, 0, 1, 0, 0, "no_lookup_matches"),
        (2, 1, 1, 1, 0, "partial_lookup_coverage"),
        (1, 1, 0, 0, 1, "matched_without_annotations"),
        (2, 2, 0, 1, 1, "partial_annotation_coverage"),
        (1, 1, 0, 1, 0, "complete_annotation_coverage"),
    ],
)
def test_coverage_states(
    reported: int, matched: int, unmatched: int, annotated: int, unannotated: int, expected: str
) -> None:
    assert (
        annotation_coverage_status(
            reported=reported,
            matched=matched,
            unmatched=unmatched,
            with_annotations=annotated,
            without_annotations=unannotated,
        )
        == expected
    )


def test_coverage_rejects_forged_status_and_overlapping_id_sets(tmp_path: Path) -> None:
    matrix = TaxaTraitMatrix([TRAIT], {"NCBITaxon:1": Counter()})
    coverage = build_profile_coverage({KEY: ["NCBITaxon:1", "NCBITaxon:2"]}, matrix)
    path = tmp_path / "coverage.tsv"
    write_profile_coverage(
        path, {KEY: replace(coverage[KEY], coverage_status="complete_annotation_coverage")}
    )
    with pytest.raises(ValueError, match="status mismatch"):
        read_profile_coverage(path)
    write_profile_coverage(
        path, {KEY: replace(coverage[KEY], unmatched_taxon_ids=("NCBITaxon:1",))}
    )
    with pytest.raises(ValueError, match="overlap"):
        read_profile_coverage(path)


def test_feature_catalog_preserves_qualified_ids_media_and_unresolved_origin(
    tmp_path: Path,
) -> None:
    metadata = _write(
        tmp_path / "metadata.tsv",
        list(FEATURE_METADATA_COLUMNS),
        [
            [TRAIT, "Trait", "", "", TRAIT, "Trait", "biolink:PhenotypicQuality"],
            [
                MEDIUM,
                "Medium",
                "METPO:2000001",
                "grows in",
                "mediadive.medium:1",
                "Medium 1",
                "biolink:ChemicalMixture",
            ],
        ],
    )
    catalog = tmp_path / "catalog.tsv"
    write_feature_catalog(metadata, catalog, include_policy_metadata=False)
    report = validate_product(
        catalog,
        "feature_catalog",
        context={"taxa_matrix": _matrix(tmp_path)},
        include_policy_metadata=False,
    )
    assert report["culture_medium_features"] == 1
    assert report["annotation_evidence"] == "direct_or_descendant_unresolved"
    assert report["feature_filter_id"] == feature_filter_id([TRAIT, MEDIUM])
    catalog.write_text(catalog.read_text().replace("culture_medium", "trait_annotation"))
    with pytest.raises(ValueError, match="MediaDive classification"):
        validate_product(catalog, "feature_catalog", include_policy_metadata=False)


def test_table_foreign_keys_and_duplicate_observations(tmp_path: Path) -> None:
    context = _sources(tmp_path)
    facts = validate_product(context["reported_taxa"], "reported_taxa", context=context)
    assert facts["rows"] == 7
    assert facts["duplicate_observations_retained"] is True
    validate_product(context["signatures"], "signature_table", context=context)
    path = context["reported_taxa"]
    path.write_text(path.read_text().replace("bsdb:1-1-1", "bsdb:1-1-unknown", 1))
    with pytest.raises(ValueError, match="missing parent"):
        validate_product(path, "reported_taxa", context=context)


def test_matrix_cross_product_vocabulary_and_zero_rows(tmp_path: Path) -> None:
    context = _sources(tmp_path)
    report = validate_product(context["taxa_matrix"], "taxa_matrix")
    assert report["all_zero_rows"] == 1
    with pytest.raises(ValueError, match="taxon vocabulary"):
        validate_product(context["taxa_matrix"], "taxa_matrix", context=context)


def test_long_taxon_profiles_are_checked_as_transposes(tmp_path: Path) -> None:
    wide = _write(
        tmp_path / "wide.tsv",
        ["study_id", "experiment_id", "inc_dec", "Taxon A", "Taxon B"],
        [[*KEY, 2, 0]],
    )
    long = _write(tmp_path / "long.tsv", ["taxon", "_".join(KEY)], [["Taxon A", 2], ["Taxon B", 0]])
    assert (
        validate_product(long, "taxon_profile_long", context={"taxon_profiles_wide": wide})[
            "layout"
        ]
        == "transpose_of_wide_not_tidy_long"
    )
    long.write_text(long.read_text().replace("Taxon B\t0", "Taxon B\t1"))
    with pytest.raises(ValueError, match="values differ"):
        validate_product(long, "taxon_profile_long", context={"taxon_profiles_wide": wide})


def test_rectangular_table_and_unique_identifiers(tmp_path: Path) -> None:
    path = tmp_path / "table.tsv"
    path.write_text("id\tid\na\tb\n")
    with pytest.raises(ValueError, match="duplicate column"):
        read_table(path)
    path.write_text("id\tvalue\na\n")
    with pytest.raises(ValueError, match="wrong row width"):
        read_table(path)


def test_taxon_profiles_count_observations_before_distinct_id_trait_aggregation(
    tmp_path: Path,
) -> None:
    context = _sources(tmp_path)
    path = _write(
        tmp_path / "taxon_wide.tsv",
        ["study_id", "experiment_id", "inc_dec", "Mouse species", "Named taxon"],
        [[*KEY, 0, 6], ["bsdb:1", "bsdb:1-2", "inc", 1, 0]],
    )
    assert (
        validate_product(path, "taxon_profile_wide", context=context)["source_aggregation_checked"]
        is True
    )
    path.write_text(path.read_text().replace("inc\t0\t6", "inc\t0\t5"))
    with pytest.raises(ValueError, match="counts differ from reported observations"):
        validate_product(path, "taxon_profile_wide", context=context)


def test_coverage_schema_rejects_non_digest_filter_identity(tmp_path: Path) -> None:
    matrix = read_taxa_trait_matrix(_matrix(tmp_path))
    coverage = build_profile_coverage({KEY: ["NCBITaxon:1"]}, matrix)
    path = tmp_path / "coverage.tsv"
    write_profile_coverage(path, {KEY: replace(coverage[KEY], trait_filter_id="unknown")})
    with pytest.raises(ValueError, match="trait_filter_id"):
        validate_product(path, "profile_coverage")
