"""Consumer annotation coverage follows its retained vocabulary, not the core total."""

from pathlib import Path
from unittest.mock import Mock

import pytest
import yaml

from mtp.count_input import read_profile_coverage
from mtp.workflow.consumers import (
    CoreBundle,
    build_consumer_coverage,
    validate_consumer_coverage,
)


@pytest.fixture
def bundle(tmp_path: Path) -> CoreBundle:
    studies = tmp_path / "studies"
    studies.mkdir()
    study = {
        "id": "bsdb:1",
        "experiments": [
            {
                "id": "bsdb:1-1",
                "host_species": "Homo sapiens",
                "sequencing_type": "16S",
                "signatures": [
                    {
                        "id": "bsdb:1-1-1",
                        "direction": "increased",
                        "taxa": [
                            {"id": "NCBITaxon:1"},
                            {"id": "NCBITaxon:1"},
                            {"id": "NCBITaxon:9"},
                        ],
                    },
                    {"id": "bsdb:1-1-2", "direction": "decreased", "taxa": [{"id": "NCBITaxon:2"}]},
                ],
            },
            {
                "id": "bsdb:1-2",
                "host_species": "Mus musculus",
                "sequencing_type": "WMS",
                "signatures": [
                    {"id": "bsdb:1-2-1", "direction": "increased", "taxa": [{"id": "NCBITaxon:3"}]},
                ],
            },
        ],
    }
    (studies / "bsdb_1.yaml").write_text(yaml.safe_dump(study), encoding="utf-8")
    matrix = tmp_path / "matrix.tsv"
    matrix.write_text(
        "taxon_id\ttaxon_name\ttaxonomic_rank\tMETPO:1\tMETPO:2|MediaDive.Medium:1\n"
        "NCBITaxon:1\tMedium only\tgenus\t0\t2\n"
        "NCBITaxon:2\tZero annotations\tspecies\t0\t0\n"
        "NCBITaxon:3\tTrait annotated\tspecies\t1\t0\n",
        encoding="utf-8",
    )
    result = Mock(spec=CoreBundle)
    result.cohort = "all"
    result.path.side_effect = {"studies": studies, "taxa_matrix": matrix}.__getitem__
    return result


def test_filtered_coverage_reclassifies_medium_only_without_losing_lookup(
    bundle: CoreBundle, tmp_path: Path
) -> None:
    all_path = build_consumer_coverage(bundle, tmp_path / "all.tsv", "include_mediadive")
    filtered_path = build_consumer_coverage(bundle, tmp_path / "filtered.tsv")
    key = ("bsdb:1", "bsdb:1-1", "inc")
    original = read_profile_coverage(all_path)[key]
    retained = read_profile_coverage(filtered_path)[key]
    assert original.matched_taxa_with_annotations == 1
    assert retained.matched_taxa_with_annotations == 0
    assert retained.matched_without_annotation_taxon_ids == ("NCBITaxon:1",)
    assert original.unmatched_taxon_ids == retained.unmatched_taxon_ids == ("NCBITaxon:9",)
    assert original.lookup_matched_taxa == retained.lookup_matched_taxa == 1
    assert retained.reported_distinct_taxa == 2
    assert retained.trait_filter_id != original.trait_filter_id
    facts = validate_consumer_coverage(bundle, filtered_path)
    assert facts["features"] == 1
    assert facts["profiles"] == 3  # Input cohort, including unpaired/nonhuman profiles.
    assert facts["coverage_scope"] == "input_cohort_before_analysis_specific_exclusions"


def test_filtered_coverage_rejects_core_all_feature_substitution(
    bundle: CoreBundle, tmp_path: Path
) -> None:
    path = build_consumer_coverage(bundle, tmp_path / "all.tsv", "include_mediadive")
    with pytest.raises(ValueError, match="retained features"):
        validate_consumer_coverage(bundle, path, "exclude_mediadive")


def test_filtered_coverage_respects_effective_input_cohort(
    bundle: CoreBundle, tmp_path: Path
) -> None:
    human = Mock(spec=CoreBundle)
    human.path = bundle.path
    human.cohort = "Homo sapiens"
    path = build_consumer_coverage(human, tmp_path / "human.tsv")
    assert validate_consumer_coverage(human, path)["profiles"] == 2
    with pytest.raises(ValueError, match="input cohort"):
        validate_consumer_coverage(bundle, path)


def test_consumer_coverage_rejects_unknown_feature_policy(
    bundle: CoreBundle, tmp_path: Path
) -> None:
    with pytest.raises(ValueError, match="Unknown consumer feature policy"):
        build_consumer_coverage(bundle, tmp_path / "coverage.tsv", "silently_drop_other_traits")
