"""Issue #1044: complete ORA calls, paired study votes and full-family correction."""

from __future__ import annotations

import csv
import io
from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from scipy.stats import hypergeom

from mtp.datamodel.workflow_pydantic import (
    AnalysisBasis,
    ArmCommandExecution,
    ArmRun,
    ArtifactInstance,
    ComparisonRun,
    ContrastPopulation,
    DirectionBasis,
    WorkflowRun,
)
from mtp.enrichment_stats import benjamini_hochberg
from mtp.multiplicity.families import rows
from mtp.multiplicity.ora import ora_comments
from mtp.multiplicity.sources import VerifiedSources, recorded_arm_output
from mtp.paths import file_sha256
from mtp.three_test_table import ora
from mtp.three_test_table.ora import CALIBRATION_STATUS, ora_table
from mtp.workflow.coverage import feature_filter_id
from multiplicity_fixture import file_snapshot, write_table

FEATURES = ("trait-a", "trait-b", "trait-c", "trait-empty")


def ora_sources(
    root: Path,
    patterns: dict[str, list[tuple[int | None, int | None]]],
    *,
    depleted: frozenset[tuple[str, str, str]] = frozenset(),
    unpaired: bool = False,
    draw_sizes: dict[tuple[str, str, str], int] | None = None,
) -> VerifiedSources:
    """Construct genuine hypergeometric exports; stub execution records, not statistics."""
    root.mkdir()
    products = {
        name: root / f"{name}.tsv"
        for name in ("ora_profiles", "taxa_matrix", "experiments", "feature_catalog")
    }
    profiles: list[dict[str, Any]] = []
    experiments = []
    for study, pairs in patterns.items():
        for number, pair in enumerate(pairs):
            experiment = f"e{number}"
            experiments.append({"study_id": study, "experiment_id": experiment})
            for direction, hit in zip(("inc", "dec"), pair, strict=True):
                profiles.append(
                    {
                        "study_id": study,
                        "experiment_id": experiment,
                        "inc_dec": direction,
                        "taxa_count": 4 if hit is not None else 0,
                        "trait-a": 4 if hit else 0,
                        "trait-b": 4 if hit is not None else 0,
                        "trait-c": 4 if hit == 0 else 0,
                        "trait-empty": 0,
                    }
                )
                if (study, experiment, direction) in depleted:
                    # This list draws the four all-zero matrix taxa, retaining N=4.
                    profiles[-1].update(dict.fromkeys(FEATURES, 0))
                if draw_sizes is not None and (study, experiment, direction) in draw_sizes:
                    profiles[-1]["taxa_count"] = draw_sizes[(study, experiment, direction)]
    if unpaired:
        profiles.append({**profiles[0], "study_id": "unpaired-study"})
    write_table(products["ora_profiles"], profiles)
    write_table(products["experiments"], experiments)
    write_table(
        products["feature_catalog"],
        [{"column_id": feature, "feature_class": "trait_annotation"} for feature in FEATURES],
    )
    write_table(
        products["taxa_matrix"],
        [
            {
                "taxon_id": f"NCBITaxon:{index}",
                "taxon_name": f"Fixture {index}",
                "taxonomic_rank": "species",
                "trait-a": int(index < 4),
                "trait-b": int(index < 96),
                "trait-c": int(4 <= index < 8),
                "trait-empty": 0,
            }
            for index in range(100)
        ],
    )
    paired_keys = {(row["study_id"], row["experiment_id"]) for row in experiments}
    eligible = [
        row
        for row in profiles
        if row["taxa_count"] and (row["study_id"], row["experiment_id"]) in paired_keys
    ]
    background = {"trait-a": 4, "trait-b": 96, "trait-c": 4}
    steps = []
    for tail, step in (("upper", "enrichment"), ("lower", "depletion")):
        path = root / f"trait_{step}_all_results.tsv"
        records: list[dict[str, Any]] = []
        for profile in eligible:
            for feature, carriers in background.items():
                observed = profile[feature]
                draws = profile["taxa_count"]
                p = float(
                    hypergeom.sf(observed - 1, 100, carriers, draws)
                    if tail == "upper"
                    else hypergeom.cdf(observed, 100, carriers, draws)
                )
                expected = draws * carriers / 100
                records.append(
                    {
                        **{key: profile[key] for key in ("study_id", "experiment_id", "inc_dec")},
                        "trait": feature,
                        "p_value": p,
                        "fdr_q": 1.0,
                        "fold": observed / expected,
                        "count": observed,
                        "expected": expected,
                        "testable_profiles": len(eligible),
                        "background_m": 100,
                        "background_k": carriers,
                        "is_significant": False,
                        "correction_family": "global",
                        "eligibility_status": "testable",
                    }
                )
        q_values = benjamini_hochberg([record["p_value"] for record in records])
        for record, q in zip(records, q_values, strict=True):
            record["fdr_q"] = q
            record["is_significant"] = q <= 0.001 and (
                record["fold"] >= 1.5 if tail == "upper" else record["fold"] <= 0.5
            )
        comments = {
            "Profiles-SHA256": file_sha256(products["ora_profiles"]),
            "Taxa-Matrix-SHA256": file_sha256(products["taxa_matrix"]),
            "Tail": tail,
            "Testable-Set": "K>0 & N>0",
            "FDR": "True",
            "FDR-Scope": "global",
            "P-Threshold": "0.001",
            "Testable-Profiles": str(len(eligible)),
            "Background-Taxa-Universe-M": "100",
            "Background-Mode": f"TaxaMatrix ({products['taxa_matrix'].name})",
            "Total-Tests-Performed": str(len(records)),
            "Min-Fold" if tail == "upper" else "Max-Fold": "1.5" if tail == "upper" else "0.5",
        }
        if records:
            write_table(path, records, comments)
        else:
            with path.open("w", encoding="utf-8", newline="") as handle:
                for key, value in comments.items():
                    handle.write(f"# {key}: {value}\n")
                csv.writer(handle, delimiter="\t").writerow(
                    [
                        "study_id",
                        "experiment_id",
                        "inc_dec",
                        "trait",
                        "p_value",
                        "fdr_q",
                        "fold",
                        "count",
                        "expected",
                        "testable_profiles",
                        "background_m",
                        "background_k",
                        "is_significant",
                        "correction_family",
                        "eligibility_status",
                    ]
                )
        steps.append(ArmCommandExecution.model_construct(id=step, outputs=[file_snapshot(path)]))
    manifests = {name: root / f"{name}.json" for name in ("core", "comparison", "ilyas")}
    for path in manifests.values():
        path.write_text('{"synthetic": true}\n', encoding="utf-8")
    core = WorkflowRun.model_construct(
        run_id="synthetic-ora",
        manifest_path=str(manifests["core"]),
        artifacts=[
            ArtifactInstance.model_construct(product_id=name, path=str(path))
            for name, path in products.items()
        ],
    )
    basis = AnalysisBasis(
        id="synthetic-ora",
        feature_catalog_sha256=file_sha256(products["feature_catalog"]),
        feature_filter_id=feature_filter_id(FEATURES),
        retained_feature_count=len(FEATURES),
        profile_unit="ora_profiles",
        aggregation_policy="binary distinct taxa",
        rollup_policy="ncbi_taxonomy_closure_rollup_off",
        denominator_policy="matched taxa",
        background_inventory_id="fixture",
        weighting_policy="equal studies",
    )
    comparison = ComparisonRun.model_construct(
        manifest_path=str(manifests["comparison"]),
        contrast_population=ContrastPopulation.all_paired,
        direction_basis=DirectionBasis.group_1_vs_group_0,
    )
    arm = ArmRun.model_construct(manifest_path=str(manifests["ilyas"]), steps=steps)
    paths = [
        *products.values(),
        *manifests.values(),
        *(root / f"trait_{name}_all_results.tsv" for name in ("enrichment", "depletion")),
    ]
    return VerifiedSources(
        comparison,
        core,
        {"ilyas": arm},
        {"ora_profiles": basis},
        {"ora_profiles": FEATURES},
        "off",
        {path.resolve(): file_sha256(path) for path in paths},
    )


def test_ora_study_votes_use_common_pairs_and_equal_study_weights(tmp_path: Path) -> None:
    sources = ora_sources(
        tmp_path / "source",
        {
            "A": [(1, 0), (0, 0)],
            "B": [(1, 0), (0, 1)],
            "C": [(0, 1)],
            "D": [(1, None)],
            "E": [(None, None)],
        },
    )
    table, inputs = ora_table(sources)
    row = table.set_index("trait").loc["trait-a"]
    assert (row.ora_enriched_inc, row.ora_enriched_dec) == (3, 2)
    assert (
        row.ora_direction_studies_up,
        row.ora_direction_studies_down,
        row.ora_direction_studies_tied,
    ) == (1, 1, 1)
    assert row.ora_direction_mean_study_delta == pytest.approx(-1 / 6)
    assert row.ora_direction_pairs_eligible == 5
    assert row.ora_direction_pairs_excluded == 2
    assert row.ora_direction_studies_eligible == 3
    assert row.ora_direction_studies_excluded == 2
    assert row.ora_direction_share_up == 0.5
    assert row.ora_direction_p == 1.0
    assert row.ora_background_m == 100
    assert row.ora_background_k == 4
    assert set(table.ora_calibration_status) == {CALIBRATION_STATUS}
    assert set(inputs) >= {sources.product("ora_profiles"), sources.product("taxa_matrix")}


def test_ora_full_family_holds_untestable_slots_and_reports_exact_interval(tmp_path: Path) -> None:
    sources = ora_sources(tmp_path / "source", {f"S{i}": [(1, 0)] for i in range(9)})
    table, _ = ora_table(sources)
    indexed = table.set_index("trait")
    row = indexed.loc["trait-a"]
    assert row.ora_direction_p == pytest.approx(1 / 256)
    assert row.ora_direction_bh_q == pytest.approx(1 / 128)
    assert row.ora_direction_by_q == pytest.approx(25 / 1536)
    assert row.ora_direction_share_up_ci95_low == pytest.approx(0.025 ** (1 / 9))
    assert row.ora_direction_share_up_ci95_high == 1.0
    assert row.ora_direction_family_size == 4
    assert indexed.loc["trait-b", "ora_direction_status"] == "all_tied"
    assert indexed.loc["trait-b", "ora_direction_vote_preference"] == "none"
    assert indexed.loc["trait-empty", "ora_direction_status"] == "no_background_taxa"
    for feature in ("trait-b", "trait-empty"):
        for column in ("ora_direction_p", "ora_direction_bh_q", "ora_direction_by_q"):
            assert pd.isna(indexed.loc[feature, column])
    assert pd.isna(indexed.loc["trait-empty", "ora_enriched_inc"])
    assert pd.isna(indexed.loc["trait-empty", "ora_enriched_dec"])
    assert indexed.loc["trait-empty", "ora_background_k"] == 0


def test_ora_correction_receives_literal_one_for_each_untestable_trait(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Check the complete correction input, before downstream handling of missing p (#1077)."""
    sources = ora_sources(tmp_path / "source", {f"S{i}": [(1, 0)] for i in range(9)})
    seen: dict[str, tuple[list[float], int | None]] = {}

    def record_correction(
        values: list[float], *, family_size: int | None = None, method: str = "bh"
    ) -> list[float]:
        assert method not in seen
        seen[method] = (list(values), family_size)
        return benjamini_hochberg(values, family_size=family_size, method=method)

    monkeypatch.setattr(ora, "benjamini_hochberg", record_correction)
    ora_table(sources)
    assert set(seen) == {"bh", "by"}
    for probabilities, family_size in seen.values():
        assert family_size == len(FEATURES)
        assert probabilities == [1 / 256, 1.0, 1 / 256, 1.0], (
            "Untestable ORA hypotheses must enter corrections at p=1"
        )


def test_ora_empty_native_family_is_explicitly_untestable(tmp_path: Path) -> None:
    sources = ora_sources(tmp_path / "source", {"S": [(None, None)]})
    table, _ = ora_table(sources)
    row = table.set_index("trait").loc["trait-a"]
    assert row.ora_direction_status == "no_eligible_pairs"
    assert row.ora_direction_pairs_eligible == 0
    assert row.ora_direction_pairs_excluded == 1
    assert pd.isna(row.ora_direction_p)
    assert pd.isna(row.ora_direction_mean_study_delta)


def test_ora_reports_list_size_balance_for_only_direction_eligible_pairs(tmp_path: Path) -> None:
    sources = ora_sources(
        tmp_path / "source",
        {"A": [(1, 0), (1, 0)], "B": [(0, 1)], "C": [(0, 0)], "D": [(1, None)]},
        unpaired=True,
        draw_sizes={
            ("A", "e0", "inc"): 8,
            ("A", "e1", "inc"): 8,
            ("B", "e0", "dec"): 8,
            ("C", "e0", "inc"): 5,
            ("C", "e0", "dec"): 5,
            ("D", "e0", "inc"): 8,
        },
    )
    table, _ = ora_table(sources)
    indexed = table.set_index("trait")
    columns = [
        "ora_direction_pairs_n_inc_gt_dec",
        "ora_direction_pairs_n_inc_lt_dec",
        "ora_direction_pairs_n_inc_eq_dec",
    ]
    for trait in ("trait-a", "trait-b", "trait-c"):
        assert [indexed.at[trait, column] for column in columns] == [2, 1, 1], (
            "List-size balance must retain increased/decreased ordering"
        )
        assert indexed.at[trait, "ora_direction_pairs_eligible"] == 4
        assert indexed.at[trait, "ora_direction_pairs_total"] == 5
        assert indexed.at[trait, "ora_direction_pairs_excluded"] == 1
    assert [indexed.at["trait-empty", column] for column in columns] == [0, 0, 0]


def test_ora_nullable_counts_preserve_integer_text_and_missing_background(tmp_path: Path) -> None:
    sources = ora_sources(tmp_path / "source", {f"S{i}": [(1, 0)] for i in range(3)})
    table, _ = ora_table(sources)
    columns = ("ora_enriched_inc", "ora_enriched_dec", "ora_depleted_inc", "ora_depleted_dec")
    assert all(str(table[column].dtype) == "Int64" for column in columns), (
        "ORA call counts must use nullable integers"
    )
    exported = {
        row["trait"]: row
        for row in csv.DictReader(io.StringIO(table.to_csv(sep="\t", index=False)), delimiter="\t")
    }
    assert exported["trait-a"]["ora_enriched_inc"] == "3"
    assert exported["trait-a"]["ora_enriched_dec"] == "0"
    assert exported["trait-c"]["ora_enriched_dec"] == "3"
    assert exported["trait-a"]["ora_depleted_inc"] == "0"
    assert exported["trait-a"]["ora_depleted_dec"] == "0"
    assert all(exported["trait-empty"][column] == "" for column in columns)


def test_ora_balanced_opposing_votes_are_mixed_not_none(tmp_path: Path) -> None:
    sources = ora_sources(tmp_path / "source", {"A": [(1, 0)], "B": [(0, 1)]})
    table, _ = ora_table(sources)
    row = table.set_index("trait").loc["trait-a"]
    assert row.ora_direction_vote_preference == "mixed"
    assert (row.ora_direction_studies_up, row.ora_direction_studies_down) == (1, 1)


def test_ora_missing_background_counts_are_blank_not_observed_zero(tmp_path: Path) -> None:
    sources = ora_sources(tmp_path / "source", {"A": [(1, 0)]})
    table, _ = ora_table(sources)
    row = table.set_index("trait").loc["trait-empty"]
    for column in ("ora_enriched_inc", "ora_enriched_dec", "ora_depleted_inc", "ora_depleted_dec"):
        assert pd.isna(row[column])
    assert row.ora_direction_vote_preference == "unavailable"


def test_ora_lower_tail_is_separate_from_enrichment_preference(tmp_path: Path) -> None:
    sources = ora_sources(
        tmp_path / "source", {"S": [(0, 0)]}, depleted=frozenset({("S", "e0", "inc")})
    )
    table, _ = ora_table(sources)
    row = table.set_index("trait").loc["trait-b"]
    assert (row.ora_depleted_inc, row.ora_depleted_dec) == (1, 0)
    assert (row.ora_enriched_inc, row.ora_enriched_dec) == (0, 0)
    assert row.ora_direction_status == "all_tied"
    assert row.ora_direction_studies_tied == 1
    assert pd.isna(row.ora_direction_p)


def test_ora_unpaired_profile_never_enters_native_counts_or_direction(tmp_path: Path) -> None:
    sources = ora_sources(tmp_path / "source", {"S": [(1, 0)]}, unpaired=True)
    table, _ = ora_table(sources)
    row = table.set_index("trait").loc["trait-a"]
    assert row.ora_enriched_inc == 1
    assert row.ora_native_testable_inc_profiles == 1
    assert row.ora_direction_studies_total == 1
    assert row.ora_direction_pairs_total == 1


@pytest.mark.parametrize("damage", ["duplicate", "missing", "extra", "q", "hit", "threshold"])
def test_ora_rejects_incomplete_or_inconsistent_native_family(tmp_path: Path, damage: str) -> None:
    sources = ora_sources(tmp_path / "source", {"S": [(1, 0)]})
    path = recorded_arm_output(
        sources.arms["ilyas"], "enrichment", "trait_enrichment_all_results.tsv"
    )
    records = list(rows(path))
    comments = ora_comments(path)
    if damage == "duplicate":
        records[-1] = records[0].copy()
    elif damage == "missing":
        records.pop()
    elif damage == "extra":
        records[-1]["trait"] = "not-declared"
    elif damage == "q":
        records[0]["fdr_q"] = "0.5"
    elif damage == "hit":
        records[0]["is_significant"] = "False"
    else:
        comments["P-Threshold"] = "0.05"
    write_table(path, records, comments)
    sources = replace(sources, recorded_hashes={})
    expected = {
        "duplicate": "duplicate coordinates",
        "missing": "truncated or has unexpected members",
        "extra": "undeclared coordinate",
        "q": "BH disagrees",
        "hit": "hit flag disagrees",
        "threshold": "P-Threshold",
    }
    with pytest.raises(ValueError, match=expected[damage]):
        ora_table(sources)


def test_ora_rejects_source_hash_change(tmp_path: Path) -> None:
    sources = ora_sources(tmp_path / "source", {"S": [(1, 0)]})
    path = recorded_arm_output(
        sources.arms["ilyas"], "enrichment", "trait_enrichment_all_results.tsv"
    )
    with path.open("a", encoding="utf-8") as handle:
        handle.write("# altered after verification\n")
    with pytest.raises(ValueError, match=r"hash|bytes"):
        ora_table(sources)


def test_ora_calibration_label_never_claims_certification() -> None:
    assert CALIBRATION_STATUS == "uncalibrated (methods/calibration.md excludes native ORA)"


def test_ora_requires_the_complete_depletion_source(tmp_path: Path) -> None:
    sources = ora_sources(tmp_path / "source", {"S": [(1, 0)]})
    path = recorded_arm_output(
        sources.arms["ilyas"], "depletion", "trait_depletion_all_results.tsv"
    )
    path.unlink()
    with pytest.raises(FileNotFoundError):
        ora_table(sources)
