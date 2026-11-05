"""Bind calibration and native-influence MT declarations to real producers (#1494/#1495)."""

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control, hypergeom, ttest_1samp

from mtp.composition_calibration import compute_composition_contrasts, run_composition_calibration
from mtp.composition_calibration.composition import load_composition_trait_columns
from mtp.contribution_inputs import Base, Key
from mtp.density_adjustment import TABLE_1_TRAITS
from mtp.native_taxon_influence import correct_shared_comparisons, ora_cells, transformed
from mtp.taxon_influence import comparison
from mtp.workflow.downstream_contracts import get_standard_analysis_spec


def test_manifest_calibration_policy_matches_full_catalog_bh_and_by(tmp_path: Path) -> None:
    focal = list(TABLE_1_TRAITS)
    retained = [*focal, "METPO:extra_constant", "METPO:extra_zero"]
    media = "mediadive.medium:fixture"
    catalog = tmp_path / "catalog.tsv"
    pd.DataFrame(
        {
            "column_id": [*retained, media],
            "object_category": ["biolink:OntologyClass"] * len(retained) + ["METPO:1004005"],
        }
    ).to_csv(catalog, sep="\t", index=False)
    rows = []
    increased = [12, 9, 11, 13, 6, 15, 10, 8]
    decreased = [2, 3, 4, 2, 3, 4, 5, 4]
    for study, (inc, dec) in enumerate(zip(increased, decreased, strict=True)):
        for direction, value in (("inc", inc), ("dec", dec)):
            other = 2 + ((study + (3 if direction == "inc" else 0)) % 6)
            rows.append(
                {
                    "study_id": f"s{study}",
                    "experiment_id": "e1",
                    "inc_dec": direction,
                    focal[0]: value,
                    focal[1]: 20 - value,
                    focal[2]: other,
                    focal[3]: 10 - other,
                    focal[4]: 2,
                    focal[5]: 0,
                    "METPO:extra_constant": 3,
                    "METPO:extra_zero": 0,
                    media: 1000 + study,
                }
            )
    profiles = pd.DataFrame(rows)
    traits = load_composition_trait_columns(catalog, profiles)
    assert traits == retained
    assert len(traits) > 6
    summary, summary_rows, replicates, null_pvalues, table1 = run_composition_calibration(
        profiles, traits, replicates=3, seed=182
    )
    assert summary["n_traits"] == len(retained)
    assert summary_rows["n_tests"].eq(len(retained)).all()
    assert replicates["n_tests"].eq(len(retained)).all()
    assert set(replicates["correction"]) == {"BH", "BY"}
    assert set(null_pvalues["trait"]) == set(retained)
    assert (
        null_pvalues.loc[null_pvalues["trait"].eq("METPO:extra_zero"), "null_p_value"].eq(1).all()
    )

    contrasts = compute_composition_contrasts(profiles, traits).study_contrasts
    observed_p = np.array(
        [
            float(ttest_1samp(contrasts[trait], 0.0).pvalue)
            if contrasts[trait].std(ddof=1) > 1e-12
            else 1.0
            for trait in retained
        ]
    )
    bh = false_discovery_control(observed_p, method="bh")
    by = false_discovery_control(observed_p, method="by")
    harmonic = sum(1 / i for i in range(1, len(retained) + 1))
    assert by == pytest.approx(np.minimum(bh * harmonic, 1.0))
    assert any(0 < a < b < 1 for a, b in zip(bh, by, strict=True))
    quoted = table1.set_index("trait_id")
    assert set(quoted.index) == set(focal)
    for index, trait in enumerate(retained):
        if trait in quoted.index:
            assert quoted.loc[trait, "observed_p"] == pytest.approx(
                observed_p[index], rel=1e-12, abs=1e-15
            )
            assert quoted.loc[trait, "observed_bh_q"] == pytest.approx(
                bh[index], rel=1e-12, abs=1e-15
            )
            assert quoted.loc[trait, "observed_by_q"] == pytest.approx(
                by[index], rel=1e-12, abs=1e-15
            )
    for row in summary_rows.itertuples():
        q = bh if row.correction == "BH" else by
        assert row.observed_calls == int((q <= row.threshold).sum())
    for replicate, group in null_pvalues.groupby("replicate"):
        p = group.set_index("trait").loc[retained, "null_p_value"].to_numpy()
        for row in replicates.loc[replicates["replicate"].eq(replicate)].itertuples():
            assert isinstance(row.correction, str)
            q = false_discovery_control(p, method=row.correction.lower())
            assert row.significant == int((q <= row.threshold).sum())
    spec = get_standard_analysis_spec("jed_composition_arcsine_sqrt_ttest_calibration_v1")
    assert (
        spec.multiple_testing_policy == "BH_and_BY_q_across_retained_non_MediaDive_catalog_traits"
    )


def test_manifest_native_influence_policy_matches_four_selected_bh_families() -> None:
    profiles: dict[Key, dict[str, int]] = {}
    for study, (inc, dec) in enumerate(
        zip([12, 8, 15, 9, 13, 7, 14, 11], [3, 10, 8, 5, 6, 9, 4, 6], strict=True)
    ):
        for direction, value in (("inc", inc), ("dec", dec)):
            profiles[(f"s{study}", "e1", direction)] = {"A": value, "B": 20 - value, "C": 0}
    base = Base(["A", "B", "C"], {"a": "A", "c": "C"}, profiles, {}, {}, {}, {}, {}, [])
    keys = list(profiles)
    policy = {
        "minimum_profile_support": 1,
        "maximum_profile_support_fraction": 1.0,
        "minimum_pairs": 3,
        "transform_multiplier": 1000,
    }
    composition, contrasts, _ = transformed(base, profiles, keys, policy)
    tested = [row for row in composition.values() if row["status"] == "tested"]
    assert {row["trait_key"] for row in tested} == {"A", "B"}
    assert all(row["family_tests"] == 2 and row["family_selected_traits"] == 2 for row in tested)
    expected = false_discovery_control([row["p"] for row in tested], method="bh")
    assert any(0 < value < 1 for value in expected)
    for row, q in zip(tested, expected, strict=True):
        assert row["q"] == pytest.approx(q, rel=1e-12, abs=1e-15)
        for threshold in (0.001, 0.05):
            assert row[f"significant_{threshold}"] == bool(q <= threshold)
    assert composition["C"]["status"] == "not_selected"
    assert composition["C"]["p"] == composition["C"]["q"] == ""

    background = Counter(
        {trait: sum(row[trait] for row in profiles.values()) for trait in base.columns}
    )
    cells, family = ora_cells(base, profiles, keys, background)
    cell_index = {
        (row["study_id"], row["experiment_id"], row["direction"], row["trait_key"]): row
        for row in cells
    }
    total = sum(background.values())
    for tail in ("enrichment", "depletion"):
        identities = [
            (key, trait)
            for key in keys
            for trait in base.columns
            if tail == "depletion" or profiles[key][trait] > 0
        ]
        assert family[f"{tail}_tests"] == len(identities)
        pvalues = [
            float(
                hypergeom.sf(
                    profiles[key][trait] - 1, total, background[trait], sum(profiles[key].values())
                )
            )
            if tail == "enrichment"
            else float(
                hypergeom.cdf(
                    profiles[key][trait], total, background[trait], sum(profiles[key].values())
                )
            )
            for key, trait in identities
        ]
        adjusted = false_discovery_control(pvalues, method="bh")
        assert any(0 < value < 1 for value in adjusted)
        for (key, trait), p, q in zip(identities, pvalues, adjusted, strict=True):
            if trait not in base.traits.values():
                continue
            row = cell_index[(*key, trait)]
            assert row[f"{tail}_p"] == pytest.approx(p, rel=1e-12, abs=1e-15)
            assert row[f"{tail}_q"] == pytest.approx(q, rel=1e-12, abs=1e-15)
            for threshold in (0.001, 0.05):
                fold_pass = row["fold"] != "" and (
                    row["fold"] >= 1.5 if tail == "enrichment" else row["fold"] <= 0.67
                )
                assert row[f"{tail}_significant_{threshold}"] == bool(q <= threshold and fold_pass)
    assert family["enrichment_tests"] < family["depletion_tests"]

    shared = [
        {
            "scenario": "baseline",
            "scope": "pooled",
            "reference": "refitted",
            **comparison(values, {key: value / 2 for key, value in values.items()}),
        }
        for values in contrasts.values()
    ]
    correct_shared_comparisons(shared, ("scenario", "scope", "reference"))
    assert {row["shared_correction"] for row in shared} == {"BH3_directional_sign_test"}
    for population in ("before", "after", "before_matched", "after_matched", "matched_change"):
        raw = [row[f"{population}_p"] for row in shared]
        qvalues = false_discovery_control([p if p != "" else 1.0 for p in raw], method="bh")
        assert any(p == "" for p in raw)
        assert any(0 < q < 1 for q in qvalues)
        for row, p, q in zip(shared, raw, qvalues, strict=True):
            if p == "":
                assert row[f"{population}_q"] == ""
            else:
                assert row[f"{population}_q"] == pytest.approx(q, rel=1e-12, abs=1e-15)
    spec = get_standard_analysis_spec("historical_native_taxon_influence_experiment_weighted_v1")
    assert (
        spec.multiple_testing_policy
        == "BH_q_within_selected_native_influence_families_per_scenario_scope_tail"
    )
