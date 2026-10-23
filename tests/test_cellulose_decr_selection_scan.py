"""The #835 selection scan: every condition triple against cellulose_decr's U6 null."""

from __future__ import annotations

import itertools
import json
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mtp.cellulose_decr.analysis import (
    CONDITIONS,
    FECES,
    TRAIT_ID,
    compartment_of,
    condition_sets,
    set_summary,
    signature_null,
    trait_taxa,
    u6_units,
)
from mtp.cellulose_decr.selection_scan import (
    ALL_TRIPLES,
    CORE,
    CORE_PRODUCTS,
    DISJOINT,
    MEAN_DISTANCE,
    NO_TAXA,
    OUTPUTS,
    OWN,
    POOL_TOO_SMALL,
    PROTOCOL_PATH,
    SHARED,
    TIER1_2,
    bin_counts,
    block_positions,
    choose_traits,
    condition_table,
    lead_and_distribution,
    lead_rank,
    load_protocol,
    main,
    prepare,
    scan,
    triple_table,
)
from mtp.paths import file_sha256
from mtp.trait_screen import ALL_TIED, TESTED

PD, HIV, CRC = (CONDITIONS[short] for short in ("PD", "HIV", "CRC"))
LEAD = (PD, HIV, CRC)
OTHER_TRAIT = "METPO:1000608"

Signature = tuple[str, str, str, str, Sequence[str], Sequence[str]]


def _frames(records: Sequence[Signature]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Profiles and taxa as ``load_population`` builds them, from paired signatures.

    Each record is (study, experiment, condition, body site, decreased taxa, increased taxa).
    """
    profiles: list[dict[str, str]] = []
    taxa: list[dict[str, str]] = []
    for study, experiment, condition, site, decreased, increased in records:
        for direction, members in (("dec", decreased), ("inc", increased)):
            profile = f"{study}|{experiment}|{direction}"
            row = {
                "study_id": study,
                "experiment_id": experiment,
                "condition": condition,
                "body_site": site,
                "compartment": compartment_of(site),
                "dir": direction,
            }
            profiles.append({"profile_id": profile, **row})
            taxa.extend(
                {"profile_id": profile, "taxon_id": taxon, **row}
                for taxon in dict.fromkeys(members)
            )
    profile_frame = pd.DataFrame(profiles).set_index("profile_id")
    columns = ["profile_id", "taxon_id", "study_id", "experiment_id", "condition"]
    taxa_frame = pd.DataFrame(taxa)[[*columns, "compartment", "dir"]]
    return profile_frame, taxa_frame


def _matrix(n_taxa: int, traits: Mapping[str, Sequence[int]]) -> pd.DataFrame:
    """A matrix of taxa t0..t{n-1}, each trait borne by the listed taxa."""
    index = pd.Index([f"t{i}" for i in range(n_taxa)], name="taxon_id")
    frame = pd.DataFrame({"taxon_name": index, "taxonomic_rank": "genus"}, index=index)
    for trait, members in traits.items():
        frame[trait] = 0
        frame.loc[[f"t{i}" for i in members], trait] = 1
    return frame


def _random_corpus(
    seed: int,
    studies: Mapping[str, int],
    *,
    n_taxa: int = 40,
    size: int = 8,
    planted: Mapping[str, Sequence[int]] | None = None,
) -> list[Signature]:
    """Stool signatures with random taxa; ``planted`` adds fixed taxa to a condition's."""
    rng = np.random.default_rng(seed)
    records: list[Signature] = []
    number = 0
    for condition, count in studies.items():
        for _ in range(count):
            number += 1
            decreased = [f"t{i}" for i in rng.choice(n_taxa, size, replace=False)]
            decreased += [f"t{i}" for i in (planted or {}).get(condition, ())]
            increased = [f"t{i}" for i in rng.choice(n_taxa, size, replace=False)]
            study = f"bsdb:{number}"
            records.append((study, f"{study}-1", condition, FECES, decreased, increased))
    return records


# --------------------------------------------------------------------------- the rule


def test_protocol_holds_the_declared_rule() -> None:
    """The parameters the code reads are the ones methods/ declares (#835)."""
    protocol = load_protocol()
    assert protocol.min_studies == 3
    assert protocol.top_n == 10
    assert protocol.permutations == 100_000
    assert protocol.seed == 20261022
    assert protocol.tolerance == 1e-12
    assert protocol.bins == 20
    assert protocol.lead_triple == LEAD
    assert protocol.lead_trait == TRAIT_ID == protocol.always_included


@pytest.mark.parametrize(
    ("section", "key", "value", "message"),
    [
        (None, "protocol", "cellulose-decr-selection-scan-v0", "but this code is"),
        ("population", "core_contrast_population", "all_paired", "another population"),
        (None, "lead_triple", [PD, HIV, "Asthma"], "another lead"),
    ],
)
def test_protocol_is_refused_when_it_changes_what_the_code_assumes(
    tmp_path: Path, section: str | None, key: str, value: object, message: str
) -> None:
    data = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    (data if section is None else data[section])[key] = value
    path = tmp_path / "protocol.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        load_protocol(path)


def test_eligibility_counts_studies_with_a_stool_comparison() -> None:
    """Three stool comparisons from one study are one cohort; other body sites, blank
    labels and fewer than three studies do not make a condition eligible."""
    records: list[Signature] = [
        *(("bsdb:1", f"bsdb:1-{i}", "Asthma", FECES, ["t0"], ["t1"]) for i in range(3)),
        *((f"bsdb:{i}", f"bsdb:{i}-1", PD, FECES, ["t0"], ["t1"]) for i in range(2, 5)),
        *((f"bsdb:{i}", f"bsdb:{i}-1", HIV, FECES, ["t0"], ["t1"]) for i in range(5, 7)),
        *((f"bsdb:{i}", f"bsdb:{i}-1", HIV, "Saliva", ["t0"], ["t1"]) for i in range(7, 9)),
        *((f"bsdb:{i}", f"bsdb:{i}-1", "", FECES, ["t0"], ["t1"]) for i in range(9, 13)),
    ]
    profiles, _ = _frames(records)
    table = condition_table(profiles, 3).set_index("condition")
    assert sorted(table.index) == ["Asthma", HIV, PD]
    counts = table.loc["Asthma", ["stool_comparisons", "stool_studies", "u6_units"]]
    assert list(counts) == [3, 1, 1]
    assert table.loc[HIV, "stool_studies"] == 2
    assert table["eligible"].to_dict() == {"Asthma": False, HIV: False, PD: True}


def test_trait_set_is_the_top_decreasing_traits_and_cellulose() -> None:
    """Decreasing traits in the screen report's order, then cellulose, flagged (#835)."""
    rows = [
        # trait, up, down, tied, p, status
        ("up_strong", 30, 2, 0, 1e-8, TESTED),
        ("down_a", 2, 30, 0, 1e-6, TESTED),
        ("down_b", 3, 29, 0, 1e-5, TESTED),
        ("down_c_more_scored", 5, 40, 0, 1e-5, TESTED),
        ("tied_out", 0, 0, 9, None, ALL_TIED),
        ("down_d", 10, 12, 0, 0.8, TESTED),
        (TRAIT_ID, 12, 14, 0, 0.84, TESTED),
    ]
    screen = pd.DataFrame(
        [
            {
                "trait_id": trait,
                "trait_label": trait,
                "studies_up": up,
                "studies_down": down,
                "studies_tied": tied,
                "studies_scored": up + down,
                "share_up": up / (up + down) if up + down else None,
                "screen_p": p,
                "screen_test_status": status,
            }
            for trait, up, down, tied, p, status in rows
        ]
    )
    chosen = choose_traits(screen, 3, TRAIT_ID)
    assert chosen["trait_id"].tolist() == ["down_a", "down_c_more_scored", "down_b", TRAIT_ID]
    assert chosen["in_top_n"].tolist() == [True, True, True, False]
    assert chosen["decreasing_rank"].tolist() == [1, 2, 3, 5]
    # Already among the top N, it is not added twice.
    assert choose_traits(screen, 5, TRAIT_ID)["trait_id"].tolist().count(TRAIT_ID) == 1


# --------------------------------------------------------------------------- U6 reuse


# cellulose_decr's NullDraws.add takes a nanmean, which warns on the draws without a
# distance that the sparse trait creates on purpose.
@pytest.mark.filterwarnings("ignore:Mean of empty slice:RuntimeWarning")
def test_scan_reuses_cellulose_decr_u6_exactly() -> None:
    """Laid out as (PD, HIV, CRC), the scan's triple is cellulose_decr's U6, draw for draw.

    The sets equal ``condition_sets``, the statistics ``set_summary``, and the null the
    summary of ``signature_null`` with the same seed and draws. Other conditions are in
    the pool but not in the triple, as in cellulose_decr. A sparse trait, borne by one
    taxon that one study of each condition reports, gives draws with empty synthetic
    sets: pairs without a distance, and draws with none.
    """
    studies = {PD: 9, HIV: 5, CRC: 7, "Asthma": 6, "Autism": 2}
    records = _random_corpus(5, studies, n_taxa=40, size=6)
    for first in (0, 9, 14):
        study, experiment, condition, site, decreased, increased = records[first]
        records[first] = (study, experiment, condition, site, [*decreased, "t40"], increased)
    # One study reports two conditions, so a study can carry two units.
    records.append(("bsdb:7", "bsdb:7-2", "Asthma", FECES, ["t1", "t2", "t30"], ["t3"]))
    profiles, taxa = _frames(records)
    sparse = "METPO:1000623"
    traits = {TRAIT_ID: range(0, 40, 2), OTHER_TRAIT: range(1, 40, 3), sparse: [40]}
    matrix = _matrix(41, traits)
    permutations, seed = 3000, 20261023
    result = scan(
        profiles,
        taxa,
        matrix,
        list(traits),
        LEAD,
        permutations=permutations,
        seed=seed,
        tolerance=1e-12,
    )
    table = triple_table(result, LEAD).set_index("trait_id")
    assert table.loc[[sparse], "u6_draws_without_distance"].gt(0).all()
    for trait in traits:
        universe = trait_taxa(matrix, trait)
        sets = condition_sets(taxa, universe)
        observed = set_summary(sets)
        units, unit_study, counts, _ = u6_units(profiles, taxa, universe)
        draws = signature_null(np.random.default_rng(seed), units, unit_study, counts, permutations)
        null = draws.summary(observed)
        row = table.loc[trait]
        assert [row["n_a"], row["n_b"], row["n_c"]] == [len(sets[s]) for s in ("PD", "HIV", "CRC")]
        assert row["three_way_core"] == observed["three_way_core"]
        assert row["three_way_union"] == observed["three_way_union"]
        assert row[MEAN_DISTANCE] == observed["mean_pairwise_jaccard_distance"]
        mean, core = null["mean_pairwise_jaccard_distance"], null["three_way_core"]
        assert row["u6_p_mean_distance"] == mean["p"]
        assert row["u6_p_core"] == core["p"]
        assert row["u6_expected_core"] == pytest.approx(core["expected"], rel=1e-12)
        # cellulose_decr averages over every draw, so a draw without a distance makes its
        # expected value missing; the scan averages over the draws that have one.
        distances = np.concatenate(draws.mean_distance)
        assert row["u6_draws_without_distance"] == np.isnan(distances).sum()
        assert row["u6_expected_mean_distance"] == pytest.approx(np.nanmean(distances), rel=1e-12)
        assert row["u6_sd_mean_distance"] == pytest.approx(np.nanstd(distances, ddof=1), rel=1e-9)
        if row["u6_draws_without_distance"] == 0:
            assert row["u6_expected_mean_distance"] == pytest.approx(mean["expected"], rel=1e-12)
            assert row["u6_z_mean_distance"] == pytest.approx(mean["z"], rel=1e-9)
        assert [row["u6_units_a"], row["u6_units_b"], row["u6_units_c"]] == [
            counts[s] for s in ("PD", "HIV", "CRC")
        ]
        assert row["u6_draws"] == permutations
        assert row["u6_draws_source"] == SHARED


def test_signature_draws_leave_u6_unchanged() -> None:
    """``signature_null`` takes its draws from ``signature_draws`` (#835) and still gives the
    draws it gave before that refactor. The values were recorded on the code before it."""
    rng = np.random.default_rng(835)
    per_study = rng.integers(1, 4, 40)
    unit_study = np.repeat(np.arange(40), per_study).astype(np.int64)
    units = rng.random((unit_study.size, 25)) < 0.2
    counts = {"PD": 9, "HIV": 5, "CRC": 7}
    draws = signature_null(np.random.default_rng(20261023), units, unit_study, counts, 12_345)
    core = np.concatenate(draws.core)
    mean = np.concatenate(draws.mean_distance)
    assert int(core.sum()) == 131273
    assert core[-5:].tolist() == [7, 12, 11, 10, 9]
    assert {short: int(np.concatenate(values).sum()) for short, values in draws.sizes.items()} == {
        "PD": 264116,
        "HIV": 201696,
        "CRC": 239651,
    }
    assert float(mean.sum()) == 4912.902420421743
    assert mean[:3].tolist() == [0.31417654808959156, 0.38166666666666665, 0.4105555555555555]


def test_every_triple_draws_distinct_studies() -> None:
    """Each condition takes its own block of the shared draws, so a triple never draws one
    study twice. Every study carries only taxa of its own, so any shared study would make
    a synthetic pair overlap; the null distance must be exactly 1 in every draw."""
    records: list[Signature] = []
    for number in range(1, 25):
        condition = ("A", "B", "C", "D", "E")[number % 5]
        study = f"bsdb:{number}"
        own = [f"t{2 * number}", f"t{2 * number + 1}"]
        records.append((study, f"{study}-1", condition, FECES, own, ["t0"]))
    profiles, taxa = _frames(records)
    matrix = _matrix(60, {TRAIT_ID: range(60)})
    conditions = ["A", "B", "C", "D", "E"]
    result = scan(
        profiles, taxa, matrix, [TRAIT_ID], conditions, permutations=500, seed=1, tolerance=1e-12
    )
    table = triple_table(result, LEAD)
    assert set(table["u6_draws_source"]) == {SHARED}
    assert (table["u6_expected_mean_distance"] == 1).all()
    # Every draw ties the observed distance of 1, and a tie counts as extreme.
    assert (table["u6_p_mean_distance"] == 1).all()
    assert (table["u6_sd_mean_distance"] == 0).all()
    assert (table["u6_expected_core"] == 0).all()


def test_overlapping_blocks_get_their_own_u6_draws() -> None:
    """When the conditions have more units than the pool has studies, the blocks wrap and
    a triple whose blocks overlap is drawn on its own, with U6's code and seed
    (seed, triple index), laid out as (a, b, c)."""
    labels = ["A", "B", "C", "D", "E"]
    records = _random_corpus(7, dict.fromkeys(labels, 3), n_taxa=30, size=5)
    # The first study of each condition also reports the next condition: 20 units, 15 studies.
    nexts = dict(zip(labels, labels[1:] + labels[:1], strict=True))
    records += [
        (study, f"{study}-2", nexts[condition], site, dec, inc)
        for study, _, condition, site, dec, inc in records[::3]
    ]
    profiles, taxa = _frames(records)
    matrix = _matrix(30, {TRAIT_ID: range(0, 30, 2)})
    seed, permutations = 11, 1500
    prepared = prepare(profiles, taxa, matrix, [TRAIT_ID], labels)
    assert prepared.units.tolist() == [4, 4, 4, 4, 4]
    result = scan(
        profiles,
        taxa,
        matrix,
        [TRAIT_ID],
        labels,
        permutations=permutations,
        seed=seed,
        tolerance=1e-12,
    )
    table = triple_table(result, LEAD)
    assert set(table["u6_draws_source"]) == {SHARED, OWN}
    assert set(table["test_status"]) == {TESTED}
    universe = trait_taxa(matrix, TRAIT_ID)
    units, unit_study, _, _ = u6_units(profiles, taxa, universe)
    for index in range(len(table)):
        row = table.iloc[index]
        triple = [row["condition_a"], row["condition_b"], row["condition_c"]]
        names = dict(zip(triple, LEAD, strict=True))
        relabelled = taxa.assign(condition=taxa["condition"].replace(names))
        observed = set_summary(condition_sets(relabelled, universe))
        if row["u6_draws_source"] == OWN:
            rng = np.random.default_rng([seed, index])
            counts = dict.fromkeys(("PD", "HIV", "CRC"), 4)
            null = signature_null(rng, units, unit_study, counts, permutations).summary(observed)
            assert row["u6_p_mean_distance"] == null["mean_pairwise_jaccard_distance"]["p"]
            assert row["u6_p_core"] == null["three_way_core"]["p"]
        assert row[MEAN_DISTANCE] == observed["mean_pairwise_jaccard_distance"]


def test_study_disjoint_reference_drops_triples_that_share_a_study() -> None:
    """A study reporting two of a triple's conditions is counted, and the declared
    sensitivity ranks only the triples whose conditions share no study."""
    labels = ["A", "B", "C", "D", "E"]
    records = _random_corpus(9, dict.fromkeys(labels, 4), n_taxa=30, size=5)
    # The first study of A and of C also reports B.
    records += [
        (study, f"{study}-2", "B", site, dec, inc)
        for study, _, _, site, dec, inc in (records[0], records[8])
    ]
    profiles, taxa = _frames(records)
    matrix = _matrix(30, {TRAIT_ID: range(30)})
    result = scan(
        profiles, taxa, matrix, [TRAIT_ID], labels, permutations=200, seed=1, tolerance=1e-12
    )
    table = triple_table(result, LEAD)
    studies = {label: {r[0] for r in records if r[2] == label} for label in labels}
    for _, row in table.iterrows():
        members = [studies[row[f"condition_{x}"]] for x in "abc"]
        shared = {s for s in set.union(*members) if sum(s in m for m in members) >= 2}
        assert row["shared_studies"] == len(shared)
    assert set(table["shared_studies"]) == {0, 1, 2}
    _, histogram, quantiles = lead_and_distribution(table, [TRAIT_ID], 20)
    disjoint = int(table["shared_studies"].eq(0).sum())
    assert quantiles[TRAIT_ID][MEAN_DISTANCE][DISJOINT]["triples"] == disjoint < len(table)
    counted = histogram[histogram["reference_set"].eq(DISJOINT)].groupby("statistic")["triples"]
    assert counted.sum().tolist() == [disjoint, disjoint]


def test_a_triple_larger_than_the_pool_is_not_tested() -> None:
    records = _random_corpus(3, {"A": 3, "B": 3, "C": 3}, n_taxa=10, size=3)
    records += [
        (study, f"{study}-{other}", other, site, dec, inc)
        for study, _, condition, site, dec, inc in list(records)
        for other in ("A", "B", "C")
        if other != condition
    ]
    profiles, taxa = _frames(records)
    matrix = _matrix(10, {TRAIT_ID: range(10)})
    result = scan(
        profiles, taxa, matrix, [TRAIT_ID], ["A", "B", "C"], permutations=50, seed=1, tolerance=0
    )
    # Each condition has a unit in all 9 studies, and the triple would need 27.
    table = triple_table(result, LEAD)
    assert table["test_status"].tolist() == [POOL_TOO_SMALL]
    assert table["u6_p_mean_distance"].isna().all()


def test_block_positions_wrap_at_the_end_of_the_order() -> None:
    blocks = block_positions([3, 2, 4], 7)
    assert [block.tolist() for block in blocks] == [[0, 1, 2], [3, 4], [5, 6, 0, 1]]


# --------------------------------------------------------------------------- the scan


def test_planted_overlap_ranks_first() -> None:
    """Three conditions whose studies all report the same rare annotated taxa decreased
    overlap more than U6 expects; the scan must rank that triple first on both statistics.

    U6 draws whole studies from the pool, planted studies included, so the planted
    studies are kept a small share of it: 100 conditions of two studies each fill the
    pool without entering any triple.
    """
    rng = np.random.default_rng(13)
    eligible = ["A", "B", "C", "D", "E", "F"]
    planted = [f"t{i}" for i in range(500, 510)]
    records: list[Signature] = []
    background = [f"X{i}" for i in range(100)]
    for number, condition in enumerate([*np.repeat(eligible, 3), *np.repeat(background, 2)]):
        decreased = [f"t{i}" for i in rng.choice(500, 4, replace=False)]
        if condition in ("B", "D", "F"):
            decreased += planted
        study = f"bsdb:{number + 1}"
        records.append((study, f"{study}-1", str(condition), FECES, decreased, ["t0"]))
    profiles, taxa = _frames(records)
    matrix = _matrix(510, {TRAIT_ID: range(510)})
    assert sorted(condition_table(profiles, 3).query("eligible")["condition"]) == eligible
    result = scan(
        profiles, taxa, matrix, [TRAIT_ID], eligible, permutations=2000, seed=4, tolerance=1e-12
    )
    table = triple_table(result, LEAD)
    assert len(table) == 20
    chosen = table["condition_a"].eq("B") & table["condition_b"].eq("D")
    planted_row = table[chosen & table["condition_c"].eq("F")].iloc[0]
    assert planted_row["rank_mean_distance"] == 1
    assert planted_row["rank_core"] == 1
    others = table.drop(index=planted_row.name)
    assert planted_row["u6_p_mean_distance"] < others["u6_p_mean_distance"].min()
    assert planted_row["u6_p_core"] < 0.01 < others["u6_p_core"].min()


def _lead_corpus() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    studies = {PD: 6, HIV: 4, CRC: 5, "Asthma": 5, "Autism": 4, "Obesity": 3}
    records = _random_corpus(21, studies, n_taxa=30, size=6)
    profiles, taxa = _frames(records)
    return profiles, taxa, _matrix(30, {TRAIT_ID: range(0, 30, 2), OTHER_TRAIT: range(1, 30, 2)})


def test_scan_is_deterministic_under_the_seed() -> None:
    profiles, taxa, matrix = _lead_corpus()
    conditions = sorted(condition_table(profiles, 3).query("eligible")["condition"])

    def run(seed: int) -> pd.DataFrame:
        result = scan(
            profiles,
            taxa,
            matrix,
            [TRAIT_ID, OTHER_TRAIT],
            conditions,
            permutations=1200,
            seed=seed,
            tolerance=1e-12,
        )
        return triple_table(result, LEAD)

    first, again, other = run(8), run(8), run(9)
    pd.testing.assert_frame_equal(first, again)
    assert not first["u6_p_mean_distance"].equals(other["u6_p_mean_distance"])


def test_lead_rank_reports_ties_as_a_range() -> None:
    p_values = pd.Series([0.01, 0.2, 0.2, 0.2, 0.5])
    assert lead_rank(p_values, 0.2) == {
        "triples": 5,
        "rank_min": 2,
        "rank_max": 4,
        "share_at_most": 0.8,
    }
    # A lead outside the reference set is placed among it as one more triple.
    assert lead_rank(pd.Series([0.01, 0.5]), 0.2, counted=False) == {
        "triples": 3,
        "rank_min": 2,
        "rank_max": 2,
        "share_at_most": 2 / 3,
    }


def test_p_value_bins_are_closed_above() -> None:
    counts = bin_counts(pd.Series([0.05, 0.0500001, 0.1, 1.0, 1e-5]), 20)
    assert counts[0] == 2
    assert counts[1] == 2
    assert counts[19] == 1
    assert sum(counts) == 5


def test_lead_triple_is_ranked_in_both_reference_sets() -> None:
    profiles, taxa, matrix = _lead_corpus()
    conditions = sorted(condition_table(profiles, 3).query("eligible")["condition"])
    result = scan(
        profiles,
        taxa,
        matrix,
        [TRAIT_ID, OTHER_TRAIT],
        conditions,
        permutations=800,
        seed=2,
        tolerance=1e-12,
    )
    triples = triple_table(result, LEAD)
    assert triples["is_lead_triple"].sum() == 2
    assert len(triples) == 2 * len(list(itertools.combinations(conditions, 3)))
    lead, histogram, quantiles = lead_and_distribution(triples, [TRAIT_ID, OTHER_TRAIT], 20)
    row = lead.set_index("trait_id").loc[TRAIT_ID]
    mine = triples[triples["trait_id"].eq(TRAIT_ID)]
    lead_p = mine.loc[mine["is_lead_triple"], "u6_p_mean_distance"].iloc[0]
    tested = mine[mine["test_status"].eq(TESTED)]
    assert (
        row[f"{MEAN_DISTANCE}_{ALL_TRIPLES}_rank_min"]
        == (tested["u6_p_mean_distance"] < lead_p).sum() + 1
    )
    assert row[f"{MEAN_DISTANCE}_{ALL_TRIPLES}_triples"] == len(tested)
    # No study reports two conditions here, so every triple is study-disjoint.
    assert row[f"{CORE}_{DISJOINT}_triples"] == len(tested)
    assert row["lead_condition_a"] == CRC
    counts = histogram[
        histogram["trait_id"].eq(TRAIT_ID)
        & histogram["statistic"].eq(MEAN_DISTANCE)
        & histogram["reference_set"].eq(ALL_TRIPLES)
    ]["triples"]
    assert counts.sum() == len(tested)
    assert quantiles[TRAIT_ID][MEAN_DISTANCE][ALL_TRIPLES]["triples"] == len(tested)


def test_triples_whose_sets_are_all_empty_are_not_tested() -> None:
    records = _random_corpus(4, {"A": 3, "B": 3, "C": 3, "D": 3}, n_taxa=10, size=3)
    profiles, taxa = _frames(records)
    # The trait is borne only by a taxon no signature reports.
    matrix = _matrix(12, {TRAIT_ID: [11]})
    result = scan(
        profiles,
        taxa,
        matrix,
        [TRAIT_ID],
        ["A", "B", "C", "D"],
        permutations=50,
        seed=1,
        tolerance=1e-12,
    )
    table = triple_table(result, LEAD)
    assert set(table["test_status"]) == {NO_TAXA}
    assert table["u6_p_core"].isna().all()
    assert table["rank_mean_distance"].isna().all()


# --------------------------------------------------------------------------- end to end

RAW_HEADER = [
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
    "Study design",
    "PMID",
    "DOI",
    "EFO ID",
    "UBERON ID",
    "Group 1 definition",
    "Group 0 sample size",
    "Group 1 sample size",
    "Description",
    "Sequencing type",
]
STUDIES = {PD: 4, HIV: 3, CRC: 3, "Asthma": 3, "Obesity": 3, "Autism": 2}


def _write_core_inputs(root: Path) -> dict[str, Path]:
    """A raw dump of paired stool case-control experiments, and a KG annotating its taxa.

    Species NCBITaxon:2001-2030 sit under one genus. Cellulose catabolism is asserted on
    the first fifteen and facultative anaerobiosis on the last twenty.
    """
    rng = np.random.default_rng(835)
    rows = []
    number = 0
    for condition, count in STUDIES.items():
        for _ in range(count):
            number += 1
            for signature, direction in ((1, "increased"), (2, "decreased")):
                species = sorted(rng.choice(np.arange(2001, 2031), 4, replace=False).tolist())
                rows.append(
                    {
                        "BSDB ID": f"bsdb:{number}/1/{signature}",
                        "Study": f"Study {number}",
                        "Title": "Synthetic selection-scan contrast",
                        "Year": "2026",
                        "Host species": "Homo sapiens",
                        "Body site": "Feces",
                        "Condition": condition,
                        "Group 0 name": "Healthy controls",
                        "Group 1 name": f"{condition} patients",
                        "Abundance in Group 1": direction,
                        "NCBI Taxonomy IDs": ";".join(map(str, species)),
                        "MetaPhlAn taxon names": ",".join(f"s__Example_{i}" for i in species),
                        "State": "Complete",
                        "Reviewer": "Fixture",
                        "Study design": "case-control",
                        "Group 1 definition": f"{condition} patients",
                        "Sequencing type": "16S",
                    }
                )
    raw = root / "full_dump.csv"
    frame = pd.DataFrame(rows).reindex(columns=RAW_HEADER).fillna("")
    raw.write_text(
        "# BugSigDB selection-scan fixture-835, License: CC0 synthetic data\n"
        + frame.to_csv(index=False, lineterminator="\n"),
        encoding="utf-8",
    )
    edges = ["id\tsubject\tpredicate\tobject\tprimary_knowledge_source"]
    nodes = ["id\tcategory\tname", "NCBITaxon:2000\tbiolink:OrganismTaxon\tExample"]
    for taxon in range(2001, 2031):
        edges.append(f"tax{taxon}\tNCBITaxon:{taxon}\tbiolink:subclass_of\tNCBITaxon:2000\tfixture")
        nodes.append(f"NCBITaxon:{taxon}\tbiolink:OrganismTaxon\tExample {taxon}")
        if taxon <= 2015:
            edges.append(f"cel{taxon}\tNCBITaxon:{taxon}\tMETPO:2000103\tGO:0030245\tfixture")
        if taxon > 2010:
            edges.append(
                f"fac{taxon}\tNCBITaxon:{taxon}\tbiolink:has_phenotype\t{OTHER_TRAIT}\tfixture"
            )
    nodes += [
        "METPO:2000103\tbiolink:OntologyClass\tcapable of",
        "GO:0030245\tbiolink:BiologicalProcess\tcellulose catabolic process",
        f"{OTHER_TRAIT}\tbiolink:PhenotypicQuality\tfacultatively anaerobic",
    ]
    (root / "merged-kg_edges.tsv").write_text("\n".join(edges) + "\n", encoding="utf-8")
    (root / "merged-kg_nodes.tsv").write_text("\n".join(nodes) + "\n", encoding="utf-8")
    return {
        "raw_dump": raw,
        "kg_edges": root / "merged-kg_edges.tsv",
        "kg_nodes": root / "merged-kg_nodes.tsv",
    }


@pytest.fixture(scope="module")
def cores(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """A Tier 1+2 core and an all-paired core of the same synthetic inputs."""
    from mtp.paths import REPO_ROOT
    from mtp.workflow.operations import OPERATIONS, validate_core_product
    from mtp.workflow.runtime import execute_workflow

    root = tmp_path_factory.mktemp("selection-scan-core")
    inputs = _write_core_inputs(root)
    runs = {}
    for name, spec, population in (
        ("tier1_2", "bugsigdb_core_subset.yaml", TIER1_2),
        ("all_paired", "bugsigdb_core.yaml", None),
    ):
        run = execute_workflow(
            REPO_ROOT / "workflows" / spec,
            inputs,
            root / name,
            OPERATIONS,
            parameter_overrides={
                "source_mode": "fixture",
                **({"contrast_population": population} if population else {}),
            },
            validator=validate_core_product,
            repository_root=REPO_ROOT,
        )
        runs[name] = Path(run.manifest_path)
    return runs


def test_main_writes_every_output_bound_to_the_core(cores: dict[str, Path], tmp_path: Path) -> None:
    """A run on a Tier 1+2 core writes the declared outputs, a provenance naming its inputs
    and rule, and a consumer record that verifies against the core."""
    from mtp.workflow.consumers import verify_consumer_record

    out = tmp_path / "scan"
    argv = ["--core-run", str(cores["tier1_2"]), "--output-dir", str(out)]
    assert main([*argv, "--permutations", "400", "--seed", "3", "--code-commit", "abc"]) == 0
    assert sorted(path.name for path in out.iterdir()) == sorted(
        [*OUTPUTS.values(), "core_consumption.json", "core_consumption.coverage.tsv"]
    )
    record = verify_consumer_record(out / "core_consumption.json")
    assert record.consumer == "mtp.cellulose_decr.selection_scan"

    conditions = pd.read_csv(out / OUTPUTS["conditions"], sep="\t")
    eligible = conditions[conditions["eligible"]]
    assert sorted(eligible["condition"]) == sorted(set(STUDIES) - {"Autism"})
    assert eligible["stool_studies"].ge(3).all()
    triples = pd.read_csv(out / OUTPUTS["triples"], sep="\t")
    traits = pd.read_csv(out / OUTPUTS["traits"], sep="\t")
    assert TRAIT_ID in set(traits["trait_id"])
    assert len(triples) == 10 * len(traits)
    assert triples["is_lead_triple"].sum() == len(traits)
    assert set(triples["u6_draws"]) <= {0, 400}

    summary = json.loads((out / OUTPUTS["summary"]).read_text(encoding="utf-8"))
    assert summary["follows_protocol"] is False
    assert summary["counts"]["eligible_conditions"] == 5
    assert summary["lead_triple"]["eligible"] == dict.fromkeys(LEAD, True)
    provenance = json.loads((out / OUTPUTS["provenance"]).read_text(encoding="utf-8"))
    assert provenance["code_commit"] == "abc"
    assert provenance["core_run"]["contrast_population"] == TIER1_2
    assert set(provenance["inputs"]) == set(CORE_PRODUCTS)
    assert provenance["protocol"]["sha256"] == file_sha256(PROTOCOL_PATH)
    assert (provenance["permutations"], provenance["seed"]) == (400, 3)

    # The same core, rule and seed give the same tables.
    again = tmp_path / "again"
    assert (
        main(
            [
                "--core-run",
                str(cores["tier1_2"]),
                "--output-dir",
                str(again),
                "--permutations",
                "400",
                "--seed",
                "3",
            ]
        )
        == 0
    )
    for name in ("conditions", "traits", "triples", "distribution", "summary"):
        assert (out / OUTPUTS[name]).read_bytes() == (again / OUTPUTS[name]).read_bytes(), name


def test_main_refuses_another_population_or_a_used_directory(
    cores: dict[str, Path], tmp_path: Path
) -> None:
    with pytest.raises(ValueError, match=TIER1_2):
        main(["--core-run", str(cores["all_paired"]), "--output-dir", str(tmp_path / "a")])
    used = tmp_path / "used"
    used.mkdir()
    (used / "old.tsv").write_text("x\n", encoding="utf-8")
    with pytest.raises(ValueError, match="new or empty"):
        main(["--core-run", str(cores["tier1_2"]), "--output-dir", str(used)])


def test_plan_reports_counts_without_writing(
    cores: dict[str, Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "plan"
    assert main(["--core-run", str(cores["tier1_2"]), "--output-dir", str(out), "--plan"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["eligible_conditions"] == 5
    assert report["triples"] == 10
    assert report["permutations"] == 100_000
    assert report["estimated_null_minutes"] > 0
    assert not out.exists()
