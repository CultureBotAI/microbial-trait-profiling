"""The #835 selection scan: every condition triple against cellulose_decr's U6 null."""

from __future__ import annotations

import dataclasses
import itertools
import json
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

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
    NO_UNITS,
    NOT_DRAWN,
    OUTCOME_COLUMNS,
    OUTPUTS,
    OWN,
    POOL_TOO_SMALL,
    PRIMARY,
    PROTOCOL_PATH,
    SENSITIVITY,
    SHARED,
    TIER1_2,
    TRIPLE_COLUMNS,
    ScanInputs,
    bin_counts,
    block_positions,
    build_population,
    choose_traits,
    closest_disjoint_triples,
    condition_table,
    evaluate,
    lead_and_distribution,
    lead_rank,
    load_core_inputs,
    load_protocol,
    main,
    prepare,
    rank_candidates,
    rollup_scope,
    scan,
    select,
    triple_table,
)
from mtp.paths import file_sha256
from mtp.trait_guilds import Profile
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
    """The parameters the code reads are the ones methods/ declares (#835), with the
    amendments of 2026-10-22 applied (#872)."""
    protocol = load_protocol()
    assert protocol.min_studies == 3
    assert protocol.top_n == 10
    assert protocol.permutations == 100_000
    assert protocol.seed == 20261022
    assert protocol.tolerance == 1e-12
    assert protocol.bins == 20
    assert protocol.lead_triple == LEAD
    assert protocol.lead_trait == TRAIT_ID == protocol.always_included
    assert (protocol.scope_primary, protocol.scope_sensitivity) == ("off", ("phenotype",))
    assert protocol.scope_absent == "off"
    data = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    # The original declaration stays as committed; the amendments carry the changes.
    assert "strain_rollup" not in data
    assert "secondary_rank" not in data
    assert [amendment["date"] for amendment in data["amendments"]] == ["2026-10-22"] * 2


def _edited(tmp_path: Path, edit: Callable[[dict[str, Any]], None]) -> Path:
    data = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    edit(data)
    path = tmp_path / "protocol.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def _set(keys: tuple[str | int, ...], value: object) -> Callable[[dict[str, Any]], None]:
    def edit(data: dict[str, Any]) -> None:
        target: Any = data
        for key in keys[:-1]:
            target = target[key]
        target[keys[-1]] = value

    return edit


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        (_set(("protocol",), "cellulose-decr-selection-scan-v0"), "but this code is"),
        (_set(("population", "core_contrast_population"), "all_paired"), "this code does"),
        (_set(("population", "stool_body_site"), "Saliva"), "stool_body_site"),
        (_set(("population", "core_products"), ["experiments"]), "core_products"),
        (_set(("statistics", "primary"), "three_way_core"), "statistics.primary"),
        (_set(("statistics", "secondary"), "pairwise"), "statistics.secondary"),
        (_set(("reference_sets",), ["all_tested_triples"]), "reference_sets"),
        (_set(("lead_triple",), [PD, HIV, "Asthma"]), "lead_triple"),
        (_set(("lead_trait",), "METPO:1000623"), "lead_trait"),
        (_set(("amendments",), []), "strain_rollup is None"),
        (_set(("amendments", 0, "sets", "strain_rollup", "primary"), "phenotype"), "strain_rollup"),
        (_set(("amendments", 1, "sets", "secondary_rank", "by"), "u6_p"), "secondary_rank"),
        (_set(("amendments", 0, "sets", "seed"), 1), "unknown fields"),
        (_set(("amendments", 1, "date"), "2026-10-21"), "date order"),
        (lambda data: data["amendments"][0].pop("reason"), "needs exactly"),
    ],
)
def test_protocol_is_refused_when_it_says_what_the_code_does_not_do(
    tmp_path: Path, edit: Callable[[dict[str, Any]], None], message: str
) -> None:
    """Every field the code relies on is checked, so an edit cannot pass as followed (#872)."""
    with pytest.raises(ValueError, match=message):
        load_protocol(_edited(tmp_path, edit))


@pytest.mark.parametrize(
    ("parameters", "expected"),
    [
        ({}, ("off", PRIMARY)),
        ({"strain_rollup": "off"}, ("off", PRIMARY)),
        ({"strain_rollup": "phenotype"}, ("phenotype", SENSITIVITY)),
        ({"strain_rollup": "phenotype", "source_mode": "pinned"}, ("phenotype", SENSITIVITY)),
    ],
)
def test_rollup_scope_reads_the_declared_scopes(
    parameters: dict[str, str], expected: tuple[str, str]
) -> None:
    """A core without the parameter predates #828 and is off (#872)."""
    assert rollup_scope(parameters, load_protocol()) == expected


def test_rollup_scope_refuses_an_undeclared_scope() -> None:
    with pytest.raises(ValueError, match="neither the declared primary"):
        rollup_scope({"strain_rollup": "full"}, load_protocol())


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
        ("even_votes", 15, 15, 0, 1.0, TESTED),
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
    # Without an always-included trait, the screen alone chooses (#854).
    assert TRAIT_ID not in set(choose_traits(screen, 3, None)["trait_id"])
    # As many studies down as up is not decreasing, however many traits are taken (#872).
    assert "even_votes" not in set(choose_traits(screen, 10, None)["trait_id"])


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
    the first fifteen and facultative anaerobiosis on the last twenty. Increased
    signatures report cellulose-only species and decreased ones mostly facultative
    species, so the screen selects facultative anaerobiosis (every study shifts down)
    while cellulose, shifting up, joins only as the always-included trait.
    """
    rng = np.random.default_rng(835)
    rows = []
    number = 0
    for condition, count in STUDIES.items():
        for _ in range(count):
            number += 1
            for signature, direction in ((1, "increased"), (2, "decreased")):
                if direction == "increased":
                    picked = rng.choice(np.arange(2001, 2011), 4, replace=False)
                else:
                    picked = np.concatenate(
                        [
                            rng.choice(np.arange(2016, 2031), 3, replace=False),
                            rng.choice(np.arange(2001, 2016), 1),
                        ]
                    )
                species = sorted(picked.tolist())
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
    # The screen selects facultative anaerobiosis; cellulose joins as always included.
    in_top = dict(zip(traits["trait_id"], traits["in_top_n"].astype(bool), strict=True))
    assert in_top == {OTHER_TRAIT: True, TRAIT_ID: False}
    assert len(triples) == 10 * len(traits)
    assert triples["is_lead_triple"].sum() == len(traits)
    assert set(triples["u6_draws"]) <= {0, 400}

    summary = json.loads((out / OUTPUTS["summary"]).read_text(encoding="utf-8"))
    assert summary["follows_protocol"] is False
    # The fixture core predates #828, so it records no scope and reads as off (#872).
    assert summary["strain_rollup"]["scope"] == "off"
    assert summary["strain_rollup"]["role"] == PRIMARY
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
    assert report["strain_rollup"]["scope"] == "off"
    assert report["eligible_conditions"] == 5
    assert report["triples"] == 10
    assert report["permutations"] == 100_000
    assert report["estimated_null_minutes"] > 0
    assert not out.exists()


# --------------------------------------------------------------------------- study splits (#854)


def _inputs(records: Sequence[Signature], matrix: pd.DataFrame) -> ScanInputs:
    """Scan inputs from signatures; each screen profile counts its taxa bearing each trait."""
    profiles, taxa = _frames(records)
    columns = tuple(c for c in matrix.columns if c not in ("taxon_name", "taxonomic_rank"))
    bearing = {trait: set(matrix.index[matrix[trait] > 0]) for trait in columns}
    members = taxa.groupby("profile_id")["taxon_id"].agg(set).to_dict()
    screen = []
    for profile_id, row in profiles.iterrows():
        found = members.get(str(profile_id), set())
        counts = {trait: len(found & bearing[trait]) for trait in columns}
        screen.append(
            Profile(
                str(row["study_id"]),
                str(row["experiment_id"]),
                str(row["dir"]),
                str(row["condition"]),
                str(row["body_site"]),
                counts,
                sum(counts.values()),
            )
        )
    return ScanInputs(profiles, taxa, matrix, tuple(screen), columns, {})


def _split_corpus() -> tuple[ScanInputs, list[str], list[str]]:
    """Eight conditions of eight studies, split by alternate studies into two halves.

    Taxon t40 bears every trait and no signature reports it; ``_tamper`` uses it.
    """
    labels = (PD, HIV, CRC, "Asthma", "Autism", "Obesity", "Psoriasis", "Gout")
    records = _random_corpus(31, dict.fromkeys(labels, 8), n_taxa=40, size=6)
    traits = {
        TRAIT_ID: [*range(0, 40, 2), 40],
        OTHER_TRAIT: [*range(1, 40, 3), 40],
        "METPO:1000623": [*range(0, 40, 5), 40],
    }
    inputs = _inputs(records, _matrix(41, traits))
    ordered = sorted(inputs.studies, key=lambda study: int(study.split(":")[1]))
    return inputs, ordered[::2], ordered[1::2]


def _tamper(inputs: ScanInputs, studies: Sequence[str]) -> ScanInputs:
    """The same population with other signature taxa and screen counts for ``studies``:
    every taxon they report becomes t40, which no other study reports."""
    taxa = inputs.taxa.copy()
    taxa.loc[taxa["study_id"].isin(studies), "taxon_id"] = "t40"
    screen = tuple(
        dataclasses.replace(
            profile,
            guild_counts={trait: count + 3 for trait, count in profile.guild_counts.items()},
            total=profile.total + 3 * len(profile.guild_counts),
        )
        if profile.study_id in studies
        else profile
        for profile in inputs.screen
    )
    return dataclasses.replace(inputs, taxa=taxa, screen=screen)


KEYS = ["condition_a", "condition_b", "condition_c", "trait_id"]


def _keys(frame: pd.DataFrame) -> list[tuple[str, str, str, str]]:
    return list(frame[KEYS].itertuples(index=False, name=None))


def test_select_reads_only_the_training_studies() -> None:
    """Eligibility, the trait screen, the sets and the U6 pool come from the training
    studies alone: changing every held-out study leaves the selection unchanged."""
    inputs, train, test = _split_corpus()
    seed = (20261022, 0, 0)
    chosen = select(inputs, train, permutations=300, seed=seed)
    assert chosen.studies == frozenset(train)
    assert chosen.conditions["stool_studies"].max() == 4
    again = select(_tamper(inputs, test), train, permutations=300, seed=seed)
    pd.testing.assert_frame_equal(chosen.triples, again.triples)
    pd.testing.assert_frame_equal(chosen.traits, again.traits)
    pd.testing.assert_frame_equal(chosen.conditions, again.conditions)
    assert np.array_equal(chosen.result.u3_taxa, again.result.u3_taxa)
    changed = select(_tamper(inputs, train), train, permutations=300, seed=seed)
    assert not changed.triples.equals(chosen.triples)


def test_evaluate_reads_only_the_held_out_studies() -> None:
    """Hypotheses chosen on the training half are tested on the held-out half alone,
    in the given order, and changing every training study leaves the tests unchanged."""
    inputs, train, test = _split_corpus()
    chosen = rank_candidates(select(inputs, train, permutations=300, seed=(1, 0, 0)).triples)
    hypotheses = _keys(chosen.head(4))
    tested = evaluate(inputs, test, hypotheses, permutations=300, seed=(1, 0, 1))
    assert tested["hypothesis"].tolist() == [0, 1, 2, 3]
    assert _keys(tested) == hypotheses
    assert set(tested["test_status"]) <= {TESTED, NO_TAXA}
    again = evaluate(_tamper(inputs, train), test, hypotheses, permutations=300, seed=(1, 0, 1))
    pd.testing.assert_frame_equal(tested, again)
    changed = evaluate(_tamper(inputs, test), test, hypotheses, permutations=300, seed=(1, 0, 1))
    assert not changed.equals(tested)


def test_evaluate_on_the_scan_layout_reproduces_the_scan() -> None:
    """Listed triples are scored one at a time with the same draws and arithmetic as the
    scan's blocks: on every study, with every eligible condition involved, the held-out
    test of a scanned triple gives the scan's values."""
    inputs, _, _ = _split_corpus()
    full = select(inputs, permutations=400, seed=7)
    rows = full.triples.iloc[::2]
    involved = set(rows[["condition_a", "condition_b", "condition_c"]].to_numpy().ravel())
    assert involved == set(full.result.conditions)
    tested = evaluate(inputs, inputs.studies, _keys(rows), permutations=400, seed=7)
    ranks = ["rank_mean_distance", "rank_core", "z_rank_mean_distance", "z_rank_core"]
    expected = rows.drop(columns=ranks).reset_index(drop=True)
    for column in ("u6_p_mean_distance", "u6_p_core", "u6_draws_source", "test_status"):
        assert tested[column].tolist() == expected[column].tolist(), column
    pd.testing.assert_frame_equal(
        tested.drop(columns=["hypothesis", *OUTCOME_COLUMNS]),
        expected,
        check_dtype=False,
        check_exact=False,
    )


def test_evaluate_keeps_a_hypothesis_it_cannot_test() -> None:
    """A held-out half without a study of one condition does not replace the hypothesis:
    it is returned as untestable, beside the testable ones."""
    inputs, _, test = _split_corpus()
    gout = set(inputs.profiles.loc[inputs.profiles["condition"].eq("Gout"), "study_id"])
    held_out = [study for study in test if study not in gout]
    hypotheses = [(PD, HIV, "Gout", TRAIT_ID), (PD, HIV, CRC, TRAIT_ID)]
    tested = evaluate(inputs, held_out, hypotheses, permutations=200, seed=3)
    assert tested["test_status"].tolist() == [NO_UNITS, TESTED]
    first = tested.iloc[0]
    assert (first["condition_a"], first["condition_b"], first["condition_c"]) == ("Gout", HIV, PD)
    assert [first["u6_units_a"], first["u6_units_b"] > 0] == [0, True]
    assert pd.isna(first["u6_p_mean_distance"])
    assert tested["is_lead_triple"].tolist() == [False, True]


def test_split_inputs_are_checked() -> None:
    inputs, _, test = _split_corpus()
    with pytest.raises(ValueError, match="three conditions"):
        evaluate(inputs, test, [(PD, PD, HIV, TRAIT_ID)], permutations=10, seed=1)
    with pytest.raises(ValueError, match="not a screen candidate"):
        evaluate(inputs, test, [(PD, HIV, CRC, "METPO:0")], permutations=10, seed=1)
    with pytest.raises(ValueError, match="not in the population"):
        inputs.restricted(["bsdb:999999"])


def test_split_seeds_are_reproducible() -> None:
    """A (seed, split, stage) seed gives the same draws every time and other draws for
    another split or stage."""
    inputs, train, test = _split_corpus()
    first = select(inputs, train, permutations=300, seed=(20261022, 4, 0))
    assert first.triples.equals(
        select(inputs, train, permutations=300, seed=(20261022, 4, 0)).triples
    )
    other = select(inputs, train, permutations=300, seed=(20261022, 5, 0))
    assert not first.triples["u6_p_mean_distance"].equals(other.triples["u6_p_mean_distance"])
    hypotheses = _keys(rank_candidates(first.triples).head(2))
    one = evaluate(inputs, test, hypotheses, permutations=300, seed=(20261022, 4, 1))
    assert one.equals(evaluate(inputs, test, hypotheses, permutations=300, seed=(20261022, 4, 1)))


def test_select_can_leave_the_trait_set_to_the_screen() -> None:
    inputs, train, _ = _split_corpus()
    declared = select(inputs, train, permutations=50, seed=1)
    screened = select(inputs, train, permutations=50, seed=1, always_included=False)
    in_top = declared.traits.loc[declared.traits["in_top_n"], "trait_id"]
    assert screened.traits["trait_id"].tolist() == in_top.tolist()


def test_select_without_traits_returns_an_empty_table() -> None:
    """A training half whose screen finds no decreasing trait selects nothing."""
    inputs, train, _ = _split_corpus()
    flat = tuple(
        dataclasses.replace(profile, guild_counts=dict.fromkeys(profile.guild_counts, 1), total=3)
        for profile in inputs.screen
    )
    empty = select(dataclasses.replace(inputs, screen=flat), train, always_included=False)
    assert empty.traits.empty
    assert empty.triples.empty
    assert list(empty.triples.columns) == list(TRIPLE_COLUMNS)


def test_triple_table_columns_are_the_declared_list() -> None:
    inputs, _, _ = _split_corpus()
    assert list(select(inputs, permutations=20, seed=1).triples.columns) == list(TRIPLE_COLUMNS)


def test_rank_candidates_orders_by_p_then_z() -> None:
    rows = pd.DataFrame(
        {
            "trait_id": ["a", "b", "c", "d", "e"],
            "condition_a": "X",
            "condition_b": "Y",
            "condition_c": "Z",
            "u6_p_mean_distance": [0.2, 0.01, 0.01, 0.5, 0.001],
            "u6_z_mean_distance": [0.0, -1.0, -3.0, 1.0, -9.0],
            "u6_p_core": [0.3, 0.02, 0.02, 0.3, 0.001],
            "u6_z_core": [0.0, 1.0, 3.0, 2.0, 9.0],
            "test_status": [TESTED, TESTED, TESTED, TESTED, NO_TAXA],
        }
    )
    assert rank_candidates(rows)["trait_id"].tolist() == ["c", "b", "a", "d"]
    ranked = rank_candidates(rows, CORE)
    assert ranked["trait_id"].tolist() == ["c", "b", "d", "a"]
    assert ranked["candidate_rank"].tolist() == [1, 2, 3, 4]


def test_load_core_inputs_reads_a_tier1_2_core_only(cores: dict[str, Path]) -> None:
    _, inputs = load_core_inputs(cores["tier1_2"])
    assert len(inputs.studies) == sum(STUDIES.values())
    assert TRAIT_ID in inputs.columns
    with pytest.raises(ValueError, match=TIER1_2):
        load_core_inputs(cores["all_paired"])


# --------------------------------------------------------------------------- review (#872)


def _wrapped_layout() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Five conditions with 4 units each in 15 studies, so the blocks wrap and overlap."""
    labels = ["A", "B", "C", "D", "E"]
    records = _random_corpus(7, dict.fromkeys(labels, 3), n_taxa=30, size=5)
    nexts = dict(zip(labels, labels[1:] + labels[:1], strict=True))
    records += [
        (study, f"{study}-2", nexts[condition], site, dec, inc)
        for study, _, condition, site, dec, inc in records[::3]
    ]
    profiles, taxa = _frames(records)
    return profiles, taxa, labels


def test_every_triple_is_drawn_on_its_own_exactly_when_its_blocks_overlap() -> None:
    """On a wrapped layout, a triple gets its own draws exactly when any pair of its
    blocks overlaps, the a-c pair included (#873)."""
    profiles, taxa, labels = _wrapped_layout()
    prepared = prepare(profiles, taxa, _matrix(30, {TRAIT_ID: range(30)}), [TRAIT_ID], labels)
    n_studies = int(prepared.unit_study.max()) + 1
    blocks = [set(block.tolist()) for block in prepared.positions]
    expected = []
    for a, b, c in prepared.triples:
        if prepared.units[[a, b, c]].sum() > n_studies:
            expected.append(NOT_DRAWN)
        elif blocks[a] & blocks[b] or blocks[a] & blocks[c] or blocks[b] & blocks[c]:
            expected.append(OWN)
        else:
            expected.append(SHARED)
    assert prepared.source.tolist() == expected
    assert {OWN, SHARED} <= set(expected)
    # (A, B, D) overlaps only in its a-c pair, A with D.
    index = [tuple(row) for row in prepared.triples.tolist()].index((0, 1, 3))
    assert not blocks[0] & blocks[1]
    assert not blocks[1] & blocks[3]
    assert blocks[0] & blocks[3]
    assert prepared.source[index] == OWN


def test_a_triple_that_fills_the_pool_is_drawn() -> None:
    """Units adding up to exactly the pool's studies still fit without overlap (#872)."""
    records = _random_corpus(2, dict.fromkeys(("A", "B", "C"), 3), n_taxa=12, size=3)
    profiles, taxa = _frames(records)
    matrix = _matrix(12, {TRAIT_ID: range(12)})
    result = scan(
        profiles, taxa, matrix, [TRAIT_ID], ["A", "B", "C"], permutations=50, seed=1, tolerance=0
    )
    table = triple_table(result, LEAD)
    assert result.pool_studies == int(result.units.sum()) == 9
    assert table["u6_draws_source"].tolist() == [SHARED]
    assert table["test_status"].tolist() == [TESTED]


def test_a_lead_that_shares_a_study_is_placed_among_the_disjoint_triples() -> None:
    """When PD and HIV share a study, the study-disjoint reference leaves the lead out,
    and its rank there counts it as one more triple (#872)."""
    studies = {PD: 4, HIV: 3, CRC: 3, "Asthma": 3, "Obesity": 3}
    records = _random_corpus(17, studies, n_taxa=30, size=5)
    study, _, _, site, dec, inc = records[0]
    records.append((study, f"{study}-2", HIV, site, dec, inc))
    profiles, taxa = _frames(records)
    matrix = _matrix(30, {TRAIT_ID: range(0, 30, 2)})
    result = scan(
        profiles, taxa, matrix, [TRAIT_ID], sorted(studies), permutations=300, seed=2, tolerance=0
    )
    triples = triple_table(result, LEAD)
    lead_row = triples[triples["is_lead_triple"]].iloc[0]
    assert lead_row["shared_studies"] == 1
    row = lead_and_distribution(triples, [TRAIT_ID], 20)[0].iloc[0]
    disjoint = triples[triples["test_status"].eq(TESTED) & triples["shared_studies"].eq(0)]
    p = lead_row["u6_p_mean_distance"]
    assert row[f"{MEAN_DISTANCE}_{DISJOINT}_triples"] == len(disjoint) + 1
    at_most = int((disjoint["u6_p_mean_distance"] <= p).sum()) + 1
    assert row[f"{MEAN_DISTANCE}_{DISJOINT}_rank_max"] == at_most


def test_z_ranks_put_more_overlap_first() -> None:
    """The amendment's secondary rank: ascending z for the mean distance and descending z
    for the core, among the tested triples whose z is defined (#872)."""
    profiles, taxa, matrix = _lead_corpus()
    conditions = sorted(condition_table(profiles, 3).query("eligible")["condition"])
    result = scan(
        profiles, taxa, matrix, [TRAIT_ID], conditions, permutations=500, seed=3, tolerance=0
    )
    triples = triple_table(result, LEAD)
    tested = triples[triples["test_status"].eq(TESTED)]
    for column, rank, ascending in (
        ("u6_z_mean_distance", "z_rank_mean_distance", True),
        ("u6_z_core", "z_rank_core", False),
    ):
        defined = tested[tested[column].notna()]
        expected = defined[column].rank(method="min", ascending=ascending).astype(int)
        assert defined[rank].astype(int).tolist() == expected.tolist()
        assert tested.loc[tested[column].isna(), rank].isna().all()
    row = lead_and_distribution(triples, [TRAIT_ID], 20)[0].iloc[0]
    z = tested["u6_z_mean_distance"].dropna()
    lead_z = tested.loc[tested["is_lead_triple"], "u6_z_mean_distance"].iloc[0]
    assert row[f"{MEAN_DISTANCE}_{ALL_TRIPLES}_z_rank_min"] == int((z < lead_z).sum()) + 1
    core_z = tested["u6_z_core"].dropna()
    lead_core = tested.loc[tested["is_lead_triple"], "u6_z_core"].iloc[0]
    assert row[f"{CORE}_{ALL_TRIPLES}_z_rank_min"] == int((core_z > lead_core).sum()) + 1


def test_evaluate_flags_explicit_outcomes() -> None:
    """A missing condition, shared studies, an empty set and an undefined null are each
    flagged on the held-out side (#854)."""
    labels = (PD, HIV, CRC, "Asthma")
    records = _random_corpus(41, dict.fromkeys(labels, 4), n_taxa=20, size=4)
    study, _, _, site, dec, inc = records[0]
    records.append((study, f"{study}-2", HIV, site, dec, inc))
    sparse = "METPO:1000623"
    records = [
        (s, e, c, site, [*dec, "t30"] if c == CRC else dec, inc)
        for s, e, c, site, dec, inc in records
    ]
    inputs = _inputs(records, _matrix(31, {TRAIT_ID: range(0, 20, 2), sparse: [30]}))
    hypotheses = [
        (PD, HIV, CRC, TRAIT_ID),
        (PD, "Asthma", CRC, TRAIT_ID),
        (PD, HIV, "Gout", TRAIT_ID),
        (PD, "Asthma", CRC, sparse),
    ]
    tested = evaluate(inputs, inputs.studies, hypotheses, permutations=200, seed=5)
    assert tested["test_status"].tolist() == [TESTED, TESTED, NO_UNITS, TESTED]
    assert tested["missing_condition"].tolist() == [False, False, True, False]
    assert tested["shares_studies"].tolist() == [True, False, True, False]
    assert tested["empty_set"].tolist() == [False, False, True, True]
    assert tested["null_undefined"].tolist() == [False, False, True, False]
    assert tested["shared_studies"].tolist() == [1, 0, 1, 0]
    assert [tested.loc[3, f"n_{x}"] for x in "abc"] == [0, 1, 0]


def test_closest_disjoint_triples_is_the_854_rule() -> None:
    """Per trait, the study-disjoint tested triple with the least observed distance, ties
    within the tolerance going to the first conditions in lexical order."""
    rows = pd.DataFrame(
        {
            "trait_id": ["x", "x", "x", "x", "y", "y", "z"],
            "condition_a": ["B", "A", "A", "A", "A", "A", "A"],
            "condition_b": ["C", "C", "B", "B", "B", "C", "B"],
            "condition_c": ["D", "D", "C", "D", "C", "D", "C"],
            MEAN_DISTANCE: [0.3, 0.3 + 1e-13, 0.1, 0.31, 0.5, 0.5, 0.2],
            "shared_studies": [0, 0, 2, 0, 0, 0, 1],
            "test_status": [TESTED] * 6 + [TESTED],
        }
    )
    chosen = closest_disjoint_triples(rows)
    assert _keys(chosen) == [("A", "C", "D", "x"), ("A", "B", "C", "y")]


def test_select_without_draws_gives_observed_statistics_only() -> None:
    """``permutations=0`` skips U6: the observed statistics are those of a drawn run and
    the U6 columns are empty, which is all #854's selection rule needs."""
    inputs, train, _ = _split_corpus()
    observed = select(inputs, train, permutations=0, seed=1, always_included=False)
    drawn = select(inputs, train, permutations=200, seed=1, always_included=False)
    same = [*KEYS, "three_way_core", "three_way_union", MEAN_DISTANCE, "shared_studies"]
    pd.testing.assert_frame_equal(observed.triples[same], drawn.triples[same])
    tested = observed.triples[observed.triples["test_status"].eq(TESTED)]
    assert tested["u6_p_mean_distance"].isna().all()
    assert tested["u6_draws"].eq(0).all()
    picked = closest_disjoint_triples(observed.triples)
    assert _keys(picked) == _keys(closest_disjoint_triples(drawn.triples))
    assert len(picked) == len(observed.traits)


def _core_tables(records: Sequence[Signature]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The core's experiments, profile keys and signature taxa for these signatures."""
    experiments = pd.DataFrame(
        [
            {
                "study_id": s,
                "experiment_id": e,
                "condition": c,
                "body_site": b,
                "host_species": "Homo sapiens",
            }
            for s, e, c, b, _, _ in records
        ]
    )
    keyed = pd.DataFrame(
        [
            {"study_id": s, "experiment_id": e, "inc_dec": d}
            for s, e, *_ in records
            for d in ("inc", "dec")
        ]
    )
    reported = pd.DataFrame(
        [
            {"study_id": s, "experiment_id": e, "direction": direction, "taxon_id": taxon}
            for s, e, _, _, dec, inc in records
            for direction, members in (("decreased", dec), ("increased", inc))
            for taxon in members
        ]
    )
    return experiments, keyed, reported


def test_build_population_takes_in_memory_tables_for_some_studies() -> None:
    """The population builds from in-memory tables and keeps whole studies when given a
    subset, as if the other studies were never there (#872, #854)."""
    records = _random_corpus(5, {PD: 4, HIV: 4, CRC: 4}, n_taxa=20, size=4)
    experiments, keyed, reported = _core_tables(records)
    matrix_taxa = pd.Index([f"t{i}" for i in range(20)])
    keep = sorted({record[0] for record in records})[::2]
    profiles, taxa = build_population(experiments, keyed, reported, matrix_taxa, keep)
    whole_profiles, whole_taxa = build_population(experiments, keyed, reported, matrix_taxa)
    assert set(profiles["study_id"]) == set(keep)
    pd.testing.assert_frame_equal(profiles, whole_profiles[whole_profiles["study_id"].isin(keep)])
    kept = whole_taxa[whole_taxa["study_id"].isin(keep)].reset_index(drop=True)
    pd.testing.assert_frame_equal(taxa.reset_index(drop=True), kept)
    unpaired = keyed[~(keyed["experiment_id"].eq(records[0][1]) & keyed["inc_dec"].eq("inc"))]
    with pytest.raises(ValueError, match="paired"):
        build_population(experiments, unpaired, reported, matrix_taxa)
