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
    CONSUMER,
    CONSUMER_RECORD,
    CORE,
    CORE_PRODUCTS,
    DISJOINT,
    MEAN_DISTANCE,
    NO_TAXA,
    NO_TRAITS,
    NO_UNITS,
    NOT_DRAWN,
    OUTCOME_COLUMNS,
    OUTPUTS,
    OWN,
    POOL_TOO_SMALL,
    PRIMARY,
    PROTOCOL_PATH,
    PROTOCOL_V1_PATH,
    SELECTED,
    SENSITIVITY,
    SHARED,
    STRAIN_ROLLUP,
    TIER1_2,
    TOO_FEW_CONDITIONS,
    TRIPLE_COLUMNS,
    VERSION,
    VERSION_V1,
    ScanInputs,
    bin_counts,
    block_positions,
    build_population,
    choose_traits,
    closest_disjoint_triples,
    condition_table,
    core_scope,
    evaluate,
    lead_and_distribution,
    lead_rank,
    load_core_inputs,
    load_inputs,
    load_protocol,
    main,
    plan,
    prepare,
    rank_candidates,
    recorded_protocol,
    scan,
    scope_record,
    scope_role,
    select,
    trait_table,
    triple_table,
)
from mtp.cellulose_decr.selection_scan_scopes import OUTPUTS as SCOPES_OUTPUTS
from mtp.cellulose_decr.selection_scan_scopes import ScanRun, load_run, pair_scopes, side_by_side
from mtp.cellulose_decr.selection_scan_scopes import main as scopes_main
from mtp.paths import file_sha256
from mtp.trait_guilds import Profile
from mtp.trait_screen import ALL_TIED, TESTED

PD, HIV, CRC = (CONDITIONS[short] for short in ("PD", "HIV", "CRC"))
LEAD = (PD, HIV, CRC)
# Tests write to temporary storage, which the scan and scopes refuse without the opt-in (#1087)
THROWAWAY = "--allow-temporary-output"
OTHER_TRAIT = "METPO:1000608"
# The v1 protocol file's hash as v2 pins it, which is the hash the records of 2026-10-23
# give; the file must stay byte-identical so that they still verify against it (#1050).
V1_SHA256 = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))["amends"]["sha256"]
V1_PATH = "methods/cellulose_decr_selection_scan_protocol.json"

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
    amendments of 2026-10-22 (#872) and the v2 amendment of 2026-10-25 (#1050) applied."""
    protocol = load_protocol()
    assert protocol.version == VERSION == "cellulose-decr-selection-scan-v2"
    assert protocol.min_studies == 3
    assert protocol.top_n == 10
    assert protocol.permutations == 100_000
    assert protocol.seed == 20261022
    assert protocol.tolerance == 1e-12
    assert protocol.bins == 20
    assert protocol.lead_triple == LEAD
    assert protocol.lead_trait == TRAIT_ID == protocol.always_included
    assert (protocol.scope_primary, protocol.scope_sensitivity) == ("phenotype", ("off",))
    assert protocol.scope_absent == "off"
    # Declared after the off and phenotype results had been seen: a post-results choice of
    # primary scope, so both scopes are always reported.
    assert protocol.results_seen == ("off", "phenotype")
    assert protocol.side_by_side


def test_v2_amends_v1_after_both_scopes_results_and_leaves_v1_as_recorded() -> None:
    """v2 is one amendment on top of v1 (#1050). The v1 file stays byte-identical, so the
    records written under it on 2026-10-23 still verify against it, and it still reads as
    it declared: off primary, phenotype the sensitivity, no results seen."""
    assert file_sha256(PROTOCOL_V1_PATH) == V1_SHA256
    v1 = load_protocol(PROTOCOL_V1_PATH)
    assert v1.version == VERSION_V1
    assert (v1.scope_primary, v1.scope_sensitivity, v1.scope_absent) == (
        "off",
        ("phenotype",),
        "off",
    )
    assert (v1.results_seen, v1.side_by_side) == ((), False)
    # Only the scope declaration differs.
    v2 = load_protocol()
    scope_fields = ("version", "scope_primary", "scope_sensitivity", "results_seen", "side_by_side")
    assert dataclasses.replace(v2, **{name: getattr(v1, name) for name in scope_fields}) == v1

    data = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    assert data["amends"] == {"protocol": VERSION_V1, "path": V1_PATH, "sha256": V1_SHA256}
    [amendment] = data["amendments"]
    assert amendment["date"] == "2026-10-25"
    assert set(amendment["sets"]) == {"strain_rollup"}
    known = amendment["known_when_amended"]
    assert known.startswith("Made after results.")
    # Both scopes' results had been seen: off in the rerun record, phenotype in a run of the
    # same recipe on 2026-10-24 that the repository does not record.
    assert "docs/reviews/2026-10-23-final-rerun-record.md" in known
    for seen in (
        "2026-10-24",
        "Tier 1+2 x phenotype",
        "outside the repository's recorded evidence",
    ):
        assert seen in known
    # Symmetric with off: the phenotype core, by its run-ID prefix, and its lead row.
    phenotype = known[known.index("Phenotype:") :]
    for seen in ("run-ec264369", "1,368 of 5,984", "p = 0.3515", "3,446", "p = 0.822"):
        assert seen in phenotype
    assert "no scan result on a phenotype core" not in known
    assert "post-results choice of primary scope" in amendment["change"]
    for issue in ("#788", "#1050"):
        assert issue in amendment["reason"]
    # The original declaration stays as committed; the amendments carry the changes.
    original = json.loads(PROTOCOL_V1_PATH.read_text(encoding="utf-8"))
    assert "strain_rollup" not in original
    assert "secondary_rank" not in original
    assert [amendment["date"] for amendment in original["amendments"]] == ["2026-10-22"] * 2


def _edited(
    tmp_path: Path, edit: Callable[[dict[str, Any]], None], source: Path = PROTOCOL_V1_PATH
) -> Path:
    data = json.loads(source.read_text(encoding="utf-8"))
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
        (_set(("protocol",), "cellulose-decr-selection-scan-v0"), "this code reads"),
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
        (_set(("amends",), {"protocol": VERSION}), "amends no earlier version"),
    ],
)
def test_protocol_is_refused_when_it_says_what_the_code_does_not_do(
    tmp_path: Path, edit: Callable[[dict[str, Any]], None], message: str
) -> None:
    """Every field of v1 the code relies on is checked, so an edit cannot pass as followed
    (#872)."""
    with pytest.raises(ValueError, match=message):
        load_protocol(_edited(tmp_path, edit))


@pytest.mark.parametrize(
    ("edit", "message"),
    [
        (_set(("amends", "sha256"), "0" * 64), "changed since"),
        (_set(("amends", "protocol"), "cellulose-decr-selection-scan-v0"), "must amend"),
        (_set(("amends", "path"), "methods/cellulose_decr_selection_scan.md"), "must amend"),
        (lambda data: data.pop("amends"), "must amend"),
        (_set(("amendments",), []), "strain_rollup is"),
        (_set(("amendments", 0, "sets", "strain_rollup", "primary"), "off"), "strain_rollup"),
        (_set(("amendments", 0, "sets", "strain_rollup", "reported"), "alone"), "strain_rollup"),
        (
            _set(("amendments", 0, "sets", "strain_rollup", "declared_after_results"), []),
            "strain_rollup",
        ),
        (
            _set(("amendments", 0, "sets", "strain_rollup", "declared_after_results"), ["off"]),
            "strain_rollup",
        ),
        (_set(("amendments", 0, "sets", "seed"), 1), "unknown fields"),
        (_set(("null",), {"permutations": 10}), "only by amendment"),
        (_set(("amendments", 0, "date"), "2026-10-21"), "date order"),
        (lambda data: data["amendments"][0].pop("known_when_amended"), "needs exactly"),
    ],
)
def test_v2_is_refused_unless_it_amends_v1_as_declared(
    tmp_path: Path, edit: Callable[[dict[str, Any]], None], message: str
) -> None:
    """v2 must name v1 by its committed path and hash, come after v1's amendments, and set
    the scope declaration this code implements: phenotype primary, off the sensitivity,
    reported side by side, declared after the off and phenotype results (#1050)."""
    with pytest.raises(ValueError, match=message):
        load_protocol(_edited(tmp_path, edit, PROTOCOL_PATH))


def test_scope_roles_follow_the_declaring_version() -> None:
    """v2: phenotype is primary and off the declared sensitivity (#1050). v1 declared the
    reverse (#872). Both refuse full."""
    for protocol, primary, sensitivity in (
        (load_protocol(), "phenotype", "off"),
        (load_protocol(PROTOCOL_V1_PATH), "off", "phenotype"),
    ):
        assert scope_role(primary, protocol) == PRIMARY
        assert scope_role(sensitivity, protocol) == SENSITIVITY
        with pytest.raises(ValueError, match="neither the declared primary"):
            scope_role("full", protocol)


def test_a_scope_names_the_scope_it_must_be_reported_beside() -> None:
    """Each scope's record names the other declared scope, so neither is read alone, and
    says that its role was declared after both scopes' results (#1050)."""
    protocol = load_protocol()
    phenotype = scope_record("phenotype", "recorded", protocol)
    assert phenotype == {
        "scope": "phenotype",
        "scope_source": "recorded",
        "role": PRIMARY,
        "declared_primary": "phenotype",
        "declared_sensitivity": ["off"],
        "declared_by": VERSION,
        "declared_after_results": ["off", "phenotype"],
        "reported": "side_by_side",
        "report_beside": ["off"],
    }
    off = scope_record("off", "absent_means", protocol)
    assert (off["role"], off["scope_source"], off["report_beside"]) == (
        SENSITIVITY,
        "absent_means",
        ["phenotype"],
    )


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
    """Tier 1+2 cores of each strain-rollup scope, and an all-paired core, of the same
    synthetic inputs. ``tier1_2_phenotype`` is the declared primary scope and
    ``tier1_2_off`` the declared sensitivity (#1050)."""
    from mtp.paths import REPO_ROOT
    from mtp.workflow.operations import OPERATIONS, validate_core_product
    from mtp.workflow.runtime import execute_workflow

    root = tmp_path_factory.mktemp("selection-scan-core")
    inputs = _write_core_inputs(root)
    runs = {}
    for name, spec, parameters in (
        (
            "tier1_2_off",
            "bugsigdb_core_subset.yaml",
            {"contrast_population": TIER1_2, "strain_rollup": "off"},
        ),
        (
            "tier1_2_phenotype",
            "bugsigdb_core_subset.yaml",
            {"contrast_population": TIER1_2, "strain_rollup": "phenotype"},
        ),
        (
            "tier1_2_full",
            "bugsigdb_core_subset.yaml",
            {"contrast_population": TIER1_2, "strain_rollup": "full"},
        ),
        ("all_paired", "bugsigdb_core.yaml", {"strain_rollup": "off"}),
    ):
        run = execute_workflow(
            REPO_ROOT / "workflows" / spec,
            inputs,
            root / name,
            OPERATIONS,
            parameter_overrides={"source_mode": "fixture", **parameters},
            validator=validate_core_product,
            repository_root=REPO_ROOT,
        )
        runs[name] = Path(run.manifest_path)
    return runs


def test_main_writes_every_output_bound_to_the_core(cores: dict[str, Path], tmp_path: Path) -> None:
    """A run on a Tier 1+2 core writes the declared outputs, a provenance naming its inputs,
    rule, protocol version and scope, and a consumer record that verifies against the core."""
    from mtp.workflow.consumers import verify_consumer_record

    out = tmp_path / "scan"
    argv = ["--core-run", str(cores["tier1_2_phenotype"]), "--output-dir", str(out), THROWAWAY]
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
    assert summary["protocol"] == VERSION
    assert summary["follows_protocol"] is False
    # The core was built with the declared primary scope, phenotype since v2 (#1050), and
    # the summary names off as the scope it must be reported beside.
    protocol = load_protocol()
    assert summary["strain_rollup"] == scope_record("phenotype", "recorded", protocol)
    assert summary["strain_rollup"]["report_beside"] == ["off"]
    assert any("declared after the off and phenotype results" in note for note in summary["notes"])
    assert summary["counts"]["eligible_conditions"] == 5
    assert summary["lead_triple"]["eligible"] == dict.fromkeys(LEAD, True)
    provenance = json.loads((out / OUTPUTS["provenance"]).read_text(encoding="utf-8"))
    assert provenance["code_commit"] == "abc"
    assert provenance["core_run"]["contrast_population"] == TIER1_2
    assert provenance["core_run"]["parameters"]["strain_rollup"] == "phenotype"
    assert set(provenance["core_run"]["sources"]) == {"raw_dump", "kg_edges", "kg_nodes"}
    assert set(provenance["inputs"]) == set(CORE_PRODUCTS)
    assert provenance["protocol"] == {
        "path": "methods/cellulose_decr_selection_scan_protocol_v2.json",
        "sha256": file_sha256(PROTOCOL_PATH),
        "version": VERSION,
        "amends": {"protocol": VERSION_V1, "path": V1_PATH, "sha256": V1_SHA256},
    }
    assert provenance["strain_rollup"] == summary["strain_rollup"]
    assert (provenance["permutations"], provenance["seed"]) == (400, 3)
    assert recorded_protocol(summary, provenance) == protocol

    # The same core, rule and seed give the same tables.
    again = tmp_path / "again"
    assert (
        main(
            [
                "--core-run",
                str(cores["tier1_2_phenotype"]),
                "--output-dir",
                str(again),
                "--permutations",
                "400",
                "--seed",
                "3",
                THROWAWAY,
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
        main(
            ["--core-run", str(cores["all_paired"]), "--output-dir", str(tmp_path / "a"), THROWAWAY]
        )
    with pytest.raises(ValueError, match="neither the declared primary"):
        main(
            [
                "--core-run",
                str(cores["tier1_2_full"]),
                "--output-dir",
                str(tmp_path / "b"),
                THROWAWAY,
            ]
        )
    used = tmp_path / "used"
    used.mkdir()
    (used / "old.tsv").write_text("x\n", encoding="utf-8")
    with pytest.raises(ValueError, match="new or empty"):
        main(["--core-run", str(cores["tier1_2_off"]), "--output-dir", str(used), THROWAWAY])


def test_plan_reports_counts_without_writing(
    cores: dict[str, Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "plan"
    assert main(["--core-run", str(cores["tier1_2_off"]), "--output-dir", str(out), "--plan"]) == 0
    report = json.loads(capsys.readouterr().out)
    # off is the declared sensitivity since v2, to be reported beside phenotype (#1050).
    assert report["protocol"] == VERSION
    assert report["strain_rollup"] == scope_record("off", "recorded", load_protocol())
    assert report["strain_rollup"]["role"] == SENSITIVITY
    assert report["eligible_conditions"] == 5
    assert report["triples"] == 10
    assert report["permutations"] == 100_000
    assert report["estimated_null_minutes"] > 0
    assert not out.exists()


def test_a_core_without_a_scope_follows_the_declared_absent_means(cores: dict[str, Path]) -> None:
    """A core built before #828 records no strain_rollup. The scan reads it as the
    protocol's absent_means, off, which v2 makes the sensitivity, and records that the
    scope came from the declaration; no default elsewhere decides it (#927, #1050)."""
    from mtp.datamodel.workflow_pydantic import WorkflowRun

    run = WorkflowRun.model_validate_json(cores["tier1_2_off"].read_text(encoding="utf-8"))
    legacy = run.model_copy(
        update={
            "parameters": [item for item in run.parameters or [] if item.name != "strain_rollup"]
        }
    )
    protocol = load_protocol()
    assert core_scope(run, protocol) == ("off", "recorded")
    assert core_scope(legacy, protocol) == ("off", "absent_means")
    assert scope_record("off", "absent_means", protocol)["role"] == SENSITIVITY
    # The protocol's declaration is what is followed.
    declared = dataclasses.replace(protocol, scope_absent="phenotype")
    assert core_scope(legacy, declared) == ("phenotype", "absent_means")
    assert core_scope(run, declared) == ("off", "recorded")


# --------------------------------------------------------------------------- both scopes (#1050)


def _scan(core: Path, out: Path, seed: int = 3) -> Path:
    argv = ["--core-run", str(core), "--output-dir", str(out), "--code-commit", "abc", THROWAWAY]
    assert main([*argv, "--permutations", "200", "--seed", str(seed)]) == 0
    return out


@pytest.fixture(scope="module")
def scans(cores: dict[str, Path], tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """A scan of each declared scope, under the same rule and code. Tests leave them as
    written; a test that changes a scan makes its own."""
    root = tmp_path_factory.mktemp("selection-scan-scopes")
    return {
        "phenotype": _scan(cores["tier1_2_phenotype"], root / "phenotype"),
        "off": _scan(cores["tier1_2_off"], root / "off"),
    }


def _scopes(out: Path, *scans: Path) -> int:
    argv = [part for scan_dir in scans for part in ("--scan", str(scan_dir))]
    return scopes_main([*argv, "--output-dir", str(out), "--code-commit", "abc", THROWAWAY])


def _read_text_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def test_scopes_report_puts_phenotype_and_off_side_by_side(
    scans: dict[str, Path], tmp_path: Path
) -> None:
    """The report holds both scopes' trait tables as written, primary first whatever the
    order given, and names the amendment that declared the roles after both scopes' results."""
    out = tmp_path / "scopes"
    assert _scopes(out, scans["off"], scans["phenotype"]) == 0
    assert sorted(path.name for path in out.iterdir()) == sorted(SCOPES_OUTPUTS.values())

    table = _read_text_table(out / SCOPES_OUTPUTS["table"])
    assert list(table.columns[:2]) == ["strain_rollup", "role"]
    assert table["strain_rollup"].drop_duplicates().tolist() == ["phenotype", "off"]
    for scope, role in (("phenotype", PRIMARY), ("off", SENSITIVITY)):
        mine = table[table["strain_rollup"].eq(scope)]
        assert set(mine["role"]) == {role}
        own = _read_text_table(scans[scope] / OUTPUTS["traits"])
        pd.testing.assert_frame_equal(
            mine.drop(columns=["strain_rollup", "role"]).reset_index(drop=True), own
        )

    report = json.loads((out / SCOPES_OUTPUTS["summary"]).read_text(encoding="utf-8"))
    assert report["protocol"] == VERSION
    assert report["follows_protocol"] is False
    assert [(entry["scope"], entry["role"]) for entry in report["scopes"]] == [
        ("phenotype", PRIMARY),
        ("off", SENSITIVITY),
    ]
    declaration = report["strain_rollup"]
    assert (declaration["declared_primary"], declaration["declared_sensitivity"]) == (
        "phenotype",
        ["off"],
    )
    assert declaration["declared_after_results"] == ["off", "phenotype"]
    assert declaration["amendment"]["date"] == "2026-10-25"
    assert declaration["amendment"]["known_when_amended"].startswith("Made after results.")
    for entry in report["scopes"]:
        scan_dir = scans[entry["scope"]]
        for name in (OUTPUTS["summary"], OUTPUTS["provenance"], OUTPUTS["traits"]):
            assert entry["files"][name] == file_sha256(scan_dir / name)
        summary = json.loads((scan_dir / OUTPUTS["summary"]).read_text(encoding="utf-8"))
        assert entry["core_run_id"] == summary["core_run"]["run_id"]
        assert entry["lead_trait"]["trait_id"] == TRAIT_ID
    assert report["provenance"]["code_commit"] == "abc"


def test_scopes_report_never_reports_one_scope_alone(
    scans: dict[str, Path], tmp_path: Path
) -> None:
    """One scan, or two of the same scope, is refused, and nothing is written."""
    for given in ((scans["phenotype"],), (scans["off"], scans["off"])):
        out = tmp_path / "refused"
        with pytest.raises(ValueError, match="one scan of each declared scope"):
            _scopes(out, *given)
        assert not out.exists()
    used = tmp_path / "used"
    used.mkdir()
    (used / "old.tsv").write_text("x\n", encoding="utf-8")
    with pytest.raises(ValueError, match="new or empty"):
        _scopes(used, scans["phenotype"], scans["off"])


def test_scopes_report_refuses_scans_that_differ_in_more_than_scope(
    cores: dict[str, Path], scans: dict[str, Path], tmp_path: Path
) -> None:
    """The two scans must share the rule, code, protocol file and core sources, and every
    file read must be one the scan bound to its core."""
    reseeded = _scan(cores["tier1_2_off"], tmp_path / "off-seed-4", seed=4)
    with pytest.raises(ValueError, match="differ in rule"):
        _scopes(tmp_path / "a", scans["phenotype"], reseeded)

    tampered = _scan(cores["tier1_2_off"], tmp_path / "off-tampered")
    traits = tampered / OUTPUTS["traits"]
    traits.write_text(traits.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Consumer output changed"):
        _scopes(tmp_path / "b", scans["phenotype"], tampered)

    copied = tmp_path / "off-copied"
    copied.mkdir()
    for path in scans["off"].iterdir():
        (copied / path.name).write_bytes(path.read_bytes())
    with pytest.raises(ValueError, match="not bound"):
        load_run(copied)

    phenotype, off = load_run(scans["phenotype"]), load_run(scans["off"])
    assert [run.scope for run in pair_scopes([off, phenotype])] == ["phenotype", "off"]
    v1_run = dataclasses.replace(off, protocol=load_protocol(PROTOCOL_V1_PATH))
    with pytest.raises(ValueError, match=f"needs both under {VERSION}"):
        pair_scopes([phenotype, v1_run])
    for section, key, value in (
        ("core_run", "sources", {"raw_dump": "0" * 64}),
        ("core_run", "git_commit", "0" * 40),
        ("core_run", "parameters", {"strain_rollup": "off", "source_mode": "pinned"}),
        (None, "scripts", {}),
        (None, "code_commit", "def"),
    ):
        provenance = json.loads(json.dumps(off.provenance))
        (provenance[section] if section else provenance)[key] = value
        with pytest.raises(ValueError, match="differ in"):
            pair_scopes([phenotype, dataclasses.replace(off, provenance=provenance)])


def test_records_verify_against_the_protocol_version_they_name(scans: dict[str, Path]) -> None:
    """A v2 record verifies against v2. A record written under v1, as on 2026-10-23, still
    verifies against the unchanged v1 file (#1050). A record whose version, file hash or
    scope role disagrees with the protocol it names is refused."""
    summary = json.loads((scans["off"] / OUTPUTS["summary"]).read_text(encoding="utf-8"))
    provenance = json.loads((scans["off"] / OUTPUTS["provenance"]).read_text(encoding="utf-8"))
    assert recorded_protocol(summary, provenance).version == VERSION

    # The fields v1's code wrote for a run on an off core.
    v1_rollup = {
        "scope": "off",
        "role": PRIMARY,
        "declared_primary": "off",
        "declared_sensitivity": ["phenotype"],
    }
    v1_protocol = {"path": V1_PATH, "sha256": V1_SHA256, "version": VERSION_V1}
    v1_summary = {**summary, "protocol": VERSION_V1, "strain_rollup": v1_rollup}
    v1_provenance = {**provenance, "protocol": v1_protocol, "strain_rollup": v1_rollup}
    assert recorded_protocol(v1_summary, v1_provenance) == load_protocol(PROTOCOL_V1_PATH)

    sensitivity = {**v1_rollup, "role": SENSITIVITY}
    for record_summary, record_provenance, message in (
        (
            {**v1_summary, "strain_rollup": sensitivity},
            {**v1_provenance, "strain_rollup": sensitivity},
            "role",
        ),
        (v1_summary, {**v1_provenance, "strain_rollup": sensitivity}, "differ"),
        (v1_summary, {**v1_provenance, "protocol": {**v1_protocol, "sha256": "0" * 64}}, "hash"),
        (v1_summary, {**v1_provenance, "protocol": {**v1_protocol, "path": "x.json"}}, "path"),
        ({**v1_summary, "protocol": VERSION}, v1_provenance, "names"),
        (
            {**summary, "protocol": "cellulose-decr-selection-scan-v0"},
            {
                **provenance,
                "protocol": {
                    **provenance["protocol"],
                    "version": "cellulose-decr-selection-scan-v0",
                },
            },
            "unknown",
        ),
        (
            {**summary, "strain_rollup": {**summary["strain_rollup"], "role": PRIMARY}},
            {**provenance, "strain_rollup": {**summary["strain_rollup"], "role": PRIMARY}},
            "role",
        ),
        # The right role but v1's declaration.
        (
            {**summary, "strain_rollup": {**summary["strain_rollup"], "declared_primary": "off"}},
            {
                **provenance,
                "strain_rollup": {**summary["strain_rollup"], "declared_primary": "off"},
            },
            "declares",
        ),
        # v2's own file and hash, but another v1 hash in its amends block.
        (
            summary,
            {
                **provenance,
                "protocol": {
                    **provenance["protocol"],
                    "amends": {**provenance["protocol"]["amends"], "sha256": "0" * 64},
                },
            },
            "do not match",
        ),
    ):
        with pytest.raises(ValueError, match=message):
            recorded_protocol(record_summary, record_provenance)


def _rebind(
    scan_dir: Path,
    core: Path,
    edit: Callable[[dict[str, Any], dict[str, Any]], None] | None = None,
    consumer: str = CONSUMER,
) -> Path:
    """Edit a scan's summary and provenance, then bind them anew to ``core`` with a fresh
    consumer record, as a forger with write access would."""
    from mtp.workflow.consumers import open_core_bundle, write_consumer_record

    paths = {name: scan_dir / OUTPUTS[name] for name in ("summary", "provenance")}
    records = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in paths.items()}
    if edit is not None:
        edit(records["summary"], records["provenance"])
    for name, path in paths.items():
        path.write_text(json.dumps(records[name], indent=1) + "\n", encoding="utf-8")
    outputs = [scan_dir / name for name in OUTPUTS.values()]
    write_consumer_record(
        open_core_bundle(core), scan_dir / CONSUMER_RECORD, consumer=consumer, outputs=outputs
    )
    return scan_dir


def test_scopes_report_checks_each_scope_against_its_bound_core(
    cores: dict[str, Path], scans: dict[str, Path], tmp_path: Path
) -> None:
    """A scan of the off core relabelled phenotype, with its consumer record rewritten, is
    refused: the scope must be the bound core's, the records must name that core, the
    record must be the scan's, and a pair needs two cores (#1079)."""
    protocol = load_protocol()

    def relabel(summary: dict[str, Any], provenance: dict[str, Any]) -> None:
        summary[STRAIN_ROLLUP] = provenance[STRAIN_ROLLUP] = scope_record(
            "phenotype", "recorded", protocol
        )

    forged = _rebind(
        _scan(cores["tier1_2_off"], tmp_path / "forged"), cores["tier1_2_off"], relabel
    )
    with pytest.raises(ValueError, match="bound core"):
        _scopes(tmp_path / "a", scans["off"], forged)
    with pytest.raises(ValueError, match="bound core"):
        load_run(forged)

    phenotype_run = json.loads(
        (scans["phenotype"] / OUTPUTS["summary"]).read_text(encoding="utf-8")
    )["core_run"]["run_id"]

    def renamed(summary: dict[str, Any], provenance: dict[str, Any]) -> None:
        summary["core_run"]["run_id"] = provenance["core_run"]["run_id"] = phenotype_run

    misnamed = _rebind(
        _scan(cores["tier1_2_off"], tmp_path / "misnamed"), cores["tier1_2_off"], renamed
    )
    with pytest.raises(ValueError, match="consumer record binds"):
        load_run(misnamed)

    def rehashed(_summary: dict[str, Any], provenance: dict[str, Any]) -> None:
        provenance["core_run"]["manifest_sha256"] = "0" * 64

    unhashed = _rebind(
        _scan(cores["tier1_2_off"], tmp_path / "unhashed"), cores["tier1_2_off"], rehashed
    )
    with pytest.raises(ValueError, match="consumer record binds"):
        load_run(unhashed)

    other = _rebind(
        _scan(cores["tier1_2_off"], tmp_path / "other"), cores["tier1_2_off"], consumer="mtp.other"
    )
    with pytest.raises(ValueError, match=f"not of {CONSUMER}"):
        load_run(other)

    phenotype, off = load_run(scans["phenotype"]), load_run(scans["off"])
    assert (phenotype.core_run_id, off.core_run_id) == (
        phenotype.summary["core_run"]["run_id"],
        off.summary["core_run"]["run_id"],
    )
    with pytest.raises(ValueError, match="same core run"):
        pair_scopes([dataclasses.replace(phenotype, core_run_id=off.core_run_id), off])


def test_main_records_an_absent_scope_end_to_end(
    cores: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A core that records no strain_rollup runs as the protocol's absent_means, and the
    summary and provenance say so. The scopes report checks that against the core itself,
    which here records off, so the run is refused there (#927, #1079)."""
    from mtp.workflow import consumers

    opened = consumers.open_core_bundle

    def legacy(manifest: Path, *args: Any, **kwargs: Any) -> Any:
        bundle = opened(manifest, *args, **kwargs)
        parameters = [item for item in bundle.run.parameters or [] if item.name != STRAIN_ROLLUP]
        run = bundle.run.model_copy(update={"parameters": parameters})
        return dataclasses.replace(bundle, run=run)

    monkeypatch.setattr(consumers, "open_core_bundle", legacy)
    out = _scan(cores["tier1_2_off"], tmp_path / "legacy")
    monkeypatch.undo()
    for name in ("summary", "provenance"):
        record = json.loads((out / OUTPUTS[name]).read_text(encoding="utf-8"))
        assert record[STRAIN_ROLLUP] == scope_record("off", "absent_means", load_protocol())
    with pytest.raises(ValueError, match="bound core"):
        load_run(out)


def _tagged_run(scope: str, role: str, tag: str) -> ScanRun:
    """A scan whose every reported value carries ``tag``, to show where each row came from."""
    traits = pd.DataFrame({"trait_id": [TRAIT_ID, OTHER_TRAIT], "tag": [f"{tag}-1", f"{tag}-2"]})
    summary = {
        STRAIN_ROLLUP: {"scope": scope, "role": role, "scope_source": "recorded"},
        "core_run": {"run_id": f"run-{tag}"},
        "follows_protocol": True,
        "counts": {"tag": tag},
        "lead_triple": {"tag": tag},
        "traits": [
            {"trait_id": OTHER_TRAIT, "tag": f"{tag}-other"},
            {"trait_id": TRAIT_ID, "tag": tag},
        ],
        "rule": {"seed": 1},
    }
    return ScanRun(
        Path(tag), summary, {}, traits, load_protocol(), {"file": tag}, core_run_id=f"run-{tag}"
    )


def test_side_by_side_keeps_each_scope_to_its_own_scan() -> None:
    """Each scope's rows, lead-trait row, counts, core and hashes are its own scan's (#1079)."""
    ordered = [_tagged_run("phenotype", PRIMARY, "ph"), _tagged_run("off", SENSITIVITY, "of")]
    table, report = side_by_side(ordered, load_protocol(), "abc")
    assert table.to_dict("records") == [
        {"strain_rollup": "phenotype", "role": PRIMARY, "trait_id": TRAIT_ID, "tag": "ph-1"},
        {"strain_rollup": "phenotype", "role": PRIMARY, "trait_id": OTHER_TRAIT, "tag": "ph-2"},
        {"strain_rollup": "off", "role": SENSITIVITY, "trait_id": TRAIT_ID, "tag": "of-1"},
        {"strain_rollup": "off", "role": SENSITIVITY, "trait_id": OTHER_TRAIT, "tag": "of-2"},
    ]
    for entry, tag in zip(report["scopes"], ("ph", "of"), strict=True):
        assert entry["lead_trait"] == {"trait_id": TRAIT_ID, "tag": tag}
        assert (entry["counts"], entry["lead_triple"], entry["files"]) == (
            {"tag": tag},
            {"tag": tag},
            {"file": tag},
        )
        assert (entry["core_run_id"], entry["directory"]) == (f"run-{tag}", tag)


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
    return ScanInputs(
        profiles,
        taxa,
        matrix,
        tuple(screen),
        columns,
        {},
        contrast_population=TIER1_2,
        strain_rollup="phenotype",
        strain_rollup_role=PRIMARY,
    )


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
    _, inputs = load_core_inputs(cores["tier1_2_off"])
    assert len(inputs.studies) == sum(STUDIES.values())
    assert TRAIT_ID in inputs.columns
    assert (inputs.strain_rollup, inputs.strain_rollup_role) == ("off", SENSITIVITY)
    _, phenotype = load_core_inputs(cores["tier1_2_phenotype"])
    assert (phenotype.strain_rollup, phenotype.strain_rollup_role) == ("phenotype", PRIMARY)
    with pytest.raises(ValueError, match="neither the declared primary"):
        load_core_inputs(cores["tier1_2_full"])
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
    # The z rank places it the same way, among the disjoint triples with a z (#907).
    z = disjoint["u6_z_mean_distance"].dropna()
    lead_z = lead_row["u6_z_mean_distance"]
    assert row[f"{MEAN_DISTANCE}_{DISJOINT}_z_triples"] == len(z) + 1
    assert row[f"{MEAN_DISTANCE}_{DISJOINT}_z_rank_max"] == int((z <= lead_z).sum()) + 1


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
    # Exact ties by default, as #854 declares: 0.3 is smaller than 0.3 + 1e-13 (#907).
    assert _keys(closest_disjoint_triples(rows)) == [("B", "C", "D", "x"), ("A", "B", "C", "y")]
    # With a tolerance the two tie, and the lexically first conditions win.
    chosen = closest_disjoint_triples(rows, tolerance=1e-12)
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
    # An unknown study is refused rather than dropped (#907).
    with pytest.raises(ValueError, match="not in the population"):
        build_population(experiments, keyed, reported, matrix_taxa, [*keep, "bsdb:999999"])


# --------------------------------------------------------------------------- empty families (#901)


def _sparse_corpus() -> ScanInputs:
    """Nine conditions of two studies each, so none has three; every decreased signature
    carries cellulose-annotated taxa and every increased one does not."""
    records: list[Signature] = []
    for number in range(1, 19):
        study = f"bsdb:{number}"
        decreased = ["t0", "t1", f"t{number % 10}"]
        records.append((study, f"{study}-1", f"C{(number - 1) // 2}", FECES, decreased, ["t10"]))
    return _inputs(records, _matrix(20, {TRAIT_ID: range(10), OTHER_TRAIT: range(10, 20)}))


def test_select_without_eligible_conditions_returns_an_empty_family() -> None:
    """Sparse labels under the default minimum: the screen selects a trait, no condition
    is eligible, and the selection is an empty family with a status, not an error (#901)."""
    inputs = _sparse_corpus()
    ordered = sorted(inputs.studies, key=lambda study: int(study.split(":")[1]))
    for studies in (ordered[::2], ordered):
        chosen = select(
            inputs, studies, permutations=19, seed=(20261022, 0, 0), always_included=False
        )
        assert chosen.traits["trait_id"].tolist() == [TRAIT_ID]
        assert not chosen.conditions["eligible"].any()
        assert chosen.triples.empty
        assert list(chosen.triples.columns) == list(TRIPLE_COLUMNS)
        assert chosen.status == TOO_FEW_CONDITIONS
        assert closest_disjoint_triples(chosen.triples).empty
        lead, histogram, _ = lead_and_distribution(chosen.triples, [TRAIT_ID], 20)
        assert lead["lead_in_scan"].tolist() == [False]
        assert histogram["triples"].sum() == 0


def test_select_with_a_minimum_no_condition_meets_returns_an_empty_family() -> None:
    """The same with a minimum no condition can meet (#901)."""
    inputs, train, _ = _split_corpus()
    protocol = dataclasses.replace(load_protocol(), min_studies=100)
    chosen = select(inputs, train, protocol=protocol, permutations=19, seed=1)
    assert not chosen.traits.empty
    assert chosen.conditions["eligible"].sum() == 0
    assert chosen.triples.empty
    assert chosen.status == TOO_FEW_CONDITIONS
    prepared = prepare(inputs.profiles, inputs.taxa, inputs.matrix, [TRAIT_ID], [])
    report = plan(chosen.conditions, chosen.traits, prepared, protocol, 100)
    assert (report["triples"], report["estimated_null_minutes"]) == (0, 0)


def test_evaluate_with_no_drawable_condition_keeps_every_hypothesis() -> None:
    """When no condition of any hypothesis has a held-out study, each hypothesis is still
    returned as untestable (#901)."""
    inputs, _, test = _split_corpus()
    asthma = set(inputs.profiles.loc[inputs.profiles["condition"].eq("Asthma"), "study_id"])
    held_out = [study for study in test if study in asthma]
    hypotheses = [(PD, HIV, CRC, TRAIT_ID), (PD, "Gout", "Lupus", OTHER_TRAIT)]
    tested = evaluate(inputs, held_out, hypotheses, permutations=19, seed=1)
    assert tested["test_status"].tolist() == [NO_UNITS, NO_UNITS]
    assert tested["missing_condition"].all()
    assert tested["null_undefined"].all()


def test_selection_status_names_why_a_family_is_empty() -> None:
    inputs, train, _ = _split_corpus()
    assert select(inputs, train, permutations=19, seed=1).status == SELECTED
    flat = tuple(
        dataclasses.replace(profile, guild_counts=dict.fromkeys(profile.guild_counts, 1), total=3)
        for profile in inputs.screen
    )
    empty = select(
        dataclasses.replace(inputs, screen=flat), train, permutations=19, always_included=False
    )
    assert empty.status == NO_TRAITS


# --------------------------------------------------------------------------- delta review (#907)


def test_select_equals_select_on_a_population_of_those_studies() -> None:
    """``select`` on a subset equals ``select`` on a population holding only that subset,
    so nothing outside it reaches the selection, the U6 pool included (#907)."""
    inputs, train, _ = _split_corpus()
    seed = (20261022, 3, 0)
    chosen = select(inputs, train, permutations=300, seed=seed)
    alone = select(inputs.restricted(train), permutations=300, seed=seed)
    assert chosen.studies == alone.studies
    pd.testing.assert_frame_equal(chosen.triples, alone.triples)
    pd.testing.assert_frame_equal(chosen.traits, alone.traits)
    pd.testing.assert_frame_equal(chosen.conditions, alone.conditions)


def test_evaluate_equals_evaluate_on_a_population_of_those_studies() -> None:
    """The same for the held-out test (#907)."""
    inputs, train, test = _split_corpus()
    ranked = rank_candidates(select(inputs, train, permutations=200, seed=1).triples)
    hypotheses = _keys(ranked.head(3))
    seed = (20261022, 3, 1)
    tested = evaluate(inputs, test, hypotheses, permutations=300, seed=seed)
    alone = evaluate(inputs.restricted(test), test, hypotheses, permutations=300, seed=seed)
    pd.testing.assert_frame_equal(tested, alone)


def test_traits_rank_by_the_lead_z_share_ascending_apart_from_the_p_share() -> None:
    """The traits rank by the share of triples at least as extreme in z as the lead's,
    smallest first, separately from the p share (#907)."""
    inputs, _, _ = _split_corpus()
    traits = [TRAIT_ID, OTHER_TRAIT, "METPO:1000623"]
    eligible = sorted(condition_table(inputs.profiles, 3).query("eligible")["condition"])
    result = scan(
        inputs.profiles,
        inputs.taxa,
        inputs.matrix,
        traits,
        eligible,
        permutations=10,
        seed=1,
        tolerance=0,
    )
    lead = pd.DataFrame(
        {
            "trait_id": traits,
            "lead_u6_p_mean_distance": [0.2, 0.4, 0.6],
            "lead_u6_p_core": [0.6, 0.4, 0.2],
            f"{MEAN_DISTANCE}_{ALL_TRIPLES}_share_at_most": [0.1, 0.5, 0.9],
            f"{MEAN_DISTANCE}_{ALL_TRIPLES}_z_share_at_most": [0.5, 0.9, 0.1],
            f"{CORE}_{ALL_TRIPLES}_share_at_most": [0.9, 0.5, 0.1],
            f"{CORE}_{ALL_TRIPLES}_z_share_at_most": [0.5, 0.1, 0.9],
        }
    )
    chosen = pd.DataFrame({"trait_id": traits})
    frame = trait_table(chosen, lead, result, TRAIT_ID)
    ranks = {
        basis: frame[f"rank_among_traits_by_{basis}"].tolist()
        for basis in (
            f"lead_share_{MEAN_DISTANCE}",
            f"lead_z_share_{MEAN_DISTANCE}",
            f"lead_p_{CORE}",
            f"lead_z_share_{CORE}",
        )
    }
    assert ranks == {
        f"lead_share_{MEAN_DISTANCE}": [1, 2, 3],
        f"lead_z_share_{MEAN_DISTANCE}": [2, 3, 1],
        f"lead_p_{CORE}": [3, 2, 1],
        f"lead_z_share_{CORE}": [2, 1, 3],
    }


def test_scan_inputs_require_the_population_and_a_declared_scope() -> None:
    """Inputs name their population and scope, and select and evaluate check the scope's
    role against the protocol (#907)."""
    inputs, train, test = _split_corpus()
    with pytest.raises(ValueError, match="population is"):
        dataclasses.replace(inputs, contrast_population="all_paired")
    with pytest.raises(ValueError, match="unknown strain-rollup role"):
        dataclasses.replace(inputs, strain_rollup_role="exploratory")
    full = dataclasses.replace(inputs, strain_rollup="full")
    mislabelled = dataclasses.replace(inputs, strain_rollup="off")
    for wrong, message in (
        (full, "neither the declared primary"),
        (mislabelled, "the protocol's declared_sensitivity"),
    ):
        with pytest.raises(ValueError, match=message):
            select(wrong, train, permutations=10, seed=1)
        with pytest.raises(ValueError, match=message):
            evaluate(wrong, test, [(PD, HIV, CRC, TRAIT_ID)], permutations=10, seed=1)
    sensitivity = dataclasses.replace(inputs, strain_rollup="off", strain_rollup_role=SENSITIVITY)
    assert select(sensitivity, train, permutations=10, seed=1).status == SELECTED


def test_scan_inputs_have_no_default_population_or_scope() -> None:
    """Inputs cannot be built without their population, scope and role, so a default
    scope cannot come back unnoticed (#927)."""
    required = {
        field.name
        for field in dataclasses.fields(ScanInputs)
        if field.default is dataclasses.MISSING and field.default_factory is dataclasses.MISSING
    }
    assert {"contrast_population", "strain_rollup", "strain_rollup_role"} <= required


def test_load_inputs_checks_the_stated_population_and_scope(cores: dict[str, Path]) -> None:
    bundle, _ = load_core_inputs(cores["tier1_2_off"])
    paths = {product: bundle.path(product) for product in CORE_PRODUCTS}
    with pytest.raises(ValueError, match="population is"):
        load_inputs(paths, contrast_population="all_paired", strain_rollup="off")
    with pytest.raises(ValueError, match="neither the declared primary"):
        load_inputs(paths, contrast_population=TIER1_2, strain_rollup="full")
    inputs = load_inputs(paths, contrast_population=TIER1_2, strain_rollup="phenotype")
    assert inputs.strain_rollup_role == PRIMARY
