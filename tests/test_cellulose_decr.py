"""cellulose_decr overlap, null-model and direction statistics (#552, #557)."""

from __future__ import annotations

import itertools
import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats
from scipy.stats import false_discovery_control
from scipy.stats.contingency import odds_ratio

from mtp.cellulose_decr.analysis import (
    ALL_CONTRASTS,
    FECES,
    INPUTS,
    NO_MSM,
    NON_FECAL_GI,
    ORAL,
    ORIENTED,
    OTHER,
    SHORT,
    STRONG,
    TRAIT_ID,
    Corpus,
    NullDraws,
    benchmark_table,
    benjamini_hochberg,
    binomial_power,
    case_control_inventory,
    catalog_table,
    check_input_drift,
    choose,
    clean,
    compartment_of,
    condition_sets,
    direct_trait_subjects,
    direction_row,
    hypergeom_pair,
    interaction_interpretation,
    interaction_table,
    label_corrections,
    load_corpus,
    merge_synonyms,
    most_common_label,
    msm_mentioned,
    net_decreased,
    paired_deltas,
    population,
    scenario_corpora,
    scenario_tables,
    set_summary,
    shared_trait_null,
    sign_tests,
    signature_null,
    stratified_exact,
    strictly_net_decreased,
    study_votes,
    subset_null,
    three_way_exact,
    u5_weights,
    u6_units,
    without_studies,
)
from mtp.contrast_orientation import (
    CASE_VS_CONTROL_SUBJECT,
    CLASSIFIER_VERSION,
    COMORBIDITY_WITHIN_DISEASE,
    PAIRED_TISSUE_AUTOLOGOUS,
    TIER_1_HEALTHY_NORMAL,
    TIER_2_CONDITION_NEGATIVE_SCREENED,
    TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
)

# --------------------------------------------------------------------------- exact tests


@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_hypergeom_pair_matches_enumeration(k: int) -> None:
    """The overlap p-value and expected distance equal a full enumeration of B."""
    a = {0, 1, 2}
    draws = [set(b) for b in itertools.combinations(range(8), 4)]
    p = sum(len(a & b) >= k for b in draws) / len(draws)
    distance = sum(1 - len(a & b) / len(a | b) for b in draws) / len(draws)
    result = hypergeom_pair(8, 3, 4, k)
    assert result["p"] == pytest.approx(p)
    assert result["expected_jaccard_distance"] == pytest.approx(distance)
    assert result["expected_intersection"] == pytest.approx(1.5)


def test_hypergeom_pair_odds_ratio() -> None:
    """k(N - a - b + k) / ((a - k)(b - k)): 3 * 4 / (1 * 2) = 6."""
    assert hypergeom_pair(10, 4, 5, 3)["odds_ratio"] == pytest.approx(6.0)


@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_three_way_exact_matches_enumeration(k: int) -> None:
    """P(|A∩B∩C| >= k) equals the share of all (B, C) draws for a fixed A."""
    a = {0, 1, 2}
    cores = [
        len(a & set(b) & set(c))
        for b in itertools.combinations(range(7), 4)
        for c in itertools.combinations(range(7), 2)
    ]
    expected = sum(core >= k for core in cores) / len(cores)
    result = three_way_exact(7, (3, 4, 2), k)
    assert result["p_3way_exact"] == pytest.approx(expected)
    # The expected core is the enumeration's mean (#632)
    assert result["expected_3way_intersection"] == pytest.approx(sum(cores) / len(cores))


# --------------------------------------------------------------------------- Monte Carlo


@pytest.mark.parametrize("weighted", [False, True])
def test_choose_marks_exactly_k_items(weighted: bool) -> None:
    rng = np.random.default_rng(1)
    weights = np.linspace(0.5, 3.0, 9) if weighted else None
    for k in (0, 1, 4, 9):
        mask = choose(rng, 200, 9, k, weights)
        assert mask.shape == (200, 9)
        assert (mask.sum(axis=1) == k).all()


def test_choose_weights_follow_successive_sampling() -> None:
    """With one pick, an item is chosen in proportion to its weight."""
    mask = choose(np.random.default_rng(2), 20_000, 2, 1, np.array([1.0, 9.0]))
    assert mask[:, 1].mean() == pytest.approx(0.9, abs=0.01)


def test_signature_null_draws_distinct_studies() -> None:
    """U6 gives each synthetic condition whole units from studies no other one uses.

    Every unit of a study carries that study's taxon, so two units of one study in
    the same draw, within or across conditions, would show up as a shared taxon.
    """
    n_studies, per_study = 6, 2
    unit_study = np.repeat(np.arange(n_studies), per_study).astype(np.int64)
    units = np.zeros((unit_study.size, n_studies + unit_study.size), dtype=bool)
    for unit, study in enumerate(unit_study):
        units[unit, study] = True
        units[unit, n_studies + unit] = True
    counts = {"PD": 2, "HIV": 2, "CRC": 2}
    draws = signature_null(np.random.default_rng(3), units, unit_study, counts, 500)
    for name, values in draws.pair_intersection.items():
        assert (np.concatenate(values) == 0).all(), name
    # Two studies per condition, each adding its own taxon and one unit's taxon.
    for short, values in draws.sizes.items():
        assert (np.concatenate(values) == 4).all(), short


def test_signature_null_draws_studies_and_units_at_random() -> None:
    """With fewer studies drawn than exist, every study and every unit within a study is
    drawn about equally often (#595).

    Unit k carries its study's taxon and k + 1 taxa of its own, so the size of a
    condition's set when it draws one unit is k + 2 and identifies the unit.
    """
    n_studies, per_study = 6, 2
    unit_study = np.repeat(np.arange(n_studies), per_study).astype(np.int64)
    n_units = unit_study.size
    units = np.zeros((n_units, n_studies + n_units * (n_units + 1) // 2), dtype=bool)
    offset = n_studies
    for unit, study in enumerate(unit_study):
        units[unit, study] = True
        units[unit, offset : offset + unit + 1] = True
        offset += unit + 1
    counts = {"PD": 1, "HIV": 1, "CRC": 1}
    draws = signature_null(np.random.default_rng(4), units, unit_study, counts, 6000)
    for name, values in draws.pair_intersection.items():
        assert (np.concatenate(values) == 0).all(), name
    sizes = np.concatenate(draws.sizes["PD"])
    shares = np.bincount(sizes - 2, minlength=n_units) / sizes.size
    assert shares == pytest.approx(np.full(n_units, 1 / n_units), abs=0.02)


def test_u4_and_within_draws_stay_inside_each_condition() -> None:
    """U4 puts one random trait on each condition's own decreased taxa, and the within
    reference draws each set from its own condition's taxa (#632)."""
    width = 9
    base = {
        short: np.isin(np.arange(width), np.arange(3 * index, 3 * index + 3))
        for index, short in enumerate(SHORT)
    }
    shared = shared_trait_null(np.random.default_rng(5), base, 4, 300)
    within = subset_null(np.random.default_rng(6), base, dict.fromkeys(SHORT, 2), 300)
    # The three conditions share no taxa, so neither null can give them one.
    for draws in (shared, within):
        for name, values in draws.pair_intersection.items():
            assert (np.concatenate(values) == 0).all(), name
    # The trait's 4 taxa are split among the three sets.
    totals = np.vstack([np.concatenate(shared.sizes[short]) for short in SHORT]).sum(axis=0)
    assert (totals == 4).all()
    for short in SHORT:
        assert (np.concatenate(within.sizes[short]) == 2).all()


def test_u5_weights_count_other_studies_plus_a_half() -> None:
    """A U3 taxon weighs the studies outside the three conditions that report it decreased
    in feces, plus 0.5, so a taxon only the three conditions report can be drawn (#632)."""
    rows = [
        ("bsdb:1", "Asthma", "dec", "a"),
        ("bsdb:2", "Asthma", "dec", "a"),
        ("bsdb:2", "Asthma", "dec", "a"),
        ("bsdb:3", "Parkinson's disease", "dec", "b"),
        ("bsdb:4", "Asthma", "inc", "b"),
        ("bsdb:5", "Asthma", "dec", "c"),
    ]
    feces_taxa = pd.DataFrame(rows, columns=["study_id", "condition", "dir", "taxon_id"])
    # c is not annotated, and d is never reported decreased, so U3 is a and b.
    assert list(u5_weights(feces_taxa, {"a", "b", "d"})) == [2.5, 0.5]


def test_size_matched_draws_keep_sizes_within_tolerance() -> None:
    """Draws are kept when every set size is within max(2, 10%) of the observed size."""
    draws = NullDraws.empty()
    draws.sizes = {short: [np.array([18, 23, 20, 22])] for short in SHORT}
    draws.core = [np.array([1, 2, 3, 4])]
    draws.mean_distance = [np.array([0.1, 0.2, 0.3, 0.4])]
    for name in draws.pair_distance:
        draws.pair_distance[name] = [np.array([0.1, 0.2, 0.3, 0.4])]
        draws.pair_intersection[name] = [np.array([1, 2, 3, 4])]
    matched, kept = draws.size_matched({"sizes": dict.fromkeys(SHORT, 20)}, 0.1)
    assert kept == 3
    assert list(np.concatenate(matched.core)) == [1, 3, 4]


def test_null_p_values_count_the_observed_value() -> None:
    """Monte Carlo p-values are (1 + draws at least as extreme) / (1 + draws)."""
    draws = NullDraws.empty()
    draws.mean_distance = [np.array([0.2, 0.5, 0.8, 0.9])]
    draws.core = [np.array([0, 1, 3, 5])]
    for name in draws.pair_distance:
        draws.pair_distance[name] = [np.array([0.1, 0.4, np.nan, 0.9])]
        draws.pair_intersection[name] = [np.array([3, 2, 0, 1])]
    draws.sizes = {short: [np.array([3, 3, 3, 3])] for short in SHORT}
    observed = {
        "jaccard_distance": dict.fromkeys(draws.pair_distance, 0.4),
        "mean_pairwise_jaccard_distance": 0.5,
        "three_way_core": 3,
    }
    summary = draws.summary(observed)
    assert summary["iterations"] == 4
    assert summary["mean_pairwise_jaccard_distance"]["p"] == pytest.approx(3 / 5)
    assert summary["three_way_core"]["p"] == pytest.approx(3 / 5)
    # A draw with an empty synthetic set has no distance and leaves that pair's p.
    assert summary["pairs"]["PD vs. HIV"]["p"] == pytest.approx(3 / 4)

    beyond_every_draw = {**observed, "mean_pairwise_jaccard_distance": 0.1, "three_way_core": 6}
    summary = draws.summary(beyond_every_draw)
    assert summary["mean_pairwise_jaccard_distance"]["p"] == pytest.approx(1 / 5)
    assert summary["three_way_core"]["p"] == pytest.approx(1 / 5)


# --------------------------------------------------------------------------- direction tests


def _deltas(values: dict[str, list[float]]) -> pd.Series:
    keys = [(study, f"{study}-{i}") for study, deltas in values.items() for i in range(len(deltas))]
    index = pd.MultiIndex.from_tuples(keys, names=["study_id", "experiment_id"])
    return pd.Series([delta for deltas in values.values() for delta in deltas], index=index)


def test_study_votes_by_mean_and_majority() -> None:
    """A study's mean and its majority of experiment signs can disagree; ties vote 0."""
    deltas = _deltas({"s1": [-0.2, -0.1, 0.5], "s2": [0.1, -0.1], "s3": [0.0], "s4": [0.3]})
    means = study_votes(deltas)
    assert means.to_dict() == pytest.approx({"s1": 0.2 / 3, "s2": 0.0, "s3": 0.0, "s4": 0.3})
    majority = study_votes(deltas, "majority")
    assert majority.to_dict() == {"s1": -1.0, "s2": 0.0, "s3": 0.0, "s4": 1.0}

    tests = sign_tests(majority)
    assert (tests["down"], tests["up"], tests["tied"]) == (1, 1, 2)
    assert tests["binom_p_two_sided"] == pytest.approx(1.0)


def test_sign_tests_leave_out_ties() -> None:
    """Nine down, none up and one tie: one-sided p = 2**-9, two-sided 2**-8."""
    tests = sign_tests(pd.Series([-0.1] * 9 + [0.0]))
    assert (tests["down"], tests["up"], tests["tied"]) == (9, 0, 1)
    assert tests["share_down_untied"] == 1.0
    assert tests["binom_p_down_one_sided"] == pytest.approx(2**-9)
    assert tests["binom_p_two_sided"] == pytest.approx(2**-8)


def test_paired_deltas_need_both_directions() -> None:
    profiles = pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s1", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e3"],
            "dir": ["inc", "dec", "dec", "inc"],
            "frac": [0.1, 0.4, 0.3, 0.2],
        }
    )
    assert paired_deltas(profiles).to_dict() == {("s1", "e1"): pytest.approx(-0.3)}
    with pytest.raises(ValueError, match="two profiles in the same direction"):
        paired_deltas(pd.concat([profiles, profiles.iloc[[0]]]))


@pytest.mark.parametrize(
    "table",
    [
        [[32, 6], [4, 11]],
        [[7, 2], [3, 9]],
        [[5, 0], [2, 4]],
        [[0, 3], [4, 2]],
    ],
)
def test_stratified_exact_one_stratum_matches_scipy(table: list[list[int]]) -> None:
    """With one stratum: Fisher's p, and scipy's conditional odds ratio and interval."""
    result = stratified_exact([table])
    expected = odds_ratio(table, kind="conditional")
    interval = expected.confidence_interval(0.95)
    assert result["p"] == pytest.approx(stats.fisher_exact(table).pvalue)
    assert result["conditional_or"] == pytest.approx(expected.statistic, rel=1e-6)
    assert result["low"] == pytest.approx(interval.low, rel=1e-6)
    assert result["high"] == pytest.approx(interval.high, rel=1e-6)
    (a, b), (c, d) = table
    if b * c:
        assert result["mantel_haenszel_or"] == pytest.approx(a * d / (b * c))


def test_stratified_exact_matches_enumeration() -> None:
    """Three strata: the p-value and the conditional MLE from the full exact null."""
    tables = [[[3, 1], [1, 2]], [[2, 2], [0, 3]], [[4, 0], [2, 2]]]
    supports = []
    observed = 0
    for (a, b), (c, d) in tables:
        total, row, column = a + b + c + d, a + b, a + c
        support = range(max(0, row - (total - column)), min(row, column) + 1)
        supports.append({x: stats.hypergeom.pmf(x, total, column, row) for x in support})
        observed += a

    def null(odds: float) -> dict[int, float]:
        weights: dict[int, float] = {}
        for combo in itertools.product(*(pmf.items() for pmf in supports)):
            total = sum(x for x, _ in combo)
            weight = math.prod(p for _, p in combo) * odds**total
            weights[total] = weights.get(total, 0.0) + weight
        scale = sum(weights.values())
        return {total: weight / scale for total, weight in weights.items()}

    central = null(1.0)
    p = sum(prob for prob in central.values() if prob <= central[observed] * (1 + 1e-7))
    result = stratified_exact(tables)
    assert result["p"] == pytest.approx(p)
    at_estimate = null(result["conditional_or"])
    assert sum(total * prob for total, prob in at_estimate.items()) == pytest.approx(observed)
    numerator = sum(a * d / (a + b + c + d) for (a, b), (c, d) in tables)
    denominator = sum(b * c / (a + b + c + d) for (a, b), (c, d) in tables)
    assert result["mantel_haenszel_or"] == pytest.approx(numerator / denominator)


def test_stratified_exact_without_informative_strata() -> None:
    """Strata whose margins fix the top-left cell carry no information."""
    result = stratified_exact([[[0, 0], [0, 0]], [[3, 0], [0, 0]]])
    assert math.isnan(result["p"])
    assert math.isnan(result["conditional_or"])
    assert math.isnan(result["mantel_haenszel_or"])
    assert stratified_exact([[[2, 0], [0, 3]]])["mantel_haenszel_or"] == math.inf


def test_benjamini_hochberg_matches_scipy_and_keeps_missing() -> None:
    p = pd.Series([0.01, np.nan, 0.04, 0.03, 0.5, 0.001], index=list("abcdef"))
    q = benjamini_hochberg(p)
    present = p.dropna()
    assert q[present.index].to_numpy() == pytest.approx(false_discovery_control(present.to_numpy()))
    assert math.isnan(q["b"])
    stepped = benjamini_hochberg(pd.Series([0.01, 0.011, 0.04, 0.5]))
    assert list(stepped) == pytest.approx([0.022, 0.022, 0.04 * 4 / 3, 0.5])


# --------------------------------------------------------------------------- inputs


@pytest.mark.parametrize(
    ("site", "compartment"),
    [
        ("Feces", FECES),
        ("feces", FECES),
        ("Feces,Blood", OTHER),
        ("Oral cavity,Feces", OTHER),
        ("Saliva", ORAL),
        ("Tongue", ORAL),
        ("Sigmoid colon", NON_FECAL_GI),
        ("Rectal mucosa", NON_FECAL_GI),
        ("Colonic mucosa", NON_FECAL_GI),
        ("Blood", OTHER),
        ("Vagina", OTHER),
        ("Unspecified", OTHER),
    ],
)
def test_compartment_of(site: str, compartment: str) -> None:
    assert compartment_of(site) == compartment


def _hashes(**digests: str) -> dict[str, dict[str, object]]:
    return {
        key: {"path": f"data/{key}.tsv", "sha256": digest, "bytes": 1}
        for key, digest in digests.items()
    }


def test_input_drift_is_refused_unless_allowed_or_repinned(tmp_path: Path) -> None:
    """Pins persist through --allow-input-drift; only --repin replaces them (#562)."""
    pinned = tmp_path / "provenance.json"
    pins = _hashes(a="aaa", b="bbb")
    pinned.write_text(json.dumps({"inputs": _hashes(a="zzz", b="bbb"), "pinned_inputs": pins}))
    assert check_input_drift(pins, pinned, allow=False, repin=False) == ([], pins)

    now = _hashes(a="aaa", b="ccc", c="ddd")
    with pytest.raises(ValueError, match="b: bbb pinned, ccc now"):
        check_input_drift(now, pinned, allow=False, repin=False)
    drifted, kept = check_input_drift(now, pinned, allow=True, repin=False)
    assert drifted == ["b"]
    assert kept == {"a": pins["a"], "b": pins["b"], "c": now["c"]}
    assert check_input_drift(now, pinned, allow=False, repin=True) == ([], now)


def test_input_pins_need_a_provenance_or_repin(tmp_path: Path) -> None:
    """A missing pin file is an error, and a provenance without pins uses its inputs."""
    now = _hashes(a="aaa")
    absent = tmp_path / "absent.json"
    with pytest.raises(FileNotFoundError, match="--repin"):
        check_input_drift(now, absent, allow=True, repin=False)
    assert check_input_drift(now, absent, allow=False, repin=True) == ([], now)
    legacy = tmp_path / "legacy.json"
    legacy.write_text(json.dumps({"inputs": _hashes(a="old")}))
    with pytest.raises(ValueError, match="a: old pinned, aaa now"):
        check_input_drift(now, legacy, allow=False, repin=False)


# Study bsdb:1 has two PD experiments and one HIV experiment, so its modal condition is
# PD. Study bsdb:2 has two colon experiments and one feces experiment, so its modal site
# is colon. Since #555 study_trait_profile_study_metadata.tsv carries each experiment's
# own labels, so the study-modal labels must come from full_dump.csv (#594).
EXPERIMENTS = [
    ("bsdb:1", "bsdb:1-1", "Parkinson's disease", "Feces", "disease_vs_healthy", "native"),
    ("bsdb:1", "bsdb:1-2", "Parkinson's disease", "Feces", "disease_vs_healthy", "reversed"),
    ("bsdb:1", "bsdb:1-3", "HIV infection", "Feces", "treatment_exposure", "not_applicable"),
    ("bsdb:2", "bsdb:2-1", "Colorectal cancer", "Colon", "disease_vs_healthy", "native"),
    ("bsdb:2", "bsdb:2-2", "Colorectal cancer", "Colon", "disease_vs_healthy", "native"),
    ("bsdb:2", "bsdb:2-3", "Colorectal cancer", "Feces", "disease_vs_healthy", "native"),
]
HALOPHILIC = "METPO:1000623"
MEDIUM = "mediadive.medium:1"


def _write(path: Path, frame: pd.DataFrame, index: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=index)


def _input_root(root: Path, classifier_version: str = CLASSIFIER_VERSION) -> dict[str, Path]:
    """The analysis inputs that load_corpus reads, for the six experiments above."""
    paths = {key: root / relative for key, relative in INPUTS.items()}
    profiles, counts, taxa, dump = [], [], [], []
    for index, (study, experiment, condition, site, *_rest) in enumerate(EXPERIMENTS):
        for signature, (inc_dec, cellulose) in enumerate((("inc", index), ("dec", 6 - index)), 1):
            other, medium = 2, 5
            profiles.append(
                {
                    "profile_id": f"{study}|{experiment}|{inc_dec}",
                    "study_id": study,
                    "experiment_id": experiment,
                    "inc_dec": inc_dec,
                    "condition": condition,
                    "body_site": site,
                    "total_trait_signal": cellulose + other,
                }
            )
            counts.append(
                {
                    "study_id": study,
                    "experiment_id": experiment,
                    "inc_dec": inc_dec,
                    TRAIT_ID: cellulose,
                    HALOPHILIC: other,
                    MEDIUM: medium,
                }
            )
            direction = "increased" if inc_dec == "inc" else "decreased"
            for taxon in ("t1", "t2"):
                taxa.append(
                    {
                        "study_id": study,
                        "experiment_id": experiment,
                        "direction": direction,
                        "taxon_id": taxon,
                    }
                )
            number = experiment.rsplit("-", 1)[1]
            dump.append(
                {
                    "BSDB ID": f"{study}/{number}/{signature}",
                    "Condition": condition,
                    "Body site": site,
                }
            )
    _write(paths["profile_metadata"], pd.DataFrame(profiles), index=False)
    _write(paths["profile_counts"], pd.DataFrame(counts), index=False)
    _write(paths["signature_taxa"], pd.DataFrame(taxa), index=False)
    paths["full_dump"].parent.mkdir(parents=True, exist_ok=True)
    paths["full_dump"].write_text(
        "# BugSigDB 2026-04-24_00:41_UTC\n" + pd.DataFrame(dump).to_csv(index=False)
    )
    columns = ["study_id", "experiment_id", "condition", "body_site"]
    experiments = pd.DataFrame([row[:4] for row in EXPERIMENTS], columns=columns)
    _write(paths["experiments"], experiments.assign(host_species="Homo sapiens"), index=False)
    manifest = pd.DataFrame(
        [(*row[:2], *row[4:]) for row in EXPERIMENTS],
        columns=["study_id", "experiment_id", "contrast_class", "orientation"],
    ).assign(
        classifier_version=classifier_version,
        contrast_type=CASE_VS_CONTROL_SUBJECT,
        control_tier=TIER_1_HEALTHY_NORMAL,
        flags="",
        group_0_name="Controls",
        group_1_name="Cases",
        group_1_definition="",
        description="",
    )
    _write(paths["contrast_manifest"], manifest, index=False)
    matrix = pd.DataFrame(
        {"taxonomic_rank": ["genus", "species"], TRAIT_ID: [1, 0], HALOPHILIC: [0, 1]},
        index=pd.Index(["t1", "t2"], name="taxon_id"),
    )
    _write(paths["taxa_trait_matrix"], matrix)
    trait_columns = pd.DataFrame(
        {
            "column_id": [TRAIT_ID, HALOPHILIC],
            "object_name": ["cellulose catabolic process", "moderately halophilic"],
        }
    )
    _write(paths["trait_columns"], trait_columns, index=False)
    return paths


@pytest.mark.parametrize("label_source", ["experiment", "study_modal"])
def test_load_corpus_label_sources(tmp_path: Path, label_source: str) -> None:
    """``experiment`` labels each profile with its own experiment, ``study_modal`` with the
    study's most common label in full_dump.csv, as the October 16 outputs did, even though
    the metadata table now holds per-experiment labels (#594)."""
    corpus = load_corpus(_input_root(tmp_path), label_source)
    profiles = corpus.profiles
    if label_source == "experiment":
        assert profiles.loc["bsdb:1|bsdb:1-3|dec", "condition"] == "HIV infection"
        assert profiles.loc["bsdb:2|bsdb:2-3|dec", "compartment"] == FECES
    else:
        assert profiles.loc["bsdb:1|bsdb:1-3|dec", "condition"] == "Parkinson's disease"
        assert profiles.loc["bsdb:2|bsdb:2-3|dec", "compartment"] == NON_FECAL_GI
    # Signature taxa carry their profile's labels.
    labelled = corpus.taxa.set_index("profile_id")
    assert (labelled["condition"] == profiles.loc[labelled.index, "condition"]).all()
    # The share is over retained columns only; the MediaDive column is left out.
    assert profiles.loc["bsdb:1|bsdb:1-1|dec", "frac"] == pytest.approx(6 / 8)


def test_label_corrections_compare_study_modal_and_experiment_labels(tmp_path: Path) -> None:
    corrections = label_corrections(_input_root(tmp_path))
    assert list(corrections["experiment_id"]) == ["bsdb:1-3", "bsdb:2-3"]
    assert list(corrections["condition_study_modal"]) == [
        "Parkinson's disease",
        "Colorectal cancer",
    ]


def test_most_common_label_ignores_blanks_and_breaks_ties_alphabetically() -> None:
    assert most_common_label(["B", "A", "B", "A", "NA", " ", "C"]) == "A"
    assert most_common_label(["nan", "", "NA"]) == "unknown"
    assert most_common_label([" Feces ", "Feces", "Colon"]) == "Feces"


def test_load_corpus_refuses_another_classifier_version(tmp_path: Path) -> None:
    paths = _input_root(tmp_path, classifier_version="contrast-orientation-v0")
    with pytest.raises(ValueError, match="contrast manifest written by"):
        load_corpus(paths, "experiment")


def test_oriented_population_flips_reversed_experiments(tmp_path: Path) -> None:
    """Only disease-vs-healthy experiments are kept, and reversed ones swap inc and dec."""
    corpus = load_corpus(_input_root(tmp_path), "experiment")
    profiles, taxa = population(corpus, ORIENTED)
    assert {"bsdb:1|bsdb:1-3|inc", "bsdb:1|bsdb:1-3|dec"}.isdisjoint(profiles.index)
    assert profiles.loc["bsdb:1|bsdb:1-1|inc", "dir"] == "inc"
    assert profiles.loc["bsdb:1|bsdb:1-2|inc", "dir"] == "dec"
    assert profiles.loc["bsdb:1|bsdb:1-2|dec", "dir"] == "inc"
    assert set(taxa["profile_id"]) <= set(profiles.index)


def test_strong_population_drops_weak_comparisons(tmp_path: Path) -> None:
    corpus = load_corpus(_input_root(tmp_path), "experiment")
    manifest = corpus.manifest.copy()
    manifest.loc[("bsdb:2", "bsdb:2-1"), "flags"] = "weak_comparison_evidence"
    corpus = replace(corpus, manifest=manifest)
    oriented, _ = population(corpus, ORIENTED)
    strong, _ = population(corpus, STRONG)
    assert "bsdb:2|bsdb:2-1|dec" in oriented.index
    assert "bsdb:2|bsdb:2-1|dec" not in strong.index
    assert len(strong) == len(oriented) - 2


def test_synonym_labels_merge_under_their_disease() -> None:
    merged = merge_synonyms(
        pd.Series(["Human immunodeficiency virus", "Colorectal carcinoma", "Asthma"])
    )
    assert list(merged) == ["HIV infection", "Colorectal cancer", "Asthma"]


def test_direct_trait_subjects_need_the_predicate_and_a_taxon(tmp_path: Path) -> None:
    edges = tmp_path / "edges.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\n"
        "e1\tNCBITaxon:1\tMETPO:2000103\tGO:0030245\n"
        "e2\tNCBITaxon:2\tMETPO:2000006\tGO:0030245\n"
        "e3\tCHEBI:1\tMETPO:2000103\tGO:0030245\n"
        "e4\tNCBITaxon:3\tMETPO:2000103\tGO:0045493\n"
    )
    assert direct_trait_subjects(edges, "METPO:2000103|GO:0030245") == {"NCBITaxon:1"}
    with pytest.raises(ValueError, match="no METPO:2000103 -> GO:9999999 edges"):
        direct_trait_subjects(edges, "METPO:2000103|GO:9999999")


def test_net_decreased_needs_more_decreased_than_increased() -> None:
    """A tie between decreased and increased studies is not net decreased."""
    rows = [
        ("t1", "dec", "a"),
        ("t1", "dec", "b"),
        ("t1", "inc", "c"),
        ("t1", "inc", "d"),
        ("t2", "dec", "a"),
        ("t2", "dec", "b"),
        ("t2", "inc", "c"),
    ]
    taxa = pd.DataFrame(rows, columns=["taxon_id", "dir", "study_id"]).assign(
        condition="Parkinson's disease", compartment=FECES
    )
    assert net_decreased(taxa, {"t1", "t2"}, "study_id")["PD"] == {"t2"}


# --------------------------------------------------------------------------- review gaps


def test_stratified_exact_ignores_uninformative_strata_beside_informative_ones() -> None:
    """An uninformative stratum's cell must not join the observed total (#632)."""
    informative = [[2, 1], [1, 2]]
    alone = stratified_exact([informative])
    mixed = stratified_exact([informative, [[1, 0], [0, 0]]])
    for key in ("p", "conditional_or", "mantel_haenszel_or"):
        assert mixed[key] == pytest.approx(alone[key])


def test_net_decreased_flags_need_more_decreased_than_increased_everywhere() -> None:
    row = {f"{short}_{direction}_studies": 3 for short in SHORT for direction in ("dec", "inc")}
    assert not strictly_net_decreased(row, "studies")
    row.update({f"{short}_dec_studies": 4 for short in SHORT})
    assert strictly_net_decreased(row, "studies")


def test_u6_units_number_the_studies_in_order() -> None:
    """signature_null relies on each study's units being contiguous (#632)."""
    studies = ["bsdb:9", "bsdb:1", "bsdb:5", "bsdb:1"]
    conditions = [
        "Parkinson's disease",
        "HIV infection",
        "Colorectal cancer",
        "Parkinson's disease",
    ]
    profiles = pd.DataFrame(
        {"study_id": studies, "condition": conditions, "compartment": FECES, "dir": "dec"}
    )
    taxa = profiles.assign(taxon_id="t1")
    _units, unit_study, _counts, frame = u6_units(profiles, taxa, {"t1"})
    assert list(frame["study_id"]) == sorted(frame["study_id"])
    assert list(unit_study) == sorted(unit_study)
    for numbering in ([1, 0, 1], [0, 2, 2]):
        # Unsorted, or sorted with a gap in the study numbers
        with pytest.raises(ValueError, match="sorted"):
            signature_null(
                np.random.default_rng(0),
                np.ones((3, 1), dtype=bool),
                np.array(numbering, dtype=np.int64),
                dict.fromkeys(SHORT, 0),
                1,
            )


def _case_control_corpus(root: Path) -> Corpus:
    """The fixture corpus with one disease-vs-healthy contrast of each kind the rules tell
    apart (#675):

    - bsdb:1-1, cases against screened controls (Tier 2);
    - bsdb:1-2, untreated cases against hospital controls (Tier 3);
    - bsdb:2-1, tumour against adjacent tissue of the same patients;
    - bsdb:2-2, cases named as treated against Tier 3 controls;
    - bsdb:2-3, a comorbidity contrast.

    bsdb:1-3 is a treatment contrast, as real manifests mark them.
    """
    corpus = load_corpus(_input_root(root), "experiment")
    manifest = corpus.manifest.copy()
    tier3 = TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    manifest.loc[("bsdb:1", "bsdb:1-1"), "control_tier"] = TIER_2_CONDITION_NEGATIVE_SCREENED
    manifest.loc[("bsdb:1", "bsdb:1-2"), "control_tier"] = tier3
    manifest.loc[("bsdb:1", "bsdb:1-3"), ["contrast_type", "control_tier"]] = [
        "treatment_exposure",
        "not_applicable",
    ]
    manifest.loc[("bsdb:2", "bsdb:2-1"), ["contrast_type", "control_tier"]] = [
        PAIRED_TISSUE_AUTOLOGOUS,
        tier3,
    ]
    manifest.loc[("bsdb:2", "bsdb:2-2"), ["control_tier", "flags"]] = [
        tier3,
        "no_group_0_definition_in_export,treated_disease_group",
    ]
    manifest.loc[("bsdb:2", "bsdb:2-3"), ["contrast_type", "control_tier"]] = [
        COMORBIDITY_WITHIN_DISEASE,
        tier3,
    ]
    return replace(corpus, manifest=manifest)


def test_case_control_rules_keep_patients_against_controls(tmp_path: Path) -> None:
    """The primary rule leaves out tissue and comorbidity contrasts; each other rule
    changes one thing, and no two keep the same experiments (#580, #675)."""
    corpus = _case_control_corpus(tmp_path)

    def kept(rule: str) -> set[str]:
        profiles, _ = population(replace(corpus, case_control=rule), ORIENTED)
        return set(profiles["experiment_id"])

    subject = {"bsdb:1-1", "bsdb:1-2", "bsdb:2-2"}
    assert kept("subject_controls") == subject
    assert kept("cases_not_named_treated") == {"bsdb:1-1", "bsdb:1-2"}
    assert kept("control_tiers_1_2") == {"bsdb:1-1"}
    assert kept("with_paired_tissue") == subject | {"bsdb:2-1"}
    assert kept("with_comorbidity") == subject | {"bsdb:2-3"}
    with pytest.raises(ValueError, match="unknown case-control rule 'tissue'"):
        population(replace(corpus, case_control="tissue"), ORIENTED)


def test_case_control_inventory_names_the_rules_that_count_each_experiment(
    tmp_path: Path,
) -> None:
    """Only disease-vs-healthy experiments that some rule sets aside are listed, so none is
    dropped silently; the treatment contrast bsdb:1-3 is not one (#675)."""
    inventory = case_control_inventory(_case_control_corpus(tmp_path))
    every = ["subject_controls", "cases_not_named_treated", "with_paired_tissue"]
    assert [(row["experiment_id"], row["counted_under"]) for row in inventory] == [
        ("bsdb:2-1", ["with_paired_tissue"]),
        ("bsdb:2-2", ["subject_controls", "with_paired_tissue", "with_comorbidity"]),
        ("bsdb:2-3", ["with_comorbidity"]),
        ("bsdb:1-2", [*every, "with_comorbidity"]),
    ]
    assert inventory[0]["contrast_type"] == PAIRED_TISSUE_AUTOLOGOUS
    assert inventory[0]["compartment"] == NON_FECAL_GI


def test_case_control_inventory_needs_every_experiment_in_the_manifest(tmp_path: Path) -> None:
    corpus = _case_control_corpus(tmp_path)
    manifest = corpus.manifest.drop(index=("bsdb:2", "bsdb:2-2"))
    with pytest.raises(ValueError, match="1 experiments have no contrast-manifest row"):
        case_control_inventory(replace(corpus, manifest=manifest))


def test_each_scenario_changes_the_labels_or_the_case_control_rule(tmp_path: Path) -> None:
    """Scenarios start from the primary corpus and change one definition (#561, #566, #580)."""
    corpus = load_corpus(_input_root(tmp_path), "experiment")
    scenarios = {slug: scenario for slug, _label, scenario in scenario_corpora(corpus)}
    assert scenarios["exact_labels"] is corpus
    assert scenarios["synonyms_merged"].condition_labels == "synonyms"
    assert {slug: scenario.case_control for slug, scenario in scenarios.items()} == {
        "exact_labels": "subject_controls",
        "synonyms_merged": "subject_controls",
        "named_treated_cases_excluded": "cases_not_named_treated",
        "control_tiers_1_2": "control_tiers_1_2",
        "paired_tissue_counted": "with_paired_tissue",
        "comorbidity_counted": "with_comorbidity",
        "without_site_miscurations": "subject_controls",
    }
    with pytest.raises(ValueError, match="primary case-control rule"):
        scenario_corpora(replace(corpus, case_control="cases_not_named_treated"))


def test_case_control_scenarios_build_their_overlap_on_disease_vs_healthy_contrasts(
    tmp_path: Path,
) -> None:
    """A case-control rule chooses among disease-vs-healthy contrasts, so its scenario's
    overlap, U3 and U6 are built on them. Setting aside the treated PD stool experiment, the
    only one reporting t3 decreased, takes t3 out of the PD set, the U3 universe and the U6
    pool. Built on all contrasts, where the rule is not applied, the row repeated the
    primary one (#769, #792)."""
    corpus = load_corpus(_input_root(tmp_path), "experiment")
    manifest = corpus.manifest.copy()
    # bsdb:1-1 compares PD cases named as treated with Tier 3 controls; bsdb:1-3 compares
    # people with HIV with controls, so every condition has a set on those contrasts.
    manifest.loc[("bsdb:1", "bsdb:1-1"), ["control_tier", "flags"]] = [
        TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
        "treated_disease_group",
    ]
    manifest.loc[("bsdb:1", "bsdb:1-3"), ["contrast_class", "orientation"]] = [
        "disease_vs_healthy",
        "native",
    ]
    taxa = corpus.taxa
    t3 = taxa[taxa["profile_id"].eq("bsdb:1|bsdb:1-1|dec") & taxa["taxon_id"].eq("t1")]
    taxa = pd.concat([taxa, t3.assign(taxon_id="t3")], ignore_index=True)
    corpus = replace(corpus, manifest=manifest, taxa=taxa)

    _, overlap, _ = scenario_tables(corpus, {"t1", "t3"}, permutations=20, seed=0)
    rows = overlap.groupby("scenario", sort=False)
    for slug in ("named_treated_cases_excluded", "control_tiers_1_2"):
        assert rows.get_group(slug)[["n_PD", "u3_N"]].to_numpy().tolist() == [[1, 1]]
    assert rows.get_group("exact_labels")[["n_PD", "u3_N"]].to_numpy().tolist() == [
        [2, 2],
        [2, 2],
    ]
    both = (ALL_CONTRASTS, ORIENTED)
    assert rows["contrast_population"].agg(tuple).to_dict() == {
        "exact_labels": both,
        "synonyms_merged": both,
        "named_treated_cases_excluded": (ORIENTED,),
        "control_tiers_1_2": (ORIENTED,),
        "paired_tissue_counted": (ORIENTED,),
        "comorbidity_counted": (ORIENTED,),
        "without_site_miscurations": both,
    }

    # U6 draws whole studies from the same disease-vs-healthy contrasts as the sets.
    cellulose = {"t1", "t3"}
    scenario = replace(corpus, case_control="cases_not_named_treated")
    profiles, oriented = population(scenario, ORIENTED)
    observed = set_summary(condition_sets(oriented, cellulose))
    units, unit_study, counts, _ = u6_units(profiles, oriented, cellulose)
    null = signature_null(np.random.default_rng(0), units, unit_study, counts, 20).summary(observed)
    mean, core = null["mean_pairwise_jaccard_distance"], null["three_way_core"]
    row = rows.get_group("named_treated_cases_excluded").iloc[0]
    assert [row[f"u6_units_{short}"] for short in SHORT] == [counts[short] for short in SHORT]
    assert [
        row["u6_expected_mean_jaccard_distance"],
        row["u6_mean_jaccard_distance_p"],
        row["u6_expected_core"],
        row["u6_core_p"],
    ] == pytest.approx([mean["expected"], mean["p"], core["expected"], core["p"]])
    # All contrasts still hold t3, so a pool drawn from them would give another U6.
    everything = u6_units(*population(scenario, ALL_CONTRASTS), cellulose)
    other = signature_null(np.random.default_rng(0), *everything[:3], 20).summary(observed)
    assert other["mean_pairwise_jaccard_distance"]["expected"] != pytest.approx(mean["expected"])


def test_case_control_scenarios_test_direction_and_interaction_on_disease_vs_healthy_contrasts(
    tmp_path: Path,
) -> None:
    """The all-contrasts population does not apply a case-control rule, so a case-control
    scenario's feces direction and Feces-vs-GI rows there repeated the primary rows. Like
    its overlap (#769), they are tested on disease-vs-healthy contrasts only, where the
    rows kept apply the rule, and the three condition-definition tables share one layout
    (#781, #831)."""
    corpus = load_corpus(_input_root(tmp_path), "experiment")
    manifest = corpus.manifest.copy()
    # bsdb:1-1 compares PD cases named as treated with Tier 3 controls; bsdb:1-3 compares
    # people with HIV with controls, so every condition has a set on those contrasts.
    manifest.loc[("bsdb:1", "bsdb:1-1"), ["control_tier", "flags"]] = [
        TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
        "treated_disease_group",
    ]
    manifest.loc[("bsdb:1", "bsdb:1-3"), ["contrast_class", "orientation"]] = [
        "disease_vs_healthy",
        "native",
    ]
    corpus = replace(corpus, manifest=manifest)

    tables = scenario_tables(corpus, {"t1"}, permutations=20, seed=0)
    both = (ALL_CONTRASTS, ORIENTED)
    for table in tables:
        populations = table.groupby("scenario", sort=False)["contrast_population"].unique()
        assert populations.map(tuple).to_dict() == {
            "exact_labels": both,
            "synonyms_merged": both,
            "named_treated_cases_excluded": (ORIENTED,),
            "control_tiers_1_2": (ORIENTED,),
            "paired_tissue_counted": (ORIENTED,),
            "comorbidity_counted": (ORIENTED,),
            "without_site_miscurations": both,
        }
    # PD, HIV, CRC and the pooled three, by experiment and by study, on each population.
    direction, _, interaction = tables
    for table in (direction, interaction):
        assert table.groupby(["scenario", "contrast_population"]).size().eq(8).all()

    # The rows kept apply the scenario's rule: setting aside cases named as treated, or
    # Tier 3 controls, drops bsdb:1-1, PD's decreased feces experiment, from both tables.
    applied = ["exact_labels", "named_treated_cases_excluded", "control_tiers_1_2"]
    for table, counts in ((direction, ["down", "up"]), (interaction, ["a_down", "a_up"])):
        rows = table[
            table["contrast_population"].eq(ORIENTED)
            & table["short_condition"].eq("PD")
            & table["unit"].eq("experiment")
        ].set_index("scenario")
        assert rows.loc[applied, counts].to_numpy().tolist() == [[1, 1], [0, 1], [0, 1]]

    # The rows dropped could only repeat the primary ones: the all-contrasts population
    # does not apply the rule and keeps every profile.
    tiers = replace(corpus, case_control="control_tiers_1_2")
    assert population(tiers, ALL_CONTRASTS)[0].equals(population(corpus, ALL_CONTRASTS)[0])


def test_direction_tests_ask_whether_decreased_signatures_carry_more() -> None:
    """One-sided tests point at the decreased side, and the minimum p is two-sided (#632)."""
    rows = []
    for index in range(8):
        study, experiment = f"bsdb:{index}", f"bsdb:{index}-1"
        rows.append((study, experiment, "dec", 0.5 + index / 100, 1.0))
        rows.append((study, experiment, "inc", 0.1, -1.0))
    profiles = pd.DataFrame(rows, columns=["study_id", "experiment_id", "dir", "frac", "z"])
    taxa = pd.DataFrame({"dir": ["dec"], "taxon_id": ["t1"]})
    row = direction_row(profiles, taxa, {"t1"}, share=0.6)
    assert row["paired_exp_down_dec_gt_inc"] == 8
    assert row["mwu_p_dec_gt_inc_one_sided"] < 0.01
    assert row["wilcoxon_p_down_one_sided"] < 0.01
    assert row["min_attainable_binom_p_two_sided"] == pytest.approx(2 * 0.5**8)


def test_binomial_power_uses_the_two_sided_rejection_region() -> None:
    expected = stats.binom.cdf(1, 10, 0.1) + stats.binom.sf(8, 10, 0.1)
    assert binomial_power(10, 0.1) == pytest.approx(expected)


def test_mean_pairwise_distance_is_a_mean() -> None:
    sets = {"PD": {"a"}, "HIV": {"a", "b"}, "CRC": {"a", "b", "c", "d"}}
    distances = [1 - 1 / 2, 1 - 1 / 4, 1 - 2 / 4]
    assert set_summary(sets)["mean_pairwise_jaccard_distance"] == pytest.approx(sum(distances) / 3)


def test_interaction_wording_follows_p_unit_and_population() -> None:
    """A row that is not significant, or not all contrasts, is not described as one (#631)."""
    significant = {
        "fisher_p_two_sided": 0.001,
        "conditional_or_ci95_low": 2.9,
        "conditional_or_ci95_high": 81.0,
        "unit": "experiment",
        "contrast_population": ALL_CONTRASTS,
    }
    text = interaction_interpretation(significant)
    assert text.startswith("Direction differs between stool and GI mucosa/lumen samples")
    assert "experiments of one study are treated as independent" in text
    assert "do not compare patients with controls" in text
    healthy = {
        "fisher_p_two_sided": 0.13,
        "conditional_or_ci95_low": 0.61,
        "conditional_or_ci95_high": 65.5,
        "unit": "study",
        "contrast_population": ORIENTED,
    }
    text = interaction_interpretation(healthy, shared=2)
    assert text.startswith("No difference in direction")
    # Study rows say how many studies sampled both compartments (#631).
    assert "; 2 studies sampled both compartments;" in text
    assert "; 1 study sampled both compartments;" in interaction_interpretation(healthy, shared=1)
    assert "sampled both" not in interaction_interpretation(healthy)
    assert "do not compare patients with controls" not in text
    assert text.endswith("not a test of anatomical specificity.")
    # The exact test and its interval can disagree near p = 0.05.
    borderline = {**healthy, "fisher_p_two_sided": 0.03, "conditional_or_ci95_low": 0.93}
    assert "the 95% interval of the odds ratio includes 1" in interaction_interpretation(borderline)
    # Boundaries: p = 0.05 is not below 0.05, and an interval ending at 1 includes it.
    at_alpha = {**healthy, "fisher_p_two_sided": 0.05}
    assert interaction_interpretation(at_alpha).startswith("No difference in direction")
    touching = {**borderline, "conditional_or_ci95_low": 1.0}
    assert "includes 1" in interaction_interpretation(touching)
    untestable = {**healthy, "fisher_p_two_sided": math.nan}
    assert interaction_interpretation(untestable).startswith("Not testable")


def test_interaction_tables_put_feces_first(tmp_path: Path) -> None:
    """Row a of each Feces-versus-other table counts the Feces experiments (#632)."""
    corpus = load_corpus(_input_root(tmp_path), "experiment")
    profiles = corpus.profiles.copy()
    # CRC: the feces experiment goes down, both colon experiments go up.
    fractions = {"bsdb:2-1": (0.6, 0.2), "bsdb:2-2": (0.6, 0.2), "bsdb:2-3": (0.2, 0.6)}
    for experiment, (inc, dec) in fractions.items():
        profiles.loc[f"bsdb:2|{experiment}|inc", "frac"] = inc
        profiles.loc[f"bsdb:2|{experiment}|dec", "frac"] = dec
    table = interaction_table(replace(corpus, profiles=profiles))
    row = table[
        table["contrast_population"].eq(ALL_CONTRASTS)
        & table["short_condition"].eq("CRC")
        & table["unit"].eq("experiment")
        & table["compartment_b"].eq(NON_FECAL_GI)
    ].iloc[0]
    assert (row["a_down"], row["a_up"], row["b_down"], row["b_up"]) == (1, 0, 0, 2)


def test_infinite_values_are_not_written_as_missing() -> None:
    """+inf and NaN were both written as null (#631)."""
    values = {"a": math.inf, "b": -math.inf, "c": math.nan, "d": 1.23456789}
    assert clean(values) == {"a": "inf", "b": "-inf", "c": None, "d": 1.23457}


def test_synonym_mode_merges_labels_when_loading(tmp_path: Path) -> None:
    """--condition-labels synonyms counts a synonym under its disease (#632)."""
    paths = _input_root(tmp_path)
    experiments = pd.read_csv(paths["experiments"], sep="\t")
    experiments["condition"] = experiments["condition"].replace("HIV infection", "HIV-1 infection")
    experiments.to_csv(paths["experiments"], sep="\t", index=False)
    exact = load_corpus(paths, "experiment").profiles
    assert exact.loc["bsdb:1|bsdb:1-3|dec", "condition"] == "HIV-1 infection"
    merged = load_corpus(paths, "experiment", "synonyms").profiles
    assert merged.loc["bsdb:1|bsdb:1-3|dec", "condition"] == "HIV infection"


def test_no_msm_population_drops_msm_mentions(tmp_path: Path) -> None:
    """msm_mentioned checks group names, group_1_definition and description (#595)."""
    corpus = load_corpus(_input_root(tmp_path), "experiment")
    manifest = corpus.manifest.copy()
    manifest.loc[("bsdb:1", "bsdb:1-1"), "description"] = "Men who have sex with men (MSM)"
    manifest.loc[("bsdb:2", "bsdb:2-1"), "group_1_definition"] = "Gay and bisexual men"
    assert list(msm_mentioned(manifest)) == [True, False, False, True, False, False]
    corpus = replace(corpus, manifest=manifest)
    oriented, _ = population(corpus, ORIENTED)
    no_msm, taxa = population(corpus, NO_MSM)
    assert {"bsdb:1|bsdb:1-1|dec", "bsdb:2|bsdb:2-1|dec"} <= set(oriented.index)
    assert {"bsdb:1|bsdb:1-1|dec", "bsdb:2|bsdb:2-1|dec"}.isdisjoint(no_msm.index)
    assert set(taxa["profile_id"]) <= set(no_msm.index)


def test_without_studies_drops_profiles_and_taxa_of_dropped_studies(tmp_path: Path) -> None:
    """without_studies removes every profile and taxon of the excluded studies (#595)."""
    corpus = load_corpus(_input_root(tmp_path), "experiment")
    trimmed = without_studies(corpus, ["bsdb:1"])
    assert set(trimmed.profiles["study_id"]) == {"bsdb:2"}
    assert set(trimmed.taxa["study_id"]) == {"bsdb:2"}


def test_catalog_table_records_rolled_up_counts_direct_edges_and_core_class() -> None:
    """catalog_table records rolled-up counts, direct edges and net-dec (#595)."""
    matrix = pd.DataFrame(
        {
            "taxon_name": ["Roseburia", "Lachnospiraceae", "Bacteroides"],
            "taxonomic_rank": ["genus", "family", "genus"],
            TRAIT_ID: [3, 17, 1],
        },
        index=pd.Index(["t1", "t2", "t3"], name="taxon_id"),
    )
    rows = [
        ("p1", "s1", "Parkinson's disease", FECES, "dec", "t1"),
        ("p2", "s1", "Parkinson's disease", FECES, "dec", "t2"),
        ("p3", "s2", "HIV infection", FECES, "dec", "t1"),
        ("p4", "s2", "HIV infection", FECES, "dec", "t2"),
        ("p5", "s3", "Colorectal cancer", FECES, "dec", "t1"),
        ("p6", "s3", "Colorectal cancer", FECES, "dec", "t2"),
        ("p7", "s4", "Colorectal cancer", FECES, "inc", "t2"),
        ("p8", "s5", "Parkinson's disease", FECES, "dec", "t3"),
    ]
    taxa = pd.DataFrame(
        rows, columns=["profile_id", "study_id", "condition", "compartment", "dir", "taxon_id"]
    )
    sets = {"PD": {"t1", "t2", "t3"}, "HIV": {"t1", "t2"}, "CRC": {"t1", "t2"}}
    catalog = catalog_table(taxa, matrix, sets, direct={"t3"}).set_index("taxon_id")
    assert list(catalog.index) == ["t2", "t1", "t3"]
    assert catalog.loc["t1", "matrix_trait_bearing_taxa_at_or_below"] == 3
    assert catalog.loc["t2", "matrix_trait_bearing_taxa_at_or_below"] == 17
    assert bool(catalog.loc["t1", "is_strictly_net_decreased_all3"]) is True
    assert bool(catalog.loc["t2", "is_strictly_net_decreased_all3"]) is False
    assert catalog.loc["t1", "core_class"] == "3-way core, genus or species rank"
    assert catalog.loc["t2", "core_class"] == "3-way core, family-to-phylum rank"
    assert catalog.loc["t3", "core_class"] == "not in 3-way core"
    assert bool(catalog.loc["t3", "has_direct_trait_edge"]) is True


def test_benchmark_table_filters_to_ten_per_condition_and_ranks_non_media() -> None:
    """benchmark_table keeps traits with >=10 dec taxa per condition and ranks non-media (#595)."""
    taxon_ids = [f"t{i}" for i in range(12)]
    matrix = pd.DataFrame(
        {
            "taxon_name": taxon_ids,
            "taxonomic_rank": ["genus"] * 12,
            "trait_high": [1] * 12,
            "mediadive.medium:1": [1] * 10 + [0, 0],
            "trait_low": [1] * 10 + [0, 0],
            "trait_too_small": [1] * 9 + [0, 0, 0],
        },
        index=pd.Index(taxon_ids, name="taxon_id"),
    )
    # PD, HIV, CRC each report 10 taxa decreased for trait_low and mediadive.medium:1;
    # CRC also reports t10 and t11 decreased for trait_high.
    rows = []
    for short, members in (
        ("Parkinson's disease", taxon_ids[:10]),
        ("HIV infection", taxon_ids[:10]),
        ("Colorectal cancer", taxon_ids[:12]),
    ):
        for taxon in members:
            rows.append((f"{short}-{taxon}", short, FECES, "dec", taxon))
    taxa = pd.DataFrame(rows, columns=["profile_id", "condition", "compartment", "dir", "taxon_id"])
    labels = pd.DataFrame(
        {
            "object_name": ["high overlap", "medium 1", "perfect overlap", "small"],
            "predicate_name": ["p", "p", "p", "p"],
        },
        index=pd.Index(
            ["trait_high", "mediadive.medium:1", "trait_low", "trait_too_small"], name="column_id"
        ),
    )
    bench = benchmark_table(taxa, matrix, labels)
    # trait_too_small has only 9 taxa per condition and is excluded.
    # mediadive.medium:1 and trait_low have identical 10-taxon sets across all 3 conditions
    # (d_J = 0); trait_high has 10 in PD/HIV and 12 in CRC (d_J > 0).
    assert list(bench["trait_id"]) == ["mediadive.medium:1", "trait_low", "trait_high"]
    assert list(bench["rank_by_convergence"]) == [1, 2, 3]
    assert list(bench["rank_excluding_growth_media"].dropna()) == [1, 2]
    assert pd.isna(
        bench.loc[
            bench["trait_id"].eq("mediadive.medium:1"), "percentile_excluding_growth_media"
        ].iloc[0]
    )
    assert bench.loc[bench["trait_id"].eq("trait_low"), "percentile_excluding_growth_media"].iloc[
        0
    ] == pytest.approx(0.5)


def test_tier1_2_catalog_scopes_taxa_to_tier1_2_oriented_contrasts(tmp_path: Path) -> None:
    """The Tier 1-2 catalog includes only taxa from Tier 1-2 disease-vs-healthy contrasts."""
    corpus = load_corpus(_input_root(tmp_path), "experiment")
    manifest = corpus.manifest.copy()
    manifest.loc[("bsdb:1", "bsdb:1-1"), "control_tier"] = TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    manifest.loc[("bsdb:1", "bsdb:1-3"), ["contrast_class", "orientation"]] = [
        "disease_vs_healthy",
        "native",
    ]
    taxa = corpus.taxa
    t3 = taxa[taxa["profile_id"].eq("bsdb:1|bsdb:1-1|dec") & taxa["taxon_id"].eq("t1")]
    taxa = pd.concat([taxa, t3.assign(taxon_id="t3")], ignore_index=True)
    corpus = replace(corpus, manifest=manifest, taxa=taxa)

    _all_p, all_taxa = population(corpus, ALL_CONTRASTS)
    all_sets = condition_sets(all_taxa, {"t1", "t3"})
    tier12_corpus = replace(corpus, case_control="control_tiers_1_2")
    _t12_p, t12_taxa = population(tier12_corpus, ORIENTED)
    t12_sets = condition_sets(t12_taxa, {"t1", "t3"})

    assert "t3" in all_sets["PD"]
    assert "t3" not in t12_sets["PD"]
