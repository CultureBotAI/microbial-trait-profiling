"""cellulose_decr overlap, null-model and direction statistics (#552, #557)."""

from __future__ import annotations

import itertools
import json
import math
import shutil
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats
from scipy.stats import false_discovery_control
from scipy.stats.contingency import odds_ratio

from mtp.cellulose_decr import figures
from mtp.cellulose_decr.analysis import (
    ANALYSIS_DIR,
    FECES,
    INPUTS,
    NON_FECAL_GI,
    ORAL,
    ORIENTED,
    OTHER,
    SHORT,
    STRONG,
    TRAIT_ID,
    NullDraws,
    benjamini_hochberg,
    check_input_drift,
    choose,
    compartment_of,
    direct_trait_subjects,
    hypergeom_pair,
    label_corrections,
    load_corpus,
    merge_synonyms,
    most_common_label,
    net_decreased,
    paired_deltas,
    population,
    sign_tests,
    signature_null,
    stratified_exact,
    study_votes,
    three_way_exact,
)
from mtp.contrast_orientation import CLASSIFIER_VERSION

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
    assert three_way_exact(7, (3, 4, 2), k)["p_3way_exact"] == pytest.approx(expected)


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


# --------------------------------------------------------------------------- figures


def test_figures_draw_from_the_committed_tables(tmp_path: Path) -> None:
    """Every figure renders from the committed tables, with a slide caption (#556).

    The file names describe the figures rather than claims the analysis rejects (#570),
    and the captions quote the tables.
    """
    for name in figures.FIGURE_INPUTS:
        shutil.copy(ANALYSIS_DIR / name, tmp_path / name)
    commit = "0123456789" * 4
    assert figures.main(["--analysis-dir", str(tmp_path), "--code-commit", commit]) == 0
    provenance = json.loads((tmp_path / "figures" / "figures_provenance.json").read_text())
    assert provenance["script"] == "src/mtp/cellulose_decr/figures.py"
    # Decks link the figures' code and tables at this commit (#556).
    assert provenance["code_commit"] == commit
    assert len(provenance["outputs"]) == 7 * 3 + 2 + 1
    for output in provenance["outputs"]:
        assert (tmp_path / output).stat().st_size > 0
    for claim in ("surprising", "convergence", "guild", "specificity", "podium"):
        assert not [output for output in provenance["outputs"] if claim in output], claim

    captions = json.loads((tmp_path / "figures" / "figure_captions.json").read_text())
    svgs = {output for output in provenance["outputs"] if output.endswith(".svg")}
    assert set(captions["figures"]) == svgs
    for caption in captions["figures"].values():
        assert all(caption[field] for field in ("title", "takeaway", "caveat"))
        assert (ANALYSIS_DIR / caption["evidence"]).exists()

    summary = json.loads((ANALYSIS_DIR / "cellulose_decr_summary.json").read_text())
    observed = summary["observed"]
    overlap = captions["figures"]["figures_graphical/fig1_overlap_and_null_models.svg"]
    assert overlap["title"].startswith(
        f"{observed['three_way_core']} of {observed['three_way_union']} "
    )
    assert overlap["title"] == captions["statements"]["overlap"]
    # PD, CRC and the pooled result keep p < 0.05 under every condition definition; HIV
    # does not (SIGNIFICANT_RESULTS_SUMMARY.md, Table A2).
    direction = captions["figures"]["figures_graphical/fig3_direction_by_body_site.svg"]
    assert direction["title"].endswith(
        "in PD, CRC and all three pooled under every condition definition; in HIV only under some"
    )
    assert "HIV 9 : 0 : 2" in direction["takeaway"]
    assert captions["statements"]["fragile"].startswith("the HIV result depends")
    assert captions["statements"]["crc_compartments"] in direction["caveat"]


# --------------------------------------------------------------------------- captions


@pytest.fixture(scope="module")
def results() -> figures.Results:
    return figures.load(ANALYSIS_DIR)


def _with_scenarios(
    results: figures.Results, edit: Callable[[pd.DataFrame], None]
) -> figures.Results:
    table = results.scenarios.copy()
    edit(table)
    return replace(results, scenarios=table)


def _healthy_studies(table: pd.DataFrame, short: str) -> pd.Series:
    return (
        table["contrast_population"].eq(ORIENTED)
        & table["short_condition"].eq(short)
        & table["unit"].eq("study")
    )


def test_direction_caption_follows_the_direction_not_only_p(results: figures.Results) -> None:
    """A condition significant the other way is named as such, not as a decrease (#589)."""

    def swap_crc(table: pd.DataFrame) -> None:
        rows = _healthy_studies(table, "CRC")
        table.loc[rows, ["down", "up"]] = table.loc[rows, ["up", "down"]].to_numpy()

    title = figures.caption_direction(_with_scenarios(results, swap_crc)).title
    assert "the increased signature does in CRC under every definition" in title
    assert "in PD and all three pooled under every condition definition" in title

    def crc_up_under_one(table: pd.DataFrame) -> None:
        swap_crc(table)
        rows = table.index[_healthy_studies(table, "CRC")]
        table.loc[rows[1:], "binom_p_two_sided"] = 0.5

    edited = _with_scenarios(results, crc_up_under_one)
    assert "the increased signature does in CRC only under some" in (
        figures.caption_direction(edited).title
    )
    assert "CRC: p from " in figures.caption_direction(edited).caveat
    fragile = figures.fragile_statement(edited)
    assert fragile.startswith("the HIV and CRC results depend")
    # Quoted beside decreased-direction claims, it must say which way CRC goes (#608).
    assert fragile.endswith(" toward the increased signature)")


def test_mixed_direction_is_named_not_only_its_p_range(results: figures.Results) -> None:
    """Small p-values in opposite directions are not a robust result (#589)."""

    def flip_one_hiv(table: pd.DataFrame) -> None:
        rows = table.index[_healthy_studies(table, "HIV") & table["binom_p_two_sided"].lt(0.05)]
        down, up = table.loc[rows[0], ["down", "up"]]
        table.loc[rows[0], ["down", "up"]] = [up, down]

    edited = _with_scenarios(results, flip_one_hiv)
    assert figures.direction_verdicts(edited)["HIV"].kind == "mixed"
    caption = figures.caption_direction(edited)
    assert "HIV changes direction with the definition" in caption.title
    assert "HIV: p from " in caption.caveat
    assert " in opposite directions across " in caption.caveat
    assert figures.fragile_statement(edited).endswith(" in opposite directions)")


def test_fragile_pooled_result_starts_its_sentence_and_is_named(
    results: figures.Results,
) -> None:
    def pooled_under_one(table: pd.DataFrame) -> None:
        rows = table.index[_healthy_studies(table, "Pooled")]
        table.loc[rows[1:], "binom_p_two_sided"] = 0.5

    edited = _with_scenarios(results, pooled_under_one)
    assert "All three pooled: p from " in figures.caption_direction(edited).caveat
    statement = figures.fragile_statement(edited)
    assert statement.startswith("the HIV and pooled results depend on the condition definition")
    assert "(study-level p HIV " in statement
    assert "; pooled " in statement


def test_direction_caption_without_significance_or_data(results: figures.Results) -> None:
    def no_significance(table: pd.DataFrame) -> None:
        table["binom_p_two_sided"] = 0.5

    title = figures.caption_direction(_with_scenarios(results, no_significance)).title
    assert title.endswith("in no condition at p < 0.05")
    assert " not " not in title

    def no_hiv(table: pd.DataFrame) -> None:
        table.drop(table.index[_healthy_studies(table, "HIV")], inplace=True)

    title = figures.caption_direction(_with_scenarios(results, no_hiv)).title
    assert title.endswith("; no data for HIV")
    assert figures.fragile_statement(_with_scenarios(results, no_hiv)) == ""


@pytest.mark.parametrize(
    ("mean_p", "core_p", "ending"),
    [
        (0.19, 0.73, "random sets of whole stool signatures overlap about as much"),
        (0.01, 0.02, "they overlap more than random sets of whole stool signatures"),
        (0.99, 0.98, "they overlap less than random sets of whole stool signatures"),
        (0.04, 0.73, "mean-distance p = 0.040 and core p = 0.73"),
    ],
)
def test_overlap_statement_reads_both_u6_endpoints(
    results: figures.Results, mean_p: float, core_p: float, ending: str
) -> None:
    summary = json.loads(json.dumps(results.summary))
    summary["null_models"]["U6"]["mean_pairwise_jaccard_distance"]["p"] = mean_p
    summary["null_models"]["U6"]["three_way_core"]["p"] = core_p
    assert figures.overlap_statement(replace(results, summary=summary)).endswith(ending)


def test_shared_taxa_caption_counts_direct_edges(results: figures.Results) -> None:
    catalog = results.catalog.copy()
    core = catalog.index[catalog["is_3way_core"].astype(bool)]
    catalog.loc[core[0], "has_direct_trait_edge"] = True
    summary = {**results.summary, "core_taxa_with_direct_trait_edge": 1}
    edited = replace(results, catalog=catalog, summary=summary)
    caption = figures.caption_shared_taxa(edited)
    assert caption.title.endswith(", 1 annotated directly")
    assert caption.takeaway.startswith("1 carries its own cellulose annotation edge")
    assert figures.direct_core_text(edited).startswith("1 of the 17 core taxa has its own")


def test_crc_statement_names_the_excluded_experiment_only_when_it_matters(
    results: figures.Results,
) -> None:
    assert "without the tumour-vs-adjacent-tissue experiment" in (
        figures.crc_compartment_statement(results)
    )
    table = results.scenario_interaction.copy()
    primary = figures.interaction_row(results, ORIENTED, "CRC", "study")
    excluded = table["scenario"].eq("crc_within_patient_excluded")
    table.loc[excluded, "fisher_p_two_sided"] = primary["fisher_p_two_sided"]
    statement = figures.crc_compartment_statement(replace(results, scenario_interaction=table))
    assert "without" not in statement
    assert "does not show anatomical specificity" in statement


def test_p_values_keep_two_significant_digits() -> None:
    assert figures.sci(0.0097) == "0.0097"
    assert figures.sci(0.0596) == "0.060"
    assert figures.sci(0.19) == "0.19"
    assert figures.sci(1.2e-5) == "1.2" + figures.TIMES + "10⁻⁵"
