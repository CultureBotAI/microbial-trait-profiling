"""Source-lineage removal and conditional, study-weighted influence summaries.

Matched-corpus Table 1-trait analog sensitivity analysis (not a reproduction of
the published Table 1 baseline).
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from statistics import mean, stdev

from scipy.stats import binomtest, t

from mtp.contribution_inputs import RANKS, Base, Graph, Key, Pair, encoded
from mtp.taxon_contributions import PROFILE_FIELDS, Ledger, Record, profile_record, write_table
from mtp.trait_composition import paired_log2_ratio_change

TIE = 1e-12


@lru_cache(maxsize=20000)
def directional_test(positive: int, negative: int) -> tuple[float, float, float]:
    test = binomtest(positive, positive + negative)
    interval = test.proportion_ci(confidence_level=0.95, method="exact")
    return float(test.pvalue), float(interval.low), float(interval.high)


@lru_cache(maxsize=5000)
def critical_t(degrees: int) -> float:
    return float(t.ppf(0.975, degrees))


def scope_keys(base: Base, scope: Record) -> list[Key]:
    return [
        key
        for key in sorted(base.profiles)
        if all(
            scope[field] == "*" or base.metadata[key[:2]][field] == scope[field]
            for field in ("condition", "body_site")
        )
    ]


def study_means(values: dict[Pair, float]) -> dict[str, float]:
    studies: defaultdict[str, list[float]] = defaultdict(list)
    for pair, value in values.items():
        studies[pair[0]].append(value)
    return {study: mean(items) for study, items in studies.items()}


def directional_counts(values: Iterable[float]) -> tuple[int, int]:
    """Count signs with the shared numerical tie tolerance."""
    items = list(values)
    return sum(value > TIE for value in items), sum(value < -TIE for value in items)


def describe(values: dict[Pair, float]) -> Record:
    """The t interval targets a mean; the exact interval targets positive share."""
    studies = study_means(values)
    items = list(studies.values())
    positive, negative = directional_counts(items)
    result: Record = {
        "pairs": len(values),
        "studies": len(items),
        "effect": "",
        "mean_ci_low": "",
        "mean_ci_high": "",
        "positive": positive,
        "negative": negative,
        "ties": len(items) - positive - negative,
        "share_ci_low": "",
        "share_ci_high": "",
        "positive_share": "",
        "p": "",
        "status": "no_eligible_pairs",
    }
    if items:
        result.update(effect=mean(items), status="insufficient_studies")
    if len(items) >= 2:
        radius = critical_t(len(items) - 1) * stdev(items) / math.sqrt(len(items))
        result.update(
            mean_ci_low=mean(items) - radius,
            mean_ci_high=mean(items) + radius,
            status="conditional_descriptive",
        )
    if positive + negative:
        p, low, high = directional_test(positive, negative)
        result.update(
            p=p,
            share_ci_low=low,
            share_ci_high=high,
            positive_share=positive / (positive + negative),
        )
    elif items:
        result["status"] = "all_tied"
    return result


def adjust(pvalues: list[float | None], by: bool = False) -> list[float | None]:
    """Keep the full predeclared family, including untestable members."""
    size = len(pvalues)
    factor = sum(1 / i for i in range(1, size + 1)) if by else 1.0
    ranked = sorted((p, i) for i, p in enumerate(pvalues) if p is not None)
    result: list[float | None] = [None] * size
    ceiling = 1.0
    for rank in range(len(ranked), 0, -1):
        p, index = ranked[rank - 1]
        ceiling = min(ceiling, p * size * factor / rank)
        result[index] = ceiling
    return result


def pair_values(
    profiles: dict[Key, dict[str, int]],
    keys: list[Key],
    trait: str,
    denominator: str = "",
    native_zero: bool = False,
    totals: dict[Key, int] | None = None,
) -> dict[Pair, float]:
    values = {}
    for key in keys:
        if key[2] != "inc":
            continue
        pair = key[:2]
        inc, dec = profiles[key], profiles[(*pair, "dec")]
        itotal, dtotal = (
            (totals[key], totals[(*pair, "dec")])
            if totals is not None
            else (sum(inc.values()), sum(dec.values()))
        )
        if denominator:
            contrast = paired_log2_ratio_change(
                inc[trait], inc[denominator], dec[trait], dec[denominator]
            )
            if contrast is not None:
                values[pair] = contrast
        elif itotal and dtotal:
            values[pair] = inc[trait] / itotal - dec[trait] / dtotal
        elif native_zero:
            values[pair] = (inc[trait] / itotal if itotal else 0) - (
                dec[trait] / dtotal if dtotal else 0
            )
    return values


def removal(
    base: Base, graph: Graph, ledger: Ledger, lineage: str
) -> tuple[dict[Key, dict[str, int]], dict[str, Counter[str]], dict[str, set[str]]]:
    """Remove source units through every reported ancestor, across the full denominator."""
    removed: dict[str, Counter[str]] = {}
    ambiguous: dict[str, set[str]] = {}
    lineage_ancestors = graph.ancestors(lineage) - {lineage}
    for target, sources in ledger.sources.items():
        vector: Counter[str] = Counter()
        for source in sources:
            if lineage in graph.ancestors(source):
                vector.update(graph.traits[source])
            elif source in lineage_ancestors:
                ambiguous.setdefault(target, set()).add(source)
        if vector:
            removed[target] = vector
    remaining = {}
    for key, baseline in base.profiles.items():
        subtract: Counter[str] = Counter()
        for target in base.reported[key]:
            subtract.update(removed.get(target, {}))
        vector_dict = {trait: value - subtract[trait] for trait, value in baseline.items()}
        if any(value < 0 for value in vector_dict.values()):
            raise ValueError(
                f"Negative remainder for {key}, {lineage}; reconstruction is inconsistent"
            )
        remaining[key] = vector_dict
    return remaining, removed, ambiguous


def source_occurrences(
    base: Base, ledger: Ledger, keys: list[Key]
) -> tuple[Counter[tuple[str, str]], dict[tuple[str, str], set[str]]]:
    """Aggregate repeated reported targets before traversing their source units."""
    targets: Counter[tuple[str, str]] = Counter()
    target_studies: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    for key in keys:
        for target in base.reported[key]:
            targets[(key[2], target)] += 1
            target_studies[(key[2], target)].add(key[0])
    counts: Counter[tuple[str, str]] = Counter()
    studies: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    for (direction, target), number in targets.items():
        for source in ledger.sources[target]:
            counts[(direction, source)] += number
            studies[(direction, source)].update(target_studies[(direction, target)])
    return counts, dict(studies)


def taxonomic_counts(
    base: Base, graph: Graph, ledger: Ledger, keys: list[Key]
) -> tuple[dict[tuple[str, str, str], Counter[str]], list[Record]]:
    counts: defaultdict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    support: defaultdict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    rank_cache = {
        s: {rank: graph.ranked(s, rank, base.ranks) for rank in RANKS}
        for s in {source for sources in ledger.sources.values() for source in sources}
    }
    six = set(base.traits.values())
    occurrences, study_support = source_occurrences(base, ledger, keys)
    for (direction, source), number in occurrences.items():
        for trait in graph.traits[source] & six:
            for rank in RANKS:
                taxon, status = rank_cache[source][rank]
                identity = taxon or status
                counts[(direction, trait, rank)][identity] += number
                support[(direction, trait, rank, identity)].update(
                    study_support[(direction, source)]
                )
            # Source-record incidence is not an additive partition of count units.
            evidence = graph.assertions.get((source, trait), [])
            identifiers = {a.primary_source for a in evidence if a.primary_source}
            for identity in identifiers or {"unresolved_source"}:
                counts[(direction, trait, "primary_source_incidence")][identity] += number
                support[(direction, trait, "primary_source_incidence", identity)].update(
                    study_support[(direction, source)]
                )
    records = []
    for (direction, trait, rank), vector in sorted(counts.items()):
        total = sum(vector.values())
        for taxon, value in sorted(vector.items()):
            records.append(
                {
                    "direction": direction,
                    "trait_key": trait,
                    "rank": rank,
                    "taxon": taxon,
                    "name": graph.names.get(taxon, base.names.get(taxon, "")),
                    "count": value,
                    "share": value / total,
                    "studies": len(support[(direction, trait, rank, taxon)]),
                    "status": "resolved"
                    if not taxon.startswith(("unresolved", "ambiguous"))
                    else taxon,
                }
            )
    return dict(counts), records


def concentration(
    counts: dict[tuple[str, str, str], Counter[str]], traits: list[str]
) -> list[Record]:
    records = []
    for direction in ("inc", "dec"):
        for trait in traits:
            for rank in (*RANKS, "primary_source_incidence"):
                vector = counts.get((direction, trait, rank), Counter())
                total = sum(vector.values())
                resolved = {
                    k: v for k, v in vector.items() if not k.startswith(("unresolved", "ambiguous"))
                }
                resolved_total = sum(resolved.values())
                records.append(
                    {
                        "direction": direction,
                        "trait_key": trait,
                        "rank": rank,
                        "total": total,
                        "resolved_total": resolved_total,
                        "unresolved_total": total - resolved_total,
                        "resolved_taxa": len(resolved),
                        "hhi_resolved": sum((v / resolved_total) ** 2 for v in resolved.values())
                        if resolved_total
                        else "",
                        "largest_share_resolved": max(resolved.values()) / resolved_total
                        if resolved_total
                        else "",
                    }
                )
    return records


def overlap(counts: dict[tuple[str, str, str], Counter[str]], traits: list[str]) -> list[Record]:
    records = []
    for direction in ("inc", "dec"):
        for rank in RANKS:
            for first, second in combinations(traits, 2):
                left, right = (
                    {
                        k: v
                        for k, v in counts.get((direction, tr, rank), {}).items()
                        if not k.startswith(("unresolved", "ambiguous"))
                    }
                    for tr in (first, second)
                )
                union, shared = left.keys() | right.keys(), left.keys() & right.keys()
                records.append(
                    {
                        "direction": direction,
                        "rank": rank,
                        "first_trait": first,
                        "second_trait": second,
                        "shared_taxa": len(shared),
                        "union_taxa": len(union),
                        "jaccard": len(shared) / len(union) if union else "",
                        "weighted_overlap": sum(
                            min(left[k] / sum(left.values()), right[k] / sum(right.values()))
                            for k in shared
                        )
                        if left and right
                        else "",
                        "status": "resolved_subset_only" if left and right else "untestable",
                    }
                )
    return records


def scenarios(base: Base, graph: Graph, ledger: Ledger, minimum: int) -> list[Record]:
    """Enumerate every supported lineage without using observed effects or cutoffs."""
    source_studies: defaultdict[str, set[str]] = defaultdict(set)
    source_counts: Counter[str] = Counter()
    occurrences, study_support = source_occurrences(base, ledger, list(base.profiles))
    for (direction, source), number in occurrences.items():
        source_studies[source].update(study_support[(direction, source)])
        source_counts[source] += number * len(graph.traits[source])
    lineage_studies: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    lineage_counts: Counter[tuple[str, str]] = Counter()
    for source, studies in source_studies.items():
        for rank in ("genus", "family"):
            lineage, _ = graph.ranked(source, rank, base.ranks)
            if lineage:
                lineage_studies[(rank, lineage)].update(studies)
                lineage_counts[(rank, lineage)] += source_counts[source]
    return [
        {
            "scenario": f"exclude_{rank}_{lineage}",
            "rank": rank,
            "lineage": lineage,
            "name": graph.names.get(lineage, base.names.get(lineage, "")),
            "studies": len(studies),
            "propagated_non_mediadive_count": lineage_counts[(rank, lineage)],
            "status": "run" if len(studies) >= minimum else "insufficient_support",
        }
        for (rank, lineage), studies in sorted(lineage_studies.items())
    ]


def trait_leading_genera(base: Base, graph: Graph, ledger: Ledger, minimum: int) -> dict[str, str]:
    """Select each Table 1 trait's supported leading genus without using effects."""
    occurrences, study_support = source_occurrences(base, ledger, list(base.profiles))
    counts: Counter[tuple[str, str]] = Counter()
    studies: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    six = set(base.traits.values())
    for (direction, source), number in occurrences.items():
        genus, _ = graph.ranked(source, "genus", base.ranks)
        if not genus:
            continue
        for trait in graph.traits[source] & six:
            counts[(trait, genus)] += number
            studies[(trait, genus)].update(study_support[(direction, source)])
    leaders = {}
    for trait in base.traits.values():
        candidates = sorted(
            (
                (lineage, count)
                for (candidate_trait, lineage), count in counts.items()
                if candidate_trait == trait and len(studies[(trait, lineage)]) >= minimum
            ),
            key=lambda item: (-item[1], item[0]),
        )
        if candidates:
            leaders[trait] = candidates[0][0]
    return leaders


def comparison(before: dict[Pair, float], after: dict[Pair, float]) -> Record:
    common = before.keys() & after.keys()
    bm = {k: before[k] for k in common}
    am = {k: after[k] for k in common}
    # Re-average both versions on the SAME experiments before differencing studies.
    delta = {k: after[k] - before[k] for k in common}
    records = {
        label: describe(values)
        for label, values in (
            ("before", before),
            ("after", after),
            ("before_matched", bm),
            ("after_matched", am),
            ("matched_change", delta),
        )
    }
    result = {
        f"{label}_{field}": value
        for label, record in records.items()
        for field, value in record.items()
    }
    result.update(
        lost_pairs=len(before.keys() - after.keys()), gained_pairs=len(after.keys() - before.keys())
    )
    result["eligibility_change"] = (
        records["after"]["effect"]
        - records["before"]["effect"]
        - records["matched_change"]["effect"]
        if before and after and common
        else ""
    )
    result["before_eligibility_component"] = (
        records["before_matched"]["effect"] - records["before"]["effect"]
        if before and common
        else ""
    )
    result["after_eligibility_component"] = (
        records["after"]["effect"] - records["after_matched"]["effect"] if after and common else ""
    )
    return result


def effect_rows(base: Base, remaining: dict[Key, dict[str, int]], keys: list[Key]) -> list[Record]:
    results = []
    before_totals = {key: sum(base.profiles[key].values()) for key in keys}
    after_totals = {key: sum(remaining[key].values()) for key in keys}
    traits = list(base.traits.values())
    valid_before = {
        key[:2]
        for key in keys
        if key[2] == "inc" and before_totals[key] > 0 and before_totals[(*key[:2], "dec")] > 0
    }
    valid_after = {
        key[:2]
        for key in keys
        if key[2] == "inc" and after_totals[key] > 0 and after_totals[(*key[:2], "dec")] > 0
    }
    for representation in ("shared_fraction", "native_noel_zero_fraction", "shared_log2_ratio"):
        features = (
            list(combinations(traits, 2))
            if representation == "shared_log2_ratio"
            else [(trait, "") for trait in traits]
        )
        family = []
        for trait, other in features:
            native_zero = representation == "native_noel_zero_fraction"
            before = pair_values(base.profiles, keys, trait, other, native_zero, before_totals)
            after = pair_values(remaining, keys, trait, other, native_zero, after_totals)
            if native_zero:
                before_cmp = {p: v for p, v in before.items() if p in valid_before}
                after_cmp = {p: v for p, v in after.items() if p in valid_after}
            else:
                before_cmp = before
                after_cmp = after
            row = dict(
                representation=representation,
                trait_key=trait,
                denominator_trait=other,
                unit="log2_ratio_contrast" if other else "fraction_difference",
                correction="BY15" if other else "BH6",
                **comparison(before_cmp, after_cmp),
            )
            family.append(row)
        if representation == "native_noel_zero_fraction":
            for row in family:
                for population in (
                    "before",
                    "after",
                    "before_matched",
                    "after_matched",
                    "matched_change",
                ):
                    if row[f"{population}_p"] == "":
                        # Preserve native p=1 for unscored features, but leave their
                        # undefined directional share/interval missing, with status.
                        row[f"{population}_p"] = 1.0
        for population in ("before", "after", "before_matched", "after_matched", "matched_change"):
            pvalues = [
                float(r[f"{population}_p"]) if r[f"{population}_p"] != "" else None for r in family
            ]
            for row, q in zip(
                family, adjust(pvalues, representation == "shared_log2_ratio"), strict=True
            ):
                row[f"{population}_q"] = q if q is not None else ""
        results.extend(family)
    return results


def pair_manifest(base: Base, remaining: dict[Key, dict[str, int]]) -> Iterator[Record]:
    for key in sorted(base.profiles):
        if key[2] != "inc":
            continue
        pair = key[:2]
        record: Record = {
            "study_id": pair[0],
            "experiment_id": pair[1],
            "condition": base.metadata[pair]["condition"],
            "body_site": base.metadata[pair]["body_site"],
            "orientation_status": "unreviewed_188",
        }
        for version, profiles in (("before", base.profiles), ("after", remaining)):
            for direction in ("inc", "dec"):
                vector = profiles[(*pair, direction)]
                record[f"{version}_{direction}_total"] = sum(vector.values())
                record[f"{version}_{direction}_zero_trait_keys"] = encoded(
                    [trait for trait in base.traits.values() if vector[trait] == 0]
                )
        yield record


def run_influence(
    base: Base, graph: Graph, ledger: Ledger, protocol: Record, output: Path
) -> Record:
    all_scenarios = scenarios(base, graph, ledger, protocol["minimum_lineage_studies"])
    selection_reasons: defaultdict[str, set[str]] = defaultdict(set)
    for rank in ("genus", "family"):
        ranked = sorted(
            (r for r in all_scenarios if r["rank"] == rank and r["status"] == "run"),
            key=lambda r: (-r["propagated_non_mediadive_count"], r["lineage"]),
        )
        for scenario in ranked[: protocol["top_lineages_per_rank"]]:
            selection_reasons[scenario["lineage"]].add(f"overall_top_{rank}")
    for scenario in all_scenarios:
        if scenario["rank"] == "genus" and scenario["name"] in protocol["predicted_genera"]:
            selection_reasons[scenario["lineage"]].add("named_prediction_genus")
            family, _ = graph.ranked(scenario["lineage"], "family", base.ranks)
            if family:
                selection_reasons[family].add("named_prediction_family")
    if protocol["include_trait_leading_genus"]:
        for trait, lineage in trait_leading_genera(
            base, graph, ledger, protocol["minimum_lineage_studies"]
        ).items():
            selection_reasons[lineage].add(f"trait_leader:{trait}")
    for scenario in all_scenarios:
        scenario["selection_reasons"] = encoded(
            sorted(selection_reasons.get(scenario["lineage"], set()))
        )
        if scenario["status"] == "run" and scenario["lineage"] not in selection_reasons:
            scenario["status"] = "outside_rule_declared_panel"
    write_table(
        output / "scenarios.tsv",
        [
            "scenario",
            "rank",
            "lineage",
            "name",
            "studies",
            "propagated_non_mediadive_count",
            "selection_reasons",
            "status",
        ],
        all_scenarios,
    )
    predictions = [
        {
            "prediction": name,
            "matched_lineages": encoded([r["lineage"] for r in all_scenarios if r["name"] == name]),
            "status": "prediction_only_not_assumed_driver",
        }
        for name in protocol["predicted_genera"]
    ]
    write_table(
        output / "predictions.tsv", ["prediction", "matched_lineages", "status"], predictions
    )
    universe: list[Record] = []
    concentrations: list[Record] = []
    overlaps: list[Record] = []
    scopes = {scope["id"]: scope_keys(base, scope) for scope in protocol["scopes"]}
    for scope, keys in scopes.items():
        counts, records = taxonomic_counts(base, graph, ledger, keys)
        universe.extend(dict(scope=scope, **r) for r in records)
        concentrations.extend(
            dict(scope=scope, **r) for r in concentration(counts, list(base.traits.values()))
        )
        overlaps.extend(dict(scope=scope, **r) for r in overlap(counts, list(base.traits.values())))
    for name, records in (
        ("contributing_universe", universe),
        ("concentration", concentrations),
        ("cross_trait_overlap", overlaps),
    ):
        write_table(output / f"{name}.tsv", list(records[0]), records)
    active = [r for r in all_scenarios if r["status"] == "run"]
    effect_count = 0
    pooled_effects: list[Record] = []

    # Sparse reported-taxon removals reconstruct every profile using the ledger links.
    def exports() -> Iterator[Record]:
        nonlocal effect_count
        for index, scenario in enumerate([{"scenario": "baseline", "lineage": ""}, *active]):
            name, lineage = scenario["scenario"], scenario["lineage"]
            print(f"Influence {index}/{len(active)}: {name}", flush=True)
            remaining, removed, ambiguous = (
                removal(base, graph, ledger, lineage) if lineage else (base.profiles, {}, {})
            )
            directory = output / name.replace(":", "_")
            directory.mkdir(exist_ok=True)
            pair_fields = [
                "study_id",
                "experiment_id",
                "condition",
                "body_site",
                "orientation_status",
                *(
                    f"{version}_{direction}_{suffix}"
                    for version in ("before", "after")
                    for direction in ("inc", "dec")
                    for suffix in ("total", "zero_trait_keys")
                ),
            ]
            write_table(
                directory / "pair_eligibility.tsv.gz", pair_fields, pair_manifest(base, remaining)
            )
            write_table(
                directory / "removed_source_counts.tsv.gz",
                ["reported_taxon", "trait_key", "removed_count"],
                (
                    {"reported_taxon": taxon, "trait_key": trait, "removed_count": value}
                    for taxon, vector in sorted(removed.items())
                    for trait, value in sorted(vector.items())
                ),
            )
            write_table(
                directory / "ambiguous_retained_sources.tsv",
                ["reported_taxon", "source_taxon", "status"],
                (
                    {
                        "reported_taxon": taxon,
                        "source_taxon": source,
                        "status": "higher_rank_direct_assertion_retained",
                    }
                    for taxon, sources in sorted(ambiguous.items())
                    for source in sorted(sources)
                ),
            )
            write_table(
                directory / "profile_counts.tsv.gz",
                [
                    *PROFILE_FIELDS,
                    "trait_key",
                    "original_count",
                    "remaining_count",
                    "original_total",
                    "remaining_total",
                    "ambiguity_flag",
                ],
                (
                    {
                        **profile_record(base, key),
                        "trait_key": trait,
                        "original_count": base.profiles[key][trait],
                        "remaining_count": vector[trait],
                        "original_total": sum(base.profiles[key].values()),
                        "remaining_total": sum(vector.values()),
                        "ambiguity_flag": any(taxon in ambiguous for taxon in base.reported[key]),
                    }
                    for key, vector in sorted(remaining.items())
                    for trait in base.traits.values()
                ),
            )
            backgrounds: list[Record] = []
            for scope, keys in scopes.items():
                for treatment, profiles in (
                    ("fixed_baseline", base.profiles),
                    ("recomputed_remaining", remaining),
                ):
                    counts = Counter(dict.fromkeys(base.columns, 0))
                    for key in keys:
                        counts.update(profiles[key])
                    total = sum(counts.values())
                    backgrounds.extend(
                        {
                            "scope": scope,
                            "treatment": treatment,
                            "trait_key": trait,
                            "count": value,
                            "full_total": total,
                            "profiles": len(keys),
                            "native_ora_status": "not_rerun_pending_Ilyas",
                        }
                        for trait, value in sorted(counts.items())
                    )
                effects = effect_rows(base, remaining, keys)
                for row in effects:
                    effect_count += 1
                    if scope == "pooled" and row["representation"] == "shared_fraction":
                        pooled_effects.append(dict(scenario=name, scope=scope, **row))
                    yield dict(scenario=name, scope=scope, **row)
            write_table(
                directory / "ora_background_counts.tsv",
                [
                    "scope",
                    "treatment",
                    "trait_key",
                    "count",
                    "full_total",
                    "profiles",
                    "native_ora_status",
                ],
                backgrounds,
            )

    example = effect_rows(base, base.profiles, [])[0]
    write_table(output / "effects.tsv.gz", ["scenario", "scope", *example], exports())
    write_table(
        output / "pooled_fraction_effects.tsv", ["scenario", "scope", *example], pooled_effects
    )
    result = {
        "scenarios": len(all_scenarios),
        "run_scenarios": len(active),
        "effect_rows": effect_count,
        "native_noel": "zero-total convention preserved; six-trait directional family",
        "native_jed": "not_rerun: transformed scale and correction family need owner integration",
        "native_ora": "not_rerun: both background count contracts exported for Ilyas",
    }
    summary = (
        "# Matched-corpus Table 1-trait analog source contributions and lineage influence\n\n"
        "> [!IMPORTANT]\n"
        "> **Matched-Corpus Table 1-Trait Analog (Issue #210)**\n"
        "> The taxon contribution and lineage influence audits reported here are a "
        "matched-corpus Table 1-trait analog sensitivity analysis, not influence on "
        "the published manuscript Table 1. They use the pinned human paired KG-Microbe "
        "profile corpus and are not a reproduction of the published Table 1 baseline "
        "(they do not reproduce published Table 1 baseline values such as facultatively "
        "aerobic $q = 2.8 \\times 10^{-10}$).\n\n"
        f"Reconciled {len(base.profiles):,} human paired profiles from "
        f"{len({k[0] for k in base.profiles}):,} studies across all "
        f"{len(base.columns)} non-MediaDive traits. "
        f"Enumerated {len(all_scenarios)} resolved genus/family lineages; "
        f"{len(active)} entered the "
        "rule-declared panel with at least three supporting studies. See manifest.json for "
        "hashes, "
        "exact scopes and protocol.\n\n"
        "The full contributing universe, concentration and all 15 cross-trait overlaps precede any "
        "display cutoff. Effects include all six raw annotation keys, both curated directions, "
        "shared fractions, all 15 log ratios and Noel's native zero-total convention. "
        "Each scenario preserves matched-pair changes and eligibility changes separately.\n\n"
        "These are provisional annotation influence diagnostics. Contrast orientation (#188) and "
        "source semantics (#189) remain unreviewed. Explicit reported ranks only resolve part of "
        "the pinned taxonomy; unresolved rank and source identity remain visible. Original source "
        "assertions cannot be counted from KG edge IDs alone. Higher-rank direct annotations are "
        "retained and flagged when their allocation to an excluded lineage is ambiguous.\n\n"
        "Native Jed transformed statistics and native ORA reruns remain for their owners; no "
        "substitute native q-values are supplied. Fixed and recomputed ORA backgrounds "
        "are distinct "
        "export contracts. Conditional intervals and separate BH6/BY15 families do not validate "
        "the full selection pipeline (#182). Loss of significance is not proof of explanation. "
        "These outputs feed #178 and the evidence limits in #179; they do not validate the "
        "Rosalind mechanisms, oxygen exposure, SCFA production, or causal ecological readouts.\n"
    )
    summary += "\nLargest contributing resolved genus, with share of **all** propagated counts:\n\n"
    summary += "| Trait key | Direction | Genus | Count share | Genus resolution coverage |\n"
    summary += "| --- | --- | --- | --- | --- |\n"
    for trait in base.traits.values():
        for direction in ("inc", "dec"):
            candidates = [
                r
                for r in universe
                if r["scope"] == "pooled"
                and r["trait_key"] == trait
                and r["direction"] == direction
                and r["rank"] == "genus"
                and r["status"] == "resolved"
            ]
            coverage = next(
                r
                for r in concentrations
                if r["scope"] == "pooled"
                and r["trait_key"] == trait
                and r["direction"] == direction
                and r["rank"] == "genus"
            )
            if candidates and coverage["total"]:
                leader = max(candidates, key=lambda r: (r["count"], r["taxon"]))
                label = trait.replace("|", "\\|")
                summary += (
                    f"| {label} | {direction} | {leader['name']} | "
                    f"{leader['count'] / coverage['total']:.1%} | "
                    f"{coverage['resolved_total'] / coverage['total']:.1%} |\n"
                )
    summary += (
        "\nPooled shared fraction effects, in percentage points. The displayed exclusion is "
        "the largest absolute matched change within the rule-declared panel; this is a "
        "descriptive "
        "ranking, not evidence that the lineage explains the association. All scenarios remain "
        "in `pooled_fraction_effects.tsv`, with uncertainty and eligibility components.\n\n"
        "| Trait key | Baseline effect | Exclusion | Matched change | "
        "Remaining effect | Lost pairs |\n"
        "| --- | --- | --- | --- | --- | --- |\n"
    )
    for trait in base.traits.values():
        candidates = [
            r
            for r in pooled_effects
            if r["trait_key"] == trait
            and r["scenario"] != "baseline"
            and r["matched_change_effect"] != ""
        ]
        if candidates:
            leader = max(candidates, key=lambda r: abs(r["matched_change_effect"]))
            label = trait.replace("|", "\\|")
            summary += (
                f"| {label} | {100 * leader['before_effect']:.4f} | {leader['scenario']} | "
                f"{100 * leader['matched_change_effect']:.4f} | "
                f"{100 * leader['after_effect']:.4f} | {leader['lost_pairs']} |\n"
            )
    (output / "influence_summary.md").write_text(summary, encoding="utf-8")
    return result


__all__ = [
    "TIE",
    "adjust",
    "comparison",
    "concentration",
    "critical_t",
    "describe",
    "directional_counts",
    "directional_test",
    "effect_rows",
    "overlap",
    "pair_manifest",
    "pair_values",
    "removal",
    "run_influence",
    "scenarios",
    "scope_keys",
    "source_occurrences",
    "study_means",
    "taxonomic_counts",
    "trait_leading_genera",
]
