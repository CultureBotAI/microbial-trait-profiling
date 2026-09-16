"""Pinned native statistical contracts on the reconciled lineage-removal corpus."""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.metadata
import json
import math
import platform
from collections import Counter
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path
from time import perf_counter

import numpy as np
from numpy.typing import NDArray
from scipy.stats import hypergeom, ttest_1samp

from bsdbng.contribution_inputs import Base, Key, Pair, digest, load_base, rows
from bsdbng.paths import DERIVED_DATA_DIR, REPO_ROOT
from bsdbng.taxon_contributions import Record, write_table
from bsdbng.taxon_influence import adjust, comparison, critical_t, directional_test, scope_keys

INPUT_NAMES = (
    "study_trait_profiles_wide.tsv",
    "experiments.tsv",
    "bugsigdb_taxa.tsv",
    "taxa_trait_matrix.tsv",
)


@lru_cache(maxsize=100000)
def null_count_interval(total: int, success: int, draws: int) -> tuple[int, int]:
    values = hypergeom.ppf([0.025, 0.975], total, success, draws)
    return int(values[0]), int(values[1])


def hypergeom_probabilities(
    triples: NDArray[np.int64],
    total: int,
    tail: str,
    workers: int = 8,
) -> NDArray[np.float64]:
    """Independent Boost-backed evaluations keep exact values and original order."""
    if tail not in {"enrichment", "depletion"}:
        raise ValueError(f"Unknown ORA tail: {tail}")

    def evaluate(block: NDArray[np.int64]) -> NDArray[np.float64]:
        values = (
            hypergeom.sf(block[:, 0] - 1, total, block[:, 1], block[:, 2])
            if tail == "enrichment"
            else hypergeom.cdf(block[:, 0], total, block[:, 1], block[:, 2])
        )
        return np.asarray(values, dtype=np.float64)

    if len(triples) < 4096 or workers == 1:
        return evaluate(triples)
    # np.unique sorts by count: contiguous chunks concentrate the expensive tail
    # in one worker. Striding balances those cases without changing cell order.
    result = np.empty(len(triples), dtype=np.float64)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for index, values in enumerate(
            pool.map(evaluate, [triples[i::workers] for i in range(workers)])
        ):
            result[index::workers] = values
    return result


def compressed_rows(path: Path) -> Iterator[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def reconstruct(base: Base, directory: Path) -> dict[Key, dict[str, int]]:
    removed: dict[str, Counter[str]] = {}
    for row in compressed_rows(directory / "removed_source_counts.tsv.gz"):
        target, trait = row["reported_taxon"], row["trait_key"]
        if trait not in base.columns or int(row["removed_count"]) <= 0:
            raise ValueError("Invalid sparse removal vector")
        vector = removed.setdefault(target, Counter())
        if trait in vector:
            raise ValueError("Duplicate sparse removal entry")
        vector[trait] = int(row["removed_count"])
    remaining = {}
    for key, original in base.profiles.items():
        subtraction: Counter[str] = Counter()
        for target in base.reported[key]:
            subtraction.update(removed.get(target, {}))
        vector_dict = {trait: value - subtraction[trait] for trait, value in original.items()}
        if min(vector_dict.values()) < 0:
            raise ValueError(f"Negative reconstructed count: {key}")
        remaining[key] = vector_dict
    expected = {(key, trait) for key in base.profiles for trait in base.traits.values()}
    for row in compressed_rows(directory / "profile_counts.tsv.gz"):
        key = (row["study_id"], row["experiment_id"], row["direction"])
        trait = row["trait_key"]
        identity = (key, trait)
        if identity not in expected:
            raise ValueError("Duplicate or unexpected scenario profile record")
        expected.remove(identity)
        actual = (
            base.profiles[key][trait],
            remaining[key][trait],
            sum(base.profiles[key].values()),
            sum(remaining[key].values()),
        )
        declared = tuple(
            int(row[f])
            for f in ("original_count", "remaining_count", "original_total", "remaining_total")
        )
        if actual != declared:
            raise ValueError(f"Scenario reconstruction mismatch: {key}, {trait}")
    if expected:
        raise ValueError("Missing scenario profile records")
    return remaining


def transformed(
    base: Base,
    profiles: dict[Key, dict[str, int]],
    keys: list[Key],
    policy: Record,
) -> tuple[dict[str, Record], dict[str, dict[Pair, float]], dict[str, float]]:
    """Mirror selection and t-test kernels, preserving the complete native family."""
    totals = {key: sum(profiles[key].values()) for key in keys}
    records: dict[str, Record] = {}
    contrasts: dict[str, dict[Pair, float]] = {}
    scales: dict[str, float] = {}
    for trait in base.columns:
        fractions = {
            key: profiles[key][trait] / totals[key] if totals[key] else 0.0 for key in keys
        }
        values = list(fractions.values())
        support = sum(profiles[key][trait] > 0 for key in keys)
        signal = sum(profiles[key][trait] for key in keys)
        var = float(np.var(values, ddof=1)) if len(values) > 1 else 0.0
        eligible = (
            support >= policy["minimum_profile_support"]
            and support <= max(2, int(len(keys) * policy["maximum_profile_support_fraction"]))
            and var > 0
        )
        logs = {
            key: math.log1p(policy["transform_multiplier"] * value)
            for key, value in fractions.items()
        }
        scale = float(np.std(list(logs.values()), ddof=0)) if logs else 0.0
        scales[trait] = scale
        # Centering cancels in a pair difference; use population SD (ddof=0).
        paired = {
            key[:2]: (logs[key] - logs[(*key[:2], "dec")]) / scale if scale else 0.0
            for key in keys
            if key[2] == "inc"
        }
        contrasts[trait] = paired
        items = list(paired.values())
        sd = float(np.std(items, ddof=1)) if len(items) > 1 else 0.0
        valid = eligible and len(items) >= policy["minimum_pairs"] and sd > 1e-8
        effect = float(np.mean(items)) if items else 0.0
        radius = critical_t(len(items) - 1) * sd / math.sqrt(len(items)) if len(items) > 1 else None
        positive, negative = sum(v > 0 for v in items), sum(v < 0 for v in items)
        share = directional_test(positive, negative) if positive + negative else None
        records[trait] = {
            "trait_key": trait,
            "selected": eligible,
            "profile_support": support,
            "total_signal": signal,
            "composition_variance": var,
            "selection_score": var * math.log1p(signal),
            "reference_sd": scale,
            "profiles": len(keys),
            "zero_total_profiles": sum(v == 0 for v in totals.values()),
            "pairs": len(items),
            "studies": len({k[0] for k in paired}),
            "effect": effect if items else "",
            "mean_ci_low": effect - radius if radius is not None else "",
            "mean_ci_high": effect + radius if radius is not None else "",
            "positive": positive,
            "negative": negative,
            "ties": len(items) - positive - negative,
            "positive_share": positive / (positive + negative) if share else "",
            "share_ci_low": share[1] if share else "",
            "share_ci_high": share[2] if share else "",
            "p": float(ttest_1samp(items, 0.0).pvalue) if valid else "",
            "q": "",
            "status": "tested" if valid else "not_selected" if not eligible else "untestable",
            "direction": "increased"
            if items and effect > 0
            else "decreased"
            if items and effect < 0
            else "null"
            if items
            else "untestable",
            "unit": "within_scope_refitted_zscore_contrast",
            "weighting": "equal_experiments",
        }
    tested = [r for r in records.values() if r["p"] != ""]
    for row, q in zip(tested, adjust([r["p"] for r in tested]), strict=True):
        row["q"] = q
    for row in records.values():
        if row["selected"] and row["p"] == "":
            row["q"] = 1.0
        row["family_tests"] = len(tested)
        row["family_selected_traits"] = sum(r["selected"] for r in records.values())
        for threshold in (0.001, 0.05):
            row[f"significant_{threshold}"] = row["status"] == "tested" and row["q"] <= threshold
    return records, contrasts, scales


def ora_cells(
    base: Base,
    profiles: dict[Key, dict[str, int]],
    keys: list[Key],
    background: Counter[str],
) -> tuple[list[Record], Record]:
    """Compute global native families before retaining the six-trait display rows."""
    total = sum(background.values())
    totals = {key: sum(profiles[key].values()) for key in keys}
    active = [key for key in keys if totals[key] > 0]
    six = set(base.traits.values())
    selected: dict[tuple[Key, str], Record] = {}
    for key in keys:
        for trait in base.traits.values():
            observed_count, success = profiles[key][trait], background[trait]
            expected = totals[key] * success / total if total else 0.0
            selected[(key, trait)] = {
                "study_id": key[0],
                "experiment_id": key[1],
                "direction": key[2],
                "trait_key": trait,
                "count": observed_count,
                "row_total": totals[key],
                "background_count": success,
                "background_total": total,
                "expected": expected,
                "fold": observed_count / expected if expected else "",
                "status": "tested"
                if expected
                else "zero_total"
                if not totals[key]
                else "zero_background_trait",
                "enrichment_p": "",
                "enrichment_q": "",
                "depletion_p": "",
                "depletion_q": "",
                "null_count_low": "",
                "null_count_high": "",
            }
    family: Record = {
        "profiles": len(keys),
        "testable_profiles": len(active),
        "background_total": total,
    }
    for tail in ("enrichment", "depletion"):
        identities = [
            (key, trait)
            for key in active
            for trait in base.columns
            if tail == "depletion" or profiles[key][trait] > 0
        ]
        observed = [profiles[key][trait] for key, trait in identities]
        successes = [background[trait] for _, trait in identities]
        draws = [totals[key] for key, _ in identities]
        if any(
            x > k or n > total or n - x > total - k
            for x, k, n in zip(observed, successes, draws, strict=True)
        ):
            raise ValueError("ORA counts do not fit the declared finite population")
        print(f"    ORA {tail}: {len(identities):,} tests", flush=True)
        started = perf_counter()
        # Repeated count/success/draw triples have identical native probabilities.
        # Expand back before BH so multiplicity still counts EVERY original test.
        triples, inverse = np.unique(
            np.asarray(list(zip(observed, successes, draws, strict=True)), dtype=np.int64).reshape(
                -1, 3
            ),
            axis=0,
            return_inverse=True,
        )
        pvalues: list[float] = (
            hypergeom_probabilities(triples, total, tail)[inverse].tolist() if identities else []
        )
        if not all(math.isfinite(p) and 0 <= p <= 1 for p in pvalues):
            raise ValueError("Nonfinite native ORA p-value")
        qvalues = adjust(list(pvalues))
        print(
            f"    {len(triples):,} unique probabilities and full BH "
            f"in {perf_counter() - started:.1f}s",
            flush=True,
        )
        family[f"{tail}_tests"] = len(identities)
        for (key, trait), p, q in zip(identities, pvalues, qvalues, strict=True):
            if trait in six:
                row = selected[(key, trait)]
                row[f"{tail}_p"], row[f"{tail}_q"] = p, q
    print("    ORA null count intervals", flush=True)
    for key in active:
        for trait in base.traits.values():
            lo, hi = null_count_interval(total, background[trait], totals[key])
            selected[(key, trait)].update(null_count_low=lo, null_count_high=hi)
    for row in selected.values():
        for tail in ("enrichment", "depletion"):
            row[f"{tail}_status"] = (
                "tested"
                if row[f"{tail}_p"] != ""
                else "zero_total"
                if not row["row_total"]
                else "zero_count_not_tested"
            )
            for threshold in (0.001, 0.05):
                row[f"{tail}_significant_{threshold}"] = (
                    row[f"{tail}_q"] != ""
                    and row[f"{tail}_q"] <= threshold
                    and row["fold"] != ""
                    and (row["fold"] >= 1.5 if tail == "enrichment" else row["fold"] <= 0.67)
                )
    return list(selected.values()), family


def fold_pairs(records: list[Record], trait: str) -> dict[Pair, float]:
    cells = {
        (r["study_id"], r["experiment_id"], r["direction"]): r["fold"]
        for r in records
        if r["trait_key"] == trait and r["fold"] != ""
    }
    return {
        key[:2]: value - cells[(*key[:2], "dec")]
        for key, value in cells.items()
        if key[2] == "inc" and (*key[:2], "dec") in cells
    }


def correct_shared_comparisons(records: list[Record], dimensions: tuple[str, ...]) -> None:
    """These sign-test families are separate from every native t/ORA family."""
    groups: dict[tuple[str, ...], list[Record]] = {}
    for row in records:
        groups.setdefault(tuple(row[field] for field in dimensions), []).append(row)
    for family in groups.values():
        for population in ("before", "after", "before_matched", "after_matched", "matched_change"):
            pvalues = [r[f"{population}_p"] if r[f"{population}_p"] != "" else None for r in family]
            for row, q in zip(family, adjust(pvalues), strict=True):
                row[f"{population}_q"] = q if q is not None else ""
                row["shared_correction"] = f"BH{len(family)}_directional_sign_test"


def composition_q_sensitivity_table(base: Base, composition: list[Record]) -> list[str]:
    """Summarize whether exclusions remove pooled native q<=0.05 membership."""
    lines = [
        "## Baseline-versus-exclusion native q-values",
        "",
        "This table directly checks whether any predeclared genus/family exclusion causes "
        "a Table 1 trait that is significant at pooled baseline to lose native q<=0.05 "
        "membership. Q-values come from the complete supported-trait family, not a six-trait "
        "recalculation. `NA` means the exclusion result was not tested or not selected.",
        "",
        "| Trait | Baseline q | Highest exclusion q | Highest-q exclusion | "
        "Exclusions losing q<=0.05 | Losing exclusions |",
        "| --- | ---: | ---: | --- | ---: | --- |",
    ]
    for slug, trait in base.traits.items():
        baseline = next(
            row
            for row in composition
            if row["scenario"] == "baseline"
            and row["scope"] == "pooled"
            and row["reference"] == "refitted"
            and row["trait_key"] == trait
        )
        candidates = [
            row
            for row in composition
            if row["scenario"] != "baseline"
            and row["scope"] == "pooled"
            and row["reference"] == "refitted"
            and row["trait_key"] == trait
        ]
        baseline_q = baseline["native_before_q"]
        baseline_significant = (
            baseline["native_before_status"] == "tested"
            and baseline_q != ""
            and float(baseline_q) <= 0.05
        )
        if candidates:
            worst = max(
                candidates,
                key=lambda row: (
                    float(row["native_after_q"]) if row["native_after_q"] != "" else math.inf
                ),
            )
            worst_q = (
                f"{float(worst['native_after_q']):.6g}"
                if worst["native_after_q"] != ""
                else f"NA ({worst['native_after_status']})"
            )
            worst_scenario = str(worst["scenario"])
        else:
            worst_q = "none"
            worst_scenario = "none"
        if baseline_significant:
            losses = [
                row
                for row in candidates
                if row["native_after_status"] != "tested"
                or row["native_after_q"] == ""
                or float(row["native_after_q"]) > 0.05
            ]
            loss_count = str(len(losses))
            loss_names = "<br>".join(str(row["scenario"]) for row in losses) or "none"
        else:
            loss_count = "n/a"
            loss_names = "baseline not significant"
        baseline_text = f"{float(baseline_q):.6g}" if baseline_q != "" else "NA"
        lines.append(
            f"| {slug} | {baseline_text} | {worst_q} | {worst_scenario} | "
            f"{loss_count} | {loss_names} |"
        )
    return lines


def native_summary(base: Base, composition: list[Record], ora: list[Record]) -> str:
    lines = [
        "# Native lineage influence completion",
        "",
        "These reruns extend the reconciled 190-v1 audit using protocol 190-native-v1. "
        "The 22 exclusion scenarios, four scopes, raw trait keys and human paired corpus "
        "are unchanged. This is matched-corpus sensitivity analysis, not a reproduction "
        "of published results with other filters. See native_manifest.json for hashes, "
        "the implementation revision, model policies and result files.",
        "",
        "## Composition",
        "",
        "Native effects are experiment-weighted increased-minus-decreased z-score contrasts. "
        "Native p-values test the mean; native BH families contain all supported variable "
        "traits, not just Table 1. Zero-total profiles remain as zero fractions. "
        "Study-weighted comparisons, sign tests and conditional intervals are separately "
        "labeled in native_composition_comparison.tsv. Refitted scales and fixed baseline "
        "SD comparisons must not be conflated.",
        "",
        "| Trait | Pooled native effect | Native p | Native q | "
        "Largest matched change, fixed baseline SD | Exclusion |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for slug, trait in base.traits.items():
        baseline = next(
            r
            for r in composition
            if r["scenario"] == "baseline" and r["scope"] == "pooled" and r["trait_key"] == trait
        )
        candidates = [
            r
            for r in composition
            if r["scenario"] != "baseline"
            and r["scope"] == "pooled"
            and r["reference"] == "fixed_baseline_sd"
            and r["trait_key"] == trait
            and r["matched_change_effect"] != ""
        ]
        leader = (
            max(candidates, key=lambda r: abs(r["matched_change_effect"])) if candidates else None
        )
        native_values = [
            f"{baseline[f'native_before_{f}']:.6g}"
            if baseline[f"native_before_{f}"] != ""
            else "untestable"
            for f in ("effect", "p", "q")
        ]
        change = f"{leader['matched_change_effect']:.6g}" if leader else "untestable"
        scenario = leader["scenario"] if leader else "none"
        lines.append(f"| {slug} | {' | '.join(native_values)} | {change} | {scenario} |")
    lines.extend(
        [
            "",
            "The displayed maximum is descriptive selection within the predeclared panel, "
            "not an explanation or causal attribution. The TSV retains all opposite, null and "
            "untestable results, eligibility components, intervals and family membership.",
            "",
            *composition_q_sensitivity_table(base, composition),
            "",
            "## ORA",
            "",
            "Native enrichment and depletion were rerun against fixed-baseline and recomputed "
            "scope backgrounds separately. The global BH families use every eligible retained "
            "trait/profile cell. Enrichment tests positive cells; depletion includes zeros in "
            "nonzero-total profiles. All six traits remain in the export, "
            "not only significant hits. "
            "Fold limits are >=1.5 for enrichment and <=0.67 for depletion.",
            "",
            "The profile table includes exact 95% central hypergeometric null count intervals: "
            "these are NOT confidence intervals for fold effects. The fold-comparison TSV "
            "instead gives study-weighted paired fold differences, conditional mean intervals, "
            "directional support and matched versus eligibility changes. It also reports "
            "native hit counts at q<=0.001 and q<=0.05 separately by curated direction.",
            "",
            "| Trait | Background | Largest absolute matched fold change | "
            "Exclusion | Lost pairs |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for slug, trait in base.traits.items():
        for treatment in ("fixed_baseline", "recomputed_remaining"):
            candidates = [
                r
                for r in ora
                if r["scenario"] != "baseline"
                and r["scope"] == "pooled"
                and r["trait_key"] == trait
                and r["background"] == treatment
                and r["matched_change_effect"] != ""
            ]
            if candidates:
                leader = max(candidates, key=lambda r: abs(r["matched_change_effect"]))
                lines.append(
                    f"| {slug} | {treatment} | {leader['matched_change_effect']:.6g} | "
                    f"{leader['scenario']} | {leader['lost_pairs']} |"
                )
    lines.extend(
        [
            "",
            "## Integration Limits",
            "",
            "For #178, use qualified trait keys, exact scope, scenario, background/reference, "
            "effect units and correction-family fields. Keep native and shared statistics "
            "separate; do not combine ORA backgrounds in a common ranking. ORA profile-count "
            "tests are not independent corroboration of study-level effects. The native "
            "composition t intervals assume independent experiments; study-weighted summaries "
            "are supplied separately, without silently changing the native model (#133).",
            "",
            "For #179, these are annotation influence diagnostics. Higher-rank direct assertions "
            "that cannot be allocated remain retained and flagged in the source ledger. Missing "
            "source identities/ranks are not repaired by statistical reruns. Disease/control "
            "orientation (#188), feature semantics (#189), source-count validity (#184) and "
            "full-pipeline calibration (#182) remain independent review requirements. No "
            "causal, oxygen-exposure, SCFA or taxonomy-independent mechanism follows from "
            "these results. Native kernel/ORA owner review and integration sign-off are still "
            "required before manuscript adoption.",
            "",
        ]
    )
    return "\n".join(lines)


def export_native(base: Base, audit: Record, protocol: Record, output: Path) -> Record:
    scopes = {s["id"]: scope_keys(base, s) for s in audit["protocol"]["scopes"]}
    scenario_names = [
        "baseline",
        *(r["scenario"] for r in rows(output / "scenarios.tsv") if r["status"] == "run"),
    ]
    backgrounds = {
        scope: Counter(
            {trait: sum(base.profiles[k][trait] for k in keys) for trait in base.columns}
        )
        for scope, keys in scopes.items()
    }
    baseline_composition = {
        scope: transformed(base, base.profiles, keys, protocol["composition"])
        for scope, keys in scopes.items()
    }
    baseline_ora: dict[str, list[Record]] = {}
    native_stats: list[Record] = []
    native_comparisons: list[Record] = []
    ora_comparisons: list[Record] = []
    families: list[Record] = []
    output_hashes: dict[str, str] = {}
    input_hashes: dict[str, str] = {"scenarios.tsv": digest(output / "scenarios.tsv")}
    for name in scenario_names:
        print(f"Native influence: {name}", flush=True)
        directory = output / name.replace(":", "_")
        for filename in ("removed_source_counts.tsv.gz", "profile_counts.tsv.gz"):
            path = directory / filename
            input_hashes[path.relative_to(output).as_posix()] = digest(path)
        remaining = reconstruct(base, directory)
        cells: list[Record] = []
        for scope, keys in scopes.items():
            print(f"  Scope {scope}: {len(keys)} profiles", flush=True)
            before, bpairs, bscales = baseline_composition[scope]
            after, apairs, ascales = transformed(base, remaining, keys, protocol["composition"])
            native_stats.extend(dict(scenario=name, scope=scope, **r) for r in after.values())
            for trait in base.traits.values():
                native = {
                    f"{stage}_{field}": value
                    for stage, records in (("before", before), ("after", after))
                    for field, value in records[trait].items()
                    if field != "trait_key"
                }
                anchored = (
                    {
                        pair: value * ascales[trait] / bscales[trait]
                        for pair, value in apairs[trait].items()
                    }
                    if bscales[trait]
                    else {}
                )
                for reference, values in (
                    ("refitted", apairs[trait]),
                    ("fixed_baseline_sd", anchored),
                ):
                    native_comparisons.append(
                        dict(
                            scenario=name,
                            scope=scope,
                            trait_key=trait,
                            reference=reference,
                            **{f"native_{k}": v for k, v in native.items()},
                            **comparison(
                                bpairs[trait] if reference == "refitted" or bscales[trait] else {},
                                values,
                            ),
                        )
                    )
            for treatment in protocol["ora"]["backgrounds"]:
                background = (
                    backgrounds[scope]
                    if treatment == "fixed_baseline"
                    else Counter(
                        {trait: sum(remaining[k][trait] for k in keys) for trait in base.columns}
                    )
                )
                if name == "baseline" and scope in baseline_ora:
                    records = baseline_ora[scope]
                    family = next(
                        r for r in families if r["scenario"] == name and r["scope"] == scope
                    )
                    family = {
                        k: v
                        for k, v in family.items()
                        if k not in ("scenario", "scope", "background")
                    }
                else:
                    records, family = ora_cells(base, remaining, keys, background)
                if name == "baseline":
                    baseline_ora[scope] = records
                families.append(dict(scenario=name, scope=scope, background=treatment, **family))
                cells.extend(dict(scope=scope, background=treatment, **r) for r in records)
                for trait in base.traits.values():
                    result = dict(
                        scenario=name,
                        scope=scope,
                        background=treatment,
                        trait_key=trait,
                        unit="fold_difference_inc_minus_dec",
                        **comparison(
                            fold_pairs(baseline_ora[scope], trait), fold_pairs(records, trait)
                        ),
                    )
                    for stage, items in (("before", baseline_ora[scope]), ("after", records)):
                        for direction in ("inc", "dec"):
                            for tail in ("enrichment", "depletion"):
                                for threshold in (0.001, 0.05):
                                    result[f"{stage}_{direction}_{tail}_hits_{threshold}"] = sum(
                                        r["trait_key"] == trait
                                        and r["direction"] == direction
                                        and r[f"{tail}_significant_{threshold}"]
                                        for r in items
                                    )
                    ora_comparisons.append(result)
        path = directory / "native_ora_profiles.tsv.gz"
        write_table(path, list(cells[0]), cells)
        output_hashes[path.relative_to(output).as_posix()] = digest(path)
    correct_shared_comparisons(native_comparisons, ("scenario", "scope", "reference"))
    correct_shared_comparisons(ora_comparisons, ("scenario", "scope", "background"))
    for filename, records in (
        ("native_composition.tsv.gz", native_stats),
        ("native_composition_comparison.tsv", native_comparisons),
        ("native_ora_fold_comparison.tsv", ora_comparisons),
        ("native_ora_families.tsv", families),
    ):
        write_table(output / filename, list(records[0]), records)
        output_hashes[filename] = digest(output / filename)
    (output / "native_summary.md").write_text(
        native_summary(base, native_comparisons, ora_comparisons), encoding="utf-8"
    )
    output_hashes["native_summary.md"] = digest(output / "native_summary.md")
    return {
        "scenarios": len(scenario_names),
        "scopes": len(scopes),
        "composition_rows": len(native_stats),
        "composition_comparisons": len(native_comparisons),
        "ora_comparisons": len(ora_comparisons),
        "ora_families": len(families),
        "input_sha256": input_hashes,
        "output_sha256": output_hashes,
    }


def run(inputs: Path, output: Path, protocol_path: Path, code_commit: str = "") -> None:
    audit = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    if not audit["completed"] or audit["mismatches"]:
        raise ValueError("A completed, exactly reconciled contribution audit is required")
    for name in INPUT_NAMES:
        if digest(inputs / name) != audit["input_sha256"][name]:
            raise ValueError(f"Audit input changed: {name}")
    base = load_base(*(inputs / name for name in INPUT_NAMES))
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    if protocol["version"] != "190-native-v1":
        raise ValueError("Unsupported native protocol")
    if len(code_commit) != 40 or any(c not in "0123456789abcdef" for c in code_commit):
        raise ValueError("The Justfile recipe must supply the full code commit")
    manifest: Record = {
        "completed": False,
        "protocol": protocol,
        "protocol_sha256": digest(protocol_path),
        "contribution_manifest_sha256": digest(output / "manifest.json"),
        "input_sha256": {name: digest(inputs / name) for name in INPUT_NAMES},
        "code_commit": code_commit,
        "runtime": {
            "python": platform.python_version(),
            "scipy": importlib.metadata.version("scipy"),
        },
        "code_sha256": {
            name: digest(Path(__file__).parent / name)
            for name in (
                "native_taxon_influence.py",
                "taxon_influence.py",
                "contribution_inputs.py",
                "taxon_contributions.py",
            )
        },
    }
    path = output / "native_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    manifest["results"] = export_native(base, audit, protocol, output)
    manifest["completed"] = True
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--inputs", type=Path, default=DERIVED_DATA_DIR)
    parser.add_argument("--output", type=Path, default=DERIVED_DATA_DIR / "table1_taxon_influence")
    parser.add_argument(
        "--protocol", type=Path, default=REPO_ROOT / "methods/table1_native_influence_protocol.json"
    )
    args = parser.parse_args()
    run(args.inputs, args.output, args.protocol, args.code_commit)


if __name__ == "__main__":
    main()
