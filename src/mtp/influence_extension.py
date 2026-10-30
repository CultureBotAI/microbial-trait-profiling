"""Additive lineage-extension and six-fraction study-LOSO diagnostics (#771).

The source population, lineage ranking and complete six-member corrections are
explicit. None of these exploratory sensitivities changes the native arm plan.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path

from mtp.contribution_inputs import Base, Graph, Key, encoded, load_base, load_graph
from mtp.paths import REPO_ROOT, file_sha256
from mtp.taxon_contributions import Ledger, Record, reconcile, write_table
from mtp.taxon_influence import (
    TIE,
    adjust,
    comparison,
    describe,
    pair_values,
    removal,
    scope_keys,
    taxonomic_counts,
)
from mtp.workflow.arms.runtime import verify_arm
from mtp.workflow.consumers import (
    open_core_bundle,
    preflight_consumer_outputs,
    write_consumer_record,
)
from mtp.workflow.output_roots import add_temporary_output_option, require_durable_output

PROTOCOL = REPO_ROOT / "methods/table1_influence_extension_protocol.json"
PROTOCOL_SHA256 = "366857b28d00228c52323fc6c7ed68848c037d09fa06ed15c8f453f553aebfd6"
PANEL_FIELDS = [
    "scope",
    "rank",
    "lineage",
    "lineages",
    "name",
    "scenario",
    "selection_kind",
    "selection_reasons",
]
FRACTION_FIELDS = ["trait_key", "family_size", *describe({}), "q_bh_six", "q_by_six"]
BASELINE_FIELDS = ["baseline_effect", "baseline_q_bh_six", "baseline_q_by_six"]
REMOVAL_FIELDS = [
    "removed_annotation_mass",
    "removed_trait_annotation_mass",
    "removal_status",
    "lineage_accounting",
]
REDUCTION_FIELDS = [
    f"{prefix}effect_reduction_{field}"
    for prefix in ("", "matched_")
    for field in ("fraction", "percent", "status")
]
MASS_FIELDS = [
    "direction",
    "trait_key",
    "population",
    "unit",
    "profile_count",
    "before_annotation_mass",
    "after_annotation_mass",
    "removed_annotation_mass",
    "before_trait_annotation_mass",
    "after_trait_annotation_mass",
    "removed_trait_annotation_mass",
    "removal_status",
]
IMPLEMENTATIONS = (
    "influence_extension.py",
    "contribution_inputs.py",
    "taxon_contributions.py",
    "taxon_influence.py",
    "taxa_trait_matrix.py",
    "trait_guilds.py",
    "workflow/operations.py",
    "workflow/consumers.py",
    "workflow/arms/runtime.py",
)


def implementation_hashes() -> dict[str, str]:
    """Hash the declared analysis producers before and after execution."""
    return {name: file_sha256(Path(__file__).parent / name) for name in IMPLEMENTATIONS}


def load_protocol() -> Record:
    """Require the reviewed rule bytes before any analysis or output publication."""
    raw = PROTOCOL.read_bytes()
    if hashlib.sha256(raw).hexdigest() != PROTOCOL_SHA256:
        raise ValueError("Lineage-extension protocol hash differs from the declared rule")
    protocol: Record = json.loads(raw)
    if protocol["version"] != "771-extension-v1":
        raise ValueError("Unknown lineage-extension protocol version")
    return protocol


def _focal_traits(base: Base) -> list[str]:
    traits = list(base.traits.values())
    if len(traits) != 6 or len(set(traits)) != 6:
        raise ValueError("The diagnostic family must contain exactly six distinct focal traits")
    return traits


def eligible_keys(profiles: Mapping[Key, Mapping[str, int]], keys: list[Key]) -> list[Key]:
    """Keep entire pairs only when both retained annotation totals are positive."""
    selected = set(keys)
    return sorted(
        key
        for key in selected
        if (*key[:2], "inc") in selected
        and (*key[:2], "dec") in selected
        and sum(profiles[(*key[:2], "inc")].values()) > 0
        and sum(profiles[(*key[:2], "dec")].values()) > 0
    )


def fraction_family(
    base: Base, profiles: dict[Key, dict[str, int]], keys: list[Key]
) -> list[Record]:
    """Recompute all six p-values before either complete-family diagnostic adjustment."""
    traits = _focal_traits(base)
    records = [
        {"trait_key": trait, "family_size": 6, **describe(pair_values(profiles, keys, trait))}
        for trait in traits
    ]
    values = [float(row["p"]) if row["p"] != "" else None for row in records]
    for row, bh, by in zip(records, adjust(values), adjust(values, by=True), strict=True):
        row["q_bh_six"] = "" if bh is None else bh
        row["q_by_six"] = "" if by is None else by
    return records


def study_loso(base: Base, keys: list[Key]) -> list[Record]:
    """Omit every experiment of one eligible study and recorrect the whole family."""
    _focal_traits(base)
    eligible = eligible_keys(base.profiles, keys)
    studies = sorted({key[0] for key in eligible})
    baseline = {row["trait_key"]: row for row in fraction_family(base, base.profiles, eligible)}
    return [
        {
            "omitted_study": study,
            **row,
            **{
                f"baseline_{field}": baseline[row["trait_key"]][field]
                for field in ("effect", "q_bh_six", "q_by_six")
            },
        }
        for study in studies
        for row in fraction_family(
            base, base.profiles, [key for key in eligible if key[0] != study]
        )
    ]


def lineage_panel(
    base: Base, graph: Graph, ledger: Ledger, protocol: Record
) -> tuple[list[Record], list[Record]]:
    """Rank each trait/direction/rank/scope independently; preserve every candidate."""
    traits = _focal_traits(base)
    top = protocol["top_per_trait_direction_rank_scope"]
    minimum = protocol["minimum_lineage_studies"]
    if isinstance(top, bool) or not isinstance(top, int) or top < 1:
        raise ValueError("top_per_trait_direction_rank_scope must be a positive integer")
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 1:
        raise ValueError("minimum_lineage_studies must be a positive integer")
    if protocol["ranks"] != ["family", "genus"] or protocol["directions"] != ["inc", "dec"]:
        raise ValueError("Extension protocol requires both family/genus and inc/dec")
    candidates: list[Record] = []
    selected: list[Record] = []
    scope_ids = [scope["id"] for scope in protocol["scopes"]]
    if len(set(scope_ids)) != len(scope_ids):
        raise ValueError("Duplicate extension scope")
    for scope in protocol["scopes"]:
        # Match the declared contributing_universe; effect eligibility is separate.
        keys = scope_keys(base, scope)
        _, records = taxonomic_counts(base, graph, ledger, keys)
        reasons: defaultdict[tuple[str, str], list[str]] = defaultdict(list)
        for trait in traits:
            for direction in protocol["directions"]:
                for rank in protocol["ranks"]:
                    group = sorted(
                        [
                            row
                            for row in records
                            if row["trait_key"] == trait
                            and row["direction"] == direction
                            and row["rank"] == rank
                        ],
                        key=lambda row: (-row["count"], row["taxon"]),
                    )
                    eligible_rank = 0
                    for row in group:
                        supported = row["status"] == "resolved" and row["studies"] >= minimum
                        if supported:
                            eligible_rank += 1
                        chosen = supported and eligible_rank <= top
                        candidates.append(
                            {
                                "scope": scope["id"],
                                **row,
                                "eligible_rank": eligible_rank if supported else "",
                                "selected": chosen,
                                "selection_status": "selected"
                                if chosen
                                else "unresolved_rank"
                                if row["status"] != "resolved"
                                else "insufficient_studies"
                                if not supported
                                else "outside_top_k",
                            }
                        )
                        if chosen:
                            reasons[(rank, row["taxon"])].append(f"{trait}:{direction}:{rank}")
        for (rank, lineage), why in sorted(reasons.items()):
            selected.append(
                {
                    "scope": scope["id"],
                    "rank": rank,
                    "lineage": lineage,
                    "lineages": encoded([lineage]),
                    "name": graph.names.get(lineage, base.names.get(lineage, "")),
                    "scenario": f"extension_exclude_{rank}_{lineage}",
                    "selection_kind": "ranked",
                    "selection_reasons": encoded(sorted(why)),
                }
            )
        for joint in protocol["joint_exclusions"]:
            lineages = sorted(set(joint["lineages"]))
            names = dict(zip(joint["lineages"], joint["names"], strict=True))
            known_ranks = {base.ranks[lineage] for lineage in lineages if base.ranks.get(lineage)}
            if known_ranks - {"family"}:
                raise ValueError("Declared joint family exclusion contains a non-family lineage")
            selected.append(
                {
                    "scope": scope["id"],
                    "rank": "family",
                    "lineage": "+".join(lineages),
                    "lineages": encoded(lineages),
                    "name": " + ".join(names[lineage] for lineage in lineages),
                    "scenario": f"extension_exclude_joint_{joint['id']}",
                    "selection_kind": "declared_joint",
                    "selection_reasons": encoded([joint["reason"]]),
                }
            )
    return candidates, selected


def joint_removal(
    base: Base, graph: Graph, ledger: Ledger, lineages: list[str]
) -> tuple[dict[Key, dict[str, int]], dict[str, set[str]]]:
    """Remove the union of source units once, including under multiple ancestry paths."""
    excluded = set(lineages)
    if not excluded:
        raise ValueError("Joint removal requires at least one lineage")
    above = set().union(*(graph.ancestors(lineage) - {lineage} for lineage in excluded))
    removed: dict[str, Counter[str]] = {}
    ambiguous: dict[str, set[str]] = {}
    for target, sources in ledger.sources.items():
        vector: Counter[str] = Counter()
        for source in sources:
            if excluded.intersection(graph.ancestors(source)):
                vector.update(graph.traits[source])
            elif source in above:
                ambiguous.setdefault(target, set()).add(source)
        removed[target] = vector
    remaining = {}
    for key, baseline in base.profiles.items():
        subtract: Counter[str] = Counter()
        for target in base.reported[key]:
            subtract.update(removed[target])
        vector_dict = {trait: value - subtract[trait] for trait, value in baseline.items()}
        if any(value < 0 for value in vector_dict.values()):
            raise ValueError(f"Negative remainder for {key}; joint reconstruction is inconsistent")
        remaining[key] = vector_dict
    return remaining, ambiguous


def _present_lineages(graph: Graph, wanted: set[str]) -> set[str]:
    """Check actual node/edge inventories; ancestor self-inclusion is not presence."""
    present = {
        lineage
        for lineage in wanted
        if lineage in graph.names or lineage in graph.parents or lineage in graph.traits
    }
    if present != wanted:
        for parents in graph.parents.values():
            present.update(wanted.intersection(parents))
            if present == wanted:
                break
    return present


def _individual_lineage_mass(
    base: Base,
    graph: Graph,
    ledger: Ledger,
    keys: list[Key],
    lineages: list[str],
    remaining: dict[Key, dict[str, int]],
) -> dict[str, Counter[str]]:
    """Audit each exclusion independently on baseline support; overlaps may be counted twice."""
    masses = {}
    for lineage in lineages:
        individual = remaining if len(lineages) == 1 else removal(base, graph, ledger, lineage)[0]
        counts: Counter[str] = Counter()
        for key in keys:
            counts.update(
                {
                    trait: base.profiles[key][trait] - individual[key][trait]
                    for trait in base.columns
                }
            )
        masses[lineage] = counts
    return masses


def lineage_effects(
    base: Base,
    graph: Graph,
    ledger: Ledger,
    selected: list[Record],
    protocol: Record,
    *,
    removed_mass_rows: list[Record] | None = None,
) -> tuple[list[Record], list[Record]]:
    """Recompute denominator and eligibility; compare identical support separately."""
    scopes = {scope["id"]: scope for scope in protocol["scopes"]}
    present = _present_lineages(
        graph, {lineage for scenario in selected for lineage in json.loads(scenario["lineages"])}
    )
    effects: list[Record] = []
    ambiguities: list[Record] = []
    for scenario in selected:
        lineages = sorted(set(json.loads(scenario["lineages"])))
        keys = eligible_keys(base.profiles, scope_keys(base, scopes[scenario["scope"]]))
        if scenario["selection_kind"] == "declared_joint":
            remaining, ambiguous = joint_removal(
                base, graph, ledger, json.loads(scenario["lineages"])
            )
        else:
            remaining, _, ambiguous = removal(base, graph, ledger, scenario["lineage"])
        individual_mass = _individual_lineage_mass(base, graph, ledger, keys, lineages, remaining)
        traits = _focal_traits(base)
        removed_total = 0
        removed_focal: Counter[str] = Counter()
        for direction in ("inc", "dec"):
            directional_keys = [key for key in keys if key[2] == direction]
            before_counts: Counter[str] = Counter()
            after_counts: Counter[str] = Counter()
            for key in directional_keys:
                before_counts.update(base.profiles[key])
                after_counts.update(remaining[key])
            before_mass = sum(before_counts.values())
            after_mass = sum(after_counts.values())
            removed_mass = before_mass - after_mass
            if removed_mass < 0:
                raise ValueError("Lineage exclusion increased annotation mass")
            removed_total += removed_mass
            for trait in base.columns:
                before_trait = before_counts[trait]
                after_trait = after_counts[trait]
                if after_trait > before_trait:
                    raise ValueError("Lineage exclusion increased a trait annotation count")
                removed_focal[trait] += before_trait - after_trait
                if removed_mass_rows is not None:
                    removed_mass_rows.append(
                        {
                            **scenario,
                            "direction": direction,
                            "trait_key": trait,
                            "population": "eligible_baseline_pairs",
                            "unit": "retained_trait_annotation_count",
                            "profile_count": len(directional_keys),
                            "before_annotation_mass": before_mass,
                            "after_annotation_mass": after_mass,
                            "removed_annotation_mass": removed_mass,
                            "before_trait_annotation_mass": before_trait,
                            "after_trait_annotation_mass": after_trait,
                            "removed_trait_annotation_mass": before_trait - after_trait,
                            "removal_status": "removed_annotation_mass"
                            if removed_mass
                            else "no_removed_annotation_mass",
                        }
                    )
        family = []
        for trait in traits:
            before = pair_values(base.profiles, keys, trait)
            after = pair_values(remaining, keys, trait)
            contrasts = comparison(before, after)
            family.append(
                {
                    **scenario,
                    "trait_key": trait,
                    "family_size": 6,
                    **contrasts,
                    "removed_annotation_mass": removed_total,
                    "removed_trait_annotation_mass": removed_focal[trait],
                    "removal_status": "removed_annotation_mass"
                    if removed_total
                    else "no_removed_annotation_mass",
                    "lineage_accounting": encoded(
                        {
                            lineage: {
                                "graph_status": "present_in_graph"
                                if lineage in present
                                else "absent_from_graph",
                                "removed_annotation_mass": sum(counts.values()),
                                "removed_trait_annotation_mass": counts[trait],
                            }
                            for lineage, counts in individual_mass.items()
                        }
                    ),
                    **effect_reduction(contrasts["before_effect"], contrasts["after_effect"]),
                    **effect_reduction(
                        contrasts["before_matched_effect"],
                        contrasts["after_matched_effect"],
                        prefix="matched_",
                    ),
                }
            )
        for population in ("before", "after", "before_matched", "after_matched", "matched_change"):
            values = [
                float(row[f"{population}_p"]) if row[f"{population}_p"] != "" else None
                for row in family
            ]
            for row, bh, by in zip(family, adjust(values), adjust(values, by=True), strict=True):
                row[f"{population}_q_bh_six"] = "" if bh is None else bh
                row[f"{population}_q_by_six"] = "" if by is None else by
        effects.extend(family)
        active_targets = {taxon for key in keys for taxon in base.reported[key]}
        ambiguities.extend(
            {
                "scope": scenario["scope"],
                "scenario": scenario["scenario"],
                "reported_taxon": target,
                "retained_source_taxon": source,
                "reason": "direct_assertion_above_excluded_lineage_retained",
            }
            for target, sources in sorted(ambiguous.items())
            if target in active_targets
            for source in sorted(sources)
        )
    return effects, ambiguities


def effect_reduction(before: float | str, after: float | str, *, prefix: str = "") -> Record:
    """Describe 1-after/before without clipping or substituting zero for undefined data."""
    fraction: float | str = ""
    if before == "":
        status = "no_before_effect"
    elif abs(float(before)) <= TIE:
        status = "zero_before_effect"
    elif after == "":
        status = "no_after_effect"
    else:
        status = "defined"
        fraction = 1 - float(after) / float(before)
    return {
        f"{prefix}effect_reduction_fraction": fraction,
        f"{prefix}effect_reduction_percent": 100 * fraction if isinstance(fraction, float) else "",
        f"{prefix}effect_reduction_status": status,
    }


def loso_summary(baseline: list[Record], omissions: list[Record]) -> list[Record]:
    """Descriptive ranges across dependent omissions; no across-omission test."""
    result = []
    for original in baseline:
        rows = [
            r
            for r in omissions
            if r["scope"] == original["scope"] and r["trait_key"] == original["trait_key"]
        ]
        effects = [float(r["effect"]) for r in rows if r["effect"] != ""]
        qvalues = [float(r["q_by_six"]) for r in rows if r["q_by_six"] != ""]
        bh_values = [float(r["q_bh_six"]) for r in rows if r["q_bh_six"] != ""]
        original_effect = original["effect"]
        result.append(
            {
                "scope": original["scope"],
                "trait_key": original["trait_key"],
                "baseline_effect": original_effect,
                "baseline_q_bh_six": original["q_bh_six"],
                "baseline_q_by_six": original["q_by_six"],
                "omissions": len(rows),
                "testable_omissions": len(qvalues),
                "minimum_effect": min(effects) if effects else "",
                "maximum_effect": max(effects) if effects else "",
                "minimum_q_bh_six": min(bh_values) if bh_values else "",
                "maximum_q_bh_six": max(bh_values) if bh_values else "",
                "minimum_q_by_six": min(qvalues) if qvalues else "",
                "maximum_q_by_six": max(qvalues) if qvalues else "",
                "mean_shift_sign_reversal_count": sum(
                    value * float(original_effect) < 0 for value in effects
                )
                if original_effect != ""
                else "",
                "interpretation": "descriptive_dependent_omission_range",
            }
        )
    return result


def native_summary_source(core_manifest: Path, noel_manifest: Path) -> Path:
    """Locate the declared native summary through a verified same-core Noel run."""
    run = verify_arm(noel_manifest)
    if run.arm != "noel" or run.core_manifest_sha256 != file_sha256(core_manifest):
        raise ValueError("Native summary must come from a verified Noel run of this exact core")
    paths = [
        Path(artifact.path)
        for step in run.steps or []
        if step.id == "native_taxon_influence"
        for artifact in step.outputs or []
        if artifact.storage_kind == "file" and Path(artifact.path).name == "native_summary.md"
    ]
    if len(paths) != 1:
        raise ValueError("Verified Noel run must declare exactly one native_summary.md")
    return paths[0]


def _display_row(row: Record) -> Record:
    """Keep TSV precision intact while making report cells readable."""
    return {
        key: format(value, ".6g") if isinstance(value, float) else value
        for key, value in row.items()
    }


def combined_summary(
    native_text: str | None,
    panel: list[Record],
    summaries: list[Record],
    effects: list[Record] | None = None,
) -> str:
    """Show the native report and clearly separate the additive dependent diagnostics."""
    baseline = (
        native_text
        if native_text is not None
        else ("Native report was not supplied. Pass --noel-run to include its verified report.\n")
    )
    text = "# Native report and additive sensitivity\n\n"
    text += "The following text is copied from the verified native report when supplied. "
    text += "Its estimators remain distinct from the study-mean fraction diagnostics below.\n\n"
    text += baseline + "\n\n## Additional lineage exclusions\n\n"
    text += "| Scope | Rank | Lineage | Name |\n|---|---|---|---|\n"
    for row in panel:
        text += (
            "| "
            + " | ".join(
                str(row[key]).replace("|", "\\|") for key in ("scope", "rank", "lineage", "name")
            )
            + " |\n"
        )
    text += "\nSee lineage_fraction_effects.tsv for every effect, matched-support comparison "
    text += "and six-family diagnostic q-value. No exclusion is hidden by its result.\n\n"
    text += "The full-data shared-fraction baseline below uses the 190-native-v2 "
    text += "shared_fraction BH6 family in effects.tsv.gz (pooled_fraction_effects.tsv "
    text += "for pooled rows), recomputed on this same core and scope. BY-six is an "
    text += "additional correction sensitivity; native all-trait composition q-values "
    text += "refer to a different family. Every baseline column below means this full-data "
    text += "shared-fraction calculation. Report numbers use six significant digits; the "
    text += "TSV retains full precision.\n\n"
    text += "| Scope | Exclusion | Trait | Baseline shift | After shift | "
    text += "Matched shift change | Lost pairs | BH-six before / after | BY-six before / after | "
    text += "Removed annotation mass | Removed focal mass | Removal status | "
    text += "Effect reduction (%) | Matched effect reduction (%) | "
    text += "Per-lineage mass (annotation / focal; graph status) |\n"
    text += "|---|---|---|---:|---:|---:|---:|---|---|---:|---:|---|---:|---:|---|\n"
    for original in effects or []:
        row = _display_row(original)
        fields = (
            row["scope"],
            row["scenario"],
            row["trait_key"],
            row["before_effect"],
            row["after_effect"],
            row["matched_change_effect"],
            row["lost_pairs"],
            f"{row['before_q_bh_six']} / {row['after_q_bh_six']}",
            f"{row['before_q_by_six']} / {row['after_q_by_six']}",
            row["removed_annotation_mass"],
            row["removed_trait_annotation_mass"],
            row["removal_status"],
            row["effect_reduction_percent"],
            row["matched_effect_reduction_percent"],
            "; ".join(
                f"{lineage}: {entry['removed_annotation_mass']} / "
                f"{entry['removed_trait_annotation_mass']} ({entry['graph_status']})"
                for lineage, entry in sorted(json.loads(row["lineage_accounting"]).items())
            ),
        )
        text += "| " + " | ".join(str(value).replace("|", "\\|") for value in fields) + " |\n"
    text += "\nReduction is 100 times (1 - after / before), using full-data study-mean "
    text += "shared-fraction effects; matched reduction uses only common eligible experiments. "
    text += "Undefined ratios are blank with a status in the TSV. Values are not clipped: "
    text += "above 100% means a sign reversal, not stronger evidence. Removed mass counts "
    text += "retained trait annotations over eligible pre-exclusion profiles, not cells or "
    text += "distinct taxa. Every retained feature is audited by direction in "
    text += "lineage_removed_annotation_mass.tsv.\n"
    text += "Individual lineage masses may overlap and must not be summed to obtain the "
    text += "joint union mass. Graph presence is separate from removed mass in this scope; "
    text += "a present lineage can have no retained annotation mass here. "
    text += "Before-effects within the 1e-12 tie tolerance have undefined reduction ratios.\n"
    text += "\n## Six-fraction study omissions\n\n"
    text += "| Scope | Trait | Omissions | Baseline shift | Baseline BH-six / BY-six | "
    text += "Mean shift range | BH-six q range | BY-six q range |\n"
    text += "|---|---|---:|---:|---|---|---|---|\n"
    for original in summaries:
        row = _display_row(original)
        label = str(row["trait_key"]).replace("|", "\\|")
        text += (
            f"| {row['scope']} | {label} | {row['omissions']} | "
            f"{row['baseline_effect']} | "
            f"{row['baseline_q_bh_six']} / {row['baseline_q_by_six']} | "
            f"{row['minimum_effect']} to {row['maximum_effect']} | "
            f"{row['minimum_q_bh_six']} to {row['maximum_q_bh_six']} | "
            f"{row['minimum_q_by_six']} to {row['maximum_q_by_six']} |\n"
        )
    text += "\nRanges summarize dependent omissions; they are not independent replications "
    text += "or across-omission error control. These are annotation-fraction differences, "
    text += "not measurements of physiological activity.\n"
    return text


def run_extension(core_manifest: Path, output: Path, noel_manifest: Path | None = None) -> Path:
    """Execute the fixed sidecar protocol on one verified core in a fresh directory.

    Scope is read from actual core execution through the public workflow API,
    never guessed from matrix dimensions or from the current primary default.
    """
    from mtp.workflow.operations import run_strain_rollup

    if output.exists():
        raise FileExistsError(f"Extension output already exists: {output}")
    bundle = open_core_bundle(core_manifest)
    record = output / "core_consumption.json"
    preflight_consumer_outputs(bundle, record, output_directory=output)
    protocol = load_protocol()
    inputs = {
        name: bundle.path(name)
        for name in (
            "count_profiles",
            "experiments",
            "reported_taxa",
            "taxa_matrix",
            "kg_edges",
            "kg_nodes",
        )
    }
    if noel_manifest is not None:
        native_root = noel_manifest.resolve().parent
        destination = output.resolve()
        if destination.is_relative_to(native_root) or native_root.is_relative_to(destination):
            raise ValueError("Extension output must be outside the immutable Noel arm bundle")
        inputs["noel_manifest"] = noel_manifest.resolve()
        inputs["native_summary"] = native_summary_source(core_manifest, noel_manifest)
    input_hashes = {name: file_sha256(path) for name, path in inputs.items()}
    native_text = (
        inputs["native_summary"].read_text(encoding="utf-8") if "native_summary" in inputs else None
    )
    protocol_sha = PROTOCOL_SHA256
    producer_hashes = implementation_hashes()
    lock_sha = file_sha256(REPO_ROOT / "uv.lock")
    started_at = datetime.now(UTC).isoformat()
    base = load_base(
        *(
            inputs[name]
            for name in ("count_profiles", "experiments", "reported_taxa", "taxa_matrix")
        )
    )
    scope = run_strain_rollup(bundle.run)
    graph = load_graph(
        inputs["kg_edges"],
        inputs["kg_nodes"],
        set(base.columns),
        set(base.traits.values()),
        strain_rollup=scope,
    )
    ledger = reconcile(base, graph)
    if ledger.mismatches:
        raise ValueError(
            f"{len(ledger.mismatches)} source/profile reconciliation mismatches; no sidecar written"
        )
    candidates, panel = lineage_panel(base, graph, ledger, protocol)
    removed_mass: list[Record] = []
    effects, ambiguities = lineage_effects(
        base, graph, ledger, panel, protocol, removed_mass_rows=removed_mass
    )
    baseline: list[Record] = []
    omissions: list[Record] = []
    for analysis_scope in protocol["scopes"]:
        keys = scope_keys(base, analysis_scope)
        baseline.extend(
            {"scope": analysis_scope["id"], **row}
            for row in fraction_family(base, base.profiles, keys)
        )
        omissions.extend({"scope": analysis_scope["id"], **row} for row in study_loso(base, keys))
    summaries = loso_summary(baseline, omissions)
    if input_hashes != {name: file_sha256(path) for name, path in inputs.items()}:
        raise ValueError("Core inputs changed during lineage-extension analysis")
    if protocol_sha != file_sha256(PROTOCOL):
        raise ValueError("Extension protocol changed during analysis")
    if producer_hashes != implementation_hashes() or lock_sha != file_sha256(REPO_ROOT / "uv.lock"):
        raise ValueError("Extension implementation changed during analysis")
    output.mkdir(parents=True, exist_ok=False)
    eligible = set(eligible_keys(base.profiles, sorted(base.profiles)))
    input_eligibility = [
        {
            **row,
            "status": "zero_retained_total"
            if row["status"] == "retained"
            and (row["study_id"], row["experiment_id"], row["direction"]) not in eligible
            else row["status"],
            "orientation_status": str(bundle.run.direction_basis),
        }
        for row in base.eligibility
    ]
    tables: dict[str, tuple[list[str], list[Record]]] = {
        "lineage_candidate_ranking.tsv": (
            [
                "scope",
                "direction",
                "trait_key",
                "rank",
                "taxon",
                "name",
                "count",
                "share",
                "studies",
                "status",
                "eligible_rank",
                "selected",
                "selection_status",
            ],
            candidates,
        ),
        "lineage_extension_panel.tsv": (PANEL_FIELDS, panel),
        "lineage_fraction_effects.tsv": (
            [
                *PANEL_FIELDS,
                "trait_key",
                "family_size",
                *comparison({}, {}),
                *REMOVAL_FIELDS,
                *REDUCTION_FIELDS,
                *(
                    f"{population}_{method}"
                    for population in (
                        "before",
                        "after",
                        "before_matched",
                        "after_matched",
                        "matched_change",
                    )
                    for method in ("q_bh_six", "q_by_six")
                ),
            ],
            effects,
        ),
        "lineage_removed_annotation_mass.tsv": ([*PANEL_FIELDS, *MASS_FIELDS], removed_mass),
        "ambiguous_retained_sources.tsv": (
            ["scope", "scenario", "reported_taxon", "retained_source_taxon", "reason"],
            ambiguities,
        ),
        "baseline_fractions.tsv": (["scope", *FRACTION_FIELDS], baseline),
        "fraction_study_loso.tsv.gz": (
            ["scope", "omitted_study", *FRACTION_FIELDS, *BASELINE_FIELDS],
            omissions,
        ),
        "fraction_loso_summary.tsv": (list(summaries[0]), summaries),
        "input_eligibility.tsv": (
            ["study_id", "experiment_id", "direction", "status", "orientation_status"],
            input_eligibility,
        ),
    }
    counts = {
        name: write_table(output / name, fields, rows) for name, (fields, rows) in tables.items()
    }
    readme = output / "README.md"
    readme.write_text(
        "# Additive lineage and study-omission sensitivity\n\n"
        f"Protocol `{protocol['version']}`; population `{bundle.run.contrast_population}`; "
        f"direction `{bundle.run.direction_basis}`; strain rollup `{scope}`.\n\n"
        "The native analysis uses the declared 190-native-v2 panel and native_summary.md. "
        "Read this sidecar alongside that run's declared-panel results. "
        "lineage_extension_panel.tsv lists the additional rule-selected exclusions; "
        "lineage_candidate_ranking.tsv retains the ranking and unresolved support.\n\n"
        f"Selected {len(panel)} scope-lineage combinations; wrote {len(omissions)} "
        "study-omission/trait rows. Each omission removes all experiments from one study "
        "and recomputes the complete six-test family. Zero-total pairs are excluded.\n\n"
        "Effects are increased-minus-decreased annotation-fraction differences. "
        "BH-six/BY-six values are conditional sensitivity diagnostics; they do not control "
        "selection of the six traits or searches across omissions, scopes and lineages. "
        "LOSO ranges are descriptive and dependent. These results concern annotation "
        "alignment, not measured microbial activity or causal disease effects.\n",
        encoding="utf-8",
    )
    outputs = [*(output / name for name in tables), readme]
    combined = output / "influence_extension_summary.md"
    combined.write_text(
        combined_summary(native_text, panel, summaries, effects),
        encoding="utf-8",
    )
    outputs.append(combined)
    manifest = output / "extension_provenance.json"
    manifest.write_text(
        json.dumps(
            {
                "version": protocol["version"],
                "completed": True,
                "started_at": started_at,
                "ended_at": datetime.now(UTC).isoformat(),
                "protocol": protocol,
                "protocol_sha256": protocol_sha,
                "core_run_id": bundle.run.run_id,
                "core_manifest_sha256": file_sha256(bundle.manifest),
                "contrast_population": bundle.run.contrast_population,
                "direction_basis": bundle.run.direction_basis,
                "strain_rollup": scope,
                "input_sha256": input_hashes,
                "retained_features": base.columns,
                "focal_features": base.traits,
                "reconciliation_mismatches": 0,
                "rows": counts,
                "code_sha256": producer_hashes,
                "uv_lock_sha256": lock_sha,
                "outputs": {path.name: file_sha256(path) for path in outputs},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    write_consumer_record(
        bundle,
        record,
        consumer="table1-influence-extension-v1",
        outputs=[*outputs, manifest],
        additional_products=("reported_taxa", "kg_edges", "kg_nodes"),
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--noel-run", type=Path, help="Include the verified same-core native report"
    )
    add_temporary_output_option(parser)
    args = parser.parse_args(argv)
    require_durable_output(parser, args.output_dir, allowed=args.allow_temporary_output)
    print(run_extension(args.core_run, args.output_dir, args.noel_run))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
