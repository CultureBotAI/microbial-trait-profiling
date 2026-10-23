"""Additive lineage-extension and six-fraction study-LOSO diagnostics (#771).

The source population, lineage ranking and complete six-member corrections are
explicit. None of these exploratory sensitivities changes the native arm plan.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path

from mtp.contribution_inputs import Base, Graph, Key, encoded, load_base, load_graph
from mtp.paths import REPO_ROOT, file_sha256
from mtp.taxon_contributions import Ledger, Record, reconcile, write_table
from mtp.taxon_influence import (
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

PROTOCOL = REPO_ROOT / "methods/table1_influence_extension_protocol.json"
PANEL_FIELDS = ["scope", "rank", "lineage", "name", "scenario", "selection_reasons"]
FRACTION_FIELDS = ["trait_key", "family_size", *describe({}), "q_bh_six", "q_by_six"]
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
    return [
        {"omitted_study": study, **row}
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
        keys = eligible_keys(base.profiles, scope_keys(base, scope))
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
                    "name": graph.names.get(lineage, base.names.get(lineage, "")),
                    "scenario": f"extension_exclude_{rank}_{lineage}",
                    "selection_reasons": encoded(sorted(why)),
                }
            )
    return candidates, selected


def lineage_effects(
    base: Base, graph: Graph, ledger: Ledger, selected: list[Record], protocol: Record
) -> tuple[list[Record], list[Record]]:
    """Recompute denominator and eligibility; compare identical support separately."""
    scopes = {scope["id"]: scope for scope in protocol["scopes"]}
    effects: list[Record] = []
    ambiguities: list[Record] = []
    for scenario in selected:
        keys = eligible_keys(base.profiles, scope_keys(base, scopes[scenario["scope"]]))
        remaining, _, ambiguous = removal(base, graph, ledger, scenario["lineage"])
        family = []
        for trait in _focal_traits(base):
            before = pair_values(base.profiles, keys, trait)
            after = pair_values(remaining, keys, trait)
            family.append(
                {**scenario, "trait_key": trait, "family_size": 6, **comparison(before, after)}
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
        original_effect = original["effect"]
        result.append(
            {
                "scope": original["scope"],
                "trait_key": original["trait_key"],
                "baseline_effect": original_effect,
                "baseline_q_by_six": original["q_by_six"],
                "omissions": len(rows),
                "testable_omissions": len(qvalues),
                "minimum_effect": min(effects) if effects else "",
                "maximum_effect": max(effects) if effects else "",
                "minimum_q_by_six": min(qvalues) if qvalues else "",
                "maximum_q_by_six": max(qvalues) if qvalues else "",
                "mean_sign_reversals": sum(value * float(original_effect) < 0 for value in effects)
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


def combined_summary(native_text: str | None, panel: list[Record], summaries: list[Record]) -> str:
    """Show baseline text and clearly separate the additive dependent diagnostics."""
    baseline = (
        native_text
        if native_text is not None
        else ("Native baseline was not supplied. Pass --noel-run to include its verified report.\n")
    )
    text = "# Declared baseline and additive sensitivity\n\n"
    text += "The following baseline text is copied from the verified native report when supplied. "
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
    text += "## Six-fraction study omissions\n\n"
    text += "| Scope | Trait | Omissions | Mean shift range | BY-six q range |\n"
    text += "|---|---|---:|---|---|\n"
    for row in summaries:
        label = str(row["trait_key"]).replace("|", "\\|")
        text += (
            f"| {row['scope']} | {label} | {row['omissions']} | "
            f"{row['minimum_effect']} to {row['maximum_effect']} | "
            f"{row['minimum_q_by_six']} to {row['maximum_q_by_six']} |\n"
        )
    text += "\nRanges summarize dependent omissions; they are not independent replications "
    text += "or across-omission error control. These are annotation-fraction differences, "
    text += "not measurements of physiological activity.\n"
    return text


def run_extension(core_manifest: Path, output: Path, noel_manifest: Path | None = None) -> Path:
    """Execute the fixed sidecar protocol on one verified core in a fresh directory.

    The scope API is supplied by #828. Scope is read from actual core execution,
    never guessed from matrix dimensions or from the current primary default.
    """
    from mtp.workflow.operations import run_strain_rollup

    if output.exists():
        raise FileExistsError(f"Extension output already exists: {output}")
    bundle = open_core_bundle(core_manifest)
    record = output / "core_consumption.json"
    preflight_consumer_outputs(bundle, record, output_directory=output)
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if protocol["version"] != "771-extension-v1":
        raise ValueError("Unknown lineage-extension protocol version")
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
    protocol_sha = file_sha256(PROTOCOL)
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
    effects, ambiguities = lineage_effects(base, graph, ledger, panel, protocol)
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
        "ambiguous_retained_sources.tsv": (
            ["scope", "scenario", "reported_taxon", "retained_source_taxon", "reason"],
            ambiguities,
        ),
        "baseline_fractions.tsv": (["scope", *FRACTION_FIELDS], baseline),
        "fraction_study_loso.tsv.gz": (["scope", "omitted_study", *FRACTION_FIELDS], omissions),
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
        "The declared 190-native-v2 panel and native_summary.md remain the baseline panel. "
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
    combined = output / "native_summary.md"
    combined.write_text(
        combined_summary(native_text, panel, summaries),
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
        bundle, record, consumer="table1-influence-extension-v1", outputs=[*outputs, manifest]
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--noel-run", type=Path, help="Include the verified same-core native report"
    )
    args = parser.parse_args(argv)
    print(run_extension(args.core_run, args.output_dir, args.noel_run))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
