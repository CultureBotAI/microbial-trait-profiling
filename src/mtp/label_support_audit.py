"""Sidecar audit of annotation carriers, sources, co-assertions and support rules (#773)."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict, deque
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from mtp import taxa_trait_matrix as matrix_api
from mtp import trait_guilds
from mtp.count_input import CountTableReader, profile_count_reader, taxa_trait_count_reader
from mtp.paths import KG_MICROBE_DIR, REPO_ROOT, file_sha256
from mtp.taxa_trait_matrix import (
    EDGES_FILENAME,
    NODES_FILENAME,
    STRAIN_PREFIX,
    SUBCLASS_OF,
    TAXON_PREFIX,
    TRAIT_PREFIX,
    make_trait_key,
    resolve_kg_file,
)
from mtp.trait_guilds import MEDIADIVE_MARKER, OXYGEN_GUILDS, TABLE1_TRAIT_LABELS, TABLE1_TRAITS
from mtp.workflow import kg_inputs

Scope = Literal["off", "phenotype"]
Row = dict[str, Any]
MADIN_SOURCE = "infores:madin_etal"
METADATA_ONLY_MADIN_SOURCES = frozenset({"gold", "patric"})
AUDIT_TRAITS = {
    term
    if term.startswith(TRAIT_PREFIX)
    else make_trait_key("METPO:2000103", term): TABLE1_TRAIT_LABELS[name].rsplit(" (", 1)[0]
    for name, terms in TABLE1_TRAITS.items()
    for term in terms
}
FOCAL_TRAITS = tuple(AUDIT_TRAITS)
AUDIT_TRAITS.update(
    {
        "METPO:1000612": "facultative oxygen preference",
        "METPO:1000602": "aerobic",
        "METPO:1000603": "anaerobic",
        "METPO:1000607": "obligately anaerobic",
        "METPO:1000611": "strictly anaerobic",
        "METPO:1000606": "obligately aerobic",
        "METPO:1000604": "microaerophilic",
        "METPO:1000681": "rod shaped",
        "METPO:1000667": "bacillus shaped",
    }
)
FOCAL_TAXA = {
    "Faecalibacterium": "NCBITaxon:216851",
    "Roseburia": "NCBITaxon:841",
    "Blautia": "NCBITaxon:572511",
    "Anaerostipes": "NCBITaxon:207244",
    "Bacteroides": "NCBITaxon:816",
    "Ruminococcus": "NCBITaxon:1263",
    "Lachnospira": "NCBITaxon:28050",
    "Streptococcus": "NCBITaxon:1301",
    "Staphylococcus": "NCBITaxon:1279",
    "Lachnospiraceae": "NCBITaxon:186803",
    "Oscillospiraceae": "NCBITaxon:216572",
}
COASSERTION_PAIRS = (
    (
        "aerobe_guild_and_anaerobe_guild",
        OXYGEN_GUILDS["aerobe"],
        OXYGEN_GUILDS["anaerobe"],
    ),
    (
        "aerobic_and_anaerobic",
        ("METPO:1000602", "METPO:1000606"),
        ("METPO:1000603", "METPO:1000607", "METPO:1000611"),
    ),
    (
        "coccus_sphere_and_rod_bacillus",
        ("METPO:1000668", "METPO:1000683"),
        ("METPO:1000681", "METPO:1000667"),
    ),
)
MADIN_RULES = {
    "METPO:1000608": ("metabolism", "facultative"),
    "METPO:1000602": ("metabolism", "aerobic"),
    "METPO:1000603": ("metabolism", "anaerobic"),
    "METPO:1000606": ("metabolism", "obligate aerobic"),
    "METPO:1000607": ("metabolism", "obligate anaerobic"),
    "METPO:1000604": ("metabolism", "microaerophilic"),
    "METPO:1000668": ("cell_shape", "coccus"),
    "METPO:1000667": ("cell_shape", "bacillus"),
    "METPO:1000623": ("range_salinity", "moderate-halophilic"),
    "METPO:2000103|GO:0030245": ("pathways", "cellulose_degradation"),
    "METPO:2000103|pathways:arsenate_reduction": ("pathways", "arsenate_reduction"),
}
PROFILE_KEYS = ["study_id", "experiment_id", "inc_dec"]


def _require_columns(reader: csv.DictReader[str], columns: set[str], path: Path) -> None:
    missing = columns - set(reader.fieldnames or [])
    if missing:
        raise ValueError(f"{path.name}: missing required columns {sorted(missing)}")


def load_kg_provenance_graph(
    edges_path: Path,
    nodes_path: Path,
    *,
    strain_rollup: Scope = "off",
) -> dict[str, Any]:
    """Read named KG columns; phenotype requires the shared scope API from #828."""
    parents: defaultdict[str, set[str]] = defaultdict(set)
    traits: dict[str, set[str]] = {}
    sources: defaultdict[str, defaultdict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    strain_taxa: defaultdict[str, set[str]] = defaultdict(set)
    strain_traits: defaultdict[str, set[str]] = defaultdict(set)
    strain_sources: defaultdict[str, defaultdict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    serialization: dict[str, object] = {}
    terms = matrix_api.read_assertion_term_metadata(nodes_path, edges_path)
    for row in kg_inputs.iter_kg_edges(edges_path, serialization):
        if "primary_knowledge_source" not in row:
            raise ValueError(f"{edges_path.name}: missing required column primary_knowledge_source")
        subject, predicate, obj = (row[key].strip() for key in ("subject", "predicate", "object"))
        if predicate == SUBCLASS_OF and obj.startswith(TAXON_PREFIX):
            if subject.startswith(TAXON_PREFIX):
                parents[subject].add(obj)
            elif strain_rollup == "phenotype" and subject.startswith(STRAIN_PREFIX):
                strain_taxa[subject].add(obj)
            continue
        if not (obj.startswith(TRAIT_PREFIX) or predicate.startswith(TRAIT_PREFIX)):
            continue
        if not subject.startswith((TAXON_PREFIX, STRAIN_PREFIX)):
            continue
        key = matrix_api.assertion_trait_key(row, terms)
        if MEDIADIVE_MARKER in key.lower():
            continue
        primary = {
            value.strip() for value in row["primary_knowledge_source"].split("|") if value.strip()
        } or {"unknown"}
        if subject.startswith(TAXON_PREFIX):
            traits.setdefault(subject, set()).add(key)
            sources[key][subject].update(primary)
        elif strain_rollup == "phenotype" and subject.startswith(STRAIN_PREFIX):
            strain_traits[subject].add(key)
            strain_sources[key][subject].update(primary)
    if not serialization["records"]:
        raise ValueError(f"{edges_path.name}: empty KG edge export")
    if strain_rollup not in ("off", "phenotype"):
        raise ValueError(f"Unsupported strain scope: {strain_rollup}")
    if strain_rollup == "phenotype":
        traits = matrix_api.scoped_taxon_traits(
            traits, dict(strain_taxa), dict(strain_traits), strain_rollup
        )
        for key, carriers in strain_sources.items():
            if matrix_api.strain_trait_rolls_up(key, strain_rollup):
                for strain, provenance in carriers.items():
                    for taxon in strain_taxa.get(strain, ()):
                        sources[key][taxon].update(provenance)
    children: defaultdict[str, set[str]] = defaultdict(set)
    for child, ancestors in parents.items():
        for parent in ancestors:
            children[parent].add(child)
    names = {}
    with nodes_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_columns(reader, {"id", "name"}, nodes_path)
        for row in reader:
            if row["id"] in traits or row["id"] in FOCAL_TAXA.values():
                names[row["id"]] = row["name"]
    return {
        "traits": traits,
        "sources": dict(sources),
        "children": dict(children),
        "names": names,
        "edge_serialization": serialization,
    }


def load_madin_subsources(path: Path | None) -> dict[str, dict[str, set[str]]]:
    """Resolve declared Madin field/value mappings; omitted input remains unknown."""
    result: defaultdict[str, defaultdict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    if path is None:
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        _require_columns(reader, {"tax_id", "data_source"}, path)
        if not any(field in (reader.fieldnames or []) for field, _ in MADIN_RULES.values()):
            raise ValueError("Madin table contains none of the declared trait fields")
        for row in reader:
            taxon = (row["tax_id"] or "").strip().removeprefix(TAXON_PREFIX)
            source = (row["data_source"] or "").strip().lower()
            if not taxon or not source:
                continue
            for key, (field, value) in MADIN_RULES.items():
                tokens = {part.strip().lower() for part in re.split(r"[;,|]", row.get(field) or "")}
                if value in tokens:
                    result[key][f"{TAXON_PREFIX}{taxon}"].add(source)
    return {key: dict(values) for key, values in result.items()}


def metadata_source_status(primary: set[str], madin: set[str], unresolved_madin: int) -> str:
    """Source labels do not establish whether a trait was measured experimentally."""
    if not primary or "unknown" in primary:
        return "unknown"
    if primary - {MADIN_SOURCE}:
        return "includes_other_sources"
    if unresolved_madin or not madin:
        return "unknown"
    return "gold_patric_only" if madin <= METADATA_ONLY_MADIN_SOURCES else "includes_other_sources"


def _read_counts(reader: CountTableReader) -> pd.DataFrame:
    records = []
    for row in reader:
        record: Row = dict(zip(reader.header, row.values, strict=True))
        record.update(zip(reader.count_columns, row.counts, strict=True))
        records.append(record)
    if not records:
        raise ValueError(f"No count rows in {reader.path.name}")
    return pd.DataFrame(records)


def bind_occurrences(
    taxa_path: Path, matrix: pd.DataFrame, profiles: pd.DataFrame
) -> tuple[pd.DataFrame, Row]:
    """Restrict to exact profile keys, then reconcile all count columns before analysis."""
    occurrences = pd.read_csv(taxa_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"study_id", "experiment_id", "direction", "taxon_id"}
    if missing := required - set(occurrences):
        raise ValueError(f"Occurrence table missing columns: {sorted(missing)}")
    occurrences["inc_dec"] = (
        occurrences["direction"]
        .str.strip()
        .str.lower()
        .map({"increased": "inc", "decreased": "dec"})
    )
    if not profiles["inc_dec"].isin(["inc", "dec"]).all():
        raise ValueError("Profiles contain unsupported directions")
    keys = set(profiles[PROFILE_KEYS].itertuples(index=False, name=None))
    before = len(occurrences)
    occurrences = occurrences.loc[
        [key in keys for key in occurrences[PROFILE_KEYS].itertuples(index=False, name=None)]
    ].copy()
    if not occurrences["taxon_id"].str.startswith(TAXON_PREFIX).all():
        raise ValueError("Bound occurrences must have explicit NCBITaxon identifiers")
    occurrences = occurrences.drop_duplicates([*PROFILE_KEYS, "taxon_id"])
    observed_keys = set(occurrences[PROFILE_KEYS].itertuples(index=False, name=None))
    if keys != observed_keys:
        raise ValueError("Occurrence table is missing a profile side; missing is not zero")
    counts = matrix.set_index("taxon_id")
    features = list(profiles.columns[3:])
    if missing := set(features) - set(counts):
        raise ValueError(f"Profile features missing from taxon matrix: {sorted(missing)}")
    bound_counts = counts.loc[counts.index.isin(occurrences["taxon_id"])]
    retained = {
        key
        for key in matrix.columns[3:]
        if MEDIADIVE_MARKER not in key.lower() and bound_counts[key].gt(0).any()
    }
    if missing := retained - set(features):
        raise ValueError(f"Profile table missing retained matrix features: {sorted(missing)}")
    expected = profiles.set_index(PROFILE_KEYS)[features]
    for key, group in occurrences.groupby(PROFILE_KEYS, sort=False):
        # Unmapped taxa remain in the reported-taxon denominator, with no matrix annotations.
        known = group["taxon_id"].isin(counts.index)
        total = counts.loc[group.loc[known, "taxon_id"], features].sum(axis=0)
        expected_row = expected.xs(key)
        if not isinstance(expected_row, pd.Series):
            raise ValueError(f"Profile keys must be unique: {key}")
        if not total.eq(expected_row).all():
            raise ValueError(
                f"Occurrence/matrix counts disagree with profile {key}; check scope/orientation"
            )
    facts = {
        "input_occurrence_rows": before,
        "bound_distinct_occurrences": len(occurrences),
        "bound_profiles": len(profiles),
        "profile_counts_reconciled": True,
        "unmapped_occurrences": int((~occurrences["taxon_id"].isin(counts.index)).sum()),
        "population": "exact supplied profile keys; no independent host or orientation claim",
    }
    return occurrences, facts


def _subtree(taxon: str, children: dict[str, set[str]]) -> set[str]:
    found: set[str] = set()
    pending = deque([taxon])
    while pending:
        node = pending.popleft()
        if node not in found:
            found.add(node)
            pending.extend(children.get(node, ()))
    return found


def build_label_support_tables(
    graph: dict[str, Any],
    matrix: pd.DataFrame,
    occurrences: pd.DataFrame,
    madin: dict[str, dict[str, set[str]]],
) -> tuple[list[Row], list[Row], list[Row]]:
    """Reconcile carrier counts and export source evidence without changing the matrix."""
    support, focal, evidence = [], [], []
    counts = matrix.set_index("taxon_id")
    features = [key for key in AUDIT_TRAITS if key in counts]
    annotated = {taxon for taxon, traits in graph["traits"].items() if traits}
    feature_carriers = {key: set(graph["sources"].get(key, {})) for key in features}
    occurrence_counts = occurrences["taxon_id"].value_counts()
    focal_names = {taxon: name for name, taxon in FOCAL_TAXA.items()}
    for taxon in sorted(set(counts.index) | set(focal_names)):
        subtree = _subtree(taxon, graph["children"])
        annotated_count = len(subtree & annotated)
        for key in features:
            carriers = subtree & feature_carriers[key]
            if taxon in counts.index and counts.loc[taxon, key] != len(carriers):
                raise ValueError(
                    f"KG carrier count disagrees with matrix: {taxon} / {key}; check strain scope"
                )
            primary: set[str] = set()
            sub_sources: set[str] = set()
            unresolved = 0
            for carrier in sorted(carriers):
                pks = set(graph["sources"][key][carrier])
                subs = madin.get(key, {}).get(carrier, set()) if MADIN_SOURCE in pks else set()
                primary.update(pks)
                sub_sources.update(subs)
                unresolved += int(MADIN_SOURCE in pks and not subs)
                if taxon in counts.index:
                    evidence.append(
                        {
                            "reported_taxon_id": taxon,
                            "trait_key": key,
                            "carrier_taxon_id": carrier,
                            "carrier_name": graph["names"].get(carrier, "unknown"),
                            "primary_sources": ";".join(sorted(pks)),
                            "madin_subsources": ";".join(sorted(subs)),
                            "madin_resolution": "resolved"
                            if subs
                            else ("unknown" if MADIN_SOURCE in pks else "not_applicable"),
                        }
                    )
            status = (
                metadata_source_status(primary, sub_sources, unresolved)
                if carriers
                else "no_annotation"
            )
            record = {
                "taxon_id": taxon,
                "taxon_name": str(counts.loc[taxon, "taxon_name"])
                if taxon in counts.index
                else focal_names[taxon],
                "taxonomic_rank": str(counts.loc[taxon, "taxonomic_rank"])
                if taxon in counts.index
                else "unknown",
                "trait_key": key,
                "trait_label": AUDIT_TRAITS[key],
                "carrier_count": len(carriers),
                "carrier_at_reported_taxon": taxon in carriers,
                "inherited_carrier_count": len(carriers - {taxon}),
                "annotated_subtree_taxa": annotated_count,
                "annotated_subtree_fraction": len(carriers) / annotated_count
                if annotated_count
                else None,
                "single_carrier_flag": len(carriers) == 1,
                "source_status": status,
                "unresolved_madin_carriers": unresolved,
                "primary_sources": ";".join(sorted(primary)),
                "madin_subsources": ";".join(sorted(sub_sources)),
                "signature_occurrences": int(occurrence_counts.get(taxon, 0)),
            }
            if taxon in counts.index and carriers:
                support.append(record)
            if taxon in focal_names:
                focal.append({**record, "in_supplied_matrix": taxon in counts.index})
    return support, focal, evidence


def build_support_threshold_sensitivity(
    support: list[Row],
    occurrences: pd.DataFrame,
    profiles: pd.DataFrame,
) -> list[Row]:
    """Descriptive rule sensitivities and nominal study sign tests; no primary calls."""
    sizes = occurrences.groupby(PROFILE_KEYS).size()
    directions = profiles.groupby(PROFILE_KEYS[:2])["inc_dec"].nunique()
    paired = directions[directions == 2].index
    if not len(paired):
        raise ValueError("No paired experiments in the supplied profiles")
    native = profiles.set_index(PROFILE_KEYS)
    retained = [key for key in native.columns if MEDIADIVE_MARKER not in key.lower()]
    total_annotations = native[retained].sum(axis=1)
    input_pairs = len(paired)
    positive_totals = total_annotations.unstack("inc_dec").reindex(paired)
    paired = paired[positive_totals[["inc", "dec"]].gt(0).all(axis=1)]
    excluded_zero_pairs = input_pairs - len(paired)
    pair_keys = pd.MultiIndex.from_frame(occurrences[PROFILE_KEYS[:2]])
    bound_taxa = set(occurrences.loc[pair_keys.isin(paired), "taxon_id"])
    results = []
    rules = (
        "native_annotation_fraction",
        "carriers_per_reported_taxon",
        "binary_presence",
        "min_two_carriers",
        "min_subtree_fraction_0.10",
        "subtree_fraction",
        "exclude_known_gold_patric_only",
    )
    for feature in FOCAL_TRAITS:
        if feature not in profiles:
            continue
        rows = {
            row["taxon_id"]: row
            for row in support
            if row["trait_key"] == feature and row["taxon_id"] in bound_taxa
        }
        for rule in rules:
            weights = {}
            for taxon, row in rows.items():
                count, fraction = row["carrier_count"], row["annotated_subtree_fraction"]
                weight = float(count)
                if rule == "binary_presence":
                    weight = 1.0
                elif rule == "min_two_carriers":
                    weight = float(count >= 2)
                elif rule == "min_subtree_fraction_0.10":
                    weight = float(fraction >= 0.10)
                elif rule == "subtree_fraction":
                    weight = float(fraction)
                elif (
                    rule == "exclude_known_gold_patric_only"
                    and row["source_status"] == "gold_patric_only"
                ):
                    weight = 0.0
                weights[taxon] = weight
            if rule == "native_annotation_fraction":
                scores = native[feature] / total_annotations.replace(0, np.nan)
                units = "fraction_of_all_retained_annotation_counts"
            else:
                weighted = occurrences.assign(
                    weight=occurrences["taxon_id"].map(weights).fillna(0.0)
                )
                scores = weighted.groupby(PROFILE_KEYS)["weight"].sum() / sizes
                units = (
                    "carriers_per_reported_taxon"
                    if rule in ("carriers_per_reported_taxon", "exclude_known_gold_patric_only")
                    else "mean_trait_weight_per_reported_taxon"
                )
            values = scores.unstack("inc_dec").reindex(paired)
            eligible = values[["inc", "dec"]].notna().all(axis=1)
            deltas = values.loc[eligible, "inc"] - values.loc[eligible, "dec"]
            study_means = deltas.groupby(level="study_id").mean()
            positive = int(study_means.gt(1e-12).sum())
            negative = int(study_means.lt(-1e-12).sum())
            active = positive + negative
            results.append(
                {
                    "trait_key": feature,
                    "rule_id": rule,
                    "effect_units": units,
                    "n_experiments": len(deltas),
                    "n_input_paired_experiments": input_pairs,
                    "n_excluded_zero_annotation_pairs": excluded_zero_pairs,
                    "comparison_basis": "paired_positive_total_retained_annotation_profiles",
                    "n_studies": len(study_means),
                    "n_positive_studies": positive,
                    "n_negative_studies": negative,
                    "n_tied_studies": len(study_means) - active,
                    "study_mean_increased_minus_decreased": float(study_means.mean())
                    if len(study_means)
                    else None,
                    "study_sign_p_value": float(binomtest(positive, active).pvalue)
                    if active
                    else None,
                    "test_status": "tested" if active else "no_nonzero_study_differences",
                    "multiple_testing": "nominal_exploratory_p_value_no_significance_call",
                    "aggregation": "study_mean_all_eligible_pairs_including_zero_differences",
                    "source_unknown_taxa": sum(
                        row["source_status"] == "unknown" for row in rows.values()
                    ),
                }
            )
    return results


def build_coassertion_census(matrix: pd.DataFrame, occurrences: pd.DataFrame) -> list[Row]:
    """Co-assertion on a subtree is not proof of contradictory physiology."""
    counts = matrix.set_index("taxon_id")
    occurrences = occurrences.loc[occurrences["taxon_id"].isin(counts.index)].copy()
    active = counts.loc[counts.index.isin(occurrences["taxon_id"])].copy()
    active["rank_group"] = active["taxonomic_rank"].where(
        active["taxonomic_rank"].isin(["species", "genus", "family"]), "other_or_unknown"
    )
    results = []
    for label, side_a, side_b in COASSERTION_PAIRS:
        missing = sorted((set(side_a) | set(side_b)) - set(counts))
        for rank in ("species", "genus", "family", "other_or_unknown", "ALL"):
            selected = active if rank == "ALL" else active.loc[active["rank_group"] == rank]
            row: Row = {
                "label_pair": label,
                "rank": rank,
                "unique_reported_taxa": len(selected),
                "status": "unavailable_features" if missing else "assessed",
                "missing_features": ";".join(missing),
            }
            if not missing:
                a, b = selected[list(side_a)].sum(axis=1), selected[list(side_b)].sum(axis=1)
                either, both = a.gt(0) | b.gt(0), a.gt(0) & b.gt(0)
                occ = occurrences.loc[occurrences["taxon_id"].isin(selected.index), "taxon_id"]
                occ_either, occ_both = int(occ.map(either).sum()), int(occ.map(both).sum())
                a_total, b_total = float(occ.map(a).sum()), float(occ.map(b).sum())
                row.update(
                    {
                        "taxa_with_either": int(either.sum()),
                        "taxa_with_both": int(both.sum()),
                        "fraction_both_of_either": float(both.sum() / either.sum())
                        if either.sum()
                        else None,
                        "occurrences_with_either": occ_either,
                        "occurrences_with_both": occ_both,
                        "fraction_occurrences_both_of_either": occ_both / occ_either
                        if occ_either
                        else None,
                        "side_a_count_mass": a_total,
                        "side_b_count_mass": b_total,
                        "side_a_coasserted_mass_fraction": float(occ.map(a.where(both, 0)).sum())
                        / a_total
                        if a_total
                        else None,
                        "side_b_coasserted_mass_fraction": float(occ.map(b.where(both, 0)).sum())
                        / b_total
                        if b_total
                        else None,
                    }
                )
            results.append(row)
    return results


def _inventory(path: Path) -> Row:
    resolved = path.resolve()
    display = (
        str(resolved.relative_to(REPO_ROOT)) if resolved.is_relative_to(REPO_ROOT) else path.name
    )
    return {"path": display, "sha256": file_sha256(path), "size_bytes": path.stat().st_size}


def run_label_support_audit(
    *,
    taxa_path: Path,
    matrix_path: Path,
    profiles_path: Path,
    edges_path: Path,
    nodes_path: Path,
    output_dir: Path,
    madin_csv: Path | None = None,
    strain_rollup: Scope = "off",
) -> Row:
    """Write sidecar products only, after all source reconciliations pass."""
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError("Choose an empty output directory for the label-support sidecar")
    inputs = {
        "taxa": taxa_path,
        "matrix": matrix_path,
        "profiles": profiles_path,
        "kg_edges": edges_path,
        "kg_nodes": nodes_path,
    }
    if madin_csv is not None:
        inputs["madin"] = madin_csv
    provenance_inputs = {name: _inventory(path) for name, path in inputs.items()}
    matrix = _read_counts(taxa_trait_count_reader(matrix_path))
    profiles = _read_counts(profile_count_reader(profiles_path))
    occurrences, scope_facts = bind_occurrences(taxa_path, matrix, profiles)
    graph = load_kg_provenance_graph(edges_path, nodes_path, strain_rollup=strain_rollup)
    madin = load_madin_subsources(madin_csv)
    support, focal, evidence = build_label_support_tables(graph, matrix, occurrences, madin)
    sensitivity = build_support_threshold_sensitivity(support, occurrences, profiles)
    census = build_coassertion_census(matrix, occurrences)
    tables = {
        "taxon_trait_label_support.tsv": (support, ["taxon_id", "trait_key", "carrier_count"]),
        "focal_taxa_trait_support.tsv": (focal, ["taxon_id", "trait_key", "carrier_count"]),
        "carrier_source_evidence.tsv": (
            evidence,
            ["reported_taxon_id", "trait_key", "carrier_taxon_id"],
        ),
        "support_threshold_sensitivity.tsv": (
            sensitivity,
            ["trait_key", "rule_id", "study_sign_p_value"],
        ),
        "coasserted_label_census.tsv": (census, ["label_pair", "rank", "status"]),
    }
    if provenance_inputs != {name: _inventory(path) for name, path in inputs.items()}:
        raise ValueError("Audit inputs changed during execution")
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, (rows, empty_columns) in tables.items():
        pd.DataFrame(rows, columns=None if rows else empty_columns).to_csv(
            output_dir / name, sep="\t", index=False
        )
    result = {
        "schema_version": "1.0",
        "created_at": datetime.now(UTC).isoformat(),
        "analysis": "label_support_sidecar",
        "strain_rollup": strain_rollup,
        "inputs": provenance_inputs,
        "scope": scope_facts,
        "madin_status": "supplied_hashed_unpinned"
        if madin_csv is not None
        else "not_supplied_source_support_unknown",
        "madin_rules": MADIN_RULES,
        "producer": _inventory(Path(__file__)),
        "matrix_policy_code": _inventory(Path(matrix_api.__file__)),
        "kg_parser_code": _inventory(Path(kg_inputs.__file__)),
        "trait_definitions_code": _inventory(Path(trait_guilds.__file__)),
        "kg_edge_serialization": graph["edge_serialization"],
        "environment_lock": _inventory(REPO_ROOT / "uv.lock"),
        "verification_status": "hashes_and_counts_reconciled_not_core_execution_verified",
        "support_cells": len(support),
        "support_table_population": "all_supplied_matrix_taxa",
        "sensitivity_population": "paired_positive_total_retained_annotation_profiles",
        "source_status_counts": dict(Counter(row["source_status"] for row in support)),
        "sensitivity_family": {
            "traits": [key for key in FOCAL_TRAITS if key in profiles],
            "tested_rules": sorted({row["rule_id"] for row in sensitivity}),
            "selection": "exploratory focal traits; nominal p-values only",
        },
        "limitations": [
            "Annotations are not measured physiology.",
            "A single carrier need not mean a single source record.",
            "Co-assertion can reflect different descendants or compatible trait semantics.",
            "Madin source labels do not establish experimental validation.",
            "Threshold and fraction rules are sensitivities, not changes to primary arm calls.",
        ],
        "outputs": {name: _inventory(output_dir / name) for name in tables},
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxa", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--kg-dir", type=Path, default=KG_MICROBE_DIR)
    parser.add_argument("--edges", type=Path)
    parser.add_argument("--nodes", type=Path)
    parser.add_argument("--madin-csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--strain-rollup", choices=("off", "phenotype"), required=True)
    args = parser.parse_args(argv)
    run_label_support_audit(
        taxa_path=args.taxa,
        matrix_path=args.matrix,
        profiles_path=args.profiles,
        edges_path=resolve_kg_file(EDGES_FILENAME, args.kg_dir, args.edges),
        nodes_path=resolve_kg_file(NODES_FILENAME, args.kg_dir, args.nodes),
        madin_csv=args.madin_csv,
        strain_rollup=args.strain_rollup,
        output_dir=args.output_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
