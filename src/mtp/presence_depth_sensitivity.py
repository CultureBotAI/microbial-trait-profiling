"""Additive presence and annotation-depth sensitivities of verified core profiles (#772)."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import sys
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

from mtp.contribution_inputs import rows
from mtp.count_input import profile_count_reader, taxa_trait_count_reader
from mtp.guild_maps import SIGN_TOLERANCE
from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES, HUMAN_HOST_SPECIES
from mtp.multiple_testing import POLICY_PATH
from mtp.paths import REPO_ROOT, file_sha256
from mtp.presence_depth_families import declare_and_correct
from mtp.trait_composition import compute_log_ratio_shifts
from mtp.trait_guilds import TABLE1_TRAITS, Profile, column_term_id
from mtp.trait_screen import NO_ELIGIBLE_STUDIES, run_screen
from mtp.workflow.consumers import (
    CoreBundle,
    open_core_bundle,
    preflight_consumer_outputs,
    verify_consumer_record,
    write_consumer_record,
)
from mtp.workflow.coverage import feature_filter_id
from mtp.workflow.downstream_contracts import build_analysis_basis
from mtp.workflow.operations import run_strain_rollup
from mtp.workflow.output_roots import add_temporary_output_option, require_durable_output
from mtp.workflow.runtime import repository_state

Key = tuple[str, str, str]
Pair = tuple[str, str]
Result = dict[str, Any]
PROTOCOL = "presence-depth-study-v2"
REPRESENTATIONS = ("raw", "presence", "prevalence")


@dataclass(frozen=True)
class ProfileEvidence:
    """Reconstructed counts and depth, including matched zero-annotation taxa."""

    key: Key
    raw: tuple[int, ...]
    presence: tuple[int, ...]
    reported_taxa: int
    matched_taxa: int
    unmatched_taxa: tuple[str, ...]
    unannotated_taxa: int
    depth_mean_total: float | None
    depth_primary: float | None
    depth_secondary: float | None
    ora_taxa_count: int | None = None


@dataclass(frozen=True)
class StudyEvidence:
    """One equally weighted observation per study and representation."""

    study_id: str
    representation: str
    experiments: int
    mean_differences: tuple[float, ...]
    depth_primary: float
    depth_secondary: float


def reconstruct_profiles(
    observed: Mapping[Key, Sequence[int]],
    reported: Mapping[Key, Sequence[str]],
    matrix: Mapping[str, Sequence[int]],
    columns: Sequence[str],
) -> list[ProfileEvidence]:
    """Deduplicate reported taxa and verify every retained raw count before inference."""
    if not columns or len(set(columns)) != len(columns):
        raise ValueError("Retained columns must be nonempty and unique")
    if set(observed) != set(reported):
        raise ValueError("Profile and reported-taxon keys do not reconstruct the same population")
    for vector in [*observed.values(), *matrix.values()]:
        if len(vector) != len(columns) or any(
            not isinstance(value, int) or value < 0 for value in vector
        ):
            raise ValueError("Reconstruction needs aligned nonnegative integer counts")
    result = []
    for key, observed_counts in sorted(observed.items()):
        if key[2] not in {"inc", "dec"}:
            raise ValueError(f"Unknown profile direction: {key}")
        taxa = sorted(set(reported[key]))
        matched = [matrix[taxon] for taxon in taxa if taxon in matrix]
        raw = tuple(sum(vector[i] for vector in matched) for i in range(len(columns)))
        presence = tuple(sum(vector[i] > 0 for vector in matched) for i in range(len(columns)))
        if raw != tuple(observed_counts):
            raise ValueError(f"Reported taxa do not reconstruct retained profile counts: {key}")
        mean_total = sum(raw) / len(matched) if matched else None
        result.append(
            ProfileEvidence(
                key,
                raw,
                presence,
                len(taxa),
                len(matched),
                tuple(taxon for taxon in taxa if taxon not in matrix),
                sum(sum(vector) == 0 for vector in matched),
                mean_total,
                sum(math.log10(1 + max(vector)) for vector in matched) / len(matched)
                if matched
                else None,
                math.log1p(mean_total) if mean_total is not None else None,
            )
        )
    return result


def make_study_observations(
    profiles: Sequence[ProfileEvidence], metadata: Mapping[Pair, Mapping[str, str]]
) -> tuple[list[StudyEvidence], list[Result]]:
    """Pair first, exclude zero totals, then average experiments once within each study."""
    by_key = {profile.key: profile for profile in profiles}
    if len(by_key) != len(profiles):
        raise ValueError("Duplicate reconstructed profile")
    pairs = {(s, e) for s, e, _ in by_key}
    if pairs - metadata.keys():
        raise ValueError("Profile has no experiment metadata")
    grouped: dict[tuple[str, str], list[tuple[np.ndarray, float, float]]] = defaultdict(list)
    ledger = []
    for pair, meta in sorted(metadata.items()):
        inc, dec = by_key.get((*pair, "inc")), by_key.get((*pair, "dec"))
        status = "included"
        if meta["host_species"] != HUMAN_HOST_SPECIES:
            status = "non_human"
        elif meta["sequencing_type"] not in CORPUS_SEQUENCING_TYPES:
            status = "unsupported_sequencing"
        elif inc is None or dec is None:
            status = "unpaired_or_absent_profiles"
        elif sum(inc.raw) == 0 or sum(dec.raw) == 0:
            status = "zero_retained_total"
        ledger.append({"study_id": pair[0], "experiment_id": pair[1], "status": status})
        if status != "included" or inc is None or dec is None:
            continue
        if (
            inc.depth_primary is None
            or dec.depth_primary is None
            or inc.depth_secondary is None
            or dec.depth_secondary is None
        ):
            raise ValueError("Positive count profile has undefined depth")
        primary = float(inc.depth_primary) - float(dec.depth_primary)
        secondary = float(inc.depth_secondary) - float(dec.depth_secondary)
        for representation in REPRESENTATIONS:
            left, left_total = profile_values(inc, representation)
            right, right_total = profile_values(dec, representation)
            difference = np.asarray(left, dtype=float) / left_total
            difference -= np.asarray(right, dtype=float) / right_total
            grouped[(pair[0], representation)].append((difference, primary, secondary))
    studies = []
    for (study, representation), observations in sorted(grouped.items()):
        studies.append(
            StudyEvidence(
                study,
                representation,
                len(observations),
                tuple(float(x) for x in np.mean([row[0] for row in observations], axis=0)),
                float(np.mean([row[1] for row in observations])),
                float(np.mean([row[2] for row in observations])),
            )
        )
    return studies, ledger


def profile_values(profile: ProfileEvidence, representation: str) -> tuple[tuple[int, ...], int]:
    """Separate annotation mass, breadth-normalized presence, and taxon prevalence."""
    if representation == "raw":
        return profile.raw, sum(profile.raw)
    if representation == "presence":
        return profile.presence, sum(profile.presence)
    if representation == "prevalence":
        denominator = profile.ora_taxa_count
        return profile.presence, profile.matched_taxa if denominator is None else denominator
    raise ValueError(f"Unknown representation: {representation}")


def native_profile_inputs(
    profiles: Sequence[ProfileEvidence],
    metadata: Mapping[Pair, Mapping[str, str]],
    columns: Sequence[str],
    ledger: Sequence[Result],
) -> dict[str, list[Profile]]:
    """Use original paired profiles for the arms' own tests, preserving study clustering."""
    eligible = {
        (row["study_id"], row["experiment_id"]) for row in ledger if row["status"] == "included"
    }
    result: dict[str, list[Profile]] = {name: [] for name in REPRESENTATIONS}
    for profile in profiles:
        pair = profile.key[:2]
        if pair not in eligible:
            continue
        meta = metadata[pair]
        for representation in REPRESENTATIONS:
            values, denominator = profile_values(profile, representation)
            result[representation].append(
                Profile(
                    *profile.key,
                    meta.get("condition", ""),
                    meta.get("body_site", ""),
                    dict(zip(columns, values, strict=True)),
                    denominator,
                )
            )
    return result


def focal_feature_map(columns: Sequence[str]) -> dict[str, str]:
    """Resolve each declared focal term exactly; absent features remain unavailable."""
    mapped = {}
    for slug, terms in TABLE1_TRAITS.items():
        candidates = [
            column for column in columns if column in terms or column_term_id(column) in terms
        ]
        if len(candidates) > 1:
            raise ValueError(f"Ambiguous focal feature: {slug}")
        if candidates:
            mapped[slug] = candidates[0]
    return mapped


def estimate_hc3(y: Sequence[float], x: Sequence[float] | None = None) -> Result:
    """OLS intercept at zero depth difference, with HC3 SE and explicit support guards."""
    values = np.asarray(y, dtype=float)
    depth = np.asarray(x, dtype=float) if x is not None else np.zeros(len(values))
    if values.ndim != 1 or depth.shape != values.shape:
        raise ValueError("Outcome and depth must be aligned one-dimensional study observations")
    if not np.isfinite(values).all() or not np.isfinite(depth).all():
        raise ValueError("Outcome and depth must be finite")
    result: Result = {
        "n_studies": len(values),
        "effect": None,
        "slope": None,
        "standard_error": None,
        "ci_low": None,
        "ci_high": None,
        "p_value": None,
        "degrees_of_freedom": None,
        "depth_min": float(min(depth)) if len(depth) else None,
        "depth_max": float(max(depth)) if len(depth) else None,
        "design_rank": None,
        "condition_number": None,
        "maximum_leverage": None,
        "status": NO_ELIGIBLE_STUDIES if not len(values) else "insufficient_study_support",
    }
    if len(values) < 3:
        return result
    varying = np.ptp(depth) > 1e-12
    if not varying and np.max(np.abs(depth)) > 1e-12:
        result["status"] = "constant_nonzero_depth"
        return result
    if min(depth) > 0 or max(depth) < 0:
        result["status"] = "zero_depth_outside_support"
        return result
    scale = float(np.max(np.abs(depth))) if varying else 1.0
    design = (
        np.column_stack([np.ones(len(values)), depth / scale])
        if varying
        else (np.ones((len(values), 1)))
    )
    rank = int(np.linalg.matrix_rank(design))
    result.update(design_rank=rank, condition_number=float(np.linalg.cond(design)))
    if rank != design.shape[1] or result["condition_number"] > 1e10:
        result["status"] = "rank_deficient_or_ill_conditioned"
        return result
    inv = np.linalg.inv(design.T @ design)
    beta = inv @ design.T @ values
    residual = values - design @ beta
    leverage = np.sum((design @ inv) * design, axis=1)
    result["maximum_leverage"] = float(max(leverage))
    if max(leverage) >= 1 - 1e-10:
        result["status"] = "unit_leverage"
        return result
    covariance = inv @ (design.T * (residual / (1 - leverage)) ** 2) @ design @ inv
    se = math.sqrt(max(0, float(covariance[0, 0])))
    effect, df = float(beta[0]), len(values) - rank
    result.update(
        effect=effect,
        slope=float(beta[1] / scale) if varying else None,
        standard_error=se,
        degrees_of_freedom=df,
    )
    if se < 1e-15:
        result["status"] = "all_zero_outcomes" if np.all(values == 0) else "zero_residual_variance"
        return result
    critical = float(stats.t.ppf(0.975, df))
    result.update(
        status="tested",
        p_value=float(2 * stats.t.sf(abs(effect / se), df)),
        ci_low=effect - critical * se,
        ci_high=effect + critical * se,
    )
    return result


def load_core_evidence(
    bundle: CoreBundle,
) -> tuple[list[str], list[ProfileEvidence], dict[Pair, Result]]:
    """Use the core's already selected/oriented taxa; never reclassify or reverse them twice."""
    catalog = list(rows(bundle.path("feature_catalog")))
    columns = [row["column_id"] for row in catalog if row["feature_class"] == "trait_annotation"]
    profiles = profile_count_reader(bundle.path("count_profiles"))
    all_columns = profiles.count_columns
    if set(all_columns) != {row["column_id"] for row in catalog}:
        raise ValueError("Core profile columns disagree with the feature catalog")
    # The profile producer controls physical column order; the catalog controls membership.
    retained = set(columns)
    columns = [column for column in all_columns if column in retained]
    indices = [all_columns.index(column) for column in columns]
    observed: dict[Key, tuple[int, ...]] = {}
    for row in profiles:
        key = (row.key[0], row.key[1], row.key[2])
        if key in observed:
            raise ValueError("Duplicate profile key")
        observed[key] = tuple(row.counts[i] for i in indices)
    matrix_reader = taxa_trait_count_reader(bundle.path("taxa_matrix"))
    if matrix_reader.count_columns != all_columns:
        raise ValueError("Core matrix columns disagree with profile columns")
    matrix: dict[str, tuple[int, ...]] = {}
    for row in matrix_reader:
        if row.key[0] in matrix:
            raise ValueError("Duplicate matrix taxon")
        matrix[row.key[0]] = tuple(row.counts[i] for i in indices)
    reported: dict[Key, list[str]] = {key: [] for key in observed}
    for occurrence in rows(bundle.path("reported_taxa")):
        direction = {"increased": "inc", "decreased": "dec"}.get(occurrence["direction"])
        if direction is None:
            raise ValueError("Unknown reported-taxon direction")
        key = (occurrence["study_id"], occurrence["experiment_id"], direction)
        # The reported-taxon product includes unsupported sequencing methods; the
        # count-profile contract explicitly restricts to 16S/WMS.
        if key in reported:
            reported[key].append(occurrence["taxon_id"])
    metadata: dict[Pair, Result] = {}
    for experiment in rows(bundle.path("experiments")):
        pair = experiment["study_id"], experiment["experiment_id"]
        if pair in metadata:
            raise ValueError("Duplicate experiment metadata")
        metadata[pair] = experiment
    evidence = reconstruct_profiles(observed, reported, matrix, columns)
    ora = profile_count_reader(bundle.path("ora_profiles"))
    if set(ora.count_columns) != {*all_columns, "taxa_count"}:
        raise ValueError("ORA profile vocabulary differs from the verified count vocabulary")
    ora_indices = [ora.count_columns.index(column) for column in columns]
    n_index = ora.count_columns.index("taxa_count")
    ora_rows: dict[Key, tuple[tuple[int, ...], int]] = {}
    for row in ora:
        key = (row.key[0], row.key[1], row.key[2])
        if key in ora_rows:
            raise ValueError("Duplicate ORA profile key")
        ora_rows[key] = tuple(row.counts[i] for i in ora_indices), row.counts[n_index]
    expected = {
        p.key for p in evidence if metadata[p.key[:2]]["host_species"] == HUMAN_HOST_SPECIES
    }
    if set(ora_rows) != expected:
        raise ValueError("Binary ORA profile population differs from reconstructed human profiles")
    verified = []
    for profile in evidence:
        if profile.key in ora_rows:
            presence, denominator = ora_rows[profile.key]
            if presence != profile.presence or denominator != profile.matched_taxa:
                raise ValueError(
                    "Binary ORA profile counts or taxa_count fail source reconstruction"
                )
            profile = replace(profile, presence=presence, ora_taxa_count=denominator)
        verified.append(profile)
    return columns, verified, metadata


def build_screen(
    studies: Sequence[StudyEvidence], columns: Sequence[str], native: Mapping[str, list[Profile]]
) -> list[Result]:
    """Reuse the native screen, retaining exact intervals and dependence diagnostics."""
    results = []
    for representation in REPRESENTATIONS:
        cohort = [row for row in studies if row.representation == representation]
        screened = run_screen(native[representation], columns)
        native_rows = {row.trait: row for row in screened.rows}
        for test in ("sign", "mean_hc3", "adjusted_primary", "adjusted_secondary"):
            family_id = f"{PROTOCOL}:{representation}:{test}"
            family_rows = []
            for index, column in enumerate(columns):
                y = [row.mean_differences[index] for row in cohort]
                positive = sum(v > SIGN_TOLERANCE for v in y)
                negative = sum(v < -SIGN_TOLERANCE for v in y)
                if test == "sign":
                    native_row = native_rows[column]
                    inference: Result = {
                        "effect": native_row.share_up,
                        "p_value": native_row.p_value,
                        "status": native_row.test_status,
                        "ci_low": native_row.ci_low,
                        "ci_high": native_row.ci_high,
                        "native_bh_q_nominal": native_row.bh_q_nominal,
                        "native_by_q_nominal": native_row.by_q_nominal,
                        "native_by_q_effective": native_row.by_q_effective,
                        "effective_family_size": screened.effective_family_size,
                        "effective_family_size_used": screened.effective_family_size_used,
                    }
                    if (positive, negative, len(y) - positive - negative) != (
                        native_row.studies_up,
                        native_row.studies_down,
                        native_row.studies_tied,
                    ):
                        raise ValueError("Native and reconstructed study sign tallies disagree")
                else:
                    depth = (
                        [row.depth_primary for row in cohort]
                        if test == "adjusted_primary"
                        else [row.depth_secondary for row in cohort]
                        if test == "adjusted_secondary"
                        else None
                    )
                    inference = estimate_hc3(y, depth)
                family_rows.append(
                    {
                        "family_id": family_id,
                        "representation": representation,
                        "test": test,
                        "feature_id": column,
                        "hypothesis_id": column,
                        "feature_ids": [column],
                        "n_studies": len(cohort),
                        "n_positive": positive,
                        "n_negative": negative,
                        "n_tied": len(y) - positive - negative,
                        "mean_difference": float(np.mean(y)) if y else None,
                        "effect_unit": "positive_study_share"
                        if test == "sign"
                        else "fraction_difference",
                        "null_effect": 0.5 if test == "sign" else 0,
                        **inference,
                    }
                )
            results.extend(family_rows)
    return results


def focal_screen(screen: Sequence[Result], columns: Sequence[str]) -> list[Result]:
    """Keep primary full-family columns beside explicitly named six-trait comparisons."""
    mapped = focal_feature_map(columns)
    families = {row["family_id"] for row in screen}
    results = []
    for family in sorted(families):
        source = {row["feature_id"]: row for row in screen if row["family_id"] == family}
        template = next(iter(source.values()))
        for slug, terms in TABLE1_TRAITS.items():
            column = mapped.get(slug)
            row = (
                dict(source[column])
                if column is not None
                else {
                    "representation": template["representation"],
                    "test": template["test"],
                    "feature_id": terms[0],
                    "feature_ids": [terms[0]],
                    "status": "missing_feature",
                    "p_value": None,
                }
            )
            row.update(
                family_id=f"{family}:focal_six",
                hypothesis_id=slug,
                focal_trait=slug,
                full_family_id=family,
                q_by=row.get("q_by"),
                q_bh=row.get("q_bh"),
                family_size=template["family_size"],
            )
            results.append(row)
    return results


def composition_screen(
    native: Mapping[str, list[Profile]], columns: Sequence[str]
) -> tuple[list[Result], list[Result]]:
    """Call Jed's native paired-ratio test and preserve all fifteen focal pair slots."""
    mapped = focal_feature_map(columns)
    summaries: list[Result] = []
    studies: list[Result] = []
    for representation in ("raw", "presence"):
        focal_profiles = [
            replace(
                p, guild_counts={slug: p.guild_counts[column] for slug, column in mapped.items()}
            )
            for p in native[representation]
        ]
        result = (
            compute_log_ratio_shifts(focal_profiles, tuple(mapped))
            if focal_profiles and len(mapped) >= 2
            else None
        )
        available = {(r.numerator, r.reference): r for r in result.summaries} if result else {}
        if result:
            studies.extend(
                {"representation": representation, **asdict(row)} for row in result.studies
            )
        for numerator, reference in combinations(TABLE1_TRAITS, 2):
            summary = available.get((numerator, reference))
            row: Result = (
                asdict(summary)
                if summary
                else {
                    "numerator": numerator,
                    "reference": reference,
                    "status": "missing_feature"
                    if numerator not in mapped or reference not in mapped
                    else NO_ELIGIBLE_STUDIES,
                    "sign_p": None,
                    "by_q": None,
                }
            )
            row.update(
                family_id=f"{PROTOCOL}:{representation}:composition",
                representation=representation,
                test="paired_log2_ratio_sign",
                hypothesis_id=f"{numerator}/{reference}",
                feature_ids=[
                    mapped.get(slug, TABLE1_TRAITS[slug][0]) for slug in (numerator, reference)
                ],
                p_value=row["sign_p"],
                native_available_pair_by_q=row.pop("by_q"),
                native_available_pair_family_size=result.family_size if result else 0,
            )
            summaries.append(row)
    return summaries, studies


def _write_tsv(path: Path, records: Iterable[Mapping[str, Any]], fields: Sequence[str]) -> None:
    """Write a stable TSV (optionally gzip) with explicit empty fields for undefined values."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        compressed = (
            gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0)
            if (path.suffix == ".gz")
            else None
        )
        with io.TextIOWrapper(compressed or raw, encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(
                {
                    key: json.dumps(value, separators=(",", ":"), allow_nan=False)
                    if isinstance(value, (list, tuple, dict))
                    else value
                    for key, value in record.items()
                }
                for record in records
            )


def implementation_sources() -> dict[str, str]:
    """Fingerprint loaded canonical mtp modules, plus the shared correction policy.

    Inventory the actual process, including modules loaded indirectly by core
    verification. A long-lived process can include additional mtp modules; this
    is a conservative inventory, not a claim that each module executed a line.
    """
    paths = {Path(__file__).resolve()}
    package_root = (REPO_ROOT / "src/mtp").resolve()
    for name, module in sys.modules.copy().items():
        if name != "mtp" and not name.startswith("mtp."):
            continue
        filename = getattr(module, "__file__", None)
        if filename is None:
            continue
        path = Path(filename).resolve()
        if not path.is_relative_to(package_root) or path.suffix != ".py":
            raise ValueError(f"Loaded mtp module is outside canonical Python sources: {name}")
        paths.add(path)
    paths.add(POLICY_PATH)
    return {path.relative_to(REPO_ROOT).as_posix(): file_sha256(path) for path in sorted(paths)}


def run_sensitivity(core_run: Path, output_dir: Path) -> Path:
    """Publish only separately named sidecars after complete input reconstruction."""
    initially_loaded = implementation_sources()
    # Guard even not-yet-imported canonical modules against changes before a lazy
    # import. Only modules actually loaded are included in published provenance.
    source_snapshot = {
        path.relative_to(REPO_ROOT).as_posix(): file_sha256(path)
        for path in (REPO_ROOT / "src/mtp").rglob("*.py")
    }
    source_snapshot.update(initially_loaded)
    bundle = open_core_bundle(core_run)
    ora_bundle = open_core_bundle(core_run, "ora_profiles")
    binding = output_dir / "core_consumption.json"
    ora_binding = output_dir / "ora_core_consumption.json"
    preflight_consumer_outputs(bundle, binding, output_directory=output_dir)
    preflight_consumer_outputs(ora_bundle, ora_binding, output_directory=output_dir)
    if output_dir.exists():
        raise FileExistsError("Sensitivity outputs require a new directory")
    rollup = run_strain_rollup(bundle.run)
    if rollup == "full":
        raise ValueError(
            "The full strain-rollup scope is not supported by the multiplicity contract; "
            "presence/depth sensitivity supports only off and phenotype"
        )
    columns, profiles, metadata = load_core_evidence(bundle)
    studies, ledger = make_study_observations(profiles, metadata)
    native = native_profile_inputs(profiles, metadata, columns, ledger)
    screen = build_screen(studies, columns, native)
    composition, composition_studies = composition_screen(native, columns)
    missing = sorted(set(TABLE1_TRAITS) - focal_feature_map(columns).keys())
    bases = [
        build_analysis_basis(
            basis_id=f"{PROTOCOL}:{representation}",
            feature_catalog_sha256=file_sha256(bundle.path("feature_catalog")),
            feature_filter_id=feature_filter_id(columns),
            retained_feature_count=len(columns),
            profile_unit="count_profiles" if representation == "raw" else "ora_profiles",
            aggregation_policy="distinct_reported_taxa; equal_experiment_means_then_equal_studies",
            rollup_policy=f"ncbi_taxonomy_closure_rollup_{rollup}",
            denominator_policy=(
                "sum_retained_annotation_counts"
                if representation == "raw"
                else "sum_retained_binary_taxon_trait_presences_breadth_weighted"
                if representation == "presence"
                else "ora_profiles_taxa_count_matched_reported_taxa"
            ),
            background_inventory_id="none_paired_profile_contrast",
            weighting_policy="equal_study",
        ).model_dump(mode="json")
        for representation in REPRESENTATIONS
    ]
    basis_map = dict(zip(REPRESENTATIONS, bases, strict=True))
    for representation in ("raw", "presence"):
        basis = build_analysis_basis(
            basis_id=f"{PROTOCOL}:{representation}:composition",
            feature_catalog_sha256=file_sha256(bundle.path("feature_catalog")),
            feature_filter_id=feature_filter_id(columns),
            retained_feature_count=len(columns),
            profile_unit="count_profiles" if representation == "raw" else "ora_profiles",
            aggregation_policy="paired_log2_ratio_difference_then_equal_study_means",
            rollup_policy=f"ncbi_taxonomy_closure_rollup_{rollup}",
            denominator_policy="positive_reference_trait_count_no_pseudocount",
            background_inventory_id="none_paired_profile_contrast",
            weighting_policy="equal_study",
        ).model_dump(mode="json")
        basis_map[f"{representation}:composition"] = basis
        bases.append(basis)
    families, corrections = declare_and_correct(screen, basis_map, bundle, ledger, rollup)
    focal = focal_screen(screen, columns)
    extra_families, extra_corrections = declare_and_correct(
        [*focal, *composition], basis_map, bundle, ledger, rollup
    )
    families.extend(extra_families)
    corrections.extend(extra_corrections)
    final_sources = implementation_sources()
    changed_sources = sorted(
        name
        for name in initially_loaded.keys() | final_sources.keys()
        if final_sources.get(name) != source_snapshot.get(name)
    )
    if changed_sources:
        raise ValueError(
            "Sensitivity implementation or policy changed during analysis: "
            + ", ".join(changed_sources)
        )
    output_dir.mkdir(parents=True)
    outputs: list[Path] = []

    def table(name: str, records: Iterable[Mapping[str, Any]], fields: Sequence[str]) -> None:
        path = output_dir / name
        _write_tsv(path, records, fields)
        outputs.append(path)

    def fields_for(records: Sequence[Result], fallback: Sequence[str]) -> list[str]:
        return list(dict.fromkeys(key for row in records for key in row)) or list(fallback)

    fields = fields_for(screen, ["family_id", "hypothesis_id", "status"])
    # DictWriter needs complete row keys only; missing inferential fields remain blank.
    table("trait_screen.tsv", screen, fields)
    table("focal_traits.tsv", focal, fields_for(focal, ["focal_trait", *fields]))
    table("composition_presence.tsv", composition, fields_for(composition, ["hypothesis_id"]))
    table(
        "composition_studies.tsv",
        composition_studies,
        fields_for(composition_studies, ["study_id"]),
    )
    table("eligibility.tsv", ledger, ["study_id", "experiment_id", "status"])
    audit = [
        {
            "study_id": p.key[0],
            "experiment_id": p.key[1],
            "inc_dec": p.key[2],
            "raw_total": sum(p.raw),
            "presence_total": sum(p.presence),
            **{
                key: value
                for key, value in asdict(p).items()
                if key not in {"key", "raw", "presence"}
            },
            "unmatched_taxa": ";".join(p.unmatched_taxa),
            "unmatched_count": len(p.unmatched_taxa),
        }
        for p in profiles
    ]
    table("profile_audit.tsv", audit, list(audit[0]) if audit else ["study_id"])
    table(
        "study_effects.tsv.gz",
        (
            {
                "study_id": row.study_id,
                "representation": row.representation,
                "experiments": row.experiments,
                "feature_id": column,
                "mean_difference": value,
                "depth_primary": row.depth_primary,
                "depth_secondary": row.depth_secondary,
            }
            for row in studies
            for column, value in zip(columns, row.mean_differences, strict=True)
        ),
        [
            "study_id",
            "representation",
            "experiments",
            "feature_id",
            "mean_difference",
            "depth_primary",
            "depth_secondary",
        ],
    )
    code_commit, code_dirty = repository_state(REPO_ROOT)
    provenance = {
        "protocol": PROTOCOL,
        "core_run_id": bundle.run.run_id,
        "core_manifest_sha256": file_sha256(bundle.manifest),
        "contrast_population": bundle.run.contrast_population,
        "direction_basis": bundle.run.direction_basis,
        "strain_rollup": rollup,
        "retained_features": len(columns),
        "feature_filter_id": feature_filter_id(columns),
        "missing_focal_traits": missing,
        "profiles": len(profiles),
        "included_studies": sum(row.representation == "raw" for row in studies),
        "eligibility_counts": {
            status: sum(row["status"] == status for row in ledger)
            for status in sorted({row["status"] for row in ledger})
        },
        "consumed_inputs": {
            name: file_sha256(bundle.path(name))
            for name in (
                "count_profiles",
                "ora_profiles",
                "reported_taxa",
                "taxa_matrix",
                "experiments",
                "feature_catalog",
            )
        },
        "implementation_sha256": final_sources["src/mtp/presence_depth_sensitivity.py"],
        "code_commit": code_commit,
        "code_dirty": code_dirty,
        "implementation_sources_sha256": final_sources,
        "methods_sha256": file_sha256(REPO_ROOT / "methods/presence_depth_sensitivity.md"),
        "dependency_lock_sha256": file_sha256(REPO_ROOT / "uv.lock"),
        "secondary_depth": "log1p(sum retained annotation counts / "
        "distinct matrix-matched reported taxa)",
        "primary_depth": "mean log10(1 + max retained annotation count per matched reported taxon)",
        "aggregation": "profile transform, inc-minus-dec difference, within-study mean",
        "sign_tolerance": SIGN_TOLERANCE,
        "interpretation": "conditional descriptive sensitivity; not a causal correction",
        "calibration": "HC3 with t(n-rank) reference; approximate, not empirically calibrated",
        "multiplicity": "BY primary and BH comparison separately per representation and test; "
        "exploratory across families",
    }
    for name, content in (
        ("analysis_bases.json", bases),
        ("families.json", families),
        ("corrections.json", corrections),
        ("provenance.json", provenance),
    ):
        path = output_dir / name
        path.write_text(json.dumps(content, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        outputs.append(path)
    readme = output_dir / "README.md"
    readme.write_text(
        "# Presence and annotation-depth sensitivity\n\n"
        f"Protocol `{PROTOCOL}`; {len(columns)} retained features; "
        f"{provenance['included_studies']} equally weighted studies.\n\n"
        "Pairing alone does not remove between-list annotation-depth differences. "
        "The adjusted intercept estimates the fraction difference at zero paired depth difference; "
        "it is a conditional association, not a causal correction. Unsupported fits have no p/q. "
        "Presence fractions retain annotation-breadth weighting; taxon prevalence uses the "
        "verified binary ORA profile k/taxa_count. Native sign-test intervals and effective-m "
        "diagnostics are retained; the composition result uses the native paired-ratio test. "
        "The focal extract retains full-screen q_by/q_bh and family_size, and reports "
        "q_by_six_trait/q_bh_six_trait with family_size_six_trait for the separate "
        "six-trait comparison. Families remain exploratory "
        "across representations and tests. See methods/presence_depth_sensitivity.md for "
        "the declared transforms, HC3 support guards and denominators.\n",
        encoding="utf-8",
    )
    outputs.append(readme)
    write_consumer_record(
        ora_bundle,
        ora_binding,
        consumer="mtp.presence_depth_sensitivity.binary",
        outputs=outputs,
        additional_products=("reported_taxa",),
    )
    verify_consumer_record(ora_binding)
    write_consumer_record(
        bundle,
        binding,
        consumer="mtp.presence_depth_sensitivity",
        outputs=[*outputs, ora_binding],
        additional_products=("reported_taxa",),
    )
    verify_consumer_record(binding)
    return binding


def main(argv: Sequence[str] | None = None) -> int:
    """Run the standalone sensitivity without altering any native arm plan or result."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    add_temporary_output_option(parser)
    args = parser.parse_args(argv)
    require_durable_output(parser, args.output_dir, allowed=args.allow_temporary_output)
    print(run_sensitivity(args.core_run, args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
