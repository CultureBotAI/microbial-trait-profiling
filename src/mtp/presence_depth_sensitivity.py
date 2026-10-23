"""Additive presence and annotation-depth sensitivities of verified core profiles (#772)."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

from mtp.contribution_inputs import rows
from mtp.count_input import profile_count_reader, taxa_trait_count_reader
from mtp.guild_maps import SIGN_TOLERANCE
from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES, HUMAN_HOST_SPECIES
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import TABLE1_TRAITS, column_term_id
from mtp.trait_screen import step_up
from mtp.workflow.consumers import (
    CoreBundle,
    open_core_bundle,
    preflight_consumer_outputs,
    verify_consumer_record,
    write_consumer_record,
)
from mtp.workflow.coverage import feature_filter_id
from mtp.workflow.downstream_contracts import build_analysis_basis

Key = tuple[str, str, str]
Pair = tuple[str, str]
Result = dict[str, Any]
PROTOCOL = "presence-depth-study-v1"


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
                math.log1p(mean_total) if mean_total is not None else None,
                sum(math.log10(1 + max(vector)) for vector in matched) / len(matched)
                if matched
                else None,
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
        for representation in ("raw", "presence"):
            left, right = getattr(inc, representation), getattr(dec, representation)
            difference = np.asarray(left, dtype=float) / sum(left)
            difference -= np.asarray(right, dtype=float) / sum(right)
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
        "status": "insufficient_studies",
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
    if list(all_columns) != [row["column_id"] for row in catalog]:
        raise ValueError("Core profile columns disagree with the ordered feature catalog")
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
    return columns, reconstruct_profiles(observed, reported, matrix, columns), metadata


def build_screen(
    studies: Sequence[StudyEvidence], columns: Sequence[str]
) -> tuple[list[Result], list[Result]]:
    """Keep every retained hypothesis; focal rows inherit full-family BY/BH values."""
    results, families = [], []
    for representation in ("raw", "presence"):
        cohort = [row for row in studies if row.representation == representation]
        for test in ("sign", "mean_hc3", "adjusted_primary", "adjusted_secondary"):
            family_id = f"{PROTOCOL}:{representation}:{test}"
            family_rows = []
            for index, column in enumerate(columns):
                y = [row.mean_differences[index] for row in cohort]
                positive = sum(v > SIGN_TOLERANCE for v in y)
                negative = sum(v < -SIGN_TOLERANCE for v in y)
                if test == "sign":
                    n = positive + negative
                    inference: Result = {
                        "effect": positive / n if n else None,
                        "p_value": float(stats.binomtest(positive, n).pvalue) if n else None,
                        "status": "tested" if n else "all_ties" if y else "insufficient_studies",
                    }
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
            held = [row["p_value"] if row["p_value"] is not None else 1.0 for row in family_rows]
            for method in ("by", "bh"):
                adjusted = step_up(held, method=method, family_size=len(columns))
                for row, q in zip(family_rows, adjusted, strict=True):
                    row[f"q_{method}"] = q if row["p_value"] is not None else None
            families.append(
                {
                    "id": family_id,
                    "analysis_basis_id": f"{PROTOCOL}:{representation}",
                    "family_size": len(columns),
                    "feature_ids": list(columns),
                    "tested": sum(p is not None for p in [r["p_value"] for r in family_rows]),
                    "untestable_policy": "p=1 for adjustment; no reported p/q",
                    "primary": "BY",
                    "comparison": "BH",
                }
            )
            results.extend(family_rows)
    return results, families


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
            writer.writerows(records)


def run_sensitivity(core_run: Path, output_dir: Path) -> Path:
    """Publish only separately named sidecars after complete input reconstruction."""
    bundle = open_core_bundle(core_run)
    binding = output_dir / "core_consumption.json"
    preflight_consumer_outputs(bundle, binding, output_directory=output_dir)
    if output_dir.exists():
        raise FileExistsError("Sensitivity outputs require a new directory")
    columns, profiles, metadata = load_core_evidence(bundle)
    studies, ledger = make_study_observations(profiles, metadata)
    screen, families = build_screen(studies, columns)
    focal: list[Result] = []
    missing = []
    for slug, terms in TABLE1_TRAITS.items():
        candidates = [column for column in columns if column_term_id(column) in terms]
        if len(candidates) > 1:
            raise ValueError(f"Ambiguous focal feature: {slug}")
        if not candidates:
            missing.append(slug)
        else:
            focal.extend(
                {"focal_trait": slug, **row} for row in screen if row["feature_id"] == candidates[0]
            )
    rollup = next(
        (p.value for p in bundle.run.parameters or [] if p.name == "strain_rollup"),
        "off (pre-scope core)",
    )
    bases = [
        build_analysis_basis(
            basis_id=f"{PROTOCOL}:{representation}",
            feature_catalog_sha256=file_sha256(bundle.path("feature_catalog")),
            feature_filter_id=feature_filter_id(columns),
            retained_feature_count=len(columns),
            profile_unit="count_profiles" if representation == "raw" else "ora_profiles",
            aggregation_policy="distinct_reported_taxa; equal_experiment_means_then_equal_studies",
            rollup_policy=rollup,
            denominator_policy=(
                "sum_retained_annotation_counts"
                if representation == "raw"
                else "sum_retained_binary_taxon_trait_presences"
            ),
            background_inventory_id="none_paired_profile_contrast",
            weighting_policy="equal_study",
        ).model_dump(mode="json")
        for representation in ("raw", "presence")
    ]
    output_dir.mkdir(parents=True)
    outputs: list[Path] = []

    def table(name: str, records: Iterable[Mapping[str, Any]], fields: Sequence[str]) -> None:
        path = output_dir / name
        _write_tsv(path, records, fields)
        outputs.append(path)

    fields = list(dict.fromkeys(key for row in screen for key in row))
    # DictWriter needs complete row keys only; missing inferential fields remain blank.
    table("trait_screen.tsv", screen, fields)
    table("focal_traits.tsv", focal, ["focal_trait", *fields])
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
                "reported_taxa",
                "taxa_matrix",
                "experiments",
                "feature_catalog",
            )
        },
        "implementation_sha256": file_sha256(Path(__file__)),
        "methods_sha256": file_sha256(REPO_ROOT / "methods/presence_depth_sensitivity.md"),
        "dependency_lock_sha256": file_sha256(REPO_ROOT / "uv.lock"),
        "primary_depth": "log1p(sum retained annotation counts / "
        "distinct matrix-matched reported taxa)",
        "secondary_depth": "mean log10(1 + max retained annotation count "
        "per matched reported taxon)",
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
        "The focal extract retains full-screen BY/BH corrections. Families remain exploratory "
        "across representations and tests. See methods/presence_depth_sensitivity.md for "
        "the declared transforms, HC3 support guards and denominators.\n",
        encoding="utf-8",
    )
    outputs.append(readme)
    write_consumer_record(
        bundle,
        binding,
        consumer="mtp.presence_depth_sensitivity",
        outputs=outputs,
        additional_products=("reported_taxa",),
    )
    verify_consumer_record(binding)
    return binding


def main(argv: Sequence[str] | None = None) -> int:
    """Run the standalone sensitivity without altering any native arm plan or result."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(run_sensitivity(args.core_run, args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
