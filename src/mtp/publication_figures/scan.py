"""Verified candidate-triple distributions, preserving historical scan protocols (#1046).

The plotted population is every tested cellulose triple, not simulated null draws.
Ranks, p-values and bins come from the saved scan; checks do not rerun inference.
"""

from __future__ import annotations

import gzip
import json
from collections.abc import Iterator
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import pandas as pd

from mtp.cellulose_decr import selection_scan as native
from mtp.cellulose_decr import selection_scan_scopes as scopes
from mtp.datamodel.workflow_pydantic import ConsumerRecord, WorkflowRun
from mtp.publication_figures.inputs import _identity, _json, _read, _require, _table
from mtp.publication_figures.plots import (
    ACCENT,
    GRAY,
    INK,
    STYLE,
    Page,
    _figure,
    _wrap,
    validate_authored_wording,
)
from mtp.trait_screen import TESTED
from mtp.workflow.consumers import verify_consumer_record
from mtp.workflow.trait_policy import resolve_trait_policy

KEYS = ["trait_id", "condition_a", "condition_b", "condition_c"]
STATS = {
    native.MEAN_DISTANCE: ("u6_p_mean_distance", "rank_mean_distance", "Mean Jaccard distance"),
    native.CORE: ("u6_p_core", "rank_core", "Three-way core"),
}


@dataclass(frozen=True)
class ScanFigureInputs:
    """Original source records, complete tables and late publication checks."""

    record: dict[str, Any]
    snapshots: dict[Path, dict[str, Any]]
    targets: dict[Path, Path]
    tables: dict[str, pd.DataFrame]
    source_roots: list[Path]
    context: dict[str, Any]
    summary: dict[str, Any]
    labels: dict[str, str]
    native_labels: dict[str, str]
    core_manifest_sha256: str
    pairing_identity: dict[str, Any] | None
    native_run: scopes.ScanRun


def _numeric(frame: pd.DataFrame, column: str, *, integer: bool = False) -> pd.Series:
    _require(column in frame, f"Scan column missing: {column}")
    values = pd.to_numeric(frame[column], errors="raise")
    _require(bool(np.isfinite(values).all()), f"Nonfinite scan value: {column}")
    if integer:
        _require(bool(((values >= 0) & (values % 1 == 0)).all()), f"Invalid scan count: {column}")
    return values


def _booleans(frame: pd.DataFrame, column: str) -> pd.Series:
    _require(column in frame, f"Scan column missing: {column}")
    values = frame[column].astype(str)
    _require(bool(values.isin(["True", "False"]).all()), f"Invalid scan boolean: {column}")
    return values.eq("True")


def _count(value: Any, description: str) -> int:
    _require(isinstance(value, int) and not isinstance(value, bool) and value >= 0, description)
    return int(value)


def _labels(metadata: pd.DataFrame, *, display: bool) -> dict[str, str]:
    _require(
        {"column_id", "column_label"}.issubset(metadata), "Scan label metadata columns missing"
    )
    _require(not metadata["column_id"].duplicated().any(), "Duplicate scan metadata trait")
    result = {
        str(row["column_id"]): str(
            (row.get("object_name") or row["column_label"]) if display else row["column_label"]
        )
        for row in metadata.to_dict("records")
    }
    _require(
        all(result) and all(value.strip() for value in result.values()), "Blank verified scan label"
    )
    return result


def _validate(bundle: ScanFigureInputs, protocol: native.Protocol) -> None:
    summary, tables = bundle.summary, bundle.tables
    triples, traits, distribution = tables["triples"], tables["traits"], tables["distribution"]
    conditions = tables["conditions"]
    _require(set(KEYS).issubset(triples), "Scan candidate keys missing")
    _require(not triples[KEYS].duplicated().any(), "Duplicate scan candidate triple")
    _require(not traits["trait_id"].duplicated().any(), "Duplicate selected trait")
    _require(set(triples["trait_id"]) == set(traits["trait_id"]), "Scan trait inventories differ")
    _require(summary["lead_trait"] == protocol.lead_trait, "Scan lead trait differs from protocol")
    _require(
        summary["lead_triple"]["conditions"] == list(protocol.lead_triple),
        "Scan lead conditions differ from protocol",
    )
    _require(
        set(STATS)
        == {summary["rule"]["primary_statistic"], summary["rule"]["secondary_statistic"]},
        "Scan statistic declaration differs",
    )
    for key, value in (
        ("min_studies", protocol.min_studies),
        ("top_n", protocol.top_n),
        ("always_included", protocol.always_included),
        ("distance_tolerance", protocol.tolerance),
    ):
        _require(summary["rule"].get(key) == value, f"Scan selection rule differs: {key}")
    _require(
        bundle.record.get("follows_protocol") is True, "Scan does not follow its declared protocol"
    )
    for key, value in (("permutations", protocol.permutations), ("seed", protocol.seed)):
        _require(
            summary["rule"][key] == bundle.record[key] == value, f"Scan {key} differs from protocol"
        )
    _require(not conditions["condition"].duplicated().any(), "Duplicate scan condition")
    eligible = conditions.loc[_booleans(conditions, "eligible")]
    _require(
        bool((_numeric(eligible, "stool_studies", integer=True) >= protocol.min_studies).all()),
        "Ineligible scan condition",
    )
    expected_keys = set(combinations(sorted(eligible["condition"].astype(str)), 3))
    counts = summary["counts"]
    _require(
        _count(counts["eligible_conditions"], "Invalid condition count") == len(eligible),
        "Scan eligible condition count differs",
    )
    _require(
        _count(counts["traits"], "Invalid trait count") == len(traits),
        "Scan selected trait count differs",
    )
    _require(
        _count(counts["triples"], "Invalid triple count") == len(expected_keys),
        "Scan candidate denominator differs",
    )
    for trait in traits["trait_id"]:
        chosen = triples[triples["trait_id"].eq(trait)]
        actual = {
            tuple(str(row[column]) for column in KEYS[1:]) for row in chosen.to_dict("records")
        }
        _require(actual == expected_keys, "Incomplete or unexpected scan candidate population")
        label = traits.loc[traits["trait_id"].eq(trait), "trait_label"].item()
        _require(
            trait in bundle.native_labels and label == bundle.native_labels[trait],
            "Scan trait label differs from bound metadata",
        )
    mine = triples[triples["trait_id"].eq(protocol.lead_trait)].copy()
    _require(not mine.empty, "Declared cellulose trait absent from scan")
    for field in (
        "n_a",
        "n_b",
        "n_c",
        "u6_units_a",
        "u6_units_b",
        "u6_units_c",
        "shared_studies",
        "three_way_core",
        "three_way_union",
        "u6_draws",
        "u6_draws_without_distance",
    ):
        _numeric(mine, field, integer=True)
    _require(
        bool(mine["test_status"].isin([TESTED, native.NO_TAXA, native.POOL_TOO_SMALL]).all()),
        "Unknown scan test status",
    )
    lead_mask = _booleans(mine, "is_lead_triple")
    expected_lead = mine[KEYS[1:]].apply(lambda row: set(row) == set(protocol.lead_triple), axis=1)
    _require(
        lead_mask.equals(expected_lead), "Scan lead-triple flag differs from declared conditions"
    )
    _require(int(lead_mask.sum()) == 1, "Scan must contain exactly one declared lead triple")
    tested = mine[mine["test_status"].eq(TESTED)]
    _require(not tested.empty, "No tested cellulose triples")
    _require(
        bool((_numeric(tested, "u6_draws", integer=True) > 0).all()), "Tested scan has no draws"
    )
    _require(
        bool(
            (
                _numeric(tested, "u6_draws_without_distance", integer=True) <= tested["u6_draws"]
            ).all()
        ),
        "Scan missing-draw support exceeds draws",
    )
    lead = mine.loc[lead_mask].iloc[0]
    _require(lead["test_status"] == TESTED, "Declared lead triple was not tested")
    trait = traits.loc[traits["trait_id"].eq(protocol.lead_trait)].iloc[0]
    _require(
        str(trait["is_lead_trait"]) == "True" and str(trait["lead_in_scan"]) == "True",
        "Scan trait callout flag differs",
    )
    for name in (*KEYS[1:], "test_status"):
        _require(
            str(trait[f"lead_{name}"]) == str(lead[name]), f"Scan trait callout differs: {name}"
        )
    for name in (native.MEAN_DISTANCE, native.CORE):
        _require(
            float(trait[f"lead_{name}"]) == float(lead[name]),
            f"Scan trait callout differs: {name}",
        )
    own = distribution[distribution["trait_id"].eq(protocol.lead_trait)]
    _require(
        len(own) == len(STATS) * len(native.REFERENCE_SETS) * protocol.bins,
        "Scan distribution inventory differs",
    )
    _require(
        set(own["statistic"]) == set(STATS)
        and set(own["reference_set"]) == set(native.REFERENCE_SETS),
        "Scan distribution reference declarations differ",
    )
    for statistic, (p_column, rank_column, _) in STATS.items():
        p = _numeric(tested, p_column)
        _require(bool(((p > 0) & (p <= 1)).all()), "Invalid nominal scan p-value")
        ranks = _numeric(tested, rank_column, integer=True)
        _require(bool(((ranks >= 1) & (ranks <= len(tested))).all()), "Invalid saved scan rank")
        # Six-significant-figure export can merge distinct exact p-values. Retain
        # stored ranks, accepting only the rank interval consistent with those bytes.
        for value, rank in zip(p, ranks, strict=True):
            _require(
                int((p < value).sum()) + 1 <= rank <= int((p <= value).sum()),
                "Saved scan rank disagrees with candidate order",
            )
        _require(
            float(trait[f"lead_{p_column}"]) == float(lead[p_column]), "Scan callout p differs"
        )
        for reference in native.REFERENCE_SETS:
            population = (
                tested
                if reference == native.ALL_TRIPLES
                else tested[tested["shared_studies"].eq(0)]
            )
            bins = own[
                own["statistic"].eq(statistic) & own["reference_set"].eq(reference)
            ].sort_values("p_above")
            left, right = _numeric(bins, "p_above"), _numeric(bins, "p_at_most")
            edges = np.linspace(0, 1, protocol.bins + 1)
            _require(
                len(bins) == protocol.bins
                and np.allclose(left, edges[:-1], rtol=0, atol=1e-14)
                and np.allclose(right, edges[1:], rtol=0, atol=1e-14),
                "Scan distribution bin edges differ",
            )
            actual_counts = _numeric(bins, "triples", integer=True).astype(int).tolist()
            p_values = _numeric(population, p_column)
            expected_counts = [
                int(((p_values > lo) & (p_values <= hi)).sum())
                for lo, hi in zip(left, right, strict=True)
            ]
            _require(
                actual_counts == expected_counts and sum(actual_counts) == len(population),
                "Scan distribution bin counts differ from candidates",
            )
            extra = 0 if bool(population["is_lead_triple"].astype(str).eq("True").any()) else 1
            prefix = f"{statistic}_{reference}_"
            denominator = len(population) + extra
            rank_min, rank_max = int(trait[prefix + "rank_min"]), int(trait[prefix + "rank_max"])
            _require(
                int(trait[prefix + "triples"]) == denominator
                and 1 <= rank_min <= rank_max <= denominator,
                "Scan callout rank denominator differs",
            )
            low = int((p_values < float(lead[p_column])).sum()) + 1
            high = int((p_values <= float(lead[p_column])).sum()) + extra
            _require(
                low <= rank_min <= rank_max <= high,
                "Scan callout rank differs from candidate population",
            )
            if reference == native.ALL_TRIPLES:
                _require(rank_min == int(lead[rank_column]), "Scan callout and triple rank differ")
            # The producer writes six significant figures, including JSON.
            _require(
                format(float(trait[prefix + "share_at_most"]), ".6g")
                == format(rank_max / denominator, ".6g"),
                "Scan callout rank share differs",
            )


def load_scan_inputs(path: Path, *, manifest_sha256: str) -> ScanFigureInputs:
    """Verify a pinned native ConsumerRecord, including historical v1, without inference."""
    path = path.expanduser().absolute()
    manifest = path / native.CONSUMER_RECORD if path.is_dir() else path
    _require(manifest.name == native.CONSUMER_RECORD, "Expected scan core_consumption.json")
    snapshots: dict[Path, dict[str, Any]] = {}
    targets: dict[Path, Path] = {}
    raw = _read(manifest, snapshots, targets, digest=manifest_sha256)
    record = ConsumerRecord.model_validate_json(raw)
    _require(record.consumer == native.CONSUMER, "Consumer is not a selection scan")
    verified = verify_consumer_record(manifest)
    _require(record == verified, "Scan consumer changed during verification")
    _read(manifest, snapshots, targets, digest=manifest_sha256)
    root = manifest.parent
    expected_names = {*native.OUTPUTS.values(), "core_consumption.coverage.tsv"}
    _require(len(record.outputs) == len(expected_names), "Scan output inventory differs")
    bound: dict[str, Any] = {}
    for output in record.outputs:
        output_path = Path(output.path)
        _require(
            output_path.name not in bound
            and output_path.resolve() == (root / output_path.name).resolve(),
            "Ambiguous or external scan output",
        )
        bound[output_path.name] = output
    _require(set(bound) == expected_names, "Scan output inventory differs")
    buffers = {}
    for name, item in bound.items():
        expected = {"sha256": item.sha256, "size_bytes": item.size_bytes}
        buffers[name] = _read(root / name, snapshots, targets, expected=expected)
        _read(Path(item.path), snapshots, targets, expected=expected)
    summary = _json(buffers[native.OUTPUTS["summary"]])
    provenance = _json(buffers[native.OUTPUTS["provenance"]])
    protocol = native.recorded_protocol(summary, provenance)
    protocol_pins = {protocol.version: provenance["protocol"]["sha256"]}
    if protocol.version != native.VERSION_V1:
        protocol_pins[native.VERSION_V1] = provenance["protocol"]["amends"]["sha256"]
    for version, pin in protocol_pins.items():
        _read(native.PROTOCOLS[version], snapshots, targets, digest=pin)
    core_path = Path(record.run_manifest_path)
    core_raw = _read(core_path, snapshots, targets, digest=record.run_manifest_sha256)
    core = WorkflowRun.model_validate_json(core_raw)
    _require(core.run_id == record.run_id, "Scan core execution differs")
    policy = resolve_trait_policy(core).policy
    scopes.check_bound_core(root, summary, provenance, core, record.run_manifest_sha256, protocol)
    _require(
        str(core.contrast_population) == native.TIER1_2
        and str(core.direction_basis) == "disease_vs_comparison",
        "Scan core population or direction differs",
    )
    for facts in (summary["core_run"], provenance["core_run"]):
        _require(
            facts["contrast_population"] == str(core.contrast_population)
            and facts["direction_basis"] == str(core.direction_basis),
            "Scan population or direction record differs",
        )
    _require(
        provenance["core_run"]["git_commit"] == core.git_commit
        and provenance["core_run"]["git_dirty"] == core.git_dirty,
        "Scan core producer identity differs",
    )
    _require(
        set(provenance["inputs"]) == set(native.CORE_PRODUCTS),
        "Scan source product inventory differs",
    )
    artifacts = core.artifacts or []
    products = dict(zip(record.product_ids, record.artifact_ids, strict=True))
    source_buffers = {}
    roots = [root, root.resolve(), core_path.parent, core_path.resolve().parent]
    for product in native.CORE_PRODUCTS:
        # The native scan binds reported_taxa in provenance and the complete
        # core manifest, but does not add it to the default ConsumerRecord list.
        # Resolve all scientific inputs through that already verified core; an
        # explicit consumer reference must still identify the same artifact.
        matches = [item for item in artifacts if item.product_id == product]
        _require(len(matches) == 1, f"Scan core product must be unique: {product}")
        artifact = matches[0]
        _require(
            product not in products or products[product] == artifact.id,
            f"Scan explicit product reference differs from core: {product}",
        )
        item = provenance["inputs"][product]
        source = core_path.parent / item["path"]
        artifact_path = core_path.parent / artifact.path
        _require(
            artifact.product_id == product and source.resolve() == artifact_path.resolve(),
            "Scan input path differs from bound product",
        )
        expected = {"sha256": artifact.sha256, "size_bytes": artifact.size_bytes}
        _require(
            item["sha256"] == artifact.sha256 and item["bytes"] == artifact.size_bytes,
            "Scan input identity differs from core artifact",
        )
        source_buffers[product] = _read(source, snapshots, targets, expected=expected)
        _read(artifact_path, snapshots, targets, expected=expected)
        roots.extend([source.parent, artifact_path.resolve().parent])
    tables = {
        name: _table(gzip.decompress(buffers[filename]) if name == "triples" else buffers[filename])
        for name, filename in native.OUTPUTS.items()
        if name not in {"summary", "provenance"}
    }
    rollup = summary["strain_rollup"]
    beside = (
        [
            scope
            for scope in (protocol.scope_primary, *protocol.scope_sensitivity)
            if scope != rollup["scope"]
        ]
        if protocol.side_by_side
        else []
    )
    if protocol.side_by_side:
        _require(
            rollup.get("report_beside") == beside and rollup.get("reported") == native.SIDE_BY_SIDE,
            "Scan companion declaration differs",
        )
    context = {
        "core_run_id": core.run_id,
        "strain_rollup": rollup["scope"],
        "contrast_population": str(core.contrast_population),
        "direction_basis": str(core.direction_basis),
        "source_commit": core.git_commit,
        "consumer_commit": provenance.get("code_commit"),
        "protocol_version": protocol.version,
        "scope_role": rollup["role"],
        "report_beside": beside,
        "historical": protocol.version == native.VERSION_V1,
        "trait_key_policy": policy,
    }
    native_run = scopes.ScanRun(
        root,
        summary,
        provenance,
        tables["traits"],
        protocol,
        {
            native.OUTPUTS[name]: _identity(buffers[native.OUTPUTS[name]])["sha256"]
            for name in scopes.READ
        },
        core.run_id,
    )
    pairing_identity = scopes.comparable(native_run) if protocol.version == native.VERSION else None
    metadata = _table(source_buffers["feature_metadata"])
    bundle = ScanFigureInputs(
        provenance,
        snapshots,
        targets,
        tables,
        list(dict.fromkeys(roots)),
        context,
        summary,
        _labels(metadata, display=True),
        _labels(metadata, display=False),
        record.run_manifest_sha256,
        pairing_identity,
        native_run,
    )
    _validate(bundle, protocol)
    recheck_scan_inputs(bundle)
    return bundle


def recheck_scan_inputs(bundle: ScanFigureInputs) -> None:
    """Check every consumed alias and byte identity again immediately before publication."""
    for path, target in bundle.targets.items():
        _require(path.resolve(strict=True) == target, "Scan input alias target changed")
    for path, expected in bundle.snapshots.items():
        _require(
            path.is_file() and _identity(path.read_bytes()) == expected,
            f"Scan input hash changed: {path.name}",
        )


def scan_audit(bundle: ScanFigureInputs) -> pd.DataFrame:
    """Every candidate/statistic, including untestable rows without plotted marks."""
    trait = bundle.summary["lead_trait"]
    chosen = bundle.tables["triples"]
    chosen = chosen[chosen["trait_id"].eq(trait)]
    rows = []
    for statistic, (p_column, rank_column, _) in STATS.items():
        for original in chosen.to_dict("records"):
            row = {str(key): value for key, value in original.items()}
            plotted = row["test_status"] == TESTED
            row.update(
                statistic=statistic,
                p_value=float(row[p_column]) if plotted else None,
                rank_min=int(row[rank_column]) if plotted else None,
                plotted=plotted,
                trait_label=bundle.labels[trait],
                reference_set=native.ALL_TRIPLES,
                strain_rollup=bundle.context["strain_rollup"],
            )
            row["row_id"] = json_row_id(statistic, row)
            rows.append(row)
    return pd.DataFrame(rows)


def json_row_id(statistic: str, row: dict[str, Any]) -> str:
    # JSON avoids collisions in condition labels containing slashes or separators.
    return json.dumps(
        [str(row["strain_rollup"]), statistic, *(str(row[key]) for key in KEYS)], ensure_ascii=False
    )


def scan_pages(bundle: ScanFigureInputs) -> Iterator[Page]:
    """Two saved p distributions with every candidate and the declared lead callout."""
    audit = scan_audit(bundle)
    plotted = audit.loc[audit["plotted"]]
    trait_id = bundle.summary["lead_trait"]
    trait = bundle.tables["traits"].loc[bundle.tables["traits"]["trait_id"].eq(trait_id)].iloc[0]
    candidate_count = int(bundle.summary["counts"]["triples"])
    tested = len(plotted) // len(STATS)
    conditions = bundle.summary["lead_triple"]["conditions"]
    abbreviations = {
        "Parkinson's disease": "PD",
        "HIV infection": "HIV",
        "Colorectal cancer": "CRC",
    }
    callout_conditions = " / ".join(abbreviations.get(name, name) for name in conditions)
    details: dict[str, Any] = {
        "candidate_triples": candidate_count,
        "tested_triples": tested,
        "untestable_triples": candidate_count - tested,
        "reference_set": native.ALL_TRIPLES,
        "protocol_version": bundle.context["protocol_version"],
        "scope_role": bundle.context["scope_role"],
        "trait_key_policy": bundle.context["trait_key_policy"],
        "report_beside": bundle.context["report_beside"],
        "inference": "saved nominal U6 p-values; no new inference",
        "panels": [],
    }
    with matplotlib.rc_context(STYLE):
        figure = _figure(7.0)
        figure.text(
            0.12,
            0.96,
            validate_authored_wording("Candidate-triple distribution"),
            fontsize=11,
            color=INK,
        )
        figure.text(0.12, 0.915, _wrap(bundle.labels[trait_id], 95), color=INK)
        figure.text(
            0.12,
            0.865,
            validate_authored_wording(
                f"{tested} tested of {candidate_count} triples; "
                f"{candidate_count - tested} untestable. Rollup: "
                f"{bundle.context['strain_rollup']} "
                f"({bundle.context['scope_role'].replace('_', ' ')})."
            ),
            fontsize=8,
        )
        figure.text(
            0.12,
            0.835,
            validate_authored_wording(f"trait keys: {bundle.context['trait_key_policy']}"),
            fontsize=8,
        )
        for index, (statistic, (_, _, title)) in enumerate(STATS.items()):
            axis = figure.add_axes((0.12, 0.52 - index * 0.31, 0.82, 0.21))
            rows = plotted[plotted["statistic"].eq(statistic)]
            bins = bundle.tables["distribution"]
            bins = bins[
                bins["trait_id"].eq(trait_id)
                & bins["statistic"].eq(statistic)
                & bins["reference_set"].eq(native.ALL_TRIPLES)
            ].sort_values("p_above")
            edges = [*bins["p_above"].astype(float), float(bins["p_at_most"].iloc[-1])]
            axis.stairs(bins["triples"].astype(int), edges, color=GRAY, linewidth=0.8)
            # Each tested row is present as a rug mark. Coincident p-values remain
            # coincident; no jitter, sampling or invented separation is introduced.
            axis.plot(
                rows["p_value"].astype(float),
                np.zeros(len(rows)),
                "|",
                color=GRAY,
                markersize=4,
                markeredgewidth=0.5,
                linestyle="none",
            )
            lead = rows.loc[rows["is_lead_triple"].astype(str).eq("True")].iloc[0]
            p_value = float(lead["p_value"])
            axis.axvline(p_value, color=ACCENT, linewidth=0.9)
            prefix = f"{statistic}_{native.ALL_TRIPLES}_"
            low, high = int(trait[prefix + "rank_min"]), int(trait[prefix + "rank_max"])
            rank = str(low) if low == high else f"{low}\u2013{high}"
            callout = f"{callout_conditions}: p = {p_value:.4g}; Rank {rank} of {tested}"
            axis.set_title(
                validate_authored_wording(title + "\n" + callout), loc="left", fontsize=8, pad=9
            )
            axis.set_xlabel(validate_authored_wording("Saved U6 nominal p-value"))
            axis.set_ylabel(validate_authored_wording("Candidate triples"))
            axis.set_xlim(0, 1)
            axis.set_ylim(bottom=0)
            axis.spines[["top", "right"]].set_visible(False)
            details["panels"].append(
                {
                    "statistic": statistic,
                    "plotted_candidates": len(rows),
                    "bin_counts": bins["triples"].astype(int).tolist(),
                    "lead_p": p_value,
                    "lead_rank_min": low,
                    "lead_rank_max": high,
                    "lead_conditions": list(bundle.summary["lead_triple"]["conditions"]),
                }
            )
        note = (
            "Stool decreased-signature annotation sets; every p-value is nominal. "
            "rank 1 = smallest nominal p-value. "
            "U6: whole-study-signature null; draws distinct studies from the stool-signature pool, "
            "takes one random condition unit per study, and splits them into three groups "
            "of the recorded sizes. "
            "Candidate triples share studies and taxa. "
        )
        note += "; ".join(
            f"{abbreviations[name]}: {name}" for name in conditions if name in abbreviations
        )
        if any(name in abbreviations for name in conditions):
            note += "."
        if bundle.context["historical"]:
            note += " Historical protocol v1; recorded scope roles retained."
        else:
            note += " Protocol v2 requires both rollup scopes together."
        figure.text(0.12, 0.025, _wrap(validate_authored_wording(note), 110), fontsize=7)
        page = Page(figure, plotted["row_id"].astype(str).tolist(), len(plotted), details)
    yield page
