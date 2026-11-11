"""Verified comparison-bound execution and full replay of study splits (#854)."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
from collections.abc import Sequence
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import scipy

from mtp.cellulose_decr import selection_scan as scanner
from mtp.cellulose_decr.study_split import (
    VERSION,
    Protocol,
    StudySplit,
    _declaration,
    _scope_role,
    run_splits,
    summarize_stability,
)
from mtp.datamodel.multiple_testing_pydantic import SourceArtifact
from mtp.multiplicity.stage import code_sha256
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import load_profiles
from mtp.trait_screen import candidate_slots
from mtp.workflow.consumers import open_core_bundle, preflight_consumer_outputs
from mtp.workflow.output_roots import add_temporary_output_option, require_durable_output
from mtp.workflow.population import ORIENTED_DIRECTION_BASIS
from mtp.workflow.runtime import repository_state

if TYPE_CHECKING:
    from mtp.multiplicity.sources import VerifiedSources

PROTOCOL_FILE = REPO_ROOT / "methods/study_split_selection_protocol.json"
METHODS_FILE = REPO_ROOT / "methods/study_split_selection.md"
PRODUCTS = (*scanner.CORE_PRODUCTS, "feature_catalog")
IMPLEMENTATIONS = (
    "src/mtp/cellulose_decr/study_split.py",
    "src/mtp/cellulose_decr/study_split_run.py",
    "src/mtp/cellulose_decr/selection_scan.py",
    "src/mtp/cellulose_decr/analysis.py",
    "src/mtp/trait_screen.py",
    "src/mtp/trait_guilds.py",
    "src/mtp/multiple_testing.py",
    "src/mtp/multiplicity/sources.py",
    "src/mtp/datamodel/multiple_testing_pydantic.py",
    "methods/multiple_testing_policy_v1.yaml",
    "methods/cellulose_decr_selection_scan_protocol.json",
    "methods/cellulose_decr_selection_scan_protocol_v2.json",
    "methods/study_split_selection_protocol.json",
    "methods/study_split_selection_protocol_v1.json",
    "methods/study_split_selection.md",
    "uv.lock",
)
STABILITY_FIELDS = (
    "hypothesis_id",
    "trait_id",
    "conditions",
    "statistic",
    "total_splits",
    "selected_splits",
    "testable_splits",
    "by_rejections_005",
    "bh_rejections_005",
    "interpretation",
)


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n"


def _inputs(comparison_path: Path) -> tuple[VerifiedSources, dict[str, Path]]:
    # #853 supplies the single authoritative AnalysisBasis binder. Do not infer
    # a replacement basis from whatever filenames happen to be present.
    from mtp.multiplicity.sources import bind_recorded_hashes, open_verified_sources

    manifest = (
        comparison_path / "comparison_run.json" if comparison_path.is_dir() else comparison_path
    )
    comparison_digest = file_sha256(manifest)
    verified = open_verified_sources(comparison_path)
    if (
        verified.core.contrast_population != scanner.TIER1_2
        or verified.core.direction_basis != ORIENTED_DIRECTION_BASIS
    ):
        raise ValueError("Study splits require the verified, disease-oriented Tier 1+2 population")
    if verified.scope not in {"off", "phenotype"}:
        raise ValueError(
            "Study splits declare phenotype primary and off sensitivity only "
            "(amended from phenotype sensitivity)"
        )
    verified = bind_recorded_hashes(verified, comparison_digest)
    paths = {f"core:{name}": verified.product(name) for name in PRODUCTS}
    paths["comparison:manifest"] = Path(verified.comparison.manifest_path).resolve()
    paths["core:manifest"] = Path(verified.comparison.core_manifest_path)
    for arm in ("jed", "noel", "ilyas"):
        paths[f"arm:{arm}:manifest"] = Path(getattr(verified.comparison, f"{arm}_manifest_path"))
    return verified, paths


def _execution_identity(protocol: Protocol) -> dict[str, Any]:
    declared = _declaration(PROTOCOL_FILE)
    defaults = Protocol()
    numeric = {
        "splits": defaults.splits,
        "seed": defaults.seed,
        "min_studies_per_condition_per_partition": defaults.min_studies,
        "top_decreasing_traits": defaults.top_traits,
        "permutations": defaults.permutations,
        "tolerance": defaults.tolerance,
    }
    if any(declared[key] != value for key, value in numeric.items()):
        raise ValueError("Committed protocol disagrees with the implemented defaults")
    if declared["direction_basis"] != ORIENTED_DIRECTION_BASIS:
        raise ValueError("Committed protocol direction basis mismatch")
    return {
        "protocol": declared,
        "parameters": asdict(protocol),
        "follows_protocol": protocol == defaults,
        # The shared collector includes transitive helpers while rejecting empty
        # inventories and unresolved index stages, without following editor locks.
        "code_sha256": code_sha256(),
        "implementation": {name: file_sha256(REPO_ROOT / name) for name in IMPLEMENTATIONS},
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }


def _recorded_inputs(paths: dict[str, Path], verified: VerifiedSources | None) -> dict[str, str]:
    """Require upstream identities for publishing, rather than adopting disk edits."""
    digests = {}
    for name, path in sorted(paths.items()):
        expected = verified.expected_digest(path) if verified is not None else None
        if expected is None:
            raise ValueError(f"Source lacks an upstream recorded hash: {name}")
        if file_sha256(path) != expected:
            raise ValueError(f"Source bytes differ from upstream recorded hash: {name}")
        digests[name] = expected
    return digests


def _provenance(
    paths: dict[str, Path], protocol: Protocol, verified: VerifiedSources | None = None
) -> dict[str, Any]:
    return {**_execution_identity(protocol), "inputs": _recorded_inputs(paths, verified)}


def _scope_binding(verified: VerifiedSources) -> dict[str, Any]:
    """Capture the verified common corpus and the scope-specific core identity."""
    core = verified.core
    declared = _declaration(PROTOCOL_FILE)
    parameters = {item.name: item.value for item in core.parameters or []}
    if len(parameters) != len(core.parameters or []):
        raise ValueError("Duplicate core parameters in study-split scope binding")
    recorded_scope = parameters.pop("strain_rollup", declared["absent_means"])
    if recorded_scope != verified.scope:
        raise ValueError("Study-split scope contradicts its verified core parameters")
    sources = {
        item.product_id: item.sha256
        for item in core.artifacts or []
        if item.product_id in scanner.CORE_SOURCES
    }
    if set(sources) != set(scanner.CORE_SOURCES) or any(not value for value in sources.values()):
        raise ValueError("Study-split scope comparison needs every recorded raw/core source")
    return {
        "core_run_id": core.run_id,
        "strain_rollup": verified.scope,
        "role": _scope_role(verified.scope),
        "comparison_contract": {
            "population": core.contrast_population,
            "direction_basis": core.direction_basis,
            "trait_key_policy": core.trait_key_policy,
            "producer": {
                "definition_sha256": core.definition_sha256,
                "code_sha256": core.code_sha256,
                "lock_sha256": core.lock_sha256,
                "runtime_sha256": hashlib.sha256(core.runtime_json.encode("utf-8")).hexdigest(),
                "git_commit": core.git_commit,
                "git_dirty": core.git_dirty,
            },
            "parameters": parameters,
            "sources": sources,
        },
    }


def _stability_table(records: Sequence[dict[str, Any]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output, fieldnames=STABILITY_FIELDS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in summarize_stability(records):
        writer.writerow({**row, "conditions": json.dumps(row["conditions"], ensure_ascii=False)})
    return output.getvalue()


def _documents(
    verified: VerifiedSources,
    paths: dict[str, Path],
    protocol: Protocol,
    recorded_at: datetime,
) -> dict[str, str]:
    features = verified.features["count_profiles"]
    matrix = scanner.read_matrix(paths["core:taxa_matrix"], features)
    profiles, taxa = scanner.load_population(
        paths["core:experiments"],
        paths["core:reported_taxa"],
        paths["core:count_profiles"],
        matrix.index,
    )
    counts, _ = load_profiles(
        paths["core:count_profiles"],
        paths["core:experiments"],
        candidate_slots(features),
        reject_missing_counts=True,
    )
    sources = [
        SourceArtifact(identity=name, sha256=digest)
        for name, digest in _recorded_inputs(paths, verified).items()
    ]
    sources.append(
        SourceArtifact(identity=f"protocol:{VERSION}", sha256=file_sha256(PROTOCOL_FILE))
    )
    records = run_splits(
        profiles,
        taxa,
        matrix,
        counts,
        features,
        verified.bases["count_profiles"],
        sources,
        protocol,
        analysis_id=verified.comparison.id,
        strain_rollup=verified.scope,
        recorded_at=recorded_at,
    )
    selected = sum(bool(record["selected_candidates"]) for record in records)
    role = _scope_role(verified.scope)
    companion = "off" if verified.scope == "phenotype" else "phenotype"
    readme = (
        "# Internal study-split selection assessment\n\n"
        f"Strain rollup `{verified.scope}`: "
        f"{role}. "
        f"Report beside the `{companion}` scope with the side-by-side report; "
        "do not choose the favorable result. This role was amended after results.\n\n"
        f"Completed {len(records)} splits; {selected} selected at least one candidate. "
        "Empty and untestable splits remain recorded.\n\n"
        "Each split selects only on its training studies and tests only on its held-out studies. "
        "Both overlap statistics form one BY/BH family. Historical leads are not forced in.\n\n"
        "Stability counts are descriptive: splits overlap and are not independent replications. "
        "There is no across-split p-value or FDR guarantee. Historical corpus inspection and "
        "U6 exchangeability remain limitations. Taxon annotation overlap does not measure "
        "patient physiology. See the recorded protocol and methods.\n"
    )
    return {
        "split_results.json": _json(records),
        "split_stability.tsv": _stability_table(records),
        "README.md": readme,
    }


def preflight_bundle_roots(output_dir: Path, paths: dict[str, Path]) -> None:
    """Keep sidecars outside every immutable input bundle, including symlink aliases."""
    destination = output_dir.resolve()
    for name, path in paths.items():
        if not name.endswith(":manifest"):
            continue
        root = path.resolve().parent
        if destination.is_relative_to(root) or root.is_relative_to(destination):
            raise ValueError(f"Output overlaps immutable upstream bundle {name}")


def execute(
    comparison_path: Path,
    output_dir: Path,
    protocol: Protocol,
    *,
    verify: bool = False,
) -> Path:
    """Write new sidecars, or verify every result by replaying from verified source inputs."""
    execution_identity = _execution_identity(protocol)
    manifest_path = output_dir / "study_split_run.json"
    saved = None
    if verify:
        saved = json.loads(manifest_path.read_text(encoding="utf-8"))
        if saved["provenance"]["parameters"] != asdict(protocol):
            raise ValueError("Replay parameters differ from the recorded run")
        timestamp = datetime.fromisoformat(saved["recorded_at"])
    else:
        if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
            raise ValueError("Study-split output must be a new or empty directory")
        timestamp = datetime.now(UTC)
    verified, paths = _inputs(comparison_path)
    preflight_bundle_roots(output_dir, paths)
    bundle = open_core_bundle(paths["core:manifest"])
    preflight_consumer_outputs(bundle, manifest_path, output_directory=output_dir)
    provenance = _provenance(paths, protocol, verified)
    scope_binding = _scope_binding(verified)
    if any(provenance[key] != value for key, value in execution_identity.items()):
        raise ValueError("Study-split implementation changed during verification")
    if saved is not None and saved["provenance"] != provenance:
        raise ValueError("Study-split sources, protocol, code or environment changed")
    if saved is None:
        commit, dirty = repository_state(REPO_ROOT)
        repository_context = {"code_commit": commit, "code_dirty": dirty}
    else:
        # Output publication or a documentation commit can change Git context
        # without changing any consumed bytes. Preserve the execution context;
        # the provenance inventory above is the replay identity (#906).
        repository_context = saved["repository_context"]
    documents = _documents(verified, paths, protocol, timestamp)
    if (
        _provenance(paths, protocol, verified) != provenance
        or _scope_binding(verified) != scope_binding
    ):
        raise ValueError("Study-split inputs or implementation changed during execution")
    output_hashes = {
        name: hashlib.sha256(content.encode("utf-8")).hexdigest()
        for name, content in sorted(documents.items())
    }
    record = {
        "protocol_version": VERSION,
        "recorded_at": timestamp.isoformat(),
        "repository_context": repository_context,
        "provenance": provenance,
        "scope_binding": scope_binding,
        "output_sha256": output_hashes,
    }
    if saved is not None:
        if saved != record:
            raise ValueError("Study-split manifest differs from full source replay")
        entries = list(output_dir.iterdir())
        if {path.name for path in entries} != {*documents, manifest_path.name} or any(
            path.is_symlink() or not path.is_file() for path in entries
        ):
            raise ValueError("Study-split output inventory differs from its declared files")
        for name, content in documents.items():
            actual = (output_dir / name).read_bytes()
            if (
                actual != content.encode("utf-8")
                or hashlib.sha256(actual).hexdigest() != output_hashes[name]
            ):
                raise ValueError(f"Study-split replay differs: {name}")
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
        for name, content in documents.items():
            (output_dir / name).write_bytes(content.encode("utf-8"))
        manifest_path.write_bytes(_json(record).encode("utf-8"))
    return manifest_path


def _load_scope_run(directory: Path, scope: str) -> dict[str, Any]:
    """Check saved bytes and their internal contracts; do not rerun the null model."""
    expected_files = {
        "study_split_run.json",
        "split_results.json",
        "split_stability.tsv",
        "README.md",
    }
    entries = list(directory.iterdir())
    if {entry.name for entry in entries} != expected_files or any(
        entry.is_symlink() or not entry.is_file() for entry in entries
    ):
        raise ValueError("Saved study-split output inventory differs from its declared files")
    buffers = {entry.name: entry.read_bytes() for entry in entries}
    saved = json.loads(buffers["study_split_run.json"])
    declared = _declaration(PROTOCOL_FILE)
    if saved["protocol_version"] != VERSION or saved["provenance"]["protocol"] != declared:
        raise ValueError("Both saved scope runs must follow the current v2 protocol")
    outputs = saved["output_sha256"]
    if set(outputs) != expected_files - {"study_split_run.json"}:
        raise ValueError("Saved study-split manifest does not bind the complete output inventory")
    for name, digest in outputs.items():
        if hashlib.sha256(buffers[name]).hexdigest() != digest:
            raise ValueError(f"Saved study-split output changed: {name}")
    binding = saved["scope_binding"]
    if binding["strain_rollup"] != scope or binding["role"] != _scope_role(scope):
        raise ValueError("Saved study-split scope or primary/sensitivity role is contradictory")
    contract = binding["comparison_contract"]
    if set(contract) != {
        "population",
        "direction_basis",
        "trait_key_policy",
        "parameters",
        "sources",
        "producer",
    } or set(contract["producer"]) != {
        "definition_sha256",
        "code_sha256",
        "lock_sha256",
        "runtime_sha256",
        "git_commit",
        "git_dirty",
    }:
        raise ValueError("Saved study-split comparison contract or producer metadata is incomplete")
    if (
        contract["population"] != scanner.TIER1_2
        or contract["direction_basis"] != ORIENTED_DIRECTION_BASIS
        or set(contract["sources"]) != set(scanner.CORE_SOURCES)
        or "strain_rollup" in contract["parameters"]
    ):
        raise ValueError(
            "Saved study-split population, direction or comparison contract is invalid"
        )
    provenance = saved["provenance"]
    parameters = Protocol(**provenance["parameters"])
    if not provenance["inputs"].get("comparison:manifest") or not provenance["implementation"]:
        raise ValueError("Saved study-split run lacks its comparison or implementation binding")
    if not provenance["environment"] or not contract["producer"]:
        raise ValueError("Saved study-split run lacks its producer or execution environment")
    records = json.loads(buffers["split_results.json"])
    if len(records) != parameters.splits:
        raise ValueError("Saved study-split count differs from the declared parameters")
    partitions = []
    for index, record in enumerate(records):
        partition = StudySplit(**record["split"])
        if partition.index != index:
            raise ValueError("Saved study-split indices are incomplete or reordered")
        partitions.append(asdict(partition))
        family = record["correction"]["family"]
        if (
            family["family_id"] != f"{VERSION}:split-{index:03}"
            or family["strain_rollup"] != scope
            or family["profile_unit"] != "count_profiles"
            or family["direction_basis"] != ORIENTED_DIRECTION_BASIS
        ):
            raise ValueError("Saved correction family contradicts its split scope or direction")
    stability = _stability_table(records).encode("utf-8")
    if stability != buffers["split_stability.tsv"]:
        raise ValueError("Saved stability table disagrees with its split result records")
    return {
        "manifest_sha256": hashlib.sha256(buffers["study_split_run.json"]).hexdigest(),
        "saved": saved,
        "records": records,
        "partitions": partitions,
        "stability": summarize_stability(records),
    }


def report_scopes(
    phenotype_dir: Path, off_dir: Path, output_dir: Path, *, verify: bool = False
) -> Path:
    """Pair both saved v2 scopes without selection, resampling or choosing a winner.

    This checks the recorded source identities and outputs, not the live upstream
    bundles. Use each run's full ``--verify`` mode for source-to-result replay.
    """
    preflight_bundle_roots(
        output_dir,
        {
            "phenotype:manifest": phenotype_dir / "study_split_run.json",
            "off:manifest": off_dir / "study_split_run.json",
        },
    )
    if (
        not verify
        and output_dir.exists()
        and (not output_dir.is_dir() or any(output_dir.iterdir()))
    ):
        raise ValueError("Study-split scopes output must be a new or empty directory")
    try:
        runs = {
            "phenotype": _load_scope_run(phenotype_dir, "phenotype"),
            "off": _load_scope_run(off_dir, "off"),
        }
        first, second = (runs[scope] for scope in ("phenotype", "off"))
        bindings = [run["saved"]["scope_binding"] for run in (first, second)]
        if (
            any(not binding["core_run_id"] for binding in bindings)
            or bindings[0]["core_run_id"] == bindings[1]["core_run_id"]
        ):
            raise ValueError("Each study-split scope must have its own recorded core run")
        if bindings[0]["comparison_contract"] != bindings[1]["comparison_contract"]:
            raise ValueError("Study-split scopes differ in their core comparison contract")
        for field in ("protocol", "parameters", "code_sha256", "implementation", "environment"):
            if first["saved"]["provenance"][field] != second["saved"]["provenance"][field]:
                raise ValueError(f"Study-split scopes differ in {field}")
        if first["partitions"] != second["partitions"]:
            raise ValueError("Study-split scopes differ in their study partitions")
    except (KeyError, TypeError) as error:
        raise ValueError("Incomplete saved study-split scope contract") from error
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=("strain_rollup", "role", *STABILITY_FIELDS),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    scopes = {}
    for scope, run in runs.items():
        saved = run["saved"]
        role = _scope_role(scope)
        for row in run["stability"]:
            writer.writerow(
                {
                    "strain_rollup": scope,
                    "role": role,
                    **row,
                    "conditions": json.dumps(row["conditions"], ensure_ascii=False),
                }
            )
        scopes[scope] = {
            "role": role,
            "core_run_id": saved["scope_binding"]["core_run_id"],
            "manifest_sha256": run["manifest_sha256"],
            "input_sha256": saved["provenance"]["inputs"],
            "output_sha256": saved["output_sha256"],
            "splits": len(run["records"]),
            "splits_with_candidates": sum(
                bool(row["selected_candidates"]) for row in run["records"]
            ),
        }
    documents = {
        "study_split_scopes.tsv": output.getvalue(),
        "README.md": (
            "# Study-split scopes, side by side\n\n"
            "`phenotype`: primary_phenotype; `off`: declared_sensitivity. Both are required. "
            "The scope amendment was made after results were known; do not select a winner.\n\n"
            "Saved stability counts are descriptive. Overlapping splits are dependent, "
            "and no across-split p-value or FDR claim is computed. Empty selections remain "
            "in the JSON split counts. Database trait annotations do not measure physiology.\n\n"
            "This report checks saved manifest/output hashes and their internal consistency. "
            "It does not rerun nulls or authenticate a coordinated rewrite of the manifests "
            "and dependent outputs. Replay each split run with --verify to check live sources.\n"
        ),
    }
    summary = {
        "protocol_version": VERSION,
        "protocol": first["saved"]["provenance"]["protocol"],
        "parameters": first["saved"]["provenance"]["parameters"],
        "comparison_contract": bindings[0]["comparison_contract"],
        "scope_order": ["phenotype", "off"],
        "scopes": scopes,
        "output_sha256": {
            name: hashlib.sha256(text.encode()).hexdigest() for name, text in documents.items()
        },
    }
    name = "study_split_scopes.json"
    documents[name] = _json(summary)
    if verify:
        entries = list(output_dir.iterdir())
        if {path.name for path in entries} != set(documents) or any(
            path.is_symlink() or not path.is_file() for path in entries
        ):
            raise ValueError("Study-split scopes report inventory differs")
        for filename, content in documents.items():
            if (output_dir / filename).read_bytes() != content.encode():
                raise ValueError(f"Study-split scopes report differs: {filename}")
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
        for filename, content in documents.items():
            with (output_dir / filename).open("xb") as stream:
                stream.write(content.encode())
    return output_dir / name


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-run", type=Path)
    parser.add_argument("--phenotype-run", type=Path, help="Saved v2 primary study-split run")
    parser.add_argument("--off-run", type=Path, help="Saved v2 sensitivity study-split run")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--splits", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--permutations", type=int, default=100000)
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Replay a single run's nulls, or check a saved paired report; no writes",
    )
    add_temporary_output_option(parser)
    args = parser.parse_args(argv)
    pairing = args.phenotype_run is not None or args.off_run is not None
    if pairing:
        if args.comparison_run is not None or args.phenotype_run is None or args.off_run is None:
            parser.error("Pairing requires --phenotype-run and --off-run, without --comparison-run")
        if (args.splits, args.seed, args.permutations) != (20, 20261023, 100000):
            parser.error("Scope pairing reads recorded parameters; do not supply numeric overrides")
    elif args.comparison_run is None:
        parser.error("Provide --comparison-run, or both saved scope runs")
    # Verification writes nothing; an execution's outputs are cited (#1087)
    if not args.verify:
        require_durable_output(parser, args.output_dir, allowed=args.allow_temporary_output)
    if pairing:
        result = report_scopes(
            args.phenotype_run, args.off_run, args.output_dir, verify=args.verify
        )
    else:
        comparison_run = args.comparison_run
        if not isinstance(comparison_run, Path):
            parser.error("Single-scope execution requires --comparison-run")
        result = execute(
            comparison_run,
            args.output_dir,
            Protocol(splits=args.splits, seed=args.seed, permutations=args.permutations),
            verify=args.verify,
        )
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
