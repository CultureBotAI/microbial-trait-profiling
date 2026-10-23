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
from mtp.cellulose_decr.study_split import Protocol, run_splits, summarize_stability
from mtp.datamodel.multiple_testing_pydantic import SourceArtifact
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import load_profiles
from mtp.trait_screen import candidate_slots
from mtp.workflow.consumers import open_core_bundle, preflight_consumer_outputs
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
    "methods/study_split_selection_protocol.json",
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
    from mtp.multiplicity.sources import open_verified_sources

    verified = open_verified_sources(comparison_path)
    if (
        verified.core.contrast_population != scanner.TIER1_2
        or verified.core.direction_basis != ORIENTED_DIRECTION_BASIS
    ):
        raise ValueError("Study splits require the verified, disease-oriented Tier 1+2 population")
    if verified.scope not in {"off", "phenotype"}:
        raise ValueError("Study splits declare off primary and phenotype sensitivity only")
    paths = {f"core:{name}": verified.product(name) for name in PRODUCTS}
    paths["comparison:manifest"] = Path(verified.comparison.manifest_path).resolve()
    paths["core:manifest"] = Path(verified.comparison.core_manifest_path)
    for arm in ("jed", "noel", "ilyas"):
        paths[f"arm:{arm}:manifest"] = Path(getattr(verified.comparison, f"{arm}_manifest_path"))
    return verified, paths


def _provenance(paths: dict[str, Path], protocol: Protocol) -> dict[str, Any]:
    declared = json.loads(PROTOCOL_FILE.read_text(encoding="utf-8"))
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
    if (declared["primary_strain_rollup"], declared["sensitivity_strain_rollup"]) != (
        "off",
        "phenotype",
    ):
        raise ValueError("Committed protocol strain-rollup roles mismatch")
    # Selection and verification use transitive helpers in the canonical package.
    # Bind the complete package so a new helper cannot silently escape provenance.
    implementation = set(IMPLEMENTATIONS)
    implementation.update(
        path.relative_to(REPO_ROOT).as_posix() for path in (REPO_ROOT / "src" / "mtp").rglob("*.py")
    )
    return {
        "protocol": declared,
        "parameters": asdict(protocol),
        "follows_protocol": protocol == defaults,
        "inputs": {name: file_sha256(path) for name, path in sorted(paths.items())},
        "implementation": {name: file_sha256(REPO_ROOT / name) for name in sorted(implementation)},
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }


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
        SourceArtifact(identity=name, sha256=file_sha256(path))
        for name, path in sorted(paths.items())
    ]
    sources.append(
        SourceArtifact(
            identity="protocol:study-split-selection-v1", sha256=file_sha256(PROTOCOL_FILE)
        )
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
    stability = summarize_stability(records)
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output, fieldnames=STABILITY_FIELDS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in stability:
        writer.writerow({**row, "conditions": json.dumps(row["conditions"], ensure_ascii=False)})
    selected = sum(bool(record["selected_candidates"]) for record in records)
    role = "primary historical-lead scope" if verified.scope == "off" else "declared sensitivity"
    readme = (
        "# Internal study-split selection assessment\n\n"
        f"Strain rollup `{verified.scope}`: "
        f"{role}. "
        "Report both scopes separately; do not choose the favorable result.\n\n"
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
        "split_stability.tsv": output.getvalue(),
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
    provenance = _provenance(paths, protocol)
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
    if _provenance(paths, protocol) != provenance:
        raise ValueError("Study-split inputs or implementation changed during execution")
    output_hashes = {
        name: hashlib.sha256(content.encode("utf-8")).hexdigest()
        for name, content in sorted(documents.items())
    }
    record = {
        "protocol_version": "study-split-selection-v1",
        "recorded_at": timestamp.isoformat(),
        "repository_context": repository_context,
        "provenance": provenance,
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


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--splits", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--permutations", type=int, default=100000)
    parser.add_argument(
        "--verify", action="store_true", help="Replay all splits and nulls; no writes"
    )
    args = parser.parse_args(argv)
    result = execute(
        args.comparison_run,
        args.output_dir,
        Protocol(splits=args.splits, seed=args.seed, permutations=args.permutations),
        verify=args.verify,
    )
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
