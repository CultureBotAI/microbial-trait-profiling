"""The selection scan's two strain-rollup scopes, side by side (#1050).

Protocol v2 (``methods/cellulose_decr_selection_scan_protocol_v2.json``) makes phenotype
the scan's primary scope and off its declared sensitivity, by an amendment made on
2026-10-25 after the results of both scopes had been seen. It is a post-results choice of
primary scope, so v2 also requires the two scopes to be reported side by side, never one
alone.

A scan (``mtp.cellulose_decr.selection_scan``) reads one verified core, so one scope, and
its summary names the other scope as the one it must be reported beside. This module reads
the two scans' output directories, verifies each against its core
(``verify_consumer_record``) and against the protocol it names (``recorded_protocol``),
checks that each scan's scope is the one its bound core records (``check_bound_core``) and
that the two read different cores and differ in nothing but the scope, and writes them
together:

- ``selection_scan_scopes.tsv``: each scan's trait table as written, primary scope first,
  with the scope and its role in front;
- ``selection_scan_scopes.json``: the declaration and the amendment that made it, and for
  each scope its role, core run, counts, lead-trait results and the hashes of the files
  read.

It runs no null and recomputes nothing. The counts are annotation counts of curated
signature taxa, not measurements of patients, cells or pathway activity
(methods/biological_interpretation.md). Every p-value is nominal.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from mtp.cellulose_decr.analysis import write_json
from mtp.cellulose_decr.selection_scan import (
    CONSUMER,
    CONSUMER_RECORD,
    PROTOCOL_PATH,
    SIDE_BY_SIDE,
    STRAIN_ROLLUP,
    Protocol,
    core_scope,
    load_protocol,
    protocol_record,
    recorded_protocol,
    repository_path,
)
from mtp.cellulose_decr.selection_scan import OUTPUTS as SCAN_OUTPUTS
from mtp.paths import file_sha256
from mtp.workflow.output_roots import add_temporary_output_option, require_durable_output

if TYPE_CHECKING:
    from mtp.datamodel.workflow_pydantic import WorkflowRun

OUTPUTS = {"table": "selection_scan_scopes.tsv", "summary": "selection_scan_scopes.json"}
# The scan outputs the report reads. Each must be one the scan bound to its core.
READ = ("summary", "provenance", "traits")


@dataclass(frozen=True)
class ScanRun:
    """One scan's outputs, verified against its core and the protocol it names.

    ``traits`` is the trait table as written (every value as text); ``files`` holds the
    SHA-256 of each file read; ``core_run_id`` is the core run its consumer record binds.
    """

    directory: Path
    summary: dict[str, Any]
    provenance: dict[str, Any]
    traits: pd.DataFrame
    protocol: Protocol
    files: dict[str, str]
    core_run_id: str

    @property
    def scope(self) -> str:
        return str(self.summary[STRAIN_ROLLUP]["scope"])

    @property
    def role(self) -> str:
        return str(self.summary[STRAIN_ROLLUP]["role"])


def check_bound_core(
    directory: Path,
    summary: Mapping[str, Any],
    provenance: Mapping[str, Any],
    run: WorkflowRun,
    manifest_sha256: str,
    protocol: Protocol,
) -> None:
    """Refuse a scan whose records do not describe the core its consumer record binds.

    The summary and provenance must name that core run and its manifest hash, and the
    recorded scope, with where it came from, must be what ``core_scope`` reads from that
    core under the recorded protocol. A scan's own summary is not trusted for its scope
    (#1079).
    """
    named = (
        summary["core_run"]["run_id"],
        provenance["core_run"]["run_id"],
        provenance["core_run"]["manifest_sha256"],
    )
    if named != (run.run_id, run.run_id, manifest_sha256):
        raise ValueError(
            f"{directory}: the records name core run {named[0]!r} and {named[1]!r} with "
            f"manifest {named[2]!r}, but its consumer record binds {run.run_id!r} with "
            f"manifest {manifest_sha256!r}"
        )
    scope, source = core_scope(run, protocol)
    rollup = summary[STRAIN_ROLLUP]
    recorded = (rollup["scope"], rollup.get("scope_source"))
    expected = (scope, source if protocol.side_by_side else None)
    if recorded != expected:
        raise ValueError(
            f"{directory}: the records give strain_rollup and its source as {recorded}, but "
            f"the bound core's are {expected}"
        )


def load_run(directory: Path) -> ScanRun:
    """Read a scan's outputs once its consumer record, its core and its protocol verify.

    The consumer record is verified end to end (the core run, and the hash of every
    output), every file read must be one of its outputs, the summary and provenance must
    agree with the protocol version they name, and they must describe the bound core and
    its scope (``check_bound_core``).
    """
    from mtp.datamodel.workflow_pydantic import WorkflowRun
    from mtp.workflow.consumers import verify_consumer_record

    record = verify_consumer_record(directory / CONSUMER_RECORD)
    if record.consumer != CONSUMER:
        raise ValueError(f"{directory} holds a record of {record.consumer}, not of {CONSUMER}")
    bound = {Path(output.path) for output in record.outputs}
    paths = {name: directory / SCAN_OUTPUTS[name] for name in READ}
    unbound = sorted(path.name for path in paths.values() if path.resolve() not in bound)
    if unbound:
        raise ValueError(f"{directory}: {unbound} are not bound to the scan's core record")
    summary = json.loads(paths["summary"].read_text(encoding="utf-8"))
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    traits = pd.read_csv(paths["traits"], sep="\t", dtype=str, keep_default_na=False)
    files = {
        path.name: file_sha256(path) for path in [*paths.values(), directory / CONSUMER_RECORD]
    }
    protocol = recorded_protocol(summary, provenance)
    # The bound core, from the manifest bytes whose hash the record holds.
    manifest = Path(record.run_manifest_path).read_bytes()
    if hashlib.sha256(manifest).hexdigest() != record.run_manifest_sha256:
        raise ValueError(f"{directory}: the bound core's manifest changed after verification")
    run = WorkflowRun.model_validate_json(manifest)
    if run.run_id != record.run_id:
        raise ValueError(f"{directory}: the bound manifest is of {run.run_id}, not {record.run_id}")
    check_bound_core(directory, summary, provenance, run, record.run_manifest_sha256, protocol)
    return ScanRun(directory, summary, provenance, traits, protocol, files, record.run_id)


def comparable(run: ScanRun) -> dict[str, Any]:
    """What two scans must share to differ in nothing but the scope."""
    core = run.provenance["core_run"]
    return {
        "protocol": run.provenance["protocol"],
        "code_commit": run.provenance["code_commit"],
        "scripts": run.provenance["scripts"],
        "rule": run.summary["rule"],
        "core population": (core["contrast_population"], core["direction_basis"]),
        "core git_commit": (core["git_commit"], core["git_dirty"]),
        "core parameters": {
            name: value for name, value in core["parameters"].items() if name != STRAIN_ROLLUP
        },
        "core sources": core["sources"],
    }


def pair_scopes(runs: Sequence[ScanRun], protocol: Protocol | None = None) -> list[ScanRun]:
    """The scans in declared order, primary first, once they can stand side by side.

    There must be exactly one scan of each declared scope, each of its own core run, all
    under the current protocol, with the same protocol file, code, rule, core population,
    core commit, core sources and core parameters other than the scope.
    """
    protocol = protocol or load_protocol()
    declared = [protocol.scope_primary, *protocol.scope_sensitivity]
    if sorted(run.scope for run in runs) != sorted(declared):
        raise ValueError(
            f"the scopes report needs one scan of each declared scope {declared}, "
            f"not {[run.scope for run in runs]}"
        )
    if len({run.core_run_id for run in runs}) != len(runs):
        raise ValueError(
            f"the scans read the same core run {[run.core_run_id for run in runs]}; each "
            "scope needs its own core"
        )
    for run in runs:
        if run.protocol.version != protocol.version:
            raise ValueError(
                f"{run.directory.name} follows {run.protocol.version}; the scopes report "
                f"needs both under {protocol.version}"
            )
    ordered = sorted(runs, key=lambda run: declared.index(run.scope))
    first, *others = (comparable(run) for run in ordered)
    for run, other in zip(ordered[1:], others, strict=True):
        for name, value in first.items():
            if other[name] != value:
                raise ValueError(f"the {ordered[0].scope} and {run.scope} scans differ in {name}")
    return ordered


def scope_amendment(path: Path = PROTOCOL_PATH) -> dict[str, Any]:
    """The amendment of ``path`` that declared the strain-rollup scopes, without its sets."""
    amendments = json.loads(path.read_text(encoding="utf-8"))["amendments"]
    latest = [amendment for amendment in amendments if STRAIN_ROLLUP in amendment["sets"]][-1]
    return {key: value for key, value in latest.items() if key != "sets"}


def side_by_side(
    ordered: Sequence[ScanRun], protocol: Protocol, code_commit: str
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """The scopes' trait tables in one table, and the report that goes with it."""
    frames = []
    for run in ordered:
        frame = run.traits.copy()
        frame.insert(0, "role", run.role)
        frame.insert(0, STRAIN_ROLLUP, run.scope)
        frames.append(frame)
    table = pd.concat(frames, ignore_index=True)
    scopes = []
    for run in ordered:
        lead = [row for row in run.summary["traits"] if row["trait_id"] == protocol.lead_trait]
        scopes.append(
            {
                "scope": run.scope,
                "role": run.role,
                "scope_source": run.summary[STRAIN_ROLLUP]["scope_source"],
                "directory": run.directory.name,
                "core_run_id": run.core_run_id,
                "follows_protocol": run.summary["follows_protocol"],
                "counts": run.summary["counts"],
                "lead_triple": run.summary["lead_triple"],
                "lead_trait": lead[0] if lead else None,
                "files": run.files,
            }
        )
    seen = " and ".join(protocol.results_seen)
    report = {
        "protocol": protocol.version,
        "follows_protocol": all(run.summary["follows_protocol"] for run in ordered),
        STRAIN_ROLLUP: {
            "declared_primary": protocol.scope_primary,
            "declared_sensitivity": list(protocol.scope_sensitivity),
            "declared_after_results": list(protocol.results_seen),
            "reported": SIDE_BY_SIDE,
            "amendment": scope_amendment(),
        },
        "rule": ordered[0].summary["rule"],
        "scopes": scopes,
        "notes": [
            f"{protocol.scope_primary} is the primary scope and "
            f"{', '.join(protocol.scope_sensitivity)} the declared sensitivity "
            f"({protocol.version}). The roles were declared after the {seen} results had "
            "been seen, a post-results choice of primary scope, so the scopes are reported "
            "together and neither alone.",
            "Each scope's rows are its scan's trait table as written; nothing is recomputed.",
            "Every p-value is nominal; no multiplicity correction is applied here.",
            "Counts describe curated signature taxa and their annotations, not physiology.",
        ],
        "provenance": {
            "generated_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "code_commit": code_commit,
            "script": repository_path(Path(__file__).resolve()),
            "script_sha256": file_sha256(Path(__file__).resolve()),
            "protocol": protocol_record(protocol.version),
        },
    }
    return table, report


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument(
        "--scan",
        type=Path,
        action="append",
        required=True,
        help="Output directory of a selection scan; give one scan of each declared scope.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="New or empty directory for the report.",
    )
    parser.add_argument(
        "--code-commit",
        default="unrecorded",
        help="Revision of this script, e.g. $(git rev-parse HEAD) with +dirty when modified.",
    )
    add_temporary_output_option(parser)
    args = parser.parse_args(argv)
    require_durable_output(parser, args.output_dir, allowed=args.allow_temporary_output)
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    out = args.output_dir
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"{out} must be a new or empty directory")
    protocol = load_protocol()
    if len(args.scan) != 1 + len(protocol.scope_sensitivity):
        declared = [protocol.scope_primary, *protocol.scope_sensitivity]
        raise ValueError(
            f"the scopes report needs one scan of each declared scope {declared}, "
            f"not {len(args.scan)} scans"
        )
    ordered = pair_scopes([load_run(directory) for directory in args.scan], protocol)
    table, report = side_by_side(ordered, protocol, args.code_commit)
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / OUTPUTS["table"], sep="\t", index=False, lineterminator="\n")
    write_json(report, out / OUTPUTS["summary"])
    print(f"wrote {len(OUTPUTS)} files to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
