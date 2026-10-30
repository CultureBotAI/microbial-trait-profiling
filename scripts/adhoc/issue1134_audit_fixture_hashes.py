"""Audit public fixture hashes before baseline admission (#1134, 2026-10-28).

Inputs: standard detect-secrets scan, preserved baseline, historical original run,
committed fixture and local Git history. Outputs: audited baseline and JSON report.
Only independently recomputed fixture provenance values are admitted; scientific
files and previous baseline entries are never changed. No network or replay.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from git import Repo

REPOSITORY = Path(__file__).resolve().parents[2]
FIXTURE = Path("tests/fixtures/trait_key_policy/historical_raw_core")
NAMESPACE = Path("/fixture/trait-policy-511c891")
SOURCE = Path("steps/identify_sources/provenance/source_identity.json")
HEX = re.compile(r"(?<![0-9a-f])[0-9a-f]{40,64}(?![0-9a-f])")


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def content_hash(value: object) -> str:
    return digest(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    )


def read_json(path: Path) -> Any:
    return json.loads(path.read_bytes())


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


class Proofs:
    """Public value proofs, without recording private original locations."""

    def __init__(self) -> None:
        self.values: dict[str, list[dict[str, str]]] = {}

    def add(self, value: str, category: str, source: str) -> None:
        self.values.setdefault(value, []).append({"category": category, "source": source})

    def match(self, expected: str, actual: str, category: str, source: str) -> None:
        require(expected == actual, f"Digest mismatch: {source}")
        self.add(expected, category, source)


def fixture_proofs(original_path: Path) -> Proofs:
    fixture = REPOSITORY / FIXTURE
    manifest = fixture / "core/run.json"
    run, original = read_json(manifest), read_json(original_path)
    generation = read_json(fixture / "generation.json")
    source, old_source = (
        read_json(fixture / "core" / SOURCE),
        read_json(original_path.parent / SOURCE),
    )
    proofs = Proofs()
    with Repo(REPOSITORY) as repository:
        producer = repository.commit(run["git_commit"])
        require(generation["producer_commit"] == producer.hexsha, "Producer commit differs")
        proofs.add(producer.hexsha, "historical_git_commit", "generation.producer_commit")
        for entry in run["code_inventory"]:
            payload = (producer.tree / entry["path"]).data_stream.read()
            require(len(payload) == entry["size_bytes"], f"Historical size: {entry['path']}")
            proofs.match(entry["sha256"], digest(payload), "historical_source_blob", entry["path"])
        for field, relative in (
            ("spec_sha256", "workflows/bugsigdb_core.yaml"),
            ("lock_sha256", "uv.lock"),
        ):
            proofs.match(
                run[field],
                digest((producer.tree / relative).data_stream.read()),
                "historical_source_blob",
                relative,
            )
    for field, value in (
        ("code_sha256", run["code_inventory"]),
        ("definition_sha256", run["definition"]),
    ):
        proofs.match(run[field], content_hash(value), "record_aggregate", f"run.{field}")
    external = {item["id"] for item in run["definition"]["products"] if item["external"]}
    runtime = json.loads(run["runtime_json"])
    runtime.pop("executable", None)
    fingerprint = {
        "definition": run["definition_sha256"],
        "spec": run["spec_sha256"],
        "code": run["code_sha256"],
        "lock": run["lock_sha256"],
        "runtime": runtime,
        "parameters": sorted((item["name"], item["value"]) for item in run["parameters"]),
        "inputs": sorted(
            (item["product_id"], item["sha256"])
            for item in run["artifacts"]
            if item["product_id"] in external
        ),
    }
    proofs.match(
        run["execution_fingerprint"],
        content_hash(fingerprint),
        "record_aggregate",
        "run.execution_fingerprint",
    )
    kinds = {item["id"]: item["storage_kind"] for item in run["definition"]["products"]}
    for artifact in run["artifacts"]:
        declared = Path(artifact["path"])
        path = (
            fixture / declared.relative_to(NAMESPACE)
            if declared.is_absolute()
            else fixture / "core" / declared
        )
        require(path.resolve().is_relative_to(fixture.resolve()), "Artifact escapes fixture")
        paths = sorted(path.rglob("*")) if kinds[artifact["product_id"]] == "collection" else [path]
        members: list[dict[str, Any]] = []
        for member in paths:
            require(not member.is_symlink(), "Fixture contains symlink")
            if member.is_dir():
                continue
            payload = member.read_bytes()
            value = digest(payload)
            proofs.add(value, "fixture_file", member.relative_to(REPOSITORY).as_posix())
            name = member.relative_to(path).as_posix() if path.is_dir() else member.name
            members.append({"path": name, "sha256": value, "size_bytes": len(payload)})
        require(members == artifact["inventory"], f"Inventory mismatch: {artifact['product_id']}")
        require(
            sum(item["size_bytes"] for item in members) == artifact["size_bytes"],
            f"Artifact size mismatch: {artifact['product_id']}",
        )
        actual = content_hash(members) if path.is_dir() else members[0]["sha256"]
        proofs.match(artifact["sha256"], str(actual), "artifact_inventory", artifact["product_id"])
    for name, item in source["consumed_files"].items():
        payload = (fixture / Path(item["path"]).relative_to(NAMESPACE)).read_bytes()
        require(len(payload) == item["size_bytes"], f"Consumed size: {name}")
        proofs.match(item["sha256"], digest(payload), "fixture_input", name)
    columns = (
        (fixture / "core/steps/taxa_matrix/derived/taxa_trait_matrix.tsv")
        .read_text()
        .splitlines()[0]
        .split("\t")[3:]
    )
    proofs.add(
        digest(json.dumps(columns, ensure_ascii=True, separators=(",", ":")).encode()),
        "feature_vocabulary",
        "taxa_trait_matrix.tsv ordered feature columns",
    )
    bindings = (
        (
            "helper_sha256",
            REPOSITORY / "scripts/adhoc/issue1134_generate_historical_core_fixture.py",
            "generator_file",
        ),
        ("normalized_manifest_sha256", manifest, "normalized_manifest"),
        ("original_manifest_sha256", original_path, "original_manifest"),
    )
    for field, path, category in bindings:
        proofs.match(generation[field], digest(path.read_bytes()), category, f"generation.{field}")
    for field, path in (
        ("before_sha256", original_path.parent / SOURCE),
        ("after_sha256", fixture / "core" / SOURCE),
    ):
        proofs.match(
            generation["source_identity_inventory_rebinding"][field],
            digest(path.read_bytes()),
            "source_identity_rebinding",
            f"generation.source_identity.{field}",
        )
    originals = {
        "run.manifest_path": original["manifest_path"],
        "run.spec_path": original["spec_path"],
        "runtime.executable": json.loads(original["runtime_json"])["executable"],
    }
    for artifact in original["artifacts"]:
        for field in ("path", "requested_path"):
            originals[f"artifact.{artifact['product_id']}.{field}"] = artifact.get(field)
    for name, item in old_source["consumed_files"].items():
        originals[f"source_identity.{name}.path"] = item["path"]
    for transformation in generation["location_transformations"]:
        field = transformation["field"]
        proofs.match(
            transformation["original_value_sha256"],
            digest(originals[field].encode()),
            "original_location_string",
            field,
        )
    matrix = generation["producer_matrix_inventory"]
    require(matrix in run["code_inventory"], "Generation matrix inventory differs")
    return proofs


def identity(entry: dict[str, Any]) -> tuple[str, str]:
    return entry["type"], entry["hashed_secret"]


def audit(
    previous: dict[str, Any], scanned: dict[str, Any], proofs: Proofs
) -> tuple[dict[str, Any], dict[str, Any]]:
    allowed = {
        str(FIXTURE / name)
        for name in (Path("core/run.json"), Path("core") / SOURCE, Path("generation.json"))
    }
    updated = copy.deepcopy(previous)
    admissions = []
    unadmitted = []
    missing: list[dict[str, str]] = []
    for filename, entries in previous["results"].items():
        found = {identity(entry) for entry in scanned["results"].get(filename, [])}
        missing.extend(
            {"filename": filename, "hashed_secret": entry["hashed_secret"]}
            for entry in entries
            if identity(entry) not in found
        )
    for filename, entries in scanned["results"].items():
        old = {identity(entry) for entry in previous["results"].get(filename, [])}
        for entry in entries:
            if identity(entry) in old:
                continue
            if filename not in allowed:
                unadmitted.append(
                    {
                        "filename": filename,
                        "type": entry["type"],
                        "hashed_secret": entry["hashed_secret"],
                    }
                )
                continue
            require(
                entry["type"] in {"Hex High Entropy String", "Secret Keyword"},
                "Unexpected detector type",
            )
            line = (REPOSITORY / filename).read_text().splitlines()[entry["line_number"] - 1]
            candidates = [
                value
                for value in HEX.findall(line)
                if hashlib.sha1(value.encode(), usedforsecurity=False).hexdigest()
                == entry["hashed_secret"]
            ]
            require(
                len(candidates) == 1,
                f"Finding is not a unique hexadecimal digest: {filename}:{entry['line_number']}",
            )
            value = candidates[0]
            require(
                value in proofs.values, f"Unproven checksum at {filename}:{entry['line_number']}"
            )
            admitted = dict(entry, is_secret=False)
            updated["results"].setdefault(filename, []).append(admitted)
            admissions.append(
                {
                    "filename": filename,
                    "line_number": entry["line_number"],
                    "hashed_secret": entry["hashed_secret"],
                    "public_value": value,
                    "proofs": proofs.values[value],
                }
            )
    for filename, entries in previous["results"].items():
        retained = {identity(entry): entry for entry in updated["results"][filename]}
        require(
            all(retained[identity(entry)] == entry for entry in entries), "Previous entry changed"
        )
    report = {
        "issue": 1134,
        "admitted_count": len(admissions),
        "admissions": admissions,
        "counts_by_file": dict(Counter(entry["filename"] for entry in admissions)),
        "counts_by_primary_proof": dict(
            Counter(entry["proofs"][0]["category"] for entry in admissions)
        ),
        "previous_entries_preserved": sum(map(len, previous["results"].values())),
        "scan_removed_entries_restored": missing,
        "unrelated_scan_addition_count": len(unadmitted),
        "unrelated_scan_additions_not_admitted": unadmitted,
    }
    return updated, report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous-baseline", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--original-manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    updated, report = audit(
        read_json(args.previous_baseline),
        read_json(args.baseline),
        fixture_proofs(args.original_manifest),
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    args.baseline.write_text(json.dumps(updated, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: value
                for key, value in report.items()
                if key not in {"admissions", "unrelated_scan_additions_not_admitted"}
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
