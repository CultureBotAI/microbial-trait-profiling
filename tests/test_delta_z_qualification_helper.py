"""Reject incomplete producer evidence in the #1115 qualification helper (#1156)."""

from __future__ import annotations

import argparse
import hashlib
import json
import runpy
from pathlib import Path
from typing import Any

import pytest
from git import Actor, Repo

CHECK_TRIAL = runpy.run_path(
    str(Path(__file__).parents[1] / "scripts/adhoc/issue1115_compare_qualified_outputs.py")
)["check_trial"]
PRODUCER_FILES = {
    "src/mtp/delta_z_derived/tables.py": b"tables = 1\n",
    "src/mtp/delta_z_derived/stage.py": b"stage = 1\n",
    "src/mtp/extra.py": b"extra = 1\n",
    "src/mtp/nested/dependency.py": b"dependency = 1\n",
    "methods/three_test_table.md": b"Parent method\n",
    "methods/delta_z_derived.md": b"Derived method\n",
    "uv.lock": b"version = 1\n",
    "pyproject.toml": b"[project]\n",
    "Justfile": b"default:\n",
}


def digest_inventory(inventory: dict[str, str]) -> str:
    raw = json.dumps(inventory, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


@pytest.fixture
def trial(tmp_path: Path) -> tuple[dict[str, Any], argparse.Namespace, Path]:
    files = {**PRODUCER_FILES, "src/other/non_producer.py": b"outside = 1\n"}
    for name, raw in files.items():
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    author = Actor("Qualification fixture", "qualification@example.invalid")
    with Repo.init(tmp_path) as repository:
        repository.index.add(list(files))
        base = repository.index.commit("base", author=author, committer=author)
        issue = repository.index.commit(
            "1115", parent_commits=[base], head=False, author=author, committer=author
        )
        policy = repository.index.commit(
            "1130", parent_commits=[base], head=False, author=author, committer=author
        )
        merged = repository.index.commit(
            "trial", parent_commits=[issue, policy], author=author, committer=author
        )
    inventory = {name: hashlib.sha256(raw).hexdigest() for name, raw in PRODUCER_FILES.items()}
    record = {
        "code_commit": merged.hexsha,
        "code_dirty": False,
        "code_inventory": inventory,
        "implementation_sha256": digest_inventory(inventory),
        "methods_sha256": inventory["methods/delta_z_derived.md"],
        "lock_sha256": inventory["uv.lock"],
    }
    args = argparse.Namespace(
        issue1115_sha=issue.hexsha, pr1130_sha=policy.hexsha, trial_sha=merged.hexsha
    )
    return record, args, tmp_path


def test_complete_trial_inventory_passes(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
) -> None:
    record, args, root = trial
    result = CHECK_TRIAL(record, args, root)
    assert result["local_trial_merge_commit"] == args.trial_sha


@pytest.mark.parametrize("corruption", ["omitted_module", "aggregate_digest", "method_digest"])
def test_incomplete_or_stale_producer_identity_rejected(
    trial: tuple[dict[str, Any], argparse.Namespace, Path], corruption: str
) -> None:
    record, args, root = trial
    if corruption == "omitted_module":
        del record["code_inventory"]["src/mtp/extra.py"]
        record["implementation_sha256"] = digest_inventory(record["code_inventory"])
    else:
        field = "implementation_sha256" if corruption == "aggregate_digest" else "methods_sha256"
        record[field] = "0" * 64
    with pytest.raises(ValueError, match=r"inventory|implementation|method"):
        CHECK_TRIAL(record, args, root)
