"""Read-only parent audit and exclusive CLI publication contracts (#1073, 2026-10-26)."""

from __future__ import annotations

import importlib
import json
import os
import subprocess
import sys
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any
from unittest.mock import Mock

import pandas as pd
import pytest

from mtp.paths import file_sha256
from test_three_test_wilcoxon_audit import Parent
from test_three_test_wilcoxon_audit import parent as parent

ROOT = Path(__file__).resolve().parents[1]
PARENT_FILES = (
    "three_test_run.json",
    "per_trait_evidence.tsv",
    "composition_pair_delta_z.tsv",
    "composition_study_delta_z.tsv",
)
CAPTURE_SOURCE = '''"""Capture actual Justfile argv for #1073 (2026-10-26).
Inputs: uv argv and the explicit capture destination; output: one JSON list.
"""
import json
import os
import sys
from pathlib import Path
Path(os.environ["MTP_AUDIT_CAPTURE"]).write_text(json.dumps(sys.argv[1:]), encoding="utf-8")
'''


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    """Keep generated command-interception helpers and audit outputs in this repo."""
    root = ROOT / ".tmp" / "wilcoxon-cli-tests"
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as directory:
        yield Path(directory)


@pytest.fixture
def cli() -> ModuleType:
    return importlib.import_module("mtp.three_test_table.wilcoxon_audit_cli")


@dataclass(frozen=True)
class Inputs:
    manifest: Path
    pin: str
    hashes: dict[Path, str]
    audit: Mock


@pytest.fixture
def inputs(tmp_path: Path, cli: ModuleType, monkeypatch: pytest.MonkeyPatch) -> Inputs:
    parent = tmp_path / "parent"
    parent.mkdir()
    for name in PARENT_FILES:
        (parent / name).write_text("{}\n" if name.endswith("json") else "trait\nexample\n")
    hashes = {path: file_sha256(path) for path in parent.iterdir()}
    audit = Mock(
        return_value=SimpleNamespace(
            table=pd.DataFrame({"trait": ["METPO:test"], "archived_p": [0.08], "new_p": [0.125]}),
            summary={
                "family_size": 1,
                "verification_scope": "pinned_parent_outputs_and_legacy_reproduction_only",
            },
            input_hashes=hashes,
        )
    )
    monkeypatch.setattr(cli, "audit_parent", audit)
    producer = {
        "code_commit": "c" * 40,
        "code_dirty": False,
        "code_inventory": {"src/mtp/example.py": "a" * 64},
        "implementation_sha256": "b" * 64,
        "methods_sha256": "d" * 64,
        "lock_sha256": "e" * 64,
        "environment": {"python": "test", "scipy": "test"},
    }
    monkeypatch.setattr(cli.stage, "_producer", lambda _requested: producer.copy())
    return Inputs(parent / PARENT_FILES[0], hashes[parent / PARENT_FILES[0]], hashes, audit)


def run(cli: ModuleType, inputs: Inputs, output: Path) -> dict[str, Any]:
    record: dict[str, Any] = cli.run_audit(
        inputs.manifest, output, expected_sha256=inputs.pin, scope="off", code_commit="c" * 40
    )
    return record


def test_actual_just_recipe_forwards_spaced_paths_and_shell_literals(tmp_path: Path) -> None:
    tools = tmp_path / "intercepted tools"
    tools.mkdir()
    writer = tools / "capture argv.py"
    writer.write_text(CAPTURE_SOURCE, encoding="utf-8")
    command = tools / "uv"
    command.write_text(
        "#!/bin/sh\n# Intercept uv for #1073 recipe validation (2026-10-26).\n"
        'exec "$MTP_AUDIT_PYTHON" "$MTP_AUDIT_WRITER" "$@"\n',
        encoding="utf-8",
    )
    command.chmod(0o755)
    capture = tmp_path / "argv.json"
    arguments = [
        "--parent",
        "parent with spaces/$(printf literal); three_test_run.json",
        "--parent-sha256",
        "a" * 64,
        "--scope",
        "off",
        "--output",
        "output with spaces/`literal`;$not_expanded",
        "--code-commit",
        "b" * 40,
    ]
    process = subprocess.run(
        ["just", "--justfile", str(ROOT / "Justfile"), "three-test-wilcoxon-audit", *arguments],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
        env={
            "PATH": os.pathsep.join([str(tools), os.environ.get("PATH", os.defpath)]),
            "MTP_AUDIT_PYTHON": sys.executable,
            "MTP_AUDIT_WRITER": str(writer),
            "MTP_AUDIT_CAPTURE": str(capture),
        },
    )
    assert process.returncode == 0, process.stdout + process.stderr
    assert json.loads(capture.read_text()) == [
        "run",
        "python",
        "-m",
        "mtp.three_test_table.wilcoxon_audit_cli",
        *arguments,
    ]


def test_cli_forwards_required_pins_and_optional_commit(
    cli: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    execute = Mock(return_value={"run_id": "audit-fixture", "summary": {"family_size": 1}})
    monkeypatch.setattr(cli, "run_audit", execute)
    assert (
        cli.main(
            [
                "--parent",
                "parent",
                "--parent-sha256",
                "a" * 64,
                "--scope",
                "phenotype",
                "--output",
                "output",
                "--code-commit",
                "b" * 40,
            ]
        )
        == 0
    )
    execute.assert_called_once_with(
        Path("parent"),
        Path("output"),
        expected_sha256="a" * 64,
        scope="phenotype",
        code_commit="b" * 40,
    )


def test_publication_records_relative_inputs_and_own_producer_without_source_writes(
    cli: ModuleType,
    inputs: Inputs,
    tmp_path: Path,
) -> None:
    before = {path: path.read_bytes() for path in inputs.hashes}
    output = tmp_path / "new parent" / "audit"
    record = run(cli, inputs, output)
    inputs.audit.assert_called_once_with(inputs.manifest, inputs.pin, scope="off")
    assert {path.name for path in output.iterdir()} == {
        "wilcoxon_audit.tsv",
        "wilcoxon_audit_run.json",
    }
    text = (output / "wilcoxon_audit_run.json").read_text()
    assert json.loads(text) == record
    assert str(tmp_path) not in text
    assert record["parent"] == {
        "manifest": "three_test_run.json",
        "sha256": inputs.pin,
        "scope": "off",
    }
    assert set(record["inputs"]) == set(PARENT_FILES)
    assert record["code_commit"] == "c" * 40
    assert "methods/three_test_wilcoxon_audit.md" in record["code_inventory"]
    assert record["outputs"]["wilcoxon_audit.tsv"]["sha256"] == file_sha256(
        output / "wilcoxon_audit.tsv"
    )
    assert {path: path.read_bytes() for path in inputs.hashes} == before


def test_audit_failure_never_creates_output_ancestors(
    cli: ModuleType, inputs: Inputs, tmp_path: Path
) -> None:
    inputs.audit.side_effect = ValueError("historical p-value reproduction failed")
    output = tmp_path / "uncreated" / "nested" / "audit"
    with pytest.raises(ValueError, match="reproduction"):
        run(cli, inputs, output)
    assert not (tmp_path / "uncreated").exists()


@pytest.mark.parametrize("overlap", ["inside", "ancestor"])
def test_source_overlap_is_rejected_before_audit(
    cli: ModuleType, inputs: Inputs, tmp_path: Path, overlap: str
) -> None:
    output = inputs.manifest.parent / "audit" if overlap == "inside" else tmp_path
    with pytest.raises((ValueError, FileExistsError)):
        run(cli, inputs, output)
    inputs.audit.assert_not_called()


def test_temporary_output_rejected_before_audit(
    cli: ModuleType, inputs: Inputs, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    temporary = tmp_path / "temporary"
    monkeypatch.setenv("TMPDIR", str(temporary))
    with pytest.raises(ValueError, match="durable"):
        run(cli, inputs, temporary / "audit")
    inputs.audit.assert_not_called()
    assert not temporary.exists()


def test_input_mutation_after_audit_prevents_any_output_parent(
    cli: ModuleType, inputs: Inputs, tmp_path: Path
) -> None:
    original = inputs.audit.return_value

    def changed(*_args: Any, **_kwargs: Any) -> Any:
        (inputs.manifest.parent / "composition_study_delta_z.tsv").write_text("changed\n")
        return original

    inputs.audit.side_effect = changed
    output = tmp_path / "uncreated" / "audit"
    with pytest.raises(ValueError, match=r"(?i)input|source"):
        run(cli, inputs, output)
    assert not output.parent.exists()


@pytest.mark.parametrize("mutation", ["source", "producer"])
def test_change_during_serialization_refuses_publication_and_cleans_staging(
    cli: ModuleType,
    inputs: Inputs,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
) -> None:
    serialize = pd.DataFrame.to_csv
    original_producer = cli.stage._producer("c" * 40)

    def changed(frame: pd.DataFrame, *args: Any, **kwargs: Any) -> Any:
        result = serialize(frame, *args, **kwargs)
        if mutation == "source":
            inputs.manifest.write_text("changed\n")
        else:
            monkeypatch.setattr(
                cli.stage,
                "_producer",
                lambda _requested: {**original_producer, "code_commit": "f" * 40},
            )
        return result

    monkeypatch.setattr(pd.DataFrame, "to_csv", changed)
    output = tmp_path / "uncreated" / "nested" / "audit"
    with pytest.raises(ValueError, match=r"(?i)input|source|producer"):
        run(cli, inputs, output)
    assert not output.exists()
    assert not (tmp_path / "uncreated").exists()
    assert not list(tmp_path.rglob(".audit-staging-*"))


def test_competing_destination_is_preserved(
    cli: ModuleType, inputs: Inputs, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    publish = cli.stage._publish_directory
    output = tmp_path / "audit"

    def competed(staging: Path, destination: Path) -> None:
        destination.mkdir()
        (destination / "owner.txt").write_text("keep me\n")
        publish(staging, destination)

    monkeypatch.setattr(cli.stage, "_publish_directory", competed)
    with pytest.raises(FileExistsError):
        run(cli, inputs, output)
    assert (output / "owner.txt").read_text() == "keep me\n"
    assert not list(tmp_path.glob(".audit-staging-*"))


@pytest.mark.parametrize("alias", ["none", "directory", "manifest"])
def test_real_parent_audit_integrates_with_atomic_cli_publication(
    cli: ModuleType,
    parent: Parent,
    tmp_path: Path,
    alias: str,
) -> None:
    pin = parent.seal()
    before = {path: path.read_bytes() for path in parent.root.iterdir()}
    supplied = parent.manifest
    if alias == "directory":
        supplied = tmp_path / "parent alias"
        supplied.symlink_to(parent.root, target_is_directory=True)
    elif alias == "manifest":
        supplied = tmp_path / "manifest alias" / "pinned.json"
        supplied.parent.mkdir()
        supplied.symlink_to(parent.manifest)
    output = tmp_path / "audit"
    record = cli.run_audit(supplied, output, expected_sha256=pin, scope="off")
    assert record["summary"]["family_size"] == len(parent.values)
    assert record["summary"]["status_changes"] == 1
    assert (
        record["summary"]["verification_scope"]
        == "pinned_parent_outputs_and_legacy_reproduction_only"
    )
    assert set(record["inputs"]) == set(PARENT_FILES)
    assert str(tmp_path) not in (output / "wilcoxon_audit_run.json").read_text()
    assert {path: path.read_bytes() for path in parent.root.iterdir()} == before


def test_middle_manifest_alias_parent_is_protected_before_audit(
    cli: ModuleType,
    inputs: Inputs,
    tmp_path: Path,
) -> None:
    middle = tmp_path / "middle" / "relay.json"
    supplied = tmp_path / "supplied" / "manifest.json"
    middle.parent.mkdir()
    supplied.parent.mkdir()
    middle.symlink_to(inputs.manifest)
    supplied.symlink_to(middle)
    output = middle.parent / "new" / "audit"
    with pytest.raises(ValueError, match="disjoint"):
        cli.run_audit(supplied, output, expected_sha256=inputs.pin, scope="off")
    inputs.audit.assert_not_called()
    assert not output.parent.exists()


def test_identical_byte_manifest_retarget_is_refused(
    cli: ModuleType,
    inputs: Inputs,
    tmp_path: Path,
) -> None:
    alias = tmp_path / "alias" / "manifest.json"
    alias.parent.mkdir()
    alias.symlink_to(inputs.manifest)
    other = tmp_path / "other" / "three_test_run.json"
    other.parent.mkdir()
    other.write_bytes(inputs.manifest.read_bytes())
    audited = inputs.audit.return_value
    audited.input_hashes = {**inputs.hashes, alias: inputs.pin}

    def retarget(*_args: Any, **_kwargs: Any) -> Any:
        alias.unlink()
        alias.symlink_to(other)
        return audited

    inputs.audit.side_effect = retarget
    output = tmp_path / "uncreated" / "audit"
    with pytest.raises(ValueError, match=r"(?i)input|source"):
        cli.run_audit(alias, output, expected_sha256=inputs.pin, scope="off")
    assert not output.parent.exists()
