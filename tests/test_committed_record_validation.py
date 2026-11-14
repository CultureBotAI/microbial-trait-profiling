"""Guard the #1316 committed-record evidence helper, without reading real runs."""

from __future__ import annotations

import hashlib
import importlib.util
import io
import sys
import tarfile
from collections import Counter
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from git import Actor, Repo
from pydantic import BaseModel


@pytest.fixture(scope="module")
def helper() -> ModuleType:
    path = Path(__file__).parents[1] / "scripts/adhoc/issue1316_validate_committed.py"
    spec = importlib.util.spec_from_file_location("_committed_validation_test_helper", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    ("path", "target"),
    [
        ("tests/fixtures/trait_key_policy/historical_raw_core/core/run.json", "WorkflowRun"),
        (
            "tests/fixtures/trait_key_policy/historical_union_core.tar.gz!core/steps/ingest/studies/bsdb_1.yaml",
            "Study",
        ),
        (
            "tests/fixtures/trait_key_policy/historical_union_core.tar.gz!core/steps/export_tables/derived/experiments.tsv",
            "ExperimentTableRow",
        ),
        ("tests/fixtures/workflow/merged-kg_nodes.tsv", "Node"),
    ],
)
def test_direct_records_receive_explicit_classes(
    helper: ModuleType, path: str, target: str
) -> None:
    assert helper.classification(path)["class"] == target


def test_unknown_fixture_is_not_silently_excluded(helper: ModuleType) -> None:
    assert helper.classification("tests/fixtures/unreviewed.xyz")["kind"] == "unknown"


def test_archive_root_generation_projects_one_producer_inventory(helper: ModuleType) -> None:
    name = "tests/fixtures/trait_key_policy/historical_union_core.tar.gz!generation.json"
    mapping = helper.classification(name)
    assert mapping["kind"] == "instance"
    assert mapping["schema"] == "workflow"
    assert mapping["class"] == "FileInventoryEntry"
    assert mapping["encoding"] == "producer_inventory"
    producer_bytes = b"fixture producer"
    inventory = {
        "path": "src/mtp/taxa_trait_matrix.py",
        "sha256": hashlib.sha256(producer_bytes).hexdigest(),
        "size_bytes": len(producer_bytes),
    }
    record = {"producer_matrix_inventory": inventory, "fixture_kind": "synthetic"}
    assert list(helper.document_records(mapping, record)) == [
        ("/producer_matrix_inventory", inventory)
    ]


def test_physical_wide_tables_do_not_fabricate_logical_fields(helper: ModuleType) -> None:
    name = (
        "tests/fixtures/trait_key_policy/historical_raw_core/"
        "core/steps/count_profiles/derived/study_trait_profiles_wide.tsv"
    )
    mapping = helper.classification(name)
    assert mapping["kind"] == "noninstance"
    assert "No missing units" in mapping["reason"]
    assert "class" not in mapping


def test_new_constraint_detected_in_already_invalid_record(helper: ModuleType) -> None:
    inherited = {"path": ["other"], "validator": "required", "schema_path": ["required"]}
    pattern = {
        "path": ["id"],
        "validator": "pattern",
        "schema_path": ["properties", "id", "pattern"],
        "message": "new refusal",
        "raw_value": "bsdb:1\n",
    }
    delta = helper.constraint_delta(
        {"valid": False, "errors": [inherited]}, {"valid": False, "errors": [inherited, pattern]}
    )
    assert delta["valid_to_invalid"] is False
    assert delta["newly_failing_constraints"] == [pattern]
    assert delta["inherited_failure_constraints"] == 1


def test_message_change_is_not_a_new_constraint(helper: ModuleType) -> None:
    old = {
        "path": ["id"],
        "validator": "pattern",
        "schema_path": ["pattern"],
        "message": "old regex",
    }
    new = dict(old, message="new regex")
    assert (
        helper.constraint_delta(
            {"valid": False, "errors": [old]}, {"valid": False, "errors": [new]}
        )["newly_failing_constraints"]
        == []
    )


def test_schema_integer_conversion_keeps_patterned_strings_and_empty_fields(
    helper: ModuleType,
) -> None:
    side = object.__new__(helper.ValidationSide)
    side.schema_documents = {
        "core_products": {
            "classes": {
                "Example": {
                    "attributes": {
                        "n": {"range": "integer"},
                        "id": {"pattern": "^id$"},
                        "label": {"range": "string"},
                    }
                }
            }
        }
    }
    row = {"n": "2", "id": "id\n", "label": ""}
    assert side.typed_table_row("core_products", "Example", row) == {
        "n": 2,
        "id": "id\n",
        "label": "",
    }
    assert row == {"n": "2", "id": "id\n", "label": ""}


class Nested(BaseModel):
    value: str = "default"
    optional: str | None = None


class Container(BaseModel):
    nested: Nested


def test_serialization_retains_nested_defaults_and_nulls(helper: ModuleType) -> None:
    value = Container(nested=Nested())
    assert (
        helper.stable_bytes(value.model_dump(mode="json"))
        == b'{"nested":{"optional":null,"value":"default"}}'
    )


def _archive(path: Path, member: tarfile.TarInfo, content: bytes) -> None:
    with tarfile.open(path, "w:gz") as archive:
        member.size = len(content)
        archive.addfile(member, io.BytesIO(content))


@pytest.mark.parametrize("name", ["../escape.json", "/absolute.json"])
def test_archive_paths_refused_without_extraction(
    helper: ModuleType, tmp_path: Path, name: str
) -> None:
    path = tmp_path / "fixture.tar.gz"
    _archive(path, tarfile.TarInfo(name), b"{}")
    with pytest.raises(ValueError, match="Unsafe or duplicate"):
        helper.archive_members(path)


def test_archive_links_refused(helper: ModuleType, tmp_path: Path) -> None:
    path = tmp_path / "fixture.tar.gz"
    member = tarfile.TarInfo("link.json")
    member.type = tarfile.SYMTYPE
    member.linkname = "target"
    _archive(path, member, b"")
    with pytest.raises(ValueError, match="not a regular file"):
        helper.archive_members(path)


def test_tsv_rejects_ragged_record(helper: ModuleType) -> None:
    with pytest.raises(ValueError, match="Ragged"):
        list(helper.table_records(io.BytesIO(b"id\tn\na\t1\textra\n"), "fixture.tsv"))


def test_nested_pattern_error_is_preserved(helper: ModuleType) -> None:
    class Error:
        def __init__(self) -> None:
            self.absolute_path = ["id"]
            self.absolute_schema_path = ["properties", "id", "pattern"]
            self.validator = "pattern"
            self.validator_value = "^id$"
            self.message = "refused"
            self.instance = "id\n"
            self.context: list[Any] = []

    parent = Error()
    parent.validator = "anyOf"
    parent.context = [Error()]
    observer = helper.AllConstraints()
    observer._collect(parent)
    assert [error["validator"] for error in observer.errors] == ["anyOf", "pattern"]
    assert observer.errors[1]["raw_value"] == "id\n"


@pytest.mark.parametrize(
    "path",
    [
        "trait-visualization/new_consumer_record.json",
        "methods/table1-native-results/new_contract.tsv.gz",
        "tests/fixtures/trait_key_policy/historical_raw_core/core/steps/export_tables/derived/new_contract.tsv",
        ".github/new_domain_record.yaml",
        "tests/fixtures/ora_figure_1505/new_run.json",
        "tests/fixtures/trait_key_policy/historical_union_core.tar.gz!surprise/ingest_log.json",
    ],
)
def test_new_paths_require_reviewed_mapping(helper: ModuleType, path: str) -> None:
    assert helper.classification(path)["kind"] == "unknown"


@pytest.mark.parametrize(
    "state", ["unchanged", "worktree_changed", "index_changed", "staged_addition", "mode_changed"]
)
def test_compare_refuses_uncommitted_record_bytes(
    helper: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, state: str
) -> None:
    root = tmp_path / "repo"
    root.mkdir()
    repo = Repo.init(root)
    path = root / "tests/fixtures/record.json"
    path.parent.mkdir(parents=True)
    path.write_text('{"value":"original"}\n')
    relative = path.relative_to(root).as_posix()
    repo.index.add([relative])
    identity = Actor("Fixture", "fixture@example.invalid")
    commit = repo.index.commit("synthetic evidence fixture", author=identity, committer=identity)
    if state == "unchanged":
        pass
    elif state == "staged_addition":
        path = root / "tests/fixtures/addition.json"
        relative = path.relative_to(root).as_posix()
        path.write_text('{"value":"new"}\n')
        repo.index.add([relative])
    elif state == "mode_changed":
        path.chmod(0o755)
    else:
        path.write_text('{"value":"changed"}\n')
        if state == "index_changed":
            repo.index.add([relative])
            # Isolate staged bytes: HEAD and the working file agree again.
            path.write_text('{"value":"original"}\n')
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    (baseline / "COMMIT").write_text(commit.hexsha)

    class NoValidation:
        def __init__(self, _baseline: Path, _sink: Any) -> None:
            self.counts: Counter[str] = Counter()
            self.files: list[dict[str, Any]] = []

        def consume(self, name: str, _handle: Any, _size: int) -> None:
            self.files.append({"path": name})

    monkeypatch.setattr(helper, "ROOT", root)
    monkeypatch.setattr(helper, "Comparison", NoValidation)
    monkeypatch.setattr(helper, "source_inventory", lambda _root: {})
    monkeypatch.setattr(helper, "tracked_candidates", lambda: ([relative], []))
    if state == "unchanged":
        assert (
            helper.compare(baseline, tmp_path / "result.json")["status"]
            == "NO_NEW_FAILURES_OR_SERIALIZATION_CHANGES"
        )
    else:
        with pytest.raises(ValueError, match=r"committed|HEAD|index|mode"):
            helper.compare(baseline, tmp_path / "result.json")


def test_staged_fixture_deletion_is_discovered_and_refused(
    helper: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    root = tmp_path / "repo"
    root.mkdir()
    repo = Repo.init(root)
    path = root / "tests/fixtures/record.json"
    path.parent.mkdir(parents=True)
    original = '{"value":"original"}\n'
    path.write_text(original)
    relative = path.relative_to(root).as_posix()
    repo.index.add([relative])
    identity = Actor("Fixture", "fixture@example.invalid")
    commit = repo.index.commit("synthetic evidence fixture", author=identity, committer=identity)
    repo.index.remove([relative], working_tree=False)
    assert path.read_text() == original
    assert (relative, 0) not in repo.index.entries
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    (baseline / "COMMIT").write_text(commit.hexsha)

    class NoValidation:
        def __init__(self, _baseline: Path, _sink: Any) -> None:
            self.counts: Counter[str] = Counter()
            self.files: list[dict[str, Any]] = []

        def consume(self, name: str, _handle: Any, _size: int) -> None:
            self.files.append({"path": name})

    monkeypatch.setattr(helper, "ROOT", root)
    monkeypatch.setattr(helper, "Comparison", NoValidation)
    monkeypatch.setattr(helper, "source_inventory", lambda _root: {})
    # Keep discovery real: an index-only listing would omit the committed record
    # and let compare() falsely report success without examining its binding.
    candidates, _ = helper.tracked_candidates()
    with pytest.raises(ValueError, match="index"):
        helper.compare(baseline, tmp_path / "result.json")
    assert relative in candidates


@pytest.mark.parametrize("kind", ["direct", "noninstance"])
@pytest.mark.parametrize(
    "payload", [b'a\tb\n"' + b"x\n" * 100 + b'"\t2\n', b"a\tb\n" + b"\t" * 200 + b"\n"]
)
def test_csv_cumulative_record_bound(
    helper: ModuleType, monkeypatch: pytest.MonkeyPatch, kind: str, payload: bytes
) -> None:
    monkeypatch.setattr(helper, "MAX_RECORD_CHARS", 64, raising=False)
    if kind == "direct":
        with pytest.raises(ValueError, match="CSV record exceeds"):
            list(helper.table_records(io.BytesIO(payload), "fixture.tsv"))
    else:
        comparison = object.__new__(helper.Comparison)
        comparison.counts = Counter()
        comparison.files = []
        with pytest.raises(ValueError, match="CSV record exceeds"):
            comparison.consume(
                "data/derived/study_trait_profiles_wide.tsv.gz", io.BytesIO(payload), len(payload)
            )


def test_small_multiline_cell_is_preserved(
    helper: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(helper, "MAX_RECORD_CHARS", 64, raising=False)
    rows = list(helper.table_records(io.BytesIO(b'a\tb\n"x\ny"\t2\n'), "fixture.tsv"))
    assert rows == [(2, {"a": "x\ny", "b": "2"})]
