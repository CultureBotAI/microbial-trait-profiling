"""Keep the historical comparison fixture reproducible without executing it (#1610)."""

from __future__ import annotations

import importlib.util
import io
import json
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import cast

import pytest

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts/adhoc/extract_comparison_specs_f0cb5b7.py"
FIXTURE = ROOT / "tests/fixtures/comparison_specs/analysis_specs_f0cb5b7.json"


@pytest.fixture(scope="module")
def extractor() -> ModuleType:
    spec = importlib.util.spec_from_file_location("comparison_fixture_extractor", HELPER)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.__file__ is not None
    assert Path(module.__file__).resolve() == HELPER.resolve()
    return module


def fixture_specs() -> list[dict[str, str]]:
    return cast(list[dict[str, str]], json.loads(FIXTURE.read_bytes()))


def literal_source(specs: list[dict[str, str]]) -> str:
    entries = [
        f"{item['id']!r}: AnalysisSpec("
        + ", ".join(f"{key}={value!r}" for key, value in item.items())
        + ")"
        for item in specs
    ]
    return "STANDARD_ANALYSIS_SPECS: dict[str, AnalysisSpec] = {" + ",\n".join(entries) + "}\n"


class TrackedStream(io.BytesIO):
    def __init__(self, content: bytes) -> None:
        super().__init__(content)
        self.read_sizes: list[int | None] = []

    def read(self, size: int | None = -1) -> bytes:
        self.read_sizes.append(size)
        return super().read(size)


class HistoricalTree:
    def __init__(self, content: bytes) -> None:
        self.stream = TrackedStream(content)
        self.paths: list[str] = []

    def __truediv__(self, path: str) -> SimpleNamespace:
        self.paths.append(path)
        return SimpleNamespace(data_stream=self.stream)


class HistoricalRepository:
    def __init__(self, content: bytes) -> None:
        self.tree = HistoricalTree(content)
        self.revisions: list[str] = []
        self.closed = False

    def commit(self, revision: str) -> SimpleNamespace:
        self.revisions.append(revision)
        return SimpleNamespace(tree=self.tree)

    def close(self) -> None:
        self.closed = True


def install_repository(
    extractor: ModuleType, monkeypatch: pytest.MonkeyPatch, root: Path, content: bytes
) -> HistoricalRepository:
    repository = HistoricalRepository(content)

    def open_repository(path: Path) -> HistoricalRepository:
        assert path == root
        return repository

    monkeypatch.setattr(extractor, "Repo", open_repository)
    return repository


def test_extracts_only_declared_literals_in_source_order_without_execution(
    extractor: ModuleType,
) -> None:
    expected = [dict(reversed(item.items())) for item in reversed(fixture_specs())]
    source = (
        "raise RuntimeError('historical code must never execute')\n"
        + literal_source(expected)
        + "unrelated = AnalysisSpec(id='not_in_the_mapping')\n"
    )
    actual = extractor.extract_specs(source)
    assert actual == expected
    assert [list(item) for item in actual] == [list(item) for item in expected]


def test_pinned_git_extraction_reproduces_fixture_bytes_and_refuses_overwrite(
    extractor: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    original = FIXTURE.read_bytes()
    repository = install_repository(
        extractor, monkeypatch, tmp_path, literal_source(fixture_specs()).encode("utf-8")
    )
    output = tmp_path / "reproduced.json"
    extractor.write_fixture(tmp_path, output)
    assert output.read_bytes() == original
    assert FIXTURE.read_bytes() == original
    assert repository.revisions == [extractor.SOURCE_COMMIT]
    assert repository.tree.paths == ["src/mtp/workflow/downstream_contracts.py"]
    assert repository.tree.stream.read_sizes == [extractor.MAX_SOURCE_BYTES + 1]
    assert repository.closed
    with pytest.raises(FileExistsError):
        extractor.write_fixture(tmp_path, output)
    assert output.read_bytes() == original


def test_git_source_read_is_bounded_and_closed_on_refusal(
    extractor: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    repository = install_repository(
        extractor, monkeypatch, tmp_path, b" " * (extractor.MAX_SOURCE_BYTES + 1)
    )
    with pytest.raises(ValueError, match=r"source.*bound"):
        extractor.read_pinned_source(tmp_path)
    assert repository.tree.stream.read_sizes == [extractor.MAX_SOURCE_BYTES + 1]
    assert repository.closed


@pytest.mark.parametrize("count", [6, 8])
def test_refuses_a_partial_or_extra_declaration_list(extractor: ModuleType, count: int) -> None:
    specs = fixture_specs()
    if count == 6:
        specs.pop()
    else:
        specs.append(dict(specs[0], id="unexpected_eighth_spec"))
    with pytest.raises(ValueError, match="seven"):
        extractor.extract_specs(literal_source(specs))


@pytest.mark.parametrize(
    "replacement",
    ["kernel=compute_kernel()", "kernel=17", "**extra_fields"],
    ids=["computed-value", "non-string-value", "keyword-expansion"],
)
def test_refuses_nonliteral_fields(extractor: ModuleType, replacement: str) -> None:
    source = literal_source(fixture_specs()).replace(
        "kernel='two_sided_exact_binomial_sign_test'", replacement, 1
    )
    with pytest.raises(ValueError, match="literal"):
        extractor.extract_specs(source)


def test_refuses_missing_fields_and_mismatched_mapping_ids(extractor: ModuleType) -> None:
    specs = fixture_specs()
    del specs[0]["weighting"]
    with pytest.raises(ValueError, match="fields"):
        extractor.extract_specs(literal_source(specs))
    source = literal_source(fixture_specs()).replace(
        "'noel_trait_direction_screen_sign_test_v1':", "'wrong_mapping_key':", 1
    )
    with pytest.raises(ValueError, match="mapping key"):
        extractor.extract_specs(source)


def test_real_pinned_blob_reproduces_the_fixture_when_the_commit_is_present(
    extractor: ModuleType, tmp_path: Path
) -> None:
    """Read the real historical blob, not a stand-in (#1721); skip on a shallow clone."""
    from git import Repo
    from git.exc import BadName, InvalidGitRepositoryError, NoSuchPathError

    try:
        commit = Repo(ROOT).commit(extractor.SOURCE_COMMIT)
    except (BadName, InvalidGitRepositoryError, NoSuchPathError, ValueError):
        pytest.skip("the pinned commit is not in this checkout (shallow clone); see the README")
    # The short SHA resolves to the commit the README names in full (#1722).
    assert commit.hexsha.startswith(extractor.SOURCE_COMMIT)
    output = tmp_path / "from_real_blob.json"
    extractor.write_fixture(ROOT, output)
    assert output.read_bytes() == FIXTURE.read_bytes()
