"""#1294: strict row fast-path parity, diagnostics, and bounded verifier proofs."""

import csv
import importlib.util
import io
import sys
from collections.abc import Sequence
from pathlib import Path
from types import ModuleType

import pytest

from mtp import count_input
from mtp.count_input import (
    COMPLETE_COUNT_MODE,
    LEGACY_SPARSE_ZERO_MODE,
    CountInputMode,
    CountReadSummary,
    CountTableReader,
    CountTableRow,
    parse_nonnegative_integer,
)

KEYS = ("study_id", "experiment_id", "inc_dec")
TRAITS = ("trait_a", "trait_b", "trait_c")
INVALID = ("", "   ", "-1", "+1", "1.5", "NaN", "٢", "²", "\uff11\uff12", "0x10", "1_000", "1e2")


def _buffer(header: Sequence[str], rows: Sequence[Sequence[str]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return stream.getvalue().encode()


def _reader(
    path: Path,
    rows: Sequence[Sequence[str]],
    *,
    mode: CountInputMode = COMPLETE_COUNT_MODE,
) -> CountTableReader:
    return CountTableReader(
        path,
        leading_columns=KEYS,
        key_columns=KEYS,
        mode=mode,
        raw_bytes=_buffer((*KEYS, *TRAITS), rows),
    )


def _forbid_parser(*_args: object, **_kwargs: object) -> int:
    raise AssertionError("per-cell parser should not be called on valid COMPLETE_COUNT_MODE row")


def _spy_parser(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, int, str, CountInputMode]]:
    calls: list[tuple[str, int, str, CountInputMode]] = []

    def spy(
        value: str,
        *,
        path: Path,
        line_number: int,
        column_name: str,
        mode: CountInputMode = COMPLETE_COUNT_MODE,
        row_key: Sequence[str] | None = None,
    ) -> int:
        calls.append((value, line_number, column_name, mode))
        return parse_nonnegative_integer(
            value,
            path=path,
            line_number=line_number,
            column_name=column_name,
            mode=mode,
            row_key=row_key,
        )

    monkeypatch.setattr(count_input, "parse_nonnegative_integer", spy)
    return calls


def test_complete_valid_rows_bypass_parser(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    raw = [
        (" s1 ", " e1 ", " inc ", "0", " 42 ", "000123"),
        ("s2", "e2", "dec", "\u00a07\u00a0", "98765432101234567890", "1"),
    ]
    reader = _reader(tmp_path / "profiles.tsv", raw)
    monkeypatch.setattr(count_input, "parse_nonnegative_integer", _forbid_parser)
    rows = list(reader)
    assert [row.counts for row in rows] == [(0, 42, 123), (7, 98765432101234567890, 1)]
    assert [row.values for row in rows] == raw
    assert [row.key for row in rows] == [("s1", "e1", "inc"), ("s2", "e2", "dec")]
    assert reader.summary == CountReadSummary(2, 3, 0, COMPLETE_COUNT_MODE)


def test_complete_rows_match_per_cell_reference(tmp_path: Path) -> None:
    raw = [("s1", "e1", "inc", "0", " 42 ", "001"), ("s2", "e2", "dec", "10000", "7", "99")]
    path = tmp_path / "profiles.tsv"
    expected = [
        CountTableRow(
            line,
            tuple(value.strip() for value in row[:3]),
            tuple(row),
            tuple(
                parse_nonnegative_integer(
                    value,
                    path=path,
                    line_number=line,
                    column_name=column,
                    row_key=row[:3],
                )
                for column, value in zip(TRAITS, row[3:], strict=True)
            ),
        )
        for line, row in enumerate(raw, start=2)
    ]
    reader = _reader(path, raw)
    assert list(reader) == expected
    assert reader.summary == CountReadSummary(len(expected), len(TRAITS), 0, COMPLETE_COUNT_MODE)
    assert list(reader) == expected
    assert reader.summary.rows == len(expected)


@pytest.mark.parametrize("bad", INVALID)
@pytest.mark.parametrize("position", [0, 1, 2])
def test_complete_invalid_cells_preserve_exact_parser_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    bad: str,
    position: int,
) -> None:
    path = tmp_path / "profiles.tsv"
    values = ["11", "22", "33"]
    values[position] = bad
    key = ("s1", "e1", "inc")
    with pytest.raises(
        ValueError, match=r"Missing count|finite nonnegative integer count"
    ) as expected:
        parse_nonnegative_integer(
            bad,
            path=path,
            line_number=2,
            column_name=TRAITS[position],
            row_key=key,
        )
    calls = _spy_parser(monkeypatch)
    reader = _reader(path, [(*key, *values)])
    # Compare strings explicitly so mutations must fail with native AssertionError.
    try:
        list(reader)
    except ValueError as actual:
        actual_message = str(actual)
    else:
        raise AssertionError("Malformed complete rows must retain the strict per-column error")
    assert actual_message == str(expected.value)
    assert calls == [
        (value, 2, column, COMPLETE_COUNT_MODE)
        for column, value in zip(TRAITS[: position + 1], values[: position + 1], strict=True)
    ]
    assert reader.summary == CountReadSummary(0, 3, 0, COMPLETE_COUNT_MODE)


def test_legacy_always_uses_per_cell_parser(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    raw = [("s1", "e1", "inc", "1", "2", "3"), ("s2", "e2", "dec", "", "  ", " 42 ")]
    calls = _spy_parser(monkeypatch)
    reader = _reader(tmp_path / "profiles.tsv", raw, mode=LEGACY_SPARSE_ZERO_MODE)
    for _ in range(2):
        calls.clear()
        assert [row.counts for row in reader] == [(1, 2, 3), (0, 0, 42)]
        assert calls == [
            (value, line, column, LEGACY_SPARSE_ZERO_MODE)
            for line, row in enumerate(raw, start=2)
            for column, value in zip(TRAITS, row[3:], strict=True)
        ]
        assert reader.summary == CountReadSummary(2, 3, 2, LEGACY_SPARSE_ZERO_MODE)


@pytest.mark.parametrize("bad", INVALID[2:])
def test_legacy_invalid_cells_preserve_exact_parser_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    bad: str,
) -> None:
    path = tmp_path / "profiles.tsv"
    key = ("s1", "e1", "inc")
    with pytest.raises(
        ValueError, match=r"Missing count|finite nonnegative integer count"
    ) as expected:
        parse_nonnegative_integer(
            bad,
            path=path,
            line_number=2,
            column_name=TRAITS[1],
            mode=LEGACY_SPARSE_ZERO_MODE,
            row_key=key,
        )
    calls = _spy_parser(monkeypatch)
    reader = _reader(path, [(*key, "", bad, "3")], mode=LEGACY_SPARSE_ZERO_MODE)
    with pytest.raises(
        ValueError, match=r"Missing count|finite nonnegative integer count"
    ) as actual:
        list(reader)
    assert str(actual.value) == str(expected.value)
    assert calls == [
        ("", 2, TRAITS[0], LEGACY_SPARSE_ZERO_MODE),
        (bad, 2, TRAITS[1], LEGACY_SPARSE_ZERO_MODE),
    ]
    assert reader.summary == CountReadSummary(0, 3, 1, LEGACY_SPARSE_ZERO_MODE)


def test_declared_count_order_preserves_first_error(tmp_path: Path) -> None:
    path = tmp_path / "profiles.tsv"
    header = (*KEYS, "metadata", *TRAITS)
    raw = [
        ("s1", "e1", "inc", "note", "1", "2", "3"),
        ("s2", "e2", "dec", "quoted\nmetadata", "bad-a", "bad-b", "bad-c"),
    ]
    reader = CountTableReader(
        path,
        leading_columns=KEYS,
        key_columns=KEYS,
        count_columns=("trait_c", "trait_a"),
        raw_bytes=_buffer(header, raw),
    )
    iterator = iter(reader)
    assert next(iterator).counts == (3, 1)
    with pytest.raises(
        ValueError, match=r"Missing count|finite nonnegative integer count"
    ) as expected:
        parse_nonnegative_integer(
            "bad-c",
            path=path,
            line_number=4,
            column_name="trait_c",
            row_key=("s2", "e2", "dec"),
        )
    with pytest.raises(
        ValueError, match=r"Missing count|finite nonnegative integer count"
    ) as actual:
        next(iterator)
    assert str(actual.value) == str(expected.value)
    assert reader.summary == CountReadSummary(1, 2, 0, COMPLETE_COUNT_MODE)


@pytest.mark.parametrize("explicit", [False, True])
def test_zero_count_columns_keep_constructor_refusal(tmp_path: Path, explicit: bool) -> None:
    # This pre-existing constructor contract is outside #1294's protected edit block.
    with pytest.raises(ValueError, match="has no count columns"):
        CountTableReader(
            tmp_path / "empty.tsv",
            leading_columns=KEYS,
            key_columns=KEYS,
            count_columns=() if explicit else None,
            raw_bytes=_buffer(KEYS, []),
        )


def test_wide_benchmark_fixture_has_parity_without_per_cell_calls(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Deterministic workload proof; the separate timing report has no flaky speed threshold.
    traits = tuple(f"trait_{i}" for i in range(317))
    raw = [
        (f"s{index}", "e1", "inc", *(str((index + i) % 101) for i in range(len(traits))))
        for index in range(64)
    ]
    path = tmp_path / "wide.tsv"
    expected = [
        tuple(
            parse_nonnegative_integer(
                value,
                path=path,
                line_number=line,
                column_name=column,
                row_key=row[:3],
            )
            for column, value in zip(traits, row[3:], strict=True)
        )
        for line, row in enumerate(raw, start=2)
    ]
    reader = CountTableReader(
        path, leading_columns=KEYS, key_columns=KEYS, raw_bytes=_buffer((*KEYS, *traits), raw)
    )
    monkeypatch.setattr(count_input, "parse_nonnegative_integer", _forbid_parser)
    assert [row.counts for row in reader] == expected
    assert reader.summary == CountReadSummary(64, 317, 0, COMPLETE_COUNT_MODE)


def _verifier() -> ModuleType:
    name = "count_input_1294_verifier_under_test"
    if name in sys.modules:
        return sys.modules[name]
    path = (
        Path(__file__).resolve().parents[1] / "scripts/adhoc/verify_1294_mutants_codex_20261111.py"
    )
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the verifier")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("placement", ["root", "nested", "inbound", "outbound"])
def test_verifier_rejects_unsafe_scratch_before_any_mkdir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    placement: str,
) -> None:
    verifier = _verifier()
    root = tmp_path / "repo"
    root.mkdir()
    external = tmp_path / "external"
    external.mkdir()
    if placement == "root":
        scratch = root
    elif placement == "nested":
        scratch = root / "not_created" / "scratch"
    elif placement == "inbound":
        link = external / "into_repo"
        link.symlink_to(root, target_is_directory=True)
        scratch = link / "scratch"
    else:
        link = root / "out_of_repo"
        link.symlink_to(external, target_is_directory=True)
        scratch = link / "scratch"
    calls: list[str] = []

    def forbid(*_args: object, **_kwargs: object) -> None:
        calls.append("filesystem mutation")
        raise AssertionError("Containment must be checked before mkdir or mkdtemp")

    monkeypatch.setattr(Path, "mkdir", forbid)
    monkeypatch.setattr(verifier.tempfile, "mkdtemp", forbid)
    with pytest.raises(ValueError, match="outside the repository"):
        verifier.verify_mutants(root, scratch, canary_only=True)
    assert calls == []


@pytest.mark.parametrize("exit_code", [0, 2, 3, 4, 5, 124])
def test_verifier_only_accepts_exact_test_failure_exit(exit_code: int) -> None:
    call = {
        "outcome": "failed",
        "native_assertion": True,
        "test_origin_ok": True,
        "assertion_in_test": True,
    }
    verifier = _verifier()
    assert verifier.classify_pytest_result(1, [call], []) is True
    assert verifier.classify_pytest_result(exit_code, [call], []) is False


@pytest.mark.parametrize("missing", ["native_assertion", "test_origin_ok", "assertion_in_test"])
def test_verifier_requires_native_bound_call_assertion(missing: str) -> None:
    call = {
        "outcome": "failed",
        "native_assertion": True,
        "test_origin_ok": True,
        "assertion_in_test": True,
    }
    call[missing] = False
    assert _verifier().classify_pytest_result(1, [call], []) is False


@pytest.mark.parametrize("phase", ["setup", "teardown", "collection"])
def test_verifier_rejects_non_call_failures(phase: str) -> None:
    call = {
        "outcome": "failed",
        "native_assertion": True,
        "test_origin_ok": True,
        "assertion_in_test": True,
    }
    assert _verifier().classify_pytest_result(1, [call], [{"phase": phase}]) is False


@pytest.mark.parametrize("case", ["internal_error", "setup", "assertion", "value_error"])
def test_verifier_real_pytest_failure_classification(tmp_path: Path, case: str) -> None:
    verifier = _verifier()
    scratch = tmp_path / "external_probe"
    scratch.mkdir()
    tests = scratch / "test_probe.py"
    if case == "internal_error":
        (scratch / "conftest.py").write_text(
            "def pytest_collection_modifyitems(items):\n"
            "    raise RuntimeError('deliberate internal error')\n",
            encoding="utf-8",
        )
        tests.write_text("def test_probe():\n    pass\n", encoding="utf-8")
    elif case == "setup":
        tests.write_text(
            "import pytest\n@pytest.fixture\ndef broken():\n"
            "    raise AssertionError('setup is not a call failure')\n"
            "def test_probe(broken):\n    pass\n",
            encoding="utf-8",
        )
    else:
        exception = "AssertionError" if case == "assertion" else "ValueError"
        tests.write_text(
            f"def test_probe():\n    raise {exception}('deliberate call failure')\n",
            encoding="utf-8",
        )
    result = verifier.run_targeted_pytest(
        scratch, "test_probe.py::test_probe", tmp_path / "proof.xml"
    )
    assert result["exit_code"] == (3 if case == "internal_error" else 1)
    assert result["canonical_unchanged"] is True
    assert isinstance(result["detail"], str)
    assert result["detail"].strip()
    assert result["killed"] is (case == "assertion")
    assert result["status"] == ("KILLED" if case == "assertion" else "SURVIVED_OR_INVALID")
    assert (tmp_path / "proof.receipt.json").stat().st_size > 0
    if case == "assertion":
        assert result["call_reports"][0]["native_assertion"] is True
    elif case == "setup":
        assert result["call_reports"] == []
        assert result["non_call_failures"]


def test_verifier_refuses_inbound_evidence_file_symlink_before_mkdir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    verifier = _verifier()
    canonical = tmp_path / "synthetic_canonical"
    canonical.mkdir()
    for relative in verifier.INPUTS:
        target = canonical / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(b"bounded fixture\n")
    external = tmp_path / "external"
    external.mkdir()
    evidence = external / "proof.xml"
    evidence.symlink_to(canonical / "must_not_create" / "proof.xml")
    calls: list[str] = []

    def forbid(*_args: object, **_kwargs: object) -> None:
        calls.append("mkdir")
        raise AssertionError("An evidence-file symlink must not redirect directory creation")

    monkeypatch.setattr(Path, "mkdir", forbid)
    with pytest.raises(ValueError, match="Evidence file must not be a symlink"):
        verifier.run_targeted_pytest(
            external,
            "test_probe.py::test_probe",
            evidence,
            repo_root=canonical,
        )
    assert calls == []
