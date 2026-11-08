"""Preserve local-path spans and connect both consumers to one source (#1759)."""

from __future__ import annotations

import importlib
import importlib.util
import re
import sys
from pathlib import Path
from types import ModuleType

import pytest

import test_no_absolute_local_paths as guard
from mtp import provenance_paths as writer

ROOT = Path(__file__).resolve().parents[1]
HOME = "/" + "Users" + "/" + "example"
DOTTED_HOME = "/" + "Users" + "/" + "first.last"
LINUX_HOME = "/" + "home" + "/" + "example"
TEMP = "/" + "private" + "/" + "var" + "/" + "folders"
WINDOWS = "C:" + "\\" + "Users" + "\\"
MIXED_WINDOWS = "q:" + "\\/" + "Users" + "/\\"
ESCAPED_WINDOWS = "C:" + "\\\\" + "Users" + "\\\\"
MatchRow = tuple[tuple[int, int], str, str | None]
CANONICAL_CASES = (
    pytest.param("a " + HOME + "/data", [((2, 16), HOME, "user_home")], id="macos"),
    pytest.param(DOTTED_HOME + "/data", [((0, 17), DOTTED_HOME, "user_home")], id="dotted-name"),
    pytest.param(LINUX_HOME, [((0, 13), LINUX_HOME, "user_home")], id="linux-name-end"),
    pytest.param(
        "/" + "System/Volumes/Data" + HOME,
        [((20, 34), HOME, "user_home")],
        id="firmlink",
    ),
    pytest.param(WINDOWS + "example", [((0, 9), WINDOWS, "user_home")], id="windows"),
    pytest.param(
        MIXED_WINDOWS + "example", [((0, 11), MIXED_WINDOWS, "user_home")], id="mixed-slashes"
    ),
    pytest.param(
        ESCAPED_WINDOWS + "example",
        [((0, 11), ESCAPED_WINDOWS, "user_home")],
        id="escaped-slashes",
    ),
    pytest.param(TEMP + "/leaf", [((0, 20), TEMP, "temp_root")], id="private-temp"),
    pytest.param(
        "/" + "var" + "/folders/leaf",
        [((0, 12), "/" + "var" + "/folders", "temp_root")],
        id="temp-firmlink",
    ),
    pytest.param(
        "/" + "private" + "/tmp/leaf",
        [((0, 12), "/" + "private" + "/tmp", "temp_root")],
        id="private-tmp",
    ),
    pytest.param(
        TEMP + HOME,
        [((0, 20), TEMP, "temp_root"), ((20, 34), HOME, "user_home")],
        id="adjacent-families",
    ),
    pytest.param(
        "tag " + HOME + "; " + LINUX_HOME,
        [((4, 18), HOME, "user_home"), ((20, 33), LINUX_HOME, "user_home")],
        id="two-homes",
    ),
)


@pytest.fixture(autouse=True)
def assigned_consumer_origins() -> None:
    assert Path(guard.__file__).resolve() == ROOT / "tests/test_no_absolute_local_paths.py"
    assert Path(writer.__file__).resolve() == ROOT / "src/mtp/provenance_paths.py"


def shared_module() -> ModuleType:
    # Missing implementation must fail an assertion during the test call, not collection.
    assert importlib.util.find_spec("mtp.local_path_patterns") is not None, (
        "The consumers need a shared local-path pattern module"
    )
    module = importlib.import_module("mtp.local_path_patterns")
    assert module.__file__ is not None
    assert Path(module.__file__).resolve() == ROOT / "src/mtp/local_path_patterns.py"
    assert callable(getattr(module, "compile_absolute_local_paths", None))
    return module


def matches(pattern: re.Pattern[str], text: str) -> list[MatchRow]:
    return [(match.span(), match.group(0), match.lastgroup) for match in pattern.finditer(text)]


@pytest.mark.parametrize(("text", "expected"), CANONICAL_CASES)
def test_existing_consumers_preserve_exact_spans_and_family_order(
    text: str, expected: list[MatchRow]
) -> None:
    assert matches(guard.ABSOLUTE_LOCAL_PATH_PATTERN, text) == expected
    for family, pattern in (
        ("user_home", guard.USER_HOME_PATTERN),
        ("temp_root", guard.MACOS_TEMP_ROOT_PATTERN),
    ):
        assert matches(pattern, text) == [
            (span, value, None) for span, value, kind in expected if kind == family
        ]
    assert matches(writer._LOCAL_PATH, text) == [(span, value, None) for span, value, _ in expected]
    with pytest.raises(ValueError, match="Provenance contains a local path or local prefix"):
        writer.assert_portable_text(text)


@pytest.mark.parametrize(("text", "expected"), CANONICAL_CASES)
def test_shared_compiler_preserves_exact_spans_and_family_order(
    text: str, expected: list[MatchRow]
) -> None:
    shared = shared_module()
    classified = shared.compile_absolute_local_paths(flags=0, classified=True)
    plain = shared.compile_absolute_local_paths(flags=re.IGNORECASE, classified=False)
    assert classified.groupindex == {"user_home": 1, "temp_root": 2}
    assert plain.groupindex == {}
    assert not classified.flags & re.IGNORECASE
    assert plain.flags & re.IGNORECASE
    assert matches(classified, text) == expected
    assert matches(plain, text) == [(span, value, None) for span, value, _ in expected]
    for family, source in (
        ("user_home", shared.USER_HOME_SOURCE),
        ("temp_root", shared.TEMP_ROOT_SOURCE),
    ):
        assert matches(re.compile(source), text) == [
            (span, value, None) for span, value, kind in expected if kind == family
        ]


@pytest.mark.parametrize(
    "text",
    [
        "/" + "users" + "/example",
        "/" + "HOME" + "/example",
        "C:" + "\\" + "uSeRs" + "\\" + "example",
        "/" + "PRIVATE" + "/VAR/FOLDERS/leaf",
        "Users" + "/example",
        "home" + "/example",
        "var" + "/folders/leaf",
        "private" + "/tmp/leaf",
    ],
)
def test_provenance_keeps_its_additional_case_and_slashless_checks(text: str) -> None:
    assert matches(guard.ABSOLUTE_LOCAL_PATH_PATTERN, text) == []
    with pytest.raises(ValueError, match="Provenance contains a local path or local prefix"):
        writer.assert_portable_text(text)


def test_shared_posix_matches_keep_their_slashless_provenance_rejection() -> None:
    """Check the two-way relation across pattern roots and name-class characters (#1835, #1871)."""
    shared = shared_module()
    absolute = shared.compile_absolute_local_paths(flags=re.IGNORECASE, classified=True)
    pattern_sources = (
        shared.USER_HOME_SOURCE,
        shared.TEMP_ROOT_SOURCE,
        writer._SLASHLESS_LOCAL_PATH.pattern,
    )
    char_classes = [
        re.compile(cls)
        for src in (shared.USER_HOME_SOURCE, writer._SLASHLESS_LOCAL_PATH.pattern)
        for cls in re.findall(r"/(\[[^\]]+\])\+", src)
    ]
    assert len(char_classes) == 2
    name_chars = [
        ch
        for code in range(33, 127)
        if (ch := chr(code)) not in "/\\" and any(cls.fullmatch(ch) for cls in char_classes)
    ]
    extracted_roots = [
        word
        for src in pattern_sources
        for word in re.findall(r"[A-Za-z0-9_.-]+", re.sub(r"\[[^\]]+\]|\\[A-Za-z]", " ", src))
    ]
    names = tuple(
        dict.fromkeys(
            [
                "example",
                f"a{''.join(name_chars)}z",
                *(f"a{ch}z" for ch in name_chars if not ch.isalnum()),
            ]
        )
    )
    home_roots = tuple(
        dict.fromkeys(
            [
                "Users",
                "home",
                "people",
                "accounts",
                "Volumes",
                "mnt",
                "media",
                "srv",
                *extracted_roots,
                *names,
            ]
        )
    )
    temp_roots = (
        ("var", "folders"),
        ("private", "var", "folders"),
        ("private", "tmp"),
        ("var", "tmp"),
        ("tmp",),
    )
    candidates = ["/" + root + "/" + name for root in home_roots for name in names]
    candidates.extend("/" + "/".join(parts) for parts in temp_roots)
    checked = 0
    converse_checked = 0
    families: set[str | None] = set()
    converse_families: set[str | None] = set()
    for candidate in candidates:
        for spelling in (candidate, candidate.upper()):
            for suffix in ("", "/", "/leaf"):
                text = spelling + suffix
                for match in absolute.finditer(text):
                    if not match.group(0).startswith("/"):
                        continue  # A drive-qualified Windows path has no leading slash to remove.
                    # Keep the original suffix: temp-root matches rely on its next slash.
                    slashless = "tag " + match.group(0)[1:] + text[match.end() :]
                    found = writer._SLASHLESS_LOCAL_PATH.search(slashless)
                    assert found is not None, (
                        f"Shared POSIX match lost after slash removal: {text!r}"
                    )
                    assert found.span() == (4, 4 + len(match.group(0)) - 1)
                    with pytest.raises(
                        ValueError, match="Provenance contains a local path or local prefix"
                    ):
                        writer.assert_portable_text(slashless)
                    checked += 1
                    families.add(match.lastgroup)

                slashless_candidate = "tag " + spelling[1:] + suffix
                for found in writer._SLASHLESS_LOCAL_PATH.finditer(slashless_candidate):
                    slashed = (
                        slashless_candidate[: found.start()]
                        + "/"
                        + slashless_candidate[found.start() :]
                    )
                    found_abs = absolute.search(slashed)
                    assert found_abs is not None, (
                        f"Slashless match lost after prepending '/': {slashless_candidate!r}"
                    )
                    assert found_abs.span() == (found.start(), found.end() + 1)
                    assert found_abs.group(0) == "/" + found.group(0)
                    converse_checked += 1
                    converse_families.add(found_abs.lastgroup)
    # Current home spellings and temp cases in both directions.
    assert checked >= 60
    assert converse_checked >= 60
    assert families == {"user_home", "temp_root"}
    assert converse_families == {"user_home", "temp_root"}


@pytest.mark.parametrize(
    "text",
    [
        "data/table.tsv",
        "someUsers/example",
        "x_Users/example",
        "var/foldersish",
        "/" + "tmp" + "/leaf",
        "var/folders",
    ],
)
def test_safe_text_controls_stay_outside_both_patterns(text: str) -> None:
    # These are text-only rules; structural private-root checks remain separate.
    assert matches(guard.ABSOLUTE_LOCAL_PATH_PATTERN, text) == []
    writer.assert_portable_text(text)


def test_consumer_grouping_and_case_flags_remain_distinct() -> None:
    assert guard.ABSOLUTE_LOCAL_PATH_PATTERN.groupindex == {"user_home": 1, "temp_root": 2}
    assert writer._LOCAL_PATH.groupindex == {}
    for pattern in (
        guard.ABSOLUTE_LOCAL_PATH_PATTERN,
        guard.USER_HOME_PATTERN,
        guard.MACOS_TEMP_ROOT_PATTERN,
    ):
        assert not pattern.flags & re.IGNORECASE
    assert writer._LOCAL_PATH.flags & re.IGNORECASE
    assert writer._SLASHLESS_LOCAL_PATH.flags & re.IGNORECASE


def test_guard_scan_keeps_adjacent_hits_on_the_same_line() -> None:
    hits = guard._scan_records((("example.json", "fixture", "\n" + TEMP + HOME),))
    assert set(hits) == {"example.json"}
    assert [(hit.line_no, hit.family, hit.matched) for hit in hits["example.json"]] == [
        (2, "temp_root", TEMP),
        (2, "user_home", HOME),
    ]


def sentinel_compiler(monkeypatch: pytest.MonkeyPatch) -> list[tuple[int, bool]]:
    shared = shared_module()
    calls: list[tuple[int, bool]] = []

    def compile_sentinel(*, flags: int, classified: bool) -> re.Pattern[str]:
        calls.append((flags, classified))
        source = (
            "(?P<user_home>shared_home_marker)|(?P<temp_root>shared_temp_marker)"
            if classified
            else "shared_home_marker|shared_temp_marker"
        )
        return re.compile(source, flags)

    monkeypatch.setattr(shared, "compile_absolute_local_paths", compile_sentinel)
    return calls


def isolated_consumer(relative: str, name: str, monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, module)
    spec.loader.exec_module(module)
    return module


def test_guard_uses_the_shared_compiler_result_for_scanning(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = sentinel_compiler(monkeypatch)
    consumer = isolated_consumer(
        "tests/test_no_absolute_local_paths.py", "_guard_pattern_connection", monkeypatch
    )
    assert calls == [(0, True)]
    hits = consumer._scan_records(
        (("example.json", "fixture", "shared_temp_marker shared_home_marker"),)
    )
    assert set(hits) == {"example.json"}
    assert [(hit.family, hit.matched) for hit in hits["example.json"]] == [
        ("temp_root", "shared_temp_marker"),
        ("user_home", "shared_home_marker"),
    ]
    assert consumer._scan_records((("example.json", "fixture", "SHARED_HOME_MARKER"),)) == {}


@pytest.mark.parametrize(
    ("source_name", "pattern_name", "marker"),
    [
        ("USER_HOME_SOURCE", "USER_HOME_PATTERN", "shared_home_source_marker"),
        ("TEMP_ROOT_SOURCE", "MACOS_TEMP_ROOT_PATTERN", "shared_temp_source_marker"),
    ],
)
def test_guard_family_patterns_use_the_shared_source_value(
    monkeypatch: pytest.MonkeyPatch, source_name: str, pattern_name: str, marker: str
) -> None:
    shared = shared_module()
    monkeypatch.setattr(shared, source_name, marker)
    consumer = isolated_consumer(
        "tests/test_no_absolute_local_paths.py", "_guard_family_source_connection", monkeypatch
    )
    pattern = getattr(consumer, pattern_name)
    assert pattern.pattern == marker
    assert matches(pattern, "prefix " + marker) == [((7, 7 + len(marker)), marker, None)]
    assert matches(pattern, marker.upper()) == []


def test_provenance_uses_the_shared_compiler_result_for_refusal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = sentinel_compiler(monkeypatch)
    consumer = isolated_consumer(
        "src/mtp/provenance_paths.py", "_provenance_pattern_connection", monkeypatch
    )
    assert calls == [(re.IGNORECASE, False)]
    for text in ("shared_home_marker", "SHARED_TEMP_MARKER"):
        with pytest.raises(ValueError, match="Provenance contains a local path or local prefix"):
            consumer.assert_portable_text(text)


def test_shared_pattern_module_and_tests_do_not_match_the_guard() -> None:
    shared_module()
    names = ("src/mtp/local_path_patterns.py", "tests/test_local_path_patterns.py")
    records = tuple((name, "fixture", (ROOT / name).read_text(encoding="utf-8")) for name in names)
    assert guard._scan_records(records) == {}


def test_shared_pattern_module_is_a_required_guard_input() -> None:
    assert "src/mtp/local_path_patterns.py" in guard.REQUIRED_SCANNED
