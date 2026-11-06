"""Guard tracked text files against machine-local user home and temporary paths (#1622)."""

from __future__ import annotations

import functools
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MAX_SCAN_BYTES = 5_000_000
BINARY_SUFFIXES: frozenset[str] = frozenset(
    {
        ".gif",
        ".gz",
        ".ico",
        ".jpeg",
        ".jpg",
        ".pdf",
        ".png",
        ".pptx",
        ".tar",
        ".woff",
        ".woff2",
        ".zip",
    }
)

# Regeneration follow-up issue for the four hash-bound provenance records.
PROVENANCE_REGENERATION_ISSUE = "#1664"


@dataclass(frozen=True)
class AllowlistEntry:
    """Explicit exemption for a hash-bound provenance record awaiting regeneration."""

    reason: str
    follow_up_issue: str


EXPECTED_ALLOWLIST_PATHS: frozenset[str] = frozenset(
    {
        "data/derived/integrated_evidence/integrated_evidence_assessments.tsv",
        "data/derived/integrated_evidence/integrated_evidence_provenance.json",
        "docs/reviews/2026-10-19-ora-provenance-refresh-audit.json",
        "docs/reviews/2026-10-19-ora-provenance-refresh-inputs.json",
    }
)

ALLOWLIST: dict[str, AllowlistEntry] = {
    "data/derived/integrated_evidence/integrated_evidence_assessments.tsv": AllowlistEntry(
        reason=(
            "Hash-bound integrated-evidence assessment sidecar recording resolved source paths "
            "from mtp.integrated_evidence._bound_assessment; coordinated regeneration tracked "
            "in #1664."
        ),
        follow_up_issue=PROVENANCE_REGENERATION_ISSUE,
    ),
    "data/derived/integrated_evidence/integrated_evidence_provenance.json": AllowlistEntry(
        reason=(
            "Hash-bound integrated-evidence provenance record storing resolved native_sources "
            "and assessment_status source paths from mtp.integrated_evidence.main; coordinated "
            "regeneration tracked in #1664."
        ),
        follow_up_issue=PROVENANCE_REGENERATION_ISSUE,
    ),
    "docs/reviews/2026-10-19-ora-provenance-refresh-audit.json": AllowlistEntry(
        reason=(
            "Hash-bound ORA provenance refresh audit output produced by "
            "scripts/adhoc/issue720_audit_ora_refresh.py from resolved manifest/input paths; "
            "coordinated regeneration tracked in #1664."
        ),
        follow_up_issue=PROVENANCE_REGENERATION_ISSUE,
    ),
    "docs/reviews/2026-10-19-ora-provenance-refresh-inputs.json": AllowlistEntry(
        reason=(
            "Hash-bound ORA provenance refresh input manifest bound by manifest_sha256 in "
            "docs/reviews/2026-10-19-ora-provenance-refresh-audit.json; coordinated "
            "regeneration tracked in #1664."
        ),
        follow_up_issue=PROVENANCE_REGENERATION_ISSUE,
    ),
}

EXPECTED_SYNTHETIC_FIXTURE_PATHS: frozenset[str] = frozenset(
    {
        "tests/test_output_roots.py",
        "tests/test_three_test_stage.py",
    }
)

SYNTHETIC_FIXTURE_FILES: dict[str, str] = {
    "tests/test_output_roots.py": "synthetic temp-root paths used to test rejection",
    "tests/test_three_test_stage.py": "synthetic temp-root paths used to test rejection",
}

USER_HOME_PATTERN: re.Pattern[str] = re.compile(
    r"(?:/Users/[A-Za-z0-9._-]+/|/home/[A-Za-z0-9._-]+/|[A-Za-z]:\\Users\\)"
)
MACOS_TEMP_ROOT_PATTERN: re.Pattern[str] = re.compile(r"(?:/private/var/folders/|/private/tmp/)")
ABSOLUTE_LOCAL_PATH_PATTERN: re.Pattern[str] = re.compile(
    r"(?:/Users/[A-Za-z0-9._-]+/|/home/[A-Za-z0-9._-]+/|[A-Za-z]:\\Users\\|"
    r"/private/var/folders/|/private/tmp/)"
)


@functools.lru_cache(maxsize=1)
def _tracked_text_files(root: Path = REPO_ROOT) -> tuple[tuple[str, str], ...]:
    """Return `(repo_relative_path, decoded_text)` for tracked non-binary files within size cap."""
    env = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}
    output = subprocess.check_output(["git", "-C", str(root), "ls-files", "-z"], env=env)
    records: list[tuple[str, str]] = []
    for raw_path in output.split(b"\0"):
        if not raw_path:
            continue
        rel_path = raw_path.decode("utf-8")
        path = root / rel_path
        if path.suffix.lower() in BINARY_SUFFIXES or not path.is_file() or path.is_symlink():
            continue
        try:
            size = path.stat().st_size
        except OSError:
            continue
        if size > MAX_SCAN_BYTES:
            continue
        data = path.read_bytes()
        if b"\0" in data[:8192]:
            continue
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            continue
        records.append((rel_path, text))
    return tuple(records)


def _scan_records(
    pattern: re.Pattern[str],
    records: tuple[tuple[str, str], ...] | None = None,
) -> dict[str, list[str]]:
    """Scan `(rel_path, text)` records for `pattern`, returning hits per relative path."""
    source_records = _tracked_text_files() if records is None else records
    hits: dict[str, list[str]] = {}
    for rel_path, text in source_records:
        file_hits: list[str] = []
        for match in pattern.finditer(text):
            line_no = text.count("\n", 0, match.start()) + 1
            line_start = text.rfind("\n", 0, match.start()) + 1
            line_end = text.find("\n", match.start())
            if line_end == -1:
                line_end = len(text)
            snippet = text[line_start:line_end].strip()
            if len(snippet) > 120:
                snippet = snippet[:117] + "..."
            file_hits.append(f"L{line_no} ({match.group(0)}): {snippet}")
        if file_hits:
            hits[rel_path] = file_hits
    return hits


def _scan_tracked_files(pattern: re.Pattern[str]) -> dict[str, list[str]]:
    """Scan tracked text files for `pattern`, returning hits per relative path."""
    return _scan_records(pattern)


def _unexpected_user_home_hits(
    records: tuple[tuple[str, str], ...] | None = None,
) -> dict[str, list[str]]:
    """Return user-home path hits outside the four-file provenance `ALLOWLIST`."""
    hits = _scan_records(USER_HOME_PATTERN, records)
    return {path: matches for path, matches in hits.items() if path not in ALLOWLIST}


def _unexpected_local_or_temp_hits(
    records: tuple[tuple[str, str], ...] | None = None,
) -> dict[str, list[str]]:
    """Return unallowlisted user-home or temp-root path hits.

    `ALLOWLIST` exempts the four hash-bound provenance files.
    `SYNTHETIC_FIXTURE_FILES` exempts only `MACOS_TEMP_ROOT_PATTERN` (never `USER_HOME_PATTERN`).
    """
    unexpected: dict[str, list[str]] = {}
    for path, matches in _unexpected_user_home_hits(records).items():
        unexpected.setdefault(path, []).extend(matches)
    temp_hits = _scan_records(MACOS_TEMP_ROOT_PATTERN, records)
    for path, matches in temp_hits.items():
        if path in ALLOWLIST or path in SYNTHETIC_FIXTURE_FILES:
            continue
        unexpected.setdefault(path, []).extend(matches)
    return unexpected


def test_allowlist_has_exact_four_provenance_records_with_reason_and_issue() -> None:
    """Allowlist must hold only the four hash-bound provenance files with reason and issue."""
    assert set(ALLOWLIST) == EXPECTED_ALLOWLIST_PATHS
    assert len(ALLOWLIST) == 4
    hits = _scan_tracked_files(ABSOLUTE_LOCAL_PATH_PATTERN)
    for rel_path, entry in ALLOWLIST.items():
        assert (REPO_ROOT / rel_path).is_file(), f"Allowlisted file missing: {rel_path}"
        assert entry.reason.strip(), f"Allowlist entry lacks reason: {rel_path}"
        assert re.fullmatch(r"#\d+", entry.follow_up_issue.strip()), (
            f"Allowlist entry lacks issue reference (#N): {rel_path}"
        )
        assert rel_path in hits, f"Stale allowlist entry (no matches remain): {rel_path}"


def test_synthetic_fixture_files_have_exact_two_test_files_with_reason() -> None:
    """Synthetic temp-root fixture exemption holds only the two rejection-test files."""
    assert set(SYNTHETIC_FIXTURE_FILES) == EXPECTED_SYNTHETIC_FIXTURE_PATHS
    assert len(SYNTHETIC_FIXTURE_FILES) == 2
    assert set(SYNTHETIC_FIXTURE_FILES).isdisjoint(ALLOWLIST)
    temp_hits = _scan_tracked_files(MACOS_TEMP_ROOT_PATTERN)
    home_hits = _scan_tracked_files(USER_HOME_PATTERN)
    for rel_path, reason in SYNTHETIC_FIXTURE_FILES.items():
        assert (REPO_ROOT / rel_path).is_file(), f"Fixture file missing: {rel_path}"
        assert reason == "synthetic temp-root paths used to test rejection"
        assert rel_path in temp_hits, f"Stale fixture exemption (no temp-root hits): {rel_path}"
        assert rel_path not in home_hits, f"Fixture file contains user-home path: {rel_path}"


def test_exemptions_do_not_widen_to_unknown_temp_files_or_user_home_in_fixtures() -> None:
    """Unknown test files with temp paths and fixture files with user-home paths must fail."""
    unknown_temp = (
        ("tests/test_unknown_temp.py", 'path = "/" "' + "private/tmp/mtp-leak/out.tsv" + '"'),
        ("tests/test_unknown_var.py", 'path = "/" "' + "private/var/folders/ab/cd/T/x" + '"'),
    )
    synthetic_records = (
        ("tests/test_unknown_temp.py", "/" + "private/tmp/mtp-leak/out.tsv"),
        ("tests/test_unknown_var.py", "/" + "private/var/folders/ab/cd/T/x"),
        ("tests/test_output_roots.py", "/" + "private/tmp/core"),
        ("tests/test_three_test_stage.py", "/" + "private/tmp/mtp-table"),
    )
    assert unknown_temp  # keep tuple non-empty
    unexpected_temp = _unexpected_local_or_temp_hits(synthetic_records)
    assert set(unexpected_temp) == {"tests/test_unknown_temp.py", "tests/test_unknown_var.py"}

    fixture_home_records = (
        ("tests/test_output_roots.py", "/" + "Users/alice/Documents/mtp/out.tsv"),
        ("tests/test_three_test_stage.py", "/" + "home/bob/work/mtp/out.tsv"),
    )
    assert set(_unexpected_user_home_hits(fixture_home_records)) == {
        "tests/test_output_roots.py",
        "tests/test_three_test_stage.py",
    }
    assert set(_unexpected_local_or_temp_hits(fixture_home_records)) == {
        "tests/test_output_roots.py",
        "tests/test_three_test_stage.py",
    }


def test_absolute_local_path_pattern_matches_all_target_families() -> None:
    """Regex must catch macOS/Linux/Windows user homes and macOS temp roots."""
    macos_home = "/" + "Users/alice/Documents/project/file.tsv"
    linux_home = "/" + "home/bob/work/project/file.tsv"
    windows_home = "C:\\" + "Users\\carol\\Documents\\file.tsv"
    macos_var_folders = "/" + "private/var/folders/ab/cd1234/T/run/out.json"
    macos_private_tmp = "/" + "private/tmp/mtp-run/out.json"

    for sample in (
        macos_home,
        linux_home,
        windows_home,
        macos_var_folders,
        macos_private_tmp,
    ):
        assert ABSOLUTE_LOCAL_PATH_PATTERN.search(sample) is not None, sample

    for safe in (
        "data/derived/ora_provenance_refresh/issue720",
        "<local checkout>/steps/ora_profiles/derived/study_trait_profiles_ora_wide.tsv",
        "<temp dir>/mtp-workflow-629",
        "/Users/",
        "/home/",
        "/private/tmp",
        "/private/var/folders",
    ):
        assert ABSOLUTE_LOCAL_PATH_PATTERN.search(safe) is None, safe


def test_tracked_files_contain_no_unallowlisted_user_home_paths() -> None:
    """No tracked text file outside the allowlist may contain a user home path."""
    unexpected = _unexpected_user_home_hits()
    assert not unexpected, "Tracked files contain unallowlisted user-home paths:\n" + "\n".join(
        f"  {path}: {matches[:3]} (total {len(matches)})"
        for path, matches in sorted(unexpected.items())
    )


def test_tracked_files_contain_no_unallowlisted_absolute_local_paths() -> None:
    """No tracked text file outside the exemptions may contain a local home or macOS temp path."""
    unexpected = _unexpected_local_or_temp_hits()
    assert not unexpected, (
        "Tracked files contain unallowlisted absolute local or temp paths:\n"
        + "\n".join(
            f"  {path}: {matches[:3]} (total {len(matches)})"
            for path, matches in sorted(unexpected.items())
        )
    )
