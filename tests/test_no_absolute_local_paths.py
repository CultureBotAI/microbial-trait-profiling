"""Guard tracked text files against machine-local user-home and temporary paths (#1622).

What the guard matches:

- user homes: macOS ``Users`` and Linux ``home`` directories followed by a user name,
  with or without a trailing slash (a JSON value can end at the name), and Windows
  homes after a drive letter with back, forward, mixed or JSON-escaped separators;
- macOS temp roots: the per-user ``var/folders`` directory, with or without the
  ``/private`` prefix (``$TMPDIR`` and ``tempfile.gettempdir()`` give it without),
  and ``tmp`` under ``/private``.

The home pattern has no left anchor on purpose (#1691): an anchor would miss a home
reached through the macOS firmlink root (``System/Volumes/Data``), a form that disk
tools print. A false positive inside a URL fails loudly, which is the safer error.

What it deliberately leaves out, tracked in #1682: a bare ``tmp`` directory at the
root (the tree uses it in synthetic examples), ``~``-relative homes and external
volumes.

No pattern or sample in this file spells out a path it matches, because the guard
scans this file too (#1674).

Every tracked path is accounted for (#1676): it is scanned, skipped for a binary
suffix, or skipped as larger than ``MAX_SCAN_BYTES`` and named in ``SKIPPED_LARGE``.
A tracked text file with NUL bytes or invalid UTF-8 fails the guard instead of being
skipped, and a file with a binary suffix must really be binary (#1692). Every index
version of a path is scanned (all conflict stages, symlink targets and files missing
from the working tree), and so is the working-tree file when it differs (#1689).

Git runs without any inherited ``GIT_*`` variable (#1687): a hook exports ``GIT_DIR``,
which would otherwise point the synthetic-repository tests at the caller's index.
"""

from __future__ import annotations

import collections
import functools
import os
import re
import subprocess
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pytest

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

# Tracked text files over MAX_SCAN_BYTES. The set must equal what the scan skips for size.
SKIPPED_LARGE: frozenset[str] = frozenset(
    {"trait-visualization/study_trait_profile_clustered_heatmap_matrix.tsv"}
)

# Files the scan must read: the four docs redacted under #1622, this guard and core files.
REQUIRED_SCANNED: frozenset[str] = frozenset(
    {
        "Justfile",
        "README.md",
        "docs/evidence-archives.md",
        "docs/reviews/2026-10-17-prompt-1-implementation-report.md",
        "docs/reviews/2026-10-19-core-workflow-implementation.md",
        "docs/reviews/2026-10-19-ora-provenance-refresh.md",
        "pyproject.toml",
        "src/mtp/__init__.py",
        "tests/test_no_absolute_local_paths.py",
    }
)

# Regeneration follow-up issue for the four allowlisted provenance records.
PROVENANCE_REGENERATION_ISSUE = "#1664"


@dataclass(frozen=True)
class AllowlistEntry:
    """Explicit exemption for a recorded provenance file awaiting regeneration."""

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
            "Run output recorded as written: the integrated-evidence assessment sidecar from "
            "mtp.integrated_evidence._bound_assessment, with resolved source paths. Editing it "
            "would misstate what the run wrote; regenerate under #1664."
        ),
        follow_up_issue=PROVENANCE_REGENERATION_ISSUE,
    ),
    "data/derived/integrated_evidence/integrated_evidence_provenance.json": AllowlistEntry(
        reason=(
            "Run output recorded as written: the integrated-evidence provenance record from "
            "mtp.integrated_evidence.main, with resolved native_sources and assessment_status "
            "source paths. Editing it would misstate what the run wrote; regenerate under #1664."
        ),
        follow_up_issue=PROVENANCE_REGENERATION_ISSUE,
    ),
    "docs/reviews/2026-10-19-ora-provenance-refresh-audit.json": AllowlistEntry(
        reason=(
            "Audit output recorded as written by scripts/adhoc/issue720_audit_ora_refresh.py, "
            "with resolved manifest and input paths. Its manifest_sha256 binds the inputs "
            "manifest, so regenerate the two together under #1664."
        ),
        follow_up_issue=PROVENANCE_REGENERATION_ISSUE,
    ),
    "docs/reviews/2026-10-19-ora-provenance-refresh-inputs.json": AllowlistEntry(
        reason=(
            "Hash-bound input manifest: its SHA-256 is the manifest_sha256 recorded in "
            "docs/reviews/2026-10-19-ora-provenance-refresh-audit.json. Regenerate the two "
            "together under #1664."
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

# Hits per family that each exempt file holds today (#1678). Counts, not hashes, so that
# detect-secrets has nothing to flag; any added or removed local path changes a count.
EXPECTED_EXEMPT_HIT_COUNTS: dict[str, dict[str, int]] = {
    "data/derived/integrated_evidence/integrated_evidence_assessments.tsv": {"user_home": 2},
    "data/derived/integrated_evidence/integrated_evidence_provenance.json": {"user_home": 8},
    "docs/reviews/2026-10-19-ora-provenance-refresh-audit.json": {"temp_root": 13, "user_home": 8},
    "docs/reviews/2026-10-19-ora-provenance-refresh-inputs.json": {"temp_root": 12, "user_home": 8},
    "tests/test_output_roots.py": {"temp_root": 3},
    "tests/test_three_test_stage.py": {"temp_root": 1},
}

# Each alternative is written so that its own source text does not match (#1674).
_USER_HOME_SOURCE = r"/(?:Users|home)/[A-Za-z0-9._-]+|[A-Za-z]:[\\/]+Users[\\/]+"
# The optional prefix only widens the reported match. The trailing slash is a lookahead,
# so a user home directly after a temp root is still found.
_TEMP_ROOT_SOURCE = r"(?:/private)?/var/(?:folders)(?=/)|/private/(?:tmp)(?=/)"

USER_HOME_PATTERN: re.Pattern[str] = re.compile(_USER_HOME_SOURCE)
MACOS_TEMP_ROOT_PATTERN: re.Pattern[str] = re.compile(_TEMP_ROOT_SOURCE)
# Both families in one pass; the named group says which family a match belongs to.
ABSOLUTE_LOCAL_PATH_PATTERN: re.Pattern[str] = re.compile(
    f"(?P<user_home>{_USER_HOME_SOURCE})|(?P<temp_root>{_TEMP_ROOT_SOURCE})"
)

USER_HOME = "user_home"
TEMP_ROOT = "temp_root"

NUL_REASON = "NUL byte in a file without a binary suffix"
UTF8_REASON = "not valid UTF-8"


@dataclass(frozen=True)
class TrackedCorpus:
    """Every tracked path, sorted into what was scanned and why the rest was not.

    `texts` holds one `(path, text)` record per distinct version of a scanned path.
    """

    texts: tuple[tuple[str, str], ...]
    binary_suffix: frozenset[str]
    large: frozenset[str]
    unscannable: dict[str, str]
    working_tree_versions: frozenset[str]

    @property
    def scanned(self) -> frozenset[str]:
        """Return the repository-relative paths whose text was scanned."""
        return frozenset(path for path, _ in self.texts)


@dataclass(frozen=True)
class Hit:
    """One local-path match in a tracked file."""

    family: str
    line_no: int
    matched: str
    snippet: str

    def __str__(self) -> str:
        """Render the hit for an assertion message."""
        return f"L{self.line_no} {self.family} ({self.matched}): {self.snippet}"


def _git_env() -> dict[str, str]:
    """Return an environment with no inherited `GIT_*` variable and no user or system config."""
    env = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    env["GIT_CONFIG_GLOBAL"] = os.devnull
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    return env


def _git(root: Path, *args: str, stdin: bytes | None = None) -> bytes:
    """Run git in `root` and return its standard output."""
    completed = subprocess.run(
        ["git", "-C", str(root), *args],
        env=_git_env(),
        input=stdin,
        capture_output=True,
        check=True,
    )
    return completed.stdout


def _index_entries(root: Path) -> list[tuple[str, str, str]]:
    """Return `(path, mode, object_id)` for every index entry, including conflict stages."""
    entries: list[tuple[str, str, str]] = []
    for record in _git(root, "ls-files", "-s", "-z").split(b"\0"):
        if not record:
            continue
        meta, _, raw_path = record.partition(b"\t")
        mode, object_id, _stage = meta.decode("ascii").split(" ")
        entries.append((raw_path.decode("utf-8"), mode, object_id))
    return entries


def _blob_sizes(root: Path, object_ids: list[str]) -> dict[str, int]:
    """Return the size of each blob, failing on a missing object."""
    if not object_ids:
        return {}
    output = _git(root, "cat-file", "--batch-check", stdin=("\n".join(object_ids) + "\n").encode())
    sizes: dict[str, int] = {}
    for line in output.decode("ascii").splitlines():
        object_id, object_type, size = line.split(" ")
        assert object_type == "blob", f"Index entry is not a blob: {line}"
        sizes[object_id] = int(size)
    return sizes


def _iter_blobs(root: Path, object_ids: list[str]) -> Iterator[tuple[str, bytes]]:
    """Yield `(object_id, content)` one blob at a time from a single `git cat-file --batch`."""
    if not object_ids:
        return
    with subprocess.Popen(
        ["git", "-C", str(root), "cat-file", "--batch"],
        env=_git_env(),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
    ) as process:
        assert process.stdin is not None
        assert process.stdout is not None
        for object_id in object_ids:
            process.stdin.write(object_id.encode("ascii") + b"\n")
            process.stdin.flush()
            header = process.stdout.readline().decode("ascii").split(" ")
            assert len(header) == 3, f"git cat-file could not read {object_id}: {header}"
            size = int(header[2])
            content = process.stdout.read(size)
            terminator = process.stdout.read(1)
            assert len(content) == size, f"Short read of {object_id}"
            assert terminator == b"\n", f"Missing record terminator after {object_id}"
            yield object_id, content
        process.stdin.close()
    assert process.returncode == 0, f"git cat-file exited {process.returncode}"


def _unreadable_reason(data: bytes) -> str | None:
    """Return why `data` cannot be scanned as text, or None if it can."""
    if b"\0" in data:
        return NUL_REASON
    try:
        data.decode("utf-8")
    except UnicodeDecodeError:
        return UTF8_REASON
    return None


def _read_tracked(root: Path, max_bytes: int) -> TrackedCorpus:
    """Read every version of every tracked path under `root` and sort it."""
    binary_suffix: set[str] = set()
    large: set[str] = set()
    unscannable: dict[str, str] = {}
    index_versions: dict[str, list[str]] = {}
    for rel_path, mode, object_id in _index_entries(root):
        if mode == "160000":
            unscannable[rel_path] = "submodule"
        elif Path(rel_path).suffix.lower() in BINARY_SUFFIXES:
            binary_suffix.add(rel_path)
        else:
            index_versions.setdefault(rel_path, []).append(object_id)

    sizes = _blob_sizes(root, sorted({oid for oids in index_versions.values() for oid in oids}))
    working_files: dict[str, Path] = {}
    for rel_path, object_ids in index_versions.items():
        path = root / rel_path
        working = path if path.is_file() and not path.is_symlink() else None
        working_size = working.stat().st_size if working is not None else 0
        if working_size > max_bytes or any(sizes[oid] > max_bytes for oid in object_ids):
            large.add(rel_path)
        elif working is not None:
            working_files[rel_path] = working

    readable = {path: oids for path, oids in index_versions.items() if path not in large}
    blobs = dict(_iter_blobs(root, sorted({oid for oids in readable.values() for oid in oids})))
    texts: list[tuple[str, str]] = []
    extra_versions: set[str] = set()
    for rel_path in sorted(readable):
        versions = [blobs[oid] for oid in readable[rel_path]]
        if rel_path in working_files:
            working_data = working_files[rel_path].read_bytes()
            if working_data not in versions:
                versions.append(working_data)
                extra_versions.add(rel_path)
        reasons = [reason for reason in map(_unreadable_reason, versions) if reason]
        if reasons:
            unscannable[rel_path] = reasons[0]
            extra_versions.discard(rel_path)
            continue
        texts.extend((rel_path, data.decode("utf-8")) for data in dict.fromkeys(versions))
    return TrackedCorpus(
        texts=tuple(texts),
        binary_suffix=frozenset(binary_suffix),
        large=frozenset(large),
        unscannable=unscannable,
        working_tree_versions=frozenset(extra_versions),
    )


def _text_files_with_binary_suffix(root: Path) -> list[str]:
    """Return tracked paths that have a binary suffix but whose index content is text."""
    paths_by_object: dict[str, list[str]] = {}
    for rel_path, mode, object_id in _index_entries(root):
        if mode != "160000" and Path(rel_path).suffix.lower() in BINARY_SUFFIXES:
            paths_by_object.setdefault(object_id, []).append(rel_path)
    found: list[str] = []
    for object_id, data in _iter_blobs(root, sorted(paths_by_object)):
        if _unreadable_reason(data) is None:
            found.extend(paths_by_object[object_id])
    return sorted(found)


@functools.lru_cache(maxsize=1)
def _tracked_corpus() -> TrackedCorpus:
    """Read the repository's tracked files once per test session."""
    return _read_tracked(REPO_ROOT, MAX_SCAN_BYTES)


def _scan_records(records: tuple[tuple[str, str], ...]) -> dict[str, tuple[Hit, ...]]:
    """Scan `(rel_path, text)` records once for both path families; versions are merged."""
    hits: dict[str, list[Hit]] = {}
    for rel_path, text in records:
        for match in ABSOLUTE_LOCAL_PATH_PATTERN.finditer(text):
            family = USER_HOME if match.group(USER_HOME) is not None else TEMP_ROOT
            line_start = text.rfind("\n", 0, match.start()) + 1
            line_end = text.find("\n", match.start())
            snippet = text[line_start : len(text) if line_end == -1 else line_end].strip()
            if len(snippet) > 120:
                snippet = snippet[:117] + "..."
            line_no = text.count("\n", 0, match.start()) + 1
            hits.setdefault(rel_path, []).append(Hit(family, line_no, match.group(0), snippet))
    return {path: tuple(path_hits) for path, path_hits in hits.items()}


@functools.lru_cache(maxsize=1)
def _tracked_hits() -> dict[str, tuple[Hit, ...]]:
    """Scan the tracked corpus once per test session."""
    return _scan_records(_tracked_corpus().texts)


def _unexpected_hits(
    hits: dict[str, tuple[Hit, ...]], families: frozenset[str]
) -> dict[str, list[Hit]]:
    """Return hits of `families` outside the exemptions.

    `ALLOWLIST` exempts the four provenance records from both families.
    `SYNTHETIC_FIXTURE_FILES` exempts only temp-root hits, never user-home hits.
    `EXPECTED_EXEMPT_HIT_COUNTS` pins how many hits each exempt file may hold.
    """
    unexpected: dict[str, list[Hit]] = {}
    for path, file_hits in hits.items():
        if path in ALLOWLIST:
            continue
        for hit in file_hits:
            if hit.family not in families:
                continue
            if hit.family == TEMP_ROOT and path in SYNTHETIC_FIXTURE_FILES:
                continue
            unexpected.setdefault(path, []).append(hit)
    return unexpected


def _exempt_count_mismatches(
    records: tuple[tuple[str, str], ...],
) -> dict[str, list[dict[str, int]]]:
    """Return, per exempt path, the family counts of each version that differs from its pin."""
    mismatches: dict[str, list[dict[str, int]]] = {}
    for rel_path, text in records:
        if rel_path not in EXPECTED_EXEMPT_HIT_COUNTS:
            continue
        version_hits = _scan_records(((rel_path, text),)).get(rel_path, ())
        counts = dict(collections.Counter(hit.family for hit in version_hits))
        if counts != EXPECTED_EXEMPT_HIT_COUNTS[rel_path]:
            mismatches.setdefault(rel_path, []).append(counts)
    return mismatches


def _format_unexpected(unexpected: dict[str, list[Hit]]) -> str:
    """Render unexpected hits, three per file, for an assertion message."""
    return "\n".join(
        f"  {path}: {[str(hit) for hit in hits[:3]]} (total {len(hits)})"
        for path, hits in sorted(unexpected.items())
    )


def _build_synthetic_repo(repo: Path) -> Path:
    """Create a git repository with one entry of each kind the reader must handle."""
    repo.mkdir()
    _git(repo, "init", "-q")
    (repo / "plain.md").write_text("no local paths here\n", encoding="utf-8")
    (repo / "sparse.md").write_text("cwd=/" + "home/bob\n", encoding="utf-8")
    (repo / "cleaned.md").write_text("cwd=/" + "home/bob\n", encoding="utf-8")
    (repo / "figure.png").write_bytes(b"\x89PNG\0binary")
    (repo / "fake.png").write_text("text under a binary suffix\n", encoding="utf-8")
    (repo / "big.tsv").write_text("x" * 30_000, encoding="utf-8")
    (repo / "big_sparse.tsv").write_text("y" * 30_000, encoding="utf-8")
    (repo / "nul.txt").write_bytes(b"text\0more")
    (repo / "late_nul.txt").write_bytes(b"a" * 9_000 + b"\0")
    (repo / "latin1.txt").write_bytes(b"A\xe9B")  # Latin-1 e-acute: not valid UTF-8
    os.symlink("/" + "Users/alice/data", repo / "link")
    _git(repo, "add", "-A")
    (repo / "sparse.md").unlink()
    (repo / "big_sparse.tsv").unlink()
    (repo / "cleaned.md").write_text("cleaned in the working tree only\n", encoding="utf-8")
    return repo


# Samples are assembled from pieces so that this file's own text does not match (#1674).
FAMILY_SAMPLES: tuple[tuple[str, str], ...] = (
    (USER_HOME, "/" + "Users/alice/Documents/project/file.tsv"),
    (USER_HOME, '{"root": "/' + 'Users/alice"}'),
    (USER_HOME, "/" + "home/bob/work/project/file.tsv"),
    (USER_HOME, "cwd=/" + "home/bob"),
    (USER_HOME, "C:\\" + "Users\\carol\\Documents\\file.tsv"),
    (USER_HOME, "C:/" + "Users/carol/Documents/file.tsv"),
    (USER_HOME, "C:/" + "Users\\carol\\file.tsv"),
    (USER_HOME, "C:\\" + "Users/carol/file.tsv"),
    (USER_HOME, '{"path": "C:\\\\' + 'Users\\\\carol\\\\file.tsv"}'),
    (TEMP_ROOT, "/" + "private/var" + "/" + "folders/ab/cd1234/T/run/out.json"),
    (TEMP_ROOT, "/" + "var/folders/ab/cd1234/T/run/out.json"),
    (TEMP_ROOT, "/" + "private/tmp/mtp-run/out.json"),
)

SAFE_SAMPLES: tuple[str, ...] = (
    "data/derived/ora_provenance_refresh/issue720",
    "<local checkout>/steps/ora_profiles/derived/study_trait_profiles_ora_wide.tsv",
    "<temp dir>/mtp-workflow-629",
    "/Users/",
    "/home/",
    "/Users/<name>/project",
    "/path/to/mtp-workflow-629/.tmp/core-production-final",
    "/absolute/path/to/a/new/ora-refresh",
    "/private/tmp",
    "/private/var/folders",
    "https://github.com/CultureBotAI/microbial-trait-profiling/issues/1622",
)

# Left out on purpose and tracked in #1682; a change here should update that issue.
OUT_OF_SCOPE_SAMPLES: tuple[str, ...] = (
    "/" + "tmp/mtp-run/out.json",
    "~/" + ".cache/uv",
    "/" + "Volumes/External/mtp/out.tsv",
)


def test_each_enforcing_pattern_matches_its_family_and_not_the_other() -> None:
    """Each enforcing pattern catches its own family; the combined scan labels it right."""
    patterns = {USER_HOME: USER_HOME_PATTERN, TEMP_ROOT: MACOS_TEMP_ROOT_PATTERN}
    for family, sample in FAMILY_SAMPLES:
        assert patterns[family].search(sample) is not None, (family, sample)
        other = TEMP_ROOT if family == USER_HOME else USER_HOME
        assert patterns[other].search(sample) is None, (other, sample)
        hits = _scan_records((("sample.txt", sample),))
        assert [hit.family for hit in hits["sample.txt"]] == [family], sample


def test_safe_and_out_of_scope_samples_match_no_pattern() -> None:
    """Placeholders, partial roots and the documented out-of-scope forms do not match."""
    for sample in SAFE_SAMPLES + OUT_OF_SCOPE_SAMPLES:
        for pattern in (USER_HOME_PATTERN, MACOS_TEMP_ROOT_PATTERN, ABSOLUTE_LOCAL_PATH_PATTERN):
            assert pattern.search(sample) is None, (pattern.pattern, sample)


def test_user_home_after_a_temp_root_is_still_reported() -> None:
    """A temp-root match must not consume the slash that starts a following user home."""
    sample = "/" + "private/tmp" + "/" + "Users/alice/out.tsv"
    hits = _scan_records((("tests/test_output_roots.py", sample),))
    assert [hit.family for hit in hits["tests/test_output_roots.py"]] == [TEMP_ROOT, USER_HOME]
    unexpected = _unexpected_hits(hits, frozenset({USER_HOME, TEMP_ROOT}))
    assert [hit.family for hit in unexpected["tests/test_output_roots.py"]] == [USER_HOME]


def test_read_tracked_accounts_for_every_kind_of_entry(tmp_path: Path) -> None:
    """Sparse, cleaned, symlinked, binary, large and undecodable entries are each accounted for."""
    repo = _build_synthetic_repo(tmp_path / "repo")

    corpus = _read_tracked(repo, max_bytes=20_000)

    assert corpus.scanned == {"plain.md", "sparse.md", "cleaned.md", "link"}
    assert corpus.working_tree_versions == {"cleaned.md"}
    assert corpus.binary_suffix == {"figure.png", "fake.png"}
    assert corpus.large == {"big.tsv", "big_sparse.tsv"}
    assert corpus.unscannable == {
        "nul.txt": NUL_REASON,
        "late_nul.txt": NUL_REASON,
        "latin1.txt": UTF8_REASON,
    }
    hits = _scan_records(corpus.texts)
    assert {path: [hit.family for hit in path_hits] for path, path_hits in hits.items()} == {
        "sparse.md": [USER_HOME],
        "cleaned.md": [USER_HOME],
        "link": [USER_HOME],
    }
    assert _text_files_with_binary_suffix(repo) == ["fake.png"]


def test_git_ignores_an_inherited_git_dir_and_index(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A hook's `GIT_DIR` and `GIT_INDEX_FILE` must not redirect the synthetic repo (#1687)."""
    caller = tmp_path / "caller"
    caller.mkdir()
    _git(caller, "init", "-q")
    (caller / "kept.md").write_text("kept\n", encoding="utf-8")
    _git(caller, "add", "kept.md")
    caller_index = caller / ".git" / "index"
    before = caller_index.read_bytes()
    monkeypatch.setenv("GIT_DIR", str(caller / ".git"))
    monkeypatch.setenv("GIT_INDEX_FILE", str(caller_index))
    monkeypatch.setenv("GIT_WORK_TREE", str(caller))

    corpus = _read_tracked(_build_synthetic_repo(tmp_path / "repo"), max_bytes=20_000)

    assert caller_index.read_bytes() == before
    assert "kept.md" not in corpus.scanned
    assert "plain.md" in corpus.scanned


def test_every_tracked_path_is_scanned_or_skipped_by_name() -> None:
    """The scan covers the whole tree; each skip is a binary file or a named large file."""
    corpus = _tracked_corpus()
    assert not corpus.unscannable, f"Tracked text files the guard cannot read: {corpus.unscannable}"
    assert corpus.large == SKIPPED_LARGE, (
        f"Files over {MAX_SCAN_BYTES} bytes differ from SKIPPED_LARGE: "
        f"{sorted(corpus.large ^ SKIPPED_LARGE)}"
    )
    missing = REQUIRED_SCANNED - corpus.scanned
    assert not missing, f"Required files were not scanned: {sorted(missing)}"
    tracked = {path for path, _, _ in _index_entries(REPO_ROOT)}
    accounted = corpus.scanned | corpus.binary_suffix | corpus.large
    assert accounted == tracked, f"Unaccounted tracked paths: {sorted(tracked ^ accounted)}"
    text_as_binary = _text_files_with_binary_suffix(REPO_ROOT)
    assert not text_as_binary, f"Text files skipped for a binary suffix: {text_as_binary}"


def test_allowlist_has_exact_four_provenance_records_with_reason_and_issue() -> None:
    """Allowlist must hold only the four provenance records, each with reason and issue."""
    assert set(ALLOWLIST) == EXPECTED_ALLOWLIST_PATHS
    assert len(ALLOWLIST) == 4
    hits = _tracked_hits()
    for rel_path, entry in ALLOWLIST.items():
        assert rel_path in _tracked_corpus().scanned, f"Allowlisted file not scanned: {rel_path}"
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
    hits = _tracked_hits()
    for rel_path, reason in SYNTHETIC_FIXTURE_FILES.items():
        assert rel_path in _tracked_corpus().scanned, f"Fixture file not scanned: {rel_path}"
        assert reason == "synthetic temp-root paths used to test rejection"
        families = {hit.family for hit in hits.get(rel_path, ())}
        assert TEMP_ROOT in families, f"Stale fixture exemption (no temp-root hits): {rel_path}"
        assert USER_HOME not in families, f"Fixture file contains user-home path: {rel_path}"


def test_exempt_files_hold_exactly_their_recorded_hit_counts() -> None:
    """Each version of an exempt file holds exactly its pinned hits per family (#1678).

    A new local path added to an exempt file changes a count and fails here, so the
    per-file exemptions cannot silently absorb new paths.
    """
    assert set(EXPECTED_EXEMPT_HIT_COUNTS) == set(ALLOWLIST) | set(SYNTHETIC_FIXTURE_FILES)
    mismatches = _exempt_count_mismatches(_tracked_corpus().texts)
    assert not mismatches, f"Exempt files whose hit counts changed: {mismatches}"
    one_more = (
        ("tests/test_three_test_stage.py", "/" + "var/folders/a/T/x\n/" + "private/tmp/y\n"),
        ("tests/test_output_roots.py", ("/" + "private/tmp/a\n") * 3),
    )
    assert _exempt_count_mismatches(one_more) == {
        "tests/test_three_test_stage.py": [{TEMP_ROOT: 2}],
    }


def test_exemptions_do_not_widen_to_unknown_temp_files_or_user_home_in_fixtures() -> None:
    """Unknown test files with temp paths and fixture files with user-home paths must fail."""
    both = frozenset({USER_HOME, TEMP_ROOT})
    synthetic_records = (
        ("tests/test_unknown_temp.py", "/" + "private/tmp/mtp-leak/out.tsv"),
        ("tests/test_unknown_var.py", "/" + "private/var" + "/" + "folders/ab/cd/T/x"),
        ("tests/test_unknown_unresolved.py", "/" + "var/folders/ab/cd/T/x"),
        ("tests/test_output_roots.py", "/" + "private/tmp/core"),
        ("tests/test_three_test_stage.py", "/" + "var/folders/ab/T/mtp-table"),
    )
    unexpected_temp = _unexpected_hits(_scan_records(synthetic_records), both)
    assert set(unexpected_temp) == {
        "tests/test_unknown_temp.py",
        "tests/test_unknown_var.py",
        "tests/test_unknown_unresolved.py",
    }

    fixture_home_records = (
        ("tests/test_output_roots.py", "/" + "Users/alice/Documents/mtp/out.tsv"),
        ("tests/test_three_test_stage.py", "/" + "home/bob/work/mtp/out.tsv"),
    )
    fixture_home_hits = _scan_records(fixture_home_records)
    for families in (frozenset({USER_HOME}), both):
        assert set(_unexpected_hits(fixture_home_hits, families)) == {
            "tests/test_output_roots.py",
            "tests/test_three_test_stage.py",
        }


def test_tracked_files_contain_no_unallowlisted_user_home_paths() -> None:
    """No tracked text file outside the allowlist may contain a user home path."""
    unexpected = _unexpected_hits(_tracked_hits(), frozenset({USER_HOME}))
    assert not unexpected, "Tracked files contain unallowlisted user-home paths:\n" + (
        _format_unexpected(unexpected)
    )


def test_tracked_files_contain_no_unallowlisted_absolute_local_paths() -> None:
    """No tracked text file outside the exemptions may contain a local home or macOS temp path."""
    unexpected = _unexpected_hits(_tracked_hits(), frozenset({USER_HOME, TEMP_ROOT}))
    assert not unexpected, "Tracked files contain unallowlisted absolute local or temp paths:\n" + (
        _format_unexpected(unexpected)
    )
