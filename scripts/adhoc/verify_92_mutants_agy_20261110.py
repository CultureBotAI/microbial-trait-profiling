"""Mutant verifier for Issue #92 download-from-release fallback (2026-11-10).

Purpose:
    Verify that the unit tests in tests/test_download_from_release_92.py kill
    >= 7 non-trivial production mutants of mtp.download release-tarball
    restoration by raising AssertionError.
Inputs:
    In-memory monkeypatched mutations of mtp.download and temporary test
    directories.
Outputs:
    Per-mutant kill status printed to stdout; non-zero exit if any mutant
    survives.
Issue: #92
Date: 2026-11-10
"""

from __future__ import annotations

import importlib.util
import tempfile
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType

import pytest

from mtp import download as dl
from mtp.paths import REPO_ROOT


def _load_test_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "test_download_from_release_92",
        REPO_ROOT / "tests/test_download_from_release_92.py",
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load tests/test_download_from_release_92.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


test_mod = _load_test_module()


@contextmanager
def _patch_attr(target: object, attr: str, replacement: object) -> Iterator[None]:
    original = getattr(target, attr)
    setattr(target, attr, replacement)
    try:
        yield
    finally:
        setattr(target, attr, original)


def _run_target_tests() -> None:
    with tempfile.TemporaryDirectory(prefix="mtp-verify-92-") as tmp_str:
        base = Path(tmp_str)
        mp = pytest.MonkeyPatch()
        try:
            test_mod.test_download_from_release_extracts_csvs_and_writes_provenance(base / "t1", mp)
        finally:
            mp.undo()

        mp = pytest.MonkeyPatch()
        try:
            test_mod.test_download_from_release_passes_explicit_tag_and_verifies_pinned_match(
                base / "t2", mp
            )
        finally:
            mp.undo()

        mp = pytest.MonkeyPatch()
        try:
            test_mod.test_download_from_release_rejects_corrupt_pinned_size_dump(base / "t3", mp)
        finally:
            mp.undo()

        for idx, bad_member in enumerate(
            ("../evil.csv", "/etc/passwd", "subdir/full_dump.csv", "unexpected.csv")
        ):
            mp = pytest.MonkeyPatch()
            try:
                test_mod.test_download_from_release_rejects_unsafe_or_unexpected_tar_members(
                    base / f"t4_{idx}", mp, bad_member
                )
            finally:
                mp.undo()

        mp = pytest.MonkeyPatch()
        try:
            test_mod.test_download_from_release_rejects_symlink_tar_member(base / "t5", mp)
        finally:
            mp.undo()

        mp = pytest.MonkeyPatch()
        try:
            test_mod.test_download_from_release_rejects_missing_full_dump_and_bad_args(
                base / "t6", mp
            )
        finally:
            mp.undo()

        mp = pytest.MonkeyPatch()
        try:
            test_mod.test_fetch_release_archive_handles_missing_gh_and_nonzero_exit(base / "t7", mp)
        finally:
            mp.undo()

        test_mod.test_justfile_declares_download_from_release_recipe()


def _expect_killed(name: str, mutate_and_run: Callable[[], None]) -> None:
    try:
        mutate_and_run()
    except AssertionError:
        print(f"KILLED: {name}")
        return
    raise RuntimeError(f"SURVIVED: {name}")


def main() -> int:
    # Baseline check: unmutated code passes all target tests
    _run_target_tests()

    # Mutant 1: SOURCE_RELEASE is recorded as "live" instead of "release"
    def _m1() -> None:
        with _patch_attr(dl, "SOURCE_RELEASE", "live"):
            _run_target_tests()

    _expect_killed("M1_wrong_source_in_provenance", _m1)

    # Mutant 2: _validate_release_tag forwards "latest" instead of normalizing to None
    def _m2() -> None:
        def _bad_validate_tag(tag: str | None) -> str | None:
            return tag

        with _patch_attr(dl, "_validate_release_tag", _bad_validate_tag):
            _run_target_tests()

    _expect_killed("M2_latest_tag_forwarded_literally", _m2)

    # Mutant 3: verify_pinned_dump is a no-op on corrupt pinned-size dump
    def _m3() -> None:
        with _patch_attr(dl, "verify_pinned_dump", lambda _b: None):
            _run_target_tests()

    _expect_killed("M3_pinned_integrity_check_skipped", _m3)

    # Mutant 4: _read_release_archive_members skips unsafe-path validation
    def _m4() -> None:
        def _lenient_members(_archive_path: Path) -> dict[str, bytes]:
            return {
                dl.FULL_DUMP_NAME: test_mod.SNAPSHOT_HEADER
                + b"BSDB ID,Study\nbsdb:1/1/1,StudyOne\n",
                "studies.csv": b"Study,State\nbsdb:1,Complete\n",
                "experiments.csv": b"Study,Experiment\nbsdb:1,bsdb:1/1\n",
                "signatures.csv": b"Study,Experiment,Signature\nbsdb:1,bsdb:1/1,bsdb:1/1/1\n",
            }

        with _patch_attr(dl, "_read_release_archive_members", _lenient_members):
            _run_target_tests()

    _expect_killed("M4_unsafe_tar_member_check_skipped", _m4)

    # Mutant 5: _infer_release_tag ignores archive tag and always returns "unknown"
    def _m5() -> None:
        with _patch_attr(dl, "_infer_release_tag", lambda _name, _req: "unknown"):
            _run_target_tests()

    _expect_killed("M5_infer_release_tag_broken", _m5)

    # Mutant 6: _validate_release_repo accepts empty or invalid repo strings
    def _m6() -> None:
        with _patch_attr(dl, "_validate_release_repo", lambda repo: repo):
            _run_target_tests()

    _expect_killed("M6_validate_release_repo_skipped", _m6)

    # Mutant 7: parse_snapshot ignores snapshot header and returns None
    def _m7() -> None:
        with _patch_attr(dl, "parse_snapshot", lambda _content: None):
            _run_target_tests()

    _expect_killed("M7_snapshot_header_ignored", _m7)

    print("All 7 mutants killed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
