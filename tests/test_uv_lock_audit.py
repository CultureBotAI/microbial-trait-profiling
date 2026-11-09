"""Keep dependency drift distinct from lock serialization metadata (#1286)."""

from __future__ import annotations

import hashlib
import json
import runpy
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from mtp.local_path_patterns import compile_absolute_local_paths
from mtp.provenance_paths import assert_portable_text, resolve_reference, verify_reference

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts/adhoc/issue1286_compare_uv_locks.py"
GUARD_PROBE_SCRIPT = REPO_ROOT / "scripts/adhoc/issue928_guard_isolation.py"
HELPER = runpy.run_path(str(SCRIPT))
GUARD_PROBE = runpy.run_path(str(GUARD_PROBE_SCRIPT))
ABSOLUTE_LOCAL_PATH = compile_absolute_local_paths(flags=0, classified=False)
LOCK = """version = 1
revision = 1
requires-python = ">=3.11"
resolution-markers = ["python_version < '3.12'", "python_version >= '3.12'"]
[[package]]
name = "example"
version = "1.0"
source = { registry = "https://example.org/simple" }
dependencies = [{name = "other", marker = "python_version < '3.12'"}]
sdist = {url = "https://example.org/source.tar.gz", hash = "sha256:original", size = 12}
wheels = [{url = "https://example.org/package.whl", hash = "sha256:wheel", size = 9}]
[[package]]
name = "example"
version = "2.0"
source = { registry = "https://example.org/simple" }
"""


def compare(tmp_path: Path, before: str, after: str) -> dict[str, Any]:
    paths = [tmp_path / "before.lock", tmp_path / "after.lock"]
    for path, content in zip(paths, (before, after), strict=True):
        path.write_text(content)
    result: dict[str, Any] = HELPER["compare_locks"](*paths)
    assert paths[0].read_text() == before
    assert paths[1].read_text() == after
    for role, path in zip(("before", "after"), paths, strict=True):
        raw = path.read_bytes()
        assert result[role]["path"] == path.name
        assert not Path(result[role]["path"]).is_absolute()
        assert result[role]["sha256"] == hashlib.sha256(raw).hexdigest()
        assert result[role]["size_bytes"] == len(raw)
    assert result.get("path_basis") == {"kind": "root-relative-v1", "root": "."}
    assert_portable_text(json.dumps(result))
    return result


def test_only_root_revision_and_artifact_upload_times_are_metadata(tmp_path: Path) -> None:
    after = LOCK.replace("revision = 1", "revision = 5")
    after = after.replace("size = 12}", 'size = 12, upload-time = "2026-01-01T00:00:00Z"}')
    after = after.replace("size = 9}", 'size = 9, upload-time = "2026-01-02T00:00:00Z"}')
    result = compare(tmp_path, LOCK, after)
    assert result["substantive_equal"] is True
    assert result["normalized_diff"] == ""
    assert result["package_version_changes"] == {}
    assert result["format_metadata"]["after"] == [
        {"path": ["revision"], "value": 5},
        {"path": ["package", 0, "sdist", "upload-time"], "value": "2026-01-01T00:00:00Z"},
        {"path": ["package", 0, "wheels", 0, "upload-time"], "value": "2026-01-02T00:00:00Z"},
    ]


@pytest.mark.parametrize(
    ("old", "new"),
    [
        ("example.org/simple", "different.example/simple"),
        ("sha256:original", "sha256:changed"),
        ("size = 12", "size = 13"),
        ("python_version < '3.12'", "python_version < '3.13'"),
        ("version = 1\n", "version = 2\n"),
        ('requires-python = ">=3.11"', 'requires-python = ">=3.12"'),
    ],
)
def test_sources_hashes_sizes_markers_and_other_fields_remain_substantive(
    tmp_path: Path, old: str, new: str
) -> None:
    result = compare(tmp_path, LOCK, LOCK.replace(old, new))
    assert result["substantive_equal"] is False
    assert result["normalized_diff"].startswith("--- before.normalized.json\n")


def test_metadata_named_fields_elsewhere_are_not_discarded(tmp_path: Path) -> None:
    before = LOCK.replace("source = { registry =", 'source = { revision = "first", registry =')
    after = before.replace('revision = "first"', 'revision = "second"')
    result = compare(tmp_path, before, after)
    assert result["substantive_equal"] is False
    before = LOCK + '\n[custom]\nupload-time = "first"\n'
    result = compare(
        tmp_path, before, before.replace('upload-time = "first"', 'upload-time = "next"')
    )
    assert result["substantive_equal"] is False


def test_array_and_package_order_are_not_normalized_away(tmp_path: Path) -> None:
    markers = "[\"python_version < '3.12'\", \"python_version >= '3.12'\"]"
    reversed_markers = "[\"python_version >= '3.12'\", \"python_version < '3.12'\"]"
    assert (
        compare(tmp_path, LOCK, LOCK.replace(markers, reversed_markers))["substantive_equal"]
        is False
    )
    header, first, second = LOCK.split("[[package]]")
    result = compare(tmp_path, LOCK, header + "[[package]]" + second + "[[package]]" + first)
    assert result["substantive_equal"] is False


def test_multiple_versions_are_reported_without_collapsing_a_package(tmp_path: Path) -> None:
    result = compare(tmp_path, LOCK, LOCK.replace('version = "2.0"', 'version = "2.1"'))
    assert result["package_version_changes"] == {
        "example": {"before": ["1.0", "2.0"], "after": ["1.0", "2.1"]}
    }


def test_cli_writes_fresh_report_and_never_overwrites(tmp_path: Path) -> None:
    before, after, output = (
        tmp_path / name for name in ("before.lock", "after.lock", "audit.json")
    )
    before.write_text(LOCK)
    after.write_text(LOCK)
    argv = ["--before", str(before), "--after", str(after), "--output", str(output)]
    assert HELPER["main"](argv) == 0
    saved = output.read_bytes()
    assert json.loads(saved)["substantive_equal"] is True
    with pytest.raises(FileExistsError) as exc_info:
        HELPER["main"](argv)
    assert str(output) not in str(exc_info.value)
    assert_portable_text(str(exc_info.value))
    assert output.read_bytes() == saved
    alias = tmp_path / "dangling.json"
    alias.symlink_to(tmp_path / "missing.json")
    with pytest.raises(FileExistsError) as symlink_exc:
        HELPER["main"]([*argv[:-1], str(alias)])
    assert str(alias) not in str(symlink_exc.value)
    assert not (tmp_path / "missing.json").exists()


def test_compare_locks_uses_portable_provenance_and_refuses_account_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#1290: lock comparator emits portable paths and refuses account-bearing paths."""
    home = tmp_path / "Users" / "residentacct"
    workspace = home / "workspace"
    workspace.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", lambda: home)
    monkeypatch.setenv("USER", "residentacct")
    monkeypatch.setenv("LOGNAME", "residentacct")

    before = workspace / "before.lock"
    after = workspace / "sub" / "after.lock"
    after.parent.mkdir()
    before.write_text(LOCK, encoding="utf-8")
    after.write_text(LOCK, encoding="utf-8")
    output = workspace / "reports" / "audit.json"
    output.parent.mkdir()

    assert (
        HELPER["main"](["--before", str(before), "--after", str(after), "--output", str(output)])
        == 0
    )
    raw_report = output.read_text(encoding="utf-8")
    assert ABSOLUTE_LOCAL_PATH.search(raw_report) is None
    assert "residentacct" not in raw_report.casefold()
    report = json.loads(raw_report)
    assert report.get("path_basis") == {"kind": "root-relative-v1", "root": ".."}
    assert report["before"]["path"] == "before.lock"
    assert report["after"]["path"] == "sub/after.lock"
    assert_portable_text(raw_report)
    assert (
        verify_reference(output, report["before"], path_basis=report["path_basis"])
        == before.resolve()
    )
    assert (
        verify_reference(output, report["after"], path_basis=report["path_basis"])
        == after.resolve()
    )
    doc = HELPER["compare_locks"].__doc__ or ""
    assert "output" in doc
    assert "path_basis" in doc
    assert "audit.json" in doc

    local_dep_lock = workspace / "local-dep.lock"
    local_dep_lock.write_text(
        LOCK
        + '\n[[package]]\nname = "local-pkg"\nversion = "0.1.0"\n'
        + f'source = {{ directory = "{home.as_posix()}/projects/local-pkg" }}\n',
        encoding="utf-8",
    )
    refused_output = workspace / "reports" / "refused.json"
    with pytest.raises(ValueError, match="Provenance contains a local path or local prefix"):
        HELPER["compare_locks"](before, local_dep_lock, output=refused_output)
    assert not refused_output.exists()

    with pytest.raises(ValueError, match="Provenance contains a local path or local prefix"):
        HELPER["main"](
            [
                "--before",
                str(before),
                "--after",
                str(local_dep_lock),
                "--output",
                str(refused_output),
            ]
        )
    assert not refused_output.exists()

    with monkeypatch.context() as mp:
        mp.setitem(
            HELPER["main"].__globals__,
            "compare_locks",
            lambda *_a, **_kw: {"leaked": f"{home.as_posix()}/unfiltered"},
        )
        with pytest.raises(ValueError, match="Provenance contains a local path or local prefix"):
            HELPER["main"](
                [
                    "--before",
                    str(before),
                    "--after",
                    str(after),
                    "--output",
                    str(refused_output),
                ]
            )
    assert not refused_output.exists()

    account_input = workspace / "residentacct-locks" / "before.lock"
    account_input.parent.mkdir()
    account_input.write_text(LOCK, encoding="utf-8")
    with pytest.raises(ValueError, match="Provenance contains a local account name"):
        HELPER["main"](
            [
                "--before",
                str(account_input),
                "--after",
                str(after),
                "--output",
                str(refused_output),
            ]
        )
    assert not refused_output.exists()

    with pytest.raises(ValueError, match=r"Provenance root must be below a home or temporary root"):
        HELPER["main"](
            [
                "--before",
                str(before),
                "--after",
                str(after),
                "--output",
                str(refused_output),
                "--provenance-root",
                str(home),
            ]
        )
    assert not refused_output.exists()


def test_guard_isolation_probe_emits_portable_paths_and_refuses_account_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#1027, #2054, #2055, #2056: guard probe portable paths, boundary scrub, and root basis."""
    home = tmp_path / "Users" / "probeacct"
    workspace = home / "workspace"
    workspace.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", lambda: home)
    monkeypatch.setenv("USER", "probeacct")
    monkeypatch.setenv("LOGNAME", "probeacct")

    # Unit checks for Probe temp_roots parent recording and scrubber branches (#2055, #2057)
    probe_unit = GUARD_PROBE["Probe"](SimpleNamespace())
    basetemp_leaf = workspace / "pytest-of-unit" / "pytest-0"
    basetemp_leaf.mkdir(parents=True)
    dummy_cfg = SimpleNamespace(
        _tmp_path_factory=SimpleNamespace(getbasetemp=lambda: basetemp_leaf)
    )
    probe_unit.pytest_configure(dummy_cfg)
    assert basetemp_leaf.resolve() in probe_unit.temp_roots
    assert basetemp_leaf.resolve().parent in probe_unit.temp_roots

    private_scrubbed = GUARD_PROBE["scrub_pytest_log"](
        "opened /opt/custom-private-test/run/artifact.tsv\n",
        output=Path("/private/opt/custom-private-test/run"),
        repo_root=REPO_ROOT,
        temp_roots=set(),
    )
    assert "<output-dir>/artifact.tsv" in private_scrubbed

    ellipsis_scrubbed = GUARD_PROBE["scrub_pytest_log"](
        'truncated = ".../custom-ellipsis-root/pytest-0/case/labels.tsv"\n',
        output=workspace / "guard-evidence",
        repo_root=REPO_ROOT,
        temp_roots={Path("/opt/custom-ellipsis-root/pytest-0")},
    )
    assert ".../<pytest-tempdir>/case/labels.tsv" in ellipsis_scrubbed, (
        "ellipsis fragment not scrubbed"
    )

    with monkeypatch.context() as mp:
        mp.chdir(workspace)
        rel_repo_scrubbed = GUARD_PROBE["scrub_pytest_log"](
            "out/tests/test_sample.py:10: in "
            "test_min_stratified_pairs_is_rejected_outside_site_stats\n"
            "output = isolated_investigation(infographic)\n"
            "Captured stdout call\n",
            output=Path("out"),
            repo_root=workspace / "out",
            temp_roots=set(),
        )
    assert "<repo-root>/tests/test_sample.py:10:" in rel_repo_scrubbed
    assert "test_min_stratified_pairs_is_rejected_outside_site_stats" in rel_repo_scrubbed
    assert "output = isolated_investigation(infographic)" in rel_repo_scrubbed
    assert "Captured stdout call" in rel_repo_scrubbed
    assert "<output-dir>side" not in rel_repo_scrubbed
    assert "std<output-dir>" not in rel_repo_scrubbed

    with monkeypatch.context() as mp:
        mp.chdir(REPO_ROOT / "scripts")
        dotdot_scrubbed = GUARD_PROBE["scrub_pytest_log"](
            "../outside/file.py .. [100%] ... passed\n",
            output=workspace / "guard-evidence",
            repo_root=REPO_ROOT,
            temp_roots=set(),
        )
    assert "../outside/file.py .. [100%] ... passed" in dotdot_scrubbed

    account_output = workspace / "probeacct-evidence"
    with pytest.raises(ValueError, match="Provenance contains a local account name"):
        GUARD_PROBE["main"](["--output-dir", str(account_output)])
    assert not account_output.exists()

    with pytest.raises(ValueError, match=r"Provenance root must be below a home or temporary root"):
        GUARD_PROBE["main"](
            ["--output-dir", str(home / "evidence"), "--provenance-root", str(home)]
        )
    assert not (home / "evidence").exists()

    with pytest.raises(ValueError, match="Source is outside the declared provenance root"):
        GUARD_PROBE["main"](
            ["--output-dir", str(workspace / "evidence"), "--provenance-root", str(workspace)]
        )
    assert not (workspace / "evidence").exists()

    with pytest.raises(ValueError, match="Provenance contains a local path or local prefix"):
        GUARD_PROBE["scrub_pytest_log"](
            f"Wrote: {home.as_posix()}/unreplaced/path.tsv\n",
            output=workspace / "guard-evidence",
            repo_root=REPO_ROOT,
            temp_roots=set(),
        )

    with monkeypatch.context() as mp:
        mp.setitem(
            GUARD_PROBE["main"].__globals__,
            "inventory",
            lambda *_args, **_kwargs: {f"{home.as_posix()}/leaked.tsv": "0" * 64},
        )
        mp.setattr(GUARD_PROBE["pytest"], "main", lambda *_args, **_kwargs: 1)
        leaked_output = workspace / "leaked-inventory"
        with pytest.raises(ValueError, match="Provenance contains a local path or local prefix"):
            GUARD_PROBE["main"](["--output-dir", str(leaked_output)])
        assert not leaked_output.exists()

    probe_output = workspace / "guard-evidence"
    assert GUARD_PROBE["main"](["--output-dir", str(probe_output)]) == 0
    with pytest.raises(FileExistsError) as exists_exc:
        GUARD_PROBE["main"](["--output-dir", str(probe_output)])
    assert str(probe_output) not in str(exists_exc.value)
    assert "probeacct" not in str(exists_exc.value).casefold()
    assert_portable_text(str(exists_exc.value))

    isolation_path = probe_output / "isolation.json"
    raw_isolation = isolation_path.read_text(encoding="utf-8")
    record = json.loads(raw_isolation)
    assert record["verified"] is True
    assert record.get("path_basis") == {"kind": "supplied-root-v1"}
    assert record.get("isolated_path_basis") == {"kind": "root-relative-v1", "root": "."}
    assert record["canonical_origin"] == "scripts/create_trait_profile_infographic.py"
    assert record["isolated_origin"] == "infographic_without_guard.py"
    assert record["checkout_outputs"] == "trait-visualization/trait-profile-investigation"
    assert not Path(record["canonical_origin"]).is_absolute()
    assert not Path(record["isolated_origin"]).is_absolute()
    assert ABSOLUTE_LOCAL_PATH.search(raw_isolation) is None
    assert "probeacct" not in raw_isolation.casefold()
    assert str(REPO_ROOT) not in raw_isolation
    assert str(probe_output) not in raw_isolation
    assert_portable_text(raw_isolation)
    assert set(record["source_before"]) == {
        "scripts/create_trait_profile_infographic.py",
        "tests/test_trait_profile_infographic.py",
        "scripts/adhoc/issue928_guard_isolation.py",
    }
    assert record["source_before"] == record["source_after"]
    for ref in (
        record["canonical_origin"],
        record["checkout_outputs"],
        *record["before"],
        *record["after"],
        *record["source_before"],
        *record["source_after"],
    ):
        resolved = resolve_reference(
            isolation_path,
            ref,
            path_basis=record["path_basis"],
            provenance_root=REPO_ROOT,
        )
        assert resolved == (REPO_ROOT / ref).resolve()
    assert (
        resolve_reference(
            isolation_path,
            record["isolated_origin"],
            path_basis=record["isolated_path_basis"],
        )
        == (probe_output / "infographic_without_guard.py").resolve()
    )
    raw_log = (probe_output / "pytest.log").read_text(encoding="utf-8")
    assert "DID NOT RAISE SystemExit" in raw_log
    assert "<pytest-tempdir>" in raw_log
    assert "test_min_stratified_pairs_is_rejected_outside_site_stats" in raw_log
    assert ABSOLUTE_LOCAL_PATH.search(raw_log) is None
    assert "probeacct" not in raw_log.casefold()
    assert str(REPO_ROOT) not in raw_log
    assert str(probe_output) not in raw_log
    assert_portable_text(raw_log)

    # Verify covering --provenance-root above checkout (#2054), relative --output-dir out,
    # and running from inside the checkout (#2055).
    covering = workspace / "covering"
    nested_repo = covering / "repo"
    for rel in (
        "scripts/create_trait_profile_infographic.py",
        "tests/test_trait_profile_infographic.py",
        "scripts/adhoc/issue928_guard_isolation.py",
    ):
        dst = nested_repo / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / rel, dst)
    nested_inv = nested_repo / "trait-visualization/trait-profile-investigation"
    nested_inv.mkdir(parents=True, exist_ok=True)
    sample_inv_file = nested_inv / "sample_output.tsv"
    sample_inv_file.write_text("trait\tcount\n", encoding="utf-8")

    basetemp_sim = home / "pytest-of-probeacct" / "pytest-0"
    basetemp_sim.mkdir(parents=True, exist_ok=True)

    def _fast_pytest_main(_args: list[str], plugins: list[Any] | None = None) -> int:
        active = (plugins or [])[0]
        active.temp_roots.update((basetemp_sim.resolve(), basetemp_sim.resolve().parent))
        active.failed.extend(["call_a", "call_b"])
        sys.stdout.write(
            f"{active.temp_roots and ''}"
            "tests/test_trait_profile_infographic.py::"
            "test_min_stratified_pairs_is_rejected_outside_site_stats[mode0]\n"
            "output = isolated_investigation(infographic)\n"
            "path = '.../probeacct/pytest-of-probeacct/pytest-0"
            "/test_min_stratified_pairs_is_r0/labels.tsv'\n"
            "E   Failed: DID NOT RAISE <class 'SystemExit'>\n"
            "Captured stdout call\n"
        )
        return 1

    with monkeypatch.context() as mp:
        mp.setitem(GUARD_PROBE["main"].__globals__, "ROOT", nested_repo.resolve())
        mp.setitem(GUARD_PROBE["main"].__globals__, "REPO_ROOT", nested_repo.resolve())
        mp.setitem(
            GUARD_PROBE["main"].__globals__,
            "__file__",
            str((nested_repo / "scripts/adhoc/issue928_guard_isolation.py").resolve()),
        )
        mp.setitem(
            GUARD_PROBE["main"].__globals__,
            "PRODUCER",
            (nested_repo / "scripts/create_trait_profile_infographic.py").resolve(),
        )
        mp.setitem(
            GUARD_PROBE["main"].__globals__,
            "TEST",
            (nested_repo / "tests/test_trait_profile_infographic.py").resolve(),
        )
        mp.setattr(GUARD_PROBE["pytest"], "main", _fast_pytest_main)
        orig_load = GUARD_PROBE["load"]

        def _load_nested(path: Path, name: str) -> Any:
            mod = orig_load(path, name)
            if name == "issue950_canonical_infographic":
                vars(mod)["INVESTIGATION_DIR"] = nested_inv.resolve()
            return mod

        mp.setitem(GUARD_PROBE["main"].__globals__, "load", _load_nested)
        cov_out = covering / "runs" / "out-covering"
        cov_out.parent.mkdir(parents=True, exist_ok=True)
        assert (
            GUARD_PROBE["main"](["--output-dir", str(cov_out), "--provenance-root", str(covering)])
            == 0
        )
        cov_iso_path = cov_out / "isolation.json"
        cov_rec = json.loads(cov_iso_path.read_text(encoding="utf-8"))
        assert cov_rec["checkout_outputs"] == "repo/trait-visualization/trait-profile-investigation"
        assert set(cov_rec["before"]) == {
            "repo/trait-visualization/trait-profile-investigation/sample_output.tsv"
        }
        for ref in (
            cov_rec["canonical_origin"],
            cov_rec["isolated_origin"],
            cov_rec["checkout_outputs"],
            *cov_rec["before"],
            *cov_rec["after"],
            *cov_rec["source_before"],
            *cov_rec["source_after"],
        ):
            resolved_cov = resolve_reference(
                cov_iso_path,
                ref,
                path_basis=cov_rec["path_basis"],
            )
            assert resolved_cov.exists(), f"{ref} resolved to missing {resolved_cov}"

        mp.chdir(covering / "runs")
        assert GUARD_PROBE["main"](["--output-dir", "out"]) == 0
        rel_out_log = (covering / "runs" / "out" / "pytest.log").read_text(encoding="utf-8")
        assert "test_min_stratified_pairs_is_rejected_outside_site_stats[mode0]" in rel_out_log
        assert "output = isolated_investigation(infographic)" in rel_out_log
        assert "Captured stdout call" in rel_out_log
        assert "<output-dir>side" not in rel_out_log
        assert "std<output-dir>" not in rel_out_log

        mp.chdir(nested_repo / "scripts")
        inside_out = covering / "runs" / "from-inside-scripts"
        assert GUARD_PROBE["main"](["--output-dir", str(inside_out)]) == 0
        inside_log = (inside_out / "pytest.log").read_text(encoding="utf-8")
        assert ".../<pytest-tempdir>/test_min_stratified_pairs_is_r0/labels.tsv" in inside_log
        assert "probeacct" not in inside_log.casefold()
