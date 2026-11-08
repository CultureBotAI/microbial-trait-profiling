"""Containment follows filesystem identity without guessing across devices (#1178).

Controlled probes use real pytest directories and override only Path.stat results;
resolution is not mocked. Native macOS probes are conditional on an actual alias.
"""

from __future__ import annotations

import errno
import os
import stat
import sys
from collections.abc import Mapping
from pathlib import Path
from types import SimpleNamespace
from typing import cast

import pytest

from mtp.workflow import output_roots


def _identity(value: os.stat_result, *, device: int, inode: int) -> os.stat_result:
    fields = list(value)
    fields[stat.ST_DEV] = device
    fields[stat.ST_INO] = inode
    return os.stat_result(fields)


def _patch_stat(
    monkeypatch: pytest.MonkeyPatch,
    overrides: Mapping[Path, os.stat_result | OSError],
) -> list[Path]:
    original = Path.stat
    observed: list[Path] = []

    def controlled(path: Path, *, follow_symlinks: bool = True) -> os.stat_result:
        observed.append(path)
        if path in overrides:
            result = overrides[path]
            if isinstance(result, OSError):
                raise result
            return result
        return original(path, follow_symlinks=follow_symlinks)

    monkeypatch.setattr(Path, "stat", controlled)
    return observed


def _directories(tmp_path: Path) -> tuple[Path, Path]:
    protected = tmp_path / "protected"
    candidate = tmp_path / "other-spelling"
    protected.mkdir()
    candidate.mkdir()
    return protected, candidate


@pytest.mark.parametrize("device", [0, 37])
def test_disjoint_identity_alias_contains_a_missing_multicomponent_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, device: int
) -> None:
    protected, alias = _directories(tmp_path)
    root_stat = _identity(protected.stat(), device=device, inode=4242)
    alias_stat = _identity(alias.stat(), device=device, inode=4242)
    assert not alias.resolve().is_relative_to(protected.resolve())
    _patch_stat(monkeypatch, {protected: root_stat, alias: alias_stat})
    output = alias / "new" / "nested" / "output"

    assert output_roots.temporary_root(output, [protected]) == protected
    assert not output.exists()


def test_identity_search_climbs_past_a_nearer_existing_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    protected, alias = _directories(tmp_path)
    nearer = alias / "existing"
    nearer.mkdir()
    root_stat = protected.stat()
    _patch_stat(
        monkeypatch,
        {alias: _identity(alias.stat(), device=root_stat.st_dev, inode=root_stat.st_ino)},
    )
    output = nearer / "missing" / "new" / "output"

    assert output_roots.temporary_root(output, [protected]) == protected
    assert not output.exists()


def test_same_inode_on_a_different_device_is_not_containment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    protected, outside = _directories(tmp_path)
    root_stat = protected.stat()
    _patch_stat(
        monkeypatch,
        {outside: _identity(outside.stat(), device=root_stat.st_dev + 1, inode=root_stat.st_ino)},
    )

    assert output_roots.temporary_root(outside / "missing" / "output", [protected]) is None


@pytest.mark.parametrize("same_identity", [False, True])
def test_casefold_equal_spellings_are_decided_by_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, same_identity: bool
) -> None:
    protected = tmp_path / "CaseSensitiveRoot"
    spelling = tmp_path / "casesensitiveroot"
    protected.mkdir()
    spelling.mkdir(exist_ok=True)
    assert str(protected).casefold() == str(spelling).casefold()
    _patch_stat(
        monkeypatch,
        {
            protected: _identity(protected.stat(), device=37, inode=4242),
            spelling: _identity(spelling.stat(), device=37, inode=4242 if same_identity else 4243),
        },
    )

    result = output_roots.temporary_root(spelling / "missing" / "output", [protected])
    assert result == (protected if same_identity else None)


def test_missing_protected_root_retains_its_lexical_boundary(tmp_path: Path) -> None:
    missing = tmp_path / "missing-protected"
    output = missing / "new" / "nested" / "output"

    assert output_roots.temporary_root(output, [missing]) == missing
    assert not missing.exists()


def test_missing_protected_root_does_not_protect_its_existing_parent(tmp_path: Path) -> None:
    missing = tmp_path / "missing-protected"
    outside = tmp_path / "sibling-output" / "new"

    assert output_roots.temporary_root(outside, [missing]) is None
    assert not missing.exists()
    assert not outside.exists()


@pytest.mark.parametrize("error_number", [errno.EACCES, errno.EIO, errno.ENOTDIR])
@pytest.mark.parametrize("allowed", [False, True])
def test_candidate_probe_errors_fail_closed_even_with_temporary_opt_in(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, error_number: int, allowed: bool
) -> None:
    protected, candidate = _directories(tmp_path)
    monkeypatch.setattr(output_roots, "system_temporary_directories", lambda: [protected])
    _patch_stat(monkeypatch, {candidate: OSError(error_number, "controlled identity probe error")})
    output = candidate / "missing" / "nested" / "output"

    with pytest.raises(ValueError, match=r"(?i)(output|identity|contain)"):
        output_roots.require_durable_output_path(output, allowed=allowed)
    assert not output.exists()


@pytest.mark.parametrize("error_number", [errno.EACCES, errno.EIO, errno.ENOTDIR])
def test_protected_root_probe_errors_are_not_silently_treated_as_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, error_number: int
) -> None:
    protected, outside = _directories(tmp_path)
    _patch_stat(monkeypatch, {protected: OSError(error_number, "controlled root probe error")})

    with pytest.raises(ValueError, match=r"(?i)(output|identity|contain)"):
        output_roots.temporary_root(outside / "new" / "output", [protected])


def _invalid_stat(original: os.stat_result, field: str, value: object) -> os.stat_result:
    values: dict[str, object] = {
        "st_mode": original.st_mode,
        "st_dev": original.st_dev,
        "st_ino": original.st_ino,
    }
    values[field] = value
    # Intentionally model a filesystem adapter that supplies unusable metadata.
    return cast(os.stat_result, SimpleNamespace(**values))


@pytest.mark.parametrize(
    ("field", "value"),
    [("st_ino", 0), ("st_ino", -1), ("st_dev", -1), ("st_ino", True), ("st_dev", False)],
)
@pytest.mark.parametrize("allowed", [False, True])
def test_invalid_candidate_identity_is_refused_even_with_temporary_opt_in(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    field: str,
    value: object,
    allowed: bool,
) -> None:
    protected, candidate = _directories(tmp_path)
    monkeypatch.setattr(output_roots, "system_temporary_directories", lambda: [protected])
    invalid = _invalid_stat(candidate.stat(), field, value)
    _patch_stat(monkeypatch, {candidate: invalid})

    with pytest.raises(ValueError, match=r"(?i)(output|identity|contain)"):
        output_roots.require_durable_output_path(candidate / "missing" / "output", allowed=allowed)


@pytest.mark.parametrize("missing_field", ["st_dev", "st_ino"])
@pytest.mark.parametrize("allowed", [False, True])
def test_missing_candidate_identity_metadata_is_a_clear_refusal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, missing_field: str, allowed: bool
) -> None:
    protected, candidate = _directories(tmp_path)
    monkeypatch.setattr(output_roots, "system_temporary_directories", lambda: [protected])
    actual = candidate.stat()
    fields = {"st_mode": actual.st_mode, "st_dev": actual.st_dev, "st_ino": actual.st_ino}
    del fields[missing_field]
    incomplete = cast(os.stat_result, SimpleNamespace(**fields))
    _patch_stat(monkeypatch, {candidate: incomplete})

    with pytest.raises(ValueError, match=r"(?i)(output|identity|contain)"):
        output_roots.require_durable_output_path(candidate / "missing" / "output", allowed=allowed)


@pytest.mark.parametrize(("field", "value"), [("st_ino", 0), ("st_dev", False)])
def test_invalid_protected_root_identity_does_not_allow_an_unproven_outside_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, field: str, value: object
) -> None:
    protected, outside = _directories(tmp_path)
    invalid = _invalid_stat(protected.stat(), field, value)
    _patch_stat(monkeypatch, {protected: invalid})

    with pytest.raises(ValueError, match=r"(?i)(output|identity|contain)"):
        output_roots.temporary_root(outside / "missing" / "output", [protected])


@pytest.mark.parametrize("error_number", [errno.EACCES, errno.EIO, errno.ENOTDIR])
@pytest.mark.parametrize("failure_stage", ["initial", "after_normalization"])
def test_strict_resolution_errors_are_not_treated_as_missing_suffixes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    error_number: int,
    failure_stage: str,
) -> None:
    output = tmp_path / "missing" / ".." / "nested" / "output"
    normalized = tmp_path / "nested" / "output"
    original = Path.resolve
    probes: list[tuple[Path, bool]] = []

    def controlled(path: Path, strict: bool = False) -> Path:
        if path == output:
            probes.append((path, strict))
            if strict:
                number = error_number if failure_stage == "initial" else errno.ENOENT
                raise OSError(number, "controlled first strict resolver probe")
            return normalized
        if path == normalized:
            probes.append((path, strict))
            if strict and failure_stage == "after_normalization":
                raise OSError(error_number, "controlled normalized strict resolver probe")
            return normalized
        return original(path, strict=strict)

    monkeypatch.setattr(Path, "resolve", controlled)
    with pytest.raises(ValueError, match=r"(?i)(output|identity|contain)") as refused:
        output_roots._resolved_path(output)
    assert isinstance(refused.value.__cause__, OSError)
    assert refused.value.__cause__.errno == error_number
    expected = [(output, True)]
    if failure_stage == "after_normalization":
        expected.extend([(output, False), (normalized, True)])
    assert probes == expected


@pytest.mark.parametrize("mode_case", ["missing", "none", "text", "float", "boolean", "file"])
def test_unusable_protected_root_mode_is_a_clear_refusal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode_case: str
) -> None:
    protected, outside = _directories(tmp_path)
    original = protected.stat()
    fields: dict[str, object] = {"st_dev": original.st_dev, "st_ino": original.st_ino}
    modes: dict[str, object] = {
        "none": None,
        "text": "directory",
        "float": 3.5,
        "boolean": True,
        "file": stat.S_IFREG | 0o644,
    }
    if mode_case != "missing":
        fields["st_mode"] = modes[mode_case]
    malformed = cast(os.stat_result, SimpleNamespace(**fields))
    _patch_stat(monkeypatch, {protected: malformed})

    with pytest.raises(ValueError, match=r"(?i)(output|identity|contain)"):
        output_roots.temporary_root(outside / "missing" / "output", [protected])


def test_shared_containment_api_is_available_without_temporary_root_defaults(
    tmp_path: Path,
) -> None:
    # This is an API-availability check, not one of the semantic red witnesses.
    containing_root = getattr(output_roots, "containing_root", None)
    assert callable(containing_root)
    protected, outside = _directories(tmp_path)

    assert containing_root(protected / "new" / "output", [protected]) == protected
    assert containing_root(outside / "new" / "output", [protected]) is None


@pytest.mark.skipif(sys.platform != "darwin", reason="Native macOS filesystem alias probe")
def test_native_case_alias_is_refused_when_the_scratch_filesystem_supports_it(
    tmp_path: Path,
) -> None:
    protected = tmp_path / "NativeCaseRoot"
    protected.mkdir()
    alias = tmp_path / "nativecaseroot"
    try:
        matching = alias.samefile(protected)
    except FileNotFoundError:
        pytest.skip("Scratch filesystem is case-sensitive")
    if not matching:
        pytest.skip("Scratch filesystem does not expose this case alias")
    protected_stat, alias_stat = protected.stat(), alias.stat()
    assert (protected_stat.st_dev, protected_stat.st_ino) == (alias_stat.st_dev, alias_stat.st_ino)

    assert output_roots.temporary_root(alias / "fresh" / "output", [protected]) is not None


@pytest.mark.skipif(sys.platform != "darwin", reason="Native macOS filesystem alias probe")
def test_native_firmlink_alias_is_refused_when_scratch_has_an_identity_alias(
    tmp_path: Path,
) -> None:
    protected = tmp_path.resolve() / "native-firmlink-root"
    protected.mkdir()
    alias = Path("/System/Volumes/Data") / protected.relative_to(Path("/"))
    try:
        matching = alias.samefile(protected)
    except (FileNotFoundError, PermissionError):
        pytest.skip("Scratch filesystem does not expose an accessible Data-volume alias")
    if not matching:
        pytest.skip("The Data-volume spelling is not the same scratch directory")
    protected_stat, alias_stat = protected.stat(), alias.stat()
    assert (protected_stat.st_dev, protected_stat.st_ino) == (alias_stat.st_dev, alias_stat.st_ino)

    assert output_roots.temporary_root(alias / "fresh" / "output", [protected]) is not None
