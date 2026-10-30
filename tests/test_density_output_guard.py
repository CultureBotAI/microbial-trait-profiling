"""Protect six-trait density outputs before declared-list I/O (#1093)."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import NoReturn

import pytest

from mtp import density_adjustment as density
from mtp.trait_lists import DeclaredTraitList

_TABLE = "table1_density_adjusted_effects.tsv"
_DECLARED = DeclaredTraitList(("METPO:1000608",), "traits.txt", "a" * 64)


def _forbidden_io(*_args: object, **_kwargs: object) -> NoReturn:
    raise AssertionError("Refused trait-list destination reached mkdir or input reading")


def _invoke(
    entrypoint: str, output: Path, monkeypatch: pytest.MonkeyPatch, *, declared: bool = True
) -> None:
    if entrypoint == "api":
        density.run_density_audit(output_dir=output, trait_list=_DECLARED if declared else None)
    else:
        arguments = ["density_adjustment", "--output-dir", str(output)]
        if declared:
            arguments.extend(["--trait-list", str(output.parent / "missing-traits.txt")])
        monkeypatch.setattr(sys, "argv", arguments)
        density.main()


def _protected_directory(tmp_path: Path, case: str, monkeypatch: pytest.MonkeyPatch) -> Path:
    if case == "arm":
        return tmp_path / "noel" / "density_adjustment"
    if case == "nested":
        return tmp_path / "jed" / "density_adjustment" / "declared" / "nested"
    if case == "alias_into_arm":
        target = tmp_path / "ilyas" / "density_adjustment"
        target.mkdir(parents=True)
        alias = tmp_path / "innocent_alias"
        alias.symlink_to(target, target_is_directory=True)
        return alias / "nested"
    if case in {"arm_alias_outward", "chained_alias_outward"}:
        target = tmp_path / "external_storage"
        target.mkdir()
        alias = tmp_path / "density_adjustment"
        alias.symlink_to(target, target_is_directory=True)
        if case == "chained_alias_outward":
            indirect = tmp_path / "innocent_alias"
            indirect.symlink_to(alias, target_is_directory=True)
            return indirect / "nested"
        return alias / "nested"
    default = tmp_path / "legacy_default_with_another_name"
    monkeypatch.setattr(density, "DEFAULT_OUTPUT_DIR", default)
    return default / "nested"


@pytest.mark.parametrize("entrypoint", ["api", "cli"])
@pytest.mark.parametrize(
    "case",
    ["arm", "nested", "alias_into_arm", "arm_alias_outward", "chained_alias_outward", "default"],
)
def test_declared_list_refuses_arm_and_default_paths_before_any_io(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    case: str,
    entrypoint: str,
) -> None:
    output = _protected_directory(tmp_path, case, monkeypatch)
    assert not output.exists()
    monkeypatch.setattr(density, "load_non_media_profile_matrix", _forbidden_io)
    monkeypatch.setattr(density, "load_trait_list", _forbidden_io)
    monkeypatch.setattr(Path, "mkdir", _forbidden_io)
    with pytest.raises(ValueError, match="own --output-dir"):
        _invoke(entrypoint, output, monkeypatch)
    assert not output.exists()


@pytest.mark.parametrize("entrypoint", ["api", "cli"])
@pytest.mark.parametrize("spelling", ["Density_Adjustment", "DENSITY_ADJUSTMENT"])
@pytest.mark.parametrize("location", ["arm", "default", "chained_alias_outward"])
def test_declared_list_refuses_case_variants_before_any_io(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    entrypoint: str,
    spelling: str,
    location: str,
) -> None:
    """Reserve case variants portably and exercise real aliases on insensitive volumes."""
    if location == "chained_alias_outward":
        target = tmp_path / "external_storage"
        target.mkdir()
        variant = tmp_path / spelling
        variant.symlink_to(target, target_is_directory=True)
        indirect = tmp_path / "innocent_alias"
        indirect.symlink_to(variant, target_is_directory=True)
        output = indirect / "declared" / "nested"
    else:
        canonical = tmp_path / "arm" / "density_adjustment"
        canonical.mkdir(parents=True)
        variant = canonical.parent / spelling
        if variant.exists():
            # On case-insensitive APFS this is the real reserved directory.
            assert variant.samefile(canonical)
        if location == "default":
            monkeypatch.setattr(density, "DEFAULT_OUTPUT_DIR", canonical)
        output = variant / "declared" / "nested"
    assert not output.exists()
    monkeypatch.setattr(density, "load_non_media_profile_matrix", _forbidden_io)
    monkeypatch.setattr(density, "load_trait_list", _forbidden_io)
    monkeypatch.setattr(Path, "mkdir", _forbidden_io)
    with pytest.raises(ValueError, match="own --output-dir"):
        _invoke(entrypoint, output, monkeypatch)
    assert not output.exists()


@pytest.mark.parametrize("entrypoint", ["api", "cli"])
@pytest.mark.parametrize(
    "prior", ["six_traits", "malformed", "empty", "symlink", "dangling", "directory"]
)
def test_declared_list_reserves_existing_density_table_without_parsing_or_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    prior: str,
    entrypoint: str,
) -> None:
    output = tmp_path / "ordinary_results"
    output.mkdir()
    table = output / _TABLE
    content = b"trait_id\traw_mean_diff\n" + b"".join(
        f"{trait}\t0.25\n".encode() for trait in density.TABLE_1_TRAITS
    )
    target = tmp_path / "existing-table.tsv"
    if prior == "directory":
        table.mkdir()
    elif prior in {"symlink", "dangling"}:
        if prior == "symlink":
            target.write_bytes(content)
        table.symlink_to(target)
    else:
        table.write_bytes(
            content
            if prior == "six_traits"
            else b"incomplete table\n"
            if prior == "malformed"
            else b""
        )
    before = table.read_bytes() if table.is_file() else None
    monkeypatch.setattr(density, "load_non_media_profile_matrix", _forbidden_io)
    monkeypatch.setattr(density, "load_trait_list", _forbidden_io)
    monkeypatch.setattr(Path, "mkdir", _forbidden_io)
    with pytest.raises(ValueError, match=r"fresh.*--output-dir"):
        _invoke(entrypoint, output, monkeypatch)
    if before is not None:
        assert table.read_bytes() == before
    elif prior == "directory":
        assert table.is_dir()
    else:
        assert table.is_symlink()
        assert not target.exists()
    assert tuple(output.iterdir()) == (table,)


class _InputReadReached(Exception):
    pass


def _stop_at_input_read(*_args: object, **_kwargs: object) -> NoReturn:
    raise _InputReadReached


@pytest.mark.parametrize("entrypoint", ["api", "cli"])
def test_default_six_trait_execution_keeps_existing_directory_and_table_behavior(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    entrypoint: str,
) -> None:
    output = tmp_path / "arm" / "density_adjustment"
    output.mkdir(parents=True)
    table = output / _TABLE
    original = b"unchanged six-trait output\n"
    table.write_bytes(original)
    monkeypatch.setattr(density, "load_non_media_profile_matrix", _stop_at_input_read)
    monkeypatch.setattr(density, "load_trait_list", _forbidden_io)
    with pytest.raises(_InputReadReached):
        _invoke(entrypoint, output, monkeypatch, declared=False)
    assert table.read_bytes() == original


@pytest.mark.parametrize("entrypoint", ["api", "cli"])
@pytest.mark.parametrize("name", ["declared_results", "density_adjustment_extra"])
def test_fresh_declared_destination_reaches_input_reading(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    entrypoint: str,
    name: str,
) -> None:
    output = tmp_path / name
    monkeypatch.setattr(density, "load_non_media_profile_matrix", _stop_at_input_read)
    monkeypatch.setattr(density, "load_trait_list", lambda _path: _DECLARED)
    with pytest.raises(_InputReadReached):
        _invoke(entrypoint, output, monkeypatch)
    assert output.is_dir()
