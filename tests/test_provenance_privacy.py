"""Private root and account-name publication boundaries (#1765/#1767/#1768)."""

from __future__ import annotations

import json
import pwd
import tempfile
from pathlib import Path
from types import ModuleType

import pytest

from mtp.paths import file_sha256
from test_provenance_paths import api as api


@pytest.fixture
def private_environment(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """Keep all identity authorities distinct so each collection source is pinned."""
    home = tmp_path / "homes" / "homeresident"
    passwd_home = tmp_path / "homes" / "pwdresident"
    temporary = tmp_path / "temporary-boundary"
    for directory in (home, passwd_home, temporary):
        directory.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", lambda: home)
    monkeypatch.setattr(tempfile, "gettempdir", lambda: str(temporary))
    monkeypatch.setattr(
        pwd,
        "getpwuid",
        lambda uid: pwd.struct_passwd(("systemperson", "x", uid, 1, "", str(passwd_home), "shell")),
    )
    for variable, value in {
        "LOGNAME": "loginperson",
        "USER": "envperson",
        "LNAME": "aliasperson",
        "USERNAME": "windowsperson",
        "TMPDIR": str(temporary),
        "TEMP": str(temporary),
        "TMP": str(temporary),
    }.items():
        monkeypatch.setenv(variable, value)
    return {"home": home, "passwd_home": passwd_home, "temporary": temporary}


@pytest.mark.parametrize("boundary", ["home", "passwd_home", "temporary"])
@pytest.mark.parametrize("ancestor", [False, True])
@pytest.mark.parametrize("reader", [False, True])
def test_private_root_or_ancestor_is_refused_even_without_sources(
    tmp_path: Path,
    private_environment: dict[str, Path],
    api: ModuleType,
    boundary: str,
    ancestor: bool,
    reader: bool,
) -> None:
    root = private_environment[boundary]
    if ancestor:
        root = root.parent
    record = tmp_path / "external" / "record.json"

    def exercise() -> None:
        if reader:
            api.resolve_provenance_root(
                record, path_basis={"kind": "supplied-root-v1"}, provenance_root=root
            )
        else:
            api.PortablePaths.create(record, [], root=root)

    with pytest.raises(ValueError, match=r"home|temp|private"):
        exercise()


@pytest.mark.parametrize("spelling", ["symlink", "dotdot"])
def test_private_root_check_uses_resolved_location(
    tmp_path: Path, private_environment: dict[str, Path], api: ModuleType, spelling: str
) -> None:
    boundary = private_environment["temporary"]
    child = boundary / "bundle"
    child.mkdir()
    if spelling == "symlink":
        root = tmp_path / "alias"
        root.symlink_to(boundary, target_is_directory=True)
    else:
        root = child / ".."
    with pytest.raises(ValueError, match=r"home|temp|private"):
        api.PortablePaths.create(tmp_path / "record.json", [], root=root)


@pytest.mark.parametrize(
    "account",
    [
        "loginperson",
        "envperson",
        "aliasperson",
        "windowsperson",
        "systemperson",
        "homeresident",
        "pwdresident",
    ],
)
def test_every_identity_authority_is_scanned_in_serialized_metadata(
    private_environment: dict[str, Path], api: ModuleType, account: str
) -> None:
    assert private_environment
    text = json.dumps({"detail": "pytest-of-" + account.upper()})
    with pytest.raises(ValueError, match="account"):
        api.assert_portable_text(text)


@pytest.mark.parametrize(
    "name", ["envperson", "-Users-envperson-work", "pytest-of-envperson", "prefixENVPERSONsuffix"]
)
@pytest.mark.parametrize("operation", ["create", "encode", "resolve"])
def test_account_names_in_relative_references_are_refused(
    private_environment: dict[str, Path],
    api: ModuleType,
    name: str,
    operation: str,
) -> None:
    root = private_environment["temporary"] / "bundle"
    root.mkdir()
    source = root / name / "source.tsv"
    source.parent.mkdir()
    source.write_text("source bytes\n")
    record = root / "record.json"

    def exercise() -> None:
        if operation == "create":
            api.PortablePaths.create(record, [source], root=root)
        elif operation == "encode":
            api.PortablePaths(record, root).encode(source)
        else:
            api.resolve_reference(
                record,
                source.relative_to(root).as_posix(),
                path_basis={"kind": "root-relative-v1", "root": "."},
            )

    with pytest.raises(ValueError, match="account"):
        exercise()


def test_root_below_account_named_scratch_directory_can_publish_clean_references(
    tmp_path: Path, private_environment: dict[str, Path], api: ModuleType
) -> None:
    root = private_environment["temporary"] / "-Users-envperson-project" / "bundle"
    root.mkdir(parents=True)
    source = root / "source.tsv"
    source.write_text("bound bytes\n")
    record = tmp_path / "elsewhere" / "record.json"
    paths = api.PortablePaths.create(record, [source], root=root)
    reference = {"path": paths.encode(source), "sha256": file_sha256(source)}
    payload = {"path_basis": paths.declaration, "source": reference}
    assert reference["path"] == "source.tsv"
    assert paths.declaration == {"kind": "supplied-root-v1"}
    api.assert_portable_text(json.dumps(payload))
    assert (
        api.verify_reference(record, reference, path_basis=paths.declaration, provenance_root=root)
        == source
    )


@pytest.mark.parametrize("kind", ["missing", "file"])
@pytest.mark.parametrize("reader", [False, True])
def test_invalid_explicit_root_reports_directory_problem_first(
    tmp_path: Path, api: ModuleType, kind: str, reader: bool
) -> None:
    root = tmp_path / "candidate"
    if kind == "file":
        root.write_text("not a directory\n")
    expected = "does not exist" if kind == "missing" else "not a directory"

    def exercise() -> None:
        if reader:
            api.resolve_provenance_root(
                tmp_path / "record.json",
                path_basis={"kind": "supplied-root-v1"},
                provenance_root=root,
            )
        else:
            api.PortablePaths.create(
                tmp_path / "record.json", [tmp_path / "outside.tsv"], root=root
            )

    with pytest.raises(ValueError, match=expected):
        exercise()


def test_supplied_root_declaration_extra_key_is_refused(tmp_path: Path, api: ModuleType) -> None:
    """N14: a valid, empty source tree isolates the exact declaration shape check."""
    with pytest.raises(ValueError, match="supplied-root declaration"):
        api.resolve_provenance_root(
            tmp_path / "record.json",
            path_basis={"kind": "supplied-root-v1", "extra": "untrusted"},
            provenance_root=tmp_path,
        )


@pytest.mark.parametrize("prefix", ["Users", "home", "uSeRs", "HOME"])
def test_other_home_prefix_is_refused_without_present_sources(api: ModuleType, prefix: str) -> None:
    """Root-relative empty records cannot weaken another account's home boundary."""
    root = Path(Path.cwd().anchor) / prefix / "otherresident"
    with pytest.raises(ValueError, match=r"home|private|prefix"):
        api.resolve_provenance_root(
            root / "record.json", path_basis={"kind": "root-relative-v1", "root": "."}
        )


@pytest.mark.parametrize("boundary", ["home", "temporary"])
def test_case_variant_private_boundary_is_refused(
    private_environment: dict[str, Path], api: ModuleType, boundary: str
) -> None:
    """Case-insensitive filesystems must not accept a differently spelled boundary."""
    root = private_environment[boundary]
    alias = root.with_name(root.name.upper())
    with pytest.raises(ValueError, match=r"home|temp|private"):
        api.resolve_provenance_root(
            alias / "record.json", path_basis={"kind": "root-relative-v1", "root": "."}
        )
