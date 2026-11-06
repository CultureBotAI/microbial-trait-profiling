"""Private root and account-name publication boundaries (#1765/#1767/#1768)."""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from types import ModuleType, SimpleNamespace

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
    password_module = ModuleType("pwd")
    monkeypatch.setattr(
        password_module,
        "getpwuid",
        lambda _uid: SimpleNamespace(pw_name="systemperson", pw_dir=str(passwd_home)),
        raising=False,
    )
    monkeypatch.setitem(sys.modules, "pwd", password_module)
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

    with pytest.raises(
        ValueError, match=r"^Provenance root must be below a home or temporary root\."
    ) as error:
        exercise()
    assert str(root) not in str(error.value)


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
    with pytest.raises(
        ValueError, match=r"^Provenance root must be below a home or temporary root\."
    ) as error:
        api.PortablePaths.create(tmp_path / "record.json", [], root=root)
    assert str(root) not in str(error.value)


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
def test_every_identity_authority_is_scanned_in_encoded_paths(
    private_environment: dict[str, Path], api: ModuleType, account: str
) -> None:
    assert private_environment
    root = private_environment["temporary"] / "bundle"
    root.mkdir()
    source = root / ("pytest-of-" + account.upper()) / "source.tsv"
    with pytest.raises(ValueError, match="account"):
        api.PortablePaths.create(root / "record.json", [source], root=root)


@pytest.mark.parametrize(
    "name", ["envperson", "-Users-envperson-work", "pytest-of-envperson", "prefix_ENVPERSON_suffix"]
)
@pytest.mark.parametrize("operation", ["create", "encode"])
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
    with pytest.raises(ValueError, match=r"^Provenance root equals a private home/temp prefix\."):
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
    with pytest.raises(
        ValueError, match=r"^Provenance root must be below a home or temporary root\."
    ):
        api.resolve_provenance_root(
            alias / "record.json", path_basis={"kind": "root-relative-v1", "root": "."}
        )


def test_account_substrings_are_not_identifying_path_tokens(
    private_environment: dict[str, Path], api: ModuleType
) -> None:
    root = private_environment["temporary"] / "bundle"
    root.mkdir()
    source = root / "prefixENVPERSONsuffix" / "source.tsv"
    try:
        paths = api.PortablePaths.create(root / "record.json", [source], root=root)
        assert paths.encode(source) == "prefixENVPERSONsuffix/source.tsv"
    except ValueError as error:
        raise AssertionError(
            "Account substrings must remain usable in ordinary path words"
        ) from error


def test_reader_account_does_not_reject_another_writers_reference(
    private_environment: dict[str, Path], api: ModuleType
) -> None:
    root = private_environment["temporary"] / "bundle"
    root.mkdir()
    source = root / "envperson" / "source.tsv"
    try:
        assert (
            api.resolve_reference(
                root / "record.json",
                "envperson/source.tsv",
                path_basis={"kind": "root-relative-v1", "root": "."},
            )
            == source
        )
    except ValueError as error:
        raise AssertionError("The reader's account must not reject a writer's reference") from error


@pytest.mark.parametrize("name", ["root", "x"])
def test_generic_and_single_character_account_tokens_are_not_identifying(
    tmp_path: Path,
    api: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    name: str,
) -> None:
    monkeypatch.setattr(api, "_account_names", lambda: {name})
    source = tmp_path / name / "input.tsv"
    try:
        paths = api.PortablePaths.create(tmp_path / "record.json", [source])
        assert paths.encode(source) == name + "/input.tsv"
    except ValueError as error:
        raise AssertionError("Generic and single-character names are exempt path tokens") from error


@pytest.mark.parametrize("name", ["al", "ana", "joan"])
@pytest.mark.parametrize(
    "template", ["-Users-{}-proj", "-home-{}-x", "pytest-of-{}", "{}_scratch", "{}.backup"]
)
@pytest.mark.parametrize("operation", ["create", "encode"])
def test_short_account_tokens_are_refused_in_scratch_paths(
    tmp_path: Path,
    api: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    name: str,
    template: str,
    operation: str,
) -> None:
    monkeypatch.setattr(api, "_account_names", lambda: {name})
    source = tmp_path / template.format(name.upper()) / "input.tsv"
    paths = api.PortablePaths(tmp_path / "record.json", tmp_path)

    def exercise() -> None:
        if operation == "create":
            api.PortablePaths.create(paths.record, [source], root=tmp_path)
        else:
            paths.encode(source)

    with pytest.raises(ValueError, match="local account name"):
        exercise()


@pytest.mark.parametrize("spelling", ["realana", "anas-drive"])
def test_each_account_token_boundary_is_required(
    tmp_path: Path,
    api: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    spelling: str,
) -> None:
    monkeypatch.setattr(api, "_account_names", lambda: {"ana"})
    source = tmp_path / spelling / "input.tsv"
    try:
        paths = api.PortablePaths.create(tmp_path / "record.json", [source])
        assert paths.encode(source) == spelling + "/input.tsv"
    except ValueError as error:
        raise AssertionError("Both account-token boundaries must match") from error


@pytest.mark.parametrize("name", ["tsv", "json"])
@pytest.mark.parametrize("compressed", [False, True])
def test_known_file_format_suffix_is_not_an_account_token(
    tmp_path: Path,
    api: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    name: str,
    compressed: bool,
) -> None:
    monkeypatch.setattr(api, "_account_names", lambda: {name})
    filename = "input." + name.upper() + (".gz" if compressed else "")
    source = tmp_path / filename
    source.write_text("fixture bytes\n")
    try:
        paths = api.PortablePaths.create(tmp_path / "record.json", [source])
        assert paths.encode(source) == filename
    except ValueError as error:
        raise AssertionError("A known file-format suffix is not an account locator") from error


@pytest.mark.parametrize("name", ["tsv", "json"])
@pytest.mark.parametrize("layout", ["parent", "stem", "directory", "unknown", "backup"])
def test_format_names_are_still_checked_outside_known_file_suffixes(
    tmp_path: Path,
    api: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    name: str,
    layout: str,
) -> None:
    monkeypatch.setattr(api, "_account_names", lambda: {name})
    if layout == "parent":
        source = tmp_path / name / "input.tsv"
        source.parent.mkdir()
        source.write_text("bytes\n")
    elif layout == "stem":
        source = tmp_path / (name + ".tsv")
        source.write_text("bytes\n")
    elif layout == "backup":
        source = tmp_path / (name + ".backup")
        source.write_text("bytes\n")
    else:
        source = tmp_path / ("input." + name)
        if layout == "directory":
            source.mkdir()
    with pytest.raises(ValueError, match="local account name"):
        api.PortablePaths.create(tmp_path / "record.json", [source], root=tmp_path)


@pytest.mark.parametrize("variable", ["TMPDIR", "TEMP", "TMP"])
def test_each_environment_temporary_boundary_is_independent(
    tmp_path: Path,
    private_environment: dict[str, Path],
    api: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    variable: str,
) -> None:
    root = tmp_path / ("independent-" + variable)
    root.mkdir()
    assert root != private_environment["temporary"]
    monkeypatch.setenv(variable, str(root))
    with pytest.raises(
        ValueError, match=r"^Provenance root must be below a home or temporary root\."
    ) as error:
        api.PortablePaths.create(root / "record.json", [], root=root)
    assert str(root) not in str(error.value)


@pytest.mark.parametrize(
    "parts",
    [
        ("tmp",),
        ("private", "tmp"),
        ("var", "tmp"),
        ("private", "var", "tmp"),
        ("var", "folders"),
        ("private", "var", "folders"),
    ],
)
def test_platform_temporary_boundaries_are_refused(
    private_environment: dict[str, Path],
    api: ModuleType,
    parts: tuple[str, ...],
) -> None:
    """Runtime authorities are elsewhere; the real platform-prefix list is untouched."""
    root = Path(Path.cwd().anchor).joinpath(*parts)
    assert root not in private_environment.values()

    def exercise() -> None:
        if root.is_dir():
            api.PortablePaths.create(root / "record.json", [], root=root)
        else:
            api.resolve_provenance_root(
                root / "record.json", path_basis={"kind": "root-relative-v1", "root": "."}
            )

    with pytest.raises(
        ValueError, match=r"^Provenance root must be below a home or temporary root\."
    ) as error:
        exercise()
    assert str(root) not in str(error.value)


@pytest.mark.parametrize("prefix", ["Users", "home"])
def test_another_account_home_cannot_start_below_the_root(
    tmp_path: Path,
    private_environment: dict[str, Path],
    api: ModuleType,
    prefix: str,
) -> None:
    root = tmp_path / "not-runtime-home" / prefix
    assert root not in private_environment.values()
    root.mkdir(parents=True)
    source = root / "anotherresident" / "source.tsv"
    with pytest.raises(
        ValueError, match=r"^Provenance root leaves a local prefix in a reference\."
    ) as error:
        api.PortablePaths.create(root / "record.json", [source], root=root)
    assert str(root) not in str(error.value)


@pytest.mark.parametrize("prefix", ["Users", "home"])
def test_reader_refuses_a_home_prefix_crossing_below_the_root(
    tmp_path: Path,
    private_environment: dict[str, Path],
    api: ModuleType,
    prefix: str,
) -> None:
    root = tmp_path / "not-runtime-home" / prefix
    assert root not in private_environment.values()
    root.mkdir(parents=True)
    with pytest.raises(
        ValueError, match=r"^Provenance root leaves a local prefix in a reference\."
    ) as error:
        api.resolve_reference(
            root / "record.json",
            "anotherresident/source.tsv",
            path_basis={"kind": "root-relative-v1", "root": "."},
        )
    assert str(root) not in str(error.value)
