"""Prove write prevention and teardown detection in a pytester checkout (#1026)."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest_plugins = ["pytester"]
CONFTEXT = Path(__file__).with_name("conftest.py")


def synthetic_checkout(pytester: pytest.Pytester, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Load the actual guard unchanged, with only its imported repository root isolated."""
    root = pytester.path
    package = root / "src" / "bsdbng"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("", encoding="utf-8")
    (package / "paths.py").write_text(
        "from pathlib import Path\nREPO_ROOT = Path(__file__).resolve().parents[2]\n",
        encoding="utf-8",
    )
    for name in ("data", "trait-visualization", "outside"):
        (root / name).mkdir()
    pytester.makeconftest(CONFTEXT.read_text(encoding="utf-8"))
    monkeypatch.setenv("PYTHONPATH", str(root / "src"))
    monkeypatch.setenv("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")
    monkeypatch.setenv("PYTHONDONTWRITEBYTECODE", "1")
    monkeypatch.delenv("PYTEST_ADDOPTS", raising=False)
    return root


WRITE_CASES = {
    "default": "def emit(path=target):\n    path.write_text('changed')\nemit()",
    "keyword-default": "def emit(*, path=target):\n    path.write_text('changed')\nemit()",
    "tuple": "outputs = (target,)\noutputs[0].write_text('changed')",
    "composite": "outputs = {'paths': [target]}\noutputs['paths'][0].write_text('changed')",
    "dynamic": "(root / tree / 'probe.txt').write_text('changed')",
    "cached-open": "with cached_open(target, 'w') as handle:\n    handle.write('changed')",
    "cached-io": "with cached_io_open(target, 'w') as handle:\n    handle.write('changed')",
    "os-open": "descriptor = os.open(target, os.O_WRONLY | os.O_TRUNC)\nos.close(descriptor)",
    "rename-in": "(root / 'incoming.txt').replace(target)",
    "rename-out": "target.rename(root / 'outgoing.txt')",
    "unlink": "target.unlink()",
    "mkdir": "(root / tree / 'new-directory').mkdir()",
}


@pytest.mark.parametrize("tree", ["data", "trait-visualization"])
@pytest.mark.parametrize("operation", WRITE_CASES)
def test_guard_refuses_before_default_composite_or_dynamic_write(
    pytester: pytest.Pytester, monkeypatch: pytest.MonkeyPatch, tree: str, operation: str
) -> None:
    root = synthetic_checkout(pytester, monkeypatch)
    target = root / tree / "probe.txt"
    original = b"original output\n"
    target.write_bytes(original)
    incoming = root / "incoming.txt"
    incoming.write_text("replacement", encoding="utf-8")
    body = WRITE_CASES[operation]
    pytester.makepyfile(
        test_attempt=(
            "import builtins, io, os\nfrom pathlib import Path\n"
            "root = Path(__file__).parent\n"
            f"tree = {tree!r}\ntarget = root / tree / 'probe.txt'\n"
            "cached_open = builtins.open\ncached_io_open = io.open\n"
            "def test_attempt():\n" + "\n".join("    " + line for line in body.splitlines())
        )
    )
    result = pytester.runpytest_subprocess("-q", "-p", "no:cacheprovider", timeout=60)
    result.assert_outcomes(failed=1)
    result.stdout.fnmatch_lines(["*tests must not write protected output*"])
    assert "ERROR at teardown" not in result.stdout.str()
    assert target.read_bytes() == original
    assert incoming.read_text(encoding="utf-8") == "replacement"
    assert not (root / "outgoing.txt").exists()
    assert not (root / tree / "new-directory").exists()


def test_pytester_prevention_and_session_teardown_backstop(
    pytester: pytest.Pytester, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = synthetic_checkout(pytester, monkeypatch)
    target = root / "trait-visualization" / "probe.txt"
    original = b"unchanged output\n"
    target.write_bytes(original)
    # First prove the new pre-write guard is active in this actual pytest session.
    pytester.makepyfile(
        test_attempt="from pathlib import Path\n"
        "def test_attempt():\n"
        "    (Path(__file__).parent / 'trait-visualization/probe.txt').write_text('blocked')\n"
    )
    blocked = pytester.runpytest_subprocess("-q", "-p", "no:cacheprovider", timeout=60)
    blocked.assert_outcomes(failed=1)
    assert target.read_bytes() == original

    # Child-process writes do not inherit Python audit hooks. Exercise the retained
    # teardown assertion through that real limitation, without disabling either guard.
    pytester.makepyfile(
        test_attempt="from pathlib import Path\nimport subprocess, sys\n"
        "def test_attempt():\n"
        "    target = Path(__file__).parent / 'trait-visualization/probe.txt'\n"
        "    subprocess.run([sys.executable, '-c', "
        "\"from pathlib import Path; import sys; Path(sys.argv[1]).write_text('injected')\", "
        "str(target)], check=True)\n"
    )
    injected = pytester.runpytest_subprocess("-q", "-p", "no:cacheprovider", timeout=60)
    injected.assert_outcomes(passed=1, errors=1)
    injected.stdout.fnmatch_lines(
        [
            "*ERROR at teardown of test_attempt*",
            "*tests changed real files under protected output roots*",
        ]
    )
    assert target.read_text(encoding="utf-8") == "injected"

    pytester.makepyfile(
        test_attempt="from pathlib import Path\n"
        "def test_attempt():\n"
        "    root = Path(__file__).parent\n"
        "    assert (root / 'trait-visualization/probe.txt').read_text() == 'injected'\n"
        "    (root / 'outside/result.txt').write_text('allowed')\n"
    )
    quiet = pytester.runpytest_subprocess("-q", "-p", "no:cacheprovider", timeout=60)
    quiet.assert_outcomes(passed=1)
    assert (root / "outside/result.txt").read_text(encoding="utf-8") == "allowed"


@pytest.mark.parametrize(
    ("tree", "operation"),
    [("data", "unlink"), ("trait-visualization", "rename"), (None, "unlink-cycle")],
    ids=["unlink-protected-target", "rename-protected-target", "unlink-cycle"],
)
def test_outside_symlink_entry_cleanup_preserves_protected_targets(
    pytester: pytest.Pytester,
    monkeypatch: pytest.MonkeyPatch,
    tree: str | None,
    operation: str,
) -> None:
    """Removing or renaming an outside link acts on its entry, not its referent."""
    root = synthetic_checkout(pytester, monkeypatch)
    protected = {
        "data": b"original data\n",
        "trait-visualization": b"original visualization\n",
    }
    for name, content in protected.items():
        (root / name / "probe.txt").write_bytes(content)
    link = root / "outside" / "entry"
    peer = root / "outside" / "cycle-peer"
    destination = root / "outside" / "moved-entry"
    if tree is None:
        link.symlink_to(peer.name)
        peer.symlink_to(link.name)
    else:
        link.symlink_to(root / tree / "probe.txt")
    action = (
        "entry.rename(root / 'outside/moved-entry')" if operation == "rename" else "entry.unlink()"
    )
    pytester.makepyfile(
        test_cleanup="from pathlib import Path\n"
        "def test_cleanup():\n"
        "    root = Path(__file__).parent\n"
        "    entry = root / 'outside/entry'\n"
        f"    {action}\n"
    )
    result = pytester.runpytest_subprocess("-q", "-p", "no:cacheprovider", timeout=60)
    result.assert_outcomes(passed=1)
    assert not link.is_symlink()
    assert not link.exists()
    if operation == "rename":
        assert destination.is_symlink()
        assert destination.readlink() == root / "trait-visualization" / "probe.txt"
    else:
        assert not destination.exists()
        assert not destination.is_symlink()
    if tree is None:
        assert peer.is_symlink()
        assert peer.readlink() == Path(link.name)
    for name, content in protected.items():
        assert (root / name / "probe.txt").read_bytes() == content


@pytest.mark.parametrize("tree", ["data", "trait-visualization"])
@pytest.mark.parametrize(
    "operation",
    [
        "write-final-link",
        "unlink-parent-link",
        "rename-parent-link",
        "rename-parent-destination",
    ],
)
def test_symlink_writes_and_protected_parent_entries_remain_guarded(
    pytester: pytest.Pytester,
    monkeypatch: pytest.MonkeyPatch,
    tree: str,
    operation: str,
) -> None:
    """Writable opens follow the final link; entry operations still resolve parents."""
    root = synthetic_checkout(pytester, monkeypatch)
    target = root / tree / "probe.txt"
    original = b"protected output\n"
    target.write_bytes(original)
    incoming = root / "outside" / "incoming-output"
    incoming_bytes = b"unprotected incoming output\n"
    incoming.write_bytes(incoming_bytes)
    link = root / "outside" / "protected-link"
    if operation == "write-final-link":
        link.symlink_to(target)
        action = "entry.write_text('changed')"
        entry = "outside/protected-link"
    else:
        link.symlink_to(root / tree, target_is_directory=True)
        entry = "outside/protected-link/probe.txt"
        if operation == "unlink-parent-link":
            action = "entry.unlink()"
        elif operation == "rename-parent-link":
            action = "entry.rename(root / 'outside/moved-output')"
        else:
            action = "(root / 'outside/incoming-output').rename(entry)"
    pytester.makepyfile(
        test_attempt="from pathlib import Path\n"
        "def test_attempt():\n"
        "    root = Path(__file__).parent\n"
        f"    entry = root / {entry!r}\n"
        f"    {action}\n"
    )
    result = pytester.runpytest_subprocess("-q", "-p", "no:cacheprovider", timeout=60)
    result.assert_outcomes(failed=1)
    result.stdout.fnmatch_lines(["*tests must not write protected output*"])
    assert "ERROR at teardown" not in result.stdout.str()
    assert link.is_symlink()
    assert target.read_bytes() == original
    assert incoming.read_bytes() == incoming_bytes
    assert not (root / "outside" / "moved-output").exists()
    assert not (root / "outside" / "moved-output").is_symlink()


@pytest.mark.parametrize("tree", ["data", "trait-visualization"])
@pytest.mark.parametrize("operation", ["unlink", "rename"])
def test_protected_root_symlink_entry_cannot_be_removed(
    pytester: pytest.Pytester,
    monkeypatch: pytest.MonkeyPatch,
    tree: str,
    operation: str,
) -> None:
    """The protected root entry stays protected when its storage is reached by a link."""
    root = synthetic_checkout(pytester, monkeypatch)
    protected = root / tree
    storage = root / "outside" / "storage"
    storage.mkdir()
    original = b"protected backing storage\n"
    (storage / "probe.txt").write_bytes(original)
    protected.rmdir()
    protected.symlink_to(storage, target_is_directory=True)
    action = (
        "entry.unlink()" if operation == "unlink" else "entry.rename(root / 'outside/moved-root')"
    )
    pytester.makepyfile(
        test_attempt="from pathlib import Path\n"
        "def test_attempt():\n"
        "    root = Path(__file__).parent\n"
        f"    entry = root / {tree!r}\n"
        f"    {action}\n"
    )
    result = pytester.runpytest_subprocess("-q", "-p", "no:cacheprovider", timeout=60)
    result.assert_outcomes(failed=1)
    result.stdout.fnmatch_lines(["*tests must not write protected output*"])
    assert "ERROR at teardown" not in result.stdout.str()
    assert protected.is_symlink()
    assert protected.readlink() == storage
    assert (storage / "probe.txt").read_bytes() == original
    assert not (root / "outside" / "moved-root").exists()
    assert not (root / "outside" / "moved-root").is_symlink()
