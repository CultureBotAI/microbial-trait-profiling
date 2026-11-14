"""Exercise the session snapshot backstop inside a synthetic checkout (#1026)."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

CONFTEXT = Path(__file__).with_name("conftest.py")
CHILD_TEST = """from pathlib import Path
import os
import subprocess
import sys
import conftest

def test_output_operation():
    root = Path(os.environ["MTP_GUARD_SYNTHETIC_ROOT"])
    action = os.environ["MTP_GUARD_ACTION"]
    tree = os.environ["MTP_GUARD_TREE"]
    target = root / tree / "probe.txt"
    if action == "roots":
        actual = getattr(conftest, "_PROTECTED_OUTPUT_DIRS", (conftest.GUARDED_DATA_DIR,))
        assert set(actual) == {root / "data", root / "trait-visualization"}
    elif action in {"create", "modify", "delete"}:
        # Audit hooks do not cross a process boundary. Keep proving that the
        # independent teardown snapshot catches all three kinds of mutation.
        script = (
            "from pathlib import Path; import sys; "
            "target = Path(sys.argv[1]); action = sys.argv[2]; "
            "target.unlink() if action == 'delete' else "
            "target.write_text(('created synthetic output' if action == 'create' else "
            "'modified synthetic output with a different size') + chr(10), encoding='utf-8')"
        )
        subprocess.run([sys.executable, "-c", script, str(target), action], check=True)
    else:
        assert action == "unchanged"
"""


def _guard_run(
    tmp_path: Path, tree: str, action: str
) -> tuple[subprocess.CompletedProcess[str], Path]:
    """Keep conftest bytes unchanged; only its imported root belongs to the child."""
    root = tmp_path / "checkout"
    tests = root / "tests"
    package = root / "src" / "bsdbng"
    tests.mkdir(parents=True)
    package.mkdir(parents=True)
    source = CONFTEXT.read_bytes()
    (tests / "conftest.py").write_bytes(source)
    (package / "__init__.py").write_text("", encoding="utf-8")
    (package / "paths.py").write_text(
        "from pathlib import Path\nREPO_ROOT = Path(__file__).resolve().parents[2]\n",
        encoding="utf-8",
    )
    for dirname in ("data", "trait-visualization", "outside"):
        (root / dirname).mkdir()
    target = root / tree / "probe.txt"
    if action in {"modify", "delete", "unchanged"}:
        target.write_text("original synthetic output\n", encoding="utf-8")
    (tests / "test_guard_child.py").write_text(CHILD_TEST, encoding="utf-8")
    executable = Path(sys.executable).with_name("pytest")
    assert executable.is_file(), "The current environment must provide its pytest executable"
    environment = {
        name: value
        for name, value in os.environ.items()
        if not name.startswith(("PYTEST_", "PYTHONPATH", "GIT_"))
    }
    environment.update(
        PYTHONPATH=str(root / "src"),
        PYTHONDONTWRITEBYTECODE="1",
        PYTEST_DISABLE_PLUGIN_AUTOLOAD="1",
        MTP_GUARD_SYNTHETIC_ROOT=str(root),
        MTP_GUARD_TREE=tree,
        MTP_GUARD_ACTION=action,
    )
    result = subprocess.run(
        [
            str(executable),
            "-q",
            "-o",
            "addopts=",
            "-p",
            "no:cacheprovider",
            str(tests / "test_guard_child.py"),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    (tmp_path / "child.stdout.txt").write_text(result.stdout, encoding="utf-8")
    (tmp_path / "child.stderr.txt").write_text(result.stderr, encoding="utf-8")
    assert CONFTEXT.read_bytes() == source
    assert (tests / "conftest.py").read_bytes() == source
    return result, target


def test_default_protected_roots_include_data_and_visualization(tmp_path: Path) -> None:
    result, _target = _guard_run(tmp_path, "trait-visualization", "roots")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "1 passed" in result.stdout


@pytest.mark.parametrize("tree", ["trait-visualization", "data"])
@pytest.mark.parametrize("action", ["create", "modify", "delete"])
def test_session_guard_refuses_protected_output_changes(
    tmp_path: Path, tree: str, action: str
) -> None:
    result, target = _guard_run(tmp_path, tree, action)
    output = result.stdout + result.stderr
    assert result.returncode == pytest.ExitCode.TESTS_FAILED, output
    assert "ERROR at teardown of test_output_operation" in output
    assert "tests changed real files under" in output
    assert tree in output
    assert "probe.txt" in output
    assert "1 passed, 1 error" in output
    assert "ERROR collecting" not in output
    if action == "delete":
        assert not target.exists()
    else:
        expected = (
            "created synthetic output\n"
            if action == "create"
            else "modified synthetic output with a different size\n"
        )
        assert target.read_text(encoding="utf-8") == expected


@pytest.mark.parametrize(
    ("tree", "action"), [("trait-visualization", "unchanged"), ("outside", "create")]
)
def test_session_guard_accepts_unchanged_or_unprotected_outputs(
    tmp_path: Path, tree: str, action: str
) -> None:
    result, target = _guard_run(tmp_path, tree, action)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "1 passed" in result.stdout
    expected = (
        "original synthetic output\n" if action == "unchanged" else "created synthetic output\n"
    )
    assert target.read_text(encoding="utf-8") == expected
