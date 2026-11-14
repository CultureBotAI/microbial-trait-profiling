"""Pin module-binding isolation and restoration in the #839 probe (#1031 item 3)."""

import importlib.util
import sys
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType

import pytest

SOURCE = (
    Path(__file__).resolve().parents[1] / "scripts/adhoc/issue839_probe_guild_layout_mutants.py"
)
SPEC = importlib.util.spec_from_file_location("issue839_binding_probe_1031", SOURCE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load the issue839 probe")
probe = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(probe)
fresh_test_modules = probe.fresh_test_modules
NAMES = ("test_trait_guilds", "test_table1_trait_maps")


@pytest.fixture(autouse=True)
def preserve_real_bindings() -> Iterator[None]:
    """Keep tests isolated even when a broken context manager fails an assertion."""
    previous = {name: sys.modules[name] for name in NAMES if name in sys.modules}
    try:
        yield
    finally:
        for name in NAMES:
            sys.modules.pop(name, None)
        sys.modules.update(previous)


def test_entry_removes_both_bindings() -> None:
    for name in NAMES:
        sys.modules[name] = ModuleType(name)
    with fresh_test_modules():
        assert all(name not in sys.modules for name in NAMES)


class InsideError(Exception):
    """Exception deliberately raised by the caller inside the context."""


@pytest.mark.parametrize(
    ("present", "raise_inside"),
    [
        pytest.param(present, raises, id=f"{label}-{'exception' if raises else 'normal'}")
        for present, label in [
            ((), "absent"),
            (NAMES[:1], "guild"),
            (NAMES[1:], "table"),
            (NAMES, "both"),
        ]
        for raises in (False, True)
    ],
)
def test_restores_bindings(present: tuple[str, ...], raise_inside: bool) -> None:
    for name in NAMES:
        sys.modules.pop(name, None)
    originals = {name: ModuleType(f"original_{name}") for name in present}
    sys.modules.update(originals)
    inner = {name: ModuleType(f"inner_{name}") for name in NAMES}
    expected_error = InsideError("caller failure")
    caught = None
    try:
        with fresh_test_modules():
            sys.modules.update(inner)
            if raise_inside:
                raise expected_error
    except InsideError as exc:
        caught = exc
    assert caught is (expected_error if raise_inside else None)
    for name in NAMES:
        if name in originals:
            assert sys.modules.get(name) is originals[name]
        else:
            assert name not in sys.modules
        assert sys.modules.get(name) is not inner[name]
