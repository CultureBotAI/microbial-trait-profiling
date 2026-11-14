"""Load a single-module mutant from exclusive scratch bytes (#1184, 2026-11-13).

Inputs: this checkout, an mtp module's replacement bytes, selected test paths.
Outputs: a new module copy outside the checkout; no canonical file writes.
For single-threaded probes only: dependencies remain canonical. Probes requiring
whole-tree provenance must use a complete private checkout instead.
"""

from __future__ import annotations

import importlib
import importlib.util
import sys
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType

import mtp
from mtp import paths


def verify_checkout_identity(repository: Path) -> None:
    """Reject an installed or foreign mtp before creating any output."""
    root = repository.resolve(strict=True)
    if (
        Path(str(mtp.__file__)).resolve() != root / "src/mtp/__init__.py"
        or paths.REPO_ROOT.resolve() != root
    ):
        raise AssertionError("Imported mtp belongs to a different checkout")


def _test_modules(test_paths: set[Path]) -> dict[str, ModuleType]:
    return {
        name: module
        for name, module in sys.modules.copy().items()
        if module is not None
        and (filename := getattr(module, "__file__", None))
        and Path(filename).resolve() in test_paths
    }


@contextmanager
def isolated_module(
    repository: Path,
    module_name: str,
    replacement: bytes,
    output: Path,
    test_paths: Sequence[Path],
) -> Iterator[ModuleType]:
    """Stage new bytes and restore module/package/test bindings even on failure.

    The caller must verify the collected tests use the yielded module, and route
    pytest caches/temporary files to scratch. Existing output is never reused.
    """
    verify_checkout_identity(repository)
    root = repository.resolve(strict=True)
    parts = module_name.split(".")
    if len(parts) < 2 or parts[0] != "mtp" or not all(part.isidentifier() for part in parts):
        raise ValueError("Expected a qualified mtp module name")
    canonical = root / "src" / Path(*parts).with_suffix(".py")
    if not canonical.is_file() or canonical.resolve() != canonical:
        raise ValueError("Expected a canonical module file without source symlinks")
    destination = output.resolve()
    if destination.is_relative_to(root) or root.is_relative_to(destination):
        raise ValueError("Isolation output must be outside the canonical checkout")
    if output.exists() or output.is_symlink():
        raise FileExistsError("Isolation output already exists; require an exclusive directory")
    selected = {path.resolve(strict=True) for path in test_paths}
    if any(not path.is_relative_to(root / "tests") for path in selected):
        raise ValueError("Selected tests must belong to the canonical checkout")
    parent = importlib.import_module(".".join(parts[:-1]))
    if Path(str(parent.__file__)).resolve().parent != canonical.parent:
        raise AssertionError("Parent package belongs to a different checkout")
    previous_modules = {module_name: sys.modules[module_name]} if module_name in sys.modules else {}
    had_attribute = hasattr(parent, parts[-1])
    previous_attribute = getattr(parent, parts[-1], None)
    previous_tests = _test_modules(selected)
    previous_bytecode = sys.dont_write_bytecode
    destination.mkdir(parents=True, exist_ok=False)
    copied = destination / canonical.name
    with copied.open("xb") as handle:
        handle.write(replacement)
    spec = importlib.util.spec_from_file_location(module_name, copied)
    if spec is None or spec.loader is None:
        raise AssertionError("Cannot load isolated module")
    module = importlib.util.module_from_spec(spec)
    try:
        sys.dont_write_bytecode = True
        for name in previous_tests:
            sys.modules.pop(name, None)
        sys.modules[module_name] = module
        setattr(parent, parts[-1], module)
        spec.loader.exec_module(module)
        if Path(str(module.__file__)).resolve() != copied or module.__spec__ is not spec:
            raise AssertionError("Loaded module origin differs from the isolated copy")
        yield module
    finally:
        for name in _test_modules(selected):
            sys.modules.pop(name, None)
        sys.modules.update(previous_tests)
        sys.modules.pop(module_name, None)
        sys.modules.update(previous_modules)
        if had_attribute:
            setattr(parent, parts[-1], previous_attribute)
        elif hasattr(parent, parts[-1]):
            delattr(parent, parts[-1])
        sys.dont_write_bytecode = previous_bytecode
