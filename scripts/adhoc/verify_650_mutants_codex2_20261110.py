"""Kill #650 guard-removal mutants in scratch copies; print verdicts, 2026-11-10."""

from __future__ import annotations

import ast
import contextlib
import importlib.util
import io
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path("scripts/build_slide_deck.py")
TEST = Path("tests/test_build_slide_deck_650.py")


class RemoveGuard(ast.NodeTransformer):
    def __init__(self, name: str) -> None:
        self.name = name
        self.removed = 0

    def visit_Call(self, node: ast.Call) -> ast.AST:
        if isinstance(node.func, ast.Name) and node.func.id == self.name:
            node.func.id = "removed_guard"
            self.removed += 1
        return self.generic_visit(node)


def main() -> None:
    spec = importlib.util.spec_from_file_location("checks_650", ROOT / TEST)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module spec for {ROOT / TEST}")
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    originals = {path: (ROOT / path).read_bytes() for path in (SCRIPT, TEST)}
    replacements = {
        "refuse_or_review": "lambda *_a: []",
        "check_results_cite_the_deck": "lambda *_a: None",
        "verified_links": "lambda links, *_a: links",
        "placed_slides": "lambda items, *_a: items",
        "check_slides": "lambda *_a: ({}, [])",
        "plate_problems": "lambda *_a: {}",
        "cellulose_problems": "lambda *_a: []",
    }
    try:
        with TemporaryDirectory() as directory:
            base = Path(directory)
            for index, name in enumerate(["baseline", *replacements]):
                root = base / str(index)
                root.mkdir()
                source = originals[SCRIPT].decode()
                if name != "baseline":
                    tree = ast.parse(source)
                    function = next(
                        n
                        for n in tree.body
                        if isinstance(n, ast.FunctionDef) and n.name == "build_deck"
                    )
                    remover = RemoveGuard(name)
                    remover.visit(function)
                    if remover.removed <= 0:
                        raise RuntimeError(f"Failed to remove guard call: {name}")
                    source = ast.unparse(tree) + f"\nremoved_guard = {replacements[name]}\n"
                target = root / SCRIPT
                target.parent.mkdir()
                target.write_text(source)
                module = checks.load_deck(target)
                try:
                    with (
                        pytest.MonkeyPatch.context() as patch,
                        contextlib.redirect_stdout(io.StringIO()),
                    ):
                        checks.exercise_build(module, root, patch)
                except AssertionError:
                    if name == "baseline":
                        raise
                    print(f"{name}: killed by AssertionError")
                else:
                    if name != "baseline":
                        raise AssertionError(f"survived: {name}")
    finally:
        for path, content in originals.items():
            if (ROOT / path).read_bytes() != content:
                raise RuntimeError(f"Worktree file modified during verification: {path}")
    print("7/7 killed; worktree sources byte-identical")


if __name__ == "__main__":
    main()
