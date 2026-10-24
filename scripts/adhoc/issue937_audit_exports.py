"""Audit canonical export declarations for #937 (2026-10-23).

Inputs: repository Python sources under src/mtp and scripts, including ignored files.
Output: JSON with source hashes and explicit cross-module public-use omissions.
This static audit follows import aliases/relative imports, not dynamic imports.
It never imports an analysis module or changes source files.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
from collections import defaultdict
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path


def declared_exports(tree: ast.Module, module: str) -> set[str] | None:
    """Read a single literal export declaration; malformed declarations fail closed."""
    declarations: list[ast.expr | None] = []
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if any(isinstance(target, ast.Name) and target.id == "__all__" for target in targets):
                declarations.append(node.value)
    if not declarations:
        return None
    if len(declarations) != 1:
        raise AssertionError(f"{module} repeats its export declaration")
    if declarations[0] is None:
        raise AssertionError(f"{module} has an uninitialized export declaration")
    try:
        values = ast.literal_eval(declarations[0])
    except (ValueError, TypeError) as error:
        raise AssertionError(f"{module} has a nonliteral export declaration") from error
    if not isinstance(values, (list, tuple)):
        raise AssertionError(f"{module} exports must be a list or tuple")
    if not all(isinstance(name, str) for name in values):
        raise AssertionError(f"{module} exports must be names")
    if len(set(values)) != len(values):
        raise AssertionError(f"Duplicate shim exports in {module}")
    return set(values)


def public_definitions(tree: ast.Module) -> set[str]:
    """Top-level public definitions/assignments, excluding incidental dependency imports."""
    names: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                names.update(bound_names(target))
    return {name for name in names if not name.startswith("_")}


def bound_names(target: ast.expr) -> set[str]:
    """Assignment to an attribute/subscript does not bind its base or index names."""
    if isinstance(target, ast.Name):
        return {target.id}
    if isinstance(target, (ast.Tuple, ast.List)):
        return {name for item in target.elts for name in bound_names(item)}
    if isinstance(target, ast.Starred):
        return bound_names(target.value)
    return set()


def module_name(path: Path, root: Path) -> str:
    relative = path.relative_to(root).with_suffix("")
    parts = relative.parts[:-1] if relative.name == "__init__" else relative.parts
    return ".".join(parts)


def import_origin(node: ast.ImportFrom, module: str, is_package: bool) -> str:
    if not node.level:
        return node.module or ""
    package = module if is_package else module.rpartition(".")[0]
    return importlib.util.resolve_name("." * node.level + (node.module or ""), package)


def function_arguments(arguments: ast.arguments) -> list[ast.arg]:
    """Exclude parameters of nested lambdas appearing inside a default expression."""
    return [
        *arguments.posonlyargs,
        *arguments.args,
        *arguments.kwonlyargs,
        *([arguments.vararg] if arguments.vararg is not None else []),
        *([arguments.kwarg] if arguments.kwarg is not None else []),
    ]


def cross_module_uses(
    sources: Mapping[str, str], packages: set[str]
) -> dict[tuple[str, str], set[str]]:
    """Index explicit imports and possible qualified uses within each lexical scope.

    Import bindings within a scope are conservative alternatives, without executing
    control flow. Nested function aliases cannot overwrite their siblings' bindings.
    """
    uses: dict[tuple[str, str], set[str]] = defaultdict(set)
    for consumer, source in sources.items():
        tree = ast.parse(source)

        def scope_nodes(nodes: Sequence[ast.AST]) -> Iterator[ast.AST]:
            def walk(node: ast.AST) -> Iterator[ast.AST]:
                yield node
                if not isinstance(
                    node,
                    (
                        ast.FunctionDef,
                        ast.AsyncFunctionDef,
                        ast.ClassDef,
                        ast.Lambda,
                        ast.ListComp,
                        ast.SetComp,
                        ast.DictComp,
                        ast.GeneratorExp,
                    ),
                ):
                    for child in ast.iter_child_nodes(node):
                        yield from walk(child)

            for statement in nodes:
                yield from walk(statement)

        def inspect_scope(
            body: Sequence[ast.AST],
            inherited: Mapping[str, set[str]],
            parameters: set[str],
            consumer_name: str,
            method_parent: Mapping[str, set[str]] | None = None,
        ) -> None:
            nodes = list(scope_nodes(body))
            aliases = dict(inherited)
            local_imports: dict[str, set[str]] = defaultdict(set)
            shadowed = set(parameters)
            for node in nodes:
                if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
                    shadowed.add(node.id)
                elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    shadowed.add(node.name)
                elif isinstance(node, ast.Import):
                    for alias in node.names:
                        name = alias.asname or alias.name.split(".")[0]
                        shadowed.add(name)
                        if alias.name in sources or alias.name.startswith("mtp."):
                            local_imports[name].add(
                                alias.name if alias.asname else alias.name.split(".")[0]
                            )
                elif isinstance(node, ast.ImportFrom):
                    origin = import_origin(node, consumer_name, consumer_name in packages)
                    for alias in node.names:
                        name = alias.asname or alias.name
                        shadowed.add(name)
                        qualified = f"{origin}.{alias.name}"
                        if qualified in sources:
                            local_imports[name].add(qualified)
                        elif origin in sources and not alias.name.startswith("_"):
                            uses[(origin, alias.name)].add(consumer_name)
            for name in shadowed:
                aliases.pop(name, None)
            aliases.update(local_imports)
            # A class namespace is a defining scope, but is not a closure for
            # methods, lambdas, comprehensions, or a nested class body.
            lexical_parent = aliases if method_parent is None else method_parent
            for node in nodes:
                if isinstance(node, ast.Attribute):
                    parts: list[str] = []
                    parent: ast.expr = node
                    while isinstance(parent, ast.Attribute):
                        parts.append(parent.attr)
                        parent = parent.value
                    if isinstance(parent, ast.Name):
                        for prefix in aliases.get(parent.id, set()):
                            qualified_parts = prefix.split(".") + list(reversed(parts))
                            origin, name = ".".join(qualified_parts[:-1]), qualified_parts[-1]
                            if origin in sources and not name.startswith("_"):
                                uses[(origin, name)].add(consumer_name)
                elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    argument_nodes = function_arguments(node.args)
                    arguments = {arg.arg for arg in argument_nodes}
                    outer_expressions = [
                        *node.decorator_list,
                        *node.args.defaults,
                        *(value for value in node.args.kw_defaults if value is not None),
                        *(arg.annotation for arg in argument_nodes if arg.annotation is not None),
                        *([node.returns] if node.returns is not None else []),
                    ]
                    inspect_scope(
                        outer_expressions,
                        aliases,
                        set(),
                        consumer_name,
                        method_parent=method_parent,
                    )
                    inspect_scope(
                        node.body,
                        lexical_parent,
                        arguments,
                        consumer_name,
                    )
                elif isinstance(node, ast.ClassDef):
                    inspect_scope(
                        [
                            *node.bases,
                            *node.decorator_list,
                            *(item.value for item in node.keywords),
                        ],
                        aliases,
                        set(),
                        consumer_name,
                        method_parent=method_parent,
                    )
                    inspect_scope(
                        node.body,
                        lexical_parent,
                        set(),
                        consumer_name,
                        method_parent=lexical_parent,
                    )
                elif isinstance(node, ast.Lambda):
                    arguments = {arg.arg for arg in function_arguments(node.args)}
                    inspect_scope(
                        [*node.args.defaults, *(v for v in node.args.kw_defaults if v is not None)],
                        aliases,
                        set(),
                        consumer_name,
                        method_parent=method_parent,
                    )
                    inspect_scope([node.body], lexical_parent, arguments, consumer_name)
                elif isinstance(node, (ast.ListComp, ast.SetComp, ast.DictComp, ast.GeneratorExp)):
                    # The first iterable is evaluated outside the comprehension's local scope.
                    inspect_scope(
                        [node.generators[0].iter],
                        aliases,
                        set(),
                        consumer_name,
                        method_parent=method_parent,
                    )
                    bindings = {name for gen in node.generators for name in bound_names(gen.target)}
                    expressions = [item for gen in node.generators for item in gen.ifs]
                    expressions.extend(gen.iter for gen in node.generators[1:])
                    if isinstance(node, ast.DictComp):
                        expressions.extend([node.key, node.value])
                    else:
                        expressions.append(node.elt)
                    inspect_scope(expressions, lexical_parent, bindings, consumer_name)

        inspect_scope(tree.body, {}, set(), consumer)
    return dict(uses)


def export_omissions(sources: Mapping[str, str], packages: set[str]) -> dict[str, list[str]]:
    """Find public canonical members used elsewhere but omitted from explicit __all__."""
    uses = cross_module_uses(sources, packages)
    missing: dict[str, list[str]] = {}
    for module, source in sources.items():
        if module != "mtp" and not module.startswith("mtp."):
            continue
        tree = ast.parse(source)
        exported = declared_exports(tree, module)
        if exported is None:
            continue
        public = public_definitions(tree)
        # Canonical re-exports count; imported third-party libraries do not.
        for node in tree.body:
            if isinstance(node, ast.Import):
                public.update(
                    alias.asname or alias.name.split(".")[0]
                    for alias in node.names
                    if alias.name in sources
                    and (alias.name == "mtp" or alias.name.startswith("mtp."))
                    and not (alias.asname or alias.name).startswith("_")
                )
            elif isinstance(node, ast.ImportFrom):
                origin = import_origin(node, module, module in packages)
                if origin in sources:
                    public.update(
                        alias.asname or alias.name
                        for alias in node.names
                        if not (alias.asname or alias.name).startswith("_")
                    )
        absent = sorted(
            name for name in public - exported if uses.get((module, name), set()) - {module}
        )
        if absent:
            missing[module] = absent
    return missing


def source_inventory(root: Path) -> tuple[dict[str, str], set[str], dict[str, str]]:
    """Read Python files directly, independent of gitignore; omit legacy mirrors as consumers."""
    sources: dict[str, str] = {}
    packages: set[str] = set()
    hashes: dict[str, str] = {}
    for directory, relative_root in [(root / "src/mtp", root / "src"), (root / "scripts", root)]:
        for path in sorted(directory.rglob("*.py")):
            name = module_name(path, relative_root)
            content = path.read_bytes()
            sources[name] = content.decode("utf-8")
            hashes[str(path.relative_to(root))] = hashlib.sha256(content).hexdigest()
            if path.name == "__init__.py":
                packages.add(name)
    return sources, packages, hashes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources, packages, hashes = source_inventory(args.root)
    omissions = export_omissions(sources, packages)
    uses = cross_module_uses(sources, packages)
    report = {
        "scope": "Explicit canonical public definitions/re-exports used by src/mtp or scripts",
        "limitations": (
            "Static imports and possible qualified attributes per lexical scope; "
            "runtime control flow and dynamic imports are not inferred"
        ),
        "source_count": len(sources),
        "source_sha256": hashes,
        "omissions": omissions,
        "consumers": {
            f"{module}.{name}": sorted(uses[(module, name)])
            for module, names in omissions.items()
            for name in names
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"source_count": len(sources), "omissions": omissions}, sort_keys=True))


if __name__ == "__main__":
    main()
