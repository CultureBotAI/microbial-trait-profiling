"""Static audit of canonical module export declarations and cross-module public usage."""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
from collections import defaultdict
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path

_MUTATING_ALL_METHODS = frozenset(
    {
        "append",
        "clear",
        "extend",
        "insert",
        "pop",
        "remove",
        "reverse",
        "sort",
        "__delitem__",
        "__iadd__",
        "__imul__",
        "__setitem__",
    }
)
_MUTATING_OPERATOR_NAMES = frozenset(
    {
        "delattr",
        "delitem",
        "iadd",
        "imul",
        "setattr",
        "setitem",
    }
)


def _root_target_name(expr: ast.expr) -> str | None:
    """Return the base identifier of a subscript/attribute chain, if any."""
    current = expr
    while isinstance(current, (ast.Subscript, ast.Attribute)):
        current = current.value
    return current.id if isinstance(current, ast.Name) else None


def _target_mutates_all(target: ast.expr) -> bool:
    """Check whether an assignment or deletion target mutates an element or attribute of __all__."""
    if isinstance(target, (ast.Subscript, ast.Attribute)):
        return _root_target_name(target) == "__all__"
    if isinstance(target, (ast.Tuple, ast.List)):
        return any(_target_mutates_all(item) for item in target.elts)
    if isinstance(target, ast.Starred):
        return _target_mutates_all(target.value)
    return False


def _pattern_bound_names(pattern: ast.pattern) -> set[str]:
    """Collect variable names bound by a structural pattern matching case."""
    names: set[str] = set()
    for node in ast.walk(pattern):
        if isinstance(node, (ast.MatchAs, ast.MatchStar)) and node.name is not None:
            names.add(node.name)
        elif isinstance(node, ast.MatchMapping) and node.rest is not None:
            names.add(node.rest)
    return names


def _call_mutates_all(node: ast.Call) -> bool:
    """Check whether a call expression mutates __all__ in place."""
    if isinstance(node.func, ast.Attribute):
        if (
            _root_target_name(node.func.value) == "__all__"
            and node.func.attr in _MUTATING_ALL_METHODS
        ):
            return True
        return bool(
            node.args
            and _root_target_name(node.args[0]) == "__all__"
            and node.func.attr in (_MUTATING_ALL_METHODS | _MUTATING_OPERATOR_NAMES)
        )
    if isinstance(node.func, ast.Name):
        return bool(
            node.args
            and _root_target_name(node.args[0]) == "__all__"
            and node.func.id in _MUTATING_OPERATOR_NAMES
        )
    return False


def _other_node_binds_all(node: ast.AST) -> bool:
    """Check whether a non-assignment AST node binds or redeclares __all__."""
    if isinstance(node, ast.NamedExpr):
        return node.target.id == "__all__"
    if isinstance(node, (ast.For, ast.AsyncFor)):
        return "__all__" in bound_names(node.target)
    if isinstance(node, (ast.With, ast.AsyncWith)):
        return any(
            item.optional_vars is not None and "__all__" in bound_names(item.optional_vars)
            for item in node.items
        )
    if isinstance(node, ast.ExceptHandler):
        return node.name == "__all__"
    if isinstance(node, ast.Match):
        return any("__all__" in _pattern_bound_names(case.pattern) for case in node.cases)
    if isinstance(node, (ast.Import, ast.ImportFrom)):
        return any((alias.asname or alias.name.split(".")[0]) == "__all__" for alias in node.names)
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
        return node.name == "__all__"
    if isinstance(node, (ast.Global, ast.Nonlocal)):
        return "__all__" in node.names
    return False


def _has_all_mutation_or_invalid_binding(tree: ast.Module, valid_declaration_ids: set[int]) -> bool:
    """Detect any in-place mutation, deletion, or non-top-level/compound binding of __all__."""
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if any(_target_mutates_all(target) for target in targets):
                return True
            if id(node) not in valid_declaration_ids and any(
                "__all__" in bound_names(target) for target in targets
            ):
                return True
        elif isinstance(node, ast.AugAssign):
            if "__all__" in bound_names(node.target) or _target_mutates_all(node.target):
                return True
        elif isinstance(node, ast.Delete):
            if any(
                "__all__" in bound_names(target) or _target_mutates_all(target)
                for target in node.targets
            ):
                return True
        elif isinstance(node, ast.Call):
            if _call_mutates_all(node):
                return True
        elif _other_node_binds_all(node):
            return True
    return False


def declared_exports(tree: ast.Module, module: str) -> set[str] | None:
    """Read a single literal export declaration; malformed declarations fail closed."""
    declarations: list[ast.expr | None] = []
    valid_declaration_ids: set[int] = set()
    for node in tree.body:
        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id == "__all__"
        ) or (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and node.target.id == "__all__"
        ):
            declarations.append(node.value)
            valid_declaration_ids.add(id(node))
    if len(declarations) > 1:
        raise AssertionError(f"{module} repeats its export declaration")
    if _has_all_mutation_or_invalid_binding(tree, valid_declaration_ids):
        raise AssertionError(f"{module} mutates or conditionally binds __all__")
    if not declarations:
        return None
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


def _is_type_checking_expr(expr: ast.expr) -> bool:
    """Return True for TYPE_CHECKING or typing.TYPE_CHECKING."""
    return (isinstance(expr, ast.Name) and expr.id == "TYPE_CHECKING") or (
        isinstance(expr, ast.Attribute)
        and isinstance(expr.value, ast.Name)
        and expr.value.id == "typing"
        and expr.attr == "TYPE_CHECKING"
    )


def _main_guard_polarity(expr: ast.expr) -> bool | None:
    """Return True for __name__ == '__main__', False for __name__ != '__main__', else None."""
    if not isinstance(expr, ast.Compare) or len(expr.ops) != 1 or len(expr.comparators) != 1:
        return None
    left, right = expr.left, expr.comparators[0]
    is_name_main = (
        isinstance(left, ast.Name)
        and left.id == "__name__"
        and isinstance(right, ast.Constant)
        and right.value == "__main__"
    ) or (
        isinstance(right, ast.Name)
        and right.id == "__name__"
        and isinstance(left, ast.Constant)
        and left.value == "__main__"
    )
    if not is_name_main:
        return None
    if isinstance(expr.ops[0], ast.Eq):
        return True
    if isinstance(expr.ops[0], ast.NotEq):
        return False
    return None


def _iter_runtime_module_statements(statements: Sequence[ast.stmt]) -> Iterator[ast.stmt]:
    """Yield module-scope statements including nested blocks, skipping static/main-only guards."""
    for stmt in statements:
        yield stmt
        if isinstance(stmt, ast.If):
            if _is_type_checking_expr(stmt.test):
                yield from _iter_runtime_module_statements(stmt.orelse)
            elif (
                isinstance(stmt.test, ast.UnaryOp)
                and isinstance(stmt.test.op, ast.Not)
                and _is_type_checking_expr(stmt.test.operand)
            ):
                yield from _iter_runtime_module_statements(stmt.body)
            else:
                main_polarity = _main_guard_polarity(stmt.test)
                if main_polarity is True:
                    yield from _iter_runtime_module_statements(stmt.orelse)
                elif main_polarity is False:
                    yield from _iter_runtime_module_statements(stmt.body)
                else:
                    yield from _iter_runtime_module_statements(stmt.body)
                    yield from _iter_runtime_module_statements(stmt.orelse)
        elif isinstance(stmt, (ast.For, ast.AsyncFor, ast.While)):
            yield from _iter_runtime_module_statements(stmt.body)
            yield from _iter_runtime_module_statements(stmt.orelse)
        elif isinstance(stmt, (ast.With, ast.AsyncWith)):
            yield from _iter_runtime_module_statements(stmt.body)
        elif isinstance(stmt, (ast.Try, ast.TryStar)):
            yield from _iter_runtime_module_statements(stmt.body)
            for handler in stmt.handlers:
                yield from _iter_runtime_module_statements(handler.body)
            yield from _iter_runtime_module_statements(stmt.orelse)
            yield from _iter_runtime_module_statements(stmt.finalbody)
        elif isinstance(stmt, ast.Match):
            for case in stmt.cases:
                yield from _iter_runtime_module_statements(case.body)


def _iter_statement_expr_nodes(stmt: ast.stmt) -> Iterator[ast.AST]:
    """Walk expressions attached directly to a statement, excluding child statements and lambdas."""

    def walk_expr(node: ast.AST) -> Iterator[ast.AST]:
        if isinstance(node, (ast.stmt, ast.Lambda)):
            return
        yield node
        for child in ast.iter_child_nodes(node):
            yield from walk_expr(child)

    for child in ast.iter_child_nodes(stmt):
        if isinstance(child, ast.match_case):
            if child.guard is not None:
                yield from walk_expr(child.guard)
        elif not isinstance(child, ast.stmt):
            yield from walk_expr(child)


def public_definitions(tree: ast.Module) -> set[str]:
    """Top-level public definitions/assignments, excluding incidental dependency imports."""
    names: set[str] = set()
    for node in _iter_runtime_module_statements(tree.body):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                names.update(bound_names(target))
        elif isinstance(node, (ast.AugAssign, ast.For, ast.AsyncFor)):
            names.update(bound_names(node.target))
        elif isinstance(node, (ast.With, ast.AsyncWith)):
            for item in node.items:
                if item.optional_vars is not None:
                    names.update(bound_names(item.optional_vars))
        elif isinstance(node, ast.Match):
            for case in node.cases:
                names.update(_pattern_bound_names(case.pattern))
        for expr_node in _iter_statement_expr_nodes(node):
            if isinstance(expr_node, ast.NamedExpr):
                names.add(expr_node.target.id)
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
        for node in _iter_runtime_module_statements(tree.body):
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


__all__ = [
    "bound_names",
    "cross_module_uses",
    "declared_exports",
    "export_omissions",
    "function_arguments",
    "import_origin",
    "main",
    "module_name",
    "public_definitions",
    "source_inventory",
]


if __name__ == "__main__":
    main()
