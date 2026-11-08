"""Guard the newly declared correction-family registry (#851).

Discovery parses every Python file below src/ and scripts/ without gitignore
filtering. It resolves import/assignment aliases, named BH/BY/FDR helpers, the
project's adjust/correct wrappers, and local forwarding wrappers to a fixed
point. Definitions of hand-written correction helpers are registered too.
Each source file is size-capped and reduced to a bounded compact summary before
the next file is parsed; alias expansion and graph traversal are bounded too.
This is an AST guard, not a proof about dynamically imported functions or
runtime-selected methods. It never imports analysis modules or runs a test.
"""

from __future__ import annotations

import ast
import re
import unicodedata
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "methods/multiple_testing_families.yaml"


@dataclass(frozen=True)
class Limits:
    """Deterministic work limits; they are not an RSS measurement or guarantee."""

    source_bytes: int = 1_048_576
    ast_nodes: int = 100_000
    ast_depth: int = 200
    files: int = 1_000
    summary_items: int = 100_000
    summary_bytes: int = 33_554_432
    name_characters: int = 4_096
    alias_steps: int = 64
    derivation_steps: int = 1_000_000
    file_derivation_steps: int = 1_000_000


DEFAULT_LIMITS = Limits()
Sources = Mapping[str, str] | Iterable[tuple[str, str]]


class DiscoveryLimitError(ValueError):
    """Refuse a scan before its declared static-work budget is exceeded."""


@dataclass
class FileDerivationBudget:
    """Share a work allowance across candidates while one file's AST is live."""

    limit: int
    work: int = 0

    def charge(self, amount: int, path: str, scope: str) -> None:
        self.work += amount
        if self.work > self.limit:
            raise DiscoveryLimitError(
                f"per-file derivation work budget exceeded: {path}:{scope} "
                f"(attempted {self.work}, limit {self.limit}). "
                "Simplify the file or review and raise Limits.file_derivation_steps."
            )


@dataclass
class SummaryBudget:
    limits: Limits
    items: int = 0
    bytes: int = 0

    def charge(self, *values: str) -> None:
        self.items += 1
        # The fixed allowance is conservative bookkeeping, not Python's sizeof.
        self.bytes += 256 + sum(len(value.encode("utf-8")) for value in values)
        if self.items > self.limits.summary_items:
            raise DiscoveryLimitError("compact summary item budget exceeded")
        if self.bytes > self.limits.summary_bytes:
            raise DiscoveryLimitError("compact summary byte budget exceeded")


PRIMITIVES = {
    "step_up",
    "step_up_array",
    "false_discovery_control",
    "benjamini_hochberg",
    "family_bh_q_values",
    "multipletests",
    "adjusted_probabilities",
    "correct_family",
    "correct_inventory",
}
NAMED_WRAPPERS = {"mtp.taxon_influence.adjust", "mtp.log_ratio_sensitivity.correct"}
HELPER_NAME = re.compile(r"(?:^|_)(?:bh|by|fdr|benjamini)(?:_|$)")
PROCEDURES = {"bh": "BH", "by": "BY", "fdr_bh": "BH", "fdr_by": "BY"}
REQUIRED = {
    "family_id",
    "area",
    "hypotheses",
    "family_definition",
    "size_source",
    "primary_procedure",
    "comparison_procedure",
    "threshold",
    "outputs",
    "consumers",
    "policy_rule",
    "deliberate_difference",
    "open_issue",
    "anchors",
}


def dotted(node: ast.AST, limits: Limits = DEFAULT_LIMITS) -> str:
    """Return a simple qualified expression, without evaluating it."""
    parts: list[str] = []
    count = 0
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        count += len(node.attr) + 1
        if count > limits.name_characters:
            raise DiscoveryLimitError("qualified name length budget exceeded")
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
    elif not parts:
        return ""
    name = ".".join(reversed(parts))
    if len(name) > limits.name_characters:
        raise DiscoveryLimitError("qualified name length budget exceeded")
    return name


@dataclass
class Unit:
    path: str
    module: str
    tree: ast.Module
    limits: Limits = DEFAULT_LIMITS
    functions: dict[str, ast.FunctionDef | ast.AsyncFunctionDef] = field(default_factory=dict)
    aliases: dict[str, dict[str, str]] = field(default_factory=dict)
    calls: list[tuple[str, ast.Call]] = field(default_factory=list)

    def resolve(self, expression: ast.AST, scope: str = "") -> str:
        name = dotted(expression, self.limits)
        scope = "" if scope == "<module>" else scope
        pieces = scope.split(".") if scope else []
        scopes = [".".join(pieces[:index]) for index in range(len(pieces), -1, -1)]
        seen: set[str] = set()
        steps = 0
        while name:
            steps += 1
            if steps > self.limits.alias_steps:
                raise DiscoveryLimitError("alias expansion step budget exceeded")
            prefix, _, suffix = name.partition(".")
            if prefix in seen:
                break
            seen.add(prefix)
            alias = next(
                (self.aliases[s][prefix] for s in scopes if prefix in self.aliases.get(s, {})), None
            )
            if alias is not None:
                if len(alias) + len(suffix) + bool(suffix) > self.limits.name_characters:
                    raise DiscoveryLimitError("alias expansion length budget exceeded")
                name = alias + (f".{suffix}" if suffix else "")
                continue
            for local_scope in scopes:
                qualified = f"{local_scope}.{name}" if local_scope else name
                if qualified in self.functions:
                    return f"{self.module}.{qualified}"
            break
        return name


class Visitor(ast.NodeVisitor):
    def __init__(self, unit: Unit) -> None:
        self.unit = unit
        self.stack: list[str] = []

    def bind(self, name: str, target: str) -> None:
        if max(len(name), len(target)) > self.unit.limits.name_characters:
            raise DiscoveryLimitError("alias binding length budget exceeded")
        self.unit.aliases.setdefault(".".join(self.stack), {})[name] = target

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        self.stack.append(node.name)
        self.generic_visit(node)
        self.stack.pop()

    def visit_FunctionDef(self, node: ast.FunctionDef | ast.AsyncFunctionDef) -> None:
        self.stack.append(node.name)
        self.unit.functions[".".join(self.stack)] = node
        self.generic_visit(node)
        self.stack.pop()

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            self.bind(
                alias.asname or alias.name.split(".")[0],
                alias.name if alias.asname else alias.name.split(".")[0],
            )

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        module = node.module or ""
        if node.level:
            parts = self.unit.module.split(".")[: -node.level]
            module = ".".join([*parts, module]).rstrip(".")
        for alias in node.names:
            self.bind(alias.asname or alias.name, f"{module}.{alias.name}")

    def visit_Assign(self, node: ast.Assign) -> None:
        value = dotted(node.value, self.unit.limits)
        if value:
            for target in node.targets:
                if isinstance(target, ast.Name):
                    self.bind(target.id, value)
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        self.unit.calls.append((".".join(self.stack) or "<module>", node))
        self.generic_visit(node)


def is_seed(name: str) -> bool:
    leaf = name.rsplit(".", 1)[-1]
    return (
        leaf in PRIMITIVES
        or leaf.strip("_") in {"bh", "by", "fdr", "correct", "adjust"}
        or name in NAMED_WRAPPERS
        or bool(
            HELPER_NAME.search(leaf)
            and any(word in leaf for word in ("adjust", "correct", "q_value", "benjamini"))
        )
    )


def owned_nodes(function: ast.FunctionDef | ast.AsyncFunctionDef) -> Iterator[ast.AST]:
    """Do not credit a nested function's return to its enclosing function."""
    pending: list[ast.AST] = list(function.body)
    while pending:
        node = pending.pop()
        if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
            continue
        yield node
        pending.extend(ast.iter_child_nodes(node))


def function_dependencies(unit: Unit, budget: SummaryBudget) -> dict[str, set[str]]:
    """Reduce a file's functions to dependency strings while its AST is in memory."""
    dependencies: dict[str, set[str]] = {}
    for name, function in unit.functions.items():
        assigned: dict[str, set[str]] = defaultdict(set)
        returned: set[str] = set()
        targets: set[str] = set()
        for node in owned_nodes(function):
            if isinstance(node, ast.Assign):
                calls = {
                    unit.resolve(child.func, name)
                    for child in ast.walk(node.value)
                    if isinstance(child, ast.Call)
                }
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        assigned[target.id].update(calls)
            if isinstance(node, ast.Return) and node.value is not None:
                for child in ast.walk(node.value):
                    if isinstance(child, ast.Call):
                        targets.add(unit.resolve(child.func, name))
                    if isinstance(child, ast.Name):
                        returned.add(child.id)
        for name_returned in returned:
            targets.update(assigned[name_returned])
        qualified = f"{unit.module}.{name}"
        for target_name in targets:
            budget.charge(qualified, target_name)
        dependencies[qualified] = targets
    return dependencies


def forwarding_names(dependencies: dict[str, set[str]]) -> set[str]:
    """Fixed point over bounded strings only; source and ASTs are already released."""
    names = {name for name in dependencies if is_seed(name)}
    reverse: dict[str, set[str]] = defaultdict(set)
    initial = set(names)
    for qualified, targets in dependencies.items():
        for target in targets:
            reverse[target].add(qualified)
            if is_seed(target):
                initial.add(target)
    pending = list(initial)
    visited: set[str] = set()
    while pending:
        target = pending.pop()
        if target in visited:
            continue
        visited.add(target)
        for caller in reverse.get(target, set()):
            if caller not in names:
                names.add(caller)
                pending.append(caller)
    return names


def literal_procedures(call: ast.Call, target: str) -> set[str]:
    """Infer literal method/by flags; dynamic policies remain explicitly unchecked."""
    leaf = target.rsplit(".", 1)[-1]
    keywords = {keyword.arg: keyword.value for keyword in call.keywords}
    method = keywords.get("method")
    if method is None and leaf == "multipletests" and len(call.args) >= 3:
        method = call.args[2]
    if method is None and target == "mtp.log_ratio_sensitivity.correct" and len(call.args) >= 3:
        method = call.args[2]
    if isinstance(method, ast.Constant) and isinstance(method.value, str):
        return {PROCEDURES.get(method.value.lower(), method.value)}
    if method is not None:
        return set()
    if target == "mtp.taxon_influence.adjust":
        by = keywords.get("by", call.args[1] if len(call.args) > 1 else ast.Constant(False))
        return {"BY" if by.value else "BH"} if isinstance(by, ast.Constant) else set()
    if leaf in {"benjamini_hochberg", "family_bh_q_values", "false_discovery_control"}:
        return {"BH"}
    if leaf == "multipletests":
        return {"holm-sidak"}
    return set()


Fact = tuple[str, str]
AliasOrigin = tuple[tuple[str, ...], str]


class DerivationEnvironment(dict[str, Fact]):
    """Arithmetic facts and separately joined possible object references."""

    def __init__(self, facts: Mapping[str, Fact] | None = None) -> None:
        super().__init__({} if facts is None else facts)
        self.aliases: dict[str, frozenset[AliasOrigin]] = {}


def visible_by_derivation(
    unit: Unit,
    scope: str,
    *,
    work_counts: dict[str, int] | None = None,
    file_budget: FileDerivationBudget | None = None,
) -> bool:
    """Recognize local BH * H_m followed by clip, with the same visible extent.

    This deliberately small recognizer is an anchor-level witness, not output
    lineage or a general Python interpreter. Opaque helpers and mutating calls
    are outside its model. Facts exist only while this one bounded AST is live.
    Writes through unrelated parameters/self can refuse a free-vector witness;
    a post-BH write into an np.float64/np.int64 allocation also refuses it. Only
    default/unshadowed built-in numeric dtypes get the fixed-size shortcut.
    These conservative refusals are intentional (#1928), not evidence that the
    source's arithmetic is incorrect. Literal BH registration is unaffected.
    Optional counters accumulate attempted work in thirteen fixed categories,
    including a charge that exceeds the limit. Reusing a dictionary accumulates
    counts but never shares the per-invocation budget. Overlapping operations
    remain separate charges; declared_names includes all special bound names.
    """
    function = unit.functions[scope]
    witnessed = False
    work = 0

    def charge(category: str, amount: int = 1) -> None:
        nonlocal work
        work += amount
        if work_counts is not None:
            work_counts[category] = work_counts.get(category, 0) + amount
        if work > unit.limits.derivation_steps:
            raise DiscoveryLimitError(
                f"derivation work budget exceeded: {unit.path}:{scope} "
                f"(attempted {work}, limit {unit.limits.derivation_steps}). "
                "Simplify the function or review and raise Limits.derivation_steps."
            )
        if file_budget is not None:
            file_budget.charge(amount, unit.path, scope)

    def copy_env(env: DerivationEnvironment) -> DerivationEnvironment:
        charge("env_copy", len(env))
        copied = DerivationEnvironment(env)
        charge("env_copy", len(env.aliases))
        copied.aliases = env.aliases.copy()
        return copied

    def local_nodes(root: ast.AST) -> Iterator[ast.AST]:
        pending = [root]
        while pending:
            charge("walk_node")
            node = pending.pop()
            yield node
            if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
                continue
            for _, value in ast.iter_fields(node):
                if isinstance(value, list):
                    # iter_child_nodes scans string entries too (e.g. globals).
                    charge("walk_list", len(value))
            pending.extend(ast.iter_child_nodes(node))

    module_writes = set()
    for node in local_nodes(unit.tree):
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store | ast.Del):
            module_writes.add(node.id)
        elif isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
            module_writes.add(node.name)

    def number(node: ast.AST, value: int) -> bool:
        return (
            isinstance(node, ast.Constant)
            and type(node.value) in {int, float}
            and node.value == value
        )

    def target(call: ast.Call, env: DerivationEnvironment) -> str:
        # A locally rebound callable is not the imported primitive.
        root = dotted(call.func, unit.limits).split(".")[0]
        return "" if root in env or root in module_writes else unit.resolve(call.func, scope)

    def extent(call: ast.Call, env: DerivationEnvironment, expected: str) -> str | None:
        if target(call, env) != expected or call.keywords or len(call.args) not in {2, 3}:
            return None
        start, stop = call.args[:2]
        if not number(start, 1) or (len(call.args) == 3 and not number(call.args[2], 1)):
            return None
        if (
            not isinstance(stop, ast.BinOp)
            or not isinstance(stop.op, ast.Add)
            or not number(stop.right, 1)
        ):
            return None
        size = expression(stop.left, env)
        return size[1] if size is not None and size[0] == "size" else None

    def harmonic(call: ast.Call, env: DerivationEnvironment) -> Fact | None:
        if len(call.args) != 1 or call.keywords:
            return None
        value = call.args[0]
        family = None
        if (
            target(call, env) == "numpy.sum"
            and isinstance(value, ast.BinOp)
            and isinstance(value.op, ast.Div)
            and number(value.left, 1)
            and isinstance(value.right, ast.Call)
        ):
            family = extent(value.right, env, "numpy.arange")
        elif (
            target(call, env) == "sum"
            and isinstance(value, ast.GeneratorExp)
            and len(value.generators) == 1
        ):
            generator = value.generators[0]
            term = value.elt
            if (
                not generator.ifs
                and not generator.is_async
                and isinstance(generator.target, ast.Name)
                and isinstance(generator.iter, ast.Call)
                and isinstance(term, ast.BinOp)
                and isinstance(term.op, ast.Div)
                and number(term.left, 1)
                and isinstance(term.right, ast.Name)
                and term.right.id == generator.target.id
            ):
                family = extent(generator.iter, env, "range")
        return ("harmonic", family) if family is not None else None

    def expression(node: ast.AST, env: DerivationEnvironment) -> Fact | None:
        charge("expression")
        if isinstance(node, ast.Name):
            return env.get(node.id, ("vector", node.id))
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
            left, right = expression(node.left, env), expression(node.right, env)
            if left and right and {left[0], right[0]} == {"BH", "harmonic"} and left[1] == right[1]:
                return ("scaled", left[1])
        if not isinstance(node, ast.Call):
            return None
        name = target(node, env)
        if name == "float" and len(node.args) == 1 and not node.keywords:
            value = expression(node.args[0], env)
            return value if value and value[0] == "harmonic" else None
        if name in {"numpy.sum", "sum"}:
            return harmonic(node, env)
        if name == "len" and len(node.args) == 1 and not node.keywords:
            value = expression(node.args[0], env)
            return ("size", value[1]) if value and value[0] in {"vector", "array"} else None
        if (
            name in {"numpy.ones", "numpy.zeros"}
            and len(node.args) == 1
            and all(
                k.arg == "dtype"
                and isinstance(k.value, ast.Name)
                and k.value.id not in env
                and k.value.id not in module_writes
                and unit.resolve(k.value, scope) in {"float", "int", "bool"}
                for k in node.keywords
            )
        ):
            # Object/unknown dtypes may hold references inserted later. Their
            # assignments use the ordinary container model, not this shortcut.
            value = expression(node.args[0], env)
            return ("array", value[1]) if value and value[0] == "size" else None
        if (
            name.rsplit(".", 1)[-1] == "benjamini_hochberg"
            and len(node.args) == 1
            and not node.keywords
        ):
            value = expression(node.args[0], env)
            return ("BH", value[1]) if value and value[0] in {"vector", "array"} else None
        if name == "numpy.clip" and len(node.args) == 3 and not node.keywords:
            value = expression(node.args[0], env)
            if (
                value
                and value[0] == "scaled"
                and number(node.args[1], 0)
                and number(node.args[2], 1)
            ):
                return ("BY", value[1])
        return None

    def charge_origins(category: str, possible: frozenset[AliasOrigin]) -> None:
        charge(category, len(possible))
        for path, origin in possible:
            charge(category, len(path) + len(origin))
            for selector in path:
                charge(category, len(selector))

    no_origins: frozenset[AliasOrigin] = frozenset()
    singleton_origins: dict[tuple[str, bool], frozenset[AliasOrigin]] = {}
    resolved_origins: dict[
        tuple[str, Fact, int, bool], tuple[frozenset[AliasOrigin], frozenset[AliasOrigin]]
    ] = {}
    merged_origins: dict[
        tuple[int, int],
        tuple[frozenset[AliasOrigin], frozenset[AliasOrigin], frozenset[AliasOrigin]],
    ] = {}
    selected_origins: dict[
        tuple[int, str], tuple[frozenset[AliasOrigin], frozenset[AliasOrigin]]
    ] = {}
    contained_origins: dict[
        tuple[int, str], tuple[frozenset[AliasOrigin], frozenset[AliasOrigin]]
    ] = {}
    subtracted_origins: dict[
        tuple[int, int],
        tuple[frozenset[AliasOrigin], frozenset[AliasOrigin], frozenset[AliasOrigin]],
    ] = {}
    wrapper_identities: dict[int, tuple[frozenset[AliasOrigin], frozenset[str]]] = {}
    matched_prefixes: dict[
        tuple[int, int],
        tuple[frozenset[AliasOrigin], frozenset[AliasOrigin], frozenset[tuple[str, ...]]],
    ] = {}
    reference_indexes: dict[
        int, tuple[frozenset[AliasOrigin], dict[str, frozenset[tuple[str, ...]]]]
    ] = {}
    partitioned_origins: dict[
        int, tuple[frozenset[AliasOrigin], frozenset[str], frozenset[str]]
    ] = {}

    def direct_origin(origin: str, *, external: bool = False) -> frozenset[AliasOrigin]:
        # Per-invocation immutable singletons share the same reference through
        # copied environments. Charge each lookup and the full construction on
        # a miss before allocating the string, tuple, set and cache entry.
        charge("expression")
        key = origin, external
        existing = singleton_origins.get(key)
        if existing is not None:
            return existing
        charge("expression", len(origin) + (len("external/") if external else 0) + 4)
        value = "external/" + origin if external else origin
        created = frozenset({((), value)})
        singleton_origins[key] = created
        return created

    def merge_origins(*groups: frozenset[AliasOrigin]) -> frozenset[AliasOrigin]:
        merged: frozenset[AliasOrigin] = frozenset()
        for group in groups:
            charge("expression")
            if not group or group is merged:
                continue
            if not merged:
                merged = group
                continue
            # Strong input references in each entry prevent identity-key reuse.
            # These operations depend only on immutable inputs, never on env.
            # Defensive lifetime guarantee (#1915, #1929): source fixtures do
            # not force allocator identity reuse, and other caches can retain
            # these inputs too. No semantic mutant kill is claimed for dropping
            # these references; passing probes do not justify dropping them.
            charge("expression")
            left_id, right_id = id(merged), id(group)
            key = (left_id, right_id) if left_id < right_id else (right_id, left_id)
            existing = merged_origins.get(key)
            if existing is not None:
                merged = existing[2]
                continue
            left, right = merged, group
            # Reuse the larger immutable input when it already contains the
            # smaller one. A real subset scan and any new union remain charged.
            if len(group) > len(merged):
                merged, group = group, merged
            charge_origins("expression", group)
            if not group <= merged:
                charge_origins("expression", merged)
                charge_origins("expression", group)
                merged = merged | group
            charge("expression", 4)
            merged_origins[key] = left, right, merged
        return merged

    def named_origins(name: str, env: DerivationEnvironment) -> frozenset[AliasOrigin]:
        possible = env.aliases.get(name, no_origins)
        fact = env.get(name, ("vector", name))
        if fact[0] != "vector":
            return possible
        is_free = name not in env
        # Read the current binding first. Only the ensuing pure computation is
        # cached; rebinding, free-name materialization and side-map updates all
        # produce a different key. No environment is retained in this cache.
        # The free bit is defensive (#1915, #1929): no reachable same-fact,
        # same-side-map free/materialized collision is established by tests.
        # Mutation-driven materialization invalidates the fact before rereading;
        # an explicit alias assignment retains external side references instead.
        charge("expression")
        key = name, fact, id(possible), is_free
        existing = resolved_origins.get(key)
        if existing is not None:
            return existing[1]
        original = possible
        possible = merge_origins(possible, direct_origin(fact[1]))
        if is_free:
            # Retain external uncertainty through local aliases/containers;
            # inspecting only the eventual write root would lose it.
            possible = merge_origins(possible, direct_origin(name, external=True))
        charge("expression", 7)
        resolved_origins[key] = original, possible
        return possible

    def select_origins(possible: frozenset[AliasOrigin], selector: str) -> frozenset[AliasOrigin]:
        charge("expression")
        if not possible:
            return possible
        key = id(possible), selector
        existing = selected_origins.get(key)
        if existing is not None:
            return existing[1]
        charge_origins("expression", possible)
        selected = frozenset(
            (path[1:], origin) for path, origin in possible if path and path[0] == selector
        )
        charge("expression", len(selector) + 3)
        selected_origins[key] = possible, selected
        return selected

    def contain_origins(possible: frozenset[AliasOrigin], selector: str) -> frozenset[AliasOrigin]:
        charge("expression")
        if not possible:
            return possible
        key = id(possible), selector
        existing = contained_origins.get(key)
        if existing is not None:
            return existing[1]
        charge_origins("expression", possible)
        charge("expression", len(possible) * (len(selector) + 1))
        contained = frozenset(((selector, *path), origin) for path, origin in possible)
        charge("expression", len(selector) + 3)
        contained_origins[key] = possible, contained
        return contained

    def possible_origins(node: ast.AST, env: DerivationEnvironment) -> frozenset[AliasOrigin]:
        """Track references, never infer arithmetic from an ambiguous value."""
        charge("expression")
        if isinstance(node, ast.Name):
            return named_origins(node.id, env)
        if isinstance(node, ast.NamedExpr):
            return possible_origins(node.value, env)
        if isinstance(node, ast.BoolOp):
            possible: frozenset[AliasOrigin] = frozenset()
            for value in node.values:
                possible = merge_origins(possible, possible_origins(value, env))
            return possible
        if isinstance(node, ast.IfExp):
            return merge_origins(
                possible_origins(node.body, env), possible_origins(node.orelse, env)
            )
        if isinstance(node, ast.List | ast.Tuple | ast.Set):
            possible = frozenset()
            for element in node.elts:
                possible = merge_origins(possible, possible_origins(element, env))
            return contain_origins(possible, "[]")
        if isinstance(node, ast.Dict):
            possible = frozenset()
            for value in node.values:
                possible = merge_origins(possible, possible_origins(value, env))
            return contain_origins(possible, "[]")
        if isinstance(node, ast.Subscript):
            possible = possible_origins(node.value, env)
            if isinstance(node.slice, ast.Slice):
                # A shallow copy owns its outer sequence, but retains references
                # to any objects inside it. Direct vector extents are not shared.
                charge_origins("expression", possible)
                return frozenset(item for item in possible if item[0])
            return select_origins(possible, "[]")
        if isinstance(node, ast.Attribute):
            charge("expression", len(node.attr) + 1)
            return select_origins(possible_origins(node.value, env), "." + node.attr)
        if isinstance(node, ast.Call):
            name = target(node, env)
            if name in {"len", "float", "sum", "numpy.sum"}:
                return frozenset()
            if (
                name in {"numpy.ones", "numpy.zeros"}
                and len(node.args) == 1
                and all(k.arg == "dtype" for k in node.keywords)
            ):
                # These allocations retain no shape/dtype object references.
                # In particular, dtype=float is not an external-vector alias.
                # Keep dispatch signatures such as like= opaque.
                return frozenset()
            if name == "typing.cast" and len(node.args) == 2 and not node.keywords:
                return possible_origins(node.args[1], env)
            if name == "contextlib.nullcontext" and len(node.args) == 1 and not node.keywords:
                return possible_origins(node.args[0], env)
            if name == "types.SimpleNamespace":
                possible = frozenset()
                for keyword in node.keywords:
                    if keyword.arg is not None:
                        charge("expression", len(keyword.arg) + 1)
                        possible = merge_origins(
                            possible,
                            contain_origins(
                                possible_origins(keyword.value, env), "." + keyword.arg
                            ),
                        )
                return possible
            copied = None
            if (
                name in {"list", "tuple", "copy.copy", "sorted"}
                and len(node.args) == 1
                and not node.keywords
            ):
                copied = node.args[0]
            elif (
                isinstance(node.func, ast.Attribute)
                and node.func.attr == "copy"
                and not node.args
                and not node.keywords
            ):
                copied = node.func.value
            if copied is not None:
                possible = possible_origins(copied, env)
                charge_origins("expression", possible)
                return frozenset(item for item in possible if item[0])
            # An opaque binding may return an argument. It does not provide an
            # arithmetic witness, and merely binding it does not invalidate one.
            possible = frozenset()
            for argument in node.args:
                possible = merge_origins(possible, possible_origins(argument, env))
            for keyword in node.keywords:
                possible = merge_origins(possible, possible_origins(keyword.value, env))
            return possible
        return frozenset()

    def root_identities(references: frozenset[AliasOrigin]) -> frozenset[str]:
        charge("write_alias_scan")
        if not references:
            return frozenset()
        existing = wrapper_identities.get(id(references))
        if existing is not None:
            return existing[1]
        charge_origins("write_alias_scan", references)
        identities = frozenset(origin for path, origin in references if not path)
        charge("write_alias_scan", 3)
        wrapper_identities[id(references)] = references, identities
        return identities

    def reference_index(
        references: frozenset[AliasOrigin],
    ) -> dict[str, frozenset[tuple[str, ...]]]:
        charge("write_alias_scan")
        existing = reference_indexes.get(id(references))
        if existing is not None:
            return existing[1]
        # Index the immutable set once rather than scanning every origin/path
        # again for each different wrapper queried against that same set.
        charge_origins("write_alias_scan", references)
        grouped: dict[str, set[tuple[str, ...]]] = {}
        for path, origin in references:
            if origin not in grouped:
                charge("write_alias_scan", 2)
                grouped[origin] = set()
            charge("write_alias_scan")
            grouped[origin].add(path)
        charge("write_alias_scan", len(grouped))
        indexed: dict[str, frozenset[tuple[str, ...]]] = {}
        for origin, paths in grouped.items():
            charge("write_alias_scan", len(paths) + 1)
            indexed[origin] = frozenset(paths)
        charge("write_alias_scan", 3)
        reference_indexes[id(references)] = references, indexed
        return indexed

    def reference_prefixes(
        root_references: frozenset[AliasOrigin],
        identities: frozenset[str],
        references: frozenset[AliasOrigin],
    ) -> frozenset[tuple[str, ...]]:
        charge("write_alias_scan")
        if not identities or not references:
            return frozenset()
        key = id(root_references), id(references)
        existing = matched_prefixes.get(key)
        if existing is not None:
            return existing[2]
        indexed = reference_index(references)
        charge("write_alias_scan", len(identities))
        matched: set[tuple[str, ...]] = set()
        for origin in identities:
            paths = indexed.get(origin, frozenset())
            charge("write_alias_scan", len(paths))
            for path in paths:
                charge("write_alias_scan", len(path))
                for selector in path:
                    charge("write_alias_scan", len(selector))
            matched.update(paths)
        charge("write_alias_scan", len(matched) + 1)
        prefixes = frozenset(matched)
        charge("write_alias_scan", 4)
        matched_prefixes[key] = root_references, references, prefixes
        return prefixes

    def bind_origins(
        name: str, possible: frozenset[AliasOrigin], env: DerivationEnvironment
    ) -> None:
        fact = env.get(name)
        if possible and fact is not None and fact[0] == "vector":
            # The arithmetic fact already retains this exact direct identity.
            # Keep the side map for additional possibilities/contained objects.
            direct = direct_origin(fact[1])
            if possible is direct:
                possible = frozenset()
            else:
                charge("assignment")
                key = id(possible), id(direct)
                existing = subtracted_origins.get(key)
                if existing is not None:
                    possible = existing[2]
                else:
                    original = possible
                    charge_origins("assignment", direct)
                    if direct <= possible:
                        charge_origins("assignment", possible)
                        charge_origins("assignment", direct)
                        possible = possible - direct
                    charge("assignment", 4)
                    subtracted_origins[key] = original, direct, possible
        charge("assignment")
        existing_value = env.aliases.get(name)
        if existing_value is possible or (existing_value is None and not possible):
            return
        if possible:
            charge("assignment", len(possible))
            env.aliases[name] = possible
        else:
            env.aliases.pop(name, None)

    def mutation_origins(node: ast.AST, env: DerivationEnvironment) -> frozenset[str]:
        possible = possible_origins(node, env)
        charge("write_alias_scan")
        if not possible:
            return frozenset()
        partition = partitioned_origins.get(id(possible))
        if partition is None:
            charge_origins("write_alias_scan", possible)
            origins = frozenset(
                origin
                for path, origin in possible
                if not path and not origin.startswith("external/")
            )
            charge_origins("write_alias_scan", possible)
            external = frozenset(
                origin.removeprefix("external/")
                for path, origin in possible
                if not path and origin.startswith("external/")
            )
            charge("write_alias_scan", 4)
            partition = possible, origins, external
            partitioned_origins[id(possible)] = partition
        _, origins, external = partition
        charge("write_alias_scan", len(origins))
        parameter_write = not origins.isdisjoint(parameter_names)
        if external or parameter_write:
            # Unknown external objects can alias incoming arguments. The final
            # syntactic Unit.aliases map cannot establish their temporal identity.
            # Only free/incoming extent keys are candidates; a local allocation's
            # position-tagged identity and scalar/array values are not aliases.
            charge("write_alias_scan", len(env))
            incoming = {
                fact[1]
                for fact in env.values()
                if fact[0] in {"vector", "size", "harmonic", "BH", "scaled", "BY"}
                and fact[1].isidentifier()
                and (external or fact[1] not in parameter_names)
            }
            charge("write_alias_scan", len(external))
            unresolved = any(not (name in env and env[name][1] != name) for name in external)
            if unresolved or parameter_write:
                charge("write_alias_scan", len(origins) + len(incoming))
                origins = origins | incoming
        return origins

    def remember_origins(origins: frozenset[str], env: DerivationEnvironment) -> None:
        charge("write_alias_scan", len(origins))
        # Only ordered traversals may insert environment/cache entries. Charge
        # the comparisons (including string lengths) and sorted-list allocation.
        # The #1929 regression observes insertion order itself. It does not
        # claim this sort changes current BY credits or real-producer charges.
        charge(
            "write_alias_scan",
            len(origins) + sum(map(len, origins)) * len(origins).bit_length(),
        )
        for origin in sorted(origins):
            if origin.isidentifier():
                env.setdefault(origin, ("vector", origin))

    def remember_mutated_origin(previous: Fact | None, env: DerivationEnvironment) -> None:
        # A free vector read may exist only through its aliases. Include its
        # original name in the mutation scan so a later read cannot recreate
        # the old extent. Position-tagged local assignments have no free origin;
        # an original name already rebound or invalidated must keep its fact.
        if previous and previous[0] == "vector" and previous[1].isidentifier():
            env.setdefault(previous[1], previous)

    def assign(
        node: ast.AST,
        value: Fact | None,
        env: DerivationEnvironment,
        possible: frozenset[AliasOrigin] = frozenset(),
    ) -> None:
        charge("assignment")
        if isinstance(node, ast.Name):
            env[node.id] = value or ("vector", f"{node.id}:{node.lineno}:{node.col_offset}")
            bind_origins(node.id, possible, env)
        elif isinstance(node, ast.Subscript | ast.Attribute):
            # Only a direct write into a known fixed-size allocation preserves
            # extent. Other writes may mutate a view of any tracked BH result.
            if (
                isinstance(node, ast.Subscript)
                and isinstance(node.value, ast.Name)
                and env.get(node.value.id, ("unknown", ""))[0] == "array"
                and not isinstance(node.ctx, ast.Del)
            ):
                return
            root = node.value
            while isinstance(root, ast.Subscript | ast.Attribute):
                root = root.value
            previous = expression(root, env)
            wrapper_aliases: dict[
                str, tuple[frozenset[AliasOrigin], frozenset[tuple[str, ...]]]
            ] = {}
            if isinstance(root, ast.Name):
                root_references = named_origins(root.id, env)
                identities = root_identities(root_references)
                # Mutating a wrapper changes its arithmetic extent, not the
                # identity shared by its aliases. Retain direct and contained
                # identity paths, so insertion reaches every known spelling
                # without invalidating an enclosing wrapper's own extent.
                charge("write_alias_scan", len(env))
                for name in env:
                    if not env.aliases.get(name):
                        # A materialized name with no side references has just
                        # its vector identity, or no references at all. Avoid
                        # indexing every unrelated singleton for every write.
                        fact = env[name]
                        charge("write_alias_scan")
                        if fact[0] != "vector" or fact[1] not in identities:
                            continue
                        charge("assignment", 3)
                        wrapper_aliases[name] = direct_origin(fact[1]), frozenset({()})
                        continue
                    references = named_origins(name, env)
                    prefixes = reference_prefixes(root_references, identities, references)
                    if prefixes:
                        wrapper_aliases[name] = references, prefixes
                if root.id not in wrapper_aliases:
                    charge("assignment")
                    wrapper_aliases[root.id] = root_references, frozenset({()})
            origins = mutation_origins(node.value, env)
            remember_mutated_origin(previous, env)
            remember_origins(origins, env)
            if isinstance(root, ast.Name):
                env[root.id] = ("unknown", root.id)
            charge("write_alias_scan", len(env))
            for alias, fact in list(env.items()):
                if (
                    fact == previous
                    or fact[0] in {"BH", "scaled", "BY"}
                    or (fact[0] == "vector" and fact[1] in origins)
                ):
                    env[alias] = ("unknown", alias)
            for alias, (references, _) in wrapper_aliases.items():
                bind_origins(alias, references, env)
            if not isinstance(node.ctx, ast.Del) and possible and isinstance(root, ast.Name):
                # The write changes the wrapper, not the inserted object. Keep
                # possible references for a later write *through* the new slot.
                selectors: list[str] = []
                child: ast.AST = node
                while isinstance(child, ast.Subscript | ast.Attribute):
                    charge(
                        "assignment", 3 if isinstance(child, ast.Subscript) else len(child.attr) + 2
                    )
                    selector = "[]" if isinstance(child, ast.Subscript) else "." + child.attr
                    selectors.append(selector)
                    child = child.value
                contained = possible
                for selector in selectors:
                    contained = contain_origins(contained, selector)
                for alias, (_, prefixes) in wrapper_aliases.items():
                    charge("assignment", len(prefixes))
                    charge("assignment", sum(map(len, prefixes)))
                    charge(
                        "assignment",
                        len(prefixes)
                        + sum(len(s) + 1 for p in prefixes for s in p) * len(prefixes).bit_length(),
                    )
                    for prefix in sorted(prefixes):
                        at_alias = contained
                        for selector in reversed(prefix):
                            at_alias = contain_origins(at_alias, selector)
                        bind_origins(
                            alias,
                            merge_origins(env.aliases.get(alias, frozenset()), at_alias),
                            env,
                        )
        elif isinstance(node, ast.Tuple | ast.List):
            elements = select_origins(possible, "[]")
            for element in node.elts:
                assign(element, None, env, elements)

    def invalidate_aliases(node: ast.AST, env: DerivationEnvironment) -> None:
        if isinstance(node, ast.Name):
            previous = env.get(node.id, ("vector", node.id))
            origins = mutation_origins(node, env)
            remember_mutated_origin(previous, env)
            remember_origins(origins, env)
            charge("alias_invalidation", len(env))
            for alias, fact in list(env.items()):
                if fact == previous or (fact[0] == "vector" and fact[1] in origins):
                    env[alias] = ("unknown", f"{alias}:{node.lineno}")
        else:
            # Augmented slot assignment first mutates the slot's referent,
            # unlike an ordinary replacement of the wrapper's slot.
            origins = mutation_origins(node, env)
            remember_origins(origins, env)
            charge("alias_invalidation", len(env))
            for alias, fact in list(env.items()):
                if fact[0] == "vector" and fact[1] in origins:
                    env[alias] = ("unknown", alias)
            assign(node, None, env)

    def intersection(
        left: DerivationEnvironment, right: DerivationEnvironment
    ) -> DerivationEnvironment:
        # Keep invalidated names present, so an earlier fact cannot reappear
        # through the fallback for an unbound parameter.
        charge("env_join", len(left) + len(right))
        joined = DerivationEnvironment(
            {
                name: left[name]
                if name in left and left.get(name) == right.get(name)
                else ("unknown", name)
                for name in dict.fromkeys((*left, *right))
            }
        )
        for name in joined:
            left_possible = left.aliases.get(name, frozenset())
            right_possible = right.aliases.get(name, frozenset())
            same_fact = left.get(name) == right.get(name)
            if same_fact and not left_possible and not right_possible:
                continue
            if not same_fact:
                if name in left:
                    left_possible = named_origins(name, left)
                if name in right:
                    right_possible = named_origins(name, right)
            if left_possible is right_possible or not right_possible:
                possible = left_possible
            elif not left_possible:
                possible = right_possible
            else:
                possible = merge_origins(left_possible, right_possible)
            if possible:
                joined.aliases[name] = possible
        return joined

    def same_aliases(left: DerivationEnvironment, right: DerivationEnvironment) -> bool:
        charge("loop_pass", len(left.aliases) + len(right.aliases))
        if len(left.aliases) != len(right.aliases):
            return False
        for name, possible in left.aliases.items():
            other = right.aliases.get(name)
            if possible is other:
                continue
            if other is None:
                return False
            charge_origins("loop_pass", possible)
            charge_origins("loop_pass", other)
            if possible != other:
                return False
        return True

    def invalidate_bindings(statement: ast.stmt, env: DerivationEnvironment) -> None:
        """Special bindings are unknown, including those inside unsupported blocks."""
        for node in local_nodes(statement):
            names: list[str] = []
            if isinstance(node, ast.NamedExpr) and isinstance(node.target, ast.Name):
                possible = possible_origins(node.value, env)
                env[node.target.id] = ("unknown", node.target.id)
                bind_origins(node.target.id, possible, env)
                names.append(node.target.id)
            elif isinstance(node, ast.Import | ast.ImportFrom):
                charge("import_names", len(node.names))
                names.extend(
                    alias.asname
                    or (alias.name.split(".")[0] if isinstance(node, ast.Import) else alias.name)
                    for alias in node.names
                )
            elif isinstance(node, ast.MatchAs | ast.MatchStar) and node.name:
                names.append(node.name)
            elif isinstance(node, ast.MatchMapping) and node.rest:
                names.append(node.rest)
            elif isinstance(node, ast.ExceptHandler) and node.name:
                names.append(node.name)
            elif isinstance(node, ast.Global | ast.Nonlocal):
                # Declaration-aware resolution of enclosing/callable bindings
                # is outside this recognizer; pin conservative rejection (#1609).
                names = node.names
            elif isinstance(node, ast.Delete):
                for item in node.targets:
                    if isinstance(item, ast.Name):
                        env.aliases.pop(item.id, None)
                        names.append(item.id)
                    else:
                        assign(item, None, env)
            # Global/nonlocal names are string fields, not child AST nodes.
            # Repeated declarations can do much work while adding only one fact.
            charge("declared_names", len(names))
            for name in names:
                env[name] = ("unknown", name)

    def terminates(body: list[ast.stmt]) -> bool:
        for statement in body:
            charge("terminates")
            if isinstance(statement, ast.Return | ast.Raise | ast.Break | ast.Continue):
                return True
            if isinstance(statement, ast.With | ast.AsyncWith) and terminates(statement.body):
                return True
            if isinstance(statement, ast.If):
                if isinstance(statement.test, ast.Constant):
                    branch = statement.body if statement.test.value else statement.orelse
                    if terminates(branch):
                        return True
                elif terminates(statement.body) and terminates(statement.orelse):
                    return True
        return False

    loop_paths: list[tuple[list[DerivationEnvironment], list[DerivationEnvironment]]] = []

    def bind_expression(
        node: ast.AST, value: ast.AST, env: DerivationEnvironment, fact: Fact | None
    ) -> None:
        if (
            isinstance(node, ast.Tuple | ast.List)
            and isinstance(value, ast.Tuple | ast.List)
            and len(node.elts) == len(value.elts)
        ):
            charge("assignment", len(node.elts))
            evaluated = [
                (target, expression(item, env), possible_origins(item, env))
                for target, item in zip(node.elts, value.elts, strict=True)
            ]
            for target_node, fact, possible in evaluated:
                assign(target_node, fact, env, possible)
        else:
            possible = (
                possible_origins(value, env)
                if fact is None or fact[0] in {"vector", "unknown"}
                else frozenset()
            )
            assign(node, fact, env, possible)

    def statements(
        body: list[ast.stmt], env: DerivationEnvironment, *, credit: bool = True
    ) -> DerivationEnvironment:
        nonlocal witnessed
        for statement in body:
            charge("statement")
            invalidate_bindings(statement, env)
            if isinstance(statement, ast.Assign | ast.AnnAssign):
                value = expression(statement.value, env) if statement.value is not None else None
                targets = (
                    statement.targets if isinstance(statement, ast.Assign) else [statement.target]
                )
                for item in targets:
                    if statement.value is not None:
                        bind_expression(item, statement.value, env, value)
                    else:
                        assign(item, value, env)
                if (
                    credit
                    and value
                    and value[0] == "BY"
                    and any(isinstance(item, ast.Name) for item in targets)
                ):
                    witnessed = True
            elif isinstance(statement, ast.AugAssign):
                invalidate_aliases(statement.target, env)
            elif isinstance(statement, ast.Return):
                if statement.value is not None:
                    value = expression(statement.value, env)
                    if credit and value and value[0] == "BY":
                        witnessed = True
                break
            elif isinstance(statement, ast.Raise | ast.Break | ast.Continue):
                if loop_paths and isinstance(statement, ast.Break | ast.Continue):
                    exits = loop_paths[-1][0 if isinstance(statement, ast.Break) else 1]
                    exits.append(copy_env(env))
                break
            elif isinstance(statement, ast.If):
                if isinstance(statement.test, ast.Constant):
                    branch = statement.body if statement.test.value else statement.orelse
                    env = statements(branch, env, credit=credit)
                    if terminates(branch):
                        break
                else:
                    env = intersection(
                        statements(statement.body, copy_env(env), credit=credit),
                        statements(statement.orelse, copy_env(env), credit=credit),
                    )
                    if terminates(statement.body) and terminates(statement.orelse):
                        break
            elif isinstance(statement, ast.For | ast.AsyncFor | ast.While):
                if (
                    isinstance(statement, ast.While)
                    and isinstance(statement.test, ast.Constant)
                    and not statement.test.value
                ):
                    env = statements(statement.orelse, env, credit=credit)
                    continue
                # The descending intersection reaches a fixed point before
                # any loop-body calculation is allowed to supply a witness.
                entry = copy_env(env)
                while True:
                    charge("loop_pass")
                    inside = copy_env(entry)
                    if isinstance(statement, ast.For | ast.AsyncFor):
                        assign(
                            statement.target,
                            None,
                            inside,
                            select_origins(possible_origins(statement.iter, inside), "[]"),
                        )
                    break_paths: list[DerivationEnvironment] = []
                    continue_paths: list[DerivationEnvironment] = []
                    loop_paths.append((break_paths, continue_paths))
                    inside = statements(statement.body, inside, credit=False)
                    loop_paths.pop()
                    updated = intersection(entry, inside)
                    for continued in continue_paths:
                        updated = intersection(updated, continued)
                    if updated == entry and same_aliases(updated, entry):
                        break
                    entry = updated
                inside = copy_env(entry)
                if isinstance(statement, ast.For | ast.AsyncFor):
                    assign(
                        statement.target,
                        None,
                        inside,
                        select_origins(possible_origins(statement.iter, inside), "[]"),
                    )
                break_paths, continue_paths = [], []
                loop_paths.append((break_paths, continue_paths))
                inside = statements(statement.body, inside, credit=credit)
                loop_paths.pop()
                normal = intersection(env, inside)
                for continued in continue_paths:
                    normal = intersection(normal, continued)
                # A break skips the loop's else clause, so break states join after it.
                env = statements(statement.orelse, normal, credit=credit)
                for exited in break_paths:
                    env = intersection(env, exited)
            elif isinstance(statement, ast.With | ast.AsyncWith):
                for context_item in statement.items:
                    if context_item.optional_vars is not None:
                        assign(
                            context_item.optional_vars,
                            None,
                            env,
                            possible_origins(context_item.context_expr, env),
                        )
                env = statements(statement.body, env, credit=credit)
                if terminates(statement.body):
                    break
            elif isinstance(statement, ast.Try | ast.Match):
                # Unsupported control flow must not lend assignments to a
                # later derivation; do not claim witnesses inside these blocks.
                for node in local_nodes(statement):
                    if isinstance(node, ast.Assign | ast.AnnAssign) and node.value is not None:
                        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                        possible = possible_origins(node.value, env)
                        for item in targets:
                            if isinstance(item, ast.Name):
                                bind_origins(
                                    item.id,
                                    merge_origins(named_origins(item.id, env), possible),
                                    env,
                                )
                    if isinstance(node, ast.AugAssign):
                        invalidate_aliases(node.target, env)
                    if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
                        possible = named_origins(node.id, env)
                        env[node.id] = ("unknown", node.id)
                        bind_origins(node.id, possible, env)
                    elif isinstance(node, ast.Subscript | ast.Attribute) and isinstance(
                        node.ctx, ast.Store | ast.Del
                    ):
                        assign(node, None, env)
            # Nested functions/classes belong to their own qualified anchors.
        return env

    args = function.args
    parameters = [*args.posonlyargs, *args.args, *args.kwonlyargs]
    parameters.extend(arg for arg in (args.vararg, args.kwarg) if arg is not None)
    parameter_names = frozenset(arg.arg for arg in parameters)
    statements(
        function.body,
        DerivationEnvironment({arg.arg: ("vector", arg.arg) for arg in parameters}),
    )
    return witnessed


@dataclass
class CompactUnit:
    """The only per-file state retained after parsing: strings and small scalars."""

    path: str
    module: str
    functions: tuple[str, ...]
    dependencies: dict[str, set[str]]
    calls: list[tuple[str, str, tuple[str, ...], int]]
    derived: tuple[str, ...]


def check_ast_budget(tree: ast.AST, limits: Limits) -> None:
    count = 0
    pending = [(tree, 0)]
    while pending:
        node, depth = pending.pop()
        count += 1
        if count > limits.ast_nodes:
            raise DiscoveryLimitError("per-file AST node budget exceeded")
        if depth > limits.ast_depth:
            raise DiscoveryLimitError("per-file AST depth budget exceeded")
        pending.extend((child, depth + 1) for child in ast.iter_child_nodes(node))


def summarize_source(path: str, source: str, budget: SummaryBudget) -> CompactUnit:
    """Parse one bounded file, then discard its source/AST before reading another."""
    limits = budget.limits
    if len(source) > limits.source_bytes or len(source.encode("utf-8")) > limits.source_bytes:
        raise DiscoveryLimitError(f"per-file source byte budget exceeded: {path}")
    module = path.removeprefix("src/").removesuffix(".py").replace("/", ".")
    module = module.removesuffix(".__init__")
    unit = Unit(path, module, ast.parse(source, filename=path), limits)
    check_ast_budget(unit.tree, limits)
    Visitor(unit).visit(unit.tree)
    budget.charge(path, "<module>")
    for function in unit.functions:
        budget.charge(path, function)
    calls = []
    bh_scopes: set[str] = set()
    clip_scopes: set[str] = set()
    for function, call in unit.calls:
        target = unit.resolve(call.func, function)
        procedures = tuple(sorted(literal_procedures(call, target)))
        budget.charge(path, function, target, *procedures)
        calls.append((function, target, procedures, call.lineno))
        if target.rsplit(".", 1)[-1] == "benjamini_hochberg":
            bh_scopes.add(function)
        if target == "numpy.clip":
            clip_scopes.add(function)
    candidates = bh_scopes & clip_scopes
    file_budget = FileDerivationBudget(limits.file_derivation_steps)
    derived = tuple(
        name
        for name in unit.functions
        if name in candidates and visible_by_derivation(unit, name, file_budget=file_budget)
    )
    for name in derived:
        budget.charge(path, name, "BY")
    return CompactUnit(
        path, module, tuple(unit.functions), function_dependencies(unit, budget), calls, derived
    )


def discover(
    sources: Sources,
    *,
    limits: Limits = DEFAULT_LIMITS,
) -> tuple[dict[tuple[str, str], dict[str, Any]], set[tuple[str, str]]]:
    budget = SummaryBudget(limits)
    summaries = []
    entries = sources.items() if isinstance(sources, Mapping) else sources
    for count, (path, source) in enumerate(entries, 1):
        if count > limits.files:
            raise DiscoveryLimitError("source file count budget exceeded")
        summaries.append(summarize_source(path, source, budget))
    dependencies = {
        name: targets for summary in summaries for name, targets in summary.dependencies.items()
    }
    names = forwarding_names(dependencies)
    anchors = {(summary.path, name) for summary in summaries for name in summary.functions}
    anchors.update((summary.path, "<module>") for summary in summaries)
    sites: dict[tuple[str, str], dict[str, Any]] = {}
    for summary in summaries:
        for function, target, inferred, line in summary.calls:
            if not is_seed(target) and target not in names:
                continue
            key = (summary.path, function)
            site = sites.setdefault(key, {"procedures": set(), "calls": [], "definition": False})
            site["procedures"].update(inferred)
            site["calls"].append({"target": target, "procedures": list(inferred), "line": line})
        for function in summary.functions:
            if is_seed(f"{summary.module}.{function}"):
                sites.setdefault((summary.path, function), {"procedures": set(), "calls": []})[
                    "definition"
                ] = True
        for function in summary.derived:
            key = (summary.path, function)
            if key in sites:
                sites[key]["derived_procedures"] = {"BY"}
    return sites, anchors


def nonempty_text(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def reference_errors(value: Any, existing: set[tuple[str, str]], label: str) -> list[str]:
    if not isinstance(value, dict) or not all(
        nonempty_text(value.get(field)) for field in ("path", "function")
    ):
        return [f"{label} needs non-empty path and function"]
    key = (value["path"], value["function"])
    return [] if key in existing else [f"{label} does not exist: {key[0]}::{key[1]}"]


def procedure_list(value: Any) -> bool:
    return isinstance(value, list) and all(
        isinstance(item, str) and item in {"BH", "BY", "none"} for item in value
    )


def family_schema_errors(
    family: Any, existing: set[tuple[str, str]], sites: Mapping[tuple[str, str], Any]
) -> list[str]:
    if not isinstance(family, dict):
        return ["family must be a mapping"]
    errors: list[str] = []
    if missing := REQUIRED - family.keys():
        errors.append(f"missing fields {sorted(missing)}")
    for name in (
        "family_id",
        "area",
        "hypotheses",
        "family_definition",
        "threshold",
        "deliberate_difference",
    ):
        if not nonempty_text(family.get(name)):
            errors.append(f"{name} must be non-empty text")
    primary = family.get("primary_procedure")
    if not isinstance(primary, str) or primary not in {"BH", "BY", "none"}:
        errors.append("primary_procedure must be BH, BY or none")
    if not procedure_list(family.get("comparison_procedure")):
        errors.append("comparison_procedure must be a list of BH, BY or none")
    issues = family.get("open_issue")
    if not isinstance(issues, list) or any(type(issue) is not int for issue in issues):
        errors.append("open_issue must be a list of integers")
    size = family.get("size_source")
    if not isinstance(size, dict):
        errors.append("size_source must be a mapping")
    else:
        if size.get("kind") not in ("executed_inventory", "selected_rows"):
            errors.append("size_source kind must be executed_inventory or selected_rows")
        for name in ("selection_rule", "count_rule"):
            if not nonempty_text(size.get(name)):
                errors.append(f"size_source {name} must be non-empty text")
        if size.get("selection_rule") == family.get("family_definition"):
            errors.append("selection_rule must explain selection, not repeat family_definition")
        if size.get("kind") == "selected_rows" and not nonempty_text(
            size.get("unfiltered_count_rule")
        ):
            errors.append("selected_rows needs an unfiltered_count_rule")
    for name in ("anchors", "consumers"):
        references = family.get(name)
        if not isinstance(references, list) or not references:
            errors.append(f"{name} must be a non-empty list")
            continue
        for reference in references:
            errors.extend(
                reference_errors(reference, existing, "anchor" if name == "anchors" else "consumer")
            )
            if (
                name == "anchors"
                and isinstance(reference, dict)
                and not procedure_list(reference.get("procedures"))
            ):
                errors.append("anchor procedures must be a list of BH, BY or none")
    if "launcher" in family and type(family["launcher"]) is not bool:
        errors.append("launcher must be a boolean")
    if not errors:
        discovered = any(
            (anchor["path"], anchor["function"]) in sites for anchor in family["anchors"]
        )
        if family.get("launcher") is True and discovered:
            errors.append("launcher: true is redundant with a discovered correction anchor")
        elif not discovered and family.get("launcher") is not True:
            errors.append("family needs a discovered correction anchor or explicit launcher: true")
    for name in ("outputs", "policy_rule"):
        value = family.get(name)
        if not isinstance(value, list) or not value:
            errors.append(f"{name} must be a non-empty list")
        elif name == "outputs":
            for output in value:
                if not isinstance(output, dict) or not all(
                    nonempty_text(output.get(key)) for key in ("file", "column")
                ):
                    errors.append("outputs items need non-empty file and column")
        elif any(
            not isinstance(rule, str) or rule not in {f"R{i}" for i in range(1, 10)}
            for rule in value
        ):
            errors.append("policy_rule items must be R1 through R9")
    return errors


def family_procedure_errors(
    family: dict[str, Any], sites: Mapping[tuple[str, str], Any]
) -> list[str]:
    """Check each family's own declarations, without borrowing from other families.

    Arithmetic declarations need a visible derivation at their own anchor.
    Dynamic calls and helper definitions do not provide an exhaustive literal set.
    """
    errors = []
    identity = family["family_id"]
    anchored = {method for anchor in family["anchors"] for method in anchor["procedures"]}
    roles = {family["primary_procedure"], *family["comparison_procedure"]} - {"none"}
    if missing := roles - anchored:
        errors.append(f"{identity}: procedures absent from own anchors: {sorted(missing)}")
    for anchor in family["anchors"]:
        key = (anchor["path"], anchor["function"])
        derived = anchor.get("derived_procedures", {})
        if not isinstance(derived, dict) or any(
            method not in {"BH", "BY"}
            or method not in anchor["procedures"]
            or not nonempty_text(reason)
            for method, reason in derived.items()
        ):
            errors.append(f"{identity}: invalid derived_procedures at {key}")
            continue
        site = sites.get(key)
        if not site or not site["calls"] or any(not call["procedures"] for call in site["calls"]):
            if derived:
                errors.append(f"{identity}: derived procedures need an all-literal site: {key}")
            continue
        literal = {method for call in site["calls"] for method in call["procedures"]}
        visible = site.get("derived_procedures", set())
        if unsupported := derived.keys() - visible:
            errors.append(
                f"{identity}: derived procedures lack visible arithmetic at {key}: "
                f"{sorted(unsupported)}"
            )
        if redundant := derived.keys() & literal:
            errors.append(f"{identity}: derived procedures already literal: {sorted(redundant)}")
        if extra := set(anchor["procedures"]) - literal - (derived.keys() & visible) - {"none"}:
            errors.append(
                f"{identity}: procedures absent from all-literal site {key}: {sorted(extra)}"
            )
    return errors


def normalized_reason(reason: str) -> str:
    normalized = " ".join(reason.casefold().split())
    end = len(normalized)
    while end and (
        normalized[end - 1].isspace() or unicodedata.category(normalized[end - 1]).startswith("P")
    ):
        end -= 1
    return normalized[:end]


def literal_count_errors(
    registry: dict[str, Any], sites: Mapping[tuple[str, str], Any], existing: set[tuple[str, str]]
) -> list[str]:
    records = registry.get("literal_procedure_counts")
    if not isinstance(records, list):
        return ["literal_procedure_counts must be a list"]
    expected: dict[tuple[str, str], dict[str, int]] = {}
    errors: list[str] = []
    for record in records:
        invalid = reference_errors(record, existing, "literal procedure count anchor")
        if invalid:
            errors.extend(invalid)
            continue
        key = (record["path"], record["function"])
        counts = record.get("counts")
        if not isinstance(counts, dict) or any(
            method not in {"BH", "BY", "none"} or type(count) is not int or count < 1
            for method, count in counts.items()
        ):
            errors.append(f"invalid literal procedure counts: {key}")
            continue
        if key in expected:
            errors.append(f"duplicate literal procedure count anchor: {key}")
        expected[key] = counts
        if key not in sites:
            errors.append(f"literal procedure count anchor is not a discovered site: {key}")
    for key, site in sorted(sites.items()):
        actual = dict(
            Counter(procedure for call in site["calls"] for procedure in call["procedures"])
        )
        if key not in expected:
            errors.append(f"missing literal procedure counts: {key[0]}::{key[1]}")
        elif actual != expected[key]:
            errors.append(
                f"literal procedure count mismatch: {key[0]}::{key[1]} "
                f"expected {expected[key]}, found {actual}"
            )
    return errors


def registry_errors(sources: Sources, registry: dict[str, Any]) -> list[str]:
    sites, existing = discover(sources)
    errors = []
    if registry.get("family_ids_are_new") is not True:
        errors.append("family_ids_are_new must explicitly declare these new identifiers")
    declared: dict[tuple[str, str], set[str]] = defaultdict(set)
    reported: dict[tuple[str, str], set[str]] = defaultdict(set)
    identities: set[str] = set()
    bh_reasons: set[str] = set()
    families = registry.get("families")
    if not isinstance(families, list) or not families:
        return [*errors, "families must be a non-empty list"]
    for family in families:
        if invalid := family_schema_errors(family, existing, sites):
            errors.extend(invalid)
            continue
        identity = family["family_id"]
        if identity in identities:
            errors.append(f"duplicate family_id {identity}")
        identities.add(identity)
        errors.extend(family_procedure_errors(family, sites))
        if family["primary_procedure"] == "BH":
            reason = normalized_reason(family["deliberate_difference"])
            if not reason or reason in bh_reasons:
                errors.append(f"duplicate BH-primary reason: {identity}")
            bh_reasons.add(reason)
        for anchor in family.get("anchors", []):
            key = (anchor["path"], anchor["function"])
            declared[key].update(anchor["procedures"])
            for field_name in ("primary_procedure", "comparison_procedure"):
                value = family.get(field_name)
                reported[key].update(value if isinstance(value, list) else [value])
    for anchor in registry.get("infrastructure", []):
        if invalid := reference_errors(anchor, existing, "infrastructure anchor"):
            errors.extend(invalid)
            continue
        if not anchor.get("reason"):
            errors.append(f"infrastructure anchor needs reason: {anchor}")
        key = (anchor["path"], anchor["function"])
        declared[key].update(anchor["procedures"])
        reported[key].update(anchor["procedures"])
    for key in sorted(declared):
        if key not in existing:
            errors.append(f"anchor does not exist: {key[0]}::{key[1]}")
        if missing := declared[key] - reported[key]:
            errors.append(
                f"family procedures omit anchor method: {key[0]}::{key[1]} {sorted(missing)}"
            )
    for key, site in sorted(sites.items()):
        if key not in declared:
            errors.append(f"unregistered correction: {key[0]}::{key[1]}")
        elif missing := site["procedures"] - declared[key]:
            errors.append(f"procedure contradicts code: {key[0]}::{key[1]} needs {sorted(missing)}")
    errors.extend(literal_count_errors(registry, sites, existing))
    return errors


def repository_sources(root: Path, limits: Limits = DEFAULT_LIMITS) -> Iterator[tuple[str, str]]:
    """Read one capped file at a time, including gitignored sources."""
    count = 0
    for directory in (root / "src", root / "scripts"):
        for path in directory.rglob("*.py"):
            if not path.is_file():
                continue
            count += 1
            if count > limits.files:
                raise DiscoveryLimitError("source file count budget exceeded")
            with path.open("rb") as stream:
                raw = stream.read(limits.source_bytes + 1)
            if len(raw) > limits.source_bytes:
                raise DiscoveryLimitError(f"per-file source byte budget exceeded: {path}")
            yield path.relative_to(root).as_posix(), raw.decode("utf-8")


def test_repository_corrections_have_registered_families() -> None:
    assert REGISTRY.is_file(), "missing correction-family registry"
    registry = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))
    assert isinstance(registry, dict)
    assert not (errors := registry_errors(repository_sources(ROOT), registry)), "\n".join(errors)


def example_registry(*anchors: tuple[str, str, list[str]]) -> dict[str, Any]:
    procedures = {method for _, _, methods in anchors for method in methods}
    primary = "BY" if "BY" in procedures else "BH" if "BH" in procedures else "none"
    return {
        "version": 1,
        "family_ids_are_new": True,
        "families": [
            {
                **dict.fromkeys(REQUIRED, "synthetic fixture"),
                "family_id": "synthetic-family",
                "primary_procedure": primary,
                "comparison_procedure": sorted(procedures - {primary}),
                "size_source": {
                    "kind": "executed_inventory",
                    "selection_rule": "Include every input hypothesis.",
                    "count_rule": "Count every input hypothesis, including untestable slots.",
                },
                "outputs": [{"file": "results.tsv", "column": "q"}],
                "consumers": [
                    {"path": path, "function": function} for path, function, _ in anchors
                ],
                "policy_rule": ["R1"],
                "open_issue": [],
                "anchors": [
                    {"path": path, "function": function, "procedures": procedures}
                    for path, function, procedures in anchors
                ],
            }
        ],
        "literal_procedure_counts": [
            {"path": path, "function": function, "counts": dict.fromkeys(procedures, 1)}
            for path, function, procedures in anchors
        ],
    }


def test_incomplete_registry_fails_at_exact_enclosing_function() -> None:
    sources = {"src/x.py": "def screen(p):\n    return false_discovery_control(p, method='by')\n"}
    errors = registry_errors(sources, example_registry())
    assert "unregistered correction: src/x.py::screen" in errors


def test_missing_anchor_fails() -> None:
    assert "anchor does not exist: src/x.py::gone" in registry_errors(
        {"src/x.py": "def screen():\n    pass\n"}, example_registry(("src/x.py", "gone", ["BH"]))
    )


@pytest.mark.parametrize("method", ["bh", "by", "fdr_bh", "fdr_by"])
def test_literal_procedure_contradiction_fails(method: str) -> None:
    sources = {"src/x.py": f"def screen(p):\n    return multipletests(p, method={method!r})\n"}
    wrong = "BY" if method.endswith("bh") else "BH"
    assert any(
        "procedure contradicts code" in error
        for error in registry_errors(sources, example_registry(("src/x.py", "screen", [wrong])))
    )


def test_import_assignment_and_local_forwarding_aliases_are_discovered() -> None:
    sources = {
        "src/x.py": (
            "from scipy.stats import false_discovery_control as fdc\n"
            "correct = fdc\n"
            "def wrapper(p):\n    q = correct(p, method='by')\n    return q\n"
            "def screen(p):\n    return wrapper(p)\n"
        ),
        "scripts/y.py": "import x as local\ndef run(p):\n    return local.wrapper(p)\n",
    }
    errors = registry_errors(sources, example_registry(("src/x.py", "wrapper", ["BY"])))
    assert "unregistered correction: src/x.py::screen" in errors
    assert "unregistered correction: scripts/y.py::run" in errors
    registry = example_registry(
        ("src/x.py", "wrapper", ["BY"]),
        ("src/x.py", "screen", ["BY"]),
        ("scripts/y.py", "run", ["BY"]),
    )
    registry["literal_procedure_counts"][1]["counts"] = {}
    registry["literal_procedure_counts"][2]["counts"] = {}
    assert not registry_errors(sources, registry)


def test_nested_qualified_function_anchor_does_not_accept_parent() -> None:
    sources = {
        "src/x.py": "class Screen:\n    def run(self, p):\n        return step_up(p, method='by')\n"
    }
    errors = registry_errors(sources, example_registry(("src/x.py", "run", ["BY"])))
    assert "unregistered correction: src/x.py::Screen.run" in errors
    assert "anchor does not exist: src/x.py::run" in errors


def test_hand_written_helper_and_ignored_script_are_not_omitted(tmp_path: Path) -> None:
    (tmp_path / "src").mkdir()
    (tmp_path / "scripts").mkdir()
    (tmp_path / ".gitignore").write_text("scripts/hidden.py\n", encoding="utf-8")
    (tmp_path / "scripts/hidden.py").write_text(
        "def benjamini_hochberg(p):\n    return sorted(p)\n"
        "def run(p):\n    return benjamini_hochberg(p)\n",
        encoding="utf-8",
    )
    errors = registry_errors(repository_sources(tmp_path), example_registry())
    assert "unregistered correction: scripts/hidden.py::benjamini_hochberg" in errors
    assert "unregistered correction: scripts/hidden.py::run" in errors


def test_dynamic_method_is_not_falsely_inferred() -> None:
    sources = {
        "src/x.py": "def screen(p, method):\n    return false_discovery_control(p, method=method)\n"
    }
    sites, _ = discover(sources)
    assert sites[("src/x.py", "screen")]["procedures"] == set()
    registry = example_registry(("src/x.py", "screen", ["BH", "BY"]))
    registry["literal_procedure_counts"][0]["counts"] = {}
    assert not registry_errors(sources, registry)


def test_top_level_procedure_change_cannot_hide_behind_correct_anchor() -> None:
    sources = {"src/x.py": "def screen(p):\n    return false_discovery_control(p, method='by')\n"}
    registry = example_registry(("src/x.py", "screen", ["BY"]))
    registry["families"][0]["primary_procedure"] = "BH"
    assert any(
        "family procedures omit anchor method" in error
        for error in registry_errors(sources, registry)
    )


def test_local_alias_does_not_shadow_import_in_another_function() -> None:
    sources = {
        "src/x.py": (
            "from scipy.stats import false_discovery_control as fdc\n"
            "def unrelated(x):\n    fdc = x.identity\n    return fdc(x)\n"
            "def screen(p):\n    return fdc(p, method='by')\n"
        )
    }
    sites, _ = discover(sources)
    assert ("src/x.py", "screen") in sites
    assert sites[("src/x.py", "screen")]["procedures"] == {"BY"}
    assert ("src/x.py", "unrelated") not in sites


def test_correction_wrapper_without_return_is_not_omitted() -> None:
    sources = {
        "src/x.py": (
            "def _correct(p):\n    p[:] = false_discovery_control(p, method='by')\n"
            "def screen(p):\n    _correct(p)\n"
        )
    }
    errors = registry_errors(sources, example_registry(("src/x.py", "_correct", ["BY"])))
    assert "unregistered correction: src/x.py::screen" in errors


def test_nested_forwarding_wrapper_is_resolved_by_qualified_name() -> None:
    sources = {
        "src/x.py": (
            "def outer(p):\n"
            "    def wrap(values):\n        return false_discovery_control(values, method='bh')\n"
            "    return wrap(p)\n"
        )
    }
    sites, _ = discover(sources)
    assert ("src/x.py", "outer.wrap") in sites
    assert ("src/x.py", "outer") in sites


def test_attribute_alias_cycle_terminates_without_hiding_other_corrections() -> None:
    sources = {
        "src/x.py": (
            "def screen(p):\n    p = p.values\n    p.items()\n"
            "    return false_discovery_control(p, method='by')\n"
        )
    }
    sites, _ = discover(sources)
    assert sites[("src/x.py", "screen")]["procedures"] == {"BY"}


@pytest.mark.parametrize(
    ("limits", "message"),
    [
        (Limits(source_bytes=20), "source byte budget"),
        (Limits(ast_nodes=3), "AST node budget"),
        (Limits(ast_depth=1), "AST depth budget"),
        (Limits(files=0), "file count budget"),
        (Limits(summary_items=0), "summary item budget"),
        (Limits(summary_bytes=1), "summary byte budget"),
        (Limits(name_characters=10), "name length budget"),
    ],
)
def test_discovery_refuses_exceeded_budgets(limits: Limits, message: str) -> None:
    source = "def screen(p):\n    return false_discovery_control(p, method='by')\n"
    with pytest.raises(DiscoveryLimitError, match=message):
        discover({"src/x.py": source}, limits=limits)


def test_long_alias_chain_refuses_before_unbounded_expansion() -> None:
    source = (
        "a = b\nb = c\nc = false_discovery_control\ndef screen(p):\n    return a(p, method='by')\n"
    )
    with pytest.raises(DiscoveryLimitError, match="alias expansion step budget"):
        discover({"src/x.py": source}, limits=Limits(alias_steps=2))


def test_alias_expansion_length_is_checked_before_concatenation() -> None:
    source = "a = lengthy\ndef screen(p):\n    return a.attribute(p)\n"
    with pytest.raises(DiscoveryLimitError, match="alias expansion length budget"):
        discover({"src/x.py": source}, limits=Limits(name_characters=12))


def test_repository_reader_is_lazy_and_refuses_oversized_files(tmp_path: Path) -> None:
    directory = tmp_path / "src"
    directory.mkdir()
    source = directory / "x.py"
    source.write_text("pass\n", encoding="utf-8")
    entries = repository_sources(tmp_path, Limits(source_bytes=10))
    # If the reader eagerly read the file on construction, this replacement
    # would not trigger the cap when the first item is actually requested.
    source.write_text("#" * 11, encoding="utf-8")
    with pytest.raises(DiscoveryLimitError, match="source byte budget"):
        next(entries)


def test_compact_file_summary_retains_no_ast_or_source() -> None:
    source = "def screen(p):\n    return false_discovery_control(p, method='by')\n"
    summary = summarize_source("src/x.py", source, SummaryBudget(DEFAULT_LIMITS))
    assert summary.functions == ("screen",)
    assert summary.calls == [("screen", "false_discovery_control", ("BY",), 2)]
    assert summary.dependencies == {"x.screen": {"false_discovery_control"}}
    assert summary.derived == ()
    assert set(vars(summary)) == {"path", "module", "functions", "dependencies", "calls", "derived"}


def test_streamed_source_input_matches_small_mapping_input() -> None:
    source = "def screen(p):\n    return false_discovery_control(p, method='by')\n"
    items = [("src/x.py", source)]
    assert discover(iter(items)) == discover(dict(items))
