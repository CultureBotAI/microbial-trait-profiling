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
from collections import defaultdict
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


DEFAULT_LIMITS = Limits()
Sources = Mapping[str, str] | Iterable[tuple[str, str]]


class DiscoveryLimitError(ValueError):
    """Refuse a scan before its declared static-work budget is exceeded."""


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


@dataclass
class CompactUnit:
    """The only per-file state retained after parsing: strings and small scalars."""

    path: str
    module: str
    functions: tuple[str, ...]
    dependencies: dict[str, set[str]]
    calls: list[tuple[str, str, tuple[str, ...], int]]


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
    for function, call in unit.calls:
        target = unit.resolve(call.func, function)
        procedures = tuple(sorted(literal_procedures(call, target)))
        budget.charge(path, function, target, *procedures)
        calls.append((function, target, procedures, call.lineno))
    return CompactUnit(
        path, module, tuple(unit.functions), function_dependencies(unit, budget), calls
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
    return sites, anchors


def registry_errors(sources: Sources, registry: dict[str, Any]) -> list[str]:
    sites, existing = discover(sources)
    errors = []
    if registry.get("family_ids_are_new") is not True:
        errors.append("family_ids_are_new must explicitly declare these new identifiers")
    declared: dict[tuple[str, str], set[str]] = defaultdict(set)
    reported: dict[tuple[str, str], set[str]] = defaultdict(set)
    identities: set[str] = set()
    for family in registry.get("families", []):
        if missing := REQUIRED - family.keys():
            errors.append(f"family {family.get('family_id')} missing {sorted(missing)}")
        identity = family.get("family_id")
        if identity in identities:
            errors.append(f"duplicate family_id {identity}")
        identities.add(identity)
        for anchor in family.get("anchors", []):
            key = (anchor["path"], anchor["function"])
            declared[key].update(anchor["procedures"])
            for field_name in ("primary_procedure", "comparison_procedure"):
                value = family.get(field_name)
                reported[key].update(value if isinstance(value, list) else [value])
    for anchor in registry.get("infrastructure", []):
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
    return {
        "version": 1,
        "family_ids_are_new": True,
        "families": [
            {
                **dict.fromkeys(REQUIRED, "synthetic fixture"),
                "family_id": "synthetic-family",
                "primary_procedure": ["BH", "BY"],
                "comparison_procedure": "none",
                "anchors": [
                    {"path": path, "function": function, "procedures": procedures}
                    for path, function, procedures in anchors
                ],
            }
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
    assert not registry_errors(
        sources,
        example_registry(
            ("src/x.py", "wrapper", ["BY"]),
            ("src/x.py", "screen", ["BY"]),
            ("scripts/y.py", "run", ["BY"]),
        ),
    )


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
    assert not registry_errors(sources, example_registry(("src/x.py", "screen", ["BH", "BY"])))


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
    assert set(vars(summary)) == {"path", "module", "functions", "dependencies", "calls"}


def test_streamed_source_input_matches_small_mapping_input() -> None:
    source = "def screen(p):\n    return false_discovery_control(p, method='by')\n"
    items = [("src/x.py", source)]
    assert discover(iter(items)) == discover(dict(items))
