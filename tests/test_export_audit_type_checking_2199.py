"""Guard classification through real definition and re-export audits (#2199)."""

from __future__ import annotations

import ast
import textwrap

import pytest

from mtp.export_audit import export_omissions, public_definitions


@pytest.mark.parametrize(
    ("before", "guard", "after", "expected"),
    [
        ("from typing import TYPE_CHECKING as TC", "TC", "", {"in_else"}),
        ("import typing as t", "t.TYPE_CHECKING", "", {"in_else"}),
        ("from typing import TYPE_CHECKING as checking", "checking", "", {"in_else"}),
        ("from typing import TYPE_CHECKING as TC", "TC and enabled", "", {"in_else"}),
        ("from typing import TYPE_CHECKING as TC", "enabled and TC", "", {"in_else"}),
        ("import typing as t", "t.TYPE_CHECKING and enabled", "", {"in_else"}),
        ("from typing import TYPE_CHECKING", "TYPE_CHECKING and enabled", "", {"in_else"}),
        ("import typing", "typing.TYPE_CHECKING and enabled", "", {"in_else"}),
        ("from typing import TYPE_CHECKING as TC", "not TC", "", {"in_body"}),
        ("import typing as t", "not t.TYPE_CHECKING", "", {"in_body"}),
        ("from typing import TYPE_CHECKING as TC", "not (TC and enabled)", "", {"in_body"}),
        ("from typing import TYPE_CHECKING as TC", "TC and (first or second)", "", {"in_else"}),
        (
            "from typing import TYPE_CHECKING as TC\nimport typing as t",
            "(TC or enabled) and t.TYPE_CHECKING",
            "",
            {"in_else"},
        ),
        (
            "from typing import TYPE_CHECKING as TC\ndef helper():\n    TC = True",
            "TC",
            "",
            {"in_else"},
        ),
        (
            "import typing as t\nclass Local:\n    t = object()",
            "t.TYPE_CHECKING",
            "",
            {"in_else"},
        ),
        # A runtime alternative must keep both possible branches in the static audit.
        ("from typing import TYPE_CHECKING as TC", "TC or enabled", "", {"in_body", "in_else"}),
        ("from typing import TYPE_CHECKING as TC", "enabled or TC", "", {"in_body", "in_else"}),
        (
            "from typing import TYPE_CHECKING as TC",
            "not TC and enabled",
            "",
            {"in_body", "in_else"},
        ),
        ("import typing as t", "not (t.TYPE_CHECKING or enabled)", "", {"in_body", "in_else"}),
        # Added aliases need actual, preceding, unambiguous module import bindings.
        ("TC = True", "TC", "", {"in_body", "in_else"}),
        ("t = object()", "t.TYPE_CHECKING", "", {"in_body", "in_else"}),
        ("from unrelated import TYPE_CHECKING as TC", "TC", "", {"in_body", "in_else"}),
        ("from .typing import TYPE_CHECKING as TC", "TC", "", {"in_body", "in_else"}),
        ("from typing import Any as TC", "TC", "", {"in_body", "in_else"}),
        (
            "def helper():\n    from typing import TYPE_CHECKING as TC",
            "TC",
            "",
            {"in_body", "in_else"},
        ),
        ("class Local:\n    import typing as t", "t.TYPE_CHECKING", "", {"in_body", "in_else"}),
        ("", "TC", "from typing import TYPE_CHECKING as TC", {"in_body", "in_else"}),
        (
            "if enabled:\n    from typing import TYPE_CHECKING as TC",
            "TC",
            "",
            {"in_body", "in_else"},
        ),
        ("from typing import TYPE_CHECKING as TC\nTC = True", "TC", "", {"in_body", "in_else"}),
        # Ambiguous aliases are conservative even if the other binding occurs later.
        ("from typing import TYPE_CHECKING as TC", "TC", "TC = True", {"in_body", "in_else"}),
        (
            "import typing as t\nt.TYPE_CHECKING = True",
            "t.TYPE_CHECKING",
            "",
            {"in_body", "in_else"},
        ),
        (
            "from typing import TYPE_CHECKING as TC",
            "(TC := enabled) and TC",
            "",
            {"in_body", "in_else"},
        ),
        (
            "from typing import TYPE_CHECKING as TC\nfor TC in values:\n    pass",
            "TC",
            "",
            {"in_body", "in_else"},
        ),
        (
            "from typing import TYPE_CHECKING as TC\nwith manager() as TC:\n    pass",
            "TC",
            "",
            {"in_body", "in_else"},
        ),
        (
            "from typing import TYPE_CHECKING as TC\n"
            "try:\n    work()\nexcept Exception as TC:\n    pass",
            "TC",
            "",
            {"in_body", "in_else"},
        ),
        (
            "from typing import TYPE_CHECKING as TC\nmatch subject:\n    case TC:\n        pass",
            "TC",
            "",
            {"in_body", "in_else"},
        ),
        ("from typing import TYPE_CHECKING as TC\ndel TC", "TC", "", {"in_body", "in_else"}),
        (
            "from typing import TYPE_CHECKING as TC\nfrom other import *",
            "TC",
            "",
            {"in_body", "in_else"},
        ),
    ],
)
def test_type_checking_guard_definitions_and_reexports(
    before: str, guard: str, after: str, expected: set[str]
) -> None:
    """Both audit paths select the right branches without executing fixture modules."""
    _assert_guard_audits(before, guard, after, expected)


def _assert_guard_audits(before: str, guard: str, after: str, expected: set[str]) -> None:
    definitions = (
        f"{before}\nif {guard}:\n    def in_body(): pass\nelse:\n    def in_else(): pass\n{after}\n"
    )
    assert public_definitions(ast.parse(definitions)) & {"in_body", "in_else"} == expected

    reexports = (
        f"{before}\nif {guard}:\n    from .origin import shared as in_body\n"
        f"else:\n    from .origin import shared as in_else\n{after}\n__all__ = []\n"
    )
    sources = {
        "mtp.origin": "shared = 1\n__all__ = ['shared']\n",
        "mtp.fixture": reexports,
        "mtp.consumer": "from mtp.fixture import in_body, in_else\n",
    }
    assert export_omissions(sources, set()) == {"mtp.fixture": sorted(expected)}


@pytest.mark.parametrize(
    ("before", "guard"),
    [
        ("import typing as t\nt.TYPE_CHECKING, other = True, 0", "t.TYPE_CHECKING"),
        ("import typing as t\n[t.TYPE_CHECKING, other] = [True, 0]", "t.TYPE_CHECKING"),
        ("import typing as t\nother, *t.TYPE_CHECKING = [0, True]", "t.TYPE_CHECKING"),
        ("from typing import TYPE_CHECKING as TC\nf = lambda value=(TC := True): None", "TC"),
        ("from typing import TYPE_CHECKING as TC\nf = lambda *, value=(TC := True): None", "TC"),
        (
            "from typing import TYPE_CHECKING as TC\nclass Mutate:\n    global TC\n    TC = True",
            "TC",
        ),
        (
            "from typing import TYPE_CHECKING as TC\n"
            "def mutate():\n    global TC\n    TC = True\nmutate()",
            "TC",
        ),
        ("import typing as t\nclass Mutate:\n    t.TYPE_CHECKING = True", "t.TYPE_CHECKING"),
        (
            "import typing as first\nimport typing as second\nfirst.TYPE_CHECKING = True",
            "second.TYPE_CHECKING",
        ),
        (
            "import typing as t\nt.TYPE_CHECKING = True\nfrom typing import TYPE_CHECKING as TC",
            "TC",
        ),
        (
            "from typing import TYPE_CHECKING as TC\n"
            "TYPE_CHECKING = True\nif TYPE_CHECKING:\n    TC = True",
            "TC",
        ),
    ],
)
def test_ambiguous_alias_mutations_keep_possible_runtime_branches(before: str, guard: str) -> None:
    """Binding collection must see compound targets, eager defaults and shared namespaces."""
    _assert_guard_audits(before, guard, "", {"in_body", "in_else"})


def test_lambda_body_local_binding_does_not_shadow_module_alias() -> None:
    """A lambda body is its own scope, unlike its immediately evaluated defaults."""
    _assert_guard_audits(
        "from typing import TYPE_CHECKING as TC\nf = lambda: (TC := True)",
        "TC",
        "",
        {"in_else"},
    )


def test_type_checking_aliases_reach_nested_module_blocks() -> None:
    """Nested runtime blocks inherit aliases; nested lexical scopes do not supply them."""
    source = textwrap.dedent(
        """\
        from typing import TYPE_CHECKING as TC
        import typing as t
        if enabled:
            if TC and supported:
                type_only = 1
            else:
                runtime_if = 1
        try:
            if not t.TYPE_CHECKING:
                runtime_try = 1
            else:
                type_only_else = 1
        except Exception:
            fallback = 1
        """
    )
    assert public_definitions(ast.parse(source)) == {"runtime_if", "runtime_try", "fallback"}


@pytest.mark.parametrize("kind", ["definition", "reexport"])
@pytest.mark.parametrize(
    ("before", "guard"),
    [
        ("import typing as t\nsetattr(t, 'TYPE_CHECKING', True)", "t.TYPE_CHECKING"),
        ("import typing as t\nt.__dict__.update(TYPE_CHECKING=True)", "t.TYPE_CHECKING"),
        (
            "import typing as t\nsetattr(t, 'TYPE_CHECKING', True)\n"
            "from typing import TYPE_CHECKING as TC",
            "TC",
        ),
        (
            "class Mutate:\n    import typing as t\n    t.TYPE_CHECKING = True\n"
            "from typing import TYPE_CHECKING as TC",
            "TC",
        ),
        (
            "def mutate():\n    import typing as t\n    t.TYPE_CHECKING = True\n"
            "mutate()\nfrom typing import TYPE_CHECKING as TC",
            "TC",
        ),
        (
            "import typing as t\nt.__dict__.update({'TYPE_CHECKING': True})",
            "t.TYPE_CHECKING",
        ),
        (
            "import typing as t\ndict.update(t.__dict__, TYPE_CHECKING=True)",
            "t.TYPE_CHECKING",
        ),
        (
            "import typing as t\nt.__dict__.__setitem__('TYPE_CHECKING', True)",
            "t.TYPE_CHECKING",
        ),
        (
            "import builtins\nimport typing as t\nbuiltins.setattr(t, 'TYPE_CHECKING', True)",
            "t.TYPE_CHECKING",
        ),
        (
            "class Mutate:\n    import typing as t\n    setattr(t, 'TYPE_CHECKING', True)\n"
            "from typing import TYPE_CHECKING as TC",
            "TC",
        ),
    ],
    ids=[
        "setattr",
        "namespace-update",
        "copied-after-mutation",
        "class-import",
        "function-import",
        "update-mapping",
        "unbound-update",
        "namespace-setitem",
        "qualified-setattr",
        "class-import-setattr",
    ],
)
def test_namespace_mutation_cannot_hide_the_only_consumed_export(
    kind: str, before: str, guard: str
) -> None:
    """F1: visible typing mutations must not make the omission audit falsely clean."""
    body = (
        "def runtime_export(): pass"
        if kind == "definition"
        else "from .origin import shared as runtime_export"
    )
    source = f"{before}\nif {guard}:\n    {body}\n__all__ = []\n"
    sources = {
        "mtp.origin": "shared = 1\n__all__ = ['shared']\n",
        "mtp.fixture": source,
        "mtp.consumer": "from mtp.fixture import runtime_export\n",
    }
    assert export_omissions(sources, set()) == {"mtp.fixture": ["runtime_export"]}
    if kind == "definition":
        assert "runtime_export" in public_definitions(ast.parse(source))


def test_namespace_deletion_invalidates_alias_inference() -> None:
    """Deletion need not produce a runtime export; its guard is no longer known false."""
    _assert_guard_audits(
        "import typing as t\ndelattr(t, 'TYPE_CHECKING')",
        "t.TYPE_CHECKING",
        "",
        {"in_body", "in_else"},
    )


@pytest.mark.parametrize(
    "before",
    [
        "def helper():\n    import typing as t",
        "class Local:\n    import typing as t",
        "unrelated = object()\nsetattr(unrelated, 'TYPE_CHECKING', True)",
        "unrelated = {}\nunrelated.update(TYPE_CHECKING=True)",
    ],
)
def test_unmutated_typing_alias_remains_type_only(before: str) -> None:
    """Local imports and mutations of unrelated objects do not poison a stable alias."""
    _assert_guard_audits(f"from typing import TYPE_CHECKING as TC\n{before}", "TC", "", {"in_else"})


@pytest.mark.parametrize("kind", ["definition", "reexport"])
@pytest.mark.parametrize(
    ("before", "guard"),
    [
        (
            "import operator\nimport typing as t\n"
            "operator.ior(t.__dict__, {'TYPE_CHECKING': True})",
            "t.TYPE_CHECKING",
        ),
        (
            "from operator import ior\nimport typing as t\n"
            "ior(t.__dict__, {'TYPE_CHECKING': True})",
            "t.TYPE_CHECKING",
        ),
        (
            "from operator import ior as combine\nimport typing as t\n"
            "combine(t.__dict__, {'TYPE_CHECKING': True})",
            "t.TYPE_CHECKING",
        ),
        ("import typing as t\nexternal_mutator(t)", "t.TYPE_CHECKING"),
        ("import typing as t\nexternal_mutator(t.__dict__)", "t.TYPE_CHECKING"),
        ("import typing as t\nt.external_mutator()", "t.TYPE_CHECKING"),
        (
            "from typing import __dict__ as ns\nns['TYPE_CHECKING'] = True\n"
            "from typing import TYPE_CHECKING as TC",
            "TC",
        ),
        (
            "class Mutate:\n    from typing import __dict__ as ns\n"
            "    ns['TYPE_CHECKING'] = True\nfrom typing import TYPE_CHECKING as TC",
            "TC",
        ),
        (
            "import operator\nimport typing as t\n"
            "operator.ior((ns := t.__dict__), {'TYPE_CHECKING': True})",
            "t.TYPE_CHECKING",
        ),
        (
            "import operator\nimport typing as t\n"
            "operator.ior(*(t.__dict__, {'TYPE_CHECKING': True}))",
            "t.TYPE_CHECKING",
        ),
        ("import typing as t\n(ns := t).TYPE_CHECKING = True", "t.TYPE_CHECKING"),
        (
            "import operator\nimport typing as t\noperator.ior(vars(t), {'TYPE_CHECKING': True})",
            "t.TYPE_CHECKING",
        ),
    ],
    ids=[
        "operator-ior",
        "imported-ior",
        "renamed-ior",
        "unknown-call-module",
        "unknown-call-dict",
        "unknown-namespace-call",
        "imported-namespace-write",
        "class-imported-namespace-write",
        "named-expression-argument",
        "unpacked-arguments",
        "named-expression-write",
        "vars-argument",
    ],
)
def test_exposed_typing_namespace_keeps_runtime_possibilities(
    kind: str, before: str, guard: str
) -> None:
    """F1b: call spelling and expression wrappers cannot imply a known-false guard."""
    body = (
        "def runtime_export(): pass"
        if kind == "definition"
        else "from .origin import shared as runtime_export"
    )
    source = f"{before}\nif {guard}:\n    {body}\n__all__ = []\n"
    sources = {
        "mtp.origin": "shared = 1\n__all__ = ['shared']\n",
        "mtp.fixture": source,
        "mtp.consumer": "from mtp.fixture import runtime_export\n",
    }
    assert export_omissions(sources, set()) == {"mtp.fixture": ["runtime_export"]}
    if kind == "definition":
        assert "runtime_export" in public_definitions(ast.parse(source))


@pytest.mark.parametrize(
    "reads",
    ["t.TYPE_CHECKING\nt.__dict__", "unrelated = object()\nexternal_mutator(unrelated)"],
)
def test_attribute_reads_and_unrelated_calls_preserve_typing_alias(reads: str) -> None:
    """Plain attribute reads do not expose the namespace to a possible mutator."""
    _assert_guard_audits(f"import typing as t\n{reads}", "t.TYPE_CHECKING", "", {"in_else"})


@pytest.mark.parametrize("kind", ["definition", "reexport"])
@pytest.mark.parametrize(
    ("before", "guard"),
    [
        ("import typing as t\nsetattr(*(), t, 'TYPE_CHECKING', True)", "t.TYPE_CHECKING"),
        (
            "from typing import __dict__ as ns\nns |= {'TYPE_CHECKING': True}\n"
            "from typing import TYPE_CHECKING as TC",
            "TC",
        ),
        ("import typing as t\nexternal_mutator(module=t)", "t.TYPE_CHECKING"),
    ],
    ids=["empty-starred-prefix", "imported-namespace-augmented-write", "keyword-exposure"],
)
def test_indirect_argument_positions_and_augmented_namespace_writes(
    kind: str, before: str, guard: str
) -> None:
    """The first runtime argument and an in-place store need no mutator-name inference."""
    body = (
        "def runtime_export(): pass"
        if kind == "definition"
        else "from .origin import shared as runtime_export"
    )
    source = f"{before}\nif {guard}:\n    {body}\n__all__ = []\n"
    sources = {
        "mtp.origin": "shared = 1\n__all__ = ['shared']\n",
        "mtp.fixture": source,
        "mtp.consumer": "from mtp.fixture import runtime_export\n",
    }
    assert export_omissions(sources, set()) == {"mtp.fixture": ["runtime_export"]}
    if kind == "definition":
        assert "runtime_export" in public_definitions(ast.parse(source))


@pytest.mark.parametrize("kind", ["definition", "reexport"])
@pytest.mark.parametrize(
    "before",
    [
        "import functools\nimport operator\nimport typing as t\n"
        "functools.reduce(operator.ior, [t.__dict__, {'TYPE_CHECKING': True}])",
        "import typing as t\n"
        "(lambda ignored, module: setattr(module, 'TYPE_CHECKING', True))(None, t)",
        "import typing as t\ndef mutate(value, module):\n"
        "    setattr(module, 'TYPE_CHECKING', value)\nmutate(True, t)",
        "import functools\nimport operator\nimport typing as t\n"
        "functools.reduce(operator.ior, *([t.__dict__, {'TYPE_CHECKING': True}],))",
        "import typing as t\nvalue = t.cast(int, setattr(t, 'TYPE_CHECKING', True))",
        "import typing as t\nT = t.TypeVar('T', bound=(setattr(t, 'TYPE_CHECKING', True), int)[1])",
        "import typing as t\nvalue = t.cast(int, 1)\nt.TYPE_CHECKING = True",
        "import typing as t\nt.TYPE_CHECKING = True\nT = t.TypeVar('T')",
        "import typing as t\nclass Local:\n"
        "    class t:\n        @staticmethod\n        def cast(*args):\n"
        "            setattr(__import__('typing'), 'TYPE_CHECKING', True)\n"
        "    t.cast(int, 1)",
        "import typing as t\nclass Fake:\n    def cast(self, *args):\n"
        "        setattr(__import__('typing'), 'TYPE_CHECKING', True)\n"
        "def mutate(t):\n    t.cast(int, 1)\nmutate(Fake())",
        "import typing as t\nclass Fake:\n    def cast(self, *args):\n"
        "        setattr(__import__('typing'), 'TYPE_CHECKING', True)\n"
        "t = Fake()\nt.cast(int, 1)",
        "import typing as t\nclass Bound:\n    def __bool__(self):\n"
        "        setattr(__import__('typing'), 'TYPE_CHECKING', True)\n"
        "        return False\nT = t.TypeVar('T', covariant=Bound())",
    ],
    ids=[
        "later-reduce-argument",
        "later-lambda-argument",
        "later-named-argument",
        "later-starred-argument",
        "nested-cast-mutation",
        "nested-typevar-mutation",
        "mutation-after-cast",
        "mutation-before-typevar",
        "class-shadowed-cast",
        "parameter-shadowed-cast",
        "rebound-cast-receiver",
        "typevar-callback",
    ],
)
def test_later_namespace_arguments_and_typing_call_mutations(kind: str, before: str) -> None:
    """F1c: all argument positions and unsafe lookalikes retain the sole runtime export."""
    # Copy after the mutation so even rebinding t cannot hide the changed flag.
    source = _fixback3_provider(kind, f"{before}\nfrom typing import TYPE_CHECKING as TC", "TC")
    sources = {
        "mtp.origin": "shared = 1\n__all__ = ['shared']\n",
        "mtp.fixture": source,
        "mtp.consumer": "from mtp.fixture import runtime_export\n",
    }
    assert export_omissions(sources, set()) == {"mtp.fixture": ["runtime_export"]}
    if kind == "definition":
        assert "runtime_export" in public_definitions(ast.parse(source))


@pytest.mark.parametrize("kind", ["definition", "reexport"])
@pytest.mark.parametrize("call", ["value = t.cast(int, 1)", "T = t.TypeVar('T')"])
@pytest.mark.parametrize("nested", [False, True], ids=["module-call", "function-call"])
def test_ordinary_typing_calls_preserve_type_only_provider_and_consumer(
    kind: str, call: str, nested: bool
) -> None:
    """F2: ordinary cast/simple TypeVar calls do not invent missing runtime exports."""
    before = "import typing as t\n"
    before += f"def helper():\n    {call}\nhelper()" if nested else call
    source = _fixback3_provider(kind, before, "t.TYPE_CHECKING")
    sources = {
        "mtp.origin": "shared = 1\n__all__ = ['shared']\n",
        "mtp.fixture": source,
        "mtp.consumer": "from typing import TYPE_CHECKING\n"
        "if TYPE_CHECKING:\n    from mtp.fixture import runtime_export\n",
    }
    assert export_omissions(sources, set()) == {}
    if kind == "definition":
        assert "runtime_export" not in public_definitions(ast.parse(source))


def _fixback3_provider(kind: str, before: str, guard: str) -> str:
    body = (
        "def runtime_export(): pass"
        if kind == "definition"
        else "from .origin import shared as runtime_export"
    )
    return f"{before}\nif {guard}:\n    {body}\n__all__ = []\n"
