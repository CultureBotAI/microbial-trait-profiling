"""Explicit conservative refusals and survivor boundaries (#1915, #1928, #1929).

The namespace-selector witness is already in review_1904. The named-cache free
bit remains defensive: no source-level counterexample for removing it is
established here. In normal invalidation, materializing a free origin is followed
by an unknown fact before another lookup; explicit alias assignments carry a
different side-reference set. No named-cache mutant kill is claimed.

Merge-cache retention is witnessed by uncached container copies (#1969). Slice
and .copy() results can lose their only strong reference without retention;
a later merge may then reuse their identity and return stale origins. Each test
requires refusal of the write and credit of its no-write control, independently
of allocator behaviour. The demonstrated retention-mutant kill depends on the
CPython allocator; the refusal contract does not.

The parameter/self writes and NumPy scalar dtype below are deliberately refused,
even when the particular runtime objects are independent. Refusal is loss of the
local arithmetic BY witness, not refusal of the source's literal BH registration.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from textwrap import indent

import pytest

import test_multiple_testing_registry as guard
import test_multiple_testing_registry_alias_forms as aliases
import test_multiple_testing_registry_review_1904 as review


@pytest.fixture(autouse=True)
def assigned_modules() -> None:
    for module in (guard, aliases, review):
        assert module.__file__ is not None
        expected = Path(__file__).with_name(f"{module.__name__}.py").resolve()
        assert Path(module.__file__).resolve() == expected


@pytest.mark.parametrize("scope", ["module", "closure"])
@pytest.mark.parametrize(
    ("operation", "output"),
    [
        ('out["m"] = m\n', "{}"),
        ("del out[0]\n", "[0.5, 0.6]"),
        ("out.m = m\n", "types.SimpleNamespace()"),
    ],
)
def test_unrelated_parameter_write_conservatively_refuses_free_extent(
    scope: review.FreeScope, operation: str, output: str
) -> None:
    """No call-site proof distinguishes a separate output object from the free p."""
    control, anchor = review.parameter_write_source(scope, "observed = out\n")
    changed, _ = review.parameter_write_source(scope, operation)
    # The existing helper models real caller aliasing; here the caller supplies
    # a distinct, type-appropriate output object. Both scopes include the call.
    assert control.count("screen(p)") == changed.count("screen(p)") == 1
    control = control.replace("screen(p)", f"screen({output})")
    changed = changed.replace("screen(p)", f"screen({output})")
    assert aliases.derived(control) == (anchor,)
    assert aliases.derived(changed) == ()


def test_self_attribute_write_conservatively_refuses_module_vector_extent() -> None:
    """self is another incoming parameter; the guard does not infer its type."""
    header = aliases.IMPORTS + f"p = {aliases.VECTOR}\nclass Screen:\n"
    for operation, expected in (("observed = self\n", ("Screen.run",)), ("self.m = m\n", ())):
        method = "def run(self):\n" + indent(aliases.EXTENT + operation + aliases.RESULT, "    ")
        text = header + indent(method, "    ") + "Screen().run()\n"
        assert aliases.derived(text) == expected


@pytest.mark.parametrize("allocator", ["ones", "zeros"])
@pytest.mark.parametrize("dtype", ["np.float64", "np.int64"])
@pytest.mark.parametrize("scope", aliases.SCOPES)
def test_numpy_scalar_dtype_post_bh_write_is_a_documented_refusal(
    allocator: str, dtype: str, scope: aliases.Scope
) -> None:
    """Only default/built-in numeric dtype yields the fixed-size array shortcut."""
    prefix = aliases.EXTENT + "q = benjamini_hochberg(p)\n"
    end = "return np.clip(q * h, 0, 1)\n"
    allocation = f"allocated = np.{allocator}(m, dtype={dtype})\n"
    control, anchor = aliases.source(scope, prefix + allocation + end)
    changed, _ = aliases.source(scope, prefix + allocation + "allocated[0] = 1\n" + end)
    builtin, _ = aliases.source(
        scope, prefix + allocation.replace(dtype, "float") + "allocated[0] = 1\n" + end
    )
    assert aliases.derived(control) == (anchor,)
    assert aliases.derived(builtin) == (anchor,)
    assert aliases.derived(changed) == ()


# Uncached container copies expose identity reuse if merge inputs are not retained.
# The refusal oracle is portable; killing a retention mutant depends on allocation.
MERGE_COPY_SETUP = (
    "zeta = [0.5]\nother = [None]\nbox1 = [zeta]\nbox2 = [p]\nw1 = box1\nw2 = box2\nw3 = other\n"
)
MERGE_COPY_FORMS = {
    "slice": "c1 = box1[:] if m else other\nc2 = box2[:] if m else other\n",
    "copy-method": "c1 = box1.copy() if m else other\nc2 = box2.copy() if m else other\n",
}


@pytest.mark.parametrize("scope", aliases.SCOPES)
@pytest.mark.parametrize("form", sorted(MERGE_COPY_FORMS))
def test_copied_container_write_after_repeated_merge_is_refused(
    scope: aliases.Scope, form: str
) -> None:
    """A later copied holder must retain its possible reference to the BH vector."""
    prefix = aliases.EXTENT + MERGE_COPY_SETUP + MERGE_COPY_FORMS[form]
    control, anchor = aliases.source(scope, prefix + aliases.RESULT)
    changed, _ = aliases.source(scope, prefix + "alias = c2[0]\ndel alias[0]\n" + aliases.RESULT)
    assert aliases.derived(control) == (anchor,)
    assert aliases.derived(changed) == ()


def _remembered_order() -> list[str]:
    """Observe real insertion, without changing a fact or forcing allocator reuse."""
    inserted: list[str] = []
    original = guard.DerivationEnvironment.setdefault

    def record(env: guard.DerivationEnvironment, key: str, default: guard.Fact) -> guard.Fact:
        if key not in env:
            inserted.append(key)
        return original(env, key, default)

    text, _ = aliases.source(
        "module",
        aliases.EXTENT
        + "alias = p if flag else (zeta if other_flag else beta)\ndel alias[0]\n"
        + aliases.RESULT,
    )
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(guard.DerivationEnvironment, "setdefault", record)
        assert aliases.derived(text) == ()
    return inserted


@pytest.mark.parametrize("seed", [str(seed) for seed in range(8)])
def test_remembered_free_origins_have_a_stable_insertion_order(seed: str) -> None:
    """An intermediate ordering contract, not a claim of changed BY or work cost."""
    assert guard.__file__ is not None
    root = Path(__file__).resolve().parents[1]
    guard_path = root / "tests" / "test_multiple_testing_registry.py"
    assert Path(guard.__file__).resolve() == guard_path
    environment = os.environ.copy()
    environment.update(
        PYTHONHASHSEED=seed,
        PYTHONDONTWRITEBYTECODE="1",
        PYTHONPATH=os.pathsep.join(str(root / part) for part in ("src", "tests")),
    )
    completed = subprocess.run(
        [sys.executable, "-B", str(Path(__file__).resolve()), "--remember-order"],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=True,
        timeout=30,
    )
    report = json.loads(completed.stdout)
    assert report["seed"] == seed
    assert report["guard_path"] == str(guard_path)
    assert report["guard_sha256"] == hashlib.sha256(guard_path.read_bytes()).hexdigest()
    assert report["inserted"] == ["beta", "p", "zeta"]


if __name__ == "__main__":
    if sys.argv[1:] != ["--remember-order"]:
        raise SystemExit("Use --remember-order to observe the insertion-order contract")
    assert guard.__file__ is not None
    path = Path(guard.__file__).resolve()
    print(
        json.dumps(
            {
                "seed": os.environ.get("PYTHONHASHSEED"),
                "guard_path": str(path),
                "guard_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "inserted": _remembered_order(),
            },
            sort_keys=True,
        )
    )
