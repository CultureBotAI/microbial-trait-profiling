"""Pin the six semantic survivors of the independent delta review (#1592)."""

from __future__ import annotations

import pytest

from test_multiple_testing_registry import registry_errors
from test_multiple_testing_registry_derived import DERIVATION_ERROR, by_registry, source

HEAD = "m = len(p)\nh = np.sum(1.0 / np.arange(1, m + 1))\nq = benjamini_hochberg(p)\n"


@pytest.mark.parametrize(
    "body",
    [
        pytest.param(
            HEAD + "for x in other:\n"
            "    if x:\n"
            "        return np.clip(q * h, 0, 1)\n"
            "    q = p\n"
            "return q\n",
            id="return-credit-ignores-flag",
        ),
        pytest.param(
            HEAD + "if True:\n    return q\nby = np.clip(q * h, 0, 1)\n",
            id="constant-if-terminator-ignored",
        ),
        pytest.param(
            HEAD + "with manager():\n    if True:\n        return q\nby = np.clip(q * h, 0, 1)\n",
            id="terminates-ignores-constant-if",
        ),
        pytest.param(
            HEAD + "try:\n    q[:] = p\nexcept ValueError:\n    pass\n"
            "return np.clip(q * h, 0, 1)\n",
            id="try-writes-not-invalidated",
        ),
        pytest.param(
            HEAD + "import q.sub\nreturn np.clip(q * h, 0, 1)\n",
            id="dotted-import-binds-full-name",
        ),
        pytest.param(
            HEAD + "out = {}\nout['by'] = np.clip(q * h, 0, 1)\nreturn out\n",
            id="assign-credit-any-target",
        ),
    ],
)
def test_delta_mutant_witness_is_not_an_established_derivation(body: str) -> None:
    """Keep dictionary stores conservatively unsupported, along with stale/dead facts."""
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors
