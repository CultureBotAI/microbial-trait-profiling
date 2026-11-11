"""Only candidate functions enter the BY derivation interpreter (#1593)."""

from __future__ import annotations

from collections import Counter
from textwrap import dedent

import pytest

import test_multiple_testing_registry as guard


def test_derivation_gate_uses_resolved_calls_owned_by_each_function(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = dedent(
        """\
        import numpy as arr
        import other
        from numpy import clip as clip_alias
        from mtp.enrichment_stats import benjamini_hochberg as bh

        def both(p):
            return arr.clip(bh(p), 0, 1)

        def bh_only(p):
            return bh(p)

        def clip_only(p):
            return arr.clip(p, 0, 1)

        def unrelated_clip(p):
            return other.clip(bh(p), 0, 1)

        def parent_bh(p):
            q = bh(p)
            def child_clip():
                return arr.clip(q, 0, 1)
            return q

        def parent_clip(p):
            q = arr.clip(p, 0, 1)
            def child_bh():
                return bh(q)
            return q

        def enclosing(p):
            def own_both():
                return clip_alias(bh(p), 0, 1)
            return own_both()

        def aliases(p):
            adjust = bh
            bound = clip_alias
            return bound(adjust(p), 0, 1)

        def empty(p):
            return p

        arr.clip(bh(p), 0, 1)
        """
    )
    entered: list[str] = []

    def record_candidate(
        unit: guard.Unit, scope: str, *, file_budget: guard.FileDerivationBudget
    ) -> bool:
        assert scope in unit.functions
        assert file_budget.limit == guard.DEFAULT_LIMITS.file_derivation_steps
        entered.append(scope)
        return False

    monkeypatch.setattr(guard, "visible_by_derivation", record_candidate)
    summary = guard.summarize_source("src/x.py", source, guard.SummaryBudget(guard.DEFAULT_LIMITS))

    assert entered == ["both", "enclosing.own_both", "aliases"]
    assert summary.derived == ()
    bh = "mtp.enrichment_stats.benjamini_hochberg"
    assert Counter((scope, target) for scope, target, _, _ in summary.calls) == Counter(
        [
            ("both", "numpy.clip"),
            ("both", bh),
            ("bh_only", bh),
            ("clip_only", "numpy.clip"),
            ("unrelated_clip", "other.clip"),
            ("unrelated_clip", bh),
            ("parent_bh", bh),
            ("parent_bh.child_clip", "numpy.clip"),
            ("parent_clip", "numpy.clip"),
            ("parent_clip.child_bh", bh),
            ("enclosing.own_both", "numpy.clip"),
            ("enclosing.own_both", bh),
            ("enclosing", "x.enclosing.own_both"),
            ("aliases", "numpy.clip"),
            ("aliases", bh),
            ("<module>", "numpy.clip"),
            ("<module>", bh),
        ]
    )
    for _, target, procedures, _ in summary.calls:
        assert procedures == (("BH",) if target == bh else ())


@pytest.mark.parametrize(
    ("imports", "bh_call", "clip_call", "sum_call", "range_call"),
    [
        (
            "import numpy as arr\nfrom mtp.enrichment_stats import benjamini_hochberg as adjust\n",
            "adjust",
            "arr.clip",
            "arr.sum",
            "arr.arange",
        ),
        (
            "from numpy import clip as bound, sum as total, arange as sequence\n"
            "from mtp.enrichment_stats import benjamini_hochberg as adjust\n",
            "adjust",
            "bound",
            "total",
            "sequence",
        ),
    ],
    ids=["module-alias", "from-import-aliases"],
)
def test_derivation_gate_preserves_supported_alias_witnesses(
    imports: str, bh_call: str, clip_call: str, sum_call: str, range_call: str
) -> None:
    source = (
        imports
        + "def screen(p):\n"
        + "    m = len(p)\n"
        + f"    h = {sum_call}(1.0 / {range_call}(1, m + 1))\n"
        + f"    q = {bh_call}(p)\n"
        + f"    return {clip_call}(q * h, 0, 1)\n"
    )
    summary = guard.summarize_source("src/x.py", source, guard.SummaryBudget(guard.DEFAULT_LIMITS))
    assert summary.derived == ("screen",)
