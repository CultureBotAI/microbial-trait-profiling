"""Deterministic float-reduction order contracts, independent of hash seeds (#2323)."""

from __future__ import annotations

import builtins
import unittest
from collections.abc import Iterable
from unittest.mock import patch

from mtp import taxon_influence as influence


class InfluenceOrderTests(unittest.TestCase):
    def test_hhi_resolved_consumes_canonical_summands(self) -> None:
        resolved = {"taxon:c": 4, "taxon:a": 1, "taxon:b": 2}
        adversarial = list(resolved)
        canonical = sorted(resolved)
        assert adversarial != canonical
        total = sum(resolved.values())
        expected = [(resolved[key] / total) ** 2 for key in canonical]
        assert expected != [(resolved[key] / total) ** 2 for key in adversarial]
        observed: list[list[float]] = []

        def recording_sum(values: Iterable[float]) -> float:
            terms = list(values)
            if any(isinstance(term, float) for term in terms):
                observed.append(terms)
            return builtins.sum(terms)

        with patch.object(influence, "sum", recording_sum, create=True):
            actual = influence._hhi_resolved(resolved, total)
        assert observed == [expected]
        assert actual == builtins.sum(expected)

    def test_weighted_overlap_consumes_canonical_summands(self) -> None:
        left = {"taxon:c": 6, "taxon:a": 1, "taxon:b": 3}
        right = {"taxon:b": 2, "taxon:c": 3, "taxon:a": 5}
        shared = ("taxon:c", "taxon:a", "taxon:b")
        canonical = sorted(left.keys() & right.keys())
        assert set(shared) == set(canonical)
        assert list(shared) != canonical
        left_total, right_total = sum(left.values()), sum(right.values())
        expected = [min(left[key] / left_total, right[key] / right_total) for key in canonical]
        assert expected != [min(left[key] / left_total, right[key] / right_total) for key in shared]
        observed: list[list[float]] = []

        def recording_sum(values: Iterable[float]) -> float:
            terms = list(values)
            if any(isinstance(term, float) for term in terms):
                observed.append(terms)
            return builtins.sum(terms)

        with patch.object(influence, "sum", recording_sum, create=True):
            actual = influence._weighted_overlap(left, right, shared)
        assert observed == [expected]
        assert actual == builtins.sum(expected)


if __name__ == "__main__":
    unittest.main()
