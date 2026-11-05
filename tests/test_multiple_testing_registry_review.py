"""PR #1438 regression cases for issues #1441 and #1442 (2026-11-05)."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

import pytest

from test_multiple_testing_registry import registry_errors

SOURCE = {"src/x.py": "def screen(p):\n    return step_up(p, method='by')\n"}


def fixture_registry(counts: dict[str, int] | None = None) -> dict[str, Any]:
    return {
        "family_ids_are_new": True,
        "families": [
            {
                "family_id": "synthetic-family",
                "area": "synthetic study screen",
                "hypotheses": "One null for each supplied trait.",
                "family_definition": "Every supplied trait, including untestable members.",
                "size_source": {
                    "kind": "executed_inventory",
                    "selection_rule": "Include every input hypothesis.",
                    "count_rule": "Count the input hypotheses, including untestable slots.",
                },
                "primary_procedure": "BY",
                "comparison_procedure": ["BH"],
                "threshold": "q <= 0.05",
                "outputs": [{"file": "results.tsv", "column": "q"}],
                "consumers": [{"path": "src/x.py", "function": "screen"}],
                "policy_rule": ["R1"],
                "deliberate_difference": "BY primary for the synthetic all-trait screen.",
                "open_issue": [],
                "anchors": [{"path": "src/x.py", "function": "screen", "procedures": ["BH", "BY"]}],
            }
        ],
        "literal_procedure_counts": [
            {
                "path": "src/x.py",
                "function": "screen",
                "counts": {"BY": 1} if counts is None else counts,
            }
        ],
    }


@pytest.mark.parametrize("boolean_flag", [False, True], ids=["M08-method", "M13-by-flag"])
def test_one_literal_flip_cannot_hide_in_mixed_function(boolean_flag: bool) -> None:
    if boolean_flag:
        source = (
            "from mtp.taxon_influence import adjust\n"
            "def screen(p):\n"
            "    a = adjust(p, by=False)\n"
            "    b = adjust(p, by=True)\n"
            "    return a, b\n"
        )
        changed = source.replace("by=True", "by=False")
        expected = {"BH": 1, "BY": 1}
    else:
        source = (
            "def screen(p):\n"
            "    a = step_up(p, method='bh')\n"
            "    b = step_up(p, method='by')\n"
            "    c = step_up(p, method='by')\n"
            "    return a, b, c\n"
        )
        changed = source.replace("method='by'", "method='bh'", 1)
        expected = {"BH": 1, "BY": 2}
    registry = fixture_registry(expected)
    assert not registry_errors({"src/x.py": source}, registry)
    errors = registry_errors({"src/x.py": changed}, registry)
    assert any("literal procedure count" in error for error in errors), errors


@pytest.mark.parametrize(
    ("field_name", "value"),
    [
        ("family_definition", ""),
        ("hypotheses", " "),
        ("threshold", None),
        ("area", ""),
        ("family_id", []),
        ("deliberate_difference", ""),
        ("primary_procedure", "Bonferroni"),
        ("primary_procedure", ["BY"]),
        ("comparison_procedure", "BH"),
        ("comparison_procedure", ["Holm"]),
        ("open_issue", [True]),
        ("open_issue", ["1442"]),
        ("open_issue", "not an issue"),
        ("anchors", []),
        ("anchors", [{"path": "src/x.py", "function": ""}]),
        ("consumers", []),
        ("consumers", [{"path": "src/x.py", "function": "gone"}]),
        ("size_source", ""),
        ("size_source", {"kind": "constant"}),
        ("size_source", {"kind": "selected_rows", "selection_rule": "finite p", "count_rule": "n"}),
        ("policy_rule", []),
        ("outputs", []),
    ],
)
def test_schema_rejects_malformed_entry(field_name: str, value: Any) -> None:
    registry = fixture_registry()
    assert not registry_errors(SOURCE, registry)
    registry["families"][0][field_name] = value
    assert registry_errors(SOURCE, registry), (field_name, value)


def test_empty_comparison_and_issue_lists_are_valid() -> None:
    registry = fixture_registry()
    registry["families"][0].update(comparison_procedure=[], open_issue=[])
    registry["families"][0]["anchors"][0]["procedures"] = ["BY"]
    assert not registry_errors(SOURCE, registry)


def test_family_needs_discovered_anchor_or_explicit_launcher() -> None:
    sources = {**SOURCE, "src/y.py": "def launch(p):\n    pass\n"}
    registry = fixture_registry()
    launcher = deepcopy(registry["families"][0])
    launcher.update(
        family_id="replicate-launcher",
        anchors=[{"path": "src/y.py", "function": "launch", "procedures": ["BY"]}],
    )
    registry["families"].append(launcher)
    assert any("discovered correction anchor" in e for e in registry_errors(sources, registry))
    launcher["launcher"] = True
    assert not registry_errors(sources, registry)


@pytest.mark.parametrize("counts", [{"BY": True}, {"BY": -1}, {"BY": 1.5}, {"Holm": 1}, []])
def test_invalid_literal_count_declaration_fails(counts: Any) -> None:
    registry = fixture_registry()
    registry["literal_procedure_counts"][0]["counts"] = counts
    assert registry_errors(SOURCE, registry)


def test_missing_or_duplicate_count_binding_fails() -> None:
    registry = fixture_registry()
    registry["literal_procedure_counts"] *= 2
    assert registry_errors(SOURCE, registry)
    registry["literal_procedure_counts"] = []
    assert registry_errors(SOURCE, registry)


def test_shared_anchor_counts_are_not_multiplied_by_families() -> None:
    registry = fixture_registry()
    other = deepcopy(registry["families"][0])
    other["family_id"] = "another-family-at-same-site"
    registry["families"].append(other)
    assert not registry_errors(SOURCE, registry)


def test_duplicate_bh_reason_and_repeated_selection_definition_are_refused() -> None:
    registry = fixture_registry()
    first = registry["families"][0]
    first.update(primary_procedure="BH", comparison_procedure=["BY"])
    second = deepcopy(first)
    second["family_id"] = "second-family"
    registry["families"].append(second)
    assert any("BH-primary reason" in e for e in registry_errors(SOURCE, registry))
    registry["families"].pop()
    first["size_source"]["selection_rule"] = first["family_definition"]
    assert any("selection_rule" in e for e in registry_errors(SOURCE, registry))
