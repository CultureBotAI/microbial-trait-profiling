"""PR1438 delta regressions for #1470, #1472 and #1473, 2026-11-05."""

from copy import deepcopy
from typing import Any

import pytest

from test_multiple_testing_registry import registry_errors
from test_multiple_testing_registry_review import SOURCE, fixture_registry


def single_method_registry() -> dict[str, Any]:
    registry = fixture_registry()
    registry["families"][0]["comparison_procedure"] = []
    registry["families"][0]["anchors"][0]["procedures"] = ["BY"]
    return registry


@pytest.mark.parametrize("role", ["primary", "comparison"])
def test_own_family_roles_cannot_borrow_another_family_anchor(role: str) -> None:
    registry = single_method_registry()
    # A different family covers BH at the shared site. The first still may not
    # claim BH without declaring it at one of its own anchors.
    source = "def screen(p):\n    return step_up(p, method='by'), step_up(p, method='bh')\n"
    other = deepcopy(registry["families"][0])
    other.update(family_id="other", primary_procedure="BH")
    other["anchors"][0]["procedures"] = ["BH"]
    registry["families"].append(other)
    registry["literal_procedure_counts"][0]["counts"]["BH"] = 1
    assert not registry_errors({"src/x.py": source}, registry)
    first = registry["families"][0]
    if role == "primary":
        first.update(primary_procedure="BH", comparison_procedure=["BY"])
    else:
        first["comparison_procedure"] = ["BH"]
    assert any("own anchors" in e for e in registry_errors({"src/x.py": source}, registry))


def test_anchor_edit_cannot_invent_a_literal_method() -> None:
    registry = single_method_registry()
    assert not registry_errors(SOURCE, registry)
    family = registry["families"][0]
    family.update(primary_procedure="BH", comparison_procedure=["BY"])
    family["anchors"][0]["procedures"].append("BH")
    assert any("all-literal" in e for e in registry_errors(SOURCE, registry))


def test_derived_by_is_narrow_and_does_not_disable_literal_counts() -> None:
    source = {"src/x.py": "def screen(p):\n    return benjamini_hochberg(p) * harmonic(p)\n"}
    registry = single_method_registry()
    family = registry["families"][0]
    family["comparison_procedure"] = ["BH"]
    anchor = family["anchors"][0]
    anchor["procedures"] = ["BH", "BY"]
    registry["literal_procedure_counts"][0]["counts"] = {"BH": 1}
    assert registry_errors(source, registry)
    anchor["derived_procedures"] = {"BY": "BH multiplied by the family harmonic sum."}
    assert not registry_errors(source, registry)
    registry["literal_procedure_counts"][0]["counts"]["BH"] = 2
    assert any("literal procedure count mismatch" in e for e in registry_errors(source, registry))


@pytest.mark.parametrize(
    "derived", [[], {"Holm": "formula"}, {"BH": "formula"}, {"BY": ""}, {"BY": "formula"}]
)
def test_invalid_derived_declarations_fail(derived: Any) -> None:
    registry = single_method_registry()
    registry["families"][0]["anchors"][0]["derived_procedures"] = derived
    assert any("derived" in e for e in registry_errors(SOURCE, registry))


def test_dynamic_calls_are_not_treated_as_exhaustive_literals() -> None:
    registry = single_method_registry()
    registry["families"][0]["comparison_procedure"] = ["BH"]
    registry["families"][0]["anchors"][0]["procedures"] = ["BH", "BY"]
    source = "def screen(p, m):\n    return step_up(p, method='by'), step_up(p, method=m)\n"
    assert not registry_errors({"src/x.py": source}, registry)


@pytest.mark.parametrize(
    "source",
    [
        "def screen(p, m):\n    return step_up(p, method=m)\n",
        "def benjamini_hochberg(p):\n    return p\n",
    ],
)
def test_derived_exemptions_require_all_literal_calls(source: str) -> None:
    registry = single_method_registry()
    family = registry["families"][0]
    function = "benjamini_hochberg" if "def benjamini_hochberg" in source else "screen"
    family["anchors"][0].update(function=function, derived_procedures={"BY": "formula"})
    family["consumers"][0]["function"] = function
    registry["literal_procedure_counts"][0].update(function=function, counts={})
    assert any("all-literal" in e for e in registry_errors({"src/x.py": source}, registry))


def test_redundant_launcher_flag_fails() -> None:
    registry = single_method_registry()
    registry["families"][0]["launcher"] = True
    assert any("launcher" in e for e in registry_errors(SOURCE, registry))


@pytest.mark.parametrize(
    "variant", [" SHARED REASON ", "shared\t reason", "shared reason...", "Shared reason —"]
)
def test_duplicate_bh_reason_normalization(variant: str) -> None:
    registry = single_method_registry()
    first = registry["families"][0]
    first.update(primary_procedure="BH", deliberate_difference="Shared reason.")
    first["anchors"][0]["procedures"] = ["BH"]
    registry["literal_procedure_counts"][0]["counts"] = {"BH": 1}
    second = deepcopy(first)
    second.update(family_id="second", deliberate_difference=variant)
    registry["families"].append(second)
    source = {"src/x.py": SOURCE["src/x.py"].replace("'by'", "'bh'")}
    assert any("duplicate BH-primary reason" in e for e in registry_errors(source, registry))


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("outputs", [None]),
        ("outputs", ["results.tsv"]),
        ("outputs", [{}]),
        ("outputs", [{"file": "", "column": "q"}]),
        ("outputs", [{"file": "results.tsv", "column": []}]),
        ("policy_rule", ["R0"]),
        ("policy_rule", ["R10"]),
        ("policy_rule", ["r1"]),
        ("policy_rule", [None]),
        ("policy_rule", [{}]),
    ],
)
def test_output_items_and_policy_rules_have_a_schema(field: str, value: Any) -> None:
    registry = single_method_registry()
    assert not registry_errors(SOURCE, registry)
    registry["families"][0][field] = value
    assert registry_errors(SOURCE, registry)
