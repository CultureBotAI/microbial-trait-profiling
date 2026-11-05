"""Loop-exit facts must reach subsequent iterations and statements (#1591)."""

import pytest

from test_multiple_testing_registry import registry_errors
from test_multiple_testing_registry_derived import DERIVATION_ERROR, by_registry, source
from test_multiple_testing_registry_flow import BH, HEAD, RESULT


@pytest.mark.parametrize("loop", ["for x in other:", "while other:"])
@pytest.mark.parametrize("exit_kind", ["continue", "break"])
@pytest.mark.parametrize("changed", [True, False], ids=["stale", "preserved"])
def test_exit_path_facts_reach_loop_head_and_following_code(
    loop: str, exit_kind: str, changed: bool
) -> None:
    body = HEAD + BH + "keep = q\n" + loop + "\n"
    if exit_kind == "continue":
        body += f"    by = {RESULT}\n"
    body += (
        "    if other:\n"
        f"        q = {'p' if changed else 'keep'}\n"
        f"        {exit_kind}\n"
        "    q = keep\n"
    )
    body += "return by\n" if exit_kind == "continue" else f"return {RESULT}\n"
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    if changed:
        assert any(DERIVATION_ERROR in error for error in errors), errors
    else:
        assert not errors


@pytest.mark.parametrize("exit_kind", ["continue", "break"])
def test_inner_loop_exits_do_not_bypass_outer_restoration(exit_kind: str) -> None:
    body = (
        HEAD + BH + "keep = q\nfor x in other:\n"
        "    for y in other:\n"
        "        q = p\n"
        f"        {exit_kind}\n"
        "    q = keep\n"
        f"    by = {RESULT}\n"
        "return by\n"
    )
    assert not registry_errors({"src/x.py": source(body)}, by_registry())


@pytest.mark.parametrize("loop", ["for x in other:", "while other:"])
@pytest.mark.parametrize("changed", [True, False], ids=["stale", "preserved"])
def test_break_skips_the_loop_else_clause(loop: str, changed: bool) -> None:
    """Python skips a loop's else on break, so the else must not overwrite break facts."""
    body = (
        HEAD + BH + "keep = q\n" + loop + "\n"
        "    if other:\n"
        f"        q = {'p' if changed else 'keep'}\n"
        "        break\n"
        "else:\n"
        "    q = keep\n"
        f"return {RESULT}\n"
    )
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    if changed:
        assert any(DERIVATION_ERROR in error for error in errors), errors
    else:
        assert not errors


@pytest.mark.parametrize("loop", ["for x in other:", "while other:"])
def test_continue_state_reaches_code_after_the_loop(loop: str) -> None:
    """The last iteration may end on continue, so its state reaches the code after the loop."""
    body = (
        HEAD + BH + "keep = q\n" + loop + "\n"
        "    if other:\n"
        "        q = p\n"
        "        continue\n"
        "    q = keep\n"
        f"return {RESULT}\n"
    )
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


def test_inner_loop_break_belongs_to_the_inner_loop() -> None:
    """A break in an inner loop exits that loop, so its state reaches the outer loop's body."""
    body = (
        HEAD + BH + "keep = q\nfor x in other:\n"
        "    for y in other:\n"
        "        if y:\n"
        "            q = p\n"
        "            break\n"
        "        q = keep\n"
        f"    by = {RESULT}\n"
        "return by\n"
    )
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors
