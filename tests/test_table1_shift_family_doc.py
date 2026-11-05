"""Bind the exported shift-family documentation to synthetic writer output (#1531)."""

from __future__ import annotations

import ast
import csv
import re
from pathlib import Path

import pytest

from mtp import table1_trait_maps
from mtp.trait_guilds import TABLE1_TRAITS, Profile

REPO_ROOT = Path(__file__).resolve().parents[1]
METHOD = REPO_ROOT / "methods/table1_trait_condition_body_site_maps.md"
FAMILY_COLUMNS = {"family_id", "family_size", "correction_method"}
STATISTIC_COLUMNS = {"binomial_p", "bh_q", "status"}
PRODUCTION_CALLERS = (
    ("table1_trait_maps", "main"),
    ("calibration", "fraction_sign_test"),
    ("calibration", "screened_fraction_sign_test"),
    ("contrast_orientation", "comparison_rows"),
    ("integrated_evidence", "fraction_sign_records"),
)


@pytest.fixture
def section() -> str:
    content = METHOD.read_text(encoding="utf-8")
    exported = content.split("### Exported correction-family columns\n", 1)[1]
    return " ".join(exported.split("\n## ", 1)[0].split())


def _written_shifts(
    output: Path, traits: tuple[str, ...]
) -> tuple[list[str], dict[str, dict[str, str]]]:
    """Four independent pairs score every trait except the final, all-tied trait."""
    profiles = [
        Profile(
            study_id=f"study:{study}",
            experiment_id=f"experiment:{study}",
            inc_dec=direction,
            condition="synthetic condition",
            body_site="Feces",
            guild_counts={trait: (3 if trait == traits[-1] else count) for trait in traits},
            total=100,
        )
        for study in range(4)
        for direction, count in (("inc", 7), ("dec", 1))
    ]
    result = table1_trait_maps.compute_direction_shifts(profiles, traits)
    shift_path, _ = table1_trait_maps.write_shift_table(result, output, "synthetic")
    with shift_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return list(reader.fieldnames), {row["trait"]: row for row in reader}


@pytest.fixture
def exported(tmp_path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    return _written_shifts(tmp_path, table1_trait_maps.TRAIT_ORDER)


def test_table1_maps_doc_describes_exported_shift_family_fields(section: str) -> None:
    """The exported-family section must carry its own provenance, not nearby prose."""
    content = METHOD.read_text(encoding="utf-8")
    assert "Outstanding: the correction family is described here but not exported" not in content
    assert "nothing in the row says which correction family produced that q" not in content
    for token in (
        *sorted(FAMILY_COLUMNS | STATISTIC_COLUMNS),
        table1_trait_maps.SHIFT_FAMILY_ID,
        table1_trait_maps.SHIFT_CORRECTION,
        "src/mtp/table1_trait_maps.py",
    ):
        assert token in section, f"Missing expected family provenance token {token!r}"


def test_documented_column_names_exist_in_real_writer_header(
    section: str, exported: tuple[list[str], dict[str, dict[str, str]]]
) -> None:
    header, rows = exported
    identifiers = set(re.findall(r'`"?([a-z][a-z_]+)"?`', section))
    noncolumns = {
        "trait_order",
        "compute_direction_shifts",
        "write_shift_table",
        table1_trait_maps.SHIFT_FAMILY_ID,
    }
    documented_columns = identifiers - noncolumns
    assert documented_columns >= FAMILY_COLUMNS | STATISTIC_COLUMNS
    assert documented_columns <= set(header)
    assert set(rows) == set(table1_trait_maps.TRAIT_ORDER)
    for row in rows.values():
        assert row["family_id"] == table1_trait_maps.SHIFT_FAMILY_ID
        assert row["correction_method"] == table1_trait_maps.SHIFT_CORRECTION


def test_documented_family_size_matches_production_list_and_export(
    section: str, exported: tuple[list[str], dict[str, dict[str, str]]]
) -> None:
    documented_sizes = re.findall(r"`family_size = (\d+)`", section)
    assert documented_sizes, "The production family size must be stated explicitly"
    expected_size = len(table1_trait_maps.TRAIT_ORDER)
    assert expected_size == len(TABLE1_TRAITS)
    assert {int(size) for size in documented_sizes} == {expected_size}
    _, rows = exported
    assert {int(row["family_size"]) for row in rows.values()} == {expected_size}


def test_unscored_trait_has_blank_p_and_q_but_remains_in_the_family(
    section: str, exported: tuple[list[str], dict[str, dict[str, str]]]
) -> None:
    assert "every trait in `trait_order` is included in the family" in section.lower()
    assert "`family_size = len(result.shifts)`" in section
    assert "exported `binomial_p` and `bh_q` cells remain blank" in section
    _, rows = exported
    unscored = rows[table1_trait_maps.TRAIT_ORDER[-1]]
    assert unscored["studies_scored"] == "0"
    assert unscored["status"] == "all_tied"
    assert unscored["binomial_p"] == unscored["bh_q"] == ""
    assert int(unscored["family_size"]) == len(rows)
    assert all(row["binomial_p"] and row["bh_q"] for row in list(rows.values())[:-1])


def test_documented_p_padding_and_bh_agree_with_synthetic_q_values(
    section: str, exported: tuple[list[str], dict[str, dict[str, str]]]
) -> None:
    padding = re.findall(r"`p = ([0-9.]+)`", section)
    assert padding
    assert {float(value) for value in padding} == {1.0}
    assert "no scored study" in section
    assert "Benjamini-Hochberg" in section
    _, rows = exported
    scored = [row for row in rows.values() if row["binomial_p"]]
    # Equal non-null p-values followed by one padded 1: BH gives p * m / (m - 1).
    p_values = {float(row["binomial_p"]) for row in scored}
    assert p_values == {2 / 2**4}
    expected_q = p_values.pop() * len(rows) / len(scored)
    assert expected_q < 1
    assert all(float(row["bh_q"]) == pytest.approx(expected_q) for row in scored)


def test_constant_family_id_also_describes_a_nonstandard_trait_list(
    section: str, tmp_path: Path
) -> None:
    assert (
        "`family_id` remains the constant `table1_six_fractions` even when the supplied "
        "trait list is not six long"
    ) in section
    traits = table1_trait_maps.TRAIT_ORDER[:2]
    assert len(traits) != len(table1_trait_maps.TRAIT_ORDER)
    _, rows = _written_shifts(tmp_path, traits)
    assert len(rows) == len(traits)
    for row in rows.values():
        assert row["family_id"] == table1_trait_maps.SHIFT_FAMILY_ID
        assert int(row["family_size"]) == len(traits)


@pytest.mark.parametrize(
    "symbol",
    ["compute_direction_shifts", "write_shift_table", "SHIFT_FAMILY_ID", "SHIFT_CORRECTION"],
)
def test_section_uses_stable_symbol_citations(section: str, symbol: str) -> None:
    assert f"`{symbol}`" in section
    assert not re.search(r"table1_trait_maps\.py:\d+", section)
    assert not re.search(r"`:\d+", section)


def test_section_identifies_the_verified_export_change(section: str) -> None:
    assert "added in commit `c86a312` (PR #241)" in section


@pytest.mark.parametrize(("module", "function"), PRODUCTION_CALLERS)
def test_named_production_callers_use_the_production_trait_list(
    section: str, module: str, function: str
) -> None:
    assert f"`{module}.{function}`" in section
    assert (
        "all use the six Table 1 traits through `TRAIT_ORDER`, explicitly or by default" in section
    )
    source = REPO_ROOT / "src/mtp" / f"{module}.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    caller = next(
        node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == function
    )
    calls = [
        node
        for node in ast.walk(caller)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "compute_direction_shifts"
    ]
    assert len(calls) == 1
    call = calls[0]
    trait_argument = (
        call.args[1]
        if len(call.args) > 1
        else next((kw.value for kw in call.keywords if kw.arg == "trait_order"), None)
    )
    if trait_argument is not None:
        assert isinstance(trait_argument, ast.Name)
        assert trait_argument.id == "TRAIT_ORDER"
    else:
        assert table1_trait_maps.compute_direction_shifts.__defaults__ == (
            table1_trait_maps.TRAIT_ORDER,
        )
