"""Keep historical table citations and study-family attribution scoped (#2061/#2062)."""

from pathlib import Path

DOC = Path(__file__).resolve().parents[1] / "methods/three_arm_statistics.md"
TABLE_LINK = (
    "[`table1_trait_condition_body_site_maps.md`](table1_trait_condition_body_site_maps.md)"
)


def directional_cell(label: str) -> str:
    rows = [line for line in DOC.read_text().splitlines() if line.startswith(f"| {label} |")]
    assert len(rows) == 1
    return rows[0].split("|")[-2].strip()


def test_typical_n_table_citation() -> None:
    cell = directional_cell("Typical n per test")
    assert f"across the 317 `off`-scope candidates ({TABLE_LINK})" in cell
    assert "trait_direction_screen.md" not in cell
    assert "table1_direction_concordance.md" not in cell
    assert "not identified ([Vintage](#vintage)). 606 to 898 for the six of Table 1" in cell
    assert "932 total" in cell
    assert "519a095b" in cell
    assert "ties excluded" in cell


def test_limitations_table_referent() -> None:
    cell = directional_cell("Main limitations")
    assert f"across the 317 `off`-scope candidates in that table ({TABLE_LINK};" in cell
    assert "that historical screen" not in cell
    assert "trait_direction_screen.md" not in cell
    assert "table1_direction_concordance.md" not in cell
    assert "[Vintage](#vintage)" in cell


def test_study_family_helper_and_correction() -> None:
    text = DOC.read_text()
    start = text.index("Jed's transformed experiment-level mean")
    paragraph = " ".join(text[start:].split("\n\n", 1)[0].split())
    assert (
        "`mt-core-native-study-t` via `study_contrasts` in "
        "[`src/mtp/native_composition.py`](../src/mtp/native_composition.py), "
        "corrected in [`scripts/create_trait_profile_infographic.py`]"
        "(../scripts/create_trait_profile_infographic.py)"
    ) in paragraph
    assert "[#133]" in paragraph
