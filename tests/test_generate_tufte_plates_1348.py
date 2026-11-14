"""Pin Plate 2's oxygen annotation label to METPO:1000606 (#1348)."""

import ast
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1] / "scripts/generate_tufte_plates.py"


def test_plate2_obligately_aerobic_label() -> None:
    tree = ast.parse(SOURCE.read_text())
    assignments = [
        node
        for node in tree.body
        if isinstance(node, ast.AnnAssign)
        and isinstance(node.target, ast.Name)
        and node.target.id == "PLATE2_TRAIT_ROWS"
    ]
    assert len(assignments) == 1
    value = assignments[0].value
    assert value is not None
    rows = ast.literal_eval(value)
    matches = [row for row in rows if row[0] == "METPO:1000606"]
    assert matches == [("METPO:1000606", "Obligately aerobic", "O₂ / Aerotolerance")]
