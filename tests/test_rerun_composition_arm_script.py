"""The composition wrapper delegates to the verified standard suite (#658).

The historical PD contributor suite and its extra pinned-KG requirement (#695)
remain outside the standard analysis/calibration scope.
"""

from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "rerun_composition_arm_jed.sh"


def test_wrapper_requires_explicit_core_and_fresh_output() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "if [[ $# -lt 2 ]]; then" in text
    assert "CORE_RUN_JSON FRESH_OUTPUT_DIR" in text
    assert 'exec just rerun-jed "$@"' in text
    assert "\nexport MPLBACKEND=Agg\n" in text
