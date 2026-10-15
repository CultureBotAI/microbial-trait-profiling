import importlib.util
import sys
from pathlib import Path

import pandas as pd

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_pd_ecological_reporter_combined.py"
)
SPEC = importlib.util.spec_from_file_location("pd_ecological_reporter_combined", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
reporter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = reporter
SPEC.loader.exec_module(reporter)


def _condition_rows(
    condition: str,
    *,
    significant_trait: str | None = None,
    paired_study_count: int = 3,
    significant_composite: bool = False,
) -> list[dict[str, object]]:
    rows = [
        {
            "condition": condition,
            "trait": trait,
            "paired_study_count": paired_study_count,
            "ci_excludes_zero": trait == significant_trait,
        }
        for trait in reporter.SIX_TRAITS
    ]
    rows.append(
        {
            "condition": condition,
            "trait": reporter.REPORTER_COMPOSITE,
            "paired_study_count": paired_study_count,
            "ci_excludes_zero": significant_composite,
        }
    )
    return rows


def test_qualifying_conditions_requires_one_of_six_individual_traits() -> None:
    summary = pd.DataFrame(
        [
            *_condition_rows(
                reporter.PD_CONDITION,
                significant_trait="arsenate reduction",
            ),
            *_condition_rows("Composite only", significant_composite=True),
            *_condition_rows(
                "Too little support",
                significant_trait="coccus shaped",
                paired_study_count=2,
            ),
            *_condition_rows("No change"),
        ]
    )

    selected = reporter.qualifying_conditions(summary)

    assert selected["condition"].tolist() == [reporter.PD_CONDITION]
    assert selected.iloc[0]["qualifying_significant_trait_count"] == 1
    assert selected.iloc[0]["qualifying_significant_traits"] == "arsenate reduction"


def test_panel_b_composite_uses_the_three_supported_positive_pd_traits() -> None:
    assert reporter.INCREASED_SIDE_TRAITS == [
        "moderately halophilic",
        "facultatively aerobic",
        "arsenate reduction",
    ]
