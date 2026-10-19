"""The autism follow-up sentence is built from the computed trait z shifts (#647)."""

import importlib.util
import sys
from pathlib import Path

import pandas as pd

MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "investigate_facultative_aerobe_condition_followup.py"
)
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("aerobe_condition_followup", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
followup = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(followup)


def profile_summary(autism_shift: float, asd_shift: float) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "condition": ["Autism", "Autism spectrum disorder", "Autoimmune disease"],
            "increased_minus_decreased_guild_z": [autism_shift, asd_shift, 0.421],
            "increased_profile_count": [8, 5, 4],
            "decreased_profile_count": [6, 2, 4],
            "study_count": [4, 3, 2],
        }
    )


def test_diverging_shifts_are_stated_with_sign_size_and_n() -> None:
    sentence = followup.autism_label_contrast_sentence(profile_summary(0.744, -1.682))

    assert "directional skews diverge" in sentence
    assert "`Autism` has a positive increased-minus-decreased trait z shift (+0.744" in sentence
    assert (
        "`Autism spectrum disorder` has a negative increased-minus-decreased trait z shift "
        "(-1.682; 5 increased and 2 decreased profiles from 3 studies)"
    ) in sentence
    assert "slightly" not in sentence


def test_same_sign_shifts_do_not_claim_divergence() -> None:
    sentence = followup.autism_label_contrast_sentence(profile_summary(0.744, 0.2))

    assert "diverge" not in sentence
    assert "`Autism spectrum disorder` has a positive" in sentence


def test_missing_label_keeps_only_the_reason_for_separate_reporting() -> None:
    summary = profile_summary(0.744, -1.682)
    summary = summary[summary["condition"] != "Autism spectrum disorder"]

    sentence = followup.autism_label_contrast_sentence(summary)

    assert sentence.endswith("reported separately rather than pooled (#128).")
