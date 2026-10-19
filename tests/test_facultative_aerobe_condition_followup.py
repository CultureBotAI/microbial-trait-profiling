"""The autism follow-up sentence is built from the computed trait z shifts (#647)."""

import importlib.util
import math
import sys
from pathlib import Path

import pandas as pd
import pytest

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


@pytest.mark.parametrize("undefined", ["Autism", "Autism spectrum disorder"])
def test_undefined_shift_is_not_called_zero_or_divergent(undefined: str) -> None:
    summary = profile_summary(0.744, -1.682)
    no_decreased = summary["condition"].eq(undefined)
    summary.loc[no_decreased, "increased_minus_decreased_guild_z"] = math.nan
    summary.loc[no_decreased, "decreased_profile_count"] = 0

    sentence = followup.autism_label_contrast_sentence(summary)

    assert (
        f"`{undefined}` has no computable increased-minus-decreased trait z shift, which needs "
        "both increased and decreased profiles ("
    ) in sentence
    assert " and 0 decreased profiles from " in sentence
    assert "diverge" not in sentence
    assert "zero" not in sentence
    assert "nan" not in sentence


def test_readme_states_the_computed_shifts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    summary = profile_summary(0.744, -1.682).assign(
        condition_group=["Autism", "Autism", "Autoimmune disease"],
        profile_count=[14, 7, 8],
        body_site_groups="Feces",
    )
    grouped = pd.DataFrame(
        columns=[
            "condition_group",
            "direction",
            "taxon_name",
            "taxonomic_rank",
            "profile_count",
            "study_count",
            "guild_trait_names",
        ]
    )
    monkeypatch.setattr(followup, "OUTPUT_DIR", tmp_path)

    followup.write_readme(
        summary,
        grouped,
        pd.DataFrame({"taxon_id": []}),
        pd.DataFrame(),
        pd.DataFrame({"pattern_name": []}),
        pd.DataFrame(),
        pd.DataFrame(),
    )

    readme = (tmp_path / "autism_autoimmune_followup.md").read_text(encoding="utf-8")
    assert followup.autism_label_contrast_sentence(summary) in readme
    assert "slightly" not in readme
