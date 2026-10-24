"""Compatibility shim forwarding bsdbng.deck_visualizations to mtp.deck_visualizations."""

import sys
from typing import Any

import mtp.deck_visualizations as _target
from mtp.deck_visualizations import (
    EXPLORATORY_LEADS,
    TABLE_1_TRAITS,
    TRAIT_SHORT_NAMES,
    clean_trait_label,
    generate_all_deck_figures,
    generate_figure1_condition_anatomy_ora_heatmap,
    generate_figure2_calibrated_ora_infographic,
    generate_figure3_lead_findings_ora_evidence,
    generate_figure4_annotation_density_pc1_landscape,
    generate_figure5_pd_crc_taxon_attribution_bubble,
    load_profile_data,
    main,
)

__all__ = [
    "EXPLORATORY_LEADS",
    "TABLE_1_TRAITS",
    "TRAIT_SHORT_NAMES",
    "clean_trait_label",
    "generate_all_deck_figures",
    "generate_figure1_condition_anatomy_ora_heatmap",
    "generate_figure2_calibrated_ora_infographic",
    "generate_figure3_lead_findings_ora_evidence",
    "generate_figure4_annotation_density_pc1_landscape",
    "generate_figure5_pd_crc_taxon_attribution_bubble",
    "load_profile_data",
    "main",
]

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.deck_visualizations", run_name="__main__")
