#!/usr/bin/env python3
"""Run ORA background stratification and annotation-density sensitivity evaluation.

Addresses Issue #118:
- Partitions study trait profiles into anatomical site strata (pooled, feces, feces_excluded).
- Evaluates hypergeometric ORA enrichment within each stratum background.
- Exports stratified trait and pair ranking tables under stable filenames:
    data/derived/trait_enrichment_traits_overall_ranked_<stratum>.tsv
    data/derived/trait_enrichment_trait_study_pairs_ranked_<stratum>.tsv
  emitting n_profiles_tested, trait_total, and body_site_group.
- Generates sensitivity comparison reports (TSV and TXT):
    data/derived/stratification_sensitivity_summary.tsv
    data/derived/stratification_sensitivity_summary.txt
  evaluating annotation-density deciles, effect heterogeneity, and support states.
"""

from __future__ import annotations

import sys

from mtp.stratification import main

if __name__ == "__main__":
    sys.exit(main())
