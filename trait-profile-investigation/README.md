# Trait Profile Investigation — Taxa Coverage Outputs

This directory contains the core output figures, data tables, and diagnostics for the **Approach A Taxa Coverage Analysis** between BugSigDB and KG-Microbe METPO traits.

The complete sub-investigation with detailed methodology, diagnostics, and documentation is maintained at:
- [`trait-visualization/trait-profile-investigation/taxa-coverage/`](../trait-visualization/trait-profile-investigation/taxa-coverage/)

---

## Output Manifest

### Visualizations:
- `taxa_coverage_histogram_rank_stratified.png` / `.svg`: Rank-stratified stacked histogram (Panel A: Absolute counts; Panel B: 100% normalized proportions).
- `taxa_coverage_occurrences_vs_unique.png` / `.svg`: Literature signature occurrences vs. unique taxa coverage disparity.
- `taxa_trait_count_distribution.png` / `.svg`: Trait degree distribution per taxon highlighting the zero-trait spike and rank rollups.

### Data Tables (TSVs):
- `taxa_coverage_by_rank.tsv`: Quantitative summary table by rank with counts, percentages, and occurrence weights.
- `taxa_kgmicrobe_coverage_details.tsv`: Complete audit table for all 5,829 unique taxa.
- `top_zero_trait_bugsigdb_taxa.tsv`: Top 100 zero-trait BugSigDB taxa with diagnosis category and suggested remediation.

---

## Key Metrics Summary

- **Total Unique BugSigDB Taxa**: 5,829 (38.69% with 0 traits, 61.31% trait-bearing).
- **Total Signature Literature Mentions**: 59,486 (**90.55% trait-bearing**, only 9.45% with 0 traits).
- **Rank Drop-off**: Species has a 45.87% zero-trait rate; higher ranks (Genus, Family, Order, Class, Phylum) maintain >72% to >79% trait-bearing rates for taxa, and >96% for literature occurrences.

---

## Regeneration Command

```bash
uv run python scripts/create_taxa_coverage_histogram.py
```
