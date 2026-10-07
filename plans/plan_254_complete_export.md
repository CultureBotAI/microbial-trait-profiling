# Plan: #254 Complete Unfiltered Export

The goal is to emit a complete per-trait table with no selection applied, making the three arms (composition, directional, ORA) joinable.

## Proposed Schema (All-Results TSV)

The analyzer will write a new file `trait_enrichment_all_results.tsv` (and `trait_depletion_all_results.tsv`) with these columns:

- `study_id`, `experiment_id`, `inc_dec`, `trait` (raw IDs)
- `p_value`, `fdr_q`
- `fold` (effect size)
- `count` (observed hits k)
- `expected` (expected hits)
- `testable_profiles` (n)
- `background_m` (M)
- `background_k` (K)
- `is_significant` (boolean, based on current thresholds)
- `correction_family` (e.g., "global" or "profile")
- `eligibility_status` (e.g., "testable", "ineligible")

## Implementation Steps

1. **Update `trait_enrichment.py` and `trait_depletion.py`**:
   - Add `--all-results-output` argument.
   - In the analyzer, instead of filtering `significant_hits` before writing, write the full `all_hits` list to the new file.
   - Add the new columns (`background_m`, `background_k`, `is_significant`, `correction_family`, `eligibility_status`) to the `all_hits` dictionaries.
   - Use the shared `enrichment_results.py` to write the metadata header (Result-Status, etc.).

2. **Wire into `Justfile`**:
   - Update `analyze` and `analyze-depletion` recipes to generate these new files.

3. **Tests**:
   - Add tests to `tests/test_enrichment_pipeline.py` verifying:
     - Row count in the all-results file equals `Total-Tests-Performed`.
     - The significant-hits file is a subset of the all-results file.
     - Hand-computed M, N, K, k values match the output.
     - Zero-discovery runs still write the full table (with `is_significant=False`).

4. **Review**:
   - Confirm filenames and column list with you.
   - Open PR with `Refs #254`.

Does this schema and plan meet your requirements for #254?
