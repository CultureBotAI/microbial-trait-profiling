# Run tables for the NCASM plate tests

Copied by `scripts/adhoc/issue2264_capture_run_table_fixture.py` from one analysis
run (see `comparison/comparison_run.json` for its population and code revision).
They are real run output, kept so the table-driven plates can be drawn and checked
without the run directory. Do not edit them; capture a new fixture instead.

| File | Source SHA-256 | Copy |
| --- | --- | --- |
| `comparison/arm1_pooled_stats/trait_profile_paired_contrast_trait_stats_pooled.tsv` | `590fed8ba72317f2d2389171e01b5f8d7a2daac02a0986d7848e2b1126ccb83d` | byte for byte |
| `arms/noel/trait_direction_screen/trait_direction_screen.tsv` | `a83428df4ead141328914d2f4865e8774693e79aa98f928edae1aaf62e052caa` | byte for byte |
| `comparison/comparison_run.json` | `76cd13f6c0b4e936e7b04294137cc1dca33c0d4a24486522d1e04740f07f16fb` | reduced to 7 fields: id, contrast_population, direction_basis, status, core_run_id, git_commit, git_dirty |
