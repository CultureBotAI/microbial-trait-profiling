# Historical comparison declarations (#1606)

`analysis_specs_f0cb5b7.json` contains every field of all seven literal
`STANDARD_ANALYSIS_SPECS` declarations in
`src/mtp/workflow/downstream_contracts.py` at commit `f0cb5b7`.

The fixture was extracted on 2026-11-05 using
`extract_issue1606_historical_specs_20261105.py`, retained in the coordination
folder's `scripts/adhoc/`. It calls the reviewer's
`pr1590_spec_strings_by_commit_20261105.py` AST extractor. Historical code is
read with `git show` and never imported. The extractor preserves declaration
order and writes the complete list, rather than patching current declarations.

The comparison test places this historical list into its small synthetic
verified comparison. It checks schema acceptance, preservation of every field,
and unchanged manifest bytes. It does not claim to verify an archived full run
or establish that the historical labels correctly describe the analyses.
