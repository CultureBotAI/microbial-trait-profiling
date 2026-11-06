# Historical comparison declarations (#1606)

`analysis_specs_f0cb5b7.json` contains every field of all seven literal
`STANDARD_ANALYSIS_SPECS` declarations in
`src/mtp/workflow/downstream_contracts.py` at commit `f0cb5b7`.

The fixture was originally extracted on 2026-11-05. The committed
[extractor](../../../scripts/adhoc/extract_comparison_specs_f0cb5b7.py)
reproduces it without importing historical code. It uses the existing GitPython
dependency to read the Git blob equivalent to
`git show f0cb5b7:src/mtp/workflow/downstream_contracts.py`, then parses the
top-level `STANDARD_ANALYSIS_SPECS` mapping with Python's standard-library AST.
It requires seven `AnalysisSpec(...)` calls with all eight fields as literal
strings, preserving declaration and keyword order. JSON uses two-space
indentation and a final newline. Comparison records themselves sort specs by id;
this fixture retains the source declaration order.

From the repository root, with dependencies installed and commit `f0cb5b7`
available locally, reproduce the fixture into a fresh temporary directory and
compare every byte:

```bash
comparison_specs_dir=$(mktemp -d)
uv run python scripts/adhoc/extract_comparison_specs_f0cb5b7.py \
  --repo . --output "$comparison_specs_dir/analysis_specs_f0cb5b7.json"
cmp tests/fixtures/comparison_specs/analysis_specs_f0cb5b7.json \
  "$comparison_specs_dir/analysis_specs_f0cb5b7.json"
```

The extractor bounds the source read to 128 KiB and refuses an existing output
path. These commands leave the committed fixture unchanged; `cmp` exits zero
without output when the reproduced bytes match.

The comparison test places this historical list into its small synthetic
verified comparison. It checks schema acceptance, preservation of every field,
and unchanged manifest bytes. It does not claim to verify an archived full run
or establish that the historical labels correctly describe the analyses.
