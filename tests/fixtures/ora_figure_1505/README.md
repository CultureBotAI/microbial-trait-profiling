# Figure 2 pre-change captures (#1505)

`separate.json` and `combined.json` were captured from the unmodified producer at
`ee819b2` before the streaming change. The inputs are created by `_inputs` in
`tests/test_ora_figure_single_read.py`. The two cases exercise a separate
depletion table and decreased rows already present in the enrichment table.

The test intercepts `Figure.write_image`, checks its output path and scale, and
compares every UTF-8 byte of `Figure.to_json()` against these captures, excluding
the fixture file's terminal newline if present. This pins
the figure data and layout passed to the image renderer without starting a
browser. It does not claim byte identity for a raster image produced by Kaleido.
