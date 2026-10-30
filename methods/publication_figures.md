# Publication figures from completed analysis tables

`just publication-figures` consumes explicit comparison, three-test and full
Delta-z-derived sidecars. Pass `--comparison-run`, `--three-test-run` and
`--derived-run` as directories or manifest files, and pin the two sidecar
manifests with `--three-test-sha256` and `--derived-sha256`. The native comparison
and completed-table verifier checks the bound core and arms; the derived reader
checks its historical producer inventory, parent binding and complete output
inventory. The consumer hashes and parses the same table buffers, then checks
source identities again before publication. It never reruns the analysis.

The required `--output` is a fresh directory outside the checkout and all source
bundles. Existing destinations, including dangling symlinks, are refused. The
shared durable-output policy applies; `--allow-temporary-output` is for disposable
tests only. The viz extra supplies Matplotlib. `--code-commit` may pin the
executing checkout; dirty development renders are explicitly recorded as dirty.

Each figure has a combined vector PDF, individual page PDFs and SVGs, PNGs at
`--dpi` (at least 300; default 300), and `provenance.json`. The record names every
input with SHA256 and size, source and consumer commits, run identifiers,
contrast population, direction basis, analysis bases and rollup. It records
page row identifiers, mark counts, selection settings, rendering environment,
and hashes of the exported files. The final `figures.json` binds the individual
provenance records. Publication uses an exclusive atomic directory rename.

## Figures and estimands

`--figure` can be repeated to select `forest`, `quadrant` or `three-test`; omitting
it renders all three. All trait names come from verified feature metadata, with
qualified identifiers alongside the labels. A label is an organism annotation,
not a measured physiological activity in a study participant. Label wording is
preserved even when an ontology label is awaiting a later terminology update.

The forest selects positive and negative tails using the persisted native ranks,
retaining rank ties at `--top-n` (default 5). Its default `--ranking
comparison_mean` reproduces the historical comparison-weighted ordering. A
comparison-weighted mean and an equal-study mean are different estimands: the
comparison mean is drawn separately from the study mean and its own 95% confidence
interval. No interval is transferred between these quantities. `--ranking
study_mean` explicitly changes selection to the persisted study-mean ranks.
Selection is descriptive; it does not define a new multiple-testing family.

The quadrant scatter uses the saved paired Delta-z values for `METPO:1000608`
and `METPO:1000607`. It shows every comparison, including repeated coordinates
and exact-zero axis points, and a second panel for literal `body_site == "Feces"`
(stool). Multi-site strings are excluded from that panel. Shared limits keep the
two panels comparable. No points are clipped or sampled. Comparisons may share
studies and controls, and the two composition-derived coordinates can be
structurally correlated; the scatter does not establish independent biological
replication or causal physiology.
Counts appear in a separate strip outside the data panels. Its key identifies
the signs of the horizontal and vertical coordinates; quadrant counts exclude
exact-axis comparisons, whose counts are shown separately.

The table figure paginates the full retained trait family, in stable identifier
order, rather than selecting significant rows. It preserves native Noel study
votes and BY q, the study-mean Delta-z estimate, interval and its own BY q, and
ORA direction counts. ORA remains an uncalibrated directional column and does not
decide a headline call. Missing or untestable values remain missing. The two
statistical families retain their recorded definitions; their q values are not
compared across tests.

## Verification and scope

`just test-publication-figures` checks input binding, complete page inventories,
repeated scatter points, rank ties, estimator-specific intervals, missing values,
wording, source protection, and publication failures. Synthetic tests precede
rendering against archived evidence. Final-size PDF/PNG inspection and the
external design-check CLI, when already installed, supplement these tests.
`scripts/adhoc/issue1046_design_check.sh` accepts a completed figure bundle,
a fresh durable audit directory and the existing design-check executable. It
checks each page at journal double-column width, validates the reports and
binds their PDF/PNG hashes back to the exports. Visual findings can override
an automated pass; a rejected canary must be corrected and checked again
before rendering the full set.
The bundle audit compares table row identities against the hash-bound
`per_trait_evidence.tsv` input and checks each page's row and mark counts. It
checks quadrant marks and population-prefixed comparison identities against
each population's declared comparison count, including axis points. The table
audit therefore needs the recorded input-table path as well as the figure
bundle. Automated design checks supplement manual inspection: a passing
report alone does not rule out text covering a data mark.

This first slice supplies the builder and these three figures. Other publication
figures and the final paper-basis rendering remain under issue #1046. Development
on a historical verified bundle does not convert it into the final paper basis.

## Slice 2: declared ranking, candidate triples and trait overlap

`just publication-figures-slice2` extends the completed-table contract with three
figures. The comparison, three-test and derived paths and SHA256 pins use the
same arguments as slice 1. Repeat `--figure` to select `rank-sensitivity`,
`candidate-triples`, or `trait-overlap`; omitting it selects all three. Each figure
exports a complete `audit.tsv`, combined/page PDFs, page SVGs and PNGs at at least
300 dpi, and input/output/producer provenance. Source buffers are hash-bound and
rechecked before the exclusive atomic publication. Core identity, population,
direction and rollup must match the parent unit; a sidecar from another unit is
refused. Scientific producers are never rerun by figure code.

The forest's default experiment-weighted selection is unchanged. The rank
sensitivity uses the saved `experiment_mean_delta_z`, `study_mean_delta_z` and
their saved ascending/descending minimum-tie ranks in `per_trait_evidence.tsv`.
It refuses missing columns. It neither averages the study ledger nor reranks
estimates. The TSV retains the whole trait family; the rendered table shows the
complete union of the positive and negative top-five sets under the two saved
rankings (or `--top-n`), including boundary ties and study-only selections.
Each rank belongs to its own mean's sign. Equal-study and equal-comparison
weighting describe different estimates; the sensitivity is descriptive and
does not introduce an inferential family.

For candidate triples, pass `--scan-run` and `--scan-sha256`, pinning the native
`core_consumption.json`. The figure consumes the saved distribution and compressed
triple tables for cellulose, retaining every tested candidate and identifying the
protocol's PD/HIV/CRC triple with its saved ranks. The panels display nominal U6
p-values for mean pairwise Jaccard distance and three-way core across candidate
triples. This is a **candidate-triple distribution**, not a distribution of null
draws. There is no simulated-null export to plot. Counts, excluded/untestable
states and the candidate reference family remain explicit. Candidate triples
can share studies and randomization draws; these panels are neither an independent
replication count nor a uniform-p-value calibration check.

Historical protocol-v1 F scans retain their version and historical scope labels.
Protocol-v2 scans require their declared companion scope, passed with
`--scan-companion-run` and `--scan-companion-sha256`. Both scopes appear in the
same figure bundle with separate identities and audit rows; the companion has
its own core run, but must match the protocol, source commit, population and
direction. The phenotype-primary designation was chosen after both scopes had
been inspected, and does not retrospectively change the F protocol.

For overlap, pass `--matrix-run`, `--matrix-sha256` (the native
`core_consumption.json` pin), and the explicit `--trait-list`. The declared list
must have nine keys in the recorded order and bind to the native output metadata.
Every one of the 36 unordered pairs is retained. Cells show the saved taxon-carrier
Jaccard fraction across all matrix taxa, with verified feature labels and IDs.
An empty union is missing, while a nonempty union with no common carrier is zero.
Diagonal cells are structural, not observations. No tenth trait is selected,
and `METPO:1000605` is not restored to the declared union list. Co-annotation
does not demonstrate co-expression, physiological activity or taxon interaction.
Historical F raw-policy outputs remain identified as such; selecting nine keys
does not convert F's matrix into the later union-policy matrix.

Development renders use F evidence only. Final P figures require approval of
the completed sidecars and a separate figure-producer pin; they are not supplied
by a successful F render.
