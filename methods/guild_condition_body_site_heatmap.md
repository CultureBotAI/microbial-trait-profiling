# Oxygen guild composition by condition and body site

Figure: `guild_condition_body_site_heatmap.png`. Produced by
`src/bsdbng/guild_maps.py`, driven by `just guild-maps`. How this analysis
differs from the composition and enrichment arms is set out in
[`three_arm_statistics.md`](three_arm_statistics.md).

Interpret all outputs using the [shared guide](biological_interpretation.md).
Guild labels summarize annotations, not measured oxygen exposure or respiration.
`inc`/`dec` refer to curated comparison directions; human host and pairing do not
verify disease/control orientation (#188). Heatmap cells are descriptive means
without cell-level confidence intervals; use the support and status tables, and
do not assign the pooled statistical test to individual cells.

## Input corpus

BugSigDB differential-abundance signatures were ingested from the full_dump.csv
of the **pinned BugSigDBExports v1.3.1 release** (Zenodo 10.5281/zenodo.19740674;
internal timestamp 2026-04-24_00:41_UTC, stamped on every output row) into
per-study YAML records, flattened to per-experiment tables, and joined to
KG-Microbe taxon-trait assertions to give a study/experiment/direction by trait
count matrix. Each row of that matrix is one BugSigDB profile: the taxa reported
as increased, or as decreased, in one experiment, summed over the METPO traits
those taxa carry (a taxonomic rollup, so one reported taxon can contribute many
assertion counts).

Two oxygen guilds were defined a priori and were not selected from the data. The
aerobe guild is METPO:1000602 (aerobic), METPO:1000604 (microaerophilic), and
METPO:1000608 (facultatively aerobic); the anaerobe guild is METPO:1000603
(anaerobic) and METPO:1000607 (obligately anaerobic). The term labels were
checked against the merged KG-Microbe node table rather than taken from prior
analysis code, and the check is re-runnable: `just verify-trait-terms` looks every
grouped identifier up in `data/kg-microbe/merged-kg_nodes.tsv` and fails on a
missing or changed label. It verifies the node table's MD5 against
`data/kg-microbe/provenance.json` **before** reading any label, and stops if they
differ, because otherwise a run against a different KG-Microbe build reports a
term as missing when the real finding is that the build is not the pinned one
(#289). `--allow-unpinned-nodes` reports the labels anyway, and says so. On the
pinned build (nodes MD5 `dd76b4656f9341db25a7c4afa96cd58e`, fetched 2026-04-21)
all twelve grouped identifiers carry their expected label, `pathways:arsenate_reduction`
included. A label
check cannot validate the grouping itself. Placing microaerophilic and
facultatively aerobic organisms in the aerobe guild is an a priori decision: both
tolerate oxygen, but neither is an obligate aerobe, and a guild of strict aerobes
alone would be a different quantity. METPO:1000608 (facultatively aerobic) is also
one of the six Table 1 traits, so this figure and the Table 1 figures share
evidence and are not independent of each other.

Three input filters were applied once, upstream, so that the figure operates on a
single corpus definition. Non-human samples were excluded on the per-experiment
host_species field, leaving Homo sapiens only. Experiments were restricted to
those carrying both an increased and a decreased signature, so every contributing
experiment supports a paired contrast. Growth-medium columns, identified by a
mediadive.medium marker in the trait column identifier, were removed from the
trait space before any total was computed. Of 5,629 profile rows read under the
`count_input_16s_wms_v1` profile policy, from 1,047 studies, none was dropped as
non-human, because the matrix holds human profiles only (#656), and 1,115
experiments were dropped as unpaired; 2,290 trait columns were MediaDive,
leaving 317 traits (the `off` scope; 344 under `phenotype`), 2,257 experiments,
4,514 profiles, and 932 studies. Every
experiment in the profile matrix resolved against the experiment metadata
table, so no rows were lost to a failed join.

The figure, its cell-status table and its McNemar then leave out every
experiment pair whose retained-annotation total is zero on a side (#783). That
leaves 55 fewer pairs: 2,202 experiments, 4,404 profiles and 922 studies. This is
the shared comparison basis of [`three_arm_statistics.md`](three_arm_statistics.md),
which the integrated evidence table's tests and the calibration already used.
The repository owner chose it for the figure on 2026-10-22. Such a pair's guild
share is undefined, not zero. Scored as 0, it lowered cell means and shaded
cells dark grey as if annotated profiles lacked the guild. It also gave this
figure's McNemar a different population from the integrated table's, so the two
reported different results for one run.

The host filter drops nothing on the matrix identified below, which is written
from human-host 16S and shotgun profiles only (#656, #742). Earlier revisions of
this note named a matrix that was not human-only, sha256
`a0c28978c5e58c153afba517f8a9fb9a50bd62d5dddefb78f583e0cb2aedf334` (2,653
header columns, 2,333 of them MediaDive, over 7,111 profile rows), and there the
filter dropped 1,482 rows (#222).

The matrix these counts were measured on is
`data/derived/study_trait_profiles_wide.tsv`, sha256
`519a095b09e179dc92fb51f9d892984a52dd100dece4dadf6c94b036a5323e6c`, the `off`-scope
`study_trait_profiles_wide.tsv.gz` committed until the final rerun at 511c891 replaced
it with the `phenotype` table (`43b2e230…`: 2,637 header columns, 344 retained),
decompressed: 2,610 header columns (3 identifier, 2,290 MediaDive, 317 retained) over 5,629 profile rows, joined
from BugSigDB `2026-04-24_00:41_UTC` and KG-Microbe release 20260120 (nodes md5
`dd76b4656f9341db25a7c4afa96cd58e`), as the download sidecar
`data/kg-microbe/provenance.json` declares it rather than as measured when the
matrix was built. The experiment metadata table is sha256
`16271237a70f61c8195d9f2acdd16f71e88f65859160f95bb734626bc1328e06` and the full
dump sha256 `af84bfca7bfe76913a20c1fedbaadbfad4e4bf78605a2cfd7d00eb0c0fdd5788`.
Naming the dump does not identify the matrix, so the run writes this record to
`corpus_identity.json` beside its outputs (#222).

The counts in this note were measured on those inputs with the code main
carries as commit 35ece98 (#808): the corpus cascade, the zero-total exclusion,
the multi-experiment share, the pooled labels, the row cuts and clusters, the
cell census, the McNemar and the rank correlations. The figure's numbers from
before #783, given for comparison, are from main's code at 65850b0 on the same
inputs.


## Cell values and axes

Each cell of the figure is a mean guild trait fraction. Within a profile, the
guild's summed trait count was divided by that profile's total count across the
retained traits — a share of taxonomically rolled-up trait assertions, not of
taxa. Because 511 of the 922 studies (55.4 per cent) contribute more than one
experiment, profile fractions were averaged within a study before being averaged
across the studies in a cell, and the companion count matrix reports the number
of contributing studies rather than experiments.

Conditions and body sites were reduced to the most frequent 60 and 20 labels
respectively, with the remainder pooled into an Other bucket. That bucket is
large and its size is reported on every run rather than left implicit: the corpus
holds 429 distinct condition labels and 147 distinct body-site labels, so keeping
the top 60 conditions still pools 369 labels and 40.4 per cent of profiles, and
keeping the top 20 body sites pools 127 labels and 13.4 per cent. The condition
tail is long and shallow, with several hundred labels attached to two profiles
each, so no practical cut makes Other small; the figure should be read with that
in mind. Note also that BugSigDB records study variables such as Diet, Age, and
Extraction protocol in the same Condition field as diseases, so the condition
axis is not purely a disease axis.

Cells supported by fewer than three studies were masked rather than drawn, so
that a single study is never rendered as a result. Rows are conditions and
columns are body sites crossed with the four guild-by-direction combinations, so
each body site contributes the two-by-two guild by direction table unrolled into
four adjacent columns.

The two guilds are drawn in separate hues, aerobe on a red-orange ramp and
anaerobe on a green-yellow ramp, each with its own colour bar, but both ramps
span the same range, from zero to the larger of the two observed maxima. Keeping
the range common means a cell in one guild is directly comparable in magnitude
to a cell in the other, while the hue still identifies which guild is being read.

## Row selection and clustering

A condition whose every cell was masked would render as a blank row and is
dropped before clustering; the dropped labels are listed on each run. A row that
keeps only dark-grey cells is retained, because "sampled, but no recorded guild
annotation signal" distinguishes mapped zero signal from unsampled cells.
On the run described here 9
conditions were dropped (Atopic asthma, Cervical cancer, Head and neck squamous
cell carcinoma, Helminthiasis, Metabolic toxicity, Phenylketonuria, Population,
Psoriasis vulgaris, Smoking behaviour measurement), leaving 52 rows.

Condition rows are ordered by hierarchical clustering of their value vectors
under Euclidean distance with Ward linkage, and the tree is cut into six groups
which are separated in the figure by blank rows; no dendrogram is drawn. Masked
cells are treated as zero for the distance calculation only, so the clustering
responds to which body sites a condition was sampled in as well as to the guild
fractions themselves. Average linkage was tried first and rejected: it chained,
leaving one dominant group and a few near-singletons. On the run described here
the cut groups the oral and salivary conditions together (Periodontitis,
Rheumatoid arthritis, Oral squamous cell carcinoma and Smoking behavior form one
cluster with Cesarean section). It puts the two cervicovaginal conditions (Human
papilloma virus infection and Cervical glandular intraepithelial neoplasia) in a
cluster of their own, so the grouping is driven by
sampling site as much as by guild fraction — which is the intended behaviour,
but means the clusters must not be read as guild-biology groupings on their own.
(An earlier draft narrated a Periodontitis/Alzheimer separation from a
pre-v1.3.1 run. Since #783 Periodontitis is back inside the top-60 cut, in the
oral cluster above, which is cluster 5 in the row-cluster table; Alzheimer's
disease is still outside it and is pooled into Other.)

## Cell states

Cells fall into four states, recorded per cell in a companion status table so the
distinction is machine-readable rather than only visual. A cell carries a value
when at least the minimum number of studies support it and at least one study has
non-zero guild signal. A cell is filled dark grey when the body site was sampled
for that condition but there is no recorded guild annotation signal under the
selected sources and mappings. Previously such a cell was drawn as a fraction
of zero and was indistinguishable from a genuinely small fraction. This state
does not establish biological absence. Since #783 every profile behind a dark-grey cell carries
some retained annotation. Before, a profile with none at all scored share 0 and
shaded cells too: 8 of this corpus's 20 dark-grey cells were backed only by
such profiles. A cell is left blank when the site was sampled by fewer
studies than the threshold, and blank when the site was never sampled for that
condition at all. On the run described here the grid is 52 condition rows (51
conditions and the pooled Other) by 80 columns, which are 20 body-site blocks: 19
sites and the pooled Other. Tongue is inside the top-20 cut, but nothing in its
block is drawn, so the block is left out. Of 4,160 cells, 340 carry a value, 302
fall below the study threshold, 10 are sampled with no guild signal, and 3,508
were never sampled.

## Statistics

The association between guild and direction was tested on studies scored once
per guild by the sign of the mean increased fraction minus the mean decreased
fraction. Over 922 studies, the aerobe guild's annotation share was higher in
the increased signatures than in the decreased ones in 520 studies and lower in
392, while the anaerobe share was higher in 402 and lower in 519. Because the
**same studies provide both signs**, the two guilds form matched pairs, not
independent samples, so the test is an **exact McNemar** — a two-sided binomial
on the discordant studies. Of the 912 studies with both signs (697 discordant),
410 have the higher aerobe share and the lower anaerobe share in the increased
signatures and 287 the reverse, giving a paired odds ratio of 1.429 with
p = 3.6e-06. That is the same test on the same population as the
integrated evidence table's McNemar row. Before the zero-total exclusion (#783)
this figure reported 410 and 284, p = 2.0e-06, over 932 studies. (An earlier draft applied Fisher's exact test to the
marginal two-by-two; that assumes independent rows, and on this data it was
anticonservative by roughly three orders of magnitude and inflated the odds
ratio by a fifth.) Counting studies with any non-zero guild signal instead of
the sign of the contrast produces a degenerate table in which nearly every study
occupies every cell; that formulation was rejected.

The paired odds ratio mitigates study-level dependence, but it does **not**
escape compositionality: both signs derive from shares of the same profile
denominator. An increase in one guild's counts can lower the other guild's share
without a change in the latter's counts. This does not force every pair of
components to correlate negatively. On this corpus 76 per cent of the scored
studies (697 of 912) carry opposite signs for the two guilds; that observation
alone does not quantify the structural contribution. The odds ratio reads as
"discordant studies more often pair a higher aerobe share with a lower anaerobe
share in the increased signatures than the reverse", not as independent
confirmation of the oxygen hypothesis. The matched-corpus Table 1-trait analog
directional-shift figure tests each trait against itself across directions, but
remains dependent on the same profile denominators. The [paired log-ratio analysis](table1_paired_log_ratios.md)
illustrates reference-explicit tests over the individual Table 1-named traits; it does
not independently validate the oxygen guild hypothesis.

Rank correlations are computed across studies, one observation per study, using
each study's mean guild fraction per direction and its paired contrast (the mean
increased minus mean decreased fraction). An earlier version correlated the 86
condition-by-body-site cells that survived masking; one study reaches many cells,
so that counted dependent cells as independent observations, and it is replaced
here by the per-study template the matched-corpus Table 1-trait analog figure already uses (#165). A study
that carries no oxygen annotation on either side has a missing share, not a
measured one, so it is excluded from the rank correlations (#289). The McNemar
drops such a study as a tie, so both statistics describe the same population.
On this corpus the five such studies were all zero-total pairs, which the
exclusion above already removes, so all 922 studies are ranked. Over those 922
studies, aerobe and anaerobe fractions correlate negatively within the same
direction, rho = -0.719 on the increased side and -0.740 on the decreased side,
and the two guilds' per-study contrasts correlate at rho = -0.737. These are
largely the structural artefact described above - two shares of one denominator -
and are not to be read as support. Studies are still not fully independent where
they share cohorts, which this computation does not model.

The McNemar test is a single test and needs no correction. The seven Spearman
p-values are reported uncorrected, and should be read as descriptive rather than
as a screened set; they are not the basis of any claim in the manuscript.

## Reproducing

All statistics use scipy (binomtest, spearmanr) rather than hand-rolled
implementations. The figure was produced by `src/bsdbng/guild_maps.py`, driven
deterministically by `just guild-maps`, alongside the value matrix, the
study-count matrix, the cell-status table, the row-cluster table, and the
statistics table as tab-separated files. The statistics table's first row and a
subtitle line state how many zero-total pairs were excluded (#783). The
integrated evidence table inventories the same cells from the same builder,
`build_figure_matrix`. A test runs both command lines on one corpus and requires
the two cell tables (#770) and the two McNemars (#783) to agree. The BugSigDB export timestamp is read
from the dump's leading comment line and is recorded in every output. PNG export
runs through kaleido 1.x, which drives a headless Chrome. For a reproducible
figure run, run `just get-chrome` once first: it downloads (~170 MB) a
version-pinned Chrome for Testing build that kaleido then prefers over any
system Chrome, which would otherwise be whatever version the machine carries
(plotly's `plotly_get_chrome` is equivalent). Without any Chrome the recipe
fails only at the figure-writing step, with every backing table already written
to disk; the corpus-cascade report prints to stdout either way. None of the
statistics depend on Chrome.

The pipeline code and this description were drafted with the assistance of a
coding agent (Claude Code) and reviewed by the author.
