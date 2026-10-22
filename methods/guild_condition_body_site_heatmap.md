# Oxygen guild composition by condition and body site

Figure: `guild_condition_body_site_heatmap.png`. Produced by
`src/bsdbng/guild_maps.py`, driven by `just guild-maps`. How this analysis
differs from the composition and enrichment arms is set out in
[`three_arm_statistics.md`](three_arm_statistics.md).

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
trait space before any total was computed. Of 7,111 profile rows read under the
`count_input_16s_wms_v1` profile policy, 1,482 were
dropped as non-human (leaving 5,629 human profiles from 1,047 studies) and 1,115
experiments as unpaired; 2,333 trait columns were
MediaDive, leaving 317 traits, 2,257 experiments, 4,514 profiles, and 932
studies. Every experiment in the profile matrix resolved against the experiment
metadata table, so no rows were lost to a failed join.

An earlier revision of this note said the v1.3.1 wide matrix is already human-only
upstream and that the host filter drops 0 rows. That is not true of the matrix
identified below, where the host filter drops 1,482 rows (#222).

The matrix these counts were measured on is
`data/derived/study_trait_profiles_wide.tsv`, sha256
`a0c28978c5e58c153afba517f8a9fb9a50bd62d5dddefb78f583e0cb2aedf334`: 2,653
header columns (3 identifier, 2,333 MediaDive, 317 retained) over 7,111 profile
rows, joined from BugSigDB `2026-04-24_00:41_UTC` and a KG-Microbe build whose
nodes are md5 `dd76b4656f9341db25a7c4afa96cd58e` and which records no release
tag. Naming the dump does not identify the matrix, so the run writes this
record to `corpus_identity.json` beside its outputs (#222).


## Cell values and axes

Each cell of the figure is a mean guild trait fraction. Within a profile, the
guild's summed trait count was divided by that profile's total count across the
retained traits — a share of taxonomically rolled-up trait assertions, not of
taxa. Because 56.2 percent of studies contribute more than one
experiment, profile fractions were averaged within a study before being averaged
across the studies in a cell, and the companion count matrix reports the number
of contributing studies rather than experiments.

Conditions and body sites were reduced to the most frequent 60 and 20 labels
respectively, with the remainder pooled into an Other bucket. That bucket is
large and its size is reported on every run rather than left implicit: the corpus
holds 451 distinct condition labels and 152 distinct body-site labels, so keeping
the top 60 conditions still pools 391 labels and 41.1 per cent of profiles, and
keeping the top 20 body sites pools 132 labels and 14.1 per cent. The condition
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
keeps only dark-grey cells is retained, because "sampled, but no taxon carried
the guild trait" is a result rather than an absence. On the run described here 9
conditions were dropped (Atopic asthma, Cervical cancer, Head and neck squamous
cell carcinoma, Helminthiasis, Metabolic toxicity, Phenylketonuria, Population,
Psoriasis vulgaris, Smoking behaviour measurement), leaving 52 rows.

Condition rows are ordered by hierarchical clustering of their value vectors
under Euclidean distance with Ward linkage, and the tree is cut into six groups
which are separated in the figure by blank rows; no dendrogram is drawn. Masked
cells are treated as zero for the distance calculation only, so the clustering
responds to which body sites a condition was sampled in as well as to the guild
fractions themselves. Average linkage was tried first and rejected: it chained,
leaving one dominant group and a few near-singletons. On the pinned v1.3.1 run
the cut groups the salivary conditions together (Rheumatoid arthritis, Oral
squamous cell carcinoma, and Smoking behavior form one cluster) and puts the two
cervicovaginal conditions (Human papilloma virus infection and Cervical
glandular intraepithelial neoplasia) in their own, so the grouping is driven by
sampling site as much as by guild fraction — which is the intended behaviour,
but means the clusters must not be read as guild-biology groupings on their own.
(An earlier draft narrated a Periodontitis/Alzheimer separation from a
pre-v1.3.1 run; those conditions no longer survive this corpus's row cuts.)

## Cell states

Cells fall into four states, recorded per cell in a companion status table so the
distinction is machine-readable rather than only visual. A cell carries a value
when at least the minimum number of studies support it and at least one study has
non-zero guild signal. A cell is filled dark grey when the body site was sampled
for that condition but no taxon in those profiles carried any trait of that
guild; previously such a cell was drawn as a fraction of zero and was
indistinguishable from a genuinely small fraction, which understated how often a
guild is simply absent. A cell is left blank when the site was sampled by fewer
studies than the threshold, and blank when the site was never sampled for that
condition at all. On the run described here the grid is 52 conditions by 84
columns; of 4,368 cells, 344 carry a value, 308 fall below the study threshold,
24 are sampled with no guild signal, and 3,692 were never sampled.

## Statistics

The association between guild and direction was tested on studies scored once
per guild by the sign of the mean increased fraction minus the mean decreased
fraction. Over 932 studies, aerobes shifted toward the increased side in 525
studies and toward the decreased side in 390, while anaerobes shifted toward the
increased side in 405 and toward the decreased side in 520. Because the **same
studies provide both signs**, the two guilds form matched pairs, not independent
samples, so the test is an **exact McNemar** — a two-sided binomial on the
discordant studies. Of the 915 studies with both signs (694 discordant), 410 shift aerobe up
while anaerobe shifts down and 284 the reverse, giving a paired odds ratio of
1.444 with p = 2.0e-06. (An earlier draft applied Fisher's exact test to the
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
studies (694 of 915) carry opposite signs for the two guilds; that observation
alone does not quantify the structural contribution. The odds ratio reads as
"the discordance runs more
often in the aerobe-up direction than the reverse", not as independent
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
here by the per-study template the matched-corpus Table 1-trait analog figure already uses (#165). Five of
the 932 studies carry no oxygen annotation on either side. Their zero is
missingness rather than a measured share, so they are excluded from the rank
correlations, which leaves 927 (#289). The McNemar already dropped them as ties,
so both statistics now describe the same population. Over those 927 studies, aerobe and anaerobe fractions
correlate negatively within the same direction, rho = -0.714 on the increased side
and -0.708 on the decreased side, and the two guilds' per-study contrasts
correlate at rho = -0.724. These are
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
statistics table as tab-separated files. The BugSigDB export timestamp is read
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
