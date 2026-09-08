# Oxygen guild composition by condition and body site

Figure: `guild_condition_body_site_heatmap.png`. Produced by
`src/bsdbng/guild_maps.py`, driven by `just guild-maps`. How this analysis
differs from the composition and enrichment arms is set out in
[`three_arm_statistics.md`](three_arm_statistics.md).

## Input corpus

BugSigDB differential-abundance signatures were ingested from the full_dump.csv
export (BugSigDB 2026-04-18_11:12_UTC) into per-study YAML records, flattened to
per-experiment tables, and joined to KG-Microbe taxon-trait assertions to give a
study/experiment/direction by trait count matrix. Each row of that matrix is one
BugSigDB profile: the taxa reported as increased, or as decreased, in one
experiment, summed over the METPO traits those taxa carry.

Two oxygen guilds were defined a priori and were not selected from the data. The
aerobe guild is METPO:1000602 (aerobic), METPO:1000604 (microaerophilic), and
METPO:1000608 (facultatively aerobic); the anaerobe guild is METPO:1000603
(anaerobic) and METPO:1000607 (obligately anaerobic). Guild membership and the
term labels were verified against the merged KG-Microbe node table rather than
taken from prior analysis code.

Three input filters were applied once, upstream, so that the figure operates on a
single corpus definition. Non-human samples were excluded on the per-experiment
host_species field, leaving Homo sapiens only. Experiments were restricted to
those carrying both an increased and a decreased signature, so every contributing
experiment supports a paired contrast. Growth-medium columns, identified by a
mediadive.medium marker in the trait column identifier, were removed from the
trait space before any total was computed. Of 12,122 profile rows read, 2,764
were dropped as non-human and 1,886 as unpaired; 2,334 of 2,651 trait columns
were MediaDive, leaving 317 traits, 3,736 experiments, and 7,472 profiles. Every
experiment in the profile matrix resolved against the experiment metadata table,
so no rows were lost to a failed join.

## Cell values and axes

Each cell of the figure is a mean guild trait fraction. Within a profile, the
guild's summed trait count was divided by that profile's total count across the
retained traits. Because 60.2 percent of studies contribute more than one
experiment, profile fractions were averaged within a study before being averaged
across the studies in a cell, and the companion count matrix reports the number
of contributing studies rather than experiments.

Conditions and body sites were reduced to the most frequent 60 and 20 labels
respectively, with the remainder pooled into an Other bucket. That bucket is
large and its size is reported on every run rather than left implicit: the corpus
holds 645 distinct condition labels and 200 distinct body-site labels, so keeping
the top 60 conditions still pools 585 labels and 46.5 per cent of profiles, and
keeping the top 20 body sites pools 180 labels and 15.3 per cent. The condition
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
the guild trait" is a result rather than an absence. On the run described here 4
conditions were dropped (Atopic asthma, Extraction protocol, Helminthiasis,
Metabolic toxicity), leaving 57 rows.

Condition rows are ordered by hierarchical clustering of their value vectors
under Euclidean distance with Ward linkage, and the tree is cut into six groups
which are separated in the figure by blank rows; no dendrogram is drawn. Masked
cells are treated as zero for the distance calculation only, so the clustering
responds to which body sites a condition was sampled in as well as to the guild
fractions themselves. Average linkage was tried first and rejected: it chained,
leaving one 52-condition group and five near-singletons. The cut assigns
Periodontitis, Smoking behavior, and Rheumatoid arthritis to one group, and
Gingivitis to a neighbouring one; these are the conditions sampled in saliva and
subgingival dental plaque rather than in feces. Alzheimer disease does not join
them. It is the only one of those conditions with faecal as well as salivary
data, and its salivary aerobe fraction is higher on the increased side, 0.146
against 0.101, while the oral conditions run the other way, for example
Periodontitis at 0.081 against 0.119. The separation reflects that difference in
sampling and direction rather than an artefact of the linkage.

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
condition at all. On the run described here the grid is 57 conditions by 84
columns; of 4,788 cells, 432 carry a value, 352 fall below the study threshold,
24 are sampled with no guild signal, and 3,980 were never sampled.

## Statistics

The association between guild and direction was tested on a two-by-two table of
studies. Each study was scored once per guild by the sign of its mean increased
fraction minus its mean decreased fraction, giving cells of guild by shift
direction. Over 1,455 studies, aerobes shifted toward the increased side in 802
studies and toward the decreased side in 619, while anaerobes shifted toward the
increased side in 641 and toward the decreased side in 794; Fisher's exact test
gives an odds ratio of 1.605 with p = 3.15e-10. Counting studies with any
non-zero guild signal instead of the sign of the contrast produces a degenerate
table in which nearly every study occupies every cell and the odds ratio is 1.01;
that formulation was rejected.

Rank correlations were computed across the 108 condition-by-body-site cells that
survived masking, using only cells where both members of a pair were present so
that the series remain aligned. Aerobe and anaerobe fractions correlate
negatively within the same direction, rho = -0.622 on the increased side and
-0.715 on the decreased side. These within-direction correlations are partly
structural: trait fractions within a profile are compositional, so a higher
aerobe share mechanically implies a lower anaerobe share. They should not be read
as independent confirmation of the oxygen hypothesis. The paired shift-direction
odds ratio above does not share that artefact, because it compares a guild
against itself across the two directions rather than comparing the two guilds
within one profile.

The Fisher test is a single test and needs no correction. The six pairwise
Spearman p-values are reported uncorrected, and should be read as descriptive
rather than as a screened set; they are not the basis of any claim in the
manuscript.

## Reproducing

All statistics use scipy (fisher_exact, spearmanr) rather than hand-rolled
implementations. The figure was produced by `src/bsdbng/guild_maps.py`, driven
deterministically by `just guild-maps`, alongside the value matrix, the
study-count matrix, the cell-status table, the row-cluster table, and the
statistics table as tab-separated files. The BugSigDB export timestamp is read
from the dump's leading comment line and is recorded in every output.

The pipeline code and this description were drafted with the assistance of a
coding agent (Claude Code) and reviewed by the author.
