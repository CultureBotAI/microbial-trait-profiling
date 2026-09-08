# Three statistics on one matrix

Three analyses in this project consume the same table,
`data/derived/study_trait_profiles_wide.tsv`, and report significance for
traits. They are not three attempts at the same test, and their p-values are not
on comparable scales. This note records what each one asks, so that a reader who
meets two of them in the same manuscript can tell why they disagree and why
agreement between them is informative rather than redundant.

The shared input has one row per study, experiment, and direction, and one column
per trait; a cell is the number of taxa in that profile carrying that trait.

## What each arm asks

**Composition arm.** Within one experiment, is a trait's share of the signal
higher on the increased side than on the decreased side? The unit of observation
is a matched study/experiment pair, and the quantity tested is a difference
between the two sides of the same pair.

**Enrichment arm.** Taking one profile on its own, does a trait appear more often
than its frequency across the whole corpus would predict? The unit is a
trait-annotation assertion, and the comparison is against a pooled background
rather than against the profile's own counterpart.

**Directional-consistency arm.** Across studies, does a trait more often shift up
than down? The unit is a study, and the quantity tested is the sign of the paired
difference, with its magnitude deliberately discarded.

## Side by side

| | Composition | Enrichment | Directional consistency |
| --- | --- | --- | --- |
| Unit of observation | matched study/experiment pair | trait-annotation assertion | study |
| Transformation | fractions, x1000, log1p, per-trait z-score | none; raw counts enter the urn | fractions, study mean per direction, sign of the difference |
| Test | one-sample t against 0, two-sided | hypergeometric upper tail, one-tailed | exact binomial against 0.5, two-sided |
| Paired | yes, by construction | no; increased and decreased tested separately | yes, then reduced to one sign per study |
| Multiple testing | Benjamini-Hochberg within each output table | Benjamini-Hochberg across every test in the run | Benjamini-Hochberg across the six Table 1 traits |
| Significance rule | none; the figure takes the top ten by p | q at most 0.001 and fold at least 1.5 | none; q reported per trait |
| Typical n per test | hundreds to thousands of pairs | tens of thousands of assertions | roughly 950 to 1,400 studies |

All three restrict to human samples and to experiments carrying both directions,
and all three drop MediaDive growth-medium columns, but they do so at different
points in the pipeline: the enrichment arm filters host species during ingest,
while the other two filter in the figure code. Their MediaDive tests also differ
by one character, matching `mediadive.medium` and `mediadive.medium:`
respectively.

## Why the p-values are not comparable

The three tests count different things, so their sample sizes differ by orders of
magnitude and their p-values cannot be ranked against one another.

The enrichment arm's urn is built from trait-annotation assertions rather than
from taxa. Its committed run reports a mean of 23,686.87 assertions per testable
profile and 556,516 tests in total, and many rows come out at a p-value and a
q-value of literally zero. At that draw size the hypergeometric has effectively
unbounded power, so the fold-change floor rather than the q-value is what selects
the reported set. This is the sampling-unit problem raised in issue #131: one
organism can enter the urn many times. Each profile is also included in its own
background, and the pooled background is restricted to paired experiments, so
"corpus frequency" means the frequency across profiles that survived the pairing
filter.

The composition arm's n is the number of matched pairs, in the hundreds to low
thousands. Its column z-scoring is cosmetic for inference: the column mean and
standard deviation cancel in an increased-minus-decreased difference, so the
t-statistic and its p-value are invariant to that step, which only sets the units
of the reported effect size. The `log1p` step is not cosmetic and does move the
p-values. Its Benjamini-Hochberg correction runs within each output table, so the
pooled, faecal, faeces-excluded, and per-body-site tables are each corrected
independently and no correction spans them.

The directional-consistency arm's n is the number of studies with a non-tied
contrast. Reducing each study to a sign throws away effect magnitude, which makes
it the most conservative of the three and immune to compositional scaling, but it
leaves no effect size beyond the share of studies shifting up. The observed
shares sit within about six points of 0.5, so its significance comes from the
number of studies rather than from the size of the shift.

## Compositionality

Trait fractions within a profile sum to a constant, so two traits compared within
the same profile are negatively correlated by construction. This affects the
guild-against-guild correlation reported alongside the oxygen guild heatmap,
which is why that correlation is described there as partly structural and is not
offered as independent support for the oxygen hypothesis. The composition arm's
paired contrast and the directional-consistency arm's sign test both compare a
trait against itself across the two directions rather than against another trait
within one profile, and so do not share the artefact. Issue #132 records the
related effect in the enrichment arm, where enrichment and depletion mirror each
other by construction.

## Vintage

All three arms currently read an unpinned devel-branch BugSigDB export. Any
comparison between them is only valid where they read the same dump, which is the
subject of issue #115. The numbers quoted in this note come from each arm's
committed outputs at the time of writing and will move when the dump is pinned.

## Provenance of this note

The composition and enrichment descriptions were read from the committed scripts
on `origin/clustered-heatmap-generation` and
`origin/enrichment-analysis-heatmap-cluster` respectively, not supplied by their
authors. They should be confirmed by those authors before any part of this note
goes into the paper. The directional-consistency description is of this branch's
own code and is documented in full in
[`table1_trait_condition_body_site_maps.md`](table1_trait_condition_body_site_maps.md)
and [`guild_condition_body_site_heatmap.md`](guild_condition_body_site_heatmap.md).

This description was drafted with the assistance of a coding agent (Claude Code)
and reviewed by the author.
