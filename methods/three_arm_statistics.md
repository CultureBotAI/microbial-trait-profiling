# Three statistics on one matrix

Three analyses in this project consume the same table,
`data/derived/study_trait_profiles_wide.tsv`, and report significance for
traits. They test different hypotheses with different sampling units. A smaller
p-value in one arm does not establish stronger biological evidence than a larger
p-value in another. This note compares what each asks and what joint
interpretation requires; agreement on shared data is not independent replication.

The shared input has one row per study, experiment, and direction, and one column
per trait; a cell is the number of taxonomically rolled-up trait assertions for
that trait in that profile — one reported taxon can contribute many counts via
its descendant nodes, so cells are not taxon counts (the same sampling-unit
inflation issue #131 records for the enrichment urn).

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
| Owner | Jed | Ilyas | Noel |
| Unit of observation | matched study/experiment pair | trait-annotation assertion | study |
| Transformation | fractions, x1000, log1p, per-trait z-score | none; raw counts enter the urn | fractions, study mean per direction, sign of the difference |
| Test | one-sample t against 0, two-sided | hypergeometric upper tail, one-tailed | exact binomial against 0.5, two-sided |
| Paired | yes, by construction | no; increased and decreased tested separately | yes, then reduced to one sign per study |
| Multiple testing | Benjamini-Hochberg within each output table | Benjamini-Hochberg across every test in the run | Benjamini-Hochberg across the six Table 1 traits |
| Significance rule | none; the figure takes the top ten by p | q at most 0.001 and fold at least 1.5 | none; q reported per trait |
| Typical n per test | 2,335 pairs, identical for every trait | tens of thousands of assertions | 616 to 917 of the 959 studies (ties excluded) |
| Strength | Uses matched contrasts and their magnitudes; screens many traits | Locates profile-specific excess relative to a declared corpus background | Gives each study one directional vote; exact test without a normality assumption for magnitudes |
| Main limitations | Repeated experiments within studies are not independent; denominator and transformation affect the estimand | Correlated annotations undermine the urn sampling model; coverage and background composition affect calls | Discards magnitude, excludes ties, and remains denominator-dependent; six traits selected on overlapping data |

The core paired analyses restrict to human samples and experiments carrying both
directions, at different points in the pipeline. The enrichment branch filters
host species at ingest; the other analyses filter downstream. The same retained
rows do not guarantee the same normalization: Ilyas and Noel remove MediaDive
before computing totals, whereas Jed's inspected paired-screen code excludes
MediaDive from trait selection but includes it in the denominator. Issue
[#174](https://github.com/CultureBotAI/bsdbng/issues/174) tracks harmonization and
the numerical rerun; the historical results below must retain their original
denominator definitions until that work is complete.

## Why the p-values are not comparable

The three tests count different things, so their sample sizes differ by orders of
magnitude and their p-values cannot be ranked against one another.

The enrichment arm's urn is built from trait-annotation assertions rather than
from taxa. Its committed run reports a mean of 23,686.87 assertions per testable
profile and 556,516 tests in total, and many rows come out at a p-value and a
q-value of numerically zero. Large assertion counts can produce extremely small
model-based p-values even for modest enrichment, making the fold floor material
to selection. Numerical zeros do not establish certainty. Calibration still
requires a suitable urn model: this is the sampling-unit problem in #131, where one
organism can enter the urn many times. Each profile is also included in its own
background, and the pooled background is restricted to paired experiments, so
"corpus frequency" means the frequency across profiles that survived the pairing
filter.

The composition arm's n is the number of matched pairs, in the hundreds to low
thousands. For a nonconstant trait, column centering cancels in an
increased-minus-decreased difference, and its common positive scale factor
cancels in the t-statistic. Thus z-scoring sets the reported effect's units while
leaving that t-test unchanged. The `log1p` step does change the tested values and can move the
p-values. Its Benjamini-Hochberg correction runs within each output table, so the
pooled, faecal, faeces-excluded, and per-body-site tables are each corrected
independently and no correction spans them.

The directional-consistency arm's n is the number of studies with a non-tied
contrast. Reducing each study to a sign discards contrast magnitude and avoids a
normality assumption for those magnitudes. It does not remove denominator
dependence, and there is no universal ranking of conservativeness or power across
these tests: they have different null hypotheses, sampling units and correction
families. The share of studies shifting up is an effect estimate for directional
consistency, not for the magnitude of a trait change. Its distance from 0.5 and
its sample size both affect significance.

## Overlap, complementarity and integration

Jed and Noel both compare increased with decreased signatures; Jed retains
contrast magnitude at the experiment level, whereas Noel tests consistency of
study-level direction. Ilyas asks whether either profile is unusual relative to
the background. Both sides can be enriched relative to that background while
their paired difference is small. A large mean difference can also coexist with
inconsistent signs across studies. Such disagreements need not be contradictions.

The methods complement one another as views of relative contrast, background
enrichment and between-study consistency. Their effects, sampling units and
nulls are incompatible with pooling raw effect values or treating their q-values
as interchangeable votes. Shared input and overlapping trait selection prevent
agreement from serving as independent validation. Harmonization must preserve
each original hypothesis as well as expose differences in normalization.

Significance retains its tested level: a body-site result does not establish a
result for every condition within it. A majority of studies containing at least
one significant profile is a descriptive support rule unless a cell-level test
has been calibrated. [#176](https://github.com/CultureBotAI/bsdbng/issues/176) and
[#177](https://github.com/CultureBotAI/bsdbng/issues/177) track the grouping,
display-selection and profile-support questions in PR #173.

Noel owns the integrated evidence table in
[#178](https://github.com/CultureBotAI/bsdbng/issues/178). Each method's result
must retain trait identity, tested condition/site and level, version/input hash,
denominator/background, selection rules, study/pair/profile counts, effect units
and uncertainty type, p/q and correction family, thresholds, and separate
eligibility, missingness and display-mask fields. BH over six selected traits is
different from BH over hundreds of traits or hundreds of thousands of profile-trait tests;
separate corrections do not control FDR for a combined "significant under any
method" set. Untested, all tied, zero annotation, masked and tested nonsignificant
results must survive joins and exports as distinguishable states.

All arms measure annotations attached to curated taxon signatures. They do not
directly measure patient-level trait abundance, pathway activity, oxygen exposure
or causation. [#179](https://github.com/CultureBotAI/bsdbng/issues/179) tracks the
shared biological-interpretation guidance for all contributors.

## Compositionality

The full set of trait fractions in a nonempty profile sums to a constant. This
constrains their covariances, but does not force every pair to correlate
negatively. For example, with counts A=10 and B=20 unchanged, increasing a third
component lowers both A's and B's fractions. Across profiles with different
third-component changes, those two fractions and their paired differences can
correlate positively even though A and B did not change their counts.

Pairing increased and decreased signatures and averaging within study do not
remove this effect. A fraction-based sign test still asks about changes relative
to the total retained annotation signal. It cannot establish an absolute change
in either the trait or its abundance in samples. The observed guild correlations
alone cannot quantify how much of their association comes from closure. Issue
#132 records related interpretation problems in enrichment/depletion.

The [paired log-ratio analysis](table1_paired_log_ratios.md), run with
`just table1-compositional-stats`, adds an explicit reference to each comparison.
For A/B it tests study-level direction consistency of paired log2(A/B) changes.
The profile total cancels, so changing unrelated components does not change that
ratio. All 15 Table 1 pairs are reported, with positive-count eligibility,
excluded experiments/studies, ties, effect magnitudes, exact binomial intervals
for directional proportions, and Benjamini-Yekutieli adjustment across the full
pair family. These are different, conditional hypotheses, not corrected versions
of the six original fraction hypotheses or independent replication.

The [pinned-corpus validation table](table1_log_ratio_validation.md) reports all
15 comparisons alongside the original six fraction directions, with source
fingerprints and the study populations retained by each ratio.

Reference choice, zero exclusions, shared annotations and selection of these six
traits remain interpretation limits. The ratio cannot distinguish an absolute
increase in its numerator from an absolute decrease in its reference. See
[Morton et al. (2019)](https://doi.org/10.1038/s41467-019-10656-5) on reference
frames and [Lovell et al. (2015)](https://doi.org/10.1371/journal.pcbi.1004075)
on associations in relative data. The tests in `tests/test_trait_composition.py`
exercise the denominator-only counterexample and the ratio invariants.

## Vintage

All three arms pin **BugSigDBExports v1.3.1** (Zenodo 10.5281/zenodo.19740674;
issues #115/#130, resolved by PR #150's SHA-verified downloader). The
directional-consistency numbers in this note are from a v1.3.1 run; the
composition and enrichment figures come from those arms' committed v1.3.1
regenerations (`3afeecd` and `c75744a` respectively). Matching the source dump is
necessary for harmonized comparisons, but not sufficient: input fingerprints,
filters, trait identities, denominators and tested populations must also align.

## Provenance of this note

The composition and enrichment descriptions were read from the committed scripts
on `origin/clustered-heatmap-generation` and
`origin/enrichment-analysis-heatmap-cluster` respectively, not supplied by their
authors. They should be confirmed by those authors before any part of this note
goes into the paper. The directional-consistency description is of this branch's
own code and is documented in full in
[`table1_trait_condition_body_site_maps.md`](table1_trait_condition_body_site_maps.md)
and [`guild_condition_body_site_heatmap.md`](guild_condition_body_site_heatmap.md).

The original note acknowledged Claude Code assistance and author review. The
subsequent statistical corrections and log-ratio extension are tracked in
[PR #181](https://github.com/CultureBotAI/bsdbng/pull/181) and require confirmation
by the method owners before manuscript use.
