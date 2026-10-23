# Three statistics on one matrix

Three analyses in this project report significance for traits from the same
curated signatures. The composition and directional-consistency arms consume
`data/derived/study_trait_profiles_wide.tsv`; the enrichment arm consumes its
binary companion, `data/derived/study_trait_profiles_ora_wide.tsv`, which records
for each profile how many distinct reported taxa bear each trait and how many taxa
it reports. They test different hypotheses with different sampling units. A
smaller p-value in one arm does not establish stronger biological evidence than a
larger p-value in another. This note compares what each asks and what joint
interpretation requires; agreement on shared data is not independent replication.

The shared input has one row per study, experiment, and direction, and one column
per trait; a cell is the number of taxonomically rolled-up trait assertions for
that trait in that profile — one reported taxon can contribute many counts via
its descendant nodes, so cells are not taxon counts (the sampling-unit inflation
that #131 recorded for the earlier, assertion-level enrichment urn).

## What each arm asks

**Composition arm.** Within one experiment, is a trait's share of the signal
higher on the increased side than on the decreased side? The unit of observation
is a matched study/experiment pair, and the quantity tested is a difference
between the two sides of the same pair.

**Enrichment arm.** Taking one profile on its own, does a trait appear more often
than its frequency across the whole corpus would predict? In the production run
(`just analyze`) the unit is a distinct reported taxon, drawn against a background
of the 5,829 BugSigDB-reported taxa in the trait matrix, and the comparison is
against that background rather than against the profile's own counterpart.
Earlier runs counted trait-annotation assertions instead (#131). The Table 1
native reruns in [`table1_native_influence.md`](table1_native_influence.md)
still follow that earlier kernel. Their `table1-native-results/native_ora_families.tsv`
draws assertion counts: 4,517 testable profiles against a background of
106,993,585 assertions. So the two runs' urns, and their test counts, differ.

**Directional-consistency arm.** Across studies, does a trait more often shift up
than down? The unit is a study, and the quantity tested is the sign of the paired
difference, with its magnitude deliberately discarded.

## Shared comparison basis

The comparison specification is `three-framework-comparison-v1`. Its paired
question is: **within a declared set of curated human signature pairs, how does
trait annotation representation differ between increased and decreased
signatures?** Background enrichment supplies a related, separately labeled
question about each signature. The comparison produces common descriptive
estimates from the shared counts and retains each framework's original test as
its own evidence record. Computing the same descriptive estimate through three
exporters is a consistency check, not three independent findings.

This section specifies the harmonized comparison to implement in
[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178). It does not certify that
the historical runs below already meet the specification. No new significance
cutoff, minimum effect, support cutoff, ranking score or top-k value is selected
here; the first comparison displays every declared candidate in manifest order.

### Corpus, features and reference

| Component | Shared definition |
| --- | --- |
| Measurement | Counts of trait-annotation assertions attached to curated taxon signatures, including the documented taxonomic roll-up. These are not patient measurements or counts of independent organisms. |
| Source manifest | Pin the profile matrix, experiment metadata, source/code versions, trait mapping and retained-column list, with hashes. Materialize the base list of `(study_id, experiment_id)` pairs and both direction records; matching row counts alone is insufficient. |
| Base population | Human experiments with exactly one increased and one decreased profile and consistent pair metadata. Reject ambiguous keys and malformed, missing or negative counts; never silently coerce them to zero. Record all filtering and available annotation-coverage information. Unknown coverage remains flagged. |
| Analysis scope | The initial comparison pools the full base population. A body-site or condition-by-body-site analysis is a separately identified scope using exact source labels and its own eligible-pair list. Missing labels remain explicit. Heatmap top-N cuts, masking and an Other display group never define the analysis population. |
| Candidate families | Keep the six traits named in Table 1 and all 15 unordered pairs as separate feature types for the matched-corpus Table 1-trait analog sensitivity analysis. Freeze the source qualifiers and ID mapping in `TABLE1_TRAITS`; orient every ratio in that declared trait order, with an explicit numerator and reference. A larger trait screen requires its own candidate manifest: the directional arm's full-trait screen declares one in [`trait_screen_manifest.json`](trait_screen_manifest.json) and is corrected in its own family, separate from the six (#255). |
| Fraction denominator | For each profile, sum every column in the declared non-MediaDive trait universe, including traits outside Table 1. That is 317 retained columns in the matrix identified below, built without BacDive strain records, with 2,333 MediaDive columns excluded; the strain rollup scope sets the universe (see below). Selecting six traits or a ratio must not redefine this total. |
| ORA background | For the harmonized baseline, fix the background to all valid profiles in the full paired human base population, using that same retained-column universe. Record its trait totals, grand total and contributing profile IDs. Do not rebuild it from a displayed stratum or a ratio's positive-count subset. Alternative backgrounds are separately versioned sensitivity scenarios. |
| Weighting | Average eligible experiment contrasts within a study, then give each eligible study equal weight within the declared scope. Pool a study's experiments directly across conditions; do not count each condition-study group as another study. Identify shared cohorts where known when assessing independence and uncertainty. |

Source validity and annotation coverage are tracked in
[#184](https://github.com/CultureBotAI/microbial-trait-profiling/issues/184). A numeric zero can be a
valid observed annotation count without proving biological absence. A profile
with zero total annotations remains a documented base record, but its fraction
and enrichment fold are undefined. If the entire ORA background total is zero,
the background is undefined; if a trait's background total is zero, its fold
enrichment is undefined. Preserve these states rather than substituting a fold.
The background remains weighted by profiles and their assertion counts, while
the descriptive comparison weights studies equally. Fixing that definition does
not validate the assertion-level urn model in
[#131](https://github.com/CultureBotAI/microbial-trait-profiling/issues/131).
This harmonized-baseline background is a proposed definition for the shared
comparison. It is not the background of the production enrichment run, which
draws distinct taxa against every taxon in the trait matrix (see
[Why the p-values are not comparable](#why-the-p-values-are-not-comparable)).

### Strain rollup scope and the trait universe

Which BacDive strain records reach the matrix is the `strain_rollup` scope
([#788](https://github.com/CultureBotAI/microbial-trait-profiling/issues/788)). The core
workflow takes it as a parameter and records it in every run. The legacy matrix builder
(`just taxa-trait-matrix`) and the Table 1 contribution ledger take it as an option, and
both default to the primary scope like the workflow. A strain's
admitted trait edges count at the NCBITaxon node it is a `biolink:subclass_of`, once per
node, so the matrix still counts trait-bearing NCBITaxon nodes. The three scopes have
distinct roles and distinct trait universes. Sizes below are for the pinned KG (edge md5
`555ecf5b…`, node md5 `dd76b465…`) and the 5,829 full-source reported taxa of the core
run at main 8b280a3.

| Scope | Role | Retained traits (m) | Retained count mass, all reported taxa | What strain records add |
| --- | --- | ---: | ---: | --- |
| `phenotype` | Primary; the default of the core workflow and the legacy builder | 344 | 1,335,025 | 27 bare METPO phenotype classes, facultatively anaerobic among them; biosafety levels and pathogen labels are excluded |
| `off` | Sensitivity; continuity with every result before #788 | 317 | 1,198,119 | Nothing |
| `full` | Exploratory; its own correction family, never headline claims | 4,224 | 3,945,417 | Every strain edge that does not assert an absence, mostly CHEBI, API-assay, EC and GO relations; every negated key is dropped |

Moving the primary from `off` to `phenotype` changes every arm at its next rerun, not only
the new traits:

- **Multiple testing.** The full-trait screen's BY and BH families grow from m = 317 to
  m = 344. The ORA family grows with the trait count, as profiles × traits, and the
  calibrations are rerun at the new m.
- **Denominators.** Every profile's fraction sums 344 retained columns. The retained count
  mass summed over reported taxa rises by 11.4% (136,906). The 27 new traits add only 6.8%
  of that rise. Most of it comes from existing traits that strain records also assert:
  mesophilic, whose count mass goes from 27,509 to 96,358, adds 50.3% of the rise, rod
  shaped 15.9% and aerobic 10.4%.
- **Existing traits.** METPO phenotypes that strain records also assert gain coverage:
  coccus shaped goes from 721 to 815 reported taxa, facultatively aerobic from 916 to 920,
  sphere shaped from 305 to 307 and moderately halophilic from 350 to 351. Cellulose
  catabolic process and arsenate reduction have no strain edges and do not change.
- **Oxygen guilds.** Aerobic goes from 1,588 to 1,746 reported taxa, anaerobic from 1,707
  to 1,884 and obligately anaerobic from 319 to 414. The aerobe guild's count mass rises
  by 16.7% and the anaerobe guild's by 9.2%. Facultatively anaerobic is in neither guild.
- **New traits.** Facultatively anaerobic enters with 712 reported taxa. The six Table 1
  trait keys, and the Table 1 families built on them, stay the same.

A result names its scope. An `off` result is comparable with earlier runs, and a `full`
result is tested only within the `full` family.

### Common paired quantities and eligibility

Let `C(e,d,t)` be the count for trait `t` in experiment `e`, direction `d`, and
let `T(e,d)` be the full non-MediaDive total. Define two comparison types:

| Comparison type | Experiment contrast, increased minus decreased | Eligible complete experiment pair | Reported effect units |
| --- | --- | --- | --- |
| Individual trait fraction | `D_F(e,t) = C(e,inc,t)/T(e,inc) - C(e,dec,t)/T(e,dec)` | Both totals positive and the trait counts known. A zero trait count is allowed. | Percentage points: multiply the study-aggregated fraction contrast by 100 for display. |
| Trait ratio A/B | `D_R(e,A,B) = log2(A_inc/B_inc) - log2(A_dec/B_dec)` | All four named trait counts positive. Exclude the whole pair if any is zero; add no pseudocount. | Change in log2(A/B), with numerator and reference retained. |

For either type, compute `D(s) = mean(D(e))` over that study's eligible
experiments in the scope. Report the mean and median of `D(s)` across studies,
plus the proportion of non-tied studies with `D(s) > 0`. Score ties before
display conversion, using two separate absolute tolerances of 1e-12: one in
fraction units for `D_F`, one in log2 units for `D_R`. These are not one shared
rule — fraction contrasts are of order 1e-3 and log2 contrasts of order 1e-1, so
the same constant sits a different distance below each scale. Tied eligible studies stay in the mean/median and coverage
counts. If all eligible studies tie, the directional proportion is undefined;
if there are no eligible pairs, both the effect and directional proportion are
undefined. An interval must identify its parameter: the existing exact 95%
binomial interval concerns the non-tied directional proportion under its study
sampling assumptions, not the mean effect or a simultaneous 15-pair statement.
Unavailable effect intervals must be labeled as unavailable.

For a given comparison row, every implementation must use the **same eligible
experiment IDs and study weights**, not merely the same studies. An excluded
experiment removes both directions. Record eligible/excluded pairs and studies,
partially observed studies and ties. Fraction analyses exclude complete
experiment pairs whenever either side's retained total is zero via
`drop_zero_total_pairs`, conforming to the shared comparison basis above and
avoiding artificial ties (#228). Historical fraction analyses that mapped
zero-total profiles to zero (such as scenario S0 in `mtp.log_ratio_sensitivity`)
are retained for historical audit only and are distinct from this baseline.

The oxygen guild figure uses the same basis (#783). The repository owner chose
this on 2026-10-22, over keeping zero-total pairs in the figure and adding them
back to the integrated McNemar. The basis covers:

- the figure's condition-by-body-site cells;
- its cell-status table;
- its aerobe/anaerobe McNemar and rank correlations;
- the integrated evidence table's figure-cell inventory and McNemar row.

A zero-total pair's guild share is undefined, so it cannot enter a cell mean or
a study sign as a measured 0. Scored that way, it lowered cell means. It also
shaded cells as sampled without guild signal: on the committed corpus, 8 of the
20 such cells were backed only by profiles with no retained annotation. And it
gave the figure's McNemar a different population from the integrated table's
and the calibration's (#730), so one run reported two McNemar results. On the
committed corpus the exclusion removes 55 of 2,257 experiment pairs and leaves
922 of 932 studies.

The fraction and ratio populations also need separate names. First retain the
full fraction-eligible analysis. For each A/B comparison, additionally calculate
the fractions of A and B on exactly that ratio's eligible pairs. Comparing these
with the full fraction results exposes population changes; comparing the matched
fractions with A/B then exposes the changed estimand. Different ratio populations
must not be collapsed into a common ordering of individual traits. A comparison
across references requires a declared common-support scenario using the same
experiment pairs for every included reference, with its coverage reported.
[#183](https://github.com/CultureBotAI/microbial-trait-profiling/issues/183) tracks this planned
sensitivity analysis. Neither fraction contrasts nor ratios establish absolute
trait abundance changes.

### How the native results attach

Each common row needs its comparison type, trait ID or oriented pair, scope,
base-manifest identifier, eligible-pair identifier and weighting rule. Method
evidence additionally retains its native hypothesis, observation unit,
transformation/background, source version, p/q and correction family. Classify
an estimate as directly comparable only when its quantity, units, population and
weighting match. Different hypotheses can supply complementary evidence; an
unresolved population or identity mismatch remains explicitly non-comparable.

Jed's transformed experiment-level mean is not the shared raw fraction effect.
His source counts supply that common effect, while his native effect and test
retain their own units until the study-dependence work in
[#133](https://github.com/CultureBotAI/microbial-trait-profiling/issues/133) is completed. Noel's
study contrasts can match the shared quantities after eligibility is reconciled.
Their tests still target directional consistency rather than a nonzero mean.

Ilyas's profile enrichment records retain their study, experiment and direction.
Define `b(t) = sum(C(t)) / sum(T)` over all base profiles and, for positive `T`
and `b(t)`, `fold(e,d,t) = (C(e,d,t)/T(e,d)) / b(t)`. Derive these folds from
validated counts, including observed trait zeros; an omitted native output row
does not establish a zero fold or a nonsignificant test.

For a descriptive summary at a common scope, compute arithmetic mean fold
enrichment within study and direction on the common fraction-eligible pairs,
then average those study means equally. Keep increased and decreased summaries
separate, and record the fixed background and undefined-fold states. Do not use
the minimum profile p-value or an any-significant-profile indicator as a shared
effect. This descriptive aggregation creates no study-, site- or condition-level
p-value. Native profile p/q values remain linked evidence at their tested level;
they must not be inserted as q-values for a pooled trait or an A/B contrast.
[#177](https://github.com/CultureBotAI/microbial-trait-profiling/issues/177) tracks formal support
aggregation. Enrichment on both sides can coexist with a small paired contrast.

## Side by side

The table describes the inspected native analyses: for enrichment, the
production run, with the Table 1 native reruns' assertion urn noted where it
differs. The last two rows identify
their contribution to the shared basis above. Historical selection settings are
reported for provenance; they are not the thresholds chosen for harmonized
comparison. The 15-ratio extension is a separate result type within Noel's
framework, described under [Compositionality](#compositionality).

| | Composition | Enrichment | Directional consistency |
| --- | --- | --- | --- |
| Owner | Jed | Ilyas | Noel |
| Unit of observation | matched study/experiment pair | distinct reported taxon in the production run; trait-annotation assertion in the Table 1 native reruns | study |
| Transformation | fractions, x1000, log1p, per-trait z-score | none; each taxon counts once per trait it bears | fractions, study mean per direction, sign of the difference |
| Test | one-sample t against 0, two-sided | hypergeometric upper tail, one-tailed | exact binomial against 0.5, two-sided |
| Paired | yes, by construction | no; increased and decreased tested separately | yes, then reduced to one sign per study |
| Native multiple testing | Benjamini-Hochberg within each output table | Benjamini-Hochberg across every test in the inspected run | Benjamini-Yekutieli across all 317 retained traits, reported at the nominal and the effective family size (#255, #265); the six-trait Benjamini-Hochberg family is retained separately for the published Table 1 numbers |
| Historical selection/display rule | no significance cutoff; the figure takes the top ten by p | q at most 0.001 and fold at least 1.5 | no display cutoff; q reported per trait |
| Typical n per test | 2,257 pairs, identical for every trait | a mean of 9.05 taxa per profile against 5,829 background taxa in the production run; tens of thousands of assertions per profile in the Table 1 native reruns | 2 to 905 of the 922 nonzero studies (932 total) across the 317 candidates, 606 to 898 for the six of Table 1 (ties excluded) |
| Strength | Uses matched contrasts and their magnitudes; screens many traits | Locates profile-specific excess relative to a declared corpus background | Gives each study one directional vote; exact test without a normality assumption for magnitudes |
| Main limitations | Repeated experiments within studies are not independent; denominator and transformation affect the estimand | Nested taxa inherit the same annotations, so draws are not independent; coverage and background composition affect calls | Discards magnitude, excludes ties, and remains denominator-dependent; support varies by orders of magnitude across the 317 candidates, so a lopsided share on few studies is not a finding |
| Contribution to the shared comparison | Supply matched source counts and native transformed effects; reconcile the denominator, then report the common study-weighted fraction/ratio quantities | Supply profile counts, fixed-background folds and native profile evidence; provide separate descriptive increased/decreased study summaries | Supply study-level direction and coverage for the shared fraction quantities and all 15 explicitly referenced ratios |
| Meaning of disagreement | A mean contrast can be large despite inconsistent study signs | Background enrichment can occur on both sides without a paired shift | Consistent small shifts can coexist with a small mean; ratio directions also depend on the named reference |

The core paired analyses restrict to human samples and experiments carrying both
directions, at different points in the pipeline. The enrichment branch filters
host species at ingest; the other analyses filter downstream. The same retained
rows do not guarantee the same normalization: Ilyas and Noel remove MediaDive
before computing totals, whereas Jed's inspected paired-screen code excludes
MediaDive from trait selection but includes it in the denominator. Issue
[#174](https://github.com/CultureBotAI/microbial-trait-profiling/issues/174) tracks harmonization and
the numerical rerun; the historical results below must retain their original
denominator definitions until that work is complete.

## Thresholds and ranks follow the shared basis

Before choosing a score or cutoff, the owners must reconcile the source manifest,
trait mapping, eligible pairs, study weights, common descriptive estimates and
missingness counts. Record which native effects match the comparison target and
which provide background context. Disagreement from a different reference,
population or hypothesis is not a failed replication. Freeze these definitions
for the comparison run; the six traits' prior selection on overlapping data still
makes the exercise exploratory.

The initial displays contain all six individual traits and all 15 ratios in
their declared manifest order, with effect, uncertainty, coverage, native p/q,
family and status fields. Subsequent rank comparison must name the score, its
direction and the common candidates with defined scores. Compute rank agreement
before significance filtering, preserve ties and report excluded candidates.
Missing scores are not bottom ranks; a defined zero effect can remain a tie even
when its test is uninformative. Rankings of different native effect types measure
prioritization, not agreement on one effect scale. Any cross-reference ranking
also needs the common-support scenario defined above.

Thresholded displays must specify a valid tested level, correction family and
effect units before setting q/effect cutoffs. Keep full source families when
joining or subsetting results; do not recalculate correction only over the
overlap or the displayed winners. Any newly declared test family requires a
separate analysis record. Equal numerical q cutoffs do not equate the hypotheses
or repair invalid p-values
([SciPy's FDR assumptions](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html)).
No top-k list or agreement count is a combined significance test.

Rank and selection stability should resample whole studies, or shared-cohort
blocks where appropriate, retaining their paired annotation vectors and rerunning
the relevant filtering, background and analysis steps. This assesses stability;
formal error-rate claims for a threshold, rank-overlap statistic or combined
decision require the method-appropriate nulls and selection-aware calibration in
[#182](https://github.com/CultureBotAI/microbial-trait-profiling/issues/182).

## Why the p-values are not comparable

The three tests count different things, so their sample sizes differ by orders of
magnitude and their p-values cannot be ranked against one another.

The enrichment arm's production urn is built from taxa. Each of the 4,514
testable paired human profiles is a draw of its distinct reported taxa (mean 9.05
per profile) from the 5,829 BugSigDB-reported taxa of `taxa_trait_matrix.tsv`; a
taxon counts as trait-bearing if the matrix marks it so, including positivity
inherited from descendants at any rank. Over the 317 non-MediaDive traits that is
1,430,938 tests, which at Benjamini-Hochberg q at most 0.001 and fold at least 1.5
call 66,064 enrichments; the lower-tail run at fold at most 0.5 calls 3
depletions (`data/derived/trait_enrichment_all_results.tsv` and
`trait_depletion_all_results.tsv` headers). An earlier assertion-level run,
with a mean of 23,686.87 assertions per testable profile and 556,516 tests,
produced many p- and q-values of numerically zero; that sampling-unit problem is
#131, now closed. The taxon urn removes the repeated counting of one organism, but
nested taxa are still not independent draws: a genus, its family and its phylum
can all be reported and all inherit the same annotation. Calibration at the
production thresholds is tracked in #551. The background is not filtered the way
the profiles are. Profiles are paired human experiments, but `M` and `K` count
all 5,829 rows of `taxa_trait_matrix.tsv`: every taxon reported in any BugSigDB
signature, from any host. 645 of those taxa occur only in non-human samples, and
862 never appear in a paired human profile. "Corpus frequency" therefore means
the share of all reported taxa that bear the trait, and each profile's own taxa
are part of its background.

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
has been calibrated. [#176](https://github.com/CultureBotAI/microbial-trait-profiling/issues/176) and
[#177](https://github.com/CultureBotAI/microbial-trait-profiling/issues/177) track the grouping,
display-selection and profile-support questions in PR #173.

Noel owns the integrated evidence table in
[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178). Each method's result
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
or causation. [#179](https://github.com/CultureBotAI/microbial-trait-profiling/issues/179) tracks the
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
15 comparisons alongside the original six fraction directions as a matched-corpus
Table 1-trait analog sensitivity analysis (not a reproduction of the published
Table 1 baseline), with source fingerprints and the study populations retained
by each ratio.

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

### The matrix the counts in this note were measured on

Naming the dump does not identify the matrix. The matrix every arm consumes is a
**join** of the BugSigDB dump with a KG-Microbe build, and the KG side is pinned
by Drive file id and checksum but carries **no release tag**, so two runs can both
be "on pinned v1.3.1" and read different matrices. That is what produced the
disagreement this section was rewritten to settle (#222): an earlier revision of
this note reported 2,290 MediaDive columns, and the run below measures 2,333.

| | |
|---|---|
| matrix | `data/derived/study_trait_profiles_wide.tsv` |
| sha256 | `4c99aa8b4e386d42cc7938ab15c452ad4278b143b7779c4eb1635f71d386f21a` |
| header columns | 2,610 = 3 identifier + **2,290 MediaDive** + **317 retained** |
| profile rows | 5,629 (`count_input_16s_wms_v1` human 16S/WMS profiles; pre-#742 run had 7,329 rows before host filtering) |
| BugSigDB vintage | `2026-04-24_00:41_UTC` (v1.3.1) |
| KG-Microbe nodes | declared md5 `dd76b4656f9341db25a7c4afa96cd58e` (declared by `data/kg-microbe/provenance.json` sidecar) |
| KG-Microbe edges | declared md5 `555ecf5ba5d89288cdc44e368375374f` (declared by `data/kg-microbe/provenance.json` sidecar) |
| KG-Microbe release | `20260120` (pinned Drive build folder in `scripts/download_kg_microbe.py`) |

Filter cascade on that matrix: 5,629 profile rows read (0 unmatched in the
metadata join), **0 rows dropped as non-human** (already restricted to human 16S
and WMS profiles under `count_input_16s_wms_v1`), leaving 3,372 human
experiments after the host filter, of which 1,115 are dropped as unpaired,
leaving **2,257 experiments, 4,514 profiles and 932 studies**. Every recipe prints
this block and writes it to a `corpus_identity.json` beside its outputs, so a
number can be traced to the file it came from rather than to a vintage string.

Two further points, because matching one field does not establish the rest:

- The 317 **retained column names** are identical between this matrix and the
  committed `study_trait_profiles_wide.tsv.gz`. That does **not** make them the
  same corpus. The `.gz` holds 12,108 rows to this matrix's 7,329, the 7,329 are a
  strict subset of them, and of those shared rows **one differs**:
  `bsdb:25487798-1 inc` disagrees on 81 of the 317 retained columns, its retained
  total 12,587 here against 12,931 there. So the populations differ by design and
  the values differ in one place.
- Earlier revisions of the per-figure notes stated that the v1.3.1 wide matrix is
  "already human-only upstream, so the host filter drops 0 rows". That is not true
  of this matrix: the host filter drops 1,500 rows. Those notes have been
  regenerated from the run described here.

**Which corpus the project should standardise on is still open** — the committed
`.gz` or a regeneration — and belongs with #115. Until it is settled, each note
states the matrix it was measured on rather than implying there is only one.

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
statistical corrections, log-ratio extension and shared comparison specification
are tracked in
[PR #181](https://github.com/CultureBotAI/microbial-trait-profiling/pull/181) and require confirmation
by the method owners before manuscript use.
