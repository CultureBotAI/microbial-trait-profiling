# Sensitivity of the paired log-ratios to eligibility, reference and population

Produced by `src/bsdbng/log_ratio_sensitivity.py`, driven by
`just table1-log-ratio-sensitivity`, which writes
`data/derived/table1_composition_sensitivity/`: every scenario row
(`sensitivity_results.tsv`), eligibility by ratio and by exact body site and
condition (`eligibility_by_ratio.tsv`, `eligibility_by_stratum.tsv`),
leave-one-study-out influence (`leave_one_study_out.tsv`), leave-one-stratum-out
influence (`leave_one_stratum_out.tsv`), the full-population
against restricted-population comparison of every matched fraction
(`population_sensitivity.tsv`), a status record for the integrated evidence table
(`sensitivity_status.json`) and a generated table of
every row (`sensitivity_report.md`). It implements
[#183](https://github.com/CultureBotAI/microbial-trait-profiling/issues/183) and the
population rules of the [shared comparison basis](three_arm_statistics.md#shared-comparison-basis).

## Why this analysis exists

The [paired log-ratio analysis](table1_paired_log_ratios.md) excludes an
experiment from ratio A/B when either trait is zero on either side, adds no
pseudocount, and averages the remaining experiments within study. That is a
deliberate complete-case rule, and it means each of the 15 ratios is measured on
its own population: 243 to 676 of the 959 studies. Comparing a ratio with the
original six fraction results therefore changes two things at once - who is
measured and what is measured. This analysis separates them.

## The declared scenario grid, and what can be said about when it was declared

The scenarios, their populations, their correction families, the 30-study minimum
for a test and the reporting threshold are in
[`table1_log_ratio_sensitivity_protocol.json`](table1_log_ratio_sensitivity_protocol.json).
Each scenario keeps its own family, no scenario is preferred, and q <= 0.05 is used
only to count rows. The same test is used throughout: study-level mean contrasts, an
exact two-sided binomial test of the non-tied studies shifting up, and an exact 95%
interval for that share (not for the mean effect).

**The claim that this was declared before the first run cannot be verified from this
repository, and is no longer asserted as a fact.** The protocol, the code, the
results and the tests were all committed together in `c70e342`, so nothing in the
history orders the protocol before the first run; the earlier
`declared_before_first_run: true` was self-attested and has been removed rather than
restated. What the available evidence does support is narrower: issue
[#183](https://github.com/CultureBotAI/microbial-trait-profiling/issues/183) names the
three axes - zero-exclusion rule, reference choice, eligible study population -
before this branch existed, so the shape of the grid predates the results; every
change made after a result was seen is listed in `protocol_amendments` with its date
and the reviewing issue; and a test enforces the scenario list against the code, so
an undeclared scenario cannot be added silently. S9 to S11 and the revised S8 were
added after review and are recorded as such.

| Scenario | Population | Family |
| --- | --- | --- |
| S0 historical fractions | all 2,335 complete pairs; a zero-total profile scores 0 | BH over 6 |
| S1 fraction baseline | 2,241 pairs with positive totals on both sides | BH over 6 |
| S2 complete-case ratios | per ratio, pairs with all four counts positive | BY over 15 |
| S3 matched fractions | each ratio's numerator and reference fractions on exactly that ratio's S2 pairs | BY over 30 |
| S4 six-trait common support | 354 pairs (239 studies) with all six traits positive | BH over 6 fractions; BY over 15 ratios |
| S5 triple common support | per trait triple, pairs with all three positive; A/B, B/C, A/C | BY over 60 |
| S6 / S7 fecal strata | S2 within exact body site Feces / outside it | BY over 15 each |
| S9 / S10 / S11 label-set strata | S2 partitioned by body-site label set: any Feces label / other lower gut with no Feces label / no gut label | BY over 15 each |
| S8 influence | S2 with each study, then each stratum of 30+ eligible studies, removed in turn | BY over 15, other 14 p held fixed; crossing at q = 0.05 |

S0 reproduces the published six-trait directions and S2 reproduces the log-ratio
module's output exactly (study counts, partially observed studies, pairs, effects
and BY q), which the tests check.

## Who each ratio keeps

Exclusion is heavy and not random. Arsenate reduction is zero on at least one side
in 1,950 of 2,335 pairs, cellulose catabolic process in 1,464, sphere shaped in
1,195, moderately halophilic in 1,308, coccus shaped in 809 and facultatively
aerobic in 769. Between 124 and 247 studies per ratio are only partially observed,
contributing some experiments but not others. The pairs a ratio keeps come from
larger signatures: across ratios the eligible pairs have a median of 27,857 to
158,832 summed annotations and 17 to 27 reported taxa, against 1,745 to 6,536
annotations and 7 to 10 reported taxa for the excluded pairs. The complete-case
populations therefore over-represent richly annotated, taxon-rich signatures, and
a ratio result does not generalise to the excluded studies. Eligibility by exact
body site and condition is in `eligibility_by_stratum.tsv`; the per-profile
annotation-coverage sidecar of #184 is not yet available and is recorded as pending.

## What changes and what does not

**Excluding zero-total pairs changes nothing material (S0 against S1).** Dropping
94 pairs, which removes 20 studies entirely, leaves every one of the six fraction directions in place,
each with q <= 0.008.

**Restriction to a ratio's pairs changes four things at once (S1 against S3).**
S1 tests 6 fractions under BH on 939 eligible studies, while S3 tests 30 matched
fractions under BY on 243 to 676 eligible studies (26% to 72% of the baseline). The
population, the eligible n, the family size and the correction all change together,
so neither a comparison of q status nor a comparison of intervals can attribute a
lost detection to any one of them. The four are reported separately below; the
row-by-row table is `population_sensitivity.tsv`.

*Effect estimates.* The restricted share differs from the baseline share by a
median of -0.037 and by at most 0.107 in absolute value.

*Interval overlap.* 29 of the 30 restricted 95% intervals overlap the baseline
interval for the same trait. Overlap does not show that the effect is unchanged: two
intervals can overlap while the estimates differ, and a wider interval overlaps more
easily. One row's interval has moved clear of its baseline: sphere shaped on the
sphere shaped/arsenate reduction pairs, 0.447 (0.384 to 0.510) on 253 studies
against 0.554 (0.520 to 0.588) on 939. Seven rows sit on the other side of 0.5 from
their baseline, but every one of the seven has 0.5 inside its restricted interval,
so none of them establishes a reversal.

*Eligible counts and correction.* Detection at q <= 0.05 is lost on 24 rows, in
two patterns:

- **Five rows lose only the correction.** They keep a raw p <= 0.05 and fail
  q <= 0.05 once the family is 30 rows under BY. Moderately halophilic on the
  arsenate reduction/moderately halophilic pairs is the clearest: share 0.588
  (0.525 to 0.649) on 255 studies, *above* its own baseline of 0.572 (0.538 to
  0.605) on 939, with p = 0.0058 and BY q = 0.0985. Facultatively aerobic on the
  facultatively aerobic/moderately halophilic pairs is the same case: 0.554 (0.513
  to 0.595) on 582 studies, p = 0.0101, q = 0.150.
- **Nineteen rows lose p <= 0.05 before any correction.** On these rows the
  eligible n falls to between a quarter and three quarters of 939, and eighteen of
  the nineteen intervals overlap their baseline.

How much of each lost detection comes from a shifted estimate, from the smaller
eligible n, and from the larger correction family has not been separated here; the
data above are compatible with any mixture of the three. The earlier phrasing -
that most fraction signals depend on the population - attributed the loss to the
population alone, which these scenarios cannot establish either, and is withdrawn.

One asymmetry survives this accounting. Cellulose catabolic process is the only
trait whose restricted share moves *further* from 0.5 on every population
(-0.027 to -0.076), and the only one that keeps q <= 0.05 on all five of its
populations (0.354 to 0.403). The other five traits together retain detection on 1
of their 25 rows (facultatively aerobic on the facultatively aerobic/coccus shaped
pairs, 0.567, q = 0.0104).

**The three significant ratios are carried by cellulose (S3 against S2).** The
ratios with BY q <= 0.05 are facultatively aerobic/cellulose (share up 0.595,
q = 0.0011), cellulose/sphere shaped (0.407, q = 0.0026) and cellulose/arsenate
reduction (0.365, q = 0.0011). On exactly the same pairs the cellulose fraction
alone shows the same direction and similar strength (0.403, 0.401 and 0.354), while
the facultatively aerobic, sphere shaped and arsenate reduction fractions do not
(0.531, 0.499 and 0.510, all q = 1). On these populations the ratio results are
therefore consistent with a relative decrease of cellulose-associated annotations
rather than with a shift in the other member of each pair.

**On one shared population the pattern is the same (S4, S5).** S4 runs two separate
families on that population - 6 fractions under BH and 15 ratios under BY - and the
counts below belong to one or the other, never to a pooled 21. On the 239 studies
with all six traits positive, cellulose has the only strong fraction shift (share
0.356, q = 5.7e-05) and moderately halophilic a weaker increase (0.590, q = 0.019),
2 of the 6 fractions; and 4 of the 15 ratios reach
q <= 0.05, those with cellulose against facultatively aerobic,
sphere shaped, arsenate reduction and moderately halophilic. Across the 20 triple
populations, all 13 rows with q <= 0.05 involve cellulose, and no ratio without
cellulose reaches q <= 0.05 in any population. Ratios without cellulose do not form
a stable ordering: their shares stay between 0.43 and 0.60 and some change side
between populations.

**Nothing is detected outside the fecal stratum, which is not the same as the signal
being fecal (S6, S7, S9, S10, S11).**

First, the label handling, because it changes what S6 and S7 mean. BugSigDB packs
several body sites into one comma-joined metadata cell, and S6 and S7 match the
whole cell. S7 is therefore "outside the exact label Feces", not "outside the gut":
of its 875 pairs, 15 (from 10 studies) carry the label Feces alongside another site
(`Feces,Blood`, `Meconium,Feces`, `Colon,Feces` and so on), and a further 105 pairs
(47 studies) come from the rest of the lower gut - Rectum, Intestine, Colon,
Colorectum, Ileum and their mucosae, plus 11 pairs labelled only Meconium. 73 of the
2,335 pairs carry more than one label. S9, S10 and S11 therefore partition the same
pairs by label-set membership instead, and S6 and S7 are kept unchanged as declared.
The full accounting is in `sensitivity_status.json` under
`body_site_label_accounting`, and the per-ratio eligible counts per stratum are in
`eligibility_by_stratum.tsv` under `stratum_kind` `gut_label_set`.

Second, the counts that are actually tested. A stratum's `studies_total` is not its
eligible count, and the earlier text quoted the former. Per ratio:

| Stratum | Studies in stratum | Eligible studies per ratio | Rows tested | Rows at q <= 0.05 |
| --- | --- | --- | --- | --- |
| S6 exact label Feces | 602 | 146 to 403 | 15 | 4 |
| S7 outside exact label Feces | 381 | 99 to 287 | 15 | 0 |
| S9 any Feces label | 609 | 147 to 407 | 15 | 4 |
| S10 other gut, no Feces label | 47 | 12 to 35 | 5 | 0 |
| S11 no gut label | 331 | 85 to 249 | 15 | 0 |

Within the fecal stratum the same four ratios reach q <= 0.021 under either
definition: facultatively aerobic/cellulose, cellulose/sphere shaped,
cellulose/arsenate reduction and cellulose/moderately halophilic. Moving the 15
multi-label fecal pairs into the stratum (S6 to S9) changes none of the four
conclusions.

Outside it, nothing is detected, but the three strata establish different amounts.
S11, away from the gut altogether, is tested, with no detection at the declared
threshold, on 331 studies and 85 to 249 eligible per ratio: all 15 rows are tested
and the largest share is 0.554. Non-detection there does not establish that the
effect is absent or equivalent to zero. S10, the
rest of the lower gut, establishes almost nothing - 47 studies, 12 to 35 eligible
per ratio, so 10 of its 15 rows fall below the 30-study minimum and are reported
descriptively with no p or q; the 5 that are tested reach shares of 0.34 to 0.70 on
30 to 35 studies. Cellulose/arsenate reduction, the strongest fecal ratio, cannot be
tested there at all. And in S7 the cellulose ratios still point the same way without
reaching the threshold - cellulose/arsenate reduction sits at 0.418 on 99 eligible
studies, against 0.331 on 146 inside Feces.

So the defensible claim is that these ratios are detected in fecal signatures and
not detected outside them, that the non-gut stratum was tested with no detection at
the declared threshold, and that the rest of the lower gut is largely untested. Whether the fecal
restriction reflects biology, the fecal stratum simply holding two thirds of the
corpus, or heterogeneous pooling of many small sites cannot be separated here.

**No single study drives a ratio; the fecal stratum does (S8).** The tested family
is the 15 S2 ratios under BY with the other 14 p-values held at their full-data
values, and a ratio crosses when that held-family q passes 0.05.

Leave-one-study-out cannot detect influence at these population sizes, and is
reported only as a bound: removing any one of a ratio's 243 to 676 eligible studies
moves its share by at most 0.0026, in *either* direction, and no single removal
changes any ratio's q status. Every study within one sign group produces the
identical share change, so the previous `most_influential_studies` column ranked
studies by their identifier rather than by influence; it has been dropped.

What it takes to cross the rule is reported instead. The three significant ratios
need 19 to 26 of their studies removed, or 9 to 12 studies to reverse sign, before
they stop crossing. The twelve non-significant ratios need 7 to 62 removals, or 4
to 33 sign reversals, to start crossing - the lowest being cellulose catabolic
process/moderately halophilic at 7 removals and 4 flips. Every ratio can be made to
cross by some number of removals, so none is unreachable.

Removing a whole stratum is a different matter. Of the 24 ratio-by-stratum removals
large enough to test - a stratum must hold at least 30 of that ratio's eligible
studies, which leaves out 4,605 smaller ratio-by-stratum cells - exactly three
cross the rule, and all three are a significant ratio with the Feces body site
removed: facultatively aerobic/cellulose (q 0.0011 to 1), cellulose/sphere shaped
(0.0026 to 1) and cellulose/arsenate reduction (0.0011 to 0.748). Removing Saliva,
colorectal cancer or Parkinson's disease leaves every ratio on the same side of the
rule. The three significant ratios therefore depend on one stratum, not on any
single study; `leave_one_stratum_out.tsv` holds every tested removal.

## Interpretation limits

The stable finding across these declared scenarios is narrow: in fecal signatures
with positive counts, cellulose-associated annotations are relatively lower on the
increased side than several other Table 1 annotations, and that relative shift
survives restriction to matched populations. It is conditional on positive
annotation counts, is measured on over-represented, richly annotated signatures,
and does not identify an absolute change in any trait. It holds in one stratum:
removing the fecal pairs is the only tested stratum removal that takes any ratio
back across the rule, the non-gut stratum was tested with no detection at the
declared threshold, and the rest of the
lower gut is too small to test. Five of the six original
fraction directions stop being detected once restricted to a ratio-eligible
population, while their estimates shift (median -0.037, at most 0.107, one interval
moving clear of its baseline), 29 of 30 intervals overlap, the eligible n falls to
between 243 and 676, and the family grows from 6 under BH to 30 under BY. The
contribution of each of these to the lost detections has not been separated. These scenarios therefore neither confirm nor refute those five
directions on the restricted populations; they establish only that the restricted
data are too few, under a correction family five times larger, to detect a shift of
that size. Agreement between scenarios is robustness to these analytical choices on
the same annotation data, not independent replication. The six traits were
selected on overlapping data, the signature directions have not been verified as
disease versus control (#188), none of these procedures has been calibrated (#182),
and annotation coverage is not yet available per profile (#184).

## Reproducing

`just table1-log-ratio-sensitivity` on the pinned BugSigDBExports v1.3.1 inputs
(export timestamp 2026-04-24_00:41_UTC). The status JSON records the code revision,
the SHA-256 of the profile matrix, experiment table, taxa table, analysis source
and protocol, the corpus filter cascade, the body-site label accounting, the
influence decision rule, and the test counts. All statistics use scipy
(`binomtest`, `false_discovery_control`).

Counts are reported **by family**, in `family_status`, and that is the field the
integrated evidence table of #178 should read. `scenario_status` is kept for
continuity but merges families where one scenario tests two: S4 runs 6 fractions
under BH and 15 ratios under BY, and its scenario-level "tested 21" is not a usable
count for either.

S2's declared 30-study minimum is a guard that `compute_log_ratio_shifts` does not
apply, so S2 reproduces that module exactly only while every ratio clears it. On
this corpus the smallest eligible count is 243, so it is inactive and the
reproduction is exact; the guard is kept because the strata of S4 to S11 do fall
below it, and a test now asserts the divergence on a fixture sized to straddle 30
rather than leaving it hidden behind `min_studies=1`.

## Corrections to the commit record

The message on `c70e342` says the scenario grid evaluates a pseudocount. It does
not: no scenario adds a pseudocount, no pseudocount exists anywhere in the code, and
the protocol's `not_done` list has always said so. The message is wrong and is left
uncorrected rather than rewriting shared history.

The protocol, code and this description were drafted with the assistance of a
coding agent (Claude Code) and reviewed by the author.
