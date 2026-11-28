# Sensitivity of the paired log-ratios to eligibility, reference and population

Produced by `src/mtp/log_ratio_sensitivity.py`, driven by
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
its own population: 241 to 682 of the 932 studies. Comparing a ratio with the
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
| S0 historical fractions | all 2,257 complete pairs; a zero-total profile scores 0 | BH over 6 |
| S1 fraction baseline | 2,214 pairs with positive totals on both sides | BH over 6 |
| S2 complete-case ratios | per ratio, pairs with all four counts positive | BY over 15 |
| S3 matched fractions | each ratio's numerator and reference fractions on exactly that ratio's S2 pairs | BY over 30 |
| S4 six-trait common support | 347 pairs (237 studies) with all six traits positive | BH over 6 fractions; BY over 15 ratios |
| S5 triple common support | per trait triple, pairs with all three positive; A/B, B/C, A/C | BY over 60 |
| S6 / S7 fecal strata | S2 within exact body site Feces / outside it | BY over 15 each |
| S9 / S10 / S11 label-set strata | S2 partitioned by body-site label set: any Feces label / other lower gut with no Feces label / no gut label | BY over 15 each |
| S8 influence | S2 with each study, then each stratum of 30+ eligible studies, removed in turn | BY over 15, other 14 p held fixed; crossing at q = 0.05 |

Registry IDs are `mt-supp-ratio-sensitivity-s0`, `-s1`, `-s2`, `-s3`, `-s5`,
`-s6`, `-s7`, `-s9`, `-s10` and `-s11` with the same prefix;
S4 has two distinct entries, `mt-supp-ratio-sensitivity-s4-fractions` and
`mt-supp-ratio-sensitivity-s4-ratios`. S8 uses
`mt-supp-ratio-influence-study` and `mt-supp-ratio-influence-stratum`.
Each scenario's q-values below use only its listed family and procedure.
Read the actual `family`, `family_size` and `correction` columns from that
execution's `sensitivity_results.tsv`; the sizes above belong to the described
historical six-trait grid. Bare p-values are nominal. There is no correction
across the collection of scenarios or omissions.

S0 reproduces the published six-trait directions and S2 reproduces the log-ratio
module's output exactly (study counts, partially observed studies, pairs, effects
and BY q), which the tests check.

## Who each ratio keeps

Exclusion is heavy and not random. Arsenate reduction is zero on at least one side
in 1,884 of 2,257 pairs, cellulose catabolic process in 1,396, sphere shaped in
1,134, moderately halophilic in 1,253, coccus shaped in 724 and facultatively
aerobic in 646. Between 123 and 240 studies per ratio are only partially observed,
contributing some experiments but not others. The pairs a ratio keeps come from
larger signatures: across ratios the eligible pairs have a median of 26,313 to
176,794 summed annotations and 17 to 27 reported taxa, against 1,949 to 7,365
annotations and 7 to 10 reported taxa for the excluded pairs. The complete-case
populations therefore over-represent richly annotated, taxon-rich signatures, and
a ratio result does not generalise to the excluded studies. Eligibility by exact
body site and condition is in `eligibility_by_stratum.tsv`; the per-profile
annotation-coverage sidecar of #184 is not yet available and is recorded as pending.

## What changes and what does not

**Excluding zero-total pairs changes nothing material (S0 against S1).** Dropping
43 pairs, which removes 9 studies entirely, leaves every one of the six fraction directions in place,
each with q <= 0.0112.

**Restriction to a ratio's pairs changes four things at once (S1 against S3).**
S1 tests 6 fractions under BH on 923 eligible studies, while S3 tests 30 matched
fractions under BY on 241 to 682 eligible studies (26.1% to 73.9% of the baseline). The
population, the eligible n, the family size and the correction all change together,
so neither a comparison of q status nor a comparison of intervals can attribute a
lost detection to any one of them. The four are reported separately below; the
row-by-row table is `population_sensitivity.tsv`.

*Effect estimates.* The restricted share differs from the baseline share by a
median of -0.040 and by at most 0.094 in absolute value.

*Interval overlap.* All 30 restricted 95% intervals overlap the baseline
interval for the same trait. Overlap does not show that the effect is unchanged: two
intervals can overlap while the estimates differ, and a wider interval overlaps more
easily. The sphere shaped interval on the
sphere shaped/arsenate reduction pairs, 0.466 (0.403 to 0.530) on 249 studies,
overlaps its baseline of 0.559 (0.525 to 0.593) on 923. Seven rows sit on the other side of 0.5 from
their baseline, but every one of the seven has 0.5 inside its restricted interval,
so none of them establishes a reversal.

*Eligible counts and correction.* Detection at q <= 0.05 is lost on 24 rows, in
two patterns:

- **Five rows lose only the correction (`lost at correction`).** They keep a
  nominal (uncorrected) p <= 0.05 and fail corrected q <= 0.05 once the family is
  30 rows under BY (so `lost at correction` compares the restricted run's nominal p
  with its corrected q). Moderately halophilic on the
  arsenate reduction/moderately halophilic pairs is the clearest: share 0.594
  (0.530 to 0.655) on 251 studies, *above* its own baseline of 0.571 (0.536 to
  0.605) on 923, with p = 0.0036 and BY q = 0.0617. Facultatively aerobic on the
  facultatively aerobic/moderately halophilic pairs is the same case: 0.553 (0.511
  to 0.594) on 574 studies, p = 0.0121, q = 0.182.
- **Nineteen rows lose nominal (uncorrected) p <= 0.05 before any correction (`lost before correction`).** On these rows the
  eligible n falls to between a quarter and three quarters of 923, and all
  nineteen intervals overlap their baseline.

How much of each lost detection comes from a shifted estimate, from the smaller
eligible n, and from the larger correction family has not been separated here; the
data above are compatible with any mixture of the three. The earlier phrasing -
that most fraction signals depend on the population - attributed the loss to the
population alone, which these scenarios cannot establish either, and is withdrawn.

One asymmetry survives this accounting. Cellulose catabolic process is the only
trait whose restricted share moves *further* from 0.5 on every population
(-0.024 to -0.081), and the only one that keeps q <= 0.05 on all five of its
populations (0.353 to 0.409). The other five traits together retain detection on 1
of their 25 rows (facultatively aerobic on the facultatively aerobic/coccus shaped
pairs, 0.574, q = 0.0030).

**The three ratios passing the S2 BY rule involve cellulose (S3 against S2).** The
ratios with BY q <= 0.05 are facultatively aerobic/cellulose (share up 0.593,
q = 0.0019), cellulose/sphere shaped (0.409, q = 0.0035) and cellulose/arsenate
reduction (0.368, q = 0.0019). On exactly the same pairs the cellulose fraction
alone shows the same direction and similar strength (0.409, 0.402 and 0.353), while
the facultatively aerobic, sphere shaped and arsenate reduction fractions do not
(0.528, 0.506 and 0.506, all q = 1). On these populations the ratio results are
therefore consistent with a relative decrease of cellulose-associated annotations
rather than with a shift in the other member of each pair.

**On one shared population the pattern is the same (S4, S5).** S4 runs two separate
families on that population - 6 fractions under BH and 15 ratios under BY - and the
counts below belong to one or the other, never to a pooled 21. On the 237 studies
with all six traits positive, cellulose has the only strong fraction shift (share
0.354, q = 5.2e-05) and moderately halophilic a weaker increase (0.586, q = 0.028),
2 of the 6 fractions; and 4 of the 15 ratios reach
q <= 0.05, those with cellulose against facultatively aerobic,
sphere shaped, arsenate reduction and moderately halophilic. Across the 20 triple
populations, all 12 rows with q <= 0.05 involve cellulose, and no ratio without
cellulose reaches q <= 0.05 in any population. Ratios without cellulose do not form
a stable ordering: their shares stay between about 0.44 and 0.58 and some change side
between populations.

**Nothing is detected outside the fecal stratum, which is not the same as the signal
being fecal (S6, S7, S9, S10, S11).**

First, the label handling, because it changes what S6 and S7 mean. BugSigDB packs
several body sites into one comma-joined metadata cell, and S6 and S7 match the
whole cell. S7 is therefore "outside the exact label Feces", not "outside the gut":
of its 839 pairs, 15 (from 10 studies) carry the label Feces alongside another site
(`Feces,Blood`, `Meconium,Feces`, `Colon,Feces` and so on), and a further 103 pairs
(46 studies) come from the rest of the lower gut - Rectum, Intestine, Colon,
Colorectum, Ileum and their mucosae, plus 11 pairs labelled only Meconium. 71 of the
2,257 pairs carry more than one label. S9, S10 and S11 therefore partition the same
pairs by label-set membership instead, and S6 and S7 are kept unchanged as declared.
The full accounting is in `sensitivity_status.json` under
`body_site_label_accounting`, and the per-ratio eligible counts per stratum are in
`eligibility_by_stratum.tsv` under `stratum_kind` `gut_label_set`.

Second, the counts that are actually tested. A stratum's `studies_total` is not its
eligible count, and the earlier text quoted the former. Per ratio:

| Stratum | Studies in stratum | Eligible studies per ratio | Rows tested | Rows at q <= 0.05 |
| --- | --- | --- | --- | --- |
| S6 exact label Feces | 588 | 144 to 407 | 15 | 3 |
| S7 outside exact label Feces | 368 | 99 to 290 | 15 | 0 |
| S9 any Feces label | 595 | 145 to 412 | 15 | 3 |
| S10 other gut, no Feces label | 46 | 12 to 36 | 5 | 0 |
| S11 no gut label | 319 | 85 to 250 | 15 | 0 |

Within the fecal stratum the same three ratios reach q <= 0.0027 under either
definition: facultatively aerobic/cellulose, cellulose/sphere shaped,
and cellulose/arsenate reduction. Cellulose/moderately halophilic has q = 0.0621
(S6) and 0.0642 (S9). Moving the 15 multi-label fecal pairs into the stratum
(S6 to S9) changes none of the three conclusions.

Outside it, nothing is detected, but the three strata establish different amounts.
S11, away from the gut altogether, is tested, with no detection at the declared
threshold, on 319 studies and 85 to 250 eligible per ratio: all 15 rows are tested
and the largest share is 0.562. Non-detection there does not establish that the
effect is absent or equivalent to zero. S10, the
rest of the lower gut, establishes almost nothing - 46 studies, 12 to 36 eligible
per ratio, so 10 of its 15 rows fall below the 30-study minimum and are reported
descriptively with no p or q; the 5 that are tested reach shares of 0.33 to 0.71 on
31 to 36 studies. Cellulose/arsenate reduction, the strongest fecal ratio, cannot be
tested there at all. And in S7 the cellulose ratios still point the same way without
reaching the threshold - cellulose/arsenate reduction sits at 0.418 on 99 eligible
studies, against 0.336 on 144 inside Feces.

So the defensible claim is that these ratios are detected in fecal signatures and
not detected outside them, that the non-gut stratum was tested with no detection at
the declared threshold, and that the rest of the lower gut is largely untested. Whether the fecal
restriction reflects biology, the fecal stratum simply holding two thirds of the
corpus, or heterogeneous pooling of many small sites cannot be separated here.

**No single study drives a ratio; the fecal stratum does (S8).** The tested family
is the 15 S2 ratios under BY with the other 14 p-values held at their full-data
values, and a ratio crosses when that held-family q passes 0.05.

Leave-one-study-out cannot detect influence at these population sizes, and is
reported only as a bound: removing any one of a ratio's 241 to 682 eligible studies
moves its share by at most 0.0027, in *either* direction, and no single removal
changes any ratio's q status. Every study within one sign group produces the
identical share change, so the previous `most_influential_studies` column ranked
studies by their identifier rather than by influence; it has been dropped.

What it takes to cross the rule is reported instead. The three significant ratios
need 17 to 23 of their studies removed, or 8 to 11 studies to reverse sign, before
they stop crossing. The twelve non-significant ratios need 10 to 61 removals, or 6
to 32 sign reversals, to start crossing - the lowest being cellulose catabolic
process/moderately halophilic at 10 removals and 6 flips. Every ratio can be made to
cross by some number of removals, so none is unreachable.

Removing a whole stratum is a different matter. Of the 59 ratio-by-stratum removals
large enough to test - a stratum must hold at least 30 of that ratio's eligible
studies, which leaves out 4,548 smaller ratio-by-stratum cells - exactly six
cross the rule: the three significant ratios with either the exact Feces body site
or any-Feces-label stratum removed. For exact Feces, facultatively aerobic/cellulose
changes from q 0.0019 to 1, cellulose/sphere shaped from 0.0035 to 1 and
cellulose/arsenate reduction from 0.0019 to 0.804. For any Feces label, the first
two again reach q = 1 and cellulose/arsenate reduction reaches 0.780. Removing Saliva,
colorectal cancer or Parkinson's disease leaves every ratio on the same side of the
rule. The three significant ratios therefore depend on the fecal stratum under either
label definition, not on any single study; `leave_one_stratum_out.tsv` holds every tested removal.

## Interpretation limits

The stable finding across these declared scenarios is narrow: in fecal signatures
with positive counts, cellulose-associated annotations are relatively lower on the
increased side than several other Table 1 annotations, and that relative shift
survives restriction to matched populations. It is conditional on positive
annotation counts, is measured on over-represented, richly annotated signatures,
and does not identify an absolute change in any trait. It holds in the fecal strata:
removing fecal pairs under either label definition is the only tested omission that takes any ratio
back across the rule, the non-gut stratum was tested with no detection at the
declared threshold, and the rest of the
lower gut is too small to test. Five of the six original
fraction directions stop being detected once restricted to a ratio-eligible
population, while their estimates shift (median -0.040, at most 0.094),
all 30 intervals overlap, the eligible n falls to
between 241 and 682, and the family grows from 6 under BH to 30 under BY. The
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
this corpus the smallest eligible count is 241, so it is inactive and the
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
