# Do these directional tests hold their error rates?

Produced by `src/bsdbng/calibration.py`, driven by `just calibration-report`, which
writes `data/derived/calibration/`: every scenario's per-test and family results
(`calibration_results.tsv`), one row per procedure and replicate
(`calibration_replicates.tsv`), the effective injection sizes
(`calibration_injection.tsv`), the declared criteria and their verdicts
(`calibration_criteria.tsv`), the shared flip assignments of #263
(`calibration_flip_assignments.tsv` and its manifest) and a status record for the
integrated evidence table (`calibration_status.json`). It implements
[#182](https://github.com/CultureBotAI/microbial-trait-profiling/issues/182) for the
procedures this repository owns.

## What is being asked

Agreement between arms, a correct binomial and a Benjamini-Hochberg step do not show
that a procedure rejects a true null at the rate it claims. This measures that
directly: make the null true, run the whole procedure, and count rejections.

| Procedure | Null | Unit | Family | Claimed rate |
| --- | --- | --- | --- | --- |
| `fraction_sign_test` | a study is equally likely to shift a trait's annotation share up or down | study | BH over the six Table 1 traits | per test 0.05; FDR 0.05 |
| `log_ratio_sign_test` | a study's paired log2(A/B) change is equally likely to be positive or negative | study | BY over the 15 pairs | per test 0.05; FDR 0.05 |
| `guild_mcnemar` | among studies whose guild signs disagree, either direction is equally likely | study | one test | 0.05 |
| `screened_fraction_sign_test` | as the first, but the three traits furthest from an even split are chosen inside the replicate and corrected over those three | study | BH over the three selected | none claimed |

The last is not a published procedure. It stands in for the step that produced the
six Table 1 traits - a selection made on this same corpus - to illustrate what
selecting and then testing the selection costs. It selects on the test's own
statistic, three of six; the published selection used other statistics over hundreds
of traits. Neither dominates the other, so the measured inflation is an illustration
of the cost, not a bound on it in either direction (#363). It is also the
**negative control**: the protocol declares in advance that it is expected to fail
the acceptance criteria, because a criterion this procedure passes is too weak to be
evidence about any procedure.

## How the null is made

Not by simulating counts, but by flipping directions on the real corpus: for a
sampled half of the units, the increased and decreased labels of every profile in
that unit are exchanged. Under each null above the two sides of a pair are
exchangeable, so relabelling leaves the distribution unchanged. Because the counts
are untouched, everything the issue asks to be retained is retained by construction:
whole annotation vectors and their correlations, unequal experiment counts per
study, zeros and uneven annotation coverage, and unequal body-site and condition
representation.

The flip unit matters, so both are reported. Flipping **whole studies** preserves the
dependence between a study's experiments and matches the unit these tests treat as
independent. Flipping **whole experiments** breaks that dependence and is the weaker
null.

**What the flip cannot test (#360).** A flip negates a unit's contrast without
changing its magnitude, so every scored count is fixed and each sign becomes a fair
coin. An exact sign test then sees exactly Binomial(n, 1/2), which is its own null,
and its per-test rejection rate is its exact size: known in closed form, at most
0.05, and independent of the data. For the three published procedures the simulated
per-test rate therefore only re-estimates that constant, under either flip unit, and
`per_test_level` is reported for them as a **sanity check** on the implementation
beside the analytic size, not as evidence that the method holds its level. It is
excluded from their overall status unless it fails. Neither flip can detect
dependence *between* studies - shared cohorts, reused datasets, overlapping
populations - because both treat studies as independent; that is listed under
`not_covered` in the protocol and is not addressed here. What the flip does test is
everything the construction leaves open: the family rates, where the dependence
between traits matters; the screened selection, which breaks the construction; and
power and the false-discovery proportion under an injected effect.

The generator is `numpy.random.default_rng([seed, sha256_64(stream), replicate])`.
Any one replicate is therefore reproducible on its own, which is what lets the same
assignments be handed to another arm (see **The shared null** below).

Power is measured by raising one trait's counts (coccus shaped) on the increased side
of half the studies by a declared multiplier. The added counts come **from the
complement**: a profile's total sums all retained trait columns (317 under `off`, 344
under `phenotype`), the eight the procedures read are disjoint subsets of them, and
the raise moves annotations in from the other 309 (336 under `phenotype`), which no
procedure reads. The total is held fixed, so every other tracked trait's
share is unchanged *exactly* - which is the only thing that makes the "untouched"
truth labels, and therefore the false-discovery proportion, mean anything. An earlier
version added the injected counts to the total, which lowered every other trait's
increased-side share and manufactured real effects on the traits it then counted as
false discoveries; `calibration_protocol.json` records that amendment.

## The independent unit, and what it costs

Labels tested inside one replicate share that replicate's resampled corpus. Pooling
them as independent Bernoulli trials understates the uncertainty, so
`calibration_replicates.tsv` keeps the cluster structure and the criteria read it:
the pooled per-test rate is a ratio estimator over replicate clusters, which reduces
to the binomial variance when a replicate contributes one test and inflates it by the
measured design effect when it contributes several correlated ones. The measured
design effects are 1.00 for the single McNemar, 1.16 for the six-trait family, 1.21
for the screened three and 1.62 for the 15 overlapping ratios - the ratios are the
most correlated, as expected.

Per-label rates keep their exact two-sided binomial interval, which is valid because
replicates are independent. Family rates are one Bernoulli per replicate, so their
bound is exact too.

## Declared criteria and their verdicts

The criteria are **amended** from `calibration-v1`, and the amendment is recorded in
[`calibration_protocol.json`](calibration_protocol.json) rather than applied quietly.
v1 passed a label when the lower bound of its exact two-sided 95% interval was at
most 0.05, which is not evidence of control: at 200 replicates it accepts 16/200, a
true rate of 0.08 clears it 57% of the time and 0.06 clears it 91%, and the screened
procedure passed it on 10 of 12 label x scenario rows. v2 pools per procedure and
passes only when the **one-sided 95% upper bound** is at most a declared tolerance of
**0.065** - the nominal 0.05 plus 0.015 of tolerated inflation. The burden of proof is
on the procedure, because an interval from a procedure whose true rate is 0.08 usually
covers 0.05 too.

The budget follows from the tolerance rather than the other way round. `guild_mcnemar`
contributes one test per replicate, so it is the least-powered procedure, and the
family criterion is one Bernoulli per replicate on the same count. 1400 null
replicates per scenario give 2800 pooled trials. The operating characteristics are
computed by **exact enumeration of the rule the verdicts use**, not by a normal
approximation; the two rules - the per-test bound and the family exact bound - give
the same numbers at this design:

| True per-test rate | What the criterion does at 2800 trials |
| --- | --- |
| 0.05 | passes with probability 0.960 |
| 0.06 | fails with probability 0.722 |
| 0.07 | fails with probability 0.997 |
| 0.08 | fails with probability 1.000 |

The budget was 1200 replicates until 2026-10-09. Enumerated exactly, both rules then
passed a true 0.05 with probability 0.925, below the protocol's requirement of 0.95;
the 0.958 reported earlier came from a normal model of a rule no longer in use
(#361). 1400 is the first round budget that meets the requirement; 1325 is the
smallest, at 0.953. The change is recorded as an amendment in the protocol.

A true 0.06 lies **inside** the declared tolerance, so where it is measured precisely
it passes by declaration, not by weakness - at the 16,800 trials the six-trait family
contributes, a true 0.06 passes 78% of the time with that family's measured design
effect of about 1.18 (85% if its six tests were independent). Treating 0.06 as
material inflation would mean declaring a tolerance of 0.055 and a budget several
times larger. That budget was not taken. The criterion is sized to separate rates at
or below 0.065 from rates at or above 0.07,
and every criterion row in `calibration_criteria.tsv` carries these four
probabilities for the design that produced it, so no verdict is reported without the
probability it would have been reached.

The pooled per-test bound is the larger of the replicate-clustered normal bound and
an exact Clopper-Pearson bound at the effective sample size (trials divided by the
measured design effect). The normal bound alone collapses to the point estimate when
no replicate rejects, because the empirical variance is then zero: a five-replicate
quick run reported `guild_mcnemar` at 0/10 with an upper bound of 0 and exported it
as a pass (#240 review, 2026-10-03). With zero rejections the design effect cannot
be measured, so each replicate counts as one trial, the fully correlated worst case;
0/10 then gives 0.259, the same bound the family criterion reports for those ten
replicates. At the declared budget the exact bound sits within 0.001 of the normal
one.

A run below the declared replicate budget cannot certify calibration. A criterion it
would pass, or would fail only because its interval is still wide, is reported as
`insufficient_budget`; it still fails a criterion whose point estimate already
exceeds the tolerance. `calibration_status.json` exports every criterion's verdict
per procedure beside an overall status: `fail` if any judged criterion failed,
otherwise `pass` only if every judged criterion passed, otherwise `not_certified`.

At 1400 null and 400 power replicates, seed 182, on the pinned v1.3.1 corpus:

| Criterion | Procedure | Verdict | Evidence |
| --- | --- | --- | --- |
| Pooled per-test upper bound at most 0.065 | fraction_sign_test | pass (**sanity check**) | 0.0465 (781/16,800), upper bound 0.0495; analytic exact size 0.0466 |
| | log_ratio_sign_test | pass (**sanity check**) | 0.0439 (1,843/42,000), upper bound 0.0460; analytic exact size 0.0432 |
| | guild_mcnemar | pass (**sanity check**) | 0.0425 (119/2,800), upper bound 0.0493; analytic exact size 0.0480 |
| | screened_fraction_sign_test | **fail** (expected) | 0.0927 (779/8,400), upper bound 0.0987 |
| Family upper bound at most 0.065 | fraction_sign_test | **pass** | 0.0454 (127/2,800), exact upper bound 0.0524 |
| | log_ratio_sign_test | **pass** | 0.0136 (38/2,800), exact upper bound 0.0178 |
| | guild_mcnemar | **pass** | 0.0425 (119/2,800), exact upper bound 0.0493 |
| | screened_fraction_sign_test | **fail** (expected) | 0.0861 (241/2,800), exact upper bound 0.0953 |
| Power reported with an interval no wider than 0.12 | fraction, log-ratio | **pass** | widest interval 0.069 and 0.099 |
| Selection inflation measured, not judged | screened | **reported** | per-test rate 0.070 to 0.122 |

The analytic exact size is the mean over each procedure's null labels at their median
scored count; the simulated rate agrees with it, as it must under this null (#360).
For `guild_mcnemar` the family and per-test rows are the same single test, so its
family pass is subject to the same construction.

The screened procedure fails both judged criteria, which is the declared expectation
and the point of including it: the criterion has teeth on real data, not only on a
fixture. Had it passed, `expectation_met` would say so in the output.

The three published procedures pass the family criterion and report power within
the declared width. That is a statement about these four procedures under these two
nulls on this corpus at this tolerance, and nothing more; in particular it says
nothing about dependence between studies.

## What the numbers say

**Under the complete null the per-test rates are the tests' exact sizes, as
constructed.** Mean per-test rejection at p <= 0.05, over 1400 replicates per scenario:

| Procedure | Study flip | Experiment flip | Highest single test |
| --- | --- | --- | --- |
| fraction_sign_test (6 tests) | 0.0462 | 0.0468 | 0.0585 (facultatively aerobic, experiment flip) |
| log_ratio_sign_test (15 tests) | 0.0427 | 0.0450 | 0.0550 (facultatively aerobic/sphere shaped, study flip) |
| guild_mcnemar (1 test) | 0.0421 | 0.0429 | 0.0429 |
| screened_fraction_sign_test (3 of 6) | 0.0924 | 0.0931 | 0.1222 (facultatively aerobic, experiment flip) |

The first three rows estimate known constants and are reported as a check on the
code. The screened row is the informative one: selection breaks the construction.

At the family level, the share of replicates with any q <= 0.05 is 0.040 and 0.051
for the six-trait BH family, 0.015 and 0.012 for the 15-pair BY family, 0.042 and
0.043 for the single McNemar, and 0.088 and 0.084 for the screened three.
Benjamini-Yekutieli is visibly conservative here, as its arbitrary-dependence
guarantee implies: it spends power to be safe under the dependence between
overlapping ratios, and the measured design effect of 1.62 is that dependence.

The tests are discrete, so the attainable p-values are a finite set, and the smallest
p *observed* is a different quantity from the resolution the test can *attain*. Both
are reported per label. Under the study flip the smallest observed null p is 2.1e-04
for the fraction family, 4.3e-05 for the ratios and 2.7e-04 for the McNemar; at the
median scored unit count the two attainable two-sided p-values straddling 0.05 are
0.0460 and 0.0540 for the fraction family, so the decision at 0.05 is about 0.008
coarse. Uniformity of p is neither expected nor required.

**Power is modest for small effects, and the nominal lift overstates the real one.**
`round` leaves a small count where it is, so 16.6% of eligible profiles are unchanged
at a nominal 1.1x and 7.0% at 1.3x; at 4x the complement caps about 30 profiles per
replicate. The effective multipliers are therefore 1.086, 1.284, 2.000 and 3.979.
Rejection at q <= 0.05 for the injected trait, 400 replicates, with its exact
interval:

| Nominal lift | Effective lift | fraction_sign_test | log_ratio_sign_test (mean of the 5 pairs naming the trait) |
| --- | --- | --- | --- |
| 1.1x | 1.086 | 0.018 [0.007, 0.036] | 0.003 (best pair 0.005) |
| 1.3x | 1.284 | 0.088 [0.062, 0.120] | 0.023 (best pair 0.073) |
| 2x | 2.000 | 0.868 [0.830, 0.899] | 0.441 (best pair 0.890) |
| 4x | 3.979 | 1.000 [0.991, 1.000] | 0.870 (best pair 1.000) |

A 28 per cent effective lift in one trait's annotation share across half of 959
studies is detected about one time in eleven. **A non-significant result from these
procedures is therefore weak evidence of no effect**, which matters wherever the
analyses report that something does not reach a threshold - the non-fecal strata, the
eleven non-cellulose ratios, and every "no cell survives correction" statement.

**The false-discovery proportion is now measured against valid truth labels.** With
the total held fixed, the five untouched traits reject at 0.051, 0.046, 0.052 and
0.051 across the four lifts - at the nominal level, where under v1's construction the
untouched traits rejected at a mean of 0.095 at 4x because the construction had given
them a real effect. The mean false-discovery proportion of the fraction family is
0.051, 0.038, 0.043 and 0.040 across the four lifts, and below 0.01 for the ratio
family. `guild_mcnemar` is untouched by the injection by construction, and its
false-discovery proportion of 0.038 to 0.055 is simply its null rate, which confirms
that the injection no longer leaks into it.

**Selecting traits on the data and testing the selection inflates the error rate.**
The screened procedure rejects a true null 0.093 of the time per test against a
nominal 0.05, its family rate is 0.086 against a claimed 0.05, and its
false-discovery proportion under an injected signal reaches 0.104. The six Table 1
traits were selected on this same corpus from hundreds of traits - by the composition
arm's paired contrast at q <= 0.001, according to the selection record on
`clustered-heatmap-generation` (`six-shared-traits/README.md`); `trait_screen.py` on
`main` describes them as the intersection of the composition and enrichment arms'
significant sets, and which account is right is open on #363. Either way the
published q-values for those six are not post-selection valid, and the measured
inflation here illustrates the cost without bounding it. This is the concrete statistical
basis for treating the six-trait results as **exploratory** rather than confirmatory,
and it cannot be repaired by correcting within the six.

## The shared null

`just calibration-flip-assignments` writes `calibration_flip_assignments.tsv`: one
row per (replicate, study) carrying the seed, the stream, the unit and a 0/1 flipped
flag, with a manifest recording the generator, the draw order, the input digests and
the table's own sha256. The stream is named for the scenario it comes from, so
replicate *i* of the artefact is replicate *i* of `complete_null_study_flip`: an arm
that reruns its own pipeline on these assignments is measured on the same null as
this arm's first 100 replicates, not merely a similar one. 100 replicates matches the
"roughly 100 permutations" of
[#263](https://github.com/CultureBotAI/microbial-trait-profiling/issues/263), where
Ilyas leads the battery; the format is a proposal, not a fait accompli.

## What this does not cover

- The enrichment arm's hypergeometric procedure. Its production urn now draws
  distinct taxa rather than trait-annotation assertions (#131), but nested taxa
  inherit the same annotations and are not independent draws. That is a model
  question a resampling study of this shape cannot settle; its calibration at the
  production thresholds is tracked in #551, and its adapter belongs to its owner.
- The composition arm's native transformed test, for the same reason. Because those
  two arms are not covered, this work is part of #182 and does not close it.
- The aggregated profile-support rule of #177, whose code is not on main.
- Calibration is not validity of interpretation. A procedure that holds its error
  rate still measures annotation shares, and none of this addresses orientation
  (#188), coverage (#184) or biological meaning (#179).

There is no untouched validation set: every result in this repository is computed on
one corpus, and the six traits were chosen on it. A retrospective split would not be
independent validation, and none is claimed.

## Reproducing

`just calibration-report` at the declared budget of 1400 null and 400 power
replicates, about twenty-five minutes, on the pinned BugSigDBExports v1.3.1 inputs.
`just calibration-report --replicates 5 --power-replicates 5` is the quick pass. It
reports every judged criterion as `insufficient_budget` and every procedure as
`not_certified`, because five replicates can neither meet the tolerance bound nor
show that a procedure exceeds it; it is a smoke test, not a calibration. (The form documented here
previously, `just calibration-report -- --replicates 5`, fails: `just` forwards the
literal `--` and argparse rejects it.)

`tests/test_calibration.py` holds the deterministic CI-sized checks. The #240 review
found five mutations that left the whole suite passing - selecting the least extreme
traits, reporting the raw p instead of the q-value, pooling only the study-flip null,
dropping the c/(c-1) factor, and changing a `<=` to `<` - so "each test was
mutation-tested" was not true of the suite as a whole (#362). Those five, and every
other `<=` on a p or q threshold, now each fail at least one test; that set of
mutations is what has been checked, not every possible defect. The tests
assert that the resampler relabels whole units and keeps the corpus paired; that the
injection leaves every non-injected trait's study contrast identical to the bit; that
the effective lift is reported rather than the nominal one; that the amended
criterion fails a true rate of 0.06 to 0.08 at the declared budget **and** fails
16/200 at the v1 budget, which the v1 rule passed; that the criterion's analytic
operating characteristics match sampled ledgers; that the family rule accounts for
Monte Carlo error; that the negative control's unmet expectation is surfaced; that
power without an interval is not reported as a pass; that zero rejections in ten
trials do not certify a zero rate; that a pass below the declared budget is withheld;
that the exported overall status keeps every failed criterion visible; that an
end-to-end quick run certifies nothing; that a replicate stream is
addressable on its own; and that the shared artefact matches the scenario it names.

The run that produced the numbers above: code commit `8ec18b1`; seed 182; 1400 + 1400 + 400 + 400 + 400 +
400 completed replicates; BugSigDB vintage `2026-04-24_00:41_UTC`;
`study_trait_profiles_wide.tsv` sha256 `6e70724e...92c6a73`; `experiments.tsv` sha256
`403afea9...f30ad230`; `full_dump.csv` sha256 `af84bfca...0fdd5788`; corpus 4,670
paired profiles over 2,335 experiments and 959 studies after the human-host and
pairing filters. The full digests, the protocol digest and the analysis-source digest
are in `calibration_status.json`, which also carries each procedure's tested level,
family, per-criterion verdicts and overall status for the integrated evidence table
of #178.

The protocol, code and this description were drafted with the assistance of a coding
agent (Claude Code) and require review by the author before manuscript use.

## Calibration-v3 eligibility amendment (#730, 2026-10-19)

The historical results above use calibration-v2. Review found that its guild
McNemar adapter retained experiment pairs with a zero annotation total on one
side, although the integrated production analysis excludes the entire pair.
On the current core, 94 pairs are affected: eligible support is 939 studies and
4,482 profiles, compared with the 959 studies and 4,670 profiles used to generate
the flip assignments. Guild discordant study counts change from 420/290 to
420/293. A v2 guild verdict therefore does not certify the production procedure.

Calibration-v3 applies the existing whole-pair exclusion inside each guild
replicate and exports its policy and support. It preserves the complete loaded
corpus for deterministic flip assignment generation, the seed 182, 1,400
replicates per null, 400 per injected scenario and the acceptance criteria.
The protocol amendment records what was already observed. Historical outputs
remain historical; a fresh full-budget execution is required for a v3 finding.

That execution completed from clean revision
`58b2ce841a88725603624007ca43d12c4a541c9d` with all declared replicates.
The [follow-up review](../docs/reviews/2026-10-19-ora-refresh-follow-up-review.md)
records its exact family-error bounds, retained historical comparison and
integration checks. The three directional procedures meet their
declared criteria; the screened negative control remains failed. Passing the
power-reporting criterion does not imply sufficient power for small signals,
and these results do not calibrate the native ORA or composition methods.
