# Do these directional tests hold their error rates?

Produced by `src/bsdbng/calibration.py`, driven by `just calibration-report`, which
writes `data/derived/calibration/`: every scenario's per-test and family results
(`calibration_results.tsv`), the declared criteria and their verdicts
(`calibration_criteria.tsv`) and a status record for the integrated evidence table
(`calibration_status.json`). It implements
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
six Table 1 traits - a selection made on this same corpus by another arm - to measure
what selecting and then testing the selection costs. Selecting three of six rather
than six of hundreds makes the measured inflation a lower bound.

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
null. Results below are near-identical under the two, which is itself informative:
on this corpus the within-study dependence is not what drives the error rate.

Power is measured by raising one trait's counts (coccus shaped) on the increased side
of half the studies by a declared multiplier, moving each profile's total with the
injected counts so the other traits' shares stay coherent. The traits left untouched
give the false-discovery proportion in the same replicates.

200 replicates per scenario, seed 182. A rate near 0.05 therefore carries a Monte
Carlo interval of roughly plus or minus 0.03; the intervals are exact binomial and are
in the results table. The tests are discrete, so the attainable p-values are a finite
set - the smallest observed under the null was 4.3e-04 - and uniformity of p is
neither expected nor required.

## Declared criteria and their verdicts

The protocol, including these criteria, was written before any result was computed
([`calibration_protocol.json`](calibration_protocol.json)).

| Criterion | Procedures | Verdict | Evidence |
| --- | --- | --- | --- |
| Per-test rejection interval includes or falls below 0.05 | fraction, log-ratio, guild | **pass** | highest lower interval bound 0.031, at cellulose catabolic process and at facultatively aerobic/cellulose |
| At most 10% of replicates give any q <= 0.05 | fraction, log-ratio | **pass** | highest rate 0.050 (fraction), 0.010 (log-ratio) |
| Power reported with its interval | fraction, log-ratio | **pass** | curve below |
| Selection inflation measured, not judged | screened | **reported** | per-test rejection rate 0.050 to 0.115 under the null |

One amendment is recorded in the protocol: the single injected multiplier of 1.3 was
replaced by a curve of 1.1, 1.3, 2 and 4 after the first run showed a 1.3x lift moves
almost no study signs. The null scenarios, the criteria and every threshold are
unchanged, and no procedure was altered in response to any result.

## What the numbers say

**Under the complete null the published procedures hold their rates.** Mean per-test
rejection at p <= 0.05, over 200 replicates:

| Procedure | Study flip | Experiment flip | Highest single test |
| --- | --- | --- | --- |
| fraction_sign_test (6 tests) | 0.042 | 0.041 | 0.060 |
| log_ratio_sign_test (15 tests) | 0.041 | 0.041 | 0.060 |
| guild_mcnemar (1 test) | 0.060 | 0.045 | 0.060 |
| screened_fraction_sign_test (3 of 6) | 0.083 | 0.081 | 0.115 |

At the family level, the share of replicates with any q <= 0.05 is 0.050 and 0.045
for the six-trait BH family, 0.010 for the 15-pair BY family, and 0.060 and 0.045 for
the single McNemar. Benjamini-Yekutieli is visibly conservative here, as its
arbitrary-dependence guarantee implies: it spends power to be safe under the
dependence between overlapping ratios.

**Power is modest for small effects.** For the injected trait, rejection at q <= 0.05:

| Injected lift | fraction_sign_test | log_ratio_sign_test (mean of the 5 pairs naming the trait) |
| --- | --- | --- |
| 1.1x | 0.020 | 0.004 |
| 1.3x | 0.095 | 0.021 |
| 2x | 0.775 | 0.426 (best pair 0.900) |
| 4x | 1.000 | 0.862 (best pair 1.000) |

A 30 per cent lift in one trait's annotation counts across half of 959 studies is
detected about one time in ten. **A non-significant result from these procedures is
therefore weak evidence of no effect**, which matters wherever the analyses report
that something does not reach a threshold - the non-fecal strata, the eleven
non-cellulose ratios, and every "no cell survives correction" statement. False
discoveries among the untouched traits stayed near nominal in the same replicates
(mean false-discovery proportion 0.03 to 0.09 for the fraction family, below 0.01 for
the ratio family).

**Selecting traits on the data and testing the selection inflates the error rate.**
The screened procedure rejects a true null 0.081 to 0.083 of the time per test
against a nominal 0.05, and its false-discovery proportion under an injected signal
reaches 0.16. The six Table 1 traits were selected on this corpus by the composition
arm's criteria over hundreds of traits, so the published q-values for those six are
not post-selection valid, and the measured inflation here is a lower bound on the
cost. This is the concrete statistical basis for treating the six-trait results as
**exploratory** rather than confirmatory, and it cannot be repaired by correcting
within the six.

## What this does not cover

- The enrichment arm's hypergeometric procedure. Its urn draws trait-annotation
  assertions rather than taxa (#131), which is a model question a resampling study
  of this shape cannot settle, and its adapter belongs to its owner.
- The composition arm's native transformed test, for the same reason.
- The aggregated profile-support rule of #177, whose code is not on main.
- Calibration is not validity of interpretation. A procedure that holds its error
  rate still measures annotation shares, and none of this addresses orientation
  (#188), coverage (#184) or biological meaning (#179).

There is no untouched validation set: every result in this repository is computed on
one corpus, and the six traits were chosen on it. A retrospective split would not be
independent validation, and none is claimed.

## Reproducing

`just calibration-report` (about two minutes at 200 replicates) on the pinned
BugSigDBExports v1.3.1 inputs; `just calibration-report -- --replicates 5` for a quick
pass. `tests/test_calibration.py` holds the deterministic CI-sized checks: the
resampler relabels whole units and keeps the corpus paired, the injection touches
only the increased side, a fixed seed reproduces a run exactly, and an
anticonservative procedure is reported as failing rather than passing. The status
JSON carries the version, seed, replicate count, each procedure's tested level and
family, and its calibration verdict, for the integrated evidence table of #178.

The protocol, code and this description were drafted with the assistance of a coding
agent (Claude Code) and require review by the author before manuscript use.
