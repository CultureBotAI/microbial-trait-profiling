# Denominator-effect review and common study-weighted quantities

Run `just table1-common-quantities`. The implementation is
`src/bsdbng/common_trait_quantities.py`; it writes the machine-readable common
table, exact pair and study manifests, denominator audit, provenance JSON, and a
generated readable report to `data/derived/table1_common_quantities/`. This is the
methods decision and comparison deliverable for
[#256](https://github.com/CultureBotAI/microbial-trait-profiling/issues/256).

## Decision

Use the denominator-harmonized, study-weighted fraction contrast as the common
six-trait descriptive baseline. Report all 15 oriented pairwise log-ratios beside
it, but for every ratio first report both constituent fractions recomputed on that
ratio's exact eligible experiment pairs. This three-row comparison separates:

1. a denominator change: the focal counts stay fixed while another retained
   component changes a fraction but not a focal/reference ratio;
2. a population change: a full-population fraction changes when restricted to the
   ratio's positive-count pairs; and
3. an estimand change: the two matched fractions are replaced by a named A/B
   comparison on the same pairs.

No reference is selected by its q-value. The common table has no p-value, q-value,
or correction family and is not independent replication. Native transformed
composition tests, ORA, and the existing fraction/ratio sign tests retain their own
hypotheses and correction families. Full-pipeline calibration remains with #182.

## Pinned comparison and denominator audit

The run pins the exact profile-matrix Git blob
`d837f449e7e967c478d10f92938923e57eec9a19` (SHA-256
`c07d589f6cc88724e43c683ebd132bcd5c7f316f36b95df3ff9acad24ebbf62f`),
the experiment table SHA-256
`403afea98d9b0a5b56939f68bcc9aace56a8eda350232a586806fe7bf30ad230`,
and BugSigDB export `2026-04-24_00:41_UTC`. The PR #197 reference is commit
`309f9b355b157ee496f844413b5dd2e275f6a3a0`; its committed denominator-parity
artifact has SHA-256
`bf6bf566a91a746b3945b271a0dce71f878f680a40b2b7d217b6c32ddc6cb8d5`.
The recipe extracts both pinned objects from Git, so an ignored local derived
matrix cannot silently change the comparison.

The matrix has 317 retained non-MediaDive columns and 2,290 excluded MediaDive
columns. The retained-column manifest and the exact source column for every Table
1 trait are emitted. Selection, support summaries, normalization, transformation,
pooled analysis, and stratified analysis in the inspected PR #197 implementation
use the full retained universe before selection. On every one of the 4,670 profiles
in the current shared base population, the current loader total equals PR #197's
`harmonized_total`; its independently computed Ilyas total also equals that value.
The legacy denominator differs for 4,495 profiles, exactly by the excluded
MediaDive total.

There is one remaining population bookkeeping discrepancy in PR #197, not a
denominator-value discrepancy. Its `paired` flag differs from the current shared
study/experiment/direction keys in 160 memberships: 80 current-base profiles are
not flagged and 80 other profiles are substituted. Those 80 substituted PR rows
also have `noel_total=NA`, even though the PR provenance reports 4,670 Noel matched
profiles. The exact keys are in `denominator_mismatches.tsv`. The common estimates
therefore use the current shared loader's 2,335 complete experiment pairs from 959
studies and do not reuse PR #197's paired-membership flag. This population
difference must be resolved in #174/#197 before calling the native composition
rerun matched to the common table.

PR #197 also retains zero-total profiles in its transformed native analysis and
weights experiment pairs directly. Those are its documented native analysis
choices, but they do not match this common contract. Here, 94 zero-total pairs are
excluded from the full fraction summaries, leaving 2,241 pairs and 939 studies; 20
studies are entirely excluded and 28 are partially observed. Transformations and
study weighting therefore remain separate from the denominator correction.

## Methods compared

| Approach | Estimand and null | Independent unit | Assumptions and zero/coverage policy | Denominator dependence diagnosed or avoided |
| --- | --- | --- | --- | --- |
| Harmonized fraction contrast | Mean/median across study means of `100 * (C_inc/T_inc - C_dec/T_dec)`, plus the proportion of non-tied study effects above zero. The common report has no tested null. | Study; eligible experiments are averaged first. | `T` is all 317 retained non-MediaDive columns. Exclude the whole pair if either total is zero; a focal zero with positive totals is valid. Blank, malformed, and negative counts fail before analysis. The exact 95% interval is only for the non-tied positive-study proportion. | Directly estimates representation relative to the retained total, so it exposes and remains susceptible to changes in unrelated retained components. It is invariant to excluded MediaDive counts. |
| Complete-case paired log-ratio A/B | Mean/median across study means of `log2(A_inc/B_inc) - log2(A_dec/B_dec)`, plus the same directional proportion. The common report has no tested null; the existing sign test remains separate. | Study; eligible experiments are averaged first. | All four named counts must be positive. Exclude the whole pair otherwise; add no pseudocount. Coverage-conditioned populations differ by ratio. | The profile total and unrelated retained components cancel, so this diagnoses closure through the shared total. It still depends on B and cannot identify whether A rose, B fell, or both changed. |
| Prespecified balance or multi-feature reference | A log contrast between geometric means of declared numerator and reference feature sets. A future test would need a predeclared null and calibration. | Study after experiment-level balance contrasts. | Requires biologically justified, fixed, non-overlapping reference sets and an explicit zero policy. Overlapping trait annotations violate a simple independent-component reading. | Cancels the shared total and spreads reliance across a reference set, but remains reference-, coverage-, selection-, and dependence-sensitive. It is a #183 sensitivity candidate, not the primary estimate here. |
| Zero-aware two-part analysis | Separately estimate the probability of observed annotation coverage and the positive-count log-ratio magnitude. Any tests would have distinct coverage and magnitude nulls. | Study, with paired coverage and positive-magnitude components. | Needs #184 to distinguish observed zero from unknown/malformed or low annotation coverage and #182 to calibrate a combined rule. A pseudocount-only analysis is not used. | Avoids silently dropping all zeros or pretending a replacement is observed, but it does not remove reference dependence or establish absolute abundance. Not yet implemented. |

The first two approaches are implemented because they answer the required common
questions with auditable eligibility. A balance is reasonable only after a
biological reference set is declared independently of the results. A zero-aware
model is premature while per-profile coverage validity is unresolved.

## Denominator-only and excluded-MediaDive invariants

Hold focal counts at A=10 and B=20 on both sides. Let an unrelated retained count
be 10 on the decreased side and 30 on the increased side. The totals are 40 and
60, so A moves from 25% to 16.7% (-8.33 percentage points) and B from 50% to
33.3% (-16.67 points), despite neither focal count changing. A/B stays 0.5 and its
paired log2 contrast is zero. Fractions therefore cannot distinguish this example
from a focal decrease without inspecting the component counts; A/B can show that A
did not move relative to B, but cannot say that either was absolutely unchanged in
an observational result.

Changing an excluded MediaDive column from zero to any positive count changes none
of the retained totals, fractions, ratio eligibility, or ratios. The focused tests
exercise both invariants, reference reversal, exact matched populations,
zero-total versus focal-zero handling, and unequal numbers of experiments per
study.

## Common results and interpretation

The six full-population fraction mean effects, in percentage points, are:
facultatively aerobic `+0.8901`, cellulose catabolic process `-0.06363`, coccus
shaped `+0.7249`, sphere shaped `+0.02134`, arsenate reduction `+0.006435`, and
moderately halophilic `+0.01861`. Their directional proportions and exact
intervals are in the generated report and machine table. These are annotation
fractions, not absolute trait abundance.

Ratio eligibility ranges from 243 to 676 studies and 358 to 1,305 experiment
pairs. For example, facultatively aerobic/cellulose uses 467 studies and 780 pairs.
On those exact pairs the mean fraction effects are `+0.5470` and `-0.07977`
percentage points, while the mean paired log2 ratio change is `+0.7906`. The full
facultatively aerobic fraction mean is `+0.8901`; moving to `+0.5470` is a
population change, while moving from the two matched fractions to `+0.7906` log2
is a change of estimand. The same decomposition is emitted for every oriented
pair.

The ratio means reproduce the existing complete-case log-ratio effects, but the
common table intentionally omits their sign-test p/q values. Broader common
support, alternative references, zero sensitivity, strata, and study influence
belong to #183 and must not be duplicated here. The output contract is ready for
#178: every row carries trait/reference identity, exact population ID, source
columns, denominator, eligibility, study/pair/partial-exclusion/tie counts, effect
units, uncertainty parameter, validity state, and calibration status.

Neither approach establishes absolute abundance. The recommendation is to retain
the full fractions as the common population-level description, use matched
fractions to diagnose selection into each ratio population, and use all 15 named
ratios as reference-explicit closure diagnostics. Biological claims require an
external load or absolute-abundance reference; statistical combination requires
the method-specific calibration in #182.
