# Sensitivity of the paired log-ratios to eligibility, reference and population

Produced by `src/bsdbng/log_ratio_sensitivity.py`, driven by
`just table1-log-ratio-sensitivity`, which writes
`data/derived/table1_composition_sensitivity/`: every scenario row
(`sensitivity_results.tsv`), eligibility by ratio and by exact body site and
condition (`eligibility_by_ratio.tsv`, `eligibility_by_stratum.tsv`),
leave-one-study-out influence (`leave_one_study_out.tsv`), a status record for the
integrated evidence table (`sensitivity_status.json`) and a generated table of
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

## Declared before the first run

The scenarios, their populations, their correction families, the 30-study minimum
for a test and the reporting threshold were written to
[`table1_log_ratio_sensitivity_protocol.json`](table1_log_ratio_sensitivity_protocol.json)
before any result was computed, and were not changed afterwards. Each scenario
keeps its own family, no scenario is preferred, and q <= 0.05 is used only to count
rows. The same test is used throughout: study-level mean contrasts, an exact
two-sided binomial test of the non-tied studies shifting up, and an exact 95%
interval for that share (not for the mean effect).

| Scenario | Population | Family |
| --- | --- | --- |
| S0 historical fractions | all 2,335 complete pairs; a zero-total profile scores 0 | BH over 6 |
| S1 fraction baseline | 2,241 pairs with positive totals on both sides | BH over 6 |
| S2 complete-case ratios | per ratio, pairs with all four counts positive | BY over 15 |
| S3 matched fractions | each ratio's numerator and reference fractions on exactly that ratio's S2 pairs | BY over 30 |
| S4 six-trait common support | 354 pairs (239 studies) with all six traits positive | BH over 6 fractions; BY over 15 ratios |
| S5 triple common support | per trait triple, pairs with all three positive; A/B, B/C, A/C | BY over 60 |
| S6 / S7 fecal strata | S2 within exact body site Feces / outside it | BY over 15 each |
| S8 leave one study out | S2 with each study removed in turn | BY over 15, other 14 p held fixed |

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

**Most fraction signals depend on the population (S1 against S3).** Restricted to
a ratio's eligible pairs, facultatively aerobic keeps its increased-side share only
on the facultatively aerobic/coccus shaped pairs (0.567, q = 0.010) and loses it on
the other four (shares 0.531 to 0.554, q >= 0.15). Coccus shaped, sphere shaped,
arsenate reduction and moderately halophilic lose theirs on every restricted
population (q >= 0.098). Cellulose catabolic process is the exception: its
decreased-side share holds on all five of its ratio populations (0.354 to 0.403,
q <= 0.003).

**The three significant ratios are carried by cellulose (S3 against S2).** The
ratios with BY q <= 0.05 are facultatively aerobic/cellulose (share up 0.595,
q = 0.0011), cellulose/sphere shaped (0.407, q = 0.0026) and cellulose/arsenate
reduction (0.365, q = 0.0011). On exactly the same pairs the cellulose fraction
alone shows the same direction and similar strength (0.403, 0.401 and 0.354), while
the facultatively aerobic, sphere shaped and arsenate reduction fractions do not
(0.531, 0.499 and 0.510, all q = 1). On these populations the ratio results are
therefore consistent with a relative decrease of cellulose-associated annotations
rather than with a shift in the other member of each pair.

**On one shared population the pattern is the same (S4, S5).** On the 239 studies
with all six traits positive, cellulose has the only strong fraction shift (share
0.356, q = 5.7e-05), moderately halophilic a weaker increase (0.590, q = 0.019), and
the significant ratios are the four with cellulose against facultatively aerobic,
sphere shaped, arsenate reduction and moderately halophilic. Across the 20 triple
populations, all 13 rows with q <= 0.05 involve cellulose, and no ratio without
cellulose reaches q <= 0.05 in any population. Ratios without cellulose do not form
a stable ordering: their shares stay between 0.43 and 0.60 and some change side
between populations.

**The signal is fecal (S6, S7).** Within Feces (602 studies; a study with experiments both in and outside Feces
counts in each stratum), facultatively
aerobic/cellulose, cellulose/sphere shaped, cellulose/arsenate reduction and
cellulose/moderately halophilic reach q <= 0.02. Outside Feces (381 studies across
all other exact body sites) no ratio does; the cellulose ratios sit at shares of
0.42 to 0.52. Whether that reflects biology, fewer eligible studies per site or
heterogeneous pooling of many sites cannot be separated here.

**No single study drives a ratio (S8).** Removing any one eligible study moves a
ratio's share up by at most about 0.003 and changes no ratio's q across 0.05.

## Interpretation limits

The stable finding across these declared scenarios is narrow: in fecal signatures
with positive counts, cellulose-associated annotations are relatively lower on the
increased side than several other Table 1 annotations, and that relative shift
survives restriction to matched populations. It is conditional on positive
annotation counts, is measured on over-represented, richly annotated signatures,
and does not identify an absolute change in any trait. Most of the six original
fraction directions do not survive restriction to ratio-eligible populations; that
is a population effect, not evidence that those directions are false in the full
corpus. Agreement between scenarios is robustness to these analytical choices on
the same annotation data, not independent replication. The six traits were
selected on overlapping data, the signature directions have not been verified as
disease versus control (#188), none of these procedures has been calibrated (#182),
and annotation coverage is not yet available per profile (#184).

## Reproducing

`just table1-log-ratio-sensitivity` on the pinned BugSigDBExports v1.3.1 inputs
(export timestamp 2026-04-24_00:41_UTC). The status JSON records the code revision,
the SHA-256 of the profile matrix, experiment table, taxa table, analysis source
and protocol, the corpus filter cascade, and the per-scenario test counts. All
statistics use scipy (`binomtest`, `false_discovery_control`).

The protocol, code and this description were drafted with the assistance of a
coding agent (Claude Code) and reviewed by the author.
