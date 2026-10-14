# The widened directional screen under the complete null (#263)

This describes `just calibration-screen`. It runs the 317-trait directional
screen of `just trait-direction-screen` (#255, #265) on copies of the corpus in
which no trait differs between the increased and decreased signatures, and
counts how many traits each correction calls. It is the directional arm's share
of the calibration battery in #263, stacked on PR #240 (the calibration harness
and shared flip artefact) and PR #276 (the screen).

## Design

The null is the directional arm's complete study-flip null from
[`calibration.md`](calibration.md), reused rather than re-implemented. For a
random half of the 959 studies, the increased and decreased labels of every
profile in the study are exchanged; counts are untouched, so annotation
correlations, unequal experiment counts per study, zeros and the body-site
imbalance survive by construction. The flips come from the same
`complete_null_study_flip` stream and seed (182) as the shared flip-assignment
artefact, so replicate *i* here flips exactly the studies that row *i* of
`calibration_flip_assignments.tsv` lists. A test checks that equality.

Each replicate runs the production screen unchanged: one vote per study per trait,
a two-sided exact binomial test against 0.5, and Benjamini-Yekutieli (BY) and
Benjamini-Hochberg (BH) across all 317 traits, plus BY at the effective family size
re-estimated from that replicate's study contrasts. Under this null every call is a
false call. The per-replicate count is therefore the number of false discoveries,
and the share of replicates with any call is the family error of the whole screen,
reported with an exact binomial interval over independent replicates.

## Result

100 replicates, seed 182, code commit `0cc6c66`, on the pinned v1.3.1 corpus
(BugSigDB vintage `2026-04-24_00:41_UTC`; 4,670 paired profiles over 2,335
experiments and 959 studies). The effective family size ranged from 112 to 114
across replicates.

| Correction | Threshold | Observed calls | Mean null calls | Most in one replicate | Replicates with any call (95% CI) | Null replicates reaching the observed count |
| --- | --- | --- | --- | --- | --- | --- |
| BY, m = 317 | q <= 0.05 | 18 | 0.00 | 0 | 0/100 (0 to 0.036) | 0 |
| BY, m = 317 | q <= 0.001 | 0 | 0.00 | 0 | 0/100 (0 to 0.036) | not applicable |
| BH, m = 317 | q <= 0.05 | 75 | 0.02 | 1 | 2/100 (0.002 to 0.070) | 0 |
| BH, m = 317 | q <= 0.001 | 7 | 0.00 | 0 | 0/100 (0 to 0.036) | 0 |
| BY, effective m | q <= 0.05 | 56 | 0.01 | 1 | 1/100 (0.0003 to 0.054) | 0 |
| BY, effective m | q <= 0.001 | 5 | 0.00 | 0 | 0/100 (0 to 0.036) | 0 |

## What this does and does not show

Under the complete null, the screen calls almost nothing: at most one trait in any
replicate under any correction, and none under BY at m = 317, the reported
correction. No null replicate comes near the observed 18, 56 or 75. The
effective-size correction, which carries no FDR guarantee (#265), also stays
within a family error of 0.054 at its upper bound on this null.

That is a statement about the complete null only. It does not show that the 18
traits are true effects, and it does not measure the false-discovery proportion
when some traits really do shift, because a partial null - some traits shifted,
others not - is not simulated here. It also says nothing about the selection that
chose the six Table 1 traits; this screen tests all 317 and selects nothing. One
hundred replicates bound the family error at 0.036 when no replicate calls
anything, which is below the 0.065 tolerance of
[`calibration_protocol.json`](calibration_protocol.json), but this screen is not
one of that protocol's declared procedures and no verdict is issued against it.

Like the calibration harness ([`calibration.md`](calibration.md)), the study-flip
null treats studies as independent. It cannot detect dependence *between*
studies (shared cohorts, overlapping populations, or reused datasets); any
resulting inflation remains out of reach of this null.

## Reproducing

`just calibration-screen` writes `screen_null_replicates.tsv` (long: seed,
stream, unit, replicate, q column, correction, family basis and size, threshold,
calls), `screen_null_summary.tsv` and `screen_null_provenance.json` (seed,
completed replicates, input digests, code commit) into
`data/derived/calibration/screen/`. `--replicates` and `--seed` change the
budget and the stream. A run takes about two minutes.

`tests/test_screen_calibration.py` checks that each correction and threshold is
counted once, that the null flips are the shared artefact's flips, that whole-study
flipping removes a perfectly consistent synthetic signal, and that the summary's
family error and interval are computed over replicates.

The code and this description were drafted with the assistance of a coding agent
(Claude Code) and require review by the author before manuscript use.
