# Directional consistency over the full retained trait space

This describes `just trait-direction-screen`, which runs the directional-consistency
arm's study sign test over every retained trait of the pinned corpus rather than over
the six traits of manuscript Table 1 (issue #255), and reports the significant set at
both the nominal and the effective family size (issue #265). The nominal BY set is the
reported result; the effective-size set is a sensitivity bound with no FDR guarantee.

## Why the family changed

The six Table 1 traits were selected as the intersection of the composition and
enrichment arms' significant sets, computed on this same corpus. Testing only those six
made the directional arm a third statistic on a feature set the other two arms had
already agreed on: its q-values were conditional on a selection the test did not model,
it could not discover a seventh trait, and its agreement with the other arms was not
independent confirmation, because the features had been chosen for having already
agreed. The six remain in the table as a flagged subset. They are no longer the family.

## What is unchanged

The statistic. For each candidate trait and each study, the mean increased-side fraction
minus the mean decreased-side fraction gives one contrast per study; its sign is one
vote. A study counts as shifted up when the contrast exceeds the tie tolerance of 1e-12,
down when it falls below the negative tolerance, and tied otherwise. Tied studies are
excluded from the test and reported separately. Each trait then gets a two-sided exact
binomial test of its non-tied studies against a half-and-half null, with an exact
Clopper-Pearson 95% interval for the share shifting up.

The study contrasts and the sign tally are computed by `study_trait_contrasts` and
`tally_direction_signs` in `src/bsdbng/table1_trait_maps.py`, which the six-trait arm
also calls. The widened screen shares that code rather than reimplementing it, so the
two cannot drift apart. Regenerating the six-trait outputs after this change produced
files byte-identical to the previous run.

The denominator is also unchanged: each profile's fraction is over the sum of every
retained non-MediaDive column, including traits outside Table 1. Selecting a candidate
never redefines that total.

## The candidate family

Every retained column of `data/derived/study_trait_profiles_wide.tsv` — all columns
except the three identifier columns and every MediaDive growth-medium column. On the
pinned v1.3.1 corpus, under the `off` scope, that is **317 candidates** (344 under the
primary `phenotype` scope), with 2,333 MediaDive columns excluded,
tested over 4,670 paired human profiles from 2,335 experiments and 959 studies.

The set is read from the matrix at run time and never written down, so a corpus change
moves the family rather than silently mismatching it. `retained_trait_columns` in
`src/bsdbng/trait_guilds.py` is the single implementation of that rule and is also what
the profile loader uses to build its denominator, so the tested universe and the
denominator universe are the same set by construction.

The candidate rule, the correction, the estimators and the reported thresholds are
declared in [`trait_screen_manifest.json`](trait_screen_manifest.json), written before
the first run. `methods/three_arm_statistics.md` requires a larger trait screen to
declare its own candidate manifest rather than reuse the frozen six-trait family.

## Correction, and two family sizes

The primary correction is **Benjamini-Yekutieli** across the full family, which is valid
under arbitrary dependence; PR #181 already applies BY to the 15-ratio family for the
same reason. Benjamini-Hochberg over the same family is reported beside it, because the
other two arms use BH and a reader comparing arms needs the comparable column.

Candidates that could not be tested are held at p = 1 inside every correction and left
blank in the output, the convention PR #181 set for the ratio family. Dropping them
would shrink the family they were declared into.

Correcting across m traits (317 under `off`, 344 under `phenotype`) implies m
independent opportunities to be wrong. It is
not obvious that the data contain that many, so the significant set is reported twice:
once at the nominal family size, and once with the same procedure run as though the
family held its effective number of tests. `step_up` in `src/bsdbng/trait_screen.py`
takes the family size as a parameter for that reason — `scipy.stats.false_discovery_control`
always corrects at the number of p-values it is handed. With the family size equal to
the number of tests it reproduces scipy exactly, which the test suite asserts against
scipy on random input for both methods. An adjusted p-value is never reported below its
own raw p-value; that floor binds only on the effective-size run.

Effective family size is estimated from the eigenvalue spectrum of the Spearman
correlation between candidates, by **Li & Ji (2005)** and **Galwey (2009)**.
**Cheverud/Nyholt is computed, reported, and excluded**: Li & Ji (2005) report that its
eigenvalue-variance estimate is overly large and gives excessively conservative results,
and Salyakina et al. (2005) found its nominal error rate off in both directions depending
on the dependence structure. On this corpus it disagrees with the other two by a wide
margin. It is excluded as unreliable, not as lenient. Here it is the largest of the three (294.4 against Li & Ji 112.0 on the
contrast basis, and 198.9 against 60.0 on the profile basis), so it would give the
strictest correction. Leaving it out makes the effective-size sensitivity run less strict,
not more. It is named in the output rather than silently omitted.
All three return the candidate count for an independent family and 1.0 for a perfectly
redundant one, which the test suite checks at both limits.

Because the substitution is monotone, the nominal significant set is always contained in
the effective one. The difference between them is candidates gained, not candidates
swapped, and the family table reports both the overlap and its size.

### Which set is the reported result

The reported significant set is **BY at the nominal m = 317** (the `off` scope; m = 344
under `phenotype`), with BH at m = 317 beside it for comparison with the other arms.
On the pinned corpus that is 18 traits at BY q <= 0.05. The effective-size run (m = 112 by Li & Ji on the study-contrast
basis, 56 traits) is a sensitivity bound only, and carries **no FDR guarantee**: Li & Ji
and Galwey were derived to set per-test thresholds for familywise error, so putting
their m into a BH or BY step-up is a heuristic, and BY at the nominal m already
controls FDR under arbitrary dependence, so shrinking m on top of it would adjust for
dependence twice. Because the count moves from 18 to 56 with the assumed dependence,
it is not a stable quantity: report it as unstable, and do not quote either number
alone. This follows the ruling on issue #265 (2026-09-28); the manuscript wording is
issue #205.

## Two correlation bases, and why they disagree

Dependence is measured on two bases and both are reported.

The **study-contrast basis** is the correlation between candidates across the 959
study-level contrasts that are actually tested. This is what drives the effective family
size used in the correction, because it is the dependence between the tested quantities.

The **profile basis** is the correlation between candidates across the 5,829 human
profiles, which is the basis quoted in issues #263 and #265. It is retained so that
those figures can be checked against this code rather than restated. They reproduce:
median Spearman +0.611 over 50,086 pairs, 63.5% above +0.5, first component 60.4% of the
variance, Li & Ji 60.0, Galwey 48.9, Cheverud/Nyholt 198.9. An integration test asserts
each of these.

The two bases disagree sharply, and that is a result rather than a discrepancy. On the
tested contrasts the median Spearman is **+0.089**, only 7.6% of pairs exceed +0.5, the
first component holds **20.6%** of the variance, and the effective family size is
**112.0** by Li & Ji rather than 60.0. Much of the between-trait correlation seen across
profiles is shared annotation depth — a profile with more annotations has more of every
trait — and taking the paired within-study contrast removes a large part of it. The
correction is therefore less aggressive than the profile figures alone would suggest,
and the "effective m of roughly 50 to 60" quoted for profiles is not the effective m of
this family.

A candidate that does not vary across the units is dropped from the correlation and
counted, because Spearman is undefined against a constant series. It stays in the tested
family; only the dependence estimate loses it. If the corpus is too small to estimate
dependence at all, the effective family size is undefined, the effective column is blank
throughout, and the family table says so rather than carrying a fabricated number.

## Missing states

Three states are kept distinct and never collapsed into a number:

| State | Meaning |
| --- | --- |
| `tested` | at least one study contributed a non-tied contrast |
| `all_tied` | every eligible study contrast fell within the tie tolerance |
| `no_eligible_studies` | no study contributed a contrast for this candidate |

An untested candidate has no share, no interval, no p and no q. It does not have a share
of 0.0 or a p-value of 1.0. On the pinned `off` corpus all 317 candidates are testable.

Every paired study contributes a contrast for every candidate, so the two untested states
arise differently. A candidate that is zero, or otherwise equal, on both sides of every
study screened is `all_tied`. If the screen were run on a subset of studies (the CLI has
no such option today), this is how a trait absent from that subset would show up.
`no_eligible_studies` cannot single out one candidate. Every candidate reads it at once
when the screened profiles hold no paired study at all; on the real pipeline that means
`load_profiles` kept no paired human study.

## Outputs

Written to `data/derived/trait_screen/` (gitignored, regenerable):

| File | Contents |
| --- | --- |
| `trait_direction_screen.tsv` | one row per candidate: identity and label, Table 1 membership, up/down/tied/scored studies, share up with its exact interval, binomial p, BH and BY q at the nominal family size, BY q at the effective one, test status, and significance flags at both thresholds |
| `trait_screen_family.tsv` | the family record: nominal and effective size, every estimator including the excluded one, both correlation bases, and the significant-set counts and overlap at each threshold |
| `trait_screen_correlation.tsv` | the eigenvalue spectrum of each basis, with per-component and cumulative variance share |
| `trait_screen_provenance.json` | input SHA-256 hashes, code commit, manifest identity and hash, corpus filter counts, and every analysis rule |
| `trait_screen_report.md` | a readable summary: the six in the widened family, and the strongest candidates by p-value |

Candidate labels come from `data/derived/taxa_trait_matrix_columns.tsv`, the column index
the pipeline emits beside the matrix it built. Without that file the screen still runs
and every candidate falls back to its CURIE, which the provenance record states.

## What this does and does not establish

It measures relative annotation representation between increased and decreased
signatures. It does not measure absolute trait abundance, pathway activity, oxygen
exposure, or causation, and agreement with another arm on this shared input is not
independent replication.

Widening the family removes the circularity of testing traits the other arms had already
chosen, but it does not address calibration (#182, #263), the assertion-level sampling
unit (#131), reference and population sensitivity (#183), or denominator harmonisation
across arms (#174). The remaining selection effect now sits in one place: the definition
of Table 1 itself.

The six traits' q-values under this family are not the six-trait q-values published in
[`table1_trait_condition_body_site_maps.md`](table1_trait_condition_body_site_maps.md),
and this screen does not change them. Revising any manuscript number is #205.
