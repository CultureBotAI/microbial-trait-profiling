# Gram-stain directional Sankey diagrams

Figures: `gram_positive_inc_sankey.png`, `gram_positive_dec_sankey.png`,
`gram_negative_inc_sankey.png` and `gram_negative_dec_sankey.png`, with the edge
table `gram_sankey_edges.tsv`, written to `data/derived/gram_sankey/`. Produced by
`src/bsdbng/gram_sankey.py`, driven by `just gram-sankey`.

## Terms and corpus

Gram stain classes were taken from KG-Microbe as METPO:1000698 (gram positive)
and METPO:1000699 (gram negative); both labels were read from the merged
KG-Microbe node table. The corpus and per-profile trait fractions come from the
shared loader used by the oxygen guild and Table 1 figures: the **pinned
BugSigDBExports v1.3.1 release** (Zenodo 10.5281/zenodo.19740674; export
timestamp 2026-04-24_00:41_UTC), human samples only, experiments carrying both an
increased and a decreased signature, and MediaDive growth-medium columns removed
before any total is computed. That leaves 317 traits, 2,335 experiments, 4,670
profiles and 959 studies. The loader rejects blank count cells rather than
reading them as zero.

A further 94 experiments are removed whole because one side has no retained
annotation at all. Such a profile has no defined fraction; scoring it as 0 would
pull its study's contrast toward zero and could create the exact-zero contrasts
that carry no direction. No experiment was removed for missing condition or body
site.

## The paired contrast, and why it does not remove the mirror artefact

Gram positive and gram negative annotations are close to exhaustive and close to
equal across KG-Microbe taxa, so counting how often each class appears on the
increased and decreased sides produces two near-mirror-image totals by
construction (issue #132). Per-direction hit counts are therefore not reported.
Instead each unit - one study within one condition and body site - contributes a
paired contrast for each class: its mean increased fraction minus its mean
decreased fraction, averaged over that unit's experiments. The contrast is
computed by the shared helper `guild_maps.grouped_direction_means`, which raises
if a unit's increased and decreased profile counts differ, so it is a paired
contrast by construction. A contrast within 1e-12 of zero is a tie.

The paired contrast does not escape the artefact; it moves it one level up. The
two classes are near-complementary shares of the same profile denominator, so a
unit's two contrasts usually point opposite ways: of 1,219 units, the
gram-positive and gram-negative contrasts have opposite signs in 957, the same
sign in 250, and at least one tie in 12. Each gram-negative diagram is therefore
largely a mirror of a gram-positive one, and the four figures are not four
independent readouts.

## Tests and what is drawn

Directionality is tested rather than asserted (issue #223). For every gram class,
condition and body site with at least three non-tied units, the numbers of units
shifting up and down are compared with an exact two-sided binomial sign test,
Benjamini-Hochberg corrected across all such flows of that gram class. For every
condition and body site, the units whose gram-positive and gram-negative
contrasts disagree are compared with an exact McNemar test (a binomial test on
the discordant units), the paired test the oxygen guild figure uses for aerobes
against anaerobes, Benjamini-Hochberg corrected across condition and body-site
pairs.

The unit is deliberately finer than the study, so it is worth stating what that
does and does not cost. A flow is one condition and body site, and a unit is one
study within one condition and body site, so **within any one flow a study votes
at most once** - each sign test's n is a count of distinct studies, not of
experiments, which is the sampling unit issue #133 asks for. The 1,219 units come
from 939 distinct studies, and the 169 studies that span more than one condition
or body site contribute a unit to each. That dependence is therefore between
flows, not inside any one of them: the individual sign tests are unaffected, and
what it weakens is the assumption of independence in the Benjamini-Hochberg
family taken across the 174 flows. Since no flow passes its sign test on this
corpus, the point is a limitation on the family rather than on any reported call.

174 flows met the three-unit floor. 153 have a net direction and 21 split evenly.
A flow is drawn once, in the figure for its net direction, with link width equal
to the net number of units (up minus down), so a six-to-five split is drawn as one
unit wide rather than as two near-equal flows. Flows whose sign test passes
q <= 0.05 are drawn dark and the rest faded. On this corpus **none passes**: the
smallest sign-test q is 0.32 (gram negative, colorectal cancer, feces, 22 units
up and 6 down), and only 2 of 174 flows reach an unadjusted p below 0.05. No
condition and body site passes the McNemar test either (smallest q 0.63; 1 of 87
pairs below an unadjusted 0.05). The diagrams therefore show where the unit
counts lean, not tested directional shifts, and each figure says so in its
subtitle.

The edge table lists every tested flow, drawn or not, with its up, down and tied
unit counts, net units, mean contrast, sign-test p and q with its family, and the
McNemar discordant counts, p and q, so every link can be checked without reading
pixels. All quantities are shares of trait annotations attached to curated taxon
signatures, not abundances, and the signature direction has not been verified as
disease versus control (issue #188).

## Reproducing

`just gram-sankey` runs the analysis deterministically and stamps the BugSigDB
export timestamp, read from the dump's leading comment line, into the edge table
and every figure. PNG export needs a headless Chrome (`just get-chrome`). The
analysis does not reuse the enrichment ranking tables, whose hypergeometric urn
draws taxon-trait assertions rather than taxa; no p-value from that procedure
enters these figures. An earlier version of this note reported 255 edges over
3,736 experiments from the unpinned 2026-04-18 export; those numbers are
superseded.

The pipeline code and this description were drafted with the assistance of a
coding agent (Claude Code) and reviewed by the author.
