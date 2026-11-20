# Gram-stain directional Sankey

Figures: `gram_positive_inc_sankey.png`, `gram_positive_dec_sankey.png`,
`gram_negative_inc_sankey.png`, `gram_negative_dec_sankey.png`, with the edge table
`gram_sankey_edges.tsv` and `gram_sankey_summary.json`. Produced by
`src/mtp/gram_sankey.py`, driven by `just gram-sankey`, written to
`data/derived/gram_sankey/` (#147, #223).

Interpret all outputs using the [shared guide](biological_interpretation.md). The
figures show trait-annotation shares of reported taxa, not abundance, cell density or
cell-wall chemistry measured in samples. `inc`/`dec` are curated comparison
directions; they are not verified as disease versus control (#188).

## Input corpus

The corpus is the one the oxygen guild figure reads
([`guild_condition_body_site_heatmap.md`](guild_condition_body_site_heatmap.md)):
`study_trait_profiles_wide.tsv` from `just study-trait-profiles`, which keeps human
hosts and 16S and shotgun (WMS) experiments only (#742), joined to `experiments.tsv`
for condition and body site. `load_profiles` keeps experiments carrying both an
increased and a decreased profile and removes MediaDive growth-medium columns before
any total is computed. Experiments with a zero retained-annotation total on a side are
removed whole (`drop_zero_total_pairs`), as are experiments with no condition or body
site. The two classes are METPO:1000698 (gram positive) and METPO:1000699 (gram
negative), each a single-term group in `GRAM_GUILDS`.

## Quantity

Gram positive and gram negative are near exhaustive and near equal across KG-Microbe,
so raw per-direction hit counts mirror each other by construction (#132). The unit is
therefore a paired contrast: for one study within one condition and body site, its mean
increased gram-class fraction minus its mean decreased fraction, averaged within the
study first through the guild figure's paired helper, which refuses a study whose
increased and decreased profile counts differ. Within one flow a study votes once, so
each sign test's n is a count of distinct studies (#133). A study that spans several
conditions or body sites contributes one unit to each, so flows are dependent with each
other through shared studies.

The contrast does not remove the mirror; it moves it one level up. Each gram-negative
diagram largely mirrors a gram-positive one, and the four figures are not four
independent readouts. `gram_sankey_summary.json` reports the agreement directly as
`complementarity`.

## Tests

A test needs at least `min_studies` (default 3) informative studies. The rule reads
support only and is fixed before testing.

- **Sign test per flow.** For each gram class and each condition -> body site flow with
  at least three non-tied studies, an exact two-sided binomial test of studies up
  against studies down. Family `gram_sankey_flow_sign:<gram class>` (registry
  `mt-supp-gram-sankey-flow-sign`): the tested flows of one gram class.
- **McNemar per condition and body site.** For each condition and body site with at
  least three studies whose two gram contrasts have opposite signs, an exact McNemar
  test (`guild_maps.exact_mcnemar_p`, the test the guild figure applies to aerobes
  against anaerobes). Family `gram_sankey_site_mcnemar` (registry
  `mt-supp-gram-sankey-site-mcnemar`).

Both families report BH q (primary, threshold q <= 0.05) and BY q beside it. Every row
of the edge table carries its family id, the selected family size, the unfiltered
inventory size, the procedure and the threshold (methods/multiple_testing_policy.md R1,
R6). Flows and sites below the support rule are written with blank p and q and are
outside the family.

## Figures

A flow is drawn once, in the direction most of its studies take, with width equal to
net studies (up minus down), so a six-to-five split is one study wide. Untested flows
and even splits are not drawn. A flow passing its sign test is drawn dark and the rest
faded. When no flow passes, the subtitle says so and states that the diagram shows
where study counts lean, not tested shifts.

## Results on the committed corpus

Run on 2026-11-19 against the committed `study_trait_profiles_wide.tsv.gz` (sha256
`a9c7300a…c18655`, BugSigDB 2026-04-24_00:41_UTC) and a local `experiments.tsv`, with
`min_studies` 3. These numbers describe that run; rerun the recipe for current values.

- 2,257 paired experiments; 43 removed for a zero total on a side; 923 studies give
  1,193 study-within-site units over 675 condition and body-site pairs. 165 studies
  appear in more than one pair.
- The two classes shift in opposite directions in 934 units, the same direction in
  246, and at least one is tied in 13.
- 85 flows per gram class meet the support rule. **No flow passes its sign test.** The
  smallest BH q is 0.32: taxa annotated as gram negative were reported increased in 22
  colorectal cancer feces studies and decreased in 6 (nominal p 0.0037). One of 170
  tested flows has a nominal p below 0.05.
- 63 condition and body-site pairs meet the McNemar support rule. **None passes.** The
  smallest BH q is 0.45 (colorectal cancer, feces: 4 studies gram positive up and
  gram negative down, 17 the reverse; nominal p 0.0072).

These results match the earlier, unmerged implementation (#148) in direction and
size, though that version ran before the corpus was narrowed to 16S and shotgun
experiments and its counts differ.
