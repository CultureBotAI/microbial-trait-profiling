# Selected-figure deck

The September 24 refresh contains 43 slides and 42 source figures. It includes
all 21 Sankeys in Jed's current manifest and the two published Ilyas Sankeys.
The three combined six-color views appear on slides 7–9; all 18 alternate and
individual-trait views appear in the appendix, slides 26–43. Every presentation
title now states a biological finding or a specific data point. Work is tracked
in #282 and PR #283, continuing the deck from #274.

Jed's sources are pinned to `205257c5ec887bdf145e2c66e12f17ee005fafa0` on
`clustered-heatmap-generation`. The artifact-only copy briefly merged in #281 was
reverted on `main`; this deck deliberately follows Jed's current analysis branch.

`selected_figures.html` is the print source. `slide_manifest.tsv` records source
revisions, paths, captions, evidence links, and image SHA-256 hashes.
`replacements.tsv` maps the September 24 slides to their September 23
predecessors. `sankey_coverage.tsv` records all 21 added views, their source TSVs,
link counts, and displayed positive-z profile membership totals. Images and
PDFs are generated local artifacts and are not duplicated in this directory.

## Changes

| Slides | Replacement |
| --- | --- |
| 1–6 | Exploratory cellulose leads, grouped associations, testing-family sensitivity, and pooled contrasts |
| 7–9 | Combined bidirectional and separate increased/decreased Sankeys, six trait colors per selected pair and direction |
| 10–17 | Study-level context, 13 selected ecological-reporter conditions, and anatomy-specific evidence limits |
| 18 | Five retained PD species; nine displayed effects, all 30 audited intervals include zero |
| 19 | Current evidence landscape: PD halophily alone survives 336 study-level trait tests |
| 20 | Eleven PD contributor taxa retained by a five-percentage-point share difference |
| 21–25 | Original Ilyas enrichment figures and both saved guild Sankeys, with revised titles |
| 26–43 | Full overlays, six-lane layouts, individual-condition maps, and twelve single-trait Sankeys |

The overview reports the source branch's four exploratory profile-level leads
and the 2,379-row individual-condition BH family. The grouped family has 747
rows. These differ from the study-level analysis on slide 19: 56 fecal
conditions × six traits, with global BH correction and one surviving trait
association (PD halophily, q = 0.04798; 45 paired studies). Its separate
264-test taxon family has no survivors. Both updated PD generators are now
committed and linked in the figure footers. Source test families are preserved;
no statistical rerun or cross-analysis evidence combination is performed.

Sankey coverage was checked against every fetched remote branch's published
PNG paths, not only the working tree. All 23 unique published Sankey PNG paths
are represented. Jed's three reduced views have exactly six trait-color links
for every selected condition/site/direction: 72 links in each direction and 144
in the combined view, spanning a union of 14 condition/site pairs. The six
colors identify facultative aerobiosis, cellulose catabolism, coccus shape,
sphere shape, arsenate reduction, and moderate halophily. Solid and dashed
links distinguish increased and decreased profiles in the combined map.

Appendix percentages divide fecal positive-z memberships by all memberships
across that figure's displayed links. These are selected descriptive views;
memberships can repeat across traits and overlapping condition groups. They
are not counts of independent patients, studies, or unique profiles. The oval
counts in reduced views deduplicate across traits within each selected
condition/site/direction. All original PNG pixels and their legends are kept.

The source for slide 12 is the pre-existing, unversioned local oxygen-guild PNG.
Its hash is recorded in the manifest; supply that original image explicitly
when reconstructing the assets. Do not substitute a newly generated figure.

## Rebuild

Fetch the source revisions listed in the manifest before reconstructing assets.
The shell recipe reads the manifest's existing figure/revision/URL/asset columns.
The following commands run from the repository root; replace the two explicit
input paths with the original guild image and an installed dedicated headless
shell.

```bash
deck_dir="$(pwd)/data/derived/selected-figure-deck"
mkdir -p "$deck_dir"
cp docs/selected-figure-deck/selected_figures.html "$deck_dir/"
just review-deck-assets docs/selected-figure-deck/slide_manifest.tsv "$deck_dir" \
  /absolute/path/to/original/guild_condition_body_site_heatmap.png
just review-deck-pdf "file://$deck_dir/selected_figures.html" \
  "$deck_dir/selected_figures.pdf" /absolute/path/to/chrome-headless-shell
```

The dedicated headless shell is used because desktop Chrome's failure to exit
after PDF export is documented in #275. The renderer must return normally with
exit code zero. The recipes restore the earlier local-only build workflow as
versioned commands; they do not suppress renderer failures or install tools.

Before publishing a PDF, compare all 42 PNG hashes with the manifest, confirm
43 pages and 42 native-resolution embedded images, check every page's source
footer, and visually inspect the changed slides. Verify all 21 source Sankey
manifest entries and the two Ilyas paths are represented. Run `just check`
before a commit. The original local deck is backed up before its replacement.
