# Selected-figure deck v2

The October 2 revision contains 31 slides and 30 original figure images.
It retains all 21 of Jed's Sankeys, including the combined six-color views on
slides 5–7. The main sequence ends on slide 13; the 18 alternate and
individual-trait Sankeys occupy slides 14–31. Work continues under issue #282
and PR #283.

## October 2 changes

Jed requested the expanded PD influence figure on slide 12 and removal of the
separate contributor view on slide 14. Slide 12 now uses
`pd_exploratory_taxa_six_traits_influence_heatmap.png`, pinned to
`205257c5ec887bdf145e2c66e12f17ee005fafa0`. It displays 30 species, genus and
family rows, including 23 matches to the curated PD literature reference.
The source tables audit all 180 taxon–trait cells: one has an
increased-direction call, six have decreased-direction calls, 167 are mixed
and six have low support. These are exploratory, non-additive removal
diagnostics, not causal effects or a calibrated discovery family.

The expanded figure preserves the literature comparison and makes the rank
difference visible: E. coli recurs in one increased and two decreased PD
profiles, while the genus Escherichia recurs in nine increased and two
decreased profiles. The prior genus-group result is linked to
[Li et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/28536926/),
DOI 10.1007/s11427-016-9001-4. Overlapping rank roll-ups, different estimands
and sparse species recurrence limit the comparison. The previous species-only
figure's statement that all 30 intervals include zero does not describe this
expanded figure and has been replaced.

The previous v2 slides 15–19 were Ilyas's saved enrichment/depletion
infographics, clustered heatmap and two guild Sankeys, pinned to
`0a0d739a70c7356ad0285154601af52a83316c94`. They are removed pending validated
regeneration. At the reviewed ORA head `e71f967`, background/sampling-unit,
complete-export and calibration findings remain open in #131, #254 and #263;
the command/export/check blockers are #295–#297. This deck revision does not
rerun ORA or replace those outputs with unsupported new results.

The six removals use the previous 37-slide v2 numbering: 14, 15, 16, 17, 18
and 19. Slides after those removals are renumbered. The complete 43-slide v1
sources in the parent directory retain the historical figures.

## Selection and provenance

The earlier v2 selection removed original v1 pages 4, 5, 11, 14, 15 and 16,
which centered on null, inconclusive or testing-family sensitivity findings.
Decreased traits and under-representation were not themselves exclusion
criteria. The cervical descriptive enrichment and study-level PD halophily
pages remain with their existing qualifications.

`page_selection.tsv` maps all 43 v1 pages to their current retention,
replacement or removal decisions. `slide_manifest.tsv` records every retained
figure's hash, source revision, evidence and prior-result reference.
`sankey_coverage.tsv` records all 21 Jed views with the current page numbers.
All retained PNGs use their original pixels.

## Rebuild and validation

Use the existing repository recipes from the root of the issue-linked branch:

```bash
deck_dir="$(pwd)/data/derived/selected-figure-deck-v2"
mkdir -p "$deck_dir"
cp docs/selected-figure-deck/v2/selected_figures.html "$deck_dir/"
just review-deck-assets docs/selected-figure-deck/v2/slide_manifest.tsv "$deck_dir" \
  /absolute/path/to/original/guild_condition_body_site_heatmap.png
just review-deck-pdf "file://$deck_dir/selected_figures.html" \
  "$deck_dir/selected_figures.pdf" /absolute/path/to/chrome-headless-shell
```

The unversioned original guild image remains an explicit input, identified by
its manifest hash. Before delivery, verify 31 pages, 30 native-resolution
figure images, all image hashes, source footers, numbering, all 21 Jed
Sankeys, the expanded slide 12 and absence of the six newly removed pages.
Inspect the rendered layout, particularly slide 12 and the appendix transition.
Run `just check` and `just precommit` before committing.

The local PDF is delivered beside v1 as `microbla-trait-profiling_v2.pdf`.
The preceding v2 is backed up before replacement; v1 is preserved.
