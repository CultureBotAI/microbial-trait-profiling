# Selected-figure deck

The September 23 refresh replaces all ten of Jed's figure images in the 22-slide
deck from issue #274. The other eleven figures and their slide content are
preserved. Work is tracked in #282.

Jed's sources are pinned to `6db7c3df16342500fe677f12c55300e251b1b1d0` on
`clustered-heatmap-generation`. The artifact-only copy briefly merged in #281 was
reverted on `main`; this deck deliberately follows Jed's current analysis branch.

`selected_figures.html` is the print source. `slide_manifest.tsv` records source
revisions, paths, captions, methods links, and image SHA-256 hashes.
`replacements.tsv` maps every replaced image to its predecessor. Images and PDFs
are generated local artifacts and are not duplicated in this directory.

## Changes

| Slides | Replacement |
| --- | --- |
| 2, 4, 5, 11 | Current individual-condition heatmaps and matching q values |
| 3 | Complete grouped heatmap, including the new dendrogram trait-pattern labels |
| 6 | Current branch's pooled paired-contrast infographic |
| 10 | Combined ecological-reporter condition comparison and taxon provenance |
| 15 | Current strict species-influence heatmap |
| 16 | Condition–direction–taxon–trait evidence network |
| 17 | Contributor map with direction differences, unresolved mass, and normalized specificity |

The overview now reports the current branch's four selected cells and the
2,379-row individual-condition BH family. These values differ from the prior
deck's denominator-harmonized branch. The source tables and their own testing
families are preserved; no statistical rerun or cross-arm evidence combination
is performed. The network and combined figure have committed methods and audit
tables, but no generating script was found in the branch's source, scripts,
tests, or Justfile, including an ignored-file search of the checkout. Their
footers therefore link the methods and state that the generator is not committed.

The source for slide 9 is the pre-existing, unversioned local oxygen-guild PNG.
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

Before publishing a PDF, compare all 21 PNG hashes with the manifest, confirm
22 pages and 21 native-resolution embedded images, check every page's source
footer, and visually inspect the changed slides. Run `just check` before a
commit. The original local deck is backed up before its replacement.
