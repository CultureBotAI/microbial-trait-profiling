# Selected-figure deck v2

Version 2 contains 37 slides and 36 original figure images. It retains all 23
Sankeys, including the combined six-color views on slides 5–7. The main sequence
ends on slide 19; the 18 alternate and individual-trait Sankeys occupy slides
20–37. The complete 43-slide v1 sources in the parent directory are unchanged.
Work continues under issue #282 and PR #283.

## Selection

The user requested removal of pages about negative results, with exceptions for
documented controversies or disagreement with prior results. Negative here means
null, inconclusive, or sensitivity-focused findings; it does not mean a trait
decreases, is under-represented, or occurs in decreased taxon signatures.

Removed original v1 pages:

| v1 slide | Reason |
| --- | --- |
| 4 | PD facultative aerobiosis loses significance in the full screening family |
| 5 | COVID-19 cellulose loses significance in the full screening family |
| 11 | No adjusted fecal condition–trait association in the broad screen |
| 14 | Oral coccus contrast remains inconclusive after correction |
| 15 | Zero adjusted associations across 18,152 condition–site tests |
| 16 | All six salivary results are nonsignificant |

The testing-family comparisons on v1 slides 4 and 5 are sensitivity analyses,
without a documented disagreement with an independent prior study. They do not
use the requested exception.

V1 slide 18 is retained as v2 slide 12 under the prior-result exception. The
figure explicitly compares sparse E. coli recurrence (one increased versus two
decreased PD profiles) with the increased Escherichia–Shigella finding reported
by [Li et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/28536926/), DOI
10.1007/s11427-016-9001-4. The title now states that comparison, and the footer
links the primary report. The notes preserve the different taxonomic ranks,
sparse recurrence, and all 30 influence intervals crossing zero. This is a
descriptive discrepancy to investigate, not a demonstrated contradiction.

The cervical slide remains because it reports descriptive profile-enrichment
support. Its title now names the coccus pattern in decreased cervical-neoplasia
signatures. The study-level PD slide remains because halophily passes its
336-test correction family. Necessary limits remain in the retained slides'
notes; source PNGs and their embedded annotations are unchanged.

`page_selection.tsv` maps all 43 v1 pages to retention/removal decisions and
new page numbers. `slide_manifest.tsv` records each retained figure's hash,
source revision, evidence, and any prior-result reference. `sankey_coverage.tsv`
keeps the complete 21-view Jed audit with the new page numbering; the two Ilyas
Sankeys are also present in the slide manifest.

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

The unversioned original guild image is the same input as in v1. All pinned
source revisions and PNG bytes are retained from v1. Before delivery, check 37
pages, 36 native-resolution embedded images, all retained image hashes, all
titles/notes/source footers, page numbering, the prior-report citation, and
absence of the six removed pages. Check for layout overflow and inspect the
revised pages. Run `just check` and `just precommit` before committing.

The local PDF is delivered beside v1 as `microbla-trait-profiling_v2.pdf`, with
separate v2 companion files. V1 is preserved byte-for-byte.
