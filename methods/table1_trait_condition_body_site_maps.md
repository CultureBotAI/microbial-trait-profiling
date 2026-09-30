# Manuscript Table 1 traits by condition and body site, and their paired directional shift

Figures: `table1_traits_condition_body_site_heatmap.png`, the four cuts
`table1_traits_focus_feces.png`, `table1_traits_focus_feces_other.png`,
`table1_traits_focus_cervical.png`, and `table1_traits_focus_oral.png`, and
`table1_traits_direction_shift.png`. Produced by
`src/bsdbng/table1_trait_maps.py`, driven by `just table1-trait-maps`. How this
analysis differs from the composition and enrichment arms is set out in
[`three_arm_statistics.md`](three_arm_statistics.md).

The figures described here repeat the oxygen guild analysis
([`guild_condition_body_site_heatmap.md`](guild_condition_body_site_heatmap.md))
over the six individual traits of manuscript Table 1 rather than over the two
pooled oxygen guilds, with all six traits combined into one figure each. The
traits are facultatively aerobic (METPO:1000608), cellulose catabolic process
(GO:0030245), coccus shaped (METPO:1000668), sphere shaped (METPO:1000683),
arsenate reduction (pathways:arsenate_reduction), and moderately halophilic
(METPO:1000623). Trait slugs, display names, and per-trait colours were taken
verbatim from six_shared_trait_selection.tsv on the composition-arm branch, so a
given trait is the same colour in both arms' figures.

## Input corpus

Inputs and corpus filters are the same code path as the guild figure and are not
re-derived here. BugSigDB differential-abundance signatures were ingested from
the full_dump.csv of the **pinned BugSigDBExports v1.3.1 release** (Zenodo
10.5281/zenodo.19740674; the dump's internal timestamp, 2026-04-24_00:41_UTC, is
stamped on every output row) into per-study YAML records, flattened to
per-experiment tables, and joined to KG-Microbe taxon-trait assertions to give a
study/experiment/direction by trait count matrix. Non-human samples are excluded
on the per-experiment host_species field (the v1.3.1 wide matrix is already
human-only upstream, so that filter drops 0 rows on this run), experiments were
restricted to those carrying both an increased and a decreased signature, and
growth-medium columns carrying a mediadive.medium marker were removed from the
trait space before any total was computed. Of 5,829 profile rows read, 1,159
experiments were dropped as unpaired; 2,290 of the trait columns were MediaDive,
leaving 317 traits, 2,335 experiments, 4,670 profiles, and 959 studies — the
same corpus as the composition arm's 2,335 pairs.

## Resolving the six traits to columns

Two of the six traits reach the profile matrix behind a METPO qualifier, as
METPO:2000103|GO:0030245 and METPO:2000103|pathways:arsenate_reduction. The
qualifier records how the term was attached rather than which term it is, so
trait lookup matches either the whole column identifier or the segment after the
final pipe. Because a term that resolved to two columns would be counted twice in
a trait sum and once in the denominator, the loader raises rather than proceeding
if any term matches more than one column. Across the full 2,610-column header of
the v1.3.1 matrix the only suffix shared by two columns is a MediaDive medium,
which is dropped before lookup runs.

## First figure: trait fraction by condition and body site

Each cell is a mean trait fraction. Within a profile, the trait's count was
divided by that profile's total count across the retained traits; profile
fractions were averaged within a study before being averaged across the studies
in a cell, and the companion count matrix reports contributing studies rather
than experiments. Rows are conditions and columns are body sites crossed with the
twelve trait-by-direction combinations, so each body site contributes the
six-by-two table unrolled into twelve adjacent columns, separated from the next
site by a rule. Cells supported by fewer than three studies were masked.

Three axis decisions differ from the guild figure, and all three follow from
fitting six traits rather than two guilds into one figure.

First, the conditions outside the top-60 cut are **excluded rather than pooled
into an Other row**. Such a row mixes hundreds of unrelated conditions, so no cell
in it reads as one condition, and at twelve columns per body site it is expensive
to draw something that cannot be interpreted. Excluding them keeps every drawn row
a single label at the cost of showing only part of the corpus, and the cost is
reported on every run: the condition axis keeps the top 60 labels and drops 391
labels carrying 41.1 per cent of profiles. Passing `--include-other-conditions`
restores the pooled row. Note that BugSigDB records study variables such as Diet,
Age, and Extraction protocol in the same Condition field as diseases, so the
condition axis is not purely a disease axis.

The body-site axis is treated differently and **does keep its Other block**, which
pools 144 labels carrying 21.8 per cent of profiles. A pooled body site still
reads as a coherent quantity, everywhere the corpus sampled outside the eight
commonest sites, in a way a pooled condition does not, and keeping it means the
figure accounts for the whole corpus along that axis. `--drop-other-body-sites`
removes it.

Second, a body site whose column block carries fewer than two drawn conditions is
dropped, because twelve columns to show a single row is not worth the width;
`--min-site-conditions` sets the threshold. On this run that removed Mouth, Oral
cavity, Subgingival dental plaque, and Throat, leaving Feces, Saliva, Uterine
cervix, Vagina, and Other.

Third, cells where a site was sampled for a condition but no taxon carried the
trait are **left blank by default rather than shaded dark grey**. The guild figure
shades them, because a guild is absent rarely enough that the shading marks an
exception. Four of the six Table 1 traits are rare enough that the shading would
cover most of the figure and read as the dominant signal. The state is still
recorded per cell in the companion status table either way, and `--mark-no-signal`
turns the shading back on. Because a condition kept solely by such cells now
renders blank, it is dropped before clustering: on this run 9 conditions were
dropped (Atopic asthma, Cervical cancer, Head and neck squamous cell carcinoma,
Helminthiasis, Metabolic toxicity, Phenylketonuria, Population, Psoriasis
vulgaris, Smoking behaviour measurement).

On this run the grid is 51 conditions by 60 columns; of 3,060 cells, 806 carry a
value, 387 fall below the study threshold, 163 are sampled with no signal for that
trait, and 1,704 were never sampled. Condition rows are ordered by hierarchical
clustering of their value vectors under Euclidean distance with Ward linkage, cut
into six groups separated by blank rows, with masked cells treated as zero for the
distance calculation only — the same procedure as the guild figure.

The guild figure draws its two guilds on one shared colour range, because a
guild's share is directly comparable to the other guild's share. The six Table 1
traits are not comparable in that way: their corpus-wide prevalence spans more
than an order of magnitude, and on a shared range the four rarer traits would be
uniformly white. Each trait is therefore drawn on its own colour range, each
colour bar carries its own maximum, and the figure states that magnitudes are
readable within a trait and not across traits. Separately, in each trait's own
range a single extreme cell was flattening the rest of the distribution to one
pale tint, so ramps are cut at the 98th percentile of that trait's observed cells
and cells above the cut saturate at the top colour.

## Second set of figures: four cuts of the same matrix

Sixty columns across one figure leaves each column too narrow to label legibly,
and the Feces block, which carries 43 of the 51 drawn conditions, sets the eye's
scale so that every other block reads as empty beside it. The combined heatmap
therefore answers what was measured where rather than showing a result. The same
matrix is cut into four smaller heatmaps, each holding conditions that were
sampled in the same places:

- `focus_feces`, the 43 conditions drawn in Feces, on that block alone.
- `focus_feces_other`, the same conditions plus those drawn only in the pooled
  Other block, 49 in total, over the two blocks. That is 49 of the combined
  figure's 51 conditions: cervical glandular intraepithelial neoplasia and oral
  squamous cell carcinoma carry no Feces or Other cell, so neither appears here.
- `focus_cervical`, cervical glandular intraepithelial neoplasia and human
  papilloma virus infection, over Uterine cervix, Vagina, and Other.
- `focus_oral`, asked for as the cluster of gingivitis, periodontitis, smoking
  behavior, and rheumatoid arthritis. On this corpus it draws two of those four.
  Gingivitis and periodontitis are both in the pinned v1.3.1 export — 2 and 7
  human paired studies respectively, 9 studies in all under those exact condition
  labels — but neither reaches the drawn matrix, for two different reasons given
  below. The figure names each with its own reason rather than drawing a smaller
  group under the same title.

That last cut is worth stating plainly, because it was proposed from a figure
built on the unpinned devel export. There, gingivitis, periodontitis, smoking
behavior, and rheumatoid arthritis clustered together, and rheumatoid arthritis
grouping with three oral conditions rather than with the other systemic diseases
was the observation the cut was meant to show. That observation cannot be checked
here, but the reason is the display rules, not the vintage, and the two conditions
fail different rules:

- **Periodontitis is cut by width.** Of 451 distinct conditions it ranks 62nd on
  20 profiles — the same profile count as the 60th-ranked condition, losing its
  place only to the deterministic alphabetical tiebreak in `_top_labels`. Raising
  `--max-conditions` to 62 draws it, which is measured, not inferred.
- **Gingivitis is cut by support, so no width draws it.** It ranks 212nd on 4
  profiles, and its two studies sit in different body sites, leaving at most one
  study per cell against `--min-studies 3`. Widening the cut to 250 or beyond
  admits the row and then drops it again as entirely blank. Only lowering
  `--min-studies`, which would weaken every other cell in the figure, would draw
  it.

So drawing all four together is not a matter of widening the cut. It needs a
corpus where gingivitis carries at least three studies in one body site.

Each cut fixes one axis and derives the other from the data: naming conditions
keeps the body sites those conditions were actually drawn in, naming body sites
keeps the conditions drawn in them, and rows or column blocks left carrying no
value are dropped either way. So a cut never asserts that a condition was studied
in a site where it was not, and never renders a blank band.

Every cell is the same size in every figure, the combined one included. A cut
holding only a few conditions is padded with blank rows rather than stretched to
fill the panel, so a cell in the two-row oral cut is the same square as a cell in
the 43-row Feces cut and the two can be compared by eye. The six colour bars lie
flat in a band under each cut rather than stacked down its side: stacked, their
height would set the figure's height, and a two-row cut would have to be padded
out with empty rows to give them room or let them collapse into hairlines.

The cuts are a redraw and not a re-analysis. They use the same cell values, the
same masking, the same row order and cluster ids from the single clustering of
the full matrix, and the same per-trait colour ranges, taken from the full matrix
rather than recomputed per cut, so a given cell is the same colour in every
figure it appears in.

What these figures do not yet do is state a condition's increased-versus-decreased
signature directly: they place the two directions side by side and leave the
comparison to the reader. The paired contrast that would state it is already
computed for the third figure, one value per study per trait, and could be drawn
per condition on a diverging scale; that is a change of quantity, not of layout,
and is not attempted here.

## Third figure: paired directional shift

There is no aerobe-against-anaerobe contrast to test here, so each trait is
tested against itself across the two directions. Every study contributes one value
per trait, the mean increased fraction minus the mean decreased fraction over that
study's paired experiments; collapsing to the study first is the correction that
issue #133 asks for, since 56.2 per cent of studies contribute more than one
experiment.
A study counts as shifted up for a trait when that difference is positive and down
when it is negative; a difference within 1e-12 of zero is treated as a tie, since
at that magnitude it is floating-point residue of summing fractions rather than a
direction (exactly one study crosses that line on this corpus). Ties are recorded
in the output table and excluded from the test rather than split between the
sides, which is why the denominator differs between traits. A trait whose every
study ties has no defined share: its share, interval, p and q are written blank
with status `all_tied`, it keeps its place in the six-trait family at p = 1, and
the figure draws no marker for it rather than a point at 0 that would read as
every study shifting down. Ties are common for the rarer traits — 343 studies for
arsenate reduction against 42 for facultatively aerobic — and almost all of them
are studies where the trait is absent from both sides. Each trait then gets a two-sided exact binomial test against a half-and-half
null with an exact binomial 95 per cent interval, and the six p-values are
Benjamini-Hochberg corrected together as one family.

One boundary case is scored rather than dropped: a profile whose taxa carry none
of the retained traits enters with fraction 0 instead of as missing. 24
experiments have such an empty decreased side against 11 with an empty increased
side, a small net push toward "up". A sensitivity run dropping every experiment
with an empty side (2,241 experiments, 939 studies) changes no conclusion — the
largest movement is facultatively aerobic's share, 0.582 to 0.581.

This figure reads the profile table directly rather than the heatmap matrix, so
it is unaffected by the axis cuts above: every study in the corpus contributes,
not only those on a drawn row.

Over 959 studies, five of the six traits sit on the increased side and one on
the decreased side. Facultatively aerobic shifts up in 534 of 917 scored studies
(share 0.582, q = 4.1e-06), moderately halophilic in 485 of 845 (0.574,
q = 5.8e-05), arsenate reduction in 346 of 616 (0.562, q = 0.0030), sphere
shaped in 467 of 844 (0.553, q = 0.0030), and coccus shaped in 499 of 915
(0.545, q = 0.0067). Cellulose catabolic process runs the other way, shifting up
in only 334 of 774 studies (0.432, q = 0.00031). These directions agree
trait-for-trait with the composition arm's paired t-test on the same matrix,
which also places cellulose catabolic process alone on the decreased side — the
same corpus and the same six traits, so this is a robustness check against
study-level dependence (issue #133), not an independent line of evidence. The
effect sizes are small: every trait sits within about eight points of 0.5, and
the significance comes from the number of studies rather than from the size of
the shift, which is worth stating explicitly next to the q-values.

The right panel reports descriptive Spearman rank correlations between traits
across those per-study contrasts. These still depend on the profile totals:
differencing and study-level averaging do not remove compositional effects.
Closure can produce positive as well as negative component correlations, so
neither sign alone establishes biological association. Coccus shaped and sphere
shaped correlate at rho = 0.45 over 959
studies, which is high enough to treat them as overlapping evidence rather than
two independent lines of it, but well short of the near-identity their names
suggest. Facultatively aerobic correlates positively with coccus shaped (0.45),
arsenate reduction (0.40), and moderately halophilic (0.38), and negatively with
cellulose catabolic process (-0.26); arsenate reduction and moderately halophilic
correlate at 0.34. A pair is omitted rather than written as a
missing value if either series is constant, since Spearman is undefined there.
The fifteen pairwise Spearman p-values are reported uncorrected, unlike the six
binomial p-values; they describe the redundancy between traits rather than
supporting a claim of their own.

The separate [paired log-ratio analysis](table1_paired_log_ratios.md) tests
reference-explicit relative contrasts for all 15 trait pairs. It reports zero
exclusions and the resulting conditional study populations. It does not replace
these historical fraction results or turn them into absolute abundance claims.

## Reproducing

All statistics use scipy (binomtest, false_discovery_control, spearmanr) rather
than hand-rolled implementations. The six figures and their six backing tables
are produced by `src/bsdbng/table1_trait_maps.py`, driven deterministically by
`just table1-trait-maps`, which shares its matrix construction, masking, row
clustering, and table writing with `src/bsdbng/guild_maps.py` rather than
duplicating them. The BugSigDB export timestamp is read from the dump's leading
comment line and recorded in every output table. PNG export runs through kaleido
1.x, which drives a headless Chrome. For a reproducible figure run, run
`just get-chrome` once first: it downloads (~170 MB) a version-pinned Chrome for
Testing build that kaleido then prefers over any system Chrome, which would
otherwise be whatever version the machine carries (plotly's `plotly_get_chrome`
is equivalent). Without any Chrome the recipe fails only at the figure-writing
step, with every backing table already written to disk; the corpus-cascade
report prints to stdout either way. None of the statistics depend on Chrome.

## Vintage

All numbers in this note are from the **pinned BugSigDBExports v1.3.1 release**
(Zenodo 10.5281/zenodo.19740674, downloaded with SHA verification by main's
`download_exports`), the same vintage both other arms pin (issues #115/#130,
resolved by PR #150). The dump's internal timestamp, 2026-04-24_00:41_UTC, is
recorded in every output table. An earlier draft of this note carried numbers
from an unpinned devel-branch export; those are superseded.

The pipeline code and this description were drafted with the assistance of a
coding agent (Claude Code) and reviewed by the author.
