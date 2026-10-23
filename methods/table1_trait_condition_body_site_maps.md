# Matched-corpus Table 1-trait analog by condition and body site, and paired directional shift

> [!IMPORTANT]
> **Matched-Corpus Table 1-Trait Analog (Issue #210)**
> The figures and numbers described here are a matched-corpus Table 1-trait analog sensitivity analysis, not a reproduction of the published Table 1 baseline. They use the pinned human paired KG-Microbe profile corpus (932 studies, 2,257 pairs under `count_input_16s_wms_v1`) and do not reproduce published Table 1 baseline values (e.g. facultatively aerobic $q = 2.8 \times 10^{-10}$).

Figures: `table1_traits_condition_body_site_heatmap.png`, the four cuts
`table1_traits_focus_feces.png`, `table1_traits_focus_feces_other.png`,
`table1_traits_focus_cervical.png`, and `table1_traits_focus_oral.png`, and
`table1_traits_direction_shift.png`. Produced by
`src/bsdbng/table1_trait_maps.py`, driven by `just table1-trait-maps`. How this
analysis differs from the composition and enrichment arms is set out in
[`three_arm_statistics.md`](three_arm_statistics.md).

The figures described here repeat the oxygen guild analysis
([`guild_condition_body_site_heatmap.md`](guild_condition_body_site_heatmap.md))
over the six individual traits named in manuscript Table 1 rather than over the two
pooled oxygen guilds, with all six traits combined into one figure each. The
traits are facultatively aerobic (METPO:1000608), cellulose catabolic process
(GO:0030245), coccus shaped (METPO:1000668), sphere shaped (METPO:1000683),
arsenate reduction (pathways:arsenate_reduction), and moderately halophilic
(METPO:1000623). Trait slugs, display names, and per-trait colours were taken
verbatim from six_shared_trait_selection.tsv on the composition-arm branch, so a
given trait is the same colour in both arms' figures.

## Input corpus

Inputs, corpus filters and population are the same code path as the guild figure
and are not re-derived here; the population rule ends this section.
BugSigDB differential-abundance signatures were ingested from
the full_dump.csv of the **pinned BugSigDBExports v1.3.1 release** (Zenodo
10.5281/zenodo.19740674; the dump's internal timestamp, 2026-04-24_00:41_UTC, is
stamped on every output row) into per-study YAML records, flattened to
per-experiment tables, and joined to KG-Microbe taxon-trait assertions to give a
study/experiment/direction by trait count matrix. Non-human samples are excluded
on the per-experiment host_species field, experiments were
restricted to those carrying both an increased and a decreased signature, and
growth-medium columns carrying a mediadive.medium marker were removed from the
trait space before any total was computed. Of 5,629 profile rows read under the
`count_input_16s_wms_v1` profile policy, from 1,047 studies, none was dropped as
non-human, because the matrix holds human profiles only (#656), and 1,115
experiments were dropped as unpaired; 2,290 trait columns were MediaDive,
leaving 317 traits, 2,257 experiments, 4,514 profiles, and 932 studies — the
same corpus as the composition arm's 2,257 pairs.

The matrix these counts were measured on is
`data/derived/study_trait_profiles_wide.tsv`, sha256
`519a095b09e179dc92fb51f9d892984a52dd100dece4dadf6c94b036a5323e6c`, the
committed `study_trait_profiles_wide.tsv.gz` decompressed: 2,610 header columns
(3 identifier, 2,290 MediaDive, 317 retained) over 5,629 profile rows, joined
from BugSigDB `2026-04-24_00:41_UTC` and KG-Microbe release 20260120 (nodes md5
`dd76b4656f9341db25a7c4afa96cd58e`), as the download sidecar
`data/kg-microbe/provenance.json` declares it. The experiment metadata table is
sha256 `16271237a70f61c8195d9f2acdd16f71e88f65859160f95bb734626bc1328e06` and the
full dump sha256 `af84bfca7bfe76913a20c1fedbaadbfad4e4bf78605a2cfd7d00eb0c0fdd5788`.
These are the guild figure's inputs too
([`guild_condition_body_site_heatmap.md`](guild_condition_body_site_heatmap.md)).
Naming the dump does not identify the matrix, so the run writes this record to
`corpus_identity.json` beside its outputs (#222).

Every figure here then runs on the shared comparison basis of
[`three_arm_statistics.md`](three_arm_statistics.md), as the guild figure does
(#783): an experiment pair whose retained-annotation total is zero on a side is
left out, both directions together. A profile with no retained annotation has
an undefined trait share, not a share of 0. Scored as 0, it lowered a cell's
mean, added a study to the cell's count, and could mark the cell as sampled
without signal. The paired directional shift has left these pairs out since
#228. The maps did not until #811: they now build their matrix with the guild
figure's `build_figure_matrix`, which drops the pairs whatever its caller passes.
Because the top-N cuts rank conditions and body sites by profile counts, the
exclusion can also change which rows and site blocks are drawn. The run prints
the number of pairs excluded, writes it as `excluded_experiment_pairs` on every
row of `table1_traits_direction_shift.tsv`, and states it in the subtitle of the
combined figure and of each cut. A cut names a requested condition that only
such pairs carried. It does not report that condition as absent from the export
or as cut by width. The map counts quoted below were measured with this
exclusion on the inputs identified above (matrix `519a095b…`, experiment table
`16271237…`), where it removes 55 pairs and leaves 922 studies. The final rerun
re-derives them from its own outputs (#908).


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
reported on every run: the condition axis keeps the top 60 labels and drops 369
labels carrying 40.4 per cent of profiles. Passing `--include-other-conditions`
restores the pooled row. Note that BugSigDB records study variables such as Diet,
Age, and Extraction protocol in the same Condition field as diseases, so the
condition axis is not purely a disease axis.

The body-site axis is treated differently and **does keep its Other block**, which
pools 139 labels carrying 21.1 per cent of profiles. A pooled body site still
reads as a coherent quantity, everywhere the corpus sampled outside the eight
commonest sites, in a way a pooled condition does not, and keeping it means the
figure accounts for the whole corpus along that axis. `--drop-other-body-sites`
removes it.

Second, a body site whose column block carries fewer than two drawn conditions is
dropped, because twelve columns to show a single row is not worth the width;
`--min-site-conditions` sets the threshold. On this run that removed Mouth, Oral
cavity, and Throat, leaving Feces, Saliva, Subgingival dental plaque, Uterine
cervix, Vagina, and Other. Before the zero-total exclusion it removed Subgingival
dental plaque too: periodontitis, which joins rheumatoid arthritis there, was
outside the condition cut (see the oral cut below).

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

On this run the grid is 51 conditions by 72 columns; of 3,672 cells, 814 carry a
value, 419 fall below the study threshold, 147 are sampled with no signal for that
trait, and 2,292 were never sampled. Before the exclusion it was 51 by 60, with
806, 399, 151 and 1,704 of 3,060 cells in those states. Condition rows are
ordered by hierarchical clustering of their value vectors under Euclidean
distance with Ward linkage, cut into six groups separated by blank rows, with
masked cells treated as zero for the distance calculation only — the same
procedure as the guild figure.

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

Seventy-two columns across one figure leave each column too narrow to label legibly,
and the Feces block, which carries 42 of the 51 drawn conditions, sets the eye's
scale so that every other block reads as empty beside it. The combined heatmap
therefore answers what was measured where rather than showing a result. The same
matrix is cut into four smaller heatmaps, each holding conditions that were
sampled in the same places:

- `focus_feces`, the 42 conditions drawn in Feces, on that block alone.
- `focus_feces_other`, the same conditions plus those drawn only in the pooled
  Other block, 48 in total, over the two blocks. That is 48 of the combined
  figure's 51 conditions: cervical glandular intraepithelial neoplasia, oral
  squamous cell carcinoma, and periodontitis carry no drawn Feces or Other cell,
  so none of them appears here.
- `focus_cervical`, cervical glandular intraepithelial neoplasia and human
  papilloma virus infection, over Uterine cervix, Vagina, and Other.
- `focus_oral`, asked for as the cluster of gingivitis, periodontitis, smoking
  behavior, and rheumatoid arthritis. On this corpus it draws three of those
  four, over Feces, Saliva, Subgingival dental plaque, and Other; before the
  zero-total exclusion it drew two, without periodontitis. Gingivitis is in the
  pinned v1.3.1 export — 2 human paired studies under that exact condition label
  — but does not reach the drawn matrix. The figure names it as outside the
  top-60 cut rather than drawing a smaller group under the same title, and the
  reason no width would draw it is given below.

That last cut is worth stating plainly, because it was proposed from a figure
built on the unpinned devel export. There, gingivitis, periodontitis, smoking
behavior, and rheumatoid arthritis clustered together, and rheumatoid arthritis
grouping with three oral conditions rather than with the other systemic diseases
was the observation the cut was meant to show. That observation cannot be checked
in full here, because gingivitis is not drawn. The reason is the display rules,
not the vintage, and the two conditions meet different rules:

- **Periodontitis sits at the edge of the width cut.** Of 429 distinct conditions
  it ranks 58th on 20 profiles, so it is drawn, and it joins rheumatoid arthritis
  in Subgingival dental plaque, which keeps that block under
  `--min-site-conditions 2`. The zero-total exclusion keeps all 20 of its
  profiles but takes three conditions ranked above it down to its count or
  below: hypertension to 18 profiles, non-alcoholic fatty liver disease to 16,
  and response to transplant to 20, which then loses the deterministic
  alphabetical tiebreak in `_top_labels`. Before the exclusion it ranked 61st of
  437 on the same 20 profiles, level with the 60th-ranked condition and behind it
  only on that tiebreak, and raising `--max-conditions` to 61 drew it and the
  same Subgingival dental plaque block.
- **Gingivitis is cut by support, so no width draws it.** It ranks 205th on 4
  profiles (206th before the exclusion), and its two studies sit in different
  body sites, leaving at most one study per cell against `--min-studies 3`.
  Widening the cut to admit it only drops the row again as entirely blank. Only
  lowering `--min-studies`, which would weaken every other cell in the figure,
  would draw it.

So drawing all four together is not a matter of widening the cut. It needs a
corpus where gingivitis carries at least three studies in one body site.

**Flagged for the rerun (#908).** This account, like the other map counts in
this note, was measured on the committed corpus before the final rerun. The
rerun must re-derive it from its own outputs rather than carry it over: it builds
its matrix under #828's strain rollup default, and periodontitis sits close
enough to the width cut that another matrix can move it back out, taking the
Subgingival dental plaque block and the cuts' contents above with it.

With `--include-other-conditions`, conditions outside the width cut are
assigned to a pooled Other row. The panel note uses the build's recorded
membership to distinguish this from omission and reports whether that row is
drawn in the particular panel. Assignment is not contribution: when a pooled
condition's body site is binned out, or its Other cells are masked, none of its
profiles reach a drawn Other cell, so the note says "assigned to" rather than
"included in". A masked pool or a condition-focused cut that excludes Other
does not display those conditions. Repeated requested labels are reported once,
in their first-request order.

Each cut fixes one axis and derives the other from the data: naming conditions
keeps the body sites those conditions were actually drawn in, naming body sites
keeps the conditions drawn in them, and rows or column blocks left carrying no
value are dropped either way. So a cut never asserts that a condition was studied
in a site where it was not, and never renders a blank band.

Every cell is the same size in every figure, the combined one included. A cut
holding only a few conditions is padded with blank rows rather than stretched to
fill the panel, so a cell in the three-row oral cut is the same square as a cell
in the 42-row Feces cut and the two can be compared by eye. The six colour bars lie
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
issue #133 asks for, since 511 of the 922 studies (55.4 per cent) contribute
more than one experiment.
A study counts as shifted up for a trait when that difference is positive and down
when it is negative; a difference within 1e-12 of zero is treated as a tie, since
at that magnitude it is floating-point residue of summing fractions rather than a
direction (exactly one study crosses that line on this corpus). Ties are recorded
in the output table and excluded from the test rather than split between the
sides, which is why the denominator differs between traits. A trait whose every
study ties has no defined share: its share, interval, p and q are written blank
with status `all_tied`, it keeps its place in the six-trait family at p = 1, and
the figure draws no marker for it rather than a point at 0 that would read as
every study shifting down. Ties are common for the rarer traits — 316 studies for
arsenate reduction against 24 for facultatively aerobic — and almost all of them
are studies where the trait is absent from both sides. Each trait then gets a
two-sided exact binomial test against a half-and-half null with an exact
binomial 95 per cent interval, and the six p-values are Benjamini-Hochberg
corrected together as one family.

One boundary case is dropped rather than scored: a profile with no retained
annotation has no defined fraction, so its whole experiment pair is left out
(#228), as it is from the maps. Of the 55 such pairs, 24 have an empty decreased
side only, 10 an empty increased side only and 21 both. Scored as 0, the first
would push their studies toward "up", the second toward "down", and the third
would add zero contrasts: a small net push toward "up". Dropping them leaves
2,202 experiments and 922 studies. Scoring them as 0 instead, as this test did
before #228, keeps all 932 studies and changes no conclusion — the largest
movement is moderately halophilic's share, 0.576 with them against 0.574
without.

This figure reads the profile table directly rather than the heatmap matrix, so
it is unaffected by the axis cuts above: every study on the shared basis
contributes, not only those on a drawn row.

Over 922 studies, five of the six traits sit on the increased side and one on
the decreased side. Facultatively aerobic shifts up in 521 of 898 scored studies
(share 0.580, q = 1.0e-05), moderately halophilic in 476 of 829 (0.574,
q = 6.6e-05), arsenate reduction in 341 of 606 (0.563, q = 0.0027), sphere
shaped in 463 of 833 (0.556, q = 0.0021), and coccus shaped in 488 of 896
(0.545, q = 0.0083). Cellulose catabolic process runs the other way, shifting up
in only 330 of 763 studies (0.433, q = 0.00043). These directions agree
trait-for-trait with the composition arm's paired t-test on the same matrix,
which also places cellulose catabolic process alone on the decreased side — the
same corpus and the same six traits, so this is a robustness check against
study-level dependence (issue #133), not an independent line of evidence and
not a reproduction of the published Table 1 baseline. The
effect sizes are small: every trait sits within about eight points of 0.5, and
the significance comes from the number of studies rather than from the size of
the shift, which is worth stating explicitly next to the q-values.

The right panel reports descriptive Spearman rank correlations between traits
across those per-study contrasts. These still depend on the profile totals:
differencing and study-level averaging do not remove compositional effects.
Closure can produce positive as well as negative component correlations, so
neither sign alone establishes biological association. Coccus shaped and sphere
shaped correlate at rho = 0.44 over 922
studies, which is high enough to treat them as overlapping evidence rather than
two independent lines of it, but well short of the near-identity their names
suggest. Facultatively aerobic correlates positively with coccus shaped (0.46),
arsenate reduction (0.39), and moderately halophilic (0.38), and negatively with
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

### Outstanding: the correction family is described here but not exported

`table1_traits_direction_shift.tsv` carries each trait's p and q and its status,
but nothing in the row says which correction family produced that q or how large
the family was. A reader joining this table to the composition arm's per-stratum
BH tables or to the enrichment arm's global FDR sees three columns called `q`
whose families differ in size by orders of magnitude — six selected traits here
against hundreds of screened traits and hundreds of thousands of profile tests
there, per #178 — and no exported row says which. The family here is Table 1's
six traits corrected together, with every trait included at p = 1 when it has no
scored study, so the family size is 6 on every run of this figure and the
untestable traits are inside it rather than dropped from it.

That is a property of the exported row, not of this figure, and the same gap
exists in the other two arms' exports. It is therefore tracked as part of the
shared result contract in
[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178) and
implemented in
[#241](https://github.com/CultureBotAI/microbial-trait-profiling/pull/241),
whose acceptance criteria already require a family identifier and its size on
every exported result. It is deliberately not added here: adding a
`correction_family` / `correction_family_size` pair to this one TSV ahead of
that contract would fix a column spelling that the integration then has to
reconcile across three producers. When #241 lands the implementing commit, this
paragraph should name it and the requirement closes.

## Reproducing

All statistics use scipy (binomtest, false_discovery_control, spearmanr) rather
than hand-rolled implementations. The six figures and their six backing tables
are produced by `src/bsdbng/table1_trait_maps.py`, driven deterministically by
`just table1-trait-maps`, which shares its zero-total pair exclusion, matrix
construction, masking, row clustering, and table writing with
`src/bsdbng/guild_maps.py` rather than duplicating them. The BugSigDB export
timestamp is read from the dump's leading comment line and recorded in every
output table. PNG export runs through kaleido
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
