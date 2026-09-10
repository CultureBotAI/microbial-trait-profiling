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
  Other block, 49 in total, over the two blocks. This is the cut that carries
  every condition the combined figure draws.
- `focus_cervical`, cervical glandular intraepithelial neoplasia and human
  papilloma virus infection, over Uterine cervix, Vagina, and Other.
- `focus_oral`, asked for as the cluster of gingivitis, periodontitis, smoking
  behavior, and rheumatoid arthritis. On this corpus it draws two of those four:
  the pinned v1.3.1 export contains no gingivitis or periodontitis study, so
  those rows do not exist to cut. The figure names them as absent rather than
  drawing a smaller group under the same title.

That last cut is worth stating plainly, because it was proposed from a figure
built on the unpinned devel export. There, gingivitis, periodontitis, smoking
behavior, and rheumatoid arthritis clustered together, and rheumatoid arthritis
grouping with three oral conditions rather than with the other systemic diseases
was the observation the cut was meant to show. It is an artefact of vintage, not
a result: v1.3.1 is a reviewed-only export and the oral studies are not in it.
Any reading of that cluster has to wait for a corpus that contains it.

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
sides, which is why the denominator differs between traits: ties are common for
the rarer traits, 343 studies for arsenate reduction against 42 for facultatively
aerobic, and almost all of them are studies where the trait is absent from both
sides. Each trait then gets a two-sided exact binomial test against a half-and-half
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

The right panel reports Spearman rank correlations between traits across those
per-study contrasts. This is deliberately not the same quantity as the guild
figure's within-profile correlation, which compares two shares of one profile
total and is therefore pushed negative by compositionality; correlating one
trait's contrast against another trait's contrast across studies does not share
that artefact. Coccus shaped and sphere shaped correlate at rho = 0.45 over 959
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

## Third layer: outlining the cells the other two arms call significant

The heatmap says how large a trait's share of a profile's signal is. It does not
say whether that share is more than chance, and this project already has two arms
that answer exactly that question in two different ways. At Marcin Joachimiak's
request (issue #172), both were run over the same corpus the figure is drawn from
and their verdicts drawn on top of it, so one figure carries all three approaches.

Neither test is invented here. The **enrichment arm** (Analysis B, Ilyas
Abdolcader, `scripts/trait_enrichment.py` on
`origin/enrichment-analysis-heatmap-cluster`) asks, one profile at a time, whether
that profile drew more assertions of a trait than a draw of its size from the
corpus pool would give: a one-tailed hypergeometric test, Benjamini-Hochberg
across every test at once, and a call at `q <= 0.001` with a fold of at least 1.5.
The **composition arm** (Analysis A, Jed Kim-Ozaeta,
`scripts/create_trait_profile_infographic.py` on
`origin/clustered-heatmap-generation`) scales trait shares per 1,000,
`log1p`-transforms them, z-scores each trait within the stratum, and tests the
increased side of a matched pair minus its decreased side against zero with a
one-sample t-test. Both are transcribed into `src/bsdbng/trait_significance.py`
with their own thresholds.

Marcin's prerequisite was to check whether either arm's committed output already
resolves significance at the level of these cells. The answer differs by arm, and
it is what decides how each one is drawn.

The enrichment arm's `trait_enrichment_results.tsv` is keyed
`(study_id, experiment_id, inc_dec, trait)`, which aggregates directly onto a cell.
Its header records `Total-Profiles-Input: 4670` and
`Total-Tests-Performed: 556516`; both match this corpus exactly, so the branch tip
already applies the human-only, paired-only, and MediaDive-excluded filters of
issues #121, #122, and #129. All six Table 1 traits are present, the two non-METPO
ones under their METPO-qualified names (`METPO:2000103|GO:0030245`,
`METPO:2000103|pathways:arsenate_reduction`), which is why a lookup on the bare
`GO:0030245` finds nothing. Re-running the test here reproduces that table's
per-trait significant-profile counts exactly - 860, 660, 864, 515, 267, and 262
for the six traits, over the same 556,516 tests - which is asserted by
`tests/test_trait_significance.py::test_reproduces_the_enrichment_arm_committed_counts`
(marked `integration`, since it reads the other branch).

The composition arm's `trait_profile_paired_contrast_trait_stats_<site>.tsv` is
keyed `(body_site_group, trait)` and pools every condition in the site, so it
cannot say which row of a block is significant. Running its test per
`(condition, body site, trait)` on this corpus shows why it was never emitted that
way: only 413 of the 1,836 trait cells reach the three matched study pairs the
test needs, and after Benjamini-Hochberg exactly one survives. That is a finding
about the corpus, not a defect in the arm - a single condition in a single body
site rarely carries enough paired studies to support a t-test.

So the two arms are drawn differently, at the grain each one actually resolves.

- **A black outline on a cell** is the enrichment arm: most of the studies behind
  that cell were individually called significant. A cell needs at least two
  studies and a strict majority of them, so one study never speaks for a cell, and
  a cell the figure masked for want of studies is never outlined - an outline on
  blank paper would read as a result rather than as a withheld cell. 62 cells are
  outlined.
- **A triangle on a column label** is the composition arm: that trait's paired
  contrast is significant across the whole body site, with the triangle pointing
  at the side carrying more of it. It says nothing about any one condition in the
  block, and the figure says so.

Two departures from the arms' committed runs are deliberate. Pairs are averaged to
one value per study before the composition t-test, which is the correction issue
#133 asks for, and the enrichment background pool is built from the figure's own
filtered corpus rather than an unfiltered one, so the test scores a profile
against the pool it belongs to. Both are stated here rather than silently applied.

The per-cell verdicts and the paired contrasts are written to
`table1_traits_cell_significance.tsv` and
`table1_traits_composition_contrasts.tsv`, so every outline in the figure can be
traced to the studies behind it.

### What the outlines show

Pooled over conditions, **all six traits shift significantly between the increased
and decreased sides of a fecal signature**, and five of the six shift the same
way. Facultatively aerobic (mean contrast +0.311, q = 2.9e-06 over 394 studies),
arsenate reduction (+0.291, q = 8.6e-05), sphere shaped (+0.284, q = 1.0e-05),
moderately halophilic (+0.251, q = 6.2e-06), and coccus shaped (+0.248,
q = 2.2e-04) all sit higher among the taxa a study reported as increased.
Cellulose catabolic process is the exception and runs the other way (-0.249,
q = 5.7e-04): it sits higher among the taxa reported as decreased. Read together
that is the oxygen hypothesis stated in one line - the organisms that rise in a
dysbiotic gut are the oxygen-tolerant ones, and the ones that fall are the fibre
degraders. Outside the gut only one contrast survives, coccus shaped in the pooled
Other block (+0.274, q = 0.014).

The cell-level outlines concentrate in the same two traits: coccus shaped
increased (17 cells) and facultatively aerobic increased (14) account for half of
the 62. Several conditions carry a coherent set of them. In feces, COVID-19 is
outlined for facultatively aerobic and coccus shaped on the increased side and
cellulose catabolic process on the decreased side, on 13 studies; gastric cancer
carries the same three-trait pattern on 4 studies and repeats facultatively
aerobic increased in saliva; diarrhea is outlined for five of the six traits, all
on the increased side except cellulose catabolism; stroke is outlined for
facultatively aerobic, coccus shaped, and arsenate reduction on the increased
side. Rheumatoid arthritis carries coccus shaped and sphere shaped on the
increased side in both feces and saliva, which is the only trait pattern that
crosses those two sites in the same condition.

Two results run against the pooled pattern and are worth keeping in view rather
than smoothing over. Chronic constipation is outlined for cellulose catabolic
process on the **increased** side over 14 studies, the opposite of the fecal
contrast, which is what slow transit favouring fibre degraders would look like.
And in the cervical conditions the direction inverts: cervical glandular
intraepithelial neoplasia is outlined for coccus shaped on the decreased side in
both uterine cervix and vagina. A vaginal community is rod-dominated at baseline,
so "cocci among the taxa that fell" is not the same claim as "cocci among the taxa
that rose" in a gut, and these two rows should not be read on the gut's terms.

Marcin's own reading of the earlier cuts was that a trait is rarely called in both
directions at once. That holds and is now countable: of the 57 condition-by-body-
site-by-trait combinations with any call, 5 carry one in both directions.

One caveat limits every outline. The enrichment arm's urn still draws trait
assertions rather than taxa, which is what issue #131 measures and expects to
shrink the significant set substantially once fixed. An outline here therefore
means "unusually assertion-rich in that trait", not "unusually taxon-rich", and
the same rollup that inflates a cell's value inflates its significance. The
outlines are a cross-arm integration of what the two arms currently say, not a
claim that survives that fix.
