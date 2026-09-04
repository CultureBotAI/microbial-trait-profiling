Manuscript Table 1 traits by condition and body site, and their paired
directional shift

The pair of figures described here repeats the oxygen guild analysis
(methods/guild_condition_body_site_heatmap.txt) over the six individual traits of
manuscript Table 1 rather than over the two pooled oxygen guilds, with all six
traits combined into one figure each. The traits are facultatively aerobic
(METPO:1000608), cellulose catabolic process (GO:0030245), coccus shaped
(METPO:1000668), sphere shaped (METPO:1000683), arsenate reduction
(pathways:arsenate_reduction), and moderately halophilic (METPO:1000623). Trait
slugs, display names, and per-trait colours were taken verbatim from
six_shared_trait_selection.tsv on the composition-arm branch, so a given trait is
the same colour in both arms' figures.

Inputs and corpus filters are the same code path as the guild figure and are not
re-derived here. BugSigDB differential-abundance signatures were ingested from
the full_dump.csv export (BugSigDB 2026-04-18_11:12_UTC) into per-study YAML
records, flattened to per-experiment tables, and joined to KG-Microbe taxon-trait
assertions to give a study/experiment/direction by trait count matrix. Non-human
samples were excluded on the per-experiment host_species field, experiments were
restricted to those carrying both an increased and a decreased signature, and
growth-medium columns carrying a mediadive.medium marker were removed from the
trait space before any total was computed. Of 12,122 profile rows read, 2,764
were dropped as non-human and 1,886 as unpaired; 2,334 of 2,651 trait columns
were MediaDive, leaving 317 traits, 3,736 experiments, and 7,472 profiles.

Two of the six traits reach the profile matrix behind a METPO qualifier, as
METPO:2000103|GO:0030245 and METPO:2000103|pathways:arsenate_reduction. The
qualifier records how the term was attached rather than which term it is, so
trait lookup now matches either the whole column identifier or the segment after
the final pipe. Because a term that resolved to two columns would be counted
twice in a trait sum and once in the denominator, the loader raises rather than
proceeding if any term matches more than one column. Across the full 2,654-column
header the only suffix shared by two columns is a MediaDive medium, which is
dropped before lookup runs.

First figure: trait fraction by condition and body site. Each cell is a mean
trait fraction. Within a profile, the trait's count was divided by that profile's
total count across the retained traits; profile fractions were averaged within a
study before being averaged across the studies in a cell, and the companion count
matrix reports contributing studies rather than experiments. Rows are conditions
and columns are body sites crossed with the twelve trait-by-direction
combinations, so each body site contributes the six-by-two table unrolled into
twelve adjacent columns, separated from the next site by a rule. Because twelve
columns per site is three times the guild figure's four, the body-site axis was
cut at the top 8 labels rather than the top 20; that pools 192 labels and 22.8
per cent of profiles into Other, against 15.3 per cent in the guild figure. The
condition axis keeps the top 60 labels as before, which still pools 585 labels
and 46.5 per cent of profiles, and BugSigDB records study variables such as Diet,
Age, and Extraction protocol in the same Condition field as diseases, so the
condition axis is not purely a disease axis. Cells supported by fewer than three
studies were masked. On this run the grid is 59 conditions by 108 columns; of
6,372 cells, 1,129 carry a value, 587 fall below the study threshold, 204 are
sampled with no signal for that trait, and 4,452 were never sampled. Two
conditions, Atopic asthma and Helminthiasis, had no surviving cell and were
dropped before clustering.

Two presentation choices differ from the guild figure and both were forced by the
traits themselves. The guild figure draws its two guilds on one shared colour
range, because a guild's share is directly comparable to the other guild's share.
The six Table 1 traits are not comparable in that way: their corpus-wide
prevalence spans more than an order of magnitude, and on a shared range the four
rarer traits would be uniformly white. Each trait is therefore drawn on its own
colour range, each colour bar carries its own maximum, and the figure states that
magnitudes are readable within a trait and not across traits. Separately, in each
trait's own range a single extreme cell was flattening the rest of the
distribution to one pale tint, so ramps are cut at the 98th percentile of that
trait's observed cells and cells above the cut saturate at the top colour. Cells
where a site was sampled for a condition but no taxon carried the trait remain
filled dark grey, distinct from a genuinely small fraction, and the four cell
states are recorded per cell in a companion status table.

Second figure: paired directional shift. There is no aerobe-against-anaerobe
contrast to test here, so each trait is tested against itself across the two
directions. Every study contributes one value per trait, the mean increased
fraction minus the mean decreased fraction over that study's paired experiments;
collapsing to the study first is the correction issue #133 asks for, since 60.2
per cent of studies contribute more than one experiment. A study counts as
shifted up for a trait when that difference is positive and down when it is
negative. Exact ties are recorded in the output table and excluded from the test
rather than split between the sides, which is why the denominator differs between
traits: ties are common for the rarer traits, 510 studies for arsenate reduction
against 56 for facultatively aerobic. Each trait then gets a two-sided exact
binomial test against a half-and-half null with an exact binomial 95 per cent
interval, and the six p-values are Benjamini-Hochberg corrected together as one
family.

Over 1,455 studies, five of the six traits sit on the increased side and one on
the decreased side. Facultatively aerobic shifts up in 791 of 1,399 scored
studies (share 0.565, q = 6.6e-06), moderately halophilic in 722 of 1,289 (0.560,
q = 5.3e-05), sphere shaped in 690 of 1,278 (0.540, q = 0.0071), arsenate
reduction in 507 of 945 (0.537, q = 0.027), and coccus shaped in 746 of 1,393
(0.536, q = 0.010). Cellulose catabolic process runs the other way, shifting up
in only 527 of 1,189 studies (0.443, q = 0.00020). These directions agree with
the composition arm's independent paired result on the same six traits, which
also places cellulose catabolic process alone on the decreased side. The effect
sizes are small: every trait sits within about six points of 0.5, and the
significance comes from the number of studies rather than from the size of the
shift, which is worth stating explicitly next to the q-values.

The right panel reports Spearman rank correlations between traits across those
per-study contrasts. This is deliberately not the same quantity as the guild
figure's within-profile correlation, which compares two shares of one profile
total and is therefore pushed negative by compositionality; correlating one
trait's contrast against another trait's contrast across studies does not share
that artefact. Coccus shaped and sphere shaped correlate at rho = 0.47 over 1,455
studies, which is high enough to treat them as overlapping evidence rather than
two independent lines of it, but well short of the near-identity their names
suggest. Facultatively aerobic correlates positively with coccus shaped (0.45),
arsenate reduction (0.42), and moderately halophilic (0.37), and negatively with
cellulose catabolic process (-0.23). A pair is omitted rather than written as a
missing value if either series is constant, since Spearman is undefined there.

All statistics use scipy (binomtest, false_discovery_control, spearmanr) rather
than hand-rolled implementations. Both figures and their five backing tables are
produced by src/bsdbng/table1_trait_maps.py, driven deterministically by
`just table1-trait-maps`, which shares its matrix construction, masking, row
clustering, and table writing with src/bsdbng/guild_maps.py rather than
duplicating them. The BugSigDB export timestamp is read from the dump's leading
comment line and recorded in every output table.

Caveat on vintage. This run reads the full_dump.csv currently in data/raw, an
unpinned devel-branch export dated 2026-04-18. The enrichment arm has since been
rerun against an updated BugSigDB release, so these figures will need
regenerating once the ingested dump is pinned and agreed; the numbers above are
expected to move. The command is a single recipe and the regeneration is cheap.

The pipeline code and this description were drafted with the assistance of a
coding agent (Claude Code) and reviewed by the author.
