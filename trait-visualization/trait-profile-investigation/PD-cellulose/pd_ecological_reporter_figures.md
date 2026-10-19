# PD Ecological-Reporter Pattern Across Fecal Conditions, and Provenance

## Terminology

An **ecological reporter** is a microbial community or trait pattern that reflects the
environment in which the microbes live. Here, the phrase means that gut taxa and their
KG-Microbe trait annotations may report changes associated with host state, intestinal
transit, diet, medication, inflammation, oxygen/redox conditions, or other ecological
pressures. It does not mean that the organisms are literal sensors, and it does not establish
which host change produced the pattern.

For each trait and study, the analysis first calculates:

```text
trait contrast = mean trait z in increased profiles - mean trait z in decreased profiles
```

The exploratory **ecological-reporter composite** is:

```text
mean(
  moderately halophilic contrast,
  facultatively aerobic contrast,
  arsenate reduction contrast
) - cellulose catabolic process contrast
```

A positive composite means the mean contrast for the three selected traits exceeds the
cellulose catabolic process contrast. It does not require opposite signs: contrasts of `+2`
and `+1`, respectively, would yield a positive composite of `+1`. The proposed two-sided pattern
requires a positive three-trait mean and a negative cellulose contrast; Panel B's upper-right
quadrant identifies that combination. The composite is a descriptive, study-level summary of
trait annotations, not a validated diagnostic score or a patient-level measurement.

## Combined Figure

![Combined figure](pd_ecological_reporter_combined.png)

The combined figure and every number in this section come from
`pd_ecological_reporter_combined_*.tsv`, which `just pd-ecological-reporter-combined` builds from
the profile inputs (#585). Each profile takes its own experiment's condition and body site, and the
trait z-scores are the six-trait z-scores of the other PD-cellulose outputs, so for the conditions
both cover, the study contrasts here equal those of the condition screen in the evidence landscape
below.

### Panel A: paired trait contrasts

Rows are the fecal condition labels that have at least three paired studies and a descriptive
95% study-level t interval excluding zero for at least one of the six individual traits. Neither
the three-trait score nor the ecological-reporter composite can qualify a condition. Columns are
the six shared traits plus the ecological-reporter composite. Each cell is the mean of paired
within-study contrasts, so studies are weighted equally rather than allowing studies with more
profiles to dominate. Red means the trait is stronger in increased profiles; blue means it is
stronger in decreased profiles. `n=` is the number of studies containing both directions. A dark
border marks the individual cells whose descriptive 95% t interval excludes zero. The orange
outline identifies the PD row.

### Panel B: PD among the qualifying conditions

Each point is one of the same qualifying fecal condition labels shown in Panel A. Conditions
with no interval-excluding-zero change among the six individual traits are not plotted. Point size
increases with paired-study count, and PD is highlighted in orange.

The x-axis is the cellulose contrast with its sign reversed. Moving right means cellulose
catabolic process is relatively stronger in decreased profiles.

The y-axis, **three-trait increased-side composite**, means:

```text
mean(
  moderately halophilic increased-minus-decreased z,
  facultatively aerobic increased-minus-decreased z,
  arsenate reduction increased-minus-decreased z
)
```

Moving upward means those three traits are, on average, stronger in increased profiles. The
upper-right quadrant therefore matches both components of the proposed ecological-reporter
pattern. PD being in that quadrant shows pattern agreement, not PD specificity: the numbers below
count the other conditions in the same region and those with larger composites.

Those three traits were chosen because, in the PD paired-study analysis, they were the three
individual traits with positive increased-minus-decreased contrasts whose descriptive 95%
intervals excluded zero: moderately halophilic, facultatively aerobic, and arsenate reduction.
Cellulose catabolic process supplied the opposite, decreased-profile component because its PD
contrast was negative and its interval also excluded zero. Coccus shaped and sphere shaped were
not included in the y-axis composite because their PD intervals crossed zero. The resulting axis
is therefore a data-derived exploratory summary of the PD pattern, not an independently validated
or prospectively specified score.

### Panel C: PD taxon and trait provenance

Rows are exact genus or species labels reported in PD fecal BugSigDB signatures; columns are the
six shared traits. A circle means the local KG-Microbe matrix annotates that taxon with the trait.
Circle size is the number of increased plus decreased study reports, so a study that reports the
taxon in both directions counts twice; the same count ranks the rows. Circle color shows the log2
ratio of increased to
decreased BugSigDB study support, with red favoring increased and blue favoring decreased. An `x`
means no local annotation for that taxon-trait cell. Taxonomic ranks are kept separate and are not
summed.

### Numbers

`just pd-ecological-reporter-combined` writes this block from the tables below; do not edit it by hand.

<!-- BEGIN numbers written by `just pd-ecological-reporter-combined` -->
- Qualifying fecal conditions in panels A and B: `16`.
- In panel B's upper-right quadrant: `13` of the `16`, PD among them. Of the other `12`, `9` have a larger ecological-reporter composite than PD and `10` a larger three-trait composite.
- PD, across `47` paired studies: cellulose contrast `-0.677` (descriptive 95% t interval `-1.139` to `-0.214`); three-trait increased-side composite `+0.459` (`+0.194` to `+0.724`); ecological-reporter composite `+1.136` (`+0.538` to `+1.734`).
<!-- END numbers written by `just pd-ecological-reporter-combined` -->

The ecological-reporter composite equals the mean contrasts for moderately halophilic, facultatively aerobic, and arsenate reduction minus the cellulose contrast; it is not a validated diagnostic score.

[Condition selection](pd_ecological_reporter_combined_condition_selection.tsv) | [Condition/trait audit table](pd_ecological_reporter_combined_condition_trait_summary.tsv) | [Study contrasts](pd_ecological_reporter_combined_study_contrasts.tsv) | [Taxon-panel cells](pd_ecological_reporter_combined_taxon_panel.tsv) | [Vector figure](pd_ecological_reporter_combined.svg)

## Condition-Taxon-Trait Evidence Landscape

![Evidence landscape](pd_condition_taxon_trait_network.png)

The regenerated figure separates inference from annotation provenance. It screened every human
fecal condition with at least three studies containing both increased and decreased BugSigDB
signatures. Exact genus and species IDs were kept separate, required at least three reporting
studies, and were tested with an exact McNemar/binomial test on discordant paired studies. The
complete eligible condition-taxon family was corrected globally with Benjamini-Hochberg.

The evidence landscape displays all eligible tests without weakening the significance standard:

- Filled points pass global `q <= 0.05` and are called FDR-significant.
- Hollow amber points have `0.05 < q <= 0.20`. They are exploratory follow-up hypotheses, not
  significant discoveries.
- Gray points have `q > 0.20` and provide the full tested background.

Panel A shows condition-taxon tests. The x-axis is the fraction of paired studies favoring
increased rather than decreased reporting, while the y-axis is `-log10(global q)`. Circles are
genera, triangles are species, and point area represents the number of reporting studies. No
condition-taxon result passes `q <= 0.05`.

Panel B shows all condition-trait tests. Its x-axis is the mean study-level
increased-minus-decreased trait contrast, its y-axis is `-log10(global q)`, and point area is the
number of paired studies. No condition-trait result passes `q <= 0.05`; the smallest global q is
Parkinson's disease to moderately halophilic (`mean contrast = +0.497`, `q = 0.0795`, `n = 47`).

The exploratory tier contains one taxon result—Parkinson's disease to *Roseburia* (`q = 0.116`)—
and five trait results: Parkinson's disease to moderately halophilic (`q = 0.0795`), Parkinson's
disease to arsenate reduction, COVID-19 to cellulose catabolic process, and colorectal cancer to
cellulose catabolic process (each `q = 0.0893`), and colorectal cancer to moderately halophilic
(`q = 0.0942`). These are shown to support follow-up planning, not as significant findings.

Panel C shows all six Parkinson's disease trait effects with their raw `p`, globally adjusted `q`,
and study count. This provides PD-specific context while retaining the correction across the full
condition-trait family.

The node and edge audit tables remain significant-only. If taxon results pass in a future data
release, KG-Microbe taxon-to-trait annotations will remain explicitly noninferential because they
do not have p-values and do not measure microbial activity.

[Statistical methods and result counts](pd_condition_taxon_trait_network_methods.md) |
[All condition-taxon tests](pd_condition_taxon_trait_network_condition_taxon_tests.tsv) |
[All condition-trait tests](pd_condition_taxon_trait_network_condition_trait_tests.tsv) |
[Condition eligibility](pd_condition_taxon_trait_network_condition_eligibility.tsv) |
[Study trait contrasts](pd_condition_taxon_trait_network_study_trait_contrasts.tsv) |
[Network nodes](pd_condition_taxon_trait_network_nodes.tsv) |
[Network edges](pd_condition_taxon_trait_network_edges.tsv) |
[Taxon selection](pd_condition_taxon_trait_network_taxon_selection.tsv) |
[Vector figure](pd_condition_taxon_trait_network.svg)

## Interpretation Boundary

These are ecological annotations of taxa found in curated differential-abundance signatures. They
are not measurements of gene expression, metabolite flux, cellulose degradation, or longitudinal
host physiology. The condition comparison is descriptive and is not an external diagnostic
validation. The figures can support a hypothesis that microbial traits act as host-state readouts;
they cannot identify the underlying host signal or establish causality.
