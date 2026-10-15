# Oral experiment sample sizes

[Histogram PNG](oral_experiment_sample_sizes.png) ·
[Editable SVG](oral_experiment_sample_sizes.svg) ·
[Vector PDF](oral_experiment_sample_sizes.pdf) ·
[Experiment audit](oral_experiment_sample_size_audit.tsv) ·
[Summary](oral_experiment_sample_size_summary.tsv)

The red line is an external reference study with **200 samples: 100 disease cases
and 100 healthy controls**. It is not added to the observed histogram.

The saved oral heatmap cohort contains 425 experiments from 137 studies. Both
group sizes are usable for 408 experiments from 136 studies; 17 experiments lack
one or both counts. Among the 408 available totals, the median is 54, the middle
50% span 30–99.25, and the range is 6–1,206. There are 354 totals below 200 (86.8%),
one equal to 200, and 53 above 200.

The left panel uses 25-sample bins through n=500; the right panel shows all 408
totals with logarithmic bins and a logarithmic x-axis. Twenty totals exceed 500.
Both panels describe the same cohort, not separate datasets.

## Cohort and denominator

Membership starts with the saved individual-condition oral heatmap rows, using
the existing oral-site OR exact-oral-diagnosis selection. It then joins the saved
profile metadata by condition and the source workflow's grouped body site. All
80 oral condition/site rows are represented; the cohort is not restricted to
significant cells or the rows retained in a particular display cut. In this
snapshot every selected experiment is sampled orally. Increased/decreased
profiles and repeated signature records are collapsed to one experiment.

This is the oral **heatmap input cohort**, which can include experiments with
only one signature direction. It is not restricted to Jed's separate matched
increased/decreased paired analysis. Signature directions describe taxonomic
changes; they are not the two participant groups.

Sample sizes come from the local pinned BugSigDBExports v1.3.1 `full_dump.csv`
(snapshot 2026-04-24_00:41_UTC), not the newer live `experiments.csv` export.
`BSDB ID` supplies the study and experiment identifiers using the ingestion
convention: `bsdb:study/experiment/signature` becomes `bsdb:study-experiment`.
The plotted n is **Group 0 sample size + Group 1 sample size**. Each source count,
group label, study design, inclusion status and experiment ID is retained in the
audit. [Input checksums](inputs.tsv) identify the three exact source files.

Missing, nonpositive, nonintegral, nonfinite or conflicting counts are excluded
and audited, never converted to zero or silently selected from conflicting
signature rows. All 425 selected experiment IDs matched the pinned dump; none
had invalid or conflicting counts. Sixteen lacked both group sizes and one
lacked group 1 size.

Experiment totals are not necessarily counts of unique people: longitudinal
comparisons can reuse participants across groups, and experiments from the same
study can share controls or participants. Existing comparisons are not all
disease-versus-healthy contrasts. The histogram compares reported sizes, not
statistical power or effective independent sample size. Study-level totals are
not summed because participant overlap cannot be resolved from these fields.

## Regenerate

With the pinned dump available locally:

```bash
just sync
just oral-experiment-sizes
```

An existing dump elsewhere can be supplied without downloading another copy:

```bash
just oral-experiment-sizes /path/to/full_dump.csv
```

The recipe writes the figure, per-experiment audit, summary and input checksums.
It does not refit heatmap effects or tests. Tracked in issue #271.
