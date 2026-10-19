# Parkinson's Disease Cellulose Figures

This folder contains the combined ecological-reporter figure, the globally FDR-tiered
condition-taxon-trait evidence landscape, vector versions, statistical methods, and
machine-readable audit tables. Panels A and B of the combined figure retain only conditions with
a descriptive 95% study-level interval excluding zero for at least one of the six individual
traits; composite scores do not qualify a condition. The evidence landscape shows all eligible
tests while separating FDR-significant results (`q <= 0.05`), exploratory results
(`0.05 < q <= 0.20`), and the remaining tested background. The node and edge tables remain
significant-only; any KG-Microbe annotation links are explicitly noninferential provenance.

See [the interpretation guide](pd_ecological_reporter_figures.md) for definitions of
**ecological reporter** and **ecological-reporter composite**, a panel-by-panel explanation, and
details of the second panel's axes.

## Outputs not yet relabelled (#555)

Profiles now take their own experiment's condition and body site. The `pd_cellulose_*`,
`pd_six_trait_direction_*`, `pd_condition_taxon_trait_network*` and
`pd_ecological_reporter_combined*` outputs were regenerated with those labels;
`just pd-ecological-reporter-combined` builds the reporter tables from the profile inputs (#585).
The following still reflect the earlier study-majority labels:

- `pd_trait_taxa_lookup.md`, `pd_six_trait_taxon_carriers.tsv`,
  `pd_six_trait_taxon_profile_evidence.tsv` and `pd_four_positive_trait_shared_taxa.tsv`, which
  count 253 PD fecal profiles instead of the current 247. `just pd-trait-taxa` stops by design
  because coccus shaped has moved to the negative side of the six-trait lollipop, leaving three
  positive-side traits instead of four.
- The `pd_species_six_traits_*`, `pd_exploratory_taxa_six_traits_*` and `pd_trait_contributors*`
  outputs, which read that evidence table.
