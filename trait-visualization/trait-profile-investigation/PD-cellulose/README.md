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

## Label sources (#555, #584)

Profiles now take their own experiment's condition and body site. These outputs were regenerated
with those labels:

- `pd_cellulose_*`, `pd_six_trait_direction_*` and `pd_condition_taxon_trait_network*` (#555).
- The PD trait-taxa outputs (#584): `pd_trait_taxa_lookup.md`, `pd_six_trait_taxon_*` and
  `pd_six_trait_taxon_evidence_provenance.json` from `just pd-trait-taxa`; the
  `pd_species_six_traits_*` and `pd_exploratory_taxa_six_traits_*` heatmaps and reference tables
  from `just pd-trait-taxa-heatmap`; and the `pd_trait_contributors*` tables and figures from
  `just pd-trait-contributors`.

Only the `pd_ecological_reporter_combined_*` source tables still reflect the earlier
study-majority labels; no script in this repository regenerates them (#585).

## PD trait-taxa outputs (#584)

`just pd-trait-taxa` traces all six lollipop traits whatever their sign, over the 247 PD fecal
profiles (128 increased and 119 decreased). Each trait's observed direction, the sign of
`increased_minus_decreased_mean_z` in `pd_six_trait_direction_summary.tsv`, travels with it as
the `trait_direction` column and as up and down arrows in the figures. Moderately halophilic,
arsenate reduction and facultatively aerobic are higher in increased profiles; coccus shaped,
sphere shaped and cellulose catabolic process are higher in decreased profiles.

`pd_six_trait_taxon_evidence_provenance.json` records the PD fecal profile set the evidence was
built from (count and sha256 of the sorted profile ids) and each trait's mean-z difference.
`just pd-trait-taxa-heatmap` and `just pd-trait-contributors` stop before writing when the
current profiles, traits, directions or lollipop values differ from the recorded ones.

The contributor tables come from the Table 1 contribution ledger (#190) on the pinned KG in
`data/kg-microbe`: 198 PD fecal profiles from human experiments with both directions. The same
generator reproduces the contributor counts, shares and selected taxa that were committed before
it existed; only the PD-versus-other-feces specificity columns change, to 247 PD and 3,357
other fecal profiles.
