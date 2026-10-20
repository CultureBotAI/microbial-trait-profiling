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

`just pd-cellulose-analysis` writes the `pd_cellulose_*` and `pd_six_trait_direction_*` outputs
and [interpretations.md](interpretations.md). It reads the individual-condition Sankey edge
tables from `just six-shared-trait-sankeys` and stops when they no longer match the profile
inputs. `just pd-ecological-reporter-combined` writes the `pd_ecological_reporter_combined_*`
tables and figure, and the numbers block of the interpretation guide.

## Label sources (#555, #584, #585)

Profiles now take their own experiment's condition and body site. All outputs in this folder
were regenerated with those labels:

- `pd_cellulose_*`, `pd_six_trait_direction_*` and `pd_condition_taxon_trait_network*` (#555).
- `pd_ecological_reporter_combined*` from `just pd-ecological-reporter-combined`, which builds
  the reporter tables from the profile inputs (#585).
- The PD trait-taxa outputs (#584): `pd_trait_taxa_lookup.md`, `pd_six_trait_taxon_*` and
  `pd_six_trait_taxon_evidence_provenance.json` from `just pd-trait-taxa`; the
  `pd_species_six_traits_*` and `pd_exploratory_taxa_six_traits_*` heatmaps and reference tables
  from `just pd-trait-taxa-heatmap`; and the `pd_trait_contributors*` tables and figures from
  `just pd-trait-contributors`.

## PD trait-taxa outputs (#584)

`just pd-trait-taxa` traces all six lollipop traits whatever their sign, over the 247 PD fecal
profiles (128 increased and 119 decreased). Each trait's observed direction, the sign of
`increased_minus_decreased_mean_z` in `pd_six_trait_direction_summary.tsv`, travels with it as
the `trait_direction` column and as up and down arrows in the figures. Moderately halophilic,
arsenate reduction and facultatively aerobic are higher in increased profiles; coccus shaped,
sphere shaped and cellulose catabolic process are higher in decreased profiles. In the influence
heatmaps an asterisk marks a column whose paired baseline Δz has the opposite sign to the
lollipop (coccus shaped and sphere shaped).

`pd_six_trait_taxon_evidence_provenance.json` records the PD fecal profile set the evidence was
built from (count and sha256 of the sorted profile ids), each trait's mean-z difference, the
sha256 of every file it read (`bugsigdb_taxa.tsv`, `taxa_trait_matrix.tsv`,
`study_trait_profiles_wide.tsv`, `study_trait_profile_study_metadata.tsv` and
`pd_six_trait_direction_summary.tsv`), and the sha256 of the evidence and carrier tables it
wrote. `just pd-trait-taxa` hashes its inputs before reading them and stops if any of them changes
while it runs. `just pd-trait-taxa-heatmap` and `just pd-trait-contributors` stop before writing
when the current profiles, traits, directions, lollipop values, input files or evidence tables
differ from the recorded ones. A trait that no PD-reported taxon carries has no rows, so only the
tables' hashes show that rows were lost. `just pd-trait-contributors` also stops when the species
heatmap's trait table no longer matches the lollipop.

After a `data/derived` change, run these steps in order:

1. `just trait-heatmaps`, which rewrites `study_trait_profile_study_metadata.tsv`. The PD scripts
   stop at the denominator-cache check, which names this command, when the metadata's totals no
   longer match the profile matrix.
2. `just pd-cellulose-analysis`, which rebuilds `pd_six_trait_direction_summary.tsv`.
   `just pd-trait-taxa` stops when the saved summary no longer matches the current profiles, and
   names this command.
3. `just pd-trait-taxa`.
4. `just pd-trait-taxa-heatmap`.
5. `just pd-trait-contributors`.

The contributor tables come from the Table 1 contribution ledger (#190) on the pinned KG in
`data/kg-microbe` (check it with `just verify-pinned-kg`, and refetch a changed file with
`just download-kg-force`, since `just download-kg` compares sizes only): 198 PD fecal profiles
from human experiments with both directions. The same
generator reproduces the contributor counts, shares and selected taxa that were committed before
it existed; only the PD-versus-other-feces specificity columns change, to 247 PD and 3,357
other fecal profiles.
