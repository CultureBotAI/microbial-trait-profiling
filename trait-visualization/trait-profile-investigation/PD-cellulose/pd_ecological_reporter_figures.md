# PD Ecological-Reporter Specificity and Provenance

![Combined figure](pd_ecological_reporter_combined.png)

The combined figure separates three evidence levels. Panel A uses paired, study-level contrasts of annotation-derived profile z-scores. Panel B places PD among all fecal condition labels with at least three paired studies (56 in this dataset). Panel C intersects exact reported PD genus/species labels with the local KG-Microbe taxa-trait matrix and colors annotations by BugSigDB direction support.

For PD, the study-weighted cellulose increased-minus-decreased contrast is `-0.651` across `45` paired studies. The exploratory reporter composite is `+1.150`. The composite equals the mean contrasts for moderately halophilic, facultatively aerobic, and arsenate reduction minus the cellulose contrast; it is not a validated diagnostic score.

[Condition/trait audit table](pd_ecological_reporter_combined_condition_trait_summary.tsv) | [Study contrasts](pd_ecological_reporter_combined_study_contrasts.tsv) | [Taxon-panel cells](pd_ecological_reporter_combined_taxon_panel.tsv) | [Vector figure](pd_ecological_reporter_combined.svg)

## Condition-Taxon-Trait Network

![Network](pd_condition_taxon_trait_network.png)

The network contains `31` nodes and `44` edges. BugSigDB provides the condition, signature direction, exact taxon report, and study support. KG-Microbe provides the taxon-to-trait annotations. Genus and species are shown together as separate exact IDs but are never summed. Taxon-to-trait links repeat the taxon's evidence support for each annotation and therefore are not conserved flows.

[Network nodes](pd_condition_taxon_trait_network_nodes.tsv) | [Network edges](pd_condition_taxon_trait_network_edges.tsv) | [Taxon selection](pd_condition_taxon_trait_network_taxon_selection.tsv) | [Vector figure](pd_condition_taxon_trait_network.svg)

## Interpretation Boundary

These are ecological annotations of taxa found in curated differential-abundance signatures. They are not measurements of gene expression, metabolite flux, cellulose degradation, or longitudinal host physiology. The condition comparison is descriptive and is not an external diagnostic validation.
