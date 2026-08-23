# All BugSigDB Condition/Body-Site Direction Maps

These maps summarize all BugSigDB study/experiment/direction profiles available in `trait-visualization/study_trait_profile_study_metadata.tsv`. They do not use a selected trait set or trait-enrichment score; they are direction-support maps for increased and decreased profiles.

Profiles summarized: 12,108
Studies summarized: 1,924

## Files

- `all_bugsigdb_condition_body_site_direction_summary.tsv`: all condition/body-site pairs, including broad labels and `Other`/`unknown` body-site groups.
- `all_bugsigdb_condition_body_site_direction_summary_v3_individual.tsv`: all condition/body-site pairs after v3 raw-label display mapping.
- `all_bugsigdb_condition_body_site_direction_summary_v3_grouped.tsv`: all condition/body-site pairs after v3 disease/body-site display rollups.
- `all_bugsigdb_increased_profile_condition_body_site_map.png`: displayed increased-profile map.
- `all_bugsigdb_increased_profile_condition_body_site_map_edges.tsv`: display-candidate edge table 
  for the increased-profile map; 24 rows are drawn.
- `all_bugsigdb_increased_profile_condition_body_site_map_v2.png`: uncapped v2 increased-profile map.
- `all_bugsigdb_increased_profile_condition_body_site_map_v2_edges.tsv`: uncapped v2 edge table for the increased-profile map.
- `all_bugsigdb_increased_profile_condition_body_site_map_v3_individual.png`: v3 increased-profile map with raw labels and grey comparator-background edges.
- `all_bugsigdb_increased_profile_condition_body_site_map_v3_individual_edges.tsv`: v3 individual edge table for the increased-profile map; 32 rows are drawn.
- `all_bugsigdb_increased_profile_condition_body_site_map_v3_grouped.png`: v3 increased-profile map with display rollups and grey comparator-background edges.
- `all_bugsigdb_increased_profile_condition_body_site_map_v3_grouped_edges.tsv`: v3 grouped edge table for the increased-profile map; 30 rows are drawn.
- `all_bugsigdb_increased_profile_condition_body_site_map_v3_grouped_condition_rollups.tsv`: raw-to-display condition-label audit table.
- `all_bugsigdb_increased_profile_condition_body_site_map_v3_grouped_body_site_rollups.tsv`: raw-to-display body-site audit table.
- `all_bugsigdb_decreased_profile_condition_body_site_map.png`: displayed decreased-profile map.
- `all_bugsigdb_decreased_profile_condition_body_site_map_edges.tsv`: display-candidate edge table 
  for the decreased-profile map; 23 rows are drawn.
- `all_bugsigdb_decreased_profile_condition_body_site_map_v2.png`: uncapped v2 decreased-profile map.
- `all_bugsigdb_decreased_profile_condition_body_site_map_v2_edges.tsv`: uncapped v2 edge table for the decreased-profile map.
- `all_bugsigdb_decreased_profile_condition_body_site_map_v3_individual.png`: v3 decreased-profile map with raw labels and grey comparator-background edges.
- `all_bugsigdb_decreased_profile_condition_body_site_map_v3_individual_edges.tsv`: v3 individual edge table for the decreased-profile map; 29 rows are drawn.
- `all_bugsigdb_decreased_profile_condition_body_site_map_v3_grouped.png`: v3 decreased-profile map with display rollups and grey comparator-background edges.
- `all_bugsigdb_decreased_profile_condition_body_site_map_v3_grouped_edges.tsv`: v3 grouped edge table for the decreased-profile map; 30 rows are drawn.
- `all_bugsigdb_decreased_profile_condition_body_site_map_v3_grouped_condition_rollups.tsv`: raw-to-display condition-label audit table.
- `all_bugsigdb_decreased_profile_condition_body_site_map_v3_grouped_body_site_rollups.tsv`: raw-to-display body-site audit table.

## Interpretation

Each row of the source profile metadata is one BugSigDB study/experiment/direction profile. Edges connect condition labels to grouped body sites. Edge width is the number of profiles in the focus direction. Edge color is the direction balance: `(focus profiles - comparator profiles) / all profiles` for that condition/body-site pair. Node labels use `n=` for focus-profile count, `s=` for unique focus-direction study count, and `focus=` for the focus-direction fraction among profiles represented by displayed edges.

The TSV preserves all sampled condition/body-site pairs. The PNG is intentionally filtered for readability: broad metadata-like condition labels and `Other`/`unknown` body-site groups are kept in the TSV but not drawn.

The `_v2` PNGs remove the readability caps on displayed edges, conditions, body-site groups, and per-body-site edge selection. They still use the same support threshold and label filters as the displayed maps.

The `_v3_individual` PNGs keep raw condition labels and raw body-site labels. The `_v3_grouped` PNGs use the same v3 display rollups as the guild Sankey maps. Disease-like condition labels are grouped when the display rule matches, but transplant/graft labels and metadata/exposure labels are kept raw. Raw body sites are rolled up to readable display groups such as `Upper airway / nasal`, `Oral cavity / dental`, `Lower airway / lung`, `Vaginal / reproductive`, and `Other body sites`.

In each v3 direction map, colored edges count the focus-direction BugSigDB profiles and parallel grey edges count comparator-direction profiles in the same condition/body-site pair for that display mode. V3 condition nodes are ordered by the colored focus-profile count reported as `n=`; `bg=` reports comparator-profile support.

The v3 individual summary contains 956 condition/body-site pairs before PNG selection. The v3 grouped summary contains 753 display condition/body-site pairs before PNG selection.
