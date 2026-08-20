# All BugSigDB Condition/Body-Site Direction Maps

These maps summarize all BugSigDB study/experiment/direction profiles available in `trait-visualization/study_trait_profile_study_metadata.tsv`. They do not use a selected trait set or trait-enrichment score; they are direction-support maps for increased and decreased profiles.

Profiles summarized: 9,344
Studies summarized: 1,646

## Files

- `all_bugsigdb_condition_body_site_direction_summary.tsv`: all condition/body-site pairs, including broad labels and `Other`/`unknown` body-site groups.
- `all_bugsigdb_increased_profile_condition_body_site_map.png`: displayed increased-profile map.
- `all_bugsigdb_increased_profile_condition_body_site_map_edges.tsv`: display-candidate edge table 
  for the increased-profile map; 24 rows are drawn.
- `all_bugsigdb_decreased_profile_condition_body_site_map.png`: displayed decreased-profile map.
- `all_bugsigdb_decreased_profile_condition_body_site_map_edges.tsv`: display-candidate edge table 
  for the decreased-profile map; 24 rows are drawn.

## Interpretation

Each row of the source profile metadata is one BugSigDB study/experiment/direction profile. Edges connect condition labels to grouped body sites. Edge width is the number of profiles in the focus direction. Edge color is the direction balance: `(focus profiles - comparator profiles) / all profiles` for that condition/body-site pair. Node labels use `n=` for focus-profile count, `s=` for unique focus-direction study count, and `focus=` for the focus-direction fraction among profiles represented by displayed edges.

The TSV preserves all sampled condition/body-site pairs. The PNG is intentionally filtered for readability: broad metadata-like condition labels and `Other`/`unknown` body-site groups are kept in the TSV but not drawn.
