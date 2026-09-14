# Methods

Prose descriptions of how each result in this repository was produced, written so
that the text can be used directly in the Methods section of the paper. Each file
names the script that produced the figure, the `just` recipe that runs it, and the
BugSigDB vintage the numbers came from.

Start with the shared [biological interpretation guide](biological_interpretation.md)
for measurement units, evidence provenance, defensible wording and caption/review
requirements across every analysis. The companion
[claim audit](biological_interpretation_audit.md) covers current prose and
generators, PR #173's overlays and PR #192's archived Rosalind discussion/review.
These rules apply to manuscript and public descriptions as well as methods.

| Figure | Produced by | Recipe | Methods |
| --- | --- | --- | --- |
| `guild_condition_body_site_heatmap.png` | `src/bsdbng/guild_maps.py` | `just guild-maps` | [`guild_condition_body_site_heatmap.md`](guild_condition_body_site_heatmap.md) |
| `table1_traits_condition_body_site_heatmap.png` | `src/bsdbng/table1_trait_maps.py` | `just table1-trait-maps` | [`table1_trait_condition_body_site_maps.md`](table1_trait_condition_body_site_maps.md) |
| `table1_traits_first_last_body_site_panels.png` | `src/bsdbng/table1_trait_maps.py` | `just table1-trait-maps` | [`table1_trait_condition_body_site_maps.md`](table1_trait_condition_body_site_maps.md) |
| `table1_traits_direction_shift.png` | `src/bsdbng/table1_trait_maps.py` | `just table1-trait-maps` | [`table1_trait_condition_body_site_maps.md`](table1_trait_condition_body_site_maps.md) |

Figures write to `data/derived/guild_maps/` and `data/derived/table1_traits/`
alongside the tab-separated tables backing them, so every number quoted in a
methods file can be checked without measuring pixels. Those outputs are
gitignored and are meant to be regenerated from the recipe rather than committed.
PNG export needs a headless Chrome (kaleido 1.x); run `just get-chrome` once —
it fetches a version-pinned Chrome for Testing build (~170 MB) that kaleido
then prefers over a system Chrome, which is what makes the renders reproducible.

[`three_arm_statistics.md`](three_arm_statistics.md) is not about one figure. It
records how the three analyses built on `study_trait_profiles_wide.tsv` differ,
and why their p-values are not on comparable scales.

The gram-stain directional Sankey diagrams live on the
`147-add-gram-stain-directional-sankey-gram-positive-gram-negative-x-increased-decreased`
branch together with their own methods file; they are not part of this branch.
