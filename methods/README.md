# Methods

Prose descriptions of how each result in this repository was produced, written so
that the text can be used directly in the Methods section of the paper. Each file
names the script that produced the figure, the `just` recipe that runs it, and the
BugSigDB vintage the numbers came from.

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
and defines the shared corpus, paired quantities and evidence levels required
before choosing comparison thresholds or ranking scores. It also explains why
the native p-values are not on comparable scales.

[`table1_paired_log_ratios.md`](table1_paired_log_ratios.md) documents
`just table1-compositional-stats`: paired log-ratio direction tests for all 15
Table 1 trait pairs, with explicit references, zero exclusions, study-level
coverage, and a separate correction family. It emits statistics and a Markdown
report without rendering figures.

[`table1_log_ratio_validation.md`](table1_log_ratio_validation.md) records the
pinned-corpus rerun of all 15 ratios alongside the original six fraction
directions, including study coverage, uncertainty and interpretation limits.

The [GPT Rosalind discussion and review](../docs/rosalind-six-trait-interpretation.md)
archives proposed biological explanations for the six traits, with a separate
evidence review and assigned follow-up issues. It records hypotheses to evaluate;
it is not an approved Results or Discussion section.

[`validation_feasibility.md`](validation_feasibility.md) documents
`just validation-evidence`: what measured inflammation, SCFA, redox, exposure and
covariate evidence is reachable, which BugSigDB experiments are stratified by a
measured marker, what the SCFA-related annotations actually assert, which external
cohorts overlap the discovery corpus, and the protocol a real validation would need.
The [second exchange](../docs/rosalind-six-trait-interpretation-response-2.md)
preserves the follow-up questions and Rosalind's revised response, including its
prediction and discriminating-evidence tables. It is archived without a separate
evidence review.

The gram-stain directional Sankey diagrams live on the
`147-add-gram-stain-directional-sankey-gram-positive-gram-negative-x-increased-decreased`
branch together with their own methods file; they are not part of this branch.
