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
| Matched-corpus Table 1-trait analog source contributions and genus/family influence TSVs | `src/mtp/taxon_contributions.py` | `just table1-taxon-influence` | [`table1_taxon_influence.md`](table1_taxon_influence.md); [measured results and cross-scale interpretation flags](table1_taxon_influence_results.md) |
| `trait_direction_screen.tsv` (no figure) | `src/bsdbng/trait_screen.py` | `just trait-direction-screen` | [`trait_direction_screen.md`](trait_direction_screen.md) |
| `trait_coupling_loto_matrix.tsv`, `trait_coupling_lco_summary.tsv`, `trait_coupling_report.txt` (no figure) | `src/mtp/coupling.py` | `just audit-trait-coupling` | [`trait_coupling_sensitivity.md`](trait_coupling_sensitivity.md) |

Figures write to `data/derived/guild_maps/` and `data/derived/table1_traits/`
alongside the tab-separated tables backing them, so every number quoted in a
methods file can be checked without measuring pixels. Those outputs are
gitignored and are meant to be regenerated from the recipe rather than committed.
PNG export needs a headless Chrome (kaleido 1.x); run `just get-chrome` once —
it fetches a version-pinned Chrome for Testing build (~170 MB) that kaleido
then prefers over a system Chrome, which is what makes the renders reproducible.

[`integrated_evidence.md`](integrated_evidence.md) documents
`just integrated-evidence`: the shared result contract in which every method keeps its
own row, family, units and states, with nothing pooled across methods and every absent
input recorded as a state rather than a blank.

[`three_arm_statistics.md`](three_arm_statistics.md) is not about one figure. It
records how the three analyses built on `study_trait_profiles_wide.tsv` differ,
and defines the shared corpus, paired quantities and evidence levels required
before choosing comparison thresholds or ranking scores. It also explains why
the native p-values are not on comparable scales.

[`contrast_orientation.md`](contrast_orientation.md) documents
`just contrast-orientation`: what each BugSigDB experiment compares, which way its
direction points, the disease-versus-healthy population with its review, and the
six Table 1 directions and 15 ratios on that population beside the historical corpus.

[`table1_paired_log_ratios.md`](table1_paired_log_ratios.md) documents
`just table1-compositional-stats`: paired log-ratio direction tests for all 15
Table 1-trait pairs as a matched-corpus analog, with explicit references, zero exclusions, study-level
coverage, and a separate correction family. It emits statistics and a Markdown
report without rendering figures.

[`calibration.md`](calibration.md) documents `just calibration-report`: what the
directional tests do when their own nulls are true, their power against a declared
injected effect, and what selecting traits on the data costs, against criteria fixed
in [`calibration_protocol.json`](calibration_protocol.json). That protocol is on its
second version; every amendment, what prompted it and what was already known when it
was made are recorded in its `amendments` array rather than edited in quietly.
`just calibration-flip-assignments` writes only the shared flip-assignment table the
other two arms would need to be measured on the same null (#263).
[`calibration_screen.md`](calibration_screen.md) documents
`just calibration-screen`: the 317-trait directional screen run under that same
study-flip null, counting the traits each correction calls when nothing is there.

[`calibration_composition.md`](calibration_composition.md) documents
`just calibrate-composition`: the paired composition arm calibrated against the
shared study-flip complete null (#263), measuring empirical Type I error control,
family-wise error rates, and Table 1 canonical trait calibration.

[`table1_log_ratio_validation.md`](table1_log_ratio_validation.md) records the
pinned-corpus rerun of all 15 ratios alongside the original six fraction
directions as a matched-corpus Table 1-trait analog sensitivity analysis (not a
reproduction of the published Table 1 baseline), including study coverage,
uncertainty and interpretation limits.

[`table1_log_ratio_sensitivity.md`](table1_log_ratio_sensitivity.md) documents
`just table1-log-ratio-sensitivity`: predeclared scenarios separating population,
estimand, reference, stratum and single-study effects on the 15 ratios, with the
protocol in [`table1_log_ratio_sensitivity_protocol.json`](table1_log_ratio_sensitivity_protocol.json).

[`cellulose_decr.md`](cellulose_decr.md) documents `just cellulose-decr` and
`just cellulose-decr-figures`, which build `data/analysis/cellulose_decr/`. That
analysis covers feces-decreased taxa annotated with cellulose catabolic process
in Parkinson's disease, HIV infection and colorectal cancer:

- their overlap across the three conditions, tested against six null models
  (U1–U6) and compared with a within-condition reference;
- the direction of the cellulose annotation shift at experiment and study level,
  under exact and synonym condition labels.

Its tables are committed and its figures are regenerated. The code is the
`mtp.cellulose_decr` package (`src/mtp/cellulose_decr/`).

The selected-figure slide deck in `data/derived/figure_review_2026-09-21/` is built
by `scripts/build_slide_deck.py` (`just figure-review-deck`), which places the
methods and Tufte plates drawn by `scripts/generate_tufte_plates.py`
(`just figure-review-plates`). `just figure-review-deck` also validates the figure
provenance YAML against its LinkML schema and prints the PDF with headless Chrome
(`just html-to-pdf`). The deck, its assets, the plates and the provenance YAML are
gitignored local outputs, and the builder reads gitignored inputs as well (the
earlier deck manifests and figures under `data/derived/`), so a fresh clone cannot
rebuild the deck (#573). Each slide shows the SVG its generator drew beside the
figure's PNG, or else the PNG; the September decks' `svg_sources/` copies are not used,
because no rerun redraws them (#661). Both builders write to the checkout that holds
them, whatever `--input-root` points at, and the deck reads the plates and the
cellulose_decr figures from that checkout too; everything else comes from the input root.
They refuse to replace an existing deck or plates without `--overwrite`; there is no other
copy, so back up `slides/`, `provenance/` and `data/derived/tufte_plates/` first. The
Table 1 numbers the cover, the plates and their slides quote are read from
`data/derived/contrast_orientation/disease_control_comparison.tsv` and the tables behind
each figure, so they follow the contrast classifier of the rerun (#579, #606). So are
Arm 2's tests and calls (the ORA tables' headers), the oxygen-guild McNemar test
(`guild_aerobe_anaerobe_correlation.tsv`), the screened sign test's rejection rate
(`calibration_criteria.tsv`, refused below its replicate budget), Plate V's PD split, Plate
II's stratum sizes and Plate VII's counts (#680, #683, #704). The fraction sign tests are
quoted with the studies left after the zero-total exclusion, the log-ratios with each
population's studies, since their eligibility is decided pair by pair.

The slides are stamped with the rerun revision (`COMMIT` in the builder, `RERUN_COMMIT` in
the plate script, which a test keeps equal) and link their committed files at it. Before
writing anything, `just figure-review-deck` checks each linked figure and evidence table,
each committed table the slide text quotes, each plate's recorded inputs and the
cellulose_decr tables against the revision their slides cite, using `git ls-tree` listings
the recipe writes. It refuses when a file differs. `--allow-unverified` builds a review
deck instead: the files that differ are not linked, their slides say so, and the cover and
every stamp read "review build" (#637). Slides whose code was never merged into main
(PRs #148, #173 and #315) stay in the deck: they link that code at the pull request's head
commit, and their stamp reads "not rerun on main" (#614). PR #173 also committed its Table 1
maps, whose significance outlines and composition triangles main's `just table1-trait-maps`
does not draw; the builder reads those five maps from that pull request's head and checks
them there, so no checkout needs them on disk (#677, #704). The build also stops, or with
`--allow-unverified` names the reason as a review build, when the plates' record lacks a
count the slides quote (#682) or the cellulose_decr tables were built on another contrast
classifier than the Table 1 tables (#679).
Gitignored figures cannot be checked against a revision; the provenance records their
SHA-256. The check reads the local repository, so push the commits that drew the plates and
cellulose figures before sharing a deck: the slides link them on GitHub. The provenance folder's `figure_provenance_and_results.json`,
`figures_manifest.tsv` and `figure_statistical_results.tsv` have no generator and are not
updated by a rebuild (#639).

To build into other folders, set them as
`just deck_slides=DIR deck_provenance=DIR figure-review-deck`, so the validation, PDF and
copy steps use them too; the recipe ends by writing `deck_build_status.json`.

The [GPT Rosalind discussion and review](../docs/rosalind-six-trait-interpretation.md)
archives proposed biological explanations for the six traits, with a separate
evidence review and assigned follow-up issues. It records hypotheses to evaluate;
it is not an approved Results or Discussion section.
The [second exchange](../docs/rosalind-six-trait-interpretation-response-2.md)
preserves the follow-up questions and Rosalind's revised response, including its
prediction and discriminating-evidence tables. It is archived without a separate
evidence review.

[`validation_feasibility.md`](validation_feasibility.md) documents
`just validation-evidence`: what measured inflammation, SCFA, redox, exposure and
covariate evidence is reachable, which BugSigDB experiments are stratified by a
measured marker, what the SCFA-related annotations actually assert, which external
cohorts overlap the discovery corpus, and the protocol a real validation would need.

[`biological_interpretation.md`](biological_interpretation.md) specifies project-wide
biological interpretation standards, distinguishing organism trait annotation from in vivo
metabolic flux, formalizing the 4 meanings of direction, providing defensible wording
rules and caption contracts, and illustrating them with worked examples.

[`biological_interpretation_audit.md`](biological_interpretation_audit.md) audits all
claims across methods, protocols, and documentation against these standards (findings
C01–C11, P01–P10, R01–R09, G01–G13), recording resolution status for each finding.

[`table1_trait_semantics_audit.md`](../docs/table1_trait_semantics_audit.md) documents
the semantic scope, source ontology mappings, approved vocabulary, and unverified
inferences for each of the six Table 1 traits.

[`coassertion_audit.md`](coassertion_audit.md) defines the scope-specific KG
node and strain co-assertion census, the generated count-bearing methods table,
and the proposed oxygen-guild carrier-exclusion sensitivity. It distinguishes
annotation unions from physiological contradiction; final counts and the
resolution decision remain pending under #846.

The gram-stain directional Sankey diagrams live on the
`147-add-gram-stain-directional-sankey-gram-positive-gram-negative-x-increased-decreased`
branch together with their own methods file; they are not part of this branch.
